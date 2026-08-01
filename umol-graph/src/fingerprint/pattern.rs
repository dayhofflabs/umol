//! RDKit `PatternFingerprint` replica (bit-exact, pinned RDKit 2026.03.3): a fixed
//! library of small SMARTS templates matched by subgraph isomorphism, each match
//! keyed by the frozen 32-bit boost hash and folded to a fixed-width [`BitFp`].
//!
//! Per match RDKit sets a *feature* bit (pattern index, then each matched atom's
//! atomic number in query order, then each pattern bond's type in pattern-edge
//! order) and a chained *count* bit (occurrence counter). Matching uses
//! `uniquify=false` — every embedding, including symmetric ones. The molecule must
//! be ground.
//!
//! Templates so far: the acyclic library plus the small-ring cycles (3–6 rings, as
//! `*`-atom cycles — a cycle pattern already forces ring membership, so no `[R]`
//! constraint is needed). The four ring-junction motifs (RDKit pattern indices
//! 9–12) are pending: they need ring-bond (`@`) matching, which depends on
//! ring-membership facts being evaluated during substructure matching.

use std::sync::LazyLock;

use umol_ast::ast::{
    AsLit, AtomId, BondId, MoleculeAst, SubstructureMatchAlgorithm, SubstructureMatchConfig,
};
use umol_ast::mol_dsl;
use umol_graph_core::{NodeId, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm};

use super::bit_fp::BitFp;
use super::feature_set::FeatureSet;
use super::featurizer::FingerprintError;
use crate::hash::gboost_combine;

/// RDKit default fold width for `PatternFingerprint`.
pub const PATTERN_FP_WIDTH: usize = 2048;

/// RDKit's occurrence-counter salt (`0xBEEF`) chained into the count bit per match.
const COUNT_SALT: u32 = 0xBEEF;

/// Bit-exact replica of RDKit's `PatternFingerprint`, folded to `width` bits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PatternFingerprinter {
    pub width: usize,
    pub match_algorithm: SubstructureMatchAlgorithm,
    pub subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    pub relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
}

impl Default for PatternFingerprinter {
    fn default() -> Self {
        Self::new()
    }
}

impl PatternFingerprinter {
    pub fn new() -> Self {
        Self {
            width: PATTERN_FP_WIDTH,
            match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        }
    }

    /// `mol` must be ground. Returns the folded pattern fingerprint.
    pub fn fingerprint(&self, mol: &MoleculeAst) -> Result<BitFp, FingerprintError> {
        if !mol.is_ground() {
            return Err(FingerprintError::NotGround);
        }
        let mut ids: Vec<u64> = Vec::new();
        for template in TEMPLATES.iter() {
            let matches = template.pattern.substructure_matches(
                mol,
                SubstructureMatchConfig {
                    match_algorithm: self.match_algorithm,
                    subgraph_isomorphism_algorithm: self.subgraph_isomorphism_algorithm,
                    relevant_cycle_algorithm: self.relevant_cycle_algorithm,
                },
            );
            let mut count_id = template.index + template.atom_count + template.bond_count;
            for embedding in &matches {
                count_id = gboost_combine(count_id, COUNT_SALT);
                ids.push(u64::from(count_id));

                let host: Vec<AtomId> = embedding
                    .atoms()
                    .matched_pairs()
                    .iter()
                    .map(|&(_, host)| AtomId::from(host))
                    .collect();
                let mut bit_id = template.index;
                for &atom in &host {
                    let atomic_number = mol
                        .atom(atom)
                        .element()
                        .as_lit()
                        .expect("ground atom")
                        .atomic_number();
                    bit_id = gboost_combine(bit_id, u32::from(atomic_number));
                }
                for &(query_a, query_b) in &template.bonds {
                    let edge = mol
                        .raw_graph()
                        .find_edge(NodeId::from(host[query_a]), NodeId::from(host[query_b]))
                        .expect("matched bond");
                    let bond = mol.bond(BondId::from(edge));
                    let bond_type = if bond.is_in_aromatic_system() {
                        12
                    } else {
                        bond.order().as_lit().expect("ground bond") as u32
                    };
                    bit_id = gboost_combine(bit_id, bond_type);
                }
                ids.push(u64::from(bit_id));
            }
        }
        FeatureSet::from_features(ids).fold(self.width)
    }
}

/// A precomputed RDKit pattern template: its 1-based `pqs[]` index, the pattern
/// molecule, the atom/bond counts, and the pattern-edge order — all consumed by the
/// keying, so they are derived once at first use rather than per fingerprint.
struct PatternTemplate {
    index: u32,
    pattern: MoleculeAst,
    atom_count: u32,
    bond_count: u32,
    bonds: Vec<(usize, usize)>,
}

/// The RDKit template library, parsed once. Atom and bond order within each
/// template mirrors RDKit's SMARTS parse so the keying is bit-identical.
static TEMPLATES: LazyLock<Vec<PatternTemplate>> = LazyLock::new(|| {
    [
        (1, mol_dsl!(r#"{:atoms ["*" "*"] :bonds [[0 1 "*"]]}"#)),
        (2, mol_dsl!(r#"{:atoms ["*" "*" "*"] :bonds [[0 1 "*"] [1 2 "*"]]}"#)),
        (3, mol_dsl!(r#"{:atoms ["*" "*" "*"] :bonds [[0 1 "*"] [1 2 "*"] [2 0 "*"]]}"#)),
        (4, mol_dsl!(r#"{:atoms ["*" "*" "*" "*"] :bonds [[0 1 "*"] [1 2 "*"] [1 3 "*"]]}"#)),
        (5, mol_dsl!(r#"{:atoms ["*" "*" "*" "*"] :bonds [[0 1 "*"] [1 2 "*"] [2 3 "*"] [3 0 "*"]]}"#)),
        (6, mol_dsl!(r#"{:atoms ["*" "*" "*" "*" "*"] :bonds [[0 1 "*"] [1 2 "*"] [2 3 "*"] [2 4 "*"]]}"#)),
        (7, mol_dsl!(r#"{:atoms ["*" "*" "*" "*" "*"] :bonds [[0 1 "*"] [1 2 "*"] [2 3 "*"] [3 4 "*"] [4 0 "*"]]}"#)),
        (8, mol_dsl!(r#"{:atoms ["*" "*" "*" "*" "*" "*"] :bonds [[0 1 "*"] [1 2 "*"] [2 3 "*"] [3 4 "*"] [4 5 "*"] [5 0 "*"]]}"#)),
        (13, mol_dsl!(r#"{:atoms ["*"] :bonds []}"#)),
    ]
    .into_iter()
    .map(|(index, pattern)| {
        let graph = pattern.raw_graph();
        let atom_count = graph.node_count() as u32;
        let bond_count = graph.edge_count() as u32;
        let bonds = graph
            .edge_ids()
            .map(|edge| {
                let [a, b] = graph.edge_endpoints(edge);
                (a.index(), b.index())
            })
            .collect();
        PatternTemplate {
            index,
            pattern,
            atom_count,
            bond_count,
            bonds,
        }
    })
    .collect()
});

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_core::ARCMATCH_DEFAULT_PATH_LENGTH;

    use super::*;
    use crate::ingest::ingest_smiles;

    const ETHANOL_BITS: &[usize] = &[
        54, 173, 217, 429, 622, 759, 778, 874, 946, 967, 1022, 1033, 1061, 1236, 1289, 1295,
    ];
    const BENZENE_BITS: &[usize] = &[
        173, 217, 230, 389, 394, 410, 429, 434, 465, 469, 523, 527, 550, 617, 622, 663, 702, 789,
        797, 898, 923, 963, 967, 972, 1003, 1007, 1022, 1033, 1061, 1155, 1164, 1165, 1182, 1236,
        1295, 1387, 1414, 1416, 1417, 1449, 1465, 1531, 1607, 1764, 1982, 2019,
    ];

    // RDKit 2026.03.3 `PatternFingerprint` on-bits (width 2048). CCO is acyclic, so
    // its bits come solely from the non-ring templates.
    #[rstest]
    fn test_pattern_fingerprinter_new() {
        assert_eq!(
            PatternFingerprinter::new(),
            PatternFingerprinter {
                width: PATTERN_FP_WIDTH,
                match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
                subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            }
        );
    }

    #[rstest]
    #[case::ethanol("CCO", ETHANOL_BITS)]
    #[case::benzene("c1ccccc1", BENZENE_BITS)]
    #[case::pyridine(
        "c1ccncc1",
        &[16, 173, 217, 222, 230, 358, 389, 394, 410, 429, 434, 465, 469, 523, 527, 550, 617, 622, 655, 663, 667, 702, 789, 797, 898, 922, 923, 963, 967, 972, 980, 1003, 1007, 1022, 1033, 1061, 1155, 1164, 1165, 1182, 1236, 1294, 1295, 1414, 1416, 1417, 1449, 1465, 1531, 1590, 1607, 1764, 1883, 1903, 1982, 2019]
    )]
    #[case::furan(
        "c1ccoc1",
        &[43,173,217,230,323,333,429,465,469,495,527,555,617,622,687,702,746,786,789,834,923,928,963,967,986,1003,1022,1033,1061,1132,1155,1165,1182,1236,1277,1289,1295,1414,1416,1417,1476,1583,1671,1722,1757,1764,1945,1980,2019]
    )]
    #[case::chlorobenzene(
        "Clc1ccccc1",
        &[173,217,230,276,285,296,343,389,394,410,429,434,465,469,474,490,512,523,527,550,571,617,622,663,673,702,789,797,865,898,923,963,967,972,1003,1007,1022,1033,1061,1092,1124,1155,1164,1165,1182,1185,1189,1236,1280,1295,1328,1364,1387,1398,1414,1416,1417,1449,1460,1465,1531,1558,1573,1607,1645,1697,1713,1764,1891,1961,1982,2017,2019]
    )]
    #[case::cyclohexane(
        "C1CCCCC1",
        &[173,217,389,394,410,429,434,465,469,523,527,550,617,622,663,702,778,797,898,923,963,967,972,1003,1007,1022,1033,1061,1155,1164,1165,1177,1182,1236,1295,1414,1416,1417,1449,1465,1531,1607,1764,1812,1982,2019]
    )]
    fn test_pattern_fingerprinter_fingerprint(#[case] smiles: &str, #[case] expected: &[usize]) {
        let mol = ingest_smiles(smiles).expect("ingest");
        let fingerprint = PatternFingerprinter::new().fingerprint(&mol).unwrap();
        let on_bits: Vec<usize> = (0..fingerprint.width())
            .filter(|&bit| fingerprint.get(bit) == Some(true))
            .collect();
        assert_eq!(on_bits, expected);
    }

    #[rstest]
    #[case::zero_width(0, FingerprintError::ZeroWidth)]
    fn test_pattern_fingerprinter_fingerprint_error(
        #[case] width: usize,
        #[case] expected: FingerprintError,
    ) {
        assert_eq!(
            PatternFingerprinter {
                width,
                ..PatternFingerprinter::new()
            }
            .fingerprint(&ingest_smiles("CCO").expect("ingest")),
            Err(expected)
        );
    }

    #[rstest]
    #[case::incidence(
        SubstructureMatchAlgorithm::Incidence,
        SubgraphIsomorphismAlgorithm::Vf2
    )]
    #[case::ullmann(
        SubstructureMatchAlgorithm::GraphAndOverlays,
        SubgraphIsomorphismAlgorithm::Ullmann
    )]
    #[case::ri(
        SubstructureMatchAlgorithm::GraphAndOverlays,
        SubgraphIsomorphismAlgorithm::Ri
    )]
    #[case::arc_match(
        SubstructureMatchAlgorithm::GraphAndOverlays,
        SubgraphIsomorphismAlgorithm::ArcMatch {
            path_length: ARCMATCH_DEFAULT_PATH_LENGTH,
        }
    )]
    #[case::vf2_rdkit(
        SubstructureMatchAlgorithm::GraphAndOverlays,
        SubgraphIsomorphismAlgorithm::Vf2Rdkit
    )]
    #[case::ray_kirsch(
        SubstructureMatchAlgorithm::GraphAndOverlays,
        SubgraphIsomorphismAlgorithm::RayKirsch
    )]
    fn test_pattern_fingerprinter_fingerprint_algorithm(
        #[case] match_algorithm: SubstructureMatchAlgorithm,
        #[case] subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    ) {
        for (smiles, expected) in [("CCO", ETHANOL_BITS), ("c1ccccc1", BENZENE_BITS)] {
            let fingerprint = PatternFingerprinter {
                width: PATTERN_FP_WIDTH,
                match_algorithm,
                subgraph_isomorphism_algorithm,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            }
            .fingerprint(&ingest_smiles(smiles).expect("ingest"))
            .unwrap();

            assert_eq!(
                (0..fingerprint.width())
                    .filter(|&bit| fingerprint.get(bit) == Some(true))
                    .collect::<Vec<_>>(),
                expected
            );
        }
    }
}
