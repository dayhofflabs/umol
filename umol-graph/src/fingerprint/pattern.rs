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

use std::ops::ControlFlow;
use std::sync::{LazyLock, OnceLock};

use umol_graph_core::{
    AutomorphismAlgorithm, NodeId, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::{
    AsLit, Molecule, PreparedSubstructureHost, PreparedSubstructurePattern,
    SubstructureMatchAlgorithm, SubstructureMatchConfig,
};
use umol_graph_ir::mol_dsl;

use super::bit_fp::BitFp;
use super::featurizer::FingerprintError;
use crate::hash::gboost_combine;

/// RDKit default fold width for `PatternFingerprint`.
pub const PATTERN_FP_WIDTH: usize = 2048;

/// RDKit's occurrence-counter salt (`0xBEEF`) chained into the count bit per match.
const COUNT_SALT: u32 = 0xBEEF;

const MATCH_ALGORITHM_COUNT: usize = 2;
const SUBISO_ALGORITHM_COUNT: usize = 4;
const RELEVANT_CYCLE_ALGORITHM_COUNT: usize = 1;
const PREPARED_TEMPLATE_SLOTS: usize =
    MATCH_ALGORITHM_COUNT * SUBISO_ALGORITHM_COUNT * RELEVANT_CYCLE_ALGORITHM_COUNT;

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
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::BitsetRi,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        }
    }

    /// `mol` must be ground. Returns the folded pattern fingerprint.
    pub fn fingerprint(&self, mol: &Molecule) -> Result<BitFp, FingerprintError> {
        if !mol.is_concrete() {
            return Err(FingerprintError::NotConcrete);
        }
        if self.width == 0 {
            return Err(FingerprintError::ZeroWidth);
        }

        let atom_numbers: Vec<u32> = mol
            .atoms()
            .iter()
            .map(|atom| {
                u32::from(
                    atom.element()
                        .as_lit()
                        .expect("ground atom")
                        .atomic_number(),
                )
            })
            .collect();
        let bond_types: Vec<u32> = mol
            .bonds()
            .iter()
            .map(|bond| {
                if bond.is_in_aromatic_system() {
                    12
                } else {
                    bond.order().as_lit().expect("ground bond") as u32
                }
            })
            .collect();

        let host = PreparedSubstructureHost::new(mol);
        let graph = mol.raw_graph();
        let templates = prepared_templates(
            self.match_algorithm,
            self.subgraph_isomorphism_algorithm,
            self.relevant_cycle_algorithm,
        );
        let mut bits = BitFp::zeros(self.width);

        for template in templates {
            let mut count_id = template.index + template.atom_count + template.bond_count;
            let _: ControlFlow<()> = template.prepared.visit_embeddings_in_host(&host, |atoms| {
                count_id = gboost_combine(count_id, COUNT_SALT);
                bits.set((count_id as usize) % self.width);

                let mut bit_id = template.index;
                for &atom in atoms {
                    bit_id = gboost_combine(bit_id, atom_numbers[atom.index()]);
                }
                for &(u, v) in &template.edges {
                    let edge = graph
                        .find_edge(NodeId::from(atoms[u]), NodeId::from(atoms[v]))
                        .expect("matched bond");
                    bit_id = gboost_combine(bit_id, bond_types[edge.index()]);
                }
                bits.set((bit_id as usize) % self.width);
                ControlFlow::Continue(())
            });
        }

        Ok(bits)
    }
}

/// A precomputed RDKit pattern template: its 1-based `pqs[]` index, the pattern
/// molecule, the atom/bond counts, and the pattern-edge endpoints in bond-id order.
struct PatternTemplate {
    index: u32,
    pattern: Molecule,
    atom_count: u32,
    bond_count: u32,
    edges: Vec<(usize, usize)>,
}

/// A pattern template precompiled for one `(match_algorithm,
/// subgraph_isomorphism_algorithm, relevant_cycle_algorithm)` combination with
/// symmetry reduction disabled (`uniquify=false`).
struct PreparedPatternTemplate {
    index: u32,
    atom_count: u32,
    bond_count: u32,
    edges: Vec<(usize, usize)>,
    prepared: PreparedSubstructurePattern,
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
        let edges = graph
            .edge_ids()
            .map(|edge| {
                let [u, v] = graph.edge_endpoints(edge);
                (u.index(), v.index())
            })
            .collect();
        PatternTemplate {
            index,
            pattern,
            atom_count,
            bond_count,
            edges,
        }
    })
    .collect()
});

static PREPARED_TEMPLATES: [OnceLock<Vec<PreparedPatternTemplate>>; PREPARED_TEMPLATE_SLOTS] =
    [const { OnceLock::new() }; PREPARED_TEMPLATE_SLOTS];

fn prepared_template_slot(
    match_algorithm: SubstructureMatchAlgorithm,
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
) -> usize {
    let match_index = match match_algorithm {
        SubstructureMatchAlgorithm::GraphAndOverlays => 0,
        SubstructureMatchAlgorithm::Incidence => 1,
    };
    let subiso_index = match subgraph_isomorphism_algorithm {
        SubgraphIsomorphismAlgorithm::Vf2Rdkit => 0,
        SubgraphIsomorphismAlgorithm::RayKirsch => 1,
        SubgraphIsomorphismAlgorithm::BitsetRi => 2,
        SubgraphIsomorphismAlgorithm::Glasgow => 3,
    };
    let cycle_index = match relevant_cycle_algorithm {
        RelevantCycleEnumerationAlgorithm::Vismara => 0,
    };
    ((match_index * SUBISO_ALGORITHM_COUNT) + subiso_index) * RELEVANT_CYCLE_ALGORITHM_COUNT
        + cycle_index
}

fn prepared_templates(
    match_algorithm: SubstructureMatchAlgorithm,
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
) -> &'static [PreparedPatternTemplate] {
    let slot = prepared_template_slot(
        match_algorithm,
        subgraph_isomorphism_algorithm,
        relevant_cycle_algorithm,
    );
    PREPARED_TEMPLATES[slot]
        .get_or_init(|| {
            let config = SubstructureMatchConfig {
                match_algorithm,
                subgraph_isomorphism_algorithm,
                relevant_cycle_algorithm,
                automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                reduce_pattern_symmetry: false,
                reduce_host_symmetry: false,
            };
            TEMPLATES
                .iter()
                .map(|template| PreparedPatternTemplate {
                    index: template.index,
                    atom_count: template.atom_count,
                    bond_count: template.bond_count,
                    edges: template.edges.clone(),
                    prepared: PreparedSubstructurePattern::new(&template.pattern, config),
                })
                .collect()
        })
        .as_slice()
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

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
    const PYRIDINE_BITS: &[usize] = &[
        16, 173, 217, 222, 230, 358, 389, 394, 410, 429, 434, 465, 469, 523, 527, 550, 617, 622,
        655, 663, 667, 702, 789, 797, 898, 922, 923, 963, 967, 972, 980, 1003, 1007, 1022, 1033,
        1061, 1155, 1164, 1165, 1182, 1236, 1294, 1295, 1414, 1416, 1417, 1449, 1465, 1531, 1590,
        1607, 1764, 1883, 1903, 1982, 2019,
    ];
    const FURAN_BITS: &[usize] = &[
        43, 173, 217, 230, 323, 333, 429, 465, 469, 495, 527, 555, 617, 622, 687, 702, 746, 786,
        789, 834, 923, 928, 963, 967, 986, 1003, 1022, 1033, 1061, 1132, 1155, 1165, 1182, 1236,
        1277, 1289, 1295, 1414, 1416, 1417, 1476, 1583, 1671, 1722, 1757, 1764, 1945, 1980, 2019,
    ];
    const CHLOROBENZENE_BITS: &[usize] = &[
        173, 217, 230, 276, 285, 296, 343, 389, 394, 410, 429, 434, 465, 469, 474, 490, 512, 523,
        527, 550, 571, 617, 622, 663, 673, 702, 789, 797, 865, 898, 923, 963, 967, 972, 1003, 1007,
        1022, 1033, 1061, 1092, 1124, 1155, 1164, 1165, 1182, 1185, 1189, 1236, 1280, 1295, 1328,
        1364, 1387, 1398, 1414, 1416, 1417, 1449, 1460, 1465, 1531, 1558, 1573, 1607, 1645, 1697,
        1713, 1764, 1891, 1961, 1982, 2017, 2019,
    ];
    const CYCLOHEXANE_BITS: &[usize] = &[
        173, 217, 389, 394, 410, 429, 434, 465, 469, 523, 527, 550, 617, 622, 663, 702, 778, 797,
        898, 923, 963, 967, 972, 1003, 1007, 1022, 1033, 1061, 1155, 1164, 1165, 1177, 1182, 1236,
        1295, 1414, 1416, 1417, 1449, 1465, 1531, 1607, 1764, 1812, 1982, 2019,
    ];

    const RDKIT_CASES: &[(&str, &[usize])] = &[
        ("CCO", ETHANOL_BITS),
        ("c1ccccc1", BENZENE_BITS),
        ("c1ccncc1", PYRIDINE_BITS),
        ("c1ccoc1", FURAN_BITS),
        ("Clc1ccccc1", CHLOROBENZENE_BITS),
        ("C1CCCCC1", CYCLOHEXANE_BITS),
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
                subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::BitsetRi,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            }
        );
    }

    #[rstest]
    #[case::ethanol("CCO", ETHANOL_BITS)]
    #[case::benzene("c1ccccc1", BENZENE_BITS)]
    #[case::pyridine("c1ccncc1", PYRIDINE_BITS)]
    #[case::furan("c1ccoc1", FURAN_BITS)]
    #[case::chlorobenzene("Clc1ccccc1", CHLOROBENZENE_BITS)]
    #[case::cyclohexane("C1CCCCC1", CYCLOHEXANE_BITS)]
    fn test_pattern_fingerprinter_fingerprint(#[case] smiles: &str, #[case] expected: &[usize]) {
        let mol = ingest_smiles(smiles).expect("ingest");
        let fingerprint = PatternFingerprinter::new().fingerprint(&mol).unwrap();
        let on_bits: Vec<usize> = (0..fingerprint.width())
            .filter(|&bit| fingerprint.get(bit) == Some(true))
            .collect();
        assert_eq!(on_bits, expected);
    }

    #[rstest]
    fn test_pattern_fingerprinter_fingerprint_empty() {
        let fingerprint = PatternFingerprinter::new()
            .fingerprint(&Molecule::new())
            .unwrap();
        assert_eq!(fingerprint.width(), PATTERN_FP_WIDTH);
        assert_eq!(fingerprint.count_ones(), 0);
    }

    #[rstest]
    #[case::zero_width(ingest_smiles("CCO").expect("ingest"), 0, FingerprintError::ZeroWidth)]
    #[case::not_concrete(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        PATTERN_FP_WIDTH,
        FingerprintError::NotConcrete
    )]
    #[case::not_concrete_before_zero_width(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        0,
        FingerprintError::NotConcrete
    )]
    fn test_pattern_fingerprinter_fingerprint_error(
        #[case] mol: Molecule,
        #[case] width: usize,
        #[case] expected: FingerprintError,
    ) {
        assert_eq!(
            PatternFingerprinter {
                width,
                ..PatternFingerprinter::new()
            }
            .fingerprint(&mol),
            Err(expected)
        );
    }

    #[rstest]
    fn test_pattern_fingerprinter_fingerprint_algorithm(
        #[values(
            SubstructureMatchAlgorithm::GraphAndOverlays,
            SubstructureMatchAlgorithm::Incidence
        )]
        match_algorithm: SubstructureMatchAlgorithm,
        #[values(
            SubgraphIsomorphismAlgorithm::Vf2Rdkit,
            SubgraphIsomorphismAlgorithm::RayKirsch,
            SubgraphIsomorphismAlgorithm::BitsetRi,
            SubgraphIsomorphismAlgorithm::Glasgow
        )]
        subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    ) {
        let fingerprinter = PatternFingerprinter {
            width: PATTERN_FP_WIDTH,
            match_algorithm,
            subgraph_isomorphism_algorithm,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        };
        for &(smiles, expected) in RDKIT_CASES {
            let fingerprint = fingerprinter
                .fingerprint(&ingest_smiles(smiles).expect("ingest"))
                .unwrap();

            assert_eq!(
                (0..fingerprint.width())
                    .filter(|&bit| fingerprint.get(bit) == Some(true))
                    .collect::<Vec<_>>(),
                expected,
                "{smiles} with {match_algorithm:?}/{subgraph_isomorphism_algorithm:?}"
            );
        }
    }
}
