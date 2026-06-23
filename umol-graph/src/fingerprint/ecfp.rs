//! ECFP featurizer (Rogers & Hahn 2010): circular refinement over Daylight atom
//! invariants, with structural (bond-set) duplicate removal.
//!
//! The molecule must be ground; the caller ([`super::Featurizer::featurize`])
//! guarantees it, so the seed reads concrete literals directly. The hash is
//! `xxh3_64` with [`ECFP_SEED`]; the paper leaves the hash unspecified, so this is
//! a frozen choice (placeholder identity, like the WL schemes).

use std::collections::{BTreeSet, HashMap};

use umol_ast::ast::{AsLit, AtomId, BondId, MoleculeAst};
use umol_graph_core::{CircularRefinementAlgorithm, EdgeId, NodeId, TraversalAlgorithm};

use super::feature_set::FeatureSet;

/// Frozen ECFP hash seed. Placeholder identity, not yet formalized.
pub const ECFP_SEED: u64 = 0xECF0_5EED_0000_0001;

/// ECFP fingerprint of `radius` iterations (diameter `2 * radius`, i.e. ECFP_{2r}).
#[derive(Clone, Copy, Debug)]
pub struct EcfpFeaturizer {
    pub radius: u32,
    pub seed: u64,
}

impl EcfpFeaturizer {
    pub fn new(radius: u32) -> Self {
        Self {
            radius,
            seed: ECFP_SEED,
        }
    }

    /// `mol` must be ground. Returns the deduplicated set of feature identifiers.
    pub fn featurize(&self, mol: &MoleculeAst) -> FeatureSet<u64> {
        let graph = mol.raw_graph();
        let node_count = graph.node_count();

        let atom_labels: Vec<u64> = (0..node_count)
            .map(|i| atom_seed(mol, AtomId(i as u32)))
            .collect();
        let bond_labels: Vec<u64> = (0..graph.edge_count())
            .map(|e| bond_seed(mol, BondId(e as u32)))
            .collect();

        let rounds = graph.circular_refine(
            |node: NodeId| atom_labels[node.index()],
            |edge: EdgeId| bond_labels[edge.index()],
            self.radius,
            CircularRefinementAlgorithm::RogersHahn { seed: self.seed },
        );

        // Round 0: initial identifiers enter directly (identifier-dedup via the set).
        let mut identifiers: Vec<u64> = rounds[0].clone();

        // Rounds >= 1: structural dedup by covered bond set, keeping the feature with
        // the smallest (round, identifier) per set (Rogers & Hahn rules 1 and 2).
        let mut kept: HashMap<Vec<u32>, (u32, u64)> = HashMap::new();
        if self.radius >= 1 {
            for atom in 0..node_count {
                let source = NodeId(atom as u32);
                let neighborhood =
                    graph.neighborhood(source, self.radius - 1, TraversalAlgorithm::Bfs);
                let mut bond_set: BTreeSet<u32> = BTreeSet::new();
                let mut shell = 0;
                for round in 1..=self.radius {
                    // The radius-`round` bond set is every edge incident to a node
                    // within distance `round - 1`; add the newly reached shell.
                    while shell < neighborhood.len() && neighborhood[shell].1 == round - 1 {
                        for neighbor in graph.neighbors(neighborhood[shell].0) {
                            bond_set.insert(neighbor.edge.index() as u32);
                        }
                        shell += 1;
                    }
                    let identifier = rounds[round as usize][source.index()];
                    let key: Vec<u32> = bond_set.iter().copied().collect();
                    kept.entry(key)
                        .and_modify(|best| {
                            if (round, identifier) < *best {
                                *best = (round, identifier);
                            }
                        })
                        .or_insert((round, identifier));
                }
            }
        }

        identifiers.extend(kept.values().map(|&(_, identifier)| identifier));
        FeatureSet::from_features(identifiers)
    }
}

/// Daylight initial invariant (Rogers & Hahn 2010): heavy-atom degree, heavy
/// valence, atomic number, isotope mass, formal charge, attached hydrogens, and
/// ring membership, bit-packed into disjoint fields (the engine rehashes it).
fn atom_seed(mol: &MoleculeAst, id: AtomId) -> u64 {
    let atom = mol.atom(id);
    let atomic_number = atom.element().as_lit().expect("ground atom").atomic_number();
    let heavy_degree = atom.heavy_atom_degree().as_lit().expect("ground atom");
    let heavy_valence = atom.heavy_atom_valence().as_lit().expect("ground atom");
    let hydrogens = atom.total_hydrogens().as_lit().expect("ground atom");
    let charge = atom.charge().as_lit().expect("ground atom");
    let mass = atom.isotope_mass().as_lit().unwrap_or(0); // Natural isotope -> 0
    let in_ring = u64::from(atom.is_in_ring());

    u64::from(atomic_number)
        | ((heavy_degree as u64 & 0xF) << 8)
        | ((heavy_valence as u64 & 0xF) << 12)
        | ((hydrogens as u64 & 0xF) << 16)
        | (u64::from(charge as u8) << 20)
        | ((u64::from(mass) & 0xFFFF) << 28)
        | (in_ring << 44)
}

/// Seed a bond from (bond order, aromatic-system membership).
fn bond_seed(mol: &MoleculeAst, id: BondId) -> u64 {
    let bond = mol.bond(id);
    let order = bond.order().as_lit().expect("ground bond");
    (order as u16 as u64) | ((bond.is_in_aromatic_system() as u64) << 16)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::mol_ground;

    use super::*;

    const BUTYRAMIDE: &str = r#"{
        :atoms ["C #h3" "C #h2" "C #h2" "C #h0" "O #h0" "N #h2"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "2"] [3 5 "1"]]
    }"#;

    // Rogers & Hahn 2010, Figure 8: butyramide feature counts per diameter.
    #[rstest]
    #[case::diameter_0(0, 5)]
    #[case::diameter_2(1, 11)]
    #[case::diameter_4(2, 14)]
    #[case::diameter_6(3, 14)]
    fn test_ecfp_featurizer_featurize_butyramide(#[case] radius: u32, #[case] expected: usize) {
        let fingerprint = EcfpFeaturizer::new(radius).featurize(&mol_ground!(BUTYRAMIDE));
        assert_eq!(fingerprint.len(), expected);
    }

    #[rstest]
    #[case::relabeled_propane(
        r#"{:atoms ["C #h3" "C #h2" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#,
        r#"{:atoms ["C #h2" "C #h3" "C #h3"] :bonds [[0 1 "1"] [0 2 "1"]]}"#
    )]
    fn test_ecfp_featurizer_featurize_order_independent(#[case] a: &str, #[case] b: &str) {
        let featurizer = EcfpFeaturizer::new(2);
        assert_eq!(
            featurizer.featurize(&mol_ground!(a)),
            featurizer.featurize(&mol_ground!(b))
        );
    }
}
