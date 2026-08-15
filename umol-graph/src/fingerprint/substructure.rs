//! Unhashed substructure fingerprint (direction A): every single atom and every
//! connected subgraph up to a bond-count bound, reduced to an exact,
//! collision-free canonical key.
//!
//! The key is the configured canonical form of the (edge-induced) subgraph, with
//! bond orders carried as label-colored vertices through edge subdivision. It
//! uses only embedding-monotone labels — element, charge, bond order — so
//! `query.is_subset(target)` is a sound substructure prescreen (no false
//! negatives). Degree, ring membership, and aromaticity are excluded: they change
//! with context and would break the subset test.

use std::ops::ControlFlow;

use umol_graph_core::{AutomorphismAlgorithm, GraphCorrespondence, SubgraphEnumerationAlgorithm};
use umol_graph_ir::ir::{AsLit, AtomId, BondId, Molecule};

use super::feature_set::FeatureSet;
use super::featurizer::FingerprintError;

/// Unhashed substructure fingerprint over single atoms and connected subgraphs of
/// ≤ `max_bonds` bonds. Connected subgraphs subsume simple paths, so paths are not
/// enumerated separately here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SubstructureFeaturizer {
    pub max_bonds: u32,
    pub subgraph_enumeration_algorithm: SubgraphEnumerationAlgorithm,
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

impl SubstructureFeaturizer {
    pub fn new(max_bonds: u32) -> Self {
        Self {
            max_bonds,
            subgraph_enumeration_algorithm: SubgraphEnumerationAlgorithm::Esu,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        }
    }

    /// `mol` must be ground. Returns the set of canonical structural keys for
    /// every atom and every connected subgraph up to `max_bonds` bonds.
    pub fn featurize(&self, mol: &Molecule) -> Result<FeatureSet<Vec<u8>>, FingerprintError> {
        if !mol.is_concrete() {
            return Err(FingerprintError::NotConcrete);
        }
        let graph = mol.raw_graph();
        let mut keys: Vec<Vec<u8>> = Vec::new();
        for node in graph.node_ids() {
            keys.push(canonical_key(
                mol,
                &graph.induced_subgraph(&[node]),
                self.automorphism_algorithm,
            ));
        }
        let _: ControlFlow<()> = graph.visit_connected_subgraphs(
            self.max_bonds,
            self.subgraph_enumeration_algorithm,
            |edges| {
                keys.push(canonical_key(
                    mol,
                    &graph.edge_induced_subgraph(edges),
                    self.automorphism_algorithm,
                ));
                ControlFlow::Continue(())
            },
        );
        Ok(FeatureSet::from_features(keys))
    }
}

/// Exact canonical key of the subgraph correspondence `sub`, with bond orders
/// carried as label-colored vertices. Two subgraphs collide iff they are
/// isomorphic as element/charge/bond-order-labeled graphs.
fn canonical_key(
    mol: &Molecule,
    sub: &GraphCorrespondence,
    automorphism_algorithm: AutomorphismAlgorithm,
) -> Vec<u8> {
    let subgraph = mol.raw_graph().extract(sub);
    subgraph.canonical_key(
        |node| {
            let id = AtomId::from(
                sub.nodes()
                    .right_of(node)
                    .expect("subgraph node maps to a host node"),
            );
            let view = mol.atom(id);
            let atomic_number = view
                .element()
                .as_lit()
                .expect("ground atom")
                .atomic_number();
            let charge = view.charge().as_lit().expect("ground atom") as i16;
            let mut color = Vec::with_capacity(4);
            color.extend_from_slice(&u16::from(atomic_number).to_le_bytes());
            color.extend_from_slice(&charge.to_le_bytes());
            color
        },
        |edge| {
            let id = BondId::from(
                sub.edges()
                    .right_of(edge)
                    .expect("subgraph edge maps to a host edge"),
            );
            let order = mol.bond(id).order().as_lit().expect("ground bond") as u16;
            order.to_le_bytes().to_vec()
        },
        automorphism_algorithm,
    )
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::mol_dsl_concrete;

    use super::*;

    const ETHANOL: &str = r#"{:atoms ["C #h3" "C #h2" "O #h1"] :bonds [[0 1 "1"] [1 2 "1"]]}"#;
    const PROPANE: &str = r#"{:atoms ["C #h3" "C #h2" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#;
    const ETHANE: &str = r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#;

    #[rstest]
    fn test_substructure_featurizer_new() {
        assert_eq!(
            SubstructureFeaturizer::new(2),
            SubstructureFeaturizer {
                max_bonds: 2,
                subgraph_enumeration_algorithm: SubgraphEnumerationAlgorithm::Esu,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        );
    }

    // Distinct ethanol features by bond bound: atoms {C, O}=2; +bonds {C-C, C-O}=4;
    // +the two-bond path {C-C-O}=5.
    #[rstest]
    #[case::bonds_0(0, 2)]
    #[case::bonds_1(1, 4)]
    #[case::bonds_2(2, 5)]
    fn test_substructure_featurizer_featurize(#[case] max_bonds: u32, #[case] expected: usize) {
        let fingerprint = SubstructureFeaturizer::new(max_bonds)
            .featurize(&mol_dsl_concrete!(ETHANOL))
            .unwrap();
        assert_eq!(fingerprint.len(), expected);
    }

    // Ethane's C-C is a substructure of propane: every ethane feature key appears
    // in propane, but propane's C-C-C path does not appear in ethane.
    #[rstest]
    fn test_substructure_featurizer_featurize_subset() {
        let featurizer = SubstructureFeaturizer::new(2);
        let ethane = featurizer.featurize(&mol_dsl_concrete!(ETHANE)).unwrap();
        let propane = featurizer.featurize(&mol_dsl_concrete!(PROPANE)).unwrap();
        assert!(ethane.is_subset(&propane));
        assert!(!propane.is_subset(&ethane));
    }

    // The same molecule under two atom numberings yields identical keys (the
    // canonical form is numbering-invariant).
    #[rstest]
    fn test_substructure_featurizer_featurize_order_independent() {
        let featurizer = SubstructureFeaturizer::new(2);
        let forward = featurizer.featurize(&mol_dsl_concrete!(ETHANOL)).unwrap();
        let relabeled = featurizer
            .featurize(&mol_dsl_concrete!(
                r#"{:atoms ["O #h1" "C #h2" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
            ))
            .unwrap();
        assert_eq!(forward, relabeled);
    }

    // A single and a double bond between carbons are distinct features: neither
    // fingerprint screens into the other.
    #[rstest]
    fn test_substructure_featurizer_featurize_bond_order() {
        let featurizer = SubstructureFeaturizer::new(1);
        let ethane = featurizer.featurize(&mol_dsl_concrete!(ETHANE)).unwrap();
        let ethene = featurizer
            .featurize(&mol_dsl_concrete!(
                r#"{:atoms ["C #h2" "C #h2"] :bonds [[0 1 "2"]]}"#
            ))
            .unwrap();
        assert!(!ethane.is_subset(&ethene));
        assert!(!ethene.is_subset(&ethane));
    }
}
