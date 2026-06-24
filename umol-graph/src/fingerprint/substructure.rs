//! Unhashed substructure fingerprint (direction A): every single atom and every
//! connected subgraph up to a bond-count bound, reduced to an exact,
//! collision-free canonical key.
//!
//! The key is the nauty canonical form of the (edge-induced) subgraph, with bond
//! orders carried as label-colored vertices (edge subdivision, since nauty colors
//! vertices only). It uses only embedding-monotone labels — element, charge, bond
//! order — so `query.is_subset(target)` is a sound substructure prescreen (no
//! false negatives). Degree, ring membership, and aromaticity are excluded: they
//! change with context and would break the subset test.

use umol_ast::ast::{AsLit, AtomId, BondId, MoleculeAst};
use umol_graph_core::{
    AutomorphismAlgorithm, EdgeId, Embedding, Graph, NodeId, SubgraphEnumerationAlgorithm,
};

use super::feature_set::FeatureSet;
use super::featurizer::FingerprintError;

/// Unhashed substructure fingerprint over single atoms and connected subgraphs of
/// ≤ `max_bonds` bonds. Connected subgraphs subsume simple paths, so paths are not
/// enumerated separately here.
#[derive(Clone, Copy, Debug)]
pub struct SubstructureFeaturizer {
    pub max_bonds: u32,
}

impl SubstructureFeaturizer {
    pub fn new(max_bonds: u32) -> Self {
        Self { max_bonds }
    }

    /// `mol` must be ground. Returns the set of canonical structural keys for
    /// every atom and every connected subgraph up to `max_bonds` bonds.
    pub fn featurize(&self, mol: &MoleculeAst) -> Result<FeatureSet<Vec<u8>>, FingerprintError> {
        if !mol.is_ground() {
            return Err(FingerprintError::NotGround);
        }
        let graph = mol.raw_graph();
        let mut keys: Vec<Vec<u8>> = Vec::new();
        for node in graph.node_ids() {
            keys.push(canonical_key(mol, &graph.induced_subgraph(&[node])));
        }
        for edges in
            graph.enumerate_connected_subgraphs(self.max_bonds, SubgraphEnumerationAlgorithm::Esu)
        {
            keys.push(canonical_key(mol, &graph.edge_induced_subgraph(&edges)));
        }
        Ok(FeatureSet::from_features(keys))
    }
}

/// Exact canonical key of `embedding`: its nauty canonical form with bond orders
/// carried as label-colored vertices. Two subgraphs collide iff they are
/// isomorphic as element/charge/bond-order-labeled graphs.
fn canonical_key(mol: &MoleculeAst, embedding: &Embedding) -> Vec<u8> {
    let subgraph = embedding.extract();
    let atom_count = subgraph.node_count();
    let bond_count = subgraph.edge_count();

    // Subdivide each bond into a vertex so nauty (vertex-color only) canonicalizes
    // bond orders: atom nodes 0..atom_count, bond nodes atom_count.. .
    let mut subdivided_edges: Vec<[u32; 2]> = Vec::with_capacity(2 * bond_count);
    for bond in 0..bond_count {
        let [a, b] = subgraph.edge_endpoints(EdgeId(bond as u32));
        let bond_node = (atom_count + bond) as u32;
        subdivided_edges.push([a.0, bond_node]);
        subdivided_edges.push([bond_node, b.0]);
    }
    let subdivided = Graph::new(atom_count + bond_count, &subdivided_edges);

    // Colors: (class, value, signed). Atoms (class 0) carry element + charge; bond
    // nodes (class 1) carry bond order. The class split keeps the two disjoint.
    let mut colors: Vec<(u8, u16, i16)> = Vec::with_capacity(atom_count + bond_count);
    for atom in 0..atom_count {
        let id = AtomId::from(embedding.host_node(NodeId(atom as u32)));
        let view = mol.atom(id);
        let atomic_number = view.element().as_lit().expect("ground atom").atomic_number();
        let charge = view.charge().as_lit().expect("ground atom") as i16;
        colors.push((0, u16::from(atomic_number), charge));
    }
    for bond in 0..bond_count {
        let id = BondId::from(embedding.host_edge(EdgeId(bond as u32)));
        let order = mol.bond(id).order().as_lit().expect("ground bond") as u16;
        colors.push((1, order, 0));
    }

    let automorphism =
        subdivided.automorphisms(|node| colors[node.index()], AutomorphismAlgorithm::Nauty);
    let canonical = automorphism.canonical_labeling();
    let mut position = vec![0u32; atom_count + bond_count];
    for (rank, node) in canonical.iter().enumerate() {
        position[node.index()] = rank as u32;
    }

    let mut key: Vec<u8> = Vec::new();
    key.extend_from_slice(&((atom_count + bond_count) as u32).to_le_bytes());
    for &node in canonical {
        let (class, value, signed) = colors[node.index()];
        key.push(class);
        key.extend_from_slice(&value.to_le_bytes());
        key.extend_from_slice(&signed.to_le_bytes());
    }
    let mut canonical_edges: Vec<(u32, u32)> = subdivided
        .edge_ids()
        .map(|edge| {
            let [u, v] = subdivided.edge_endpoints(edge);
            let (u, v) = (position[u.index()], position[v.index()]);
            (u.min(v), u.max(v))
        })
        .collect();
    canonical_edges.sort_unstable();
    key.extend_from_slice(&(canonical_edges.len() as u32).to_le_bytes());
    for (u, v) in canonical_edges {
        key.extend_from_slice(&u.to_le_bytes());
        key.extend_from_slice(&v.to_le_bytes());
    }
    key
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::mol_ground;

    use super::*;

    const ETHANOL: &str = r#"{:atoms ["C #h3" "C #h2" "O #h1"] :bonds [[0 1 "1"] [1 2 "1"]]}"#;
    const PROPANE: &str = r#"{:atoms ["C #h3" "C #h2" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#;
    const ETHANE: &str = r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#;

    // Distinct ethanol features by bond bound: atoms {C, O}=2; +bonds {C-C, C-O}=4;
    // +the two-bond path {C-C-O}=5.
    #[rstest]
    #[case::bonds_0(0, 2)]
    #[case::bonds_1(1, 4)]
    #[case::bonds_2(2, 5)]
    fn test_substructure_featurizer_featurize(#[case] max_bonds: u32, #[case] expected: usize) {
        let fingerprint = SubstructureFeaturizer::new(max_bonds)
            .featurize(&mol_ground!(ETHANOL))
            .unwrap();
        assert_eq!(fingerprint.len(), expected);
    }

    // Ethane's C-C is a substructure of propane: every ethane feature key appears
    // in propane, but propane's C-C-C path does not appear in ethane.
    #[rstest]
    fn test_substructure_featurizer_featurize_subset() {
        let featurizer = SubstructureFeaturizer::new(2);
        let ethane = featurizer.featurize(&mol_ground!(ETHANE)).unwrap();
        let propane = featurizer.featurize(&mol_ground!(PROPANE)).unwrap();
        assert!(ethane.is_subset(&propane));
        assert!(!propane.is_subset(&ethane));
    }

    // The same molecule under two atom numberings yields identical keys (the
    // canonical form is numbering-invariant).
    #[rstest]
    fn test_substructure_featurizer_featurize_order_independent() {
        let featurizer = SubstructureFeaturizer::new(2);
        let forward = featurizer.featurize(&mol_ground!(ETHANOL)).unwrap();
        let relabeled = featurizer
            .featurize(&mol_ground!(
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
        let ethane = featurizer.featurize(&mol_ground!(ETHANE)).unwrap();
        let ethene = featurizer
            .featurize(&mol_ground!(r#"{:atoms ["C #h2" "C #h2"] :bonds [[0 1 "2"]]}"#))
            .unwrap();
        assert!(!ethane.is_subset(&ethene));
        assert!(!ethene.is_subset(&ethane));
    }
}
