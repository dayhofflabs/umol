//! Resolution conformance tests

#[cfg(feature = "conformance")]
#[path = "resolution/mod.rs"]
pub mod resolution_suite;

/// Tests verifying that `StableGraph` and `MoleculeBuilder` permit self-loops and
/// parallel edges at the data-structure level, without any topology validation.
///
/// These confirm that topology checking must be a separate pass on a fully-constructed
/// builder, not a constraint enforced by the underlying graph type.
#[cfg(test)]
mod topology_graph_behavior {
    use petgraph::prelude::Undirected;
    use petgraph::stable_graph::StableGraph;
    use umol_shared::element::Element;
    use umol_graph::graph_ir::atom_pattern::AtomPattern;
    use umol_graph::graph_ir::bond_pattern::BondPattern;
    use umol_graph::graph_ir::molecule_builder::MoleculeBuilder;

    #[test]
    fn stable_graph_allows_self_loop() {
        let mut g: StableGraph<(), (), Undirected, u32> = StableGraph::default();
        let n = g.add_node(());
        let e = g.add_edge(n, n, ());
        assert!(g.edge_endpoints(e) == Some((n, n)));
    }

    #[test]
    fn stable_graph_allows_parallel_edges() {
        let mut g: StableGraph<(), (), Undirected, u32> = StableGraph::default();
        let a = g.add_node(());
        let b = g.add_node(());
        let e1 = g.add_edge(a, b, ());
        let e2 = g.add_edge(a, b, ());
        assert_ne!(e1, e2);
        assert_eq!(g.edge_count(), 2);
    }

    #[test]
    fn molecule_builder_allows_self_loop_bond() {
        let mut builder = MoleculeBuilder::new();
        let atom = builder.add_atom(AtomPattern::new(Element::C));
        let bond_idx = builder.add_bond_unchecked(atom, atom, BondPattern::new(1));
        assert_eq!(builder.bond_count(), 1);
        let (a, b) = builder.bond_atom_indices(bond_idx).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn molecule_builder_allows_parallel_bonds() {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(AtomPattern::new(Element::C));
        let b = builder.add_atom(AtomPattern::new(Element::C));
        let e1 = builder.add_bond_unchecked(a, b, BondPattern::new(1));
        let e2 = builder.add_bond_unchecked(a, b, BondPattern::new(2));
        assert_ne!(e1, e2);
        assert_eq!(builder.bond_count(), 2);
    }
}
