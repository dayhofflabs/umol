//! Graph view: typed adapter over the underlying `Graph`.

use std::ops::ControlFlow;

use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismGroupOrder, AutomorphismOutput,
    BiconnectedComponentsAlgorithm, BipartiteMaximumMatchingAlgorithm,
    ConnectedComponentsAlgorithm, EdgeId, GeneralMaximumMatchingAlgorithm, Graph,
    MatchingEnumerationAlgorithm, MaximumIndependentSetAlgorithm, NodeId, NonBipartiteGraphError,
    PerfectMatchingAlgorithm, RelevantCycleEnumerationAlgorithm, ShortestCycleAlgorithm,
    SimpleCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};

use super::super::id::{AtomId, BondId};
use super::super::matching::BondMatching;

/// AtomId/BondId-typed adapter over the underlying `Graph`. Holds the
/// pure-graph algorithms (connectivity, cycles, matchings, isomorphisms)
/// without exposing graph-core's `NodeId` / `EdgeId` types in the public
/// API. Construct via `Molecule::graph()`.
#[derive(Clone, Copy)]
pub struct GraphView<'a> {
    graph: &'a Graph,
}

impl<'a> GraphView<'a> {
    pub(crate) fn new(graph: &'a Graph) -> Self {
        Self { graph }
    }

    pub fn degree(&self, atom: AtomId) -> usize {
        self.graph.degree(NodeId::from(atom))
    }

    pub fn enumerate_connected_components(
        &self,
        alg: ConnectedComponentsAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .enumerate_connected_components(alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn enumerate_biconnected_components(
        &self,
        alg: BiconnectedComponentsAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .enumerate_biconnected_components(alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn shortest_cycle_through_bond(
        &self,
        bond: BondId,
        alg: ShortestCycleAlgorithm,
    ) -> Option<usize> {
        self.graph
            .shortest_cycle_through_edge(EdgeId::from(bond), alg)
    }

    pub fn shortest_cycle_through_atom(
        &self,
        atom: AtomId,
        alg: ShortestCycleAlgorithm,
    ) -> Option<usize> {
        self.graph
            .shortest_cycle_through_node(NodeId::from(atom), alg)
    }

    /// Visits each simple cycle's atoms until traversal completes or the visitor
    /// returns [`ControlFlow::Break`]. The slice borrows visit state and is only
    /// valid for the duration of the call.
    pub fn visit_simple_cycles<B>(
        &self,
        max_size: usize,
        alg: SimpleCycleEnumerationAlgorithm,
        mut visitor: impl FnMut(&[AtomId]) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        let mut scratch: Vec<AtomId> = Vec::new();
        self.graph.visit_simple_cycles(max_size, alg, |cycle| {
            scratch.clear();
            scratch.extend(cycle.nodes().iter().copied().map(AtomId::from));
            visitor(&scratch)
        })
    }

    pub fn enumerate_simple_cycles(
        &self,
        max_size: usize,
        alg: SimpleCycleEnumerationAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        let mut cycles = Vec::new();
        let _: ControlFlow<()> = self.visit_simple_cycles(max_size, alg, |cycle| {
            cycles.push(cycle.to_vec());
            ControlFlow::Continue(())
        });
        cycles
    }

    /// Visits each relevant cycle's atoms until traversal completes or the visitor
    /// returns [`ControlFlow::Break`]. The slice borrows visit state and is only
    /// valid for the duration of the call.
    pub fn visit_relevant_cycles<B>(
        &self,
        max_size: usize,
        alg: RelevantCycleEnumerationAlgorithm,
        mut visitor: impl FnMut(&[AtomId]) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        let mut scratch: Vec<AtomId> = Vec::new();
        self.graph.visit_relevant_cycles(max_size, alg, |cycle| {
            scratch.clear();
            scratch.extend(cycle.nodes().iter().copied().map(AtomId::from));
            visitor(&scratch)
        })
    }

    pub fn enumerate_relevant_cycles(
        &self,
        max_size: usize,
        alg: RelevantCycleEnumerationAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        let mut cycles = Vec::new();
        let _: ControlFlow<()> = self.visit_relevant_cycles(max_size, alg, |cycle| {
            cycles.push(cycle.to_vec());
            ControlFlow::Continue(())
        });
        cycles
    }

    pub fn maximum_independent_set(&self, alg: MaximumIndependentSetAlgorithm) -> Vec<AtomId> {
        self.graph
            .maximum_independent_set(alg)
            .into_iter()
            .map(AtomId::from)
            .collect()
    }

    pub fn bipartite_maximum_matching(
        &self,
        node_order: &[AtomId],
        algorithm: BipartiteMaximumMatchingAlgorithm,
    ) -> Result<BondMatching, NonBipartiteGraphError> {
        let nodes: Vec<NodeId> = node_order.iter().copied().map(NodeId::from).collect();
        self.graph
            .bipartite_maximum_matching(&nodes, algorithm)
            .map(BondMatching)
    }

    pub fn general_maximum_matching(
        &self,
        node_order: &[AtomId],
        algorithm: GeneralMaximumMatchingAlgorithm,
    ) -> BondMatching {
        let nodes: Vec<NodeId> = node_order.iter().copied().map(NodeId::from).collect();
        BondMatching(self.graph.general_maximum_matching(&nodes, algorithm))
    }

    pub fn bipartite_maximum_matching_or_general(
        &self,
        node_order: &[AtomId],
        bipartite_algorithm: BipartiteMaximumMatchingAlgorithm,
        general_algorithm: GeneralMaximumMatchingAlgorithm,
    ) -> BondMatching {
        let nodes: Vec<NodeId> = node_order.iter().copied().map(NodeId::from).collect();
        BondMatching(self.graph.bipartite_maximum_matching_or_general(
            &nodes,
            bipartite_algorithm,
            general_algorithm,
        ))
    }

    pub fn perfect_matching(
        &self,
        node_order: &[AtomId],
        alg: PerfectMatchingAlgorithm,
    ) -> Option<BondMatching> {
        let nodes: Vec<NodeId> = node_order.iter().copied().map(NodeId::from).collect();
        self.graph.perfect_matching(&nodes, alg).map(BondMatching)
    }

    /// Visits every perfect matching until traversal completes or the visitor
    /// returns [`ControlFlow::Break`].
    pub fn visit_perfect_matchings<B>(
        &self,
        alg: MatchingEnumerationAlgorithm,
        mut visitor: impl FnMut(BondMatching) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        self.graph
            .visit_perfect_matchings(alg, |matching| visitor(BondMatching(matching)))
    }

    pub fn enumerate_perfect_matchings(
        &self,
        alg: MatchingEnumerationAlgorithm,
    ) -> Vec<BondMatching> {
        let mut matchings = Vec::new();
        let _: ControlFlow<()> = self.visit_perfect_matchings(alg, |matching| {
            matchings.push(matching);
            ControlFlow::Continue(())
        });
        matchings
    }

    /// Visits every maximum matching until traversal completes or the visitor
    /// returns [`ControlFlow::Break`].
    pub fn visit_maximum_matchings<B>(
        &self,
        alg: MatchingEnumerationAlgorithm,
        mut visitor: impl FnMut(BondMatching) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        self.graph
            .visit_maximum_matchings(alg, |matching| visitor(BondMatching(matching)))
    }

    pub fn enumerate_maximum_matchings(
        &self,
        alg: MatchingEnumerationAlgorithm,
    ) -> Vec<BondMatching> {
        let mut matchings = Vec::new();
        let _: ControlFlow<()> = self.visit_maximum_matchings(alg, |matching| {
            matchings.push(matching);
            ControlFlow::Continue(())
        });
        matchings
    }

    pub fn automorphisms<C: Ord + Copy>(
        &self,
        atom_color: impl Fn(AtomId) -> C,
        alg: AutomorphismAlgorithm,
    ) -> AtomAutomorphism {
        AtomAutomorphism(
            self.graph
                .automorphisms(|n| atom_color(AtomId::from(n)), alg),
        )
    }

    /// Visits each occurrence of `query` within `self` (the host) as a
    /// query→host `AtomId` embedding (query atom index → host atom) until
    /// traversal completes or the visitor returns [`ControlFlow::Break`]. The
    /// slice borrows visit state and is only valid for the duration of the
    /// call. `atom_match`/`bond_match` receive `(query, host)` — i.e.
    /// `(pattern, target)`, matching the `pattern.matches(target)` convention.
    pub fn visit_subgraph_isomorphisms<B>(
        &self,
        query: &GraphView<'_>,
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
        mut visitor: impl FnMut(&[AtomId]) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        let mut scratch: Vec<AtomId> = Vec::new();
        self.graph.visit_subgraph_isomorphisms(
            query.graph,
            &mut |query_node, host_node| {
                atom_match(AtomId::from(query_node), AtomId::from(host_node))
            },
            &mut |query_edge, host_edge| {
                bond_match(BondId::from(query_edge), BondId::from(host_edge))
            },
            alg,
            |embedding| {
                scratch.clear();
                scratch.extend(embedding.iter().copied().map(AtomId::from));
                visitor(&scratch)
            },
        )
    }

    /// Occurrences of `query` within `self` (the host), one query→host `AtomId`
    /// vector per occurrence, collected from
    /// [`visit_subgraph_isomorphisms`](Self::visit_subgraph_isomorphisms).
    pub fn enumerate_subgraph_isomorphisms(
        &self,
        query: &GraphView<'_>,
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        let mut occurrences = Vec::new();
        let _: ControlFlow<()> =
            self.visit_subgraph_isomorphisms(query, atom_match, bond_match, alg, |embedding| {
                occurrences.push(embedding.to_vec());
                ControlFlow::Continue(())
            });
        occurrences
    }

    /// Like [`visit_subgraph_isomorphisms`](Self::visit_subgraph_isomorphisms) with query atom
    /// `anchor.0` pinned to host atom `anchor.1`. Closures receive `(query, host)`.
    pub fn visit_subgraph_isomorphisms_at<B>(
        &self,
        query: &GraphView<'_>,
        anchor: (AtomId, AtomId),
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
        mut visitor: impl FnMut(&[AtomId]) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        let mut scratch: Vec<AtomId> = Vec::new();
        self.graph.visit_subgraph_isomorphisms_at(
            query.graph,
            (NodeId::from(anchor.0), NodeId::from(anchor.1)),
            &mut |query_node, host_node| {
                atom_match(AtomId::from(query_node), AtomId::from(host_node))
            },
            &mut |query_edge, host_edge| {
                bond_match(BondId::from(query_edge), BondId::from(host_edge))
            },
            alg,
            |embedding| {
                scratch.clear();
                scratch.extend(embedding.iter().copied().map(AtomId::from));
                visitor(&scratch)
            },
        )
    }

    /// Like [`enumerate_subgraph_isomorphisms`](Self::enumerate_subgraph_isomorphisms) with query atom
    /// `anchor.0` pinned to host atom `anchor.1`, collected from
    /// [`visit_subgraph_isomorphisms_at`](Self::visit_subgraph_isomorphisms_at).
    pub fn enumerate_subgraph_isomorphisms_at(
        &self,
        query: &GraphView<'_>,
        anchor: (AtomId, AtomId),
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        let mut occurrences = Vec::new();
        let _: ControlFlow<()> = self.visit_subgraph_isomorphisms_at(
            query,
            anchor,
            atom_match,
            bond_match,
            alg,
            |embedding| {
                occurrences.push(embedding.to_vec());
                ControlFlow::Continue(())
            },
        );
        occurrences
    }
}

/// Atom-level wrapper over `umol_graph_core::AutomorphismOutput` — the result of
/// a graph automorphism query. Indexes the permutation in terms of `AtomId`
/// rather than raw `NodeId`.
#[derive(Clone, Debug)]
pub struct AtomAutomorphism(AutomorphismOutput);

impl AtomAutomorphism {
    pub fn atom_count(&self) -> usize {
        self.0.node_count()
    }

    pub fn orbit_count(&self) -> usize {
        self.0.orbit_count()
    }

    pub fn orbit_of(&self, atom: AtomId) -> AtomId {
        AtomId::from(self.0.orbit_of(NodeId::from(atom)))
    }

    pub fn same_orbit(&self, a: AtomId, b: AtomId) -> bool {
        self.0.same_orbit(NodeId::from(a), NodeId::from(b))
    }

    pub fn canonical_labels(&self) -> Vec<AtomId> {
        self.0
            .canonical_labels()
            .iter()
            .map(|&n| AtomId::from(n))
            .collect()
    }

    pub fn group_order(&self) -> AutomorphismGroupOrder {
        self.0.group_order()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::id::AtomId;
    use crate::ir::molecule::{Molecule, MoleculeEntries};

    #[fixture]
    fn hexagon() -> AtomAutomorphism {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(4), AtomId(5), BondForm::from_order(1)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
            ],
            ..Default::default()
        })
        .graph()
        .automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty)
    }

    #[fixture]
    fn chain_3() -> AtomAutomorphism {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
            ],
            ..Default::default()
        })
        .graph()
        .automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty)
    }

    #[rstest]
    #[case::hexagon(hexagon(), 6)]
    #[case::chain_3(chain_3(), 3)]
    fn test_atom_automorphism_atom_count(
        #[case] automorphism: AtomAutomorphism,
        #[case] expected: usize,
    ) {
        assert_eq!(automorphism.atom_count(), expected);
    }

    #[rstest]
    #[case::hexagon(hexagon(), 1)]
    #[case::chain_3(chain_3(), 2)]
    fn test_atom_automorphism_orbit_count(
        #[case] automorphism: AtomAutomorphism,
        #[case] expected: usize,
    ) {
        assert_eq!(automorphism.orbit_count(), expected);
    }

    #[rstest]
    #[case::hexagon(hexagon(), AtomId(3), AtomId(0))]
    #[case::chain_endpoint(chain_3(), AtomId(2), AtomId(0))]
    #[case::chain_center(chain_3(), AtomId(1), AtomId(1))]
    fn test_atom_automorphism_orbit_of(
        #[case] automorphism: AtomAutomorphism,
        #[case] atom: AtomId,
        #[case] expected: AtomId,
    ) {
        assert_eq!(automorphism.orbit_of(atom), expected);
    }

    #[rstest]
    #[case::hexagon_equivalent(hexagon(), AtomId(0), AtomId(3), true)]
    #[case::chain_endpoints_equivalent(chain_3(), AtomId(0), AtomId(2), true)]
    #[case::chain_endpoint_vs_middle(chain_3(), AtomId(0), AtomId(1), false)]
    fn test_atom_automorphism_same_orbit(
        #[case] automorphism: AtomAutomorphism,
        #[case] first: AtomId,
        #[case] second: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(automorphism.same_orbit(first, second), expected);
    }

    #[rstest]
    #[case::hexagon(hexagon(), (0..6).map(|index| AtomId(index as u32)).collect::<Vec<_>>())]
    fn test_atom_automorphism_canonical_labels(
        #[case] automorphism: AtomAutomorphism,
        #[case] expected: Vec<AtomId>,
    ) {
        let mut sorted = automorphism.canonical_labels();
        sorted.sort_unstable();
        assert_eq!(sorted, expected);
    }

    #[rstest]
    #[case::hexagon(hexagon(), AutomorphismGroupOrder::Exact(12))]
    #[case::chain_3(chain_3(), AutomorphismGroupOrder::Exact(2))]
    fn test_atom_automorphism_group_order(
        #[case] automorphism: AtomAutomorphism,
        #[case] expected: AutomorphismGroupOrder,
    ) {
        assert_eq!(automorphism.group_order(), expected);
    }
}
