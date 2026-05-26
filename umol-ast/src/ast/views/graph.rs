//! Graph view: typed adapter over the underlying `Graph`.

use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm, EdgeId, Graph, MatchingEnumerationAlgorithm,
    MaxIndependentSetAlgorithm, MaxMatchingAlgorithm, NodeId, PerfectMatchingAlgorithm,
    ShortestCycleAlgorithm, SubgraphIsomorphismAlgorithm,
};

use super::super::automorphism::AtomAutomorphism;
use super::super::ids::{AtomId, BondId};
use super::super::matching::BondMatching;

/// AtomId/BondId-typed adapter over the underlying `Graph`. Holds the
/// pure-graph algorithms (connectivity, cycles, matchings, isomorphisms)
/// without exposing graph-core's `NodeId` / `EdgeId` types in the public
/// API. Construct via `MoleculeAst::graph()`.
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

    pub fn connected_components(&self, alg: ConnectedComponentsAlgorithm) -> Vec<Vec<AtomId>> {
        self.graph
            .connected_components(alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn biconnected_components(&self, alg: BiconnectedComponentsAlgorithm) -> Vec<Vec<AtomId>> {
        self.graph
            .biconnected_components(alg)
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

    pub fn enumerate_cycles(
        &self,
        max_size: usize,
        alg: CycleEnumerationAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .enumerate_cycles(max_size, alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn maximum_independent_set(&self, alg: MaxIndependentSetAlgorithm) -> Vec<AtomId> {
        self.graph
            .maximum_independent_set(alg)
            .into_iter()
            .map(AtomId::from)
            .collect()
    }

    pub fn maximum_matching(&self, alg: MaxMatchingAlgorithm) -> BondMatching {
        BondMatching(self.graph.maximum_matching(alg))
    }

    pub fn perfect_matching(
        &self,
        node_order: &[AtomId],
        alg: PerfectMatchingAlgorithm,
    ) -> Option<BondMatching> {
        let nodes: Vec<NodeId> = node_order.iter().copied().map(NodeId::from).collect();
        self.graph.perfect_matching(&nodes, alg).map(BondMatching)
    }

    pub fn enumerate_perfect_matchings(
        &self,
        alg: MatchingEnumerationAlgorithm,
    ) -> Vec<BondMatching> {
        self.graph
            .enumerate_perfect_matchings(alg)
            .into_iter()
            .map(BondMatching)
            .collect()
    }

    pub fn enumerate_maximum_matchings(
        &self,
        alg: MatchingEnumerationAlgorithm,
    ) -> Vec<BondMatching> {
        self.graph
            .enumerate_maximum_matchings(alg)
            .into_iter()
            .map(BondMatching)
            .collect()
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

    pub fn subgraph_isomorphisms(
        &self,
        query: &GraphView<'_>,
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .subgraph_isomorphisms(
                query.graph,
                &mut |tn, qn| atom_match(AtomId::from(tn), AtomId::from(qn)),
                &mut |te, qe| bond_match(BondId::from(te), BondId::from(qe)),
                alg,
            )
            .into_iter()
            .map(|m| m.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn subgraph_isomorphisms_at(
        &self,
        query: &GraphView<'_>,
        anchor: (AtomId, AtomId),
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .subgraph_isomorphisms_at(
                query.graph,
                (NodeId::from(anchor.0), NodeId::from(anchor.1)),
                &mut |tn, qn| atom_match(AtomId::from(tn), AtomId::from(qn)),
                &mut |te, qe| bond_match(BondId::from(te), BondId::from(qe)),
                alg,
            )
            .into_iter()
            .map(|m| m.into_iter().map(AtomId::from).collect())
            .collect()
    }
}
