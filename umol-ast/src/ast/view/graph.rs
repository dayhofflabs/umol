//! Graph view: typed adapter over the underlying `Graph`.

use umol_graph_core::{
    AutoGroupOrder, Automorphism, AutomorphismAlgorithm, BiconnectedComponentsAlgorithm,
    ConnectedComponentsAlgorithm, CycleEnumerationAlgorithm, EdgeId, Graph,
    MatchingEnumerationAlgorithm, MaxIndependentSetAlgorithm, MaxMatchingAlgorithm, NodeId,
    PerfectMatchingAlgorithm, ShortestCycleAlgorithm, SubgraphIsomorphismAlgorithm,
};

use super::super::id::{AtomId, BondId};
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

/// Atom-level wrapper over `umol_graph_core::Automorphism` — the result of
/// [`GraphView::automorphisms`]. Indexes the permutation in terms of `AtomId`
/// rather than raw `NodeId`.
#[derive(Clone, Debug)]
pub struct AtomAutomorphism(Automorphism);

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

    pub fn canonical_labeling(&self) -> Vec<AtomId> {
        self.0
            .canonical_labeling()
            .iter()
            .map(|&n| AtomId::from(n))
            .collect()
    }

    pub fn auto_group_order(&self) -> AutoGroupOrder {
        self.0.auto_group_order()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::id::AtomId;
    use crate::ast::molecule::MoleculeAst;

    fn ring(n: usize) -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); n];
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomId(i as u32),
                    AtomId(((i + 1) % n) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::from_atoms_and_bonds(atoms, bonds)
    }

    fn chain(n: usize) -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); n];
        let bonds: Vec<_> = (0..n.saturating_sub(1))
            .map(|i| {
                (
                    AtomId(i as u32),
                    AtomId((i + 1) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::from_atoms_and_bonds(atoms, bonds)
    }

    #[fixture]
    fn hexagon() -> AtomAutomorphism {
        ring(6)
            .graph()
            .automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty)
    }

    #[fixture]
    fn chain_3() -> AtomAutomorphism {
        chain(3)
            .graph()
            .automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty)
    }

    #[rstest]
    #[case::hexagon_wraps_six(hexagon(), 6, 1)]
    #[case::chain_3_wraps_three(chain_3(), 3, 2)]
    fn test_atom_automorphism_counts(
        #[case] auto: AtomAutomorphism,
        #[case] expected_atom_count: usize,
        #[case] expected_orbit_count: usize,
    ) {
        assert_eq!(auto.atom_count(), expected_atom_count);
        assert_eq!(auto.orbit_count(), expected_orbit_count);
    }

    #[rstest]
    fn test_atom_automorphism_orbit_of_agrees_within_orbit(hexagon: AtomAutomorphism) {
        // In a hexagon every atom is in the same orbit, so `orbit_of` must
        // return the same representative for every atom.
        let rep = hexagon.orbit_of(AtomId(0));
        for i in 1..6 {
            assert_eq!(hexagon.orbit_of(AtomId(i)), rep);
        }
    }

    #[rstest]
    fn test_atom_automorphism_orbit_of_distinguishes_non_equivalent(chain_3: AtomAutomorphism) {
        // In C-C-C the center atom sits in a different orbit from the endpoints.
        assert_ne!(chain_3.orbit_of(AtomId(0)), chain_3.orbit_of(AtomId(1)));
        assert_eq!(chain_3.orbit_of(AtomId(0)), chain_3.orbit_of(AtomId(2)));
    }

    #[rstest]
    #[case::hexagon_equivalent(hexagon(), AtomId(0), AtomId(3), true)]
    #[case::chain_endpoints_equivalent(chain_3(), AtomId(0), AtomId(2), true)]
    #[case::chain_endpoint_vs_middle(chain_3(), AtomId(0), AtomId(1), false)]
    fn test_atom_automorphism_same_orbit(
        #[case] auto: AtomAutomorphism,
        #[case] a: AtomId,
        #[case] b: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(auto.same_orbit(a, b), expected);
    }

    #[rstest]
    fn test_atom_automorphism_canonical_labeling_is_permutation(hexagon: AtomAutomorphism) {
        let labeling = hexagon.canonical_labeling();
        assert_eq!(labeling.len(), 6);
        let mut sorted = labeling.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..6).map(|i| AtomId(i as u32)).collect::<Vec<_>>());
    }

    #[rstest]
    fn test_atom_automorphism_auto_group_order_hexagon(hexagon: AtomAutomorphism) {
        // D6 has 12 automorphisms on the vertex set.
        assert_eq!(hexagon.auto_group_order(), AutoGroupOrder::Exact(12));
    }

    #[rstest]
    fn test_atom_automorphism_auto_group_order_chain(chain_3: AtomAutomorphism) {
        // Linear C-C-C: only identity and the endpoint-swap reflection.
        assert_eq!(chain_3.auto_group_order(), AutoGroupOrder::Exact(2));
    }
}
