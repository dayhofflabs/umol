//! Graph view: typed adapter over the underlying `Graph`.

use std::ops::ControlFlow;

use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismGroup, AutomorphismGroupOrder,
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

    /// Compute the atom automorphism group of this graph under `atom_color` and
    /// `bond_color` using `alg`.
    ///
    /// Skips canonical-labeling construction when only symmetry orbits, group
    /// order, and generators are needed.
    ///
    /// # Semantic properties
    ///
    /// Orbits, orbit count, group order, and the generated subgroup depend only
    /// on the graph and the equivalence classes induced by `atom_color` and
    /// `bond_color`, and agree across [`AutomorphismAlgorithm`] variants.
    pub fn automorphisms<AC: Ord, BC: Ord>(
        &self,
        atom_color: impl Fn(AtomId) -> AC,
        bond_color: impl Fn(BondId) -> BC,
        alg: AutomorphismAlgorithm,
    ) -> AtomAutomorphismGroup {
        AtomAutomorphismGroup(self.graph.automorphisms(
            |node| atom_color(AtomId::from(node)),
            |edge| bond_color(BondId::from(edge)),
            alg,
        ))
    }

    /// Compute a canonical atom ordering and automorphism group of this graph
    /// under `atom_color` and `bond_color` using `alg`.
    ///
    /// # Semantic properties
    ///
    /// For a fixed [`AutomorphismAlgorithm`] and color ordering, the returned
    /// labeling induces an identical canonical colored adjacency matrix across
    /// all atom renumberings of isomorphic graphs.
    pub fn canonical_labeling<AC: Ord, BC: Ord>(
        &self,
        atom_color: impl Fn(AtomId) -> AC,
        bond_color: impl Fn(BondId) -> BC,
        alg: AutomorphismAlgorithm,
    ) -> AtomCanonicalLabeling {
        let labeling = self.graph.canonical_labeling(
            |node| atom_color(AtomId::from(node)),
            |edge| bond_color(BondId::from(edge)),
            alg,
        );
        let labels = labeling
            .labels()
            .iter()
            .copied()
            .map(AtomId::from)
            .collect();
        AtomCanonicalLabeling {
            labels,
            automorphisms: AtomAutomorphismGroup(labeling.into_automorphisms()),
        }
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

    /// Returns whether at least one subgraph isomorphism of `query` exists within `self` (the
    /// host). `atom_match`/`bond_match` receive `(query, host)`.
    pub fn has_subgraph_isomorphism(
        &self,
        query: &GraphView<'_>,
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> bool {
        self.graph.has_subgraph_isomorphism(
            query.graph,
            &mut |query_node, host_node| {
                atom_match(AtomId::from(query_node), AtomId::from(host_node))
            },
            &mut |query_edge, host_edge| {
                bond_match(BondId::from(query_edge), BondId::from(host_edge))
            },
            alg,
        )
    }

    /// Returns the first subgraph isomorphism of `query` within `self` (the host) as a
    /// query→host `AtomId` vector, or `None` if no match exists. `atom_match`/`bond_match`
    /// receive `(query, host)`.
    pub fn first_subgraph_isomorphism(
        &self,
        query: &GraphView<'_>,
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Option<Vec<AtomId>> {
        match self.graph.visit_subgraph_isomorphisms(
            query.graph,
            &mut |query_node, host_node| {
                atom_match(AtomId::from(query_node), AtomId::from(host_node))
            },
            &mut |query_edge, host_edge| {
                bond_match(BondId::from(query_edge), BondId::from(host_edge))
            },
            alg,
            |embedding| ControlFlow::Break(embedding.iter().copied().map(AtomId::from).collect()),
        ) {
            ControlFlow::Break(embedding) => Some(embedding),
            ControlFlow::Continue(()) => None,
        }
    }

    /// Counts every subgraph isomorphism of `query` within `self` (the host) without
    /// materializing embeddings. `atom_match`/`bond_match` receive `(query, host)`.
    pub fn count_subgraph_isomorphisms(
        &self,
        query: &GraphView<'_>,
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> u64 {
        self.graph.count_subgraph_isomorphisms(
            query.graph,
            &mut |query_node, host_node| {
                atom_match(AtomId::from(query_node), AtomId::from(host_node))
            },
            &mut |query_edge, host_edge| {
                bond_match(BondId::from(query_edge), BondId::from(host_edge))
            },
            alg,
        )
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

    /// Returns whether at least one subgraph isomorphism mapping `anchor.0` to `anchor.1` exists
    /// within `self` (the host). `atom_match`/`bond_match` receive `(query, host)`.
    pub fn has_subgraph_isomorphism_at(
        &self,
        query: &GraphView<'_>,
        anchor: (AtomId, AtomId),
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> bool {
        self.graph.has_subgraph_isomorphism_at(
            query.graph,
            (NodeId::from(anchor.0), NodeId::from(anchor.1)),
            &mut |query_node, host_node| {
                atom_match(AtomId::from(query_node), AtomId::from(host_node))
            },
            &mut |query_edge, host_edge| {
                bond_match(BondId::from(query_edge), BondId::from(host_edge))
            },
            alg,
        )
    }

    /// Returns the first subgraph isomorphism mapping `anchor.0` to `anchor.1` within `self` (the
    /// host) as a query→host `AtomId` vector, or `None` if no anchored match exists.
    /// `atom_match`/`bond_match` receive `(query, host)`.
    pub fn first_subgraph_isomorphism_at(
        &self,
        query: &GraphView<'_>,
        anchor: (AtomId, AtomId),
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Option<Vec<AtomId>> {
        match self.graph.visit_subgraph_isomorphisms_at(
            query.graph,
            (NodeId::from(anchor.0), NodeId::from(anchor.1)),
            &mut |query_node, host_node| {
                atom_match(AtomId::from(query_node), AtomId::from(host_node))
            },
            &mut |query_edge, host_edge| {
                bond_match(BondId::from(query_edge), BondId::from(host_edge))
            },
            alg,
            |embedding| ControlFlow::Break(embedding.iter().copied().map(AtomId::from).collect()),
        ) {
            ControlFlow::Break(embedding) => Some(embedding),
            ControlFlow::Continue(()) => None,
        }
    }

    /// Counts every subgraph isomorphism mapping `anchor.0` to `anchor.1` within `self` (the
    /// host) without materializing embeddings. `atom_match`/`bond_match` receive `(query, host)`.
    pub fn count_subgraph_isomorphisms_at(
        &self,
        query: &GraphView<'_>,
        anchor: (AtomId, AtomId),
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> u64 {
        self.graph.count_subgraph_isomorphisms_at(
            query.graph,
            (NodeId::from(anchor.0), NodeId::from(anchor.1)),
            &mut |query_node, host_node| {
                atom_match(AtomId::from(query_node), AtomId::from(host_node))
            },
            &mut |query_edge, host_edge| {
                bond_match(BondId::from(query_edge), BondId::from(host_edge))
            },
            alg,
        )
    }
}

/// Atom-level symmetry group wrapping [`AutomorphismGroup`], issued by
/// `GraphView::automorphisms` or [`AtomCanonicalLabeling`].
///
/// Every instance preserves the invariants of [`AutomorphismGroup`] in the
/// molecule's [`AtomId`] space: each orbit representative is the minimum
/// [`AtomId`] in its orbit, [`Self::orbit_count`] equals the number of distinct
/// orbits, and every permutation in [`Self::generators`] is an atom bijection
/// preserving adjacency, atom colors, and bond colors.
///
/// Queries taking an [`AtomId`] assume the identifier belongs to the molecule
/// that produced this value.
///
/// # Semantic properties
///
/// - Orbits, [`Self::orbit_count`], [`Self::group_order`], and the generated
///   subgroup are invariant across [`AutomorphismAlgorithm`] variants.
/// - The specific generating basis returned by [`Self::generators`] is
///   algorithm-dependent.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomAutomorphismGroup(AutomorphismGroup);

impl AtomAutomorphismGroup {
    /// Number of atoms in the graph that produced this group.
    pub fn atom_count(&self) -> usize {
        self.0.node_count()
    }

    /// Number of distinct atom orbits.
    pub fn orbit_count(&self) -> usize {
        self.0.orbit_count()
    }

    /// Minimum [`AtomId`] in the orbit containing `atom`.
    ///
    /// # Panics
    ///
    /// Panics if `atom.0 as usize >= self.atom_count()`.
    pub fn orbit_of(&self, atom: AtomId) -> AtomId {
        AtomId::from(self.0.orbit_of(NodeId::from(atom)))
    }

    /// Whether `first` and `second` lie in the same automorphism orbit.
    ///
    /// # Panics
    ///
    /// Panics if `first.0 as usize >= self.atom_count()` or
    /// `second.0 as usize >= self.atom_count()`.
    pub fn same_orbit(&self, first: AtomId, second: AtomId) -> bool {
        self.0.same_orbit(NodeId::from(first), NodeId::from(second))
    }

    /// Order of the automorphism group.
    pub fn group_order(&self) -> AutomorphismGroupOrder {
        self.0.group_order()
    }

    /// A generating set of the automorphism group, each generator stored as a
    /// permutation image over `0..atom_count`.
    pub fn generators(&self) -> Vec<Vec<AtomId>> {
        self.0
            .generators()
            .iter()
            .map(|generator| generator.iter().copied().map(AtomId::from).collect())
            .collect()
    }
}

/// Canonical atom ordering paired with [`AtomAutomorphismGroup`], issued by
/// `GraphView::canonical_labeling`.
///
/// [`Self::labels`] has length `self.automorphisms().atom_count()` and is a
/// bijection over `0..atom_count`, where `labels()[rank]` is the source
/// [`AtomId`] assigned canonical position `rank`.
///
/// # Semantic properties
///
/// - For a fixed [`AutomorphismAlgorithm`] and color ordering, [`Self::labels`]
///   induces an identical canonical colored adjacency matrix across all atom
///   renumberings of isomorphic graphs.
/// - Tie-breaking between equivalent canonical orderings is
///   algorithm-dependent.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomCanonicalLabeling {
    labels: Vec<AtomId>,
    automorphisms: AtomAutomorphismGroup,
}

impl AtomCanonicalLabeling {
    /// Canonical atom order: `labels()[rank]` is the source [`AtomId`] assigned
    /// canonical position `rank`.
    pub fn labels(&self) -> &[AtomId] {
        &self.labels
    }

    /// Atom automorphism group computed alongside this canonical labeling.
    pub fn automorphisms(&self) -> &AtomAutomorphismGroup {
        &self.automorphisms
    }

    /// Consume this labeling and return its [`AtomAutomorphismGroup`].
    pub fn into_automorphisms(self) -> AtomAutomorphismGroup {
        self.automorphisms
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use AutomorphismAlgorithm::{Nauty, TrailIr};

    use super::*;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::id::AtomId;
    use crate::ir::molecule::{Molecule, MoleculeEntries};

    #[fixture]
    fn hexagon() -> Molecule {
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
    }

    #[fixture]
    fn cyclohexene() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(4), AtomId(5), BondForm::from_order(1)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
            ],
            ..Default::default()
        })
    }

    #[fixture]
    fn chain_3() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
            ],
            ..Default::default()
        })
    }

    fn bond_order_key(molecule: &Molecule, bond: BondId) -> u8 {
        if molecule.bond(bond).attributes().order == BondForm::from_order(2).order {
            2
        } else {
            1
        }
    }

    #[rstest]
    #[case::hexagon(hexagon(), 6)]
    #[case::cyclohexene(cyclohexene(), 6)]
    #[case::chain_3(chain_3(), 3)]
    fn test_atom_automorphism_group_atom_count(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] molecule: Molecule,
        #[case] expected: usize,
    ) {
        let group =
            molecule
                .graph()
                .automorphisms(|_| 0u8, |bond| bond_order_key(&molecule, bond), alg);
        assert_eq!(group.atom_count(), expected);
    }

    #[rstest]
    #[case::hexagon(hexagon(), 1)]
    #[case::cyclohexene(cyclohexene(), 3)]
    #[case::chain_3(chain_3(), 2)]
    fn test_atom_automorphism_group_orbit_count(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] molecule: Molecule,
        #[case] expected: usize,
    ) {
        let group =
            molecule
                .graph()
                .automorphisms(|_| 0u8, |bond| bond_order_key(&molecule, bond), alg);
        assert_eq!(group.orbit_count(), expected);
    }

    #[rstest]
    #[case::hexagon(hexagon(), AtomId(3), AtomId(0))]
    #[case::cyclohexene_alkene(cyclohexene(), AtomId(1), AtomId(0))]
    #[case::cyclohexene_allylic(cyclohexene(), AtomId(5), AtomId(2))]
    #[case::cyclohexene_homoallylic(cyclohexene(), AtomId(4), AtomId(3))]
    #[case::chain_endpoint(chain_3(), AtomId(2), AtomId(0))]
    #[case::chain_center(chain_3(), AtomId(1), AtomId(1))]
    fn test_atom_automorphism_group_orbit_of(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: AtomId,
    ) {
        let group =
            molecule
                .graph()
                .automorphisms(|_| 0u8, |bond| bond_order_key(&molecule, bond), alg);
        assert_eq!(group.orbit_of(atom), expected);
    }

    #[rstest]
    #[case::hexagon_equivalent(hexagon(), AtomId(0), AtomId(3), true)]
    #[case::cyclohexene_alkene_equivalent(cyclohexene(), AtomId(0), AtomId(1), true)]
    #[case::cyclohexene_alkene_vs_allylic(cyclohexene(), AtomId(0), AtomId(2), false)]
    #[case::chain_endpoints_equivalent(chain_3(), AtomId(0), AtomId(2), true)]
    #[case::chain_endpoint_vs_middle(chain_3(), AtomId(0), AtomId(1), false)]
    fn test_atom_automorphism_group_same_orbit(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] molecule: Molecule,
        #[case] first: AtomId,
        #[case] second: AtomId,
        #[case] expected: bool,
    ) {
        let group =
            molecule
                .graph()
                .automorphisms(|_| 0u8, |bond| bond_order_key(&molecule, bond), alg);
        assert_eq!(group.same_orbit(first, second), expected);
    }

    #[rstest]
    #[case::hexagon(hexagon(), AutomorphismGroupOrder::Exact(12))]
    #[case::cyclohexene(cyclohexene(), AutomorphismGroupOrder::Exact(2))]
    #[case::chain_3(chain_3(), AutomorphismGroupOrder::Exact(2))]
    fn test_atom_automorphism_group_group_order(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] molecule: Molecule,
        #[case] expected: AutomorphismGroupOrder,
    ) {
        let group =
            molecule
                .graph()
                .automorphisms(|_| 0u8, |bond| bond_order_key(&molecule, bond), alg);
        assert_eq!(group.group_order(), expected);
    }

    #[rstest]
    #[case::chain_3(chain_3(), vec![vec![AtomId(2), AtomId(1), AtomId(0)]])]
    #[case::cyclohexene(
        cyclohexene(),
        vec![vec![AtomId(1), AtomId(0), AtomId(5), AtomId(4), AtomId(3), AtomId(2)]]
    )]
    fn test_atom_automorphism_group_generators(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] molecule: Molecule,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        let group =
            molecule
                .graph()
                .automorphisms(|_| 0u8, |bond| bond_order_key(&molecule, bond), alg);
        assert_eq!(group.generators(), expected);
    }

    #[rstest]
    #[should_panic]
    #[case::orbit_of_oob(AtomId(3), None)]
    #[should_panic]
    #[case::same_orbit_first_oob(AtomId(3), Some(AtomId(0)))]
    #[should_panic]
    #[case::same_orbit_second_oob(AtomId(0), Some(AtomId(3)))]
    fn test_atom_automorphism_group_out_of_bounds_panic(
        #[case] first: AtomId,
        #[case] second: Option<AtomId>,
    ) {
        let group = chain_3().graph().automorphisms(|_| 0u8, |_| 0u8, TrailIr);
        match second {
            None => {
                let _ = group.orbit_of(first);
            }
            Some(second) => {
                let _ = group.same_orbit(first, second);
            }
        }
    }

    #[rstest]
    #[case::hexagon(hexagon(), 1, AutomorphismGroupOrder::Exact(12))]
    #[case::cyclohexene(cyclohexene(), 3, AutomorphismGroupOrder::Exact(2))]
    #[case::chain_3(chain_3(), 2, AutomorphismGroupOrder::Exact(2))]
    fn test_atom_canonical_labeling(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] molecule: Molecule,
        #[case] expected_orbits: usize,
        #[case] expected_order: AutomorphismGroupOrder,
    ) {
        let labeling = molecule.graph().canonical_labeling(
            |_| 0u8,
            |bond| bond_order_key(&molecule, bond),
            alg,
        );
        let mut sorted = labeling.labels().to_vec();
        sorted.sort_unstable();
        let expected_atoms: Vec<AtomId> = molecule.atoms().ids().collect();
        assert_eq!(sorted, expected_atoms);
        assert_eq!(labeling.automorphisms().orbit_count(), expected_orbits);
        assert_eq!(labeling.automorphisms().group_order(), expected_order);

        let expected_group = labeling.automorphisms().clone();
        assert_eq!(labeling.into_automorphisms(), expected_group);
    }

    #[rstest]
    #[case::hexagon_chain_3(hexagon(), chain_3(), (AtomId(1), AtomId(0)), 12, 2)]
    #[case::cyclohexene_chain_3_bond_orders(cyclohexene(), chain_3(), (AtomId(1), AtomId(3)), 8, 2)]
    #[case::cyclohexene_chain_3_anchored_double_bond_endpoint(cyclohexene(), chain_3(), (AtomId(1), AtomId(0)), 8, 0)]
    #[case::chain_3_in_hexagon_no_match(chain_3(), hexagon(), (AtomId(0), AtomId(0)), 0, 0)]
    fn test_graph_view_subgraph_isomorphism_modes(
        #[values(
            SubgraphIsomorphismAlgorithm::Vf2Rdkit,
            SubgraphIsomorphismAlgorithm::RayKirsch,
            SubgraphIsomorphismAlgorithm::BitsetRi,
            SubgraphIsomorphismAlgorithm::Glasgow
        )]
        alg: SubgraphIsomorphismAlgorithm,
        #[case] host: Molecule,
        #[case] query: Molecule,
        #[case] anchor: (AtomId, AtomId),
        #[case] expected_count: u64,
        #[case] expected_anchored_count: u64,
    ) {
        let host_view = host.graph();
        let query_view = query.graph();

        let all = host_view.enumerate_subgraph_isomorphisms(
            &query_view,
            &mut |_, _| true,
            &mut |qb, hb| bond_order_key(&query, qb) == bond_order_key(&host, hb),
            alg,
        );
        let count = host_view.count_subgraph_isomorphisms(
            &query_view,
            &mut |_, _| true,
            &mut |qb, hb| bond_order_key(&query, qb) == bond_order_key(&host, hb),
            alg,
        );
        let has = host_view.has_subgraph_isomorphism(
            &query_view,
            &mut |_, _| true,
            &mut |qb, hb| bond_order_key(&query, qb) == bond_order_key(&host, hb),
            alg,
        );
        let first = host_view.first_subgraph_isomorphism(
            &query_view,
            &mut |_, _| true,
            &mut |qb, hb| bond_order_key(&query, qb) == bond_order_key(&host, hb),
            alg,
        );
        assert_eq!(all.len() as u64, expected_count);
        assert_eq!(count, expected_count);
        assert_eq!(has, expected_count > 0);
        assert_eq!(first.is_some(), expected_count > 0);
        if let Some(first_emb) = first {
            assert!(all.contains(&first_emb));
        }

        let all_at = host_view.enumerate_subgraph_isomorphisms_at(
            &query_view,
            anchor,
            &mut |_, _| true,
            &mut |qb, hb| bond_order_key(&query, qb) == bond_order_key(&host, hb),
            alg,
        );
        let count_at = host_view.count_subgraph_isomorphisms_at(
            &query_view,
            anchor,
            &mut |_, _| true,
            &mut |qb, hb| bond_order_key(&query, qb) == bond_order_key(&host, hb),
            alg,
        );
        let has_at = host_view.has_subgraph_isomorphism_at(
            &query_view,
            anchor,
            &mut |_, _| true,
            &mut |qb, hb| bond_order_key(&query, qb) == bond_order_key(&host, hb),
            alg,
        );
        let first_at = host_view.first_subgraph_isomorphism_at(
            &query_view,
            anchor,
            &mut |_, _| true,
            &mut |qb, hb| bond_order_key(&query, qb) == bond_order_key(&host, hb),
            alg,
        );
        assert_eq!(all_at.len() as u64, expected_anchored_count);
        assert_eq!(count_at, expected_anchored_count);
        assert_eq!(has_at, expected_anchored_count > 0);
        assert_eq!(first_at.is_some(), expected_anchored_count > 0);
        if let Some(first_emb) = first_at {
            assert_eq!(first_emb[anchor.0.index()], anchor.1);
            assert!(all_at.contains(&first_emb));
        }
    }
}
