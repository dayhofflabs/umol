//! Cycle queries, enumeration, bases, and ring-family decomposition.
//!
//! Current operations provide breadth-first shortest-cycle queries; bounded
//! Read--Tarjan simple-cycle visitation and collection; Vismara relevant-cycle
//! visitation and collection; Horton minimum cycle bases; and Kolodzik Unique
//! Ring Families. Read--Tarjan and Vismara expose strict simple-graph,
//! explicit multigraph-fallback, and total dispatch operations. See [Read and
//! Tarjan (1975)](https://doi.org/10.1002/net.1975.5.3.237),
//! [Horton (1987)](https://doi.org/10.1137/0216026),
//! [Vismara (1997)](https://doi.org/10.37236/1294), and
//! [Kolodzik, Urbaczek, and Rarey
//! (2012)](https://doi.org/10.1021/ci200629w).

use std::collections::{HashSet, VecDeque};
use std::error::Error;
use std::fmt;
use std::ops::ControlFlow;
use std::slice::Iter;

use num_bigint::BigUint;

use crate::graph::{EdgeId, Graph, NodeId, SubdividedGraph, SubdivisionNodeSource};

mod basis;
mod relevant;
mod simple;
mod urf;

use self::basis::minimum_cycle_basis_horton;

/// Error returned when a graph operation requires a simple graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NonSimpleGraphError;

impl fmt::Display for NonSimpleGraphError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "graph is not simple")
    }
}

impl Error for NonSimpleGraphError {}

/// An undirected cycle represented by corresponding node and edge sequences.
///
/// Edge `i` connects node `i` to node `(i + 1) % len`. Rotation and reversal
/// do not affect equality because cycles are normalized when constructed.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Cycle {
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
}

impl Cycle {
    fn normalized(graph: &Graph, nodes: Vec<NodeId>, edges: Vec<EdgeId>) -> Self {
        assert!(!nodes.is_empty(), "a cycle must contain a node");
        assert_eq!(
            nodes.len(),
            edges.len(),
            "a cycle must have one edge per node"
        );

        let node_count = nodes.len();
        let distinct_nodes: HashSet<_> = nodes.iter().copied().collect();
        let distinct_edges: HashSet<_> = edges.iter().copied().collect();
        assert_eq!(
            distinct_nodes.len(),
            node_count,
            "cycle nodes must be distinct"
        );
        assert_eq!(
            distinct_edges.len(),
            node_count,
            "cycle edges must be distinct"
        );

        for index in 0..node_count {
            let first = nodes[index];
            let second = nodes[(index + 1) % node_count];
            assert!(graph.contains_node(first), "cycle node is not in the graph");
            assert!(
                graph.contains_node(second),
                "cycle node is not in the graph"
            );
            let edge = edges[index];
            assert!(graph.contains_edge(edge), "cycle edge is not in the graph");
            let [source, target] = graph.edge_endpoints(edge);
            assert!(
                (source == first && target == second) || (source == second && target == first),
                "cycle edge does not connect consecutive nodes"
            );
        }

        let start = nodes
            .iter()
            .enumerate()
            .min_by_key(|(_, node)| *node)
            .map(|(index, _)| index)
            .expect("a non-empty cycle has a minimum node");

        let mut forward_nodes = Vec::with_capacity(node_count);
        let mut forward_edges = Vec::with_capacity(node_count);
        let mut reverse_nodes = Vec::with_capacity(node_count);
        let mut reverse_edges = Vec::with_capacity(node_count);
        for offset in 0..node_count {
            forward_nodes.push(nodes[(start + offset) % node_count]);
            forward_edges.push(edges[(start + offset) % node_count]);
            reverse_nodes.push(nodes[(start + node_count - offset) % node_count]);
            reverse_edges.push(edges[(start + node_count - offset - 1) % node_count]);
        }

        let (nodes, edges) = if (&reverse_nodes, &reverse_edges) < (&forward_nodes, &forward_edges)
        {
            (reverse_nodes, reverse_edges)
        } else {
            (forward_nodes, forward_edges)
        };
        Self { nodes, edges }
    }

    fn map_edges(&self, source: &Graph, edge_sources: &[EdgeId]) -> Self {
        let edges = self
            .edges
            .iter()
            .map(|edge| edge_sources[edge.index()])
            .collect();
        Self::normalized(source, self.nodes.clone(), edges)
    }

    fn map_subdivision(
        &self,
        source: &Graph,
        subdivision: &SubdividedGraph,
        edge_sources: &[EdgeId],
    ) -> Self {
        let mut nodes = Vec::with_capacity(self.length() / 2);
        let mut edges = Vec::with_capacity(self.length() / 2);
        for &node in &self.nodes {
            match subdivision.node_source(node) {
                SubdivisionNodeSource::Node(node) => nodes.push(node),
                SubdivisionNodeSource::Edge(edge) => {
                    edges.push(edge_sources[edge.index()]);
                }
            }
        }
        Self::normalized(source, nodes, edges)
    }

    pub fn nodes(&self) -> &[NodeId] {
        &self.nodes
    }

    pub fn edges(&self) -> &[EdgeId] {
        &self.edges
    }

    pub fn length(&self) -> usize {
        self.nodes.len()
    }
}

/// One minimum-total-length basis of the graph's binary cycle space.
///
/// Cycles use source graph identifiers. For graphs with parallel edges, the
/// internal subdivision does not affect the reported cycle lengths.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MinimumCycleBasis {
    cycles: Vec<Cycle>,
    total_length: usize,
}

impl MinimumCycleBasis {
    /// Number of independent cycles in the basis.
    pub fn dimension(&self) -> usize {
        self.cycles.len()
    }

    /// Sum of basis-cycle lengths in source-edge units.
    pub fn total_length(&self) -> usize {
        self.total_length
    }

    /// Iterate over the selected basis cycles.
    pub fn iter(&self) -> Iter<'_, Cycle> {
        self.cycles.iter()
    }
}

/// Exact number of relevant cycles represented by a Unique Ring Family.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct RelevantCycleCount(pub BigUint);

/// Index of a family within a [`UniqueRingFamilies`] result.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct UniqueRingFamilyId(pub u32);

impl UniqueRingFamilyId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// One equivalence class in a Unique Ring Family decomposition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UniqueRingFamily {
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
    weight: usize,
    relevant_cycle_count: RelevantCycleCount,
}

impl UniqueRingFamily {
    /// Source graph nodes occurring in at least one cycle of the family.
    pub fn nodes(&self) -> &[NodeId] {
        &self.nodes
    }

    /// Source graph edges occurring in at least one cycle of the family.
    pub fn edges(&self) -> &[EdgeId] {
        &self.edges
    }

    /// Common cycle length in source-edge units.
    pub fn weight(&self) -> usize {
        self.weight
    }

    /// Exact number of relevant cycles represented by the family.
    pub fn relevant_cycle_count(&self) -> &RelevantCycleCount {
        &self.relevant_cycle_count
    }
}

/// Polynomial representation of the graph's Unique Ring Families.
///
/// Relevant cycles are retained as compact family state and are expanded only
/// by [`UniqueRingFamilies::visit_relevant_cycles`].
#[derive(Debug)]
pub struct UniqueRingFamilies {
    families: Vec<UniqueRingFamily>,
    node_to_families: Vec<Vec<UniqueRingFamilyId>>,
    edge_to_families: Vec<Vec<UniqueRingFamilyId>>,
    decomposition: urf::UrfDecomposition,
}

impl UniqueRingFamilies {
    /// Number of Unique Ring Families.
    pub fn count(&self) -> usize {
        self.families.len()
    }

    /// Iterate over all valid family identifiers.
    pub fn ids(&self) -> impl Iterator<Item = UniqueRingFamilyId> {
        (0..self.families.len()).map(|index| UniqueRingFamilyId(index as u32))
    }

    /// Iterate over the family descriptors.
    pub fn iter(&self) -> Iter<'_, UniqueRingFamily> {
        self.families.iter()
    }

    /// Return the family identified by `id`.
    pub fn get(&self, id: UniqueRingFamilyId) -> Option<&UniqueRingFamily> {
        self.families.get(id.index())
    }

    /// Return the families containing `node`.
    pub fn families_containing_node(&self, node: NodeId) -> &[UniqueRingFamilyId] {
        self.node_to_families
            .get(node.index())
            .map_or(&[], Vec::as_slice)
    }

    /// Return the families containing `edge`.
    pub fn families_containing_edge(&self, edge: EdgeId) -> &[UniqueRingFamilyId] {
        self.edge_to_families
            .get(edge.index())
            .map_or(&[], Vec::as_slice)
    }

    /// Visit the relevant cycles in one family until traversal completes or
    /// the visitor returns [`ControlFlow::Break`].
    ///
    /// # Panics
    ///
    /// Panics when `family` is not a valid identifier in this result.
    pub fn visit_relevant_cycles<B>(
        &self,
        family: UniqueRingFamilyId,
        visitor: impl FnMut(Cycle) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        self.decomposition.visit(family, visitor)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ShortestCycleAlgorithm {
    Bfs,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SimpleCycleEnumerationAlgorithm {
    ReadTarjan,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RelevantCycleEnumerationAlgorithm {
    Vismara,
}

/// Algorithm used to select a minimum cycle basis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MinimumCycleBasisAlgorithm {
    Horton,
}

/// Algorithm used to decompose a graph into Unique Ring Families.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UniqueRingFamilyAlgorithm {
    Kolodzik,
}

impl Graph {
    pub fn shortest_cycle_through_edge(
        &self,
        edge: EdgeId,
        alg: ShortestCycleAlgorithm,
    ) -> Option<usize> {
        match alg {
            ShortestCycleAlgorithm::Bfs => self.shortest_cycle_through_edge_bfs(edge),
        }
    }

    pub fn shortest_cycle_through_node(
        &self,
        node: NodeId,
        alg: ShortestCycleAlgorithm,
    ) -> Option<usize> {
        match alg {
            ShortestCycleAlgorithm::Bfs => self.shortest_cycle_through_node_bfs(node),
        }
    }

    /// Visits every simple cycle using the direct simple-graph path.
    ///
    /// Returns [`NonSimpleGraphError`] before invoking `visitor` when the graph
    /// contains a loop or parallel edges.
    pub fn try_visit_simple_cycles<B, F>(
        &self,
        max_cycle_size: usize,
        alg: SimpleCycleEnumerationAlgorithm,
        mut visitor: F,
    ) -> Result<ControlFlow<B>, NonSimpleGraphError>
    where
        F: FnMut(Cycle) -> ControlFlow<B>,
    {
        if !self.is_simple() {
            return Err(NonSimpleGraphError);
        }
        Ok(match alg {
            SimpleCycleEnumerationAlgorithm::ReadTarjan => {
                self.visit_simple_cycles_read_tarjan(max_cycle_size, &mut visitor)
            }
        })
    }

    /// Visits every simple cycle using the multigraph-capable fallback path.
    pub fn visit_simple_cycles_fallback<B, F>(
        &self,
        max_cycle_size: usize,
        alg: SimpleCycleEnumerationAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(Cycle) -> ControlFlow<B>,
    {
        match alg {
            SimpleCycleEnumerationAlgorithm::ReadTarjan => {
                self.visit_simple_cycles_read_tarjan(max_cycle_size, &mut visitor)
            }
        }
    }

    /// Visits every simple cycle until traversal completes or the visitor
    /// returns [`ControlFlow::Break`].
    ///
    /// Traversal uses the direct path for simple graphs and the
    /// multigraph-capable fallback otherwise. It is deterministic for a fixed
    /// graph representation, but its order is not a canonical ordering
    /// contract.
    pub fn visit_simple_cycles<B, F>(
        &self,
        max_cycle_size: usize,
        alg: SimpleCycleEnumerationAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(Cycle) -> ControlFlow<B>,
    {
        if self.is_simple() {
            match alg {
                SimpleCycleEnumerationAlgorithm::ReadTarjan => {
                    self.visit_simple_cycles_read_tarjan(max_cycle_size, &mut visitor)
                }
            }
        } else {
            self.visit_simple_cycles_fallback(max_cycle_size, alg, visitor)
        }
    }

    /// Collects every simple cycle using the direct simple-graph path.
    pub fn try_enumerate_simple_cycles(
        &self,
        max_cycle_size: usize,
        alg: SimpleCycleEnumerationAlgorithm,
    ) -> Result<Vec<Cycle>, NonSimpleGraphError> {
        let mut cycles = Vec::new();
        let _: ControlFlow<()> = self.try_visit_simple_cycles(max_cycle_size, alg, |cycle| {
            cycles.push(cycle);
            ControlFlow::Continue(())
        })?;
        Ok(cycles)
    }

    /// Collects every simple cycle using the multigraph-capable fallback path.
    pub fn enumerate_simple_cycles_fallback(
        &self,
        max_cycle_size: usize,
        alg: SimpleCycleEnumerationAlgorithm,
    ) -> Vec<Cycle> {
        let mut cycles = Vec::new();
        let _: ControlFlow<()> = self.visit_simple_cycles_fallback(max_cycle_size, alg, |cycle| {
            cycles.push(cycle);
            ControlFlow::Continue(())
        });
        cycles
    }

    /// Collects every simple cycle selected by `alg`.
    pub fn enumerate_simple_cycles(
        &self,
        max_cycle_size: usize,
        alg: SimpleCycleEnumerationAlgorithm,
    ) -> Vec<Cycle> {
        let mut cycles = Vec::new();
        let _: ControlFlow<()> = self.visit_simple_cycles(max_cycle_size, alg, |cycle| {
            cycles.push(cycle);
            ControlFlow::Continue(())
        });
        cycles
    }

    /// Visits relevant cycles using the direct simple-graph path.
    ///
    /// Returns [`NonSimpleGraphError`] before invoking `visitor` when the graph
    /// contains a loop or parallel edges.
    pub fn try_visit_relevant_cycles<B, F>(
        &self,
        max_cycle_size: usize,
        alg: RelevantCycleEnumerationAlgorithm,
        mut visitor: F,
    ) -> Result<ControlFlow<B>, NonSimpleGraphError>
    where
        F: FnMut(Cycle) -> ControlFlow<B>,
    {
        if !self.is_simple() {
            return Err(NonSimpleGraphError);
        }
        Ok(match alg {
            RelevantCycleEnumerationAlgorithm::Vismara => {
                self::relevant::visit_relevant_cycles_vismara(self, max_cycle_size, &mut visitor)
            }
        })
    }

    /// Visits relevant cycles using the multigraph-capable fallback path.
    ///
    /// Cycles retain source graph identifiers after loop extraction,
    /// subdivision, and projection.
    pub fn visit_relevant_cycles_fallback<B, F>(
        &self,
        max_cycle_size: usize,
        alg: RelevantCycleEnumerationAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(Cycle) -> ControlFlow<B>,
    {
        match alg {
            RelevantCycleEnumerationAlgorithm::Vismara => {
                self::relevant::visit_relevant_cycles_vismara_fallback(
                    self,
                    max_cycle_size,
                    &mut visitor,
                )
            }
        }
    }

    /// Visits relevant cycles until traversal completes or the visitor returns
    /// [`ControlFlow::Break`].
    ///
    /// Traversal uses the direct path for simple graphs and the
    /// multigraph-capable fallback otherwise. It is deterministic for a fixed
    /// graph representation, but its order is not a canonical ordering
    /// contract. Cycles always use source graph identifiers, including for
    /// loops and parallel edges. Only cycles with at most `max_cycle_size`
    /// source edges are visited.
    pub fn visit_relevant_cycles<B, F>(
        &self,
        max_cycle_size: usize,
        alg: RelevantCycleEnumerationAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(Cycle) -> ControlFlow<B>,
    {
        if self.is_simple() {
            match alg {
                RelevantCycleEnumerationAlgorithm::Vismara => {
                    self::relevant::visit_relevant_cycles_vismara(
                        self,
                        max_cycle_size,
                        &mut visitor,
                    )
                }
            }
        } else {
            self.visit_relevant_cycles_fallback(max_cycle_size, alg, visitor)
        }
    }

    /// Collects relevant cycles using the direct simple-graph path.
    pub fn try_enumerate_relevant_cycles(
        &self,
        max_cycle_size: usize,
        alg: RelevantCycleEnumerationAlgorithm,
    ) -> Result<Vec<Cycle>, NonSimpleGraphError> {
        let mut cycles = Vec::new();
        let _: ControlFlow<()> = self.try_visit_relevant_cycles(max_cycle_size, alg, |cycle| {
            cycles.push(cycle);
            ControlFlow::Continue(())
        })?;
        Ok(cycles)
    }

    /// Collects relevant cycles using the multigraph-capable fallback path.
    pub fn enumerate_relevant_cycles_fallback(
        &self,
        max_cycle_size: usize,
        alg: RelevantCycleEnumerationAlgorithm,
    ) -> Vec<Cycle> {
        let mut cycles = Vec::new();
        let _: ControlFlow<()> =
            self.visit_relevant_cycles_fallback(max_cycle_size, alg, |cycle| {
                cycles.push(cycle);
                ControlFlow::Continue(())
            });
        cycles
    }

    /// Collects the relevant cycles selected by `alg`.
    ///
    /// This is the collecting counterpart of [`Graph::visit_relevant_cycles`].
    pub fn enumerate_relevant_cycles(
        &self,
        max_cycle_size: usize,
        alg: RelevantCycleEnumerationAlgorithm,
    ) -> Vec<Cycle> {
        let mut cycles = Vec::new();
        let _: ControlFlow<()> = self.visit_relevant_cycles(max_cycle_size, alg, |cycle| {
            cycles.push(cycle);
            ControlFlow::Continue(())
        });
        cycles
    }

    /// Select one minimum-total-length basis of the binary cycle space.
    ///
    /// The returned basis is not canonically ordered. Self-loops contribute
    /// one-edge basis cycles, and parallel-edge cycles retain their source
    /// `EdgeId` identities.
    pub fn minimum_cycle_basis(&self, alg: MinimumCycleBasisAlgorithm) -> MinimumCycleBasis {
        match alg {
            MinimumCycleBasisAlgorithm::Horton => minimum_cycle_basis_horton(self),
        }
    }

    /// Decompose the graph into Unique Ring Families.
    pub fn unique_ring_families(&self, alg: UniqueRingFamilyAlgorithm) -> UniqueRingFamilies {
        match alg {
            UniqueRingFamilyAlgorithm::Kolodzik => urf::unique_ring_families_kolodzik(self),
        }
    }

    // BFS with edge exclusion. O(V+E).
    fn shortest_cycle_through_edge_bfs(&self, edge: EdgeId) -> Option<usize> {
        let [u, v] = self.edge_endpoints(edge);
        if u == v {
            return Some(1);
        }

        let mut dist = vec![u32::MAX; self.node_bound()];
        let mut queue = VecDeque::new();
        dist[u.index()] = 0;
        queue.push_back(u);
        while let Some(current) = queue.pop_front() {
            for nbr in self.neighbors(current) {
                if nbr.edge == edge {
                    continue;
                }
                if dist[nbr.node.index()] == u32::MAX {
                    dist[nbr.node.index()] = dist[current.index()] + 1;
                    if nbr.node == v {
                        return Some(dist[v.index()] as usize + 1);
                    }
                    queue.push_back(nbr.node);
                }
            }
        }
        None
    }

    // Min over incident edges. O(deg * (V+E)).
    fn shortest_cycle_through_node_bfs(&self, node: NodeId) -> Option<usize> {
        self.neighbors(node)
            .iter()
            .filter_map(|nbr| self.shortest_cycle_through_edge_bfs(nbr.edge))
            .min()
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::ops::ControlFlow;

    use num_bigint::BigUint;
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::MinimumCycleBasisAlgorithm::Horton;
    use super::RelevantCycleEnumerationAlgorithm::Vismara;
    use super::ShortestCycleAlgorithm::Bfs;
    use super::SimpleCycleEnumerationAlgorithm::ReadTarjan;
    use super::UniqueRingFamilyAlgorithm::Kolodzik;
    use super::{Cycle, NonSimpleGraphError, UniqueRingFamilyId};
    use crate::graph::{EdgeId, Graph, NodeId};

    type ExpectedUniqueRingFamily = (Vec<NodeId>, Vec<EdgeId>, usize, u64, Vec<Cycle>);

    #[rstest]
    #[case::self_loop(
        Graph::new(1, &[[0, 0]]),
        vec![NodeId(0)],
        vec![EdgeId(0)],
        Cycle { nodes: vec![NodeId(0)], edges: vec![EdgeId(0)] },
    )]
    #[case::digon(
        Graph::new(2, &[[0, 1], [0, 1]]),
        vec![NodeId(1), NodeId(0)],
        vec![EdgeId(1), EdgeId(0)],
        Cycle {
            nodes: vec![NodeId(0), NodeId(1)],
            edges: vec![EdgeId(0), EdgeId(1)],
        },
    )]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
        vec![NodeId(0), NodeId(1), NodeId(2)],
        vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        },
    )]
    #[case::rotated(
        Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
        vec![NodeId(1), NodeId(2), NodeId(0)],
        vec![EdgeId(1), EdgeId(2), EdgeId(0)],
        Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        },
    )]
    #[case::reversed(
        Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
        vec![NodeId(0), NodeId(2), NodeId(1)],
        vec![EdgeId(2), EdgeId(1), EdgeId(0)],
        Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        },
    )]
    #[case::parallel_first(
        Graph::new(2, &[[0, 1], [0, 1], [0, 1]]),
        vec![NodeId(0), NodeId(1)],
        vec![EdgeId(1), EdgeId(0)],
        Cycle {
            nodes: vec![NodeId(0), NodeId(1)],
            edges: vec![EdgeId(0), EdgeId(1)],
        },
    )]
    #[case::parallel_second(
        Graph::new(2, &[[0, 1], [0, 1], [0, 1]]),
        vec![NodeId(0), NodeId(1)],
        vec![EdgeId(2), EdgeId(0)],
        Cycle {
            nodes: vec![NodeId(0), NodeId(1)],
            edges: vec![EdgeId(0), EdgeId(2)],
        },
    )]
    fn test_cycle_normalized(
        #[case] graph: Graph,
        #[case] nodes: Vec<NodeId>,
        #[case] edges: Vec<EdgeId>,
        #[case] expected: Cycle,
    ) {
        assert_eq!(Cycle::normalized(&graph, nodes, edges), expected);
    }

    #[rstest]
    #[case::triangle(
        Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(3), EdgeId(4), EdgeId(5)],
        },
        vec![NodeId(0), NodeId(1), NodeId(2)],
        vec![EdgeId(3), EdgeId(4), EdgeId(5)],
        3,
    )]
    fn test_cycle_accessors(
        #[case] cycle: Cycle,
        #[case] expected_nodes: Vec<NodeId>,
        #[case] expected_edges: Vec<EdgeId>,
        #[case] expected_len: usize,
    ) {
        assert_eq!(cycle.nodes(), expected_nodes.as_slice());
        assert_eq!(cycle.edges(), expected_edges.as_slice());
        assert_eq!(cycle.length(), expected_len);
    }

    #[rstest]
    #[case::rotation_reversal(
        Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
        vec![
            (
                vec![NodeId(1), NodeId(2), NodeId(0)],
                vec![EdgeId(1), EdgeId(2), EdgeId(0)],
            ),
            (
                vec![NodeId(0), NodeId(2), NodeId(1)],
                vec![EdgeId(2), EdgeId(1), EdgeId(0)],
            ),
        ],
        HashSet::from([Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }]),
    )]
    #[case::parallel_edges(
        Graph::new(2, &[[0, 1], [0, 1], [0, 1]]),
        vec![
            (
                vec![NodeId(0), NodeId(1)],
                vec![EdgeId(1), EdgeId(0)],
            ),
            (
                vec![NodeId(0), NodeId(1)],
                vec![EdgeId(2), EdgeId(0)],
            ),
        ],
        HashSet::from([
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0), EdgeId(1)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0), EdgeId(2)],
            },
        ]),
    )]
    fn test_cycle_hash(
        #[case] graph: Graph,
        #[case] paths: Vec<(Vec<NodeId>, Vec<EdgeId>)>,
        #[case] expected: HashSet<Cycle>,
    ) {
        let actual: HashSet<Cycle> = paths
            .into_iter()
            .map(|(nodes, edges)| Cycle::normalized(&graph, nodes, edges))
            .collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::self_loop(1, vec![[0, 0]], EdgeId(0), Some(1))]
    #[case::digon(2, vec![[0, 1], [0, 1]], EdgeId(0), Some(2))]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], EdgeId(0), Some(3))]
    #[case::square(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]], EdgeId(0), Some(4))]
    #[case::bridge(4, vec![[0, 1], [1, 2], [2, 3]], EdgeId(1), None)]
    #[case::naphthalene_shared(10, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [3, 6],
        [6, 7], [7, 8], [8, 9], [9, 4]], EdgeId(3), Some(6))]
    #[case::naphthalene_outer(10, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [3, 6],
        [6, 7], [7, 8], [8, 9], [9, 4]], EdgeId(0), Some(6))]
    fn test_graph_shortest_cycle_through_edge(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] edge: EdgeId,
        #[case] expected: Option<usize>,
    ) {
        let graph = Graph::new(node_count, &edges);
        assert_eq!(graph.shortest_cycle_through_edge(edge, Bfs), expected);
    }

    #[rstest]
    #[case::self_loop(1, vec![[0, 0]], NodeId(0), Some(1))]
    #[case::digon(2, vec![[0, 1], [0, 1]], NodeId(0), Some(2))]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], NodeId(0), Some(3))]
    #[case::pendant(3, vec![[0, 1], [1, 2]], NodeId(0), None)]
    #[case::naphthalene_bridgehead(10, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [3, 6],
        [6, 7], [7, 8], [8, 9], [9, 4]], NodeId(3), Some(6))]
    fn test_graph_shortest_cycle_through_node(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] node: NodeId,
        #[case] expected: Option<usize>,
    ) {
        let graph = Graph::new(node_count, &edges);
        assert_eq!(graph.shortest_cycle_through_node(node, Bfs), expected);
    }

    #[rstest]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    fn test_graph_try_visit_simple_cycles(#[case] graph: Graph, #[case] expected: Vec<Cycle>) {
        let mut visited = Vec::new();
        let result = graph.try_visit_simple_cycles(usize::MAX, ReadTarjan, |cycle| {
            visited.push(cycle);
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, Ok(ControlFlow::Continue(())));
        assert_eq!(visited, expected);
    }

    #[rstest]
    #[case::loop_edge(Graph::new(1, &[[0, 0]]))]
    #[case::parallel_edges(Graph::new(2, &[[0, 1], [0, 1]]))]
    fn test_graph_try_visit_simple_cycles_error(#[case] graph: Graph) {
        let mut visited = Vec::new();
        let result = graph.try_visit_simple_cycles(usize::MAX, ReadTarjan, |cycle| {
            visited.push(cycle);
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, Err(NonSimpleGraphError));
        assert_eq!(visited, Vec::<Cycle>::new());
    }

    #[rstest]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    fn test_graph_try_visit_simple_cycles_break(
        #[case] graph: Graph,
        #[case] expected: Vec<Cycle>,
    ) {
        let mut visited = Vec::new();
        let result = graph.try_visit_simple_cycles(usize::MAX, ReadTarjan, |cycle| {
            visited.push(cycle);
            ControlFlow::Break(visited.len())
        });

        assert_eq!(result, Ok(ControlFlow::Break(1)));
        assert_eq!(visited, expected);
    }

    #[rstest]
    #[case::bounded_multigraph(
        Graph::new(3, &[[0, 0], [0, 1], [0, 1], [1, 2], [0, 2]]),
        2,
        vec![
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(0)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1), EdgeId(2)],
            },
        ],
    )]
    fn test_graph_visit_simple_cycles_fallback(
        #[case] graph: Graph,
        #[case] max_cycle_size: usize,
        #[case] expected: Vec<Cycle>,
    ) {
        let mut visited = Vec::new();
        let result = graph.visit_simple_cycles_fallback(max_cycle_size, ReadTarjan, |cycle| {
            visited.push(cycle);
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(visited, expected);
    }

    #[rstest]
    #[case::first(
        Graph::new(3, &[[0, 0], [0, 1], [0, 1], [1, 2], [0, 2]]),
        vec![Cycle {
            nodes: vec![NodeId(0)],
            edges: vec![EdgeId(0)],
        }],
    )]
    fn test_graph_visit_simple_cycles_fallback_break(
        #[case] graph: Graph,
        #[case] expected: Vec<Cycle>,
    ) {
        let mut visited = Vec::new();
        let result = graph.visit_simple_cycles_fallback(usize::MAX, ReadTarjan, |cycle| {
            visited.push(cycle);
            ControlFlow::Break(visited.len())
        });

        assert_eq!(result, ControlFlow::Break(1));
        assert_eq!(visited, expected);
    }

    #[rstest]
    fn test_graph_visit_simple_cycles() {
        let graph = Graph::new(6, &[[0, 1], [1, 2], [2, 0], [3, 4], [3, 4], [5, 5]]);
        let mut visited = Vec::new();
        let result = graph.visit_simple_cycles(usize::MAX, ReadTarjan, |cycle| {
            visited.push(cycle);
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(
            visited,
            vec![
                Cycle {
                    nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                    edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                },
                Cycle {
                    nodes: vec![NodeId(3), NodeId(4)],
                    edges: vec![EdgeId(3), EdgeId(4)],
                },
                Cycle {
                    nodes: vec![NodeId(5)],
                    edges: vec![EdgeId(5)],
                },
            ]
        );
    }

    #[rstest]
    #[case::first(
        1,
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    #[case::prefix(
        2,
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            },
            Cycle {
                nodes: vec![NodeId(3), NodeId(4)],
                edges: vec![EdgeId(3), EdgeId(4)],
            },
        ],
    )]
    fn test_graph_visit_simple_cycles_break(
        #[case] stop_after: usize,
        #[case] expected: Vec<Cycle>,
    ) {
        let graph = Graph::new(6, &[[0, 1], [1, 2], [2, 0], [3, 4], [3, 4], [5, 5]]);
        let mut visited = Vec::new();
        let result = graph.visit_simple_cycles(usize::MAX, ReadTarjan, |cycle| {
            visited.push(cycle);
            if visited.len() == stop_after {
                ControlFlow::Break(visited.len())
            } else {
                ControlFlow::Continue(())
            }
        });

        assert_eq!(result, ControlFlow::Break(stop_after));
        assert_eq!(visited, expected);
    }

    #[rstest]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    fn test_graph_try_enumerate_simple_cycles(#[case] graph: Graph, #[case] expected: Vec<Cycle>) {
        assert_eq!(
            graph.try_enumerate_simple_cycles(usize::MAX, ReadTarjan),
            Ok(expected),
        );
    }

    #[rstest]
    #[case::loop_edge(Graph::new(1, &[[0, 0]]))]
    #[case::parallel_edges(Graph::new(2, &[[0, 1], [0, 1]]))]
    fn test_graph_try_enumerate_simple_cycles_error(#[case] graph: Graph) {
        assert_eq!(
            graph.try_enumerate_simple_cycles(usize::MAX, ReadTarjan),
            Err(NonSimpleGraphError),
        );
    }

    #[rstest]
    #[case::simple(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        usize::MAX,
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    #[case::multigraph(
        Graph::new(2, &[[0, 0], [0, 1], [0, 1]]),
        usize::MAX,
        vec![
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(0)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1), EdgeId(2)],
            },
        ],
    )]
    fn test_graph_enumerate_simple_cycles_fallback(
        #[case] graph: Graph,
        #[case] max_cycle_size: usize,
        #[case] expected: Vec<Cycle>,
    ) {
        assert_eq!(
            graph.enumerate_simple_cycles_fallback(max_cycle_size, ReadTarjan),
            expected,
        );
    }

    #[rstest]
    #[case::zero_bound(
        Graph::new(1, &[[0, 0]]),
        0,
        vec![],
    )]
    #[case::one_bound(
        Graph::new(3, &[[0, 0], [0, 1], [0, 1], [1, 2], [2, 0]]),
        1,
        vec![Cycle {
            nodes: vec![NodeId(0)],
            edges: vec![EdgeId(0)],
        }],
    )]
    #[case::two_bound(
        Graph::new(3, &[[0, 0], [0, 1], [0, 1], [1, 2], [2, 0]]),
        2,
        vec![
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(0)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1), EdgeId(2)],
            },
        ],
    )]
    #[case::disconnected(
        Graph::new(
            8,
            &[
                [0, 1], [1, 2], [2, 0], [2, 3],
                [4, 5], [5, 6], [6, 7], [7, 4],
            ],
        ),
        usize::MAX,
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            },
            Cycle {
                nodes: vec![NodeId(4), NodeId(5), NodeId(6), NodeId(7)],
                edges: vec![EdgeId(4), EdgeId(5), EdgeId(6), EdgeId(7)],
            },
        ],
    )]
    #[case::fused_bridge(
        Graph::new(
            6,
            &[
                [0, 1], [1, 2], [2, 0],
                [2, 3], [3, 4], [4, 2], [4, 5],
            ],
        ),
        usize::MAX,
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            },
            Cycle {
                nodes: vec![NodeId(2), NodeId(3), NodeId(4)],
                edges: vec![EdgeId(3), EdgeId(4), EdgeId(5)],
            },
        ],
    )]
    #[case::parallel_alternatives(
        Graph::new(2, &[[0, 1], [0, 1], [0, 1]]),
        usize::MAX,
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0), EdgeId(1)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0), EdgeId(2)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1), EdgeId(2)],
            },
        ],
    )]
    fn test_graph_enumerate_simple_cycles(
        #[case] graph: Graph,
        #[case] max_cycle_size: usize,
        #[case] expected: Vec<Cycle>,
    ) {
        assert_eq!(
            graph.enumerate_simple_cycles(max_cycle_size, ReadTarjan),
            expected
        );
    }

    #[rstest]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    fn test_graph_try_visit_relevant_cycles(#[case] graph: Graph, #[case] expected: Vec<Cycle>) {
        let mut visited = Vec::new();
        let result = graph.try_visit_relevant_cycles(usize::MAX, Vismara, |cycle| {
            visited.push(cycle);
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, Ok(ControlFlow::Continue(())));
        assert_eq!(visited, expected);
    }

    #[rstest]
    #[case::loop_edge(Graph::new(1, &[[0, 0]]))]
    #[case::parallel_edges(Graph::new(2, &[[0, 1], [0, 1]]))]
    fn test_graph_try_visit_relevant_cycles_error(#[case] graph: Graph) {
        let mut visited = Vec::new();
        let result = graph.try_visit_relevant_cycles(usize::MAX, Vismara, |cycle| {
            visited.push(cycle);
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, Err(NonSimpleGraphError));
        assert_eq!(visited, Vec::<Cycle>::new());
    }

    #[rstest]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    fn test_graph_try_visit_relevant_cycles_break(
        #[case] graph: Graph,
        #[case] expected: Vec<Cycle>,
    ) {
        let mut visited = Vec::new();
        let result = graph.try_visit_relevant_cycles(usize::MAX, Vismara, |cycle| {
            visited.push(cycle);
            ControlFlow::Break(visited.len())
        });

        assert_eq!(result, Ok(ControlFlow::Break(1)));
        assert_eq!(visited, expected);
    }

    #[rstest]
    #[case::bounded_multigraph(
        Graph::new(3, &[[0, 0], [0, 1], [1, 2], [0, 2], [0, 1]]),
        2,
        vec![
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(0)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1), EdgeId(4)],
            },
        ],
    )]
    fn test_graph_visit_relevant_cycles_fallback(
        #[case] graph: Graph,
        #[case] max_cycle_size: usize,
        #[case] expected: Vec<Cycle>,
    ) {
        let mut visited = Vec::new();
        let result = graph.visit_relevant_cycles_fallback(max_cycle_size, Vismara, |cycle| {
            visited.push(cycle);
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(visited, expected);
    }

    #[rstest]
    #[case::first(
        Graph::new(3, &[[0, 0], [0, 1], [1, 2], [0, 2], [0, 1]]),
        vec![Cycle {
            nodes: vec![NodeId(0)],
            edges: vec![EdgeId(0)],
        }],
    )]
    fn test_graph_visit_relevant_cycles_fallback_break(
        #[case] graph: Graph,
        #[case] expected: Vec<Cycle>,
    ) {
        let mut visited = Vec::new();
        let result = graph.visit_relevant_cycles_fallback(usize::MAX, Vismara, |cycle| {
            visited.push(cycle);
            ControlFlow::Break(visited.len())
        });

        assert_eq!(result, ControlFlow::Break(1));
        assert_eq!(visited, expected);
    }

    #[rstest]
    fn test_graph_visit_relevant_cycles() {
        let graph = Graph::new(3, &[[0, 0], [0, 1], [1, 2], [0, 2], [0, 1]]);
        let mut visited = Vec::new();
        let result = graph.visit_relevant_cycles(usize::MAX, Vismara, |cycle| {
            visited.push(cycle);
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(
            visited,
            vec![
                Cycle {
                    nodes: vec![NodeId(0)],
                    edges: vec![EdgeId(0)],
                },
                Cycle {
                    nodes: vec![NodeId(0), NodeId(1)],
                    edges: vec![EdgeId(1), EdgeId(4)],
                },
                Cycle {
                    nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                    edges: vec![EdgeId(1), EdgeId(2), EdgeId(3)],
                },
                Cycle {
                    nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                    edges: vec![EdgeId(4), EdgeId(2), EdgeId(3)],
                },
            ]
        );
    }

    #[rstest]
    #[case::first(
        1,
        vec![Cycle {
            nodes: vec![NodeId(0)],
            edges: vec![EdgeId(0)],
        }],
    )]
    #[case::second(
        2,
        vec![
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(0)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1), EdgeId(4)],
            },
        ],
    )]
    fn test_graph_visit_relevant_cycles_break(
        #[case] stop_after: usize,
        #[case] expected: Vec<Cycle>,
    ) {
        let graph = Graph::new(3, &[[0, 0], [0, 1], [1, 2], [0, 2], [0, 1]]);
        let mut visited = Vec::new();
        let result = graph.visit_relevant_cycles(usize::MAX, Vismara, |cycle| {
            visited.push(cycle);
            if visited.len() == stop_after {
                ControlFlow::Break(visited.len())
            } else {
                ControlFlow::Continue(())
            }
        });

        assert_eq!(result, ControlFlow::Break(stop_after));
        assert_eq!(visited, expected);
    }

    #[rstest]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    fn test_graph_try_enumerate_relevant_cycles(
        #[case] graph: Graph,
        #[case] expected: Vec<Cycle>,
    ) {
        assert_eq!(
            graph.try_enumerate_relevant_cycles(usize::MAX, Vismara),
            Ok(expected),
        );
    }

    #[rstest]
    #[case::loop_edge(Graph::new(1, &[[0, 0]]))]
    #[case::parallel_edges(Graph::new(2, &[[0, 1], [0, 1]]))]
    fn test_graph_try_enumerate_relevant_cycles_error(#[case] graph: Graph) {
        assert_eq!(
            graph.try_enumerate_relevant_cycles(usize::MAX, Vismara),
            Err(NonSimpleGraphError),
        );
    }

    #[rstest]
    #[case::simple(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        usize::MAX,
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    #[case::multigraph(
        Graph::new(2, &[[0, 0], [0, 1], [0, 1]]),
        usize::MAX,
        vec![
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(0)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1), EdgeId(2)],
            },
        ],
    )]
    fn test_graph_enumerate_relevant_cycles_fallback(
        #[case] graph: Graph,
        #[case] max_cycle_size: usize,
        #[case] expected: Vec<Cycle>,
    ) {
        assert_eq!(
            graph.enumerate_relevant_cycles_fallback(max_cycle_size, Vismara),
            expected,
        );
    }

    #[rstest]
    #[case::direct(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        usize::MAX,
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    #[case::fallback(
        Graph::new(5, &[[0, 1], [1, 2], [0, 2], [3, 4], [3, 4]]),
        usize::MAX,
        vec![
            Cycle {
                nodes: vec![NodeId(3), NodeId(4)],
                edges: vec![EdgeId(3), EdgeId(4)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            },
        ],
    )]
    #[case::loops(
        Graph::new(2, &[[0, 0], [0, 0], [1, 1]]),
        usize::MAX,
        vec![
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(0)],
            },
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(1)],
            },
            Cycle {
                nodes: vec![NodeId(1)],
                edges: vec![EdgeId(2)],
            },
        ],
    )]
    #[case::parallel(
        Graph::new(2, &[[0, 1], [0, 1], [0, 1]]),
        usize::MAX,
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0), EdgeId(1)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0), EdgeId(2)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1), EdgeId(2)],
            },
        ],
    )]
    #[case::longer_parallel(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2], [0, 1]]),
        usize::MAX,
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0), EdgeId(3)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                edges: vec![EdgeId(3), EdgeId(1), EdgeId(2)],
            },
        ],
    )]
    #[case::mixed(
        Graph::new(
            7,
            &[
                [0, 0],
                [1, 2], [2, 3], [3, 4], [1, 4],
                [5, 6], [5, 6], [5, 6],
            ],
        ),
        usize::MAX,
        vec![
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(0)],
            },
            Cycle {
                nodes: vec![NodeId(5), NodeId(6)],
                edges: vec![EdgeId(5), EdgeId(6)],
            },
            Cycle {
                nodes: vec![NodeId(5), NodeId(6)],
                edges: vec![EdgeId(5), EdgeId(7)],
            },
            Cycle {
                nodes: vec![NodeId(5), NodeId(6)],
                edges: vec![EdgeId(6), EdgeId(7)],
            },
            Cycle {
                nodes: vec![NodeId(1), NodeId(2), NodeId(3), NodeId(4)],
                edges: vec![EdgeId(1), EdgeId(2), EdgeId(3), EdgeId(4)],
            },
        ],
    )]
    #[case::bounded(
        Graph::new(3, &[[0, 0], [0, 1], [1, 2], [0, 2], [0, 1]]),
        2,
        vec![
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(0)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1), EdgeId(4)],
            },
        ],
    )]
    fn test_graph_enumerate_relevant_cycles(
        #[case] graph: Graph,
        #[case] max_cycle_size: usize,
        #[case] expected: Vec<Cycle>,
    ) {
        assert_eq!(
            graph.enumerate_relevant_cycles(max_cycle_size, Vismara),
            expected
        );
    }

    #[rstest]
    #[case::forest(
        Graph::new(4, &[[0, 1], [1, 2], [1, 3]]),
        vec![],
        0,
    )]
    #[case::disconnected(
        Graph::new(6, &[[0, 1], [1, 2], [0, 2], [3, 4], [4, 5], [3, 5]]),
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            },
            Cycle {
                nodes: vec![NodeId(3), NodeId(4), NodeId(5)],
                edges: vec![EdgeId(3), EdgeId(4), EdgeId(5)],
            },
        ],
        6,
    )]
    #[case::tied(
        Graph::new(
            4,
            &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]],
        ),
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                edges: vec![EdgeId(0), EdgeId(3), EdgeId(1)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(3)],
                edges: vec![EdgeId(0), EdgeId(4), EdgeId(2)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(2), NodeId(3)],
                edges: vec![EdgeId(1), EdgeId(5), EdgeId(2)],
            },
        ],
        9,
    )]
    #[case::loop_and_triangle(
        Graph::new(3, &[[0, 0], [0, 1], [1, 2], [0, 2]]),
        vec![
            Cycle {
                nodes: vec![NodeId(0)],
                edges: vec![EdgeId(0)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                edges: vec![EdgeId(1), EdgeId(2), EdgeId(3)],
            },
        ],
        4,
    )]
    #[case::parallel_pair(
        Graph::new(2, &[[0, 1], [0, 1]]),
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1)],
            edges: vec![EdgeId(0), EdgeId(1)],
        }],
        2,
    )]
    #[case::parallel_triple(
        Graph::new(2, &[[0, 1], [0, 1], [0, 1]]),
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0), EdgeId(1)],
            },
            Cycle {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0), EdgeId(2)],
            },
        ],
        4,
    )]
    fn test_graph_minimum_cycle_basis(
        #[case] graph: Graph,
        #[case] expected: Vec<Cycle>,
        #[case] total_length: usize,
    ) {
        let basis = graph.minimum_cycle_basis(Horton);
        assert_eq!(basis.iter().cloned().collect::<Vec<_>>(), expected);
        assert_eq!(basis.dimension(), expected.len());
        assert_eq!(basis.total_length(), total_length);
    }

    #[rstest]
    #[case::independent(
        Graph::new(
            6,
            &[[0, 1], [1, 2], [0, 2], [3, 4], [4, 5], [3, 5]],
        ),
        vec![
            (
                vec![NodeId(0), NodeId(1), NodeId(2)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                3,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                    edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                }],
            ),
            (
                vec![NodeId(3), NodeId(4), NodeId(5)],
                vec![EdgeId(3), EdgeId(4), EdgeId(5)],
                3,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(3), NodeId(4), NodeId(5)],
                    edges: vec![EdgeId(3), EdgeId(4), EdgeId(5)],
                }],
            ),
        ],
    )]
    #[case::fused(
        Graph::new(
            4,
            &[[0, 1], [1, 2], [0, 2], [1, 3], [2, 3]],
        ),
        vec![
            (
                vec![NodeId(0), NodeId(1), NodeId(2)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                3,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                    edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                }],
            ),
            (
                vec![NodeId(1), NodeId(2), NodeId(3)],
                vec![EdgeId(1), EdgeId(3), EdgeId(4)],
                3,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(1), NodeId(2), NodeId(3)],
                    edges: vec![EdgeId(1), EdgeId(4), EdgeId(3)],
                }],
            ),
        ],
    )]
    #[case::bridged(
        Graph::new(
            6,
            &[
                [0, 1], [1, 2], [0, 2],
                [3, 4], [4, 5], [3, 5],
                [2, 3],
            ],
        ),
        vec![
            (
                vec![NodeId(0), NodeId(1), NodeId(2)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                3,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                    edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                }],
            ),
            (
                vec![NodeId(3), NodeId(4), NodeId(5)],
                vec![EdgeId(3), EdgeId(4), EdgeId(5)],
                3,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(3), NodeId(4), NodeId(5)],
                    edges: vec![EdgeId(3), EdgeId(4), EdgeId(5)],
                }],
            ),
        ],
    )]
    #[case::symmetric(
        Graph::new(
            6,
            &[[0, 1], [1, 4], [0, 2], [2, 4], [0, 3], [3, 5], [4, 5]],
        ),
        vec![
            (
                vec![NodeId(0), NodeId(1), NodeId(2), NodeId(4)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
                4,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(0), NodeId(1), NodeId(4), NodeId(2)],
                    edges: vec![EdgeId(0), EdgeId(1), EdgeId(3), EdgeId(2)],
                }],
            ),
            (
                vec![
                    NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4), NodeId(5),
                ],
                vec![
                    EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3), EdgeId(4), EdgeId(5),
                    EdgeId(6),
                ],
                5,
                2,
                vec![
                    Cycle {
                        nodes: vec![
                            NodeId(0), NodeId(1), NodeId(4), NodeId(5), NodeId(3),
                        ],
                        edges: vec![
                            EdgeId(0), EdgeId(1), EdgeId(6), EdgeId(5), EdgeId(4),
                        ],
                    },
                    Cycle {
                        nodes: vec![
                            NodeId(0), NodeId(2), NodeId(4), NodeId(5), NodeId(3),
                        ],
                        edges: vec![
                            EdgeId(2), EdgeId(3), EdgeId(6), EdgeId(5), EdgeId(4),
                        ],
                    },
                ],
            ),
        ],
    )]
    #[case::quotient_relation(
        Graph::new(
            6,
            &[
                [0, 3], [0, 4], [1, 4], [3, 4],
                [0, 5], [1, 5], [2, 5], [3, 5],
            ],
        ),
        vec![
            (
                vec![NodeId(0), NodeId(3), NodeId(4)],
                vec![EdgeId(0), EdgeId(1), EdgeId(3)],
                3,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(0), NodeId(3), NodeId(4)],
                    edges: vec![EdgeId(0), EdgeId(3), EdgeId(1)],
                }],
            ),
            (
                vec![NodeId(0), NodeId(3), NodeId(5)],
                vec![EdgeId(0), EdgeId(4), EdgeId(7)],
                3,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(0), NodeId(3), NodeId(5)],
                    edges: vec![EdgeId(0), EdgeId(7), EdgeId(4)],
                }],
            ),
            (
                vec![
                    NodeId(0), NodeId(1), NodeId(3), NodeId(4), NodeId(5),
                ],
                vec![
                    EdgeId(1), EdgeId(2), EdgeId(3), EdgeId(4), EdgeId(5), EdgeId(7),
                ],
                4,
                2,
                vec![
                    Cycle {
                        nodes: vec![NodeId(0), NodeId(4), NodeId(1), NodeId(5)],
                        edges: vec![EdgeId(1), EdgeId(2), EdgeId(5), EdgeId(4)],
                    },
                    Cycle {
                        nodes: vec![NodeId(1), NodeId(4), NodeId(3), NodeId(5)],
                        edges: vec![EdgeId(2), EdgeId(3), EdgeId(7), EdgeId(5)],
                    },
                ],
            ),
        ],
    )]
    #[case::loops(
        Graph::new(2, &[[0, 0], [1, 1]]),
        vec![
            (
                vec![NodeId(0)],
                vec![EdgeId(0)],
                1,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(0)],
                    edges: vec![EdgeId(0)],
                }],
            ),
            (
                vec![NodeId(1)],
                vec![EdgeId(1)],
                1,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(1)],
                    edges: vec![EdgeId(1)],
                }],
            ),
        ],
    )]
    #[case::parallel(
        Graph::new(2, &[[0, 1], [0, 1], [0, 1]]),
        vec![
            (
                vec![NodeId(0), NodeId(1)],
                vec![EdgeId(0), EdgeId(1)],
                2,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(0), NodeId(1)],
                    edges: vec![EdgeId(0), EdgeId(1)],
                }],
            ),
            (
                vec![NodeId(0), NodeId(1)],
                vec![EdgeId(0), EdgeId(2)],
                2,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(0), NodeId(1)],
                    edges: vec![EdgeId(0), EdgeId(2)],
                }],
            ),
            (
                vec![NodeId(0), NodeId(1)],
                vec![EdgeId(1), EdgeId(2)],
                2,
                1,
                vec![Cycle {
                    nodes: vec![NodeId(0), NodeId(1)],
                    edges: vec![EdgeId(1), EdgeId(2)],
                }],
            ),
        ],
    )]
    fn test_graph_unique_ring_families(
        #[case] graph: Graph,
        #[case] expected: Vec<ExpectedUniqueRingFamily>,
    ) {
        let result = graph.unique_ring_families(Kolodzik);
        assert_eq!(result.count(), expected.len());

        for (id, (nodes, edges, weight, count, expected_cycles)) in result.ids().zip(expected) {
            let family = result.get(id).expect("a returned family id must be valid");
            let mut cycles = Vec::new();
            let flow = result.visit_relevant_cycles(id, |cycle| {
                cycles.push(cycle);
                ControlFlow::<()>::Continue(())
            });

            assert_eq!(family.nodes(), nodes);
            assert_eq!(family.edges(), edges);
            assert_eq!(family.weight(), weight);
            assert_eq!(&family.relevant_cycle_count().0, &BigUint::from(count));
            assert_eq!(flow, ControlFlow::Continue(()));
            assert_eq!(cycles, expected_cycles);
        }
    }

    #[rstest]
    #[case::first(
        1,
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(4), NodeId(5), NodeId(3)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(6), EdgeId(5), EdgeId(4)],
        }],
    )]
    fn test_unique_ring_families_visit_relevant_cycles(
        #[case] stop_after: usize,
        #[case] expected: Vec<Cycle>,
    ) {
        let graph = Graph::new(6, &[[0, 1], [1, 4], [0, 2], [2, 4], [0, 3], [3, 5], [4, 5]]);
        let families = graph.unique_ring_families(Kolodzik);
        let mut cycles = Vec::new();
        let flow = families.visit_relevant_cycles(UniqueRingFamilyId(1), |cycle| {
            cycles.push(cycle);
            if cycles.len() == stop_after {
                ControlFlow::Break(cycles.len())
            } else {
                ControlFlow::Continue(())
            }
        });

        assert_eq!(flow, ControlFlow::Break(stop_after));
        assert_eq!(cycles, expected);
    }
}
