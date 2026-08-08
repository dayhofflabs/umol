//! Compact shortest-path and cycle-family state for relevant cycles.
//!
//! The current implementation builds shortest-path DAGs and Vismara cycle
//! prototypes, then expands the selected families during visitation. See
//! [Vismara, *Union of all the Minimum Cycle Bases of a Graph*
//! (1997)](https://doi.org/10.37236/1294).

use std::cmp::Ordering;
use std::collections::{HashSet, VecDeque};
use std::ops::ControlFlow;

use num_bigint::BigUint;

use super::basis::{CycleVectorBasis, EdgeVector};
use super::Cycle;
use crate::algorithms::connectivity::BiconnectedComponentsAlgorithm;
use crate::graph::{EdgeId, Graph, Neighbor, NodeId};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct ShortestPath {
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
}

impl ShortestPath {
    pub(super) fn common_prefix_len(&self, other: &Self) -> usize {
        self.nodes
            .iter()
            .zip(&other.nodes)
            .take_while(|(left, right)| left == right)
            .count()
    }

    pub(super) fn cycle_with(
        &self,
        graph: &Graph,
        other: &Self,
        closing_edge: EdgeId,
    ) -> Option<Cycle> {
        let shared_len = self.common_prefix_len(other);
        if shared_len == 0 {
            return None;
        }

        let first_tail: HashSet<_> = self.nodes[shared_len..].iter().copied().collect();
        if other.nodes[shared_len..]
            .iter()
            .any(|node| first_tail.contains(node))
        {
            return None;
        }

        let first = *self.nodes.last()?;
        let second = *other.nodes.last()?;
        let endpoints = graph.edge_endpoints(closing_edge);
        if !((endpoints == [first, second]) || (endpoints == [second, first])) {
            return None;
        }

        let mut nodes = self.nodes[shared_len - 1..].to_vec();
        nodes.extend(other.nodes[shared_len..].iter().rev());
        if nodes.len() < 3 {
            return None;
        }

        let mut edges = self.edges[shared_len - 1..].to_vec();
        edges.push(closing_edge);
        edges.extend(other.edges[shared_len - 1..].iter().rev());
        if edges.iter().copied().collect::<HashSet<_>>().len() != edges.len() {
            return None;
        }

        Some(Cycle::normalized(graph, nodes, edges))
    }

    fn cycle_with_middle(
        &self,
        graph: &Graph,
        other: &Self,
        middle: NodeId,
        first_edge: EdgeId,
        second_edge: EdgeId,
    ) -> Option<Cycle> {
        let shared_len = self.common_prefix_len(other);
        if shared_len == 0 {
            return None;
        }

        let first_tail: HashSet<_> = self.nodes[shared_len..].iter().copied().collect();
        if other.nodes[shared_len..]
            .iter()
            .any(|node| first_tail.contains(node))
            || first_tail.contains(&middle)
            || other.nodes[shared_len..].contains(&middle)
        {
            return None;
        }

        let first = *self.nodes.last()?;
        let second = *other.nodes.last()?;
        if graph.edge_endpoints(first_edge) != [first.min(middle), first.max(middle)]
            || graph.edge_endpoints(second_edge) != [second.min(middle), second.max(middle)]
        {
            return None;
        }

        let mut nodes = self.nodes[shared_len - 1..].to_vec();
        nodes.push(middle);
        nodes.extend(other.nodes[shared_len..].iter().rev());

        let mut edges = self.edges[shared_len - 1..].to_vec();
        edges.push(first_edge);
        edges.push(second_edge);
        edges.extend(other.edges[shared_len - 1..].iter().rev());
        if edges.iter().copied().collect::<HashSet<_>>().len() != edges.len() {
            return None;
        }

        Some(Cycle::normalized(graph, nodes, edges))
    }
}

#[derive(Clone, Debug)]
pub(super) struct ShortestPathDag {
    root: NodeId,
    distances: Vec<Option<usize>>,
    predecessors: Vec<Vec<Neighbor>>,
    included: Vec<bool>,
}

impl ShortestPathDag {
    pub(super) fn new(graph: &Graph, root: NodeId) -> Self {
        let mut distances = vec![None; graph.node_count()];
        let mut predecessors = vec![Vec::new(); graph.node_count()];
        let mut queue = VecDeque::from([root]);
        distances[root.index()] = Some(0);

        while let Some(current) = queue.pop_front() {
            let next_distance = distances[current.index()].expect("queued node has a distance") + 1;
            for neighbor in graph.neighbors(current) {
                match distances[neighbor.node.index()] {
                    None => {
                        distances[neighbor.node.index()] = Some(next_distance);
                        predecessors[neighbor.node.index()].push(Neighbor {
                            node: current,
                            edge: neighbor.edge,
                        });
                        queue.push_back(neighbor.node);
                    }
                    Some(distance) if distance == next_distance => {
                        predecessors[neighbor.node.index()].push(Neighbor {
                            node: current,
                            edge: neighbor.edge,
                        });
                    }
                    Some(_) => {}
                }
            }
        }

        for alternatives in &mut predecessors {
            alternatives.sort_unstable_by_key(|neighbor| (neighbor.node, neighbor.edge));
        }
        Self {
            root,
            included: distances.iter().map(Option::is_some).collect(),
            distances,
            predecessors,
        }
    }

    fn vismara(graph: &Graph, root: NodeId, component: &[NodeId], degrees: &[usize]) -> Self {
        let mut allowed = vec![false; graph.node_bound()];
        for &node in component {
            allowed[node.index()] = true;
        }
        let distances = bfs_distances(graph, root, &allowed);

        let root_key = (degrees[root.index()], root);
        let mut preceding = allowed.clone();
        for &node in component {
            preceding[node.index()] = node == root || (degrees[node.index()], node) < root_key;
        }
        let ordered_distances = bfs_distances(graph, root, &preceding);
        let mut included = vec![false; graph.node_bound()];
        for &node in component {
            included[node.index()] =
                node != root && ordered_distances[node.index()] == distances[node.index()];
        }
        included[root.index()] = true;

        let mut predecessors = vec![Vec::new(); graph.node_bound()];
        for &node in component {
            if node == root || !included[node.index()] {
                continue;
            }
            let distance = distances[node.index()].expect("included node has a distance");
            for neighbor in graph.neighbors(node) {
                if included[neighbor.node.index()]
                    && distances[neighbor.node.index()].is_some_and(|other| other + 1 == distance)
                {
                    predecessors[node.index()].push(*neighbor);
                }
            }
            predecessors[node.index()]
                .sort_unstable_by_key(|neighbor| (neighbor.node, neighbor.edge));
        }

        Self {
            root,
            distances,
            predecessors,
            included,
        }
    }

    pub(super) fn distance(&self, node: NodeId) -> Option<usize> {
        self.distances[node.index()]
    }

    pub(super) fn path_to(&self, target: NodeId) -> Option<ShortestPath> {
        if !self.included[target.index()] {
            return None;
        }
        let mut nodes = vec![target];
        let mut edges = Vec::new();
        let mut current = target;
        while current != self.root {
            let predecessor = self.predecessors[current.index()].first()?;
            nodes.push(predecessor.node);
            edges.push(predecessor.edge);
            current = predecessor.node;
        }
        nodes.reverse();
        edges.reverse();
        Some(ShortestPath { nodes, edges })
    }

    fn visit_paths_to<B>(
        &self,
        target: NodeId,
        mut visitor: impl FnMut(&ShortestPath) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        if !self.included[target.index()] {
            return ControlFlow::Continue(());
        }
        let mut nodes = vec![target];
        let mut edges = Vec::new();
        self.visit_reconstructed(target, &mut nodes, &mut edges, &mut visitor)
    }

    fn path_count(&self, target: NodeId) -> BigUint {
        if !self.included[target.index()] {
            return BigUint::from(0_u8);
        }

        let mut nodes = self
            .included
            .iter()
            .enumerate()
            .filter_map(|(node, &included)| included.then_some(NodeId(node as u32)))
            .collect::<Vec<_>>();
        nodes.sort_unstable_by_key(|node| self.distances[node.index()]);

        let mut counts = vec![BigUint::from(0_u8); self.included.len()];
        counts[self.root.index()] = BigUint::from(1_u8);
        for node in nodes {
            if node == self.root {
                continue;
            }
            counts[node.index()] = self.predecessors[node.index()]
                .iter()
                .map(|predecessor| &counts[predecessor.node.index()])
                .sum();
        }
        counts[target.index()].clone()
    }

    fn path_union(&self, target: NodeId) -> (Vec<NodeId>, Vec<EdgeId>) {
        if !self.included[target.index()] {
            return (Vec::new(), Vec::new());
        }

        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        let mut visited = vec![false; self.included.len()];
        let mut stack = vec![target];
        visited[target.index()] = true;
        while let Some(node) = stack.pop() {
            nodes.push(node);
            for predecessor in &self.predecessors[node.index()] {
                edges.push(predecessor.edge);
                if !visited[predecessor.node.index()] {
                    visited[predecessor.node.index()] = true;
                    stack.push(predecessor.node);
                }
            }
        }
        nodes.sort_unstable();
        edges.sort_unstable();
        edges.dedup();
        (nodes, edges)
    }

    fn visit_reconstructed<B>(
        &self,
        current: NodeId,
        nodes: &mut Vec<NodeId>,
        edges: &mut Vec<EdgeId>,
        visitor: &mut impl FnMut(&ShortestPath) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        if current == self.root {
            return visitor(&ShortestPath {
                nodes: nodes.iter().rev().copied().collect(),
                edges: edges.iter().rev().copied().collect(),
            });
        }

        for predecessor in &self.predecessors[current.index()] {
            nodes.push(predecessor.node);
            edges.push(predecessor.edge);
            if let ControlFlow::Break(value) =
                self.visit_reconstructed(predecessor.node, nodes, edges, visitor)
            {
                return ControlFlow::Break(value);
            }
            edges.pop();
            nodes.pop();
        }
        ControlFlow::Continue(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FamilyConnector {
    Odd(EdgeId),
    Even {
        middle: NodeId,
        first_edge: EdgeId,
        second_edge: EdgeId,
    },
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct RelevantCycleFamily {
    dag: usize,
    first: NodeId,
    second: NodeId,
    connector: FamilyConnector,
    prototype: Cycle,
}

impl RelevantCycleFamily {
    pub(super) fn weight(&self) -> usize {
        self.prototype.length()
    }

    pub(super) fn prototype(&self) -> &Cycle {
        &self.prototype
    }

    pub(super) fn cycle_count(&self, dag: &ShortestPathDag) -> BigUint {
        dag.path_count(self.first) * dag.path_count(self.second)
    }

    pub(super) fn union(&self, dag: &ShortestPathDag) -> (Vec<NodeId>, Vec<EdgeId>) {
        let (mut nodes, mut edges) = dag.path_union(self.first);
        let (other_nodes, other_edges) = dag.path_union(self.second);
        nodes.extend(other_nodes);
        edges.extend(other_edges);
        match self.connector {
            FamilyConnector::Odd(edge) => edges.push(edge),
            FamilyConnector::Even {
                middle,
                first_edge,
                second_edge,
            } => {
                nodes.push(middle);
                edges.extend([first_edge, second_edge]);
            }
        }
        nodes.sort_unstable();
        nodes.dedup();
        edges.sort_unstable();
        edges.dedup();
        (nodes, edges)
    }

    pub(super) fn visit_cycles<B>(
        &self,
        graph: &Graph,
        dag: &ShortestPathDag,
        visitor: &mut impl FnMut(Cycle) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        dag.visit_paths_to(self.first, |first| {
            dag.visit_paths_to(self.second, |second| {
                let cycle = match self.connector {
                    FamilyConnector::Odd(edge) => first.cycle_with(graph, second, edge),
                    FamilyConnector::Even {
                        middle,
                        first_edge,
                        second_edge,
                    } => first.cycle_with_middle(graph, second, middle, first_edge, second_edge),
                }
                .expect("paths in a Vismara family form a simple cycle");
                visitor(cycle)
            })
        })
    }
}

#[derive(Debug)]
pub(super) struct RelevantCycleAnalysis {
    dags: Vec<ShortestPathDag>,
    families: Vec<RelevantCycleFamily>,
}

pub(super) fn visit_relevant_cycles_vismara<B>(
    graph: &Graph,
    max_cycle_size: usize,
    visitor: &mut impl FnMut(Cycle) -> ControlFlow<B>,
) -> ControlFlow<B> {
    debug_assert!(graph.is_simple(), "direct Vismara input must be simple");
    if max_cycle_size < 3 {
        return ControlFlow::Continue(());
    }
    RelevantCycleAnalysis::new(graph).visit_cycles(graph, max_cycle_size, visitor)
}

pub(super) fn visit_relevant_cycles_vismara_fallback<B>(
    source: &Graph,
    max_cycle_size: usize,
    visitor: &mut impl FnMut(Cycle) -> ControlFlow<B>,
) -> ControlFlow<B> {
    let mut loops = Vec::new();
    let mut loopless_edges = Vec::new();
    let mut edge_sources = Vec::new();

    for edge in source.edge_ids() {
        let [first, second] = source.edge_endpoints(edge);
        if first == second {
            loops.push(Cycle::normalized(source, vec![first], vec![edge]));
            continue;
        }
        loopless_edges.push([first.0, second.0]);
        edge_sources.push(edge);
    }

    if max_cycle_size >= 1 {
        for cycle in loops {
            if let ControlFlow::Break(value) = visitor(cycle) {
                return ControlFlow::Break(value);
            }
        }
    }
    if max_cycle_size < 2 {
        return ControlFlow::Continue(());
    }

    let loopless = Graph::new(source.node_count(), &loopless_edges);
    let subdivision = loopless.subdivide_edges();
    RelevantCycleAnalysis::new(subdivision.graph()).visit_cycles(
        subdivision.graph(),
        max_cycle_size.saturating_mul(2),
        |cycle| visitor(cycle.map_subdivision(source, &subdivision, &edge_sources)),
    )
}

impl RelevantCycleAnalysis {
    pub(super) fn new(graph: &Graph) -> Self {
        let mut dags = Vec::new();
        let mut candidates = Vec::new();
        // Vismara analyzes relevant cycles independently within biconnected
        // components. Tarjan supplies that decomposition as fixed preprocessing,
        // not as an independent relevant-cycle choice.
        let mut components = graph.enumerate_biconnected_components(BiconnectedComponentsAlgorithm::Tarjan);
        components.sort();

        for component in components {
            let mut degrees = vec![0; graph.node_bound()];
            for &node in &component {
                degrees[node.index()] = graph
                    .neighbors(node)
                    .iter()
                    .filter(|neighbor| component.binary_search(&neighbor.node).is_ok())
                    .count();
            }

            for &root in &component {
                let dag = ShortestPathDag::vismara(graph, root, &component, &degrees);
                let dag_index = dags.len();

                for &middle in &component {
                    if middle == root || !dag.included[middle.index()] {
                        continue;
                    }

                    let first_path = dag
                        .path_to(middle)
                        .expect("included Vismara node has a representative path");
                    let mut previous = Vec::new();
                    for neighbor in graph.neighbors(middle) {
                        if !dag.included[neighbor.node.index()] {
                            continue;
                        }
                        let Some(neighbor_distance) = dag.distance(neighbor.node) else {
                            continue;
                        };
                        let middle_distance = dag
                            .distance(middle)
                            .expect("included Vismara node has a distance");

                        if neighbor_distance + 1 == middle_distance {
                            previous.push(*neighbor);
                        } else if neighbor_distance == middle_distance
                            && (degrees[neighbor.node.index()], neighbor.node)
                                < (degrees[middle.index()], middle)
                        {
                            let second_path = dag
                                .path_to(neighbor.node)
                                .expect("included Vismara neighbor has a representative path");
                            if first_path.common_prefix_len(&second_path) == 1 {
                                candidates.push(RelevantCycleFamily {
                                    dag: dag_index,
                                    first: middle,
                                    second: neighbor.node,
                                    connector: FamilyConnector::Odd(neighbor.edge),
                                    prototype: first_path
                                        .cycle_with(graph, &second_path, neighbor.edge)
                                        .expect("Vismara odd prototype is a simple cycle"),
                                });
                            }
                        }
                    }

                    for first in 0..previous.len() {
                        for second in first + 1..previous.len() {
                            let first_path = dag
                                .path_to(previous[first].node)
                                .expect("included Vismara predecessor has a representative path");
                            let second_path = dag
                                .path_to(previous[second].node)
                                .expect("included Vismara predecessor has a representative path");
                            if first_path.common_prefix_len(&second_path) == 1 {
                                candidates.push(RelevantCycleFamily {
                                    dag: dag_index,
                                    first: previous[first].node,
                                    second: previous[second].node,
                                    connector: FamilyConnector::Even {
                                        middle,
                                        first_edge: previous[first].edge,
                                        second_edge: previous[second].edge,
                                    },
                                    prototype: first_path
                                        .cycle_with_middle(
                                            graph,
                                            &second_path,
                                            middle,
                                            previous[first].edge,
                                            previous[second].edge,
                                        )
                                        .expect("Vismara even prototype is a simple cycle"),
                                });
                            }
                        }
                    }
                }

                dags.push(dag);
            }
        }

        candidates.sort_by(compare_families);
        let mut basis = CycleVectorBasis::new(graph.edge_count());
        let mut families = Vec::new();
        let mut start = 0;
        while start < candidates.len() {
            let weight = candidates[start].weight();
            let end = candidates[start..]
                .iter()
                .position(|family| family.weight() != weight)
                .map_or(candidates.len(), |offset| start + offset);
            let relevant = candidates[start..end]
                .iter()
                .filter(|family| {
                    basis.is_independent(EdgeVector::from_cycle(
                        graph.edge_count(),
                        &family.prototype,
                    ))
                })
                .cloned()
                .collect::<Vec<_>>();
            for family in &relevant {
                basis.insert(EdgeVector::from_cycle(
                    graph.edge_count(),
                    &family.prototype,
                ));
            }
            families.extend(relevant);
            start = end;
        }

        Self { dags, families }
    }

    pub(super) fn visit_cycles<B>(
        &self,
        graph: &Graph,
        max_cycle_size: usize,
        mut visitor: impl FnMut(Cycle) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        for family in &self.families {
            if family.weight() > max_cycle_size {
                continue;
            }
            if let ControlFlow::Break(value) =
                family.visit_cycles(graph, &self.dags[family.dag], &mut visitor)
            {
                return ControlFlow::Break(value);
            }
        }
        ControlFlow::Continue(())
    }

    pub(super) fn families(&self) -> &[RelevantCycleFamily] {
        &self.families
    }

    pub(super) fn family_dag(&self, family: &RelevantCycleFamily) -> &ShortestPathDag {
        &self.dags[family.dag]
    }
}

fn bfs_distances(graph: &Graph, root: NodeId, allowed: &[bool]) -> Vec<Option<usize>> {
    let mut distances = vec![None; graph.node_bound()];
    let mut queue = VecDeque::from([root]);
    distances[root.index()] = Some(0);
    while let Some(current) = queue.pop_front() {
        let next_distance = distances[current.index()].expect("queued node has a distance") + 1;
        for neighbor in graph.neighbors(current) {
            if allowed[neighbor.node.index()] && distances[neighbor.node.index()].is_none() {
                distances[neighbor.node.index()] = Some(next_distance);
                queue.push_back(neighbor.node);
            }
        }
    }
    distances
}

fn compare_families(first: &RelevantCycleFamily, second: &RelevantCycleFamily) -> Ordering {
    first
        .weight()
        .cmp(&second.weight())
        .then_with(|| first.prototype.nodes().cmp(second.prototype.nodes()))
        .then_with(|| first.prototype.edges().cmp(second.prototype.edges()))
        .then_with(|| first.dag.cmp(&second.dag))
        .then_with(|| first.first.cmp(&second.first))
        .then_with(|| first.second.cmp(&second.second))
}

#[cfg(test)]
mod tests {
    use std::ops::ControlFlow;

    use rstest::rstest;

    use super::{RelevantCycleAnalysis, ShortestPath, ShortestPathDag};
    use crate::algorithms::cycles::Cycle;
    use crate::graph::{EdgeId, Graph, NodeId};

    #[rstest]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(1)],
            edges: vec![EdgeId(0)],
        },
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(2)],
            edges: vec![EdgeId(2)],
        },
        EdgeId(1),
        Some(Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }),
    )]
    #[case::shared_prefix(
        Graph::new(4, &[[0, 1], [1, 2], [1, 3], [2, 3]]),
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1)],
        },
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(1), NodeId(3)],
            edges: vec![EdgeId(0), EdgeId(2)],
        },
        EdgeId(3),
        Some(Cycle {
            nodes: vec![NodeId(1), NodeId(2), NodeId(3)],
            edges: vec![EdgeId(1), EdgeId(3), EdgeId(2)],
        }),
    )]
    #[case::intersecting_tails(
        Graph::new(5, &[[0, 1], [1, 2], [0, 3], [3, 1], [1, 4], [2, 4]]),
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1)],
        },
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(3), NodeId(1), NodeId(4)],
            edges: vec![EdgeId(2), EdgeId(3), EdgeId(4)],
        },
        EdgeId(5),
        None,
    )]
    fn test_shortest_path_cycle_with(
        #[case] graph: Graph,
        #[case] path: ShortestPath,
        #[case] other: ShortestPath,
        #[case] closing_edge: EdgeId,
        #[case] expected: Option<Cycle>,
    ) {
        assert_eq!(path.cycle_with(&graph, &other, closing_edge), expected);
    }

    #[rstest]
    #[case::square(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [0, 3]]),
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(1)],
            edges: vec![EdgeId(0)],
        },
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(3)],
            edges: vec![EdgeId(3)],
        },
        NodeId(2),
        EdgeId(1),
        EdgeId(2),
        Some(Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
        }),
    )]
    fn test_shortest_path_cycle_with_middle(
        #[case] graph: Graph,
        #[case] path: ShortestPath,
        #[case] other: ShortestPath,
        #[case] middle: NodeId,
        #[case] first_edge: EdgeId,
        #[case] second_edge: EdgeId,
        #[case] expected: Option<Cycle>,
    ) {
        assert_eq!(
            path.cycle_with_middle(&graph, &other, middle, first_edge, second_edge),
            expected
        );
    }

    #[rstest]
    #[case::connected(
        Graph::new(5, &[[0, 1], [0, 2], [1, 3], [2, 3]]),
        NodeId(0),
        vec![Some(0), Some(1), Some(1), Some(2), None],
    )]
    #[case::parallel(
        Graph::new(2, &[[0, 1], [0, 1]]),
        NodeId(0),
        vec![Some(0), Some(1)],
    )]
    fn test_shortest_path_dag_distance(
        #[case] graph: Graph,
        #[case] root: NodeId,
        #[case] expected: Vec<Option<usize>>,
    ) {
        let dag = ShortestPathDag::new(&graph, root);
        let actual = graph
            .node_ids()
            .map(|node| dag.distance(node))
            .collect::<Vec<_>>();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::selected(
        Graph::new(4, &[[0, 1], [0, 2], [1, 3], [2, 3]]),
        NodeId(0),
        NodeId(3),
        Some(ShortestPath {
            nodes: vec![NodeId(0), NodeId(1), NodeId(3)],
            edges: vec![EdgeId(0), EdgeId(2)],
        }),
    )]
    #[case::unreachable(Graph::new(2, &[]), NodeId(0), NodeId(1), None)]
    fn test_shortest_path_dag_path_to(
        #[case] graph: Graph,
        #[case] root: NodeId,
        #[case] target: NodeId,
        #[case] expected: Option<ShortestPath>,
    ) {
        assert_eq!(ShortestPathDag::new(&graph, root).path_to(target), expected);
    }

    #[rstest]
    #[case::alternatives(
        Graph::new(4, &[[0, 1], [0, 2], [1, 3], [2, 3]]),
        NodeId(0),
        NodeId(3),
        vec![
            ShortestPath {
                nodes: vec![NodeId(0), NodeId(1), NodeId(3)],
                edges: vec![EdgeId(0), EdgeId(2)],
            },
            ShortestPath {
                nodes: vec![NodeId(0), NodeId(2), NodeId(3)],
                edges: vec![EdgeId(1), EdgeId(3)],
            },
        ],
    )]
    #[case::parallel(
        Graph::new(2, &[[0, 1], [0, 1]]),
        NodeId(0),
        NodeId(1),
        vec![
            ShortestPath {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0)],
            },
            ShortestPath {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1)],
            },
        ],
    )]
    #[case::unreachable(
        Graph::new(2, &[]),
        NodeId(0),
        NodeId(1),
        vec![],
    )]
    fn test_shortest_path_dag_visit_paths_to(
        #[case] graph: Graph,
        #[case] root: NodeId,
        #[case] target: NodeId,
        #[case] expected: Vec<ShortestPath>,
    ) {
        let mut paths = Vec::new();
        let result = ShortestPathDag::new(&graph, root).visit_paths_to(target, |path| {
            paths.push(path.clone());
            ControlFlow::<()>::Continue(())
        });
        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(paths, expected);
    }

    #[rstest]
    #[case::first(
        1,
        vec![ShortestPath {
            nodes: vec![NodeId(0), NodeId(1), NodeId(3)],
            edges: vec![EdgeId(0), EdgeId(2)],
        }],
    )]
    fn test_shortest_path_dag_visit_paths_to_break(
        #[case] stop_after: usize,
        #[case] expected: Vec<ShortestPath>,
    ) {
        let graph = Graph::new(4, &[[0, 1], [0, 2], [1, 3], [2, 3]]);
        let mut paths = Vec::new();
        let result = ShortestPathDag::new(&graph, NodeId(0)).visit_paths_to(NodeId(3), |path| {
            paths.push(path.clone());
            if paths.len() == stop_after {
                ControlFlow::Break(paths.len())
            } else {
                ControlFlow::Continue(())
            }
        });
        assert_eq!(result, ControlFlow::Break(stop_after));
        assert_eq!(paths, expected);
    }

    #[rstest]
    #[case::odd(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }],
    )]
    #[case::even(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [0, 3]]),
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
        }],
    )]
    #[case::unequal_theta(Graph::new(
        6,
        &[
            [0, 1], [1, 4], [0, 2], [2, 4], [0, 3], [3, 5], [4, 5],
        ],
    ),
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(4), NodeId(2)],
                edges: vec![EdgeId(0), EdgeId(1), EdgeId(3), EdgeId(2)],
            },
            Cycle {
                nodes: vec![
                    NodeId(0), NodeId(1), NodeId(4), NodeId(5), NodeId(3),
                ],
                edges: vec![
                    EdgeId(0), EdgeId(1), EdgeId(6), EdgeId(5), EdgeId(4),
                ],
            },
        ],
    )]
    #[case::fused(
        Graph::new(
            10,
            &[
                [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [0, 5],
                [3, 6], [6, 7], [7, 8], [8, 9], [4, 9],
            ],
        ),
        vec![
            Cycle {
                nodes: vec![
                    NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4), NodeId(5),
                ],
                edges: vec![
                    EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3), EdgeId(4), EdgeId(5),
                ],
            },
            Cycle {
                nodes: vec![
                    NodeId(3), NodeId(4), NodeId(9), NodeId(8), NodeId(7), NodeId(6),
                ],
                edges: vec![
                    EdgeId(3), EdgeId(10), EdgeId(9), EdgeId(8), EdgeId(7), EdgeId(6),
                ],
            },
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
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
                edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            },
            Cycle {
                nodes: vec![NodeId(3), NodeId(4), NodeId(5)],
                edges: vec![EdgeId(3), EdgeId(4), EdgeId(5)],
            },
        ],
    )]
    fn test_relevant_cycle_analysis_new(#[case] graph: Graph, #[case] expected: Vec<Cycle>) {
        let analysis = RelevantCycleAnalysis::new(&graph);
        let prototypes = analysis
            .families
            .iter()
            .map(|family| family.prototype.clone())
            .collect::<Vec<_>>();
        assert_eq!(prototypes, expected);
    }

    #[rstest]
    #[case::multiple_shortest_paths(
        Graph::new(
            6,
            &[
                [0, 1], [1, 4], [0, 2], [2, 4], [0, 3], [3, 5], [4, 5],
            ],
        ),
        vec![
            Cycle {
                nodes: vec![NodeId(0), NodeId(1), NodeId(4), NodeId(2)],
                edges: vec![EdgeId(0), EdgeId(1), EdgeId(3), EdgeId(2)],
            },
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
    )]
    fn test_relevant_cycle_analysis_visit_cycles(
        #[case] graph: Graph,
        #[case] expected: Vec<Cycle>,
    ) {
        let analysis = RelevantCycleAnalysis::new(&graph);
        let mut cycles = Vec::new();
        let result = analysis.visit_cycles(&graph, usize::MAX, |cycle| {
            cycles.push(cycle);
            ControlFlow::<()>::Continue(())
        });
        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(cycles, expected);
    }

    #[rstest]
    #[case::first(
        1,
        vec![Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(4), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(3), EdgeId(2)],
        }],
    )]
    fn test_relevant_cycle_analysis_visit_cycles_break(
        #[case] stop_after: usize,
        #[case] expected: Vec<Cycle>,
    ) {
        let graph = Graph::new(6, &[[0, 1], [1, 4], [0, 2], [2, 4], [0, 3], [3, 5], [4, 5]]);
        let analysis = RelevantCycleAnalysis::new(&graph);
        let mut cycles = Vec::new();
        let result = analysis.visit_cycles(&graph, usize::MAX, |cycle| {
            cycles.push(cycle);
            if cycles.len() == stop_after {
                ControlFlow::Break(cycles.len())
            } else {
                ControlFlow::Continue(())
            }
        });
        assert_eq!(result, ControlFlow::Break(stop_after));
        assert_eq!(cycles, expected);
    }
}
