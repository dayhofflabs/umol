//! Maximum and perfect matching.

use std::collections::{HashSet, VecDeque};
use std::error::Error;
use std::fmt;

#[cfg(test)]
use super::enumeration::{MatchingEnumerationAlgorithm, MatchingSearchState};
use crate::algorithms::bipartite::BipartitionAlgorithm;
use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaximumMatchingAlgorithm {
    Edmonds,
    /// Bipartite-only.
    HopcroftKarp,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaximumMatchingError {
    NonBipartite,
}

impl fmt::Display for MaximumMatchingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonBipartite => write!(formatter, "Hopcroft-Karp requires a bipartite graph"),
        }
    }
}

impl Error for MaximumMatchingError {}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PerfectMatchingAlgorithm {
    BacktrackingDfs,
}

/// A matching: a set of edges with no shared endpoints.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matching {
    pub(super) edges: Vec<EdgeId>,
    pub(super) mate: Vec<Option<NodeId>>,
}

impl Matching {
    fn from_mates(graph: &Graph, mates: &[i32]) -> Self {
        let mut edges = Vec::new();
        let mut represented = vec![false; graph.node_bound()];
        for eid in graph.edge_ids() {
            let [a, b] = graph.edge_endpoints(eid);
            if mates[a.index()] == b.0 as i32 && !represented[a.index()] && !represented[b.index()]
            {
                edges.push(eid);
                represented[a.index()] = true;
                represented[b.index()] = true;
            }
        }
        edges.sort_unstable();
        let mate_opt = mates
            .iter()
            .map(|&m| if m >= 0 { Some(NodeId(m as u32)) } else { None })
            .collect();
        Self {
            edges,
            mate: mate_opt,
        }
    }

    pub fn edges(&self) -> &[EdgeId] {
        &self.edges
    }

    pub fn size(&self) -> usize {
        self.edges.len()
    }

    pub fn is_perfect(&self, node_count: usize) -> bool {
        self.edges.len() * 2 == node_count
    }

    pub fn mate(&self, node: NodeId) -> Option<NodeId> {
        self.mate.get(node.index()).copied().flatten()
    }

    pub fn is_matched(&self, node: NodeId) -> bool {
        self.mate(node).is_some()
    }
}

impl Graph {
    /// Returns a maximum matching, using `node_order` as the deterministic vertex traversal
    /// priority. `node_order` must contain every graph node exactly once. Neighbor ties retain the
    /// graph's stable adjacency order. Hopcroft-Karp returns
    /// [`MaximumMatchingError::NonBipartite`] when the graph is not bipartite.
    pub fn maximum_matching(
        &self,
        node_order: &[NodeId],
        alg: MaximumMatchingAlgorithm,
    ) -> Result<Matching, MaximumMatchingError> {
        debug_assert_eq!(
            node_order.len(),
            self.node_count(),
            "maximum-matching node order must contain every graph node exactly once",
        );
        debug_assert_eq!(
            node_order.iter().copied().collect::<HashSet<_>>(),
            self.node_ids().collect::<HashSet<_>>(),
            "maximum-matching node order must contain every graph node exactly once",
        );
        match alg {
            MaximumMatchingAlgorithm::Edmonds => Ok(self.maximum_matching_edmonds(node_order)),
            MaximumMatchingAlgorithm::HopcroftKarp => {
                self.maximum_matching_hopcroft_karp(node_order)
            }
        }
    }

    /// Returns a perfect matching (covers every vertex) when one exists in the
    /// node order given, else `None`. The order controls determinism: the same
    /// graph + same order always yields the same matching.
    pub fn perfect_matching(
        &self,
        node_order: &[NodeId],
        alg: PerfectMatchingAlgorithm,
    ) -> Option<Matching> {
        match alg {
            PerfectMatchingAlgorithm::BacktrackingDfs => {
                self.perfect_matching_backtracking_dfs(node_order)
            }
        }
    }

    // Hopcroft-Karp maximum matching. BFS constructs shortest augmenting-path
    // layers and the following DFS pass augments a maximal set of disjoint
    // paths in that layer graph, giving O(E·√V).
    fn maximum_matching_hopcroft_karp(
        &self,
        node_order: &[NodeId],
    ) -> Result<Matching, MaximumMatchingError> {
        let n = self.node_count();
        // Hopcroft-Karp needs a two-coloring before it can build alternating
        // layers. BFS supplies that coloring in O(V+E) as fixed preprocessing
        // for this implementation, not as a separate matching choice.
        let colors = self
            .bipartition(BipartitionAlgorithm::Bfs)
            .ok_or(MaximumMatchingError::NonBipartite)?;

        let mut mate = vec![-1i32; n];
        let mut distance = vec![usize::MAX; n];
        loop {
            let shortest = hopcroft_karp_bfs(self, node_order, &colors, &mate, &mut distance);
            if shortest == usize::MAX {
                break;
            }
            for &node in node_order {
                if !colors[node.index()] && mate[node.index()] < 0 {
                    hopcroft_karp_dfs(self, node, &mut mate, &mut distance, shortest);
                }
            }
        }

        Ok(Matching::from_mates(self, &mate))
    }

    // Greedy DFS with backtracking. Walks `node_order` left-to-right; at each
    // unmatched vertex tries each unmatched neighbor in adjacency order;
    // backtracks on dead ends. Returns `Some` if every vertex in `node_order`
    // can be paired, else `None`. Deterministic given a stable `node_order`
    // and the graph's stable neighbor list.
    fn perfect_matching_backtracking_dfs(&self, node_order: &[NodeId]) -> Option<Matching> {
        let n = self.node_count();
        if node_order.is_empty() {
            return Some(Matching {
                edges: Vec::new(),
                mate: Vec::new(),
            });
        }
        if !node_order.len().is_multiple_of(2) {
            return None;
        }
        let mut mate = vec![-1i32; n];
        if backtrack_pair(self, node_order, 0, &mut mate) {
            Some(Matching::from_mates(self, &mate))
        } else {
            None
        }
    }

    // Edmonds 1965, Gabow simplification 1976. Ref impl: cp-algorithms.com.
    pub(super) fn maximum_matching_edmonds(&self, node_order: &[NodeId]) -> Matching {
        let n = self.node_count();
        if n == 0 {
            return Matching {
                edges: Vec::new(),
                mate: Vec::new(),
            };
        }

        let mut mate = vec![-1i32; n];

        for &first in node_order {
            if mate[first.index()] >= 0 {
                continue;
            }
            if let Some(second) = self
                .neighbors(first)
                .iter()
                .map(|neighbor| neighbor.node)
                .find(|second| mate[second.index()] < 0)
            {
                mate[first.index()] = second.0 as i32;
                mate[second.index()] = first.0 as i32;
            }
        }

        for &node in node_order {
            if mate[node.index()] < 0 {
                augment_from(self, &mut mate, node.index());
            }
        }

        Matching::from_mates(self, &mate)
    }
}

fn hopcroft_karp_bfs(
    graph: &Graph,
    node_order: &[NodeId],
    colors: &[bool],
    mate: &[i32],
    distance: &mut [usize],
) -> usize {
    distance.fill(usize::MAX);
    let mut queue = VecDeque::new();
    for &node in node_order {
        if !colors[node.index()] && mate[node.index()] < 0 {
            distance[node.index()] = 0;
            queue.push_back(node);
        }
    }

    let mut shortest = usize::MAX;
    while let Some(u) = queue.pop_front() {
        let next_distance = distance[u.index()] + 1;
        if next_distance > shortest {
            continue;
        }
        for nbr in graph.neighbors(u) {
            let v = nbr.node;
            let paired = mate[v.index()];
            if paired < 0 {
                shortest = shortest.min(next_distance);
            } else if next_distance < shortest && distance[paired as usize] == usize::MAX {
                distance[paired as usize] = next_distance;
                queue.push_back(NodeId(paired as u32));
            }
        }
    }
    shortest
}

fn hopcroft_karp_dfs(
    graph: &Graph,
    u: NodeId,
    mate: &mut [i32],
    distance: &mut [usize],
    shortest: usize,
) -> bool {
    let next_distance = distance[u.index()] + 1;
    for nbr in graph.neighbors(u) {
        let v = nbr.node;
        let paired = mate[v.index()];
        let extends_path = if paired < 0 {
            next_distance == shortest
        } else {
            next_distance < shortest
                && distance[paired as usize] == next_distance
                && hopcroft_karp_dfs(graph, NodeId(paired as u32), mate, distance, shortest)
        };
        if extends_path {
            mate[u.index()] = v.0 as i32;
            mate[v.index()] = u.0 as i32;
            return true;
        }
    }
    distance[u.index()] = usize::MAX;
    false
}

// Backtracking helper for `perfect_matching_backtracking_dfs`. Returns true
// when every vertex in `order[idx..]` has been paired (recursing into
// `idx + 1` after each successful pair).
fn backtrack_pair(graph: &Graph, order: &[NodeId], idx: usize, mate: &mut [i32]) -> bool {
    let mut i = idx;
    while i < order.len() && mate[order[i].index()] >= 0 {
        i += 1;
    }
    if i >= order.len() {
        return true;
    }
    let v = order[i];
    for nbr in graph.neighbors(v) {
        let u = nbr.node;
        if mate[u.index()] >= 0 {
            continue;
        }
        mate[v.index()] = u.0 as i32;
        mate[u.index()] = v.0 as i32;
        if backtrack_pair(graph, order, i + 1, mate) {
            return true;
        }
        mate[v.index()] = -1;
        mate[u.index()] = -1;
    }
    false
}

fn augment_from(graph: &Graph, mate: &mut [i32], root: usize) -> bool {
    let n = graph.node_bound();
    let mut used = vec![false; n];
    let mut p = vec![-1i32; n];
    let mut base: Vec<usize> = (0..n).collect();

    used[root] = true;
    let mut queue = VecDeque::new();
    queue.push_back(root);

    while let Some(v) = queue.pop_front() {
        for nbr in graph.neighbors(NodeId(v as u32)) {
            let to = nbr.node.index();

            if base[v] == base[to] || mate[v] == to as i32 {
                continue;
            }

            if to == root || (mate[to] >= 0 && p[mate[to] as usize] >= 0) {
                let curbase = blossom_lca(&base, mate, &p, root, v, to);
                let mut in_blossom = vec![false; n];
                mark_blossom_path(&mut in_blossom, &base, mate, &mut p, v, curbase, to);
                mark_blossom_path(&mut in_blossom, &base, mate, &mut p, to, curbase, v);
                for i in 0..n {
                    if in_blossom[base[i]] {
                        base[i] = curbase;
                        if !used[i] {
                            used[i] = true;
                            queue.push_back(i);
                        }
                    }
                }
            } else if p[to] < 0 {
                p[to] = v as i32;
                if mate[to] < 0 {
                    let mut u = to as i32;
                    while u >= 0 {
                        let pv = p[u as usize];
                        let ppv = mate[pv as usize];
                        mate[u as usize] = pv;
                        mate[pv as usize] = u;
                        u = ppv;
                    }
                    return true;
                }
                let matched = mate[to] as usize;
                used[matched] = true;
                queue.push_back(matched);
            }
        }
    }

    false
}

fn blossom_lca(base: &[usize], mate: &[i32], p: &[i32], root: usize, a: usize, b: usize) -> usize {
    let n = base.len();
    let mut visited = vec![false; n];
    let mut a = base[a];
    let mut b = base[b];
    loop {
        visited[a] = true;
        if a == root || mate[a] < 0 {
            break;
        }
        a = base[p[mate[a] as usize] as usize];
    }
    loop {
        if visited[b] {
            return b;
        }
        b = base[p[mate[b] as usize] as usize];
    }
}

fn mark_blossom_path(
    in_blossom: &mut [bool],
    base: &[usize],
    mate: &[i32],
    p: &mut [i32],
    mut v: usize,
    b: usize,
    mut child: usize,
) {
    while base[v] != b {
        in_blossom[base[v]] = true;
        in_blossom[base[mate[v] as usize]] = true;
        p[v] = child as i32;
        child = mate[v] as usize;
        v = p[mate[v] as usize] as usize;
    }
}

#[cfg(test)]
mod tests {
    use std::ops::ControlFlow;
    use std::time::{Duration, Instant};

    use pretty_assertions::assert_eq;
    use proptest::prelude::*;
    use proptest::test_runner::{Config, TestRunner};
    use rstest::*;

    use super::MatchingEnumerationAlgorithm::BranchAndBound;
    use super::MaximumMatchingAlgorithm::{Edmonds, HopcroftKarp};
    use super::PerfectMatchingAlgorithm::BacktrackingDfs;
    use super::{Matching, MatchingSearchState, MaximumMatchingAlgorithm, MaximumMatchingError};
    use crate::graph::{EdgeId, Graph, NodeId};

    #[rstest]
    #[case::empty(0, vec![], 0, true)]
    #[case::single_edge(2, vec![[0, 1]], 1, true)]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], 1, false)]
    #[case::square(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]], 2, true)]
    #[case::path_4(4, vec![[0, 1], [1, 2], [2, 3]], 2, true)]
    #[case::path_5(5, vec![[0, 1], [1, 2], [2, 3], [3, 4]], 2, false)]
    #[case::k4(4, vec![[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], 2, true)]
    #[case::parallel_edges(2, vec![[0, 1], [0, 1]], 1, true)]
    #[case::petersen(
        10,
        vec![
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 0],
            [5, 7], [7, 9], [9, 6], [6, 8], [8, 5],
            [0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
        ],
        5,
        true,
    )]
    fn test_graph_maximum_matching(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected_size: usize,
        #[case] expected_perfect: bool,
    ) {
        let g = Graph::new(node_count, &edges);
        let node_order: Vec<NodeId> = g.node_ids().collect();
        let m = g.maximum_matching(&node_order, Edmonds).unwrap();
        assert_eq!(m.size(), expected_size, "matching size");
        assert_eq!(m.is_perfect(node_count), expected_perfect, "is_perfect");
        assert_matching_valid(&g, &m);
    }

    #[rstest]
    #[case::single_edge(2, vec![[0, 1]], 1, true)]
    #[case::square(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]], 2, true)]
    #[case::path_4(4, vec![[0, 1], [1, 2], [2, 3]], 2, true)]
    #[case::path_5(5, vec![[0, 1], [1, 2], [2, 3], [3, 4]], 2, false)]
    #[case::hexagon(6, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]], 3, true)]
    #[case::k_2_3(
        5,
        vec![[0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4]],
        2,
        false,
    )]
    fn test_graph_maximum_matching_hopcroft_karp(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected_size: usize,
        #[case] expected_perfect: bool,
    ) {
        let g = Graph::new(node_count, &edges);
        let node_order: Vec<NodeId> = g.node_ids().collect();
        let m = g.maximum_matching(&node_order, HopcroftKarp).unwrap();
        assert_eq!(m.size(), expected_size, "matching size");
        assert_eq!(m.is_perfect(node_count), expected_perfect, "is_perfect");
        assert_matching_valid(&g, &m);
    }

    #[rstest]
    #[case::square(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]])]
    #[case::hexagon(6, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]])]
    #[case::bipartite_4_by_4(
        8,
        vec![
            [0, 4], [0, 5], [1, 4], [1, 6],
            [2, 5], [2, 7], [3, 6], [3, 7],
        ],
    )]
    fn test_graph_maximum_matching_cross_algorithm(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
    ) {
        let graph = Graph::new(node_count, &edges);
        let node_order: Vec<NodeId> = graph.node_ids().collect();
        let hopcroft_karp = graph.maximum_matching(&node_order, HopcroftKarp).unwrap();
        let edmonds = graph.maximum_matching(&node_order, Edmonds).unwrap();
        assert_eq!(hopcroft_karp.size(), edmonds.size());
    }

    #[rstest]
    fn test_graph_maximum_matching_property() {
        const PROPERTY_CASES: u32 = 128;

        let strategy = (0_usize..=7, 0_usize..=7).prop_flat_map(|(left, right)| {
            let possible_edges: Vec<_> = (0..left as u32)
                .flat_map(|first| {
                    (left as u32..(left + right) as u32).map(move |second| [first, second])
                })
                .collect();
            (
                Just(left + right),
                Just(possible_edges.clone()),
                prop::collection::vec(any::<bool>(), possible_edges.len()),
            )
        });
        let mut runner = TestRunner::new(Config {
            cases: PROPERTY_CASES,
            ..Config::default()
        });

        runner
            .run(&strategy, |(node_count, possible_edges, present)| {
                let edges: Vec<_> = possible_edges
                    .into_iter()
                    .zip(present)
                    .filter_map(|(edge, present)| present.then_some(edge))
                    .collect();
                let graph = Graph::new(node_count, &edges);
                let node_order: Vec<NodeId> = graph.node_ids().collect();
                let hopcroft_karp = graph
                    .maximum_matching(&node_order, HopcroftKarp)
                    .expect("generated graph is bipartite");
                let edmonds = graph
                    .maximum_matching(&node_order, Edmonds)
                    .expect("Edmonds maximum matching is infallible");

                prop_assert_eq!(hopcroft_karp.size(), edmonds.size());
                Ok(())
            })
            .unwrap();
    }

    #[rstest]
    #[case::triangle(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]))]
    fn test_graph_maximum_matching_error(#[case] graph: Graph) {
        let node_order: Vec<NodeId> = graph.node_ids().collect();
        assert_eq!(
            graph.maximum_matching(&node_order, HopcroftKarp),
            Err(MaximumMatchingError::NonBipartite),
        );
    }

    #[rstest]
    #[case::edmonds(Edmonds)]
    #[case::hopcroft_karp(HopcroftKarp)]
    fn test_graph_maximum_matching_node_order(#[case] algorithm: MaximumMatchingAlgorithm) {
        let graph = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]);
        let first_order = [NodeId(0), NodeId(1), NodeId(2), NodeId(3)];
        let second_order = [NodeId(2), NodeId(0), NodeId(1), NodeId(3)];

        let first = graph.maximum_matching(&first_order, algorithm).unwrap();
        let second = graph.maximum_matching(&second_order, algorithm).unwrap();

        assert_eq!(first.edges(), &[EdgeId(0), EdgeId(2)]);
        assert_eq!(second.edges(), &[EdgeId(1), EdgeId(3)]);
        assert_eq!(graph.maximum_matching(&first_order, algorithm), Ok(first),);
        assert_eq!(graph.maximum_matching(&second_order, algorithm), Ok(second),);
    }

    #[rstest]
    #[case::single_edge(2, vec![[0, 1]], vec![NodeId(0), NodeId(1)], true, 1)]
    #[case::square(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]], vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)], true, 2)]
    #[case::hexagon(
        6,
        vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]],
        vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4), NodeId(5)],
        true,
        3,
    )]
    #[case::k4(
        4,
        vec![[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]],
        vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
        true,
        2,
    )]
    #[case::triangle_no_perfect(
        3,
        vec![[0, 1], [1, 2], [0, 2]],
        vec![NodeId(0), NodeId(1), NodeId(2)],
        false,
        0,
    )]
    #[case::path_5_no_perfect(
        5,
        vec![[0, 1], [1, 2], [2, 3], [3, 4]],
        vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4)],
        false,
        0,
    )]
    #[case::two_disconnected_edges(4, vec![[0, 1], [2, 3]], vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)], true, 2)]
    fn test_graph_perfect_matching(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] node_order: Vec<NodeId>,
        #[case] expected_some: bool,
        #[case] expected_size: usize,
    ) {
        let g = Graph::new(node_count, &edges);
        let m = g.perfect_matching(&node_order, BacktrackingDfs);
        assert_eq!(m.is_some(), expected_some, "Some-ness");
        if let Some(m) = m {
            assert_eq!(m.size(), expected_size);
            assert!(m.is_perfect(node_count));
            assert_matching_valid(&g, &m);
        }
    }

    #[rstest]
    fn test_graph_perfect_matching_determinism() {
        let g = Graph::new(4, &[[0, 1], [0, 2], [1, 3], [2, 3]]);
        let order = vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)];
        let m1 = g.perfect_matching(&order, BacktrackingDfs).unwrap();
        let m2 = g.perfect_matching(&order, BacktrackingDfs).unwrap();
        assert_eq!(m1.edges(), m2.edges());
    }

    #[rstest]
    #[case::single_edge(2, vec![[0, 1]], 1)]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], 0)]
    #[case::square(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]], 2)]
    #[case::c6(6, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]], 2)]
    #[case::k4(4, vec![[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], 3)]
    #[case::greedy_counterexample(4, vec![[1, 2], [0, 1], [2, 3]], 1)]
    #[case::parallel_edges(2, vec![[0, 1], [0, 1]], 2)]
    fn test_graph_enumerate_perfect_matchings(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected_count: usize,
    ) {
        let g = Graph::new(node_count, &edges);
        let matchings = g.enumerate_perfect_matchings(BranchAndBound);
        assert_eq!(
            matchings.len(),
            expected_count,
            "number of perfect matchings"
        );
        for m in &matchings {
            assert!(m.is_perfect(node_count));
            assert_matching_valid(&g, m);
        }
        assert_all_distinct(&matchings);
    }

    #[rstest]
    fn test_graph_enumerate_perfect_matchings_traversal() {
        let graph = Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]);

        assert_eq!(
            graph.enumerate_perfect_matchings(BranchAndBound),
            vec![
                Matching {
                    edges: vec![EdgeId(0), EdgeId(5)],
                    mate: vec![
                        Some(NodeId(1)),
                        Some(NodeId(0)),
                        Some(NodeId(3)),
                        Some(NodeId(2)),
                    ],
                },
                Matching {
                    edges: vec![EdgeId(1), EdgeId(4)],
                    mate: vec![
                        Some(NodeId(2)),
                        Some(NodeId(3)),
                        Some(NodeId(0)),
                        Some(NodeId(1)),
                    ],
                },
                Matching {
                    edges: vec![EdgeId(2), EdgeId(3)],
                    mate: vec![
                        Some(NodeId(3)),
                        Some(NodeId(2)),
                        Some(NodeId(1)),
                        Some(NodeId(0)),
                    ],
                },
            ]
        );
    }

    #[rstest]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], 3)]
    #[case::path_3(3, vec![[0, 1], [1, 2]], 2)]
    #[case::single_edge(2, vec![[0, 1]], 1)]
    #[case::parallel_edges(2, vec![[0, 1], [0, 1]], 2)]
    fn test_graph_enumerate_maximum_matchings(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected_count: usize,
    ) {
        let g = Graph::new(node_count, &edges);
        let node_order: Vec<NodeId> = g.node_ids().collect();
        let initial = g.maximum_matching(&node_order, Edmonds).unwrap();
        let target = initial.size();
        let matchings = g.enumerate_maximum_matchings(BranchAndBound);
        assert_eq!(
            matchings.len(),
            expected_count,
            "number of maximum matchings"
        );
        for m in &matchings {
            assert_eq!(m.size(), target);
            assert_matching_valid(&g, m);
        }
        assert_all_distinct(&matchings);
    }

    #[rstest]
    fn test_graph_enumerate_matchings_exhaustive() {
        const EXHAUSTIVE_NODE_BOUND: usize = 5;
        const EXHAUSTIVE_TIME_BUDGET: Duration = Duration::from_secs(30);

        let started = Instant::now();
        for node_count in 0..=EXHAUSTIVE_NODE_BOUND {
            let potential_edges: Vec<_> = (0..node_count as u32)
                .flat_map(|first| (first + 1..node_count as u32).map(move |second| [first, second]))
                .collect();
            for graph_mask in 0_usize..(1 << potential_edges.len()) {
                let edges: Vec<_> = potential_edges
                    .iter()
                    .enumerate()
                    .filter_map(|(index, &edge)| (graph_mask & (1 << index) != 0).then_some(edge))
                    .collect();
                let graph = Graph::new(node_count, &edges);

                let perfect = graph.enumerate_perfect_matchings(BranchAndBound);
                let mut perfect_edges: Vec<_> = perfect
                    .iter()
                    .map(|matching| matching.edges().to_vec())
                    .collect();
                perfect_edges.sort_unstable();
                assert_eq!(perfect_edges, exhaustive_perfect_matchings(&graph));
                for matching in &perfect {
                    assert!(matching.is_perfect(node_count));
                    assert_matching_valid(&graph, matching);
                }
                assert_all_distinct(&perfect);

                let maximum = graph.enumerate_maximum_matchings(BranchAndBound);
                let mut maximum_edges: Vec<_> = maximum
                    .iter()
                    .map(|matching| matching.edges().to_vec())
                    .collect();
                maximum_edges.sort_unstable();
                let expected_maximum = exhaustive_maximum_matchings(&graph);
                assert_eq!(maximum_edges, expected_maximum);
                let maximum_size = expected_maximum.first().map_or(0, Vec::len);
                for matching in &maximum {
                    assert_eq!(matching.size(), maximum_size);
                    assert_matching_valid(&graph, matching);
                }
                assert_all_distinct(&maximum);
            }
        }
        assert!(
            started.elapsed() <= EXHAUSTIVE_TIME_BUDGET,
            "exhaustive enumeration check exceeded {EXHAUSTIVE_TIME_BUDGET:?}",
        );
    }

    #[rstest]
    #[case::mixed(
        6,
        vec![
            [0, 3], [0, 4], [1, 3], [1, 5], [2, 4], [2, 5],
            [0, 1], [4, 5],
        ],
    )]
    #[case::non_bipartite(
        7,
        vec![
            [0, 1], [1, 2], [2, 0], [2, 3], [3, 4], [4, 5], [5, 6],
            [6, 3], [0, 6], [1, 4],
        ],
    )]
    fn test_graph_enumerate_matchings_edge_order(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
    ) {
        let forward = Graph::new(node_count, &edges);
        let reverse = Graph::new(node_count, &edges.iter().rev().copied().collect::<Vec<_>>());
        let canonicalize = |graph: &Graph, matchings: Vec<Matching>| {
            let mut canonical: Vec<_> = matchings
                .into_iter()
                .map(|matching| {
                    let mut selected: Vec<_> = matching
                        .edges()
                        .iter()
                        .map(|&edge| {
                            let [first, second] = graph.edge_endpoints(edge);
                            [first.0.min(second.0), first.0.max(second.0)]
                        })
                        .collect();
                    selected.sort_unstable();
                    selected
                })
                .collect();
            canonical.sort_unstable();
            canonical
        };

        assert_eq!(
            canonicalize(
                &forward,
                forward.enumerate_perfect_matchings(BranchAndBound)
            ),
            canonicalize(
                &reverse,
                reverse.enumerate_perfect_matchings(BranchAndBound)
            )
        );
        assert_eq!(
            canonicalize(
                &forward,
                forward.enumerate_maximum_matchings(BranchAndBound)
            ),
            canonicalize(
                &reverse,
                reverse.enumerate_maximum_matchings(BranchAndBound)
            )
        );
    }

    #[rstest]
    fn test_graph_enumerate_matchings_property() {
        const PROPERTY_CASES: u32 = 32;

        let strategy = (6_usize..=7).prop_flat_map(|node_count| {
            let pairs: Vec<_> = (0..node_count as u32)
                .flat_map(|first| (first + 1..node_count as u32).map(move |second| [first, second]))
                .collect();
            (
                Just(node_count),
                Just(pairs.clone()),
                prop::collection::vec(prop_oneof![3 => Just(false), 1 => Just(true)], pairs.len()),
            )
        });
        let mut runner = TestRunner::new(Config {
            cases: PROPERTY_CASES,
            ..Config::default()
        });

        runner
            .run(&strategy, |(node_count, pairs, present)| {
                let edges: Vec<_> = pairs
                    .into_iter()
                    .zip(present)
                    .filter_map(|(edge, present)| present.then_some(edge))
                    .collect();
                let forward = Graph::new(node_count, &edges);
                let reverse =
                    Graph::new(node_count, &edges.iter().rev().copied().collect::<Vec<_>>());
                let canonicalize = |graph: &Graph, matchings: Vec<Matching>| {
                    let mut canonical: Vec<_> = matchings
                        .into_iter()
                        .map(|matching| {
                            let mut selected: Vec<_> = matching
                                .edges()
                                .iter()
                                .map(|&edge| {
                                    let [first, second] = graph.edge_endpoints(edge);
                                    [first.0.min(second.0), first.0.max(second.0)]
                                })
                                .collect();
                            selected.sort_unstable();
                            selected
                        })
                        .collect();
                    canonical.sort_unstable();
                    canonical
                };

                prop_assert_eq!(
                    canonicalize(
                        &forward,
                        forward.enumerate_perfect_matchings(BranchAndBound)
                    ),
                    canonicalize(
                        &reverse,
                        reverse.enumerate_perfect_matchings(BranchAndBound)
                    )
                );
                prop_assert_eq!(
                    canonicalize(
                        &forward,
                        forward.enumerate_maximum_matchings(BranchAndBound)
                    ),
                    canonicalize(
                        &reverse,
                        reverse.enumerate_maximum_matchings(BranchAndBound)
                    )
                );
                Ok(())
            })
            .unwrap();
    }

    #[rstest]
    #[case::zero(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), vec![])]
    #[case::one(Graph::new(2, &[[0, 1]]), vec![vec![EdgeId(0)]])]
    #[case::full(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![vec![EdgeId(0), EdgeId(2)], vec![EdgeId(1), EdgeId(3)]],
    )]
    fn test_graph_visit_perfect_matchings(
        #[case] graph: Graph,
        #[case] expected: Vec<Vec<EdgeId>>,
    ) {
        let mut visited = Vec::new();
        let result = graph.visit_perfect_matchings(BranchAndBound, |matching| {
            visited.push(matching.edges().to_vec());
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(visited, expected);
    }

    #[rstest]
    #[case::first(1)]
    #[case::prefix(2)]
    fn test_graph_visit_perfect_matchings_break(#[case] stop_after: usize) {
        let graph = Graph::new(6, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [4, 5]]);
        let mut visited = Vec::new();
        let result = graph.visit_perfect_matchings(BranchAndBound, |matching| {
            visited.push(matching.edges().to_vec());
            if visited.len() == stop_after {
                ControlFlow::Break(visited.len())
            } else {
                ControlFlow::Continue(())
            }
        });

        assert_eq!(result, ControlFlow::Break(stop_after));
        assert_eq!(visited.len(), stop_after);
    }

    #[rstest]
    fn test_graph_visit_perfect_matchings_equivalence() {
        let graph = Graph::new(6, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [4, 5]]);
        let expected = graph.enumerate_perfect_matchings(BranchAndBound);
        let mut visited = Vec::new();
        let result = graph.visit_perfect_matchings(BranchAndBound, |matching| {
            visited.push(matching);
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(visited, expected);
    }

    #[rstest]
    fn test_graph_visit_maximum_matchings() {
        let graph = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let mut visited = Vec::new();
        let result = graph.visit_maximum_matchings(BranchAndBound, |matching| {
            visited.push(matching.edges().to_vec());
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(
            visited,
            vec![vec![EdgeId(0)], vec![EdgeId(1)], vec![EdgeId(2)]]
        );
    }

    #[rstest]
    #[case::first(1)]
    #[case::prefix(2)]
    fn test_graph_visit_maximum_matchings_break(#[case] stop_after: usize) {
        let graph = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let mut count = 0;
        let result = graph.visit_maximum_matchings(BranchAndBound, |matching| {
            count += 1;
            assert_eq!(matching.size(), 1);
            if count == stop_after {
                ControlFlow::Break(count)
            } else {
                ControlFlow::Continue(())
            }
        });

        assert_eq!(result, ControlFlow::Break(stop_after));
        assert_eq!(count, stop_after);
    }

    #[rstest]
    fn test_graph_visit_maximum_matchings_equivalence() {
        let graph = Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]);
        let expected = graph.enumerate_maximum_matchings(BranchAndBound);
        let mut visited = Vec::new();
        let result = graph.visit_maximum_matchings(BranchAndBound, |matching| {
            visited.push(matching);
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(visited, expected);
    }

    #[rstest]
    fn test_graph_visit_maximum_matchings_retention() {
        let graph = Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]);
        let mut count = 0;
        let mut last_size = 0;
        let result = graph.visit_maximum_matchings(BranchAndBound, |matching| {
            count += 1;
            last_size = matching.size();
            ControlFlow::<()>::Continue(())
        });

        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(count, 3);
        assert_eq!(last_size, 2);
    }

    fn exhaustive_matchings(graph: &Graph) -> Vec<Vec<EdgeId>> {
        let edge_count = graph.edge_count();
        let subset_count = 1_usize
            .checked_shl(edge_count as u32)
            .expect("test oracle graph fits a usize edge mask");
        let mut matchings = Vec::new();

        for mask in 0..subset_count {
            let mut matched = vec![false; graph.node_bound()];
            let mut edges = Vec::new();
            let mut valid = true;
            for edge in graph.edge_ids() {
                if mask & (1 << edge.index()) == 0 {
                    continue;
                }
                let [first, second] = graph.edge_endpoints(edge);
                if matched[first.index()] || matched[second.index()] {
                    valid = false;
                    break;
                }
                matched[first.index()] = true;
                matched[second.index()] = true;
                edges.push(edge);
            }
            if valid {
                matchings.push(edges);
            }
        }

        matchings.sort_unstable();
        matchings
    }

    fn exhaustive_exact_matchings(graph: &Graph, size: usize) -> Vec<Vec<EdgeId>> {
        exhaustive_matchings(graph)
            .into_iter()
            .filter(|matching| matching.len() == size)
            .collect()
    }

    fn exhaustive_perfect_matchings(graph: &Graph) -> Vec<Vec<EdgeId>> {
        exhaustive_exact_matchings(graph, graph.node_count() / 2)
            .into_iter()
            .filter(|matching| matching.len() * 2 == graph.node_count())
            .collect()
    }

    fn exhaustive_maximum_matchings(graph: &Graph) -> Vec<Vec<EdgeId>> {
        let matchings = exhaustive_matchings(graph);
        let maximum_size = matchings.iter().map(Vec::len).max().unwrap_or(0);
        matchings
            .into_iter()
            .filter(|matching| matching.len() == maximum_size)
            .collect()
    }

    fn greedy_matching_size(graph: &Graph) -> usize {
        let mut matched = vec![false; graph.node_bound()];
        let mut size = 0;
        for edge in graph.edge_ids() {
            let [first, second] = graph.edge_endpoints(edge);
            if !matched[first.index()] && !matched[second.index()] {
                matched[first.index()] = true;
                matched[second.index()] = true;
                size += 1;
            }
        }
        size
    }

    #[rstest]
    #[case::empty(
        Graph::default(),
        0,
        vec![vec![]],
        vec![vec![]],
        vec![vec![]],
        vec![vec![]],
    )]
    #[case::isolated(
        Graph::new(2, &[]),
        0,
        vec![vec![]],
        vec![vec![]],
        vec![],
        vec![vec![]],
    )]
    #[case::path_3(
        Graph::new(3, &[[0, 1], [1, 2]]),
        1,
        vec![vec![], vec![EdgeId(0)], vec![EdgeId(1)]],
        vec![vec![EdgeId(0)], vec![EdgeId(1)]],
        vec![],
        vec![vec![EdgeId(0)], vec![EdgeId(1)]],
    )]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        1,
        vec![vec![], vec![EdgeId(0)], vec![EdgeId(1)], vec![EdgeId(2)]],
        vec![vec![EdgeId(0)], vec![EdgeId(1)], vec![EdgeId(2)]],
        vec![],
        vec![vec![EdgeId(0)], vec![EdgeId(1)], vec![EdgeId(2)]],
    )]
    #[case::square(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        2,
        vec![
            vec![],
            vec![EdgeId(0)],
            vec![EdgeId(0), EdgeId(2)],
            vec![EdgeId(1)],
            vec![EdgeId(1), EdgeId(3)],
            vec![EdgeId(2)],
            vec![EdgeId(3)],
        ],
        vec![vec![EdgeId(0), EdgeId(2)], vec![EdgeId(1), EdgeId(3)]],
        vec![vec![EdgeId(0), EdgeId(2)], vec![EdgeId(1), EdgeId(3)]],
        vec![vec![EdgeId(0), EdgeId(2)], vec![EdgeId(1), EdgeId(3)]],
    )]
    #[case::disconnected_edges(
        Graph::new(4, &[[0, 1], [2, 3]]),
        1,
        vec![
            vec![],
            vec![EdgeId(0)],
            vec![EdgeId(0), EdgeId(1)],
            vec![EdgeId(1)],
        ],
        vec![vec![EdgeId(0)], vec![EdgeId(1)]],
        vec![vec![EdgeId(0), EdgeId(1)]],
        vec![vec![EdgeId(0), EdgeId(1)]],
    )]
    fn test_exhaustive_matching_sets(
        #[case] graph: Graph,
        #[case] exact_size: usize,
        #[case] expected_all: Vec<Vec<EdgeId>>,
        #[case] expected_exact: Vec<Vec<EdgeId>>,
        #[case] expected_perfect: Vec<Vec<EdgeId>>,
        #[case] expected_maximum: Vec<Vec<EdgeId>>,
    ) {
        assert_eq!(exhaustive_matchings(&graph), expected_all);
        assert_eq!(
            exhaustive_exact_matchings(&graph, exact_size),
            expected_exact
        );
        assert_eq!(exhaustive_perfect_matchings(&graph), expected_perfect);
        assert_eq!(exhaustive_maximum_matchings(&graph), expected_maximum);
    }

    #[rstest]
    #[case::middle_edge_first(Graph::new(4, &[[1, 2], [0, 1], [2, 3]]), 1, 2)]
    fn test_greedy_matching_cardinality_bound(
        #[case] graph: Graph,
        #[case] expected_greedy: usize,
        #[case] expected_maximum: usize,
    ) {
        assert_eq!(greedy_matching_size(&graph), expected_greedy);
        assert_eq!(
            exhaustive_maximum_matchings(&graph)[0].len(),
            expected_maximum
        );
    }

    #[rstest]
    fn test_matching_search_state_include() {
        let graph = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [0, 3]]);
        let mut state = MatchingSearchState::new(&graph);
        let initial = state.clone();

        let undo = state.include(EdgeId(1));
        assert_eq!(state.included, vec![false, true, false, false]);
        assert_eq!(state.excluded, vec![true, false, true, false]);
        assert_eq!(state.covered, vec![false, true, true, false]);
        assert_eq!(state.included_size, 1);

        state.undo_include(undo);
        assert_eq!(state, initial);
    }

    #[rstest]
    fn test_matching_search_state_exclude() {
        let graph = Graph::new(3, &[[0, 1], [1, 2]]);
        let mut state = MatchingSearchState::new(&graph);
        let initial = state.clone();

        let undo = state.exclude(EdgeId(1));
        assert_eq!(state.included, vec![false, false]);
        assert_eq!(state.excluded, vec![false, true]);
        assert_eq!(state.covered, vec![false, false, false]);
        assert_eq!(state.included_size, 0);

        state.undo_exclude(undo);
        assert_eq!(state, initial);
    }

    #[rstest]
    fn test_matching_search_state_residual_graph() {
        let graph = Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [2, 5]]);
        let mut state = MatchingSearchState::new(&graph);
        state.include(EdgeId(0));
        state.exclude(EdgeId(3));

        let (residual, correspondence) = state.residual_graph();

        assert_eq!(residual, Graph::new(4, &[[0, 1], [2, 3], [0, 3]]));
        assert_eq!(
            correspondence.nodes().mates(),
            &[
                (NodeId(0), NodeId(2)),
                (NodeId(1), NodeId(3)),
                (NodeId(2), NodeId(4)),
                (NodeId(3), NodeId(5)),
            ]
        );
        assert_eq!(
            correspondence.edges().mates(),
            &[
                (EdgeId(0), EdgeId(2)),
                (EdgeId(1), EdgeId(4)),
                (EdgeId(2), EdgeId(5)),
            ]
        );
        assert_eq!(
            correspondence.nodes().right_exposed(),
            vec![NodeId(0), NodeId(1)]
        );
        assert_eq!(
            correspondence.edges().right_exposed(),
            vec![EdgeId(0), EdgeId(1), EdgeId(3)]
        );
    }

    #[rstest]
    fn test_matching_search_state_residual_graph_exhaustive() {
        let potential_edges = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];

        for graph_mask in 0_usize..(1 << potential_edges.len()) {
            let edges: Vec<_> = potential_edges
                .iter()
                .enumerate()
                .filter_map(|(index, &edge)| (graph_mask & (1 << index) != 0).then_some(edge))
                .collect();
            let graph = Graph::new(4, &edges);
            let state_count = 3_usize.pow(graph.edge_count() as u32);

            for encoded_state in 0..state_count {
                let mut state = MatchingSearchState::new(&graph);
                let mut encoding = encoded_state;
                let mut valid = true;
                for edge in graph.edge_ids() {
                    match encoding % 3 {
                        1 => {
                            let [first, second] = graph.edge_endpoints(edge);
                            if state.covered[first.index()] || state.covered[second.index()] {
                                valid = false;
                                break;
                            }
                            state.included[edge.index()] = true;
                            state.covered[first.index()] = true;
                            state.covered[second.index()] = true;
                            state.included_size += 1;
                        }
                        2 => state.excluded[edge.index()] = true,
                        _ => {}
                    }
                    encoding /= 3;
                }
                if !valid {
                    continue;
                }
                for edge in graph.edge_ids() {
                    if state.included[edge.index()] {
                        continue;
                    }
                    let [first, second] = graph.edge_endpoints(edge);
                    if (state.covered[first.index()] || state.covered[second.index()])
                        && !state.excluded[edge.index()]
                    {
                        valid = false;
                        break;
                    }
                }
                if !valid {
                    continue;
                }

                let expected_nodes: Vec<_> = graph
                    .node_ids()
                    .filter(|node| !state.covered[node.index()])
                    .collect();
                let mut original_to_residual = vec![None; graph.node_count()];
                for (index, &original) in expected_nodes.iter().enumerate() {
                    original_to_residual[original.index()] = Some(NodeId(index as u32));
                }
                let expected_edges: Vec<_> = graph
                    .edge_ids()
                    .filter(|edge| {
                        let [first, second] = graph.edge_endpoints(*edge);
                        !state.excluded[edge.index()]
                            && !state.covered[first.index()]
                            && !state.covered[second.index()]
                    })
                    .collect();
                let expected_endpoints: Vec<_> = expected_edges
                    .iter()
                    .map(|&edge| {
                        let [first, second] = graph.edge_endpoints(edge);
                        [
                            original_to_residual[first.index()].unwrap().0,
                            original_to_residual[second.index()].unwrap().0,
                        ]
                    })
                    .collect();

                let (residual, correspondence) = state.residual_graph();
                assert_eq!(
                    residual,
                    Graph::new(expected_nodes.len(), &expected_endpoints)
                );
                assert_eq!(
                    correspondence.nodes().mates(),
                    expected_nodes
                        .iter()
                        .enumerate()
                        .map(|(index, &original)| (NodeId(index as u32), original))
                        .collect::<Vec<_>>()
                );
                assert_eq!(
                    correspondence.edges().mates(),
                    expected_edges
                        .iter()
                        .enumerate()
                        .map(|(index, &original)| (EdgeId(index as u32), original))
                        .collect::<Vec<_>>()
                );
            }
        }
    }

    #[rstest]
    #[case::greedy_counterexample(Graph::new(4, &[[1, 2], [0, 1], [2, 3]]), 2, true)]
    #[case::odd_capacity(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), 2, false)]
    #[case::zero_target(Graph::new(2, &[[0, 1]]), 0, true)]
    fn test_matching_search_state_can_extend_to(
        #[case] graph: Graph,
        #[case] target_size: usize,
        #[case] expected: bool,
    ) {
        let state = MatchingSearchState::new(&graph);
        assert_eq!(state.can_extend_to(target_size), expected);
    }

    #[rstest]
    fn test_matching_search_state_can_extend_to_greedy_counterexample() {
        let graph = Graph::new(4, &[[1, 2], [0, 1], [2, 3]]);
        let state = MatchingSearchState::new(&graph);
        let (residual, _) = state.residual_graph();

        assert_eq!(greedy_matching_size(&residual), 1);
        assert_eq!(exhaustive_maximum_matchings(&residual)[0].len(), 2);
        assert!(state.can_extend_to(2));
    }

    #[rstest]
    fn test_matching_search_state_can_extend_to_exhaustive() {
        const EXHAUSTIVE_NODE_BOUND: usize = 5;
        const EXHAUSTIVE_TIME_BUDGET: Duration = Duration::from_secs(60);

        let started = Instant::now();
        for node_count in 0..=EXHAUSTIVE_NODE_BOUND {
            let potential_edges: Vec<_> = (0..node_count as u32)
                .flat_map(|first| (first + 1..node_count as u32).map(move |second| [first, second]))
                .collect();
            let state_count = 4_usize.pow(potential_edges.len() as u32);

            for encoded_state in 0..state_count {
                let mut encoding = encoded_state;
                let mut edges = Vec::new();
                let mut statuses = Vec::new();
                for &edge in &potential_edges {
                    let status = encoding % 4;
                    if status != 0 {
                        edges.push(edge);
                        statuses.push(status);
                    }
                    encoding /= 4;
                }

                let graph = Graph::new(node_count, &edges);
                let mut state = MatchingSearchState::new(&graph);
                let mut valid_state = true;
                for (index, &status) in statuses.iter().enumerate() {
                    let edge = EdgeId(index as u32);
                    match status {
                        2 => {
                            let [first, second] = graph.edge_endpoints(edge);
                            if state.covered[first.index()] || state.covered[second.index()] {
                                valid_state = false;
                                break;
                            }
                            state.included[edge.index()] = true;
                            state.covered[first.index()] = true;
                            state.covered[second.index()] = true;
                            state.included_size += 1;
                        }
                        3 => state.excluded[edge.index()] = true,
                        _ => {}
                    }
                }
                if !valid_state {
                    continue;
                }
                for edge in graph.edge_ids() {
                    if state.included[edge.index()] {
                        continue;
                    }
                    let [first, second] = graph.edge_endpoints(edge);
                    if (state.covered[first.index()] || state.covered[second.index()])
                        && !state.excluded[edge.index()]
                    {
                        valid_state = false;
                        break;
                    }
                }
                if !valid_state {
                    continue;
                }

                let mut reachable_sizes = vec![false; node_count / 2 + 1];
                for subset in 0_usize..(1 << graph.edge_count()) {
                    let mut matched = vec![false; node_count];
                    let mut size = 0;
                    let mut valid_matching = true;
                    for edge in graph.edge_ids() {
                        let selected = subset & (1 << edge.index()) != 0;
                        if state.included[edge.index()] && !selected
                            || state.excluded[edge.index()] && selected
                        {
                            valid_matching = false;
                            break;
                        }
                        if !selected {
                            continue;
                        }
                        let [first, second] = graph.edge_endpoints(edge);
                        if matched[first.index()] || matched[second.index()] {
                            valid_matching = false;
                            break;
                        }
                        matched[first.index()] = true;
                        matched[second.index()] = true;
                        size += 1;
                    }
                    if valid_matching {
                        reachable_sizes[size] = true;
                    }
                }

                for target_size in 0..=node_count / 2 + 1 {
                    let expected = reachable_sizes.get(target_size).copied().unwrap_or(false);
                    assert_eq!(
                        state.can_extend_to(target_size),
                        expected,
                        "node_count={node_count}, encoded_state={encoded_state}, target={target_size}",
                    );
                }
            }
        }
        assert!(
            started.elapsed() <= EXHAUSTIVE_TIME_BUDGET,
            "exhaustive extension check exceeded {EXHAUSTIVE_TIME_BUDGET:?}",
        );
    }

    fn assert_matching_valid(graph: &Graph, matching: &Matching) {
        let mut matched = vec![false; graph.node_count()];
        for &eid in matching.edges() {
            let [a, b] = graph.edge_endpoints(eid);
            assert!(!matched[a.index()], "node {} matched twice", a.0);
            assert!(!matched[b.index()], "node {} matched twice", b.0);
            matched[a.index()] = true;
            matched[b.index()] = true;
        }
    }

    fn assert_all_distinct(matchings: &[Matching]) {
        for i in 0..matchings.len() {
            for j in (i + 1)..matchings.len() {
                assert_ne!(
                    matchings[i].edges(),
                    matchings[j].edges(),
                    "matchings {} and {} are identical",
                    i,
                    j
                );
            }
        }
    }
}
