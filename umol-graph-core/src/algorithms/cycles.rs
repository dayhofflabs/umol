//! Relevant cycle enumeration and shortest-cycle queries.

use std::collections::{HashSet, VecDeque};

use crate::algorithms::bcc::BiconnectedComponentsAlgorithm;
use crate::graph::{EdgeId, Graph, NodeId};

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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ShortestCycleAlgorithm {
    Bfs,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CycleEnumerationAlgorithm {
    Vismara,
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

    pub fn enumerate_cycles(
        &self,
        max_cycle_size: usize,
        alg: CycleEnumerationAlgorithm,
    ) -> Vec<Vec<NodeId>> {
        match alg {
            CycleEnumerationAlgorithm::Vismara => self.enumerate_cycles_vismara(max_cycle_size),
        }
    }

    // BFS with edge exclusion. O(V+E).
    fn shortest_cycle_through_edge_bfs(&self, edge: EdgeId) -> Option<usize> {
        let [u, v] = self.edge_endpoints(edge);
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

    // Vismara 1997. Ref impl: CDK InitialCycles.java. O(V * (V+E)) per BCC.
    fn enumerate_cycles_vismara(&self, max_cycle_size: usize) -> Vec<Vec<NodeId>> {
        if max_cycle_size < 3 || self.node_count() < 3 {
            return Vec::new();
        }

        let mut seen: HashSet<Vec<NodeId>> = HashSet::new();
        let mut result = Vec::new();

        // Vismara enumerates relevant cycles independently within each
        // biconnected component. Tarjan supplies the O(V+E) decomposition
        // fixed by this implementation.
        for component in self.biconnected_components(BiconnectedComponentsAlgorithm::Tarjan) {
            let subgraph = self.induced_subgraph(&component);
            let sub_graph = self.extract(&subgraph);
            let cycles = relevant_cycles_in_bcc(&sub_graph, max_cycle_size);
            for cycle in cycles {
                let mapped: Vec<NodeId> = cycle
                    .iter()
                    .map(|&sub| {
                        subgraph
                            .nodes()
                            .right_of(sub)
                            .expect("subgraph node maps to a host node")
                    })
                    .collect();
                let normalized = normalize_cycle(self, &mapped);
                if seen.insert(normalized.clone()) {
                    result.push(normalized);
                }
            }
        }

        result.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        result
    }
}

struct ShortestPathTree {
    dist: Vec<u32>,
    parents: Vec<Vec<NodeId>>,
}

impl ShortestPathTree {
    fn bfs(graph: &Graph, root: NodeId) -> Self {
        let n = graph.node_bound();
        let mut dist = vec![u32::MAX; n];
        let mut parents: Vec<Vec<NodeId>> = vec![Vec::new(); n];
        let mut queue = VecDeque::new();

        dist[root.index()] = 0;
        queue.push_back(root);

        while let Some(current) = queue.pop_front() {
            let d = dist[current.index()];
            for nbr in graph.neighbors(current) {
                let nd = d + 1;
                if nd < dist[nbr.node.index()] {
                    dist[nbr.node.index()] = nd;
                    parents[nbr.node.index()] = vec![current];
                    queue.push_back(nbr.node);
                } else if nd == dist[nbr.node.index()] {
                    parents[nbr.node.index()].push(current);
                }
            }
        }

        Self { dist, parents }
    }

    fn all_paths_to(&self, root: NodeId, target: NodeId) -> Vec<Vec<NodeId>> {
        let mut result = Vec::new();
        let mut path = vec![target];
        self.reconstruct(root, target, &mut path, &mut result);
        result
    }

    fn reconstruct(
        &self,
        root: NodeId,
        current: NodeId,
        path: &mut Vec<NodeId>,
        result: &mut Vec<Vec<NodeId>>,
    ) {
        if current == root {
            let mut full = path.clone();
            full.reverse();
            result.push(full);
            return;
        }
        for &p in &self.parents[current.index()] {
            path.push(p);
            self.reconstruct(root, p, path, result);
            path.pop();
        }
    }
}

// Vismara 1997 "Unions of all the minimum cycle bases of a graph".
// Ref impl: CDK InitialCycles.java (Algorithm 1).
//
// For each vertex r in the BCC, build a shortest-path tree. For each
// non-tree edge (p,q) relative to r, construct prototype cycles:
//   - Odd prototype (dist[r][p] + dist[r][q] + 1 is odd):
//     union of shortest paths r→p and r→q, plus edge (p,q).
//   - Even prototype (dist[r][p] + dist[r][q] + 1 is even):
//     for each vertex z on both shortest paths at the split point,
//     construct cycle through z.
// Relevance filter: cycle of length L is relevant iff L equals the
// shortest cycle through at least one of its edges.
fn relevant_cycles_in_bcc(graph: &Graph, max_cycle_size: usize) -> Vec<Vec<NodeId>> {
    let n = graph.node_bound();
    if n < 3 {
        return Vec::new();
    }

    let nodes: Vec<NodeId> = graph.node_ids().collect();
    let trees: Vec<ShortestPathTree> = nodes
        .iter()
        .map(|&v| ShortestPathTree::bfs(graph, v))
        .collect();

    let shortest_edge: Vec<Option<usize>> = graph
        .edge_ids()
        .map(|eid| graph.shortest_cycle_through_edge_bfs(eid))
        .collect();

    let mut seen: HashSet<Vec<NodeId>> = HashSet::new();
    let mut cycles = Vec::new();

    for (ri, &root) in nodes.iter().enumerate() {
        let tree = &trees[ri];
        for eid in graph.edge_ids() {
            let [p, q] = graph.edge_endpoints(eid);
            let dp = tree.dist[p.index()];
            let dq = tree.dist[q.index()];
            if dp == u32::MAX || dq == u32::MAX {
                continue;
            }

            let cycle_len = dp as usize + dq as usize + 1;
            if cycle_len < 3 || cycle_len > max_cycle_size {
                continue;
            }

            let is_odd = cycle_len % 2 == 1;

            if is_odd {
                // Odd prototype: paths from root to p and q must meet only at root.
                // The edge (p,q) closes the cycle.
                if dp + dq + 1 != cycle_len as u32 {
                    continue;
                }
                let paths_to_p = tree.all_paths_to(root, p);
                let paths_to_q = tree.all_paths_to(root, q);
                for pp in &paths_to_p {
                    for pq in &paths_to_q {
                        if let Some(c) = join_odd_cycle(pp, pq) {
                            if c.len() <= max_cycle_size && is_relevant(graph, &c, &shortest_edge) {
                                let norm = normalize_cycle(graph, &c);
                                if seen.insert(norm.clone()) {
                                    cycles.push(norm);
                                }
                            }
                        }
                    }
                }
            } else {
                // Even prototype: the paths from root to p and root to q must diverge
                // at the same depth. We need vertices z at distance (cycle_len/2 - 1)
                // from root that appear on shortest paths to both p and q.
                let half = cycle_len / 2;
                if dp as usize != half || dq as usize != half {
                    // For even cycles, both endpoints must be at distance half from root
                    // (since cycle_len = dp + dq + 1 and cycle_len is even, this means
                    //  dp + dq is odd, so dp != dq. The edge (p,q) contributes 1.)
                    // Actually for even prototypes: dp + dq + 1 = even means dp + dq is odd.
                    // We need the paths to share a prefix from root and diverge.
                    // Enumerate through the "even" construction: paths r→p and r→q that
                    // share nodes only up to a split point, then diverge.
                    let paths_to_p = tree.all_paths_to(root, p);
                    let paths_to_q = tree.all_paths_to(root, q);
                    for pp in &paths_to_p {
                        for pq in &paths_to_q {
                            if let Some(c) = join_even_cycle(pp, pq) {
                                if c.len() <= max_cycle_size
                                    && is_relevant(graph, &c, &shortest_edge)
                                {
                                    let norm = normalize_cycle(graph, &c);
                                    if seen.insert(norm.clone()) {
                                        cycles.push(norm);
                                    }
                                }
                            }
                        }
                    }
                    continue;
                }

                let paths_to_p = tree.all_paths_to(root, p);
                let paths_to_q = tree.all_paths_to(root, q);
                for pp in &paths_to_p {
                    for pq in &paths_to_q {
                        if let Some(c) = join_even_cycle(pp, pq) {
                            if c.len() <= max_cycle_size && is_relevant(graph, &c, &shortest_edge) {
                                let norm = normalize_cycle(graph, &c);
                                if seen.insert(norm.clone()) {
                                    cycles.push(norm);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    cycles
}

fn join_odd_cycle(path_p: &[NodeId], path_q: &[NodeId]) -> Option<Vec<NodeId>> {
    if path_p.len() < 2 || path_q.len() < 2 {
        return None;
    }
    debug_assert_eq!(path_p[0], path_q[0]);
    let p_set: HashSet<NodeId> = path_p.iter().copied().collect();
    for &v in &path_q[1..] {
        if p_set.contains(&v) {
            return None;
        }
    }
    // Cycle: path_p forward (root to p), then path_q reversed (q to root), skip root at join
    let mut cycle = path_p.to_vec();
    for &v in path_q.iter().rev() {
        if v != path_p[0] {
            cycle.push(v);
        }
    }
    Some(cycle)
}

fn join_even_cycle(path_p: &[NodeId], path_q: &[NodeId]) -> Option<Vec<NodeId>> {
    if path_p.len() < 2 || path_q.len() < 2 {
        return None;
    }
    debug_assert_eq!(path_p[0], path_q[0]);

    // Find the split point: last shared node on both paths from root
    let shared_len = path_p
        .iter()
        .zip(path_q.iter())
        .take_while(|(a, b)| a == b)
        .count();

    if shared_len == 0 {
        return None;
    }

    // After the split, the two paths must be completely disjoint
    let tail_p: HashSet<NodeId> = path_p[shared_len..].iter().copied().collect();
    for &v in &path_q[shared_len..] {
        if tail_p.contains(&v) {
            return None;
        }
    }

    // Cycle: path from split to p, then back from q to split
    // The split vertex is path_p[shared_len - 1]
    let mut cycle: Vec<NodeId> = path_p[shared_len - 1..].to_vec();
    for &v in path_q[shared_len..].iter().rev() {
        cycle.push(v);
    }

    if cycle.len() < 3 {
        return None;
    }

    Some(cycle)
}

fn is_relevant(graph: &Graph, cycle: &[NodeId], shortest_edge: &[Option<usize>]) -> bool {
    let len = cycle.len();
    for i in 0..len {
        let a = cycle[i];
        let b = cycle[(i + 1) % len];
        if let Some(eid) = graph.find_edge(a, b) {
            if let Some(shortest) = shortest_edge[eid.index()] {
                if shortest == len {
                    return true;
                }
            }
        }
    }
    false
}

fn normalize_cycle(graph: &Graph, nodes: &[NodeId]) -> Vec<NodeId> {
    let edges = nodes
        .iter()
        .copied()
        .zip(nodes.iter().copied().cycle().skip(1))
        .take(nodes.len())
        .map(|(first, second)| {
            graph
                .find_edge(first, second)
                .expect("consecutive cycle nodes must share an edge")
        })
        .collect();
    Cycle::normalized(graph, nodes.to_vec(), edges).nodes
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::Cycle;
    use super::CycleEnumerationAlgorithm::Vismara;
    use super::ShortestCycleAlgorithm::Bfs;
    use crate::graph::{EdgeId, Graph, NodeId};

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

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    #[rstest]
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
        let g = Graph::new(node_count, &edges);
        assert_eq!(g.shortest_cycle_through_edge(edge, Bfs), expected);
    }

    #[rstest]
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
        let g = Graph::new(node_count, &edges);
        assert_eq!(g.shortest_cycle_through_node(node, Bfs), expected);
    }

    #[rstest]
    #[case::hexagon(6, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]], 6,
        vec![vec![n(0), n(1), n(2), n(3), n(4), n(5)]])]
    #[case::two_fused_triangles(5, vec![[0, 1], [1, 2], [0, 2], [2, 3], [3, 4], [2, 4]], 5,
        vec![vec![n(0), n(1), n(2)], vec![n(2), n(3), n(4)]])]
    #[case::k4(4, vec![[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], 4,
        vec![vec![n(0), n(1), n(2)], vec![n(0), n(1), n(3)], vec![n(0), n(2), n(3)], vec![n(1), n(2), n(3)]])]
    #[case::pentagon_cutoff(5, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]], 4, vec![])]
    #[case::pentagon(5, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]], 5, vec![vec![n(0), n(1), n(2), n(3), n(4)]])]
    #[case::acyclic(3, vec![[0, 1], [1, 2]], 10, vec![])]
    #[case::naphthalene(10,
        vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [3, 6], [6, 7], [7, 8], [8, 9], [9, 4]],
        10,
        vec![vec![n(0), n(1), n(2), n(3), n(4), n(5)],
             vec![n(3), n(4), n(9), n(8), n(7), n(6)]])]
    #[case::cube(8,
        vec![[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]],
        4,
        vec![vec![n(0), n(1), n(2), n(3)], vec![n(0), n(1), n(5), n(4)], vec![n(0), n(3), n(7), n(4)],
             vec![n(1), n(2), n(6), n(5)], vec![n(2), n(3), n(7), n(6)], vec![n(4), n(5), n(6), n(7)]])]
    #[case::prism(6,
        vec![[0, 1], [1, 2], [0, 2], [3, 4], [4, 5], [3, 5], [0, 3], [1, 4], [2, 5]],
        6,
        vec![vec![n(0), n(1), n(2)], vec![n(3), n(4), n(5)],
             vec![n(0), n(1), n(4), n(3)], vec![n(0), n(2), n(5), n(3)], vec![n(1), n(2), n(5), n(4)]])]
    fn test_graph_enumerate_cycles(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] max_size: usize,
        #[case] expected: Vec<Vec<NodeId>>,
    ) {
        let g = Graph::new(node_count, &edges);
        assert_eq!(g.enumerate_cycles(max_size, Vismara), expected);
    }
}
