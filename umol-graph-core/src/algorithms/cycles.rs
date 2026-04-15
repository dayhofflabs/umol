use std::collections::HashSet;

use crate::graph::{Graph, NodeId};

impl<N, E> Graph<N, E> {
    pub fn enumerate_simple_cycles(&self, max_cycle_size: usize) -> Vec<Vec<NodeId>> {
        if max_cycle_size < 3 || self.node_count() < 3 {
            return Vec::new();
        }

        let mut raw_cycles: Vec<Vec<NodeId>> = Vec::new();
        for start in self.node_ids() {
            let mut path = vec![start];
            let mut visited = HashSet::new();
            visited.insert(start);
            self.dfs_cycles(
                start,
                start,
                &mut path,
                &mut visited,
                max_cycle_size,
                &mut raw_cycles,
            );
        }

        let mut seen: HashSet<Vec<NodeId>> = HashSet::new();
        let mut result = Vec::new();
        for cycle in raw_cycles {
            let normalized = normalize_cycle(&cycle);
            if seen.insert(normalized.clone()) {
                result.push(normalized);
            }
        }

        result.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        result
    }

    fn dfs_cycles(
        &self,
        start: NodeId,
        current: NodeId,
        path: &mut Vec<NodeId>,
        visited: &mut HashSet<NodeId>,
        max_cycle_size: usize,
        cycles: &mut Vec<Vec<NodeId>>,
    ) {
        if path.len() > max_cycle_size {
            return;
        }

        for neighbor in self.neighbors(current) {
            let next = neighbor.node;
            if next == start && path.len() >= 3 {
                cycles.push(path.clone());
            } else if next > start && !visited.contains(&next) && path.len() < max_cycle_size {
                visited.insert(next);
                path.push(next);
                self.dfs_cycles(start, next, path, visited, max_cycle_size, cycles);
                path.pop();
                visited.remove(&next);
            }
        }
    }
}

fn normalize_cycle(cycle: &[NodeId]) -> Vec<NodeId> {
    let n = cycle.len();
    debug_assert!(n >= 3);

    let min_pos = cycle
        .iter()
        .enumerate()
        .min_by_key(|&(_, id)| id)
        .expect("non-empty cycle")
        .0;

    let mut rotated = Vec::with_capacity(n);
    for i in 0..n {
        rotated.push(cycle[(min_pos + i) % n]);
    }

    if n > 1 && rotated[1] > rotated[n - 1] {
        rotated[1..].reverse();
    }

    rotated
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use crate::graph::{Graph, NodeId};

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    #[rstest]
    #[case::hexagon(
        6,
        vec![(0, 1, ()), (1, 2, ()), (2, 3, ()), (3, 4, ()), (4, 5, ()), (5, 0, ())],
        6,
        vec![vec![n(0), n(1), n(2), n(3), n(4), n(5)]]
    )]
    #[case::two_fused_triangles(
        5,
        vec![(0, 1, ()), (1, 2, ()), (0, 2, ()), (2, 3, ()), (3, 4, ()), (2, 4, ())],
        5,
        vec![vec![n(0), n(1), n(2)], vec![n(2), n(3), n(4)]]
    )]
    fn test_graph_enumerate_simple_cycles(
        #[case] node_count: usize,
        #[case] edges: Vec<(u32, u32, ())>,
        #[case] max_size: usize,
        #[case] expected: Vec<Vec<NodeId>>,
    ) {
        let g = Graph::<(), _>::from_edges(node_count, edges);
        assert_eq!(g.enumerate_simple_cycles(max_size), expected);
    }

    #[test]
    fn test_graph_enumerate_simple_cycles_naphthalene() {
        let g = Graph::<(), ()>::from_edges(
            10,
            vec![
                (0, 1, ()),
                (1, 2, ()),
                (2, 3, ()),
                (3, 4, ()),
                (4, 5, ()),
                (5, 0, ()),
                (3, 6, ()),
                (6, 7, ()),
                (7, 8, ()),
                (8, 9, ()),
                (9, 4, ()),
            ],
        );
        let cycles = g.enumerate_simple_cycles(10);
        assert_eq!(cycles.len(), 3);
    }

    #[test]
    fn test_graph_enumerate_simple_cycles_max_size_cutoff() {
        let g = Graph::<(), ()>::from_edges(
            5,
            vec![(0, 1, ()), (1, 2, ()), (2, 3, ()), (3, 4, ()), (4, 0, ())],
        );
        assert!(g.enumerate_simple_cycles(4).is_empty());
        assert_eq!(
            g.enumerate_simple_cycles(5),
            vec![vec![n(0), n(1), n(2), n(3), n(4)]]
        );
    }

    #[test]
    fn test_graph_enumerate_simple_cycles_empty() {
        let g = Graph::<(), ()>::from_edges(3, vec![(0, 1, ()), (1, 2, ())]);
        assert!(g.enumerate_simple_cycles(10).is_empty());
    }
}
