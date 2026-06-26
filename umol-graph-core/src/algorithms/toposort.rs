//! Deterministic topological sort of a `DiGraph`.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::digraph::DiGraph;
use crate::graph::NodeId;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TopologicalSortAlgorithm {
    Kahn,
}

impl DiGraph {
    /// A unique topological order via Kahn's algorithm draining a ready set ordered
    /// by `key` (ties broken by `NodeId`), so the result is a deterministic function
    /// of the graph and `key`. `None` if the graph contains a cycle.
    pub fn topological_order<K: Ord>(
        &self,
        key: impl Fn(NodeId) -> K,
        alg: TopologicalSortAlgorithm,
    ) -> Option<Vec<NodeId>> {
        match alg {
            TopologicalSortAlgorithm::Kahn => self.topological_order_kahn(key),
        }
    }

    fn topological_order_kahn<K: Ord>(&self, key: impl Fn(NodeId) -> K) -> Option<Vec<NodeId>> {
        let mut in_degree = vec![0usize; self.node_count()];
        for node in self.node_ids() {
            for &target in self.successors(node) {
                in_degree[target.index()] += 1;
            }
        }

        let mut ready: BinaryHeap<Reverse<(K, u32)>> = self
            .node_ids()
            .filter(|&node| in_degree[node.index()] == 0)
            .map(|node| Reverse((key(node), node.0)))
            .collect();

        let mut order = Vec::with_capacity(self.node_count());
        while let Some(Reverse((_, raw))) = ready.pop() {
            let node = NodeId(raw);
            order.push(node);
            for &target in self.successors(node) {
                in_degree[target.index()] -= 1;
                if in_degree[target.index()] == 0 {
                    ready.push(Reverse((key(target), target.0)));
                }
            }
        }

        (order.len() == self.node_count()).then_some(order)
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::TopologicalSortAlgorithm::Kahn;
    use super::*;

    #[rstest]
    #[case::linear(3, vec![[0, 1], [1, 2]], Some(vec![NodeId(0), NodeId(1), NodeId(2)]))]
    #[case::diamond(
        4,
        vec![[0, 1], [0, 2], [1, 3], [2, 3]],
        Some(vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)])
    )]
    #[case::independent(3, vec![], Some(vec![NodeId(0), NodeId(1), NodeId(2)]))]
    #[case::cycle(3, vec![[0, 1], [1, 2], [2, 0]], None)]
    #[case::self_loop(1, vec![[0, 0]], None)]
    fn test_digraph_topological_order(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected: Option<Vec<NodeId>>,
    ) {
        let g = DiGraph::new(node_count, &edges);
        assert_eq!(g.topological_order(|node| node.0, Kahn), expected);
    }
}
