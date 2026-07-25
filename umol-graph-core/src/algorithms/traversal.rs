//! Distance-limited graph neighborhoods.
//!
//! [`Graph::neighborhood`] is currently the sole traversal operation. It uses
//! breadth-first search to return nodes through a requested distance; this
//! module does not yet expose general traversal events or trees.

use std::collections::VecDeque;

use crate::graph::{Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TraversalAlgorithm {
    Bfs,
}

impl Graph {
    /// Nodes within `max_depth` edges of `source` (inclusive) as `(node, distance)`
    /// pairs, in nondecreasing distance order. `source` appears as `(source, 0)`.
    /// The result is proportional to the neighborhood, not the whole graph.
    pub fn neighborhood(
        &self,
        source: NodeId,
        max_depth: u32,
        algorithm: TraversalAlgorithm,
    ) -> Vec<(NodeId, u32)> {
        match algorithm {
            TraversalAlgorithm::Bfs => self.neighborhood_bfs(source, max_depth),
        }
    }

    fn neighborhood_bfs(&self, source: NodeId, max_depth: u32) -> Vec<(NodeId, u32)> {
        let mut visited = vec![false; self.node_count()];
        visited[source.index()] = true;
        let mut reached = vec![(source, 0)];
        let mut queue = VecDeque::from([(source, 0u32)]);
        while let Some((node, depth)) = queue.pop_front() {
            if depth >= max_depth {
                continue;
            }
            for neighbor in self.neighbors(node) {
                if !visited[neighbor.node.index()] {
                    visited[neighbor.node.index()] = true;
                    reached.push((neighbor.node, depth + 1));
                    queue.push_back((neighbor.node, depth + 1));
                }
            }
        }
        reached
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::path_full(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), 0, 3,
        vec![(NodeId(0), 0), (NodeId(1), 1), (NodeId(2), 2), (NodeId(3), 3)]
    )]
    #[case::path_bounded(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), 0, 1,
        vec![(NodeId(0), 0), (NodeId(1), 1)]
    )]
    #[case::path_from_middle(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), 2, 3,
        vec![(NodeId(0), 2), (NodeId(1), 1), (NodeId(2), 0), (NodeId(3), 1)]
    )]
    #[case::cycle(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]), 0, 3,
        vec![(NodeId(0), 0), (NodeId(1), 1), (NodeId(2), 2), (NodeId(3), 1)]
    )]
    #[case::disconnected(
        Graph::new(3, &[[0, 1]]), 0, 3,
        vec![(NodeId(0), 0), (NodeId(1), 1)]
    )]
    fn test_graph_neighborhood(
        #[case] graph: Graph,
        #[case] source: u32,
        #[case] max_depth: u32,
        #[case] expected: Vec<(NodeId, u32)>,
    ) {
        let mut actual = graph.neighborhood(NodeId(source), max_depth, TraversalAlgorithm::Bfs);
        actual.sort_by_key(|&(node, _)| node.0);
        assert_eq!(actual, expected);
    }
}
