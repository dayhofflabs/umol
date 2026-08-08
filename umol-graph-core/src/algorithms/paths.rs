//! Bounded visitation and collection of simple paths as edge sequences.
//!
//! [`Graph::visit_paths`] currently uses depth-first search with a
//! maximum-length bound and emits each undirected simple path once; the eager
//! form collects the same traversal.

use std::ops::ControlFlow;

use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathEnumerationAlgorithm {
    Dfs,
}

impl Graph {
    /// Visits each simple path (no repeated node) of 1..=`max_length` edges,
    /// each undirected path once as the ordered edge sequence from one
    /// endpoint, until traversal completes or the visitor returns
    /// [`ControlFlow::Break`]. The slice borrows the search's working buffer
    /// and is only valid for the duration of the call. Traversal is
    /// deterministic for a fixed graph representation, but its order is not a
    /// canonical ordering contract.
    pub fn visit_paths<B, F>(
        &self,
        max_length: u32,
        algorithm: PathEnumerationAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(&[EdgeId]) -> ControlFlow<B>,
    {
        match algorithm {
            PathEnumerationAlgorithm::Dfs => self.visit_paths_dfs(max_length, &mut visitor),
        }
    }

    /// Collects every simple path as an owned edge sequence by collecting
    /// [`Graph::visit_paths`].
    pub fn enumerate_paths(
        &self,
        max_length: u32,
        algorithm: PathEnumerationAlgorithm,
    ) -> Vec<Vec<EdgeId>> {
        let mut result = Vec::new();
        let _: ControlFlow<()> = self.visit_paths(max_length, algorithm, |edges| {
            result.push(edges.to_vec());
            ControlFlow::Continue(())
        });
        result
    }

    fn visit_paths_dfs<B, F>(&self, max_length: u32, visitor: &mut F) -> ControlFlow<B>
    where
        F: FnMut(&[EdgeId]) -> ControlFlow<B>,
    {
        if max_length == 0 {
            return ControlFlow::Continue(());
        }
        let mut visited = vec![false; self.node_count()];
        let mut edges: Vec<EdgeId> = Vec::with_capacity(max_length as usize);
        for start in self.node_ids() {
            visited[start.index()] = true;
            let result = self.path_dfs(start, start, max_length, &mut visited, &mut edges, visitor);
            visited[start.index()] = false;
            if let ControlFlow::Break(value) = result {
                return ControlFlow::Break(value);
            }
        }
        ControlFlow::Continue(())
    }

    fn path_dfs<B, F>(
        &self,
        start: NodeId,
        current: NodeId,
        max_length: u32,
        visited: &mut [bool],
        edges: &mut Vec<EdgeId>,
        visitor: &mut F,
    ) -> ControlFlow<B>
    where
        F: FnMut(&[EdgeId]) -> ControlFlow<B>,
    {
        if edges.len() == max_length as usize {
            return ControlFlow::Continue(());
        }
        for neighbor in self.neighbors(current) {
            if visited[neighbor.node.index()] {
                continue;
            }
            visited[neighbor.node.index()] = true;
            edges.push(neighbor.edge);
            // Canonical orientation: emit each undirected path from its lower endpoint.
            let emitted = if start.0 < neighbor.node.0 {
                visitor(edges)
            } else {
                ControlFlow::Continue(())
            };
            let result = match emitted {
                ControlFlow::Continue(()) => {
                    self.path_dfs(start, neighbor.node, max_length, visited, edges, visitor)
                }
                ControlFlow::Break(value) => ControlFlow::Break(value),
            };
            edges.pop();
            visited[neighbor.node.index()] = false;
            if let ControlFlow::Break(value) = result {
                return ControlFlow::Break(value);
            }
        }
        ControlFlow::Continue(())
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    fn sorted_sets(mut sets: Vec<Vec<EdgeId>>) -> Vec<Vec<u32>> {
        let mut out: Vec<Vec<u32>> = sets
            .iter_mut()
            .map(|s| {
                let mut ids: Vec<u32> = s.iter().map(|e| e.0).collect();
                ids.sort_unstable();
                ids
            })
            .collect();
        out.sort_unstable();
        out
    }

    #[rstest]
    #[case::path3_len1(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), 1,
        vec![vec![0], vec![1], vec![2]]
    )]
    #[case::path3_len3(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), 3,
        vec![vec![0], vec![1], vec![2], vec![0, 1], vec![1, 2], vec![0, 1, 2]]
    )]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), 3,
        vec![vec![0], vec![1], vec![2], vec![0, 1], vec![0, 2], vec![1, 2]]
    )]
    #[case::star(
        Graph::new(4, &[[0, 1], [0, 2], [0, 3]]), 2,
        vec![vec![0], vec![1], vec![2], vec![0, 1], vec![0, 2], vec![1, 2]]
    )]
    fn test_graph_visit_paths(
        #[case] graph: Graph,
        #[case] max_length: u32,
        #[case] expected: Vec<Vec<u32>>,
    ) {
        let mut collected: Vec<Vec<EdgeId>> = Vec::new();
        let flow: ControlFlow<()> =
            graph.visit_paths(max_length, PathEnumerationAlgorithm::Dfs, |edges| {
                collected.push(edges.to_vec());
                ControlFlow::Continue(())
            });
        assert_eq!(flow, ControlFlow::Continue(()));
        let mut expected = expected;
        expected.sort_unstable();
        assert_eq!(sorted_sets(collected), expected);
    }

    #[rstest]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), 3,
        vec![vec![0], vec![1], vec![2], vec![0, 1], vec![0, 2], vec![1, 2]]
    )]
    fn test_graph_visit_paths_termination(
        #[case] graph: Graph,
        #[case] max_length: u32,
        #[case] expected: Vec<Vec<u32>>,
    ) {
        let flow = graph.visit_paths(max_length, PathEnumerationAlgorithm::Dfs, |edges| {
            ControlFlow::Break(edges.to_vec())
        });
        let ControlFlow::Break(first) = flow else {
            panic!("expected Break on first emission");
        };
        let mut ids: Vec<u32> = first.iter().map(|e| e.0).collect();
        ids.sort_unstable();
        assert!(expected.contains(&ids), "invalid path {ids:?}");
    }

    #[rstest]
    #[case::path3_len1(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), 1,
        vec![vec![0], vec![1], vec![2]]
    )]
    #[case::path3_len3(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), 3,
        vec![vec![0], vec![1], vec![2], vec![0, 1], vec![1, 2], vec![0, 1, 2]]
    )]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), 3,
        vec![vec![0], vec![1], vec![2], vec![0, 1], vec![0, 2], vec![1, 2]]
    )]
    #[case::star(
        Graph::new(4, &[[0, 1], [0, 2], [0, 3]]), 2,
        vec![vec![0], vec![1], vec![2], vec![0, 1], vec![0, 2], vec![1, 2]]
    )]
    fn test_graph_enumerate_paths(
        #[case] graph: Graph,
        #[case] max_length: u32,
        #[case] expected: Vec<Vec<u32>>,
    ) {
        let paths = graph.enumerate_paths(max_length, PathEnumerationAlgorithm::Dfs);
        let mut expected = expected;
        expected.sort_unstable();
        assert_eq!(sorted_sets(paths), expected);
    }
}
