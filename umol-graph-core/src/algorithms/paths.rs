//! Bounded enumeration of simple paths as edge sequences.

use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathEnumerationAlgorithm {
    Dfs,
}

impl Graph {
    /// Simple paths (no repeated node) of 1..=`max_length` edges, each undirected
    /// path once. Each path is the ordered edge sequence from one endpoint.
    pub fn enumerate_paths(
        &self,
        max_length: u32,
        algorithm: PathEnumerationAlgorithm,
    ) -> Vec<Vec<EdgeId>> {
        match algorithm {
            PathEnumerationAlgorithm::Dfs => self.enumerate_paths_dfs(max_length),
        }
    }

    fn enumerate_paths_dfs(&self, max_length: u32) -> Vec<Vec<EdgeId>> {
        let mut result = Vec::new();
        if max_length == 0 {
            return result;
        }
        let mut visited = vec![false; self.node_count()];
        let mut edges: Vec<EdgeId> = Vec::with_capacity(max_length as usize);
        for start in self.node_ids() {
            visited[start.index()] = true;
            self.path_dfs(
                start,
                start,
                max_length,
                &mut visited,
                &mut edges,
                &mut result,
            );
            visited[start.index()] = false;
        }
        result
    }

    fn path_dfs(
        &self,
        start: NodeId,
        current: NodeId,
        max_length: u32,
        visited: &mut [bool],
        edges: &mut Vec<EdgeId>,
        result: &mut Vec<Vec<EdgeId>>,
    ) {
        if edges.len() == max_length as usize {
            return;
        }
        for neighbor in self.neighbors(current) {
            if visited[neighbor.node.index()] {
                continue;
            }
            visited[neighbor.node.index()] = true;
            edges.push(neighbor.edge);
            // Canonical orientation: emit each undirected path from its lower endpoint.
            if start.0 < neighbor.node.0 {
                result.push(edges.clone());
            }
            self.path_dfs(start, neighbor.node, max_length, visited, edges, result);
            edges.pop();
            visited[neighbor.node.index()] = false;
        }
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
