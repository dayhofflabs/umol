//! Bounded enumeration of simple paths and connected subgraphs, as edge sets.
//!
//! Both return `Vec<Vec<EdgeId>>` — a feature is an edge set, since a ring and its
//! spanning path differ only in edges. Sizes are measured in edges (bonds).

use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathEnumerationAlgorithm {
    Dfs,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubgraphEnumerationAlgorithm {
    Esu,
}

impl Graph {
    /// Simple paths (no repeated node) of 1..=`max_length` edges, each undirected
    /// path once. Each path is the ordered edge sequence from one endpoint.
    pub fn enumerate_paths(
        &self,
        max_length: usize,
        algorithm: PathEnumerationAlgorithm,
    ) -> Vec<Vec<EdgeId>> {
        match algorithm {
            PathEnumerationAlgorithm::Dfs => self.enumerate_paths_dfs(max_length),
        }
    }

    /// Connected subgraphs (possibly branched) of 1..=`max_size` edges, each once.
    /// Subgraphs are connected edge sets — a subset of bonds, not the node-induced
    /// subgraph, so chords are excluded unless their bond is in the set.
    pub fn enumerate_connected_subgraphs(
        &self,
        max_size: usize,
        algorithm: SubgraphEnumerationAlgorithm,
    ) -> Vec<Vec<EdgeId>> {
        match algorithm {
            SubgraphEnumerationAlgorithm::Esu => self.enumerate_connected_subgraphs_esu(max_size),
        }
    }

    fn enumerate_paths_dfs(&self, max_length: usize) -> Vec<Vec<EdgeId>> {
        let mut result = Vec::new();
        if max_length == 0 {
            return result;
        }
        let mut visited = vec![false; self.node_count()];
        let mut edges: Vec<EdgeId> = Vec::with_capacity(max_length);
        for start in self.node_ids() {
            visited[start.index()] = true;
            self.path_dfs(start, start, max_length, &mut visited, &mut edges, &mut result);
            visited[start.index()] = false;
        }
        result
    }

    fn path_dfs(
        &self,
        start: NodeId,
        current: NodeId,
        max_length: usize,
        visited: &mut [bool],
        edges: &mut Vec<EdgeId>,
        result: &mut Vec<Vec<EdgeId>>,
    ) {
        if edges.len() == max_length {
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

    fn enumerate_connected_subgraphs_esu(&self, max_size: usize) -> Vec<Vec<EdgeId>> {
        let mut result = Vec::new();
        if max_size == 0 {
            return result;
        }
        let edge_adjacency = self.edge_adjacency();
        for root in self.edge_ids() {
            let mut subgraph = vec![root];
            let mut extension: Vec<EdgeId> = edge_adjacency[root.index()]
                .iter()
                .copied()
                .filter(|candidate| candidate.0 > root.0)
                .collect();
            result.push(vec![root]);
            self.esu_extend(
                &edge_adjacency,
                root,
                &mut subgraph,
                &mut extension,
                max_size,
                &mut result,
            );
        }
        result
    }

    fn esu_extend(
        &self,
        edge_adjacency: &[Vec<EdgeId>],
        root: EdgeId,
        subgraph: &mut Vec<EdgeId>,
        extension: &mut Vec<EdgeId>,
        max_size: usize,
        result: &mut Vec<Vec<EdgeId>>,
    ) {
        if subgraph.len() == max_size {
            return;
        }
        while let Some(candidate) = extension.pop() {
            let mut next_extension = extension.clone();
            // Exclusive neighborhood: edges adjacent to `candidate`, above `root`, not
            // already in or adjacent to the subgraph.
            for &edge in &edge_adjacency[candidate.index()] {
                if edge.0 > root.0
                    && !subgraph.contains(&edge)
                    && !edge_adjacent_to_set(edge, subgraph, edge_adjacency)
                    && !next_extension.contains(&edge)
                {
                    next_extension.push(edge);
                }
            }
            subgraph.push(candidate);
            result.push(subgraph.clone());
            self.esu_extend(
                edge_adjacency,
                root,
                subgraph,
                &mut next_extension,
                max_size,
                result,
            );
            subgraph.pop();
        }
    }
}

/// Whether `edge` shares an endpoint with any edge in `set`.
fn edge_adjacent_to_set(edge: EdgeId, set: &[EdgeId], edge_adjacency: &[Vec<EdgeId>]) -> bool {
    edge_adjacency[edge.index()]
        .iter()
        .any(|adjacent| set.contains(adjacent))
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
        #[case] max_length: usize,
        #[case] expected: Vec<Vec<u32>>,
    ) {
        let paths = graph.enumerate_paths(max_length, PathEnumerationAlgorithm::Dfs);
        let mut expected = expected;
        expected.sort_unstable();
        assert_eq!(sorted_sets(paths), expected);
    }

    #[rstest]
    #[case::path3_len3(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), 3,
        vec![vec![0], vec![1], vec![2], vec![0, 1], vec![1, 2], vec![0, 1, 2]]
    )]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), 3,
        vec![vec![0], vec![1], vec![2], vec![0, 1], vec![0, 2], vec![1, 2], vec![0, 1, 2]]
    )]
    #[case::star(
        Graph::new(4, &[[0, 1], [0, 2], [0, 3]]), 3,
        vec![
            vec![0], vec![1], vec![2], vec![0, 1], vec![0, 2], vec![1, 2], vec![0, 1, 2]
        ]
    )]
    fn test_graph_enumerate_connected_subgraphs(
        #[case] graph: Graph,
        #[case] max_size: usize,
        #[case] expected: Vec<Vec<u32>>,
    ) {
        let subgraphs =
            graph.enumerate_connected_subgraphs(max_size, SubgraphEnumerationAlgorithm::Esu);
        let mut expected = expected;
        expected.sort_unstable();
        assert_eq!(sorted_sets(subgraphs), expected);
    }
}
