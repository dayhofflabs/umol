//! Bounded enumeration of connected subgraphs as edge sets.
//!
//! [`Graph::enumerate_connected_subgraphs`] currently uses an edge-oriented
//! adaptation of ESU to emit each connected edge set once. See
//! [Wernicke, *A Faster Algorithm for Detecting Network Motifs*
//! (2005)](https://doi.org/10.1007/11557067_14).

use crate::graph::{EdgeId, Graph};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubgraphEnumerationAlgorithm {
    Esu,
}

impl Graph {
    /// Connected subgraphs (possibly branched) of 1..=`max_size` edges, each once.
    /// Subgraphs are connected edge sets — a subset of bonds, not the node-induced
    /// subgraph, so chords are excluded unless their bond is in the set.
    pub fn enumerate_connected_subgraphs(
        &self,
        max_size: u32,
        algorithm: SubgraphEnumerationAlgorithm,
    ) -> Vec<Vec<EdgeId>> {
        match algorithm {
            SubgraphEnumerationAlgorithm::Esu => self.enumerate_connected_subgraphs_esu(max_size),
        }
    }

    fn enumerate_connected_subgraphs_esu(&self, max_size: u32) -> Vec<Vec<EdgeId>> {
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
        max_size: u32,
        result: &mut Vec<Vec<EdgeId>>,
    ) {
        if subgraph.len() == max_size as usize {
            return;
        }
        while let Some(candidate) = extension.pop() {
            let mut next_extension = extension.clone();
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
        #[case] max_size: u32,
        #[case] expected: Vec<Vec<u32>>,
    ) {
        let subgraphs =
            graph.enumerate_connected_subgraphs(max_size, SubgraphEnumerationAlgorithm::Esu);
        let mut expected = expected;
        expected.sort_unstable();
        assert_eq!(sorted_sets(subgraphs), expected);
    }
}
