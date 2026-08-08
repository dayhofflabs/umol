//! Bounded visitation and collection of connected subgraphs as edge sets.
//!
//! [`Graph::visit_connected_subgraphs`] currently uses an edge-oriented
//! adaptation of ESU to emit each connected edge set once; the eager form
//! collects the same traversal. See
//! [Wernicke, *A Faster Algorithm for Detecting Network Motifs*
//! (2005)](https://doi.org/10.1007/11557067_14).

use std::ops::ControlFlow;

use crate::graph::{EdgeId, Graph};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubgraphEnumerationAlgorithm {
    Esu,
}

impl Graph {
    /// Visits each connected subgraph (possibly branched) of 1..=`max_size`
    /// edges once, until traversal completes or the visitor returns
    /// [`ControlFlow::Break`]. Subgraphs are connected edge sets — a subset of
    /// bonds, not the node-induced subgraph, so chords are excluded unless
    /// their bond is in the set. The slice borrows the search's working buffer
    /// and is only valid for the duration of the call. Traversal is
    /// deterministic for a fixed graph representation, but its order is not a
    /// canonical ordering contract.
    pub fn visit_connected_subgraphs<B, F>(
        &self,
        max_size: u32,
        algorithm: SubgraphEnumerationAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(&[EdgeId]) -> ControlFlow<B>,
    {
        match algorithm {
            SubgraphEnumerationAlgorithm::Esu => {
                self.visit_connected_subgraphs_esu(max_size, &mut visitor)
            }
        }
    }

    /// Collects every connected subgraph as an owned edge set by collecting
    /// [`Graph::visit_connected_subgraphs`].
    pub fn enumerate_connected_subgraphs(
        &self,
        max_size: u32,
        algorithm: SubgraphEnumerationAlgorithm,
    ) -> Vec<Vec<EdgeId>> {
        let mut result = Vec::new();
        let _: ControlFlow<()> = self.visit_connected_subgraphs(max_size, algorithm, |subgraph| {
            result.push(subgraph.to_vec());
            ControlFlow::Continue(())
        });
        result
    }

    fn visit_connected_subgraphs_esu<B, F>(&self, max_size: u32, visitor: &mut F) -> ControlFlow<B>
    where
        F: FnMut(&[EdgeId]) -> ControlFlow<B>,
    {
        if max_size == 0 {
            return ControlFlow::Continue(());
        }
        let edge_adjacency = self.edge_adjacency();
        for root in self.edge_ids() {
            let mut subgraph = vec![root];
            let mut extension: Vec<EdgeId> = edge_adjacency[root.index()]
                .iter()
                .copied()
                .filter(|candidate| candidate.0 > root.0)
                .collect();
            if let ControlFlow::Break(value) = visitor(&subgraph) {
                return ControlFlow::Break(value);
            }
            if let ControlFlow::Break(value) = self.esu_extend(
                &edge_adjacency,
                root,
                &mut subgraph,
                &mut extension,
                max_size,
                visitor,
            ) {
                return ControlFlow::Break(value);
            }
        }
        ControlFlow::Continue(())
    }

    fn esu_extend<B, F>(
        &self,
        edge_adjacency: &[Vec<EdgeId>],
        root: EdgeId,
        subgraph: &mut Vec<EdgeId>,
        extension: &mut Vec<EdgeId>,
        max_size: u32,
        visitor: &mut F,
    ) -> ControlFlow<B>
    where
        F: FnMut(&[EdgeId]) -> ControlFlow<B>,
    {
        if subgraph.len() == max_size as usize {
            return ControlFlow::Continue(());
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
            let result = match visitor(subgraph) {
                ControlFlow::Continue(()) => self.esu_extend(
                    edge_adjacency,
                    root,
                    subgraph,
                    &mut next_extension,
                    max_size,
                    visitor,
                ),
                ControlFlow::Break(value) => ControlFlow::Break(value),
            };
            subgraph.pop();
            if let ControlFlow::Break(value) = result {
                return ControlFlow::Break(value);
            }
        }
        ControlFlow::Continue(())
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
    fn test_graph_visit_connected_subgraphs(
        #[case] graph: Graph,
        #[case] max_size: u32,
        #[case] expected: Vec<Vec<u32>>,
    ) {
        let mut collected: Vec<Vec<EdgeId>> = Vec::new();
        let flow: ControlFlow<()> =
            graph.visit_connected_subgraphs(max_size, SubgraphEnumerationAlgorithm::Esu, |edges| {
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
        vec![vec![0], vec![1], vec![2], vec![0, 1], vec![0, 2], vec![1, 2], vec![0, 1, 2]]
    )]
    fn test_graph_visit_connected_subgraphs_termination(
        #[case] graph: Graph,
        #[case] max_size: u32,
        #[case] expected: Vec<Vec<u32>>,
    ) {
        let flow =
            graph.visit_connected_subgraphs(max_size, SubgraphEnumerationAlgorithm::Esu, |edges| {
                ControlFlow::Break(edges.to_vec())
            });
        let ControlFlow::Break(first) = flow else {
            panic!("expected Break on first emission");
        };
        let mut ids: Vec<u32> = first.iter().map(|e| e.0).collect();
        ids.sort_unstable();
        assert!(expected.contains(&ids), "invalid subgraph {ids:?}");
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
