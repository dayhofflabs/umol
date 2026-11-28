//! CSR simple directed graph (out-adjacency only).
//!
//! Stores per-node successor lists in compressed sparse row form; node data lives
//! externally, indexed by `NodeId`. Minimal by design — predecessors, edge identity,
//! and copy-on-write are added when a directed algorithm needs them.

use crate::graph::GraphError;
use crate::id::NodeId;

/// Simple directed graph stored as compressed sparse row (CSR).
///
/// Self-loops and duplicate directed edges are rejected; opposite-direction arcs
/// `[a, b]` and `[b, a]` are distinct and permitted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DiGraph {
    offsets: Vec<u32>,
    targets: Vec<NodeId>,
    node_count: usize,
}

impl DiGraph {
    /// Build from directed `[from, to]` edges. Successor order within a node follows
    /// edge input order.
    ///
    /// # Panics
    ///
    /// Panics if an endpoint is outside `0..node_count`, if an edge is a self-loop, or if
    /// duplicate directed edges appear. Use [`Self::try_new`] for untrusted input.
    pub fn new(node_count: usize, edges: &[[u32; 2]]) -> Self {
        Self::try_new(node_count, edges)
            .unwrap_or_else(|error| panic!("invalid directed graph: {error}"))
    }

    /// Checked form of [`Self::new`].
    ///
    /// Successor order within each node follows edge input order. Opposite-direction edges
    /// `[a, b]` and `[b, a]` are distinct and permitted.
    ///
    /// # Errors
    ///
    /// Returns [`GraphError::NodeOutOfBounds`] when an endpoint is `>= node_count`,
    /// [`GraphError::SelfLoop`] when `from == to`, or [`GraphError::ParallelEdge`] when the
    /// same directed pair `[from, to]` appears more than once.
    pub fn try_new(node_count: usize, edges: &[[u32; 2]]) -> Result<Self, GraphError> {
        let mut out_degree = vec![0u32; node_count];
        for &[from, to] in edges {
            if from as usize >= node_count {
                return Err(GraphError::NodeOutOfBounds {
                    node: from,
                    node_count,
                });
            }
            if to as usize >= node_count {
                return Err(GraphError::NodeOutOfBounds {
                    node: to,
                    node_count,
                });
            }
            if from == to {
                return Err(GraphError::SelfLoop { node: from });
            }
            out_degree[from as usize] += 1;
        }
        let mut offsets = vec![0u32; node_count + 1];
        for i in 0..node_count {
            offsets[i + 1] = offsets[i] + out_degree[i];
        }
        let mut targets = vec![NodeId(0); edges.len()];
        let mut cursor: Vec<u32> = offsets[..node_count].to_vec();
        for &[from, to] in edges {
            let pos = cursor[from as usize] as usize;
            targets[pos] = NodeId(to);
            cursor[from as usize] += 1;
        }
        let mut scratch = Vec::new();
        for i in 0..node_count {
            let start = offsets[i] as usize;
            let end = offsets[i + 1] as usize;
            let slice = &targets[start..end];
            if slice.len() <= 1 {
                continue;
            }
            if slice.len() <= 16 {
                for (idx, &target) in slice.iter().enumerate() {
                    if slice[idx + 1..].contains(&target) {
                        return Err(GraphError::ParallelEdge {
                            endpoints: [i as u32, target.0],
                        });
                    }
                }
            } else {
                scratch.clear();
                scratch.extend_from_slice(slice);
                scratch.sort_unstable();
                for pair in scratch.windows(2) {
                    if pair[0] == pair[1] {
                        return Err(GraphError::ParallelEdge {
                            endpoints: [i as u32, pair[0].0],
                        });
                    }
                }
            }
        }
        Ok(Self {
            offsets,
            targets,
            node_count,
        })
    }

    pub fn node_count(&self) -> usize {
        self.node_count
    }

    pub fn node_ids(&self) -> impl ExactSizeIterator<Item = NodeId> {
        (0..self.node_count as u32).map(NodeId)
    }

    /// The out-neighbors of `id`.
    pub fn successors(&self, id: NodeId) -> &[NodeId] {
        let start = self.offsets[id.index()] as usize;
        let end = self.offsets[id.index() + 1] as usize;
        &self.targets[start..end]
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use rstest::*;

    use super::*;

    fn assert_exact_size<T>(mut iterator: impl ExactSizeIterator<Item = T>, expected: Vec<T>)
    where
        T: Debug + PartialEq,
    {
        assert_eq!(iterator.len(), expected.len());
        assert_eq!(iterator.size_hint(), (expected.len(), Some(expected.len())));
        while let Some(expected_item) = expected.get(expected.len() - iterator.len()) {
            let previous = iterator.len();
            assert_eq!(iterator.next().as_ref(), Some(expected_item));
            let remaining = iterator.len();
            assert_eq!(remaining, previous - 1);
            assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
        }
        assert_eq!(iterator.next(), None);
        assert_eq!(iterator.len(), 0);
        assert_eq!(iterator.size_hint(), (0, Some(0)));
    }

    #[rstest]
    #[case::empty(DiGraph::new(0, &[]), vec![])]
    #[case::populated(
        DiGraph::new(3, &[[0, 1], [0, 2], [1, 2]]),
        vec![NodeId(0), NodeId(1), NodeId(2)],
    )]
    fn test_digraph_node_ids(#[case] graph: DiGraph, #[case] expected: Vec<NodeId>) {
        assert_exact_size(graph.node_ids(), expected);
    }

    #[rstest]
    #[case(NodeId(0), vec![NodeId(1), NodeId(2)])]
    #[case(NodeId(1), vec![NodeId(2)])]
    #[case(NodeId(2), vec![])]
    fn test_digraph_successors(#[case] node: NodeId, #[case] expected: Vec<NodeId>) {
        let g = DiGraph::new(3, &[[0, 1], [0, 2], [1, 2]]);
        assert_eq!(g.successors(node), expected);
    }

    #[rstest]
    fn test_digraph_try_new_opposite_arcs() {
        let graph = DiGraph::try_new(2, &[[0, 1], [1, 0]]).unwrap();
        assert_eq!(graph.successors(NodeId(0)), &[NodeId(1)]);
        assert_eq!(graph.successors(NodeId(1)), &[NodeId(0)]);
    }

    #[rstest]
    fn test_digraph_try_new_high_out_degree() {
        let edges: Vec<[u32; 2]> = (1..20).map(|target| [0, target]).collect();
        let graph = DiGraph::try_new(20, &edges).unwrap();
        assert_eq!(
            graph.successors(NodeId(0)),
            &(1..20).map(NodeId).collect::<Vec<_>>()
        );

        let mut duplicate = edges;
        duplicate.push([0, 10]);
        assert_eq!(
            DiGraph::try_new(20, &duplicate),
            Err(GraphError::ParallelEdge { endpoints: [0, 10] })
        );
    }

    #[rstest]
    #[case::from_out_of_bounds(2, &[[2, 0]], GraphError::NodeOutOfBounds { node: 2, node_count: 2 })]
    #[case::to_out_of_bounds(2, &[[0, 2]], GraphError::NodeOutOfBounds { node: 2, node_count: 2 })]
    #[case::out_of_bounds_loop(2, &[[2, 2]], GraphError::NodeOutOfBounds { node: 2, node_count: 2 })]
    #[case::self_loop(2, &[[1, 1]], GraphError::SelfLoop { node: 1 })]
    #[case::parallel_edge(3, &[[0, 2], [0, 1], [0, 2]], GraphError::ParallelEdge { endpoints: [0, 2] })]
    fn test_digraph_try_new_error(
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
        #[case] expected: GraphError,
    ) {
        assert_eq!(DiGraph::try_new(node_count, edges), Err(expected));
    }

    #[rstest]
    #[case::from_out_of_bounds(2, &[[2, 0]])]
    #[case::to_out_of_bounds(2, &[[0, 2]])]
    #[case::self_loop(2, &[[0, 0]])]
    #[case::parallel_edge(2, &[[0, 1], [0, 1]])]
    #[should_panic]
    fn test_digraph_new_error(#[case] node_count: usize, #[case] edges: &[[u32; 2]]) {
        DiGraph::new(node_count, edges);
    }
}
