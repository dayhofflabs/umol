//! Connected components.
//!
//! Unrestricted breadth-first traversal supplies sorted component node sets.

use std::ops::ControlFlow;

use crate::algorithms::traversal::BreadthFirstEvent;
use crate::graph::{Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConnectedComponentsAlgorithm {
    Bfs,
}

impl Graph {
    /// Visits each connected component as a borrowed, sorted node slice.
    ///
    /// # Semantic properties
    ///
    /// On completion, the emitted components partition all nodes into maximal
    /// connected sets. Each slice is sorted by node id; components are ordered
    /// by their least node id. Loops and parallel edges do not affect membership.
    /// Empty graphs emit no components.
    ///
    /// Bfs collects one component through unrestricted breadth-first traversal
    /// and visits it at FinishTree. A visitor's [`ControlFlow::Break`] returns
    /// immediately, before another component is explored. Callback panics or
    /// nontermination remain caller-owned. The slice is borrowed for the call;
    /// copy it to retain the component.
    pub fn visit_connected_components<B, V>(
        &self,
        algorithm: ConnectedComponentsAlgorithm,
        mut visitor: V,
    ) -> ControlFlow<B>
    where
        V: FnMut(&[NodeId]) -> ControlFlow<B>,
    {
        match algorithm {
            ConnectedComponentsAlgorithm::Bfs => {
                let mut component = Vec::new();
                self.visit_breadth_first(self.node_ids(), None, |event| {
                    match event {
                        BreadthFirstEvent::Discover { node, .. } => component.push(node),
                        BreadthFirstEvent::FinishTree { .. } => {
                            component.sort_unstable();
                            visitor(&component)?;
                            component.clear();
                        }
                        BreadthFirstEvent::Finish { .. } => (),
                    }
                    ControlFlow::Continue(())
                })
            }
        }
    }

    /// Collects connected components as sorted node vectors.
    ///
    /// # Semantic properties
    ///
    /// The result is exactly the sequence emitted by
    /// [`Graph::visit_connected_components`] with the same algorithm and a
    /// visitor that always continues, preserving component and member order.
    pub fn enumerate_connected_components(
        &self,
        alg: ConnectedComponentsAlgorithm,
    ) -> Vec<Vec<NodeId>> {
        let mut components = Vec::new();
        let _: ControlFlow<()> = self.visit_connected_components(alg, |component| {
            components.push(component.to_vec());
            ControlFlow::Continue(())
        });
        components
    }
}

#[cfg(test)]
mod tests {
    use std::ops::ControlFlow;

    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::ConnectedComponentsAlgorithm::Bfs;
    use crate::graph::{Graph, NodeId};

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    #[rstest]
    #[case::empty(Graph::new(0, &[]), vec![])]
    #[case::isolated(Graph::new(3, &[]), vec![vec![NodeId(0)], vec![NodeId(1)], vec![NodeId(2)]])]
    #[case::interleaved(
        Graph::new(6, &[[0, 4], [4, 2], [2, 1], [3, 5]]),
        vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(4)], vec![NodeId(3), NodeId(5)]]
    )]
    #[case::loops_parallel(
        Graph::new(4, &[[0, 0], [0, 2], [0, 2], [1, 3], [3, 3]]),
        vec![vec![NodeId(0), NodeId(2)], vec![NodeId(1), NodeId(3)]]
    )]
    fn test_graph_visit_connected_components(
        #[case] graph: Graph,
        #[case] expected: Vec<Vec<NodeId>>,
    ) {
        let mut actual = Vec::new();
        let result = graph.visit_connected_components(Bfs, |component| {
            actual.push(component.to_vec());
            ControlFlow::<()>::Continue(())
        });
        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::first(0)]
    #[case::middle(1)]
    #[case::last(2)]
    fn test_graph_visit_connected_components_break(#[case] stop: usize) {
        let graph = Graph::new(7, &[[0, 4], [4, 2], [1, 5], [3, 6]]);
        let expected = [
            vec![NodeId(0), NodeId(2), NodeId(4)],
            vec![NodeId(1), NodeId(5)],
            vec![NodeId(3), NodeId(6)],
        ];
        let mut actual = Vec::new();
        let result = graph.visit_connected_components(Bfs, |component| {
            actual.push(component.to_vec());
            if actual.len() == stop + 1 {
                ControlFlow::Break(component.to_vec())
            } else {
                ControlFlow::Continue(())
            }
        });
        assert_eq!(result, ControlFlow::Break(expected[stop].clone()));
        assert_eq!(actual, expected[..=stop]);
    }

    #[rstest]
    #[case::empty(0, vec![], vec![])]
    #[case::isolated(3, vec![], vec![vec![n(0)], vec![n(1)], vec![n(2)]])]
    #[case::single_edge(2, vec![[0, 1]], vec![vec![n(0), n(1)]])]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], vec![vec![n(0), n(1), n(2)]])]
    #[case::two_components(4, vec![[0, 1], [2, 3]], vec![vec![n(0), n(1)], vec![n(2), n(3)]])]
    fn test_graph_enumerate_connected_components(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected: Vec<Vec<NodeId>>,
    ) {
        let g = Graph::new(node_count, &edges);
        assert_eq!(g.enumerate_connected_components(Bfs), expected);
    }
}
