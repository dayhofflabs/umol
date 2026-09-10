//! Graph traversal events and distance-limited neighborhoods.

use std::collections::VecDeque;
use std::ops::ControlFlow;

use crate::graph::{EdgeId, Graph, Neighbor, NodeId};

/// An event in an undirected depth-first traversal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DepthFirstEvent {
    /// First encounter with a node. `parent` names the parent node and tree edge;
    /// it is `None` for a traversal root.
    Discover {
        node: NodeId,
        parent: Option<Neighbor>,
    },
    /// First encounter with an edge whose endpoints have already been discovered.
    /// Each non-tree edge, including a loop, is reported once by edge identity.
    NonTreeEdge {
        from: NodeId,
        to: NodeId,
        edge: EdgeId,
    },
    /// All neighbors and descendants of a node have been processed.
    Finish { node: NodeId },
    /// The traversal tree rooted at `root` has finished.
    FinishTree { root: NodeId },
}

/// Visits undirected connectivity in depth-first order using an explicit stack.
///
/// # Semantic properties
///
/// `roots` supplies ordered candidate roots; already reached nodes are skipped.
/// Each node's neighbor iterator is suspended while its descendants are visited,
/// preserving iterator order. On completion, each reached node has one nested
/// discovery/finish pair, each non-tree edge is reported once, and each completed
/// root is followed by [`DepthFirstEvent::FinishTree`].
///
/// Connectivity must remain fixed during traversal. Every node id must be below
/// `node_bound`, every edge id below `edge_bound`, and each undirected edge must
/// have a stable identity shared by its incidences. Parallel edges have distinct
/// identities; loops may have repeated incidences. Filtered connectivity can
/// retain sparse original ids by supplying the original bounds.
///
/// Inconsistent inputs have no traversal correctness guarantee, but do not cause
/// internal indexing panics. There is no connectivity validation pass. Panics or
/// nontermination in supplied callbacks and iterators remain caller-owned.
///
/// A visitor's [`ControlFlow::Break`] is returned immediately, without further
/// callback or iterator calls or synthetic finish events. State uses space
/// proportional to the bounds plus the active path and its suspended iterators.
pub fn visit_depth_first<R, N, I, V, B>(
    node_bound: usize,
    edge_bound: usize,
    roots: R,
    neighbors: N,
    mut visitor: V,
) -> ControlFlow<B>
where
    R: IntoIterator<Item = NodeId>,
    N: Fn(NodeId) -> I,
    I: Iterator<Item = Neighbor>,
    V: FnMut(DepthFirstEvent) -> ControlFlow<B>,
{
    let mut visited_nodes = vec![false; node_bound];
    let mut visited_edges = vec![false; edge_bound];
    let mut stack = Vec::new();
    for root in roots {
        let Some(visited) = visited_nodes.get_mut(root.index()) else {
            continue;
        };
        if *visited {
            continue;
        }
        *visited = true;
        visitor(DepthFirstEvent::Discover {
            node: root,
            parent: None,
        })?;
        stack.push((root, neighbors(root)));
        while let Some((node, pending)) = stack.last_mut() {
            let from = *node;
            if let Some(Neighbor { node: to, edge }) = pending.next() {
                let (Some(node_visited), Some(edge_visited)) = (
                    visited_nodes.get_mut(to.index()),
                    visited_edges.get_mut(edge.index()),
                ) else {
                    continue;
                };
                if *edge_visited {
                    continue;
                }
                *edge_visited = true;
                if *node_visited {
                    visitor(DepthFirstEvent::NonTreeEdge { from, to, edge })?;
                } else {
                    *node_visited = true;
                    visitor(DepthFirstEvent::Discover {
                        node: to,
                        parent: Some(Neighbor { node: from, edge }),
                    })?;
                    stack.push((to, neighbors(to)));
                }
            } else {
                stack.pop();
                visitor(DepthFirstEvent::Finish { node: from })?;
            }
        }
        visitor(DepthFirstEvent::FinishTree { root })?;
    }
    ControlFlow::Continue(())
}

/// An event in an undirected breadth-first traversal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BreadthFirstEvent {
    /// First encounter with a node, at its distance from the current root.
    /// `parent` names the parent node and first reaching edge; roots have no parent.
    Discover {
        node: NodeId,
        parent: Option<Neighbor>,
        depth: usize,
    },
    /// The node's permitted neighbor expansion has finished, including at a depth limit.
    Finish { node: NodeId, depth: usize },
    /// The traversal queue for `root` is empty.
    FinishTree { root: NodeId },
}

/// Visits undirected connectivity in breadth-first order using a FIFO queue.
///
/// # Semantic properties
///
/// Candidate `roots` are processed sequentially in supplied order, skipping
/// already reached nodes. Neighbors are examined in iterator order. Nodes are
/// marked and discovered when enqueued; their parent is the first reaching
/// node/edge pair. Each reached node receives one Discover and one Finish.
/// Finish follows its neighbor expansion; FinishTree follows an empty queue,
/// before the next candidate root is considered.
///
/// `max_depth` suppresses neighbor expansion at the limit, including the neighbor
/// callback, but nodes there still receive Discover and Finish. `None` is
/// unrestricted. A limited tree need not cover a whole connected component.
/// Depths are shortest distances from the current root in the graph excluding
/// nodes reached by previous trees. This is sequential traversal, not multi-source
/// BFS; all trees share visitation state.
///
/// Connectivity must remain fixed, with node ids below `node_bound` and stable
/// edge identities shared by each undirected edge's incidences. Loops and parallel
/// edges are supported. Edge ids are carried into parent events without indexing
/// edge state, so no edge bound is required.
///
/// Inconsistent inputs have no correctness guarantee but do not cause internal
/// indexing panics. There is no validation pass. Panics and nontermination in
/// supplied callbacks or iterators remain caller-owned. A visitor's
/// [`ControlFlow::Break`] returns immediately without further callbacks, iterator
/// calls, or synthetic finish events. State uses space proportional to the node
/// bound plus the queue.
pub fn visit_breadth_first<R, N, I, V, B>(
    node_bound: usize,
    roots: R,
    max_depth: Option<usize>,
    neighbors: N,
    mut visitor: V,
) -> ControlFlow<B>
where
    R: IntoIterator<Item = NodeId>,
    N: Fn(NodeId) -> I,
    I: Iterator<Item = Neighbor>,
    V: FnMut(BreadthFirstEvent) -> ControlFlow<B>,
{
    let mut visited_nodes = vec![false; node_bound];
    let mut queue = VecDeque::new();
    let max_depth = max_depth.unwrap_or(usize::MAX);
    for root in roots {
        let Some(visited) = visited_nodes.get_mut(root.index()) else {
            continue;
        };
        if *visited {
            continue;
        }
        *visited = true;
        queue.push_back((root, 0));
        visitor(BreadthFirstEvent::Discover {
            node: root,
            parent: None,
            depth: 0,
        })?;
        while let Some((node, depth)) = queue.pop_front() {
            if depth < max_depth {
                for Neighbor { node: to, edge } in neighbors(node) {
                    let Some(visited) = visited_nodes.get_mut(to.index()) else {
                        continue;
                    };
                    if *visited {
                        continue;
                    }
                    *visited = true;
                    queue.push_back((to, depth + 1));
                    visitor(BreadthFirstEvent::Discover {
                        node: to,
                        parent: Some(Neighbor { node, edge }),
                        depth: depth + 1,
                    })?;
                }
            }
            visitor(BreadthFirstEvent::Finish { node, depth })?;
        }
        visitor(BreadthFirstEvent::FinishTree { root })?;
    }
    ControlFlow::Continue(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TraversalAlgorithm {
    Bfs,
}

impl Graph {
    /// Visits candidate roots in the supplied order and neighbors in CSR order.
    ///
    /// Pass [`Graph::node_ids`] to visit every connected component. Already
    /// reached roots are skipped. See [`visit_depth_first`] for event ordering
    /// and immediate [`ControlFlow::Break`] behavior.
    pub fn visit_depth_first<R, V, B>(&self, roots: R, visitor: V) -> ControlFlow<B>
    where
        R: IntoIterator<Item = NodeId>,
        V: FnMut(DepthFirstEvent) -> ControlFlow<B>,
    {
        visit_depth_first(
            self.node_count(),
            self.edge_count(),
            roots,
            |node| self.neighbors(node).iter().copied(),
            visitor,
        )
    }

    /// Visits each candidate root's tree in breadth-first order using CSR neighbors.
    ///
    /// Pass [`Graph::node_ids`] and no depth limit to visit all components.
    /// See [`visit_breadth_first`] for depth limits, shared visitation across
    /// roots, event ordering, and immediate [`ControlFlow::Break`] behavior.
    pub fn visit_breadth_first<R, V, B>(
        &self,
        roots: R,
        max_depth: Option<usize>,
        visitor: V,
    ) -> ControlFlow<B>
    where
        R: IntoIterator<Item = NodeId>,
        V: FnMut(BreadthFirstEvent) -> ControlFlow<B>,
    {
        visit_breadth_first(
            self.node_count(),
            roots,
            max_depth,
            |node| self.neighbors(node).iter().copied(),
            visitor,
        )
    }

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
    use std::cell::{Cell, RefCell};
    use std::iter;

    use rstest::rstest;

    use super::*;
    use DepthFirstEvent::{Discover, Finish, FinishTree, NonTreeEdge};

    #[rstest]
    #[case::neighbor_order(
        3, 6, vec![NodeId(0)],
        vec![
            vec![Neighbor { node: NodeId(2), edge: EdgeId(5) }, Neighbor { node: NodeId(1), edge: EdgeId(2) }],
            vec![Neighbor { node: NodeId(0), edge: EdgeId(2) }],
            vec![Neighbor { node: NodeId(0), edge: EdgeId(5) }],
        ],
        vec![
            Discover { node: NodeId(0), parent: None },
            Discover { node: NodeId(2), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(5) }) },
            Finish { node: NodeId(2) },
            Discover { node: NodeId(1), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(2) }) },
            Finish { node: NodeId(1) }, Finish { node: NodeId(0) }, FinishTree { root: NodeId(0) },
        ]
    )]
    fn test_visit_depth_first(
        #[case] node_bound: usize,
        #[case] edge_bound: usize,
        #[case] roots: Vec<NodeId>,
        #[case] adjacency: Vec<Vec<Neighbor>>,
        #[case] expected: Vec<DepthFirstEvent>,
    ) {
        let mut actual = Vec::new();
        let result = visit_depth_first(
            node_bound,
            edge_bound,
            roots,
            |node| adjacency[node.index()].iter().copied(),
            |event| {
                actual.push(event);
                ControlFlow::<()>::Continue(())
            },
        );
        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::discover(0)]
    #[case::child_discover(1)]
    #[case::non_tree_edge(3)]
    #[case::finish(4)]
    #[case::finish_tree(7)]
    fn test_visit_depth_first_break(#[case] stop: usize) {
        let graph = Graph::new(4, &[[0, 1], [0, 2], [1, 2]]);
        let expected = [
            Discover {
                node: NodeId(0),
                parent: None,
            },
            Discover {
                node: NodeId(1),
                parent: Some(Neighbor {
                    node: NodeId(0),
                    edge: EdgeId(0),
                }),
            },
            Discover {
                node: NodeId(2),
                parent: Some(Neighbor {
                    node: NodeId(1),
                    edge: EdgeId(2),
                }),
            },
            NonTreeEdge {
                from: NodeId(2),
                to: NodeId(0),
                edge: EdgeId(1),
            },
            Finish { node: NodeId(2) },
            Finish { node: NodeId(1) },
            Finish { node: NodeId(0) },
            FinishTree { root: NodeId(0) },
        ];
        let stopped = Cell::new(false);
        let mut candidates = graph.node_ids();
        let roots = iter::from_fn(|| {
            assert!(!stopped.get());
            candidates.next()
        });
        let mut actual = Vec::new();
        let result = visit_depth_first(
            graph.node_count(),
            graph.edge_count(),
            roots,
            |node| {
                assert!(!stopped.get());
                let mut pending = graph.neighbors(node).iter().copied();
                let stopped = &stopped;
                iter::from_fn(move || {
                    assert!(!stopped.get());
                    pending.next()
                })
            },
            |event| {
                assert!(!stopped.get());
                actual.push(event);
                if actual.len() == stop + 1 {
                    stopped.set(true);
                    ControlFlow::Break(event)
                } else {
                    ControlFlow::Continue(())
                }
            },
        );
        assert_eq!(result, ControlFlow::Break(expected[stop]));
        assert_eq!(actual, expected[..=stop]);
    }

    #[rstest]
    #[case::neighbor_order(
        vec![
            vec![Neighbor { node: NodeId(2), edge: EdgeId(u32::MAX) }, Neighbor { node: NodeId(1), edge: EdgeId(7) }],
            vec![Neighbor { node: NodeId(0), edge: EdgeId(7) }],
            vec![Neighbor { node: NodeId(0), edge: EdgeId(u32::MAX) }],
        ],
        vec![
            BreadthFirstEvent::Discover { node: NodeId(0), parent: None, depth: 0 },
            BreadthFirstEvent::Discover { node: NodeId(2), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(u32::MAX) }), depth: 1 },
            BreadthFirstEvent::Discover { node: NodeId(1), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(7) }), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(0), depth: 0 },
            BreadthFirstEvent::Finish { node: NodeId(2), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(1), depth: 1 },
            BreadthFirstEvent::FinishTree { root: NodeId(0) },
        ]
    )]
    fn test_visit_breadth_first(
        #[case] adjacency: Vec<Vec<Neighbor>>,
        #[case] expected: Vec<BreadthFirstEvent>,
    ) {
        let mut actual = Vec::new();
        let result = visit_breadth_first(
            adjacency.len(),
            [NodeId(0)],
            None,
            |node| adjacency[node.index()].iter().copied(),
            |event| {
                actual.push(event);
                ControlFlow::<()>::Continue(())
            },
        );
        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::discover(0)]
    #[case::child_discover(1)]
    #[case::root_finish(3)]
    #[case::child_finish(4)]
    #[case::finish_tree(6)]
    fn test_visit_breadth_first_break(#[case] stop: usize) {
        let graph = Graph::new(4, &[[0, 1], [0, 2], [1, 2]]);
        let expected = [
            BreadthFirstEvent::Discover {
                node: NodeId(0),
                parent: None,
                depth: 0,
            },
            BreadthFirstEvent::Discover {
                node: NodeId(1),
                parent: Some(Neighbor {
                    node: NodeId(0),
                    edge: EdgeId(0),
                }),
                depth: 1,
            },
            BreadthFirstEvent::Discover {
                node: NodeId(2),
                parent: Some(Neighbor {
                    node: NodeId(0),
                    edge: EdgeId(1),
                }),
                depth: 1,
            },
            BreadthFirstEvent::Finish {
                node: NodeId(0),
                depth: 0,
            },
            BreadthFirstEvent::Finish {
                node: NodeId(1),
                depth: 1,
            },
            BreadthFirstEvent::Finish {
                node: NodeId(2),
                depth: 1,
            },
            BreadthFirstEvent::FinishTree { root: NodeId(0) },
        ];
        let stopped = Cell::new(false);
        let mut candidates = graph.node_ids();
        let roots = iter::from_fn(|| {
            assert!(!stopped.get());
            candidates.next()
        });
        let mut actual = Vec::new();
        let result = visit_breadth_first(
            graph.node_count(),
            roots,
            None,
            |node| {
                assert!(!stopped.get());
                let mut pending = graph.neighbors(node).iter().copied();
                let stopped = &stopped;
                iter::from_fn(move || {
                    assert!(!stopped.get());
                    pending.next()
                })
            },
            |event| {
                assert!(!stopped.get());
                actual.push(event);
                if actual.len() == stop + 1 {
                    stopped.set(true);
                    ControlFlow::Break(event)
                } else {
                    ControlFlow::Continue(())
                }
            },
        );
        assert_eq!(result, ControlFlow::Break(expected[stop]));
        assert_eq!(actual, expected[..=stop]);
    }

    #[rstest]
    #[case::zero(Some(0), vec![(NodeId(0), 0)], vec![])]
    #[case::one(Some(1), vec![(NodeId(0), 0), (NodeId(1), 1)], vec![NodeId(0)])]
    #[case::exact(Some(2), vec![(NodeId(0), 0), (NodeId(1), 1), (NodeId(2), 2)], vec![NodeId(0), NodeId(1)])]
    #[case::unrestricted(None, vec![(NodeId(0), 0), (NodeId(1), 1), (NodeId(2), 2)], vec![NodeId(0), NodeId(1), NodeId(2)])]
    #[case::maximum(Some(usize::MAX), vec![(NodeId(0), 0), (NodeId(1), 1), (NodeId(2), 2)], vec![NodeId(0), NodeId(1), NodeId(2)])]
    fn test_visit_breadth_first_depth(
        #[case] max_depth: Option<usize>,
        #[case] expected: Vec<(NodeId, usize)>,
        #[case] expected_expansions: Vec<NodeId>,
    ) {
        let graph = Graph::new(3, &[[0, 1], [1, 2]]);
        let expansions = RefCell::new(Vec::new());
        let mut discoveries = Vec::new();
        let mut finishes = Vec::new();
        let mut trees = Vec::new();
        let result = visit_breadth_first(
            graph.node_count(),
            [NodeId(0)],
            max_depth,
            |node| {
                expansions.borrow_mut().push(node);
                graph.neighbors(node).iter().copied()
            },
            |event| {
                match event {
                    BreadthFirstEvent::Discover { node, depth, .. } => {
                        discoveries.push((node, depth))
                    }
                    BreadthFirstEvent::Finish { node, depth } => finishes.push((node, depth)),
                    BreadthFirstEvent::FinishTree { root } => trees.push(root),
                }
                ControlFlow::<()>::Continue(())
            },
        );
        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(discoveries, expected);
        assert_eq!(finishes, expected);
        assert_eq!(trees, vec![NodeId(0)]);
        assert_eq!(expansions.into_inner(), expected_expansions);
    }

    #[rstest]
    #[case::empty(Graph::new(0, &[]), vec![], vec![])]
    #[case::no_roots(Graph::new(2, &[[0, 1]]), vec![], vec![])]
    #[case::candidate_order(
        Graph::new(3, &[]), vec![NodeId(2), NodeId(0), NodeId(2)],
        vec![Discover { node: NodeId(2), parent: None }, Finish { node: NodeId(2) }, FinishTree { root: NodeId(2) },
             Discover { node: NodeId(0), parent: None }, Finish { node: NodeId(0) }, FinishTree { root: NodeId(0) }]
    )]
    #[case::cycle(
        Graph::new(3, &[[0, 1], [0, 2], [1, 2]]), vec![NodeId(0), NodeId(1), NodeId(2)],
        vec![
            Discover { node: NodeId(0), parent: None },
            Discover { node: NodeId(1), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(0) }) },
            Discover { node: NodeId(2), parent: Some(Neighbor { node: NodeId(1), edge: EdgeId(2) }) },
            NonTreeEdge { from: NodeId(2), to: NodeId(0), edge: EdgeId(1) },
            Finish { node: NodeId(2) }, Finish { node: NodeId(1) }, Finish { node: NodeId(0) }, FinishTree { root: NodeId(0) },
        ]
    )]
    #[case::loops_parallel_isolated(
        Graph::new(3, &[[0, 0], [0, 1], [0, 1]]), vec![NodeId(0), NodeId(1), NodeId(2)],
        vec![
            Discover { node: NodeId(0), parent: None },
            NonTreeEdge { from: NodeId(0), to: NodeId(0), edge: EdgeId(0) },
            Discover { node: NodeId(1), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(1) }) },
            NonTreeEdge { from: NodeId(1), to: NodeId(0), edge: EdgeId(2) },
            Finish { node: NodeId(1) }, Finish { node: NodeId(0) }, FinishTree { root: NodeId(0) },
            Discover { node: NodeId(2), parent: None }, Finish { node: NodeId(2) }, FinishTree { root: NodeId(2) },
        ]
    )]
    fn test_graph_visit_depth_first(
        #[case] graph: Graph,
        #[case] roots: Vec<NodeId>,
        #[case] expected: Vec<DepthFirstEvent>,
    ) {
        let mut actual = Vec::new();
        let result = graph.visit_depth_first(roots, |event| {
            actual.push(event);
            ControlFlow::<()>::Continue(())
        });
        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::path(100_000)]
    fn test_graph_visit_depth_first_depth(#[case] count: u32) {
        let edges: Vec<_> = (1..count).map(|node| [node - 1, node]).collect();
        let graph = Graph::new(count as usize, &edges);
        let mut expected = (0..count)
            .map(|node| Discover {
                node: NodeId(node),
                parent: node.checked_sub(1).map(|parent| Neighbor {
                    node: NodeId(parent),
                    edge: EdgeId(parent),
                }),
            })
            .chain((0..count).rev().map(|node| Finish { node: NodeId(node) }))
            .chain([FinishTree { root: NodeId(0) }]);
        let result = graph.visit_depth_first(graph.node_ids(), |event| {
            assert_eq!(Some(event), expected.next());
            ControlFlow::<()>::Continue(())
        });
        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(expected.next(), None);
    }

    #[rstest]
    #[case::empty(Graph::new(0, &[]), vec![], None, vec![])]
    #[case::no_roots(Graph::new(2, &[[0, 1]]), vec![], None, vec![])]
    #[case::cycle(
        Graph::new(3, &[[0, 1], [0, 2], [1, 2]]), vec![NodeId(0), NodeId(1)], None,
        vec![
            BreadthFirstEvent::Discover { node: NodeId(0), parent: None, depth: 0 },
            BreadthFirstEvent::Discover { node: NodeId(1), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(0) }), depth: 1 },
            BreadthFirstEvent::Discover { node: NodeId(2), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(1) }), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(0), depth: 0 },
            BreadthFirstEvent::Finish { node: NodeId(1), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(2), depth: 1 },
            BreadthFirstEvent::FinishTree { root: NodeId(0) },
        ]
    )]
    #[case::loops_parallel_isolated(
        Graph::new(3, &[[0, 0], [0, 1], [0, 1]]), vec![NodeId(2), NodeId(0), NodeId(1), NodeId(2)], None,
        vec![
            BreadthFirstEvent::Discover { node: NodeId(2), parent: None, depth: 0 },
            BreadthFirstEvent::Finish { node: NodeId(2), depth: 0 },
            BreadthFirstEvent::FinishTree { root: NodeId(2) },
            BreadthFirstEvent::Discover { node: NodeId(0), parent: None, depth: 0 },
            BreadthFirstEvent::Discover { node: NodeId(1), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(1) }), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(0), depth: 0 },
            BreadthFirstEvent::Finish { node: NodeId(1), depth: 1 },
            BreadthFirstEvent::FinishTree { root: NodeId(0) },
        ]
    )]
    #[case::limited_trees(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)], Some(1),
        vec![
            BreadthFirstEvent::Discover { node: NodeId(0), parent: None, depth: 0 },
            BreadthFirstEvent::Discover { node: NodeId(1), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(0) }), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(0), depth: 0 },
            BreadthFirstEvent::Finish { node: NodeId(1), depth: 1 },
            BreadthFirstEvent::FinishTree { root: NodeId(0) },
            BreadthFirstEvent::Discover { node: NodeId(2), parent: None, depth: 0 },
            BreadthFirstEvent::Discover { node: NodeId(3), parent: Some(Neighbor { node: NodeId(2), edge: EdgeId(2) }), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(2), depth: 0 },
            BreadthFirstEvent::Finish { node: NodeId(3), depth: 1 },
            BreadthFirstEvent::FinishTree { root: NodeId(2) },
        ]
    )]
    #[case::limited_barrier(
        Graph::new(7, &[[0, 1], [1, 2], [2, 3], [4, 2], [4, 5], [5, 6], [6, 3]]), vec![NodeId(0), NodeId(4)], Some(2),
        vec![
            BreadthFirstEvent::Discover { node: NodeId(0), parent: None, depth: 0 },
            BreadthFirstEvent::Discover { node: NodeId(1), parent: Some(Neighbor { node: NodeId(0), edge: EdgeId(0) }), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(0), depth: 0 },
            BreadthFirstEvent::Discover { node: NodeId(2), parent: Some(Neighbor { node: NodeId(1), edge: EdgeId(1) }), depth: 2 },
            BreadthFirstEvent::Finish { node: NodeId(1), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(2), depth: 2 },
            BreadthFirstEvent::FinishTree { root: NodeId(0) },
            BreadthFirstEvent::Discover { node: NodeId(4), parent: None, depth: 0 },
            BreadthFirstEvent::Discover { node: NodeId(5), parent: Some(Neighbor { node: NodeId(4), edge: EdgeId(4) }), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(4), depth: 0 },
            BreadthFirstEvent::Discover { node: NodeId(6), parent: Some(Neighbor { node: NodeId(5), edge: EdgeId(5) }), depth: 2 },
            BreadthFirstEvent::Finish { node: NodeId(5), depth: 1 },
            BreadthFirstEvent::Finish { node: NodeId(6), depth: 2 },
            BreadthFirstEvent::FinishTree { root: NodeId(4) },
        ]
    )]
    fn test_graph_visit_breadth_first(
        #[case] graph: Graph,
        #[case] roots: Vec<NodeId>,
        #[case] max_depth: Option<usize>,
        #[case] expected: Vec<BreadthFirstEvent>,
    ) {
        let mut actual = Vec::new();
        let result = graph.visit_breadth_first(roots, max_depth, |event| {
            actual.push(event);
            ControlFlow::<()>::Continue(())
        });
        assert_eq!(result, ControlFlow::Continue(()));
        assert_eq!(actual, expected);
    }

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
