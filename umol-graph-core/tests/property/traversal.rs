//! Depth-first event order is checked against recursive DFS using sets over
//! generated undirected multigraphs of at most eight nodes and twenty edges.
//! Ordered candidate roots include subsets and duplicates; callback adjacency
//! can reverse incidence order and retain sparse edge ids. Graph's CSR adapter
//! is checked separately. Early breaks must produce the reference prefix.
//! Raw finite callback tables exercise panic freedom without asserting a useful
//! traversal for inconsistent connectivity.
//!
//! Breadth-first events use a layer-by-layer reference over the same bounded
//! multigraph domain. Independent edge relaxation checks shortest distances and
//! reached sets for sequential roots, excluding previous trees. Depth limits,
//! reversed incidence order, sparse edge ids, and early breaks are included.
//! Event collectors must equal both the complete visitor sequence and its
//! definition-level reference for the same roots and depth limit.

use std::collections::{BTreeMap, HashSet};
use std::ops::ControlFlow;

use proptest::prelude::*;
use rstest::rstest;
use umol_graph_core::{
    visit_breadth_first, visit_depth_first, BreadthFirstEvent, DepthFirstEvent, EdgeId, Neighbor,
    NodeId,
};

use super::strategy::graph_with_edge_multiset;

fn reference_depth_first(adjacency: &[Vec<Neighbor>], roots: &[NodeId]) -> Vec<DepthFirstEvent> {
    let mut nodes = HashSet::new();
    let mut edges = HashSet::new();
    let mut events = Vec::new();
    for &root in roots {
        if nodes.contains(&root) {
            continue;
        }
        reference_visit(adjacency, root, None, &mut nodes, &mut edges, &mut events);
        events.push(DepthFirstEvent::FinishTree { root });
    }
    events
}

fn reference_visit(
    adjacency: &[Vec<Neighbor>],
    node: NodeId,
    parent: Option<Neighbor>,
    nodes: &mut HashSet<NodeId>,
    edges: &mut HashSet<EdgeId>,
    events: &mut Vec<DepthFirstEvent>,
) {
    nodes.insert(node);
    events.push(DepthFirstEvent::Discover { node, parent });
    for &Neighbor { node: to, edge } in &adjacency[node.index()] {
        if !edges.insert(edge) {
            continue;
        }
        if nodes.contains(&to) {
            events.push(DepthFirstEvent::NonTreeEdge {
                from: node,
                to,
                edge,
            });
        } else {
            reference_visit(
                adjacency,
                to,
                Some(Neighbor { node, edge }),
                nodes,
                edges,
                events,
            );
        }
    }
    events.push(DepthFirstEvent::Finish { node });
}

fn reference_breadth_first(
    adjacency: &[Vec<Neighbor>],
    roots: &[NodeId],
    max_depth: Option<usize>,
) -> Vec<BreadthFirstEvent> {
    let mut visited = HashSet::new();
    let mut events = Vec::new();
    for &root in roots {
        if !visited.insert(root) {
            continue;
        }
        events.push(BreadthFirstEvent::Discover {
            node: root,
            parent: None,
            depth: 0,
        });
        let mut layer = vec![root];
        for depth in 0..adjacency.len() {
            let mut next_layer = Vec::new();
            for node in layer {
                if max_depth.is_none_or(|limit| depth < limit) {
                    for &Neighbor { node: to, edge } in &adjacency[node.index()] {
                        if visited.insert(to) {
                            next_layer.push(to);
                            events.push(BreadthFirstEvent::Discover {
                                node: to,
                                parent: Some(Neighbor { node, edge }),
                                depth: depth + 1,
                            });
                        }
                    }
                }
                events.push(BreadthFirstEvent::Finish { node, depth });
            }
            if next_layer.is_empty() {
                break;
            }
            layer = next_layer;
        }
        events.push(BreadthFirstEvent::FinishTree { root });
    }
    events
}

fn reference_distances(
    node_count: usize,
    edges: &[[u32; 2]],
    root: NodeId,
    excluded: &HashSet<NodeId>,
) -> Vec<Option<usize>> {
    let mut distances = vec![None; node_count];
    distances[root.index()] = Some(0);
    for _ in 1..node_count {
        let previous = distances.clone();
        for &[a, b] in edges {
            if excluded.contains(&NodeId(a)) || excluded.contains(&NodeId(b)) {
                continue;
            }
            for (from, to) in [(a as usize, b as usize), (b as usize, a as usize)] {
                if let Some(distance) = previous[from] {
                    distances[to] =
                        Some(distances[to].map_or(distance + 1, |old| old.min(distance + 1)));
                }
            }
        }
    }
    distances
}

#[rstest]
fn test_visit_depth_first() {
    proptest!(|(
        (graph, edges) in graph_with_edge_multiset(8, 20),
        candidates in prop::collection::vec(0u32..8, 0..16),
        reverse in any::<bool>(),
        stop in 0usize..64,
    )| {
        let mut adjacency = vec![Vec::new(); graph.node_count()];
        for (index, &[a, b]) in edges.iter().enumerate() {
            let edge = EdgeId((2 * index + 1) as u32);
            adjacency[a as usize].push(Neighbor { node: NodeId(b), edge });
            adjacency[b as usize].push(Neighbor { node: NodeId(a), edge });
        }
        if reverse {
            for neighbors in &mut adjacency {
                neighbors.reverse();
            }
        }
        let roots: Vec<_> = candidates.into_iter()
            .filter(|&node| (node as usize) < graph.node_count()).map(NodeId).collect();
        let expected = reference_depth_first(&adjacency, &roots);
        let mut actual = Vec::new();
        let result = visit_depth_first(
            graph.node_count(), 2 * edges.len(), roots.iter().copied(),
            |node| adjacency[node.index()].iter().copied(),
            |event| {
                actual.push(event);
                ControlFlow::<()>::Continue(())
            },
        );
        prop_assert_eq!(result, ControlFlow::Continue(()));
        prop_assert_eq!(actual, expected.clone());

        let mut prefix = Vec::new();
        let result = visit_depth_first(
            graph.node_count(), 2 * edges.len(), roots,
            |node| adjacency[node.index()].iter().copied(),
            |event| {
                prefix.push(event);
                if prefix.len() == stop + 1 {
                    ControlFlow::Break(event)
                } else {
                    ControlFlow::Continue(())
                }
            },
        );
        prop_assert_eq!(result, expected.get(stop).copied().map_or(
            ControlFlow::Continue(()), ControlFlow::Break));
        prop_assert_eq!(prefix, expected[..(stop + 1).min(expected.len())].to_vec());
    });
}

#[rstest]
fn test_visit_depth_first_partial() {
    proptest!(|(
        node_bound in 0usize..8,
        edge_bound in 0usize..8,
        candidates in prop::collection::vec(prop_oneof![0u32..12, Just(u32::MAX)], 0..16),
        table in prop::collection::vec(
            prop::collection::vec(
                (prop_oneof![0u32..12, Just(u32::MAX)], prop_oneof![0u32..12, Just(u32::MAX)]),
                0..16,
            ), 0..8,
        ),
    )| {
        let result = visit_depth_first(
            node_bound, edge_bound, candidates.into_iter().map(NodeId),
            |node| table.get(node.index()).into_iter().flatten()
                .map(|&(node, edge)| Neighbor { node: NodeId(node), edge: EdgeId(edge) }),
            |_| ControlFlow::<()>::Continue(()),
        );
        prop_assert_eq!(result, ControlFlow::Continue(()));
    });
}

#[rstest]
fn test_graph_visit_depth_first() {
    proptest!(|(
        (graph, _) in graph_with_edge_multiset(8, 20),
        candidates in prop::collection::vec(0u32..8, 0..16),
    )| {
        let roots: Vec<_> = candidates.into_iter()
            .filter(|&node| (node as usize) < graph.node_count()).map(NodeId).collect();
        let adjacency: Vec<_> = graph.node_ids().map(|node| graph.neighbors(node).to_vec()).collect();
        let expected = reference_depth_first(&adjacency, &roots);
        let mut actual = Vec::new();
        let result = graph.visit_depth_first(roots, |event| {
            actual.push(event);
            ControlFlow::<()>::Continue(())
        });
        prop_assert_eq!(result, ControlFlow::Continue(()));
        prop_assert_eq!(actual, expected);
    });
}

#[rstest]
fn test_graph_enumerate_depth_first_events() {
    proptest!(|(
        (graph, _) in graph_with_edge_multiset(8, 20),
        candidates in prop::collection::vec(0u32..8, 0..16),
    )| {
        let roots: Vec<_> = candidates.into_iter()
            .filter(|&node| (node as usize) < graph.node_count()).map(NodeId).collect();
        let adjacency: Vec<_> = graph.node_ids().map(|node| graph.neighbors(node).to_vec()).collect();
        let expected = reference_depth_first(&adjacency, &roots);
        let mut visited = Vec::new();
        let result = graph.visit_depth_first(roots.iter().copied(), |event| {
            visited.push(event);
            ControlFlow::<()>::Continue(())
        });
        let actual = graph.enumerate_depth_first_events(roots);
        prop_assert_eq!(result, ControlFlow::Continue(()));
        prop_assert_eq!(&actual, &expected);
        prop_assert_eq!(actual, visited);
    });
}

#[rstest]
fn test_visit_breadth_first() {
    proptest!(|(
        (graph, edges) in graph_with_edge_multiset(8, 20),
        candidates in prop::collection::vec(0u32..8, 0..16),
        reverse in any::<bool>(),
        max_depth in prop::option::of(prop_oneof![0usize..10, Just(usize::MAX)]),
        stop in 0usize..32,
    )| {
        let mut adjacency = vec![Vec::new(); graph.node_count()];
        for (index, &[a, b]) in edges.iter().enumerate() {
            let edge = EdgeId((2 * index + 1) as u32);
            adjacency[a as usize].push(Neighbor { node: NodeId(b), edge });
            adjacency[b as usize].push(Neighbor { node: NodeId(a), edge });
        }
        if reverse {
            for neighbors in &mut adjacency {
                neighbors.reverse();
            }
        }
        let roots: Vec<_> = candidates.into_iter()
            .filter(|&node| (node as usize) < graph.node_count()).map(NodeId).collect();
        let expected = reference_breadth_first(&adjacency, &roots, max_depth);
        let mut actual = Vec::new();
        let result = visit_breadth_first(
            graph.node_count(), roots.iter().copied(), max_depth,
            |node| adjacency[node.index()].iter().copied(),
            |event| {
                actual.push(event);
                ControlFlow::<()>::Continue(())
            },
        );
        prop_assert_eq!(result, ControlFlow::Continue(()));
        prop_assert_eq!(actual, expected.clone());

        let mut prefix = Vec::new();
        let result = visit_breadth_first(
            graph.node_count(), roots, max_depth,
            |node| adjacency[node.index()].iter().copied(),
            |event| {
                prefix.push(event);
                if prefix.len() == stop + 1 {
                    ControlFlow::Break(event)
                } else {
                    ControlFlow::Continue(())
                }
            },
        );
        prop_assert_eq!(result, expected.get(stop).copied().map_or(
            ControlFlow::Continue(()), ControlFlow::Break));
        prop_assert_eq!(prefix, expected[..(stop + 1).min(expected.len())].to_vec());
    });
}

#[rstest]
fn test_visit_breadth_first_partial() {
    proptest!(|(
        node_bound in 0usize..8,
        max_depth in prop::option::of(prop_oneof![0usize..10, Just(usize::MAX)]),
        candidates in prop::collection::vec(prop_oneof![0u32..12, Just(u32::MAX)], 0..16),
        table in prop::collection::vec(
            prop::collection::vec(
                (prop_oneof![0u32..12, Just(u32::MAX)], any::<u32>()), 0..16,
            ), 0..8,
        ),
    )| {
        let result = visit_breadth_first(
            node_bound, candidates.into_iter().map(NodeId), max_depth,
            |node| table.get(node.index()).into_iter().flatten()
                .map(|&(node, edge)| Neighbor { node: NodeId(node), edge: EdgeId(edge) }),
            |_| ControlFlow::<()>::Continue(()),
        );
        prop_assert_eq!(result, ControlFlow::Continue(()));
    });
}

#[rstest]
fn test_graph_visit_breadth_first() {
    proptest!(|(
        (graph, edges) in graph_with_edge_multiset(8, 20),
        candidates in prop::collection::vec(0u32..8, 0..16),
        max_depth in prop::option::of(prop_oneof![0usize..10, Just(usize::MAX)]),
    )| {
        let roots: Vec<_> = candidates.into_iter()
            .filter(|&node| (node as usize) < graph.node_count()).map(NodeId).collect();
        let mut excluded = HashSet::new();
        let mut expected_trees = Vec::new();
        for &root in &roots {
            if excluded.contains(&root) {
                continue;
            }
            let reached: BTreeMap<_, _> = reference_distances(graph.node_count(), &edges, root, &excluded)
                .into_iter().enumerate().filter_map(|(node, distance)| {
                    distance.filter(|&depth| max_depth.is_none_or(|limit| depth <= limit))
                        .map(|depth| (NodeId(node as u32), depth))
                }).collect();
            excluded.extend(reached.keys().copied());
            expected_trees.push((root, reached));
        }
        let adjacency: Vec<_> = graph.node_ids().map(|node| graph.neighbors(node).to_vec()).collect();
        let expected = reference_breadth_first(&adjacency, &roots, max_depth);
        let mut actual = Vec::new();
        let result = graph.visit_breadth_first(roots, max_depth, |event| {
            actual.push(event);
            ControlFlow::<()>::Continue(())
        });
        prop_assert_eq!(result, ControlFlow::Continue(()));
        prop_assert_eq!(&actual, &expected);
        let mut actual_trees = Vec::new();
        let mut reached = BTreeMap::new();
        for event in actual {
            match event {
                BreadthFirstEvent::Discover { node, depth, .. } => {
                    prop_assert_eq!(reached.insert(node, depth), None);
                }
                BreadthFirstEvent::FinishTree { root } => {
                    actual_trees.push((root, reached));
                    reached = BTreeMap::new();
                }
                BreadthFirstEvent::Finish { .. } => (),
            }
        }
        prop_assert_eq!(actual_trees, expected_trees);
    });
}

#[rstest]
fn test_graph_enumerate_breadth_first_events() {
    proptest!(|(
        (graph, _) in graph_with_edge_multiset(8, 20),
        candidates in prop::collection::vec(0u32..8, 0..16),
        max_depth in prop::option::of(prop_oneof![0usize..10, Just(usize::MAX)]),
    )| {
        let roots: Vec<_> = candidates.into_iter()
            .filter(|&node| (node as usize) < graph.node_count()).map(NodeId).collect();
        let adjacency: Vec<_> = graph.node_ids().map(|node| graph.neighbors(node).to_vec()).collect();
        let expected = reference_breadth_first(&adjacency, &roots, max_depth);
        let mut visited = Vec::new();
        let result = graph.visit_breadth_first(roots.iter().copied(), max_depth, |event| {
            visited.push(event);
            ControlFlow::<()>::Continue(())
        });
        let actual = graph.enumerate_breadth_first_events(roots, max_depth);
        prop_assert_eq!(result, ControlFlow::Continue(()));
        prop_assert_eq!(&actual, &expected);
        prop_assert_eq!(actual, visited);
    });
}
