//! Depth-first event order is checked against recursive DFS using sets over
//! generated undirected multigraphs of at most eight nodes and twenty edges.
//! Ordered candidate roots include subsets and duplicates; callback adjacency
//! can reverse incidence order and retain sparse edge ids. Graph's CSR adapter
//! is checked separately. Early breaks must produce the reference prefix.
//! Raw finite callback tables exercise panic freedom without asserting a useful
//! traversal for inconsistent connectivity.

use std::collections::HashSet;
use std::ops::ControlFlow;

use proptest::prelude::*;
use rstest::rstest;
use umol_graph_core::{visit_depth_first, DepthFirstEvent, EdgeId, Neighbor, NodeId};

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
