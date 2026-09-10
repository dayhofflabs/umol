//! Component visitation and enumeration must partition nodes into maximal
//! connected sets, with sorted members and component order by least node id.
//! Generated multigraphs have at most eight nodes and twenty edges. The reference
//! computes reachability by repeatedly closing a node set under edge incidence,
//! independently of the production BFS. Early termination must preserve prefixes.

use std::collections::BTreeSet;
use std::ops::ControlFlow;

use proptest::prelude::*;
use rstest::rstest;
use umol_graph_core::{ConnectedComponentsAlgorithm, NodeId};

use super::strategy::graph_with_edge_multiset;

fn reference_components(node_count: usize, edges: &[[u32; 2]]) -> Vec<Vec<NodeId>> {
    let mut reached = BTreeSet::new();
    let mut components = Vec::new();
    for root in (0..node_count).map(NodeId::from) {
        if reached.contains(&root) {
            continue;
        }
        let mut component = BTreeSet::from([root]);
        loop {
            let previous_size = component.len();
            for &[a, b] in edges {
                if component.contains(&NodeId(a)) || component.contains(&NodeId(b)) {
                    component.extend([NodeId(a), NodeId(b)]);
                }
            }
            if component.len() == previous_size {
                break;
            }
        }
        reached.extend(component.iter().copied());
        components.push(component.into_iter().collect());
    }
    components
}

#[rstest]
fn test_graph_visit_connected_components() {
    proptest!(|(
        (graph, edges) in graph_with_edge_multiset(8, 20),
        stop in 0usize..10,
    )| {
        let expected = reference_components(graph.node_count(), &edges);
        let mut actual = Vec::new();
        let result = graph.visit_connected_components(ConnectedComponentsAlgorithm::Bfs, |component| {
            actual.push(component.to_vec());
            ControlFlow::<()>::Continue(())
        });
        prop_assert_eq!(result, ControlFlow::Continue(()));
        prop_assert_eq!(&actual, &expected);
        prop_assert_eq!(graph.enumerate_connected_components(ConnectedComponentsAlgorithm::Bfs), actual);

        let mut prefix = Vec::new();
        let result = graph.visit_connected_components(ConnectedComponentsAlgorithm::Bfs, |component| {
            prefix.push(component.to_vec());
            if prefix.len() == stop + 1 {
                ControlFlow::Break(component.to_vec())
            } else {
                ControlFlow::Continue(())
            }
        });
        prop_assert_eq!(result, expected.get(stop).cloned().map_or(
            ControlFlow::Continue(()), ControlFlow::Break));
        prop_assert_eq!(prefix, expected[..(stop + 1).min(expected.len())].to_vec());
    });
}
