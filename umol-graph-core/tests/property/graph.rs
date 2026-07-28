//! Properties of the core graph representation and its predicates.

use std::collections::HashSet;

use proptest::prelude::*;
use umol_graph_core::{EdgeId, NodeId};

use super::strategy::graph_with_edge_multiset;

proptest! {
    #[test]
    fn test_graph_is_simple(
        (graph, edges) in graph_with_edge_multiset(8, 16),
    ) {
        let mut endpoint_pairs = HashSet::with_capacity(edges.len());
        let expected = edges.into_iter().all(|[first, second]| {
            let endpoints = if first <= second {
                [first, second]
            } else {
                [second, first]
            };
            first != second && endpoint_pairs.insert(endpoints)
        });

        prop_assert_eq!(graph.is_simple(), expected);
    }

    #[test]
    fn test_graph_id_iterators_exact_size(
        (graph, _) in graph_with_edge_multiset(8, 16),
        prefix in any::<usize>(),
    ) {
        let node_count = graph.node_count();
        let node_prefix = prefix.min(node_count);
        let mut nodes = graph.node_ids();
        prop_assert_eq!(nodes.len(), node_count);
        prop_assert_eq!(nodes.size_hint(), (node_count, Some(node_count)));
        for index in 0..node_prefix {
            prop_assert_eq!(nodes.next(), Some(NodeId(index as u32)));
            let remaining = node_count - index - 1;
            prop_assert_eq!(nodes.len(), remaining);
            prop_assert_eq!(nodes.size_hint(), (remaining, Some(remaining)));
        }
        prop_assert_eq!(
            nodes.collect::<Vec<_>>(),
            (node_prefix..node_count)
                .map(|index| NodeId(index as u32))
                .collect::<Vec<_>>(),
        );

        let edge_count = graph.edge_count();
        let edge_prefix = prefix.min(edge_count);
        let mut edges = graph.edge_ids();
        prop_assert_eq!(edges.len(), edge_count);
        prop_assert_eq!(edges.size_hint(), (edge_count, Some(edge_count)));
        for index in 0..edge_prefix {
            prop_assert_eq!(edges.next(), Some(EdgeId(index as u32)));
            let remaining = edge_count - index - 1;
            prop_assert_eq!(edges.len(), remaining);
            prop_assert_eq!(edges.size_hint(), (remaining, Some(remaining)));
        }
        prop_assert_eq!(
            edges.collect::<Vec<_>>(),
            (edge_prefix..edge_count)
                .map(|index| EdgeId(index as u32))
                .collect::<Vec<_>>(),
        );
    }
}
