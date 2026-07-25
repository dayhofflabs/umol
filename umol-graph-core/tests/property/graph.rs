//! Properties of the core graph representation and its predicates.

use std::collections::HashSet;

use proptest::prelude::*;

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
}
