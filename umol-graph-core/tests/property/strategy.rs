//! Generated graph inputs used by property tests.

use std::collections::HashSet;

use proptest::prelude::*;
use umol_graph_core::Graph;

pub(super) fn simple_graph(max_nodes: usize, max_edges: usize) -> impl Strategy<Value = Graph> {
    graph_with_simple_edges(max_nodes, max_edges).prop_map(|(graph, _)| graph)
}

pub(super) fn graph_with_simple_edges(
    max_nodes: usize,
    max_edges: usize,
) -> impl Strategy<Value = (Graph, Vec<[u32; 2]>)> {
    (
        0..=max_nodes,
        prop::collection::vec((0..max_nodes as u32, 0..max_nodes as u32), 0..=max_edges),
    )
        .prop_map(|(node_count, endpoints)| {
            let mut seen = HashSet::new();
            let edges: Vec<[u32; 2]> = endpoints
                .into_iter()
                .filter(|&(first, second)| {
                    first < node_count as u32 && second < node_count as u32 && first != second
                })
                .filter(|&(first, second)| seen.insert((first.min(second), first.max(second))))
                .map(|(first, second)| [first, second])
                .collect();
            (Graph::new(node_count, &edges), edges)
        })
}
