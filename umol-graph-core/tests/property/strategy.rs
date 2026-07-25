//! Generated graph inputs used by property tests.

use proptest::prelude::*;
use umol_graph_core::Graph;

pub(super) fn multigraph(max_nodes: usize, max_edges: usize) -> impl Strategy<Value = Graph> {
    graph_with_edge_multiset(max_nodes, max_edges).prop_map(|(graph, _)| graph)
}

pub(super) fn graph_with_edge_multiset(
    max_nodes: usize,
    max_edges: usize,
) -> impl Strategy<Value = (Graph, Vec<[u32; 2]>)> {
    (
        0..=max_nodes,
        prop::collection::vec((0..max_nodes as u32, 0..max_nodes as u32), 0..=max_edges),
    )
        .prop_map(|(node_count, endpoints)| {
            let edges: Vec<[u32; 2]> = endpoints
                .into_iter()
                .filter(|&(first, second)| first < node_count as u32 && second < node_count as u32)
                .map(|(first, second)| [first, second])
                .collect();
            (Graph::new(node_count, &edges), edges)
        })
}
