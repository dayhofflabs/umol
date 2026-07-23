#[path = "cycles/oracle.rs"]
mod oracle;

use proptest::prelude::*;
use umol_graph_core::{Graph, SubdivisionNodeSource};

use self::oracle::{relevant_cycles, unique_ring_families};

/// Random multigraph with loops and repeated endpoint pairs.
fn graph(max_nodes: usize, max_edges: usize) -> impl Strategy<Value = Graph> {
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
            Graph::new(node_count, &edges)
        })
}

proptest! {
    #[test]
    fn test_unique_ring_families(graph in graph(5, 5)) {
        let relevant = relevant_cycles(&graph);
        let families = unique_ring_families(&graph);

        for family in &families {
            prop_assert!(!family.is_empty());
            let weight = family[0].len();
            prop_assert!(family.iter().all(|cycle| cycle.len() == weight));
        }

        let mut partition: Vec<_> = families.into_iter().flatten().collect();
        partition.sort_by(|left, right| {
            left.len()
                .cmp(&right.len())
                .then_with(|| left.cmp(right))
        });
        prop_assert_eq!(partition, relevant);
    }

    #[test]
    fn test_graph_subdivide_edges(source in graph(8, 12)) {
        let subdivision = source.subdivide_edges();
        let graph = subdivision.graph();

        prop_assert_eq!(
            graph.node_count(),
            source.node_count() + source.edge_count()
        );
        prop_assert_eq!(graph.edge_count(), 2 * source.edge_count());

        for node in graph.node_ids() {
            let node_source = subdivision.node_source(node);
            prop_assert_eq!(subdivision.node_of(node_source), node);
        }

        for incidence in graph.edge_ids() {
            let edge = subdivision.edge_source(incidence);
            prop_assert!(source.contains_edge(edge));
            prop_assert!(subdivision.incidence_edges_of(edge).contains(&incidence));
        }

        for edge in source.edge_ids() {
            let [first, second] = source.edge_endpoints(edge);
            let inserted = subdivision.node_of(SubdivisionNodeSource::Edge(edge));
            let [first_incidence, second_incidence] = subdivision.incidence_edges_of(edge);

            prop_assert_eq!(
                subdivision.node_source(inserted),
                SubdivisionNodeSource::Edge(edge)
            );
            prop_assert_eq!(
                graph.edge_endpoints(first_incidence),
                [first, inserted]
            );
            prop_assert_eq!(
                graph.edge_endpoints(second_incidence),
                [second, inserted]
            );
            prop_assert_eq!(subdivision.edge_source(first_incidence), edge);
            prop_assert_eq!(subdivision.edge_source(second_incidence), edge);
        }
    }
}
