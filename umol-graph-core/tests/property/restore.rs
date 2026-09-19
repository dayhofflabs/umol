//! Storage restoration laws through the public API.
//!
//! [`Graph::restore`] is compared with an independently enumerated endpoint/adjacency model
//! over multigraphs with up to 16 nodes and 40 edges. A separate producer roundtrip checks
//! compatibility with cascading removal, including retained shared graphs. Generated malformed
//! inputs independently vary compactions and saved entries to check freedom from panics.

use proptest::prelude::*;
use umol_graph_core::{Compaction, EdgeId, Graph, GraphCompaction, Neighbor, NodeId};

proptest! {
    #[test]
    fn test_graph_restore(
        node_count in 0usize..17,
        raw_edges in prop::collection::vec((0u32..16, 0u32..16, any::<bool>()), 0..41),
        node_removals in prop::collection::vec(any::<bool>(), 16),
        rotation in any::<usize>(),
    ) {
        let edges: Vec<_> = raw_edges.into_iter().filter(|_| node_count > 0)
            .map(|(a, b, remove)| ([a % node_count as u32, b % node_count as u32], remove))
            .collect();
        let surviving_nodes: Vec<_> = (0..node_count)
            .filter(|&id| !node_removals[id]).map(NodeId::from).collect();
        let mut removed = Vec::new();
        let mut compacted_edges = Vec::new();
        for (id, &(nodes, requested)) in edges.iter().enumerate() {
            if requested || nodes.iter().any(|node| !surviving_nodes.contains(&NodeId(*node))) {
                removed.push((EdgeId::from(id), nodes.map(NodeId)));
            } else {
                compacted_edges.push(nodes.map(|node| {
                    surviving_nodes.iter().position(|id| *id == NodeId(node)).unwrap() as u32
                }));
            }
        }
        let compaction = GraphCompaction::new(
            Compaction::new(node_count, (0..node_count).filter(|&id| node_removals[id])
                .map(NodeId::from).collect()).unwrap(),
            Compaction::new(edges.len(), removed.iter().map(|&(id, _)| id).collect()).unwrap(),
        );
        removed.reverse();
        if !removed.is_empty() {
            let offset = rotation % removed.len();
            removed.rotate_left(offset);
        }
        let mut graph = Graph::new(surviving_nodes.len(), &compacted_edges);
        graph.restore(&compaction, &removed);
        prop_assert_eq!(graph.node_count(), node_count);
        prop_assert_eq!(graph.edge_count(), edges.len());
        let mut expected_neighbors = vec![vec![]; node_count];
        for (id, &([a, b], _)) in edges.iter().enumerate() {
            let edge = EdgeId::from(id);
            prop_assert_eq!(graph.edge_endpoints(edge), [NodeId(a.min(b)), NodeId(a.max(b))]);
            expected_neighbors[a as usize].push(Neighbor { node: NodeId(b), edge });
            expected_neighbors[b as usize].push(Neighbor { node: NodeId(a), edge });
        }
        for (id, mut expected) in expected_neighbors.into_iter().enumerate() {
            expected.sort_unstable_by_key(|neighbor| (neighbor.node, neighbor.edge));
            let mut actual = graph.neighbors(NodeId::from(id)).to_vec();
            actual.sort_unstable_by_key(|neighbor| (neighbor.node, neighbor.edge));
            prop_assert_eq!(actual, expected);
        }
    }

    #[test]
    fn test_graph_restore_roundtrip(
        node_count in 0usize..17,
        raw_edges in prop::collection::vec((0u32..16, 0u32..16, any::<bool>()), 0..41),
        node_removals in prop::collection::vec(any::<bool>(), 16),
    ) {
        let edges: Vec<_> = raw_edges.into_iter().filter(|_| node_count > 0)
            .map(|(a, b, remove)| ([a % node_count as u32, b % node_count as u32], remove))
            .collect();
        let endpoints: Vec<_> = edges.iter().map(|&(nodes, _)| nodes).collect();
        let original = Graph::new(node_count, &endpoints);
        let nodes: Vec<_> = (0..node_count).filter(|&id| node_removals[id])
            .map(NodeId::from).collect();
        let removed_edges: Vec<_> = edges.iter().enumerate().filter(|(_, (_, remove))| *remove)
            .map(|(id, _)| EdgeId::from(id)).collect();
        let removed: Vec<_> = edges.iter().enumerate().rev()
            .filter(|(_, (pair, requested))| *requested || pair.iter().any(|&n| nodes.contains(&NodeId(n))))
            .map(|(id, &(pair, _))| (EdgeId::from(id), pair.map(NodeId))).collect();
        let mut graph = original.clone();
        let compaction = graph.tracked_remove_cascading(&nodes, &removed_edges);
        let before = Graph::new(graph.node_count(), &graph.edge_ids()
            .map(|edge| graph.edge_endpoints(edge).map(|node| node.0)).collect::<Vec<_>>());
        let retained = graph.clone();
        graph.restore(&compaction, &removed);
        prop_assert_eq!(graph, original);
        prop_assert_eq!(retained, before);
    }

    #[test]
    fn test_graph_restore_malformed(
        node_count in 0usize..17,
        raw_edges in prop::collection::vec((0u32..16, 0u32..16), 0..41),
        source_nodes in 0usize..25,
        source_edges in 0usize..49,
        node_removals in prop::collection::vec(any::<bool>(), 24),
        edge_removals in prop::collection::vec(any::<bool>(), 48),
        saved in prop::collection::vec((0u32..56, 0u32..32, 0u32..32), 0..57),
        shared in any::<bool>(),
    ) {
        let endpoints: Vec<_> = raw_edges.into_iter().filter(|_| node_count > 0)
            .map(|(a, b)| [a % node_count as u32, b % node_count as u32]).collect();
        let mut graph = Graph::new(node_count, &endpoints);
        let retained = shared.then(|| graph.clone());
        let compaction = GraphCompaction::new(
            Compaction::new(source_nodes, (0..source_nodes).filter(|&i| node_removals[i])
                .map(NodeId::from).collect()).unwrap(),
            Compaction::new(source_edges, (0..source_edges).filter(|&i| edge_removals[i])
                .map(EdgeId::from).collect()).unwrap(),
        );
        let removed: Vec<_> = saved.into_iter()
            .map(|(id, a, b)| (EdgeId(id), [NodeId(a), NodeId(b)])).collect();
        graph.restore(&compaction, &removed);
        drop(retained);
    }
}
