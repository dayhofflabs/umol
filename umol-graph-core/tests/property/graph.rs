//! Properties of the core graph representation and its predicates.

use std::collections::HashSet;

use proptest::prelude::*;
use umol_graph_core::{Correspondence, EdgeId, Graph, GraphCorrespondence, Neighbor, NodeId};

use super::strategy::graph_with_edge_multiset;

proptest! {
    /// The selected edge IDs equal the endpoint-membership definition, including multigraphs.
    /// Reversing or repeating the node sequence preserves that edge set without repeating IDs.
    #[test]
    fn test_graph_induced_edges(
        (graph, edges) in graph_with_edge_multiset(16, 40),
        nodes in prop::collection::vec(0u32..16, 0..32),
    ) {
        let mut nodes: Vec<NodeId> = nodes.into_iter()
            .filter(|&node| (node as usize) < graph.node_count())
            .map(NodeId).collect();
        let expected: Vec<EdgeId> = edges.iter().enumerate()
            .filter(|(_, [a, b])| nodes.contains(&NodeId(*a)) && nodes.contains(&NodeId(*b)))
            .map(|(index, _)| EdgeId::from(index)).collect();

        let mut actual: Vec<_> = graph.induced_edges(&nodes).collect();
        actual.sort_unstable();
        prop_assert_eq!(&actual, &expected);

        nodes.reverse();
        nodes.extend_from_within(..);
        let mut actual: Vec<_> = graph.induced_edges(&nodes).collect();
        actual.sort_unstable();
        prop_assert_eq!(actual, expected);
    }

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

    /// Bulk node addition preserves the original topology and appends isolated nodes.
    #[test]
    fn test_graph_extend_nodes(
        (original, edges) in graph_with_edge_multiset(16, 40),
        count in 0usize..16,
        shared in any::<bool>(),
    ) {
        let mut graph = if shared { original.clone() } else { Graph::new(original.node_count(), &edges) };
        let ids = graph.extend_nodes(count);
        prop_assert_eq!(ids.len(), count);
        prop_assert_eq!(ids.collect::<Vec<_>>(), (original.node_count()..original.node_count() + count).map(NodeId::from).collect::<Vec<_>>());
        prop_assert_eq!(graph, Graph::new(original.node_count() + count, &edges));
        prop_assert_eq!(&original, &Graph::new(original.node_count(), &edges));
    }

    /// Bulk edge addition agrees with constructing the complete ordered edge list.
    #[test]
    fn test_graph_extend_edges(
        (original, mut edges) in graph_with_edge_multiset(16, 40),
        added in prop::collection::vec((0u32..16, 0u32..16), 0..32),
        shared in any::<bool>(),
    ) {
        let mut graph = if shared { original.clone() } else { Graph::new(original.node_count(), &edges) };
        let added: Vec<_> = added.into_iter()
            .filter(|&(a, b)| (a as usize) < original.node_count() && (b as usize) < original.node_count())
            .map(|(a, b)| [NodeId(a), NodeId(b)]).collect();
        let ids = graph.extend_edges(&added);
        prop_assert_eq!(ids.len(), added.len());
        prop_assert_eq!(ids.collect::<Vec<_>>(), (edges.len()..edges.len() + added.len()).map(EdgeId::from).collect::<Vec<_>>());
        prop_assert_eq!(&original, &Graph::new(original.node_count(), &edges));
        edges.extend(added.iter().map(|&[a, b]| [a.0, b.0]));
        prop_assert_eq!(graph, Graph::new(original.node_count(), &edges));
    }

    /// Combined extension preserves edge order and both incidences of each edge, including loops.
    /// Neighbors are derived directly from endpoints, independently of CSR construction.
    #[test]
    fn test_graph_extend(
        (original, mut edges) in graph_with_edge_multiset(16, 40),
        count in 0usize..16,
        added in prop::collection::vec((0u32..32, 0u32..32), 0..32),
        shared in any::<bool>(),
    ) {
        let mut graph = if shared { original.clone() } else { Graph::new(original.node_count(), &edges) };
        let node_count = original.node_count() + count;
        let added: Vec<_> = added.into_iter()
            .filter(|&(a, b)| (a as usize) < node_count && (b as usize) < node_count)
            .map(|(a, b)| [NodeId(a), NodeId(b)]).collect();
        let (nodes, edge_ids) = graph.extend(count, &added);
        prop_assert_eq!(nodes.len(), count);
        prop_assert_eq!(edge_ids.len(), added.len());
        prop_assert_eq!(nodes.collect::<Vec<_>>(), (original.node_count()..node_count).map(NodeId::from).collect::<Vec<_>>());
        prop_assert_eq!(edge_ids.collect::<Vec<_>>(), (edges.len()..edges.len() + added.len()).map(EdgeId::from).collect::<Vec<_>>());
        prop_assert_eq!(&original, &Graph::new(original.node_count(), &edges));
        edges.extend(added.iter().map(|&[a, b]| [a.0, b.0]));
        prop_assert_eq!(graph.node_count(), node_count);
        prop_assert_eq!(graph.edge_count(), edges.len());
        for (index, &[a, b]) in edges.iter().enumerate() {
            prop_assert_eq!(graph.edge_endpoints(EdgeId::from(index)), [NodeId(a.min(b)), NodeId(a.max(b))]);
        }
        for node in 0..node_count as u32 {
            let mut expected = Vec::new();
            for (index, &[a, b]) in edges.iter().enumerate() {
                if a == node { expected.push(Neighbor { node: NodeId(b), edge: EdgeId::from(index) }); }
                if b == node { expected.push(Neighbor { node: NodeId(a), edge: EdgeId::from(index) }); }
            }
            expected.sort_unstable_by_key(|neighbor| (neighbor.node, neighbor.edge));
            let mut actual = graph.neighbors(NodeId(node)).to_vec();
            prop_assert!(actual.windows(2).all(|pair| pair[0].node <= pair[1].node));
            actual.sort_unstable_by_key(|neighbor| (neighbor.node, neighbor.edge));
            prop_assert_eq!(actual, expected);
        }
    }

    #[test]
    fn test_graph_tracked_remove_cascading(
        (graph, edges) in graph_with_edge_multiset(8, 16),
        node_mask in any::<u8>(),
        edge_mask in any::<u16>(),
    ) {
        let removed_nodes = (0..graph.node_count())
            .filter(|&idx| node_mask & (1 << idx) != 0)
            .map(NodeId::from)
            .collect::<Vec<_>>();
        let removed_edges = (0..graph.edge_count())
            .filter(|&idx| edge_mask & (1 << idx) != 0)
            .map(EdgeId::from)
            .collect::<Vec<_>>();
        let surviving_nodes = (0..graph.node_count())
            .map(NodeId::from)
            .filter(|id| !removed_nodes.contains(id))
            .collect::<Vec<_>>();
        let surviving_edges = edges.iter().enumerate()
            .filter(|(idx, [a, b])| {
                !removed_edges.contains(&EdgeId::from(*idx))
                    && !removed_nodes.contains(&NodeId(*a))
                    && !removed_nodes.contains(&NodeId(*b))
            })
            .collect::<Vec<_>>();
        let expected_edges = surviving_edges.iter()
            .map(|(_, [a, b])| {
                [
                    surviving_nodes.iter().position(|id| id == &NodeId(*a)).unwrap() as u32,
                    surviving_nodes.iter().position(|id| id == &NodeId(*b)).unwrap() as u32,
                ]
            })
            .collect::<Vec<_>>();
        let expected = Graph::new(surviving_nodes.len(), &expected_edges);
        let expected_correspondence = GraphCorrespondence::new(
            Correspondence::new(
                surviving_nodes.iter().enumerate()
                    .map(|(idx, &old)| (old, NodeId::from(idx)))
                    .collect(),
                graph.node_count(),
                surviving_nodes.len(),
            ).unwrap(),
            Correspondence::new(
                surviving_edges.iter().enumerate()
                    .map(|(idx, (old, _))| (EdgeId::from(*old), EdgeId::from(idx)))
                    .collect(),
                graph.edge_count(),
                surviving_edges.len(),
            ).unwrap(),
        );

        let mut plain = graph.clone();
        let mut witnessed = graph.clone();
        plain.remove_cascading(&removed_nodes, &removed_edges);
        let compaction = witnessed.tracked_remove_cascading(&removed_nodes, &removed_edges);
        prop_assert_eq!(&plain, &expected);
        prop_assert_eq!(&witnessed, &expected);
        prop_assert_eq!(GraphCorrespondence::from(&compaction), expected_correspondence);

        let dangling = edges.iter().enumerate().any(|(idx, [a, b])| {
            (removed_nodes.contains(&NodeId(*a)) || removed_nodes.contains(&NodeId(*b)))
                && !removed_edges.contains(&EdgeId::from(idx))
        });
        let mut plain = graph.clone();
        let mut witnessed = graph.clone();
        prop_assert_eq!(plain.try_remove(&removed_nodes, &removed_edges), (!dangling).then_some(()));
        prop_assert_eq!(
            witnessed.try_tracked_remove(&removed_nodes, &removed_edges),
            (!dangling).then_some(compaction),
        );
        let expected = if dangling { graph } else { expected };
        prop_assert_eq!(plain, expected.clone());
        prop_assert_eq!(witnessed, expected);
    }
}
