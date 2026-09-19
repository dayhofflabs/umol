//! Storage restoration laws through the public API.
//!
//! [`Graph::restore`] is compared with an independently enumerated endpoint/adjacency model
//! over multigraphs with up to 16 nodes and 40 edges. A separate producer roundtrip checks
//! compatibility with cascading removal, including retained shared graphs. Generated malformed
//! inputs independently vary compactions and saved entries to check freedom from panics.
//!
//! [`FixedRelationSet::restore`] is compared with independent rows and incidence scans over
//! up to 32 three-participant rows, including reordered saved entries and changed survivor
//! payloads. Separate roundtrips exercise the removal producer. [`FixedRelationSet::restore_participants`]
//! uses structured participants with node, edge, both, or neither reference and extra labels;
//! independently enumerated survivors define inverse translation over eight-id domains.
//! Separate compaction roundtrips check producer compatibility. Malformed-input properties
//! require only freedom from panics for both operations.
//!
//! VarRelationSet exercises the same restoration laws over heterogeneous rows of zero to eight
//! participants, checking packed row boundaries as well as complete values and incidence.
//!
//! FixedFixedBirelationSet uses factors of arity two and three. Independent expectations
//! retain their boundaries and scan union incidence, including references shared across factors.

use proptest::prelude::*;
use umol_graph_core::{
    Compaction, EdgeId, FixedFixedBirelationSet, FixedRelationSet, Graph, GraphCompaction,
    GraphCorrespondence, GraphRemapping, Neighbor, NodeId, ParticipantRefs, RelationId,
    RelationParticipant, VarRelationSet,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct References {
    node: Option<NodeId>,
    edge: Option<EdgeId>,
    label: u8,
}

impl RelationParticipant for References {
    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self> {
        Some(Self {
            node: match self.node {
                Some(id) => Some(id.try_map(correspondence)?),
                None => None,
            },
            edge: match self.edge {
                Some(id) => Some(id.try_map(correspondence)?),
                None => None,
            },
            ..self
        })
    }

    fn remap(self, remapping: &GraphRemapping) -> Self {
        Self {
            node: self.node.map(|id| id.remap(remapping)),
            edge: self.edge.map(|id| id.remap(remapping)),
            ..self
        }
    }

    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        Some(Self {
            node: match self.node {
                Some(id) => Some(id.compact(compaction)?),
                None => None,
            },
            edge: match self.edge {
                Some(id) => Some(id.compact(compaction)?),
                None => None,
            },
            ..self
        })
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
        Self {
            node: self.node.map(|id| id.uncompact(compaction)),
            edge: self.edge.map(|id| id.uncompact(compaction)),
            ..self
        }
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: self.node,
            edge: self.edge,
        }
    }
}

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

    #[test]
    fn test_fixed_relation_set_restore(
        rows in prop::collection::vec((prop::array::uniform3(0u32..16), any::<bool>()), 0..33),
        order in prop::collection::vec(any::<u32>(), 32),
    ) {
        let ids: Vec<_> = rows.iter().enumerate().filter(|(_, (_, remove))| *remove)
            .map(|(id, _)| RelationId::from(id)).collect();
        let compaction = Compaction::new(rows.len(), ids).unwrap();
        let expected: Vec<_> = rows.iter().enumerate()
            .map(|(id, (row, remove))| (row.map(NodeId), if *remove { id } else { id + rows.len() })).collect();
        let survivors = expected.iter().zip(&rows).filter(|(_, (_, remove))| !remove)
            .map(|(entry, _)| *entry).collect();
        let mut saved: Vec<_> = expected.iter().zip(&rows).enumerate().rev()
            .filter(|(_, (_, (_, remove)))| *remove)
            .map(|(id, ((row, data), _))| (RelationId::from(id), *row, *data)).collect();
        saved.sort_unstable_by_key(|(id, _, _)| (order[id.index()], *id));
        let mut relations = FixedRelationSet::new(survivors);
        relations.restore(&compaction, saved);
        for node in 0..17 {
            let incidence: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (row, _))| row.contains(&NodeId(node)))
                .map(|(id, _)| RelationId::from(id)).collect();
            prop_assert_eq!(relations.incident_to_node(NodeId(node)), incidence);
            prop_assert_eq!(relations.incident_to_edge(EdgeId(node)), &[]);
        }
        prop_assert_eq!(relations.into_entries(), expected);
    }

    #[test]
    fn test_fixed_relation_set_restore_roundtrip(
        rows in prop::collection::vec((prop::array::uniform3(0u32..16), any::<bool>()), 0..33),
    ) {
        let original = FixedRelationSet::new(rows.iter().enumerate()
            .map(|(id, (row, _))| (row.map(EdgeId), id)).collect());
        let ids: Vec<_> = rows.iter().enumerate().filter(|(_, (_, remove))| *remove)
            .map(|(id, _)| RelationId::from(id)).collect();
        let saved = rows.iter().enumerate().rev().filter(|(_, (_, remove))| *remove)
            .map(|(id, (row, _))| (RelationId::from(id), row.map(EdgeId), id)).collect();
        let mut relations = original.clone();
        let compaction = relations.tracked_remove(&ids);
        relations.restore(&compaction, saved);
        for edge in 0..17 {
            let incidence: Vec<_> = rows.iter().enumerate()
                .filter(|(_, (row, _))| row.contains(&edge))
                .map(|(id, _)| RelationId::from(id)).collect();
            prop_assert_eq!(relations.incident_to_edge(EdgeId(edge)), incidence);
        }
        prop_assert_eq!(relations, original);
    }

    #[test]
    fn test_fixed_relation_set_restore_malformed(
        rows in prop::collection::vec(prop::array::uniform2(any::<u32>()), 0..33),
        count in 0usize..41,
        removals in prop::collection::vec(any::<bool>(), 40),
        saved in prop::collection::vec((0u32..48, prop::array::uniform2(any::<u32>())), 0..41),
    ) {
        let mut relations = FixedRelationSet::new(rows.into_iter().enumerate()
            .map(|(id, row)| (row.map(NodeId), id)).collect());
        let compaction = Compaction::new(count, (0..count).filter(|&id| removals[id])
            .map(RelationId::from).collect()).unwrap();
        relations.restore(&compaction, saved.into_iter().enumerate()
            .map(|(data, (id, row))| (RelationId(id), row.map(NodeId), data)).collect());
    }

    #[test]
    fn test_fixed_relation_set_restore_participants(
        nodes in 0usize..9,
        edges in 0usize..9,
        node_removals in prop::collection::vec(any::<bool>(), 8),
        edge_removals in prop::collection::vec(any::<bool>(), 8),
        rows in prop::collection::vec(prop::array::uniform2((0usize..8, 0usize..8, any::<bool>(), any::<bool>(), any::<u8>())), 0..33),
    ) {
        let surviving_nodes: Vec<_> = (0..nodes).filter(|&id| !node_removals[id]).map(NodeId::from).collect();
        let surviving_edges: Vec<_> = (0..edges).filter(|&id| !edge_removals[id]).map(EdgeId::from).collect();
        let compaction = GraphCompaction::new(
            Compaction::new(nodes, (0..nodes).filter(|&id| node_removals[id]).map(NodeId::from).collect()).unwrap(),
            Compaction::new(edges, (0..edges).filter(|&id| edge_removals[id]).map(EdgeId::from).collect()).unwrap(),
        );
        let entries: Vec<_> = rows.iter().enumerate().map(|(id, row)| (row.map(|(n, e, has_n, has_e, label)| References {
            node: (has_n && !surviving_nodes.is_empty()).then(|| NodeId::from(n % surviving_nodes.len())),
            edge: (has_e && !surviving_edges.is_empty()).then(|| EdgeId::from(e % surviving_edges.len())),
            label,
        }), id)).collect();
        let expected: Vec<_> = entries.iter().map(|(row, data)| (row.map(|p| References {
            node: p.node.map(|id| surviving_nodes[id.index()]),
            edge: p.edge.map(|id| surviving_edges[id.index()]),
            ..p
        }), *data)).collect();
        let mut relations = FixedRelationSet::new(entries);
        relations.restore_participants(&compaction);
        for id in 0..9 {
            let node_incidence: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (row, _))| row.iter().any(|p| p.node == Some(NodeId(id))))
                .map(|(id, _)| RelationId::from(id)).collect();
            let edge_incidence: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (row, _))| row.iter().any(|p| p.edge == Some(EdgeId(id))))
                .map(|(id, _)| RelationId::from(id)).collect();
            prop_assert_eq!(relations.incident_to_node(NodeId(id)), node_incidence);
            prop_assert_eq!(relations.incident_to_edge(EdgeId(id)), edge_incidence);
        }
        prop_assert_eq!(relations.into_entries(), expected);
    }

    #[test]
    fn test_fixed_relation_set_restore_participants_roundtrip(
        rows in prop::collection::vec(prop::array::uniform2((0u32..8, 0u32..8, any::<bool>(), any::<bool>(), any::<u8>())), 0..33),
        node_removals in prop::collection::vec(any::<bool>(), 8),
        edge_removals in prop::collection::vec(any::<bool>(), 8),
    ) {
        let entries: Vec<_> = rows.into_iter().enumerate().map(|(id, row)| (row.map(|(n, e, has_n, has_e, label)| References {
            node: has_n.then_some(NodeId(n)), edge: has_e.then_some(EdgeId(e)), label,
        }), id)).collect();
        let expected: Vec<_> = entries.iter().filter(|(row, _)| row.iter().all(|p|
            p.node.is_none_or(|id| !node_removals[id.index()]) && p.edge.is_none_or(|id| !edge_removals[id.index()])
        )).copied().collect();
        let compaction = GraphCompaction::new(
            Compaction::new(8, (0..8).filter(|&id| node_removals[id]).map(NodeId::from).collect()).unwrap(),
            Compaction::new(8, (0..8).filter(|&id| edge_removals[id]).map(EdgeId::from).collect()).unwrap(),
        );
        let mut relations = FixedRelationSet::new(entries).compact(&compaction);
        relations.restore_participants(&compaction);
        for id in 0..9 {
            let nodes: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (row, _))| row.iter().any(|p| p.node == Some(NodeId(id))))
                .map(|(id, _)| RelationId::from(id)).collect();
            let edges: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (row, _))| row.iter().any(|p| p.edge == Some(EdgeId(id))))
                .map(|(id, _)| RelationId::from(id)).collect();
            prop_assert_eq!(relations.incident_to_node(NodeId(id)), nodes);
            prop_assert_eq!(relations.incident_to_edge(EdgeId(id)), edges);
        }
        prop_assert_eq!(relations.into_entries(), expected);
    }

    #[test]
    fn test_fixed_relation_set_restore_participants_malformed(
        rows in prop::collection::vec(prop::array::uniform2((prop::option::of(0u32..16), prop::option::of(0u32..16), any::<u8>())), 0..33),
        nodes in 0usize..9,
        edges in 0usize..9,
        node_removals in prop::collection::vec(any::<bool>(), 8),
        edge_removals in prop::collection::vec(any::<bool>(), 8),
    ) {
        let mut relations = FixedRelationSet::new(rows.into_iter().enumerate().map(|(id, row)|
            (row.map(|(node, edge, label)| References { node: node.map(NodeId), edge: edge.map(EdgeId), label }), id)
        ).collect());
        relations.restore_participants(&GraphCompaction::new(
            Compaction::new(nodes, (0..nodes).filter(|&id| node_removals[id]).map(NodeId::from).collect()).unwrap(),
            Compaction::new(edges, (0..edges).filter(|&id| edge_removals[id]).map(EdgeId::from).collect()).unwrap(),
        ));
    }

    #[test]
    fn test_var_relation_set_restore(
        rows in prop::collection::vec((prop::collection::vec(0u32..16, 0..9), any::<bool>()), 0..33),
        order in prop::collection::vec(any::<u32>(), 32),
    ) {
        let ids: Vec<_> = rows.iter().enumerate().filter(|(_, (_, remove))| *remove)
            .map(|(id, _)| RelationId::from(id)).collect();
        let compaction = Compaction::new(rows.len(), ids).unwrap();
        let expected: Vec<_> = rows.iter().enumerate()
            .map(|(id, (row, remove))| (row.iter().copied().map(NodeId).collect::<Vec<_>>(), if *remove { id } else { id + rows.len() })).collect();
        let survivors = expected.iter().zip(&rows).filter(|(_, (_, remove))| !remove)
            .map(|(entry, _)| entry.clone()).collect();
        let mut saved: Vec<_> = expected.iter().zip(&rows).enumerate().rev()
            .filter(|(_, (_, (_, remove)))| *remove)
            .map(|(id, ((row, data), _))| (RelationId::from(id), row.clone(), *data)).collect();
        saved.sort_unstable_by_key(|(id, _, _)| (order[id.index()], *id));
        let mut relations = VarRelationSet::new(survivors);
        relations.restore(&compaction, saved);
        for node in 0..17 {
            let incidence: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (row, _))| row.contains(&NodeId(node)))
                .map(|(id, _)| RelationId::from(id)).collect();
            prop_assert_eq!(relations.incident_to_node(NodeId(node)), incidence);
            prop_assert_eq!(relations.incident_to_edge(EdgeId(node)), &[]);
        }
        prop_assert_eq!(relations.into_entries(), expected);
    }

    #[test]
    fn test_var_relation_set_restore_roundtrip(
        rows in prop::collection::vec((prop::collection::vec(0u32..16, 0..9), any::<bool>()), 0..33),
    ) {
        let original = VarRelationSet::new(rows.iter().enumerate()
            .map(|(id, (row, _))| (row.iter().copied().map(EdgeId).collect::<Vec<_>>(), id)).collect());
        let ids: Vec<_> = rows.iter().enumerate().filter(|(_, (_, remove))| *remove)
            .map(|(id, _)| RelationId::from(id)).collect();
        let saved = rows.iter().enumerate().rev().filter(|(_, (_, remove))| *remove)
            .map(|(id, (row, _))| (RelationId::from(id), row.iter().copied().map(EdgeId).collect::<Vec<_>>(), id)).collect();
        let mut relations = original.clone();
        let compaction = relations.tracked_remove(&ids);
        relations.restore(&compaction, saved);
        for edge in 0..17 {
            let incidence: Vec<_> = rows.iter().enumerate()
                .filter(|(_, (row, _))| row.contains(&edge))
                .map(|(id, _)| RelationId::from(id)).collect();
            prop_assert_eq!(relations.incident_to_edge(EdgeId(edge)), incidence);
        }
        prop_assert_eq!(relations, original);
    }

    #[test]
    fn test_var_relation_set_restore_malformed(
        rows in prop::collection::vec(prop::collection::vec(any::<u32>(), 0..9), 0..33),
        count in 0usize..41,
        removals in prop::collection::vec(any::<bool>(), 40),
        saved in prop::collection::vec((0u32..48, prop::collection::vec(any::<u32>(), 0..9)), 0..41),
    ) {
        let mut relations = VarRelationSet::new(rows.into_iter().enumerate()
            .map(|(id, row)| (row.iter().copied().map(NodeId).collect::<Vec<_>>(), id)).collect());
        let compaction = Compaction::new(count, (0..count).filter(|&id| removals[id])
            .map(RelationId::from).collect()).unwrap();
        relations.restore(&compaction, saved.into_iter().enumerate()
            .map(|(data, (id, row))| (RelationId(id), row.iter().copied().map(NodeId).collect::<Vec<_>>(), data)).collect());
    }

    #[test]
    fn test_var_relation_set_restore_participants(
        nodes in 0usize..9,
        edges in 0usize..9,
        node_removals in prop::collection::vec(any::<bool>(), 8),
        edge_removals in prop::collection::vec(any::<bool>(), 8),
        rows in prop::collection::vec(prop::collection::vec((0usize..8, 0usize..8, any::<bool>(), any::<bool>(), any::<u8>()), 0..9), 0..33),
    ) {
        let surviving_nodes: Vec<_> = (0..nodes).filter(|&id| !node_removals[id]).map(NodeId::from).collect();
        let surviving_edges: Vec<_> = (0..edges).filter(|&id| !edge_removals[id]).map(EdgeId::from).collect();
        let compaction = GraphCompaction::new(
            Compaction::new(nodes, (0..nodes).filter(|&id| node_removals[id]).map(NodeId::from).collect()).unwrap(),
            Compaction::new(edges, (0..edges).filter(|&id| edge_removals[id]).map(EdgeId::from).collect()).unwrap(),
        );
        let entries: Vec<_> = rows.iter().enumerate().map(|(id, row)| (row.iter().copied().map(|(n, e, has_n, has_e, label)| References {
            node: (has_n && !surviving_nodes.is_empty()).then(|| NodeId::from(n % surviving_nodes.len())),
            edge: (has_e && !surviving_edges.is_empty()).then(|| EdgeId::from(e % surviving_edges.len())),
            label,
        }).collect::<Vec<_>>(), id)).collect();
        let expected: Vec<_> = entries.iter().map(|(row, data)| (row.iter().copied().map(|p| References {
            node: p.node.map(|id| surviving_nodes[id.index()]),
            edge: p.edge.map(|id| surviving_edges[id.index()]),
            ..p
        }).collect::<Vec<_>>(), *data)).collect();
        let mut relations = VarRelationSet::new(entries);
        relations.restore_participants(&compaction);
        for id in 0..9 {
            let node_incidence: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (row, _))| row.iter().any(|p| p.node == Some(NodeId(id))))
                .map(|(id, _)| RelationId::from(id)).collect();
            let edge_incidence: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (row, _))| row.iter().any(|p| p.edge == Some(EdgeId(id))))
                .map(|(id, _)| RelationId::from(id)).collect();
            prop_assert_eq!(relations.incident_to_node(NodeId(id)), node_incidence);
            prop_assert_eq!(relations.incident_to_edge(EdgeId(id)), edge_incidence);
        }
        prop_assert_eq!(relations.into_entries(), expected);
    }

    #[test]
    fn test_var_relation_set_restore_participants_roundtrip(
        rows in prop::collection::vec(prop::collection::vec((0u32..8, 0u32..8, any::<bool>(), any::<bool>(), any::<u8>()), 0..9), 0..33),
        node_removals in prop::collection::vec(any::<bool>(), 8),
        edge_removals in prop::collection::vec(any::<bool>(), 8),
    ) {
        let entries: Vec<_> = rows.into_iter().enumerate().map(|(id, row)| (row.iter().copied().map(|(n, e, has_n, has_e, label)| References {
            node: has_n.then_some(NodeId(n)), edge: has_e.then_some(EdgeId(e)), label,
        }).collect::<Vec<_>>(), id)).collect();
        let expected: Vec<_> = entries.iter().filter(|(row, _)| row.iter().all(|p|
            p.node.is_none_or(|id| !node_removals[id.index()]) && p.edge.is_none_or(|id| !edge_removals[id.index()])
        )).cloned().collect();
        let compaction = GraphCompaction::new(
            Compaction::new(8, (0..8).filter(|&id| node_removals[id]).map(NodeId::from).collect()).unwrap(),
            Compaction::new(8, (0..8).filter(|&id| edge_removals[id]).map(EdgeId::from).collect()).unwrap(),
        );
        let mut relations = VarRelationSet::new(entries).compact(&compaction);
        relations.restore_participants(&compaction);
        for id in 0..9 {
            let nodes: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (row, _))| row.iter().any(|p| p.node == Some(NodeId(id))))
                .map(|(id, _)| RelationId::from(id)).collect();
            let edges: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (row, _))| row.iter().any(|p| p.edge == Some(EdgeId(id))))
                .map(|(id, _)| RelationId::from(id)).collect();
            prop_assert_eq!(relations.incident_to_node(NodeId(id)), nodes);
            prop_assert_eq!(relations.incident_to_edge(EdgeId(id)), edges);
        }
        prop_assert_eq!(relations.into_entries(), expected);
    }

    #[test]
    fn test_var_relation_set_restore_participants_malformed(
        rows in prop::collection::vec(prop::collection::vec((prop::option::of(0u32..16), prop::option::of(0u32..16), any::<u8>()), 0..9), 0..33),
        nodes in 0usize..9,
        edges in 0usize..9,
        node_removals in prop::collection::vec(any::<bool>(), 8),
        edge_removals in prop::collection::vec(any::<bool>(), 8),
    ) {
        let mut relations = VarRelationSet::new(rows.into_iter().enumerate().map(|(id, row)|
            (row.iter().copied().map(|(node, edge, label)| References { node: node.map(NodeId), edge: edge.map(EdgeId), label }).collect::<Vec<_>>(), id)
        ).collect());
        relations.restore_participants(&GraphCompaction::new(
            Compaction::new(nodes, (0..nodes).filter(|&id| node_removals[id]).map(NodeId::from).collect()).unwrap(),
            Compaction::new(edges, (0..edges).filter(|&id| edge_removals[id]).map(EdgeId::from).collect()).unwrap(),
        ));
    }

    #[test]
    fn test_fixed_fixed_birelation_set_restore(
        rows in prop::collection::vec((prop::array::uniform5(0u32..16), any::<bool>()), 0..33),
        order in prop::collection::vec(any::<u32>(), 32),
    ) {
        let ids = rows.iter().enumerate().filter(|(_, (_, remove))| *remove)
            .map(|(id, _)| RelationId::from(id)).collect();
        let compaction = Compaction::new(rows.len(), ids).unwrap();
        let expected: Vec<_> = rows.iter().enumerate().map(|(id, (row, remove))| (
            [NodeId(row[0]), NodeId(row[1])], [NodeId(row[2]), NodeId(row[3]), NodeId(row[4])],
            if *remove { id } else { id + rows.len() })).collect();
        let survivors = expected.iter().zip(&rows).filter(|(_, (_, remove))| !remove).map(|(entry, _)| *entry).collect();
        let mut saved: Vec<_> = expected.iter().zip(&rows).enumerate().rev().filter(|(_, (_, (_, remove)))| *remove)
            .map(|(id, ((first, second, data), _))| (RelationId::from(id), *first, *second, *data)).collect();
        saved.sort_unstable_by_key(|(id, _, _, _)| (order[id.index()], *id));
        let mut relations = FixedFixedBirelationSet::new(survivors);
        relations.restore(&compaction, saved);
        for node in 0..17 {
            let incidence: Vec<_> = expected.iter().enumerate()
                .filter(|(_, (first, second, _))| first.contains(&NodeId(node)) || second.contains(&NodeId(node)))
                .map(|(id, _)| RelationId::from(id)).collect();
            prop_assert_eq!(relations.incident_to_node(NodeId(node)), incidence);
            prop_assert_eq!(relations.incident_to_edge(EdgeId(node)), &[]);
        }
        prop_assert_eq!(relations.into_entries(), expected);
    }

    #[test]
    fn test_fixed_fixed_birelation_set_restore_roundtrip(
        rows in prop::collection::vec((prop::array::uniform5(0u32..16), any::<bool>()), 0..33),
    ) {
        let expected: Vec<_> = rows.iter().enumerate().map(|(id, (row, _))| (
            [NodeId(row[0]), NodeId(row[1])], [EdgeId(row[2]), EdgeId(row[3]), EdgeId(row[4])], id)).collect();
        let original = FixedFixedBirelationSet::new(expected.clone());
        let ids: Vec<_> = rows.iter().enumerate().filter(|(_, (_, remove))| *remove).map(|(id, _)| RelationId::from(id)).collect();
        let saved = expected.iter().zip(&rows).enumerate().rev().filter(|(_, (_, (_, remove)))| *remove)
            .map(|(id, ((first, second, data), _))| (RelationId::from(id), *first, *second, *data)).collect();
        let mut relations = original.clone();
        let compaction = relations.tracked_remove(&ids);
        relations.restore(&compaction, saved);
        for id in 0..17 {
            let nodes: Vec<_> = expected.iter().enumerate().filter(|(_, (row, _, _))| row.contains(&NodeId(id))).map(|(i, _)| RelationId::from(i)).collect();
            let edges: Vec<_> = expected.iter().enumerate().filter(|(_, (_, row, _))| row.contains(&EdgeId(id))).map(|(i, _)| RelationId::from(i)).collect();
            prop_assert_eq!(relations.incident_to_node(NodeId(id)), nodes);
            prop_assert_eq!(relations.incident_to_edge(EdgeId(id)), edges);
        }
        prop_assert_eq!(relations, original);
    }

    #[test]
    fn test_fixed_fixed_birelation_set_restore_malformed(
        rows in prop::collection::vec(prop::array::uniform5(any::<u32>()), 0..33),
        count in 0usize..41,
        removals in prop::collection::vec(any::<bool>(), 40),
        saved in prop::collection::vec((0u32..48, prop::array::uniform5(any::<u32>())), 0..41),
    ) {
        let mut relations = FixedFixedBirelationSet::new(rows.into_iter().enumerate().map(|(id, row)| (
            [NodeId(row[0]), NodeId(row[1])], [EdgeId(row[2]), EdgeId(row[3]), EdgeId(row[4])], id)).collect());
        let compaction = Compaction::new(count, (0..count).filter(|&id| removals[id]).map(RelationId::from).collect()).unwrap();
        relations.restore(&compaction, saved.into_iter().enumerate().map(|(data, (id, row))| (
            RelationId(id), [NodeId(row[0]), NodeId(row[1])], [EdgeId(row[2]), EdgeId(row[3]), EdgeId(row[4])], data)).collect());
    }

    #[test]
    fn test_fixed_fixed_birelation_set_restore_participants(
        nodes in 0usize..9,
        edges in 0usize..9,
        node_removals in prop::collection::vec(any::<bool>(), 8),
        edge_removals in prop::collection::vec(any::<bool>(), 8),
        rows in prop::collection::vec(prop::array::uniform5((0usize..8, 0usize..8, any::<bool>(), any::<bool>(), any::<u8>())), 0..33),
    ) {
        let surviving_nodes: Vec<_> = (0..nodes).filter(|&id| !node_removals[id]).map(NodeId::from).collect();
        let surviving_edges: Vec<_> = (0..edges).filter(|&id| !edge_removals[id]).map(EdgeId::from).collect();
        let compaction = GraphCompaction::new(
            Compaction::new(nodes, (0..nodes).filter(|&id| node_removals[id]).map(NodeId::from).collect()).unwrap(),
            Compaction::new(edges, (0..edges).filter(|&id| edge_removals[id]).map(EdgeId::from).collect()).unwrap(),
        );
        let entries: Vec<_> = rows.iter().enumerate().map(|(id, row)| {
            let row = row.map(|(n, e, has_n, has_e, label)| References {
                node: (has_n && !surviving_nodes.is_empty()).then(|| NodeId::from(n % surviving_nodes.len())),
                edge: (has_e && !surviving_edges.is_empty()).then(|| EdgeId::from(e % surviving_edges.len())), label,
            });
            ([row[0], row[1]], [row[2], row[3], row[4]], id)
        }).collect();
        let expected: Vec<_> = entries.iter().map(|(first, second, data)| (
            first.map(|p| References { node: p.node.map(|id| surviving_nodes[id.index()]), edge: p.edge.map(|id| surviving_edges[id.index()]), ..p }),
            second.map(|p| References { node: p.node.map(|id| surviving_nodes[id.index()]), edge: p.edge.map(|id| surviving_edges[id.index()]), ..p }), *data)).collect();
        let mut relations = FixedFixedBirelationSet::new(entries);
        relations.restore_participants(&compaction);
        for id in 0..9 {
            let nodes: Vec<_> = expected.iter().enumerate().filter(|(_, (first, second, _))|
                first.iter().chain(second).any(|p| p.node == Some(NodeId(id)))).map(|(i, _)| RelationId::from(i)).collect();
            let edges: Vec<_> = expected.iter().enumerate().filter(|(_, (first, second, _))|
                first.iter().chain(second).any(|p| p.edge == Some(EdgeId(id)))).map(|(i, _)| RelationId::from(i)).collect();
            prop_assert_eq!(relations.incident_to_node(NodeId(id)), nodes);
            prop_assert_eq!(relations.incident_to_edge(EdgeId(id)), edges);
        }
        prop_assert_eq!(relations.into_entries(), expected);
    }

    #[test]
    fn test_fixed_fixed_birelation_set_restore_participants_roundtrip(
        rows in prop::collection::vec(prop::array::uniform5((0u32..8, 0u32..8, any::<bool>(), any::<bool>(), any::<u8>())), 0..33),
        node_removals in prop::collection::vec(any::<bool>(), 8),
        edge_removals in prop::collection::vec(any::<bool>(), 8),
    ) {
        let entries: Vec<_> = rows.into_iter().enumerate().map(|(id, row)| {
            let row = row.map(|(n, e, has_n, has_e, label)| References { node: has_n.then_some(NodeId(n)), edge: has_e.then_some(EdgeId(e)), label });
            ([row[0], row[1]], [row[2], row[3], row[4]], id)
        }).collect();
        let expected: Vec<_> = entries.iter().filter(|(first, second, _)| first.iter().chain(second).all(|p|
            p.node.is_none_or(|id| !node_removals[id.index()]) && p.edge.is_none_or(|id| !edge_removals[id.index()]))).copied().collect();
        let compaction = GraphCompaction::new(
            Compaction::new(8, (0..8).filter(|&id| node_removals[id]).map(NodeId::from).collect()).unwrap(),
            Compaction::new(8, (0..8).filter(|&id| edge_removals[id]).map(EdgeId::from).collect()).unwrap(),
        );
        let mut relations = FixedFixedBirelationSet::new(entries).compact(&compaction);
        relations.restore_participants(&compaction);
        for id in 0..9 {
            let nodes: Vec<_> = expected.iter().enumerate().filter(|(_, (first, second, _))|
                first.iter().chain(second).any(|p| p.node == Some(NodeId(id)))).map(|(i, _)| RelationId::from(i)).collect();
            let edges: Vec<_> = expected.iter().enumerate().filter(|(_, (first, second, _))|
                first.iter().chain(second).any(|p| p.edge == Some(EdgeId(id)))).map(|(i, _)| RelationId::from(i)).collect();
            prop_assert_eq!(relations.incident_to_node(NodeId(id)), nodes);
            prop_assert_eq!(relations.incident_to_edge(EdgeId(id)), edges);
        }
        prop_assert_eq!(relations.into_entries(), expected);
    }

    #[test]
    fn test_fixed_fixed_birelation_set_restore_participants_malformed(
        rows in prop::collection::vec(prop::array::uniform5((prop::option::of(0u32..16), prop::option::of(0u32..16), any::<u8>())), 0..33),
        nodes in 0usize..9,
        edges in 0usize..9,
        node_removals in prop::collection::vec(any::<bool>(), 8),
        edge_removals in prop::collection::vec(any::<bool>(), 8),
    ) {
        let mut relations = FixedFixedBirelationSet::new(rows.into_iter().enumerate().map(|(id, row)| {
            let row = row.map(|(node, edge, label)| References { node: node.map(NodeId), edge: edge.map(EdgeId), label });
            ([row[0], row[1]], [row[2], row[3], row[4]], id)
        }).collect());
        relations.restore_participants(&GraphCompaction::new(
            Compaction::new(nodes, (0..nodes).filter(|&id| node_removals[id]).map(NodeId::from).collect()).unwrap(),
            Compaction::new(edges, (0..edges).filter(|&id| edge_removals[id]).map(EdgeId::from).collect()).unwrap(),
        ));
    }
}
