//! Derived union incidence for relation storage.
//!
//! Node and edge keys occupy separate sorted arrays, paired with relation-id arrays.
//! Repeated references within a row contribute one entry per key; references from both
//! factors of a birelation contribute to the same index. Equal-key entries retain ascending
//! relation-id order, so anchored searches can select the first matching row.

use std::ops::Sub;

use super::participant::ParticipantRefs;
use crate::compact::{Compaction, GraphCompaction};
use crate::id::{EdgeId, NodeId, RelationId};

/// Union incidence index: a node → relations pair and an edge → relations pair,
/// each stored as parallel key and relation-id vectors. Participants self-route via `refs()`, so a
/// set with only node participants leaves the edge half empty, and vice versa.
#[derive(Clone, Debug, Default)]
pub(super) struct Incidence {
    node_keys: Vec<NodeId>,
    node_rels: Vec<RelationId>,
    edge_keys: Vec<EdgeId>,
    edge_rels: Vec<RelationId>,
}

impl Incidence {
    /// Build incidence for dense relation ids in `0..relation_count`.
    ///
    /// `fill(i, out)` pushes every participant's refs for row `i` into an empty buffer.
    /// Neither graph membership nor participant equality is checked. Callback panics propagate.
    pub(super) fn build(
        relation_count: usize,
        mut fill: impl FnMut(usize, &mut Vec<ParticipantRefs>),
    ) -> Self {
        let mut node_entries: Vec<(NodeId, RelationId)> = Vec::new();
        let mut edge_entries: Vec<(EdgeId, RelationId)> = Vec::new();
        let mut refs: Vec<ParticipantRefs> = Vec::new();
        let mut nodes: Vec<NodeId> = Vec::new();
        let mut edges: Vec<EdgeId> = Vec::new();
        for i in 0..relation_count {
            let rid = RelationId(i as u32);
            refs.clear();
            fill(i, &mut refs);
            nodes.clear();
            edges.clear();
            for r in &refs {
                if let Some(node) = r.node {
                    nodes.push(node);
                }
                if let Some(edge) = r.edge {
                    edges.push(edge);
                }
            }
            nodes.sort_unstable();
            nodes.dedup();
            edges.sort_unstable();
            edges.dedup();
            node_entries.extend(nodes.iter().map(|&n| (n, rid)));
            edge_entries.extend(edges.iter().map(|&e| (e, rid)));
        }
        node_entries.sort_by_key(|&(k, _)| k);
        edge_entries.sort_by_key(|&(k, _)| k);
        Self {
            node_keys: node_entries.iter().map(|&(k, _)| k).collect(),
            node_rels: node_entries.iter().map(|&(_, r)| r).collect(),
            edge_keys: edge_entries.iter().map(|&(k, _)| k).collect(),
            edge_rels: edge_entries.iter().map(|&(_, r)| r).collect(),
        }
    }

    /// Filter and renumber rows in place; monotone survivor ids retain equal-key order.
    pub(super) fn remove(&mut self, compaction: &Compaction<RelationId>) {
        remove_entries(
            &mut self.node_keys,
            &mut self.node_rels,
            compaction.removed(),
        );
        remove_entries(
            &mut self.edge_keys,
            &mut self.edge_rels,
            compaction.removed(),
        );
    }

    /// Translate sorted keys in place after successful participant compaction.
    /// All keys have surviving images; their order and relation ids are preserved.
    pub(super) fn compact_participants(&mut self, topology: &GraphCompaction) {
        compact_keys(&mut self.node_keys, topology.nodes().removed());
        compact_keys(&mut self.edge_keys, topology.edges().removed());
    }

    /// Rows referencing `node`, in ascending order; empty if absent.
    pub(super) fn incident_to_node(&self, node: NodeId) -> &[RelationId] {
        let start = self.node_keys.partition_point(|n| *n < node);
        let end = start + self.node_keys[start..].partition_point(|n| *n <= node);
        &self.node_rels[start..end]
    }

    /// Rows referencing `edge`, in ascending order; empty if absent.
    pub(super) fn incident_to_edge(&self, edge: EdgeId) -> &[RelationId] {
        let start = self.edge_keys.partition_point(|e| *e < edge);
        let end = start + self.edge_keys[start..].partition_point(|e| *e <= edge);
        &self.edge_rels[start..end]
    }

    /// Whether any row references `node`.
    pub(super) fn has_incident_to_node(&self, node: NodeId) -> bool {
        self.node_keys.binary_search(&node).is_ok()
    }

    /// Whether any row references `edge`.
    pub(super) fn has_incident_to_edge(&self, edge: EdgeId) -> bool {
        self.edge_keys.binary_search(&edge).is_ok()
    }
}

fn remove_entries<K: Copy>(
    keys: &mut Vec<K>,
    relations: &mut Vec<RelationId>,
    removed: &[RelationId],
) {
    let mut dst = 0;
    for src in 0..relations.len() {
        let relation = relations[src];
        if let Err(offset) = removed.binary_search(&relation) {
            keys[dst] = keys[src];
            relations[dst] = relation - offset;
            dst += 1;
        }
    }
    keys.truncate(dst);
    relations.truncate(dst);
}

fn compact_keys<K: Copy + Ord + Sub<usize, Output = K>>(keys: &mut [K], removed: &[K]) {
    if removed.is_empty() {
        return;
    }
    let mut offset = 0;
    for key in keys {
        while offset < removed.len() && removed[offset] < *key {
            offset += 1;
        }
        *key = *key - offset;
    }
}
