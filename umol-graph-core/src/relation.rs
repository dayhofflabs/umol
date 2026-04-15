//! Relation sets: N-ary relations over graph nodes with CSR incidence.
//!
//! `FixedRelationSet<R, N>` stores relations of compile-time-known arity
//! (e.g. binary dative bonds). `VarRelationSet<R>` stores variable-arity
//! relations (e.g. aromatic systems). Both use sorted parallel arrays for
//! incidence lookup, avoiding per-node offset tables.

use crate::graph::NodeId;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RelationId(pub u32);

impl RelationId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

fn build_incidence<'a>(
    relation_count: usize,
    participants_of: impl Fn(usize) -> &'a [NodeId],
) -> (Vec<NodeId>, Vec<RelationId>) {
    let mut entries: Vec<(NodeId, RelationId)> = Vec::new();
    for i in 0..relation_count {
        let rid = RelationId(i as u32);
        for &node in participants_of(i) {
            entries.push((node, rid));
        }
    }
    entries.sort_by_key(|&(node, _)| node);

    let nodes = entries.iter().map(|&(n, _)| n).collect();
    let rels = entries.iter().map(|&(_, r)| r).collect();
    (nodes, rels)
}

/// Fixed-arity relation set. Each relation connects exactly N nodes.
///
/// Flat CSR storage: participants are `Vec<[NodeId; N]>`, incidence is
/// a flat array with offset table. No heap allocations per node or
/// per relation.
#[derive(Clone, Debug)]
pub struct FixedRelationSet<R, const N: usize> {
    participants: Vec<[NodeId; N]>,
    data: Vec<R>,
    incidence_nodes: Vec<NodeId>,
    incidence_rels: Vec<RelationId>,
}

impl<R: PartialEq, const N: usize> PartialEq for FixedRelationSet<R, N> {
    fn eq(&self, other: &Self) -> bool {
        self.participants == other.participants && self.data == other.data
    }
}

impl<R: Eq, const N: usize> Eq for FixedRelationSet<R, N> {}

impl<R, const N: usize> FixedRelationSet<R, N> {
    pub fn new(entries: Vec<([NodeId; N], R)>) -> Self {
        let mut participants = Vec::with_capacity(entries.len());
        let mut data = Vec::with_capacity(entries.len());
        for (p, d) in entries {
            participants.push(p);
            data.push(d);
        }

        let (incidence_nodes, incidence_rels) =
            build_incidence(participants.len(), |i| &participants[i]);

        Self {
            participants,
            data,
            incidence_nodes,
            incidence_rels,
        }
    }

    pub fn relation_count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &R {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut R {
        &mut self.data[id.index()]
    }

    pub fn participants(&self, id: RelationId) -> &[NodeId; N] {
        &self.participants[id.index()]
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        let start = self.incidence_nodes.partition_point(|n| *n < node);
        let end = start + self.incidence_nodes[start..].partition_point(|n| *n <= node);
        &self.incidence_rels[start..end]
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence_nodes.binary_search(&node).is_ok()
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

}

impl<R, const N: usize> Default for FixedRelationSet<R, N> {
    fn default() -> Self {
        Self {
            participants: Vec::new(),
            data: Vec::new(),
            incidence_nodes: Vec::new(),
            incidence_rels: Vec::new(),
        }
    }
}

/// Variable-arity relation set. Each relation connects an arbitrary
/// number of nodes.
///
/// Flat CSR storage: participant ranges via offset table, incidence
/// via a second offset table. No heap allocations per node or per
/// relation.
#[derive(Clone, Debug)]
pub struct VarRelationSet<R> {
    offsets: Vec<u32>,
    participants: Vec<NodeId>,
    data: Vec<R>,
    incidence_nodes: Vec<NodeId>,
    incidence_rels: Vec<RelationId>,
}

impl<R: PartialEq> PartialEq for VarRelationSet<R> {
    fn eq(&self, other: &Self) -> bool {
        self.offsets == other.offsets
            && self.participants == other.participants
            && self.data == other.data
    }
}

impl<R: Eq> Eq for VarRelationSet<R> {}

impl<R> VarRelationSet<R> {
    pub fn new(entries: Vec<(Vec<NodeId>, R)>) -> Self {
        let relation_count = entries.len();
        let mut offsets = Vec::with_capacity(relation_count + 1);
        offsets.push(0);

        let total_participants: usize = entries.iter().map(|(p, _)| p.len()).sum();
        let mut participants = Vec::with_capacity(total_participants);
        let mut data = Vec::with_capacity(relation_count);

        for (p, d) in entries {
            participants.extend_from_slice(&p);
            offsets.push(participants.len() as u32);
            data.push(d);
        }

        let (incidence_nodes, incidence_rels) =
            build_incidence(relation_count, |i| {
                let start = offsets[i] as usize;
                let end = offsets[i + 1] as usize;
                &participants[start..end]
            });

        Self {
            offsets,
            participants,
            data,
            incidence_nodes,
            incidence_rels,
        }
    }

    pub fn relation_count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &R {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut R {
        &mut self.data[id.index()]
    }

    pub fn participants(&self, id: RelationId) -> &[NodeId] {
        let start = self.offsets[id.index()] as usize;
        let end = self.offsets[id.index() + 1] as usize;
        &self.participants[start..end]
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        let start = self.incidence_nodes.partition_point(|n| *n < node);
        let end = start + self.incidence_nodes[start..].partition_point(|n| *n <= node);
        &self.incidence_rels[start..end]
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence_nodes.binary_search(&node).is_ok()
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

}

impl<R> Default for VarRelationSet<R> {
    fn default() -> Self {
        Self {
            offsets: vec![0],
            participants: Vec::new(),
            data: Vec::new(),
            incidence_nodes: Vec::new(),
            incidence_rels: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use super::*;

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    #[test]
    fn test_fixed_relation_set_new() {
        let rs: FixedRelationSet<&str, 2> = FixedRelationSet::new(
            vec![([n(0), n(1)], "dative"), ([n(1), n(2)], "noncov")],
        );
        assert_eq!(rs.relation_count(), 2);
        assert_eq!(rs.data(RelationId(0)), &"dative");
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(2)]);
    }

    #[test]
    fn test_fixed_relation_set_incidence() {
        let rs: FixedRelationSet<(), 2> = FixedRelationSet::new(vec![
            ([n(0), n(1)], ()),
            ([n(0), n(2)], ()),
            ([n(2), n(3)], ()),
        ]);
        assert_eq!(rs.incident(n(0)), &[RelationId(0), RelationId(1)]);
        assert_eq!(rs.incident(n(1)), &[RelationId(0)]);
        assert_eq!(rs.incident(n(2)), &[RelationId(1), RelationId(2)]);
        assert_eq!(rs.incident(n(3)), &[RelationId(2)]);
        assert!(rs.has_incident(n(0)));
        assert!(!rs.has_incident(n(5)));
    }

    #[test]
    fn test_fixed_relation_set_data_mut() {
        let mut rs: FixedRelationSet<i32, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[test]
    fn test_fixed_relation_set_default() {
        let rs = FixedRelationSet::<(), 2>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[test]
    fn test_fixed_relation_set_relation_ids() {
        let rs: FixedRelationSet<(), 2> = FixedRelationSet::new(
            vec![([n(0), n(1)], ()), ([n(1), n(2)], ())],
        );
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[test]
    fn test_var_relation_set_new() {
        let rs: VarRelationSet<&str> = VarRelationSet::new(vec![
            (vec![n(0), n(1), n(2), n(3), n(4), n(5)], "benzene"),
        ]);
        assert_eq!(rs.relation_count(), 1);
        assert_eq!(rs.data(RelationId(0)), &"benzene");
        assert_eq!(
            rs.participants(RelationId(0)),
            &[n(0), n(1), n(2), n(3), n(4), n(5)]
        );
    }

    #[test]
    fn test_var_relation_set_incidence() {
        let rs: VarRelationSet<()> = VarRelationSet::new(vec![
            (vec![n(0), n(1), n(2)], ()),
            (vec![n(2), n(3), n(4)], ()),
        ]);
        assert_eq!(rs.incident(n(0)), &[RelationId(0)]);
        assert_eq!(rs.incident(n(2)), &[RelationId(0), RelationId(1)]);
        assert_eq!(rs.incident(n(4)), &[RelationId(1)]);
        assert!(rs.has_incident(n(0)));
        assert!(!rs.has_incident(n(7)));
    }

    #[test]
    fn test_var_relation_set_variable_arity() {
        let rs: VarRelationSet<&str> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], "pair"),
            (vec![n(2), n(3), n(4), n(5)], "quad"),
        ]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(2), n(3), n(4), n(5)]);
    }

    #[test]
    fn test_var_relation_set_data_mut() {
        let mut rs: VarRelationSet<i32> =
            VarRelationSet::new(vec![(vec![n(0), n(1), n(2)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[test]
    fn test_var_relation_set_default() {
        let rs = VarRelationSet::<()>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[test]
    fn test_var_relation_set_relation_ids() {
        let rs: VarRelationSet<()> = VarRelationSet::new(
            vec![(vec![n(0), n(1)], ()), (vec![n(1), n(2)], ())],
        );
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }
}
