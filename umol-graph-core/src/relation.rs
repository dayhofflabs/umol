//! Relation sets: N-ary relations over graph nodes with CSR incidence.
//!
//! `FixedRelationSet<P, O, D, N>` stores relations of compile-time-known arity
//! (e.g. binary dative bonds); `VarRelationSet<P, O, D>` stores variable-arity
//! relations (e.g. aromatic systems). Participants are typed `P`
//! (`RelationParticipant`); the factor ordering `O` (`Unordered`/`Ordered`)
//! controls canonicalization. Both use sorted parallel arrays for incidence.

use std::hash::Hash;
use std::marker::PhantomData;

use crate::graph::{EdgeId, NodeId, Remapping};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RelationId(pub u32);

impl RelationId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Canonicalization of a relation factor's participants, applied on
/// construction and after a remap relabels them. `Unordered` sorts (membership
/// is the datum); `Ordered` preserves input order (position is the datum, e.g.
/// a stereo configuration or a bond direction).
pub trait FactorOrdering {
    fn canonicalize<P: Ord>(participants: &mut [P]);
}

/// Set-valued factor: participants are canonicalized by sorting.
#[derive(Clone, Copy, Debug)]
pub struct Unordered;

/// Positional factor: canonicalization is a no-op, input order is preserved.
#[derive(Clone, Copy, Debug)]
pub struct Ordered;

impl FactorOrdering for Unordered {
    fn canonicalize<P: Ord>(participants: &mut [P]) {
        participants.sort_unstable();
    }
}

impl FactorOrdering for Ordered {
    fn canonicalize<P: Ord>(_participants: &mut [P]) {}
}

/// The id-space contents of a participant, surfaced for the incidence index.
/// At most one ref per space today (a node or an edge); a future port type
/// could fill both.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ParticipantRefs {
    pub node: Option<NodeId>,
    pub edge: Option<EdgeId>,
}

/// A value that can occupy a relation factor: routes through a `Remapping` in
/// both directions and exposes its node/edge refs for incidence. One impl per
/// concrete id type — dispatch is static, since a factor is homogeneous.
pub trait RelationParticipant: Copy + Ord + Hash {
    fn remap(self, remapping: &Remapping) -> Option<Self>;
    fn unmap(self, remapping: &Remapping) -> Self;
    fn refs(self) -> ParticipantRefs;
}

impl RelationParticipant for NodeId {
    fn remap(self, remapping: &Remapping) -> Option<Self> {
        remapping.map_node(self)
    }

    fn unmap(self, remapping: &Remapping) -> Self {
        remapping.unmap_node(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: Some(self),
            edge: None,
        }
    }
}

impl RelationParticipant for EdgeId {
    fn remap(self, remapping: &Remapping) -> Option<Self> {
        remapping.map_edge(self)
    }

    fn unmap(self, remapping: &Remapping) -> Self {
        remapping.unmap_edge(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: None,
            edge: Some(self),
        }
    }
}

fn build_incidence<'a, P: RelationParticipant + 'a>(
    relation_count: usize,
    participants_of: impl Fn(usize) -> &'a [P],
) -> (Vec<NodeId>, Vec<RelationId>) {
    let mut entries: Vec<(NodeId, RelationId)> = Vec::new();
    for i in 0..relation_count {
        let rid = RelationId(i as u32);
        for &p in participants_of(i) {
            if let Some(node) = p.refs().node {
                entries.push((node, rid));
            }
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
pub struct FixedRelationSet<P, O, D, const N: usize> {
    participants: Vec<[P; N]>,
    data: Vec<D>,
    incidence_nodes: Vec<NodeId>,
    incidence_rels: Vec<RelationId>,
    _ordering: PhantomData<O>,
}

impl<P: PartialEq, O, D: PartialEq, const N: usize> PartialEq for FixedRelationSet<P, O, D, N> {
    fn eq(&self, other: &Self) -> bool {
        self.participants == other.participants && self.data == other.data
    }
}

impl<P: Eq, O, D: Eq, const N: usize> Eq for FixedRelationSet<P, O, D, N> {}

impl<P: RelationParticipant, O: FactorOrdering, D, const N: usize> FixedRelationSet<P, O, D, N> {
    pub fn new(entries: Vec<([P; N], D)>) -> Self {
        let mut participants = Vec::with_capacity(entries.len());
        let mut data = Vec::with_capacity(entries.len());
        for (mut p, d) in entries {
            O::canonicalize(&mut p);
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
            _ordering: PhantomData,
        }
    }

    pub fn relation_count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    pub fn data_iter_mut(&mut self) -> impl Iterator<Item = &mut D> {
        self.data.iter_mut()
    }

    pub fn participants(&self, id: RelationId) -> &[P; N] {
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

    pub fn apply_remapping(&self, remapping: &Remapping) -> Self
    where
        D: Clone,
    {
        let entries: Vec<([P; N], D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let parts: Option<Vec<P>> = self
                    .participants(rid)
                    .iter()
                    .map(|&p| p.remap(remapping))
                    .collect();
                let parts: [P; N] = parts?.try_into().ok()?;
                Some((parts, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }
}

impl<P, O, D, const N: usize> Default for FixedRelationSet<P, O, D, N> {
    fn default() -> Self {
        Self {
            participants: Vec::new(),
            data: Vec::new(),
            incidence_nodes: Vec::new(),
            incidence_rels: Vec::new(),
            _ordering: PhantomData,
        }
    }
}

/// Variable-arity relation set. Each relation connects an arbitrary
/// number of nodes. Participants are sorted by `NodeId` on construction.
///
/// Flat CSR storage: participant ranges via offset table, incidence
/// via a second offset table. No heap allocations per node or per
/// relation.
#[derive(Clone, Debug)]
pub struct VarRelationSet<P, O, D> {
    offsets: Vec<u32>,
    participants: Vec<P>,
    data: Vec<D>,
    incidence_nodes: Vec<NodeId>,
    incidence_rels: Vec<RelationId>,
    _ordering: PhantomData<O>,
}

impl<P: PartialEq, O, D: PartialEq> PartialEq for VarRelationSet<P, O, D> {
    fn eq(&self, other: &Self) -> bool {
        self.offsets == other.offsets
            && self.participants == other.participants
            && self.data == other.data
    }
}

impl<P: Eq, O, D: Eq> Eq for VarRelationSet<P, O, D> {}

impl<P: RelationParticipant, O: FactorOrdering, D> VarRelationSet<P, O, D> {
    pub fn new(entries: Vec<(Vec<P>, D)>) -> Self {
        let relation_count = entries.len();
        let mut offsets = Vec::with_capacity(relation_count + 1);
        offsets.push(0);

        let total_participants: usize = entries.iter().map(|(p, _)| p.len()).sum();
        let mut participants = Vec::with_capacity(total_participants);
        let mut data = Vec::with_capacity(relation_count);

        for (mut p, d) in entries {
            O::canonicalize(&mut p);
            participants.extend_from_slice(&p);
            offsets.push(participants.len() as u32);
            data.push(d);
        }

        let (incidence_nodes, incidence_rels) = build_incidence(relation_count, |i| {
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
            _ordering: PhantomData,
        }
    }

    pub fn relation_count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    pub fn data_iter_mut(&mut self) -> impl Iterator<Item = &mut D> {
        self.data.iter_mut()
    }

    pub fn participants(&self, id: RelationId) -> &[P] {
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

    pub fn apply_remapping(&self, remapping: &Remapping) -> Self
    where
        D: Clone,
    {
        let entries: Vec<(Vec<P>, D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let parts: Option<Vec<P>> = self
                    .participants(rid)
                    .iter()
                    .map(|&p| p.remap(remapping))
                    .collect();
                Some((parts?, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }
}

impl<P, O, D> Default for VarRelationSet<P, O, D> {
    fn default() -> Self {
        Self {
            offsets: vec![0],
            participants: Vec::new(),
            data: Vec::new(),
            incidence_nodes: Vec::new(),
            incidence_rels: Vec::new(),
            _ordering: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    #[test]
    fn test_fixed_relation_set_new() {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], "dative"), ([n(1), n(2)], "noncov")]);
        assert_eq!(rs.relation_count(), 2);
        assert_eq!(rs.data(RelationId(0)), &"dative");
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(2)]);
    }

    #[test]
    fn test_fixed_relation_set_participants_sorted() {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(2), n(0)], "a"), ([n(3), n(1)], "b")]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(2)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(3)]);
    }

    #[test]
    fn test_fixed_relation_set_participants_ordered() {
        let rs: FixedRelationSet<NodeId, Ordered, &str, 2> =
            FixedRelationSet::new(vec![([n(2), n(0)], "a"), ([n(3), n(1)], "b")]);
        assert_eq!(rs.participants(RelationId(0)), &[n(2), n(0)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(3), n(1)]);
    }

    #[test]
    fn test_fixed_relation_set_apply_remapping() {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(0), n(2)], "keep"), ([n(1), n(3)], "drop")]);
        let remapping = Remapping {
            removed_nodes: vec![1],
            removed_edges: vec![],
        };
        let out = rs.apply_remapping(&remapping);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[test]
    fn test_fixed_relation_set_incidence() {
        let rs: FixedRelationSet<NodeId, Unordered, (), 2> = FixedRelationSet::new(vec![
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
        let mut rs: FixedRelationSet<NodeId, Unordered, i32, 2> = FixedRelationSet::new(vec![([n(0), n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[test]
    fn test_fixed_relation_set_data_iter_mut() {
        let mut rs: FixedRelationSet<NodeId, Unordered, i32, 2> = FixedRelationSet::new(vec![
            ([n(0), n(1)], 1),
            ([n(1), n(2)], 2),
            ([n(2), n(3)], 3),
        ]);
        for d in rs.data_iter_mut() {
            *d *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
        assert_eq!(rs.data(RelationId(2)), &30);
    }

    #[test]
    fn test_fixed_relation_set_default() {
        let rs = FixedRelationSet::<NodeId, Unordered, (), 2>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[test]
    fn test_fixed_relation_set_relation_ids() {
        let rs: FixedRelationSet<NodeId, Unordered, (), 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], ()), ([n(1), n(2)], ())]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[test]
    fn test_var_relation_set_new() {
        let rs: VarRelationSet<NodeId, Unordered, &str> =
            VarRelationSet::new(vec![(vec![n(0), n(1), n(2), n(3), n(4), n(5)], "benzene")]);
        assert_eq!(rs.relation_count(), 1);
        assert_eq!(rs.data(RelationId(0)), &"benzene");
        assert_eq!(
            rs.participants(RelationId(0)),
            &[n(0), n(1), n(2), n(3), n(4), n(5)]
        );
    }

    #[test]
    fn test_var_relation_set_participants_sorted() {
        let rs: VarRelationSet<NodeId, Unordered, ()> = VarRelationSet::new(vec![
            (vec![n(5), n(2), n(0), n(3)], ()),
            (vec![n(4), n(1)], ()),
        ]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(2), n(3), n(5)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(4)]);
    }

    #[test]
    fn test_var_relation_set_participants_ordered() {
        let rs: VarRelationSet<NodeId, Ordered, ()> = VarRelationSet::new(vec![
            (vec![n(5), n(2), n(0), n(3)], ()),
            (vec![n(4), n(1)], ()),
        ]);
        assert_eq!(rs.participants(RelationId(0)), &[n(5), n(2), n(0), n(3)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(4), n(1)]);
        assert_eq!(rs.incident(n(0)), &[RelationId(0)]);
        assert_eq!(rs.incident(n(4)), &[RelationId(1)]);
    }

    #[test]
    fn test_var_relation_set_apply_remapping() {
        let rs: VarRelationSet<NodeId, Unordered, &str> = VarRelationSet::new(vec![
            (vec![n(0), n(2), n(4)], "keep"),
            (vec![n(1), n(3)], "drop"),
        ]);
        let remapping = Remapping {
            removed_nodes: vec![1],
            removed_edges: vec![],
        };
        let out = rs.apply_remapping(&remapping);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1), n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[test]
    fn test_var_relation_set_incidence() {
        let rs: VarRelationSet<NodeId, Unordered, ()> = VarRelationSet::new(vec![
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
        let rs: VarRelationSet<NodeId, Unordered, &str> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], "pair"),
            (vec![n(2), n(3), n(4), n(5)], "quad"),
        ]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(2), n(3), n(4), n(5)]);
    }

    #[test]
    fn test_var_relation_set_data_mut() {
        let mut rs: VarRelationSet<NodeId, Unordered, i32> = VarRelationSet::new(vec![(vec![n(0), n(1), n(2)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[test]
    fn test_var_relation_set_data_iter_mut() {
        let mut rs: VarRelationSet<NodeId, Unordered, i32> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], 1),
            (vec![n(2), n(3), n(4)], 2),
            (vec![n(5)], 3),
        ]);
        for d in rs.data_iter_mut() {
            *d *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
        assert_eq!(rs.data(RelationId(2)), &30);
    }

    #[test]
    fn test_var_relation_set_default() {
        let rs = VarRelationSet::<NodeId, Unordered, ()>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[test]
    fn test_var_relation_set_relation_ids() {
        let rs: VarRelationSet<NodeId, Unordered, ()> =
            VarRelationSet::new(vec![(vec![n(0), n(1)], ()), (vec![n(1), n(2)], ())]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    #[case::already_sorted(vec![0, 1, 2], vec![0, 1, 2])]
    #[case::reversed(vec![2, 1, 0], vec![0, 1, 2])]
    #[case::shuffled(vec![2, 0, 3, 1], vec![0, 1, 2, 3])]
    fn test_unordered_canonicalize(#[case] mut input: Vec<i32>, #[case] expected: Vec<i32>) {
        Unordered::canonicalize(&mut input);
        assert_eq!(input, expected);
    }

    #[rstest]
    #[case::already_sorted(vec![0, 1, 2])]
    #[case::reversed(vec![2, 1, 0])]
    #[case::shuffled(vec![2, 0, 3, 1])]
    fn test_ordered_canonicalize(#[case] input: Vec<i32>) {
        let mut actual = input.clone();
        Ordered::canonicalize(&mut actual);
        assert_eq!(actual, input);
    }

    #[rstest]
    #[case::before_removed(NodeId(0), Some(NodeId(0)))]
    #[case::removed(NodeId(1), None)]
    #[case::after_removed(NodeId(2), Some(NodeId(1)))]
    fn test_node_id_remap(#[case] id: NodeId, #[case] expected: Option<NodeId>) {
        let remapping = Remapping {
            removed_nodes: vec![1],
            removed_edges: vec![],
        };
        assert_eq!(id.remap(&remapping), expected);
    }

    #[rstest]
    #[case::before_gap(NodeId(0), NodeId(0))]
    #[case::after_gap(NodeId(1), NodeId(2))]
    fn test_node_id_unmap(#[case] id: NodeId, #[case] expected: NodeId) {
        let remapping = Remapping {
            removed_nodes: vec![1],
            removed_edges: vec![],
        };
        assert_eq!(id.unmap(&remapping), expected);
    }

    #[test]
    fn test_node_id_refs() {
        assert_eq!(
            NodeId(3).refs(),
            ParticipantRefs {
                node: Some(NodeId(3)),
                edge: None,
            }
        );
    }

    #[rstest]
    #[case::removed(EdgeId(0), None)]
    #[case::after_removed(EdgeId(2), Some(EdgeId(1)))]
    fn test_edge_id_remap(#[case] id: EdgeId, #[case] expected: Option<EdgeId>) {
        let remapping = Remapping {
            removed_nodes: vec![],
            removed_edges: vec![0],
        };
        assert_eq!(id.remap(&remapping), expected);
    }

    #[test]
    fn test_edge_id_refs() {
        assert_eq!(
            EdgeId(2).refs(),
            ParticipantRefs {
                node: None,
                edge: Some(EdgeId(2)),
            }
        );
    }
}
