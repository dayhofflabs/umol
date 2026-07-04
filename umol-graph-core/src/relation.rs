//! Relation sets: N-ary relations over typed participants (`NodeId`, `EdgeId`,
//! or external type implementing `RelationParticipant`), each carrying a shared
//! union incidence index (a node index and an edge index) routed from every
//! participant's `refs()`.
//! `FixedRelationSet<P, O, D, N>` stores relations of compile-time-known arity,
//! `VarRelationSet<P, O, D>` stores variable-arity relations. Participants are
//! typed `P` (`RelationParticipant`); the factor ordering `O` (`Unordered`/`Ordered`)
//! controls canonicalization. `FixedFixedBirelationSet`, `FixedVarBirelationSet`,
//! and `VarVarBirelationSet` relate two factors, each with its own participant
//! type, ordering, and arity. The union incidence spans both factors, so a relation
//! is reachable from any of its participants regardless of id-space.

use std::hash::Hash;
use std::marker::PhantomData;

use crate::correspondence::Correspondence;
use crate::graph::{Compaction, EdgeId, NodeId, Remapping};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RelationId(pub u32);

impl RelationId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<usize> for RelationId {
    fn from(index: usize) -> Self {
        Self(index as u32)
    }
}

/// The glued relation set and its two relation-level coprojections — the result of a same-space
/// relation-set [`pushout`](FixedRelationSet::pushout).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelationPushout<S> {
    pub object: S,
    /// `self` relation → object relation (identity — the object keeps `self`'s relation ids).
    pub left: Correspondence<RelationId>,
    /// `right` relation → object relation (a coincidence folds onto its `self` partner; the rest are
    /// appended after `self`).
    pub right: Correspondence<RelationId>,
}

/// The two coprojections of a relation pushout: `self` is the identity prefix `0..self_count`,
/// `right` follows `right_map`. Both over the object relation space of size `object_count`.
fn relation_pushout<S>(
    object: S,
    self_count: usize,
    object_count: usize,
    right_map: Vec<RelationId>,
) -> RelationPushout<S> {
    let left: Vec<RelationId> = (0..self_count).map(RelationId::from).collect();
    RelationPushout {
        object,
        left: Correspondence::from_images(&left, object_count),
        right: Correspondence::from_images(&right_map, object_count),
    }
}

/// Position of a participant within a single relation's tuple — local to the
/// relation (frame-relative), distinct from the global `NodeId`/`EdgeId`/`RelationId`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ParticipantPosition(pub u32);

impl ParticipantPosition {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Canonicalization of a relation factor's participants, applied on
/// construction and after a remap relabels them. `Unordered` sorts (membership
/// is the datum); `Ordered` preserves input order (position is the datum).
pub trait FactorOrdering {
    fn canonicalize<P: Ord>(participants: &mut [P]);
    /// Canonicalize `participants` in place and return the position permutation
    /// `σ`: the new position `i` holds the participant from old position `σ[i]`.
    fn canonicalize_positions<P: Ord + Copy>(participants: &mut [P]) -> Vec<ParticipantPosition>;
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

    fn canonicalize_positions<P: Ord + Copy>(participants: &mut [P]) -> Vec<ParticipantPosition> {
        let mut order: Vec<usize> = (0..participants.len()).collect();
        order.sort_by(|&a, &b| participants[a].cmp(&participants[b]));
        let sorted: Vec<P> = order.iter().map(|&i| participants[i]).collect();
        participants.copy_from_slice(&sorted);
        order
            .into_iter()
            .map(|i| ParticipantPosition(i as u32))
            .collect()
    }
}

impl FactorOrdering for Ordered {
    fn canonicalize<P: Ord>(_participants: &mut [P]) {}

    fn canonicalize_positions<P: Ord + Copy>(participants: &mut [P]) -> Vec<ParticipantPosition> {
        (0..participants.len() as u32)
            .map(ParticipantPosition)
            .collect()
    }
}

/// The id-space contents of a participant, surfaced for the incidence index.
/// At most one ref per space today (a node or an edge); a future port type
/// could fill both.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ParticipantRefs {
    pub node: Option<NodeId>,
    pub edge: Option<EdgeId>,
}

/// A single node or edge to route a participant through the incidence index — the resolved,
/// exactly-one form of `ParticipantRefs`, used to narrow `find_by_participants` candidates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ParticipantAnchor {
    Node(NodeId),
    Edge(EdgeId),
}

/// A value that can occupy a relation factor: routes through a `Compaction`
/// (removal/compaction, both directions) and a `Remapping` (general relabel,
/// forward), and exposes its node/edge refs for incidence. One impl per concrete
/// id type — dispatch is static, since a factor is homogeneous.
pub trait RelationParticipant: Copy + Ord + Hash {
    fn compact(self, compaction: &Compaction) -> Option<Self>;
    fn uncompact(self, compaction: &Compaction) -> Self;
    fn remap(self, remapping: &Remapping) -> Self;
    fn refs(self) -> ParticipantRefs;

    /// The node or edge to route this participant through the incidence index, if any — narrows
    /// `find_by_participants` candidates (`None` falls back to a linear scan).
    fn anchor(self) -> Option<ParticipantAnchor>;
}

impl RelationParticipant for NodeId {
    fn compact(self, compaction: &Compaction) -> Option<Self> {
        compaction.compact_node(self)
    }

    fn uncompact(self, compaction: &Compaction) -> Self {
        compaction.uncompact_node(self)
    }

    fn remap(self, remapping: &Remapping) -> Self {
        remapping.map_node(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: Some(self),
            edge: None,
        }
    }

    fn anchor(self) -> Option<ParticipantAnchor> {
        Some(ParticipantAnchor::Node(self))
    }
}

impl RelationParticipant for EdgeId {
    fn compact(self, compaction: &Compaction) -> Option<Self> {
        compaction.compact_edge(self)
    }

    fn uncompact(self, compaction: &Compaction) -> Self {
        compaction.uncompact_edge(self)
    }

    fn remap(self, remapping: &Remapping) -> Self {
        remapping.map_edge(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: None,
            edge: Some(self),
        }
    }

    fn anchor(self) -> Option<ParticipantAnchor> {
        Some(ParticipantAnchor::Edge(self))
    }
}

/// Union incidence index: a node → relations pair and an edge → relations pair,
/// each a sorted `(keys, rels)` slice. Participants self-route via `refs()`, so a
/// set with only node participants leaves the edge half empty, and vice versa.
#[derive(Clone, Debug, Default)]
struct Incidence {
    node_keys: Vec<NodeId>,
    node_rels: Vec<RelationId>,
    edge_keys: Vec<EdgeId>,
    edge_rels: Vec<RelationId>,
}

impl Incidence {
    /// `fill(i, out)` pushes every participant's `refs()` for relation `i`.
    fn build(
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

    fn incident(&self, node: NodeId) -> &[RelationId] {
        let start = self.node_keys.partition_point(|n| *n < node);
        let end = start + self.node_keys[start..].partition_point(|n| *n <= node);
        &self.node_rels[start..end]
    }

    fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        let start = self.edge_keys.partition_point(|e| *e < edge);
        let end = start + self.edge_keys[start..].partition_point(|e| *e <= edge);
        &self.edge_rels[start..end]
    }

    fn has_incident(&self, node: NodeId) -> bool {
        self.node_keys.binary_search(&node).is_ok()
    }

    fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.edge_keys.binary_search(&edge).is_ok()
    }
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
    incidence: Incidence,
    _ordering: PhantomData<O>,
}

impl<P: PartialEq, O, D: PartialEq, const N: usize> PartialEq for FixedRelationSet<P, O, D, N> {
    fn eq(&self, other: &Self) -> bool {
        self.participants == other.participants && self.data == other.data
    }
}

impl<P: Eq, O, D: Eq, const N: usize> Eq for FixedRelationSet<P, O, D, N> {}

/// Relabel a factor's participants through a general `Remapping` and re-canonicalize,
/// returning the new participants and the position permutation `σ` (`σ[new] = old`).
/// Total — no participant is dropped.
fn remap_factor<P, O>(
    participants: &[P],
    remapping: &Remapping,
) -> (Vec<P>, Vec<ParticipantPosition>)
where
    P: RelationParticipant,
    O: FactorOrdering,
{
    let mut relabeled: Vec<P> = participants.iter().map(|&p| p.remap(remapping)).collect();
    let positions = O::canonicalize_positions(&mut relabeled);
    (relabeled, positions)
}

/// Multiset equality of a stored factor slice against a pre-sorted query — sort a copy of
/// `stored` and compare (the stored `Ordered` frame is left intact). Matches on identity (the
/// participant multiset), independent of the factor's ordering marker.
fn participants_match<P: RelationParticipant>(stored: &[P], sorted_query: &[P]) -> bool {
    if stored.len() != sorted_query.len() {
        return false;
    }
    let mut sorted_stored: Vec<P> = stored.to_vec();
    sorted_stored.sort_unstable();
    sorted_stored.as_slice() == sorted_query
}

impl<P: RelationParticipant, O: FactorOrdering, D, const N: usize> FixedRelationSet<P, O, D, N> {
    pub fn new(entries: Vec<([P; N], D)>) -> Self {
        let mut participants = Vec::with_capacity(entries.len());
        let mut data = Vec::with_capacity(entries.len());
        for (mut p, d) in entries {
            O::canonicalize(&mut p);
            participants.push(p);
            data.push(d);
        }

        let incidence = Incidence::build(participants.len(), |i, out| {
            out.extend(participants[i].iter().map(|p| p.refs()));
        });

        Self {
            participants,
            data,
            incidence,
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

    /// Id of the relation whose participants equal `query` as a multiset (order-independent),
    /// if any. §4.1 uniqueness ⇒ at most one hit.
    pub fn find_by_participants(&self, query: &[P]) -> Option<RelationId> {
        let mut sorted_query: Vec<P> = query.to_vec();
        sorted_query.sort_unstable();
        let matches = |id: RelationId| participants_match(self.participants(id), &sorted_query);
        match query.iter().find_map(|p| p.anchor()) {
            Some(ParticipantAnchor::Node(node)) => {
                self.incident(node).iter().copied().find(|&id| matches(id))
            }
            Some(ParticipantAnchor::Edge(edge)) => self
                .incident_edge(edge)
                .iter()
                .copied()
                .find(|&id| matches(id)),
            None => self.relation_ids().find(|&id| matches(id)),
        }
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_edge(edge)
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    pub fn apply_compaction(&self, compaction: &Compaction) -> Self
    where
        D: Clone,
    {
        let entries: Vec<([P; N], D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let parts: Option<Vec<P>> = self
                    .participants(rid)
                    .iter()
                    .map(|&p| p.compact(compaction))
                    .collect();
                let parts: [P; N] = parts?.try_into().ok()?;
                Some((parts, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }

    pub fn apply_remapping(&self, remapping: &Remapping) -> (Self, Vec<ParticipantPosition>)
    where
        D: Clone,
    {
        let mut positions = Vec::new();
        let entries: Vec<([P; N], D)> = (0..self.relation_count())
            .map(|i| {
                let rid = RelationId(i as u32);
                let (sorted, sigma) = remap_factor::<P, O>(self.participants(rid), remapping);
                positions.extend(sigma);
                let parts: [P; N] = sorted
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                (parts, self.data(rid).clone())
            })
            .collect();
        (Self::new(entries), positions)
    }

    /// Glue `self` and `right`, both **already in the same participant id-space**, identifying
    /// coinciding relations (equal participants) — the same-space relation pushout. `combine` merges
    /// the data of a coincidence (`None` = ⊥ ⇒ the whole glue is inadmissible ⇒ `None`); every other
    /// relation is carried. `self`'s ids are the identity prefix of the object, `right`'s
    /// non-coinciding relations are appended. The caller brings both sides into the common space with
    /// [`apply_remapping`](Self::apply_remapping) and re-indexes its `D` first.
    pub fn pushout(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<([P; N], D)> = self
            .relation_ids()
            .map(|id| (*self.participants(id), self.data(id).clone()))
            .collect();
        let self_count = entries.len();
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.relation_count());
        for id in right.relation_ids() {
            match self.find_by_participants(right.participants(id)) {
                Some(hit) => {
                    let merged = combine(&entries[hit.index()].1, right.data(id))?;
                    entries[hit.index()].1 = merged;
                    right_map.push(hit);
                }
                None => {
                    right_map.push(RelationId(entries.len() as u32));
                    entries.push((*right.participants(id), right.data(id).clone()));
                }
            }
        }
        let object_count = entries.len();
        Some(relation_pushout(
            Self::new(entries),
            self_count,
            object_count,
            right_map,
        ))
    }
}

impl<P, O, D, const N: usize> Default for FixedRelationSet<P, O, D, N> {
    fn default() -> Self {
        Self {
            participants: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
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
    incidence: Incidence,
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

        let incidence = Incidence::build(relation_count, |i, out| {
            let start = offsets[i] as usize;
            let end = offsets[i + 1] as usize;
            out.extend(participants[start..end].iter().map(|p| p.refs()));
        });

        Self {
            offsets,
            participants,
            data,
            incidence,
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

    /// Id of the relation whose participants equal `query` as a multiset (order-independent),
    /// if any. §4.1 uniqueness ⇒ at most one hit.
    pub fn find_by_participants(&self, query: &[P]) -> Option<RelationId> {
        let mut sorted_query: Vec<P> = query.to_vec();
        sorted_query.sort_unstable();
        let matches = |id: RelationId| participants_match(self.participants(id), &sorted_query);
        match query.iter().find_map(|p| p.anchor()) {
            Some(ParticipantAnchor::Node(node)) => {
                self.incident(node).iter().copied().find(|&id| matches(id))
            }
            Some(ParticipantAnchor::Edge(edge)) => self
                .incident_edge(edge)
                .iter()
                .copied()
                .find(|&id| matches(id)),
            None => self.relation_ids().find(|&id| matches(id)),
        }
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_edge(edge)
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    pub fn apply_compaction(&self, compaction: &Compaction) -> Self
    where
        D: Clone,
    {
        let entries: Vec<(Vec<P>, D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let parts: Option<Vec<P>> = self
                    .participants(rid)
                    .iter()
                    .map(|&p| p.compact(compaction))
                    .collect();
                Some((parts?, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }

    pub fn apply_remapping(&self, remapping: &Remapping) -> (Self, Vec<ParticipantPosition>)
    where
        D: Clone,
    {
        let mut positions = Vec::new();
        let entries: Vec<(Vec<P>, D)> = (0..self.relation_count())
            .map(|i| {
                let rid = RelationId(i as u32);
                let (sorted, sigma) = remap_factor::<P, O>(self.participants(rid), remapping);
                positions.extend(sigma);
                (sorted, self.data(rid).clone())
            })
            .collect();
        (Self::new(entries), positions)
    }

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`].
    pub fn pushout(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<(Vec<P>, D)> = self
            .relation_ids()
            .map(|id| (self.participants(id).to_vec(), self.data(id).clone()))
            .collect();
        let self_count = entries.len();
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.relation_count());
        for id in right.relation_ids() {
            match self.find_by_participants(right.participants(id)) {
                Some(hit) => {
                    let merged = combine(&entries[hit.index()].1, right.data(id))?;
                    entries[hit.index()].1 = merged;
                    right_map.push(hit);
                }
                None => {
                    right_map.push(RelationId(entries.len() as u32));
                    entries.push((right.participants(id).to_vec(), right.data(id).clone()));
                }
            }
        }
        let object_count = entries.len();
        Some(relation_pushout(
            Self::new(entries),
            self_count,
            object_count,
            right_map,
        ))
    }
}

impl<P, O, D> Default for VarRelationSet<P, O, D> {
    fn default() -> Self {
        Self {
            offsets: vec![0],
            participants: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
            _ordering: PhantomData,
        }
    }
}

/// Birelation with two fixed-arity factors.
#[derive(Clone, Debug)]
pub struct FixedFixedBirelationSet<L1, O1, const N1: usize, L2, O2, const N2: usize, D> {
    participants_1: Vec<[L1; N1]>,
    participants_2: Vec<[L2; N2]>,
    data: Vec<D>,
    incidence: Incidence,
    _ordering: PhantomData<(O1, O2)>,
}

impl<L1, O1, const N1: usize, L2, O2, const N2: usize, D> PartialEq
    for FixedFixedBirelationSet<L1, O1, N1, L2, O2, N2, D>
where
    L1: PartialEq,
    L2: PartialEq,
    D: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.participants_1 == other.participants_1
            && self.participants_2 == other.participants_2
            && self.data == other.data
    }
}

impl<L1, O1, const N1: usize, L2, O2, const N2: usize, D> Eq
    for FixedFixedBirelationSet<L1, O1, N1, L2, O2, N2, D>
where
    L1: Eq,
    L2: Eq,
    D: Eq,
{
}

impl<L1, O1, const N1: usize, L2, O2, const N2: usize, D>
    FixedFixedBirelationSet<L1, O1, N1, L2, O2, N2, D>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
{
    pub fn new(entries: Vec<([L1; N1], [L2; N2], D)>) -> Self {
        let relation_count = entries.len();
        let mut participants_1 = Vec::with_capacity(relation_count);
        let mut participants_2 = Vec::with_capacity(relation_count);
        let mut data = Vec::with_capacity(relation_count);
        for (mut l1, mut l2, d) in entries {
            O1::canonicalize(&mut l1);
            O2::canonicalize(&mut l2);
            participants_1.push(l1);
            participants_2.push(l2);
            data.push(d);
        }
        let incidence = Incidence::build(relation_count, |i, out| {
            out.extend(participants_1[i].iter().map(|p| p.refs()));
            out.extend(participants_2[i].iter().map(|p| p.refs()));
        });
        Self {
            participants_1,
            participants_2,
            data,
            incidence,
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

    pub fn participants_1(&self, id: RelationId) -> &[L1; N1] {
        &self.participants_1[id.index()]
    }

    pub fn participants_2(&self, id: RelationId) -> &[L2; N2] {
        &self.participants_2[id.index()]
    }

    /// Id of the relation whose two factors equal `query_1` / `query_2` as multisets
    /// (order-independent per factor), if any. §4.1 uniqueness ⇒ at most one hit.
    pub fn find_by_participants(&self, query_1: &[L1], query_2: &[L2]) -> Option<RelationId> {
        let mut sorted_1: Vec<L1> = query_1.to_vec();
        sorted_1.sort_unstable();
        let mut sorted_2: Vec<L2> = query_2.to_vec();
        sorted_2.sort_unstable();
        let matches = |id: RelationId| {
            participants_match(self.participants_1(id), &sorted_1)
                && participants_match(self.participants_2(id), &sorted_2)
        };
        match query_1
            .iter()
            .find_map(|p| p.anchor())
            .or_else(|| query_2.iter().find_map(|p| p.anchor()))
        {
            Some(ParticipantAnchor::Node(node)) => {
                self.incident(node).iter().copied().find(|&id| matches(id))
            }
            Some(ParticipantAnchor::Edge(edge)) => self
                .incident_edge(edge)
                .iter()
                .copied()
                .find(|&id| matches(id)),
            None => self.relation_ids().find(|&id| matches(id)),
        }
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_edge(edge)
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    pub fn apply_compaction(&self, compaction: &Compaction) -> Self
    where
        D: Clone,
    {
        let entries: Vec<([L1; N1], [L2; N2], D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let f1: Option<Vec<L1>> = self
                    .participants_1(rid)
                    .iter()
                    .map(|&p| p.compact(compaction))
                    .collect();
                let f1: [L1; N1] = f1?.try_into().ok()?;
                let f2: Option<Vec<L2>> = self
                    .participants_2(rid)
                    .iter()
                    .map(|&p| p.compact(compaction))
                    .collect();
                let f2: [L2; N2] = f2?.try_into().ok()?;
                Some((f1, f2, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }

    pub fn apply_remapping(
        &self,
        remapping: &Remapping,
    ) -> (Self, Vec<ParticipantPosition>, Vec<ParticipantPosition>)
    where
        D: Clone,
    {
        let mut positions_1 = Vec::new();
        let mut positions_2 = Vec::new();
        let entries: Vec<([L1; N1], [L2; N2], D)> = (0..self.relation_count())
            .map(|i| {
                let rid = RelationId(i as u32);
                let (s1, sigma1) = remap_factor::<L1, O1>(self.participants_1(rid), remapping);
                let (s2, sigma2) = remap_factor::<L2, O2>(self.participants_2(rid), remapping);
                positions_1.extend(sigma1);
                positions_2.extend(sigma2);
                let f1: [L1; N1] = s1
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                let f2: [L2; N2] = s2
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                (f1, f2, self.data(rid).clone())
            })
            .collect();
        (Self::new(entries), positions_1, positions_2)
    }
}

impl<L1, O1, const N1: usize, L2, O2, const N2: usize, D> Default
    for FixedFixedBirelationSet<L1, O1, N1, L2, O2, N2, D>
{
    fn default() -> Self {
        Self {
            participants_1: Vec::new(),
            participants_2: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
            _ordering: PhantomData,
        }
    }
}

/// Birelation with a fixed-arity factor 1 and a variable-arity factor 2. Each factor
/// has its own participant type and ordering; the union incidence spans both.
#[derive(Clone, Debug)]
pub struct FixedVarBirelationSet<L1, O1, const N1: usize, L2, O2, D> {
    participants_1: Vec<[L1; N1]>,
    f2_offsets: Vec<u32>,
    participants_2: Vec<L2>,
    data: Vec<D>,
    incidence: Incidence,
    _ordering: PhantomData<(O1, O2)>,
}

impl<L1, O1, const N1: usize, L2, O2, D> PartialEq for FixedVarBirelationSet<L1, O1, N1, L2, O2, D>
where
    L1: PartialEq,
    L2: PartialEq,
    D: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.participants_1 == other.participants_1
            && self.f2_offsets == other.f2_offsets
            && self.participants_2 == other.participants_2
            && self.data == other.data
    }
}

impl<L1, O1, const N1: usize, L2, O2, D> Eq for FixedVarBirelationSet<L1, O1, N1, L2, O2, D>
where
    L1: Eq,
    L2: Eq,
    D: Eq,
{
}

impl<L1, O1, const N1: usize, L2, O2, D> FixedVarBirelationSet<L1, O1, N1, L2, O2, D>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
{
    pub fn new(entries: Vec<([L1; N1], Vec<L2>, D)>) -> Self {
        let relation_count = entries.len();
        let mut participants_1 = Vec::with_capacity(relation_count);
        let mut f2_offsets = Vec::with_capacity(relation_count + 1);
        f2_offsets.push(0);
        let mut participants_2 = Vec::new();
        let mut data = Vec::with_capacity(relation_count);
        for (mut l1, mut l2, d) in entries {
            O1::canonicalize(&mut l1);
            participants_1.push(l1);
            O2::canonicalize(&mut l2);
            participants_2.extend_from_slice(&l2);
            f2_offsets.push(participants_2.len() as u32);
            data.push(d);
        }
        let incidence = Incidence::build(relation_count, |i, out| {
            out.extend(participants_1[i].iter().map(|p| p.refs()));
            let start = f2_offsets[i] as usize;
            let end = f2_offsets[i + 1] as usize;
            out.extend(participants_2[start..end].iter().map(|p| p.refs()));
        });
        Self {
            participants_1,
            f2_offsets,
            participants_2,
            data,
            incidence,
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

    pub fn participants_1(&self, id: RelationId) -> &[L1; N1] {
        &self.participants_1[id.index()]
    }

    pub fn participants_2(&self, id: RelationId) -> &[L2] {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        &self.participants_2[start..end]
    }

    /// Id of the relation whose two factors equal `query_1` / `query_2` as multisets
    /// (order-independent per factor), if any. §4.1 uniqueness ⇒ at most one hit.
    pub fn find_by_participants(&self, query_1: &[L1], query_2: &[L2]) -> Option<RelationId> {
        let mut sorted_1: Vec<L1> = query_1.to_vec();
        sorted_1.sort_unstable();
        let mut sorted_2: Vec<L2> = query_2.to_vec();
        sorted_2.sort_unstable();
        let matches = |id: RelationId| {
            participants_match(self.participants_1(id), &sorted_1)
                && participants_match(self.participants_2(id), &sorted_2)
        };
        match query_1
            .iter()
            .find_map(|p| p.anchor())
            .or_else(|| query_2.iter().find_map(|p| p.anchor()))
        {
            Some(ParticipantAnchor::Node(node)) => {
                self.incident(node).iter().copied().find(|&id| matches(id))
            }
            Some(ParticipantAnchor::Edge(edge)) => self
                .incident_edge(edge)
                .iter()
                .copied()
                .find(|&id| matches(id)),
            None => self.relation_ids().find(|&id| matches(id)),
        }
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_edge(edge)
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    pub fn apply_compaction(&self, compaction: &Compaction) -> Self
    where
        D: Clone,
    {
        let entries: Vec<([L1; N1], Vec<L2>, D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let f1: Option<Vec<L1>> = self
                    .participants_1(rid)
                    .iter()
                    .map(|&p| p.compact(compaction))
                    .collect();
                let f1: [L1; N1] = f1?.try_into().ok()?;
                let f2: Option<Vec<L2>> = self
                    .participants_2(rid)
                    .iter()
                    .map(|&p| p.compact(compaction))
                    .collect();
                Some((f1, f2?, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }

    pub fn apply_remapping(
        &self,
        remapping: &Remapping,
    ) -> (Self, Vec<ParticipantPosition>, Vec<ParticipantPosition>)
    where
        D: Clone,
    {
        let mut positions_1 = Vec::new();
        let mut positions_2 = Vec::new();
        let entries: Vec<([L1; N1], Vec<L2>, D)> = (0..self.relation_count())
            .map(|i| {
                let rid = RelationId(i as u32);
                let (s1, sigma1) = remap_factor::<L1, O1>(self.participants_1(rid), remapping);
                let (s2, sigma2) = remap_factor::<L2, O2>(self.participants_2(rid), remapping);
                positions_1.extend(sigma1);
                positions_2.extend(sigma2);
                let f1: [L1; N1] = s1
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                (f1, s2, self.data(rid).clone())
            })
            .collect();
        (Self::new(entries), positions_1, positions_2)
    }

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`]. Coincidence is equality of
    /// **both** factors' participants.
    pub fn pushout(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<([L1; N1], Vec<L2>, D)> = self
            .relation_ids()
            .map(|id| {
                (
                    *self.participants_1(id),
                    self.participants_2(id).to_vec(),
                    self.data(id).clone(),
                )
            })
            .collect();
        let self_count = entries.len();
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.relation_count());
        for id in right.relation_ids() {
            match self.find_by_participants(right.participants_1(id), right.participants_2(id)) {
                Some(hit) => {
                    let merged = combine(&entries[hit.index()].2, right.data(id))?;
                    entries[hit.index()].2 = merged;
                    right_map.push(hit);
                }
                None => {
                    right_map.push(RelationId(entries.len() as u32));
                    entries.push((
                        *right.participants_1(id),
                        right.participants_2(id).to_vec(),
                        right.data(id).clone(),
                    ));
                }
            }
        }
        let object_count = entries.len();
        Some(relation_pushout(
            Self::new(entries),
            self_count,
            object_count,
            right_map,
        ))
    }
}

impl<L1, O1, const N1: usize, L2, O2, D> Default for FixedVarBirelationSet<L1, O1, N1, L2, O2, D> {
    fn default() -> Self {
        Self {
            participants_1: Vec::new(),
            f2_offsets: vec![0],
            participants_2: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
            _ordering: PhantomData,
        }
    }
}

/// Birelation with two variable-arity factors.
#[derive(Clone, Debug)]
pub struct VarVarBirelationSet<L1, O1, L2, O2, D> {
    f1_offsets: Vec<u32>,
    participants_1: Vec<L1>,
    f2_offsets: Vec<u32>,
    participants_2: Vec<L2>,
    data: Vec<D>,
    incidence: Incidence,
    _ordering: PhantomData<(O1, O2)>,
}

impl<L1, O1, L2, O2, D> PartialEq for VarVarBirelationSet<L1, O1, L2, O2, D>
where
    L1: PartialEq,
    L2: PartialEq,
    D: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.f1_offsets == other.f1_offsets
            && self.participants_1 == other.participants_1
            && self.f2_offsets == other.f2_offsets
            && self.participants_2 == other.participants_2
            && self.data == other.data
    }
}

impl<L1, O1, L2, O2, D> Eq for VarVarBirelationSet<L1, O1, L2, O2, D>
where
    L1: Eq,
    L2: Eq,
    D: Eq,
{
}

impl<L1, O1, L2, O2, D> VarVarBirelationSet<L1, O1, L2, O2, D>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
{
    pub fn new(entries: Vec<(Vec<L1>, Vec<L2>, D)>) -> Self {
        let relation_count = entries.len();
        let mut f1_offsets = Vec::with_capacity(relation_count + 1);
        f1_offsets.push(0);
        let mut participants_1 = Vec::new();
        let mut f2_offsets = Vec::with_capacity(relation_count + 1);
        f2_offsets.push(0);
        let mut participants_2 = Vec::new();
        let mut data = Vec::with_capacity(relation_count);
        for (mut l1, mut l2, d) in entries {
            O1::canonicalize(&mut l1);
            participants_1.extend_from_slice(&l1);
            f1_offsets.push(participants_1.len() as u32);
            O2::canonicalize(&mut l2);
            participants_2.extend_from_slice(&l2);
            f2_offsets.push(participants_2.len() as u32);
            data.push(d);
        }
        let incidence = Incidence::build(relation_count, |i, out| {
            let s1 = f1_offsets[i] as usize;
            let e1 = f1_offsets[i + 1] as usize;
            out.extend(participants_1[s1..e1].iter().map(|p| p.refs()));
            let s2 = f2_offsets[i] as usize;
            let e2 = f2_offsets[i + 1] as usize;
            out.extend(participants_2[s2..e2].iter().map(|p| p.refs()));
        });
        Self {
            f1_offsets,
            participants_1,
            f2_offsets,
            participants_2,
            data,
            incidence,
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

    pub fn participants_1(&self, id: RelationId) -> &[L1] {
        let start = self.f1_offsets[id.index()] as usize;
        let end = self.f1_offsets[id.index() + 1] as usize;
        &self.participants_1[start..end]
    }

    pub fn participants_2(&self, id: RelationId) -> &[L2] {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        &self.participants_2[start..end]
    }

    /// Id of the relation whose two factors equal `query_1` / `query_2` as multisets
    /// (order-independent per factor), if any. §4.1 uniqueness ⇒ at most one hit.
    pub fn find_by_participants(&self, query_1: &[L1], query_2: &[L2]) -> Option<RelationId> {
        let mut sorted_1: Vec<L1> = query_1.to_vec();
        sorted_1.sort_unstable();
        let mut sorted_2: Vec<L2> = query_2.to_vec();
        sorted_2.sort_unstable();
        let matches = |id: RelationId| {
            participants_match(self.participants_1(id), &sorted_1)
                && participants_match(self.participants_2(id), &sorted_2)
        };
        match query_1
            .iter()
            .find_map(|p| p.anchor())
            .or_else(|| query_2.iter().find_map(|p| p.anchor()))
        {
            Some(ParticipantAnchor::Node(node)) => {
                self.incident(node).iter().copied().find(|&id| matches(id))
            }
            Some(ParticipantAnchor::Edge(edge)) => self
                .incident_edge(edge)
                .iter()
                .copied()
                .find(|&id| matches(id)),
            None => self.relation_ids().find(|&id| matches(id)),
        }
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_edge(edge)
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    pub fn apply_compaction(&self, compaction: &Compaction) -> Self
    where
        D: Clone,
    {
        let entries: Vec<(Vec<L1>, Vec<L2>, D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let f1: Option<Vec<L1>> = self
                    .participants_1(rid)
                    .iter()
                    .map(|&p| p.compact(compaction))
                    .collect();
                let f2: Option<Vec<L2>> = self
                    .participants_2(rid)
                    .iter()
                    .map(|&p| p.compact(compaction))
                    .collect();
                Some((f1?, f2?, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }

    pub fn apply_remapping(
        &self,
        remapping: &Remapping,
    ) -> (Self, Vec<ParticipantPosition>, Vec<ParticipantPosition>)
    where
        D: Clone,
    {
        let mut positions_1 = Vec::new();
        let mut positions_2 = Vec::new();
        let entries: Vec<(Vec<L1>, Vec<L2>, D)> = (0..self.relation_count())
            .map(|i| {
                let rid = RelationId(i as u32);
                let (s1, sigma1) = remap_factor::<L1, O1>(self.participants_1(rid), remapping);
                let (s2, sigma2) = remap_factor::<L2, O2>(self.participants_2(rid), remapping);
                positions_1.extend(sigma1);
                positions_2.extend(sigma2);
                (s1, s2, self.data(rid).clone())
            })
            .collect();
        (Self::new(entries), positions_1, positions_2)
    }
}

impl<L1, O1, L2, O2, D> Default for VarVarBirelationSet<L1, O1, L2, O2, D> {
    fn default() -> Self {
        Self {
            f1_offsets: vec![0],
            participants_1: Vec::new(),
            f2_offsets: vec![0],
            participants_2: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
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
    fn test_node_id_compact(#[case] id: NodeId, #[case] expected: Option<NodeId>) {
        let compaction = Compaction::new(vec![1], vec![]);
        assert_eq!(id.compact(&compaction), expected);
    }

    #[rstest]
    #[case::before_gap(NodeId(0), NodeId(0))]
    #[case::after_gap(NodeId(1), NodeId(2))]
    fn test_node_id_unmap(#[case] id: NodeId, #[case] expected: NodeId) {
        let compaction = Compaction::new(vec![1], vec![]);
        assert_eq!(id.uncompact(&compaction), expected);
    }

    #[rstest]
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
    fn test_edge_id_compact(#[case] id: EdgeId, #[case] expected: Option<EdgeId>) {
        let compaction = Compaction::new(vec![], vec![0]);
        assert_eq!(id.compact(&compaction), expected);
    }

    #[rstest]
    #[case::before_gap(EdgeId(0), EdgeId(0))]
    #[case::after_gap(EdgeId(1), EdgeId(2))]
    fn test_edge_id_unmap(#[case] id: EdgeId, #[case] expected: EdgeId) {
        let compaction = Compaction::new(vec![], vec![1]);
        assert_eq!(id.uncompact(&compaction), expected);
    }

    #[rstest]
    fn test_edge_id_refs() {
        assert_eq!(
            EdgeId(2).refs(),
            ParticipantRefs {
                node: None,
                edge: Some(EdgeId(2)),
            }
        );
    }

    #[rstest]
    fn test_fixed_relation_set_new() {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], "dative"), ([n(1), n(2)], "noncov")]);
        assert_eq!(rs.relation_count(), 2);
        assert_eq!(rs.data(RelationId(0)), &"dative");
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(2)]);
    }

    #[rstest]
    fn test_fixed_relation_set_data_mut() {
        let mut rs: FixedRelationSet<NodeId, Unordered, i32, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
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

    #[rstest]
    fn test_fixed_relation_set_participants_sorted() {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(2), n(0)], "a"), ([n(3), n(1)], "b")]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(2)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(3)]);
    }

    #[rstest]
    fn test_fixed_relation_set_participants_ordered() {
        let rs: FixedRelationSet<NodeId, Ordered, &str, 2> =
            FixedRelationSet::new(vec![([n(2), n(0)], "a"), ([n(3), n(1)], "b")]);
        assert_eq!(rs.participants(RelationId(0)), &[n(2), n(0)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(3), n(1)]);
    }

    #[rstest]
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

    #[rstest]
    #[case::first(RelationId(0), true)]
    #[case::last(RelationId(1), true)]
    #[case::out_of_range(RelationId(2), false)]
    fn test_fixed_relation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
        let rs: FixedRelationSet<NodeId, Unordered, (), 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], ()), ([n(1), n(2)], ())]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_fixed_relation_set_relation_ids() {
        let rs: FixedRelationSet<NodeId, Unordered, (), 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], ()), ([n(1), n(2)], ())]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_fixed_relation_set_apply_compaction() {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(0), n(2)], "keep"), ([n(1), n(3)], "drop")]);
        let compaction = Compaction::new(vec![1], vec![]);
        let out = rs.apply_compaction(&compaction);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_fixed_relation_set_apply_remapping() {
        // swap the pair: [0,1] relabeled to [1,0], re-sorted to [0,1]; σ = [1,0]
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], "x")]);
        let remapping = Remapping::new(vec![n(1), n(0)], vec![]);
        let (out, positions) = rs.apply_remapping(&remapping);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.data(RelationId(0)), &"x");
        assert_eq!(
            positions,
            vec![ParticipantPosition(1), ParticipantPosition(0)]
        );
    }

    #[rstest]
    fn test_fixed_relation_set_default() {
        let rs = FixedRelationSet::<NodeId, Unordered, (), 2>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
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

    #[rstest]
    fn test_var_relation_set_data_mut() {
        let mut rs: VarRelationSet<NodeId, Unordered, i32> =
            VarRelationSet::new(vec![(vec![n(0), n(1), n(2)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
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

    #[rstest]
    fn test_var_relation_set_participants_sorted() {
        let rs: VarRelationSet<NodeId, Unordered, ()> = VarRelationSet::new(vec![
            (vec![n(5), n(2), n(0), n(3)], ()),
            (vec![n(4), n(1)], ()),
        ]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(2), n(3), n(5)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(4)]);
    }

    #[rstest]
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

    #[rstest]
    fn test_var_relation_set_variable_arity() {
        let rs: VarRelationSet<NodeId, Unordered, &str> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], "pair"),
            (vec![n(2), n(3), n(4), n(5)], "quad"),
        ]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(2), n(3), n(4), n(5)]);
    }

    #[rstest]
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

    #[rstest]
    fn test_var_relation_set_edge_incidence() {
        let rs: VarRelationSet<EdgeId, Unordered, &str> = VarRelationSet::new(vec![
            (vec![EdgeId(0), EdgeId(2)], "a"),
            (vec![EdgeId(1), EdgeId(2)], "b"),
        ]);
        assert_eq!(rs.incident_edge(EdgeId(2)), &[RelationId(0), RelationId(1)]);
        assert_eq!(rs.incident_edge(EdgeId(0)), &[RelationId(0)]);
        assert!(rs.has_incident_edge(EdgeId(2)));
        assert!(!rs.has_incident_edge(EdgeId(5)));
        assert!(rs.incident(n(0)).is_empty());
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    #[case::first(RelationId(0), true)]
    #[case::out_of_range(RelationId(1), false)]
    fn test_var_relation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
        let rs: VarRelationSet<NodeId, Unordered, ()> =
            VarRelationSet::new(vec![(vec![n(0), n(1)], ())]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_var_relation_set_relation_ids() {
        let rs: VarRelationSet<NodeId, Unordered, ()> =
            VarRelationSet::new(vec![(vec![n(0), n(1)], ()), (vec![n(1), n(2)], ())]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_var_relation_set_apply_compaction() {
        let rs: VarRelationSet<NodeId, Unordered, &str> = VarRelationSet::new(vec![
            (vec![n(0), n(2), n(4)], "keep"),
            (vec![n(1), n(3)], "drop"),
        ]);
        let compaction = Compaction::new(vec![1], vec![]);
        let out = rs.apply_compaction(&compaction);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1), n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_var_relation_set_apply_remapping() {
        // relabel 0→2, 1→0, 2→1: [0,1,2] relabeled to [2,0,1], re-sorted to [0,1,2]; σ = [1,2,0]
        let rs: VarRelationSet<NodeId, Unordered, &str> =
            VarRelationSet::new(vec![(vec![n(0), n(1), n(2)], "x")]);
        let remapping = Remapping::new(vec![n(2), n(0), n(1)], vec![]);
        let (out, positions) = rs.apply_remapping(&remapping);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1), n(2)]);
        assert_eq!(out.data(RelationId(0)), &"x");
        assert_eq!(
            positions,
            vec![
                ParticipantPosition(1),
                ParticipantPosition(2),
                ParticipantPosition(0)
            ]
        );
    }

    #[rstest]
    fn test_var_relation_set_default() {
        let rs = VarRelationSet::<NodeId, Unordered, ()>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_new() {
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 2, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(2), n(1)], "x")]);
        assert_eq!(rs.relation_count(), 1);
        assert_eq!(rs.participants_1(RelationId(0)), &[n(0)]);
        assert_eq!(rs.participants_2(RelationId(0)), &[n(1), n(2)]);
        assert_eq!(rs.data(RelationId(0)), &"x");
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_data_mut() {
        let mut rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 1, i32> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_data_iter_mut() {
        let mut rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 1, i32> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], 1), ([n(2)], [n(3)], 2)]);
        for d in rs.data_iter_mut() {
            *d *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_incidence() {
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, EdgeId, Unordered, 1, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [EdgeId(7)], "x")]);
        assert_eq!(rs.incident(n(0)), &[RelationId(0)]);
        assert_eq!(rs.incident_edge(EdgeId(7)), &[RelationId(0)]);
        assert!(rs.has_incident(n(0)));
        assert!(rs.has_incident_edge(EdgeId(7)));
        assert!(!rs.has_incident(n(5)));
    }

    #[rstest]
    #[case::first(RelationId(0), true)]
    #[case::out_of_range(RelationId(1), false)]
    fn test_fixed_fixed_birelation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 1, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], "x")]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_relation_ids() {
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 1, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], "a"), ([n(2)], [n(3)], "b")]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_apply_compaction() {
        // dropped relation loses a factor-1 participant
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 2, &str> =
            FixedFixedBirelationSet::new(vec![
                ([n(0)], [n(2), n(4)], "keep"),
                ([n(1)], [n(5), n(6)], "drop"),
            ]);
        let compaction = Compaction::new(vec![1], vec![]);
        let out = rs.apply_compaction(&compaction);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(1), n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_apply_remapping() {
        // factor-1 [0,1] swaps (σ₁ = [1,0]); factor-2 [2] is fixed (σ₂ = [0])
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 2, NodeId, Unordered, 1, &str> =
            FixedFixedBirelationSet::new(vec![([n(0), n(1)], [n(2)], "x")]);
        let remapping = Remapping::new(vec![n(1), n(0), n(2)], vec![]);
        let (out, p1, p2) = rs.apply_remapping(&remapping);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(2)]);
        assert_eq!(p1, vec![ParticipantPosition(1), ParticipantPosition(0)]);
        assert_eq!(p2, vec![ParticipantPosition(0)]);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_default() {
        let rs =
            FixedFixedBirelationSet::<NodeId, Unordered, 1, NodeId, Unordered, 1, ()>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_fixed_var_birelation_set_new() {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1), n(2), n(3)], "ct")]);
        assert_eq!(rs.relation_count(), 1);
        assert_eq!(rs.participants_1(RelationId(0)), &[EdgeId(0)]);
        assert_eq!(rs.participants_2(RelationId(0)), &[n(1), n(2), n(3)]);
        assert_eq!(rs.data(RelationId(0)), &"ct");
    }

    #[rstest]
    fn test_fixed_var_birelation_set_data_mut() {
        let mut rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, i32> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_data_iter_mut() {
        let mut rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, i32> =
            FixedVarBirelationSet::new(vec![
                ([EdgeId(0)], vec![n(1)], 1),
                ([EdgeId(1)], vec![n(2)], 2),
            ]);
        for d in rs.data_iter_mut() {
            *d *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_incidence() {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1), n(2)], "ct")]);
        assert_eq!(rs.incident_edge(EdgeId(0)), &[RelationId(0)]);
        assert_eq!(rs.incident(n(2)), &[RelationId(0)]);
        assert!(rs.has_incident_edge(EdgeId(0)));
        assert!(rs.has_incident(n(1)));
        assert!(rs.incident(n(0)).is_empty());
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    #[case::first(RelationId(0), true)]
    #[case::out_of_range(RelationId(1), false)]
    fn test_fixed_var_birelation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1)], "ct")]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_relation_ids() {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![
                ([EdgeId(0)], vec![n(1)], "a"),
                ([EdgeId(1)], vec![n(2)], "b"),
            ]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_apply_compaction() {
        let rs: FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![
                ([n(0)], vec![n(2), n(4)], "keep"),
                ([n(5)], vec![n(1), n(3)], "drop"),
            ]);
        let compaction = Compaction::new(vec![1], vec![]);
        let out = rs.apply_compaction(&compaction);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(1), n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_fixed_var_birelation_set_apply_remapping() {
        // factor-1 [0] fixed (σ₁ = [0]); factor-2 [1,2] relabeled to [3,1], re-sorted to [1,3] (σ₂ = [1,0])
        let rs: FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, &str> =
            FixedVarBirelationSet::new(vec![([n(0)], vec![n(1), n(2)], "x")]);
        let remapping = Remapping::new(vec![n(0), n(3), n(1)], vec![]);
        let (out, p1, p2) = rs.apply_remapping(&remapping);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(1), n(3)]);
        assert_eq!(p1, vec![ParticipantPosition(0)]);
        assert_eq!(p2, vec![ParticipantPosition(1), ParticipantPosition(0)]);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_default() {
        let rs = FixedVarBirelationSet::<EdgeId, Ordered, 1, NodeId, Ordered, ()>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident_edge(EdgeId(0)));
    }

    #[rstest]
    fn test_var_var_birelation_set_new() {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![(vec![n(0), n(1)], vec![EdgeId(5)], "y")]);
        assert_eq!(rs.relation_count(), 1);
        assert_eq!(rs.participants_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants_2(RelationId(0)), &[EdgeId(5)]);
        assert_eq!(rs.data(RelationId(0)), &"y");
    }

    #[rstest]
    fn test_var_var_birelation_set_data_mut() {
        let mut rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, i32> =
            VarVarBirelationSet::new(vec![(vec![n(0)], vec![EdgeId(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_var_var_birelation_set_data_iter_mut() {
        let mut rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, i32> =
            VarVarBirelationSet::new(vec![
                (vec![n(0)], vec![EdgeId(1)], 1),
                (vec![n(2)], vec![EdgeId(3)], 2),
            ]);
        for d in rs.data_iter_mut() {
            *d *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_var_var_birelation_set_incidence() {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![(vec![n(0), n(1)], vec![EdgeId(5)], "y")]);
        assert_eq!(rs.incident(n(1)), &[RelationId(0)]);
        assert_eq!(rs.incident_edge(EdgeId(5)), &[RelationId(0)]);
        assert!(rs.has_incident(n(0)));
        assert!(rs.has_incident_edge(EdgeId(5)));
        assert!(!rs.has_incident_edge(EdgeId(0)));
    }

    #[rstest]
    #[case::first(RelationId(0), true)]
    #[case::out_of_range(RelationId(1), false)]
    fn test_var_var_birelation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![(vec![n(0)], vec![EdgeId(1)], "y")]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_var_var_birelation_set_relation_ids() {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![
                (vec![n(0)], vec![EdgeId(1)], "a"),
                (vec![n(2)], vec![EdgeId(3)], "b"),
            ]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_var_var_birelation_set_apply_compaction() {
        // dropped relation loses a factor-2 participant
        let rs: VarVarBirelationSet<NodeId, Unordered, NodeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![
                (vec![n(0), n(2)], vec![n(4)], "keep"),
                (vec![n(5)], vec![n(1)], "drop"),
            ]);
        let compaction = Compaction::new(vec![1], vec![]);
        let out = rs.apply_compaction(&compaction);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_var_var_birelation_set_apply_remapping() {
        // both factors swap: factor-1 [0,1]→[1,0]→[0,1] (σ₁ = [1,0]); factor-2 [2,3]→[3,2]→[2,3] (σ₂ = [1,0])
        let rs: VarVarBirelationSet<NodeId, Unordered, NodeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![(vec![n(0), n(1)], vec![n(2), n(3)], "x")]);
        let remapping = Remapping::new(vec![n(1), n(0), n(3), n(2)], vec![]);
        let (out, p1, p2) = rs.apply_remapping(&remapping);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(2), n(3)]);
        assert_eq!(p1, vec![ParticipantPosition(1), ParticipantPosition(0)]);
        assert_eq!(p2, vec![ParticipantPosition(1), ParticipantPosition(0)]);
    }

    #[rstest]
    fn test_var_var_birelation_set_default() {
        let rs = VarVarBirelationSet::<NodeId, Unordered, EdgeId, Unordered, ()>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    #[case::exact(vec![n(0), n(1)], Some(RelationId(0)))]
    #[case::reordered(vec![n(1), n(0)], Some(RelationId(0)))]
    #[case::second(vec![n(2), n(3)], Some(RelationId(1)))]
    #[case::absent(vec![n(0), n(3)], None)]
    #[case::wrong_arity(vec![n(0)], None)]
    fn test_fixed_relation_set_find_by_participants(
        #[case] query: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        let rs: FixedRelationSet<NodeId, Unordered, (), 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], ()), ([n(2), n(3)], ())]);
        assert_eq!(rs.find_by_participants(&query), expected);
    }

    #[rstest]
    #[case::exact(vec![n(0), n(1), n(2)], Some(RelationId(0)))]
    #[case::reordered(vec![n(2), n(0), n(1)], Some(RelationId(0)))]
    #[case::second(vec![n(3), n(4)], Some(RelationId(1)))]
    #[case::subset(vec![n(0), n(1)], None)]
    #[case::superset(vec![n(0), n(1), n(2), n(3)], None)]
    fn test_var_relation_set_find_by_participants(
        #[case] query: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        let rs: VarRelationSet<NodeId, Unordered, ()> =
            VarRelationSet::new(vec![(vec![n(0), n(1), n(2)], ()), (vec![n(3), n(4)], ())]);
        assert_eq!(rs.find_by_participants(&query), expected);
    }

    #[rstest]
    #[case::exact(vec![n(0), n(1)], vec![n(2)], Some(RelationId(0)))]
    #[case::reordered_factor(vec![n(1), n(0)], vec![n(2)], Some(RelationId(0)))]
    #[case::second(vec![n(3), n(4)], vec![n(5)], Some(RelationId(1)))]
    #[case::absent(vec![n(0), n(1)], vec![n(9)], None)]
    fn test_fixed_fixed_birelation_set_find_by_participants(
        #[case] query_1: Vec<NodeId>,
        #[case] query_2: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 2, NodeId, Unordered, 1, ()> =
            FixedFixedBirelationSet::new(vec![
                ([n(0), n(1)], [n(2)], ()),
                ([n(3), n(4)], [n(5)], ()),
            ]);
        assert_eq!(rs.find_by_participants(&query_1, &query_2), expected);
    }

    #[rstest]
    #[case::exact(vec![n(0)], vec![n(1)], Some(RelationId(0)))]
    #[case::role_swap(vec![n(1)], vec![n(0)], None)]
    #[case::multiset_reordered(vec![n(3)], vec![n(5), n(4), n(4)], Some(RelationId(1)))]
    #[case::wrong_multiplicity(vec![n(3)], vec![n(4), n(5)], None)]
    #[case::absent(vec![n(0)], vec![n(2)], None)]
    fn test_fixed_var_birelation_set_find_by_participants(
        #[case] query_1: Vec<NodeId>,
        #[case] query_2: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        // Factor2 `Ordered` (a coset frame) yet matched as a multiset: duplicate `n(4)` and the
        // role-swap case exercise the key semantics.
        let rs: FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Ordered, ()> =
            FixedVarBirelationSet::new(vec![
                ([n(0)], vec![n(1)], ()),
                ([n(3)], vec![n(4), n(4), n(5)], ()),
            ]);
        assert_eq!(rs.find_by_participants(&query_1, &query_2), expected);
    }

    #[rstest]
    #[case::exact(vec![n(1), n(2)], Some(RelationId(0)))]
    #[case::reordered(vec![n(2), n(1)], Some(RelationId(0)))]
    #[case::absent(vec![n(1), n(3)], None)]
    fn test_fixed_var_birelation_set_find_by_participants_edge_anchor(
        #[case] ligands: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        // Stereo-bond-like: factor1 is an `EdgeId` site, so the anchor routes through `incident_edge`.
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, ()> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1), n(2)], ())]);
        assert_eq!(rs.find_by_participants(&[EdgeId(0)], &ligands), expected);
    }

    #[rstest]
    #[case::exact(vec![n(0), n(1)], vec![n(2), n(3)], Some(RelationId(0)))]
    #[case::reordered(vec![n(1), n(0)], vec![n(3), n(2)], Some(RelationId(0)))]
    #[case::role_swap(vec![n(2), n(3)], vec![n(0), n(1)], None)]
    #[case::absent(vec![n(0), n(1)], vec![n(2), n(9)], None)]
    fn test_var_var_birelation_set_find_by_participants(
        #[case] query_1: Vec<NodeId>,
        #[case] query_2: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        let rs: VarVarBirelationSet<NodeId, Unordered, NodeId, Unordered, ()> =
            VarVarBirelationSet::new(vec![
                (vec![n(0), n(1)], vec![n(2), n(3)], ()),
                (vec![n(4)], vec![n(5)], ()),
            ]);
        assert_eq!(rs.find_by_participants(&query_1, &query_2), expected);
    }

    #[rstest]
    fn test_fixed_relation_set_pushout() {
        // same-space glue: self {01}=10 {23}=20 ; right {01}=5 (coincides) {45}=30 (new); combine=sum.
        let left = FixedRelationSet::<NodeId, Unordered, i32, 2>::new(vec![
            ([n(0), n(1)], 10),
            ([n(2), n(3)], 20),
        ]);
        let right = FixedRelationSet::<NodeId, Unordered, i32, 2>::new(vec![
            ([n(0), n(1)], 5),
            ([n(4), n(5)], 30),
        ]);
        let glue = left.pushout(&right, |a, b| Some(a + b)).expect("no ⊥");
        assert_eq!(glue.object.relation_count(), 3);
        let value = |q: [NodeId; 2]| {
            glue.object
                .find_by_participants(&q)
                .map(|id| *glue.object.data(id))
        };
        assert_eq!(value([n(0), n(1)]), Some(15)); // coincidence combined
        assert_eq!(value([n(2), n(3)]), Some(20)); // self-only carried
        assert_eq!(value([n(4), n(5)]), Some(30)); // right-only appended
        assert_eq!(glue.left.right_of(RelationId(0)), Some(RelationId(0))); // self identity
        assert_eq!(glue.right.right_of(RelationId(0)), Some(RelationId(0))); // right {01} folds onto self
        assert_eq!(glue.right.right_of(RelationId(1)), Some(RelationId(2))); // right {45} appended
    }

    #[rstest]
    fn test_fixed_relation_set_pushout_bottom() {
        // combine returns ⊥ on the coincidence → the whole glue is inadmissible.
        let left = FixedRelationSet::<NodeId, Unordered, i32, 2>::new(vec![([n(0), n(1)], 10)]);
        let right = FixedRelationSet::<NodeId, Unordered, i32, 2>::new(vec![([n(0), n(1)], 5)]);
        assert_eq!(left.pushout(&right, |_, _| None), None);
    }

    #[rstest]
    fn test_var_relation_set_pushout() {
        let left =
            VarRelationSet::<NodeId, Unordered, i32>::new(vec![(vec![n(0), n(1), n(2)], 10)]);
        let right = VarRelationSet::<NodeId, Unordered, i32>::new(vec![
            (vec![n(0), n(1), n(2)], 5),
            (vec![n(3), n(4)], 20),
        ]);
        let glue = left.pushout(&right, |a, b| Some(a + b)).expect("no ⊥");
        assert_eq!(glue.object.relation_count(), 2);
        assert_eq!(
            glue.object
                .find_by_participants(&[n(0), n(1), n(2)])
                .map(|id| *glue.object.data(id)),
            Some(15),
        );
        assert_eq!(glue.right.right_of(RelationId(1)), Some(RelationId(1))); // {34} appended
    }

    #[rstest]
    fn test_fixed_var_birelation_set_pushout() {
        // coincidence requires *both* factors equal (site + members).
        let left =
            FixedVarBirelationSet::<NodeId, Ordered, 1, NodeId, Unordered, i32>::new(vec![(
                [n(0)],
                vec![n(1), n(2)],
                10,
            )]);
        let right = FixedVarBirelationSet::<NodeId, Ordered, 1, NodeId, Unordered, i32>::new(vec![
            ([n(0)], vec![n(1), n(2)], 5),
            ([n(3)], vec![n(4), n(5)], 20),
        ]);
        let glue = left.pushout(&right, |a, b| Some(a + b)).expect("no ⊥");
        assert_eq!(glue.object.relation_count(), 2);
        assert_eq!(
            glue.object
                .find_by_participants(&[n(0)], &[n(1), n(2)])
                .map(|id| *glue.object.data(id)),
            Some(15),
        );
        assert_eq!(glue.right.right_of(RelationId(1)), Some(RelationId(1))); // appended
    }
}
