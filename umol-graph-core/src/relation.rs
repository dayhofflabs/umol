//! Relation sets: N-ary relations over typed participants (`NodeId`, `EdgeId`,
//! or external type implementing `RelationParticipant`), each carrying a shared
//! union incidence index (a node index and an edge index) routed from every
//! participant's `refs()`.
//! `FixedRelationSet<P, D, N>` stores relations of compile-time-known arity,
//! `VarRelationSet<P, D>` stores variable-arity relations. Participants are
//! typed `P` (`RelationParticipant`); the factor ordering `O` (`Unordered`/`Ordered`)
//! controls canonicalization. `FixedFixedBirelationSet`, `FixedVarBirelationSet`,
//! and `VarVarBirelationSet` relate two factors, each with its own participant
//! type, ordering, and arity. The union incidence spans both factors, so a relation
//! is reachable from any of its participants regardless of id-space.

use std::hash::{Hash, Hasher};
use std::ops::{Add, Sub};

use crate::compact::{Compaction, GraphCompaction};
use crate::correspondence::{Correspondence, GraphCorrespondence};
use crate::graph::{EdgeId, NodeId};
use crate::remap::GraphRemapping;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RelationId(pub u32);

impl RelationId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<RelationId> for usize {
    fn from(id: RelationId) -> Self {
        id.0 as usize
    }
}

impl From<usize> for RelationId {
    fn from(index: usize) -> Self {
        Self(index as u32)
    }
}

impl Add<usize> for RelationId {
    type Output = Self;

    fn add(self, offset: usize) -> Self {
        Self(self.0 + offset as u32)
    }
}

impl Sub<usize> for RelationId {
    type Output = Self;

    fn sub(self, offset: usize) -> Self {
        Self(self.0 - offset as u32)
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

/// The shared-relation intersection and its two projections — the result of a same-space relation
/// [`pullback`](FixedRelationSet::pullback).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelationPullback<S> {
    pub object: S,
    /// object relation → `self` relation.
    pub left: Correspondence<RelationId>,
    /// object relation → `right` relation.
    pub right: Correspondence<RelationId>,
}

/// The two projections of a relation pullback, each mapping a shared relation to its original.
fn relation_pullback<S>(
    object: S,
    left_images: Vec<RelationId>,
    right_images: Vec<RelationId>,
    self_count: usize,
    right_count: usize,
) -> RelationPullback<S> {
    RelationPullback {
        object,
        left: Correspondence::from_images(&left_images, self_count),
        right: Correspondence::from_images(&right_images, right_count),
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

/// Reorder `participants` in place so that `new[i] = old[order[i]]`.
///
/// Panics unless `order` is a permutation of `0..participants.len()`. Under that condition the
/// participant multiset is unchanged, so an incidence index built over these participants stays
/// valid and is left untouched.
fn permute_participants<P: Copy>(participants: &mut [P], order: &[ParticipantPosition]) {
    assert_eq!(
        order.len(),
        participants.len(),
        "permute: order length must equal the relation's arity"
    );
    let mut seen = vec![false; participants.len()];
    for position in order {
        let index = position.index();
        assert!(
            index < participants.len(),
            "permute: position {index} is outside 0..{}",
            participants.len()
        );
        assert!(!seen[index], "permute: position {index} is repeated");
        seen[index] = true;
    }
    let permuted: Vec<P> = order
        .iter()
        .map(|position| participants[position.index()])
        .collect();
    participants.copy_from_slice(&permuted);
}

/// The id-space contents of a participant, surfaced for the incidence index.
/// At most one ref per space today (a node or an edge); a future port type
/// could fill both.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ParticipantRefs {
    pub node: Option<NodeId>,
    pub edge: Option<EdgeId>,
}
/// A value that can occupy a relation factor: supports compaction, remapping, and
/// correspondence-based id transport, and exposes its node/edge refs for incidence. One impl per concrete
/// id type — dispatch is static, since a factor is homogeneous.
pub trait RelationParticipant: Copy + Ord + Hash {
    fn compact(self, compaction: &GraphCompaction) -> Option<Self>;
    fn uncompact(self, compaction: &GraphCompaction) -> Self;

    /// Relabel every referenced id through a correspondence, preserving other participant data.
    ///
    /// # Panics
    /// Panics when a referenced node or edge has no image.
    fn map(self, correspondence: &GraphCorrespondence) -> Self {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Relabel every referenced id, or return `None` if any reference has no image.
    /// Unused correspondence entries need not be matched. Every id looked up must be
    /// reported by [`refs`](Self::refs); all other participant data must be preserved.
    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self>;

    /// Relabel this participant through `remapping`.
    ///
    /// Every node or edge id read from `remapping` must be reported by [`refs`](Self::refs), so
    /// checked relation-set remapping can establish coverage before calling this method.
    fn remap(self, remapping: &GraphRemapping) -> Self;

    /// Return every graph id used to represent this participant.
    fn refs(self) -> ParticipantRefs;
}

impl RelationParticipant for NodeId {
    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self> {
        correspondence.nodes().right_of(self)
    }
    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        compaction.compact_node(self)
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
        compaction.uncompact_node(self)
    }

    fn remap(self, remapping: &GraphRemapping) -> Self {
        remapping.map_node(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: Some(self),
            edge: None,
        }
    }
}

impl RelationParticipant for EdgeId {
    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self> {
        correspondence.edges().right_of(self)
    }
    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        compaction.compact_edge(self)
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
        compaction.uncompact_edge(self)
    }

    fn remap(self, remapping: &GraphRemapping) -> Self {
        remapping.map_edge(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: None,
            edge: Some(self),
        }
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
pub struct FixedRelationSet<P, D, const N: usize> {
    participants: Vec<[P; N]>,
    data: Vec<D>,
    incidence: Incidence,
}

impl<P: PartialEq, D: PartialEq, const N: usize> PartialEq for FixedRelationSet<P, D, N> {
    fn eq(&self, other: &Self) -> bool {
        self.participants == other.participants && self.data == other.data
    }
}

impl<P: Eq, D: Eq, const N: usize> Eq for FixedRelationSet<P, D, N> {}

impl<P: Hash, D: Hash, const N: usize> Hash for FixedRelationSet<P, D, N> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.participants.hash(state);
        self.data.hash(state);
    }
}

/// Whether every graph id these participants reference has an image under `remapping` — the
/// precondition [`FixedRelationSet::try_remap`] and its peers check before relabelling.
fn remappable_under<P>(participants: &[P], remapping: &GraphRemapping) -> bool
where
    P: RelationParticipant,
{
    participants.iter().all(|participant| {
        let refs = participant.refs();
        refs.node
            .is_none_or(|node| remapping.try_map_node(node).is_some())
            && refs
                .edge
                .is_none_or(|edge| remapping.try_map_edge(edge).is_some())
    })
}

/// Multiset equality of a relation's stored participants against `query`, which the caller has
/// already sorted (hoisted out of its candidate scan). The stored frame is left intact.
///
/// Matches on identity — the participant multiset — independent of the factor's ordering marker.
fn participants_match<P: RelationParticipant>(participants: &[P], query: &[P]) -> bool {
    if participants.len() != query.len() {
        return false;
    }
    let mut sorted: Vec<P> = participants.to_vec();
    sorted.sort_unstable();
    sorted == query
}

impl<P: RelationParticipant, D, const N: usize> FixedRelationSet<P, D, N> {
    pub fn new(entries: Vec<([P; N], D)>) -> Self {
        let mut participants = Vec::with_capacity(entries.len());
        let mut data = Vec::with_capacity(entries.len());
        for (p, d) in entries {
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
        }
    }

    /// Consume the set into its canonical stored entries, in relation-id order.
    pub fn into_entries(self) -> Vec<([P; N], D)> {
        self.participants.into_iter().zip(self.data).collect()
    }

    pub fn count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    /// Every relation as `(id, participants, payload)` in relation-id order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (RelationId, &[P; N], &D)> {
        self.participants
            .iter()
            .zip(&self.data)
            .enumerate()
            .map(|(index, (participants, data))| (RelationId(index as u32), participants, data))
    }

    /// Every relation as `(id, participants, payload)` in relation-id order, the payload mutable.
    ///
    /// Participants stay immutable: changing them would invalidate the incidence index, which
    /// [`permute_with`](Self::permute_with) is the one operation allowed to leave intact.
    pub fn iter_mut(&mut self) -> impl ExactSizeIterator<Item = (RelationId, &[P; N], &mut D)> {
        let participants = &self.participants;
        self.data
            .iter_mut()
            .enumerate()
            .map(move |(index, data)| (RelationId(index as u32), &participants[index], data))
    }

    pub fn participants(&self, id: RelationId) -> &[P; N] {
        &self.participants[id.index()]
    }

    /// Reorder relation `id`'s participants so that `new[i] = old[order[i]]`, leaving the payload
    /// untouched.
    ///
    /// The multiset is unchanged, so incidence answers identically and is not rebuilt. Panics
    /// unless `order` is a permutation of `0..arity`.
    pub fn permute_with(&mut self, id: RelationId, order: &[ParticipantPosition]) {
        permute_participants(self.participants[id.index()].as_mut_slice(), order);
    }

    /// Id of the relation coinciding with `query` — the one whose participants equal it as a
    /// multiset, in any order.
    ///
    /// This is the identity question, not a lookup: the participant multiset is the relation's
    /// identity and the stored sequence is only the frame its payload is expressed in, so two
    /// entries presenting the same participants differently coincide. It is what
    /// [`pushout`](Self::pushout) and [`pullback`](Self::pullback) join on. Naming an entity by a
    /// subset of its constituents is a different question with a different key, and belongs to the
    /// caller that knows the key. §4.1 uniqueness ⇒ at most one hit.
    pub fn coincident(&self, node: NodeId, query: &[P]) -> Option<RelationId> {
        self.coincident_in(self.incident(node), query)
    }

    /// This is the identity question, not a lookup: the participant multiset is the relation's
    /// identity and the stored sequence is only the frame its payload is expressed in, so two
    /// entries presenting the same participants differently coincide. It is what
    /// [`pushout`](Self::pushout) and [`pullback`](Self::pullback) join on. Naming an entity by a
    /// subset of its constituents is a different question with a different key, and belongs to the
    /// caller that knows the key. §4.1 uniqueness ⇒ at most one hit.
    ///
    /// `edge` narrows the scan to the edge incidence index. The node-indexed peer is
    /// [`coincident`](Self::coincident).
    pub fn coincident_edge(&self, edge: EdgeId, query: &[P]) -> Option<RelationId> {
        self.coincident_in(self.incident_edge(edge), query)
    }

    /// Whether relation `id` coincides with `query` — the known-id sibling of
    /// [`coincident`](Self::coincident), which searches for it instead.
    ///
    /// `pushout` and `pullback` apply this to a supplied pairing before gluing on it. A caller that
    /// already holds the id and needs identity established — because a frame-invariant payload
    /// carries without reading either frame — asks here rather than deriving the comparison again.
    pub fn is_coincident(&self, id: RelationId, query: &[P]) -> bool {
        self.coincident_in(&[id], query).is_some()
    }

    fn coincident_in(&self, candidates: &[RelationId], query: &[P]) -> Option<RelationId> {
        let mut sorted_query: Vec<P> = query.to_vec();
        sorted_query.sort_unstable();
        candidates
            .iter()
            .copied()
            .find(|&id| participants_match(self.participants(id), &sorted_query))
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

    pub fn ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    /// Compact participant ids and drop relations containing a removed participant.
    ///
    /// Produces the same set as [`Self::tracked_compact`], without returning
    /// the relation-id compaction.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: Clone,
    {
        self.tracked_compact(compaction).0
    }

    /// Compact participant ids, dropping every relation that contains a removed participant, and
    /// report which relation ids the drop consumed.
    ///
    /// The returned compaction moves this set's own ids, so a caller holding relation ids can
    /// carry them across the removal without a second traversal.
    pub fn tracked_compact(&self, compaction: &GraphCompaction) -> (Self, Compaction<RelationId>)
    where
        D: Clone,
    {
        let mut removed = Vec::new();
        let mut entries: Vec<([P; N], D)> = Vec::with_capacity(self.count());
        for i in 0..self.count() {
            let rid = RelationId(i as u32);
            let parts: Option<[P; N]> = self
                .participants(rid)
                .iter()
                .map(|&p| p.compact(compaction))
                .collect::<Option<Vec<P>>>()
                .and_then(|parts| parts.try_into().ok());
            match parts {
                Some(parts) => entries.push((parts, self.data(rid).clone())),
                None => removed.push(rid),
            }
        }
        (
            Self::new(entries),
            Compaction::new(self.count(), removed)
                .expect("removed relations belong to the source set"),
        )
    }

    /// Relabel every participant, preserving rows, participant order, and payloads.
    ///
    /// # Semantic properties
    ///
    /// Each positional payload item remains attached to the participant whose id is relabeled.
    ///
    /// # Panics
    ///
    /// Panics when a participant lies outside the remapping's corresponding source range.
    pub fn remap(&self, remapping: &GraphRemapping) -> Self
    where
        D: Clone,
    {
        self.map_participants(|participant| Some(participant.remap(remapping)))
            .expect("remapping transport supplies every participant")
    }

    /// Relabel participant ids through a correspondence without changing rows, frames, or payloads.
    ///
    /// # Panics
    /// Panics when any referenced node or edge has no image.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self
    where
        D: Clone,
    {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Relabel participant ids, returning `None` when any reference has no image.
    ///
    /// Only referenced ids require images; unrelated source entries may be unmatched.
    ///
    /// # Semantic properties
    /// Row ids, participant positions, and payloads are preserved. Identity mapping is exact;
    /// sequential covered mappings agree with their correspondence composition.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self>
    where
        D: Clone,
    {
        self.map_participants(|participant| participant.try_map(correspondence))
    }

    fn map_participants(&self, mut map_1: impl FnMut(P) -> Option<P>) -> Option<Self>
    where
        D: Clone,
    {
        let entries = self
            .ids()
            .map(|id| {
                let parts_1: [P; N] = self
                    .participants(id)
                    .iter()
                    .copied()
                    .map(&mut map_1)
                    .collect::<Option<Vec<_>>>()?
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                Some((parts_1, self.data(id).clone()))
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Self::new(entries))
    }

    /// Relabel every participant, returning `None` when the remapping does not cover the set.
    pub fn try_remap(&self, remapping: &GraphRemapping) -> Option<Self>
    where
        D: Clone,
    {
        self.ids()
            .all(|id| remappable_under(self.participants(id), remapping))
            .then(|| self.remap(remapping))
    }

    /// Glue `self` and `right`, both **already in the same participant id-space**, identifying
    /// coinciding relations (equal participants) — the same-space relation pushout. `combine` merges
    /// the data of a coincidence (`None` = ⊥ ⇒ the whole glue is inadmissible ⇒ `None`); every other
    /// relation is carried. `self`'s ids are the identity prefix of the object, `right`'s
    /// non-coinciding relations are appended. The caller brings both sides, including positional
    /// data, into the common space with [`map`](Self::map) or [`remap`](Self::remap) first.
    pub fn pushout(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[P]) -> Option<RelationId>,
        mut combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<([P; N], D)> = self
            .ids()
            .map(|id| (*self.participants(id), self.data(id).clone()))
            .collect();
        let self_count = entries.len();
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.count());
        for id in right.ids() {
            match coincident(self, right.participants(id)) {
                Some(hit) => {
                    let merged = combine(
                        (self.participants(hit), self.data(hit)),
                        (right.participants(id), right.data(id)),
                    )?;
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

    /// Same-space relation pullback — the shared relations (coinciding participants), data
    /// combined by `combine` (`None` = ⊥ ⇒ inadmissible); non-coinciding relations are dropped.
    /// Its two projections map each shared relation to its `self` / `right` original. Same-space
    /// contract as [`FixedRelationSet::pushout`].
    pub fn pullback(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[P]) -> Option<RelationId>,
        mut combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
    ) -> Option<RelationPullback<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<([P; N], D)> = Vec::new();
        let mut left_images: Vec<RelationId> = Vec::new();
        let mut right_images: Vec<RelationId> = Vec::new();
        for id in self.ids() {
            if let Some(hit) = coincident(right, self.participants(id)) {
                let merged = combine(
                    (self.participants(id), self.data(id)),
                    (right.participants(hit), right.data(hit)),
                )?;
                entries.push((*self.participants(id), merged));
                left_images.push(id);
                right_images.push(hit);
            }
        }
        Some(relation_pullback(
            Self::new(entries),
            left_images,
            right_images,
            self.count(),
            right.count(),
        ))
    }
}

impl<P, D, const N: usize> Default for FixedRelationSet<P, D, N> {
    fn default() -> Self {
        Self {
            participants: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
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
pub struct VarRelationSet<P, D> {
    offsets: Vec<u32>,
    participants: Vec<P>,
    data: Vec<D>,
    incidence: Incidence,
}

impl<P: PartialEq, D: PartialEq> PartialEq for VarRelationSet<P, D> {
    fn eq(&self, other: &Self) -> bool {
        self.offsets == other.offsets
            && self.participants == other.participants
            && self.data == other.data
    }
}

impl<P: Eq, D: Eq> Eq for VarRelationSet<P, D> {}

impl<P: Hash, D: Hash> Hash for VarRelationSet<P, D> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.offsets.hash(state);
        self.participants.hash(state);
        self.data.hash(state);
    }
}

impl<P: RelationParticipant, D> VarRelationSet<P, D> {
    pub fn new(entries: Vec<(Vec<P>, D)>) -> Self
    where
        D: Clone,
    {
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
        }
    }

    /// Consume the set into its canonical stored entries, in relation-id order.
    pub fn into_entries(self) -> Vec<(Vec<P>, D)> {
        let Self {
            offsets,
            participants,
            data,
            ..
        } = self;
        let lengths = offsets
            .windows(2)
            .map(|range| (range[1] - range[0]) as usize);
        let mut participants = participants.into_iter();
        lengths
            .zip(data)
            .map(|(len, data)| (participants.by_ref().take(len).collect(), data))
            .collect()
    }

    pub fn count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    /// Every relation as `(id, participants, payload)` in relation-id order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (RelationId, &[P], &D)> {
        let offsets = &self.offsets;
        let participants = &self.participants;
        self.data.iter().enumerate().map(move |(index, data)| {
            let start = offsets[index] as usize;
            let end = offsets[index + 1] as usize;
            (RelationId(index as u32), &participants[start..end], data)
        })
    }

    /// Every relation as `(id, participants, payload)` in relation-id order, the payload mutable.
    ///
    /// Participants stay immutable: changing them would invalidate the incidence index, which
    /// [`permute_with`](Self::permute_with) is the one operation allowed to leave intact.
    pub fn iter_mut(&mut self) -> impl ExactSizeIterator<Item = (RelationId, &[P], &mut D)> {
        let offsets = &self.offsets;
        let participants = &self.participants;
        self.data.iter_mut().enumerate().map(move |(index, data)| {
            let start = offsets[index] as usize;
            let end = offsets[index + 1] as usize;
            (RelationId(index as u32), &participants[start..end], data)
        })
    }

    pub fn participants(&self, id: RelationId) -> &[P] {
        let start = self.offsets[id.index()] as usize;
        let end = self.offsets[id.index() + 1] as usize;
        &self.participants[start..end]
    }

    /// Reorder relation `id`'s participants so that `new[i] = old[order[i]]`, leaving the payload
    /// untouched.
    ///
    /// The multiset is unchanged, so incidence answers identically and is not rebuilt. Panics
    /// unless `order` is a permutation of `0..arity`.
    pub fn permute_with(&mut self, id: RelationId, order: &[ParticipantPosition]) {
        let start = self.offsets[id.index()] as usize;
        let end = self.offsets[id.index() + 1] as usize;
        permute_participants(&mut self.participants[start..end], order);
    }

    /// Id of the relation coinciding with `query` — the one whose participants equal it as a
    /// multiset, in any order.
    ///
    /// This is the identity question, not a lookup: the participant multiset is the relation's
    /// identity and the stored sequence is only the frame its payload is expressed in, so two
    /// entries presenting the same participants differently coincide. It is what
    /// [`pushout`](Self::pushout) and [`pullback`](Self::pullback) join on. Naming an entity by a
    /// subset of its constituents is a different question with a different key, and belongs to the
    /// caller that knows the key. §4.1 uniqueness ⇒ at most one hit.
    pub fn coincident(&self, node: NodeId, query: &[P]) -> Option<RelationId> {
        self.coincident_in(self.incident(node), query)
    }

    /// This is the identity question, not a lookup: the participant multiset is the relation's
    /// identity and the stored sequence is only the frame its payload is expressed in, so two
    /// entries presenting the same participants differently coincide. It is what
    /// [`pushout`](Self::pushout) and [`pullback`](Self::pullback) join on. Naming an entity by a
    /// subset of its constituents is a different question with a different key, and belongs to the
    /// caller that knows the key. §4.1 uniqueness ⇒ at most one hit.
    ///
    /// `edge` narrows the scan to the edge incidence index. The node-indexed peer is
    /// [`coincident`](Self::coincident).
    pub fn coincident_edge(&self, edge: EdgeId, query: &[P]) -> Option<RelationId> {
        self.coincident_in(self.incident_edge(edge), query)
    }

    /// Whether relation `id` coincides with these participants — the known-id sibling of
    /// [`coincident`](Self::coincident), which searches for it instead.
    ///
    /// `pushout` and `pullback` apply this to a supplied pairing before gluing on it. A caller that
    /// already holds the id and needs identity established — because a frame-invariant payload
    /// carries without reading either frame — asks here rather than deriving the comparison again.
    pub fn is_coincident(&self, id: RelationId, query: &[P]) -> bool {
        self.coincident_in(&[id], query).is_some()
    }

    fn coincident_in(&self, candidates: &[RelationId], query: &[P]) -> Option<RelationId> {
        let mut sorted_query: Vec<P> = query.to_vec();
        sorted_query.sort_unstable();
        candidates
            .iter()
            .copied()
            .find(|&id| participants_match(self.participants(id), &sorted_query))
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

    pub fn ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    /// Compact participant ids and drop relations containing a removed participant.
    ///
    /// Produces the same set as [`Self::tracked_compact`], without returning
    /// the relation-id compaction.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: Clone,
    {
        self.tracked_compact(compaction).0
    }

    /// Compact participant ids, dropping every relation that contains a removed participant, and
    /// report which relation ids the drop consumed.
    ///
    /// The returned compaction moves this set's own ids, so a caller holding relation ids can
    /// carry them across the removal without a second traversal.
    pub fn tracked_compact(&self, compaction: &GraphCompaction) -> (Self, Compaction<RelationId>)
    where
        D: Clone,
    {
        let mut removed = Vec::new();
        let mut entries: Vec<(Vec<P>, D)> = Vec::with_capacity(self.count());
        for i in 0..self.count() {
            let rid = RelationId(i as u32);
            let parts: Option<Vec<P>> = self
                .participants(rid)
                .iter()
                .map(|&p| p.compact(compaction))
                .collect();
            match parts {
                Some(parts) => entries.push((parts, self.data(rid).clone())),
                None => removed.push(rid),
            }
        }
        (
            Self::new(entries),
            Compaction::new(self.count(), removed)
                .expect("removed relations belong to the source set"),
        )
    }

    /// Relabel every participant, preserving rows, participant order, and payloads.
    ///
    /// # Semantic properties
    ///
    /// Each positional payload item remains attached to the participant whose id is relabeled.
    ///
    /// # Panics
    ///
    /// Panics when a participant lies outside the remapping's corresponding source range.
    pub fn remap(&self, remapping: &GraphRemapping) -> Self
    where
        D: Clone,
    {
        self.map_participants(|participant| Some(participant.remap(remapping)))
            .expect("remapping transport supplies every participant")
    }

    /// Relabel participant ids through a correspondence without changing rows, frames, or payloads.
    ///
    /// # Panics
    /// Panics when any referenced node or edge has no image.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self
    where
        D: Clone,
    {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Relabel participant ids, returning `None` when any reference has no image.
    ///
    /// Only referenced ids require images; unrelated source entries may be unmatched.
    ///
    /// # Semantic properties
    /// Row ids, participant positions, and payloads are preserved. Identity mapping is exact;
    /// sequential covered mappings agree with their correspondence composition.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self>
    where
        D: Clone,
    {
        self.map_participants(|participant| participant.try_map(correspondence))
    }

    fn map_participants(&self, mut map_1: impl FnMut(P) -> Option<P>) -> Option<Self>
    where
        D: Clone,
    {
        let entries = self
            .ids()
            .map(|id| {
                let parts_1: Vec<P> = self
                    .participants(id)
                    .iter()
                    .copied()
                    .map(&mut map_1)
                    .collect::<Option<Vec<_>>>()?;
                Some((parts_1, self.data(id).clone()))
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Self::new(entries))
    }

    /// Relabel every participant, returning `None` when the remapping does not cover the set.
    pub fn try_remap(&self, remapping: &GraphRemapping) -> Option<Self>
    where
        D: Clone,
    {
        self.ids()
            .all(|id| remappable_under(self.participants(id), remapping))
            .then(|| self.remap(remapping))
    }

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`].
    pub fn pushout(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[P]) -> Option<RelationId>,
        mut combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<(Vec<P>, D)> = self
            .ids()
            .map(|id| (self.participants(id).to_vec(), self.data(id).clone()))
            .collect();
        let self_count = entries.len();
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.count());
        for id in right.ids() {
            match coincident(self, right.participants(id)) {
                Some(hit) => {
                    let merged = combine(
                        (self.participants(hit), self.data(hit)),
                        (right.participants(id), right.data(id)),
                    )?;
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

    /// Same-space relation pullback — see [`FixedRelationSet::pullback`].
    pub fn pullback(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[P]) -> Option<RelationId>,
        mut combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
    ) -> Option<RelationPullback<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<(Vec<P>, D)> = Vec::new();
        let mut left_images: Vec<RelationId> = Vec::new();
        let mut right_images: Vec<RelationId> = Vec::new();
        for id in self.ids() {
            if let Some(hit) = coincident(right, self.participants(id)) {
                let merged = combine(
                    (self.participants(id), self.data(id)),
                    (right.participants(hit), right.data(hit)),
                )?;
                entries.push((self.participants(id).to_vec(), merged));
                left_images.push(id);
                right_images.push(hit);
            }
        }
        Some(relation_pullback(
            Self::new(entries),
            left_images,
            right_images,
            self.count(),
            right.count(),
        ))
    }
}

impl<P, D> Default for VarRelationSet<P, D> {
    fn default() -> Self {
        Self {
            offsets: vec![0],
            participants: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
        }
    }
}

/// Birelation with two fixed-arity factors.
#[derive(Clone, Debug)]
pub struct FixedFixedBirelationSet<L1, const N1: usize, L2, const N2: usize, D> {
    participants_1: Vec<[L1; N1]>,
    participants_2: Vec<[L2; N2]>,
    data: Vec<D>,
    incidence: Incidence,
}

impl<L1, const N1: usize, L2, const N2: usize, D> PartialEq
    for FixedFixedBirelationSet<L1, N1, L2, N2, D>
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

impl<L1, const N1: usize, L2, const N2: usize, D> Eq for FixedFixedBirelationSet<L1, N1, L2, N2, D>
where
    L1: Eq,
    L2: Eq,
    D: Eq,
{
}

impl<L1, const N1: usize, L2, const N2: usize, D> FixedFixedBirelationSet<L1, N1, L2, N2, D>
where
    L1: RelationParticipant,
    L2: RelationParticipant,
{
    pub fn new(entries: Vec<([L1; N1], [L2; N2], D)>) -> Self
    where
        D: Clone,
    {
        let relation_count = entries.len();
        let mut participants_1 = Vec::with_capacity(relation_count);
        let mut participants_2 = Vec::with_capacity(relation_count);
        let mut data = Vec::with_capacity(relation_count);
        for (l1, l2, d) in entries {
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
        }
    }

    /// Consume the set into its canonical stored entries, in relation-id order.
    pub fn into_entries(self) -> Vec<([L1; N1], [L2; N2], D)> {
        let Self {
            participants_1,
            participants_2,
            data,
            ..
        } = self;
        participants_1
            .into_iter()
            .zip(participants_2)
            .zip(data)
            .map(|((participants_1, participants_2), data)| (participants_1, participants_2, data))
            .collect()
    }

    pub fn count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    /// Every relation as `(id, first-factor participants, second-factor participants, payload)` in
    /// relation-id order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (RelationId, &[L1; N1], &[L2; N2], &D)> {
        let participants_1 = &self.participants_1;
        let participants_2 = &self.participants_2;

        self.data.iter().enumerate().map(move |(index, data)| {
            (
                RelationId(index as u32),
                &participants_1[index],
                &participants_2[index],
                data,
            )
        })
    }

    /// Every relation as `(id, first-factor participants, second-factor participants, payload)` in
    /// relation-id order, the payload mutable.
    ///
    /// Participants stay immutable: changing them would invalidate the incidence index, which
    /// [`permute_1_with`](Self::permute_1_with) and [`permute_2_with`](Self::permute_2_with) are
    /// the operations allowed to leave intact.
    pub fn iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = (RelationId, &[L1; N1], &[L2; N2], &mut D)> {
        let participants_1 = &self.participants_1;
        let participants_2 = &self.participants_2;

        self.data.iter_mut().enumerate().map(move |(index, data)| {
            (
                RelationId(index as u32),
                &participants_1[index],
                &participants_2[index],
                data,
            )
        })
    }

    pub fn participants_1(&self, id: RelationId) -> &[L1; N1] {
        &self.participants_1[id.index()]
    }

    pub fn participants_2(&self, id: RelationId) -> &[L2; N2] {
        &self.participants_2[id.index()]
    }

    /// Reorder relation `id`'s first-factor participants so that `new[i] = old[order[i]]`,
    /// leaving the second factor and the payload untouched.
    ///
    /// The multiset is unchanged, so incidence answers identically and is not rebuilt. Panics
    /// unless `order` is a permutation of `0..arity`.
    pub fn permute_1_with(&mut self, id: RelationId, order: &[ParticipantPosition]) {
        permute_participants(self.participants_1[id.index()].as_mut_slice(), order);
    }

    /// Reorder relation `id`'s second-factor participants so that `new[i] = old[order[i]]`,
    /// leaving the first factor and the payload untouched.
    ///
    /// The multiset is unchanged, so incidence answers identically and is not rebuilt. Panics
    /// unless `order` is a permutation of `0..arity`.
    pub fn permute_2_with(&mut self, id: RelationId, order: &[ParticipantPosition]) {
        permute_participants(self.participants_2[id.index()].as_mut_slice(), order);
    }

    /// Id of the relation coinciding with `query_1` / `query_2` — the one whose factors equal them
    /// as multisets, in any order.
    ///
    /// This is the identity question, not a lookup: the participant multiset is the relation's
    /// identity and the stored sequence is only the frame its payload is expressed in, so two
    /// entries presenting the same participants differently coincide. It is what
    /// [`pushout`](Self::pushout) and [`pullback`](Self::pullback) join on. Naming an entity by a
    /// subset of its constituents is a different question with a different key, and belongs to the
    /// caller that knows the key. §4.1 uniqueness ⇒ at most one hit.
    pub fn coincident(&self, node: NodeId, query_1: &[L1], query_2: &[L2]) -> Option<RelationId> {
        self.coincident_in(self.incident(node), query_1, query_2)
    }

    /// Edge-indexed peer of [`coincident`](Self::coincident).
    pub fn coincident_edge(
        &self,
        edge: EdgeId,
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<RelationId> {
        self.coincident_in(self.incident_edge(edge), query_1, query_2)
    }

    /// Whether relation `id` coincides with `query_1` / `query_2` — the known-id sibling of
    /// [`coincident`](Self::coincident), which searches for it instead.
    ///
    /// `pushout` and `pullback` apply this to a supplied pairing before gluing on it. A caller that
    /// already holds the id and needs identity established — because a frame-invariant payload
    /// carries without reading either frame — asks here rather than deriving the comparison again.
    pub fn is_coincident(&self, id: RelationId, query_1: &[L1], query_2: &[L2]) -> bool {
        self.coincident_in(&[id], query_1, query_2).is_some()
    }

    fn coincident_in(
        &self,
        candidates: &[RelationId],
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<RelationId> {
        let mut sorted_1: Vec<L1> = query_1.to_vec();
        sorted_1.sort_unstable();
        let mut sorted_2: Vec<L2> = query_2.to_vec();
        sorted_2.sort_unstable();
        candidates.iter().copied().find(|&id| {
            participants_match(self.participants_1(id), &sorted_1)
                && participants_match(self.participants_2(id), &sorted_2)
        })
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

    pub fn ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    /// Compact participant ids and drop relations containing a removed participant.
    ///
    /// Produces the same set as [`Self::tracked_compact`], without returning
    /// the relation-id compaction.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: Clone,
    {
        self.tracked_compact(compaction).0
    }

    /// Compact participant ids, dropping every relation that contains a removed participant, and
    /// report which relation ids the drop consumed.
    ///
    /// The returned compaction moves this set's own ids, so a caller holding relation ids can
    /// carry them across the removal without a second traversal.
    pub fn tracked_compact(&self, compaction: &GraphCompaction) -> (Self, Compaction<RelationId>)
    where
        D: Clone,
    {
        let mut removed = Vec::new();
        let mut entries: Vec<([L1; N1], [L2; N2], D)> = Vec::with_capacity(self.count());
        for i in 0..self.count() {
            let rid = RelationId(i as u32);
            let f1: Option<[L1; N1]> = self
                .participants_1(rid)
                .iter()
                .map(|&p| p.compact(compaction))
                .collect::<Option<Vec<L1>>>()
                .and_then(|parts| parts.try_into().ok());
            let f2: Option<[L2; N2]> = self
                .participants_2(rid)
                .iter()
                .map(|&p| p.compact(compaction))
                .collect::<Option<Vec<L2>>>()
                .and_then(|parts| parts.try_into().ok());
            match (f1, f2) {
                (Some(f1), Some(f2)) => entries.push((f1, f2, self.data(rid).clone())),
                _ => removed.push(rid),
            }
        }
        (
            Self::new(entries),
            Compaction::new(self.count(), removed)
                .expect("removed relations belong to the source set"),
        )
    }

    /// Relabel every participant, preserving rows, participant order, and payloads.
    ///
    /// # Semantic properties
    ///
    /// In both factors, each positional payload item remains attached to the participant whose id
    /// is relabeled.
    ///
    /// # Panics
    ///
    /// Panics when a participant lies outside the remapping's corresponding source range.
    pub fn remap(&self, remapping: &GraphRemapping) -> Self
    where
        D: Clone,
    {
        self.map_participants(
            |participant| Some(participant.remap(remapping)),
            |participant| Some(participant.remap(remapping)),
        )
        .expect("remapping transport supplies every participant")
    }

    /// Relabel participant ids through a correspondence without changing rows, frames, or payloads.
    ///
    /// # Panics
    /// Panics when any referenced node or edge has no image.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self
    where
        D: Clone,
    {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Relabel participant ids, returning `None` when any reference has no image.
    ///
    /// Only referenced ids require images; unrelated source entries may be unmatched.
    ///
    /// # Semantic properties
    /// Row ids, participant positions, and payloads are preserved. Identity mapping is exact;
    /// sequential covered mappings agree with their correspondence composition.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self>
    where
        D: Clone,
    {
        self.map_participants(
            |participant| participant.try_map(correspondence),
            |participant| participant.try_map(correspondence),
        )
    }

    fn map_participants(
        &self,
        mut map_1: impl FnMut(L1) -> Option<L1>,
        mut map_2: impl FnMut(L2) -> Option<L2>,
    ) -> Option<Self>
    where
        D: Clone,
    {
        let entries = self
            .ids()
            .map(|id| {
                let parts_1: [L1; N1] = self
                    .participants_1(id)
                    .iter()
                    .copied()
                    .map(&mut map_1)
                    .collect::<Option<Vec<_>>>()?
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                let parts_2: [L2; N2] = self
                    .participants_2(id)
                    .iter()
                    .copied()
                    .map(&mut map_2)
                    .collect::<Option<Vec<_>>>()?
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                Some((parts_1, parts_2, self.data(id).clone()))
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Self::new(entries))
    }

    /// Relabel every participant, returning `None` when the remapping does not cover either factor.
    pub fn try_remap(&self, remapping: &GraphRemapping) -> Option<Self>
    where
        D: Clone,
    {
        self.ids()
            .all(|id| {
                remappable_under(self.participants_1(id), remapping)
                    && remappable_under(self.participants_2(id), remapping)
            })
            .then(|| self.remap(remapping))
    }

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`]. Coincidence is equality of
    /// both factors' participants.
    pub fn pushout(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        mut combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<([L1; N1], [L2; N2], D)> = self
            .ids()
            .map(|id| {
                (
                    *self.participants_1(id),
                    *self.participants_2(id),
                    self.data(id).clone(),
                )
            })
            .collect();
        let self_count = entries.len();
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.count());
        for id in right.ids() {
            match coincident(self, right.participants_1(id), right.participants_2(id)) {
                Some(hit) => {
                    let merged = combine(
                        (
                            self.participants_1(hit),
                            self.participants_2(hit),
                            self.data(hit),
                        ),
                        (
                            right.participants_1(id),
                            right.participants_2(id),
                            right.data(id),
                        ),
                    )?;
                    entries[hit.index()].2 = merged;
                    right_map.push(hit);
                }
                None => {
                    right_map.push(RelationId(entries.len() as u32));
                    entries.push((
                        *right.participants_1(id),
                        *right.participants_2(id),
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

    /// Same-space relation pullback — see [`FixedRelationSet::pullback`].
    pub fn pullback(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        mut combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<RelationPullback<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<([L1; N1], [L2; N2], D)> = Vec::new();
        let mut left_images: Vec<RelationId> = Vec::new();
        let mut right_images: Vec<RelationId> = Vec::new();
        for id in self.ids() {
            if let Some(hit) = coincident(right, self.participants_1(id), self.participants_2(id)) {
                let merged = combine(
                    (
                        self.participants_1(id),
                        self.participants_2(id),
                        self.data(id),
                    ),
                    (
                        right.participants_1(hit),
                        right.participants_2(hit),
                        right.data(hit),
                    ),
                )?;
                entries.push((*self.participants_1(id), *self.participants_2(id), merged));
                left_images.push(id);
                right_images.push(hit);
            }
        }
        Some(relation_pullback(
            Self::new(entries),
            left_images,
            right_images,
            self.count(),
            right.count(),
        ))
    }
}

impl<L1, const N1: usize, L2, const N2: usize, D> Default
    for FixedFixedBirelationSet<L1, N1, L2, N2, D>
{
    fn default() -> Self {
        Self {
            participants_1: Vec::new(),
            participants_2: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
        }
    }
}

/// Birelation with a fixed-arity factor 1 and a variable-arity factor 2. Each factor
/// has its own participant type and ordering; the union incidence spans both.
#[derive(Clone, Debug)]
pub struct FixedVarBirelationSet<L1, const N1: usize, L2, D> {
    participants_1: Vec<[L1; N1]>,
    f2_offsets: Vec<u32>,
    participants_2: Vec<L2>,
    data: Vec<D>,
    incidence: Incidence,
}

impl<L1, const N1: usize, L2, D> PartialEq for FixedVarBirelationSet<L1, N1, L2, D>
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

impl<L1, const N1: usize, L2, D> Eq for FixedVarBirelationSet<L1, N1, L2, D>
where
    L1: Eq,
    L2: Eq,
    D: Eq,
{
}

impl<L1, const N1: usize, L2, D> Hash for FixedVarBirelationSet<L1, N1, L2, D>
where
    L1: Hash,
    L2: Hash,
    D: Hash,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.participants_1.hash(state);
        self.f2_offsets.hash(state);
        self.participants_2.hash(state);
        self.data.hash(state);
    }
}

impl<L1, const N1: usize, L2, D> FixedVarBirelationSet<L1, N1, L2, D>
where
    L1: RelationParticipant,
    L2: RelationParticipant,
{
    pub fn new(entries: Vec<([L1; N1], Vec<L2>, D)>) -> Self
    where
        D: Clone,
    {
        let relation_count = entries.len();
        let mut participants_1 = Vec::with_capacity(relation_count);
        let mut f2_offsets = Vec::with_capacity(relation_count + 1);
        f2_offsets.push(0);
        let mut participants_2 = Vec::new();
        let mut data = Vec::with_capacity(relation_count);
        for (l1, l2, d) in entries {
            participants_1.push(l1);
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
        }
    }

    /// Consume the set into its canonical stored entries, in relation-id order.
    pub fn into_entries(self) -> Vec<([L1; N1], Vec<L2>, D)> {
        let Self {
            participants_1,
            f2_offsets,
            participants_2,
            data,
            ..
        } = self;
        let lengths = f2_offsets
            .windows(2)
            .map(|range| (range[1] - range[0]) as usize);
        let mut participants_2 = participants_2.into_iter();
        participants_1
            .into_iter()
            .zip(lengths)
            .zip(data)
            .map(|((participants_1, len), data)| {
                (
                    participants_1,
                    participants_2.by_ref().take(len).collect(),
                    data,
                )
            })
            .collect()
    }

    pub fn count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    /// Every relation as `(id, first-factor participants, second-factor participants, payload)` in
    /// relation-id order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (RelationId, &[L1; N1], &[L2], &D)> {
        let participants_1 = &self.participants_1;
        let participants_2 = &self.participants_2;
        let f2_offsets = &self.f2_offsets;
        self.data.iter().enumerate().map(move |(index, data)| {
            (
                RelationId(index as u32),
                &participants_1[index],
                &participants_2[f2_offsets[index] as usize..f2_offsets[index + 1] as usize],
                data,
            )
        })
    }

    /// Every relation as `(id, first-factor participants, second-factor participants, payload)` in
    /// relation-id order, the payload mutable.
    ///
    /// Participants stay immutable: changing them would invalidate the incidence index, which
    /// [`permute_1_with`](Self::permute_1_with) and [`permute_2_with`](Self::permute_2_with) are
    /// the operations allowed to leave intact.
    pub fn iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = (RelationId, &[L1; N1], &[L2], &mut D)> {
        let participants_1 = &self.participants_1;
        let participants_2 = &self.participants_2;
        let f2_offsets = &self.f2_offsets;
        self.data.iter_mut().enumerate().map(move |(index, data)| {
            (
                RelationId(index as u32),
                &participants_1[index],
                &participants_2[f2_offsets[index] as usize..f2_offsets[index + 1] as usize],
                data,
            )
        })
    }

    pub fn participants_1(&self, id: RelationId) -> &[L1; N1] {
        &self.participants_1[id.index()]
    }

    /// Reorder relation `id`'s first-factor participants so that `new[i] = old[order[i]]`,
    /// leaving the second factor and the payload untouched.
    ///
    /// The multiset is unchanged, so incidence answers identically and is not rebuilt. Panics
    /// unless `order` is a permutation of `0..arity`.
    pub fn permute_1_with(&mut self, id: RelationId, order: &[ParticipantPosition]) {
        permute_participants(self.participants_1[id.index()].as_mut_slice(), order);
    }

    pub fn participants_2(&self, id: RelationId) -> &[L2] {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        &self.participants_2[start..end]
    }

    /// Reorder relation `id`'s second-factor participants so that `new[i] = old[order[i]]`,
    /// leaving the first factor and the payload untouched.
    ///
    /// The multiset is unchanged, so incidence answers identically and is not rebuilt. Panics
    /// unless `order` is a permutation of `0..arity`.
    pub fn permute_2_with(&mut self, id: RelationId, order: &[ParticipantPosition]) {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        permute_participants(&mut self.participants_2[start..end], order);
    }

    /// Id of the relation coinciding with `query_1` / `query_2` — the one whose factors equal them
    /// as multisets, in any order.
    ///
    /// This is the identity question, not a lookup: the participant multiset is the relation's
    /// identity and the stored sequence is only the frame its payload is expressed in, so two
    /// entries presenting the same participants differently coincide. It is what
    /// [`pushout`](Self::pushout) and [`pullback`](Self::pullback) join on. Naming an entity by a
    /// subset of its constituents is a different question with a different key, and belongs to the
    /// caller that knows the key. §4.1 uniqueness ⇒ at most one hit.
    pub fn coincident(&self, node: NodeId, query_1: &[L1], query_2: &[L2]) -> Option<RelationId> {
        self.coincident_in(self.incident(node), query_1, query_2)
    }

    /// Edge-indexed peer of [`coincident`](Self::coincident).
    pub fn coincident_edge(
        &self,
        edge: EdgeId,
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<RelationId> {
        self.coincident_in(self.incident_edge(edge), query_1, query_2)
    }

    /// Whether relation `id` coincides with these participants — the known-id sibling of
    /// [`coincident`](Self::coincident), which searches for it instead.
    ///
    /// `pushout` and `pullback` apply this to a supplied pairing before gluing on it. A caller that
    /// already holds the id and needs identity established — because a frame-invariant payload
    /// carries without reading either frame — asks here rather than deriving the comparison again.
    pub fn is_coincident(&self, id: RelationId, query_1: &[L1], query_2: &[L2]) -> bool {
        self.coincident_in(&[id], query_1, query_2).is_some()
    }

    fn coincident_in(
        &self,
        candidates: &[RelationId],
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<RelationId> {
        let mut sorted_1: Vec<L1> = query_1.to_vec();
        sorted_1.sort_unstable();
        let mut sorted_2: Vec<L2> = query_2.to_vec();
        sorted_2.sort_unstable();
        candidates.iter().copied().find(|&id| {
            participants_match(self.participants_1(id), &sorted_1)
                && participants_match(self.participants_2(id), &sorted_2)
        })
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

    pub fn ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    /// Compact participant ids and drop relations containing a removed participant.
    ///
    /// Produces the same set as [`Self::tracked_compact`], without returning
    /// the relation-id compaction.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: Clone,
    {
        self.tracked_compact(compaction).0
    }

    /// Compact participant ids, dropping every relation that contains a removed participant, and
    /// report which relation ids the drop consumed.
    ///
    /// The returned compaction moves this set's own ids, so a caller holding relation ids can
    /// carry them across the removal without a second traversal.
    pub fn tracked_compact(&self, compaction: &GraphCompaction) -> (Self, Compaction<RelationId>)
    where
        D: Clone,
    {
        let mut removed = Vec::new();
        let mut entries: Vec<([L1; N1], Vec<L2>, D)> = Vec::with_capacity(self.count());
        for i in 0..self.count() {
            let rid = RelationId(i as u32);
            let f1: Option<[L1; N1]> = self
                .participants_1(rid)
                .iter()
                .map(|&p| p.compact(compaction))
                .collect::<Option<Vec<L1>>>()
                .and_then(|parts| parts.try_into().ok());
            let f2: Option<Vec<L2>> = self
                .participants_2(rid)
                .iter()
                .map(|&p| p.compact(compaction))
                .collect();
            match (f1, f2) {
                (Some(f1), Some(f2)) => entries.push((f1, f2, self.data(rid).clone())),
                _ => removed.push(rid),
            }
        }
        (
            Self::new(entries),
            Compaction::new(self.count(), removed)
                .expect("removed relations belong to the source set"),
        )
    }

    /// Relabel every participant, preserving rows, participant order, and payloads.
    ///
    /// # Semantic properties
    ///
    /// In both factors, each positional payload item remains attached to the participant whose id
    /// is relabeled.
    ///
    /// # Panics
    ///
    /// Panics when a participant lies outside the remapping's corresponding source range.
    pub fn remap(&self, remapping: &GraphRemapping) -> Self
    where
        D: Clone,
    {
        self.map_participants(
            |participant| Some(participant.remap(remapping)),
            |participant| Some(participant.remap(remapping)),
        )
        .expect("remapping transport supplies every participant")
    }

    /// Relabel participant ids through a correspondence without changing rows, frames, or payloads.
    ///
    /// # Panics
    /// Panics when any referenced node or edge has no image.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self
    where
        D: Clone,
    {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Relabel participant ids, returning `None` when any reference has no image.
    ///
    /// Only referenced ids require images; unrelated source entries may be unmatched.
    ///
    /// # Semantic properties
    /// Row ids, participant positions, and payloads are preserved. Identity mapping is exact;
    /// sequential covered mappings agree with their correspondence composition.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self>
    where
        D: Clone,
    {
        self.map_participants(
            |participant| participant.try_map(correspondence),
            |participant| participant.try_map(correspondence),
        )
    }

    fn map_participants(
        &self,
        mut map_1: impl FnMut(L1) -> Option<L1>,
        mut map_2: impl FnMut(L2) -> Option<L2>,
    ) -> Option<Self>
    where
        D: Clone,
    {
        let entries = self
            .ids()
            .map(|id| {
                let parts_1: [L1; N1] = self
                    .participants_1(id)
                    .iter()
                    .copied()
                    .map(&mut map_1)
                    .collect::<Option<Vec<_>>>()?
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                let parts_2: Vec<L2> = self
                    .participants_2(id)
                    .iter()
                    .copied()
                    .map(&mut map_2)
                    .collect::<Option<Vec<_>>>()?;
                Some((parts_1, parts_2, self.data(id).clone()))
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Self::new(entries))
    }

    /// Relabel every participant, returning `None` when the remapping does not cover either factor.
    pub fn try_remap(&self, remapping: &GraphRemapping) -> Option<Self>
    where
        D: Clone,
    {
        self.ids()
            .all(|id| {
                remappable_under(self.participants_1(id), remapping)
                    && remappable_under(self.participants_2(id), remapping)
            })
            .then(|| self.remap(remapping))
    }

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`]. Coincidence is equality of
    /// **both** factors' participants.
    pub fn pushout(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        mut combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<([L1; N1], Vec<L2>, D)> = self
            .ids()
            .map(|id| {
                (
                    *self.participants_1(id),
                    self.participants_2(id).to_vec(),
                    self.data(id).clone(),
                )
            })
            .collect();
        let self_count = entries.len();
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.count());
        for id in right.ids() {
            match coincident(self, right.participants_1(id), right.participants_2(id)) {
                Some(hit) => {
                    let merged = combine(
                        (
                            self.participants_1(hit),
                            self.participants_2(hit),
                            self.data(hit),
                        ),
                        (
                            right.participants_1(id),
                            right.participants_2(id),
                            right.data(id),
                        ),
                    )?;
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

    /// Same-space relation pullback — see [`FixedRelationSet::pullback`].
    pub fn pullback(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        mut combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<RelationPullback<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<([L1; N1], Vec<L2>, D)> = Vec::new();
        let mut left_images: Vec<RelationId> = Vec::new();
        let mut right_images: Vec<RelationId> = Vec::new();
        for id in self.ids() {
            if let Some(hit) = coincident(right, self.participants_1(id), self.participants_2(id)) {
                let merged = combine(
                    (
                        self.participants_1(id),
                        self.participants_2(id),
                        self.data(id),
                    ),
                    (
                        right.participants_1(hit),
                        right.participants_2(hit),
                        right.data(hit),
                    ),
                )?;
                entries.push((
                    *self.participants_1(id),
                    self.participants_2(id).to_vec(),
                    merged,
                ));
                left_images.push(id);
                right_images.push(hit);
            }
        }
        Some(relation_pullback(
            Self::new(entries),
            left_images,
            right_images,
            self.count(),
            right.count(),
        ))
    }
}

impl<L1, const N1: usize, L2, D> Default for FixedVarBirelationSet<L1, N1, L2, D> {
    fn default() -> Self {
        Self {
            participants_1: Vec::new(),
            f2_offsets: vec![0],
            participants_2: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
        }
    }
}

/// Birelation with two variable-arity factors.
#[derive(Clone, Debug)]
pub struct VarVarBirelationSet<L1, L2, D> {
    f1_offsets: Vec<u32>,
    participants_1: Vec<L1>,
    f2_offsets: Vec<u32>,
    participants_2: Vec<L2>,
    data: Vec<D>,
    incidence: Incidence,
}

impl<L1, L2, D> PartialEq for VarVarBirelationSet<L1, L2, D>
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

impl<L1, L2, D> Eq for VarVarBirelationSet<L1, L2, D>
where
    L1: Eq,
    L2: Eq,
    D: Eq,
{
}

impl<L1, L2, D> VarVarBirelationSet<L1, L2, D>
where
    L1: RelationParticipant,
    L2: RelationParticipant,
{
    pub fn new(entries: Vec<(Vec<L1>, Vec<L2>, D)>) -> Self
    where
        D: Clone,
    {
        let relation_count = entries.len();
        let mut f1_offsets = Vec::with_capacity(relation_count + 1);
        f1_offsets.push(0);
        let mut participants_1 = Vec::new();
        let mut f2_offsets = Vec::with_capacity(relation_count + 1);
        f2_offsets.push(0);
        let mut participants_2 = Vec::new();
        let mut data = Vec::with_capacity(relation_count);
        for (l1, l2, d) in entries {
            participants_1.extend_from_slice(&l1);
            f1_offsets.push(participants_1.len() as u32);
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
        }
    }

    /// Consume the set into its canonical stored entries, in relation-id order.
    pub fn into_entries(self) -> Vec<(Vec<L1>, Vec<L2>, D)> {
        let Self {
            f1_offsets,
            participants_1,
            f2_offsets,
            participants_2,
            data,
            ..
        } = self;
        let lengths_1 = f1_offsets
            .windows(2)
            .map(|range| (range[1] - range[0]) as usize);
        let lengths_2 = f2_offsets
            .windows(2)
            .map(|range| (range[1] - range[0]) as usize);
        let mut participants_1 = participants_1.into_iter();
        let mut participants_2 = participants_2.into_iter();
        lengths_1
            .zip(lengths_2)
            .zip(data)
            .map(|((len_1, len_2), data)| {
                (
                    participants_1.by_ref().take(len_1).collect(),
                    participants_2.by_ref().take(len_2).collect(),
                    data,
                )
            })
            .collect()
    }

    pub fn count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    /// Every relation as `(id, first-factor participants, second-factor participants, payload)` in
    /// relation-id order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (RelationId, &[L1], &[L2], &D)> {
        let participants_1 = &self.participants_1;
        let participants_2 = &self.participants_2;
        let f1_offsets = &self.f1_offsets;
        let f2_offsets = &self.f2_offsets;
        self.data.iter().enumerate().map(move |(index, data)| {
            (
                RelationId(index as u32),
                &participants_1[f1_offsets[index] as usize..f1_offsets[index + 1] as usize],
                &participants_2[f2_offsets[index] as usize..f2_offsets[index + 1] as usize],
                data,
            )
        })
    }

    /// Every relation as `(id, first-factor participants, second-factor participants, payload)` in
    /// relation-id order, the payload mutable.
    ///
    /// Participants stay immutable: changing them would invalidate the incidence index, which
    /// [`permute_1_with`](Self::permute_1_with) and [`permute_2_with`](Self::permute_2_with) are
    /// the operations allowed to leave intact.
    pub fn iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = (RelationId, &[L1], &[L2], &mut D)> {
        let participants_1 = &self.participants_1;
        let participants_2 = &self.participants_2;
        let f1_offsets = &self.f1_offsets;
        let f2_offsets = &self.f2_offsets;
        self.data.iter_mut().enumerate().map(move |(index, data)| {
            (
                RelationId(index as u32),
                &participants_1[f1_offsets[index] as usize..f1_offsets[index + 1] as usize],
                &participants_2[f2_offsets[index] as usize..f2_offsets[index + 1] as usize],
                data,
            )
        })
    }

    pub fn participants_1(&self, id: RelationId) -> &[L1] {
        let start = self.f1_offsets[id.index()] as usize;
        let end = self.f1_offsets[id.index() + 1] as usize;
        &self.participants_1[start..end]
    }

    /// Reorder relation `id`'s first-factor participants so that `new[i] = old[order[i]]`,
    /// leaving the second factor and the payload untouched.
    ///
    /// The multiset is unchanged, so incidence answers identically and is not rebuilt. Panics
    /// unless `order` is a permutation of `0..arity`.
    pub fn permute_1_with(&mut self, id: RelationId, order: &[ParticipantPosition]) {
        let start = self.f1_offsets[id.index()] as usize;
        let end = self.f1_offsets[id.index() + 1] as usize;
        permute_participants(&mut self.participants_1[start..end], order);
    }

    pub fn participants_2(&self, id: RelationId) -> &[L2] {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        &self.participants_2[start..end]
    }

    /// Reorder relation `id`'s second-factor participants so that `new[i] = old[order[i]]`,
    /// leaving the first factor and the payload untouched.
    ///
    /// The multiset is unchanged, so incidence answers identically and is not rebuilt. Panics
    /// unless `order` is a permutation of `0..arity`.
    pub fn permute_2_with(&mut self, id: RelationId, order: &[ParticipantPosition]) {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        permute_participants(&mut self.participants_2[start..end], order);
    }

    /// Id of the relation coinciding with `query_1` / `query_2` — the one whose factors equal them
    /// as multisets, in any order.
    ///
    /// This is the identity question, not a lookup: the participant multiset is the relation's
    /// identity and the stored sequence is only the frame its payload is expressed in, so two
    /// entries presenting the same participants differently coincide. It is what
    /// [`pushout`](Self::pushout) and [`pullback`](Self::pullback) join on. Naming an entity by a
    /// subset of its constituents is a different question with a different key, and belongs to the
    /// caller that knows the key. §4.1 uniqueness ⇒ at most one hit.
    pub fn coincident(&self, node: NodeId, query_1: &[L1], query_2: &[L2]) -> Option<RelationId> {
        self.coincident_in(self.incident(node), query_1, query_2)
    }

    /// Edge-indexed peer of [`coincident`](Self::coincident).
    pub fn coincident_edge(
        &self,
        edge: EdgeId,
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<RelationId> {
        self.coincident_in(self.incident_edge(edge), query_1, query_2)
    }

    /// Whether relation `id` coincides with these participants — the known-id sibling of
    /// [`coincident`](Self::coincident), which searches for it instead.
    ///
    /// `pushout` and `pullback` apply this to a supplied pairing before gluing on it. A caller that
    /// already holds the id and needs identity established — because a frame-invariant payload
    /// carries without reading either frame — asks here rather than deriving the comparison again.
    pub fn is_coincident(&self, id: RelationId, query_1: &[L1], query_2: &[L2]) -> bool {
        self.coincident_in(&[id], query_1, query_2).is_some()
    }

    fn coincident_in(
        &self,
        candidates: &[RelationId],
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<RelationId> {
        let mut sorted_1: Vec<L1> = query_1.to_vec();
        sorted_1.sort_unstable();
        let mut sorted_2: Vec<L2> = query_2.to_vec();
        sorted_2.sort_unstable();
        candidates.iter().copied().find(|&id| {
            participants_match(self.participants_1(id), &sorted_1)
                && participants_match(self.participants_2(id), &sorted_2)
        })
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

    pub fn ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    /// Compact participant ids and drop relations containing a removed participant.
    ///
    /// Produces the same set as [`Self::tracked_compact`], without returning
    /// the relation-id compaction.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: Clone,
    {
        self.tracked_compact(compaction).0
    }

    /// Compact participant ids, dropping every relation that contains a removed participant, and
    /// report which relation ids the drop consumed.
    ///
    /// The returned compaction moves this set's own ids, so a caller holding relation ids can
    /// carry them across the removal without a second traversal.
    pub fn tracked_compact(&self, compaction: &GraphCompaction) -> (Self, Compaction<RelationId>)
    where
        D: Clone,
    {
        let mut removed = Vec::new();
        let mut entries: Vec<(Vec<L1>, Vec<L2>, D)> = Vec::with_capacity(self.count());
        for i in 0..self.count() {
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
            match (f1, f2) {
                (Some(f1), Some(f2)) => entries.push((f1, f2, self.data(rid).clone())),
                _ => removed.push(rid),
            }
        }
        (
            Self::new(entries),
            Compaction::new(self.count(), removed)
                .expect("removed relations belong to the source set"),
        )
    }

    /// Relabel every participant, preserving rows, participant order, and payloads.
    ///
    /// # Semantic properties
    ///
    /// In both factors, each positional payload item remains attached to the participant whose id
    /// is relabeled.
    ///
    /// # Panics
    ///
    /// Panics when a participant lies outside the remapping's corresponding source range.
    pub fn remap(&self, remapping: &GraphRemapping) -> Self
    where
        D: Clone,
    {
        self.map_participants(
            |participant| Some(participant.remap(remapping)),
            |participant| Some(participant.remap(remapping)),
        )
        .expect("remapping transport supplies every participant")
    }

    /// Relabel participant ids through a correspondence without changing rows, frames, or payloads.
    ///
    /// # Panics
    /// Panics when any referenced node or edge has no image.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self
    where
        D: Clone,
    {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Relabel participant ids, returning `None` when any reference has no image.
    ///
    /// Only referenced ids require images; unrelated source entries may be unmatched.
    ///
    /// # Semantic properties
    /// Row ids, participant positions, and payloads are preserved. Identity mapping is exact;
    /// sequential covered mappings agree with their correspondence composition.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self>
    where
        D: Clone,
    {
        self.map_participants(
            |participant| participant.try_map(correspondence),
            |participant| participant.try_map(correspondence),
        )
    }

    fn map_participants(
        &self,
        mut map_1: impl FnMut(L1) -> Option<L1>,
        mut map_2: impl FnMut(L2) -> Option<L2>,
    ) -> Option<Self>
    where
        D: Clone,
    {
        let entries = self
            .ids()
            .map(|id| {
                let parts_1: Vec<L1> = self
                    .participants_1(id)
                    .iter()
                    .copied()
                    .map(&mut map_1)
                    .collect::<Option<Vec<_>>>()?;
                let parts_2: Vec<L2> = self
                    .participants_2(id)
                    .iter()
                    .copied()
                    .map(&mut map_2)
                    .collect::<Option<Vec<_>>>()?;
                Some((parts_1, parts_2, self.data(id).clone()))
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Self::new(entries))
    }

    /// Relabel every participant, returning `None` when the remapping does not cover either factor.
    pub fn try_remap(&self, remapping: &GraphRemapping) -> Option<Self>
    where
        D: Clone,
    {
        self.ids()
            .all(|id| {
                remappable_under(self.participants_1(id), remapping)
                    && remappable_under(self.participants_2(id), remapping)
            })
            .then(|| self.remap(remapping))
    }

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`]. Coincidence is equality of
    /// both factors' participants.
    pub fn pushout(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        mut combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<(Vec<L1>, Vec<L2>, D)> = self
            .ids()
            .map(|id| {
                (
                    self.participants_1(id).to_vec(),
                    self.participants_2(id).to_vec(),
                    self.data(id).clone(),
                )
            })
            .collect();
        let self_count = entries.len();
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.count());
        for id in right.ids() {
            match coincident(self, right.participants_1(id), right.participants_2(id)) {
                Some(hit) => {
                    let merged = combine(
                        (
                            self.participants_1(hit),
                            self.participants_2(hit),
                            self.data(hit),
                        ),
                        (
                            right.participants_1(id),
                            right.participants_2(id),
                            right.data(id),
                        ),
                    )?;
                    entries[hit.index()].2 = merged;
                    right_map.push(hit);
                }
                None => {
                    right_map.push(RelationId(entries.len() as u32));
                    entries.push((
                        right.participants_1(id).to_vec(),
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

    /// Same-space relation pullback — see [`FixedRelationSet::pullback`].
    pub fn pullback(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        mut combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<RelationPullback<Self>>
    where
        D: Clone,
    {
        let mut entries: Vec<(Vec<L1>, Vec<L2>, D)> = Vec::new();
        let mut left_images: Vec<RelationId> = Vec::new();
        let mut right_images: Vec<RelationId> = Vec::new();
        for id in self.ids() {
            if let Some(hit) = coincident(right, self.participants_1(id), self.participants_2(id)) {
                let merged = combine(
                    (
                        self.participants_1(id),
                        self.participants_2(id),
                        self.data(id),
                    ),
                    (
                        right.participants_1(hit),
                        right.participants_2(hit),
                        right.data(hit),
                    ),
                )?;
                entries.push((
                    self.participants_1(id).to_vec(),
                    self.participants_2(id).to_vec(),
                    merged,
                ));
                left_images.push(id);
                right_images.push(hit);
            }
        }
        Some(relation_pullback(
            Self::new(entries),
            left_images,
            right_images,
            self.count(),
            right.count(),
        ))
    }
}

impl<L1, L2, D> Default for VarVarBirelationSet<L1, L2, D> {
    fn default() -> Self {
        Self {
            f1_offsets: vec![0],
            participants_1: Vec::new(),
            f2_offsets: vec![0],
            participants_2: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;
    use std::hash::{DefaultHasher, Hash, Hasher};

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::correspondence::GraphCorrespondence;
    use crate::graph::Graph;
    use crate::remap::Remapping;

    #[fixture]
    fn participant_correspondence() -> GraphCorrespondence {
        GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(5)), (NodeId(2), NodeId(1))], 4, 6)
                .unwrap(),
            Correspondence::new(vec![(EdgeId(0), EdgeId(6)), (EdgeId(2), EdgeId(3))], 4, 7)
                .unwrap(),
        )
    }

    #[rstest]
    #[case::mapped(NodeId(2), Some(NodeId(1)))]
    #[case::unmatched(NodeId(1), None)]
    #[case::outside(NodeId(4), None)]
    fn test_node_id_try_map(
        participant_correspondence: GraphCorrespondence,
        #[case] id: NodeId,
        #[case] expected: Option<NodeId>,
    ) {
        assert_eq!(id.try_map(&participant_correspondence), expected);
    }

    #[rstest]
    #[case::mapped(EdgeId(2), Some(EdgeId(3)))]
    #[case::unmatched(EdgeId(1), None)]
    #[case::outside(EdgeId(4), None)]
    fn test_edge_id_try_map(
        participant_correspondence: GraphCorrespondence,
        #[case] id: EdgeId,
        #[case] expected: Option<EdgeId>,
    ) {
        assert_eq!(id.try_map(&participant_correspondence), expected);
    }

    fn hash<T: Hash>(value: &T) -> u64 {
        let mut hasher = DefaultHasher::new();
        value.hash(&mut hasher);
        hasher.finish()
    }

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct PositionLabels(Vec<u32>);

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct BiPositionLabels {
        factor_1: Vec<u32>,
        factor_2: Vec<u32>,
    }

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    fn assert_exact_size<T>(mut iterator: impl ExactSizeIterator<Item = T>, expected: Vec<T>)
    where
        T: Debug + PartialEq,
    {
        assert_eq!(iterator.len(), expected.len());
        assert_eq!(iterator.size_hint(), (expected.len(), Some(expected.len())));
        while let Some(expected_item) = expected.get(expected.len() - iterator.len()) {
            let previous = iterator.len();
            assert_eq!(iterator.next().as_ref(), Some(expected_item));
            let remaining = iterator.len();
            assert_eq!(remaining, previous - 1);
            assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
        }
        assert_eq!(iterator.next(), None);
        assert_eq!(iterator.len(), 0);
        assert_eq!(iterator.size_hint(), (0, Some(0)));
    }

    #[rstest]
    #[case::before_removed(NodeId(0), Some(NodeId(0)))]
    #[case::removed(NodeId(1), None)]
    #[case::after_removed(NodeId(2), Some(NodeId(1)))]
    fn test_node_id_compact(#[case] id: NodeId, #[case] expected: Option<NodeId>) {
        let compaction = GraphCompaction::new(
            Compaction::new(3, vec![NodeId(1)]).unwrap(),
            Compaction::identity(0),
        );
        assert_eq!(id.compact(&compaction), expected);
    }

    #[rstest]
    #[case::before_gap(NodeId(0), NodeId(0))]
    #[case::after_gap(NodeId(1), NodeId(2))]
    fn test_node_id_uncompact(#[case] id: NodeId, #[case] expected: NodeId) {
        let compaction = GraphCompaction::new(
            Compaction::new(3, vec![NodeId(1)]).unwrap(),
            Compaction::identity(0),
        );
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
        let compaction = GraphCompaction::new(
            Compaction::identity(0),
            Compaction::new(3, vec![EdgeId(0)]).unwrap(),
        );
        assert_eq!(id.compact(&compaction), expected);
    }

    #[rstest]
    #[case::before_gap(EdgeId(0), EdgeId(0))]
    #[case::after_gap(EdgeId(1), EdgeId(2))]
    fn test_edge_id_uncompact(#[case] id: EdgeId, #[case] expected: EdgeId) {
        let compaction = GraphCompaction::new(
            Compaction::identity(0),
            Compaction::new(3, vec![EdgeId(1)]).unwrap(),
        );
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
        let rs: FixedRelationSet<NodeId, &str, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], "dative"), ([n(1), n(2)], "noncov")]);
        assert_eq!(rs.count(), 2);
        assert_eq!(rs.data(RelationId(0)), &"dative");
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(2)]);
    }

    #[rstest]
    fn test_fixed_relation_set_hash() {
        let entries = vec![([n(2), n(0)], "first"), ([n(3), n(1)], "second")];
        let left: FixedRelationSet<NodeId, &str, 2> = FixedRelationSet::new(entries.clone());
        let right: FixedRelationSet<NodeId, &str, 2> = FixedRelationSet::new(entries);
        assert_eq!(left, right);
        assert_eq!(hash(&left), hash(&right));
    }

    #[rstest]
    #[case::roundtrip(
        vec![([n(2), n(0)], "first"), ([n(3), n(1)], "second")],
    )]
    fn test_fixed_relation_set_into_entries(#[case] entries: Vec<([NodeId; 2], &str)>) {
        let rs: FixedRelationSet<NodeId, &str, 2> = FixedRelationSet::new(entries.clone());
        assert_eq!(rs.into_entries(), entries);
    }

    #[rstest]
    fn test_fixed_relation_set_data_mut() {
        let mut rs: FixedRelationSet<NodeId, i32, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    /// What frame a pushout coincidence comes out in, and whether the two sides reaching `combine`
    /// are in the same frame. Recorded so the answer is read off a run rather than argued.
    ///
    /// Both sides carry the same two participants in opposite order with a position-labelled
    /// payload. `combine` records what it was handed and returns the left payload untouched, so the
    /// test observes the inputs rather than any realignment.
    #[rstest]
    #[case::both_frames_kept(
        VarRelationSet::<NodeId, PositionLabels>::new(
            vec![(vec![n(7), n(3)], PositionLabels(vec![70, 30]))]),
        VarRelationSet::<NodeId, PositionLabels>::new(
            vec![(vec![n(3), n(7)], PositionLabels(vec![30, 70]))]),
        [n(7), n(3)],
        PositionLabels(vec![70, 30]),
        PositionLabels(vec![30, 70]),
    )]
    fn test_var_relation_set_pushout_coincidence_frame(
        #[case] left: VarRelationSet<NodeId, PositionLabels>,
        #[case] right: VarRelationSet<NodeId, PositionLabels>,
        #[case] expected_frame: [NodeId; 2],
        #[case] expected_left_seen: PositionLabels,
        #[case] expected_right_seen: PositionLabels,
    ) {
        let mut seen = None;
        let merged = left
            .pushout(
                &right,
                |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
                |(_, a), (_, b)| {
                    seen = Some((a.clone(), b.clone()));
                    Some(a.clone())
                },
            )
            .expect("combine never rejects here");
        let object = merged.object;

        assert_eq!(object.count(), 1, "the two entries coincide");
        assert_eq!(
            object.participants(RelationId(0)),
            expected_frame.as_slice()
        );
        assert_eq!(
            seen,
            Some((expected_left_seen, expected_right_seen)),
            "the two payloads as `combine` received them"
        );
    }

    /// Two right entries coinciding with one left entry is rejected, not merged: the right
    /// coprojection must be injective, and `Correspondence` asserts that.
    ///
    /// This is why `pushout` may read the left payload out of its own output buffer — the entry can
    /// never be merged into twice. `pullback` reads the same payload from the source instead, and
    /// the two agree because the case that would separate them cannot be constructed.
    #[rstest]
    #[should_panic(expected = "correspondence images must be unique")]
    fn test_var_relation_set_pushout_repeated_coincidence() {
        let left: VarRelationSet<NodeId, PositionLabels> =
            VarRelationSet::new(vec![(vec![n(0), n(1)], PositionLabels(vec![1, 1]))]);
        let right: VarRelationSet<NodeId, PositionLabels> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], PositionLabels(vec![2, 2])),
            (vec![n(0), n(1)], PositionLabels(vec![4, 4])),
        ]);

        left.pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
            |(_, a), (_, b)| {
                Some(PositionLabels(
                    a.0.iter().zip(&b.0).map(|(x, y)| x + y).collect(),
                ))
            },
        );
    }

    #[rstest]
    fn test_fixed_relation_set_iter() {
        let empty = FixedRelationSet::<NodeId, i32, 2>::default();
        assert_eq!(empty.iter().collect::<Vec<_>>(), vec![]);

        let rs: FixedRelationSet<NodeId, i32, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], 1), ([n(1), n(2)], 2)]);
        assert_eq!(rs.iter().len(), 2);
        assert_eq!(
            rs.iter().collect::<Vec<_>>(),
            vec![
                (RelationId(0), &[n(0), n(1)], &1),
                (RelationId(1), &[n(1), n(2)], &2),
            ],
        );
    }

    #[rstest]
    fn test_fixed_relation_set_iter_mut() {
        let mut empty = FixedRelationSet::<NodeId, i32, 2>::default();
        assert_eq!(empty.iter_mut().len(), 0);

        let mut rs: FixedRelationSet<NodeId, i32, 2> = FixedRelationSet::new(vec![
            ([n(0), n(1)], 1),
            ([n(1), n(2)], 2),
            ([n(2), n(3)], 3),
        ]);
        assert_eq!(rs.iter_mut().len(), 3);
        for (id, participants, data) in rs.iter_mut() {
            assert_eq!(participants[0], n(id.index() as u32));
            *data *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
        assert_eq!(rs.data(RelationId(2)), &30);
    }

    #[rstest]
    fn test_fixed_relation_set_participants_ordered() {
        let rs: FixedRelationSet<NodeId, &str, 2> =
            FixedRelationSet::new(vec![([n(2), n(0)], "a"), ([n(3), n(1)], "b")]);
        assert_eq!(rs.participants(RelationId(0)), &[n(2), n(0)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(3), n(1)]);
    }

    #[rstest]
    #[case::rotation(vec![ParticipantPosition(2), ParticipantPosition(0), ParticipantPosition(1)], [n(2), n(0), n(1)])]
    #[case::transposition(vec![ParticipantPosition(1), ParticipantPosition(0), ParticipantPosition(2)], [n(1), n(0), n(2)])]
    #[case::reversal(vec![ParticipantPosition(2), ParticipantPosition(1), ParticipantPosition(0)], [n(2), n(1), n(0)])]
    fn test_fixed_relation_set_permute_with(
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: [NodeId; 3],
    ) {
        let mut rs: FixedRelationSet<NodeId, &str, 3> =
            FixedRelationSet::new(vec![([n(0), n(1), n(2)], "a"), ([n(3), n(4), n(5)], "b")]);
        let incidence_before: Vec<Vec<RelationId>> =
            (0..6).map(|i| rs.incident(n(i)).to_vec()).collect();

        rs.permute_with(RelationId(0), &order);

        assert_eq!(rs.participants(RelationId(0)), &expected);
        assert_eq!(rs.participants(RelationId(1)), &[n(3), n(4), n(5)]);
        assert_eq!(rs.data(RelationId(0)), &"a");
        assert_eq!(rs.data(RelationId(1)), &"b");
        let incidence_after: Vec<Vec<RelationId>> =
            (0..6).map(|i| rs.incident(n(i)).to_vec()).collect();
        assert_eq!(incidence_after, incidence_before);
    }

    #[rstest]
    #[case::ordered_factor(FixedRelationSet::<NodeId, &str, 3>::new(vec![
        ([n(2), n(0), n(1)], "a"),
    ]))]
    #[case::unordered_factor(FixedRelationSet::<NodeId, &str, 3>::new(vec![
        ([n(2), n(0), n(1)], "a"),
    ]))]
    fn test_fixed_relation_set_permute_with_identity(
        #[case] input: FixedRelationSet<NodeId, &'static str, 3>,
    ) {
        let mut permuted = input.clone();
        permuted.permute_with(
            RelationId(0),
            &[
                ParticipantPosition(0),
                ParticipantPosition(1),
                ParticipantPosition(2),
            ],
        );
        assert_eq!(permuted, input);
    }

    #[rstest]
    #[case::order_too_short(vec![ParticipantPosition(0), ParticipantPosition(1)])]
    #[case::order_too_long(vec![ParticipantPosition(0), ParticipantPosition(1), ParticipantPosition(2), ParticipantPosition(0)])]
    #[case::position_out_of_range(vec![ParticipantPosition(0), ParticipantPosition(1), ParticipantPosition(3)])]
    #[case::position_repeated(vec![ParticipantPosition(0), ParticipantPosition(0), ParticipantPosition(1)])]
    #[should_panic(expected = "permute")]
    fn test_fixed_relation_set_permute_with_error(#[case] order: Vec<ParticipantPosition>) {
        let mut rs: FixedRelationSet<NodeId, &str, 3> =
            FixedRelationSet::new(vec![([n(0), n(1), n(2)], "a")]);
        rs.permute_with(RelationId(0), &order);
    }

    #[rstest]
    fn test_fixed_relation_set_incidence() {
        let rs: FixedRelationSet<NodeId, (), 2> = FixedRelationSet::new(vec![
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
        let rs: FixedRelationSet<NodeId, (), 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], ()), ([n(1), n(2)], ())]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_fixed_relation_set_relation_ids() {
        assert_exact_size(FixedRelationSet::<NodeId, (), 2>::default().ids(), vec![]);
        let rs: FixedRelationSet<NodeId, (), 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], ()), ([n(1), n(2)], ())]);
        assert_exact_size(rs.ids(), vec![RelationId(0), RelationId(1)]);
    }

    #[fixture]
    fn fixed_relation_set_compaction_input() -> FixedRelationSet<NodeId, &'static str, 2> {
        FixedRelationSet::new(vec![
            ([NodeId(0), NodeId(2)], "keep"),
            ([NodeId(1), NodeId(3)], "drop"),
        ])
    }

    #[rstest]
    #[case::partial(
        vec![NodeId(1)],
        FixedRelationSet::new(vec![([NodeId(0), NodeId(1)], "keep")]),
        vec![RelationId(1)],
    )]
    #[case::all(
        vec![NodeId(0), NodeId(1)],
        FixedRelationSet::default(),
        vec![RelationId(0), RelationId(1)],
    )]
    fn test_fixed_relation_set_tracked_compact(
        fixed_relation_set_compaction_input: FixedRelationSet<NodeId, &'static str, 2>,
        #[case] removed_nodes: Vec<NodeId>,
        #[case] expected: FixedRelationSet<NodeId, &'static str, 2>,
        #[case] removed_relations: Vec<RelationId>,
    ) {
        let input = fixed_relation_set_compaction_input;
        let compaction = GraphCompaction::new(
            Compaction::new(4, removed_nodes).unwrap(),
            Compaction::identity(0),
        );
        let (output, witness) = input.tracked_compact(&compaction);
        assert_eq!(input.compact(&compaction), expected);
        assert_eq!(output, expected);
        assert_eq!(
            witness,
            Compaction::new(2, removed_relations.clone()).unwrap()
        );
        let survivors = (0..2)
            .map(RelationId)
            .filter(|id| !removed_relations.contains(id))
            .collect::<Vec<_>>();
        for (idx, &old) in survivors.iter().enumerate() {
            assert_eq!(witness.compact(old), Some(RelationId::from(idx)));
        }
    }

    #[rstest]
    #[case::empty(FixedRelationSet::default())]
    #[case::rows(
        FixedRelationSet::new(vec![([NodeId(0), NodeId(2)], "keep"), ([NodeId(1), NodeId(3)], "drop")]),
    )]
    fn test_fixed_relation_set_compact_identity(
        #[case] input: FixedRelationSet<NodeId, &'static str, 2>,
    ) {
        let compaction = GraphCompaction::new(Compaction::identity(4), Compaction::identity(0));
        assert_eq!(input.compact(&compaction), input);
        assert_eq!(
            input.tracked_compact(&compaction),
            (input.clone(), Compaction::identity(input.count())),
        );
    }

    #[rstest]
    #[case::rows(FixedRelationSet::new(vec![([NodeId(2), NodeId(0)], vec![7, 11]), ([NodeId(2), NodeId(0)], vec![13, 17])]),
        FixedRelationSet::new(vec![([NodeId(1), NodeId(5)], vec![7, 11]), ([NodeId(1), NodeId(5)], vec![13, 17])]))]
    fn test_fixed_relation_set_map(
        participant_correspondence: GraphCorrespondence,
        #[case] input: FixedRelationSet<NodeId, Vec<u32>, 2>,
        #[case] expected: FixedRelationSet<NodeId, Vec<u32>, 2>,
    ) {
        assert_eq!(input.map(&participant_correspondence), expected);
        assert_eq!(
            input.try_map(&participant_correspondence),
            Some(expected.clone())
        );
        let reverse = GraphCorrespondence::new(
            participant_correspondence.nodes().reverse(),
            participant_correspondence.edges().reverse(),
        );
        assert_eq!(expected.map(&reverse), input);
        let composed = participant_correspondence.compose(&reverse).unwrap();
        assert_eq!(input.map(&composed), input);
        assert_eq!(
            expected.incident(NodeId(1)),
            &[RelationId(0), RelationId(1)]
        );
    }

    #[rstest]
    #[case::empty(FixedRelationSet::new(vec![]))]
    #[case::rows(FixedRelationSet::new(vec![([NodeId(2), NodeId(0)], vec![7, 11]), ([NodeId(2), NodeId(0)], vec![13, 17])]))]
    fn test_fixed_relation_set_map_identity(#[case] input: FixedRelationSet<NodeId, Vec<u32>, 2>) {
        let identity = GraphCorrespondence::new(
            Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2), NodeId(3)], 4),
            Correspondence::from_images(&[EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)], 4),
        );
        assert_eq!(input.try_map(&identity), Some(input.clone()));
        assert_eq!(input.map(&identity), input);
    }

    #[rstest]
    #[case::missing_node(1)]
    #[case::outside_node(4)]
    fn test_fixed_relation_set_try_map_error(
        participant_correspondence: GraphCorrespondence,
        #[case] node: u32,
    ) {
        let input: FixedRelationSet<NodeId, Vec<u32>, 2> = FixedRelationSet::new(vec![
            ([NodeId(2), NodeId(0)], vec![7, 11]),
            ([NodeId(node), NodeId(0)], vec![13, 17]),
        ]);
        assert_eq!(input.try_map(&participant_correspondence), None);
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_fixed_relation_set_map_error(participant_correspondence: GraphCorrespondence) {
        let node = 1;

        let input: FixedRelationSet<NodeId, Vec<u32>, 2> =
            FixedRelationSet::new(vec![([NodeId(node), NodeId(0)], vec![7, 11])]);
        input.map(&participant_correspondence);
    }

    #[rstest]
    fn test_fixed_relation_set_remap() {
        let rs: FixedRelationSet<NodeId, PositionLabels, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], PositionLabels(vec![10, 11]))]);
        let remapping = GraphRemapping::new(
            Remapping::new(vec![n(1), n(0)]).expect("permutation images"),
            Remapping::new(vec![]).expect("permutation images"),
        );
        let out = rs.remap(&remapping);
        assert_eq!(out.participants(RelationId(0)), &[n(1), n(0)]);
        assert_eq!(out.data(RelationId(0)), &PositionLabels(vec![10, 11]));
    }

    #[rstest]
    #[case::covered(vec![n(1), n(0)], true)]
    #[case::uncovered_node(vec![n(0)], false)]
    fn test_fixed_relation_set_try_remap(#[case] nodes: Vec<NodeId>, #[case] covered: bool) {
        let rs: FixedRelationSet<NodeId, PositionLabels, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], PositionLabels(vec![10, 11]))]);
        let remapping = GraphRemapping::new(
            Remapping::new(nodes).expect("permutation images"),
            Remapping::new(vec![]).expect("permutation images"),
        );
        let expected = covered.then(|| rs.remap(&remapping));
        assert_eq!(rs.try_remap(&remapping), expected);
    }

    #[rstest]
    fn test_fixed_relation_set_default() {
        let rs = FixedRelationSet::<NodeId, (), 2>::default();
        assert_eq!(rs.count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_var_relation_set_new() {
        let rs: VarRelationSet<NodeId, &str> =
            VarRelationSet::new(vec![(vec![n(0), n(1), n(2), n(3), n(4), n(5)], "benzene")]);
        assert_eq!(rs.count(), 1);
        assert_eq!(rs.data(RelationId(0)), &"benzene");
        assert_eq!(
            rs.participants(RelationId(0)),
            &[n(0), n(1), n(2), n(3), n(4), n(5)]
        );
    }

    #[rstest]
    fn test_var_relation_set_hash() {
        let entries = vec![
            (vec![n(2), n(0)], "first"),
            (vec![n(4), n(3), n(1)], "second"),
        ];
        let left: VarRelationSet<NodeId, &str> = VarRelationSet::new(entries.clone());
        let right: VarRelationSet<NodeId, &str> = VarRelationSet::new(entries);
        assert_eq!(left, right);
        assert_eq!(hash(&left), hash(&right));
    }

    #[rstest]
    #[case::roundtrip(
        vec![
            (vec![n(2), n(0)], "first"),
            (vec![n(4), n(3), n(1)], "second"),
        ],
    )]
    fn test_var_relation_set_into_entries(#[case] entries: Vec<(Vec<NodeId>, &str)>) {
        let rs: VarRelationSet<NodeId, &str> = VarRelationSet::new(entries.clone());
        assert_eq!(rs.into_entries(), entries);
    }

    #[rstest]
    fn test_var_relation_set_data_mut() {
        let mut rs: VarRelationSet<NodeId, i32> =
            VarRelationSet::new(vec![(vec![n(0), n(1), n(2)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_var_relation_set_iter() {
        let empty = VarRelationSet::<NodeId, i32>::default();
        assert_eq!(empty.iter().collect::<Vec<_>>(), vec![]);

        let rs: VarRelationSet<NodeId, i32> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], 1),
            (vec![n(2), n(3), n(4)], 2),
            (vec![n(5)], 3),
        ]);
        assert_eq!(rs.iter().len(), 3);
        assert_eq!(
            rs.iter().collect::<Vec<_>>(),
            vec![
                (RelationId(0), [n(0), n(1)].as_slice(), &1),
                (RelationId(1), [n(2), n(3), n(4)].as_slice(), &2),
                (RelationId(2), [n(5)].as_slice(), &3),
            ],
        );
    }

    #[rstest]
    fn test_var_relation_set_iter_mut() {
        let mut empty = VarRelationSet::<NodeId, i32>::default();
        assert_eq!(empty.iter_mut().len(), 0);

        let mut rs: VarRelationSet<NodeId, i32> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], 1),
            (vec![n(2), n(3), n(4)], 2),
            (vec![n(5)], 3),
        ]);
        let arities: Vec<usize> = rs
            .iter_mut()
            .map(|(_, participants, data)| {
                *data *= 10;
                participants.len()
            })
            .collect();
        assert_eq!(arities, vec![2, 3, 1]);
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
        assert_eq!(rs.data(RelationId(2)), &30);
    }

    #[rstest]
    fn test_var_relation_set_participants_ordered() {
        let rs: VarRelationSet<NodeId, ()> = VarRelationSet::new(vec![
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
        let rs: VarRelationSet<NodeId, &str> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], "pair"),
            (vec![n(2), n(3), n(4), n(5)], "quad"),
        ]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(2), n(3), n(4), n(5)]);
    }

    #[rstest]
    #[case::second_relation(RelationId(1), vec![ParticipantPosition(2), ParticipantPosition(0), ParticipantPosition(1)],
        vec![n(4), n(2), n(3)])]
    #[case::last_relation(RelationId(2), vec![ParticipantPosition(1), ParticipantPosition(0)], vec![n(6), n(5)])]
    fn test_var_relation_set_permute_with(
        #[case] id: RelationId,
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: Vec<NodeId>,
    ) {
        let mut rs: VarRelationSet<NodeId, &str> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], "a"),
            (vec![n(2), n(3), n(4)], "b"),
            (vec![n(5), n(6)], "c"),
        ]);
        let incidence_before: Vec<Vec<RelationId>> =
            (0..7).map(|i| rs.incident(n(i)).to_vec()).collect();

        rs.permute_with(id, &order);

        assert_eq!(rs.participants(id), expected.as_slice());
        assert_eq!(rs.data(id), &["a", "b", "c"][id.index()]);
        for other in [RelationId(0), RelationId(1), RelationId(2)] {
            if other != id {
                let stored: Vec<Vec<NodeId>> =
                    vec![vec![n(0), n(1)], vec![n(2), n(3), n(4)], vec![n(5), n(6)]];
                assert_eq!(rs.participants(other), stored[other.index()].as_slice());
            }
        }
        let incidence_after: Vec<Vec<RelationId>> =
            (0..7).map(|i| rs.incident(n(i)).to_vec()).collect();
        assert_eq!(incidence_after, incidence_before);
    }

    #[rstest]
    fn test_var_relation_set_permute_with_identity() {
        let input: VarRelationSet<NodeId, &str> =
            VarRelationSet::new(vec![(vec![n(2), n(0), n(1)], "a"), (vec![n(4), n(3)], "b")]);
        let mut permuted = input.clone();
        permuted.permute_with(
            RelationId(0),
            &[
                ParticipantPosition(0),
                ParticipantPosition(1),
                ParticipantPosition(2),
            ],
        );
        assert_eq!(permuted, input);
    }

    #[rstest]
    #[case::order_too_short(vec![ParticipantPosition(0), ParticipantPosition(1)])]
    #[case::position_out_of_range(vec![ParticipantPosition(0), ParticipantPosition(1), ParticipantPosition(9)])]
    #[case::position_repeated(vec![ParticipantPosition(1), ParticipantPosition(1), ParticipantPosition(0)])]
    #[should_panic(expected = "permute")]
    fn test_var_relation_set_permute_with_error(#[case] order: Vec<ParticipantPosition>) {
        let mut rs: VarRelationSet<NodeId, &str> =
            VarRelationSet::new(vec![(vec![n(0), n(1)], "a"), (vec![n(2), n(3), n(4)], "b")]);
        rs.permute_with(RelationId(1), &order);
    }

    #[rstest]
    fn test_var_relation_set_incidence() {
        let rs: VarRelationSet<NodeId, ()> = VarRelationSet::new(vec![
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
        let rs: VarRelationSet<EdgeId, &str> = VarRelationSet::new(vec![
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
        let rs: VarRelationSet<NodeId, ()> = VarRelationSet::new(vec![(vec![n(0), n(1)], ())]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_var_relation_set_relation_ids() {
        assert_exact_size(VarRelationSet::<NodeId, ()>::default().ids(), vec![]);
        let rs: VarRelationSet<NodeId, ()> =
            VarRelationSet::new(vec![(vec![n(0), n(1)], ()), (vec![n(1), n(2)], ())]);
        assert_exact_size(rs.ids(), vec![RelationId(0), RelationId(1)]);
    }

    #[fixture]
    fn var_relation_set_compaction_input() -> VarRelationSet<NodeId, &'static str> {
        VarRelationSet::new(vec![
            (vec![NodeId(0), NodeId(2), NodeId(4)], "keep"),
            (vec![NodeId(1), NodeId(3)], "drop"),
        ])
    }

    #[rstest]
    #[case::partial(
        vec![NodeId(1)],
        VarRelationSet::new(vec![(vec![NodeId(0), NodeId(1), NodeId(3)], "keep")]),
        vec![RelationId(1)],
    )]
    #[case::all(
        vec![NodeId(0), NodeId(1)],
        VarRelationSet::default(),
        vec![RelationId(0), RelationId(1)],
    )]
    fn test_var_relation_set_tracked_compact(
        var_relation_set_compaction_input: VarRelationSet<NodeId, &'static str>,
        #[case] removed_nodes: Vec<NodeId>,
        #[case] expected: VarRelationSet<NodeId, &'static str>,
        #[case] removed_relations: Vec<RelationId>,
    ) {
        let input = var_relation_set_compaction_input;
        let compaction = GraphCompaction::new(
            Compaction::new(5, removed_nodes).unwrap(),
            Compaction::identity(0),
        );
        let (output, witness) = input.tracked_compact(&compaction);
        assert_eq!(input.compact(&compaction), expected);
        assert_eq!(output, expected);
        assert_eq!(
            witness,
            Compaction::new(2, removed_relations.clone()).unwrap()
        );
        let survivors = (0..2)
            .map(RelationId)
            .filter(|id| !removed_relations.contains(id))
            .collect::<Vec<_>>();
        for (idx, &old) in survivors.iter().enumerate() {
            assert_eq!(witness.compact(old), Some(RelationId::from(idx)));
        }
    }

    #[rstest]
    #[case::empty(VarRelationSet::default())]
    #[case::rows(
        VarRelationSet::new(vec![(vec![NodeId(0), NodeId(2), NodeId(4)], "keep"), (vec![NodeId(1), NodeId(3)], "drop")]),
    )]
    fn test_var_relation_set_compact_identity(#[case] input: VarRelationSet<NodeId, &'static str>) {
        let compaction = GraphCompaction::new(Compaction::identity(5), Compaction::identity(0));
        assert_eq!(input.compact(&compaction), input);
        assert_eq!(
            input.tracked_compact(&compaction),
            (input.clone(), Compaction::identity(input.count())),
        );
    }

    #[rstest]
    #[case::rows(VarRelationSet::new(vec![(vec![NodeId(2), NodeId(0)], vec![7, 11]), (vec![NodeId(2), NodeId(0)], vec![13, 17])]),
        VarRelationSet::new(vec![(vec![NodeId(1), NodeId(5)], vec![7, 11]), (vec![NodeId(1), NodeId(5)], vec![13, 17])]))]
    fn test_var_relation_set_map(
        participant_correspondence: GraphCorrespondence,
        #[case] input: VarRelationSet<NodeId, Vec<u32>>,
        #[case] expected: VarRelationSet<NodeId, Vec<u32>>,
    ) {
        assert_eq!(input.map(&participant_correspondence), expected);
        assert_eq!(
            input.try_map(&participant_correspondence),
            Some(expected.clone())
        );
        let reverse = GraphCorrespondence::new(
            participant_correspondence.nodes().reverse(),
            participant_correspondence.edges().reverse(),
        );
        assert_eq!(expected.map(&reverse), input);
        let composed = participant_correspondence.compose(&reverse).unwrap();
        assert_eq!(input.map(&composed), input);
        assert_eq!(
            expected.incident(NodeId(1)),
            &[RelationId(0), RelationId(1)]
        );
    }

    #[rstest]
    #[case::empty(VarRelationSet::new(vec![]))]
    #[case::rows(VarRelationSet::new(vec![(vec![NodeId(2), NodeId(0)], vec![7, 11]), (vec![NodeId(2), NodeId(0)], vec![13, 17])]))]
    fn test_var_relation_set_map_identity(#[case] input: VarRelationSet<NodeId, Vec<u32>>) {
        let identity = GraphCorrespondence::new(
            Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2), NodeId(3)], 4),
            Correspondence::from_images(&[EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)], 4),
        );
        assert_eq!(input.try_map(&identity), Some(input.clone()));
        assert_eq!(input.map(&identity), input);
    }

    #[rstest]
    #[case::missing_node(1)]
    #[case::outside_node(4)]
    fn test_var_relation_set_try_map_error(
        participant_correspondence: GraphCorrespondence,
        #[case] node: u32,
    ) {
        let input: VarRelationSet<NodeId, Vec<u32>> = VarRelationSet::new(vec![
            (vec![NodeId(2), NodeId(0)], vec![7, 11]),
            (vec![NodeId(node), NodeId(0)], vec![13, 17]),
        ]);
        assert_eq!(input.try_map(&participant_correspondence), None);
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_var_relation_set_map_error(participant_correspondence: GraphCorrespondence) {
        let node = 1;

        let input: VarRelationSet<NodeId, Vec<u32>> =
            VarRelationSet::new(vec![(vec![NodeId(node), NodeId(0)], vec![7, 11])]);
        input.map(&participant_correspondence);
    }

    #[rstest]
    fn test_var_relation_set_remap() {
        let rs: VarRelationSet<EdgeId, PositionLabels> = VarRelationSet::new(vec![(
            vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            PositionLabels(vec![20, 21, 22]),
        )]);
        let remapping = GraphRemapping::new(
            Remapping::new(vec![]).expect("permutation images"),
            Remapping::new(vec![EdgeId(2), EdgeId(0), EdgeId(1)]).expect("permutation images"),
        );
        let out = rs.remap(&remapping);
        assert_eq!(
            out.participants(RelationId(0)),
            &[EdgeId(2), EdgeId(0), EdgeId(1)]
        );
        assert_eq!(out.data(RelationId(0)), &PositionLabels(vec![20, 21, 22]));
    }

    #[rstest]
    #[case::covered(vec![EdgeId(2), EdgeId(0), EdgeId(1)], true)]
    #[case::uncovered_edge(vec![EdgeId(1), EdgeId(0)], false)]
    fn test_var_relation_set_try_remap(#[case] edges: Vec<EdgeId>, #[case] covered: bool) {
        let rs: VarRelationSet<EdgeId, PositionLabels> = VarRelationSet::new(vec![(
            vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            PositionLabels(vec![20, 21, 22]),
        )]);
        let remapping = GraphRemapping::new(
            Remapping::new(vec![]).expect("permutation images"),
            Remapping::new(edges).expect("permutation images"),
        );
        let expected = covered.then(|| rs.remap(&remapping));
        assert_eq!(rs.try_remap(&remapping), expected);
    }

    #[rstest]
    fn test_var_relation_set_default() {
        let rs = VarRelationSet::<NodeId, ()>::default();
        assert_eq!(rs.count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_new() {
        let rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 2, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(2), n(1)], "x")]);
        assert_eq!(rs.count(), 1);
        assert_eq!(rs.participants_1(RelationId(0)), &[n(0)]);
        assert_eq!(rs.participants_2(RelationId(0)), &[n(2), n(1)]);
        assert_eq!(rs.data(RelationId(0)), &"x");
    }

    #[rstest]
    #[case::roundtrip(
        vec![
            ([n(2)], [n(4), n(1)], "first"),
            ([n(3)], [n(5), n(0)], "second"),
        ],
    )]
    fn test_fixed_fixed_birelation_set_into_entries(
        #[case] entries: Vec<([NodeId; 1], [NodeId; 2], &str)>,
    ) {
        let rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 2, &str> =
            FixedFixedBirelationSet::new(entries.clone());
        assert_eq!(rs.into_entries(), entries);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_data_mut() {
        let mut rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, i32> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_iter() {
        let empty = FixedFixedBirelationSet::<NodeId, 1, NodeId, 1, i32>::default();
        assert_eq!(empty.iter().collect::<Vec<_>>(), vec![]);

        let rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, i32> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], 1), ([n(2)], [n(3)], 2)]);
        assert_eq!(rs.iter().len(), 2);
        assert_eq!(
            rs.iter().collect::<Vec<_>>(),
            vec![
                (RelationId(0), &[n(0)], &[n(1)], &1),
                (RelationId(1), &[n(2)], &[n(3)], &2),
            ],
        );
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_iter_mut() {
        let mut empty = FixedFixedBirelationSet::<NodeId, 1, NodeId, 1, i32>::default();
        assert_eq!(empty.iter_mut().len(), 0);

        let mut rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, i32> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], 1), ([n(2)], [n(3)], 2)]);
        for (_, first, second, data) in rs.iter_mut() {
            assert_eq!(second[0].0, first[0].0 + 1);
            *data *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_permute_1_with() {
        let mut rs: FixedFixedBirelationSet<NodeId, 3, EdgeId, 2, &str> =
            FixedFixedBirelationSet::new(vec![([n(0), n(1), n(2)], [EdgeId(7), EdgeId(8)], "a")]);
        rs.permute_1_with(
            RelationId(0),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );
        assert_eq!(rs.participants_1(RelationId(0)), &[n(2), n(0), n(1)]);
        assert_eq!(rs.participants_2(RelationId(0)), &[EdgeId(7), EdgeId(8)]);
        assert_eq!(rs.data(RelationId(0)), &"a");
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_permute_2_with() {
        let mut rs: FixedFixedBirelationSet<NodeId, 3, EdgeId, 2, &str> =
            FixedFixedBirelationSet::new(vec![([n(0), n(1), n(2)], [EdgeId(7), EdgeId(8)], "a")]);
        rs.permute_2_with(
            RelationId(0),
            &[ParticipantPosition(1), ParticipantPosition(0)],
        );
        assert_eq!(rs.participants_1(RelationId(0)), &[n(0), n(1), n(2)]);
        assert_eq!(rs.participants_2(RelationId(0)), &[EdgeId(8), EdgeId(7)]);
        assert_eq!(rs.data(RelationId(0)), &"a");
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_incidence() {
        let rs: FixedFixedBirelationSet<NodeId, 1, EdgeId, 1, &str> =
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
        let rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], "x")]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_relation_ids() {
        assert_exact_size(
            FixedFixedBirelationSet::<NodeId, 1, NodeId, 1, &str>::default().ids(),
            vec![],
        );
        let rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], "a"), ([n(2)], [n(3)], "b")]);
        assert_exact_size(rs.ids(), vec![RelationId(0), RelationId(1)]);
    }

    #[fixture]
    fn fixed_fixed_birelation_set_compaction_input(
    ) -> FixedFixedBirelationSet<NodeId, 1, NodeId, 2, &'static str> {
        FixedFixedBirelationSet::new(vec![
            ([NodeId(0)], [NodeId(2), NodeId(4)], "keep"),
            ([NodeId(1)], [NodeId(5), NodeId(6)], "drop"),
        ])
    }

    #[rstest]
    #[case::partial(
        vec![NodeId(1)],
        FixedFixedBirelationSet::new(vec![([NodeId(0)], [NodeId(1), NodeId(3)], "keep")]),
        vec![RelationId(1)],
    )]
    #[case::all(
        vec![NodeId(0), NodeId(1)],
        FixedFixedBirelationSet::default(),
        vec![RelationId(0), RelationId(1)],
    )]
    fn test_fixed_fixed_birelation_set_tracked_compact(
        fixed_fixed_birelation_set_compaction_input: FixedFixedBirelationSet<
            NodeId,
            1,
            NodeId,
            2,
            &'static str,
        >,
        #[case] removed_nodes: Vec<NodeId>,
        #[case] expected: FixedFixedBirelationSet<NodeId, 1, NodeId, 2, &'static str>,
        #[case] removed_relations: Vec<RelationId>,
    ) {
        let input = fixed_fixed_birelation_set_compaction_input;
        let compaction = GraphCompaction::new(
            Compaction::new(7, removed_nodes).unwrap(),
            Compaction::identity(0),
        );
        let (output, witness) = input.tracked_compact(&compaction);
        assert_eq!(input.compact(&compaction), expected);
        assert_eq!(output, expected);
        assert_eq!(
            witness,
            Compaction::new(2, removed_relations.clone()).unwrap()
        );
        let survivors = (0..2)
            .map(RelationId)
            .filter(|id| !removed_relations.contains(id))
            .collect::<Vec<_>>();
        for (idx, &old) in survivors.iter().enumerate() {
            assert_eq!(witness.compact(old), Some(RelationId::from(idx)));
        }
    }

    #[rstest]
    #[case::empty(FixedFixedBirelationSet::default())]
    #[case::rows(
        FixedFixedBirelationSet::new(vec![([NodeId(0)], [NodeId(2), NodeId(4)], "keep"), ([NodeId(1)], [NodeId(5), NodeId(6)], "drop")]),
    )]
    fn test_fixed_fixed_birelation_set_compact_identity(
        #[case] input: FixedFixedBirelationSet<NodeId, 1, NodeId, 2, &'static str>,
    ) {
        let compaction = GraphCompaction::new(Compaction::identity(7), Compaction::identity(0));
        assert_eq!(input.compact(&compaction), input);
        assert_eq!(
            input.tracked_compact(&compaction),
            (input.clone(), Compaction::identity(input.count())),
        );
    }

    #[rstest]
    #[case::rows(FixedFixedBirelationSet::new(vec![([EdgeId(2), EdgeId(0)], [NodeId(2), NodeId(0)], vec![7, 11]), ([EdgeId(2), EdgeId(0)], [NodeId(2), NodeId(0)], vec![13, 17])]),
        FixedFixedBirelationSet::new(vec![([EdgeId(3), EdgeId(6)], [NodeId(1), NodeId(5)], vec![7, 11]), ([EdgeId(3), EdgeId(6)], [NodeId(1), NodeId(5)], vec![13, 17])]))]
    fn test_fixed_fixed_birelation_set_map(
        participant_correspondence: GraphCorrespondence,
        #[case] input: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, Vec<u32>>,
        #[case] expected: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, Vec<u32>>,
    ) {
        assert_eq!(input.map(&participant_correspondence), expected);
        assert_eq!(
            input.try_map(&participant_correspondence),
            Some(expected.clone())
        );
        let reverse = GraphCorrespondence::new(
            participant_correspondence.nodes().reverse(),
            participant_correspondence.edges().reverse(),
        );
        assert_eq!(expected.map(&reverse), input);
        let composed = participant_correspondence.compose(&reverse).unwrap();
        assert_eq!(input.map(&composed), input);
        assert_eq!(
            expected.incident(NodeId(1)),
            &[RelationId(0), RelationId(1)]
        );
        assert_eq!(
            expected.incident_edge(EdgeId(3)),
            expected.incident(NodeId(1))
        );
    }

    #[rstest]
    #[case::empty(FixedFixedBirelationSet::new(vec![]))]
    #[case::rows(FixedFixedBirelationSet::new(vec![([EdgeId(2), EdgeId(0)], [NodeId(2), NodeId(0)], vec![7, 11]), ([EdgeId(2), EdgeId(0)], [NodeId(2), NodeId(0)], vec![13, 17])]))]
    fn test_fixed_fixed_birelation_set_map_identity(
        #[case] input: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, Vec<u32>>,
    ) {
        let identity = GraphCorrespondence::new(
            Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2), NodeId(3)], 4),
            Correspondence::from_images(&[EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)], 4),
        );
        assert_eq!(input.try_map(&identity), Some(input.clone()));
        assert_eq!(input.map(&identity), input);
    }

    #[rstest]
    #[case::missing_node(1, 2)]
    #[case::outside_node(4, 2)]
    #[case::missing_edge(2, 1)]
    #[case::outside_edge(2, 4)]
    fn test_fixed_fixed_birelation_set_try_map_error(
        participant_correspondence: GraphCorrespondence,
        #[case] node: u32,
        #[case] edge: u32,
    ) {
        let input: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, Vec<u32>> =
            FixedFixedBirelationSet::new(vec![
                ([EdgeId(2), EdgeId(0)], [NodeId(2), NodeId(0)], vec![7, 11]),
                (
                    [EdgeId(edge), EdgeId(0)],
                    [NodeId(node), NodeId(0)],
                    vec![13, 17],
                ),
            ]);
        assert_eq!(input.try_map(&participant_correspondence), None);
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_fixed_fixed_birelation_set_map_error(participant_correspondence: GraphCorrespondence) {
        let node = 1;
        let edge = 2;
        let input: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, Vec<u32>> =
            FixedFixedBirelationSet::new(vec![(
                [EdgeId(edge), EdgeId(0)],
                [NodeId(node), NodeId(0)],
                vec![7, 11],
            )]);
        input.map(&participant_correspondence);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_remap() {
        let rs: FixedFixedBirelationSet<NodeId, 2, EdgeId, 2, BiPositionLabels> =
            FixedFixedBirelationSet::new(vec![(
                [n(0), n(1)],
                [EdgeId(0), EdgeId(1)],
                BiPositionLabels {
                    factor_1: vec![10, 11],
                    factor_2: vec![20, 21],
                },
            )]);
        let remapping = GraphRemapping::new(
            Remapping::new(vec![n(1), n(0)]).expect("permutation images"),
            Remapping::new(vec![EdgeId(1), EdgeId(0)]).expect("permutation images"),
        );
        let out = rs.remap(&remapping);
        assert_eq!(out.participants_1(RelationId(0)), &[n(1), n(0)]);
        assert_eq!(out.participants_2(RelationId(0)), &[EdgeId(1), EdgeId(0)]);
        assert_eq!(
            out.data(RelationId(0)),
            &BiPositionLabels {
                factor_1: vec![10, 11],
                factor_2: vec![20, 21],
            }
        );
    }

    #[rstest]
    #[case::covered(vec![n(1), n(0)], vec![EdgeId(1), EdgeId(0)], true)]
    #[case::uncovered_node(vec![n(0)], vec![EdgeId(1), EdgeId(0)], false)]
    #[case::uncovered_edge(vec![n(1), n(0)], vec![EdgeId(0)], false)]
    fn test_fixed_fixed_birelation_set_try_remap(
        #[case] nodes: Vec<NodeId>,
        #[case] edges: Vec<EdgeId>,
        #[case] covered: bool,
    ) {
        let rs: FixedFixedBirelationSet<NodeId, 2, EdgeId, 2, BiPositionLabels> =
            FixedFixedBirelationSet::new(vec![(
                [n(0), n(1)],
                [EdgeId(0), EdgeId(1)],
                BiPositionLabels {
                    factor_1: vec![10, 11],
                    factor_2: vec![20, 21],
                },
            )]);
        let remapping = GraphRemapping::new(
            Remapping::new(nodes).expect("permutation images"),
            Remapping::new(edges).expect("permutation images"),
        );
        let expected = covered.then(|| rs.remap(&remapping));
        assert_eq!(rs.try_remap(&remapping), expected);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_default() {
        let rs = FixedFixedBirelationSet::<NodeId, 1, NodeId, 1, ()>::default();
        assert_eq!(rs.count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_fixed_var_birelation_set_new() {
        let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1), n(2), n(3)], "ct")]);
        assert_eq!(rs.count(), 1);
        assert_eq!(rs.participants_1(RelationId(0)), &[EdgeId(0)]);
        assert_eq!(rs.participants_2(RelationId(0)), &[n(1), n(2), n(3)]);
        assert_eq!(rs.data(RelationId(0)), &"ct");
    }

    #[rstest]
    fn test_fixed_var_birelation_set_hash() {
        let entries = vec![
            ([EdgeId(2)], vec![n(3), n(1)], "first"),
            ([EdgeId(4)], vec![n(5), n(0)], "second"),
        ];
        let left: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> =
            FixedVarBirelationSet::new(entries.clone());
        let right: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> =
            FixedVarBirelationSet::new(entries);
        assert_eq!(left, right);
        assert_eq!(hash(&left), hash(&right));
    }

    #[rstest]
    #[case::roundtrip(
        vec![
            ([EdgeId(2)], vec![n(3), n(1)], "first"),
            ([EdgeId(4)], vec![n(5), n(0)], "second"),
        ],
    )]
    fn test_fixed_var_birelation_set_into_entries(
        #[case] entries: Vec<([EdgeId; 1], Vec<NodeId>, &str)>,
    ) {
        let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> =
            FixedVarBirelationSet::new(entries.clone());
        assert_eq!(rs.into_entries(), entries);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_data_mut() {
        let mut rs: FixedVarBirelationSet<EdgeId, 1, NodeId, i32> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_iter() {
        let empty = FixedVarBirelationSet::<EdgeId, 1, NodeId, i32>::default();
        assert_eq!(empty.iter().collect::<Vec<_>>(), vec![]);

        let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, i32> = FixedVarBirelationSet::new(vec![
            ([EdgeId(0)], vec![n(1), n(3)], 1),
            ([EdgeId(1)], vec![n(2)], 2),
        ]);
        assert_eq!(rs.iter().len(), 2);
        assert_eq!(
            rs.iter().collect::<Vec<_>>(),
            vec![
                (RelationId(0), &[EdgeId(0)], [n(1), n(3)].as_slice(), &1),
                (RelationId(1), &[EdgeId(1)], [n(2)].as_slice(), &2),
            ],
        );
    }

    #[rstest]
    fn test_fixed_var_birelation_set_iter_mut() {
        let mut empty = FixedVarBirelationSet::<EdgeId, 1, NodeId, i32>::default();
        assert_eq!(empty.iter_mut().len(), 0);

        let mut rs: FixedVarBirelationSet<EdgeId, 1, NodeId, i32> =
            FixedVarBirelationSet::new(vec![
                ([EdgeId(0)], vec![n(1), n(3)], 1),
                ([EdgeId(1)], vec![n(2)], 2),
            ]);
        let arities: Vec<usize> = rs
            .iter_mut()
            .map(|(_, _, ligands, data)| {
                *data *= 10;
                ligands.len()
            })
            .collect();
        assert_eq!(arities, vec![2, 1]);
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_permute_1_with() {
        let mut rs: FixedVarBirelationSet<EdgeId, 2, NodeId, &str> =
            FixedVarBirelationSet::new(vec![
                ([EdgeId(4), EdgeId(5)], vec![n(0), n(1), n(2)], "a"),
                ([EdgeId(6), EdgeId(7)], vec![n(3), n(4)], "b"),
            ]);
        rs.permute_1_with(
            RelationId(1),
            &[ParticipantPosition(1), ParticipantPosition(0)],
        );
        assert_eq!(rs.participants_1(RelationId(0)), &[EdgeId(4), EdgeId(5)]);
        assert_eq!(rs.participants_1(RelationId(1)), &[EdgeId(7), EdgeId(6)]);
        assert_eq!(rs.participants_2(RelationId(1)), &[n(3), n(4)]);
        assert_eq!(rs.data(RelationId(1)), &"b");
    }

    #[rstest]
    fn test_fixed_var_birelation_set_permute_2_with() {
        let mut rs: FixedVarBirelationSet<EdgeId, 2, NodeId, &str> =
            FixedVarBirelationSet::new(vec![
                ([EdgeId(4), EdgeId(5)], vec![n(0), n(1), n(2)], "a"),
                ([EdgeId(6), EdgeId(7)], vec![n(3), n(4)], "b"),
            ]);
        let incidence_before: Vec<Vec<RelationId>> =
            (0..5).map(|i| rs.incident(n(i)).to_vec()).collect();

        rs.permute_2_with(
            RelationId(0),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );

        assert_eq!(rs.participants_2(RelationId(0)), &[n(2), n(0), n(1)]);
        assert_eq!(rs.participants_2(RelationId(1)), &[n(3), n(4)]);
        assert_eq!(rs.participants_1(RelationId(0)), &[EdgeId(4), EdgeId(5)]);
        assert_eq!(rs.data(RelationId(0)), &"a");
        let incidence_after: Vec<Vec<RelationId>> =
            (0..5).map(|i| rs.incident(n(i)).to_vec()).collect();
        assert_eq!(incidence_after, incidence_before);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_permute_1_with_identity() {
        let input: FixedVarBirelationSet<EdgeId, 2, NodeId, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(5), EdgeId(4)], vec![n(2), n(0), n(1)], "a")]);
        let mut permuted = input.clone();
        permuted.permute_1_with(
            RelationId(0),
            &[ParticipantPosition(0), ParticipantPosition(1)],
        );
        permuted.permute_2_with(
            RelationId(0),
            &[
                ParticipantPosition(0),
                ParticipantPosition(1),
                ParticipantPosition(2),
            ],
        );
        assert_eq!(permuted, input);
    }

    #[rstest]
    #[case::order_too_short(vec![ParticipantPosition(0), ParticipantPosition(1)])]
    #[case::position_out_of_range(vec![ParticipantPosition(0), ParticipantPosition(1), ParticipantPosition(3)])]
    #[case::position_repeated(vec![ParticipantPosition(2), ParticipantPosition(2), ParticipantPosition(0)])]
    #[should_panic(expected = "permute")]
    fn test_fixed_var_birelation_set_permute_2_with_error(#[case] order: Vec<ParticipantPosition>) {
        let mut rs: FixedVarBirelationSet<EdgeId, 2, NodeId, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(4), EdgeId(5)], vec![n(0), n(1), n(2)], "a")]);
        rs.permute_2_with(RelationId(0), &order);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_incidence() {
        let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> =
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
        let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1)], "ct")]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_relation_ids() {
        assert_exact_size(
            FixedVarBirelationSet::<EdgeId, 1, NodeId, &str>::default().ids(),
            vec![],
        );
        let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> = FixedVarBirelationSet::new(vec![
            ([EdgeId(0)], vec![n(1)], "a"),
            ([EdgeId(1)], vec![n(2)], "b"),
        ]);
        assert_exact_size(rs.ids(), vec![RelationId(0), RelationId(1)]);
    }

    #[fixture]
    fn fixed_var_birelation_set_compaction_input(
    ) -> FixedVarBirelationSet<NodeId, 1, NodeId, &'static str> {
        FixedVarBirelationSet::new(vec![
            ([NodeId(0)], vec![NodeId(2), NodeId(4)], "keep"),
            ([NodeId(5)], vec![NodeId(1), NodeId(3)], "drop"),
        ])
    }

    #[rstest]
    #[case::partial(
        vec![NodeId(1)],
        FixedVarBirelationSet::new(vec![([NodeId(0)], vec![NodeId(1), NodeId(3)], "keep")]),
        vec![RelationId(1)],
    )]
    #[case::all(
        vec![NodeId(0), NodeId(1)],
        FixedVarBirelationSet::default(),
        vec![RelationId(0), RelationId(1)],
    )]
    fn test_fixed_var_birelation_set_tracked_compact(
        fixed_var_birelation_set_compaction_input: FixedVarBirelationSet<
            NodeId,
            1,
            NodeId,
            &'static str,
        >,
        #[case] removed_nodes: Vec<NodeId>,
        #[case] expected: FixedVarBirelationSet<NodeId, 1, NodeId, &'static str>,
        #[case] removed_relations: Vec<RelationId>,
    ) {
        let input = fixed_var_birelation_set_compaction_input;
        let compaction = GraphCompaction::new(
            Compaction::new(6, removed_nodes).unwrap(),
            Compaction::identity(0),
        );
        let (output, witness) = input.tracked_compact(&compaction);
        assert_eq!(input.compact(&compaction), expected);
        assert_eq!(output, expected);
        assert_eq!(
            witness,
            Compaction::new(2, removed_relations.clone()).unwrap()
        );
        let survivors = (0..2)
            .map(RelationId)
            .filter(|id| !removed_relations.contains(id))
            .collect::<Vec<_>>();
        for (idx, &old) in survivors.iter().enumerate() {
            assert_eq!(witness.compact(old), Some(RelationId::from(idx)));
        }
    }

    #[rstest]
    #[case::empty(FixedVarBirelationSet::default())]
    #[case::rows(
        FixedVarBirelationSet::new(vec![([NodeId(0)], vec![NodeId(2), NodeId(4)], "keep"), ([NodeId(5)], vec![NodeId(1), NodeId(3)], "drop")]),
    )]
    fn test_fixed_var_birelation_set_compact_identity(
        #[case] input: FixedVarBirelationSet<NodeId, 1, NodeId, &'static str>,
    ) {
        let compaction = GraphCompaction::new(Compaction::identity(6), Compaction::identity(0));
        assert_eq!(input.compact(&compaction), input);
        assert_eq!(
            input.tracked_compact(&compaction),
            (input.clone(), Compaction::identity(input.count())),
        );
    }

    #[rstest]
    #[case::rows(FixedVarBirelationSet::new(vec![([EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![7, 11]), ([EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![13, 17])]),
        FixedVarBirelationSet::new(vec![([EdgeId(3), EdgeId(6)], vec![NodeId(1), NodeId(5)], vec![7, 11]), ([EdgeId(3), EdgeId(6)], vec![NodeId(1), NodeId(5)], vec![13, 17])]))]
    fn test_fixed_var_birelation_set_map(
        participant_correspondence: GraphCorrespondence,
        #[case] input: FixedVarBirelationSet<EdgeId, 2, NodeId, Vec<u32>>,
        #[case] expected: FixedVarBirelationSet<EdgeId, 2, NodeId, Vec<u32>>,
    ) {
        assert_eq!(input.map(&participant_correspondence), expected);
        assert_eq!(
            input.try_map(&participant_correspondence),
            Some(expected.clone())
        );
        let reverse = GraphCorrespondence::new(
            participant_correspondence.nodes().reverse(),
            participant_correspondence.edges().reverse(),
        );
        assert_eq!(expected.map(&reverse), input);
        let composed = participant_correspondence.compose(&reverse).unwrap();
        assert_eq!(input.map(&composed), input);
        assert_eq!(
            expected.incident(NodeId(1)),
            &[RelationId(0), RelationId(1)]
        );
        assert_eq!(
            expected.incident_edge(EdgeId(3)),
            expected.incident(NodeId(1))
        );
    }

    #[rstest]
    #[case::empty(FixedVarBirelationSet::new(vec![]))]
    #[case::rows(FixedVarBirelationSet::new(vec![([EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![7, 11]), ([EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![13, 17])]))]
    fn test_fixed_var_birelation_set_map_identity(
        #[case] input: FixedVarBirelationSet<EdgeId, 2, NodeId, Vec<u32>>,
    ) {
        let identity = GraphCorrespondence::new(
            Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2), NodeId(3)], 4),
            Correspondence::from_images(&[EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)], 4),
        );
        assert_eq!(input.try_map(&identity), Some(input.clone()));
        assert_eq!(input.map(&identity), input);
    }

    #[rstest]
    #[case::missing_node(1, 2)]
    #[case::outside_node(4, 2)]
    #[case::missing_edge(2, 1)]
    #[case::outside_edge(2, 4)]
    fn test_fixed_var_birelation_set_try_map_error(
        participant_correspondence: GraphCorrespondence,
        #[case] node: u32,
        #[case] edge: u32,
    ) {
        let input: FixedVarBirelationSet<EdgeId, 2, NodeId, Vec<u32>> =
            FixedVarBirelationSet::new(vec![
                (
                    [EdgeId(2), EdgeId(0)],
                    vec![NodeId(2), NodeId(0)],
                    vec![7, 11],
                ),
                (
                    [EdgeId(edge), EdgeId(0)],
                    vec![NodeId(node), NodeId(0)],
                    vec![13, 17],
                ),
            ]);
        assert_eq!(input.try_map(&participant_correspondence), None);
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_fixed_var_birelation_set_map_error(participant_correspondence: GraphCorrespondence) {
        let node = 1;
        let edge = 2;
        let input: FixedVarBirelationSet<EdgeId, 2, NodeId, Vec<u32>> =
            FixedVarBirelationSet::new(vec![(
                [EdgeId(edge), EdgeId(0)],
                vec![NodeId(node), NodeId(0)],
                vec![7, 11],
            )]);
        input.map(&participant_correspondence);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_remap() {
        let rs: FixedVarBirelationSet<EdgeId, 2, NodeId, BiPositionLabels> =
            FixedVarBirelationSet::new(vec![(
                [EdgeId(0), EdgeId(1)],
                vec![n(0), n(1), n(2)],
                BiPositionLabels {
                    factor_1: vec![30, 31],
                    factor_2: vec![40, 41, 42],
                },
            )]);
        let remapping = GraphRemapping::new(
            Remapping::new(vec![n(2), n(0), n(1)]).expect("permutation images"),
            Remapping::new(vec![EdgeId(2), EdgeId(0), EdgeId(1)]).expect("permutation images"),
        );
        let out = rs.remap(&remapping);
        assert_eq!(out.participants_1(RelationId(0)), &[EdgeId(2), EdgeId(0)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(2), n(0), n(1)]);
        assert_eq!(
            out.data(RelationId(0)),
            &BiPositionLabels {
                factor_1: vec![30, 31],
                factor_2: vec![40, 41, 42],
            }
        );
    }

    #[rstest]
    #[case::forward(vec![NodeId(0), NodeId(1)], vec![NodeId(1), NodeId(2)])]
    #[case::reversed(vec![NodeId(1), NodeId(0)], vec![NodeId(2), NodeId(1)])]
    fn test_fixed_var_birelation_set_map_pushout(
        #[case] participants: Vec<NodeId>,
        #[case] expected_participants: Vec<NodeId>,
    ) {
        let left = Graph::new(2, &[[0, 1]]);
        let right = Graph::new(2, &[[0, 1]]);
        let overlap = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(1), NodeId(0))], 2, 2).unwrap(),
            Correspondence::new(vec![], 1, 1).unwrap(),
        );
        let pushout = left.pushout(&right, &overlap);

        let relations: FixedVarBirelationSet<EdgeId, 1, NodeId, Vec<u32>> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], participants, vec![7, 11])]);
        let expected =
            FixedVarBirelationSet::new(vec![([EdgeId(1)], expected_participants, vec![7, 11])]);
        assert_eq!(relations.map(&pushout.right), expected);
    }

    #[rstest]
    #[case::covered(
        vec![n(2), n(0), n(1)],
        vec![EdgeId(2), EdgeId(0), EdgeId(1)],
        true,
    )]
    #[case::uncovered_node(vec![n(1), n(0)], vec![EdgeId(1), EdgeId(0)], false)]
    #[case::uncovered_edge(vec![n(2), n(0), n(1)], vec![EdgeId(0)], false)]
    fn test_fixed_var_birelation_set_try_remap(
        #[case] nodes: Vec<NodeId>,
        #[case] edges: Vec<EdgeId>,
        #[case] covered: bool,
    ) {
        let rs: FixedVarBirelationSet<EdgeId, 2, NodeId, BiPositionLabels> =
            FixedVarBirelationSet::new(vec![(
                [EdgeId(0), EdgeId(1)],
                vec![n(0), n(1), n(2)],
                BiPositionLabels {
                    factor_1: vec![30, 31],
                    factor_2: vec![40, 41, 42],
                },
            )]);
        let remapping = GraphRemapping::new(
            Remapping::new(nodes).expect("permutation images"),
            Remapping::new(edges).expect("permutation images"),
        );
        let expected = covered.then(|| rs.remap(&remapping));
        assert_eq!(rs.try_remap(&remapping), expected);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_default() {
        let rs = FixedVarBirelationSet::<EdgeId, 1, NodeId, ()>::default();
        assert_eq!(rs.count(), 0);
        assert!(!rs.has_incident_edge(EdgeId(0)));
    }

    #[rstest]
    fn test_var_var_birelation_set_new() {
        let rs: VarVarBirelationSet<NodeId, EdgeId, &str> =
            VarVarBirelationSet::new(vec![(vec![n(0), n(1)], vec![EdgeId(5)], "y")]);
        assert_eq!(rs.count(), 1);
        assert_eq!(rs.participants_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants_2(RelationId(0)), &[EdgeId(5)]);
        assert_eq!(rs.data(RelationId(0)), &"y");
    }

    #[rstest]
    #[case::roundtrip(
        vec![
            (vec![n(2), n(0)], vec![EdgeId(4), EdgeId(1)], "first"),
            (vec![n(5), n(3)], vec![EdgeId(6), EdgeId(2)], "second"),
        ],
    )]
    fn test_var_var_birelation_set_into_entries(
        #[case] entries: Vec<(Vec<NodeId>, Vec<EdgeId>, &str)>,
    ) {
        let rs: VarVarBirelationSet<NodeId, EdgeId, &str> =
            VarVarBirelationSet::new(entries.clone());
        assert_eq!(rs.into_entries(), entries);
    }

    #[rstest]
    fn test_var_var_birelation_set_data_mut() {
        let mut rs: VarVarBirelationSet<NodeId, EdgeId, i32> =
            VarVarBirelationSet::new(vec![(vec![n(0)], vec![EdgeId(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_var_var_birelation_set_iter() {
        let empty = VarVarBirelationSet::<NodeId, EdgeId, i32>::default();
        assert_eq!(empty.iter().collect::<Vec<_>>(), vec![]);

        let rs: VarVarBirelationSet<NodeId, EdgeId, i32> = VarVarBirelationSet::new(vec![
            (vec![n(0), n(4)], vec![EdgeId(1)], 1),
            (vec![n(2)], vec![EdgeId(3)], 2),
        ]);
        assert_eq!(rs.iter().len(), 2);
        assert_eq!(
            rs.iter().collect::<Vec<_>>(),
            vec![
                (
                    RelationId(0),
                    [n(0), n(4)].as_slice(),
                    [EdgeId(1)].as_slice(),
                    &1,
                ),
                (RelationId(1), [n(2)].as_slice(), [EdgeId(3)].as_slice(), &2),
            ],
        );
    }

    #[rstest]
    fn test_var_var_birelation_set_iter_mut() {
        let mut empty = VarVarBirelationSet::<NodeId, EdgeId, i32>::default();
        assert_eq!(empty.iter_mut().len(), 0);

        let mut rs: VarVarBirelationSet<NodeId, EdgeId, i32> = VarVarBirelationSet::new(vec![
            (vec![n(0), n(4)], vec![EdgeId(1)], 1),
            (vec![n(2)], vec![EdgeId(3)], 2),
        ]);
        let arities: Vec<(usize, usize)> = rs
            .iter_mut()
            .map(|(_, first, second, data)| {
                *data *= 10;
                (first.len(), second.len())
            })
            .collect();
        assert_eq!(arities, vec![(2, 1), (1, 1)]);
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_var_var_birelation_set_permute_1_with() {
        let mut rs: VarVarBirelationSet<NodeId, EdgeId, &str> = VarVarBirelationSet::new(vec![
            (vec![n(0), n(1)], vec![EdgeId(9)], "a"),
            (vec![n(2), n(3), n(4)], vec![EdgeId(7), EdgeId(8)], "b"),
        ]);
        rs.permute_1_with(
            RelationId(1),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );
        assert_eq!(rs.participants_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants_1(RelationId(1)), &[n(4), n(2), n(3)]);
        assert_eq!(rs.participants_2(RelationId(1)), &[EdgeId(7), EdgeId(8)]);
        assert_eq!(rs.data(RelationId(1)), &"b");
    }

    #[rstest]
    fn test_var_var_birelation_set_permute_2_with() {
        let mut rs: VarVarBirelationSet<NodeId, EdgeId, &str> = VarVarBirelationSet::new(vec![
            (vec![n(0), n(1)], vec![EdgeId(9)], "a"),
            (vec![n(2), n(3), n(4)], vec![EdgeId(7), EdgeId(8)], "b"),
        ]);
        rs.permute_2_with(
            RelationId(1),
            &[ParticipantPosition(1), ParticipantPosition(0)],
        );
        assert_eq!(rs.participants_2(RelationId(0)), &[EdgeId(9)]);
        assert_eq!(rs.participants_2(RelationId(1)), &[EdgeId(8), EdgeId(7)]);
        assert_eq!(rs.participants_1(RelationId(1)), &[n(2), n(3), n(4)]);
        assert_eq!(rs.data(RelationId(1)), &"b");
    }

    #[rstest]
    fn test_var_var_birelation_set_incidence() {
        let rs: VarVarBirelationSet<NodeId, EdgeId, &str> =
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
        let rs: VarVarBirelationSet<NodeId, EdgeId, &str> =
            VarVarBirelationSet::new(vec![(vec![n(0)], vec![EdgeId(1)], "y")]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_var_var_birelation_set_relation_ids() {
        assert_exact_size(
            VarVarBirelationSet::<NodeId, EdgeId, &str>::default().ids(),
            vec![],
        );
        let rs: VarVarBirelationSet<NodeId, EdgeId, &str> = VarVarBirelationSet::new(vec![
            (vec![n(0)], vec![EdgeId(1)], "a"),
            (vec![n(2)], vec![EdgeId(3)], "b"),
        ]);
        assert_exact_size(rs.ids(), vec![RelationId(0), RelationId(1)]);
    }

    #[fixture]
    fn var_var_birelation_set_compaction_input() -> VarVarBirelationSet<NodeId, NodeId, &'static str>
    {
        VarVarBirelationSet::new(vec![
            (vec![NodeId(0), NodeId(2)], vec![NodeId(4)], "keep"),
            (vec![NodeId(5)], vec![NodeId(1)], "drop"),
        ])
    }

    #[rstest]
    #[case::partial(
        vec![NodeId(1)],
        VarVarBirelationSet::new(vec![(vec![NodeId(0), NodeId(1)], vec![NodeId(3)], "keep")]),
        vec![RelationId(1)],
    )]
    #[case::all(
        vec![NodeId(0), NodeId(1)],
        VarVarBirelationSet::default(),
        vec![RelationId(0), RelationId(1)],
    )]
    fn test_var_var_birelation_set_tracked_compact(
        var_var_birelation_set_compaction_input: VarVarBirelationSet<NodeId, NodeId, &'static str>,
        #[case] removed_nodes: Vec<NodeId>,
        #[case] expected: VarVarBirelationSet<NodeId, NodeId, &'static str>,
        #[case] removed_relations: Vec<RelationId>,
    ) {
        let input = var_var_birelation_set_compaction_input;
        let compaction = GraphCompaction::new(
            Compaction::new(6, removed_nodes).unwrap(),
            Compaction::identity(0),
        );
        let (output, witness) = input.tracked_compact(&compaction);
        assert_eq!(input.compact(&compaction), expected);
        assert_eq!(output, expected);
        assert_eq!(
            witness,
            Compaction::new(2, removed_relations.clone()).unwrap()
        );
        let survivors = (0..2)
            .map(RelationId)
            .filter(|id| !removed_relations.contains(id))
            .collect::<Vec<_>>();
        for (idx, &old) in survivors.iter().enumerate() {
            assert_eq!(witness.compact(old), Some(RelationId::from(idx)));
        }
    }

    #[rstest]
    #[case::empty(VarVarBirelationSet::default())]
    #[case::rows(
        VarVarBirelationSet::new(vec![(vec![NodeId(0), NodeId(2)], vec![NodeId(4)], "keep"), (vec![NodeId(5)], vec![NodeId(1)], "drop")]),
    )]
    fn test_var_var_birelation_set_compact_identity(
        #[case] input: VarVarBirelationSet<NodeId, NodeId, &'static str>,
    ) {
        let compaction = GraphCompaction::new(Compaction::identity(6), Compaction::identity(0));
        assert_eq!(input.compact(&compaction), input);
        assert_eq!(
            input.tracked_compact(&compaction),
            (input.clone(), Compaction::identity(input.count())),
        );
    }

    #[rstest]
    #[case::rows(VarVarBirelationSet::new(vec![(vec![EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![7, 11]), (vec![EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![13, 17])]),
        VarVarBirelationSet::new(vec![(vec![EdgeId(3), EdgeId(6)], vec![NodeId(1), NodeId(5)], vec![7, 11]), (vec![EdgeId(3), EdgeId(6)], vec![NodeId(1), NodeId(5)], vec![13, 17])]))]
    fn test_var_var_birelation_set_map(
        participant_correspondence: GraphCorrespondence,
        #[case] input: VarVarBirelationSet<EdgeId, NodeId, Vec<u32>>,
        #[case] expected: VarVarBirelationSet<EdgeId, NodeId, Vec<u32>>,
    ) {
        assert_eq!(input.map(&participant_correspondence), expected);
        assert_eq!(
            input.try_map(&participant_correspondence),
            Some(expected.clone())
        );
        let reverse = GraphCorrespondence::new(
            participant_correspondence.nodes().reverse(),
            participant_correspondence.edges().reverse(),
        );
        assert_eq!(expected.map(&reverse), input);
        let composed = participant_correspondence.compose(&reverse).unwrap();
        assert_eq!(input.map(&composed), input);
        assert_eq!(
            expected.incident(NodeId(1)),
            &[RelationId(0), RelationId(1)]
        );
        assert_eq!(
            expected.incident_edge(EdgeId(3)),
            expected.incident(NodeId(1))
        );
    }

    #[rstest]
    #[case::empty(VarVarBirelationSet::new(vec![]))]
    #[case::rows(VarVarBirelationSet::new(vec![(vec![EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![7, 11]), (vec![EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![13, 17])]))]
    fn test_var_var_birelation_set_map_identity(
        #[case] input: VarVarBirelationSet<EdgeId, NodeId, Vec<u32>>,
    ) {
        let identity = GraphCorrespondence::new(
            Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2), NodeId(3)], 4),
            Correspondence::from_images(&[EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)], 4),
        );
        assert_eq!(input.try_map(&identity), Some(input.clone()));
        assert_eq!(input.map(&identity), input);
    }

    #[rstest]
    #[case::missing_node(1, 2)]
    #[case::outside_node(4, 2)]
    #[case::missing_edge(2, 1)]
    #[case::outside_edge(2, 4)]
    fn test_var_var_birelation_set_try_map_error(
        participant_correspondence: GraphCorrespondence,
        #[case] node: u32,
        #[case] edge: u32,
    ) {
        let input: VarVarBirelationSet<EdgeId, NodeId, Vec<u32>> = VarVarBirelationSet::new(vec![
            (
                vec![EdgeId(2), EdgeId(0)],
                vec![NodeId(2), NodeId(0)],
                vec![7, 11],
            ),
            (
                vec![EdgeId(edge), EdgeId(0)],
                vec![NodeId(node), NodeId(0)],
                vec![13, 17],
            ),
        ]);
        assert_eq!(input.try_map(&participant_correspondence), None);
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_var_var_birelation_set_map_error(participant_correspondence: GraphCorrespondence) {
        let node = 1;
        let edge = 2;
        let input: VarVarBirelationSet<EdgeId, NodeId, Vec<u32>> =
            VarVarBirelationSet::new(vec![(
                vec![EdgeId(edge), EdgeId(0)],
                vec![NodeId(node), NodeId(0)],
                vec![7, 11],
            )]);
        input.map(&participant_correspondence);
    }

    #[rstest]
    fn test_var_var_birelation_set_remap() {
        let rs: VarVarBirelationSet<NodeId, EdgeId, BiPositionLabels> =
            VarVarBirelationSet::new(vec![(
                vec![n(0), n(1)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                BiPositionLabels {
                    factor_1: vec![50, 51],
                    factor_2: vec![60, 61, 62],
                },
            )]);
        let remapping = GraphRemapping::new(
            Remapping::new(vec![n(1), n(0)]).expect("permutation images"),
            Remapping::new(vec![EdgeId(2), EdgeId(0), EdgeId(1)]).expect("permutation images"),
        );
        let out = rs.remap(&remapping);
        assert_eq!(out.participants_1(RelationId(0)), &[n(1), n(0)]);
        assert_eq!(
            out.participants_2(RelationId(0)),
            &[EdgeId(2), EdgeId(0), EdgeId(1)]
        );
        assert_eq!(
            out.data(RelationId(0)),
            &BiPositionLabels {
                factor_1: vec![50, 51],
                factor_2: vec![60, 61, 62],
            }
        );
    }

    #[rstest]
    #[case::covered(
        vec![n(1), n(0)],
        vec![EdgeId(2), EdgeId(0), EdgeId(1)],
        true,
    )]
    #[case::uncovered_node(vec![n(0)], vec![EdgeId(2), EdgeId(0), EdgeId(1)], false)]
    #[case::uncovered_edge(vec![n(1), n(0)], vec![EdgeId(1), EdgeId(0)], false)]
    fn test_var_var_birelation_set_try_remap(
        #[case] nodes: Vec<NodeId>,
        #[case] edges: Vec<EdgeId>,
        #[case] covered: bool,
    ) {
        let rs: VarVarBirelationSet<NodeId, EdgeId, BiPositionLabels> =
            VarVarBirelationSet::new(vec![(
                vec![n(0), n(1)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                BiPositionLabels {
                    factor_1: vec![50, 51],
                    factor_2: vec![60, 61, 62],
                },
            )]);
        let remapping = GraphRemapping::new(
            Remapping::new(nodes).expect("permutation images"),
            Remapping::new(edges).expect("permutation images"),
        );
        let expected = covered.then(|| rs.remap(&remapping));
        assert_eq!(rs.try_remap(&remapping), expected);
    }

    #[rstest]
    fn test_var_var_birelation_set_default() {
        let rs = VarVarBirelationSet::<NodeId, EdgeId, ()>::default();
        assert_eq!(rs.count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    #[case::exact(vec![n(0), n(1)], Some(RelationId(0)))]
    #[case::reordered(vec![n(1), n(0)], Some(RelationId(0)))]
    #[case::second(vec![n(2), n(3)], Some(RelationId(1)))]
    #[case::absent(vec![n(0), n(3)], None)]
    #[case::wrong_arity(vec![n(0)], None)]
    fn test_fixed_relation_set_coincident(
        #[case] query: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        let rs: FixedRelationSet<NodeId, (), 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], ()), ([n(2), n(3)], ())]);
        assert_eq!(
            query
                .first()
                .and_then(|&anchor| rs.coincident(anchor, &query)),
            expected,
        );
    }

    #[rstest]
    #[case::exact(vec![n(0), n(1), n(2)], Some(RelationId(0)))]
    #[case::reordered(vec![n(2), n(0), n(1)], Some(RelationId(0)))]
    #[case::second(vec![n(3), n(4)], Some(RelationId(1)))]
    #[case::subset(vec![n(0), n(1)], None)]
    #[case::superset(vec![n(0), n(1), n(2), n(3)], None)]
    fn test_var_relation_set_coincident(
        #[case] query: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        let rs: VarRelationSet<NodeId, ()> =
            VarRelationSet::new(vec![(vec![n(0), n(1), n(2)], ()), (vec![n(3), n(4)], ())]);
        assert_eq!(
            query
                .first()
                .and_then(|&anchor| rs.coincident(anchor, &query)),
            expected,
        );
    }

    #[rstest]
    #[case::exact(vec![n(0), n(1)], vec![n(2)], Some(RelationId(0)))]
    #[case::reordered_factor(vec![n(1), n(0)], vec![n(2)], Some(RelationId(0)))]
    #[case::second(vec![n(3), n(4)], vec![n(5)], Some(RelationId(1)))]
    #[case::absent(vec![n(0), n(1)], vec![n(9)], None)]
    fn test_fixed_fixed_birelation_set_coincident(
        #[case] query_1: Vec<NodeId>,
        #[case] query_2: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        let rs: FixedFixedBirelationSet<NodeId, 2, NodeId, 1, ()> =
            FixedFixedBirelationSet::new(vec![
                ([n(0), n(1)], [n(2)], ()),
                ([n(3), n(4)], [n(5)], ()),
            ]);
        assert_eq!(
            query_1
                .first()
                .and_then(|&anchor| rs.coincident(anchor, &query_1, &query_2)),
            expected,
        );
    }

    #[rstest]
    #[case::exact(vec![n(0)], vec![n(1)], Some(RelationId(0)))]
    #[case::role_swap(vec![n(1)], vec![n(0)], None)]
    #[case::multiset_reordered(vec![n(3)], vec![n(5), n(4), n(4)], Some(RelationId(1)))]
    #[case::wrong_multiplicity(vec![n(3)], vec![n(4), n(5)], None)]
    #[case::absent(vec![n(0)], vec![n(2)], None)]
    fn test_fixed_var_birelation_set_coincident(
        #[case] query_1: Vec<NodeId>,
        #[case] query_2: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        // Factor2 `Ordered` (a coset frame) yet matched as a multiset: duplicate `n(4)` and the
        // role-swap case exercise the key semantics.
        let rs: FixedVarBirelationSet<NodeId, 1, NodeId, ()> = FixedVarBirelationSet::new(vec![
            ([n(0)], vec![n(1)], ()),
            ([n(3)], vec![n(4), n(4), n(5)], ()),
        ]);
        assert_eq!(
            query_1
                .first()
                .and_then(|&anchor| rs.coincident(anchor, &query_1, &query_2)),
            expected,
        );
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
        let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, ()> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1), n(2)], ())]);
        assert_eq!(
            rs.coincident_edge(EdgeId(0), &[EdgeId(0)], &ligands),
            expected,
        );
    }

    #[rstest]
    #[case::exact(vec![n(0), n(1)], vec![n(2), n(3)], Some(RelationId(0)))]
    #[case::reordered(vec![n(1), n(0)], vec![n(3), n(2)], Some(RelationId(0)))]
    #[case::role_swap(vec![n(2), n(3)], vec![n(0), n(1)], None)]
    #[case::absent(vec![n(0), n(1)], vec![n(2), n(9)], None)]
    fn test_var_var_birelation_set_coincident(
        #[case] query_1: Vec<NodeId>,
        #[case] query_2: Vec<NodeId>,
        #[case] expected: Option<RelationId>,
    ) {
        let rs: VarVarBirelationSet<NodeId, NodeId, ()> = VarVarBirelationSet::new(vec![
            (vec![n(0), n(1)], vec![n(2), n(3)], ()),
            (vec![n(4)], vec![n(5)], ()),
        ]);
        assert_eq!(
            query_1
                .first()
                .and_then(|&anchor| rs.coincident(anchor, &query_1, &query_2)),
            expected,
        );
    }

    #[rstest]
    fn test_fixed_relation_set_pushout() {
        // same-space glue: self {01}=10 {23}=20 ; right {01}=5 (coincides) {45}=30 (new); combine=sum.
        let left =
            FixedRelationSet::<NodeId, i32, 2>::new(vec![([n(0), n(1)], 10), ([n(2), n(3)], 20)]);
        let right =
            FixedRelationSet::<NodeId, i32, 2>::new(vec![([n(0), n(1)], 5), ([n(4), n(5)], 30)]);
        let glue = left
            .pushout(
                &right,
                |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
                |(_, a), (_, b)| Some(a + b),
            )
            .expect("no ⊥");
        assert_eq!(glue.object.count(), 3);
        let value = |q: [NodeId; 2]| {
            glue.object
                .coincident(q[0], &q)
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
        let left = FixedRelationSet::<NodeId, i32, 2>::new(vec![([n(0), n(1)], 10)]);
        let right = FixedRelationSet::<NodeId, i32, 2>::new(vec![([n(0), n(1)], 5)]);
        assert_eq!(
            left.pushout(
                &right,
                |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
                |_, _| None
            ),
            None
        );
    }

    #[rstest]
    fn test_var_relation_set_pushout() {
        let left = VarRelationSet::<NodeId, i32>::new(vec![(vec![n(0), n(1), n(2)], 10)]);
        let right = VarRelationSet::<NodeId, i32>::new(vec![
            (vec![n(0), n(1), n(2)], 5),
            (vec![n(3), n(4)], 20),
        ]);
        let glue = left
            .pushout(
                &right,
                |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
                |(_, a), (_, b)| Some(a + b),
            )
            .expect("no ⊥");
        assert_eq!(glue.object.count(), 2);
        assert_eq!(
            glue.object
                .coincident(n(0), &[n(0), n(1), n(2)])
                .map(|id| *glue.object.data(id)),
            Some(15),
        );
        assert_eq!(glue.right.right_of(RelationId(1)), Some(RelationId(1))); // {34} appended
    }

    #[rstest]
    fn test_fixed_var_birelation_set_pushout() {
        // coincidence requires *both* factors equal (site + members).
        let left = FixedVarBirelationSet::<NodeId, 1, NodeId, i32>::new(vec![(
            [n(0)],
            vec![n(1), n(2)],
            10,
        )]);
        let right = FixedVarBirelationSet::<NodeId, 1, NodeId, i32>::new(vec![
            ([n(0)], vec![n(1), n(2)], 5),
            ([n(3)], vec![n(4), n(5)], 20),
        ]);
        let glue = left
            .pushout(
                &right,
                |set: &_, q1: &[NodeId], q2: &_| {
                    q1.first().and_then(|&n| set.coincident(n, q1, q2))
                },
                |(_, _, a), (_, _, b)| Some(a + b),
            )
            .expect("no ⊥");
        assert_eq!(glue.object.count(), 2);
        assert_eq!(
            glue.object
                .coincident(n(0), &[n(0)], &[n(1), n(2)])
                .map(|id| *glue.object.data(id)),
            Some(15),
        );
        assert_eq!(glue.right.right_of(RelationId(1)), Some(RelationId(1))); // appended
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_pushout() {
        let left = FixedFixedBirelationSet::<NodeId, 1, NodeId, 2, i32>::new(vec![(
            [n(0)],
            [n(1), n(2)],
            10,
        )]);
        let right = FixedFixedBirelationSet::<NodeId, 1, NodeId, 2, i32>::new(vec![
            ([n(0)], [n(1), n(2)], 5),
            ([n(3)], [n(4), n(5)], 20),
        ]);
        let glue = left
            .pushout(
                &right,
                |set: &_, q1: &[NodeId], q2: &_| {
                    q1.first().and_then(|&n| set.coincident(n, q1, q2))
                },
                |(_, _, a), (_, _, b)| Some(a + b),
            )
            .expect("no ⊥");
        assert_eq!(glue.object.count(), 2);
        assert_eq!(
            glue.object
                .coincident(n(0), &[n(0)], &[n(1), n(2)])
                .map(|id| *glue.object.data(id)),
            Some(15),
        );
        assert_eq!(glue.right.right_of(RelationId(1)), Some(RelationId(1)));
    }

    #[rstest]
    fn test_var_var_birelation_set_pushout() {
        let left = VarVarBirelationSet::<NodeId, NodeId, i32>::new(vec![(
            vec![n(0), n(1)],
            vec![n(2), n(3)],
            10,
        )]);
        let right = VarVarBirelationSet::<NodeId, NodeId, i32>::new(vec![
            (vec![n(0), n(1)], vec![n(2), n(3)], 5),
            (vec![n(4)], vec![n(5)], 20),
        ]);
        let glue = left
            .pushout(
                &right,
                |set: &_, q1: &[NodeId], q2: &_| {
                    q1.first().and_then(|&n| set.coincident(n, q1, q2))
                },
                |(_, _, a), (_, _, b)| Some(a + b),
            )
            .expect("no ⊥");
        assert_eq!(glue.object.count(), 2);
        assert_eq!(
            glue.object
                .coincident(n(0), &[n(0), n(1)], &[n(2), n(3)])
                .map(|id| *glue.object.data(id)),
            Some(15),
        );
        assert_eq!(glue.right.right_of(RelationId(1)), Some(RelationId(1)));
    }

    #[rstest]
    fn test_fixed_relation_set_pullback() {
        // intersection: self {01}=10 {23}=20 ; right {01}=5 {45}=30 — only {01} is shared.
        let left =
            FixedRelationSet::<NodeId, i32, 2>::new(vec![([n(0), n(1)], 10), ([n(2), n(3)], 20)]);
        let right =
            FixedRelationSet::<NodeId, i32, 2>::new(vec![([n(0), n(1)], 5), ([n(4), n(5)], 30)]);
        let pb = left
            .pullback(
                &right,
                |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
                |(_, a), (_, b)| Some(a + b),
            )
            .expect("no ⊥");
        assert_eq!(pb.object.count(), 1); // self-only and right-only dropped
        assert_eq!(*pb.object.data(RelationId(0)), 15);
        assert_eq!(pb.left.right_of(RelationId(0)), Some(RelationId(0))); // → self {01}
        assert_eq!(pb.right.right_of(RelationId(0)), Some(RelationId(0))); // → right {01}
    }

    #[rstest]
    fn test_fixed_relation_set_pullback_bottom() {
        let left = FixedRelationSet::<NodeId, i32, 2>::new(vec![([n(0), n(1)], 10)]);
        let right = FixedRelationSet::<NodeId, i32, 2>::new(vec![([n(0), n(1)], 5)]);
        assert_eq!(
            left.pullback(
                &right,
                |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
                |_, _| None
            ),
            None
        );
    }

    #[rstest]
    fn test_fixed_var_birelation_set_pullback() {
        let left = FixedVarBirelationSet::<NodeId, 1, NodeId, i32>::new(vec![
            ([n(0)], vec![n(1), n(2)], 10),
            ([n(3)], vec![n(4)], 20),
        ]);
        let right = FixedVarBirelationSet::<NodeId, 1, NodeId, i32>::new(vec![(
            [n(0)],
            vec![n(1), n(2)],
            5,
        )]);
        let pb = left
            .pullback(
                &right,
                |set: &_, q1: &[NodeId], q2: &_| {
                    q1.first().and_then(|&n| set.coincident(n, q1, q2))
                },
                |(_, _, a), (_, _, b)| Some(a + b),
            )
            .expect("no ⊥");
        assert_eq!(pb.object.count(), 1); // only the shared ([0],[1,2])
        assert_eq!(*pb.object.data(RelationId(0)), 15);
    }
}
