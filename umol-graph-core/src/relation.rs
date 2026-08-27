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

use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

use crate::compaction::GraphCompaction;
use crate::correspondence::Correspondence;
use crate::graph::{EdgeId, NodeId};
use crate::remapping::Remapping;

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

/// A relation payload's coupling to participant order — the payload-side mirror of
/// [`RelationParticipant`] (which couples a participant to the id space). A relation set
/// canonicalizes its participants on construction and remap, so any position-indexed payload (e.g.
/// per-member electron counts) must follow that reorder via `on_permutation`. This is purely the
/// *structural* (position) coupling; value equivalence is a separate concern (in the consuming
/// crate, `on_permutation` composes with a canonical value equality to give the full framed compare).
pub trait RelationData {
    /// Reindex a position-indexed payload by `order` (the σ from `canonicalize_positions` /
    /// [`participant_permutation`](FixedRelationSet::participant_permutation)). A payload with no
    /// positional content is a no-op — but the impl is required, so the decision is explicit at every
    /// payload type.
    fn on_permutation(&mut self, order: &[ParticipantPosition]);

    /// `true` when `on_permutation` cannot change `self` (no positional content, or it is already
    /// wildcard) — a guard that lets consumers skip the reindex work. Conservative default `false`;
    /// override to `true` (or a value-dependent test) where the payload is provably invariant.
    fn is_permutation_invariant(&self) -> bool {
        false
    }
}

/// Two-factor analog of [`RelationData`] for birelation sets: `on_permutation` takes a permutation
/// per factor.
pub trait BiRelationData {
    fn on_permutation(&mut self, order_1: &[ParticipantPosition], order_2: &[ParticipantPosition]);

    fn is_permutation_invariant(&self) -> bool {
        false
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

/// A value that can occupy a relation factor: routes through a `GraphCompaction`
/// (removal/compaction, both directions) and a `Remapping` (general relabel,
/// forward), and exposes its node/edge refs for incidence. One impl per concrete
/// id type — dispatch is static, since a factor is homogeneous.
pub trait RelationParticipant: Copy + Ord + Hash {
    fn compact(self, compaction: &GraphCompaction) -> Option<Self>;
    fn uncompact(self, compaction: &GraphCompaction) -> Self;

    /// Relabel this participant through `remapping`.
    ///
    /// Every node or edge id read from `remapping` must be reported by [`refs`](Self::refs), so
    /// checked relation-set remapping can establish coverage before calling this method.
    fn remap(self, remapping: &Remapping) -> Self;

    /// Return every graph id used to represent this participant.
    fn refs(self) -> ParticipantRefs;

    /// The node or edge to route this participant through the incidence index, if any — narrows
    /// `find_by_participants` candidates (`None` falls back to a linear scan).
    fn anchor(self) -> Option<ParticipantAnchor>;
}

impl RelationParticipant for NodeId {
    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        compaction.compact_node(self)
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
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
    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        compaction.compact_edge(self)
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
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

impl<P: Hash, O, D: Hash, const N: usize> Hash for FixedRelationSet<P, O, D, N> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.participants.hash(state);
        self.data.hash(state);
    }
}

/// Relabel a factor's participants through a general `Remapping`, preserving their stored order.
/// The owning relation-set constructor canonicalizes the result and transports positional data.
fn remap_factor<P>(participants: &[P], remapping: &Remapping) -> Vec<P>
where
    P: RelationParticipant,
{
    participants.iter().map(|&p| p.remap(remapping)).collect()
}

fn participants_are_covered_by<P>(participants: &[P], remapping: &Remapping) -> bool
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
    pub fn new(entries: Vec<([P; N], D)>) -> Self
    where
        D: RelationData + Clone,
    {
        let mut participants = Vec::with_capacity(entries.len());
        let mut data = Vec::with_capacity(entries.len());
        for (mut p, mut d) in entries {
            let sigma = O::canonicalize_positions(&mut p);
            d.on_permutation(&sigma);
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

    pub fn data_iter_mut(&mut self) -> impl ExactSizeIterator<Item = &mut D> {
        self.data.iter_mut()
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

    /// The permutation σ reindexing `query` into relation `id`'s stored participant frame, or `None`
    /// if their participants differ (up to this factor's ordering). The σ-keeping, known-id sibling of
    /// [`find_by_participants`](Self::find_by_participants): the structural half of a relation compare,
    /// which the caller completes by reindexing the payload with σ and comparing values.
    pub fn participant_permutation(
        &self,
        id: RelationId,
        query: &[P],
    ) -> Option<Vec<ParticipantPosition>> {
        let mut canonical = query.to_vec();
        let sigma = O::canonicalize_positions(&mut canonical);
        (self.participants(id).as_slice() == canonical.as_slice()).then_some(sigma)
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

    pub fn relation_ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    /// Compact participant ids, dropping every relation that contains a removed participant.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: RelationData + Clone,
    {
        let entries: Vec<([P; N], D)> = (0..self.count())
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

    /// Relabel every participant and transport positional data into canonical participant order.
    ///
    /// # Semantic properties
    ///
    /// Each positional payload item remains attached to the participant whose id is relabeled.
    ///
    /// # Panics
    ///
    /// Panics when a participant lies outside the remapping's corresponding source range.
    pub fn remap(&self, remapping: &Remapping) -> Self
    where
        D: RelationData + Clone,
    {
        let entries: Vec<([P; N], D)> = (0..self.count())
            .map(|i| {
                let rid = RelationId(i as u32);
                let parts: [P; N] = remap_factor(self.participants(rid), remapping)
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                (parts, self.data(rid).clone())
            })
            .collect();
        Self::new(entries)
    }

    /// Relabel every participant, returning `None` when the remapping does not cover the set.
    pub fn try_remap(&self, remapping: &Remapping) -> Option<Self>
    where
        D: RelationData + Clone,
    {
        self.relation_ids()
            .all(|id| participants_are_covered_by(self.participants(id), remapping))
            .then(|| self.remap(remapping))
    }

    /// Glue `self` and `right`, both **already in the same participant id-space**, identifying
    /// coinciding relations (equal participants) — the same-space relation pushout. `combine` merges
    /// the data of a coincidence (`None` = ⊥ ⇒ the whole glue is inadmissible ⇒ `None`); every other
    /// relation is carried. `self`'s ids are the identity prefix of the object, `right`'s
    /// non-coinciding relations are appended. The caller brings both sides, including positional
    /// data, into the common space with [`remap`](Self::remap) first.
    pub fn pushout(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: RelationData + Clone,
    {
        let mut entries: Vec<([P; N], D)> = self
            .relation_ids()
            .map(|id| (*self.participants(id), self.data(id).clone()))
            .collect();
        let self_count = entries.len();
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.count());
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

    /// Same-space relation pullback — the shared relations (coinciding participants), data
    /// combined by `combine` (`None` = ⊥ ⇒ inadmissible); non-coinciding relations are dropped.
    /// Its two projections map each shared relation to its `self` / `right` original. Same-space
    /// contract as [`FixedRelationSet::pushout`].
    pub fn pullback(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPullback<Self>>
    where
        D: RelationData + Clone,
    {
        let mut entries: Vec<([P; N], D)> = Vec::new();
        let mut left_images: Vec<RelationId> = Vec::new();
        let mut right_images: Vec<RelationId> = Vec::new();
        for id in self.relation_ids() {
            if let Some(hit) = right.find_by_participants(self.participants(id)) {
                let merged = combine(self.data(id), right.data(hit))?;
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

impl<P: Hash, O, D: Hash> Hash for VarRelationSet<P, O, D> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.offsets.hash(state);
        self.participants.hash(state);
        self.data.hash(state);
    }
}

impl<P: RelationParticipant, O: FactorOrdering, D> VarRelationSet<P, O, D> {
    pub fn new(entries: Vec<(Vec<P>, D)>) -> Self
    where
        D: RelationData + Clone,
    {
        let relation_count = entries.len();
        let mut offsets = Vec::with_capacity(relation_count + 1);
        offsets.push(0);

        let total_participants: usize = entries.iter().map(|(p, _)| p.len()).sum();
        let mut participants = Vec::with_capacity(total_participants);
        let mut data = Vec::with_capacity(relation_count);

        for (mut p, mut d) in entries {
            let sigma = O::canonicalize_positions(&mut p);
            d.on_permutation(&sigma);
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

    pub fn data_iter_mut(&mut self) -> impl ExactSizeIterator<Item = &mut D> {
        self.data.iter_mut()
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

    /// The permutation σ reindexing `query` into relation `id`'s stored participant frame, or `None`
    /// if their participants differ (up to this factor's ordering). The σ-keeping, known-id sibling of
    /// [`find_by_participants`](Self::find_by_participants).
    pub fn participant_permutation(
        &self,
        id: RelationId,
        query: &[P],
    ) -> Option<Vec<ParticipantPosition>> {
        let mut canonical = query.to_vec();
        let sigma = O::canonicalize_positions(&mut canonical);
        (self.participants(id) == canonical.as_slice()).then_some(sigma)
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

    pub fn relation_ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    /// Compact participant ids, dropping every relation that contains a removed participant.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: RelationData + Clone,
    {
        let entries: Vec<(Vec<P>, D)> = (0..self.count())
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

    /// Relabel every participant and transport positional data into canonical participant order.
    ///
    /// # Semantic properties
    ///
    /// Each positional payload item remains attached to the participant whose id is relabeled.
    ///
    /// # Panics
    ///
    /// Panics when a participant lies outside the remapping's corresponding source range.
    pub fn remap(&self, remapping: &Remapping) -> Self
    where
        D: RelationData + Clone,
    {
        let entries: Vec<(Vec<P>, D)> = (0..self.count())
            .map(|i| {
                let rid = RelationId(i as u32);
                (
                    remap_factor(self.participants(rid), remapping),
                    self.data(rid).clone(),
                )
            })
            .collect();
        Self::new(entries)
    }

    /// Relabel every participant, returning `None` when the remapping does not cover the set.
    pub fn try_remap(&self, remapping: &Remapping) -> Option<Self>
    where
        D: RelationData + Clone,
    {
        self.relation_ids()
            .all(|id| participants_are_covered_by(self.participants(id), remapping))
            .then(|| self.remap(remapping))
    }

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`].
    pub fn pushout(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: RelationData + Clone,
    {
        let mut entries: Vec<(Vec<P>, D)> = self
            .relation_ids()
            .map(|id| (self.participants(id).to_vec(), self.data(id).clone()))
            .collect();
        let self_count = entries.len();
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.count());
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

    /// Same-space relation pullback — see [`FixedRelationSet::pullback`].
    pub fn pullback(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPullback<Self>>
    where
        D: RelationData + Clone,
    {
        let mut entries: Vec<(Vec<P>, D)> = Vec::new();
        let mut left_images: Vec<RelationId> = Vec::new();
        let mut right_images: Vec<RelationId> = Vec::new();
        for id in self.relation_ids() {
            if let Some(hit) = right.find_by_participants(self.participants(id)) {
                let merged = combine(self.data(id), right.data(hit))?;
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
    pub fn new(entries: Vec<([L1; N1], [L2; N2], D)>) -> Self
    where
        D: BiRelationData + Clone,
    {
        let relation_count = entries.len();
        let mut participants_1 = Vec::with_capacity(relation_count);
        let mut participants_2 = Vec::with_capacity(relation_count);
        let mut data = Vec::with_capacity(relation_count);
        for (mut l1, mut l2, mut d) in entries {
            let s1 = O1::canonicalize_positions(&mut l1);
            let s2 = O2::canonicalize_positions(&mut l2);
            d.on_permutation(&s1, &s2);
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

    pub fn data_iter_mut(&mut self) -> impl ExactSizeIterator<Item = &mut D> {
        self.data.iter_mut()
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

    /// The per-factor permutations (σ₁, σ₂) reindexing `query_1` / `query_2` into relation `id`'s
    /// stored participant frame, or `None` if either factor's participants differ (up to its
    /// ordering). The σ-keeping, known-id sibling of [`find_by_participants`](Self::find_by_participants).
    #[allow(clippy::type_complexity)]
    pub fn participant_permutation(
        &self,
        id: RelationId,
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<(Vec<ParticipantPosition>, Vec<ParticipantPosition>)> {
        let mut canonical_1 = query_1.to_vec();
        let s1 = O1::canonicalize_positions(&mut canonical_1);
        let mut canonical_2 = query_2.to_vec();
        let s2 = O2::canonicalize_positions(&mut canonical_2);
        (self.participants_1(id).as_slice() == canonical_1.as_slice()
            && self.participants_2(id).as_slice() == canonical_2.as_slice())
        .then_some((s1, s2))
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

    pub fn relation_ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    /// Compact participant ids, dropping every relation that contains a removed participant.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: BiRelationData + Clone,
    {
        let entries: Vec<([L1; N1], [L2; N2], D)> = (0..self.count())
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

    /// Relabel every participant and transport positional data into canonical participant order.
    ///
    /// # Semantic properties
    ///
    /// In both factors, each positional payload item remains attached to the participant whose id
    /// is relabeled.
    ///
    /// # Panics
    ///
    /// Panics when a participant lies outside the remapping's corresponding source range.
    pub fn remap(&self, remapping: &Remapping) -> Self
    where
        D: BiRelationData + Clone,
    {
        let entries: Vec<([L1; N1], [L2; N2], D)> = (0..self.count())
            .map(|i| {
                let rid = RelationId(i as u32);
                let f1: [L1; N1] = remap_factor(self.participants_1(rid), remapping)
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                let f2: [L2; N2] = remap_factor(self.participants_2(rid), remapping)
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                (f1, f2, self.data(rid).clone())
            })
            .collect();
        Self::new(entries)
    }

    /// Relabel every participant, returning `None` when the remapping does not cover either factor.
    pub fn try_remap(&self, remapping: &Remapping) -> Option<Self>
    where
        D: BiRelationData + Clone,
    {
        self.relation_ids()
            .all(|id| {
                participants_are_covered_by(self.participants_1(id), remapping)
                    && participants_are_covered_by(self.participants_2(id), remapping)
            })
            .then(|| self.remap(remapping))
    }

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`]. Coincidence is equality of
    /// both factors' participants.
    pub fn pushout(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: BiRelationData + Clone,
    {
        let mut entries: Vec<([L1; N1], [L2; N2], D)> = self
            .relation_ids()
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
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPullback<Self>>
    where
        D: BiRelationData + Clone,
    {
        let mut entries: Vec<([L1; N1], [L2; N2], D)> = Vec::new();
        let mut left_images: Vec<RelationId> = Vec::new();
        let mut right_images: Vec<RelationId> = Vec::new();
        for id in self.relation_ids() {
            if let Some(hit) =
                right.find_by_participants(self.participants_1(id), self.participants_2(id))
            {
                let merged = combine(self.data(id), right.data(hit))?;
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

impl<L1, O1, const N1: usize, L2, O2, D> Hash for FixedVarBirelationSet<L1, O1, N1, L2, O2, D>
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

impl<L1, O1, const N1: usize, L2, O2, D> FixedVarBirelationSet<L1, O1, N1, L2, O2, D>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
{
    pub fn new(entries: Vec<([L1; N1], Vec<L2>, D)>) -> Self
    where
        D: BiRelationData + Clone,
    {
        let relation_count = entries.len();
        let mut participants_1 = Vec::with_capacity(relation_count);
        let mut f2_offsets = Vec::with_capacity(relation_count + 1);
        f2_offsets.push(0);
        let mut participants_2 = Vec::new();
        let mut data = Vec::with_capacity(relation_count);
        for (mut l1, mut l2, mut d) in entries {
            let s1 = O1::canonicalize_positions(&mut l1);
            let s2 = O2::canonicalize_positions(&mut l2);
            d.on_permutation(&s1, &s2);
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
            _ordering: PhantomData,
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

    pub fn data_iter_mut(&mut self) -> impl ExactSizeIterator<Item = &mut D> {
        self.data.iter_mut()
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

    /// The per-factor permutations (σ₁, σ₂) reindexing `query_1` / `query_2` into relation `id`'s
    /// stored participant frame, or `None` if either factor's participants differ (up to its
    /// ordering). The σ-keeping, known-id sibling of [`find_by_participants`](Self::find_by_participants).
    #[allow(clippy::type_complexity)]
    pub fn participant_permutation(
        &self,
        id: RelationId,
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<(Vec<ParticipantPosition>, Vec<ParticipantPosition>)> {
        let mut canonical_1 = query_1.to_vec();
        let s1 = O1::canonicalize_positions(&mut canonical_1);
        let mut canonical_2 = query_2.to_vec();
        let s2 = O2::canonicalize_positions(&mut canonical_2);
        (self.participants_1(id).as_slice() == canonical_1.as_slice()
            && self.participants_2(id) == canonical_2.as_slice())
        .then_some((s1, s2))
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

    pub fn relation_ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    /// Compact participant ids, dropping every relation that contains a removed participant.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: BiRelationData + Clone,
    {
        let entries: Vec<([L1; N1], Vec<L2>, D)> = (0..self.count())
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

    /// Relabel every participant and transport positional data into canonical participant order.
    ///
    /// # Semantic properties
    ///
    /// In both factors, each positional payload item remains attached to the participant whose id
    /// is relabeled.
    ///
    /// # Panics
    ///
    /// Panics when a participant lies outside the remapping's corresponding source range.
    pub fn remap(&self, remapping: &Remapping) -> Self
    where
        D: BiRelationData + Clone,
    {
        let entries: Vec<([L1; N1], Vec<L2>, D)> = (0..self.count())
            .map(|i| {
                let rid = RelationId(i as u32);
                let f1: [L1; N1] = remap_factor(self.participants_1(rid), remapping)
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("factor arity preserved"));
                (
                    f1,
                    remap_factor(self.participants_2(rid), remapping),
                    self.data(rid).clone(),
                )
            })
            .collect();
        Self::new(entries)
    }

    /// Relabel every participant, returning `None` when the remapping does not cover either factor.
    pub fn try_remap(&self, remapping: &Remapping) -> Option<Self>
    where
        D: BiRelationData + Clone,
    {
        self.relation_ids()
            .all(|id| {
                participants_are_covered_by(self.participants_1(id), remapping)
                    && participants_are_covered_by(self.participants_2(id), remapping)
            })
            .then(|| self.remap(remapping))
    }

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`]. Coincidence is equality of
    /// **both** factors' participants.
    pub fn pushout(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: BiRelationData + Clone,
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
        let mut right_map: Vec<RelationId> = Vec::with_capacity(right.count());
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

    /// Same-space relation pullback — see [`FixedRelationSet::pullback`].
    pub fn pullback(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPullback<Self>>
    where
        D: BiRelationData + Clone,
    {
        let mut entries: Vec<([L1; N1], Vec<L2>, D)> = Vec::new();
        let mut left_images: Vec<RelationId> = Vec::new();
        let mut right_images: Vec<RelationId> = Vec::new();
        for id in self.relation_ids() {
            if let Some(hit) =
                right.find_by_participants(self.participants_1(id), self.participants_2(id))
            {
                let merged = combine(self.data(id), right.data(hit))?;
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
    pub fn new(entries: Vec<(Vec<L1>, Vec<L2>, D)>) -> Self
    where
        D: BiRelationData + Clone,
    {
        let relation_count = entries.len();
        let mut f1_offsets = Vec::with_capacity(relation_count + 1);
        f1_offsets.push(0);
        let mut participants_1 = Vec::new();
        let mut f2_offsets = Vec::with_capacity(relation_count + 1);
        f2_offsets.push(0);
        let mut participants_2 = Vec::new();
        let mut data = Vec::with_capacity(relation_count);
        for (mut l1, mut l2, mut d) in entries {
            let s1 = O1::canonicalize_positions(&mut l1);
            let s2 = O2::canonicalize_positions(&mut l2);
            d.on_permutation(&s1, &s2);
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
            _ordering: PhantomData,
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

    pub fn data_iter_mut(&mut self) -> impl ExactSizeIterator<Item = &mut D> {
        self.data.iter_mut()
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

    /// The per-factor permutations (σ₁, σ₂) reindexing `query_1` / `query_2` into relation `id`'s
    /// stored participant frame, or `None` if either factor's participants differ (up to its
    /// ordering). The σ-keeping, known-id sibling of [`find_by_participants`](Self::find_by_participants).
    #[allow(clippy::type_complexity)]
    pub fn participant_permutation(
        &self,
        id: RelationId,
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<(Vec<ParticipantPosition>, Vec<ParticipantPosition>)> {
        let mut canonical_1 = query_1.to_vec();
        let s1 = O1::canonicalize_positions(&mut canonical_1);
        let mut canonical_2 = query_2.to_vec();
        let s2 = O2::canonicalize_positions(&mut canonical_2);
        (self.participants_1(id) == canonical_1.as_slice()
            && self.participants_2(id) == canonical_2.as_slice())
        .then_some((s1, s2))
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

    pub fn relation_ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    /// Compact participant ids, dropping every relation that contains a removed participant.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: BiRelationData + Clone,
    {
        let entries: Vec<(Vec<L1>, Vec<L2>, D)> = (0..self.count())
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

    /// Relabel every participant and transport positional data into canonical participant order.
    ///
    /// # Semantic properties
    ///
    /// In both factors, each positional payload item remains attached to the participant whose id
    /// is relabeled.
    ///
    /// # Panics
    ///
    /// Panics when a participant lies outside the remapping's corresponding source range.
    pub fn remap(&self, remapping: &Remapping) -> Self
    where
        D: BiRelationData + Clone,
    {
        let entries: Vec<(Vec<L1>, Vec<L2>, D)> = (0..self.count())
            .map(|i| {
                let rid = RelationId(i as u32);
                (
                    remap_factor(self.participants_1(rid), remapping),
                    remap_factor(self.participants_2(rid), remapping),
                    self.data(rid).clone(),
                )
            })
            .collect();
        Self::new(entries)
    }

    /// Relabel every participant, returning `None` when the remapping does not cover either factor.
    pub fn try_remap(&self, remapping: &Remapping) -> Option<Self>
    where
        D: BiRelationData + Clone,
    {
        self.relation_ids()
            .all(|id| {
                participants_are_covered_by(self.participants_1(id), remapping)
                    && participants_are_covered_by(self.participants_2(id), remapping)
            })
            .then(|| self.remap(remapping))
    }

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`]. Coincidence is equality of
    /// both factors' participants.
    pub fn pushout(
        &self,
        right: &Self,
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPushout<Self>>
    where
        D: BiRelationData + Clone,
    {
        let mut entries: Vec<(Vec<L1>, Vec<L2>, D)> = self
            .relation_ids()
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
        mut combine: impl FnMut(&D, &D) -> Option<D>,
    ) -> Option<RelationPullback<Self>>
    where
        D: BiRelationData + Clone,
    {
        let mut entries: Vec<(Vec<L1>, Vec<L2>, D)> = Vec::new();
        let mut left_images: Vec<RelationId> = Vec::new();
        let mut right_images: Vec<RelationId> = Vec::new();
        for id in self.relation_ids() {
            if let Some(hit) =
                right.find_by_participants(self.participants_1(id), self.participants_2(id))
            {
                let merged = combine(self.data(id), right.data(hit))?;
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
    use std::fmt::Debug;
    use std::hash::{DefaultHasher, Hash, Hasher};

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn hash<T: Hash>(value: &T) -> u64 {
        let mut hasher = DefaultHasher::new();
        value.hash(&mut hasher);
        hasher.finish()
    }

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct PositionLabels(Vec<u32>);

    impl RelationData for PositionLabels {
        fn on_permutation(&mut self, order: &[ParticipantPosition]) {
            let previous = self.0.clone();
            self.0 = order
                .iter()
                .map(|position| previous[position.index()])
                .collect();
        }
    }

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct BiPositionLabels {
        factor_1: Vec<u32>,
        factor_2: Vec<u32>,
    }

    impl BiRelationData for BiPositionLabels {
        fn on_permutation(
            &mut self,
            order_1: &[ParticipantPosition],
            order_2: &[ParticipantPosition],
        ) {
            let previous_1 = self.factor_1.clone();
            self.factor_1 = order_1
                .iter()
                .map(|position| previous_1[position.index()])
                .collect();
            let previous_2 = self.factor_2.clone();
            self.factor_2 = order_2
                .iter()
                .map(|position| previous_2[position.index()])
                .collect();
        }
    }

    impl RelationData for &str {
        fn on_permutation(&mut self, _: &[ParticipantPosition]) {}
    }
    impl BiRelationData for &str {
        fn on_permutation(&mut self, _: &[ParticipantPosition], _: &[ParticipantPosition]) {}
    }
    impl RelationData for i32 {
        fn on_permutation(&mut self, _: &[ParticipantPosition]) {}
    }
    impl BiRelationData for i32 {
        fn on_permutation(&mut self, _: &[ParticipantPosition], _: &[ParticipantPosition]) {}
    }
    impl RelationData for () {
        fn on_permutation(&mut self, _: &[ParticipantPosition]) {}
    }
    impl BiRelationData for () {
        fn on_permutation(&mut self, _: &[ParticipantPosition], _: &[ParticipantPosition]) {}
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

    fn assert_data_iter_mut<'a>(
        mut iterator: impl ExactSizeIterator<Item = &'a mut i32>,
        expected: &[i32],
    ) {
        assert_eq!(iterator.len(), expected.len());
        assert_eq!(iterator.size_hint(), (expected.len(), Some(expected.len())));
        for expected_item in expected {
            let previous = iterator.len();
            let item = iterator.next().expect("expected another data item");
            assert_eq!(*item, *expected_item);
            *item *= 10;
            let remaining = iterator.len();
            assert_eq!(remaining, previous - 1);
            assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
        }
        assert_eq!(iterator.next(), None);
        assert_eq!(iterator.len(), 0);
        assert_eq!(iterator.size_hint(), (0, Some(0)));
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
        let compaction = GraphCompaction::new(vec![NodeId(1)], vec![]);
        assert_eq!(id.compact(&compaction), expected);
    }

    #[rstest]
    #[case::before_gap(NodeId(0), NodeId(0))]
    #[case::after_gap(NodeId(1), NodeId(2))]
    fn test_node_id_unmap(#[case] id: NodeId, #[case] expected: NodeId) {
        let compaction = GraphCompaction::new(vec![NodeId(1)], vec![]);
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
        let compaction = GraphCompaction::new(vec![], vec![EdgeId(0)]);
        assert_eq!(id.compact(&compaction), expected);
    }

    #[rstest]
    #[case::before_gap(EdgeId(0), EdgeId(0))]
    #[case::after_gap(EdgeId(1), EdgeId(2))]
    fn test_edge_id_unmap(#[case] id: EdgeId, #[case] expected: EdgeId) {
        let compaction = GraphCompaction::new(vec![], vec![EdgeId(1)]);
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
        assert_eq!(rs.count(), 2);
        assert_eq!(rs.data(RelationId(0)), &"dative");
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(2)]);
    }

    #[rstest]
    fn test_fixed_relation_set_hash() {
        let entries = vec![([n(2), n(0)], "first"), ([n(3), n(1)], "second")];
        let left: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(entries.clone());
        let right: FixedRelationSet<NodeId, Unordered, &str, 2> = FixedRelationSet::new(entries);
        assert_eq!(left, right);
        assert_eq!(hash(&left), hash(&right));
    }

    #[rstest]
    #[case::canonical_entries(
        vec![([n(2), n(0)], "first"), ([n(3), n(1)], "second")],
        vec![([n(0), n(2)], "first"), ([n(1), n(3)], "second")],
    )]
    fn test_fixed_relation_set_into_entries(
        #[case] entries: Vec<([NodeId; 2], &str)>,
        #[case] expected: Vec<([NodeId; 2], &str)>,
    ) {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> = FixedRelationSet::new(entries);
        assert_eq!(rs.into_entries(), expected);
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
        let mut empty = FixedRelationSet::<NodeId, Unordered, i32, 2>::default();
        assert_data_iter_mut(empty.data_iter_mut(), &[]);

        let mut rs: FixedRelationSet<NodeId, Unordered, i32, 2> = FixedRelationSet::new(vec![
            ([n(0), n(1)], 1),
            ([n(1), n(2)], 2),
            ([n(2), n(3)], 3),
        ]);
        assert_data_iter_mut(rs.data_iter_mut(), &[1, 2, 3]);
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
    #[case::rotation(vec![ParticipantPosition(2), ParticipantPosition(0), ParticipantPosition(1)], [n(2), n(0), n(1)])]
    #[case::transposition(vec![ParticipantPosition(1), ParticipantPosition(0), ParticipantPosition(2)], [n(1), n(0), n(2)])]
    #[case::reversal(vec![ParticipantPosition(2), ParticipantPosition(1), ParticipantPosition(0)], [n(2), n(1), n(0)])]
    fn test_fixed_relation_set_permute_with(
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: [NodeId; 3],
    ) {
        let mut rs: FixedRelationSet<NodeId, Ordered, &str, 3> =
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
    #[case::ordered_factor(FixedRelationSet::<NodeId, Ordered, &str, 3>::new(vec![
        ([n(2), n(0), n(1)], "a"),
    ]))]
    #[case::unordered_factor(FixedRelationSet::<NodeId, Unordered, &str, 3>::new(vec![
        ([n(2), n(0), n(1)], "a"),
    ]))]
    fn test_fixed_relation_set_permute_with_identity<O: FactorOrdering + Clone + Debug>(
        #[case] input: FixedRelationSet<NodeId, O, &'static str, 3>,
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
        let mut rs: FixedRelationSet<NodeId, Ordered, &str, 3> =
            FixedRelationSet::new(vec![([n(0), n(1), n(2)], "a")]);
        rs.permute_with(RelationId(0), &order);
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
        assert_exact_size(
            FixedRelationSet::<NodeId, Unordered, (), 2>::default().relation_ids(),
            vec![],
        );
        let rs: FixedRelationSet<NodeId, Unordered, (), 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], ()), ([n(1), n(2)], ())]);
        assert_exact_size(rs.relation_ids(), vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_fixed_relation_set_compact() {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(0), n(2)], "keep"), ([n(1), n(3)], "drop")]);
        let compaction = GraphCompaction::new(vec![NodeId(1)], vec![]);
        let out = rs.compact(&compaction);
        assert_eq!(out.count(), 1);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_fixed_relation_set_remap() {
        let rs: FixedRelationSet<NodeId, Unordered, PositionLabels, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], PositionLabels(vec![10, 11]))]);
        let remapping = Remapping::new(vec![n(1), n(0)], vec![]);
        let out = rs.remap(&remapping);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.data(RelationId(0)), &PositionLabels(vec![11, 10]));
    }

    #[rstest]
    #[case::covered(vec![n(1), n(0)], true)]
    #[case::uncovered_node(vec![n(0)], false)]
    fn test_fixed_relation_set_try_remap(#[case] nodes: Vec<NodeId>, #[case] covered: bool) {
        let rs: FixedRelationSet<NodeId, Unordered, PositionLabels, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], PositionLabels(vec![10, 11]))]);
        let remapping = Remapping::new(nodes, vec![]);
        let expected = covered.then(|| rs.remap(&remapping));
        assert_eq!(rs.try_remap(&remapping), expected);
    }

    #[rstest]
    fn test_fixed_relation_set_default() {
        let rs = FixedRelationSet::<NodeId, Unordered, (), 2>::default();
        assert_eq!(rs.count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_var_relation_set_new() {
        let rs: VarRelationSet<NodeId, Unordered, &str> =
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
        let left: VarRelationSet<NodeId, Unordered, &str> = VarRelationSet::new(entries.clone());
        let right: VarRelationSet<NodeId, Unordered, &str> = VarRelationSet::new(entries);
        assert_eq!(left, right);
        assert_eq!(hash(&left), hash(&right));
    }

    #[rstest]
    #[case::canonical_entries(
        vec![
            (vec![n(2), n(0)], "first"),
            (vec![n(4), n(3), n(1)], "second"),
        ],
        vec![
            (vec![n(0), n(2)], "first"),
            (vec![n(1), n(3), n(4)], "second"),
        ],
    )]
    fn test_var_relation_set_into_entries(
        #[case] entries: Vec<(Vec<NodeId>, &str)>,
        #[case] expected: Vec<(Vec<NodeId>, &str)>,
    ) {
        let rs: VarRelationSet<NodeId, Unordered, &str> = VarRelationSet::new(entries);
        assert_eq!(rs.into_entries(), expected);
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
        let mut empty = VarRelationSet::<NodeId, Unordered, i32>::default();
        assert_data_iter_mut(empty.data_iter_mut(), &[]);

        let mut rs: VarRelationSet<NodeId, Unordered, i32> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], 1),
            (vec![n(2), n(3), n(4)], 2),
            (vec![n(5)], 3),
        ]);
        assert_data_iter_mut(rs.data_iter_mut(), &[1, 2, 3]);
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
    #[case::second_relation(RelationId(1), vec![ParticipantPosition(2), ParticipantPosition(0), ParticipantPosition(1)],
        vec![n(4), n(2), n(3)])]
    #[case::last_relation(RelationId(2), vec![ParticipantPosition(1), ParticipantPosition(0)], vec![n(6), n(5)])]
    fn test_var_relation_set_permute_with(
        #[case] id: RelationId,
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: Vec<NodeId>,
    ) {
        let mut rs: VarRelationSet<NodeId, Ordered, &str> = VarRelationSet::new(vec![
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
        let input: VarRelationSet<NodeId, Ordered, &str> =
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
        let mut rs: VarRelationSet<NodeId, Ordered, &str> =
            VarRelationSet::new(vec![(vec![n(0), n(1)], "a"), (vec![n(2), n(3), n(4)], "b")]);
        rs.permute_with(RelationId(1), &order);
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
        assert_exact_size(
            VarRelationSet::<NodeId, Unordered, ()>::default().relation_ids(),
            vec![],
        );
        let rs: VarRelationSet<NodeId, Unordered, ()> =
            VarRelationSet::new(vec![(vec![n(0), n(1)], ()), (vec![n(1), n(2)], ())]);
        assert_exact_size(rs.relation_ids(), vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_var_relation_set_compact() {
        let rs: VarRelationSet<NodeId, Unordered, &str> = VarRelationSet::new(vec![
            (vec![n(0), n(2), n(4)], "keep"),
            (vec![n(1), n(3)], "drop"),
        ]);
        let compaction = GraphCompaction::new(vec![NodeId(1)], vec![]);
        let out = rs.compact(&compaction);
        assert_eq!(out.count(), 1);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1), n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_var_relation_set_remap() {
        let rs: VarRelationSet<EdgeId, Ordered, PositionLabels> = VarRelationSet::new(vec![(
            vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            PositionLabels(vec![20, 21, 22]),
        )]);
        let remapping = Remapping::new(vec![], vec![EdgeId(2), EdgeId(0), EdgeId(1)]);
        let out = rs.remap(&remapping);
        assert_eq!(
            out.participants(RelationId(0)),
            &[EdgeId(2), EdgeId(0), EdgeId(1)]
        );
        assert_eq!(out.data(RelationId(0)), &PositionLabels(vec![20, 21, 22]));
    }

    #[rstest]
    #[case::covered(vec![EdgeId(2), EdgeId(0), EdgeId(1)], true)]
    #[case::uncovered_edge(vec![EdgeId(2), EdgeId(0)], false)]
    fn test_var_relation_set_try_remap(#[case] edges: Vec<EdgeId>, #[case] covered: bool) {
        let rs: VarRelationSet<EdgeId, Ordered, PositionLabels> = VarRelationSet::new(vec![(
            vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            PositionLabels(vec![20, 21, 22]),
        )]);
        let remapping = Remapping::new(vec![], edges);
        let expected = covered.then(|| rs.remap(&remapping));
        assert_eq!(rs.try_remap(&remapping), expected);
    }

    #[rstest]
    fn test_var_relation_set_default() {
        let rs = VarRelationSet::<NodeId, Unordered, ()>::default();
        assert_eq!(rs.count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_new() {
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 2, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(2), n(1)], "x")]);
        assert_eq!(rs.count(), 1);
        assert_eq!(rs.participants_1(RelationId(0)), &[n(0)]);
        assert_eq!(rs.participants_2(RelationId(0)), &[n(1), n(2)]);
        assert_eq!(rs.data(RelationId(0)), &"x");
    }

    #[rstest]
    #[case::canonical_entries(
        vec![
            ([n(2)], [n(4), n(1)], "first"),
            ([n(3)], [n(5), n(0)], "second"),
        ],
        vec![
            ([n(2)], [n(1), n(4)], "first"),
            ([n(3)], [n(0), n(5)], "second"),
        ],
    )]
    fn test_fixed_fixed_birelation_set_into_entries(
        #[case] entries: Vec<([NodeId; 1], [NodeId; 2], &str)>,
        #[case] expected: Vec<([NodeId; 1], [NodeId; 2], &str)>,
    ) {
        let rs: FixedFixedBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, 2, &str> =
            FixedFixedBirelationSet::new(entries);
        assert_eq!(rs.into_entries(), expected);
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
        let mut empty =
            FixedFixedBirelationSet::<NodeId, Unordered, 1, NodeId, Unordered, 1, i32>::default();
        assert_data_iter_mut(empty.data_iter_mut(), &[]);

        let mut rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 1, i32> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], 1), ([n(2)], [n(3)], 2)]);
        assert_data_iter_mut(rs.data_iter_mut(), &[1, 2]);
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_permute_1_with() {
        let mut rs: FixedFixedBirelationSet<NodeId, Ordered, 3, EdgeId, Ordered, 2, &str> =
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
        let mut rs: FixedFixedBirelationSet<NodeId, Ordered, 3, EdgeId, Ordered, 2, &str> =
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
        assert_exact_size(
            FixedFixedBirelationSet::<NodeId, Unordered, 1, NodeId, Unordered, 1, &str>::default()
                .relation_ids(),
            vec![],
        );
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 1, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], "a"), ([n(2)], [n(3)], "b")]);
        assert_exact_size(rs.relation_ids(), vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_compact() {
        // dropped relation loses a factor-1 participant
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 2, &str> =
            FixedFixedBirelationSet::new(vec![
                ([n(0)], [n(2), n(4)], "keep"),
                ([n(1)], [n(5), n(6)], "drop"),
            ]);
        let compaction = GraphCompaction::new(vec![NodeId(1)], vec![]);
        let out = rs.compact(&compaction);
        assert_eq!(out.count(), 1);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(1), n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_remap() {
        let rs: FixedFixedBirelationSet<
            NodeId,
            Unordered,
            2,
            EdgeId,
            Unordered,
            2,
            BiPositionLabels,
        > = FixedFixedBirelationSet::new(vec![(
            [n(0), n(1)],
            [EdgeId(0), EdgeId(1)],
            BiPositionLabels {
                factor_1: vec![10, 11],
                factor_2: vec![20, 21],
            },
        )]);
        let remapping = Remapping::new(vec![n(1), n(0)], vec![EdgeId(1), EdgeId(0)]);
        let out = rs.remap(&remapping);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.participants_2(RelationId(0)), &[EdgeId(0), EdgeId(1)]);
        assert_eq!(
            out.data(RelationId(0)),
            &BiPositionLabels {
                factor_1: vec![11, 10],
                factor_2: vec![21, 20],
            }
        );
    }

    #[rstest]
    #[case::covered(vec![n(1), n(0)], vec![EdgeId(1), EdgeId(0)], true)]
    #[case::uncovered_node(vec![n(0)], vec![EdgeId(1), EdgeId(0)], false)]
    #[case::uncovered_edge(vec![n(1), n(0)], vec![EdgeId(1)], false)]
    fn test_fixed_fixed_birelation_set_try_remap(
        #[case] nodes: Vec<NodeId>,
        #[case] edges: Vec<EdgeId>,
        #[case] covered: bool,
    ) {
        let rs: FixedFixedBirelationSet<
            NodeId,
            Unordered,
            2,
            EdgeId,
            Unordered,
            2,
            BiPositionLabels,
        > = FixedFixedBirelationSet::new(vec![(
            [n(0), n(1)],
            [EdgeId(0), EdgeId(1)],
            BiPositionLabels {
                factor_1: vec![10, 11],
                factor_2: vec![20, 21],
            },
        )]);
        let remapping = Remapping::new(nodes, edges);
        let expected = covered.then(|| rs.remap(&remapping));
        assert_eq!(rs.try_remap(&remapping), expected);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_default() {
        let rs =
            FixedFixedBirelationSet::<NodeId, Unordered, 1, NodeId, Unordered, 1, ()>::default();
        assert_eq!(rs.count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_fixed_var_birelation_set_new() {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, &str> =
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
        let left: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Unordered, &str> =
            FixedVarBirelationSet::new(entries.clone());
        let right: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Unordered, &str> =
            FixedVarBirelationSet::new(entries);
        assert_eq!(left, right);
        assert_eq!(hash(&left), hash(&right));
    }

    #[rstest]
    #[case::canonical_entries(
        vec![
            ([EdgeId(2)], vec![n(3), n(1)], "first"),
            ([EdgeId(4)], vec![n(5), n(0)], "second"),
        ],
        vec![
            ([EdgeId(2)], vec![n(1), n(3)], "first"),
            ([EdgeId(4)], vec![n(0), n(5)], "second"),
        ],
    )]
    fn test_fixed_var_birelation_set_into_entries(
        #[case] entries: Vec<([EdgeId; 1], Vec<NodeId>, &str)>,
        #[case] expected: Vec<([EdgeId; 1], Vec<NodeId>, &str)>,
    ) {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Unordered, &str> =
            FixedVarBirelationSet::new(entries);
        assert_eq!(rs.into_entries(), expected);
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
        let mut empty =
            FixedVarBirelationSet::<EdgeId, Ordered, 1, NodeId, Ordered, i32>::default();
        assert_data_iter_mut(empty.data_iter_mut(), &[]);

        let mut rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, i32> =
            FixedVarBirelationSet::new(vec![
                ([EdgeId(0)], vec![n(1)], 1),
                ([EdgeId(1)], vec![n(2)], 2),
            ]);
        assert_data_iter_mut(rs.data_iter_mut(), &[1, 2]);
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_permute_1_with() {
        let mut rs: FixedVarBirelationSet<EdgeId, Ordered, 2, NodeId, Ordered, &str> =
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
        let mut rs: FixedVarBirelationSet<EdgeId, Ordered, 2, NodeId, Ordered, &str> =
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
    fn test_fixed_var_birelation_set_permute_with_identity() {
        let input: FixedVarBirelationSet<EdgeId, Ordered, 2, NodeId, Ordered, &str> =
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
        let mut rs: FixedVarBirelationSet<EdgeId, Ordered, 2, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(4), EdgeId(5)], vec![n(0), n(1), n(2)], "a")]);
        rs.permute_2_with(RelationId(0), &order);
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
        assert_exact_size(
            FixedVarBirelationSet::<EdgeId, Ordered, 1, NodeId, Ordered, &str>::default()
                .relation_ids(),
            vec![],
        );
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![
                ([EdgeId(0)], vec![n(1)], "a"),
                ([EdgeId(1)], vec![n(2)], "b"),
            ]);
        assert_exact_size(rs.relation_ids(), vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_compact() {
        let rs: FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![
                ([n(0)], vec![n(2), n(4)], "keep"),
                ([n(5)], vec![n(1), n(3)], "drop"),
            ]);
        let compaction = GraphCompaction::new(vec![NodeId(1)], vec![]);
        let out = rs.compact(&compaction);
        assert_eq!(out.count(), 1);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(1), n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_fixed_var_birelation_set_remap() {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 2, NodeId, Unordered, BiPositionLabels> =
            FixedVarBirelationSet::new(vec![(
                [EdgeId(0), EdgeId(1)],
                vec![n(0), n(1), n(2)],
                BiPositionLabels {
                    factor_1: vec![30, 31],
                    factor_2: vec![40, 41, 42],
                },
            )]);
        let remapping = Remapping::new(vec![n(2), n(0), n(1)], vec![EdgeId(2), EdgeId(0)]);
        let out = rs.remap(&remapping);
        assert_eq!(out.participants_1(RelationId(0)), &[EdgeId(2), EdgeId(0)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(0), n(1), n(2)]);
        assert_eq!(
            out.data(RelationId(0)),
            &BiPositionLabels {
                factor_1: vec![30, 31],
                factor_2: vec![41, 42, 40],
            }
        );
    }

    #[rstest]
    #[case::covered(
        vec![n(2), n(0), n(1)],
        vec![EdgeId(2), EdgeId(0)],
        true,
    )]
    #[case::uncovered_node(vec![n(2), n(0)], vec![EdgeId(2), EdgeId(0)], false)]
    #[case::uncovered_edge(vec![n(2), n(0), n(1)], vec![EdgeId(2)], false)]
    fn test_fixed_var_birelation_set_try_remap(
        #[case] nodes: Vec<NodeId>,
        #[case] edges: Vec<EdgeId>,
        #[case] covered: bool,
    ) {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 2, NodeId, Unordered, BiPositionLabels> =
            FixedVarBirelationSet::new(vec![(
                [EdgeId(0), EdgeId(1)],
                vec![n(0), n(1), n(2)],
                BiPositionLabels {
                    factor_1: vec![30, 31],
                    factor_2: vec![40, 41, 42],
                },
            )]);
        let remapping = Remapping::new(nodes, edges);
        let expected = covered.then(|| rs.remap(&remapping));
        assert_eq!(rs.try_remap(&remapping), expected);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_default() {
        let rs = FixedVarBirelationSet::<EdgeId, Ordered, 1, NodeId, Ordered, ()>::default();
        assert_eq!(rs.count(), 0);
        assert!(!rs.has_incident_edge(EdgeId(0)));
    }

    #[rstest]
    fn test_var_var_birelation_set_new() {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![(vec![n(0), n(1)], vec![EdgeId(5)], "y")]);
        assert_eq!(rs.count(), 1);
        assert_eq!(rs.participants_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants_2(RelationId(0)), &[EdgeId(5)]);
        assert_eq!(rs.data(RelationId(0)), &"y");
    }

    #[rstest]
    #[case::canonical_entries(
        vec![
            (vec![n(2), n(0)], vec![EdgeId(4), EdgeId(1)], "first"),
            (vec![n(5), n(3)], vec![EdgeId(6), EdgeId(2)], "second"),
        ],
        vec![
            (vec![n(0), n(2)], vec![EdgeId(1), EdgeId(4)], "first"),
            (vec![n(3), n(5)], vec![EdgeId(2), EdgeId(6)], "second"),
        ],
    )]
    fn test_var_var_birelation_set_into_entries(
        #[case] entries: Vec<(Vec<NodeId>, Vec<EdgeId>, &str)>,
        #[case] expected: Vec<(Vec<NodeId>, Vec<EdgeId>, &str)>,
    ) {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(entries);
        assert_eq!(rs.into_entries(), expected);
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
        let mut empty = VarVarBirelationSet::<NodeId, Unordered, EdgeId, Unordered, i32>::default();
        assert_data_iter_mut(empty.data_iter_mut(), &[]);

        let mut rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, i32> =
            VarVarBirelationSet::new(vec![
                (vec![n(0)], vec![EdgeId(1)], 1),
                (vec![n(2)], vec![EdgeId(3)], 2),
            ]);
        assert_data_iter_mut(rs.data_iter_mut(), &[1, 2]);
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_var_var_birelation_set_permute_1_with() {
        let mut rs: VarVarBirelationSet<NodeId, Ordered, EdgeId, Ordered, &str> =
            VarVarBirelationSet::new(vec![
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
        let mut rs: VarVarBirelationSet<NodeId, Ordered, EdgeId, Ordered, &str> =
            VarVarBirelationSet::new(vec![
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
        assert_exact_size(
            VarVarBirelationSet::<NodeId, Unordered, EdgeId, Unordered, &str>::default()
                .relation_ids(),
            vec![],
        );
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![
                (vec![n(0)], vec![EdgeId(1)], "a"),
                (vec![n(2)], vec![EdgeId(3)], "b"),
            ]);
        assert_exact_size(rs.relation_ids(), vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_var_var_birelation_set_compact() {
        // dropped relation loses a factor-2 participant
        let rs: VarVarBirelationSet<NodeId, Unordered, NodeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![
                (vec![n(0), n(2)], vec![n(4)], "keep"),
                (vec![n(5)], vec![n(1)], "drop"),
            ]);
        let compaction = GraphCompaction::new(vec![NodeId(1)], vec![]);
        let out = rs.compact(&compaction);
        assert_eq!(out.count(), 1);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.participants_2(RelationId(0)), &[n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_var_var_birelation_set_remap() {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Ordered, BiPositionLabels> =
            VarVarBirelationSet::new(vec![(
                vec![n(0), n(1)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                BiPositionLabels {
                    factor_1: vec![50, 51],
                    factor_2: vec![60, 61, 62],
                },
            )]);
        let remapping = Remapping::new(vec![n(1), n(0)], vec![EdgeId(2), EdgeId(0), EdgeId(1)]);
        let out = rs.remap(&remapping);
        assert_eq!(out.participants_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(
            out.participants_2(RelationId(0)),
            &[EdgeId(2), EdgeId(0), EdgeId(1)]
        );
        assert_eq!(
            out.data(RelationId(0)),
            &BiPositionLabels {
                factor_1: vec![51, 50],
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
    #[case::uncovered_node(vec![n(1)], vec![EdgeId(2), EdgeId(0), EdgeId(1)], false)]
    #[case::uncovered_edge(vec![n(1), n(0)], vec![EdgeId(2), EdgeId(0)], false)]
    fn test_var_var_birelation_set_try_remap(
        #[case] nodes: Vec<NodeId>,
        #[case] edges: Vec<EdgeId>,
        #[case] covered: bool,
    ) {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Ordered, BiPositionLabels> =
            VarVarBirelationSet::new(vec![(
                vec![n(0), n(1)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2)],
                BiPositionLabels {
                    factor_1: vec![50, 51],
                    factor_2: vec![60, 61, 62],
                },
            )]);
        let remapping = Remapping::new(nodes, edges);
        let expected = covered.then(|| rs.remap(&remapping));
        assert_eq!(rs.try_remap(&remapping), expected);
    }

    #[rstest]
    fn test_var_var_birelation_set_default() {
        let rs = VarVarBirelationSet::<NodeId, Unordered, EdgeId, Unordered, ()>::default();
        assert_eq!(rs.count(), 0);
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
        assert_eq!(glue.object.count(), 3);
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
        assert_eq!(glue.object.count(), 2);
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
        assert_eq!(glue.object.count(), 2);
        assert_eq!(
            glue.object
                .find_by_participants(&[n(0)], &[n(1), n(2)])
                .map(|id| *glue.object.data(id)),
            Some(15),
        );
        assert_eq!(glue.right.right_of(RelationId(1)), Some(RelationId(1))); // appended
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_pushout() {
        let left = FixedFixedBirelationSet::<NodeId, Ordered, 1, NodeId, Unordered, 2, i32>::new(
            vec![([n(0)], [n(1), n(2)], 10)],
        );
        let right =
            FixedFixedBirelationSet::<NodeId, Ordered, 1, NodeId, Unordered, 2, i32>::new(vec![
                ([n(0)], [n(1), n(2)], 5),
                ([n(3)], [n(4), n(5)], 20),
            ]);
        let glue = left.pushout(&right, |a, b| Some(a + b)).expect("no ⊥");
        assert_eq!(glue.object.count(), 2);
        assert_eq!(
            glue.object
                .find_by_participants(&[n(0)], &[n(1), n(2)])
                .map(|id| *glue.object.data(id)),
            Some(15),
        );
        assert_eq!(glue.right.right_of(RelationId(1)), Some(RelationId(1)));
    }

    #[rstest]
    fn test_var_var_birelation_set_pushout() {
        let left = VarVarBirelationSet::<NodeId, Unordered, NodeId, Unordered, i32>::new(vec![(
            vec![n(0), n(1)],
            vec![n(2), n(3)],
            10,
        )]);
        let right = VarVarBirelationSet::<NodeId, Unordered, NodeId, Unordered, i32>::new(vec![
            (vec![n(0), n(1)], vec![n(2), n(3)], 5),
            (vec![n(4)], vec![n(5)], 20),
        ]);
        let glue = left.pushout(&right, |a, b| Some(a + b)).expect("no ⊥");
        assert_eq!(glue.object.count(), 2);
        assert_eq!(
            glue.object
                .find_by_participants(&[n(0), n(1)], &[n(2), n(3)])
                .map(|id| *glue.object.data(id)),
            Some(15),
        );
        assert_eq!(glue.right.right_of(RelationId(1)), Some(RelationId(1)));
    }

    #[rstest]
    fn test_fixed_relation_set_pullback() {
        // intersection: self {01}=10 {23}=20 ; right {01}=5 {45}=30 — only {01} is shared.
        let left = FixedRelationSet::<NodeId, Unordered, i32, 2>::new(vec![
            ([n(0), n(1)], 10),
            ([n(2), n(3)], 20),
        ]);
        let right = FixedRelationSet::<NodeId, Unordered, i32, 2>::new(vec![
            ([n(0), n(1)], 5),
            ([n(4), n(5)], 30),
        ]);
        let pb = left.pullback(&right, |a, b| Some(a + b)).expect("no ⊥");
        assert_eq!(pb.object.count(), 1); // self-only and right-only dropped
        assert_eq!(*pb.object.data(RelationId(0)), 15);
        assert_eq!(pb.left.right_of(RelationId(0)), Some(RelationId(0))); // → self {01}
        assert_eq!(pb.right.right_of(RelationId(0)), Some(RelationId(0))); // → right {01}
    }

    #[rstest]
    fn test_fixed_relation_set_pullback_bottom() {
        let left = FixedRelationSet::<NodeId, Unordered, i32, 2>::new(vec![([n(0), n(1)], 10)]);
        let right = FixedRelationSet::<NodeId, Unordered, i32, 2>::new(vec![([n(0), n(1)], 5)]);
        assert_eq!(left.pullback(&right, |_, _| None), None);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_pullback() {
        let left = FixedVarBirelationSet::<NodeId, Ordered, 1, NodeId, Unordered, i32>::new(vec![
            ([n(0)], vec![n(1), n(2)], 10),
            ([n(3)], vec![n(4)], 20),
        ]);
        let right =
            FixedVarBirelationSet::<NodeId, Ordered, 1, NodeId, Unordered, i32>::new(vec![(
                [n(0)],
                vec![n(1), n(2)],
                5,
            )]);
        let pb = left.pullback(&right, |a, b| Some(a + b)).expect("no ⊥");
        assert_eq!(pb.object.count(), 1); // only the shared ([0],[1,2])
        assert_eq!(*pb.object.data(RelationId(0)), 15);
    }
}
