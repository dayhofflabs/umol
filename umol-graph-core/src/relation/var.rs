use std::hash::{Hash, Hasher};

use super::incidence::Incidence;
use super::participant::{ParticipantPosition, RelationParticipant};
use super::{
    participants_match, permute_participants, relation_pullback, relation_pushout,
    remappable_under, RelationId, RelationPullbackCorrespondence, RelationPushoutCorrespondence,
};
use crate::compact::{Compaction, GraphCompaction};
use crate::correspondence::GraphCorrespondence;
use crate::graph::{EdgeId, NodeId};
use crate::remap::GraphRemapping;

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

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`](crate::FixedRelationSet::pushout).
    pub fn pushout(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[P]) -> Option<RelationId>,
        combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
    ) -> Option<Self>
    where
        D: Clone,
    {
        self.tracked_pushout(right, coincident, combine)
            .map(|(object, _)| object)
    }

    /// Glue relation sets and return both input-to-result mappings with the result.
    ///
    /// Has the same result and failure behavior as [`Self::pushout`].
    pub fn tracked_pushout(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[P]) -> Option<RelationId>,
        mut combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
    ) -> Option<(Self, RelationPushoutCorrespondence)>
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

    /// Same-space relation pullback — see [`FixedRelationSet::pullback`](crate::FixedRelationSet::pullback).
    pub fn pullback(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[P]) -> Option<RelationId>,
        combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
    ) -> Option<Self>
    where
        D: Clone,
    {
        self.tracked_pullback(right, coincident, combine)
            .map(|(object, _)| object)
    }

    /// Return the shared relation set and its two result-to-input projections.
    ///
    /// Has the same result and failure behavior as [`Self::pullback`].
    pub fn tracked_pullback(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[P]) -> Option<RelationId>,
        mut combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
    ) -> Option<(Self, RelationPullbackCorrespondence)>
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
