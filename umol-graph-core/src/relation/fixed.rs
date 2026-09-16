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

    /// Replace every participant of relation `id`, preserving the supplied order and duplicates.
    ///
    /// As in [`Self::new`], participant references are indexed through
    /// [`RelationParticipant::refs`] without checking membership in an external graph.
    /// The fixed arity is enforced by the array type; zero arity is permitted.
    /// This rebuilds incidence over the entire collection.
    ///
    /// # Semantic properties
    ///
    /// Relation ids, row count, all payloads, and all other rows remain unchanged.
    /// Incidence contains each relation exactly once for every node or edge referenced by
    /// any of its final participants. The payload is not interpreted or transported.
    /// Replacing a row with its existing sequence leaves the set unchanged.
    /// These laws are exercised against a row model in `tests/property/relation.rs`.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn replace_participants(&mut self, id: RelationId, participants: [P; N]) {
        self.participants[id.index()] = participants;
        self.incidence = Incidence::build(self.count(), |i, out| {
            out.extend(
                self.participants[i]
                    .iter()
                    .map(|participant| participant.refs()),
            );
        });
    }

    /// Replace one participant at a factor-local position in relation `id`.
    ///
    /// Delegates to [`Self::replace_participants`] with the edited array, sharing its
    /// admission and preservation contract and its full incidence rebuild.
    ///
    /// # Semantic properties
    ///
    /// Equivalent to replacing the whole array after changing only `position`.
    /// Every other participant position remains unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set or `position` is outside `0..N`.
    /// Every position is invalid when `N` is zero.
    pub fn replace_participant(
        &mut self,
        id: RelationId,
        position: ParticipantPosition,
        participant: P,
    ) {
        let mut participants = self.participants[id.index()];
        participants[position.index()] = participant;
        self.replace_participants(id, participants);
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
    /// contract as [`FixedRelationSet::pushout`](crate::FixedRelationSet::pushout).
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
