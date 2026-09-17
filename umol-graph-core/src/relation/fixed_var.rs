//! A fixed first factor and variable second factor: [FixedVarBirelationSet].

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

/// Relations with a fixed-arity first factor and a variable-arity second factor.
///
/// Rows retain input order, participant multiplicity, and payloads. Coinciding rows are
/// permitted. Participant references are not checked against an external graph.
/// The first factor stores arrays; the second uses a flat buffer with row offsets.
/// Either factor may be empty.
/// Relation ids use u32 indices. Variable-factor offsets also use u32.
///
/// # Semantic properties
///
/// - `new(entries).into_entries() == entries` for representable sizes.
/// - Equality compares stored sequences and payloads, including their order.
/// - Incidence lists each relation once per referenced node or edge, in relation-id order.
///   The union spans both factors.
/// - Coincidence compares participant multisets in each factor; it ignores stored order
///   and payloads, but observes multiplicity and complete participant values.
///
/// Construction/query, mutation, and transport laws are exercised through public APIs in
/// `tests/property/relation.rs`. Transport laws require conforming
/// [`RelationParticipant`] implementations.
#[derive(Clone, Debug)]
pub struct FixedVarBirelationSet<L1, const N1: usize, L2, D> {
    participants_1: Vec<[L1; N1]>,
    f2_offsets: Vec<u32>,
    participants_2: Vec<L2>,
    data: Vec<D>,
    incidence: Incidence,
}

impl<L1, const N1: usize, L2, D> FixedVarBirelationSet<L1, N1, L2, D>
where
    L1: RelationParticipant,
    L2: RelationParticipant,
{
    /// Store entries in order and build their union incidence index.
    ///
    /// Preserves the supplied sequence in each factor, including repeated values.
    /// Node/edge references come from [`RelationParticipant::refs`]; graph membership is
    /// not validated. Payloads are stored without interpretation.
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

    /// Consume the set into its stored entries, in relation-id order.
    ///
    /// Preserves participant order and multiplicity and returns the current payloads.
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

    /// Return the number of stored relations.
    pub fn count(&self) -> usize {
        self.data.len()
    }

    /// Whether `id` lies in this set's dense relation-id range.
    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    /// Iterate over `0..count()` in ascending relation-id order.
    pub fn ids(&self) -> impl ExactSizeIterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
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
    /// Only payloads are mutable; participant storage and its derived index stay intact.
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

    /// Borrow relation `id`'s stored first-factor sequence.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn participants_1(&self, id: RelationId) -> &[L1; N1] {
        &self.participants_1[id.index()]
    }

    /// Borrow relation `id`'s stored second-factor sequence.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn participants_2(&self, id: RelationId) -> &[L2] {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        &self.participants_2[start..end]
    }

    /// Borrow the payload of relation `id`.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    /// Mutably borrow the payload of relation `id`, leaving participants and incidence intact.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    /// Relations referencing `node`, once each in ascending relation-id order.
    ///
    /// Returns an empty slice when no participant references `node`. References in
    /// either factor contribute to the same union index.
    pub fn incident_to_node(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident_to_node(node)
    }

    /// Relations referencing `edge`, once each in ascending relation-id order.
    ///
    /// Returns an empty slice when no participant references `edge`. References in
    /// either factor contribute to the same union index.
    pub fn incident_to_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_to_edge(edge)
    }

    /// Whether any participant in either factor references `node`.
    pub fn has_incident_to_node(&self, node: NodeId) -> bool {
        self.incidence.has_incident_to_node(node)
    }

    /// Whether any participant in either factor references `edge`.
    pub fn has_incident_to_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_to_edge(edge)
    }

    /// Find the first relation incident with `node` whose factors match the queries.
    ///
    /// Compares complete participant multisets in each factor, ignoring stored order and
    /// payloads. Returns the smallest matching relation id, or `None` if the incidence
    /// list contains no match. Duplicate matching rows are permitted. An unreferenced
    /// anchor yields `None`, even if a matching row exists elsewhere.
    pub fn coincident_to_node(
        &self,
        node: NodeId,
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<RelationId> {
        self.coincident_among(self.incident_to_node(node), query_1, query_2)
    }

    /// Find the first relation incident with `edge` whose factors match the queries.
    ///
    /// Uses the same multiset comparison and first-match rule as [`Self::coincident_to_node`].
    /// Returns `None` when the edge incidence list contains no matching row.
    pub fn coincident_to_edge(
        &self,
        edge: EdgeId,
        query_1: &[L1],
        query_2: &[L2],
    ) -> Option<RelationId> {
        self.coincident_among(self.incident_to_edge(edge), query_1, query_2)
    }

    /// Whether relation `id` matches the query multisets.
    ///
    /// Compares complete participant values and multiplicities in each factor, without
    /// requiring a node or edge anchor. Stored order and payloads do not affect the result.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn is_coincident(&self, id: RelationId, query_1: &[L1], query_2: &[L2]) -> bool {
        self.coincident_among(&[id], query_1, query_2).is_some()
    }

    /// Return the first candidate whose stored participant multisets match both queries.
    fn coincident_among(
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

    /// Reorder relation `id`'s first-factor participants so that `new[i] = old[order[i]]`.
    ///
    /// # Semantic properties
    ///
    /// The selected multiset, the other factor, all payloads, and all other rows are unchanged.
    /// Incidence is unchanged and is not rebuilt. Identity leaves the set equal to itself;
    /// applying a permutation followed by its inverse recovers the original sequence.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set or `order` is not a permutation of the selected
    /// factor's positions (wrong length, repeated position, or out-of-range position).
    pub fn permute_participants_1(&mut self, id: RelationId, order: &[ParticipantPosition]) {
        permute_participants(self.participants_1[id.index()].as_mut_slice(), order);
    }

    /// Reorder relation `id`'s second-factor participants so that `new[i] = old[order[i]]`.
    ///
    /// # Semantic properties
    ///
    /// The selected multiset, the other factor, all payloads, and all other rows are unchanged.
    /// Incidence is unchanged and is not rebuilt. Identity leaves the set equal to itself;
    /// applying a permutation followed by its inverse recovers the original sequence.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set or `order` is not a permutation of the selected
    /// factor's positions (wrong length, repeated position, or out-of-range position).
    pub fn permute_participants_2(&mut self, id: RelationId, order: &[ParticipantPosition]) {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        permute_participants(&mut self.participants_2[start..end], order);
    }

    /// Replace both participant factors of relation `id`, preserving their supplied sequences.
    ///
    /// As in [`Self::new`], references are indexed through [`RelationParticipant::refs`]
    /// without checking membership in an external graph. The array enforces first-factor
    /// arity; either factor may be empty. Resizing the second factor moves later participants
    /// and adjusts their row offsets. Incidence is rebuilt over the entire collection once,
    /// after both factors have been replaced.
    ///
    /// # Semantic properties
    ///
    /// Relation ids, row count, all payloads, and all other rows remain unchanged.
    /// Participant order and duplicates are preserved. Incidence lists each relation once
    /// for every node or edge referenced by either final factor; a reference remains indexed
    /// while any participant in either factor retains it. Payloads are neither interpreted
    /// nor transported when lengths or positions change. Replacing both factors with their
    /// existing sequences preserves the set.
    /// These laws are exercised against a row model in `tests/property/relation.rs`.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn replace_participants(
        &mut self,
        id: RelationId,
        participants_1: [L1; N1],
        participants_2: &[L2],
    ) {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        let old_len = end - start;
        self.participants_1[id.index()] = participants_1;
        self.participants_2
            .splice(start..end, participants_2.iter().copied());
        for offset in &mut self.f2_offsets[id.index() + 1..] {
            *offset = (*offset as usize - old_len + participants_2.len()) as u32;
        }
        self.incidence = Incidence::build(self.count(), |i, out| {
            out.extend(self.participants_1[i].iter().map(|p| p.refs()));
            out.extend(
                self.participants_2(RelationId::from(i))
                    .iter()
                    .map(|p| p.refs()),
            );
        });
    }

    /// Replace the first participant factor of relation `id`, preserving the second factor.
    ///
    /// Delegates to [`Self::replace_participants`] with a copy of the current second factor,
    /// sharing its admission, preservation, and incidence contract. References retained by
    /// the second factor remain indexed.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn replace_participants_1(&mut self, id: RelationId, participants: [L1; N1]) {
        let second = self.participants_2(id).to_vec();
        self.replace_participants(id, participants, &second);
    }

    /// Replace the second participant factor of relation `id`, preserving the first factor.
    ///
    /// Delegates to [`Self::replace_participants`] with the current first factor, sharing
    /// its admission, preservation, and incidence contract. The replacement may be empty.
    /// References retained by the first factor remain indexed.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn replace_participants_2(&mut self, id: RelationId, participants: &[L2]) {
        self.replace_participants(id, self.participants_1[id.index()], participants);
    }

    /// Replace one participant at a position local to relation `id`'s first factor.
    ///
    /// Delegates to [`Self::replace_participants_1`] with the edited array.
    ///
    /// # Semantic properties
    ///
    /// Equivalent to whole-factor replacement after changing only `position`.
    /// Every other position and the complete second factor remain unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set or `position` is outside `0..N1`.
    /// Every position is invalid when `N1` is zero.
    pub fn replace_participant_1(
        &mut self,
        id: RelationId,
        position: ParticipantPosition,
        participant: L1,
    ) {
        let mut participants = self.participants_1[id.index()];
        participants[position.index()] = participant;
        self.replace_participants_1(id, participants);
    }

    /// Replace one participant at a position local to relation `id`'s second factor.
    ///
    /// Delegates to [`Self::replace_participants_2`] with the edited sequence.
    ///
    /// # Semantic properties
    ///
    /// Equivalent to whole-factor replacement after changing only `position`.
    /// The factor length, every other position, and the complete first factor remain unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set or `position` is outside the second factor.
    /// Every position is invalid for an empty second factor.
    pub fn replace_participant_2(
        &mut self,
        id: RelationId,
        position: ParticipantPosition,
        participant: L2,
    ) {
        let mut participants = self.participants_2(id).to_vec();
        participants[position.index()] = participant;
        self.replace_participants_2(id, &participants);
    }

    /// Insert a second-factor participant before `position`, or append at the factor length.
    ///
    /// Delegates to [`Self::replace_participants_2`] with the edited sequence, sharing its
    /// admission and preservation contract and its full incidence rebuild.
    ///
    /// # Semantic properties
    ///
    /// The second factor grows by one and existing participants retain their relative order.
    /// The first factor is unchanged. Removing the inserted position restores the original set.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set or `position` exceeds the second factor's length.
    pub fn insert_participant_2(
        &mut self,
        id: RelationId,
        position: ParticipantPosition,
        participant: L2,
    ) {
        let mut participants = self.participants_2(id).to_vec();
        participants.insert(position.index(), participant);
        self.replace_participants_2(id, &participants);
    }

    /// Remove the participant at a position local to relation `id`'s second factor.
    ///
    /// Delegates to [`Self::replace_participants_2`] with the edited sequence, sharing its
    /// admission and preservation contract and its full incidence rebuild.
    ///
    /// # Semantic properties
    ///
    /// The second factor shrinks by one and surviving participants retain their relative order.
    /// Removing its sole participant leaves an empty second factor with the same id and payload.
    /// The first factor is unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set or `position` is outside the second factor.
    /// Every position is invalid for an empty second factor.
    pub fn remove_participant_2(&mut self, id: RelationId, position: ParticipantPosition) {
        let mut participants = self.participants_2(id).to_vec();
        participants.remove(position.index());
        self.replace_participants_2(id, &participants);
    }

    /// Relabel participants through a correspondence, preserving rows, positions, and payloads.
    ///
    /// The asserted peer of [`Self::try_map`]. No payload or frame normalization is performed.
    ///
    /// # Panics
    ///
    /// Panics if any participant's [`RelationParticipant::try_map`] returns `None`,
    /// including when a referenced node or edge has no image.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self
    where
        D: Clone,
    {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Relabel participants, returning `None` if any participant cannot be mapped.
    ///
    /// Delegates admission to [`RelationParticipant::try_map`]; a missing image for any
    /// referenced node or edge rejects the whole operation. Unused ids need no image.
    ///
    /// # Semantic properties
    ///
    /// For conforming participants, row ids, positions, and payloads are preserved.
    /// An identity correspondence covering the references leaves the set equal to itself;
    /// sequential covered mappings agree with mapping through their composition.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self>
    where
        D: Clone,
    {
        self.map_participants(
            |participant| participant.try_map(correspondence),
            |participant| participant.try_map(correspondence),
        )
    }

    /// Relabel participants through a remapping, preserving rows, positions, and payloads.
    ///
    /// # Semantic properties
    ///
    /// For conforming participants and a covering remapping, identity leaves the set equal
    /// to itself and a remapping followed by its inverse recovers the original set.
    /// The payload stays in the supplied participant frame.
    ///
    /// # Panics
    ///
    /// Panics if a referenced node or edge is outside the remapping's source range.
    /// Participant remapping panics propagate.
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

    /// Relabel participants if the remapping covers every reported reference.
    ///
    /// Returns `None` if any node or edge reported by [`RelationParticipant::refs`] is
    /// outside the corresponding source range. Otherwise returns [`Self::remap`]'s result.
    /// Unused ids do not require coverage.
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

    /// Map every factor while retaining row order and payloads; reject the set on any `None`.
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

    /// Compact participant ids and discard rows containing a participant that cannot survive.
    ///
    /// A row is dropped if any [`RelationParticipant::compact`] call returns `None`.
    /// For node/edge participants this includes removed or out-of-source-range references.
    /// Returns the same set as [`Self::tracked_compact`], without its relation-id compaction.
    pub fn compact(&self, compaction: &GraphCompaction) -> Self
    where
        D: Clone,
    {
        self.tracked_compact(compaction).0
    }

    /// Compact participants and return the induced compaction of this set's relation ids.
    ///
    /// Drops a row if any participant's [`RelationParticipant::compact`] returns `None`.
    /// The compaction's source count is this set's original row count.
    ///
    /// # Semantic properties
    ///
    /// Surviving rows retain their relative order, factor positions, and cloned payloads;
    /// only participant ids and dense relation ids change. For conforming participants,
    /// a covering identity compaction preserves the set exactly. Plain and tracked
    /// compaction produce equal sets.
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

    /// Glue two relation sets in the same participant id space using caller-selected pairings.
    ///
    /// Calls `coincident` on the original `self` for each right row. The callback must
    /// identify a matching row, or return `None` for an unmatched row. Matches must be
    /// injective; storage does not verify participant equality or payload compatibility.
    /// `combine` receives both original factors and payloads, without frame alignment.
    /// Returns `None` as soon as `combine` rejects a pair; otherwise returns the glued set.
    ///
    /// # Semantic properties
    ///
    /// For valid pairings, `self` keeps its row ids and factors. Matched payloads are
    /// replaced by `combine`'s result. Unmatched right rows append in right-row order,
    /// retaining their factors and payloads. Given identical callback outcomes, plain and tracked pushout
    /// produce equal sets.
    ///
    /// # Panics
    ///
    /// Panics if the callback returns an id outside `self`, or returns the same id for
    /// multiple right rows and all combinations succeed. Callback panics propagate.
    pub fn pushout(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<Self>
    where
        D: Clone,
    {
        self.tracked_pushout(right, coincident, combine)
            .map(|(object, _)| object)
    }

    /// Perform [`Self::pushout`] and return its relation correspondences.
    ///
    /// Has the same result and callback contract as the plain operation, including `None`
    /// on payload-combination rejection. Both correspondences cover their inputs and share the
    /// result count as their target count.
    ///
    /// # Panics
    ///
    /// Panics under the callback conditions documented for [`Self::pushout`].
    pub fn tracked_pushout(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        mut combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<(Self, RelationPushoutCorrespondence)>
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

    /// Retain caller-paired rows from two sets in the same participant id space.
    ///
    /// Calls `coincident` on `right` for each left row; `None` omits that row. Matching
    /// must be injective and identify equal participant multisets in each factor; storage
    /// does not validate the pairing. `combine` receives both original factors and
    /// payloads, without frame alignment. Returns `None` as soon as `combine` rejects a pair.
    ///
    /// # Semantic properties
    ///
    /// For valid pairings, the result retains matched left rows in left-row order, with
    /// their factors, dense new ids, and the combined payloads. No matches yields
    /// `Some` of an empty set. Given identical callback outcomes, plain and tracked
    /// pullback produce equal sets.
    ///
    /// # Panics
    ///
    /// Panics if the callback returns an id outside `right`, or returns the same id for
    /// multiple left rows and all combinations succeed. Callback panics propagate.
    pub fn pullback(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<Self>
    where
        D: Clone,
    {
        self.tracked_pullback(right, coincident, combine)
            .map(|(object, _)| object)
    }

    /// Perform [`Self::pullback`] and return its relation correspondences.
    ///
    /// Has the same result and callback contract as the plain operation, including `None`
    /// on payload-combination rejection. Both correspondences cover the result and share its
    /// count as their source count.
    ///
    /// # Panics
    ///
    /// Panics under the callback conditions documented for [`Self::pullback`].
    pub fn tracked_pullback(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        mut combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<(Self, RelationPullbackCorrespondence)>
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
