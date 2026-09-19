//! Array-backed fixed-arity relation storage: [FixedRelationSet].

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

/// Fixed-arity relations over typed participants, with opaque payloads.
///
/// Rows retain input order, participant multiplicity, and payloads. Coinciding rows are
/// permitted. Participant references are not checked against an external graph.
/// Participant rows are stored as arrays; fixed arity may be zero.
/// Relation ids use u32 indices.
///
/// # Semantic properties
///
/// - `new(entries).into_entries() == entries` for representable sizes.
/// - Equality compares stored sequences and payloads, including their order.
/// - Incidence lists each relation once per referenced node or edge, in relation-id order.
/// - Coincidence compares participant multisets in the factor; it ignores stored order
///   and payloads, but observes multiplicity and complete participant values.
///
/// Construction/query, mutation, and transport laws are exercised through public APIs in
/// `tests/property/relation.rs`; restoration laws in `tests/property/restore.rs`.
/// Transport and restoration laws require conforming
/// [`RelationParticipant`] implementations.
#[derive(Clone, Debug)]
pub struct FixedRelationSet<P, D, const N: usize> {
    participants: Vec<[P; N]>,
    data: Vec<D>,
    incidence: Incidence,
}

impl<P: RelationParticipant, D, const N: usize> FixedRelationSet<P, D, N> {
    /// Store entries in order and build their union incidence index.
    ///
    /// Preserves the supplied sequence in the factor, including repeated values.
    /// Node/edge references come from [`RelationParticipant::refs`]; graph membership is
    /// not validated. Payloads are stored without interpretation.
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

    /// Consume the set into its stored entries, in relation-id order.
    ///
    /// Preserves participant order and multiplicity and returns the current payloads.
    pub fn into_entries(self) -> Vec<([P; N], D)> {
        self.participants.into_iter().zip(self.data).collect()
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
    /// Only payloads are mutable; participant storage and its derived index stay intact.
    pub fn iter_mut(&mut self) -> impl ExactSizeIterator<Item = (RelationId, &[P; N], &mut D)> {
        let participants = &self.participants;
        self.data
            .iter_mut()
            .enumerate()
            .map(move |(index, data)| (RelationId(index as u32), &participants[index], data))
    }

    /// Borrow the stored participant sequence of relation `id`.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn participants(&self, id: RelationId) -> &[P; N] {
        &self.participants[id.index()]
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
    /// Returns an empty slice when no participant references `node`.
    pub fn incident_to_node(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident_to_node(node)
    }

    /// Relations referencing `edge`, once each in ascending relation-id order.
    ///
    /// Returns an empty slice when no participant references `edge`.
    pub fn incident_to_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_to_edge(edge)
    }

    /// Whether any participant references `node`.
    pub fn has_incident_to_node(&self, node: NodeId) -> bool {
        self.incidence.has_incident_to_node(node)
    }

    /// Whether any participant references `edge`.
    pub fn has_incident_to_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_to_edge(edge)
    }

    /// Find the first relation incident with `node` whose participants match `query`.
    ///
    /// Compares complete participant multisets in the factor, ignoring stored order and
    /// payloads. Returns the smallest matching relation id, or `None` if the incidence
    /// list contains no match. Duplicate matching rows are permitted. An unreferenced
    /// anchor yields `None`, even if a matching row exists elsewhere.
    pub fn coincident_to_node(&self, node: NodeId, query: &[P]) -> Option<RelationId> {
        self.coincident_among(self.incident_to_node(node), query)
    }

    /// Find the first relation incident with `edge` whose participants match `query`.
    ///
    /// Uses the same multiset comparison and first-match rule as [`Self::coincident_to_node`].
    /// Returns `None` when the edge incidence list contains no matching row.
    pub fn coincident_to_edge(&self, edge: EdgeId, query: &[P]) -> Option<RelationId> {
        self.coincident_among(self.incident_to_edge(edge), query)
    }

    /// Whether relation `id` matches the query multiset.
    ///
    /// Compares complete participant values and multiplicities in the factor, without
    /// requiring a node or edge anchor. Stored order and payloads do not affect the result.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set.
    pub fn is_coincident(&self, id: RelationId, query: &[P]) -> bool {
        self.coincident_among(&[id], query).is_some()
    }

    /// Return the first candidate whose stored participant multiset matches the query.
    fn coincident_among(&self, candidates: &[RelationId], query: &[P]) -> Option<RelationId> {
        let mut sorted_query: Vec<P> = query.to_vec();
        sorted_query.sort_unstable();
        candidates
            .iter()
            .copied()
            .find(|&id| participants_match(self.participants(id), &sorted_query))
    }

    /// Append a relation and return its id, equal to the previous row count.
    ///
    /// As in [`Self::new`], order and repeated participants are preserved, coinciding rows
    /// are permitted, and references are indexed without checking external graph membership.
    /// Fixed arity is enforced by the array type; zero arity is permitted.
    /// The payload is moved into storage. Incidence is rebuilt over the entire collection.
    ///
    /// # Semantic properties
    ///
    /// Existing ids, participants, and payloads remain unchanged. Incidence lists the new
    /// relation once for each referenced node or edge. Removing the appended relation
    /// immediately afterward restores the previous set. These laws are exercised against
    /// a row model in `tests/property/relation.rs`.
    pub fn add(&mut self, participants: [P; N], data: D) -> RelationId {
        let id = RelationId(self.count() as u32);
        self.participants.push(participants);
        self.data.push(data);
        self.incidence = Incidence::build(self.count(), |i, out| {
            out.extend(
                self.participants[i]
                    .iter()
                    .map(|participant| participant.refs()),
            );
        });
        id
    }

    /// Remove whole relations, preserving survivor order and making their ids dense.
    ///
    /// Delegates to [`Self::tracked_remove`] and discards the compaction. Input ids refer
    /// to the pre-removal set; their order and repetitions do not matter. Empty input
    /// leaves the set unchanged. Participant references and surviving payloads are unchanged.
    ///
    /// # Panics
    ///
    /// Panics before mutation if any supplied id is outside the set.
    pub fn remove(&mut self, ids: &[RelationId]) {
        self.tracked_remove(ids);
    }

    /// Remove whole relations and return the old-to-new survivor compaction.
    ///
    /// Input ids refer to the pre-removal set; their order and repetitions do not matter.
    /// The compaction records the original count and removed ids, not removed payloads.
    /// A nonempty removal rebuilds incidence over the surviving rows once.
    ///
    /// # Semantic properties
    ///
    /// Survivors retain their relative order, participant sequences, and payloads without
    /// cloning. Only relation ids change; referenced nodes and edges are not removed or
    /// relabeled. Incidence uses the new dense relation ids. Empty input returns identity
    /// over the original count and leaves the set unchanged. Discarding the compaction
    /// produces the same state as [`Self::remove`]. These laws are exercised against
    /// a row model in `tests/property/relation.rs`.
    ///
    /// # Panics
    ///
    /// Panics before mutation if any supplied id is outside the set.
    pub fn tracked_remove(&mut self, ids: &[RelationId]) -> Compaction<RelationId> {
        let compaction = Compaction::new(self.count(), ids.to_vec())
            .expect("removed relations belong to the source set");
        if ids.is_empty() {
            return compaction;
        }
        let mut dst = 0;
        for src in 0..self.count() {
            if compaction
                .removed()
                .binary_search(&RelationId(src as u32))
                .is_err()
            {
                self.participants.swap(dst, src);
                self.data.swap(dst, src);
                dst += 1;
            }
        }
        self.participants.truncate(dst);
        self.data.truncate(dst);
        self.incidence = Incidence::build(self.count(), |i, out| {
            out.extend(
                self.participants[i]
                    .iter()
                    .map(|participant| participant.refs()),
            );
        });
        compaction
    }

    /// Undo row removal using its compaction and original saved entries.
    ///
    /// Saved entries carry original relation ids and may arrive in any order. Their
    /// participants already use the desired reference space and are not translated.
    /// Payloads move without cloning. Storage appends the saved rows, reorders both
    /// columns in place using a temporary destination array, and rebuilds incidence.
    ///
    /// The compaction and saved entries must come from the removal being undone.
    /// With conforming [`RelationParticipant`] implementations, manipulated or mismatched
    /// undo data does not panic; the resulting set is unspecified.
    ///
    /// # Semantic properties
    ///
    /// Matching removal and restoration recover the original rows, ids, and incidence.
    /// Surviving participant sequences and current payloads are preserved. Reordering
    /// saved entries does not affect the result. Matching identity is a no-op, including
    /// for zero arity. Public-API properties in `tests/property/restore.rs` check these
    /// laws against independent rows and incidence scans, plus freedom from panics.
    pub fn restore(
        &mut self,
        compaction: &Compaction<RelationId>,
        removed: Vec<(RelationId, [P; N], D)>,
    ) {
        if compaction.removed().is_empty() {
            return;
        }
        let Some(count) = self.count().checked_add(removed.len()) else {
            return;
        };
        if compaction.source_count() > u32::MAX as usize || count > u32::MAX as usize {
            return;
        }
        let mut destinations = Vec::new();
        if destinations.try_reserve_exact(count).is_err()
            || self.participants.try_reserve(removed.len()).is_err()
            || self.data.try_reserve(removed.len()).is_err()
        {
            return;
        }
        let mut removed_ids = compaction.removed().iter().peekable();
        let mut original = 0;
        for _ in 0..self.count() {
            while removed_ids.peek().is_some_and(|id| id.index() == original) {
                removed_ids.next();
                original += 1;
            }
            if original >= compaction.source_count() {
                return;
            }
            destinations.push(original);
            original += 1;
        }
        for (id, participants, data) in removed {
            destinations.push(id.index());
            self.participants.push(participants);
            self.data.push(data);
        }
        for source in 0..count {
            while destinations[source] != source {
                let target = destinations[source];
                if target >= count || destinations[target] == target {
                    break;
                }
                self.participants.swap(source, target);
                self.data.swap(source, target);
                destinations.swap(source, target);
            }
        }
        self.incidence = Incidence::build(self.count(), |i, out| {
            out.extend(self.participants[i].iter().map(|p| p.refs()));
        });
    }

    /// Undo graph-reference compaction in the surviving participant rows.
    ///
    /// Expands node and edge references in place and rebuilds incidence. Payloads are
    /// not moved or cloned. Restore participant references before restoring removed
    /// rows, whose saved participants already use original graph ids.
    ///
    /// The graph compaction must match the participant compaction being undone.
    /// With conforming [`RelationParticipant`] implementations, manipulated or mismatched
    /// undo data does not panic; the resulting set is unspecified.
    ///
    /// # Semantic properties
    ///
    /// Matching compaction and restoration recover surviving participant values,
    /// including their order, multiplicity, and non-reference data. Relation ids,
    /// arity, and payloads are unchanged. Matching identity and empty storage are
    /// no-ops. Public-API properties in `tests/property/restore.rs` check independent
    /// inverse-reference expectations and incidence, roundtrips, and freedom from panics.
    pub fn restore_participants(&mut self, compaction: &GraphCompaction) {
        if compaction.nodes().removed().is_empty() && compaction.edges().removed().is_empty() {
            return;
        }
        if compaction.nodes().source_count().saturating_sub(1) > u32::MAX as usize
            || compaction.edges().source_count().saturating_sub(1) > u32::MAX as usize
        {
            return;
        }
        for row in &mut self.participants {
            for participant in row {
                let refs = participant.refs();
                if refs
                    .node
                    .is_some_and(|id| id.index() >= compaction.nodes().result_count())
                    || refs
                        .edge
                        .is_some_and(|id| id.index() >= compaction.edges().result_count())
                {
                    continue;
                }
                *participant = participant.uncompact(compaction);
            }
        }
        self.incidence = Incidence::build(self.count(), |i, out| {
            out.extend(self.participants[i].iter().map(|p| p.refs()));
        });
    }

    /// Reorder relation `id`'s participants so that `new[i] = old[order[i]]`.
    ///
    /// # Semantic properties
    ///
    /// The selected multiset, all payloads, and all other rows are unchanged.
    /// Incidence is unchanged and is not rebuilt. Identity leaves the set equal to itself;
    /// applying a permutation followed by its inverse recovers the original sequence.
    ///
    /// # Panics
    ///
    /// Panics if `id` is outside the set or `order` is not a permutation of the selected
    /// factor's positions (wrong length, repeated position, or out-of-range position).
    pub fn permute_participants(&mut self, id: RelationId, order: &[ParticipantPosition]) {
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
        self.map_participants(|participant| participant.try_map(correspondence))
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
        self.map_participants(|participant| Some(participant.remap(remapping)))
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
            .all(|id| remappable_under(self.participants(id), remapping))
            .then(|| self.remap(remapping))
    }

    /// Map every factor while retaining row order and payloads; reject the set on any `None`.
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

    /// Glue two relation sets in the same participant id space using caller-selected pairings.
    ///
    /// Calls `coincident` on the original `self` for each right row. The callback must
    /// identify a matching row, or return `None` for an unmatched row. Matches must be
    /// injective; storage does not verify participant equality or payload compatibility.
    /// `combine` receives both original participant sequences and payloads, without frame alignment.
    /// Returns `None` as soon as `combine` rejects a pair; otherwise returns the glued set.
    ///
    /// # Semantic properties
    ///
    /// For valid pairings, `self` keeps its row ids and participant sequence. Matched payloads are
    /// replaced by `combine`'s result. Unmatched right rows append in right-row order,
    /// retaining their participant sequence and payloads. Given identical callback outcomes,
    /// plain and tracked pushout produce equal sets.
    ///
    /// # Panics
    ///
    /// Panics if the callback returns an id outside `self`, or returns the same id for
    /// multiple right rows and all combinations succeed. Callback panics propagate.
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

    /// Retain caller-paired rows from two sets in the same participant id space.
    ///
    /// Calls `coincident` on `right` for each left row; `None` omits that row. Matching
    /// must be injective and identify equal participant multisets in the factor; storage
    /// does not validate the pairing. `combine` receives both original participant sequences and
    /// payloads, without frame alignment. Returns `None` as soon as `combine` rejects a pair.
    ///
    /// # Semantic properties
    ///
    /// For valid pairings, the result retains matched left rows in left-row order, with
    /// their participant sequence, dense new ids, and the combined payloads. No matches yields
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
        coincident: impl Fn(&Self, &[P]) -> Option<RelationId>,
        combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
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
