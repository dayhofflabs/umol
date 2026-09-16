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

    /// Same-space relation pushout — see [`FixedRelationSet::pushout`](crate::FixedRelationSet::pushout). Coincidence is equality of
    /// both factors' participants.
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

    /// Glue relation sets and return both input-to-result mappings with the result.
    ///
    /// Has the same result and failure behavior as [`Self::pushout`].
    pub fn tracked_pushout(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        mut combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<(Self, RelationPushoutCorrespondence)>
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

    /// Same-space relation pullback — see [`FixedRelationSet::pullback`](crate::FixedRelationSet::pullback).
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

    /// Return the shared relation set and its two result-to-input projections.
    ///
    /// Has the same result and failure behavior as [`Self::pullback`].
    pub fn tracked_pullback(
        &self,
        right: &Self,
        coincident: impl Fn(&Self, &[L1], &[L2]) -> Option<RelationId>,
        mut combine: impl FnMut((&[L1], &[L2], &D), (&[L1], &[L2], &D)) -> Option<D>,
    ) -> Option<(Self, RelationPullbackCorrespondence)>
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
