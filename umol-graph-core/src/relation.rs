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

use std::ops::{Add, Sub};

pub use self::fixed::FixedRelationSet;
pub use self::fixed_fixed::FixedFixedBirelationSet;
pub use self::fixed_var::FixedVarBirelationSet;
pub use self::participant::{ParticipantPosition, ParticipantRefs, RelationParticipant};
pub use self::var::VarRelationSet;
pub use self::var_var::VarVarBirelationSet;
use crate::correspondence::Correspondence;
use crate::remap::GraphRemapping;

mod fixed;
mod fixed_fixed;
mod fixed_var;
mod incidence;
mod participant;
mod var;
mod var_var;

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

/// The two input-to-result correspondences of a same-space relation-set pushout.
///
/// Operation-produced components have equal target counts and cover their respective
/// inputs. Public fields may be assembled independently; agreement with a particular
/// resulting relation set is contextual.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelationPushoutCorrespondence {
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
) -> (S, RelationPushoutCorrespondence) {
    let left: Vec<RelationId> = (0..self_count).map(RelationId::from).collect();
    (
        object,
        RelationPushoutCorrespondence {
            left: Correspondence::from_images(&left, object_count),
            right: Correspondence::from_images(&right_map, object_count),
        },
    )
}

/// The two result-to-input projections of a same-space relation-set pullback.
///
/// Operation-produced components cover the result and have equal source counts. Public fields may be
/// assembled independently; agreement with the result and input sets is contextual.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelationPullbackCorrespondence {
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
) -> (S, RelationPullbackCorrespondence) {
    (
        object,
        RelationPullbackCorrespondence {
            left: Correspondence::from_images(&left_images, self_count),
            right: Correspondence::from_images(&right_images, right_count),
        },
    )
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

#[cfg(test)]
mod tests;
