//! Error types for molecule graph-IR operations.

use thiserror::Error;

use super::entity::{Entity, EntityKind};
use super::id::AtomId;
use super::molecule::transact::TransactionError;
use super::validate::DpoContradiction;

/// Unsatisfiable operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[error("reached a contradiction")]
pub struct Contradiction;

/// No join: elements have no least upper bound (meet-semilattice).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[error("no join: elements have no least upper bound")]
pub struct NoJoin;

/// Error from applying a reaction onto a host molecule (`Reaction::apply_at`).
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ApplyError {
    /// DPO gluing condition is violated.
    #[error("dangling edge at deleted host atom {host_atom}")]
    Dangling { host_atom: AtomId },
    /// The reaction's deltas are inconsistent (normalization failed).
    #[error("inconsistent reaction deltas")]
    Inconsistent,
    /// Structural conflict: a parallel bond, overlapping systems, two stereo centers on one site, etc.
    #[error("applied product has a structural conflict")]
    StructuralConflict,
    /// The lowered edit transaction failed against the host.
    #[error("apply transaction failed: {0}")]
    Transaction(#[from] TransactionError),
    /// The supplied or matcher-produced correspondence does not map an entity consistently into
    /// the host.
    #[error("application correspondence does not map {entity:?} consistently into the host")]
    CorrespondenceMismatch { entity: Entity },
    /// A matched stereo entity's rule and host ligand frames are not orderings of the same set.
    #[error("matched stereo frame differs for {entity:?}")]
    StereoFrameMismatch { entity: Entity },
    /// A condition established by reaction integrity validation or matching did not hold during
    /// lowering.
    #[error("application reached an internal invariant failure")]
    InternalInvariant,
}

impl ApplyError {
    /// Whether this failure rejects only the current pattern embedding.
    pub fn is_match_rejection(&self) -> bool {
        matches!(self, Self::Dangling { .. } | Self::StructuralConflict)
    }
}

/// A reaction or host does not satisfy the structural preconditions for application.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ApplyPreconditionError {
    /// The reaction's delta sequence cannot be normalized.
    #[error("inconsistent reaction deltas")]
    InconsistentReaction,
    /// The reaction left-hand side violates an entity-structure invariant required by application.
    #[error("reaction lhs violates the application invariant for {kind}")]
    ReactionStructureInvariant { kind: EntityKind },
    /// The reaction violates its rule-local DPO invariant.
    #[error("invalid reaction: {0}")]
    ReactionDpo(DpoContradiction),
    /// The host violates an entity-structure invariant required by application.
    #[error("host violates the application invariant for {kind}")]
    HostStructureInvariant { kind: EntityKind },
    /// A delta or molecule constraint references an entity unavailable on the reaction LHS or among
    /// the entities created by the reaction.
    #[error("reaction references unavailable entity {entity:?}")]
    InvalidReactionReference { entity: Entity },
    /// A delta's endpoints, participants, site, or ligands disagree with its identified LHS entity.
    #[error("reaction incidence does not match lhs entity {entity:?}")]
    ReactionIncidenceMismatch { entity: Entity },
}

impl From<Contradiction> for ApplyError {
    fn from(_: Contradiction) -> Self {
        Self::Inconsistent
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::ir::id::{BondId, StereoAtomId};

    #[rstest]
    #[case::dangling(ApplyError::Dangling { host_atom: AtomId(3) }, true)]
    #[case::structural_conflict(ApplyError::StructuralConflict, true)]
    #[case::inconsistent(ApplyError::Inconsistent, false)]
    #[case::transaction(ApplyError::Transaction(TransactionError::OldStateMismatch), false)]
    #[case::correspondence(
        ApplyError::CorrespondenceMismatch { entity: Entity::Bond(BondId(2)) },
        false,
    )]
    #[case::stereo_frame(
        ApplyError::StereoFrameMismatch {
            entity: Entity::StereoAtom(StereoAtomId(1)),
        },
        false,
    )]
    #[case::internal(ApplyError::InternalInvariant, false)]
    fn test_apply_error_is_match_rejection(#[case] error: ApplyError, #[case] expected: bool) {
        assert_eq!(error.is_match_rejection(), expected);
    }

    #[rstest]
    #[case::correspondence(
        ApplyError::CorrespondenceMismatch { entity: Entity::Atom(AtomId(3)) },
        "application correspondence does not map Atom(AtomId(3)) consistently into the host",
    )]
    #[case::stereo_frame(
        ApplyError::StereoFrameMismatch {
            entity: Entity::StereoAtom(StereoAtomId(1)),
        },
        "matched stereo frame differs for StereoAtom(StereoAtomId(1))",
    )]
    #[case::internal(
        ApplyError::InternalInvariant,
        "application reached an internal invariant failure"
    )]
    fn test_apply_error_display(#[case] error: ApplyError, #[case] expected: &str) {
        assert_eq!(error.to_string(), expected);
    }

    #[rstest]
    #[case::invalid_reference(
        ApplyPreconditionError::InvalidReactionReference {
            entity: Entity::Atom(AtomId(3)),
        },
        "reaction references unavailable entity Atom(AtomId(3))",
    )]
    #[case::incidence(
        ApplyPreconditionError::ReactionIncidenceMismatch {
            entity: Entity::Bond(BondId(2)),
        },
        "reaction incidence does not match lhs entity Bond(BondId(2))",
    )]
    fn test_apply_precondition_error_display(
        #[case] error: ApplyPreconditionError,
        #[case] expected: &str,
    ) {
        assert_eq!(error.to_string(), expected);
    }
}
