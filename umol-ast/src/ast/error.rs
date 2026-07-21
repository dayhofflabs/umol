//! Error types for molecule AST operations.

use thiserror::Error;

use super::id::AtomId;
use super::molecule::transact::TransactionError;
use super::validate::{DpoContradiction, EntityStructureContradiction};

/// Unsatisfiable operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[error("reached a contradiction")]
pub struct Contradiction;

/// No join: elements have no least upper bound (meet-semilattice).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[error("no join: elements have no least upper bound")]
pub struct NoJoin;

/// Error from applying a reaction onto a host molecule (`ReactionAst::apply_at`).
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ApplyError {
    /// DPO gluing condition is violated.
    #[error("dangling edge at deleted host atom {host_atom}")]
    Dangling { host_atom: AtomId },
    /// The reaction's deltas are inconsistent (canonicalization failed).
    #[error("inconsistent reaction deltas")]
    Inconsistent,
    /// Structural conflict: a parallel bond, overlapping systems, two stereo centers on one site, etc.
    #[error("applied product has a structural conflict")]
    StructuralConflict,
    /// The lowered edit transaction failed against the host.
    #[error("apply transaction failed: {0}")]
    Transaction(#[from] TransactionError),
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
    /// The reaction's delta sequence cannot be canonicalized.
    #[error("inconsistent reaction deltas")]
    InconsistentReaction,
    /// The reaction left-hand side is structurally invalid.
    #[error("invalid reaction lhs: {0}")]
    ReactionStructure(EntityStructureContradiction),
    /// The reaction violates its rule-local DPO invariant.
    #[error("invalid reaction: {0}")]
    ReactionDpo(DpoContradiction),
    /// The host molecule is structurally invalid.
    #[error("invalid host: {0}")]
    HostStructure(EntityStructureContradiction),
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

    #[rstest]
    #[case::dangling(ApplyError::Dangling { host_atom: AtomId(3) }, true)]
    #[case::structural_conflict(ApplyError::StructuralConflict, true)]
    #[case::inconsistent(ApplyError::Inconsistent, false)]
    #[case::transaction(ApplyError::Transaction(TransactionError::OldStateMismatch), false)]
    fn test_apply_error_is_match_rejection(#[case] error: ApplyError, #[case] expected: bool) {
        assert_eq!(error.is_match_rejection(), expected);
    }
}
