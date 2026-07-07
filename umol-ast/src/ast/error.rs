//! Error types for molecule AST operations.

use thiserror::Error;

use super::id::AtomId;
use super::molecule::transact::TransactionError;

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

impl From<Contradiction> for ApplyError {
    fn from(_: Contradiction) -> Self {
        Self::Inconsistent
    }
}
