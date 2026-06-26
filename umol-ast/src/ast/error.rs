//! Error types for molecule AST operations.

use thiserror::Error;

use super::id::AtomId;
use super::molecule::transact::TransactionError;

/// Signal that a value is unsatisfiable — no admissible assignment remains.
/// Raised by fallible canonicalization/construction (e.g. an empty set);
/// `Lattice::meet` surfaces the same condition as `None`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[error("reached a contradiction")]
pub struct Contradiction;

/// Error from applying a reaction onto a host molecule (`ReactionAst::apply_at`).
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ApplyError {
    /// A deleted host atom retains an incident bond the rule does not also delete — the DPO
    /// gluing condition is violated, so the match is not a valid application site.
    #[error("dangling edge at deleted host atom {host_atom}")]
    Dangling { host_atom: AtomId },
    /// The reaction's deltas are inconsistent (canonicalization failed).
    #[error("inconsistent reaction deltas")]
    Inconsistent,
    /// The lowered edit transaction failed against the host.
    #[error("apply transaction failed: {0}")]
    Transaction(#[from] TransactionError),
}

impl From<Contradiction> for ApplyError {
    fn from(_: Contradiction) -> Self {
        Self::Inconsistent
    }
}
