//! Error types for molecule AST operations.

use thiserror::Error;

/// Signal that a value is unsatisfiable — no admissible assignment remains.
/// Raised by fallible canonicalization/construction (e.g. an empty set);
/// `Lattice::meet` surfaces the same condition as `None`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[error("reached a contradiction")]
pub struct Contradiction;
