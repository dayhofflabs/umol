//! Error types for molecule AST operations.

use thiserror::Error;

use super::idx::AtomId;

/// Error raised by `MoleculeAst::rewrite` when the L / R / assignment
/// triple is not self-consistent.
#[derive(Clone, Debug, PartialEq, Error)]
pub enum RewriteError {
    #[error("dangling edge: atom {atom} -> neighbor {neighbor} not in L")]
    DanglingEdge { atom: AtomId, neighbor: AtomId },
    #[error("dangling relation at atom {atom}")]
    DanglingRelation { atom: AtomId },
    #[error("LHS atom {0} not in assignment")]
    UnmappedLhsAtom(AtomId),
    #[error("assignment atom {0} not in target")]
    UnmappedAssignmentAtom(AtomId),
}

/// Error raised by `Expr::evaluate` / `Expr::evaluate_bool` when the
/// expression cannot be reduced (unbound variable, division by zero, or
/// arithmetic-vs-boolean domain mismatch).
#[derive(Clone, Debug, PartialEq, Error)]
pub enum EvaluationError {
    #[error("Unbound variable: {0}")]
    UnboundVariable(String),
    #[error("Division by zero")]
    DivisionByZero,
    #[error("Type mismatch")]
    TypeMismatch,
}
