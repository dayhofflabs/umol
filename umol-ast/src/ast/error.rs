//! Error types for molecule AST operations.

use thiserror::Error;

use super::idx::AtomIdx;

/// Error raised by `MoleculeAst::rewrite` when the L / R / assignment
/// triple is not self-consistent.
#[derive(Clone, Debug, PartialEq, Error)]
pub enum RewriteError {
    #[error("dangling edge: atom {atom} -> neighbor {neighbor} not in L")]
    DanglingEdge { atom: AtomIdx, neighbor: AtomIdx },
    #[error("dangling relation at atom {atom}")]
    DanglingRelation { atom: AtomIdx },
    #[error("LHS atom {0} not in assignment")]
    UnmappedLhsAtom(AtomIdx),
    #[error("assignment atom {0} not in target")]
    UnmappedAssignmentAtom(AtomIdx),
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
