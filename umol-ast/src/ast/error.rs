//! Error types for molecule AST operations.

use thiserror::Error;

use super::constraint::joint_domain::JointVar;
use super::ids::AtomId;

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

/// Error raised by `JointDomainAst::from_ints` (and sibling constructors) when
/// the input does not satisfy the joint-domain invariants: vars and tuples
/// must each have ≥ 2 entries, tuple arities must match the vars list, and
/// vars must be unique.
#[derive(Clone, Debug, PartialEq, Error)]
pub enum JointDomainError {
    #[error("vars must have at least 2 entries, got {0}")]
    TooFewVars(usize),
    #[error("tuples must have at least 2 distinct entries, got {0}")]
    TooFewTuples(usize),
    #[error("tuple {tuple_index} has {tuple_len} values but vars has {vars_len}")]
    ArityMismatch {
        tuple_index: usize,
        tuple_len: usize,
        vars_len: usize,
    },
    #[error("duplicate var {0:?} in vars list")]
    DuplicateVar(JointVar),
}

/// Signal raised by `Lattice::saturate` (and `saturate_atom` for `AtomAst`)
/// when cross-field constraint propagation reveals that no admissible value
/// assignment remains. `Lattice::meet` converts this to `None` at the
/// boundary so it propagates through the standard `Option<Self>` contract.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[error("constraint propagation reached a contradiction")]
pub struct Contradiction;
