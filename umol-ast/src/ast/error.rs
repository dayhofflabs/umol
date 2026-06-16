//! Error types for molecule AST operations.

use thiserror::Error;

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

/// Signal that a value is unsatisfiable — no admissible assignment remains.
/// Raised by fallible canonicalization/construction (e.g. an empty set);
/// `Lattice::meet` surfaces the same condition as `None`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[error("reached a contradiction")]
pub struct Contradiction;
