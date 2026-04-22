//! Error types for molecule AST operations.

use std::error::Error;
use std::fmt;

use super::idx::AtomIdx;

#[derive(Debug)]
pub enum RewriteError {
    DanglingEdge { atom: AtomIdx, neighbor: AtomIdx },
    DanglingRelation { atom: AtomIdx },
    UnmappedLhsAtom(AtomIdx),
    UnmappedAssignmentAtom(AtomIdx),
}

impl fmt::Display for RewriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DanglingEdge { atom, neighbor } => {
                write!(
                    f,
                    "dangling edge: atom {} -> neighbor {} not in L",
                    atom.0, neighbor.0
                )
            }
            Self::DanglingRelation { atom } => {
                write!(f, "dangling relation at atom {}", atom.0)
            }
            Self::UnmappedLhsAtom(a) => write!(f, "LHS atom {} not in assignment", a.0),
            Self::UnmappedAssignmentAtom(a) => {
                write!(f, "assignment atom {} not in target", a.0)
            }
        }
    }
}

impl Error for RewriteError {}
