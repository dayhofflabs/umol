use thiserror::Error;
use umol_data::Element;

use super::molecule::{AtomIndex, BondIndex};
use crate::diagnostics::Diagnostic;

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ResolutionError {
    #[error("Atom index already exists: {0:?}")]
    DuplicateAtom(AtomIndex),
    #[error("Atom not found: {0:?}")]
    AtomNotFound(AtomIndex),
    #[error("Bond not found: {0:?}")]
    BondNotFound(BondIndex),
    #[error("Duplicate bond between atoms {0:?} and {1:?}")]
    DuplicateBond(AtomIndex, AtomIndex),
    #[error("Self-loop bond detected on atom {0:?}")]
    SelfLoop(AtomIndex),
    #[error("Invalid atom specification: {0}")]
    InvalidAtomSpec(String),
    #[error("Invalid bond specification: {0}")]
    InvalidBondSpec(String),
    #[error("Valence violation for element {0:?}: {1}")]
    ValenceViolation(Element, String),
    #[error("Unsupported element in GraphIR: {0:?}")]
    UnsupportedElement(Element),
    #[error("Specification parsing error: {0}")]
    SpecParseError(String),
}

impl From<ResolutionError> for Diagnostic {
    fn from(_error: ResolutionError) -> Self {
        todo!()
    }
}
