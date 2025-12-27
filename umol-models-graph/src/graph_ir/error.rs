use crate::diagnostics::{Diagnostic, DiagnosticKind, Severity};
use crate::graph_ir::{AtomIndex, BondIndex};
use thiserror::Error;
use umol_data::Element;

#[derive(Debug, Error, Clone, PartialEq)]
pub enum GraphError {
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
    #[error("Conversion from TableIR to GraphIR failed: {0}")]
    ConversionFailed(String),
    #[error("Unsupported element in GraphIR: {0:?}")]
    UnsupportedElement(Element),
    #[error("Specification parsing error: {0}")]
    SpecParseError(String),
}

impl From<GraphError> for Diagnostic {
    fn from(error: GraphError) -> Self {
        use DiagnosticKind::*;
        let (kind, details) = match error {
            GraphError::DuplicateAtom(_) => (Unknown, Some(error.to_string())),
            GraphError::AtomNotFound(_) => (Unknown, Some(error.to_string())),
            GraphError::BondNotFound(_) => (Unknown, Some(error.to_string())),
            GraphError::DuplicateBond(_, _) => (GraphTopologyParallelEdges, Some(error.to_string())),
            GraphError::SelfLoop(_) => (GraphTopologySelfLoopRing, Some(error.to_string())),
            GraphError::InvalidAtomSpec(ref s) => (Unknown, Some(s.clone())),
            GraphError::InvalidBondSpec(ref s) => (Unknown, Some(s.clone())),
            GraphError::ValenceViolation(_, ref s) => {
                (GraphValenceOutOfElementRange, Some(s.clone()))
            }
            GraphError::ConversionFailed(ref s) => (GraphConversionUnknown, Some(s.clone())),
            GraphError::UnsupportedElement(_) => (Unknown, Some(error.to_string())),
            GraphError::SpecParseError(ref s) => (Unknown, Some(s.clone())),
        };
        Diagnostic {
            kind,
            category: kind.category(),
            severity: Severity::Error,
            span: None,
            details,
        }
    }
}

impl From<GraphError> for umol::error::ParseError {
    fn from(error: GraphError) -> Self {
        umol::error::ParseError::Format(Box::new(error))
    }
}

impl From<GraphError> for umol::Error {
    fn from(error: GraphError) -> Self {
        umol::Error::Parse(error.into())
    }
}
