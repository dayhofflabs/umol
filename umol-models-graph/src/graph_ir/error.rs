use thiserror::Error;
use umol_data::Element;

use crate::diagnostics::Diagnostic;

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ResolutionError {
    #[error("Invalid atom specification: {0}")]
    InvalidAtomSpec(String),
    #[error("Invalid bond specification: {0}")]
    InvalidBondSpec(String),
    #[error("Valence violation for element {0:?}: {1}")]
    ValenceViolation(Element, String),
    #[error("No valence match for {0}")]
    ValenceNoMatch(String),
    #[error("Valence ambiguous for {0}")]
    ValenceAmbiguous(String),
}

impl From<ResolutionError> for Diagnostic {
    fn from(_error: ResolutionError) -> Self {
        todo!()
    }
}
