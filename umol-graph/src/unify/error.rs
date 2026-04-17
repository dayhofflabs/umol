//! Error types for unification.

use thiserror::Error;

use crate::diagnostics::Diagnostic;
use crate::unify::aromaticity::AromaticityError;

impl From<ResolutionError> for Diagnostic {
    fn from(_error: ResolutionError) -> Self {
        todo!()
    }
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ResolutionError {
    #[error("resolution underdetermined")]
    Underdetermined,
    #[error("resolution contradictory")]
    Contradictory,
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ValidationError {
    #[error("validation underdetermined")]
    Underdetermined,
    #[error("validation contradictory")]
    Contradictory,
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ConfigError {
    #[error("Invalid atom type registry: {0}")]
    InvalidAtomTypeRegistry(String),
    #[error("Invalid valence table: {0}")]
    InvalidValenceTable(String),
}
