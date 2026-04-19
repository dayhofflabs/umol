//! Error types for unification.

use thiserror::Error;

use crate::diagnostics::Diagnostic;
use crate::ops::aromaticity::AromaticityError;

impl From<ResolutionError> for Diagnostic {
    fn from(_error: ResolutionError) -> Self {
        todo!()
    }
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ResolutionError {
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
}

#[derive(Debug, Error, Clone, PartialEq)]
pub enum ConfigError {
    #[error("Invalid atom type registry: {0}")]
    InvalidAtomTypeRegistry(String),
    #[error("Invalid valence table: {0}")]
    InvalidValenceTable(String),
}
