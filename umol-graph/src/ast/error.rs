//! AST errors.

use thiserror::Error;

#[derive(Clone, Debug, PartialEq, Error)]
#[error("molecule AST is not fully ground")]
pub struct GroundError;

#[derive(Clone, Debug, PartialEq, Error)]
pub enum LoweringError {
    #[error(transparent)]
    Ground(#[from] GroundError),
    #[error("{0}")]
    Custom(String),
}

impl From<String> for LoweringError {
    fn from(value: String) -> Self {
        Self::Custom(value)
    }
}

impl From<&str> for LoweringError {
    fn from(value: &str) -> Self {
        Self::Custom(value.to_string())
    }
}
