//! Error types for umol-data.

use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum DataError {
    #[error("invalid element: {0}")]
    InvalidElement(String),
    #[error("invalid isotope: {0}")]
    InvalidIsotope(String),
    #[error("invalid occupation: {0}")]
    InvalidOccupation(String),
}
