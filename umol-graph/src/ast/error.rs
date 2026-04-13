//! AST lowering errors.

use thiserror::Error;
use umol_shared::error::SpinStateError;

#[derive(Clone, Debug, PartialEq, Error)]
#[error("molecule AST is not fully ground")]
pub struct GroundError;

#[derive(Clone, Debug, PartialEq, Error)]
pub enum LoweringError {
    #[error("non-ground value for field '{field}'")]
    NonGround { field: &'static str },
    #[error("value {value} out of range for field '{field}'")]
    OutOfRange { field: &'static str, value: i64 },
    #[error("field '{field}' is required but not present")]
    MissingField { field: &'static str },
    #[error("invalid spin multiplicity: {0}")]
    InvalidMultiplicity(u8),
    #[error("incompatible spin state: {0}")]
    SpinState(#[from] SpinStateError),
    #[error("invalid atom spec: {0}")]
    Atom(String),
    #[error("unknown atom label: {0}")]
    UnknownLabel(String),
    #[error("invalid molecule spec: {0}")]
    Molecule(String),
}
