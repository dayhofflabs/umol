//! AST lowering and evaluation errors.

use thiserror::Error;
use umol_data::SpinStateError;

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

#[derive(Clone, Debug, PartialEq, Error)]
pub enum EvaluationError {
    #[error("Unbound variable: {0}")]
    UnboundVariable(String),
    #[error("Division by zero")]
    DivisionByZero,
    #[error("Type mismatch")]
    TypeMismatch,
}
