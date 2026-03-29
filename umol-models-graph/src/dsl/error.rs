//! Domain errors for DSL parsing.

use nom::error::{ErrorKind as NomErrorKind, ParseError as NomParseError};
use thiserror::Error;
use umol_data::SpinStateError;

#[derive(Clone, Debug, PartialEq, Error)]
pub enum ParseError {
    #[error("Invalid number")]
    InvalidNumber,
    #[error("Invalid bond order: {0}")]
    InvalidBondOrder(String),
    #[error("Unknown bond predicate: {0}")]
    UnknownBondPredicate(String),
    #[error("Duplicate {0} bond predicate")]
    DuplicateBondPredicate(String),
    #[error("Invalid bond data {0}")]
    InvalidBondData(String),
    #[error("Incomplete input")]
    Incomplete,
    #[error("Trailing input: {0:?}")]
    TrailingInput(String),
    #[error("Invalid value: {0}")]
    InvalidValue(String),
    #[error("Invalid atom element: {0}")]
    InvalidElement(String),
    #[error("Unknown atom predicate: {0}")]
    UnknownAtomPredicate(String),
    #[error("Duplicate {0} atom predicate")]
    DuplicateAtomPredicate(String),
    #[error("Nom error: {0:?}")]
    NomError(NomErrorKind),
    #[error("EDN parse error: {0}")]
    EdnParse(String),
    #[error("Missing required key: {0}")]
    MissingKey(String),
    #[error("expected {expected} for :{field}")]
    WrongFieldType {
        field: String,
        expected: String,
    },
    #[error("invalid atom DSL: {0}")]
    InvalidAtomDsl(String),
    #[error("invalid bond DSL: {0}")]
    InvalidBondDsl(String),
    #[error("invalid bond entry: expected map-based {{[:id keyword] :a :b :bond}} or vector-based [a b bond-spec]")]
    InvalidBond,
    #[error("Duplicate structural id: {0}")]
    DuplicateId(String),
    #[error("Unknown atom endpoint: {0}")]
    InvalidAtomIndex(String),
    #[error("Unknown alias: {0}")]
    UnknownAlias(String),
    #[error("Invalid spin state: {0}")]
    InvalidSpinState(#[from] SpinStateError),
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

impl<I> NomParseError<I> for ParseError {
    fn from_error_kind(_input: I, kind: NomErrorKind) -> Self {
        ParseError::NomError(kind)
    }

    fn append(_input: I, _kind: NomErrorKind, other: Self) -> Self {
        other
    }
}

#[derive(Clone, Debug, PartialEq, Error)]
pub enum LoweringError {}
