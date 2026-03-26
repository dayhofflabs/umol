//! Domain errors for DSL parsing.

use nom::error::{ErrorKind as NomErrorKind, ParseError as NomParseError};
use thiserror::Error;

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
    #[error("Nom error: {0:?}")]
    NomError(NomErrorKind),
}

impl<I> NomParseError<I> for ParseError {
    fn from_error_kind(_input: I, kind: NomErrorKind) -> Self {
        ParseError::NomError(kind)
    }

    fn append(_input: I, _kind: NomErrorKind, other: Self) -> Self {
        other
    }
}
