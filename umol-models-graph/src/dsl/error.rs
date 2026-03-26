//! Domain errors for DSL parsing.

use nom::error::{ErrorKind as NomErrorKind, ParseError as NomParseError};
use thiserror::Error;

#[derive(Clone, Debug, PartialEq, Error)]
pub enum ParseError {
    #[error("Invalid number")]
    InvalidNumber,
    #[error("Invalid bond order")]
    InvalidBondOrder,
    #[error("Unknown bond predicate")]
    UnknownBondPredicate,
    #[error("Duplicate {0} bond predicate")]
    DuplicateBondPredicate(String),
    #[error("Trailing content in bond string")]
    TrailingContent,
    #[error("Incomplete input")]
    Incomplete,
    #[error("Invalid value DSL: {0}")]
    InvalidValueDsl(String),
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
