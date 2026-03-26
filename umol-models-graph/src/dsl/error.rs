//! Domain errors for DSL parsing.

use nom::error::Error as NomError;
use nom::Err;
use thiserror::Error;

#[derive(Clone, Debug, PartialEq, Error)]
pub enum ParseError {
    #[error("Invalid value expression {0}")]
    InvalidValueExpr(String),
    #[error("Incomplete input")]
    Incomplete,
}

impl ParseError {
    pub(crate) fn value_from_nom(err: Err<NomError<&str>>) -> Self {
        match err {
            Err::Error(e) | Err::Failure(e) => ParseError::InvalidValueExpr(e.input.to_string()),
            Err::Incomplete(_) => ParseError::Incomplete,
        }
    }
}
