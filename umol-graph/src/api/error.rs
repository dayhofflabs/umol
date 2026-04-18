//! Errors for umol-graph API

use thiserror::Error;
use crate::dsl::error::ParseError as DslParseError;
use crate::ast::error::GroundError;

#[derive(Clone, Debug, PartialEq, Error)]
pub enum MoleculeEdnError {
    #[error(transparent)]
    Parse(#[from] DslParseError),
    #[error(transparent)]
    NotGround(#[from] GroundError),
}

