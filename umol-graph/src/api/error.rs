//! Errors for umol-graph API

use thiserror::Error;
use crate::ast::error::{GroundError, LoweringError};
use crate::dsl::error::ParseError as DslParseError;

#[derive(Clone, Debug, PartialEq, Error)]
pub enum MoleculeEdnError {
    #[error(transparent)]
    Parse(#[from] DslParseError),
    #[error(transparent)]
    Lowering(#[from] LoweringError),
    #[error(transparent)]
    NotGround(#[from] GroundError),
}

