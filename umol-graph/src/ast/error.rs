//! AST errors.

use thiserror::Error;

#[derive(Clone, Debug, PartialEq, Error)]
#[error("molecule AST is not fully ground")]
pub struct GroundError;
