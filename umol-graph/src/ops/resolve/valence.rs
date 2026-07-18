//! Valence resolver. Dispatches between atom-typing and counts strategies
//! defined in [`crate::ops::valence`].

use thiserror::Error;
use umol_ast::ast::{MoleculeAst, TransactionError};
use umol_utils::solution::Solution;

use crate::ops::model::ValenceModel;
use crate::ops::valence::{AtomTypingError, AtomTypingValence, CountsError, CountsValence};

#[derive(Clone, Debug)]
pub enum ValenceResolver<'a> {
    AtomTyping(AtomTypingValence<'a>),
    Counts(CountsValence<'a>),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceContradiction {
    #[error(transparent)]
    AtomTyping(#[from] AtomTypingError),
    #[error(transparent)]
    Counts(#[from] CountsError),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceError {
    #[error(transparent)]
    Transaction(#[from] TransactionError),
}

impl<'a> ValenceResolver<'a> {
    pub fn new(model: &'a ValenceModel) -> Self {
        match model {
            ValenceModel::AtomTyping(m) => Self::AtomTyping(AtomTypingValence::new(m)),
            ValenceModel::Counts(m) => Self::Counts(CountsValence::new(m)),
        }
    }

    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), ValenceContradiction>, ValenceError> {
        let outcome = match self {
            Self::AtomTyping(r) => r.resolve(ast).map_err(ValenceContradiction::from),
            Self::Counts(r) => match r.resolve(ast)? {
                Solution::Determined(()) | Solution::Underdetermined(()) => Ok(()),
                Solution::Contradictory(c) => Err(ValenceContradiction::from(c)),
            },
        };
        match outcome {
            Ok(()) => Ok(Solution::Determined(())),
            Err(c) => Ok(Solution::Contradictory(c)),
        }
    }
}
