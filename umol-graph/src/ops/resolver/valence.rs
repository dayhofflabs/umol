//! Valence resolver. Dispatches between atom-typing and counts strategies
//! defined in [`crate::ops::valence`].

use thiserror::Error;
use umol_ast::ast::MoleculeAst;

use crate::ops::model::ValenceModel;
use crate::ops::solution::Solution;
use crate::ops::valence::{AtomTypingError, AtomTypingValence, CountsError, CountsValence};

#[derive(Clone, Debug)]
pub enum ValenceResolver {
    AtomTyping(AtomTypingValence),
    Counts(CountsValence),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceContradiction {
    #[error(transparent)]
    AtomTyping(#[from] AtomTypingError),
    #[error(transparent)]
    Counts(#[from] CountsError),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceError {}

impl ValenceResolver {
    pub fn new(model: &ValenceModel) -> Self {
        match model {
            ValenceModel::AtomTyping { registry } => {
                Self::AtomTyping(AtomTypingValence::new(registry.clone()))
            }
            ValenceModel::Counts { table } => Self::Counts(CountsValence::new(table.clone())),
        }
    }

    /// Reports completion or contradiction; the global ground-status verdict
    /// is the composite `Resolver`'s job, not this sub-resolver's.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), ValenceContradiction>, ValenceError> {
        let outcome = match self {
            Self::AtomTyping(r) => r.resolve(ast).map_err(ValenceContradiction::from),
            Self::Counts(r) => r.resolve(ast).map_err(ValenceContradiction::from),
        };
        match outcome {
            Ok(()) => Ok(Solution::Determined(())),
            Err(c) => Ok(Solution::Contradictory(c)),
        }
    }
}
