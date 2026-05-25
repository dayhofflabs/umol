//! Valence resolver. Dispatches between atom-typing and counts strategies
//! defined in [`crate::ops::valence`].

use thiserror::Error;
use umol_ast::ast::MoleculeAst;

use crate::ops::model::ValenceModel;
use crate::ops::solution::Solution;
use crate::ops::valence::{
    AtomTypingError, AtomTypingValenceResolver, CountsError, CountsValenceResolver,
};

#[derive(Clone, Debug)]
pub enum ValenceResolver {
    AtomTyping(AtomTypingValenceResolver),
    Counts(CountsValenceResolver),
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
            ValenceModel::AtomTyping {
                registry,
                normal_valence,
            } => Self::AtomTyping(AtomTypingValenceResolver::new(
                registry.clone(),
                normal_valence.clone(),
            )),
            ValenceModel::Counts {
                table,
                normal_valence,
                allow_implicit_hydrogens,
            } => Self::Counts(CountsValenceResolver::new(
                table.clone(),
                normal_valence.clone(),
                *allow_implicit_hydrogens,
            )),
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
