//! Valence resolver and its supporting data.

pub mod atom_typing;
pub mod counts;
pub mod registry;
mod shared;
pub mod table;

pub use atom_typing::{AtomTypingError, AtomTypingValenceResolver};
pub use counts::{CountsError, CountsValenceResolver};
pub use registry::AtomTypeRegistry;
pub use table::{ValenceEntry, ValenceTable};

use thiserror::Error;
use umol_ast::ast::MoleculeAst;

use crate::ops::config::ValenceModel;
use crate::ops::solution::Solution;

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
            ValenceModel::AtomTyping { registry } => {
                Self::AtomTyping(AtomTypingValenceResolver::new(registry.clone()))
            }
            ValenceModel::Counts {
                table,
                allow_implicit_hydrogens,
            } => Self::Counts(CountsValenceResolver::new(
                table.clone(),
                *allow_implicit_hydrogens,
            )),
        }
    }

    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), ValenceContradiction>, ValenceError> {
        let outcome = match self {
            Self::AtomTyping(r) => r.resolve(ast).map_err(ValenceContradiction::from),
            Self::Counts(r) => r.resolve(ast).map_err(ValenceContradiction::from),
        };

        match outcome {
            Ok(()) => {
                let all_ground = ast.atoms().iter().all(|v| v.data.is_ground());
                Ok(if all_ground {
                    Solution::Determined(())
                } else {
                    Solution::Underdetermined(())
                })
            }
            Err(c) => Ok(Solution::Contradictory(c)),
        }
    }
}
