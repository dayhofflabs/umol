//! Tier-3 valence conformance validator: the read-only twin of `ValenceResolver`.
//! Dispatches on `ValenceModel` and folds each engine's per-atom classification
//! over the atoms, surfacing the first mismatch as `Contradictory`.

use thiserror::Error;
use umol_ast::ast::MoleculeAst;
use umol_utils::solution::Solution;

use crate::ops::model::ValenceModel;
use crate::ops::valence::{AtomTypingMismatch, AtomTypingValence, CountsMismatch, CountsValence};

#[derive(Clone, Debug)]
pub enum ValenceConformanceValidator<'a> {
    AtomTyping(AtomTypingValence<'a>),
    Counts(CountsValence<'a>),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceConformanceContradiction {
    #[error(transparent)]
    AtomTyping(#[from] AtomTypingMismatch),
    #[error(transparent)]
    Counts(#[from] CountsMismatch),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceConformanceError {}

impl<'a> ValenceConformanceValidator<'a> {
    pub fn new(model: &'a ValenceModel) -> Self {
        match model {
            ValenceModel::AtomTyping(m) => Self::AtomTyping(AtomTypingValence::new(m)),
            ValenceModel::Counts(m) => Self::Counts(CountsValence::new(m)),
        }
    }

    pub fn validate(
        &self,
        ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), ValenceConformanceContradiction>, ValenceConformanceError> {
        let ast = ast.as_ref();
        let mut any_undetermined = false;
        for id in ast.atoms().ids() {
            let outcome = match self {
                Self::AtomTyping(engine) => engine
                    .classify_molecule_atom(ast, id)
                    .map_contradiction(ValenceConformanceContradiction::from),
                Self::Counts(engine) => engine
                    .classify_molecule_atom(ast, id)
                    .map_contradiction(ValenceConformanceContradiction::from),
            };
            match outcome {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
            }
        }
        Ok(if any_undetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::rstest;
    use umol_ast::mol_dsl_ground;

    use super::*;
    use crate::ops::model::{AtomTypingModel, CountsModel};
    use crate::ops::valence::{AtomTypeRegistry, ValenceTable};

    #[rstest]
    #[case::counts(ValenceModel::Counts(CountsModel {
        table: Cow::Borrowed(ValenceTable::default_table()),
    }))]
    #[case::atom_typing(ValenceModel::AtomTyping(AtomTypingModel {
        registry: Cow::Borrowed(AtomTypeRegistry::default_registry()),
    }))]
    fn test_valence_conformance_validator_validate(#[case] model: ValenceModel) {
        let molecule = mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#);
        let result = ValenceConformanceValidator::new(&model)
            .validate(&molecule)
            .unwrap();
        assert_eq!(result, Solution::Determined(()));
    }
}
