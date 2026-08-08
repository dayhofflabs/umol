//! Tier-3 valence conformance validator: the read-only twin of `ValenceResolver`.
//! Dispatches on `ValenceModel` and folds each engine's per-atom classification
//! over the atoms, surfacing the first mismatch as `Contradictory`.

use thiserror::Error;
use umol_graph_ir::ir::MoleculeAst;
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
            ValenceModel::AtomTyping { registry } => {
                Self::AtomTyping(AtomTypingValence::new(registry.as_ref()))
            }
            ValenceModel::Counts { table } => Self::Counts(CountsValence::new(table.as_ref())),
        }
    }

    pub fn validate(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<(), ValenceConformanceContradiction>, ValenceConformanceError> {
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
    use umol_graph_ir::{atom_dsl, mol_dsl_ground};

    use super::*;
    use crate::ops::valence::{AtomTypeRegistry, ValenceTable};

    #[rstest]
    fn test_valence_conformance_validator_new() {
        let counts = ValenceModel::Counts {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let atom_typing = ValenceModel::AtomTyping {
            registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h4")])),
        };

        assert!(matches!(
            ValenceConformanceValidator::new(&counts),
            ValenceConformanceValidator::Counts(_)
        ));
        assert!(matches!(
            ValenceConformanceValidator::new(&atom_typing),
            ValenceConformanceValidator::AtomTyping(_)
        ));
    }

    #[rstest]
    #[case::counts(ValenceModel::Counts {
        table: Cow::Borrowed(ValenceTable::default_table()),
    })]
    #[case::atom_typing(ValenceModel::AtomTyping {
        registry: Cow::Borrowed(AtomTypeRegistry::default_registry()),
    })]
    fn test_valence_conformance_validator_validate(#[case] model: ValenceModel) {
        let molecule = mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#);
        let result = ValenceConformanceValidator::new(&model)
            .validate(&molecule)
            .unwrap();
        assert_eq!(result, Solution::Determined(()));
    }
}
