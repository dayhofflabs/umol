//! Valence resolver. Dispatches between atom-typing and counts strategies
//! defined in [`crate::ops::valence`].

use thiserror::Error;
use umol_ast::ast::{Edit, MoleculeAst, TransactionError};
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
            ValenceModel::AtomTyping { registry } => {
                Self::AtomTyping(AtomTypingValence::new(registry.as_ref()))
            }
            ValenceModel::Counts { table } => Self::Counts(CountsValence::new(table.as_ref())),
        }
    }

    /// Construct the selected valence model's complete edit plan.
    ///
    /// A non-literal element makes the whole plan underdetermined and yields
    /// no edits.
    pub fn plan(&self, ast: &MoleculeAst) -> Solution<Vec<Edit>, ValenceContradiction> {
        match self {
            Self::AtomTyping(resolver) => resolver
                .plan(ast)
                .map_contradiction(ValenceContradiction::from),
            Self::Counts(resolver) => resolver
                .plan(ast)
                .map_contradiction(ValenceContradiction::from),
        }
    }

    /// Plan and atomically apply the selected valence model.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), ValenceContradiction>, ValenceError> {
        let edits = match self.plan(ast) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let mut editor = ast.edit();
        editor.transact(edits)?;
        *ast = editor.build();
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::rstest;
    use umol_ast::ast::{AtomFieldChange, AtomHandle, AtomId, IsotopeMassAst};
    use umol_ast::{atom_dsl, mol_dsl};
    use umol_chem::element::Element;

    use super::*;
    use crate::ops::valence::{AtomTypeRegistry, ValenceTable};

    #[rstest]
    fn test_valence_resolver_new() {
        let counts = ValenceModel::Counts {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let atom_typing = ValenceModel::AtomTyping {
            registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h4")])),
        };

        assert!(matches!(
            ValenceResolver::new(&counts),
            ValenceResolver::Counts(_)
        ));
        assert!(matches!(
            ValenceResolver::new(&atom_typing),
            ValenceResolver::AtomTyping(_)
        ));
    }

    #[rstest]
    #[case::counts(ValenceModel::Counts {
        table: Cow::Borrowed(ValenceTable::default_table()),
    })]
    #[case::atom_typing(ValenceModel::AtomTyping {
        registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        )])),
    })]
    fn test_valence_resolver_plan(#[case] model: ValenceModel) {
        let molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s#v0#a!"]}"#);
        assert_eq!(
            ValenceResolver::new(&model).plan(&molecule),
            Solution::Determined(vec![Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::IsotopeMass {
                    old: IsotopeMassAst::Undetermined,
                    new: IsotopeMassAst::Natural,
                },
            }])
        );
    }

    #[rstest]
    #[case::counts(ValenceModel::Counts {
        table: Cow::Borrowed(ValenceTable::default_table()),
    })]
    #[case::atom_typing(ValenceModel::AtomTyping {
        registry: Cow::Borrowed(AtomTypeRegistry::default_registry()),
    })]
    fn test_valence_resolver_plan_partial(#[case] model: ValenceModel) {
        let molecule = mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#);
        assert_eq!(
            ValenceResolver::new(&model).plan(&molecule),
            Solution::Underdetermined(Vec::new())
        );
    }

    #[rstest]
    #[case::counts(
        ValenceModel::Counts {
            table: Cow::Borrowed(ValenceTable::default_table()),
        },
        mol_dsl!(r#"{:atoms ["C#c0" "Fe#c0#h0#a+"]}"#),
        ValenceContradiction::Counts(CountsError::UndeterminedAromaticValence)
    )]
    #[case::atom_typing(
        ValenceModel::AtomTyping {
            registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h4")])),
        },
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0#h3"]}"#),
        ValenceContradiction::AtomTyping(AtomTypingError::NoMatch {
            atom_id: AtomId(1),
            element: Element::C,
            charge: Some(0),
        })
    )]
    fn test_valence_resolver_plan_error(
        #[case] model: ValenceModel,
        #[case] molecule: MoleculeAst,
        #[case] expected: ValenceContradiction,
    ) {
        assert_eq!(
            ValenceResolver::new(&model).plan(&molecule),
            Solution::Contradictory(expected)
        );
    }

    #[rstest]
    #[case::counts(ValenceModel::Counts {
        table: Cow::Borrowed(ValenceTable::default_table()),
    })]
    #[case::atom_typing(ValenceModel::AtomTyping {
        registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        )])),
    })]
    fn test_valence_resolver_resolve(#[case] model: ValenceModel) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s#v0#a!"]}"#);
        assert_eq!(
            ValenceResolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(
            molecule,
            mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#a!"]}"#)
        );
    }

    #[rstest]
    #[case::counts(ValenceModel::Counts {
        table: Cow::Borrowed(ValenceTable::default_table()),
    })]
    #[case::atom_typing(ValenceModel::AtomTyping {
        registry: Cow::Borrowed(AtomTypeRegistry::default_registry()),
    })]
    fn test_valence_resolver_resolve_partial(#[case] model: ValenceModel) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#);
        let original = molecule.clone();
        assert_eq!(
            ValenceResolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::counts(
        ValenceModel::Counts {
            table: Cow::Borrowed(ValenceTable::default_table()),
        },
        mol_dsl!(r#"{:atoms ["C#c0" "Fe#c0#h0#a+"]}"#),
        ValenceContradiction::Counts(CountsError::UndeterminedAromaticValence)
    )]
    #[case::atom_typing(
        ValenceModel::AtomTyping {
            registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h4")])),
        },
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0#h3"]}"#),
        ValenceContradiction::AtomTyping(AtomTypingError::NoMatch {
            atom_id: AtomId(1),
            element: Element::C,
            charge: Some(0),
        })
    )]
    fn test_valence_resolver_resolve_error(
        #[case] model: ValenceModel,
        #[case] mut molecule: MoleculeAst,
        #[case] expected: ValenceContradiction,
    ) {
        let original = molecule.clone();
        assert_eq!(
            ValenceResolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
    }
}
