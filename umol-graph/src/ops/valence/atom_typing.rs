//! Atom-typing valence resolver: narrows each atom against a registry of
//! `AtomAst` patterns keyed by element and (optionally) charge.

use thiserror::Error;
use umol_ast::ast::{AsLit, AtomAst, AtomId, Lattice, MoleculeAst};
use umol_shared::element::Element;

use super::compare::compare_valence_preference;
use crate::ops::model::AtomTypingModel;

#[derive(Clone, Debug)]
pub struct AtomTypingValence<'a> {
    model: &'a AtomTypingModel,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AtomTypingError {
    #[error("no atom-typing match for {atom_id:?} (element {element}, charge {charge:?})")]
    NoMatch {
        atom_id: AtomId,
        element: Element,
        charge: Option<i8>,
    },
}

impl<'a> AtomTypingValence<'a> {
    pub fn new(model: &'a AtomTypingModel) -> Self {
        Self { model }
    }

    /// Iterates atoms, narrowing each non-ground atom against the registry.
    /// Returns `Err` on the first atom that has zero matching patterns.
    /// Multiple matches are resolved via [`compare_valence_preference`].
    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), AtomTypingError> {
        for id in ast.atoms().ids() {
            let atom = ast.atom(id);
            if atom.is_ground() {
                continue;
            }

            let Some(element) = atom.element().as_lit() else {
                continue;
            };

            // Topological derived predicates should return Lit.
            let constraints = atom.derive_constraints();
            let pattern = atom.ast.clone().with_constraints(constraints);
            let charge = pattern.charge.as_lit().map(|n| n as i8);
            let candidates: Vec<&AtomAst> = self
                .model
                .registry
                .lookup(element, charge)
                .iter()
                .filter(|entry| pattern.matches(entry))
                .collect();

            match candidates.len() {
                0 => {
                    return Err(AtomTypingError::NoMatch {
                        atom_id: id,
                        element,
                        charge,
                    });
                }
                1 => {
                    let cand = candidates[0];
                    let atom_mut = ast.atom_mut(id).ast;
                    atom_mut.narrow_from(cand);
                }
                _ => {
                    let best = candidates
                        .into_iter()
                        .max_by(|a, b| compare_valence_preference(a, b))
                        .unwrap();
                    let atom_mut = ast.atom_mut(id).ast;
                    atom_mut.narrow_from(best);
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::*;
    use umol_ast::{mol, mol_ground};

    use super::*;
    use crate::ops::valence::AtomTypeRegistry;
    use crate::registry;

    #[rstest]
    #[case::default_registry(Cow::Borrowed(AtomTypeRegistry::default_registry()))]
    #[case::empty_registry(Cow::Owned(AtomTypeRegistry::new()))]
    fn test_atom_typing_valence_resolve_identity(#[case] registry: Cow<'static, AtomTypeRegistry>) {
        let molecule = mol_ground!(r#"{:atoms ["C #h4"] :bonds []}"#);
        let model = AtomTypingModel { registry };
        let mut resolved = molecule.clone();
        AtomTypingValence::new(&model)
            .resolve(&mut resolved)
            .unwrap();
        assert_eq!(resolved, molecule);
    }

    #[rstest]
    fn test_atom_typing_valence_resolve_error() {
        let model = AtomTypingModel {
            registry: Cow::Owned(registry!["C#c0#h4#n0#u0"]),
        };
        let mut molecule = mol!(r#"{:atoms ["C #h3" "Cl"] :bonds [[0 1 "1"]]}"#);
        let err = AtomTypingValence::new(&model)
            .resolve(&mut molecule)
            .unwrap_err();
        assert_eq!(
            err,
            AtomTypingError::NoMatch {
                atom_id: AtomId(0),
                element: Element::C,
                charge: None,
            }
        );
    }
}
