//! Atom-typing valence resolver: narrows each atom against a registry of
//! `AtomAst` patterns keyed by element and (optionally) charge.

use thiserror::Error;
use umol_ast::ast::{AsLit, AtomAst, AtomId, Lattice, MoleculeAst};
use umol_chem::element::Element;
use umol_utils::solution::Solution;

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

/// A ground atom that no registry pattern admits.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
#[error("no atom-typing pattern: element {element}, charge {charge:?}")]
pub struct AtomTypingMismatch {
    pub element: Element,
    pub charge: Option<i8>,
}

impl<'a> AtomTypingValence<'a> {
    pub fn new(model: &'a AtomTypingModel) -> Self {
        Self { model }
    }

    /// Iterates atoms, narrowing each non-ground atom against the registry.
    /// Returns `Err` on the first atom that has zero matching patterns.
    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), AtomTypingError> {
        for id in ast.atoms().ids() {
            self.resolve_molecule_atom(ast, id)?;
        }
        Ok(())
    }

    /// Narrows one non-ground atom against the registry. The topology-derived
    /// predicates are set on the atom in place (resolution enriches it); matching
    /// then runs against the atom directly, and the winning registry pattern is
    /// merged via `narrow_from`. Multiple matches resolve via
    /// [`compare_valence_preference`].
    fn resolve_molecule_atom(
        &self,
        ast: &mut MoleculeAst,
        id: AtomId,
    ) -> Result<(), AtomTypingError> {
        let (element, charge, derived) = {
            let atom = ast.atom(id);
            if atom.is_ground() {
                return Ok(());
            }
            let Some(element) = atom.element().as_lit() else {
                return Ok(());
            };
            let charge = atom.charge().as_lit().map(|n| n as i8);
            (element, charge, atom.derive_constraints(false))
        };
        ast.atom_mut(id).ast.constraints.extend(derived);

        let candidates: Vec<&AtomAst> = {
            let atom = ast.atom(id).ast;
            self.model
                .registry
                .lookup(element, charge)
                .iter()
                .filter(|entry| atom.is_compatible(entry))
                .collect()
        };
        let best = match candidates.len() {
            0 => {
                return Err(AtomTypingError::NoMatch {
                    atom_id: id,
                    element,
                    charge,
                });
            }
            1 => candidates[0],
            _ => candidates
                .into_iter()
                .max_by(|a, b| compare_valence_preference(a, b))
                .unwrap(),
        };
        ast.atom_mut(id).ast.narrow_from(best);
        Ok(())
    }

    /// Read-only conformance check for a resolved atom: `Determined` if some
    /// registry pattern admits it, `Contradictory` if none does,
    /// `Underdetermined` if the atom is not ground.
    pub fn conforms_molecule_atom(
        &self,
        ast: &MoleculeAst,
        atom_id: AtomId,
    ) -> Solution<(), AtomTypingMismatch> {
        let atom = ast.atom(atom_id);
        if !atom.is_ground() {
            return Solution::Underdetermined(());
        }
        let Some(element) = atom.element().as_lit() else {
            return Solution::Underdetermined(());
        };
        let constraints = atom.derive_constraints(true);
        let pattern = atom.ast.clone().with_constraints(constraints);
        let charge = pattern.charge.as_lit().map(|n| n as i8);
        let conforms = self
            .model
            .registry
            .lookup(element, charge)
            .iter()
            .any(|entry| pattern.is_compatible(entry));
        if conforms {
            Solution::Determined(())
        } else {
            Solution::Contradictory(AtomTypingMismatch { element, charge })
        }
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::*;
    use umol_ast::ast::Constraints;
    use umol_ast::{atom, mol, mol_ground};

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

    #[rstest]
    #[case::carbon_conforms("C#i=#c0#h4#n0#u0#s#v0#a!", Solution::Determined(()))]
    #[case::wrong_carbon(
        "C#i=#c0#h3#n0#u#s2#v0#a!",
        Solution::Contradictory(AtomTypingMismatch {
            element: Element::C,
            charge: Some(0),
        })
    )]
    #[case::not_ground("C", Solution::Underdetermined(()))]
    fn test_atom_typing_valence_conforms_atom(
        #[case] input: &str,
        #[case] expected: Solution<(), AtomTypingMismatch>,
    ) {
        let model = AtomTypingModel {
            registry: Cow::Owned(registry!["C#c0#h4#n0#u0#s#v0#a!"]),
        };
        let resolver = AtomTypingValence::new(&model);
        let molecule = MoleculeAst::from_parts(
            vec![atom!(input)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        );
        assert_eq!(resolver.conforms_molecule_atom(&molecule, AtomId(0)), expected);
    }
}
