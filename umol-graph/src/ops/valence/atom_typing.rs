//! Atom-typing valence resolver: plans registry-driven narrowing for each atom
//! against `AtomForm` patterns keyed by element and (optionally) charge.

use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AsLit, AtomForm, AtomHandle, AtomId, Edits, Lattice, MoleculeAst, TransactionError,
};
use umol_utils::solution::Solution;

use super::compare::compare_valence_preference;
use super::AtomTypeRegistry;

#[derive(Clone, Debug)]
pub struct AtomTypingValence<'a> {
    registry: &'a AtomTypeRegistry,
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
    pub fn new(registry: &'a AtomTypeRegistry) -> Self {
        Self { registry }
    }

    /// Construct the complete edit plan without mutating `ast`.
    ///
    /// A non-literal element makes the whole plan underdetermined and yields
    /// no edits.
    pub fn plan(&self, ast: &MoleculeAst) -> Solution<Edits, AtomTypingError> {
        for atom in ast.atoms().iter() {
            if atom.element().as_lit().is_none() {
                return Solution::Underdetermined(Edits::new());
            }
        }

        let mut edits = Edits::new();
        for id in ast.atoms().ids() {
            let selected = match self.resolve_molecule_atom(ast, id) {
                Ok(Some(selected)) => selected,
                Ok(None) => continue,
                Err(contradiction) => return Solution::Contradictory(contradiction),
            };
            let current = ast.atom(id).ast;
            let update = current.difference_to(&selected);
            edits.update_atom(AtomHandle::Id(id), current, &update);
        }
        Solution::Determined(edits)
    }

    /// Plan and atomically apply atom-typing valence resolution.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), AtomTypingError>, TransactionError> {
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

    /// Compute the selected atom without mutating the molecule.
    fn resolve_molecule_atom(
        &self,
        ast: &MoleculeAst,
        id: AtomId,
    ) -> Result<Option<AtomForm>, AtomTypingError> {
        let atom = ast.atom(id);
        if atom.is_ground() {
            return Ok(None);
        }
        let Some(element) = atom.element().as_lit() else {
            return Ok(None);
        };
        let charge = atom.charge().as_lit().map(|n| n as i8);
        let mut selected = atom.ast.clone();
        selected.constraints.extend(atom.derive_constraints(false));
        let best =
            self.select_candidate(&selected, element, charge)
                .ok_or(AtomTypingError::NoMatch {
                    atom_id: id,
                    element,
                    charge,
                })?;
        selected.narrow_from(best);
        Ok(Some(selected))
    }

    /// Classify a molecule atom against the registry: `Determined` if some
    /// pattern admits it, `Contradictory` if none does, and `Underdetermined`
    /// if the atom is not ground.
    pub fn classify_molecule_atom(
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
        let admitted = self.select_candidate(&pattern, element, charge).is_some();
        if admitted {
            Solution::Determined(())
        } else {
            Solution::Contradictory(AtomTypingMismatch { element, charge })
        }
    }

    fn select_candidate(
        &self,
        pattern: &AtomForm,
        element: Element,
        charge: Option<i8>,
    ) -> Option<&AtomForm> {
        self.registry
            .lookup(element, charge)
            .iter()
            .filter(|entry| pattern.is_compatible(entry))
            .max_by(|a, b| compare_valence_preference(a, b))
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::{fixture, rstest};
    use umol_graph_ir::ir::{
        AtomConstraintForm, AtomFieldChange, Edit, Edits, MoleculeEntries, NumForm,
    };
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_ground};

    use super::*;
    use crate::ops::valence::AtomTypeRegistry;

    #[fixture]
    fn atom_type_registry() -> AtomTypeRegistry {
        AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h4")])
    }

    #[rstest]
    fn test_atom_typing_valence_plan(atom_type_registry: AtomTypeRegistry) {
        let resolver = AtomTypingValence::new(&atom_type_registry);
        let molecule = mol_dsl!(r#"{:atoms ["C#c0#D1"]}"#);
        assert_eq!(
            resolver.plan(&molecule),
            Solution::Determined(Edits::from_iter([
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::ImplicitHydrogens {
                        old: NumForm::Undetermined,
                        new: NumForm::Lit(4),
                    },
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(0)),
                    old: None,
                    new: Some(AtomConstraintForm::valence(0_i64)),
                },
            ]))
        );
    }

    #[rstest]
    #[case::default_registry(Cow::Borrowed(AtomTypeRegistry::default_registry()))]
    #[case::empty_registry(Cow::Owned(AtomTypeRegistry::new()))]
    fn test_atom_typing_valence_plan_identity(#[case] registry: Cow<'static, AtomTypeRegistry>) {
        let molecule = mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#);
        assert_eq!(
            AtomTypingValence::new(registry.as_ref()).plan(&molecule),
            Solution::Determined(Edits::new())
        );
    }

    #[rstest]
    #[case::later_undetermined_element(mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#))]
    fn test_atom_typing_valence_plan_partial(
        atom_type_registry: AtomTypeRegistry,
        #[case] molecule: MoleculeAst,
    ) {
        assert_eq!(
            AtomTypingValence::new(&atom_type_registry).plan(&molecule),
            Solution::Underdetermined(Edits::new())
        );
    }

    #[rstest]
    #[case::no_match(
        mol_dsl!(r#"{:atoms ["C#c0#h3"]}"#),
        AtomTypingError::NoMatch {
            atom_id: AtomId(0),
            element: Element::C,
            charge: Some(0),
        }
    )]
    fn test_atom_typing_valence_plan_error(
        atom_type_registry: AtomTypeRegistry,
        #[case] molecule: MoleculeAst,
        #[case] expected: AtomTypingError,
    ) {
        assert_eq!(
            AtomTypingValence::new(&atom_type_registry).plan(&molecule),
            Solution::Contradictory(expected)
        );
    }

    #[rstest]
    fn test_atom_typing_valence_resolve(atom_type_registry: AtomTypeRegistry) {
        let resolver = AtomTypingValence::new(&atom_type_registry);
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#D1"]}"#);
        assert_eq!(
            resolver.resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#c0#h4#v0#D1"]}"#));
    }

    #[rstest]
    #[case::later_undetermined_element(mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#))]
    fn test_atom_typing_valence_resolve_partial(
        atom_type_registry: AtomTypeRegistry,
        #[case] mut molecule: MoleculeAst,
    ) {
        let original = molecule.clone();
        assert_eq!(
            AtomTypingValence::new(&atom_type_registry).resolve(&mut molecule),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::later_atom_contradiction(
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0#h3"]}"#),
        AtomTypingError::NoMatch {
            atom_id: AtomId(1),
            element: Element::C,
            charge: Some(0),
        }
    )]
    fn test_atom_typing_valence_resolve_error(
        atom_type_registry: AtomTypeRegistry,
        #[case] mut molecule: MoleculeAst,
        #[case] expected: AtomTypingError,
    ) {
        let original = molecule.clone();
        assert_eq!(
            AtomTypingValence::new(&atom_type_registry).resolve(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
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
    fn test_atom_typing_valence_classify_molecule_atom(
        atom_type_registry: AtomTypeRegistry,
        #[case] input: &str,
        #[case] expected: Solution<(), AtomTypingMismatch>,
    ) {
        let resolver = AtomTypingValence::new(&atom_type_registry);
        let molecule = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![atom_dsl!(input)],
            ..Default::default()
        });
        assert_eq!(
            resolver.classify_molecule_atom(&molecule, AtomId(0)),
            expected
        );
    }
}
