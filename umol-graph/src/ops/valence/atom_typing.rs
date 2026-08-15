//! Atom-typing valence resolver: admits registry rows for each atom through
//! its constraints view; singleton admissions become edits, plural admissions
//! become completions.

use smallvec::SmallVec;
use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AromaticValenceForm, AsLit, AtomConstraintForm, AtomConstraintKey, AtomForm, AtomId, Lattice,
    Molecule, NumForm,
};
use umol_utils::solution::Solution;

use super::{AtomCompletions, AtomTypeRegistry};

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

    /// Admission: determine candidate sets for each atom under resolution,
    /// Underdetermined if any atom is non-literal, empty if no atoms are admitted.
    pub fn admit(&self, molecule: &Molecule) -> Solution<AtomCompletions, AtomTypingError> {
        for atom in molecule.atoms().iter() {
            if atom.element().as_lit().is_none() {
                return Solution::Underdetermined(AtomCompletions::new());
            }
        }

        let mut completions = AtomCompletions::new();
        for id in molecule.atoms().ids() {
            match self.admitted_completions(molecule, id) {
                Ok(Some(admitted)) => completions.insert(id, admitted),
                Ok(None) => {}
                Err(contradiction) => return Solution::Contradictory(contradiction),
            }
        }
        Solution::Determined(completions)
    }

    /// Determine admitted completions, compatible with the atom's form and constraints.
    /// `None` when the atom is ground or its element is not literal.
    fn admitted_completions(
        &self,
        molecule: &Molecule,
        id: AtomId,
    ) -> Result<Option<SmallVec<[AtomForm; 1]>>, AtomTypingError> {
        let atom = molecule.atom(id);
        let ground_contribution_open = atom.is_ground()
            && matches!(
                atom.constraints()
                    .asserted_complete(AtomConstraintKey::AromaticValence),
                Some(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(NumForm::Undetermined)
                ))
            )
            && !atom.is_in_aromatic_system();
        if atom.is_ground() && !ground_contribution_open {
            return Ok(None);
        }
        let Some(element) = atom.element().as_lit() else {
            return Ok(None);
        };
        let charge = atom.charge().as_lit().map(|n| n as i8);
        let constraints = atom.constraints();
        let admitted: SmallVec<[AtomForm; 1]> = self
            .registry
            .lookup(element, charge)
            .iter()
            .filter(|row| {
                // A field-ground atom is under resolution only for its open
                // aromatic contribution: only aromatic rows complete it.
                if ground_contribution_open
                    && !row
                        .constraints
                        .aromatic_valence()
                        .is_some_and(|a| a.is_aromatic())
                {
                    return false;
                }
                atom.attributes.is_compatible(row)
                    && row.constraints.iter().all(|entry| {
                        let key = entry.key();
                        let derived = match key {
                            AtomConstraintKey::RingDegree
                            | AtomConstraintKey::RingValence
                            | AtomConstraintKey::RingMembership(_) => None,
                            _ => constraints.derived(key),
                        };
                        // The resolution reading: present evidence as in the
                        // open-world check, but a key with neither side reads
                        // the closed-world assertion — absence is actual
                        // absence, so an unmarked atom admits no aromatic row.
                        let host = match (constraints.asserted(key), derived) {
                            (Some(asserted), Some(derived)) => match asserted.meet(&derived) {
                                Some(host) => host,
                                None => return false,
                            },
                            (Some(asserted), None) => asserted.clone(),
                            (None, Some(derived)) => derived,
                            (None, None) => match constraints.asserted_complete(key) {
                                Some(closed) => closed,
                                None => return true,
                            },
                        };
                        entry.is_compatible(&host)
                    })
            })
            .map(|row| {
                atom.attributes
                    .meet(row)
                    .expect("admission implies the meet exists")
            })
            .collect();
        if admitted.is_empty() {
            return Err(AtomTypingError::NoMatch {
                atom_id: id,
                element,
                charge,
            });
        }
        Ok(Some(admitted))
    }

    /// Classify a molecule atom against the registry: `Determined` if some
    /// pattern admits it, `Contradictory` if none does, and `Underdetermined`
    /// if the atom is not ground.
    ///
    /// Admission here reads the closure (`derived_complete`) per row key —
    /// the conformance reading for a ground atom — composed from the view's
    /// keyed core; ring keys compare against the asserted side only (no ring
    /// context is built).
    pub fn classify_molecule_atom(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
    ) -> Solution<(), AtomTypingMismatch> {
        let atom = molecule.atom(atom_id);
        if !atom.is_ground() {
            return Solution::Underdetermined(());
        }
        let Some(element) = atom.element().as_lit() else {
            return Solution::Underdetermined(());
        };
        let charge = atom.charge().as_lit().map(|n| n as i8);
        let constraints = atom.constraints();
        let admitted = self.registry.lookup(element, charge).iter().any(|row| {
            atom.attributes.is_compatible(row)
                && row.constraints.iter().all(|entry| {
                    let key = entry.key();
                    let derived = match key {
                        AtomConstraintKey::RingDegree
                        | AtomConstraintKey::RingValence
                        | AtomConstraintKey::RingMembership(_) => None,
                        _ => constraints.derived_complete(key),
                    };
                    let host = match (constraints.asserted(key), derived) {
                        (Some(asserted), Some(derived)) => match asserted.meet(&derived) {
                            Some(host) => host,
                            None => return false,
                        },
                        (Some(asserted), None) => asserted.clone(),
                        (None, Some(derived)) => derived,
                        (None, None) => return true,
                    };
                    entry.is_compatible(&host)
                })
        });
        if admitted {
            Solution::Determined(())
        } else {
            Solution::Contradictory(AtomTypingMismatch { element, charge })
        }
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::{fixture, rstest};
    use smallvec::smallvec;
    use umol_graph_ir::ir::MoleculeEntries;
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_ground};

    use super::*;
    use crate::ops::valence::AtomTypeRegistry;

    #[fixture]
    fn atom_type_registry() -> AtomTypeRegistry {
        AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h4")])
    }

    #[fixture]
    fn plural_registry() -> AtomTypeRegistry {
        AtomTypeRegistry::from_atoms([
            atom_dsl!("C#c0#h4"),
            atom_dsl!("N#c0#h0#a1"),
            atom_dsl!("N#c0#h1#a2"),
        ])
    }

    #[rstest]
    fn test_atom_typing_valence_admit(atom_type_registry: AtomTypeRegistry) {
        let resolver = AtomTypingValence::new(&atom_type_registry);
        let molecule = mol_dsl!(r#"{:atoms ["C#c0#D1"]}"#);
        let mut expected = AtomCompletions::new();
        expected.insert(AtomId(0), smallvec![atom_dsl!("C#c0#h4#D1")]);
        assert_eq!(resolver.admit(&molecule), Solution::Determined(expected));
    }

    #[rstest]
    fn test_atom_typing_valence_admit_unmarked(plural_registry: AtomTypeRegistry) {
        // Closed-world: an unmarked atom admits no aromatic row; with only
        // aromatic rows for its element, admission is a contradiction.
        let resolver = AtomTypingValence::new(&plural_registry);
        let molecule = mol_dsl!(r#"{:atoms ["N#c0"]}"#);
        assert_eq!(
            resolver.admit(&molecule),
            Solution::Contradictory(AtomTypingError::NoMatch {
                atom_id: AtomId(0),
                element: Element::N,
                charge: Some(0),
            })
        );
    }

    #[rstest]
    fn test_atom_typing_valence_admit_ground_evidenced(plural_registry: AtomTypeRegistry) {
        let resolver = AtomTypingValence::new(&plural_registry);
        let molecule = mol_dsl_ground!(r#"{:atoms ["N#h1" "N#h1"] :bonds [[0 1 "1#a"]]}"#);
        let Solution::Determined(completions) = resolver.admit(&molecule) else {
            panic!("ground-evidenced admission did not determine");
        };
        let admitted: Vec<String> = completions
            .iter()
            .flat_map(|(_, forms)| forms.iter().map(ToString::to_string))
            .collect();
        assert_eq!(
            admitted,
            vec!["N#i=#c0#h#n0#u0#s#a2", "N#i=#c0#h#n0#u0#s#a2"]
        );
    }

    #[rstest]
    fn test_atom_typing_valence_admit_plural(plural_registry: AtomTypeRegistry) {
        let resolver = AtomTypingValence::new(&plural_registry);
        let molecule = mol_dsl!(r#"{:atoms ["C#c0" "N#c0#a+"]}"#);
        let mut expected = AtomCompletions::new();
        expected.insert(AtomId(0), smallvec![atom_dsl!("C#c0#h4")]);
        expected.insert(
            AtomId(1),
            smallvec![atom_dsl!("N#c0#h0#a1"), atom_dsl!("N#c0#h1#a2")],
        );
        assert_eq!(resolver.admit(&molecule), Solution::Determined(expected));
    }

    #[rstest]
    #[case::topology_key_admits(
        AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h3#D1")]),
        Solution::Determined({
            let mut expected = AtomCompletions::new();
            expected.insert(AtomId(0), smallvec![atom_dsl!("C#c0#h3#D1")]);
            expected
        })
    )]
    #[case::topology_key_rejects(
        AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h2#D2")]),
        Solution::Contradictory(AtomTypingError::NoMatch {
            atom_id: AtomId(0),
            element: Element::C,
            charge: Some(0),
        })
    )]
    fn test_atom_typing_valence_admit_admission(
        #[case] registry: AtomTypeRegistry,
        #[case] expected: Solution<AtomCompletions, AtomTypingError>,
    ) {
        let molecule = mol_dsl!(r#"{:atoms ["C#c0" "C#i=#c0#h3#n0#u0#s"] :bonds [[0 1 "1"]]}"#);
        assert_eq!(AtomTypingValence::new(&registry).admit(&molecule), expected);
    }

    #[rstest]
    #[case::default_registry(Cow::Borrowed(AtomTypeRegistry::default_registry()))]
    #[case::empty_registry(Cow::Owned(AtomTypeRegistry::new()))]
    fn test_atom_typing_valence_admit_identity(#[case] registry: Cow<'static, AtomTypeRegistry>) {
        let molecule = mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#);
        assert_eq!(
            AtomTypingValence::new(registry.as_ref()).admit(&molecule),
            Solution::Determined(AtomCompletions::new())
        );
    }

    #[rstest]
    #[case::later_undetermined_element(mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#))]
    fn test_atom_typing_valence_admit_partial(
        atom_type_registry: AtomTypeRegistry,
        #[case] molecule: Molecule,
    ) {
        assert_eq!(
            AtomTypingValence::new(&atom_type_registry).admit(&molecule),
            Solution::Underdetermined(AtomCompletions::new())
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
    #[case::later_atom_contradiction(
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0#h3"]}"#),
        AtomTypingError::NoMatch {
            atom_id: AtomId(1),
            element: Element::C,
            charge: Some(0),
        }
    )]
    fn test_atom_typing_valence_admit_error(
        atom_type_registry: AtomTypeRegistry,
        #[case] molecule: Molecule,
        #[case] expected: AtomTypingError,
    ) {
        assert_eq!(
            AtomTypingValence::new(&atom_type_registry).admit(&molecule),
            Solution::Contradictory(expected)
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
    fn test_atom_typing_valence_classify_molecule_atom(
        atom_type_registry: AtomTypeRegistry,
        #[case] input: &str,
        #[case] expected: Solution<(), AtomTypingMismatch>,
    ) {
        let resolver = AtomTypingValence::new(&atom_type_registry);
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom_dsl!(input)],
            ..Default::default()
        });
        assert_eq!(
            resolver.classify_molecule_atom(&molecule, AtomId(0)),
            expected
        );
    }
}
