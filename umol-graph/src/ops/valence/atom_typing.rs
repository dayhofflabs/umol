//! Atom-typing valence resolver: admits registry rows for each atom through
//! its constraints view; singleton admissions become edits, plural admissions
//! become completions.

use smallvec::SmallVec;
use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_ir::ir::{AsLit, AtomConstraintKey, AtomForm, AtomId, Lattice, Molecule};
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

    /// Admission: every atom under resolution gets its candidate set —
    /// one disjunct per admitted registry row — and no edits are produced.
    /// A non-literal element makes the whole admission underdetermined and
    /// empty; plurality is state, not a verdict.
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

    /// The admitted completions of an atom: one disjunct per registry row
    /// surviving admission — field compatibility plus the constraints view's
    /// `is_compatible` — each the meet of the atom's form with the row.
    /// `None` when the atom is ground or its element is not literal.
    fn admitted_completions(
        &self,
        molecule: &Molecule,
        id: AtomId,
    ) -> Result<Option<SmallVec<[AtomForm; 1]>>, AtomTypingError> {
        let atom = molecule.atom(id);
        if atom.is_ground() {
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
                atom.attributes.is_compatible(row) && constraints.is_compatible(&row.constraints)
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
    fn test_atom_typing_valence_admit_plural(plural_registry: AtomTypeRegistry) {
        let resolver = AtomTypingValence::new(&plural_registry);
        let molecule = mol_dsl!(r#"{:atoms ["C#c0" "N#c0"]}"#);
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
