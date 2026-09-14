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
use crate::ops::model::ValenceTieBreak;

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
    /// Exactly equal post-meet forms occur once, in first-occurrence order.
    /// Registry patterns leave isotope information unchanged, including unresolved forms.
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
        let mut admitted = SmallVec::<[AtomForm; 1]>::new();
        for candidate in self
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
        {
            if !admitted.contains(&candidate) {
                admitted.push(candidate);
            }
        }
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

    /// Infer an implicit-H count under the valence policy, without using stored H or lone pairs.
    ///
    /// Filters element/charge lookup results by localized valence, aromatic contribution,
    /// and retained unpaired-electron count/multiplicity. Other atom constraints do not enter
    /// this calculation. Stored valence and aromatic assertions take precedence; absent values
    /// are derived from bonds and aromatic systems. Returns None for an absent atom, missing
    /// concrete lookup evidence, no admissible row, a non-literal admitted H count, or distinct
    /// H counts under Strict.
    /// MostSaturated returns the greatest admitted H count.
    ///
    /// # Semantic properties
    ///
    /// Does not modify the molecule or registry. Reordering or repeating rows does not change
    /// the result. Strict requires agreement on H, not on complete atom forms. The caller
    /// decides whether the inferred count permits eliding the stored count.
    #[cfg_attr(not(test), expect(dead_code))]
    pub(crate) fn infer_implicit_hydrogens(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
        policy: ValenceTieBreak,
    ) -> Option<i64> {
        let atom = molecule.atoms().get(atom_id)?;
        let element = atom.element().as_lit()?;
        let charge = i8::try_from(atom.charge().as_lit()?).ok()?;
        let constraints = atom.constraints();
        let valence = NumForm::Lit(match constraints.valence() {
            Some(asserted) => asserted.as_lit()?,
            None => atom.valence().as_lit()?,
        });
        let aromatic = match constraints.aromatic_valence() {
            Some(asserted) => AromaticValenceForm::from(asserted.as_lit()?),
            None if atom.is_in_aromatic_system() => {
                AromaticValenceForm::aromatic(atom.aromatic_valence().as_lit()?)
            }
            None => AromaticValenceForm::NotAromatic,
        };
        let mut hydrogens = None;
        for row in self.registry.lookup(element, Some(charge)) {
            if !atom
                .unpaired_electrons()
                .is_compatible(&row.unpaired_electrons)
                || row
                    .constraints
                    .valence()
                    .is_some_and(|v| !v.is_compatible(&valence))
                || row
                    .constraints
                    .aromatic_valence()
                    .is_some_and(|a| !a.is_compatible(&aromatic))
            {
                continue;
            }
            let count = row.implicit_hydrogens.as_lit()?;
            hydrogens = Some(match (hydrogens, policy) {
                (None, _) => count,
                (Some(previous), ValenceTieBreak::Strict) if previous != count => return None,
                (Some(previous), _) => count.max(previous),
            });
        }
        hydrogens
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::{fixture, rstest};
    use smallvec::smallvec;
    use umol_graph_ir::ir::{IsotopeMassForm, MoleculeEntries};
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_concrete};

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
    #[case::overlapping(vec![atom_dsl!("C#c0#h*#n0#u0#s"), atom_dsl!("C#c0#h4#n0#u0#s")], smallvec![atom_dsl!("C#i=#c0#h4#n0#u0#s")])]
    #[case::reversed(vec![atom_dsl!("C#c0#h4#n0#u0#s"), atom_dsl!("C#c0#h*#n0#u0#s")], smallvec![atom_dsl!("C#i=#c0#h4#n0#u0#s")])]
    #[case::distinct_constraints(vec![atom_dsl!("C#i*#c0#h4#n0#u0#s"), atom_dsl!("C#i*#c0#h4#n0#u0#s#v0")], smallvec![atom_dsl!("C#i=#c0#h4#n0#u0#s"), atom_dsl!("C#i=#c0#h4#n0#u0#s#v0")])]
    fn test_atom_typing_valence_admit_overlap(
        #[case] rows: Vec<AtomForm>,
        #[case] expected: SmallVec<[AtomForm; 1]>,
    ) {
        let registry = AtomTypeRegistry::from_atoms(rows);
        let molecule = mol_dsl!(r#"{:atoms ["C#i=#c0#h4"]}"#);
        assert_eq!(
            AtomTypingValence::new(&registry).admit(&molecule),
            Solution::Determined(AtomCompletions::from_iter([(AtomId(0), expected)]))
        );
    }

    #[rstest]
    #[case::undetermined(IsotopeMassForm::Undetermined)]
    #[case::natural(IsotopeMassForm::Natural)]
    #[case::mass(IsotopeMassForm::Lit(13))]
    #[case::set(IsotopeMassForm::lit_set([12, 13]))]
    #[case::variable(IsotopeMassForm::var("mass"))]
    fn test_atom_typing_valence_admit_isotope(#[case] isotope: IsotopeMassForm) {
        let source = AtomForm {
            isotope_mass: isotope.clone(),
            ..atom_dsl!("C#c0#h4")
        };
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![source],
            ..Default::default()
        });
        let expected = AtomForm {
            isotope_mass: isotope,
            ..atom_dsl!("C#c0#h4#n0#u0#s#v0#d0#t0#a!#m!")
        };
        assert_eq!(
            AtomTypingValence::new(AtomTypeRegistry::default_registry()).admit(&molecule),
            Solution::Determined(AtomCompletions::from_iter([(
                AtomId(0),
                smallvec![expected]
            )])),
        );
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
        let molecule = mol_dsl_concrete!(r#"{:atoms ["N#h1" "N#h1"] :bonds [[0 1 "1#a"]]}"#);
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
        let molecule = mol_dsl_concrete!(r#"{:atoms ["C #h4"] :bonds []}"#);
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
    #[case::singleton(vec!["C#c0#h4#n0#u0#s#v0#a!"], r#"{:atoms ["C#c0#h4#n0#u0#s"]}"#, Some(4), Some(4))]
    #[case::duplicate(vec!["C#c0#h4#v0#a!", "C#c0#h4#v0#a!"], r#"{:atoms ["C#c0#h4#n0#u0#s"]}"#, Some(4), Some(4))]
    #[case::same_h(vec!["C#c0#h2#n0", "C#c0#h2#n1"], r#"{:atoms ["C#c0#h2#n1#u0#s"]}"#, Some(2), Some(2))]
    #[case::plural(vec!["C#c0#h2#n1", "C#c0#h4#n0"], r#"{:atoms ["C#c0#h2#n1#u0#s"]}"#, None, Some(4))]
    #[case::reversed(vec!["C#c0#h4#n0", "C#c0#h2#n1"], r#"{:atoms ["C#c0#h2#n1#u0#s"]}"#, None, Some(4))]
    #[case::different_stored_h(vec!["C#c0#h4"], r#"{:atoms ["C#c0#h2#n1#u0#s"]}"#, Some(4), Some(4))]
    #[case::zero(vec!["C#c0#h0"], r#"{:atoms ["C#c0#h0#n2#u0#s"]}"#, Some(0), Some(0))]
    #[case::element(vec!["N#c0#h3", "C#c0#h4"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, Some(4), Some(4))]
    #[case::charge(vec!["C#c+#h3", "C#c0#h4"], r#"{:atoms ["C#c+#h3#u0#s"]}"#, Some(3), Some(3))]
    #[case::unpaired(vec!["C#c0#h3#u1#s2", "C#c0#h4#u0#s"], r#"{:atoms ["C#c0#h3#u1#s2"]}"#, Some(3), Some(3))]
    #[case::multiplicity(vec!["C#c0#h2#u2#s3", "C#c0#h0#u2#s"], r#"{:atoms ["C#c0#h2#u2#s3"]}"#, Some(2), Some(2))]
    #[case::valence(vec!["C#c0#h4#v0", "C#c0#h3#v1"], r#"{:atoms ["C#c0#h3#u0#s" "C"] :bonds [[0 1 "1"]]}"#, Some(3), Some(3))]
    #[case::aromatic(vec!["C#c0#h0#v0#a0", "C#c0#h1#v0#a1", "C#c0#h2#v0#a2"], r#"{:atoms ["C#c0#h1#u0#s#a1"]}"#, Some(1), Some(1))]
    #[case::both_valences(vec!["C#c0#h4#v0#a1", "C#c0#h2#v1#a2", "C#c0#h1#v1#a1"], r#"{:atoms ["C#c0#h1#u0#s#a1" "C"] :bonds [[0 1 "1"]]}"#, Some(1), Some(1))]
    #[case::not_aromatic(vec!["C#c0#h1#a1", "C#c0#h4#a!"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, Some(4), Some(4))]
    #[case::other_constraints(vec!["C#c0#h4#x2#y3#V8#H6#D7", "C#c0#h4#x0#y0#V4#H4#D4"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, Some(4), Some(4))]
    #[case::empty(vec![], r#"{:atoms ["C#c0#h4#u0#s"]}"#, None, None)]
    #[case::no_match(vec!["C#c+#h3"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, None, None)]
    #[case::nonliteral_h(vec!["C#c0#h*", "C#c0#h4"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, None, None)]
    #[case::filtered_nonliteral_h(vec!["C#c0#h*#v1", "C#c0#h4#v0"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, Some(4), Some(4))]
    #[case::open_charge(vec!["C#c0#h4"], r#"{:atoms ["C#h4#u0#s"]}"#, None, None)]
    #[case::open_aromatic(vec!["C#c0#h1#a1"], r#"{:atoms ["C#c0#h1#u0#s#a*"]}"#, None, None)]
    #[case::stored_valence(vec!["C#c0#h4#v0", "C#c0#h3#v1"], r#"{:atoms ["C#c0#h4#u0#s#v1"]}"#, Some(3), Some(3))]
    #[case::open_valence(vec!["C#c0#h4#v0"], r#"{:atoms ["C#c0#h4#u0#s#v*"]}"#, None, None)]
    #[case::aromatic_system(
        vec!["C#c0#h1#v2#a1", "C#c0#h0#v2#a2"],
        r#"{:atoms ["C#c0#h1#u0#s" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#,
        Some(1), Some(1)
    )]
    #[case::stored_aromatic(
        vec!["C#c0#h1#v2#a1", "C#c0#h0#v2#a2"],
        r#"{:atoms ["C#c0#h1#u0#s#a2" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#,
        Some(0), Some(0)
    )]
    #[case::aromatic_zero(vec!["C#c0#h2#a0", "C#c0#h4#a!"], r#"{:atoms ["C#c0#h2#u0#s#a0"]}"#, Some(2), Some(2))]
    fn test_atom_typing_valence_infer_implicit_hydrogens(
        #[case] rows: Vec<&str>,
        #[case] input: &str,
        #[case] strict: Option<i64>,
        #[case] most_saturated: Option<i64>,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
    ) {
        let registry = AtomTypeRegistry::from_atoms(rows.into_iter().map(|row| atom_dsl!(row)));
        let molecule = mol_dsl!(input);
        let original_molecule = molecule.clone();
        let original_registry = registry.clone();
        let expected = match policy {
            ValenceTieBreak::Strict => strict,
            ValenceTieBreak::MostSaturated => most_saturated,
        };
        assert_eq!(
            AtomTypingValence::new(&registry).infer_implicit_hydrogens(
                &molecule,
                AtomId(0),
                policy
            ),
            expected,
        );
        assert_eq!(molecule, original_molecule);
        assert_eq!(registry, original_registry);
        assert_eq!(registry.content_hash(), original_registry.content_hash());
    }

    #[rstest]
    #[case::methane(r#"{:atoms ["C#c0#h4#n0#u0#s"]}"#, None, Some(4))]
    #[case::singlet_methylene(r#"{:atoms ["C#c0#h2#n1#u0#s"]}"#, None, Some(4))]
    #[case::methyl(r#"{:atoms ["C#c0#h3#n0#u1#s2"]}"#, None, Some(3))]
    #[case::water(r#"{:atoms ["O#c0#h2#n2#u0#s"]}"#, Some(2), Some(2))]
    #[case::aromatic_carbon(
        r#"{:atoms ["C#c0#h1#n0#u0#s#a1" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"]]}"#,
        Some(1),
        Some(1)
    )]
    #[case::pyrrole_nitrogen(
        r#"{:atoms ["N#c0#h1#n0#u0#s#a2" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"]]}"#,
        Some(1),
        Some(1)
    )]
    fn test_atom_typing_valence_infer_implicit_hydrogens_default_registry(
        #[case] input: &str,
        #[case] strict: Option<i64>,
        #[case] most_saturated: Option<i64>,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
    ) {
        let molecule = mol_dsl!(input);
        let expected = match policy {
            ValenceTieBreak::Strict => strict,
            ValenceTieBreak::MostSaturated => most_saturated,
        };
        assert_eq!(
            AtomTypingValence::new(AtomTypeRegistry::default_registry()).infer_implicit_hydrogens(
                &molecule,
                AtomId(0),
                policy
            ),
            expected,
        );
    }

    #[rstest]
    #[case::carbon_conforms("C#i=#c0#h4#n0#u0#s#v0#a!", Solution::Determined(()))]
    #[case::carbon_isotope("C#i13#c0#h4#n0#u0#s#v0#a!", Solution::Determined(()))]
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
