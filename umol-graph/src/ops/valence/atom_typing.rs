//! Atom-typing valence resolver: admits registry rows for each atom through
//! its constraints view; singleton admissions become edits, plural admissions
//! become completions.

use smallvec::SmallVec;
use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AromaticValenceForm, AsLit, AtomConstraintForm, AtomConstraintKey, AtomForm, AtomId, AtomView,
    Lattice, Molecule, NumForm, RingConfig, RingModel,
};
use umol_utils::solution::Solution;

use super::{AtomCompletions, AtomTypeRegistry};
use crate::ops::model::ValenceTieBreak;
use crate::ops::resolve::WildcardPolicy;

/// Admits atom completions from a borrowed registry under the selected wildcard policy.
#[derive(Clone, Debug)]
pub struct AtomTypingValence<'a> {
    registry: &'a AtomTypeRegistry,
    wildcards: WildcardPolicy,
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
    /// Uses strict wildcard handling during admission.
    pub fn new(registry: &'a AtomTypeRegistry) -> Self {
        Self::with_config(registry, WildcardPolicy::Strict)
    }

    /// Borrows the registry and sets the wildcard policy for admission.
    pub fn with_config(registry: &'a AtomTypeRegistry, wildcards: WildcardPolicy) -> Self {
        Self {
            registry,
            wildcards,
        }
    }

    /// Determines candidate sets for atoms that need a registry completion.
    /// Admit skips atoms with an undetermined element. Every other atom requires a
    /// literal element. A required completion that depends on an open bond order
    /// yields Underdetermined under Admit. Underdetermined carries an empty set of completions.
    /// Exactly equal post-meet forms occur once, in first-occurrence order.
    /// Registry patterns leave isotope information unchanged, including unresolved forms.
    pub fn admit(&self, molecule: &Molecule) -> Solution<AtomCompletions, AtomTypingError> {
        for atom in molecule.atoms().iter() {
            if !atom.element().is_ground()
                && (self.wildcards == WildcardPolicy::Strict || !atom.attributes().is_wildcard())
            {
                return Solution::Underdetermined(AtomCompletions::new());
            }
        }

        let mut completions = AtomCompletions::new();
        let mut underdetermined = false;
        for atom in molecule.atoms().iter() {
            let Some(element) = atom.element().as_lit() else {
                continue;
            };
            match self.admitted_completions(molecule, atom, element) {
                Solution::Determined(Some(admitted)) => completions.insert(atom.id(), admitted),
                Solution::Determined(None) => {}
                Solution::Underdetermined(_) if self.wildcards == WildcardPolicy::Admit => {
                    underdetermined = true;
                }
                Solution::Underdetermined(_) => {
                    return Solution::Underdetermined(AtomCompletions::new())
                }
                Solution::Contradictory(contradiction) => {
                    return Solution::Contradictory(contradiction)
                }
            }
        }
        if underdetermined {
            Solution::Underdetermined(AtomCompletions::new())
        } else {
            Solution::Determined(completions)
        }
    }

    /// Determine admitted completions, compatible with the atom's form and constraints.
    /// A determined `None` means the atom needs no candidate. Lookup uses its complete valence reading.
    fn admitted_completions(
        &self,
        molecule: &Molecule,
        atom: AtomView<'_>,
        element: Element,
    ) -> Solution<Option<SmallVec<[AtomForm; 1]>>, AtomTypingError> {
        let charge = atom.charge().as_lit().map(|n| n as i8);
        let mismatch = AtomTypingError::NoMatch {
            atom_id: atom.id(),
            element,
            charge,
        };
        let Some(aromatic) = atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::AromaticValence)
        else {
            return Solution::Contradictory(mismatch);
        };
        let is_ground = atom.is_ground();
        let contribution_open = matches!(&aromatic,
            AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(value)) if !value.is_ground()
        );
        let ground_contribution_open = is_ground && contribution_open;
        if is_ground && !ground_contribution_open {
            return Solution::Determined(None);
        }
        let open_order_outcome = || {
            if self.wildcards == WildcardPolicy::Admit
                && !contribution_open
                && atom.charge().is_ground()
                && atom.implicit_hydrogens().is_ground()
                && atom.lone_pairs().is_ground()
                && atom.unpaired_electrons().is_ground()
            {
                Solution::Determined(None)
            } else {
                Solution::Underdetermined(None)
            }
        };
        let Some(AtomConstraintForm::Valence(valence)) = atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::Valence)
        else {
            return Solution::Contradictory(mismatch);
        };
        let Some(valence) = valence.as_lit() else {
            return open_order_outcome();
        };
        let rows = self.registry.lookup(element, valence as u8);
        let rings = rows
            .iter()
            .any(|row| {
                row.constraints.iter().any(|c| {
                    matches!(
                        c.key(),
                        AtomConstraintKey::RingDegree
                            | AtomConstraintKey::RingValence
                            | AtomConstraintKey::RingMembership(_)
                    )
                })
            })
            .then(|| {
                molecule
                    .rings(RingModel::default(), RingConfig::default())
                    .into_ring_set()
            });
        let constraints = match &rings {
            Some(rings) => atom.constraints().with_rings(rings),
            None => atom.constraints(),
        };
        let mut readings = SmallVec::<[(AtomConstraintKey, Option<AtomConstraintForm>); 5]>::new();
        readings.push((AtomConstraintKey::AromaticValence, Some(aromatic)));
        readings.push((
            AtomConstraintKey::Valence,
            Some(AtomConstraintForm::valence(valence)),
        ));
        let attributes_first = !atom.implicit_hydrogens().is_undetermined();
        let check_dative = self.wildcards == WildcardPolicy::Admit && atom.is_in_dative_bond();
        let mut admitted = SmallVec::<[AtomForm; 1]>::new();
        for row in rows {
            if ground_contribution_open && !row.constraints.aromatic_valence().is_aromatic() {
                continue;
            }
            if !atom.charge().is_compatible(&row.charge)
                || (attributes_first && !atom.attributes().is_compatible(row))
            {
                continue;
            }
            let mut needs_dative_order = false;
            let compatible = row.constraints.iter().all(|entry| {
                let key = entry.key();
                let index = match readings.iter().position(|(stored, _)| *stored == key) {
                    Some(index) => index,
                    None => {
                        let index = readings.len();
                        readings.push((key, constraints.evidenced_complete(key)));
                        index
                    }
                };
                let Some(host) = &readings[index].1 else {
                    return false;
                };
                if check_dative
                    && matches!(
                        key,
                        AtomConstraintKey::DonatedPairs | AtomConstraintKey::AcceptedPairs
                    )
                    && !host.is_ground()
                    && !entry.matches(host)
                {
                    needs_dative_order = true;
                }
                entry.is_compatible(host)
            });
            if !compatible || (!attributes_first && !atom.attributes().is_compatible(row)) {
                continue;
            }
            if needs_dative_order {
                return open_order_outcome();
            }
            let candidate = atom
                .attributes()
                .meet(row)
                .expect("admission implies the meet exists");
            if !admitted.contains(&candidate) {
                admitted.push(candidate);
            }
        }
        if admitted.is_empty() {
            return Solution::Contradictory(mismatch);
        }
        Solution::Determined(Some(admitted))
    }

    /// Classify a molecule atom against the registry: `Determined` if some
    /// pattern admits it, `Contradictory` if none does, and `Underdetermined`
    /// if the atom is not ground or its complete valence is non-literal.
    ///
    /// Admission compares each row key with `evidenced_complete`; contradictory
    /// evidence admits no row. Ring context is built when a row requests it.
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
        let Some(AtomConstraintForm::Valence(valence)) = atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::Valence)
        else {
            return Solution::Contradictory(AtomTypingMismatch { element, charge });
        };
        let Some(valence) = valence.as_lit() else {
            return Solution::Underdetermined(());
        };
        let rows = self.registry.lookup(element, valence as u8);
        let rings = rows
            .iter()
            .any(|row| {
                row.constraints.iter().any(|c| {
                    matches!(
                        c.key(),
                        AtomConstraintKey::RingDegree
                            | AtomConstraintKey::RingValence
                            | AtomConstraintKey::RingMembership(_)
                    )
                })
            })
            .then(|| {
                molecule
                    .rings(RingModel::default(), RingConfig::default())
                    .into_ring_set()
            });
        let constraints = match &rings {
            Some(rings) => atom.constraints().with_rings(rings),
            None => atom.constraints(),
        };
        let admitted = rows.iter().any(|row| {
            atom.attributes().is_compatible(row)
                && row.constraints.iter().all(|entry| {
                    constraints
                        .evidenced_complete(entry.key())
                        .is_some_and(|host| entry.is_compatible(&host))
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
    /// Filters element/valence lookup results by charge, aromatic contribution,
    /// and retained unpaired-electron count/multiplicity. Other atom constraints do not enter
    /// this calculation. Valence and aromatic evidence use the complete constraint reading,
    /// meeting assertions with relation evidence. Returns None for an absent atom, missing
    /// concrete lookup evidence, no admissible row, a non-literal admitted H count, or distinct
    /// H counts under Strict.
    /// MostSaturated returns the greatest admitted H count.
    ///
    /// # Semantic properties
    ///
    /// Does not modify the molecule or registry. Reordering or repeating rows does not change
    /// the result. Strict requires agreement on H, not on complete atom forms. The caller
    /// decides whether the inferred count permits eliding the stored count.
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
        let AtomConstraintForm::Valence(valence) =
            constraints.evidenced_complete(AtomConstraintKey::Valence)?
        else {
            unreachable!()
        };
        let valence = valence.as_lit()? as u8;
        let AtomConstraintForm::AromaticValence(aromatic) =
            constraints.evidenced_complete(AtomConstraintKey::AromaticValence)?
        else {
            unreachable!()
        };
        let aromatic = AromaticValenceForm::from(aromatic.as_lit()?);
        let mut hydrogens = None;
        for row in self.registry.lookup(element, valence) {
            if row.charge != NumForm::Lit(i64::from(charge))
                || !atom
                    .unpaired_electrons()
                    .is_compatible(&row.unpaired_electrons)
                || !row.constraints.aromatic_valence().is_compatible(&aromatic)
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
    use std::ptr;

    use rstest::{fixture, rstest};
    use smallvec::smallvec;
    use umol_graph_ir::ir::{
        BondForm, IsotopeMassForm, MoleculeEntries, RingMembershipForm, RingScope, TryIntoIr,
    };
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_concrete};
    use umol_io::smiles::Smiles;

    use super::*;
    use crate::ops::valence::AtomTypeRegistry;

    #[fixture]
    fn atom_type_registry() -> AtomTypeRegistry {
        AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h4#v0")])
    }

    #[fixture]
    fn plural_registry() -> AtomTypeRegistry {
        AtomTypeRegistry::from_atoms([
            atom_dsl!("C#c0#h4#v0"),
            atom_dsl!("N#c0#h0#a1#v0"),
            atom_dsl!("N#c0#h1#a2#v0"),
            atom_dsl!("N#c0#h0#a1#v1"),
            atom_dsl!("N#c0#h1#a2#v1"),
        ])
    }

    #[rstest]
    fn test_atom_typing_valence_new(atom_type_registry: AtomTypeRegistry) {
        let resolver = AtomTypingValence::new(&atom_type_registry);
        assert!(ptr::eq(resolver.registry, &atom_type_registry));
        assert_eq!(resolver.wildcards, WildcardPolicy::Strict);
        assert_eq!(
            resolver.admit(&mol_dsl!(r#"{:atoms ["*"]}"#)),
            Solution::Underdetermined(AtomCompletions::new())
        );
    }

    #[rstest]
    fn test_atom_typing_valence_with_config(
        atom_type_registry: AtomTypeRegistry,
        #[values(WildcardPolicy::Strict, WildcardPolicy::Admit)] wildcards: WildcardPolicy,
    ) {
        let resolver = AtomTypingValence::with_config(&atom_type_registry, wildcards);
        assert!(ptr::eq(resolver.registry, &atom_type_registry));
        assert_eq!(resolver.wildcards, wildcards);
    }

    #[rstest]
    #[case::degree(
        &atom_type_registry(), WildcardPolicy::Strict, mol_dsl!(r#"{:atoms ["C#c0#D1"]}"#),
        AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("C#c0#h4#D1#v0")])]),
    )]
    #[case::branched(
        AtomTypeRegistry::default_registry(), WildcardPolicy::Admit,
        (&Smiles::parse("*C(O)C").unwrap().into_table_ir()).try_into_ir(&()).unwrap(),
        AtomCompletions::from_iter([
            (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h1#n0#u0#s#v3#d0#t0#a!#m!")]),
            (AtomId(2), smallvec![atom_dsl!("O#i=#c0#h1#n2#u0#s#v1#d0#t0#a!#m!")]),
            (AtomId(3), smallvec![
                atom_dsl!("C#i=#c0#h3#n0#u0#s#v1#d0#t0#a!#m!"),
                atom_dsl!("C#i=#c0#h1#n1#u0#s#v1#d0#t0#a!#m!"),
            ]),
        ]),
    )]
    #[case::wildcard(
        &atom_type_registry(), WildcardPolicy::Admit, mol_dsl!(r#"{:atoms ["*"]}"#),
        AtomCompletions::new(),
    )]
    #[case::complete_localized(
        &atom_type_registry(), WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["C#c0#h3#n0#u0#s" "*"] :bonds [[0 1 "*"]]}"#),
        AtomCompletions::new(),
    )]
    #[case::complete_dative(
        AtomTypeRegistry::default_registry(), WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["N#c0#h3#n1#u0#s" "*"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "*"}]}"#),
        AtomCompletions::new(),
    )]
    #[case::unrestricted_dative(
        &AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h2#v0#d1#t1"), atom_dsl!("C#c0#h4#v0")]),
        WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["C#c0" "*"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "*"}]}"#),
        AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("C#c0#h4#v0")])]),
    )]
    fn test_atom_typing_valence_admit(
        #[case] registry: &AtomTypeRegistry,
        #[case] wildcards: WildcardPolicy,
        #[case] molecule: Molecule,
        #[case] expected: AtomCompletions,
    ) {
        let resolver = AtomTypingValence::with_config(registry, wildcards);
        assert_eq!(resolver.admit(&molecule), Solution::Determined(expected));
    }

    #[rstest]
    #[case::overlapping(vec![atom_dsl!("C#c0#h*#n0#u0#s#v0"), atom_dsl!("C#c0#h4#n0#u0#s#v0")], smallvec![atom_dsl!("C#i=#c0#h4#n0#u0#s#v0")])]
    #[case::reversed(vec![atom_dsl!("C#c0#h4#n0#u0#s#v0"), atom_dsl!("C#c0#h*#n0#u0#s#v0")], smallvec![atom_dsl!("C#i=#c0#h4#n0#u0#s#v0")])]
    #[case::distinct_constraints(vec![atom_dsl!("C#i*#c0#h4#n0#u0#s#v0"), atom_dsl!("C#i*#c0#h4#n0#u0#s#v0#D0")], smallvec![atom_dsl!("C#i=#c0#h4#n0#u0#s#v0"), atom_dsl!("C#i=#c0#h4#n0#u0#s#v0#D0")])]
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
            vec!["N#i=#c0#h#n0#u0#s#v#a2", "N#i=#c0#h#n0#u0#s#v#a2"]
        );
    }

    #[rstest]
    fn test_atom_typing_valence_admit_plural(plural_registry: AtomTypeRegistry) {
        let resolver = AtomTypingValence::new(&plural_registry);
        let molecule = mol_dsl!(r#"{:atoms ["C#c0" "N#c0#a+"]}"#);
        let mut expected = AtomCompletions::new();
        expected.insert(AtomId(0), smallvec![atom_dsl!("C#c0#h4#v0")]);
        expected.insert(
            AtomId(1),
            smallvec![atom_dsl!("N#c0#h0#a1#v0"), atom_dsl!("N#c0#h1#a2#v0")],
        );
        assert_eq!(resolver.admit(&molecule), Solution::Determined(expected));
    }

    #[rstest]
    #[case::topology_key_admits(
        AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h3#D1#v1")]),
        Solution::Determined({
            let mut expected = AtomCompletions::new();
            expected.insert(AtomId(0), smallvec![atom_dsl!("C#c0#h3#D1#v1")]);
            expected
        })
    )]
    #[case::topology_key_rejects(
        AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h2#D2#v1")]),
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
    #[case::repeated_key(NumForm::Lit(1), vec!["C#c0#h3#v2", "C#c0#h3#v1", "C#c0#h3#v1"], vec!["C#c0#h3#v1"])]
    #[case::singleton(NumForm::lit_set([1]), vec!["C#c0#h3#v2", "C#c0#h3#v1"], vec!["C#c0#h3#v1"])]
    #[case::bottom(NumForm::lit_set([]), vec!["C#c0#h3#v2", "C#c0#h3#v1"], vec![])]
    #[case::multiple_keys(NumForm::Lit(1), vec!["C#c0#h3#v1#t1", "C#c0#h3#v1#t0"], vec!["C#c0#h3#v1#t0"])]
    fn test_atom_typing_valence_admit_readings(
        #[case] order: NumForm,
        #[case] rows: Vec<&str>,
        #[case] expected: Vec<&str>,
    ) {
        let registry = AtomTypeRegistry::from_atoms(rows.into_iter().map(|row| atom_dsl!(row)));
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom_dsl!("C#c0"), atom_dsl!("C#i=#c0#h3#n0#u0#s")],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::new(order))],
            ..Default::default()
        });
        let expected = if expected.is_empty() {
            Solution::Contradictory(AtomTypingError::NoMatch {
                atom_id: AtomId(0),
                element: Element::C,
                charge: Some(0),
            })
        } else {
            Solution::Determined(AtomCompletions::from_iter([(
                AtomId(0),
                expected.into_iter().map(|row| atom_dsl!(row)).collect(),
            )]))
        };
        assert_eq!(AtomTypingValence::new(&registry).admit(&molecule), expected);
    }

    #[rstest]
    fn test_atom_typing_valence_admit_ring_readings() {
        let row = atom_dsl!("C#c0#h2#n0#u0#s").with_constraints([
            AtomConstraintForm::Valence(NumForm::Lit(2)),
            AtomConstraintForm::DonatedPairs(NumForm::Lit(0)),
            AtomConstraintForm::AcceptedPairs(NumForm::Lit(0)),
            AtomConstraintForm::Degree(NumForm::Lit(2)),
            AtomConstraintForm::TotalDegree(NumForm::Lit(4)),
            AtomConstraintForm::RingMembership(RingMembershipForm::new(RingScope::Size(3), 1)),
            AtomConstraintForm::RingMembership(RingMembershipForm::new(RingScope::Size(4), 0)),
        ]);
        let rejected = row
            .clone()
            .with_constraint(AtomConstraintForm::RingMembership(RingMembershipForm::new(
                RingScope::Size(3),
                0,
            )));
        let registry = AtomTypeRegistry::from_atoms([rejected, row.clone()]);
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#
        );
        assert_eq!(
            AtomTypingValence::new(&registry).admit(&molecule),
            Solution::Determined(AtomCompletions::from_iter([(AtomId(0), smallvec![row])])),
        );
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
    #[case::complement(mol_dsl!(r#"{:atoms ["!{C,N}"]}"#))]
    #[case::variable(mol_dsl!(r#"{:atoms ["?e"]}"#))]
    #[case::bond_order(mol_dsl!(r#"{:atoms ["C#c0" "C#c0"] :bonds [[0 1 "*"]]}"#))]
    fn test_atom_typing_valence_admit_partial(
        atom_type_registry: AtomTypeRegistry,
        #[values(WildcardPolicy::Strict, WildcardPolicy::Admit)] wildcards: WildcardPolicy,
        #[case] molecule: Molecule,
    ) {
        assert_eq!(
            AtomTypingValence::with_config(&atom_type_registry, wildcards).admit(&molecule),
            Solution::Underdetermined(AtomCompletions::new())
        );
    }

    #[rstest]
    #[case::donor(mol_dsl!(r#"{:atoms ["N#c0" "*"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "*"}]}"#))]
    #[case::acceptor(mol_dsl!(r#"{:atoms ["*" "B#c0"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "*"}]}"#))]
    fn test_atom_typing_valence_admit_dative(#[case] molecule: Molecule) {
        let resolver = AtomTypingValence::with_config(
            AtomTypeRegistry::default_registry(),
            WildcardPolicy::Admit,
        );
        assert_eq!(
            resolver.admit(&molecule),
            Solution::Underdetermined(AtomCompletions::new())
        );
    }

    #[rstest]
    #[case::no_match(
        WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["C#c0#h3"]}"#),
        AtomTypingError::NoMatch {
            atom_id: AtomId(0),
            element: Element::C,
            charge: Some(0),
        }
    )]
    #[case::later_atom_contradiction(
        WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0#h3"]}"#),
        AtomTypingError::NoMatch {
            atom_id: AtomId(1),
            element: Element::C,
            charge: Some(0),
        }
    )]
    #[case::ground_aromatic_conflict(
        WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["C#c0#h3#n0#u0#s#a!" "C"] :bonds [[0 1 "1#a"]]}"#),
        AtomTypingError::NoMatch { atom_id: AtomId(0), element: Element::C, charge: Some(0) }
    )]
    #[case::mismatch_after_open_order(
        WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["C#c0" "*" "N#c0"] :bonds [[0 1 "*"]]}"#),
        AtomTypingError::NoMatch { atom_id: AtomId(2), element: Element::N, charge: Some(0) },
    )]
    fn test_atom_typing_valence_admit_error(
        atom_type_registry: AtomTypeRegistry,
        #[case] wildcards: WildcardPolicy,
        #[case] molecule: Molecule,
        #[case] expected: AtomTypingError,
    ) {
        assert_eq!(
            AtomTypingValence::with_config(&atom_type_registry, wildcards).admit(&molecule),
            Solution::Contradictory(expected)
        );
    }

    #[rstest]
    #[case::singleton(vec!["C#c0#h4#n0#u0#s#v0#a!"], r#"{:atoms ["C#c0#h4#n0#u0#s"]}"#, Some(4), Some(4))]
    #[case::duplicate(vec!["C#c0#h4#v0#a!", "C#c0#h4#v0#a!"], r#"{:atoms ["C#c0#h4#n0#u0#s"]}"#, Some(4), Some(4))]
    #[case::same_h(vec!["C#c0#h2#n0#v0", "C#c0#h2#n1#v0"], r#"{:atoms ["C#c0#h2#n1#u0#s"]}"#, Some(2), Some(2))]
    #[case::plural(vec!["C#c0#h2#n1#v0", "C#c0#h4#n0#v0"], r#"{:atoms ["C#c0#h2#n1#u0#s"]}"#, None, Some(4))]
    #[case::reversed(vec!["C#c0#h4#n0#v0", "C#c0#h2#n1#v0"], r#"{:atoms ["C#c0#h2#n1#u0#s"]}"#, None, Some(4))]
    #[case::different_stored_h(vec!["C#c0#h4#v0"], r#"{:atoms ["C#c0#h2#n1#u0#s"]}"#, Some(4), Some(4))]
    #[case::zero(vec!["C#c0#h0#v0"], r#"{:atoms ["C#c0#h0#n2#u0#s"]}"#, Some(0), Some(0))]
    #[case::element(vec!["N#c0#h3#v0", "C#c0#h4#v0"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, Some(4), Some(4))]
    #[case::charge(vec!["C#c+#h3#v0", "C#c0#h4#v0"], r#"{:atoms ["C#c+#h3#u0#s"]}"#, Some(3), Some(3))]
    #[case::unpaired(vec!["C#c0#h3#u1#s2#v0", "C#c0#h4#u0#s#v0"], r#"{:atoms ["C#c0#h3#u1#s2"]}"#, Some(3), Some(3))]
    #[case::multiplicity(vec!["C#c0#h2#u2#s3#v0", "C#c0#h0#u2#s#v0"], r#"{:atoms ["C#c0#h2#u2#s3"]}"#, Some(2), Some(2))]
    #[case::valence(vec!["C#c0#h4#v0", "C#c0#h3#v1"], r#"{:atoms ["C#c0#h3#u0#s" "C"] :bonds [[0 1 "1"]]}"#, Some(3), Some(3))]
    #[case::aromatic(vec!["C#c0#h0#v0#a0", "C#c0#h1#v0#a1", "C#c0#h2#v0#a2"], r#"{:atoms ["C#c0#h1#u0#s#a1"]}"#, Some(1), Some(1))]
    #[case::both_valences(vec!["C#c0#h4#v0#a1", "C#c0#h2#v1#a2", "C#c0#h1#v1#a1"], r#"{:atoms ["C#c0#h1#u0#s#a1" "C"] :bonds [[0 1 "1"]]}"#, Some(1), Some(1))]
    #[case::not_aromatic(vec!["C#c0#h1#a1#v0", "C#c0#h4#a!#v0"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, Some(4), Some(4))]
    #[case::other_constraints(vec!["C#c0#h4#x2#y3#V8#H6#D7#v0", "C#c0#h4#x0#y0#V4#H4#D4#v0"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, Some(4), Some(4))]
    #[case::empty(vec![], r#"{:atoms ["C#c0#h4#u0#s"]}"#, None, None)]
    #[case::no_match(vec!["C#c+#h3#v0"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, None, None)]
    #[case::nonliteral_h(vec!["C#c0#h*#v0", "C#c0#h4#v0"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, None, None)]
    #[case::filtered_nonliteral_h(vec!["C#c0#h*#v1", "C#c0#h4#v0"], r#"{:atoms ["C#c0#h4#u0#s"]}"#, Some(4), Some(4))]
    #[case::open_charge(vec!["C#c0#h4#v0"], r#"{:atoms ["C#h4#u0#s"]}"#, None, None)]
    #[case::open_aromatic(vec!["C#c0#h1#a1#v0"], r#"{:atoms ["C#c0#h1#u0#s#a*"]}"#, None, None)]
    #[case::stored_valence(vec!["C#c0#h4#v0", "C#c0#h3#v1"], r#"{:atoms ["C#c0#h4#u0#s#v1"]}"#, None, None)]
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
        None, None
    )]
    #[case::aromatic_zero(vec!["C#c0#h2#a0#v0", "C#c0#h4#a!#v0"], r#"{:atoms ["C#c0#h2#u0#s#a0"]}"#, Some(2), Some(2))]
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
