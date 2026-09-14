//! Counts valence resolver: enumerates candidate states from the first table
//! covalence ≥ `v` (targets sorted smallest to largest), splitting
//! `covalence − v` between implicit H and aromatic covalence, then assigning
//! lone pairs and unpaired electrons from the nonbonding budget. Literals
//! constrain each step; singleton candidate sets become edits, plural sets
//! become completions.

use smallvec::SmallVec;
use thiserror::Error;
use umol_chem::element::Element;
use umol_chem::spin::{SpinState, UnpairedElectrons};
#[cfg(test)]
use umol_graph_ir::ir::MoleculeEntries;
use umol_graph_ir::ir::{
    aromatic_covalence, AromaticValence, AromaticValenceForm, AsLit, AtomConstraintForm,
    AtomConstraintKey, AtomConstraintsForm, AtomForm, AtomId, Lattice, Molecule, NumForm,
    UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

use super::{AtomCompletions, ValenceEntry, ValenceTable};
use crate::ops::model::ValenceTieBreak;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum CountsError {
    #[error("no matching valence state")]
    NoMatch,
    #[error("element out of scope: no valence table entry")]
    InvalidElement,
    #[error("aromatic valence unspecified (#a+): no valence table entry")]
    UndeterminedAromaticValence,
}

/// Atom that no valence-table state admits.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
#[error("no valence-table state: element {element}, charge {charge}, valence {valence}")]
pub struct CountsMismatch {
    pub element: Element,
    pub charge: i64,
    pub valence: i64,
}

#[derive(Clone, Copy, Debug)]
struct CountsInput {
    valence: i64,
    accepted_pairs: i64,
    is_aromatic: bool,
}

impl CountsInput {
    fn for_molecule_atom(molecule: &Molecule, atom_id: AtomId) -> Self {
        let atom = molecule.atom(atom_id);
        let constraints = atom.constraints();
        Self {
            valence: atom.valence().as_lit().unwrap_or(0),
            accepted_pairs: atom.accepted_pairs().as_lit().unwrap_or(0),
            // Aromatic evidence: membership in a stored system (a relation
            // fact, read derived) or the closed-world assertion reading,
            // which merges both mark placements and closes to `NotAromatic`.
            is_aromatic: matches!(
                constraints.derived(AtomConstraintKey::AromaticValence),
                Some(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(_)
                ))
            ) || matches!(
                constraints.asserted_complete(AtomConstraintKey::AromaticValence),
                Some(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(_)
                ))
            ),
        }
    }
}

#[derive(Clone, Debug)]
pub struct CountsValence<'a> {
    table: &'a ValenceTable,
}

impl<'a> CountsValence<'a> {
    pub fn new(table: &'a ValenceTable) -> Self {
        Self { table }
    }

    /// Admission: every atom under resolution gets its candidate set —
    /// the counts enumeration per atom — and no edits are produced. A
    /// non-literal element makes the whole admission underdetermined and
    /// empty; plurality is state, not a solution.
    /// Isotope information is preserved without defaulting or selection.
    pub fn admit(&self, molecule: &Molecule) -> Solution<AtomCompletions, CountsError> {
        for atom in molecule.atoms().iter() {
            if atom.element().as_lit().is_none() {
                return Solution::Underdetermined(AtomCompletions::new());
            }
        }

        let mut completions = AtomCompletions::new();
        for id in molecule.atoms().ids() {
            match self.admitted_completions(molecule, id) {
                Ok(Some(candidates)) => completions.insert(id, candidates),
                Ok(None) => {}
                Err(contradiction) => return Solution::Contradictory(contradiction),
            }
        }
        Solution::Determined(completions)
    }

    fn admitted_completions(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
    ) -> Result<Option<SmallVec<[AtomForm; 1]>>, CountsError> {
        let atom = molecule.atom(atom_id);
        if atom.is_ground() {
            // A field-ground atom whose aromatic evidence is present but
            // whose contribution is undetermined is still under resolution:
            // the system entity needs its contribution. Membership in a
            // stored system determines the contribution, so those atoms are
            // complete.
            let contribution_open = matches!(
                atom.constraints()
                    .asserted_complete(AtomConstraintKey::AromaticValence),
                Some(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(NumForm::Undetermined)
                ))
            ) && !atom.is_in_aromatic_system();
            if !contribution_open {
                return Ok(None);
            }
            let element = atom.element().as_lit().expect("ground element");
            let charge = atom.charge().as_lit().expect("ground charge");
            let accepted_pairs = atom.accepted_pairs().as_lit().unwrap_or(0);
            let Some(entry) = element
                .shift((2 * accepted_pairs - charge) as i8)
                .and_then(|shifted| self.table.entry(shifted))
            else {
                return Ok(None);
            };
            if entry.aromatic_valences.is_empty() {
                return Ok(None);
            }
            let mut candidates = SmallVec::new();
            for &aromatic_valence in &entry.aromatic_valences {
                let mut candidate = atom.attributes.clone();
                candidate
                    .constraints
                    .set(AtomConstraintForm::aromatic_valence(
                        AromaticValenceForm::aromatic(NumForm::Lit(i64::from(aromatic_valence))),
                    ));
                candidates.push(candidate);
            }
            return Ok(Some(candidates));
        }
        if atom.element().is_undetermined() {
            return Ok(None);
        };
        if atom.charge().is_undetermined() {
            return Ok(None);
        };

        if atom.valence().as_lit().is_none() {
            return Ok(None);
        }
        let input = CountsInput::for_molecule_atom(molecule, atom_id);
        let candidates = self.candidate_states(atom.attributes, input)?;
        Ok(Some(candidates))
    }

    /// Every candidate state admitted by the table and the atom's literals,
    /// in enumeration order (implicit hydrogens ascending, then the table's
    /// aromatic valences). The fallback aromatic-valence list is consulted
    /// only when the primary list admits no candidate.
    fn candidate_states(
        &self,
        atom: &AtomForm,
        input: CountsInput,
    ) -> Result<SmallVec<[AtomForm; 1]>, CountsError> {
        let CountsInput {
            valence,
            accepted_pairs,
            is_aromatic,
        } = input;
        let element = atom.element.as_lit().unwrap();
        let charge = atom.charge.as_lit().unwrap();

        let entry = element
            .shift((2 * accepted_pairs - charge) as i8)
            .and_then(|shifted| self.table.entry(shifted));

        let aromatic_constraint = atom
            .constraints
            .aromatic_valence()
            .unwrap_or(&AromaticValenceForm::Undetermined);
        if entry.is_none()
            && matches!(
                aromatic_constraint,
                AromaticValenceForm::Aromatic(NumForm::Undetermined)
            )
        {
            return Err(CountsError::UndeterminedAromaticValence);
        }

        // Bonding budget is next-largest saturation target - valence.
        // Above largest saturation target, bonding budget is zero.
        let bonding_budget = match entry {
            Some(entry) if atom.implicit_hydrogens.as_lit().is_none() => Some(
                entry
                    .target_covalences
                    .iter()
                    .find_map(|&c| {
                        let c = i64::from(c);
                        (c >= valence).then(|| c - valence)
                    })
                    .unwrap_or(0),
            ),
            _ => None,
        };

        let hydrogen_values = candidate_implicit_hydrogens(
            &atom.implicit_hydrogens,
            bonding_budget,
            entry.is_none(),
        )?;

        let mut candidates = SmallVec::new();
        for aromatic_values in candidate_aromatic_valences(aromatic_constraint, is_aromatic, entry)
        {
            for &implicit_hydrogens in &hydrogen_values {
                if !atom
                    .implicit_hydrogens
                    .matches(&NumForm::Lit(implicit_hydrogens))
                {
                    continue;
                }
                for &aromatic_valence in &aromatic_values {
                    if !aromatic_constraint.matches_value(aromatic_valence) {
                        continue;
                    }
                    if let Some(b) = bonding_budget {
                        if implicit_hydrogens + aromatic_covalence(aromatic_valence) > b {
                            continue;
                        }
                    }
                    let electron_budget = i64::from(element.valence_electrons()) - charge;
                    let nonbonding =
                        electron_budget - valence - aromatic_valence - implicit_hydrogens;
                    if nonbonding < 0 {
                        continue;
                    }
                    let Some((lone_pairs, unpaired_electrons)) =
                        derive_lone_pairs_and_unpaired_electrons(atom, element, nonbonding)
                    else {
                        continue;
                    };
                    let Some(multiplicity) =
                        derive_multiplicity(&atom.unpaired_electrons, unpaired_electrons)
                    else {
                        continue;
                    };
                    let derived = derive_atom(
                        implicit_hydrogens,
                        unpaired_electrons,
                        multiplicity,
                        lone_pairs,
                        valence,
                        is_aromatic,
                        aromatic_valence,
                    );
                    if let Some(candidate) = atom.meet(&derived) {
                        candidates.push(candidate);
                    }
                }
            }
            if !candidates.is_empty() {
                break;
            }
        }

        if candidates.is_empty() {
            return Err(CountsError::NoMatch);
        }
        Ok(candidates)
    }

    /// Classify molecule atom (including ground atoms) against valence table:
    /// - `Determined` if some state admits it.
    /// - `Contradictory` if no consistent state exists.
    /// - `Underdetermined` if atom is not ground.
    pub fn classify_molecule_atom(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
    ) -> Solution<(), CountsMismatch> {
        let atom = molecule.atom(atom_id);
        if !atom.is_ground() {
            return Solution::Underdetermined(());
        }
        let Some(element) = atom.element().as_lit() else {
            return Solution::Underdetermined(());
        };
        let charge = atom.charge().as_lit().unwrap_or(0);
        let input = CountsInput::for_molecule_atom(molecule, atom_id);
        match self.candidate_states(atom.attributes, input) {
            Ok(_) => Solution::Determined(()),
            Err(_) => Solution::Contradictory(CountsMismatch {
                element,
                charge,
                valence: input.valence,
            }),
        }
    }

    /// Infer implicit H under the valence policy, ignoring stored H and lone pairs.
    /// Stored #v and #a take precedence over derivation from bonds and aromatic systems.
    /// Uses the first table covalence at or above #v and retains charge/spin evidence.
    /// Returns None for missing concrete evidence, no admissible H count, or distinct H
    /// counts under Strict. MostSaturated selects the greatest admissible count.
    /// Neither the molecule nor the table is modified.
    pub fn infer_implicit_hydrogens(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
        policy: ValenceTieBreak,
    ) -> Option<i64> {
        let atom = molecule.atoms().get(atom_id)?;
        let element = atom.element().as_lit()?;
        let charge = atom.charge().as_lit()?;
        let constraints = atom.constraints();
        let valence = match constraints.valence() {
            Some(asserted) => asserted.as_lit()?,
            None => atom.valence().as_lit()?,
        };
        let aromatic = match constraints.aromatic_valence() {
            Some(asserted) => asserted.as_lit()?.valence_count(),
            None => atom.aromatic_valence().as_lit()?,
        };
        let accepted_pairs = atom.accepted_pairs().as_lit()?;
        let shifted_number = i64::from(element.atomic_number())
            .checked_add(accepted_pairs.checked_mul(2)?)?
            .checked_sub(charge)?;
        let entry = u8::try_from(shifted_number)
            .ok()
            .and_then(Element::from_atomic_number)
            .and_then(|shifted| self.table.entry(shifted));
        let bonding_budget = entry
            .and_then(|entry| {
                entry
                    .target_covalences
                    .iter()
                    .map(|&c| i64::from(c))
                    .find(|&c| c >= valence)
            })
            .map(|c| c.checked_sub(valence))
            .unwrap_or(Some(0))?;
        let nonbonding_budget = i64::from(element.valence_electrons())
            .checked_sub(charge)?
            .checked_sub(valence)?
            .checked_sub(aromatic)?;
        let spin = atom.unpaired_electrons();
        let unpaired_count = spin.count.as_lit();
        let max_lone_pairs = i64::from(element.valence_capacity()) / 2;
        let mut hydrogens = None;
        for h in 0..=bonding_budget {
            if entry.is_some() && h.checked_add(aromatic_covalence(aromatic))? > bonding_budget {
                continue;
            }
            let nonbonding = nonbonding_budget.checked_sub(h)?;
            if nonbonding < 0 {
                continue;
            }
            let unpaired = unpaired_count.unwrap_or(nonbonding % 2);
            let remaining = nonbonding.checked_sub(unpaired)?;
            if remaining < 0
                || remaining % 2 != 0
                || remaining / 2 > max_lone_pairs
                || !spin.count.matches(&NumForm::Lit(unpaired))
                || derive_multiplicity(spin, unpaired).is_none()
            {
                continue;
            }
            hydrogens = Some(match (hydrogens, policy) {
                (Some(previous), ValenceTieBreak::Strict) if previous != h => return None,
                _ => h,
            });
        }
        hydrogens
    }
}

fn candidate_implicit_hydrogens(
    implicit_hydrogens: &NumForm,
    bonding_budget: Option<i64>,
    no_entry: bool,
) -> Result<Vec<i64>, CountsError> {
    if let Some(h) = implicit_hydrogens.as_lit() {
        return Ok(vec![h]);
    }
    if no_entry {
        if implicit_hydrogens.is_undetermined() {
            return Ok(vec![0]);
        }
        return Err(CountsError::NoMatch);
    }
    let b = bonding_budget.ok_or(CountsError::NoMatch)?;
    Ok((0..=b).collect())
}

/// Aromatic-valence lists in trial order: the primary `aromatic_valences`
/// list, then `fallback_aromatic_valences`, consulted only when the primary
/// list admits no candidate.
fn candidate_aromatic_valences(
    aromatic: &AromaticValenceForm,
    is_aromatic: bool,
    entry: Option<&ValenceEntry>,
) -> Vec<Vec<i64>> {
    match aromatic.as_lit().map(AromaticValence::valence_count) {
        Some(a) => vec![vec![a]],
        None => match entry {
            Some(entry) if is_aromatic => {
                [&entry.aromatic_valences, &entry.fallback_aromatic_valences]
                    .into_iter()
                    .filter(|list| !list.is_empty())
                    .map(|list| list.iter().map(|&a| i64::from(a)).collect())
                    .collect()
            }
            _ => vec![vec![0]],
        },
    }
}

fn derive_lone_pairs_and_unpaired_electrons(
    atom: &AtomForm,
    element: Element,
    nonbonding: i64,
) -> Option<(i64, i64)> {
    let max_lone_pairs = i64::from(element.valence_capacity()) / 2;
    match (
        atom.lone_pairs.as_lit(),
        atom.unpaired_electrons.count.as_lit(),
    ) {
        (Some(lone_pairs), Some(unpaired_electrons)) => {
            if unpaired_electrons + 2 * lone_pairs == nonbonding {
                Some((lone_pairs, unpaired_electrons))
            } else {
                None
            }
        }
        (Some(lone_pairs), None) => {
            let unpaired_electrons = nonbonding - 2 * lone_pairs;
            if unpaired_electrons < 0
                || !atom
                    .unpaired_electrons
                    .count
                    .matches(&NumForm::Lit(unpaired_electrons))
            {
                return None;
            }
            Some((lone_pairs, unpaired_electrons))
        }
        (None, Some(unpaired_electrons)) => {
            let remaining = nonbonding - unpaired_electrons;
            if remaining < 0 || remaining % 2 != 0 {
                return None;
            }
            let lone_pairs = remaining / 2;
            if lone_pairs > max_lone_pairs || !atom.lone_pairs.matches(&NumForm::Lit(lone_pairs)) {
                return None;
            }
            Some((lone_pairs, unpaired_electrons))
        }
        (None, None) => {
            let unpaired_electrons = nonbonding % 2;
            let lone_pairs = (nonbonding - unpaired_electrons) / 2;
            if lone_pairs > max_lone_pairs {
                return None;
            }
            Some((lone_pairs, unpaired_electrons))
        }
    }
}

fn derive_atom(
    implicit_hydrogens: i64,
    unpaired_electrons: i64,
    multiplicity: i64,
    lone_pairs: i64,
    valence: i64,
    is_aromatic: bool,
    aromatic_valence: i64,
) -> AtomForm {
    AtomForm {
        implicit_hydrogens: NumForm::Lit(implicit_hydrogens),
        lone_pairs: NumForm::Lit(lone_pairs),
        unpaired_electrons: UnpairedElectronsForm {
            count: NumForm::Lit(unpaired_electrons),
            multiplicity: NumForm::Lit(multiplicity),
        },
        constraints: AtomConstraintsForm::from_iter([
            AtomConstraintForm::Valence(NumForm::Lit(valence)),
            AtomConstraintForm::AromaticValence(if is_aromatic {
                AromaticValenceForm::Aromatic(NumForm::Lit(aromatic_valence))
            } else {
                AromaticValenceForm::NotAromatic
            }),
        ]),
        ..Default::default()
    }
}

fn derive_multiplicity(unpaired_electrons: &UnpairedElectronsForm, count: i64) -> Option<i64> {
    let multiplicity = match unpaired_electrons.multiplicity {
        NumForm::Lit(multiplicity) => multiplicity,
        NumForm::Undetermined => count.checked_add(1)?,
        _ => return None,
    };
    SpinState::try_from(UnpairedElectrons {
        count,
        multiplicity,
    })
    .ok()?;
    Some(multiplicity)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use smallvec::smallvec;
    use umol_graph_ir::ir::{AromaticSystemForm, IsotopeMassForm};
    use umol_graph_ir::{atom_dsl, mol_dsl};

    use super::*;
    use crate::ops::valence::ValenceTable;

    #[rustfmt::skip]
    #[rstest]
    #[case::explicit_triplet(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) },
        2,
        Some(3),
    )]
    #[case::explicit_open_shell_singlet(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(1) },
        2,
        Some(1),
    )]
    #[case::incompatible(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(2) },
        2,
        None,
    )]
    #[case::negative_multiplicity(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(-1) },
        2,
        None,
    )]
    #[case::derived(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Undetermined },
        2,
        Some(3),
    )]
    #[case::negative_count(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Undetermined },
        -1,
        None,
    )]
    #[case::pattern(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::lit_set([1, 3]) },
        2,
        None,
    )]
    #[case::overflow(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Undetermined },
        i64::MAX,
        None,
    )]
    fn test_derive_multiplicity(
        #[case] unpaired_electrons: UnpairedElectronsForm,
        #[case] count: i64,
        #[case] expected: Option<i64>,
    ) {
        assert_eq!(derive_multiplicity(&unpaired_electrons, count), expected);
    }

    #[rstest]
    #[case::singleton_methane(
        mol_dsl!(r#"{:atoms ["C#c0#h4"]}"#),
        Solution::Determined({
            let mut completions = AtomCompletions::new();
            completions.insert(
                AtomId(0),
                smallvec![atom_dsl!("C#c0#h4#n0#u0#s#v0#a!")],
            );
            completions
        })
    )]
    #[case::plural_bare_carbon(
        mol_dsl!(r#"{:atoms ["C#c0"]}"#),
        Solution::Determined({
            let mut completions = AtomCompletions::new();
            completions.insert(
                AtomId(0),
                smallvec![
                    atom_dsl!("C#c0#h0#n2#u0#s#v0#a!"),
                    atom_dsl!("C#c0#h#n#u#s2#v0#a!"),
                    atom_dsl!("C#c0#h2#n#u0#s#v0#a!"),
                    atom_dsl!("C#c0#h3#n0#u#s2#v0#a!"),
                    atom_dsl!("C#c0#h4#n0#u0#s#v0#a!"),
                ],
            );
            completions
        })
    )]
    #[case::water_singletons(
        mol_dsl!(r#"{:atoms ["O #c0" "H #c0" "H #c0"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        Solution::Determined({
            let mut completions = AtomCompletions::new();
            completions.insert(
                AtomId(0),
                smallvec![atom_dsl!("O#c0#h0#n2#u0#s#v2#a!")],
            );
            for atom in 1..3 {
                completions.insert(
                    AtomId(atom),
                    smallvec![atom_dsl!("H#c0#h0#n0#u0#s#v#a!")],
                );
            }
            completions
        })
    )]
    #[case::ethane_carbons_plural(
        mol_dsl!(r#"{:atoms ["C #c0" "C #c0"] :bonds [[0 1 "1"]]}"#),
        Solution::Determined({
            let mut completions = AtomCompletions::new();
            let disjuncts: SmallVec<[AtomForm; 1]> = smallvec![
                atom_dsl!("C#c0#h0#n#u#s2#v#a!"),
                atom_dsl!("C#c0#h#n#u0#s#v#a!"),
                atom_dsl!("C#c0#h2#n0#u#s2#v#a!"),
                atom_dsl!("C#c0#h3#n0#u0#s#v#a!"),
            ];
            completions.insert(AtomId(0), disjuncts.clone());
            completions.insert(AtomId(1), disjuncts);
            completions
        })
    )]
    fn test_counts_valence_admit(
        #[case] molecule: Molecule,
        #[case] expected: Solution<AtomCompletions, CountsError>,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        assert_eq!(resolver.admit(&molecule), expected);
    }

    #[rstest]
    #[case::undetermined(IsotopeMassForm::Undetermined)]
    #[case::natural(IsotopeMassForm::Natural)]
    #[case::mass(IsotopeMassForm::Lit(13))]
    #[case::set(IsotopeMassForm::lit_set([12, 13]))]
    #[case::variable(IsotopeMassForm::var("mass"))]
    fn test_counts_valence_admit_isotope(#[case] isotope: IsotopeMassForm) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm {
                isotope_mass: isotope.clone(),
                ..atom_dsl!("C#c0#h4")
            }],
            ..Default::default()
        });
        assert_eq!(
            CountsValence::new(ValenceTable::default_table()).admit(&molecule),
            Solution::Determined(AtomCompletions::from_iter([(
                AtomId(0),
                smallvec![AtomForm {
                    isotope_mass: isotope,
                    ..atom_dsl!("C#c0#h4#n0#u0#s#v0#a!")
                }]
            )]))
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::bond_marked_ring(
        mol_dsl!(r#"{
            :atoms ["C#i=#c0#h1#n0#u0#s" "C#i=#c0#h1#n0#u0#s" "C#i=#c0#h1#n0#u0#s"
                    "C#i=#c0#h1#n0#u0#s" "C#i=#c0#h1#n0#u0#s" "C#i=#c0#h1#n0#u0#s"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"]
                    [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"]]}"#),
        vec!["C#i=#c0#h#n0#u0#s#a"; 6]
    )]
    fn test_counts_valence_admit_ground_evidenced(
        #[case] molecule: Molecule,
        #[case] expected: Vec<&str>,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        let Solution::Determined(completions) = resolver.admit(&molecule) else {
            panic!("ground-evidenced admission did not determine");
        };
        let admitted: Vec<String> = completions
            .iter()
            .flat_map(|(_, forms)| forms.iter().map(ToString::to_string))
            .collect();
        assert_eq!(admitted, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::contribution_asserted(mol_dsl!(r#"{:atoms ["C#i=#c0#h1#n0#u0#s#a2"]}"#))]
    #[case::in_system(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom_dsl!("C#i=#c0#h1#n0#u0#s")],
            aromatic: vec![(vec![AtomId(0)], AromaticSystemForm::from_electrons(vec![1]))],
            ..Default::default()
        })
    )]
    #[case::unmarked(mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#))]
    fn test_counts_valence_admit_ground_complete(#[case] molecule: Molecule) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        assert_eq!(
            resolver.admit(&molecule),
            Solution::Determined(AtomCompletions::new())
        );
    }

    #[rstest]
    fn test_counts_valence_admit_identity() {
        let resolver = CountsValence::new(ValenceTable::default_table());
        let molecule = mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#a!"]}"#);
        assert_eq!(
            resolver.admit(&molecule),
            Solution::Determined(AtomCompletions::new())
        );
    }

    #[rstest]
    #[case::later_undetermined_element(mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#))]
    fn test_counts_valence_admit_partial(#[case] molecule: Molecule) {
        assert_eq!(
            CountsValence::new(ValenceTable::default_table()).admit(&molecule),
            Solution::Underdetermined(AtomCompletions::new())
        );
    }

    #[rstest]
    #[case::later_atom_contradiction(mol_dsl!(r#"{:atoms ["C#c0#h4" "Fe#c0#h0#a+"]}"#), CountsError::UndeterminedAromaticValence)]
    fn test_counts_valence_admit_error(#[case] molecule: Molecule, #[case] expected: CountsError) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        assert_eq!(resolver.admit(&molecule), Solution::Contradictory(expected));
    }

    #[rstest]
    #[case::methane(r#"{:atoms ["C#c0#h4#n0#u0#s"]}"#, None, Some(4))]
    #[case::different_stored_h(r#"{:atoms ["C#c0#h2#n1#u0#s"]}"#, None, Some(4))]
    #[case::unresolved_h_lp(r#"{:atoms ["C#c0#u0#s"]}"#, None, Some(4))]
    #[case::methyl(r#"{:atoms ["C#c0#h3#u1#s2"]}"#, None, Some(3))]
    #[case::triplet(r#"{:atoms ["C#c0#h2#u2#s3"]}"#, None, Some(2))]
    #[case::open_shell_singlet(r#"{:atoms ["C#c0#h2#u2#s"]}"#, None, Some(2))]
    #[case::incompatible_spin(r#"{:atoms ["C#c0#u2#s2"]}"#, None, None)]
    #[case::open_spin(r#"{:atoms ["C#c0"]}"#, None, Some(4))]
    #[case::water(r#"{:atoms ["O#c0#h2#n2#u0#s"]}"#, None, Some(2))]
    #[case::fluorine(r#"{:atoms ["F#c0#h1#n3#u0#s"]}"#, Some(1), Some(1))]
    #[case::anion(r#"{:atoms ["C#c-#h3#u0#s"]}"#, None, Some(3))]
    #[case::cation(r#"{:atoms ["C#c+#h3#u0#s"]}"#, None, Some(3))]
    #[case::derived_valence(r#"{:atoms ["C#c0#h3#u0#s" "C"] :bonds [[0 1 "1"]]}"#, None, Some(3))]
    #[case::stored_valence(r#"{:atoms ["C#c0#h1#u0#s#v3"]}"#, Some(1), Some(1))]
    #[case::other_constraints(r#"{:atoms ["C#c0#h4#u0#s#x2#y3#V8#H6#D7"]}"#, None, Some(4))]
    #[case::aromatic_carbon(r#"{:atoms ["C#c0#u0#s#v2#a1"]}"#, Some(1), Some(1))]
    #[case::pyridine_nitrogen(r#"{:atoms ["N#c0#u0#s#v2#a1"]}"#, Some(0), Some(0))]
    #[case::pyrrole_nitrogen(r#"{:atoms ["N#c0#u0#s#v2#a2"]}"#, Some(1), Some(1))]
    #[case::aromatic_zero(r#"{:atoms ["C#c0#u0#s#v4#a0"]}"#, Some(0), Some(0))]
    #[case::aromatic_zero_h(r#"{:atoms ["C#c0#u0#s#v2#a0"]}"#, None, Some(2))]
    #[case::aromatic_budget(r#"{:atoms ["C#c0#u0#s#v4#a1"]}"#, None, None)]
    #[case::first_target(r#"{:atoms ["S#c0#u0#s#v0"]}"#, None, Some(2))]
    #[case::at_first_target(r#"{:atoms ["S#c0#u0#s#v2"]}"#, Some(0), Some(0))]
    #[case::between_targets(r#"{:atoms ["S#c0#u0#s#v3"]}"#, Some(1), Some(1))]
    #[case::at_second_target(r#"{:atoms ["S#c0#u0#s#v4"]}"#, Some(0), Some(0))]
    #[case::last_target(r#"{:atoms ["S#c0#u0#s#v5"]}"#, Some(1), Some(1))]
    #[case::above_targets(r#"{:atoms ["Cl#c0#u0#s#v3"]}"#, Some(0), Some(0))]
    #[case::no_electron_budget(r#"{:atoms ["C#c0#u0#s#v5"]}"#, None, None)]
    #[case::no_table_entry(r#"{:atoms ["Fe#c0#u0#s"]}"#, Some(0), Some(0))]
    #[case::lone_pair_capacity(r#"{:atoms ["O#c-4#u0#s"]}"#, Some(2), Some(2))]
    #[case::excess_lone_pairs(r#"{:atoms ["O#c-8#u0#s"]}"#, None, None)]
    #[case::accepted_pairs(
        r#"{:atoms ["C#c0#u0#s" "N"] :dative-bonds [{:donors [1] :acceptor 0 :attrs "1"}]}"#,
        None,
        Some(2)
    )]
    #[case::open_charge(r#"{:atoms ["C#u0#s"]}"#, None, None)]
    #[case::open_valence(r#"{:atoms ["C#c0#u0#s#v*"]}"#, None, None)]
    #[case::open_aromatic(r#"{:atoms ["C#c0#u0#s#a+"]}"#, None, None)]
    #[case::aromatic_system(
        r#"{:atoms ["C#c0#h1#u0#s" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#,
        Some(1),
        Some(1)
    )]
    #[case::stored_aromatic(
        r#"{:atoms ["C#c0#h1#u0#s#a0" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#,
        None,
        Some(2)
    )]
    fn test_counts_valence_infer_implicit_hydrogens(
        #[case] input: &str,
        #[case] strict: Option<i64>,
        #[case] most_saturated: Option<i64>,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
    ) {
        let molecule = mol_dsl!(input);
        let original = molecule.clone();
        let table = ValenceTable::smiles_table();
        let original_table = table.clone();
        let expected = match policy {
            ValenceTieBreak::Strict => strict,
            ValenceTieBreak::MostSaturated => most_saturated,
        };
        assert_eq!(
            CountsValence::new(table).infer_implicit_hydrogens(&molecule, AtomId(0), policy),
            expected,
        );
        assert_eq!(molecule, original);
        assert_eq!(table, &original_table);
        assert_eq!(table.content_hash(), original_table.content_hash());
    }

    #[rstest]
    #[case::smiles(ValenceTable::smiles_table(), None, None)]
    #[case::mdl(ValenceTable::mdl_table(), Some(1), Some(1))]
    #[case::general(ValenceTable::default_table(), Some(1), Some(1))]
    fn test_counts_valence_infer_implicit_hydrogens_table(
        #[case] table: &ValenceTable,
        #[case] strict: Option<i64>,
        #[case] most_saturated: Option<i64>,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
    ) {
        let molecule = mol_dsl!(r#"{:atoms ["Cl#c0#u0#s#v2"]}"#);
        let expected = match policy {
            ValenceTieBreak::Strict => strict,
            ValenceTieBreak::MostSaturated => most_saturated,
        };
        assert_eq!(
            CountsValence::new(table).infer_implicit_hydrogens(&molecule, AtomId(0), policy),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::methane_h("C#c0#h4", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["C#c0#h4#n0#u0#s#v0#a!"])]
    #[case::methane_h_inferred("C#c0#h*", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec![
        "C#c0#h0#n2#u0#s#v0#a!",
        "C#c0#h#n#u#s2#v0#a!",
        "C#c0#h2#n#u0#s#v0#a!",
        "C#c0#h3#n0#u#s2#v0#a!",
        "C#c0#h4#n0#u0#s#v0#a!",
    ])]
    #[case::ammonia("N#c0#h3", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["N#c0#h3#n#u0#s#v0#a!"])]
    #[case::water("O#c0#h2", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["O#c0#h2#n2#u0#s#v0#a!"])]
    #[case::methyl_radical("C#c0#h3", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["C#c0#h3#n0#u#s2#v0#a!"])]
    #[case::methyl_radical_h_inferred("C#c0#u", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec![
        "C#c0#h#n#u#s2#v0#a!",
        "C#c0#h3#n0#u#s2#v0#a!",
    ])]
    #[case::methyl_anion("C#c-1#h3", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["C#c-#h3#n#u0#s#v0#a!"])]
    #[case::hydroxyl_radical("O#c0#h1", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["O#c0#h#n2#u#s2#v0#a!"])]
    #[case::fluoride("F#c-1#h0", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["F#c-#h0#n4#u0#s#v0#a!"])]
    #[case::fluorine_atom("F#c0#h0", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["F#c0#h0#n3#u#s2#v0#a!"])]
    #[case::magnesium_atom("Mg#c0#h0", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["Mg#c0#h0#n#u0#s#v0#a!"])]
    #[case::ethane_carbon("C#c0#v1", CountsInput { valence: 1, accepted_pairs: 0, is_aromatic: false }, vec![
        "C#c0#h0#n#u#s2#v#a!",
        "C#c0#h#n#u0#s#v#a!",
        "C#c0#h2#n0#u#s2#v#a!",
        "C#c0#h3#n0#u0#s#v#a!",
    ])]
    #[case::methylene_carbon("C#c0#v2", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: false }, vec![
        "C#c0#h0#n#u0#s#v2#a!",
        "C#c0#h#n0#u#s2#v2#a!",
        "C#c0#h2#n0#u0#s#v2#a!",
    ])]
    #[case::methine_carbon("C#c0#v3", CountsInput { valence: 3, accepted_pairs: 0, is_aromatic: false }, vec![
        "C#c0#h0#n0#u#s2#v3#a!",
        "C#c0#h#n0#u0#s#v3#a!",
    ])]
    #[case::amine_nitrogen("N#c0#v1", CountsInput { valence: 1, accepted_pairs: 0, is_aromatic: false }, vec![
        "N#c0#h0#n2#u0#s#v#a!",
        "N#c0#h#n#u#s2#v#a!",
        "N#c0#h2#n#u0#s#v#a!",
    ])]
    #[case::alcohol_oxygen("O#c0#v1", CountsInput { valence: 1, accepted_pairs: 0, is_aromatic: false }, vec![
        "O#c0#h0#n2#u#s2#v#a!",
        "O#c0#h#n2#u0#s#v#a!",
    ])]
    #[case::benzene_carbon("C#c0#v2#h1#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec!["C#c0#h#n0#u0#s#v2#a"])]
    #[case::benzene_carbon_h_inferred("C#c0#v2#h*#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec![
        "C#c0#h0#n0#u#s2#v2#a",
        "C#c0#h#n0#u0#s#v2#a",
    ])]
    #[case::fused_aromatic_carbon_h_inferred("C#c0#v3#h*#a+", CountsInput { valence: 3, accepted_pairs: 0, is_aromatic: true }, vec!["C#c0#h0#n0#u0#s#v3#a"])]
    #[case::aromatic_carbon_unpaired("C#c0#v2#h*#u1#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec!["C#c0#h0#n0#u#s2#v2#a"])]
    #[case::pyridine_nitrogen("N#c0#v2#h0#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec![
        "N#c0#h0#n#u0#s#v2#a",
        "N#c0#h0#n0#u#s2#v2#a2",
    ])]
    #[case::pyrrole_nitrogen("N#c0#v2#h1#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec![
        "N#c0#h#n0#u#s2#v2#a",
        "N#c0#h#n0#u0#s#v2#a2",
    ])]
    #[case::furan_oxygen("O#c0#v2#h0#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec!["O#c0#h0#n#u0#s#v2#a2"])]
    #[case::furan_oxygen_h_inferred("O#c0#v2#h*#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec!["O#c0#h0#n#u0#s#v2#a2"])]
    #[case::borazine_boron("B#c0#v2#h1#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec!["B#c0#h#n0#u0#s#v2#a0"])]
    #[case::cyclopentadienyl_carbanion("C#c-1#v2#h1#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec![
        "C#c-#h#n0#u#s2#v2#a",
        "C#c-#h#n0#u0#s#v2#a2",
    ])]
    #[case::tropylium_carbocation("C#c1#v2#h1#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec!["C#c+#h#n0#u0#s#v2#a0"])]
    #[case::exocyclic_carbonyl_carbon("C#c0#v4#h*#a+", CountsInput { valence: 4, accepted_pairs: 0, is_aromatic: true }, vec!["C#c0#h0#n0#u0#s#v4#a0"])]
    #[case::sp3_carbon_h_pinned("C#c0#v2#h2#a+", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: true }, vec!["C#c0#h2#n0#u0#s#v2#a0"])]
    #[case::iron_out_of_table("Fe#c0#h0", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["Fe#c0#h0#n4#u0#s#v0#a!"])]
    fn test_counts_valence_candidate_states(
        #[case] input: &str,
        #[case] counts_input: CountsInput,
        #[case] expected: Vec<&str>,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        let candidates = resolver
            .candidate_states(&atom_dsl!(input), counts_input)
            .unwrap();
        assert_eq!(
            candidates
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>(),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::chloronium("Cl#c1", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: false }, vec!["Cl#c+#h0#n2#u0#s#v2#a!"])]
    #[case::chlorine_trifluoride("Cl#c0", CountsInput { valence: 3, accepted_pairs: 0, is_aromatic: false }, vec!["Cl#c0#h0#n2#u0#s#v3#a!"])]
    #[case::divalent_fluorine("F#c0", CountsInput { valence: 2, accepted_pairs: 0, is_aromatic: false }, vec!["F#c0#h0#n2#u#s2#v2#a!"])]
    fn test_counts_valence_candidate_states_saturated(
        #[case] input: &str,
        #[case] counts_input: CountsInput,
        #[case] expected: Vec<&str>,
    ) {
        let resolver = CountsValence::new(ValenceTable::smiles_table());
        let candidates = resolver
            .candidate_states(&atom_dsl!(input), counts_input)
            .unwrap();
        assert_eq!(
            candidates
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>(),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_aromatic_out_of_table("Fe#c0#h0#a+", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: true }, CountsError::UndeterminedAromaticValence)]
    #[case::over_valence("C#c0", CountsInput { valence: 5, accepted_pairs: 0, is_aromatic: false }, CountsError::NoMatch)]
    fn test_counts_valence_candidate_states_error(
        #[case] input: &str,
        #[case] counts_input: CountsInput,
        #[case] expected: CountsError,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        assert_eq!(
            resolver.candidate_states(&atom_dsl!(input), counts_input),
            Err(expected)
        );
    }

    #[rstest]
    #[case::methane_conforms("C#i=#c0#h4#n0#u0#s#v0#a!", Solution::Determined(()))]
    #[case::excess_lone_pairs("C#i=#c0#h4#n2#u0#s#v0#a!", Solution::Contradictory(CountsMismatch { element: Element::C, charge: 0, valence: 0, }))]
    #[case::not_ground("C", Solution::Underdetermined(()))]
    fn test_counts_valence_classify_molecule_atom(
        #[case] input: &str,
        #[case] expected: Solution<(), CountsMismatch>,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table());
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
