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
    AtomConstraintKey, AtomConstraintsForm, AtomForm, AtomHandle, AtomId, AtomView, Edits,
    IsotopeMassForm, Lattice, Molecule, NumForm, TransactionError, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

use super::{AtomCompletions, ValenceTable};

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
    fn for_molecule_atom(atom: AtomView<'_>) -> Self {
        let constraints = atom.constraints();
        Self {
            valence: atom.valence().as_lit().unwrap_or(0),
            accepted_pairs: atom.accepted_pairs().as_lit().unwrap_or(0),
            is_aromatic: matches!(
                constraints.derived(AtomConstraintKey::AromaticValence),
                Some(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(_)
                ))
            ) || constraints
                .aromatic_valence()
                .is_some_and(|a| a.is_aromatic()),
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

    /// Construct the complete plan without mutating `molecule`: singleton
    /// candidate sets produce edits, plural sets produce completions.
    ///
    /// A non-literal element makes the whole plan underdetermined and yields
    /// nothing; a plural candidate set makes it underdetermined and yields
    /// both the singleton edits and the completions.
    pub fn plan(&self, molecule: &Molecule) -> Solution<(Edits, AtomCompletions), CountsError> {
        for atom in molecule.atoms().iter() {
            if atom.element().as_lit().is_none() {
                return Solution::Underdetermined((Edits::new(), AtomCompletions::new()));
            }
        }

        let mut edits = Edits::new();
        let mut completions = AtomCompletions::new();
        for id in molecule.atoms().ids() {
            let candidates = match self.resolve_molecule_atom(molecule, id) {
                Ok(Some(candidates)) => candidates,
                Ok(None) => continue,
                Err(contradiction) => return Solution::Contradictory(contradiction),
            };
            if let [selected] = candidates.as_slice() {
                let current = molecule.atom(id).attributes;
                // The constraint channel holds assertions only: a committed
                // singleton narrows fields; the enumeration's #v/#a values
                // are candidate bookkeeping and are never written back.
                let mut selected = selected.clone();
                selected.constraints = current.constraints.clone();
                let update = current.difference_to(&selected);
                edits.update_atom(AtomHandle::Id(id), current, &update);
            } else {
                completions.insert(id, candidates);
            }
        }
        if completions.is_empty() {
            Solution::Determined((edits, completions))
        } else {
            Solution::Underdetermined((edits, completions))
        }
    }

    /// Plan and atomically apply counts-valence resolution. Singleton edits
    /// are committed even under an underdetermined verdict, whose payload is
    /// the plural survivors.
    pub fn resolve(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<AtomCompletions, CountsError>, TransactionError> {
        let (edits, completions, determined) = match self.plan(molecule) {
            Solution::Determined((edits, completions)) => (edits, completions, true),
            Solution::Underdetermined((edits, completions)) => (edits, completions, false),
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let mut editor = molecule.edit();
        editor.transact(edits)?;
        *molecule = editor.build();
        Ok(if determined {
            Solution::Determined(completions)
        } else {
            Solution::Underdetermined(completions)
        })
    }

    fn resolve_molecule_atom(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
    ) -> Result<Option<SmallVec<[AtomForm; 1]>>, CountsError> {
        let atom = molecule.atom(atom_id);
        if atom.is_ground() {
            return Ok(None);
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
        let input = CountsInput::for_molecule_atom(atom);
        let mut candidates = self.select_candidates(atom.attributes, input)?;
        for candidate in &mut candidates {
            if candidate.isotope_mass.is_undetermined() {
                candidate.isotope_mass = IsotopeMassForm::Natural;
            }
        }
        Ok(Some(candidates))
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
        let input = CountsInput::for_molecule_atom(atom);
        match self.select_candidates(atom.attributes, input) {
            Ok(_) => Solution::Determined(()),
            Err(_) => Solution::Contradictory(CountsMismatch {
                element,
                charge,
                valence: input.valence,
            }),
        }
    }

    /// Every candidate state admitted by the table and the atom's literals,
    /// in enumeration order (implicit hydrogens ascending, then the table's
    /// aromatic valences).
    fn select_candidates(
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

        let bonding_budget = match entry {
            Some(entry) if atom.implicit_hydrogens.as_lit().is_none() => Some(
                entry
                    .target_covalences
                    .iter()
                    .map(|&c| i64::from(c))
                    .find(|&c| c >= valence)
                    .map(|c| c - valence)
                    .ok_or(CountsError::NoMatch)?,
            ),
            _ => None,
        };

        let aromatic_values = candidate_aromatic_valences(
            aromatic_constraint,
            is_aromatic,
            entry.map(|e| e.aromatic_valences.as_slice()),
        );

        let mut candidates = SmallVec::new();
        for implicit_hydrogens in
            candidate_implicit_hydrogens(&atom.implicit_hydrogens, bonding_budget, entry.is_none())?
        {
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
                let nonbonding = electron_budget - valence - aromatic_valence - implicit_hydrogens;
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

        if candidates.is_empty() {
            return Err(CountsError::NoMatch);
        }
        Ok(candidates)
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

fn candidate_aromatic_valences(
    aromatic: &AromaticValenceForm,
    is_aromatic: bool,
    table: Option<&[u8]>,
) -> Vec<i64> {
    match aromatic.as_lit().map(AromaticValence::valence_count) {
        Some(a) => vec![a],
        None => match table {
            Some(table) if is_aromatic => table.iter().map(|&a| i64::from(a)).collect(),
            _ => vec![0],
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
    use umol_graph_ir::ir::{AtomFieldChange, Edit};
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
        Solution::Determined((
            Edits::from_iter([
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::IsotopeMass {
                        old: IsotopeMassForm::Undetermined,
                        new: IsotopeMassForm::Natural,
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::LonePairs {
                        old: NumForm::Undetermined,
                        new: NumForm::Lit(0),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::UnpairedElectrons {
                        old: UnpairedElectronsForm::default(),
                        new: UnpairedElectronsForm::from((0_u8, 1_u8)),
                    },
                },
            ]),
            AtomCompletions::new()
        ))
    )]
    #[case::plural_bare_carbon(
        mol_dsl!(r#"{:atoms ["C#c0"]}"#),
        Solution::Underdetermined((Edits::new(), {
            let mut completions = AtomCompletions::new();
            completions.insert(
                AtomId(0),
                smallvec![
                    atom_dsl!("C#i=#c0#h0#n2#u0#s#v0#a!"),
                    atom_dsl!("C#i=#c0#h#n#u#s2#v0#a!"),
                    atom_dsl!("C#i=#c0#h2#n#u0#s#v0#a!"),
                    atom_dsl!("C#i=#c0#h3#n0#u#s2#v0#a!"),
                    atom_dsl!("C#i=#c0#h4#n0#u0#s#v0#a!"),
                ],
            );
            completions
        }))
    )]
    fn test_counts_valence_plan(
        #[case] molecule: Molecule,
        #[case] expected: Solution<(Edits, AtomCompletions), CountsError>,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        assert_eq!(resolver.plan(&molecule), expected);
    }

    #[rstest]
    fn test_counts_valence_plan_identity() {
        let resolver = CountsValence::new(ValenceTable::default_table());
        let molecule = mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#a!"]}"#);
        assert_eq!(
            resolver.plan(&molecule),
            Solution::Determined((Edits::new(), AtomCompletions::new()))
        );
    }

    #[rstest]
    #[case::later_undetermined_element(mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#))]
    fn test_counts_valence_plan_partial(#[case] molecule: Molecule) {
        assert_eq!(
            CountsValence::new(ValenceTable::default_table()).plan(&molecule),
            Solution::Underdetermined((Edits::new(), AtomCompletions::new()))
        );
    }

    #[rstest]
    #[case::later_atom_contradiction(mol_dsl!(r#"{:atoms ["C#c0#h4" "Fe#c0#h0#a+"]}"#), CountsError::UndeterminedAromaticValence)]
    fn test_counts_valence_plan_error(#[case] molecule: Molecule, #[case] expected: CountsError) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        assert_eq!(resolver.plan(&molecule), Solution::Contradictory(expected));
    }

    #[rstest]
    #[case::water_singletons(
        mol_dsl!(r#"{:atoms ["O #c0" "H #c0" "H #c0"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        Solution::Determined(AtomCompletions::new()),
        vec!["O#i=#c0#h0#n2#u0#s", "H#i=#c0#h0#n0#u0#s", "H#i=#c0#h0#n0#u0#s"]
    )]
    #[case::ethane_carbons_plural(
        mol_dsl!(r#"{:atoms ["C #c0" "C #c0"] :bonds [[0 1 "1"]]}"#),
        Solution::Underdetermined({
            let mut completions = AtomCompletions::new();
            let disjuncts: SmallVec<[AtomForm; 1]> = smallvec![
                atom_dsl!("C#i=#c0#h0#n#u#s2#v#a!"),
                atom_dsl!("C#i=#c0#h#n#u0#s#v#a!"),
                atom_dsl!("C#i=#c0#h2#n0#u#s2#v#a!"),
                atom_dsl!("C#i=#c0#h3#n0#u0#s#v#a!"),
            ];
            completions.insert(AtomId(0), disjuncts.clone());
            completions.insert(AtomId(1), disjuncts);
            completions
        }),
        vec!["C#c0", "C#c0"]
    )]
    #[case::benzene_carbons_plural(
        mol_dsl!(r#"{:atoms ["C #c0" "C #c0" "C #c0" "C #c0" "C #c0" "C #c0"] :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"]]}"#),
        Solution::Underdetermined({
            let mut completions = AtomCompletions::new();
            let disjuncts: SmallVec<[AtomForm; 1]> = smallvec![
                atom_dsl!("C#i=#c0#h0#n0#u#s2#v2#a"),
                atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a"),
            ];
            for atom in 0..6 {
                completions.insert(AtomId(atom), disjuncts.clone());
            }
            completions
        }),
        vec!["C#c0", "C#c0", "C#c0", "C#c0", "C#c0", "C#c0"]
    )]
    fn test_counts_valence_resolve(
        #[case] mut molecule: Molecule,
        #[case] expected: Solution<AtomCompletions, CountsError>,
        #[case] expected_atoms: Vec<&str>,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        assert_eq!(resolver.resolve(&mut molecule), Ok(expected));
        for (atom, expected_atom) in molecule.atoms().iter().zip(expected_atoms) {
            assert_eq!(atom.attributes.to_string(), expected_atom);
        }
    }

    #[rstest]
    #[case::later_undetermined_element(mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#))]
    fn test_counts_valence_resolve_partial(#[case] mut molecule: Molecule) {
        let original = molecule.clone();
        assert_eq!(
            CountsValence::new(ValenceTable::default_table()).resolve(&mut molecule),
            Ok(Solution::Underdetermined(AtomCompletions::new()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::later_atom_contradiction(mol_dsl!(r#"{:atoms ["C#c0#h4" "Fe#c0#h0#a+"]}"#), CountsError::UndeterminedAromaticValence)]
    fn test_counts_valence_resolve_error(
        #[case] mut molecule: Molecule,
        #[case] expected: CountsError,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        let original = molecule.clone();
        assert_eq!(
            resolver.resolve(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
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
    #[case::iron_out_of_table("Fe#c0#h0", CountsInput { valence: 0, accepted_pairs: 0, is_aromatic: false }, vec!["Fe#c0#h0#n4#u0#s#v0#a!"])]
    fn test_counts_valence_select_candidates(
        #[case] input: &str,
        #[case] counts_input: CountsInput,
        #[case] expected: Vec<&str>,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        let candidates = resolver
            .select_candidates(&atom_dsl!(input), counts_input)
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
    fn test_counts_valence_select_candidates_error(
        #[case] input: &str,
        #[case] counts_input: CountsInput,
        #[case] expected: CountsError,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table());
        assert_eq!(
            resolver.select_candidates(&atom_dsl!(input), counts_input),
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
