//! Counts valence resolver: candidate selection picks the first table
//! covalence ≥ `v` (targets sorted smallest to largest), splits `covalence − v`
//! between implicit H and aromatic increment, then assigns lone pairs and unpaired
//! electrons from the nonbonding budget. Literals constrain each step.

use thiserror::Error;
#[cfg(test)]
use umol_ast::ast::MoleculeParts;
use umol_ast::ast::{
    aromatic_increment, AromaticValenceAst, AsLit, AtomAst, AtomConstraintAst, AtomConstraintsAst,
    AtomHandle, AtomId, AtomView, BooleanAst, Edit, IsotopeMassAst, Lattice, MoleculeAst,
    SpinStateAst, TransactionError, ValueAst,
};
use umol_chem::element::Element;
use umol_chem::spin::{SpinMultiplicity, SpinState};
use umol_utils::solution::Solution;

use crate::ops::model::CountsModel;

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
    fn for_atom(atom: &AtomAst) -> Self {
        Self {
            valence: atom
                .constraints
                .valence()
                .unwrap_or(&ValueAst::Undetermined)
                .as_lit_or(0),
            accepted_pairs: atom
                .constraints
                .accepted_pairs()
                .unwrap_or(&ValueAst::Undetermined)
                .as_lit_or(0),
            is_aromatic: atom
                .constraints
                .aromatic_valence()
                .is_some_and(|a| a.is_aromatic()),
        }
    }

    fn for_molecule_atom(atom: AtomView<'_>) -> Self {
        Self {
            valence: atom.valence().as_lit_or(0),
            accepted_pairs: atom.accepted_pairs().as_lit_or(0),
            is_aromatic: atom.is_in_aromatic_system()
                || atom
                    .neighbors()
                    .any(|n| matches!(n.bond().constraints().aromatic(), BooleanAst::Lit(true)))
                || atom
                    .constraints()
                    .aromatic_valence()
                    .is_some_and(|a| a.is_aromatic()),
        }
    }
}

#[derive(Clone, Debug)]
pub struct CountsValence<'a> {
    model: &'a CountsModel,
}

impl<'a> CountsValence<'a> {
    pub fn new(model: &'a CountsModel) -> Self {
        Self { model }
    }

    /// Construct the complete edit plan without mutating `ast`.
    ///
    /// A non-literal element makes the whole plan underdetermined and yields
    /// no edits.
    pub fn plan(&self, ast: &MoleculeAst) -> Solution<Vec<Edit>, CountsError> {
        for atom in ast.atoms().iter() {
            if atom.element().as_lit().is_none() {
                return Solution::Underdetermined(Vec::new());
            }
        }

        let mut edits = Vec::new();
        for id in ast.atoms().ids() {
            let selected = match self.resolve_molecule_atom(ast, id) {
                Ok(Some(selected)) => selected,
                Ok(None) => continue,
                Err(contradiction) => return Solution::Contradictory(contradiction),
            };
            let current = ast.atom(id).ast;
            let update = current.difference_to(&selected);
            edits.extend(Edit::for_atom_update(AtomHandle::Id(id), current, &update));
        }
        Solution::Determined(edits)
    }

    /// Plan and atomically apply counts-valence resolution.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), CountsError>, TransactionError> {
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

    fn resolve_molecule_atom(
        &self,
        ast: &MoleculeAst,
        atom_id: AtomId,
    ) -> Result<Option<AtomAst>, CountsError> {
        let atom = ast.atom(atom_id);
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
        let mut selected = self.select_candidate(atom.ast, input)?;
        if selected.isotope_mass.is_undetermined() {
            selected.isotope_mass = IsotopeMassAst::Natural;
        }
        Ok(Some(selected))
    }

    pub fn resolve_atom(&self, ast: &mut AtomAst) -> Result<(), CountsError> {
        if ast.is_ground() {
            return Ok(());
        }
        if ast.element.is_undetermined() {
            return Ok(());
        };
        if ast.charge.is_undetermined() {
            return Ok(());
        };

        let input = CountsInput::for_atom(ast);
        *ast = self.select_candidate(ast, input)?;
        Ok(())
    }

    /// Classify molecule atom (including ground atoms) against valence table:
    /// - `Determined` if some state admits it.
    /// - `Contradictory` if no consistent state exists.
    /// - `Underdetermined` if atom is not ground.
    pub fn classify_molecule_atom(
        &self,
        ast: &MoleculeAst,
        atom_id: AtomId,
    ) -> Solution<(), CountsMismatch> {
        let atom = ast.atom(atom_id);
        if !atom.is_ground() {
            return Solution::Underdetermined(());
        }
        let Some(element) = atom.element().as_lit() else {
            return Solution::Underdetermined(());
        };
        let charge = atom.charge().as_lit_or(0);
        let input = CountsInput::for_molecule_atom(atom);
        match self.select_candidate(atom.ast, input) {
            Ok(_) => Solution::Determined(()),
            Err(_) => Solution::Contradictory(CountsMismatch {
                element,
                charge,
                valence: input.valence,
            }),
        }
    }

    fn select_candidate(&self, atom: &AtomAst, input: CountsInput) -> Result<AtomAst, CountsError> {
        let CountsInput {
            valence,
            accepted_pairs,
            is_aromatic,
        } = input;
        let element = atom.element.as_lit().unwrap();
        let charge = atom.charge.as_lit().unwrap();

        let entry = element
            .shift((2 * accepted_pairs - charge) as i8)
            .and_then(|shifted| self.model.table.entry(shifted));

        let aromatic_constraint = atom
            .constraints
            .aromatic_valence()
            .unwrap_or(&AromaticValenceAst::Undetermined);
        if entry.is_none()
            && matches!(
                aromatic_constraint,
                AromaticValenceAst::Aromatic(ValueAst::Undetermined)
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

        let mut candidates = Vec::new();
        for implicit_hydrogens in
            candidate_implicit_hydrogens(&atom.implicit_hydrogens, bonding_budget, entry.is_none())?
        {
            if !atom
                .implicit_hydrogens
                .matches(&ValueAst::Lit(implicit_hydrogens))
            {
                continue;
            }
            for &aromatic_valence in &aromatic_values {
                if !aromatic_constraint.matches_value(aromatic_valence) {
                    continue;
                }
                if let Some(b) = bonding_budget {
                    if implicit_hydrogens + aromatic_increment(aromatic_valence) > b {
                        continue;
                    }
                }
                let electron_budget = i64::from(element.valence_electrons()) - charge;
                let nonbonding = electron_budget - valence - aromatic_valence - implicit_hydrogens;
                if nonbonding < 0 {
                    continue;
                }
                let Some((lone_pairs, unpaired)) =
                    derive_lone_pairs_and_unpaired(atom, element, nonbonding)
                else {
                    continue;
                };
                let Some(multiplicity) = derive_multiplicity(&atom.spin, unpaired) else {
                    continue;
                };
                let derived = derive_atom(
                    implicit_hydrogens,
                    unpaired,
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

        let best = candidates
            .into_iter()
            .max_by(super::compare::compare_valence_preference)
            .ok_or(CountsError::NoMatch)?;
        Ok(best)
    }
}

fn candidate_implicit_hydrogens(
    implicit_hydrogens: &ValueAst,
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
    aromatic: &AromaticValenceAst,
    is_aromatic: bool,
    table: Option<&[u8]>,
) -> Vec<i64> {
    match aromatic.as_lit() {
        Some(a) => vec![a],
        None => match table {
            Some(table) if is_aromatic => table.iter().map(|&a| i64::from(a)).collect(),
            _ => vec![0],
        },
    }
}

fn derive_lone_pairs_and_unpaired(
    atom: &AtomAst,
    element: Element,
    nonbonding: i64,
) -> Option<(i64, i64)> {
    let max_lone_pairs = i64::from(element.valence_capacity()) / 2;
    match (atom.lone_pairs.as_lit(), atom.spin.unpaired.as_lit()) {
        (Some(lone_pairs), Some(unpaired)) => {
            if unpaired + 2 * lone_pairs == nonbonding {
                Some((lone_pairs, unpaired))
            } else {
                None
            }
        }
        (Some(lone_pairs), None) => {
            let unpaired = nonbonding - 2 * lone_pairs;
            if unpaired < 0 || !atom.spin.unpaired.matches(&ValueAst::Lit(unpaired)) {
                return None;
            }
            Some((lone_pairs, unpaired))
        }
        (None, Some(unpaired)) => {
            let remaining = nonbonding - unpaired;
            if remaining < 0 || remaining % 2 != 0 {
                return None;
            }
            let lone_pairs = remaining / 2;
            if lone_pairs > max_lone_pairs || !atom.lone_pairs.matches(&ValueAst::Lit(lone_pairs)) {
                return None;
            }
            Some((lone_pairs, unpaired))
        }
        (None, None) => {
            let unpaired = nonbonding % 2;
            let lone_pairs = (nonbonding - unpaired) / 2;
            if lone_pairs > max_lone_pairs {
                return None;
            }
            Some((lone_pairs, unpaired))
        }
    }
}

fn derive_atom(
    implicit_hydrogens: i64,
    unpaired: i64,
    multiplicity: i64,
    lone_pairs: i64,
    valence: i64,
    is_aromatic: bool,
    aromatic_valence: i64,
) -> AtomAst {
    AtomAst {
        implicit_hydrogens: ValueAst::Lit(implicit_hydrogens),
        lone_pairs: ValueAst::Lit(lone_pairs),
        spin: SpinStateAst {
            unpaired: ValueAst::Lit(unpaired),
            multiplicity: ValueAst::Lit(multiplicity),
        },
        constraints: AtomConstraintsAst::from_iter([
            AtomConstraintAst::Valence(ValueAst::Lit(valence)),
            AtomConstraintAst::AromaticValence(if is_aromatic {
                AromaticValenceAst::Aromatic(ValueAst::Lit(aromatic_valence))
            } else {
                AromaticValenceAst::NotAromatic
            }),
        ]),
        ..Default::default()
    }
}

fn derive_multiplicity(spin: &SpinStateAst, unpaired: i64) -> Option<i64> {
    let unpaired = unpaired as u8;
    match spin.multiplicity {
        ValueAst::Lit(m) => {
            let multiplicity = SpinMultiplicity::try_from(m as u8).ok()?;
            SpinState::are_compatible(unpaired, multiplicity).then_some(m)
        }
        ValueAst::Undetermined => Some(i64::from(unpaired) + 1),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::rstest;
    use umol_ast::ast::AtomFieldChange;
    use umol_ast::{atom_dsl, mol_dsl};

    use super::*;
    use crate::ops::valence::ValenceTable;

    #[rstest]
    fn test_counts_valence_plan() {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let resolver = CountsValence::new(&model);
        let molecule = mol_dsl!(r#"{:atoms ["C#c0"]}"#);
        assert_eq!(
            resolver.plan(&molecule),
            Solution::Determined(vec![
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::IsotopeMass {
                        old: IsotopeMassAst::Undetermined,
                        new: IsotopeMassAst::Natural,
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::ImplicitHydrogens {
                        old: ValueAst::Undetermined,
                        new: ValueAst::Lit(4),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::LonePairs {
                        old: ValueAst::Undetermined,
                        new: ValueAst::Lit(0),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::Spin {
                        old: SpinStateAst::default(),
                        new: SpinStateAst::from((0_u8, 1_u8)),
                    },
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(0)),
                    old: None,
                    new: Some(AtomConstraintAst::valence(0_i64)),
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(0)),
                    old: None,
                    new: Some(AtomConstraintAst::aromatic_valence(
                        AromaticValenceAst::NotAromatic,
                    )),
                },
            ])
        );
    }

    #[rstest]
    fn test_counts_valence_plan_identity() {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let resolver = CountsValence::new(&model);
        let molecule = mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#a!"]}"#);
        assert_eq!(resolver.plan(&molecule), Solution::Determined(Vec::new()));
    }

    #[rstest]
    #[case::later_undetermined_element(mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#))]
    fn test_counts_valence_plan_partial(#[case] molecule: MoleculeAst) {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        assert_eq!(
            CountsValence::new(&model).plan(&molecule),
            Solution::Underdetermined(Vec::new())
        );
    }

    #[rstest]
    #[case::later_atom_contradiction(mol_dsl!(r#"{:atoms ["C#c0" "Fe#c0#h0#a+"]}"#), CountsError::UndeterminedAromaticValence)]
    fn test_counts_valence_plan_error(
        #[case] molecule: MoleculeAst,
        #[case] expected: CountsError,
    ) {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let resolver = CountsValence::new(&model);
        assert_eq!(resolver.plan(&molecule), Solution::Contradictory(expected));
    }

    #[rstest]
    #[case::ethane_carbon(mol_dsl!(r#"{:atoms ["C #c0" "C #c0"] :bonds [[0 1 "1"]]}"#), 0, "C#i=#c0#h3#n0#u0#s#v#a!")]
    #[case::water_oxygen(mol_dsl!(r#"{:atoms ["O #c0" "H #c0" "H #c0"] :bonds [[0 1 "1"] [0 2 "1"]]}"#), 0, "O#i=#c0#h0#n2#u0#s#v2#a!")]
    #[case::benzene_ring(mol_dsl!( r#"{:atoms ["C #c0" "C #c0" "C #c0" "C #c0" "C #c0" "C #c0"] :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"]]}"#), 0, "C#i=#c0#h#n0#u0#s#v2#a")]
    fn test_counts_valence_resolve(
        #[case] mut molecule: MoleculeAst,
        #[case] atom_id: u32,
        #[case] expected: &str,
    ) {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let resolver = CountsValence::new(&model);
        assert_eq!(
            resolver.resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule.atom(AtomId(atom_id)).ast.to_string(), expected);
    }

    #[rstest]
    #[case::later_undetermined_element(mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#))]
    fn test_counts_valence_resolve_partial(#[case] mut molecule: MoleculeAst) {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let original = molecule.clone();
        assert_eq!(
            CountsValence::new(&model).resolve(&mut molecule),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::later_atom_contradiction(mol_dsl!(r#"{:atoms ["C#c0" "Fe#c0#h0#a+"]}"#), CountsError::UndeterminedAromaticValence)]
    fn test_counts_valence_resolve_error(
        #[case] mut molecule: MoleculeAst,
        #[case] expected: CountsError,
    ) {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let resolver = CountsValence::new(&model);
        let original = molecule.clone();
        assert_eq!(
            resolver.resolve(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::methane_h("C#c0#h4", "C#c0#h4#n0#u0#s#v0#a!")]
    #[case::methane_h_inferred("C#c0#h*", "C#c0#h4#n0#u0#s#v0#a!")]
    #[case::ammonia("N#c0#h3", "N#c0#h3#n#u0#s#v0#a!")]
    #[case::water("O#c0#h2", "O#c0#h2#n2#u0#s#v0#a!")]
    #[case::methyl_radical("C#c0#h3", "C#c0#h3#n0#u#s2#v0#a!")]
    #[case::methyl_radical_h_inferred("C#c0#u", "C#c0#h3#n0#u#s2#v0#a!")]
    #[case::methyl_anion("C#c-1#h3", "C#c-#h3#n#u0#s#v0#a!")]
    #[case::hydroxyl_radical("O#c0#h1", "O#c0#h#n2#u#s2#v0#a!")]
    #[case::fluoride("F#c-1#h0", "F#c-#h0#n4#u0#s#v0#a!")]
    #[case::fluorine_atom("F#c0#h0", "F#c0#h0#n3#u#s2#v0#a!")]
    #[case::magnesium_atom("Mg#c0#h0", "Mg#c0#h0#n#u0#s#v0#a!")]
    #[case::ethane_carbon("C#c0#v1", "C#c0#h3#n0#u0#s#v#a!")]
    #[case::methylene_carbon("C#c0#v2", "C#c0#h2#n0#u0#s#v2#a!")]
    #[case::methine_carbon("C#c0#v3", "C#c0#h#n0#u0#s#v3#a!")]
    #[case::amine_nitrogen("N#c0#v1", "N#c0#h2#n#u0#s#v#a!")]
    #[case::alcohol_oxygen("O#c0#v1", "O#c0#h#n2#u0#s#v#a!")]
    #[case::benzene_carbon("C#c0#v2#h1#a+", "C#c0#h#n0#u0#s#v2#a")]
    #[case::benzene_carbon_h_inferred("C#c0#v2#h*#a+", "C#c0#h#n0#u0#s#v2#a")]
    #[case::fused_aromatic_carbon_h_inferred("C#c0#v3#h*#a+", "C#c0#h0#n0#u0#s#v3#a")]
    #[case::aromatic_carbon_unpaired("C#c0#v2#h*#u1#a+", "C#c0#h0#n0#u#s2#v2#a")]
    #[case::pyridine_nitrogen("N#c0#v2#h0#a+", "N#c0#h0#n#u0#s#v2#a")]
    #[case::pyrrole_nitrogen("N#c0#v2#h1#a+", "N#c0#h#n0#u0#s#v2#a2")]
    #[case::furan_oxygen("O#c0#v2#h0#a+", "O#c0#h0#n#u0#s#v2#a2")]
    #[case::furan_oxygen_h_inferred("O#c0#v2#h*#a+", "O#c0#h0#n#u0#s#v2#a2")]
    #[case::borazine_boron("B#c0#v2#h1#a+", "B#c0#h#n0#u0#s#v2#a0")]
    #[case::cyclopentadienyl_carbanion("C#c-1#v2#h1#a+", "C#c-#h#n0#u0#s#v2#a2")]
    #[case::tropylium_carbocation("C#c1#v2#h1#a+", "C#c+#h#n0#u0#s#v2#a0")]
    #[case::iron_out_of_table("Fe#c0#h0", "Fe#c0#h0#n4#u0#s#v0#a!")]
    fn test_counts_valence_resolve_atom(#[case] input: &str, #[case] expected: &str) {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let resolver = CountsValence::new(&model);
        let mut atom = atom_dsl!(input);
        resolver.resolve_atom(&mut atom).unwrap();
        assert_eq!(atom.to_string(), expected);
    }

    #[rstest]
    #[case::undetermined_aromatic_out_of_table(atom_dsl!("Fe#c0#h0#a+"), Err(CountsError::UndeterminedAromaticValence))]
    fn test_counts_valence_resolve_atom_error(
        #[case] mut atom: AtomAst,
        #[case] expected: Result<(), CountsError>,
    ) {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let resolver = CountsValence::new(&model);
        assert_eq!(resolver.resolve_atom(&mut atom), expected);
    }

    #[rstest]
    #[case::methane_conforms("C#i=#c0#h4#n0#u0#s#v0#a!", Solution::Determined(()))]
    #[case::excess_lone_pairs("C#i=#c0#h4#n2#u0#s#v0#a!", Solution::Contradictory(CountsMismatch { element: Element::C, charge: 0, valence: 0, }))]
    #[case::not_ground("C", Solution::Underdetermined(()))]
    fn test_counts_valence_classify_molecule_atom(
        #[case] input: &str,
        #[case] expected: Solution<(), CountsMismatch>,
    ) {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let resolver = CountsValence::new(&model);
        let molecule = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![atom_dsl!(input)],
            ..Default::default()
        });
        assert_eq!(
            resolver.classify_molecule_atom(&molecule, AtomId(0)),
            expected
        );
    }
}
