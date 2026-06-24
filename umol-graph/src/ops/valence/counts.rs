//! Counts valence resolver: per-atom electron conservation. `derive_fields`
//! picks the first table covalence ≥ `v` (targets sorted smallest to largest),
//! splits `covalence − v` between
//! implicit H and aromatic increment, then assigns lone pairs and spin from
//! the nonbonding budget. Literals constrain each step.

use thiserror::Error;
use umol_ast::ast::{
    aromatic_increment, AromaticValenceAst, AsLit, AtomAst, AtomConstraint, AtomConstraints,
    AtomId, IsotopeMassAst, Lattice, MoleculeAst, SpinStateAst, ValueAst,
};
use umol_chem::element::Element;
use umol_chem::spin::{SpinMultiplicity, SpinState};
use umol_utils::solution::Solution;

use crate::ops::model::CountsModel;

#[derive(Clone, Debug)]
pub struct CountsValence<'a> {
    model: &'a CountsModel,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum CountsError {
    #[error("no matching valence state")]
    NoMatch,
    #[error("element out of scope: no valence table entry")]
    InvalidElement,
    #[error("aromatic valence unspecified (#a+): no valence table entry")]
    UndeterminedAromaticValence,
    #[error("undetermined element")]
    UndeterminedElement,
}

/// A ground atom that no valence-table state admits.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
#[error("no valence-table state: element {element}, charge {charge}, valence {valence}")]
pub struct CountsMismatch {
    pub element: Element,
    pub charge: i64,
    pub valence: i64,
}

impl<'a> CountsValence<'a> {
    pub fn new(model: &'a CountsModel) -> Self {
        Self { model }
    }

    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), CountsError> {
        for atom in ast.atoms().iter() {
            if atom.element().as_lit().is_none() {
                return Err(CountsError::UndeterminedElement);
            }
        }

        for id in ast.atoms().ids() {
            self.resolve_molecule_atom(ast, id)?;
        }
        Ok(())
    }

    fn resolve_molecule_atom(
        &self,
        ast: &mut MoleculeAst,
        atom_id: AtomId,
    ) -> Result<(), CountsError> {
        let atom = ast.atom(atom_id);
        if atom.is_ground() {
            return Ok(());
        }
        if atom.element().is_undetermined() {
            return Ok(());
        };
        if atom.charge().is_undetermined() {
            return Ok(());
        };

        let Some(valence) = atom.valence().as_lit() else {
            return Ok(());
        };
        let accepted_pairs = atom.accepted_pairs().as_lit_or(0);
        let is_aromatic = atom.is_in_aromatic_system()
            || atom.neighbors().any(|n| n.bond().constraints().aromatic())
            || atom
                .constraints()
                .aromatic_valence()
                .is_some_and(|a| a.is_aromatic());

        let derived = self.derive_fields(atom.ast, valence, accepted_pairs, is_aromatic)?;
        let mut resolved = atom.ast.meet(&derived).ok_or(CountsError::NoMatch)?;
        if resolved.isotope_mass.is_undetermined() {
            resolved.isotope_mass = IsotopeMassAst::Natural;
        }
        *ast.atom_mut(atom_id).ast = resolved;
        Ok(())
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
        let valence = ast
            .constraints
            .valence()
            .unwrap_or(&ValueAst::Undetermined)
            .as_lit_or(0);
        let accepted_pairs = ast
            .constraints
            .accepted_pairs()
            .unwrap_or(&ValueAst::Undetermined)
            .as_lit_or(0);
        let is_aromatic = ast
            .constraints
            .aromatic_valence()
            .is_some_and(|a| a.is_aromatic());

        let derived = self.derive_fields(ast, valence, accepted_pairs, is_aromatic)?;
        *ast = ast.meet(&derived).ok_or(CountsError::NoMatch)?;
        Ok(())
    }

    /// Read-only conformance check for a resolved atom: `Determined` if some
    /// valence-table state admits it, `Contradictory` if none does,
    /// `Underdetermined` if the atom is not ground. Reuses `derive_fields`.
    /// Unlike `resolve_atom` it does not skip ground atoms — that skip is why a
    /// table-violating ground atom otherwise passes resolution unchecked.
    pub fn conforms_molecule_atom(
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
        let valence = atom.valence().as_lit_or(0);
        let accepted_pairs = atom.accepted_pairs().as_lit_or(0);
        let is_aromatic = atom.is_in_aromatic_system()
            || atom.neighbors().any(|n| n.bond().constraints().aromatic())
            || atom
                .constraints()
                .aromatic_valence()
                .is_some_and(|a| a.is_aromatic());
        match self.derive_fields(atom.ast, valence, accepted_pairs, is_aromatic) {
            Ok(_) => Solution::Determined(()),
            Err(_) => Solution::Contradictory(CountsMismatch {
                element,
                charge,
                valence,
            }),
        }
    }

    fn derive_fields(
        &self,
        atom: &AtomAst,
        valence: i64,
        accepted_pairs: i64,
        is_aromatic: bool,
    ) -> Result<AtomAst, CountsError> {
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

        let aromatic_values = aromatic_valence_values(
            aromatic_constraint,
            is_aromatic,
            entry.map(|e| e.aromatic_valences.as_slice()),
        );

        let mut candidates = Vec::new();
        for implicit_hydrogens in
            implicit_hydrogen_values(&atom.implicit_hydrogens, bonding_budget, entry.is_none())?
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
                    derive_lone_pairs_unpaired(atom, element, nonbonding)
                else {
                    continue;
                };
                let Some(multiplicity) = derive_multiplicity(&atom.spin, unpaired) else {
                    continue;
                };
                let updates = make_updates(
                    implicit_hydrogens,
                    unpaired,
                    multiplicity,
                    lone_pairs,
                    valence,
                    is_aromatic,
                    aromatic_valence,
                );
                if let Some(candidate) = atom.meet(&updates) {
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

fn implicit_hydrogen_values(
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

fn aromatic_valence_values(
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

fn derive_lone_pairs_unpaired(
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

fn make_updates(
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
        constraints: AtomConstraints::from_iter([
            AtomConstraint::Valence(ValueAst::Lit(valence)),
            AtomConstraint::AromaticValence(if is_aromatic {
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
            let multiplicity = SpinMultiplicity::from_repr(m as u8)?;
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
    use umol_ast::ast::Constraints;
    use umol_ast::{atom, mol};

    use super::*;
    use crate::ops::valence::ValenceTable;

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
        let mut atom = atom!(input);
        resolver.resolve_atom(&mut atom).unwrap();
        assert_eq!(atom.to_string(), expected);
    }

    #[rstest]
    #[case::ethane_carbon(
        mol!(r#"{:atoms ["C #c0" "C #c0"] :bonds [[0 1 "1"]]}"#),
        0,
        "C#i=#c0#h3#n0#u0#s#v#a!"
    )]
    #[case::water_oxygen(
        mol!(r#"{:atoms ["O #c0" "H #c0" "H #c0"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        0,
        "O#i=#c0#h0#n2#u0#s#v2#a!"
    )]
    #[case::benzene_ring(
        mol!(
            r#"{:atoms ["C #c0" "C #c0" "C #c0" "C #c0" "C #c0" "C #c0"]
               :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"]]}"#
        ),
        0,
        "C#i=#c0#h#n0#u0#s#v2#a"
    )]
    fn test_counts_valence_resolve_molecule_atom(
        #[case] mut molecule: MoleculeAst,
        #[case] atom_id: u32,
        #[case] expected: &str,
    ) {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let resolver = CountsValence::new(&model);
        resolver.resolve(&mut molecule).unwrap();
        assert_eq!(molecule.atom(AtomId(atom_id)).ast.to_string(), expected);
    }

    #[rstest]
    #[case::undetermined_aromatic_out_of_table(
        atom!("Fe#c0#h0#a+"),
        Err(CountsError::UndeterminedAromaticValence)
    )]
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
    #[case::excess_lone_pairs(
        "C#i=#c0#h4#n2#u0#s#v0#a!",
        Solution::Contradictory(CountsMismatch {
            element: Element::C,
            charge: 0,
            valence: 0,
        })
    )]
    #[case::not_ground("C", Solution::Underdetermined(()))]
    fn test_counts_valence_conforms_atom(
        #[case] input: &str,
        #[case] expected: Solution<(), CountsMismatch>,
    ) {
        let model = CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        };
        let resolver = CountsValence::new(&model);
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
        assert_eq!(
            resolver.conforms_molecule_atom(&molecule, AtomId(0)),
            expected
        );
    }
}
