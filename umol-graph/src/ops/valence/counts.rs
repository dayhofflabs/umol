//! Counts valence resolver: per-atom electron conservation. `derive_fields`
//! builds the resolved atom type (implicit H, lone pairs, spin, aromatic
//! valence) from the per-atom inputs; `resolve_atom` and `resolve_molecule_atom`
//! gather those inputs — from constraints and from topology respectively — and
//! apply the synthesized type with `meet`. Implicit hydrogens saturate to the
//! lowest admissible `covalence_set` entry, the aromatic valence is computed
//! from the covalence, and spin is closed-shell (fewest unpaired, then most
//! lone pairs).

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AsLit, AtomAst, AtomConstraint, AtomConstraints, AtomId, Lattice,
    MoleculeAst, ValueAst,
};

use crate::ops::valence::table::ValenceTable;

#[derive(Clone, Debug)]
pub struct CountsValence {
    pub table: ValenceTable,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum CountsError {
    #[error("no matching valence state")]
    NoMatch,
    #[error("element out of scope: no covalence set")]
    InvalidElement,
    #[error("undetermined element")]
    UndeterminedElement,
    #[error("aromatic valence requires implicit hydrogens")]
    AromaticNoImplicitHydrogens,
}

impl CountsValence {
    pub fn new(table: ValenceTable) -> Self {
        Self { table }
    }

    /// Resolve all atoms in molecule. Errors if element is undetermined or not in valence table.
    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), CountsError> {
        for atom in ast.atoms().iter() {
            let Some(element) = atom.element().as_lit() else {
                return Err(CountsError::UndeterminedElement);
            };
            if self.table.entry(element).is_none() {
                return Err(CountsError::InvalidElement);
            }
        }

        for i in 0..ast.atoms().count() as u32 {
            self.resolve_molecule_atom(ast, AtomId(i))?;
        }
        Ok(())
    }

    /// Resolve atom in molecule: valence and accepted pairs are read from topology.
    /// Atom is considered aromatic if in an aromatic system *or* has aromatic-valence constraint.
    fn resolve_molecule_atom(
        &self,
        ast: &mut MoleculeAst,
        atom_id: AtomId,
    ) -> Result<(), CountsError> {
        let atom = ast.atom(atom_id);
        if atom.ast.is_ground() {
            return Ok(());
        }
        if atom.element().is_undetermined() {
            return Ok(());
        };
        if atom.charge().is_undetermined() {
            return Ok(());
        };

        // Defer if the localized valence is not determinate yet.
        let Some(valence) = atom.valence().as_lit() else {
            return Ok(());
        };
        let accepted_pairs = atom.accepted_pairs().as_lit().unwrap_or(0);
        let is_aromatic = atom.is_in_aromatic_system()
            || AromaticValenceAst::aromatic(ValueAst::Undetermined)
                .matches(&atom.constraints().aromatic_valence());

        let derived = self.derive_fields(atom.ast, valence, accepted_pairs, is_aromatic)?;
        let resolved = atom.ast.meet(&derived).ok_or(CountsError::NoMatch)?;
        *ast.atom_mut(atom_id).ast = resolved;
        Ok(())
    }

    /// Resolve atom: valence, accepted pairs, aromatic valence are read from constraints
    /// (absent constraints default to zero / non-aromatic).
    /// Errors if element is undetermined or not in valence table or if charge is undetermined.
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
        let valence = ast.constraints.valence().as_lit().unwrap_or(0);
        let accepted_pairs = ast.constraints.accepted_pairs().as_lit().unwrap_or(0);
        let is_aromatic = AromaticValenceAst::aromatic(ValueAst::Undetermined)
            .matches(&ast.constraints.aromatic_valence());

        let derived = self.derive_fields(ast, valence, accepted_pairs, is_aromatic)?;
        *ast = ast.meet(&derived).ok_or(CountsError::NoMatch)?;
        Ok(())
    }

    /// Build atom type updates (implicit hydrogens, lone pairs, spin, and constraints).
    /// Sets derived valence and aromatic valence constraints.
    /// Returns `None` when the atom is underdetermined. Caller should call meet with the derived values.
    fn derive_fields(
        &self,
        atom: &AtomAst,
        valence: i64,
        accepted_pairs: i64,
        is_aromatic: bool,
    ) -> Result<AtomAst, CountsError> {
        let element = atom.element.as_lit().unwrap();
        let charge = atom.charge.as_lit().unwrap();

        // covalence_set of the charge/dative-shifted (isoelectronic) element.
        // Error if shift target or table data is missing.
        let covalence_set = match element
            .shift((2 * accepted_pairs - charge) as i8)
            .and_then(|shifted| self.table.entry(shifted))
        {
            Some(entry) if !entry.covalence_set.is_empty() => &entry.covalence_set,
            _ => return Err(CountsError::InvalidElement),
        };

        let mut derived = AtomAst::default();

        let implicit_hydrogens = match atom.implicit_hydrogens.as_lit() {
            Some(v) => v,
            None => {
                if is_aromatic {
                    return Err(CountsError::AromaticNoImplicitHydrogens);
                }
                let v = derive_implicit_hydrogens(covalence_set, valence)?;
                derived.implicit_hydrogens = ValueAst::lit(v);
                v
            }
        };

        let electrons = i64::from(element.valence_electrons());
        let valence_capacity = i64::from(element.valence_capacity());

        let aromatic_valence = derive_aromatic_valence(
            covalence_set,
            electrons,
            charge,
            implicit_hydrogens,
            valence,
            is_aromatic,
        )?;

        let unpaired = match (atom.lone_pairs.as_lit(), atom.spin.unpaired.as_lit()) {
            (Some(_), Some(u)) => u,
            (input_lone_pairs, input_unpaired) => {
                let (computed_lone_pairs, computed_unpaired) = derive_lone_pairs_unpaired(
                    valence_capacity,
                    electrons,
                    charge,
                    implicit_hydrogens,
                    valence,
                    aromatic_valence,
                )?;

                if input_lone_pairs.is_none() {
                    derived.lone_pairs = ValueAst::lit(computed_lone_pairs);
                }
                if input_unpaired.is_none() {
                    derived.spin.unpaired = ValueAst::lit(computed_unpaired);
                }

                input_unpaired.unwrap_or(computed_unpaired)
            }
        };

        match atom.spin.multiplicity.as_lit() {
            Some(_) => {}
            None => {
                derived.spin.multiplicity = ValueAst::lit(unpaired + 1);
            }
        }

        derived.constraints = AtomConstraints::from_iter([
            AtomConstraint::Valence(ValueAst::Lit(valence)),
            AtomConstraint::AromaticValence(if is_aromatic {
                AromaticValenceAst::Aromatic(ValueAst::Lit(aromatic_valence))
            } else {
                AromaticValenceAst::NotAromatic
            }),
        ]);
        Ok(derived)
    }
}

// Derive implicit hydrogens from lowest covalence exceeding derived valence
#[inline]
fn derive_implicit_hydrogens(covalence_set: &[u8], valence: i64) -> Result<i64, CountsError> {
    covalence_set
        .iter()
        .find_map(|&c| {
            let i = i64::from(c) - valence;
            if i >= 0 {
                Some(i)
            } else {
                None
            }
        })
        .ok_or(CountsError::NoMatch)
}

// Derive aromatic valence from the covalence via aromatic increment.
// `increment = covalence − (v+h)`.
// increment = 1 -> aromatic valence = 1
// increment = 0 -> aromatic valence = 0 or 2
// Distinguish aromatic valence = 0 and 2 from remaining electron count in aromatic system.
// Aromatic valence = 0 in acceptors (if total_localized_valence = covalence)
// Aromatic valence = 2 in donors (if total localized valence < covalence)
// Aromatic valence = 1 in other aromatic atoms
#[inline]
fn derive_aromatic_valence(
    covalence_set: &[u8],
    electrons: i64,
    charge: i64,
    implicit_hydrogens: i64,
    valence: i64,
    is_aromatic: bool,
) -> Result<i64, CountsError> {
    let total_localized_valence = valence + implicit_hydrogens;
    if is_aromatic {
        if covalence_set.contains(&(total_localized_valence as u8)) {
            match electrons - charge - total_localized_valence {
                0 => Ok(0),
                _ => Ok(2),
            }
        } else if covalence_set.contains(&((total_localized_valence + 1) as u8)) {
            Ok(1)
        } else {
            Err(CountsError::NoMatch)
        }
    } else {
        Ok(0)
    }
}

// Derive nonbonding electrons from per-atom conservation: e − c = v + h + a + u + 2n.
// Spit nonbonding electrons by max-n principle
#[inline]
fn derive_lone_pairs_unpaired(
    valence_capacity: i64,
    electrons: i64,
    charge: i64,
    implicit_hydrogens: i64,
    valence: i64,
    aromatic_valence: i64,
) -> Result<(i64, i64), CountsError> {
    let nonbonding = electrons - charge - implicit_hydrogens - valence - aromatic_valence;
    if nonbonding < 0 {
        return Err(CountsError::NoMatch);
    }
    let unpaired = nonbonding % 2;
    if unpaired > nonbonding || (nonbonding - unpaired) % 2 != 0 {
        return Err(CountsError::NoMatch);
    }
    let lone_pairs = (nonbonding - unpaired) / 2;
    if lone_pairs > valence_capacity / 2 {
        Err(CountsError::NoMatch)
    } else {
        Ok((lone_pairs, unpaired))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::{atom, mol};

    use super::*;
    use crate::valence_table;

    #[rstest]
    #[case::methane_h("C#c0#h4", "C#c0#h4#n0#u0#s#v0#a!")]
    #[case::methane_h_inferred("C#c0#h*", "C#c0#h4#n0#u0#s#v0#a!")]
    #[case::ammonia("N#c0#h3", "N#c0#h3#n#u0#s#v0#a!")]
    #[case::water("O#c0#h2", "O#c0#h2#n2#u0#s#v0#a!")]
    #[case::methyl_radical("C#c0#h3", "C#c0#h3#n0#u#s2#v0#a!")]
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
    #[case::pyridine_nitrogen("N#c0#v2#h0#a+", "N#c0#h0#n#u0#s#v2#a")]
    #[case::pyrrole_nitrogen("N#c0#v2#h1#a+", "N#c0#h#n0#u0#s#v2#a2")]
    #[case::furan_oxygen("O#c0#v2#h0#a+", "O#c0#h0#n#u0#s#v2#a2")]
    #[case::borazine_boron("B#c0#v2#h1#a+", "B#c0#h#n0#u0#s#v2#a0")]
    #[case::cyclopentadienyl_carbanion("C#c-1#v2#h1#a+", "C#c-#h#n0#u0#s#v2#a2")]
    #[case::tropylium_carbocation("C#c1#v2#h1#a+", "C#c+#h#n0#u0#s#v2#a0")]
    fn test_counts_valence_resolve_atom(#[case] input: &str, #[case] expected: &str) {
        let resolver = CountsValence::new(ValenceTable::default_table().clone());
        let mut atom = atom!(input);
        resolver.resolve_atom(&mut atom).unwrap();
        assert_eq!(atom.to_string(), expected);
    }

    #[rstest]
    #[case::ethane_carbon(
        mol!(r#"{:atoms ["C #c0" "C #c0"] :bonds [[0 1 "1"]]}"#),
        0,
        "C#c0#h3#n0#u0#s#v#a!"
    )]
    #[case::water_oxygen(
        mol!(r#"{:atoms ["O #c0" "H #c0" "H #c0"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        0,
        "O#c0#h0#n2#u0#s#v2#a!"
    )]
    fn test_counts_valence_resolve_molecule_atom(
        #[case] mut molecule: MoleculeAst,
        #[case] atom_id: u32,
        #[case] expected: &str,
    ) {
        let resolver = CountsValence::new(ValenceTable::default_table().clone());
        resolver.resolve(&mut molecule).unwrap();
        assert_eq!(molecule.atom(AtomId(atom_id)).ast.to_string(), expected);
    }

    #[test]
    fn test_counts_valence_resolve_atom_error() {
        let resolver = CountsValence::new(valence_table! { C => [] });
        let mut atom = atom!("C#c0");
        assert_eq!(
            resolver.resolve_atom(&mut atom),
            Err(CountsError::InvalidElement)
        );
    }
}
