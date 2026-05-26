//! Counts valence resolver: per-atom electron conservation. `synthesize`
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
    MoleculeAst, SpinStateAst, ValueAst,
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
}

impl CountsValence {
    pub fn new(table: ValenceTable) -> Self {
        Self { table }
    }

    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), CountsError> {
        // Every atom must name an element the table covers before resolving.
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

    /// Resolve an atom in a molecule: localized valence and accepted dative
    /// pairs come from topology, and the atom counts as aromatic if it is in an
    /// aromatic system *or* carries an aromatic-valence constraint. The
    /// synthesized type is applied with `meet` (a conflict with a pinned field
    /// is a `NoMatch`).
    fn resolve_molecule_atom(
        &self,
        ast: &mut MoleculeAst,
        atom_id: AtomId,
    ) -> Result<(), CountsError> {
        let atom = ast.atom(atom_id);
        if atom.ast.is_ground() {
            return Ok(());
        }
        // Defer if the localized valence is not determinate yet.
        let Some(valence) = atom.valence().as_lit() else {
            return Ok(());
        };
        let accepted_pairs = atom.accepted_pairs().as_lit().unwrap_or(0);
        let is_aromatic = atom.is_in_aromatic_system()
            || AromaticValenceAst::aromatic(ValueAst::Undetermined)
                .matches(&atom.constraints().aromatic_valence());

        let Some(synthesized) = self.synthesize(atom.ast, valence, accepted_pairs, is_aromatic)?
        else {
            return Ok(());
        };
        let resolved = atom.ast.meet(&synthesized).ok_or(CountsError::NoMatch)?;
        *ast.atom_mut(atom_id).ast = resolved;
        Ok(())
    }

    /// Resolve a standalone atom: localized valence (`#v`), accepted dative
    /// pairs (`#t`), and aromatic membership (`#a`) are read from the atom's
    /// own constraints (absent ones default to zero / non-aromatic). The
    /// synthesized type is applied with `meet`.
    pub fn resolve_atom(&self, ast: &mut AtomAst) -> Result<(), CountsError> {
        if ast.is_ground() {
            return Ok(());
        }
        let valence = ast.constraints.valence().as_lit().unwrap_or(0);
        let accepted_pairs = ast.constraints.accepted_pairs().as_lit().unwrap_or(0);
        let is_aromatic = AromaticValenceAst::aromatic(ValueAst::Undetermined)
            .matches(&ast.constraints.aromatic_valence());

        if let Some(synthesized) = self.synthesize(ast, valence, accepted_pairs, is_aromatic)? {
            *ast = ast.meet(&synthesized).ok_or(CountsError::NoMatch)?;
        }
        Ok(())
    }

    /// Build the resolved atom type — implicit hydrogens, lone pairs, spin, and
    /// aromatic valence — from the atom's element/charge/spin and the supplied
    /// localized valence, accepted dative pairs, and aromaticity. Returns
    /// `None` when the atom is underdetermined (no `Lit` element or charge, or
    /// aromatic with an unknown H count). The returned atom sets only the
    /// computed fields; the caller `meet`s it onto the target.
    fn synthesize(
        &self,
        atom: &AtomAst,
        valence: i64,
        accepted_pairs: i64,
        is_aromatic: bool,
    ) -> Result<Option<AtomAst>, CountsError> {
        let Some(element) = atom.element.as_lit() else {
            return Ok(None);
        };
        // Charge must be given; it is never defaulted to zero.
        let Some(charge) = atom.charge.as_lit() else {
            return Ok(None);
        };

        // covalence_set of the charge/dative-shifted (isoelectronic) element.
        // No shift target, no table entry, or an empty set: out of scope.
        let covalence_set = match element
            .shift((2 * accepted_pairs - charge) as i8)
            .and_then(|shifted| self.table.entry(shifted))
        {
            Some(entry) if !entry.covalence_set.is_empty() => &entry.covalence_set,
            _ => return Err(CountsError::InvalidElement),
        };

        // Implicit hydrogens: given, else saturate to the lowest covalence ≥ v.
        let implicit_hydrogens = match atom.implicit_hydrogens.as_lit() {
            Some(h) => h,
            None if is_aromatic => return Ok(None), // need a known H count to determine `a`
            None => match covalence_set.iter().filter(|&&c| i64::from(c) >= valence).min() {
                Some(&c) => i64::from(c) - valence,
                None => return Err(CountsError::NoMatch),
            },
        };
        let total_localized_valence = valence + implicit_hydrogens;
        let electrons = i64::from(element.valence_electrons());

        // Aromatic valence `a` from the covalence: `ai = covalence − (v+h)`.
        // `ai = 1` ⇒ `a = 1`; `ai = 0` splits acceptor (no spare electrons)
        // from lone-pair donor by the leftover electron count.
        let aromatic_valence = if is_aromatic {
            if covalence_set.contains(&(total_localized_valence as u8)) {
                if electrons - charge - total_localized_valence == 0 {
                    0
                } else {
                    2
                }
            } else if covalence_set.contains(&((total_localized_valence + 1) as u8)) {
                1
            } else {
                return Err(CountsError::NoMatch);
            }
        } else {
            0
        };

        // Per-atom conservation: e − c = v + h + a + u + 2n.
        let nonbonding = electrons - charge - total_localized_valence - aromatic_valence;
        if nonbonding < 0 {
            return Err(CountsError::NoMatch);
        }
        let unpaired = atom.spin.unpaired.as_lit().unwrap_or(nonbonding % 2);
        if unpaired > nonbonding || (nonbonding - unpaired) % 2 != 0 {
            return Err(CountsError::NoMatch);
        }
        let lone_pairs = (nonbonding - unpaired) / 2;
        if lone_pairs > i64::from(element.valence_capacity() / 2) {
            return Err(CountsError::NoMatch);
        }
        let multiplicity = atom.spin.multiplicity.as_lit().unwrap_or(unpaired + 1);

        let constraints = AtomConstraints::from_iter([
            AtomConstraint::Valence(ValueAst::Lit(valence)),
            AtomConstraint::AromaticValence(if is_aromatic {
                AromaticValenceAst::Aromatic(ValueAst::Lit(aromatic_valence))
            } else {
                AromaticValenceAst::NotAromatic
            }),
        ]);
        Ok(Some(AtomAst {
            implicit_hydrogens: ValueAst::Lit(implicit_hydrogens),
            lone_pairs: ValueAst::Lit(lone_pairs),
            spin: SpinStateAst {
                unpaired: ValueAst::Lit(unpaired),
                multiplicity: ValueAst::Lit(multiplicity),
            },
            constraints,
            ..AtomAst::default()
        }))
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
