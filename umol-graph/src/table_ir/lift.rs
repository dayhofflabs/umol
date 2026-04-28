//! TableIR → `umol_ast::MoleculeAst` lift.
//!
//! Implements `IntoAst<MoleculeAst> for &Molecule` (and the per-atom and
//! per-bond analogues), mirroring the DSL → AST conversion pattern in
//! `umol-ast`. Field-by-field copy with no defaulting beyond the literal
//! mapping of `Option<T>` → `Lit(t) | Undetermined`.

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AtomAst, AtomConstraint, AtomIdx, BondAst, BondConstraint, Constraints,
    ElementAst, ImplicitHydrogensAst, IntoAst, IsotopeAst, MoleculeAst, MulticenterBondAst,
    NoncovalentBondAst, NoncovalentBondKind, SpinStateAst, ValueAst,
};
use umol_shared::spin::SpinState;

use crate::table_ir::atom::{Atom as TableAtom, ImplicitHydrogens as TableImplicitH};
use crate::table_ir::bond::{
    Bond as TableBond, BondDonation as TableBondDonation, BondNoncovalent as TableNoncovalent,
    BondOrder as TableBondOrder,
};
use crate::table_ir::Molecule as TableMolecule;

/// Failure modes for the TableIR → MoleculeAst lift.
///
/// Empty for now: every TableIR field maps mechanically to a MoleculeAst
/// equivalent. The variant slot exists so that future strictly-checked lift
/// paths (e.g., asserting bond endpoint validity, rejecting unsupported
/// `BondOrder::Any`) can land without changing the public signature.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum LiftError {}

impl IntoAst<MoleculeAst> for &TableMolecule {
    type Ctx = ();
    type Error = LiftError;

    fn into_ast(self, ctx: &Self::Ctx) -> Result<MoleculeAst, LiftError> {
        let atoms: Vec<AtomAst> = self
            .atoms
            .iter()
            .map(|a| a.into_ast(ctx))
            .collect::<Result<_, _>>()?;

        let mut regular = Vec::new();
        let mut dative: Vec<(Vec<AtomIdx>, AtomIdx, umol_ast::ast::DativeBondAst)> = Vec::new();
        let mut noncovalent = Vec::new();
        for b in &self.bonds {
            let a_idx = AtomIdx(b.atoms.first());
            let b_idx = AtomIdx(b.atoms.second());
            if let Some(kind) = b.noncovalent.map(noncovalent_kind) {
                noncovalent.push((a_idx, b_idx, NoncovalentBondAst::from_kind(kind)));
            } else if let Some(donation) = b.donation {
                let (donor, acceptor) = match donation {
                    TableBondDonation::Donating => (a_idx, b_idx),
                    TableBondDonation::Accepting => (b_idx, a_idx),
                    _ => {
                        regular.push((a_idx, b_idx, b.into_ast(ctx)?));
                        continue;
                    }
                };
                let dative_bond =
                    umol_ast::ast::DativeBondAst::new(lift_bond_order(b.order));
                dative.push((vec![donor], acceptor, dative_bond));
            } else {
                regular.push((a_idx, b_idx, b.into_ast(ctx)?));
            }
        }

        let multicenter: Vec<(Vec<AtomIdx>, MulticenterBondAst)> = self
            .multicenter_bonds
            .iter()
            .map(|mc| {
                let mut seen = std::collections::HashSet::new();
                let atoms: Vec<AtomIdx> = mc
                    .all_atoms()
                    .into_iter()
                    .filter(|a| seen.insert(*a))
                    .map(AtomIdx)
                    .collect();
                let n = atoms.len();
                (
                    atoms,
                    MulticenterBondAst::new(
                        vec![ValueAst::Undetermined; n],
                        ValueAst::Undetermined,
                        SpinStateAst::default(),
                    ),
                )
            })
            .collect();

        // Per-atom constraints from TableIR (aromatic flag, asserted
        // valence) live inline on each `AtomAst.constraints`. The
        // molecule-scope `Constraints` container stays empty here; lift adds
        // entity-scope constraints onto the atoms themselves below.
        let mut atoms = atoms;
        for (i, a) in self.atoms.iter().enumerate() {
            if a.aromatic == Some(true) {
                atoms[i]
                    .constraints
                    .add(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(
                        ValueAst::Undetermined,
                    )));
            }
            if let Some(v) = a.valence {
                atoms[i]
                    .constraints
                    .add(AtomConstraint::Valence(ValueAst::Lit(v as i64)));
            }
        }
        let constraints = Constraints::new();

        Ok(MoleculeAst::new(
            atoms,
            regular,
            dative,
            vec![],
            multicenter,
            noncovalent,
            constraints,
        ))
    }
}

impl IntoAst<AtomAst> for &TableAtom {
    type Ctx = ();
    type Error = LiftError;

    fn into_ast(self, _ctx: &Self::Ctx) -> Result<AtomAst, LiftError> {
        Ok(AtomAst {
            element: ElementAst::Lit(self.element),
            isotope_mass: match self.isotope_mass {
                Some(m) => IsotopeAst::Lit(m as i64),
                None => IsotopeAst::Undetermined,
            },
            charge: match self.charge {
                Some(c) => ValueAst::Lit(c as i64),
                None => ValueAst::Undetermined,
            },
            implicit_hydrogens: match &self.implicit_hydrogens {
                Some(TableImplicitH::Normal) => ImplicitHydrogensAst::Normal,
                Some(TableImplicitH::Hydrogens(n)) => ImplicitHydrogensAst::Lit(*n as i64),
                None => ImplicitHydrogensAst::Undetermined,
            },
            lone_pairs: match self.lone_pairs {
                Some(n) => ValueAst::Lit(n as i64),
                None => ValueAst::Undetermined,
            },
            spin: lift_atom_spin(self),
            constraints: Default::default(),
        })
    }
}

impl IntoAst<BondAst> for &TableBond {
    type Ctx = ();
    type Error = LiftError;

    fn into_ast(self, _ctx: &Self::Ctx) -> Result<BondAst, LiftError> {
        let mut bond = BondAst::new(lift_bond_order(self.order));
        bond.charge = match self.charge {
            Some(c) => ValueAst::Lit(c as i64),
            None => ValueAst::Undetermined,
        };
        bond.spin = lift_bond_spin(self);
        if matches!(self.order, TableBondOrder::Aromatic) {
            bond.constraints.add(BondConstraint::Aromatic);
        }
        Ok(bond)
    }
}

fn lift_atom_spin(atom: &TableAtom) -> SpinStateAst {
    match (atom.unpaired_electrons, atom.multiplicity) {
        (Some(u), Some(m)) => match SpinState::try_new(u, m) {
            Ok(s) => SpinStateAst::from_state(s),
            Err(_) => SpinStateAst::from_values(
                ValueAst::Lit(u as i64),
                ValueAst::Lit(m.multiplicity() as i64),
            ),
        },
        (Some(u), None) => {
            SpinStateAst::from_values(ValueAst::Lit(u as i64), ValueAst::Undetermined)
        }
        (None, Some(m)) => SpinStateAst::from_values(
            ValueAst::Undetermined,
            ValueAst::Lit(m.multiplicity() as i64),
        ),
        (None, None) => SpinStateAst::default(),
    }
}

fn lift_bond_spin(bond: &TableBond) -> SpinStateAst {
    match (bond.unpaired_electrons, bond.multiplicity) {
        (Some(u), Some(m)) => match SpinState::try_new(u, m) {
            Ok(s) => SpinStateAst::from_state(s),
            Err(_) => SpinStateAst::from_values(
                ValueAst::Lit(u as i64),
                ValueAst::Lit(m.multiplicity() as i64),
            ),
        },
        (Some(u), None) => {
            SpinStateAst::from_values(ValueAst::Lit(u as i64), ValueAst::Undetermined)
        }
        (None, Some(m)) => SpinStateAst::from_values(
            ValueAst::Undetermined,
            ValueAst::Lit(m.multiplicity() as i64),
        ),
        (None, None) => SpinStateAst::default(),
    }
}

fn lift_bond_order(order: TableBondOrder) -> ValueAst {
    match order {
        TableBondOrder::Zero => ValueAst::Lit(0),
        TableBondOrder::Single => ValueAst::Lit(1),
        TableBondOrder::Double => ValueAst::Lit(2),
        TableBondOrder::Triple => ValueAst::Lit(3),
        TableBondOrder::Quadruple => ValueAst::Lit(4),
        TableBondOrder::Quintuple => ValueAst::Lit(5),
        TableBondOrder::Sextuple => ValueAst::Lit(6),
        // Definite-aromatic: σ-order is 1 by Kekulé convention; the aromatic
        // flag is added separately as `BondConstraint::Aromatic`. Renders as
        // `1#a`.
        TableBondOrder::Aromatic => ValueAst::Lit(1),
        // Fuzzy orders: no concrete σ-order can be assigned; lift to
        // `Undetermined`. Aromatic-flag setting (where applicable) is left
        // off — the chemistry of these is too ambiguous for the lift.
        TableBondOrder::SingleOrDouble
        | TableBondOrder::SingleOrAromatic
        | TableBondOrder::DoubleOrAromatic
        | TableBondOrder::Any => ValueAst::Undetermined,
    }
}

fn noncovalent_kind(kind: TableNoncovalent) -> NoncovalentBondKind {
    match kind {
        TableNoncovalent::Hydrogen => NoncovalentBondKind::HydrogenBond,
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::table_ir::atom::Atom as TableAtom;
    use crate::table_ir::bond::{Bond as TableBond, BondOrder as TableBondOrder};
    use crate::table_ir::Molecule as TableMolecule;

    fn methane_table() -> TableMolecule {
        let mut atom = TableAtom::from_element(Element::C);
        atom.implicit_hydrogens = Some(TableImplicitH::Hydrogens(4));
        let mut mol = TableMolecule::empty();
        mol.atoms.push(atom);
        mol
    }

    #[rstest]
    fn test_table_molecule_into_ast_methane() {
        let mol = methane_table();
        let ast: MoleculeAst = (&mol).into_ast(&()).unwrap();
        assert_eq!(ast.atoms().count(), 1);
        let atom = ast.atom(AtomIdx(0)).data;
        assert_eq!(atom.element, ElementAst::Lit(Element::C));
        assert!(matches!(atom.implicit_hydrogens, ImplicitHydrogensAst::Lit(4)));
    }

    #[rstest]
    fn test_table_molecule_into_ast_bond_order() {
        let mut mol = TableMolecule::empty();
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.bonds.push(TableBond::new(0, 1, TableBondOrder::Double));
        let ast: MoleculeAst = (&mol).into_ast(&()).unwrap();
        let bond = ast.bond(umol_ast::ast::BondIdx(0)).data;
        assert!(matches!(bond.order, ValueAst::Lit(2)));
    }

    #[rstest]
    fn test_table_molecule_into_ast_aromatic_bond() {
        let mut mol = TableMolecule::empty();
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.bonds
            .push(TableBond::new(0, 1, TableBondOrder::Aromatic));
        let ast: MoleculeAst = (&mol).into_ast(&()).unwrap();
        let bond = ast.bond(umol_ast::ast::BondIdx(0)).data;
        // Definite-aromatic lifts to Kekulé σ-order 1 plus the Aromatic
        // constraint, rendering as `1#a`.
        assert!(matches!(bond.order, ValueAst::Lit(1)));
        assert!(bond
            .constraints
            .iter()
            .any(|c| matches!(c, BondConstraint::Aromatic)));
    }
}
