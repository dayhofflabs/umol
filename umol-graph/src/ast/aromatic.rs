//! Aromatic system graph container.

use umol_ast::ast::aromatic::AromaticSystemAst;
use umol_ast::ast::value::ValueAst;

use super::bond::{coerce_spin, release_spin};
use super::config::{BondAstConfig, NumericMode};
use super::constraint::{AromaticValenceConstraint, AtomConstraint, BondConstraint};
use super::{AtomIdx, BondIdx};

pub fn coerce_aromatic_system(ast: &mut AromaticSystemAst, cfg: &BondAstConfig) {
    if matches!(ast.charge, ValueAst::Undetermined) {
        ast.charge = match cfg.charge_mode {
            NumericMode::Zero => ValueAst::Lit(0),
            NumericMode::Required => ValueAst::Undetermined,
        };
    }
    coerce_spin(&mut ast.spin, cfg);
}

pub fn release_aromatic_system(ast: &mut AromaticSystemAst, cfg: &BondAstConfig) {
    if matches!(
        (&cfg.charge_mode, &ast.charge),
        (NumericMode::Zero, ValueAst::Lit(0))
    ) {
        ast.charge = ValueAst::Undetermined;
    }
    release_spin(&mut ast.spin, cfg);
}

/// Aromatic system ground term
///
/// Container for the per-system base attributes (charge, spin, electrons) plus
/// the per-participant constraints that lift into the molecule. `atom_constraints`
/// is parallel to `atoms` (one constraint per participant, typically
/// `AromaticValence`); `bond_constraints` is parallel to `bonds` (typically
/// `Aromatic`). Atoms are stored sorted and deduplicated.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AromaticSystem {
    atoms: Vec<AtomIdx>,
    bonds: Vec<BondIdx>,
    ast: AromaticSystemAst,
    atom_constraints: Vec<AtomConstraint>,
    bond_constraints: Vec<BondConstraint>,
}

impl AromaticSystem {
    pub fn new(
        atoms: Vec<AtomIdx>,
        bonds: Vec<BondIdx>,
        ast: AromaticSystemAst,
        atom_constraints: Vec<AtomConstraint>,
        bond_constraints: Vec<BondConstraint>,
    ) -> Self {
        debug_assert_eq!(atoms.len(), atom_constraints.len());
        debug_assert_eq!(bonds.len(), bond_constraints.len());
        Self {
            atoms,
            bonds,
            ast,
            atom_constraints,
            bond_constraints,
        }
    }

    pub fn atoms(&self) -> &[AtomIdx] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[BondIdx] {
        &self.bonds
    }

    pub fn ast(&self) -> &AromaticSystemAst {
        &self.ast
    }

    pub fn atom_constraints(&self) -> &[AtomConstraint] {
        &self.atom_constraints
    }

    pub fn bond_constraints(&self) -> &[BondConstraint] {
        &self.bond_constraints
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Sum of π-electron contributions from the attached `AromaticValence`
    /// constraints. Non-`Value(Lit)` entries contribute 0.
    pub fn electron_count(&self) -> u8 {
        self.atom_constraints
            .iter()
            .map(|c| match c {
                AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(
                    ValueAst::Lit(n),
                )) => *n as u8,
                _ => 0,
            })
            .sum()
    }

    pub fn contains_atom(&self, atom: AtomIdx) -> bool {
        self.atoms.binary_search(&atom).is_ok()
    }

    pub fn contains_bond(&self, bond: BondIdx) -> bool {
        self.bonds.binary_search(&bond).is_ok()
    }
}
