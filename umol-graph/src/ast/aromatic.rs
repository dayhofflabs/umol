//! Aromatic system AST structures.

use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use super::constraint::{AromaticValenceConstraint, AtomConstraint, BondConstraint};
use super::{AtomIdx, BondIdx};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticSystemAst {
    pub charge: ValueAst,
    pub spin: SpinStateAst,
}

impl AromaticSystemAst {
    pub fn is_ground(&self) -> bool {
        self.charge.is_ground() && self.spin.is_ground()
    }
}

/// Aromatic system ground term
///
/// Container for the per-system base attributes (charge, spin) plus the
/// per-participant constraints that lift into the molecule. `atom_constraints`
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
