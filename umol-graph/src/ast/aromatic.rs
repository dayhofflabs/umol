//! Aromatic system AST container.
//!
//! `AromaticSystem` packages a perceived delocalized π system for lifting into
//! a `MoleculeAst`: the participant indices, the per-system base attributes
//! (`AromaticSystemAst`), and the per-participant constraints (typically
//! `AromaticValence` for atoms and `Aromatic` for bonds).

use std::collections::HashSet;

use umol_graph_core::NodeId;

use crate::ast::constraint::{AtomConstraint, BondConstraint};
use crate::ast::molecule::{AromaticSystemAst, MoleculeAst};
use crate::ast::{AtomIdx, BondIdx};

/// A perceived delocalized π system.
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
                AtomConstraint::AromaticValence(
                    crate::ast::constraint::AromaticValenceConstraint::Value(
                        umol_shared::value_ast::ValueAst::Lit(n),
                    ),
                ) => *n as u8,
                _ => 0,
            })
            .sum()
    }

    pub fn contains_atom(&self, atom: AtomIdx) -> bool {
        self.atoms.binary_search(&atom).is_ok()
    }
}

/// Bonds of the induced subgraph over `atoms`: edges whose both endpoints
/// are in the set. Result is sorted by `BondIdx`.
pub fn induced_bonds(ast: &MoleculeAst, atoms: &[AtomIdx]) -> Vec<BondIdx> {
    let set: HashSet<AtomIdx> = atoms.iter().copied().collect();
    let graph = ast.graph();
    let mut seen: HashSet<BondIdx> = HashSet::new();
    let mut bonds: Vec<BondIdx> = Vec::new();
    for &a in atoms {
        for n in graph.neighbors(NodeId::from(a)) {
            let other = AtomIdx::from(n.node);
            if set.contains(&other) {
                let bond = BondIdx::from(n.edge);
                if seen.insert(bond) {
                    bonds.push(bond);
                }
            }
        }
    }
    bonds.sort_unstable();
    bonds
}
