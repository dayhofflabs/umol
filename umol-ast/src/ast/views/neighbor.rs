//! Neighbor view: yielded by `MoleculeAst::neighbors`.

use super::super::idx::{AtomId, BondId};
use super::super::molecule::MoleculeAst;
use super::atom::AtomView;
use super::bond::BondView;

/// Neighbor-side view of a bond: the atom on the other end and the bond
/// connecting it. Full `AtomView` / `BondView` lookups are deferred to
/// `atom()` / `bond()`.
#[derive(Clone, Copy, Debug)]
pub struct NeighborView<'a> {
    atom_id: AtomId,
    bond_id: BondId,
    molecule: &'a MoleculeAst,
}

impl<'a> NeighborView<'a> {
    pub(crate) fn new(atom_id: AtomId, bond_id: BondId, molecule: &'a MoleculeAst) -> Self {
        Self {
            atom_id,
            bond_id,
            molecule,
        }
    }

    pub fn atom_id(&self) -> AtomId {
        self.atom_id
    }

    pub fn bond_id(&self) -> BondId {
        self.bond_id
    }

    pub fn atom(&self) -> AtomView<'a> {
        self.molecule.atom(self.atom_id)
    }

    pub fn bond(&self) -> BondView<'a> {
        self.molecule.bond(self.bond_id)
    }
}
