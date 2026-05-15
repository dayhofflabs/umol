//! Neighbor view: yielded by `MoleculeAst::neighbors`.

use super::super::bond::BondAst;
use super::super::idx::{AtomId, BondId};
use super::super::molecule::MoleculeAst;

/// Neighbor-side view of a bond: the atom on the other end (`atom`), the
/// bond index, the bond data, and the parent `MoleculeAst` for navigation
/// to the neighbor's full atom view. Yielded by `MoleculeAst::neighbors`.
#[derive(Clone, Copy, Debug)]
pub struct NeighborView<'a> {
    pub bond: BondId,
    pub atom: AtomId,
    pub ast: &'a BondAst,
    #[allow(dead_code)]
    molecule: &'a MoleculeAst,
}

impl<'a> NeighborView<'a> {
    pub(crate) fn new(
        bond: BondId,
        atom: AtomId,
        ast: &'a BondAst,
        molecule: &'a MoleculeAst,
    ) -> Self {
        Self {
            bond,
            atom,
            ast,
            molecule,
        }
    }
}
