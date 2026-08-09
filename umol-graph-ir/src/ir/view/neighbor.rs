//! Neighbor view: yielded by `MoleculeAst::neighbors`.

use super::super::id::{AtomId, BondId};
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use crate::ir::aromatic::AromaticSystemAst;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::dative::DativeBondAst;
    use crate::ir::id::{AtomId, BondId};
    use crate::ir::molecule::{MoleculeAst, MoleculeEntries};
    use crate::ir::multicenter::MulticenterBondAst;
    use crate::ir::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};

    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(3), DativeBondAst::from_order(1))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemAst::default(),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondAst::default(),
            )],
            noncovalent: vec![(
                AtomId(0),
                AtomId(3),
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        })
    }

    #[rstest]
    fn test_neighbor_view_fields(molecule: MoleculeAst) {
        let collected: Vec<(BondId, AtomId, BondForm)> = molecule
            .neighbors(AtomId(2))
            .map(|n| (n.bond_id(), n.atom_id(), n.bond().ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(1), AtomId(1), BondForm::from_order(2)),
                (BondId(2), AtomId(3), BondForm::from_order(1)),
            ],
        );
    }
}
