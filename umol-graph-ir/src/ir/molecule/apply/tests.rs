use rstest::*;
use umol_chem::element::Element;

use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::atom::AtomForm;
use crate::ir::bond::BondForm;
use crate::ir::dative::DativeBondForm;
use crate::ir::id::{AtomId, BondId};
use crate::ir::ligand::{StereoLigand, StereoLigandKind};
use crate::ir::molecule::MoleculeEntries;
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
use crate::ir::stereo::{StereoAtomForm, StereoBondForm, StereoKind};

#[fixture]
fn removal_entries() -> MoleculeEntries {
    MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        dative: vec![
            (vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1)),
            (vec![AtomId(2)], AtomId(3), DativeBondForm::from_order(1)),
        ],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(2), AtomId(3)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(1)],
                MulticenterBondForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(2), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![2, 1]),
            ),
        ],
        noncovalent: vec![
            (
                [AtomId(0), AtomId(1)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                [AtomId(2), AtomId(3)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::default(),
            ),
            (
                AtomId(2),
                vec![
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::default(),
            ),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 1_u32),
            ),
            (
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 1_u32),
            ),
        ],
        constraints: Default::default(),
    }
}

#[fixture]
fn stereo_bond_entries(mut removal_entries: MoleculeEntries) -> MoleculeEntries {
    removal_entries
        .bonds
        .push((AtomId(0), AtomId(2), BondForm::from_order(1)));
    for (_, _, attributes) in &mut removal_entries.stereo_bonds {
        *attributes = StereoBondForm::default();
    }
    removal_entries
}

mod additions;
mod entity_constraints;
mod fields;
mod handles;
mod molecule_constraints;
mod removals;
mod undo;
