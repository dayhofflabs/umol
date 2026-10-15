use rstest::*;
use umol_chem::element::Element;

use super::super::ApplicationState;
use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::atom::AtomForm;
use crate::ir::bond::BondForm;
use crate::ir::dative::DativeBondForm;
use crate::ir::edit::{
    AddAromaticSystem, AddBond, AddMulticenterBond, AddStereoAtom, AddStereoBond,
    AddedAromaticSystem, AddedAtom, AddedBond, AddedDativeBond, AddedMulticenterBond,
    AddedNoncovalentBond, AddedStereoAtom, AddedStereoBond, AromaticSystemHandle, AtomHandle,
    BondHandle, DativeBondHandle, Edit, MulticenterBondHandle, NoncovalentBondHandle,
    StereoAtomHandle, StereoBondHandle, Undo,
};
use crate::ir::entity::EntityKind;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::ligand::{StereoLigand, StereoLigandKind};
use crate::ir::molecule::MoleculeEntries;
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
use crate::ir::stereo::{StereoAtomForm, StereoBondForm, StereoKind};
use crate::ir::{Molecule, TransactionError};

#[rstest]
fn test_molecule_apply_edit_additions(#[values(false, true)] journaled: bool) {
    let initial = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edits = [
        Edit::AddAtoms {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::F),
            ],
        },
        Edit::AddBonds {
            bonds: vec![
                AddBond {
                    endpoints: [AtomHandle::Id(AtomId(1)), AtomHandle::New(0)],
                    attributes: BondForm::from_order(2),
                },
                AddBond {
                    endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::New(1)],
                    attributes: BondForm::from_order(1),
                },
                AddBond {
                    endpoints: [AtomHandle::New(0), AtomHandle::New(1)],
                    attributes: BondForm::from_order(1),
                },
            ],
        },
        Edit::AddDativeBond {
            donors: vec![AtomHandle::New(1), AtomHandle::Id(AtomId(0))],
            acceptor: AtomHandle::New(0),
            attributes: DativeBondForm::from_order(1),
        },
        Edit::AddAromaticSystem(Box::new(AddAromaticSystem {
            atoms: vec![
                AtomHandle::New(0),
                AtomHandle::Id(AtomId(1)),
                AtomHandle::Id(AtomId(0)),
            ],
            attributes: AromaticSystemForm::from_electrons(vec![2, 1, 1]),
        })),
        Edit::AddMulticenterBond(Box::new(AddMulticenterBond {
            atoms: vec![AtomHandle::New(1), AtomHandle::New(0)],
            attributes: MulticenterBondForm::from_electrons(vec![1, 2]),
        })),
        Edit::AddNoncovalentBond {
            atoms: [AtomHandle::New(1), AtomHandle::Id(AtomId(0))],
            attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        },
        Edit::AddStereoAtom(Box::new(AddStereoAtom {
            site: AtomHandle::Id(AtomId(0)),
            ligands: vec![
                (AtomHandle::New(1), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                (
                    AtomHandle::Id(AtomId(0)),
                    StereoLigandKind::ImplicitHydrogen,
                ),
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
            ],
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        })),
        Edit::AddStereoBond(Box::new(AddStereoBond {
            site: BondHandle::New(0),
            ligands: vec![
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                (
                    AtomHandle::Id(AtomId(1)),
                    StereoLigandKind::ImplicitHydrogen,
                ),
                (AtomHandle::New(1), StereoLigandKind::Atom),
                (AtomHandle::New(0), StereoLigandKind::ImplicitHydrogen),
            ],
            attributes: StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        })),
    ];
    let mut undos = Vec::new();
    for edit in edits {
        if journaled {
            undos.push(
                molecule
                    .apply_edit_with_undo(edit, &mut state)
                    .unwrap()
                    .unwrap(),
            );
        } else {
            molecule.apply_edit(edit, &mut state).unwrap();
        }
    }

    let stereo_atom_ligands = vec![
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
    ];
    let stereo_bond_ligands = vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
    ];
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::F),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        dative: vec![(
            vec![AtomId(3), AtomId(0)],
            AtomId(2),
            DativeBondForm::from_order(1),
        )],
        aromatic: vec![(
            vec![AtomId(2), AtomId(1), AtomId(0)],
            AromaticSystemForm::from_electrons(vec![2, 1, 1]),
        )],
        multicenter: vec![(
            vec![AtomId(3), AtomId(2)],
            MulticenterBondForm::from_electrons(vec![1, 2]),
        )],
        noncovalent: vec![(
            [AtomId(3), AtomId(0)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(0),
            stereo_atom_ligands.clone(),
            StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        )],
        stereo_bonds: vec![(
            BondId(1),
            stereo_bond_ligands.clone(),
            StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        )],
        ..Default::default()
    });
    assert_eq!(molecule, expected);
    assert_eq!(state.atom(AtomHandle::Id(AtomId(1))), Ok(AtomId(1)));
    assert_eq!(state.atom(AtomHandle::New(0)), Ok(AtomId(2)));
    assert_eq!(state.atom(AtomHandle::New(1)), Ok(AtomId(3)));
    assert_eq!(state.bond(BondHandle::Id(BondId(0))), Ok(BondId(0)));
    assert_eq!(state.bond(BondHandle::New(0)), Ok(BondId(1)));
    assert_eq!(state.bond(BondHandle::New(2)), Ok(BondId(3)));
    assert_eq!(
        state.dative_bond(DativeBondHandle::New(0)),
        Ok(DativeBondId(0))
    );
    assert_eq!(
        state.aromatic_system(AromaticSystemHandle::New(0)),
        Ok(AromaticSystemId(0))
    );
    assert_eq!(
        state.multicenter_bond(MulticenterBondHandle::New(0)),
        Ok(MulticenterBondId(0))
    );
    assert_eq!(
        state.noncovalent_bond(NoncovalentBondHandle::New(0)),
        Ok(NoncovalentBondId(0))
    );
    assert_eq!(
        state.stereo_atom(StereoAtomHandle::New(0)),
        Ok(StereoAtomId(0))
    );
    assert_eq!(
        state.stereo_bond(StereoBondHandle::New(0)),
        Ok(StereoBondId(0))
    );

    if journaled {
        assert_eq!(
            undos,
            vec![
                Undo::RemoveAddedTopology {
                    atoms: vec![
                        AddedAtom {
                            id: AtomId(2),
                            attributes: AtomForm::from_element(Element::O)
                        },
                        AddedAtom {
                            id: AtomId(3),
                            attributes: AtomForm::from_element(Element::F)
                        },
                    ],
                    bonds: vec![],
                },
                Undo::RemoveAddedTopology {
                    atoms: vec![],
                    bonds: vec![
                        AddedBond {
                            id: BondId(1),
                            endpoints: [AtomId(1), AtomId(2)],
                            attributes: BondForm::from_order(2)
                        },
                        AddedBond {
                            id: BondId(2),
                            endpoints: [AtomId(0), AtomId(3)],
                            attributes: BondForm::from_order(1)
                        },
                        AddedBond {
                            id: BondId(3),
                            endpoints: [AtomId(2), AtomId(3)],
                            attributes: BondForm::from_order(1)
                        },
                    ],
                },
                Undo::RemoveAddedDativeBond(AddedDativeBond {
                    id: DativeBondId(0),
                    donors: vec![AtomId(3), AtomId(0)],
                    acceptor: AtomId(2),
                    attributes: DativeBondForm::from_order(1)
                }),
                Undo::RemoveAddedAromaticSystem(Box::new(AddedAromaticSystem {
                    id: AromaticSystemId(0),
                    atoms: vec![AtomId(2), AtomId(1), AtomId(0)],
                    attributes: AromaticSystemForm::from_electrons(vec![2, 1, 1])
                })),
                Undo::RemoveAddedMulticenterBond(Box::new(AddedMulticenterBond {
                    id: MulticenterBondId(0),
                    atoms: vec![AtomId(3), AtomId(2)],
                    attributes: MulticenterBondForm::from_electrons(vec![1, 2])
                })),
                Undo::RemoveAddedNoncovalentBond(AddedNoncovalentBond {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(3), AtomId(0)],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)
                }),
                Undo::RemoveAddedStereoAtom(Box::new(AddedStereoAtom {
                    id: StereoAtomId(0),
                    site: AtomId(0),
                    ligands: stereo_atom_ligands,
                    attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32)
                })),
                Undo::RemoveAddedStereoBond(Box::new(AddedStereoBond {
                    id: StereoBondId(0),
                    site: BondId(1),
                    ligands: stereo_bond_ligands,
                    attributes: StereoBondForm::new(StereoKind::CisTrans, 1_u32)
                })),
            ]
        );
        let mut restored = molecule.clone();
        for undo in undos.iter().rev().cloned() {
            restored.apply_undo(undo);
        }
        assert_eq!(restored, initial);
        for undo in undos.into_iter().rev() {
            molecule.apply_undo(undo);
        }
        assert_eq!(molecule, initial);
    }
}

#[rstest]
#[case::first(0)]
#[case::middle(1)]
#[case::last(2)]
fn test_molecule_apply_edit_add_bonds_error(
    #[case] invalid_position: usize,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        ..Default::default()
    });
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let mut bonds = vec![
        AddBond {
            endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
            attributes: BondForm::from_order(1),
        },
        AddBond {
            endpoints: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
            attributes: BondForm::from_order(2),
        },
        AddBond {
            endpoints: [AtomHandle::Id(AtomId(2)), AtomHandle::Id(AtomId(0))],
            attributes: BondForm::from_order(1),
        },
    ];
    bonds[invalid_position].endpoints[1] = AtomHandle::Id(AtomId(3));
    let edit = Edit::AddBonds { bonds };
    let error = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).unwrap_err()
    } else {
        molecule.apply_edit(edit, &mut state).unwrap_err()
    };

    assert_eq!(
        error,
        TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 3,
            count: 3
        }
    );
    assert_eq!(molecule, initial);
    assert_eq!(
        state.bond(BondHandle::New(0)),
        Err(TransactionError::HandleOutOfRange {
            kind: EntityKind::Bond,
            index: 0,
            count: 0
        })
    );
}
