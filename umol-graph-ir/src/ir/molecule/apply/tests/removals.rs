use rstest::*;
use umol_chem::element::Element;

use super::super::ApplicationState;
use super::removal_entries;
use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::atom::AtomForm;
use crate::ir::bond::BondForm;
use crate::ir::constraint::{
    AromaticSystemConstraintForm, AtomConstraintForm, Constraint, DativeBondConstraintForm,
    MulticenterBondConstraintForm, NoncovalentBondConstraintForm, RelationalConstraint,
    StereoAtomConstraintForm, StereoBondConstraintForm, StereogenicityForm,
};
use crate::ir::dative::DativeBondForm;
use crate::ir::edit::{
    AromaticSystemHandle, AtomHandle, BondHandle, CascadedConstraints, DativeBondHandle, Edit,
    ModifiedConstraint, MulticenterBondHandle, NoncovalentBondHandle, RemovedConstraint,
    StereoAtomHandle, StereoBondHandle, Undo,
};
use crate::ir::electrons::ElectronCountsForm;
use crate::ir::entity::EntityKind;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::ligand::{StereoLigand, StereoLigandKind};
use crate::ir::molecule::MoleculeEntries;
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
use crate::ir::num::NumForm;
use crate::ir::stereo::{StereoAtomForm, StereoBondForm, StereoKind};
use crate::ir::{Molecule, TransactionError};

#[rstest]
fn test_molecule_apply_edit_remove_topology(
    mut removal_entries: MoleculeEntries,
    #[values(false, true)] journaled: bool,
) {
    let removed_constraint = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1));
    let old_constraint = Constraint::StereoBond(
        StereoBondId(1),
        StereoKind::CisTrans,
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    );
    let new_constraint = Constraint::StereoBond(
        StereoBondId(0),
        StereoKind::CisTrans,
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    );
    removal_entries.constraints = vec![
        removed_constraint.clone(),
        old_constraint.clone(),
        removed_constraint.clone(),
    ]
    .into();
    let initial = Molecule::from_entries(removal_entries);
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    molecule
        .apply_edit(
            Edit::AddAtoms {
                atoms: vec![AtomForm::from_element(Element::F)],
            },
            &mut state,
        )
        .unwrap();
    let before_removal = molecule.clone();
    let edit = Edit::RemoveTopology {
        atoms: vec![AtomHandle::Id(AtomId(0))],
        bonds: vec![],
    };
    let undo = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).unwrap()
    } else {
        molecule.apply_edit(edit, &mut state).unwrap();
        None
    };
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
        ],
        bonds: vec![(AtomId(1), AtomId(2), BondForm::from_order(1))],
        dative: vec![(vec![AtomId(1)], AtomId(2), DativeBondForm::from_order(1))],
        aromatic: vec![(
            vec![AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        multicenter: vec![(
            vec![AtomId(1), AtomId(2)],
            MulticenterBondForm::from_electrons(vec![2, 1]),
        )],
        noncovalent: vec![(
            [AtomId(1), AtomId(2)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(1),
            vec![
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
            ],
            StereoAtomForm::default(),
        )],
        stereo_bonds: vec![(
            BondId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
            ],
            StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        )],
        constraints: vec![new_constraint.clone()].into(),
    });
    assert_eq!(molecule, expected);
    assert_eq!(
        state.atom(AtomHandle::Id(AtomId(0))),
        Err(TransactionError::HandleRemoved {
            kind: EntityKind::Atom,
            index: 0,
        })
    );
    assert_eq!(state.atom(AtomHandle::Id(AtomId(3))), Ok(AtomId(2)));
    assert_eq!(state.atom(AtomHandle::New(0)), Ok(AtomId(3)));
    assert_eq!(state.bond(BondHandle::Id(BondId(1))), Ok(BondId(0)));
    assert_eq!(
        state.dative_bond(DativeBondHandle::Id(DativeBondId(1))),
        Ok(DativeBondId(0))
    );
    assert_eq!(
        state.aromatic_system(AromaticSystemHandle::Id(AromaticSystemId(1))),
        Ok(AromaticSystemId(0))
    );
    assert_eq!(
        state.multicenter_bond(MulticenterBondHandle::Id(MulticenterBondId(1))),
        Ok(MulticenterBondId(0))
    );
    assert_eq!(
        state.noncovalent_bond(NoncovalentBondHandle::Id(NoncovalentBondId(1))),
        Ok(NoncovalentBondId(0))
    );
    assert_eq!(
        state.stereo_atom(StereoAtomHandle::Id(StereoAtomId(1))),
        Ok(StereoAtomId(0))
    );
    assert_eq!(
        state.stereo_bond(StereoBondHandle::Id(StereoBondId(1))),
        Ok(StereoBondId(0))
    );
    if journaled {
        let undo = undo.unwrap();
        let Undo::RestoreRemovedTopology(removed) = &undo else {
            panic!("expected topology restoration");
        };
        assert_eq!(
            &removed.cascade,
            &CascadedConstraints {
                removed: vec![
                    RemovedConstraint {
                        position: 0,
                        constraint: removed_constraint.clone()
                    },
                    RemovedConstraint {
                        position: 2,
                        constraint: removed_constraint
                    },
                ],
                modified: vec![ModifiedConstraint {
                    position: 1,
                    old: old_constraint,
                    new: new_constraint
                }],
            }
        );
        let mut restored = molecule.clone();
        restored.apply_undo(undo.clone());
        assert_eq!(restored, before_removal);
        molecule.apply_undo(undo);
        assert_eq!(molecule, before_removal);
    }
}

#[rstest]
#[case::identity(vec![])]
#[case::site(vec![BondHandle::Id(BondId(0))])]
fn test_molecule_apply_edit_remove_bonds(
    mut removal_entries: MoleculeEntries,
    #[case] bonds: Vec<BondHandle>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.stereo_atoms.clear();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    if !bonds.is_empty() {
        removal_entries.bonds.remove(0);
        removal_entries.stereo_bonds.remove(0);
        removal_entries.stereo_bonds[0].0 = BondId(0);
    }
    let expected = Molecule::from_entries(removal_entries);
    let edit = Edit::RemoveTopology {
        atoms: vec![],
        bonds,
    };
    let undo = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).unwrap()
    } else {
        molecule.apply_edit(edit, &mut state).unwrap();
        None
    };
    assert_eq!(molecule, expected);
    if journaled {
        molecule.apply_undo(undo.unwrap());
        assert_eq!(molecule, initial);
    }
}

#[rstest]
#[case::dative(EntityKind::DativeBond)]
#[case::aromatic(EntityKind::AromaticSystem)]
#[case::multicenter(EntityKind::MulticenterBond)]
#[case::noncovalent(EntityKind::NoncovalentBond)]
#[case::stereo_atom(EntityKind::StereoAtom)]
#[case::stereo_bond(EntityKind::StereoBond)]
fn test_molecule_apply_edit_remove_overlays(
    mut removal_entries: MoleculeEntries,
    #[case] kind: EntityKind,
    #[values(false, true)] journaled: bool,
) {
    let (edit, removed_constraint, old_constraint) = match kind {
        EntityKind::DativeBond => (
            Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(0))],
                    AtomHandle::Id(AtomId(1)),
                    DativeBondForm::from_order(1),
                )],
            },
            Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(false)),
            Constraint::DativeBond(DativeBondId(1), DativeBondConstraintForm::aromatic(false)),
        ),
        EntityKind::AromaticSystem => (
            Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
                    AromaticSystemForm::from_electrons(vec![2, 1]),
                )],
            },
            Constraint::AromaticSystem(
                AromaticSystemId(0),
                AromaticSystemConstraintForm::electron_count(3),
            ),
            Constraint::AromaticSystem(
                AromaticSystemId(1),
                AromaticSystemConstraintForm::electron_count(3),
            ),
        ),
        EntityKind::MulticenterBond => (
            Edit::RemoveMulticenterBonds {
                removes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
                    MulticenterBondForm::from_electrons(vec![2, 1]),
                )],
            },
            Constraint::MulticenterBond(
                MulticenterBondId(0),
                MulticenterBondConstraintForm::ElectronCount(3.into()),
            ),
            Constraint::MulticenterBond(
                MulticenterBondId(1),
                MulticenterBondConstraintForm::ElectronCount(3.into()),
            ),
        ),
        EntityKind::NoncovalentBond => (
            Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
            },
            Constraint::NoncovalentBond(
                NoncovalentBondId(0),
                NoncovalentBondConstraintForm::intramolecular(true),
            ),
            Constraint::NoncovalentBond(
                NoncovalentBondId(1),
                NoncovalentBondConstraintForm::intramolecular(true),
            ),
        ),
        EntityKind::StereoAtom => (
            Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                        (
                            AtomHandle::Id(AtomId(0)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::default(),
                )],
            },
            Constraint::Relational(RelationalConstraint::StereoAtomAllLigands {
                stereo_atom: StereoAtomId(0),
                predicate: Box::new(AtomConstraintForm::degree(1)),
            }),
            Constraint::Relational(RelationalConstraint::StereoAtomAllLigands {
                stereo_atom: StereoAtomId(1),
                predicate: Box::new(AtomConstraintForm::degree(1)),
            }),
        ),
        EntityKind::StereoBond => (
            Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                        (
                            AtomHandle::Id(AtomId(0)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::LonePair),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, 0_u32),
                )],
            },
            Constraint::StereoBond(
                StereoBondId(0),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::StereoBond(
                StereoBondId(1),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
        ),
        EntityKind::Atom | EntityKind::Bond => unreachable!(),
    };
    let unchanged = Constraint::Atom(AtomId(3), AtomConstraintForm::degree(1));
    removal_entries.constraints = vec![
        removed_constraint.clone(),
        unchanged.clone(),
        old_constraint.clone(),
        removed_constraint.clone(),
    ]
    .into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let undo = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).unwrap()
    } else {
        molecule.apply_edit(edit, &mut state).unwrap();
        None
    };
    match kind {
        EntityKind::DativeBond => {
            removal_entries.dative.remove(0);
        }
        EntityKind::AromaticSystem => {
            removal_entries.aromatic.remove(0);
        }
        EntityKind::MulticenterBond => {
            removal_entries.multicenter.remove(0);
        }
        EntityKind::NoncovalentBond => {
            removal_entries.noncovalent.remove(0);
        }
        EntityKind::StereoAtom => {
            removal_entries.stereo_atoms.remove(0);
        }
        EntityKind::StereoBond => {
            removal_entries.stereo_bonds.remove(0);
        }
        EntityKind::Atom | EntityKind::Bond => unreachable!(),
    }
    removal_entries.constraints = vec![unchanged, removed_constraint.clone()].into();
    assert_eq!(molecule, Molecule::from_entries(removal_entries));
    assert_eq!(state.atom(AtomHandle::Id(AtomId(3))), Ok(AtomId(3)));
    assert_eq!(state.bond(BondHandle::Id(BondId(1))), Ok(BondId(1)));
    for (current_kind, removed, survivor) in [
        (
            EntityKind::DativeBond,
            state
                .dative_bond(DativeBondHandle::Id(DativeBondId(0)))
                .map(|id| id.index()),
            state
                .dative_bond(DativeBondHandle::Id(DativeBondId(1)))
                .map(|id| id.index()),
        ),
        (
            EntityKind::AromaticSystem,
            state
                .aromatic_system(AromaticSystemHandle::Id(AromaticSystemId(0)))
                .map(|id| id.index()),
            state
                .aromatic_system(AromaticSystemHandle::Id(AromaticSystemId(1)))
                .map(|id| id.index()),
        ),
        (
            EntityKind::MulticenterBond,
            state
                .multicenter_bond(MulticenterBondHandle::Id(MulticenterBondId(0)))
                .map(|id| id.index()),
            state
                .multicenter_bond(MulticenterBondHandle::Id(MulticenterBondId(1)))
                .map(|id| id.index()),
        ),
        (
            EntityKind::NoncovalentBond,
            state
                .noncovalent_bond(NoncovalentBondHandle::Id(NoncovalentBondId(0)))
                .map(|id| id.index()),
            state
                .noncovalent_bond(NoncovalentBondHandle::Id(NoncovalentBondId(1)))
                .map(|id| id.index()),
        ),
        (
            EntityKind::StereoAtom,
            state
                .stereo_atom(StereoAtomHandle::Id(StereoAtomId(0)))
                .map(|id| id.index()),
            state
                .stereo_atom(StereoAtomHandle::Id(StereoAtomId(1)))
                .map(|id| id.index()),
        ),
        (
            EntityKind::StereoBond,
            state
                .stereo_bond(StereoBondHandle::Id(StereoBondId(0)))
                .map(|id| id.index()),
            state
                .stereo_bond(StereoBondHandle::Id(StereoBondId(1)))
                .map(|id| id.index()),
        ),
    ] {
        if current_kind == kind {
            assert_eq!(
                removed,
                Err(TransactionError::HandleRemoved { kind, index: 0 })
            );
            assert_eq!(survivor, Ok(0));
        } else {
            assert_eq!(removed, Ok(0));
            assert_eq!(survivor, Ok(1));
        }
    }
    if journaled {
        let undo = undo.unwrap();
        let cascade = match &undo {
            Undo::RestoreRemovedDativeBonds(removed) => &removed.cascade,
            Undo::RestoreRemovedAromaticSystems(removed) => &removed.cascade,
            Undo::RestoreRemovedMulticenterBonds(removed) => &removed.cascade,
            Undo::RestoreRemovedNoncovalentBonds(removed) => &removed.cascade,
            Undo::RestoreRemovedStereoAtoms(removed) => &removed.cascade,
            Undo::RestoreRemovedStereoBonds(removed) => &removed.cascade,
            _ => panic!("expected overlay restoration"),
        };
        assert_eq!(
            cascade,
            &CascadedConstraints {
                removed: vec![
                    RemovedConstraint {
                        position: 0,
                        constraint: removed_constraint.clone()
                    },
                    RemovedConstraint {
                        position: 3,
                        constraint: removed_constraint.clone()
                    },
                ],
                modified: vec![ModifiedConstraint {
                    position: 2,
                    old: old_constraint,
                    new: removed_constraint
                }],
            }
        );
        let mut restored = molecule.clone();
        restored.apply_undo(undo.clone());
        assert_eq!(restored, initial);
        molecule.apply_undo(undo);
        assert_eq!(molecule, initial);
    }
}

#[rstest]
#[case::dative_bond(EntityKind::DativeBond)]
#[case::aromatic_system(EntityKind::AromaticSystem)]
#[case::multicenter_bond(EntityKind::MulticenterBond)]
#[case::noncovalent_bond(EntityKind::NoncovalentBond)]
#[case::stereo_atom(EntityKind::StereoAtom)]
#[case::stereo_bond(EntityKind::StereoBond)]
fn test_molecule_apply_edit_remove_incidence(
    mut removal_entries: MoleculeEntries,
    #[case] kind: EntityKind,
    #[values(false, true)] mismatched: bool,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.dative[0].2.order = NumForm::lit_set([]);
    removal_entries.aromatic[0].1.charge = NumForm::lit_set([]);
    removal_entries.multicenter[0].1.charge = NumForm::lit_set([]);
    removal_entries.stereo_atoms[0]
        .2
        .constraints
        .set(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::LitSet(Default::default()),
        ));
    removal_entries.stereo_bonds[0]
        .2
        .constraints
        .set(StereoBondConstraintForm::Stereogenicity(
            StereogenicityForm::LitSet(Default::default()),
        ));
    let initial = Molecule::from_entries(removal_entries.clone());
    let edit = match kind {
        EntityKind::DativeBond => {
            let (donors, acceptor, attributes) = removal_entries.dative.remove(0);
            Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    donors.into_iter().map(AtomHandle::Id).collect(),
                    AtomHandle::Id(if mismatched { AtomId(3) } else { acceptor }),
                    attributes,
                )],
            }
        }
        EntityKind::AromaticSystem => {
            let (mut atoms, mut attributes) = removal_entries.aromatic.remove(0);
            atoms.reverse();
            attributes.electrons = ElectronCountsForm::from(vec![2, 1]);
            if mismatched {
                atoms[0] = AtomId(2);
            }
            Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    atoms.into_iter().map(AtomHandle::Id).collect(),
                    attributes,
                )],
            }
        }
        EntityKind::MulticenterBond => {
            let (mut atoms, mut attributes) = removal_entries.multicenter.remove(0);
            atoms.reverse();
            attributes.electrons = ElectronCountsForm::from(vec![2, 1]);
            if mismatched {
                atoms[0] = AtomId(2);
            }
            Edit::RemoveMulticenterBonds {
                removes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    atoms.into_iter().map(AtomHandle::Id).collect(),
                    attributes,
                )],
            }
        }
        EntityKind::NoncovalentBond => {
            let (mut atoms, attributes) = removal_entries.noncovalent.remove(0);
            atoms.reverse();
            if mismatched {
                atoms[0] = AtomId(2);
            }
            Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    atoms.map(AtomHandle::Id),
                    attributes,
                )],
            }
        }
        EntityKind::StereoAtom => {
            let (site, mut ligands, attributes) = removal_entries.stereo_atoms.remove(0);
            ligands.reverse();
            Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(if mismatched { AtomId(2) } else { site }),
                    ligands
                        .into_iter()
                        .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                        .collect(),
                    attributes,
                )],
            }
        }
        EntityKind::StereoBond => {
            let (site, mut ligands, attributes) = removal_entries.stereo_bonds.remove(0);
            if mismatched {
                ligands.swap(1, 2);
            } else {
                ligands.rotate_left(2);
            }
            Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(site),
                    ligands
                        .into_iter()
                        .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                        .collect(),
                    attributes,
                )],
            }
        }
        EntityKind::Atom | EntityKind::Bond => unreachable!(),
    };
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let result = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state)
    } else {
        molecule.apply_edit(edit, &mut state).map(|()| None)
    };
    if mismatched {
        assert_eq!(result, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, initial);
    } else {
        let undo = result.unwrap();
        assert_eq!(molecule, Molecule::from_entries(removal_entries));
        if journaled {
            molecule.apply_undo(undo.unwrap());
            assert_eq!(molecule, initial);
        }
    }
}

#[rstest]
#[case::duplicate_atom(Edit::RemoveTopology {
    atoms: vec![AtomHandle::Id(AtomId(0)); 2], bonds: vec![],
}, TransactionError::DuplicateRemoval { kind: EntityKind::Atom })]
#[case::duplicate_bond(Edit::RemoveTopology {
    atoms: vec![], bonds: vec![BondHandle::Id(BondId(0)); 2],
}, TransactionError::DuplicateRemoval { kind: EntityKind::Bond })]
#[case::missing_atom(Edit::RemoveTopology {
    atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(4))], bonds: vec![],
}, TransactionError::HandleOutOfRange { kind: EntityKind::Atom, index: 4, count: 4 })]
#[case::missing_bond(Edit::RemoveTopology {
    atoms: vec![AtomHandle::Id(AtomId(0))], bonds: vec![BondHandle::Id(BondId(2))],
}, TransactionError::HandleOutOfRange { kind: EntityKind::Bond, index: 2, count: 2 })]
#[case::duplicate_dative_bond(Edit::RemoveDativeBonds { removes: vec![(
    DativeBondHandle::Id(DativeBondId(0)), vec![AtomHandle::Id(AtomId(0))],
    AtomHandle::Id(AtomId(1)), DativeBondForm::from_order(1),
); 2] }, TransactionError::DuplicateRemoval { kind: EntityKind::DativeBond })]
#[case::dative_bond_old_state(Edit::RemoveDativeBonds { removes: vec![
    (DativeBondHandle::Id(DativeBondId(0)), vec![AtomHandle::Id(AtomId(0))],
     AtomHandle::Id(AtomId(1)), DativeBondForm::from_order(1)),
    (DativeBondHandle::Id(DativeBondId(1)), vec![AtomHandle::Id(AtomId(2))],
     AtomHandle::Id(AtomId(3)), DativeBondForm::from_order(2)),
] }, TransactionError::OldStateMismatch)]
fn test_molecule_apply_edit_remove_error(
    removal_entries: MoleculeEntries,
    #[case] edit: Edit,
    #[case] expected: TransactionError,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries);
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let result = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(result, Err(expected));
    assert_eq!(molecule, initial);
    assert_eq!(state.atoms.initial, None);
    assert_eq!(state.bonds.initial, None);
    assert_eq!(state.dative_bonds.initial, None);
    assert_eq!(state.aromatic_systems.initial, None);
    assert_eq!(state.multicenter_bonds.initial, None);
    assert_eq!(state.noncovalent_bonds.initial, None);
    assert_eq!(state.stereo_atoms.initial, None);
    assert_eq!(state.stereo_bonds.initial, None);
}
