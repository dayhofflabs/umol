use rstest::*;
use umol_graph_core::Compaction;

use super::super::{ApplicationState, HandleTable};
use super::stereo_bond_entries;
use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::atom::AtomForm;
use crate::ir::bond::BondForm;
use crate::ir::constraint::{AtomConstraintForm, AtomConstraintKey};
use crate::ir::dative::DativeBondForm;
use crate::ir::edit::{AtomHandle, Edit, Edits};
use crate::ir::electrons::ElectronCountsForm;
use crate::ir::entity::EntityKind;
use crate::ir::id::AtomId;
use crate::ir::molecule::MoleculeEntries;
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm};
use crate::ir::num::NumForm;
use crate::ir::{
    AromaticSystemConstraintForm, AromaticSystemConstraintKey, AromaticSystemHandle,
    AromaticSystemId, BondConstraintForm, BondConstraintKey, BondHandle, BondId,
    DativeBondConstraintForm, DativeBondConstraintKey, DativeBondHandle, DativeBondId, Molecule,
    MulticenterBondConstraintForm, MulticenterBondConstraintKey, MulticenterBondHandle,
    MulticenterBondId, NoncovalentBondConstraintForm, NoncovalentBondConstraintKey,
    NoncovalentBondHandle, NoncovalentBondId, RollbackMethod, StereoAtomConstraintForm,
    StereoAtomConstraintKey, StereoAtomForm, StereoAtomHandle, StereoAtomId,
    StereoBondConstraintForm, StereoBondConstraintKey, StereoBondForm, StereoBondHandle,
    StereoBondId, StereoConfigurationForm, StereoKind, StereoLigand, StereoLigandKind,
    TransactionError,
};

#[rstest]
#[case::first(3, 0, Ok(AtomId(0)))]
#[case::last(3, 2, Ok(AtomId(2)))]
#[case::outside(3, 3, Err(TransactionError::HandleOutOfRange {
    kind: EntityKind::Atom, index: 3, count: 3,
}))]
#[case::empty(0, 0, Err(TransactionError::HandleOutOfRange {
    kind: EntityKind::Atom, index: 0, count: 0,
}))]
fn test_handle_table_initial(
    #[case] count: usize,
    #[case] index: usize,
    #[case] expected: Result<AtomId, TransactionError>,
) {
    let table = HandleTable::new(count);

    assert_eq!(table.initial(EntityKind::Atom, index), expected);
    assert_eq!(table.initial, None);
}

#[rstest]
#[case::first(0, Ok(AtomId(3)))]
#[case::second(1, Ok(AtomId(4)))]
#[case::outside(2, Err(TransactionError::HandleOutOfRange {
    kind: EntityKind::Atom, index: 2, count: 2,
}))]
fn test_handle_table_created(
    #[case] index: usize,
    #[case] expected: Result<AtomId, TransactionError>,
) {
    let mut table = HandleTable::new(3);
    table.push(AtomId(3));
    table.push(AtomId(4));

    assert_eq!(table.created(EntityKind::Atom, index), expected);
    assert_eq!(table.initial, None);
}

#[rstest]
#[case::identity(vec![], vec![Some(AtomId(0)), Some(AtomId(1)), Some(AtomId(2))], vec![Some(AtomId(3)), Some(AtomId(4))])]
#[case::initial_and_created(vec![AtomId(1), AtomId(3)], vec![Some(AtomId(0)), None, Some(AtomId(1))], vec![None, Some(AtomId(2))])]
#[case::all(vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)], vec![None, None, None], vec![None, None])]
fn test_handle_table_compact(
    #[case] removed: Vec<AtomId>,
    #[case] initial: Vec<Option<AtomId>>,
    #[case] created: Vec<Option<AtomId>>,
) {
    let mut table = HandleTable::new(3);
    table.push(AtomId(3));
    table.push(AtomId(4));
    let compaction = Compaction::new(5, removed).unwrap();

    table.compact(|id| compaction.compact(id));

    assert_eq!(table.initial_count, 3);
    assert_eq!(table.initial, Some(initial));
    assert_eq!(table.created, created);
}

#[rstest]
fn test_handle_table_compact_sequence() {
    let mut table = HandleTable::new(3);
    table.push(AtomId(3));
    table.push(AtomId(4));
    let first = Compaction::new(5, vec![AtomId(1), AtomId(3)]).unwrap();
    table.compact(|id| first.compact(id));
    table.push(AtomId(3));
    let second = Compaction::new(4, vec![AtomId(0), AtomId(2)]).unwrap();
    table.compact(|id| second.compact(id));

    assert_eq!(table.initial_count, 3);
    assert_eq!(table.initial, Some(vec![None, None, Some(AtomId(0))]));
    assert_eq!(table.created, vec![None, None, Some(AtomId(1))]);
    assert_eq!(table.initial(EntityKind::Atom, 2), Ok(AtomId(0)));
    assert_eq!(table.created(EntityKind::Atom, 2), Ok(AtomId(1)));
    assert_eq!(
        table.initial(EntityKind::Atom, 1),
        Err(TransactionError::HandleRemoved {
            kind: EntityKind::Atom,
            index: 1
        }),
    );
    assert_eq!(
        table.created(EntityKind::Atom, 1),
        Err(TransactionError::HandleRemoved {
            kind: EntityKind::Atom,
            index: 1
        }),
    );
    assert_eq!(
        table.initial(EntityKind::Atom, 3),
        Err(TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 3,
            count: 3
        }),
    );
    assert_eq!(
        table.created(EntityKind::Atom, 3),
        Err(TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 3,
            count: 3
        }),
    );
}

#[rstest]
fn test_molecule_apply_atom_batches_removed(
    #[values(false, true)] created: bool,
    #[values(false, true)] constraint: bool,
    #[values(false, true)] journaled: bool,
) {
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 2],
        ..Default::default()
    });
    let mut prefix = Edits::new();
    let removed = if created {
        prefix.add_atom(AtomForm::default())
    } else {
        AtomHandle::Id(AtomId(1))
    };
    prefix.remove_atom(removed.clone());
    let batch = if constraint {
        Edit::ModifyAtomConstraints {
            key: AtomConstraintKey::Valence,
            changes: vec![
                (
                    AtomHandle::Id(AtomId(0)),
                    None,
                    Some(AtomConstraintForm::valence(1)),
                ),
                (removed, None, None),
            ],
        }
    } else {
        Edit::ModifyAtomCharges {
            changes: vec![
                (
                    AtomHandle::Id(AtomId(0)),
                    NumForm::Undetermined,
                    NumForm::Lit(1),
                ),
                (removed, NumForm::Undetermined, NumForm::Lit(1)),
            ],
        }
    };
    let error = TransactionError::HandleRemoved {
        kind: EntityKind::Atom,
        index: if created { 0 } else { 1 },
    };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    for edit in prefix.clone() {
        molecule.apply_edit(edit, &mut state).unwrap();
    }
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut atoms = vec![AtomForm::default(); if created { 2 } else { 1 }];
    if !journaled {
        if constraint {
            atoms[0].constraints = vec![AtomConstraintForm::valence(1)].into();
        } else {
            atoms[0].charge = NumForm::Lit(1);
        }
    }
    assert_eq!(
        molecule,
        Molecule::from_entries(MoleculeEntries {
            atoms,
            ..Default::default()
        })
    );
    prefix.push(batch);
    let mut molecule = source.clone();
    assert_eq!(
        molecule.transact([prefix], RollbackMethod::Journal),
        Err(error.into())
    );
    assert_eq!(molecule, source);
}

#[rstest]
fn test_molecule_apply_bond_batches_removed(
    #[values(false, true)] created: bool,
    #[values(false, true)] constraint: bool,
    #[values(false, true)] journaled: bool,
) {
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::default()),
            (AtomId(1), AtomId(2), BondForm::default()),
        ],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let mut prefix = Edits::new();
    let removed = if created {
        prefix.add_bond(
            AtomHandle::Id(AtomId(2)),
            AtomHandle::Id(AtomId(3)),
            BondForm::default(),
        )
    } else {
        BondHandle::Id(BondId(1))
    };
    prefix.remove_bond(removed.clone());
    let batch = if constraint {
        Edit::ModifyBondConstraints {
            key: BondConstraintKey::Aromatic,
            changes: vec![
                (
                    BondHandle::Id(BondId(0)),
                    None,
                    Some(BondConstraintForm::aromatic(true)),
                ),
                (removed, None, None),
            ],
        }
    } else {
        Edit::ModifyBondOrders {
            changes: vec![
                (
                    BondHandle::Id(BondId(0)),
                    NumForm::Undetermined,
                    NumForm::Lit(2),
                ),
                (removed, NumForm::Undetermined, NumForm::Lit(1)),
            ],
        }
    };
    let error = TransactionError::HandleRemoved {
        kind: EntityKind::Bond,
        index: if created { 0 } else { 1 },
    };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    for edit in prefix.clone() {
        molecule.apply_edit(edit, &mut state).unwrap();
    }
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !created {
        expected.bonds.pop();
    }
    if !journaled {
        if constraint {
            expected.bonds[0].2.constraints = vec![BondConstraintForm::aromatic(true)].into();
        } else {
            expected.bonds[0].2.order = NumForm::Lit(2);
        }
    }
    assert_eq!(molecule, Molecule::from_entries(expected));
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        let edits = Edits::from_iter(prefix.clone().into_iter().chain([batch.clone()]));
        assert_eq!(
            molecule.transact([edits], recovery),
            Err(error.clone().into())
        );
        assert_eq!(molecule, source);
    }
}

#[rstest]
fn test_molecule_apply_dative_bond_batches_removed(
    #[values(false, true)] created: bool,
    #[values(false, true)] constraint: bool,
    #[values(false, true)] journaled: bool,
) {
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 4],
        dative: vec![
            (
                vec![AtomId(2), AtomId(0)],
                AtomId(1),
                DativeBondForm::default(),
            ),
            (
                vec![AtomId(3), AtomId(1)],
                AtomId(2),
                DativeBondForm::default(),
            ),
        ],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let mut prefix = Edits::new();
    let removed = if created {
        prefix.add_dative_bond(
            vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(2))],
            AtomHandle::Id(AtomId(3)),
            DativeBondForm::default(),
        )
    } else {
        DativeBondHandle::Id(DativeBondId(0))
    };
    let removal = if created {
        (
            vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(2))],
            AtomHandle::Id(AtomId(3)),
            DativeBondForm::default(),
        )
    } else {
        (
            vec![AtomHandle::Id(AtomId(2)), AtomHandle::Id(AtomId(0))],
            AtomHandle::Id(AtomId(1)),
            DativeBondForm::default(),
        )
    };
    prefix.remove_dative_bonds(vec![(removed.clone(), removal.0, removal.1, removal.2)]);
    let key = DativeBondConstraintKey::Aromatic;
    let form = DativeBondConstraintForm::undetermined(key);
    let batch = if constraint {
        Edit::ModifyDativeBondConstraints {
            key,
            changes: vec![
                (
                    DativeBondHandle::Id(DativeBondId(1)),
                    None,
                    Some(form.clone()),
                ),
                (removed, None, None),
            ],
        }
    } else {
        Edit::ModifyDativeBondOrders {
            changes: vec![
                (
                    DativeBondHandle::Id(DativeBondId(1)),
                    NumForm::Undetermined,
                    NumForm::Lit(1),
                ),
                (
                    removed,
                    NumForm::Undetermined,
                    NumForm::LitSet(Box::new([1, 2].into())),
                ),
            ],
        }
    };
    let error = TransactionError::HandleRemoved {
        kind: EntityKind::DativeBond,
        index: 0,
    };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    for edit in prefix.clone() {
        molecule.apply_edit(edit, &mut state).unwrap();
    }
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !created {
        expected.dative.remove(0);
    }
    if !journaled {
        let attributes = &mut expected.dative[usize::from(created)].2;
        if constraint {
            attributes.constraints = vec![form].into();
        } else {
            attributes.order = NumForm::Lit(1);
        }
    }
    assert_eq!(molecule, Molecule::from_entries(expected));
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        let edits = Edits::from_iter(prefix.clone().into_iter().chain([batch.clone()]));
        assert_eq!(
            molecule.transact([edits], recovery),
            Err(error.clone().into())
        );
        assert_eq!(molecule, source);
    }
}

#[rstest]
fn test_molecule_apply_aromatic_system_batches_removed(
    #[values(false, true)] created: bool,
    #[values(false, true)] constraint: bool,
    #[values(false, true)] journaled: bool,
) {
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 9],
        aromatic: vec![
            (
                vec![AtomId(2), AtomId(0), AtomId(1)],
                AromaticSystemForm::default(),
            ),
            (
                vec![AtomId(5), AtomId(3), AtomId(4)],
                AromaticSystemForm::default(),
            ),
        ],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let mut prefix = Edits::new();
    let removed = if created {
        prefix.add_aromatic_system(
            vec![
                AtomHandle::Id(AtomId(8)),
                AtomHandle::Id(AtomId(6)),
                AtomHandle::Id(AtomId(7)),
            ],
            AromaticSystemForm::default(),
        )
    } else {
        AromaticSystemHandle::Id(AromaticSystemId(0))
    };
    let removal = if created {
        (
            vec![
                AtomHandle::Id(AtomId(8)),
                AtomHandle::Id(AtomId(6)),
                AtomHandle::Id(AtomId(7)),
            ],
            AromaticSystemForm::default(),
        )
    } else {
        (
            vec![
                AtomHandle::Id(AtomId(2)),
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
            ],
            AromaticSystemForm::default(),
        )
    };
    prefix.remove_aromatic_systems(vec![(removed.clone(), removal.0, removal.1)]);
    let key = AromaticSystemConstraintKey::ElectronCount;
    let form = AromaticSystemConstraintForm::undetermined(key);
    let batch = if constraint {
        Edit::ModifyAromaticSystemConstraints {
            key,
            changes: vec![
                (
                    AromaticSystemHandle::Id(AromaticSystemId(1)),
                    None,
                    Some(form.clone()),
                ),
                (removed, None, None),
            ],
        }
    } else {
        Edit::ModifyAromaticSystemElectrons {
            changes: vec![
                (
                    AromaticSystemHandle::Id(AromaticSystemId(1)),
                    ElectronCountsForm::Undetermined,
                    ElectronCountsForm::Lit(vec![1, 2]),
                ),
                (
                    removed,
                    ElectronCountsForm::Undetermined,
                    ElectronCountsForm::Lit(vec![]),
                ),
            ],
        }
    };
    let error = TransactionError::HandleRemoved {
        kind: EntityKind::AromaticSystem,
        index: 0,
    };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    for edit in prefix.clone() {
        molecule.apply_edit(edit, &mut state).unwrap();
    }
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !created {
        expected.aromatic.remove(0);
    }
    if !journaled {
        let attributes = &mut expected.aromatic[usize::from(created)].1;
        if constraint {
            attributes.constraints = vec![form].into();
        } else {
            attributes.electrons = ElectronCountsForm::Lit(vec![1, 2]);
        }
    }
    assert_eq!(molecule, Molecule::from_entries(expected));
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        let edits = Edits::from_iter(prefix.clone().into_iter().chain([batch.clone()]));
        assert_eq!(
            molecule.transact([edits], recovery),
            Err(error.clone().into())
        );
        assert_eq!(molecule, source);
    }
}

#[rstest]
fn test_molecule_apply_multicenter_bond_batches_removed(
    #[values(false, true)] created: bool,
    #[values(false, true)] constraint: bool,
    #[values(false, true)] journaled: bool,
) {
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 4],
        multicenter: vec![
            (
                vec![AtomId(2), AtomId(0), AtomId(1)],
                MulticenterBondForm::default(),
            ),
            (
                vec![AtomId(3), AtomId(1), AtomId(2)],
                MulticenterBondForm::default(),
            ),
        ],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let mut prefix = Edits::new();
    let removed = if created {
        prefix.add_multicenter_bond(
            vec![
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(2)),
                AtomHandle::Id(AtomId(3)),
            ],
            MulticenterBondForm::default(),
        )
    } else {
        MulticenterBondHandle::Id(MulticenterBondId(0))
    };
    let removal = if created {
        (
            vec![
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(2)),
                AtomHandle::Id(AtomId(3)),
            ],
            MulticenterBondForm::default(),
        )
    } else {
        (
            vec![
                AtomHandle::Id(AtomId(2)),
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
            ],
            MulticenterBondForm::default(),
        )
    };
    prefix.remove_multicenter_bonds(vec![(removed.clone(), removal.0, removal.1)]);
    let key = MulticenterBondConstraintKey::ElectronCount;
    let form = MulticenterBondConstraintForm::undetermined(key);
    let batch = if constraint {
        Edit::ModifyMulticenterBondConstraints {
            key,
            changes: vec![
                (
                    MulticenterBondHandle::Id(MulticenterBondId(1)),
                    None,
                    Some(form.clone()),
                ),
                (removed, None, None),
            ],
        }
    } else {
        Edit::ModifyMulticenterBondElectrons {
            changes: vec![
                (
                    MulticenterBondHandle::Id(MulticenterBondId(1)),
                    ElectronCountsForm::Undetermined,
                    ElectronCountsForm::Lit(vec![1, 2]),
                ),
                (
                    removed,
                    ElectronCountsForm::Undetermined,
                    ElectronCountsForm::Lit(vec![]),
                ),
            ],
        }
    };
    let error = TransactionError::HandleRemoved {
        kind: EntityKind::MulticenterBond,
        index: 0,
    };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    for edit in prefix.clone() {
        molecule.apply_edit(edit, &mut state).unwrap();
    }
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !created {
        expected.multicenter.remove(0);
    }
    if !journaled {
        let attributes = &mut expected.multicenter[usize::from(created)].1;
        if constraint {
            attributes.constraints = vec![form].into();
        } else {
            attributes.electrons = ElectronCountsForm::Lit(vec![1, 2]);
        }
    }
    assert_eq!(molecule, Molecule::from_entries(expected));
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        let edits = Edits::from_iter(prefix.clone().into_iter().chain([batch.clone()]));
        assert_eq!(
            molecule.transact([edits], recovery),
            Err(error.clone().into())
        );
        assert_eq!(molecule, source);
    }
}

#[rstest]
fn test_molecule_apply_noncovalent_bond_batches_removed(
    #[values(false, true)] created: bool,
    #[values(false, true)] constraint: bool,
    #[values(false, true)] journaled: bool,
) {
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 6],
        noncovalent: vec![
            ([AtomId(1), AtomId(0)], NoncovalentBondForm::default()),
            ([AtomId(3), AtomId(2)], NoncovalentBondForm::default()),
        ],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let mut prefix = Edits::new();
    let removed = if created {
        prefix.add_noncovalent_bond(
            [AtomHandle::Id(AtomId(5)), AtomHandle::Id(AtomId(4))],
            NoncovalentBondForm::default(),
        )
    } else {
        NoncovalentBondHandle::Id(NoncovalentBondId(0))
    };
    let removal = if created {
        (
            [AtomHandle::Id(AtomId(5)), AtomHandle::Id(AtomId(4))],
            NoncovalentBondForm::default(),
        )
    } else {
        (
            [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
            NoncovalentBondForm::default(),
        )
    };
    prefix.remove_noncovalent_bonds(vec![(removed.clone(), removal.0, removal.1)]);
    let key = NoncovalentBondConstraintKey::Intramolecular;
    let form = NoncovalentBondConstraintForm::undetermined(key);
    let batch = if constraint {
        Edit::ModifyNoncovalentBondConstraints {
            key,
            changes: vec![
                (
                    NoncovalentBondHandle::Id(NoncovalentBondId(1)),
                    None,
                    Some(form.clone()),
                ),
                (removed, None, None),
            ],
        }
    } else {
        Edit::ModifyNoncovalentBondKinds {
            changes: vec![
                (
                    NoncovalentBondHandle::Id(NoncovalentBondId(1)),
                    NoncovalentBondKindForm::Undetermined,
                    NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                ),
                (
                    removed,
                    NoncovalentBondKindForm::Undetermined,
                    NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                ),
            ],
        }
    };
    let error = TransactionError::HandleRemoved {
        kind: EntityKind::NoncovalentBond,
        index: 0,
    };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    for edit in prefix.clone() {
        molecule.apply_edit(edit, &mut state).unwrap();
    }
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !created {
        expected.noncovalent.remove(0);
    }
    if !journaled {
        let attributes = &mut expected.noncovalent[usize::from(created)].1;
        if constraint {
            attributes.constraints = vec![form].into();
        } else {
            attributes.kind = NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond);
        }
    }
    assert_eq!(molecule, Molecule::from_entries(expected));
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        let edits = Edits::from_iter(prefix.clone().into_iter().chain([batch.clone()]));
        assert_eq!(
            molecule.transact([edits], recovery),
            Err(error.clone().into())
        );
        assert_eq!(molecule, source);
    }
}

#[rstest]
fn test_molecule_apply_stereo_atom_batches_removed(
    #[values(false, true)] created: bool,
    #[values(false, true)] constraint: bool,
    #[values(false, true)] journaled: bool,
) {
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 3],
        stereo_atoms: vec![
            (AtomId(0), vec![], StereoAtomForm::default()),
            (AtomId(1), vec![], StereoAtomForm::default()),
        ],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let mut prefix = Edits::new();
    let (removed, site) = if created {
        (
            prefix.add_stereo_atom(AtomHandle::Id(AtomId(2)), vec![], StereoAtomForm::default()),
            AtomHandle::Id(AtomId(2)),
        )
    } else {
        (
            StereoAtomHandle::Id(StereoAtomId(0)),
            AtomHandle::Id(AtomId(0)),
        )
    };
    prefix.remove_stereo_atoms(vec![(
        removed.clone(),
        site,
        vec![],
        StereoAtomForm::default(),
    )]);
    let key = StereoAtomConstraintKey::Stereogenicity;
    let form = StereoAtomConstraintForm::undetermined(key);
    let configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32);
    let batch = if constraint {
        Edit::ModifyStereoAtomConstraints {
            kind: None,
            key,
            changes: vec![
                (
                    StereoAtomHandle::Id(StereoAtomId(1)),
                    None,
                    Some(form.clone()),
                ),
                (removed, None, None),
            ],
        }
    } else {
        Edit::ModifyStereoAtomConfigurations {
            changes: vec![
                (
                    StereoAtomHandle::Id(StereoAtomId(1)),
                    StereoConfigurationForm::Undetermined,
                    configuration.clone(),
                ),
                (
                    removed,
                    StereoConfigurationForm::Undetermined,
                    configuration.clone(),
                ),
            ],
        }
    };
    let error = TransactionError::HandleRemoved {
        kind: EntityKind::StereoAtom,
        index: 0,
    };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    for edit in prefix.clone() {
        molecule.apply_edit(edit, &mut state).unwrap();
    }
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !created {
        expected.stereo_atoms.remove(0);
    }
    if !journaled {
        let attributes = &mut expected.stereo_atoms[usize::from(created)].2;
        if constraint {
            attributes.constraints = vec![form].into();
        } else {
            attributes.configuration = configuration;
        }
    }
    assert_eq!(molecule, Molecule::from_entries(expected));
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut restored = source.clone();
        let edits = Edits::from_iter(prefix.clone().into_iter().chain([batch.clone()]));
        assert_eq!(
            restored.transact([edits], recovery),
            Err(error.clone().into())
        );
        assert_eq!(restored, source);
    }
}

#[rstest]
fn test_molecule_apply_stereo_bond_batches_removed(
    stereo_bond_entries: MoleculeEntries,
    #[values(false, true)] created: bool,
    #[values(false, true)] constraint: bool,
    #[values(false, true)] journaled: bool,
) {
    let entries = stereo_bond_entries;
    let ligands = vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
        StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
    ];
    let source = Molecule::from_entries(entries.clone());
    let mut prefix = Edits::new();
    let (removed, site) = if created {
        (
            prefix.add_stereo_bond(
                BondHandle::Id(BondId(2)),
                ligands
                    .iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                StereoBondForm::default(),
            ),
            BondHandle::Id(BondId(2)),
        )
    } else {
        (
            StereoBondHandle::Id(StereoBondId(0)),
            BondHandle::Id(BondId(0)),
        )
    };
    let ligands = if created {
        ligands
    } else {
        entries.stereo_bonds[0].1.clone()
    };
    prefix.remove_stereo_bonds(vec![(
        removed.clone(),
        site,
        ligands
            .iter()
            .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
            .collect(),
        StereoBondForm::default(),
    )]);
    let key = StereoBondConstraintKey::Stereogenicity;
    let form = StereoBondConstraintForm::undetermined(key);
    let configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32);
    let batch = if constraint {
        Edit::ModifyStereoBondConstraints {
            kind: None,
            key,
            changes: vec![
                (
                    StereoBondHandle::Id(StereoBondId(1)),
                    None,
                    Some(form.clone()),
                ),
                (removed, None, None),
            ],
        }
    } else {
        Edit::ModifyStereoBondConfigurations {
            changes: vec![
                (
                    StereoBondHandle::Id(StereoBondId(1)),
                    StereoConfigurationForm::Undetermined,
                    configuration.clone(),
                ),
                (
                    removed,
                    StereoConfigurationForm::Undetermined,
                    configuration.clone(),
                ),
            ],
        }
    };
    let error = TransactionError::HandleRemoved {
        kind: EntityKind::StereoBond,
        index: 0,
    };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    for edit in prefix.clone() {
        molecule.apply_edit(edit, &mut state).unwrap();
    }
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !created {
        expected.stereo_bonds.remove(0);
    }
    if !journaled {
        let attributes = &mut expected.stereo_bonds[usize::from(created)].2;
        if constraint {
            attributes.constraints = vec![form].into();
        } else {
            attributes.configuration = configuration;
        }
    }
    assert_eq!(molecule, Molecule::from_entries(expected));
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut restored = source.clone();
        let edits = Edits::from_iter(prefix.clone().into_iter().chain([batch.clone()]));
        assert_eq!(
            restored.transact([edits], recovery),
            Err(error.clone().into())
        );
        assert_eq!(restored, source);
    }
}
