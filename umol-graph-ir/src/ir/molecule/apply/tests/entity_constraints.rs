use rstest::*;

use super::super::ApplicationState;
use super::removal_entries;
use crate::ir::atom::AtomForm;
use crate::ir::bond::BondForm;
use crate::ir::constraint::{
    AromaticSystemConstraintForm, AromaticSystemConstraintKey, AtomConstraintForm,
    AtomConstraintKey, BondConstraintForm, BondConstraintKey, DativeBondConstraintForm,
    DativeBondConstraintKey, MulticenterBondConstraintForm, MulticenterBondConstraintKey,
    NoncovalentBondConstraintForm, NoncovalentBondConstraintKey, StereoAtomConstraintForm,
    StereoAtomConstraintKey, StereoBondConstraintForm, StereoBondConstraintKey, StereogenicityForm,
};
use crate::ir::edit::{
    AromaticSystemHandle, AtomFieldChange, AtomHandle, BondHandle, DativeBondHandle, Edit, Edits,
    MulticenterBondHandle, NoncovalentBondHandle, StereoAtomHandle, StereoBondHandle, Undo,
};
use crate::ir::entity::EntityKind;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::molecule::MoleculeEntries;
use crate::ir::num::NumForm;
use crate::ir::stereo::{CisTransStereo, StereoKind, Stereogenicity};
use crate::ir::traits::Normalize;
use crate::ir::{
    Molecule, MoleculeApplyError, RingScope, RollbackMethod, Transaction, TransactionError,
};

#[rstest]
#[case(AtomConstraintKey::Valence)]
#[case(AtomConstraintKey::DonatedPairs)]
#[case(AtomConstraintKey::AcceptedPairs)]
#[case(AtomConstraintKey::AromaticValence)]
#[case(AtomConstraintKey::MulticenterValence)]
#[case(AtomConstraintKey::TetrahedralStereo)]
#[case(AtomConstraintKey::Degree)]
#[case(AtomConstraintKey::TotalDegree)]
#[case(AtomConstraintKey::TotalValence)]
#[case(AtomConstraintKey::RingDegree)]
#[case(AtomConstraintKey::RingValence)]
#[case(AtomConstraintKey::TotalHydrogens)]
#[case(AtomConstraintKey::RingMembership(RingScope::All))]
#[case(AtomConstraintKey::RingMembership(RingScope::Size(5)))]
fn test_molecule_apply_atom_constraint_batches(#[case] key: AtomConstraintKey) {
    let form = AtomConstraintForm::undetermined(key);
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default()],
        ..Default::default()
    });
    let rows = vec![
        (AtomHandle::Id(AtomId(0)), None, Some(form.clone())),
        (AtomHandle::New(0), None, Some(form.clone())),
        (
            AtomHandle::Id(AtomId(0)),
            Some(form.clone()),
            Some(form.clone()),
        ),
        (AtomHandle::Id(AtomId(0)), None, None),
        (AtomHandle::Id(AtomId(0)), Some(form.clone()), None),
    ];
    let mut edits = Edits::new();
    edits.add_atom(AtomForm::default());
    edits.modify_atom_constraints(key, rows.clone());
    let mut singles = Edits::new();
    singles.add_atom(AtomForm::default());
    for (id, old, new) in rows {
        singles.push(Edit::ModifyAtomConstraint { id, old, new });
    }
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::default(),
            AtomForm {
                constraints: vec![form].into(),
                ..Default::default()
            },
        ],
        ..Default::default()
    });
    let actual = source.clone().tracked_apply(edits.clone()).unwrap();
    assert_eq!(actual.0, expected);
    assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
    let mut molecule = source.clone();
    Transaction::run(
        &mut molecule,
        RollbackMethod::Journal,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits)?;
            transaction.rollback();
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(molecule, source);
}

#[rstest]
#[case(Edit::ModifyAtomConstraints { key: AtomConstraintKey::RingMembership(RingScope::Size(5)), changes: vec![] })]
fn test_molecule_apply_atom_constraint_batches_identity(#[case] edit: Edit) {
    let source = Molecule::new();
    assert_eq!(
        source
            .clone()
            .apply(Edits::from_iter([edit.clone()]))
            .unwrap(),
        source
    );
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(edit, &mut state)
        .unwrap()
        .unwrap();
    molecule.apply_undo(undo);
    assert_eq!(molecule, source);
}

#[rstest]
fn test_molecule_apply_atom_constraint_batches_normalized() {
    let accepted = AtomConstraintForm::valence(NumForm::LitSet(Box::new([1].into())));
    let mut molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm {
            constraints: vec![AtomConstraintForm::valence(1)].into(),
            ..Default::default()
        }],
        ..Default::default()
    });
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(
            Edit::ModifyAtomConstraints {
                key: AtomConstraintKey::Valence,
                changes: vec![(AtomHandle::Id(AtomId(0)), Some(accepted.clone()), None)],
            },
            &mut state,
        )
        .unwrap()
        .unwrap();
    assert_eq!(
        undo,
        Undo::RestoreAtomConstraints {
            key: AtomConstraintKey::Valence,
            constraints: vec![(AtomId(0), Some(accepted.clone()))]
        }
    );
    molecule.apply_undo(undo);
    assert_eq!(
        molecule,
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm {
                constraints: vec![accepted].into(),
                ..Default::default()
            }],
            ..Default::default()
        })
    );
}

#[rstest]
#[case::old_key(Some(AtomConstraintForm::degree(1)), None)]
#[case::new_key(None, Some(AtomConstraintForm::degree(1)))]
#[case::old_state(Some(AtomConstraintForm::valence(1)), None)]
#[case::ring_scope(None, Some(AtomConstraintForm::ring_membership(RingScope::Size(6), 1)))]
fn test_molecule_apply_atom_constraint_batches_error(
    #[case] old: Option<AtomConstraintForm>,
    #[case] new: Option<AtomConstraintForm>,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] journaled: bool,
) {
    let key = if matches!(&new, Some(AtomConstraintForm::RingMembership(_))) {
        AtomConstraintKey::RingMembership(RingScope::Size(5))
    } else {
        AtomConstraintKey::Valence
    };
    let form = AtomConstraintForm::undetermined(key);
    let mut rows: Vec<_> = (0..3)
        .map(|id| (AtomHandle::Id(AtomId(id)), None, Some(form.clone())))
        .collect();
    rows[failure].1 = old;
    rows[failure].2 = new;
    let batch = Edit::ModifyAtomConstraints { key, changes: rows };
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 3],
        ..Default::default()
    });
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(TransactionError::OldStateMismatch));
    let mut atoms = vec![AtomForm::default(); 3];
    if !journaled {
        for atom in &mut atoms[..failure] {
            atom.constraints = vec![form.clone()].into();
        }
    }
    assert_eq!(
        molecule,
        Molecule::from_entries(MoleculeEntries {
            atoms,
            ..Default::default()
        })
    );
    let prefix = Edit::ModifyAtomField {
        id: AtomHandle::Id(AtomId(0)),
        change: AtomFieldChange::Charge {
            old: NumForm::Undetermined,
            new: NumForm::Lit(1),
        },
    };
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        assert_eq!(
            molecule.transact(
                [Edits::from_iter([prefix.clone(), batch.clone()])],
                recovery
            ),
            Err(TransactionError::OldStateMismatch.into())
        );
        assert_eq!(molecule, source);
    }
}

#[rstest]
#[case::addition(None, None, Some(AtomConstraintForm::degree(3)))]
#[case::modification(
    Some(AtomConstraintForm::degree(2)),
    Some(AtomConstraintForm::degree(2)),
    Some(AtomConstraintForm::degree(3))
)]
#[case::removal(
    Some(AtomConstraintForm::degree(2)),
    Some(AtomConstraintForm::degree(2)),
    None
)]
#[case::absent(None, None, None)]
#[case::equivalent_old(Some(AtomConstraintForm::degree(2)), Some(AtomConstraintForm::degree(NumForm::lit_set([2]))), Some(AtomConstraintForm::degree(3)))]
fn test_molecule_apply_edit_modify_atom_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<AtomConstraintForm>,
    #[case] old: Option<AtomConstraintForm>,
    #[case] new: Option<AtomConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.atoms[0].constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyAtomConstraint {
        id: AtomHandle::Id(AtomId(0)),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreAtomConstraint {
            id: AtomId(0),
            key: AtomConstraintKey::Degree,
            constraint: old,
        })
    };
    let undo = if journaled {
        let undo = molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap();
        assert_eq!(undo, expected_undo);
        undo
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.atoms[0].constraints = new.into_iter().collect::<Vec<_>>().into();
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    if expected_undo.is_some() {
        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);
    }
    if let Some(undo) = undo {
        molecule.apply_undo(undo);
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::key(AtomConstraintForm::degree(2), AtomConstraintForm::valence(2))]
#[case::value(AtomConstraintForm::degree(3), AtomConstraintForm::degree(4))]
fn test_molecule_apply_edit_modify_atom_constraint_error(
    mut removal_entries: MoleculeEntries,
    #[case] old: AtomConstraintForm,
    #[case] new: AtomConstraintForm,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.atoms[0].constraints = vec![AtomConstraintForm::degree(2)].into();
    let initial = Molecule::from_entries(removal_entries);
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyAtomConstraint {
        id: AtomHandle::Id(AtomId(0)),
        old: Some(old),
        new: Some(new),
    };
    let result = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(result, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, initial);
}

#[rstest]
#[case(
    BondConstraintForm::aromatic(false),
    BondConstraintForm::aromatic(true)
)]
#[case(
    BondConstraintForm::undetermined(BondConstraintKey::CisTransStereo),
    BondConstraintForm::cis_trans_stereo(CisTransStereo::Stereo(0))
)]
#[case(
    BondConstraintForm::ring_membership(RingScope::All, 1),
    BondConstraintForm::ring_membership(RingScope::All, 2)
)]
#[case(
    BondConstraintForm::ring_membership(RingScope::Size(5), 1),
    BondConstraintForm::ring_membership(RingScope::Size(5), 2)
)]
fn test_molecule_apply_bond_constraint_batches(
    #[case] old: BondConstraintForm,
    #[case] new: BondConstraintForm,
) {
    let key = old.key();
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 3],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::default())],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let rows = vec![
        (BondHandle::Id(BondId(0)), None, Some(old.clone())),
        (BondHandle::New(0), None, Some(new.clone())),
        (BondHandle::Id(BondId(0)), Some(old), Some(new.clone())),
        (BondHandle::Id(BondId(0)), None, None),
        (BondHandle::Id(BondId(0)), Some(new.clone()), None),
    ];
    let mut prefix = Edits::new();
    prefix.add_bond(
        AtomHandle::Id(AtomId(1)),
        AtomHandle::Id(AtomId(2)),
        BondForm::default(),
    );
    let mut batches = prefix.clone();
    batches.modify_bond_constraints(key, rows.clone());
    let mut singles = prefix;
    for (id, old, new) in rows {
        singles.push(Edit::ModifyBondConstraint { id, old, new });
    }
    let mut expected = entries;
    expected.bonds.push((
        AtomId(1),
        AtomId(2),
        BondForm {
            constraints: vec![new].into(),
            ..Default::default()
        },
    ));
    let expected = Molecule::from_entries(expected);
    let actual = source.clone().tracked_apply(batches.clone()).unwrap();
    assert_eq!(actual.0, expected);
    assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        Transaction::run(
            &mut molecule,
            recovery,
            |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batches.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(molecule, source);
    }
}

#[rstest]
#[case(BondConstraintKey::Aromatic)]
#[case(BondConstraintKey::CisTransStereo)]
#[case(BondConstraintKey::RingMembership(RingScope::All))]
#[case(BondConstraintKey::RingMembership(RingScope::Size(5)))]
fn test_molecule_apply_bond_constraint_batches_identity(#[case] key: BondConstraintKey) {
    let edit = Edit::ModifyBondConstraints {
        key,
        changes: vec![],
    };
    let mut molecule = Molecule::new();
    let mut state = ApplicationState::new(&molecule);
    molecule.apply_edit(edit.clone(), &mut state).unwrap();
    let undo = molecule
        .apply_edit_with_undo(edit, &mut state)
        .unwrap()
        .unwrap();
    assert_eq!(
        undo,
        Undo::RestoreBondConstraints {
            key,
            constraints: vec![]
        }
    );
    molecule.apply_undo(undo);
    assert_eq!(molecule, Molecule::new());
}

#[rstest]
#[case::old_key(Some(BondConstraintForm::aromatic(true)), None)]
#[case::new_key(None, Some(BondConstraintForm::aromatic(true)))]
#[case::old_state(Some(BondConstraintForm::ring_membership(RingScope::Size(5), 1)), None)]
#[case::ring_scope(None, Some(BondConstraintForm::ring_membership(RingScope::Size(6), 1)))]
fn test_molecule_apply_bond_constraint_batches_error(
    #[case] old: Option<BondConstraintForm>,
    #[case] new: Option<BondConstraintForm>,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] journaled: bool,
) {
    let key = BondConstraintKey::RingMembership(RingScope::Size(5));
    let form = BondConstraintForm::undetermined(key);
    let mut rows: Vec<_> = (0..3)
        .map(|i| (BondHandle::Id(BondId(i)), None, Some(form.clone())))
        .collect();
    rows[failure].1 = old;
    rows[failure].2 = new;
    let batch = Edit::ModifyBondConstraints { key, changes: rows };
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 4],
        bonds: (0..3)
            .map(|i| (AtomId(i), AtomId(i + 1), BondForm::default()))
            .collect(),
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(TransactionError::OldStateMismatch));
    let mut expected = entries;
    if !journaled {
        for (_, _, attributes) in &mut expected.bonds[..failure] {
            attributes.constraints = vec![form.clone()].into();
        }
    }
    assert_eq!(molecule, Molecule::from_entries(expected));
    let prefix = Edit::ModifyAtomField {
        id: AtomHandle::Id(AtomId(0)),
        change: AtomFieldChange::Charge {
            old: NumForm::Undetermined,
            new: NumForm::Lit(1),
        },
    };
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        assert_eq!(
            molecule.transact(
                [Edits::from_iter([prefix.clone(), batch.clone()])],
                recovery
            ),
            Err(TransactionError::OldStateMismatch.into())
        );
        assert_eq!(molecule, source);
    }
}

#[rstest]
#[case::addition(None, None, Some(BondConstraintForm::aromatic(false)))]
#[case::modification(
    Some(BondConstraintForm::aromatic(true)),
    Some(BondConstraintForm::aromatic(true)),
    Some(BondConstraintForm::aromatic(false))
)]
#[case::removal(
    Some(BondConstraintForm::aromatic(true)),
    Some(BondConstraintForm::aromatic(true)),
    None
)]
#[case::absent(None, None, None)]
fn test_molecule_apply_edit_modify_bond_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<BondConstraintForm>,
    #[case] old: Option<BondConstraintForm>,
    #[case] new: Option<BondConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.bonds[0].2.constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyBondConstraint {
        id: BondHandle::Id(BondId(0)),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreBondConstraint {
            id: BondId(0),
            key: BondConstraintKey::Aromatic,
            constraint: old,
        })
    };
    let undo = if journaled {
        let undo = molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap();
        assert_eq!(undo, expected_undo);
        undo
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.bonds[0].2.constraints = new.into_iter().collect::<Vec<_>>().into();
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    if expected_undo.is_some() {
        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);
    }
    if let Some(undo) = undo {
        molecule.apply_undo(undo);
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::addition(None, None, Some(DativeBondConstraintForm::aromatic(false)))]
#[case::modification(
    Some(DativeBondConstraintForm::aromatic(true)),
    Some(DativeBondConstraintForm::aromatic(true)),
    Some(DativeBondConstraintForm::aromatic(false))
)]
#[case::removal(
    Some(DativeBondConstraintForm::aromatic(true)),
    Some(DativeBondConstraintForm::aromatic(true)),
    None
)]
#[case::absent(None, None, None)]
fn test_molecule_apply_edit_modify_dative_bond_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<DativeBondConstraintForm>,
    #[case] old: Option<DativeBondConstraintForm>,
    #[case] new: Option<DativeBondConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.dative[0].2.constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyDativeBondConstraint {
        id: DativeBondHandle::Id(DativeBondId(0)),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreDativeBondConstraint {
            id: DativeBondId(0),
            key: DativeBondConstraintKey::Aromatic,
            constraint: old,
        })
    };
    let undo = if journaled {
        let undo = molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap();
        assert_eq!(undo, expected_undo);
        undo
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.dative[0].2.constraints = new.into_iter().collect::<Vec<_>>().into();
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    if expected_undo.is_some() {
        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);
    }
    if let Some(undo) = undo {
        molecule.apply_undo(undo);
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::addition(None, None, Some(AromaticSystemConstraintForm::electron_count(4)))]
#[case::modification(
    Some(AromaticSystemConstraintForm::electron_count(2)),
    Some(AromaticSystemConstraintForm::electron_count(2)),
    Some(AromaticSystemConstraintForm::electron_count(4))
)]
#[case::removal(
    Some(AromaticSystemConstraintForm::electron_count(2)),
    Some(AromaticSystemConstraintForm::electron_count(2)),
    None
)]
#[case::absent(None, None, None)]
#[case::equivalent_old(Some(AromaticSystemConstraintForm::electron_count(2)), Some(AromaticSystemConstraintForm::electron_count(NumForm::lit_set([2]))), Some(AromaticSystemConstraintForm::electron_count(4)))]
fn test_molecule_apply_edit_modify_aromatic_system_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<AromaticSystemConstraintForm>,
    #[case] old: Option<AromaticSystemConstraintForm>,
    #[case] new: Option<AromaticSystemConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.aromatic[0].1.constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyAromaticSystemConstraint {
        id: AromaticSystemHandle::Id(AromaticSystemId(0)),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreAromaticSystemConstraint {
            id: AromaticSystemId(0),
            key: AromaticSystemConstraintKey::ElectronCount,
            constraint: old,
        })
    };
    let undo = if journaled {
        let undo = molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap();
        assert_eq!(undo, expected_undo);
        undo
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.aromatic[0].1.constraints = new.into_iter().collect::<Vec<_>>().into();
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    if expected_undo.is_some() {
        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);
    }
    if let Some(undo) = undo {
        molecule.apply_undo(undo);
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::addition(None, None, Some(MulticenterBondConstraintForm::electron_count(4)))]
#[case::modification(
    Some(MulticenterBondConstraintForm::electron_count(2)),
    Some(MulticenterBondConstraintForm::electron_count(2)),
    Some(MulticenterBondConstraintForm::electron_count(4))
)]
#[case::removal(
    Some(MulticenterBondConstraintForm::electron_count(2)),
    Some(MulticenterBondConstraintForm::electron_count(2)),
    None
)]
#[case::absent(None, None, None)]
#[case::equivalent_old(Some(MulticenterBondConstraintForm::electron_count(2)), Some(MulticenterBondConstraintForm::electron_count(NumForm::lit_set([2]))), Some(MulticenterBondConstraintForm::electron_count(4)))]
fn test_molecule_apply_edit_modify_multicenter_bond_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<MulticenterBondConstraintForm>,
    #[case] old: Option<MulticenterBondConstraintForm>,
    #[case] new: Option<MulticenterBondConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.multicenter[0].1.constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyMulticenterBondConstraint {
        id: MulticenterBondHandle::Id(MulticenterBondId(0)),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreMulticenterBondConstraint {
            id: MulticenterBondId(0),
            key: MulticenterBondConstraintKey::ElectronCount,
            constraint: old,
        })
    };
    let undo = if journaled {
        let undo = molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap();
        assert_eq!(undo, expected_undo);
        undo
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.multicenter[0].1.constraints = new.into_iter().collect::<Vec<_>>().into();
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    if expected_undo.is_some() {
        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);
    }
    if let Some(undo) = undo {
        molecule.apply_undo(undo);
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::addition(None, None, Some(NoncovalentBondConstraintForm::intramolecular(false)))]
#[case::modification(
    Some(NoncovalentBondConstraintForm::intramolecular(true)),
    Some(NoncovalentBondConstraintForm::intramolecular(true)),
    Some(NoncovalentBondConstraintForm::intramolecular(false))
)]
#[case::removal(
    Some(NoncovalentBondConstraintForm::intramolecular(true)),
    Some(NoncovalentBondConstraintForm::intramolecular(true)),
    None
)]
#[case::absent(None, None, None)]
fn test_molecule_apply_edit_modify_noncovalent_bond_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<NoncovalentBondConstraintForm>,
    #[case] old: Option<NoncovalentBondConstraintForm>,
    #[case] new: Option<NoncovalentBondConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.noncovalent[0].1.constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyNoncovalentBondConstraint {
        id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreNoncovalentBondConstraint {
            id: NoncovalentBondId(0),
            key: NoncovalentBondConstraintKey::Intramolecular,
            constraint: old,
        })
    };
    let undo = if journaled {
        let undo = molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap();
        assert_eq!(undo, expected_undo);
        undo
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.noncovalent[0].1.constraints = new.into_iter().collect::<Vec<_>>().into();
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    if expected_undo.is_some() {
        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);
    }
    if let Some(undo) = undo {
        molecule.apply_undo(undo);
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::addition(
    None,
    None,
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic
    )))
)]
#[case::modification(
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic
    )))
)]
#[case::removal(
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    None
)]
#[case::absent(None, None, None)]
fn test_molecule_apply_edit_modify_stereo_atom_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<StereoAtomConstraintForm>,
    #[case] old: Option<StereoAtomConstraintForm>,
    #[case] new: Option<StereoAtomConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.stereo_atoms[0].2.constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyStereoAtomConstraint {
        id: StereoAtomHandle::Id(StereoAtomId(0)),
        kind: Some(StereoKind::CisTrans),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreStereoAtomConstraint {
            id: StereoAtomId(0),
            key: StereoAtomConstraintKey::Stereogenicity,
            constraint: old,
        })
    };
    let undo = if journaled {
        let undo = molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap();
        assert_eq!(undo, expected_undo);
        undo
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.stereo_atoms[0].2.constraints = new.into_iter().collect::<Vec<_>>().into();
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    if expected_undo.is_some() {
        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);
    }
    if let Some(undo) = undo {
        molecule.apply_undo(undo);
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::addition(
    None,
    None,
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic
    )))
)]
#[case::modification(
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic
    )))
)]
#[case::removal(
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    None
)]
#[case::absent(None, None, None)]
fn test_molecule_apply_edit_modify_stereo_bond_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<StereoBondConstraintForm>,
    #[case] old: Option<StereoBondConstraintForm>,
    #[case] new: Option<StereoBondConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.stereo_bonds[0].2.constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyStereoBondConstraint {
        id: StereoBondHandle::Id(StereoBondId(0)),
        kind: Some(StereoKind::Tetrahedral),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreStereoBondConstraint {
            id: StereoBondId(0),
            key: StereoBondConstraintKey::Stereogenicity,
            constraint: old,
        })
    };
    let undo = if journaled {
        let undo = molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap();
        assert_eq!(undo, expected_undo);
        undo
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.stereo_bonds[0].2.constraints = new.into_iter().collect::<Vec<_>>().into();
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    if expected_undo.is_some() {
        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);
    }
    if let Some(undo) = undo {
        molecule.apply_undo(undo);
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::atom(Edit::ModifyAtomConstraint {
    id: AtomHandle::Id(AtomId(0)), old: None, new: None,
}, EntityKind::Atom)]
#[case::bond(Edit::ModifyBondConstraint {
    id: BondHandle::Id(BondId(0)), old: None, new: None,
}, EntityKind::Bond)]
#[case::dative_bond(Edit::ModifyDativeBondConstraint {
    id: DativeBondHandle::Id(DativeBondId(0)), old: None, new: None,
}, EntityKind::DativeBond)]
#[case::aromatic_system(Edit::ModifyAromaticSystemConstraint {
    id: AromaticSystemHandle::Id(AromaticSystemId(0)), old: None, new: None,
}, EntityKind::AromaticSystem)]
#[case::multicenter_bond(Edit::ModifyMulticenterBondConstraint {
    id: MulticenterBondHandle::Id(MulticenterBondId(0)), old: None, new: None,
}, EntityKind::MulticenterBond)]
#[case::noncovalent_bond(Edit::ModifyNoncovalentBondConstraint {
    id: NoncovalentBondHandle::Id(NoncovalentBondId(0)), old: None, new: None,
}, EntityKind::NoncovalentBond)]
#[case::stereo_atom(Edit::ModifyStereoAtomConstraint {
    id: StereoAtomHandle::Id(StereoAtomId(0)), kind: None, old: None, new: None,
}, EntityKind::StereoAtom)]
#[case::stereo_bond(Edit::ModifyStereoBondConstraint {
    id: StereoBondHandle::Id(StereoBondId(0)), kind: None, old: None, new: None,
}, EntityKind::StereoBond)]
fn test_molecule_apply_edit_constraint_target_error(
    #[case] edit: Edit,
    #[case] kind: EntityKind,
    #[values(false, true)] journaled: bool,
) {
    let mut molecule = Molecule::default();
    let mut state = ApplicationState::new(&molecule);
    let result = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(
        result,
        Err(TransactionError::HandleOutOfRange {
            kind,
            index: 0,
            count: 0
        })
    );
    assert_eq!(molecule, Molecule::default());
}
