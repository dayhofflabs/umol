use super::*;

#[rstest]
#[case(MulticenterBondConstraintForm::electron_count(1), MulticenterBondConstraintForm::electron_count(NumForm::LitSet(Box::new([1, 2].into()))))]
fn test_molecule_apply_multicenter_bond_constraint_batches(
    #[case] old: MulticenterBondConstraintForm,
    #[case] new: MulticenterBondConstraintForm,
) {
    let key = old.key();
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 4],
        multicenter: vec![(
            vec![AtomId(2), AtomId(0), AtomId(1)],
            MulticenterBondForm::default(),
        )],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let rows = vec![
        (
            MulticenterBondHandle::Id(MulticenterBondId(0)),
            None,
            Some(old.clone()),
        ),
        (MulticenterBondHandle::New(0), None, Some(new.clone())),
        (
            MulticenterBondHandle::Id(MulticenterBondId(0)),
            Some(old),
            Some(new.clone()),
        ),
        (MulticenterBondHandle::Id(MulticenterBondId(0)), None, None),
        (
            MulticenterBondHandle::Id(MulticenterBondId(0)),
            Some(new.clone()),
            None,
        ),
    ];
    let mut prefix = Edits::new();
    prefix.add_multicenter_bond(
        vec![
            AtomHandle::Id(AtomId(3)),
            AtomHandle::Id(AtomId(1)),
            AtomHandle::Id(AtomId(2)),
        ],
        MulticenterBondForm::default(),
    );
    let mut batches = prefix.clone();
    batches.modify_multicenter_bond_constraints(key, rows.clone());
    let mut singles = prefix;
    for (id, old, new) in rows {
        singles.push(Edit::ModifyMulticenterBondConstraint { id, old, new });
    }
    let mut expected = entries;
    expected.multicenter.push((
        vec![AtomId(3), AtomId(1), AtomId(2)],
        MulticenterBondForm {
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
#[case(MulticenterBondConstraintKey::ElectronCount)]
fn test_molecule_apply_multicenter_bond_constraint_batches_identity(
    #[case] key: MulticenterBondConstraintKey,
) {
    let edit = Edit::ModifyMulticenterBondConstraints {
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
        Undo::RestoreMulticenterBondConstraints {
            key,
            constraints: vec![]
        }
    );
    molecule.apply_undo(undo);
    assert_eq!(molecule, Molecule::new());
}

#[rstest]
#[case::old_state(
    Some(MulticenterBondConstraintForm::undetermined(
        MulticenterBondConstraintKey::ElectronCount
    )),
    None
)]
fn test_molecule_apply_multicenter_bond_constraint_batches_error(
    #[case] old: Option<MulticenterBondConstraintForm>,
    #[case] new: Option<MulticenterBondConstraintForm>,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    let key = MulticenterBondConstraintKey::ElectronCount;
    let form = MulticenterBondConstraintForm::undetermined(key);
    let mut rows: Vec<_> = (0..3)
        .map(|i| {
            (
                MulticenterBondHandle::Id(MulticenterBondId(i)),
                None,
                Some(form.clone()),
            )
        })
        .collect();
    rows[failure].1 = old;
    rows[failure].2 = new;
    if missing {
        rows[failure].0 = MulticenterBondHandle::Id(MulticenterBondId(3));
    }
    let batch = Edit::ModifyMulticenterBondConstraints { key, changes: rows };
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
            (
                vec![AtomId(0), AtomId(2), AtomId(3)],
                MulticenterBondForm::default(),
            ),
        ],
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
    let error = if missing {
        TransactionError::HandleOutOfRange {
            kind: EntityKind::MulticenterBond,
            index: 3,
            count: 3,
        }
    } else {
        TransactionError::OldStateMismatch
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !journaled {
        for entry in &mut expected.multicenter[..failure] {
            entry.1.constraints = vec![form.clone()].into();
        }
    }
    assert_eq!(molecule, Molecule::from_entries(expected));
    let prefix = Edit::ModifyAtomCharges {
        changes: vec![(
            AtomHandle::Id(AtomId(0)),
            NumForm::Undetermined,
            NumForm::Lit(1),
        )],
    };
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        assert_eq!(
            molecule.transact(
                [Edits::from_iter([prefix.clone(), batch.clone()])],
                recovery
            ),
            Err(error.clone().into())
        );
        assert_eq!(molecule, source);
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
