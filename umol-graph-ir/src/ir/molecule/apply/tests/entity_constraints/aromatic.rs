use super::*;

#[rstest]
#[case(AromaticSystemConstraintForm::electron_count(1), AromaticSystemConstraintForm::electron_count(NumForm::LitSet(Box::new([1, 2].into()))))]
fn test_molecule_apply_aromatic_system_constraint_batches(
    #[case] old: AromaticSystemConstraintForm,
    #[case] new: AromaticSystemConstraintForm,
) {
    let key = old.key();
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 9],
        aromatic: vec![(
            vec![AtomId(2), AtomId(0), AtomId(1)],
            AromaticSystemForm::default(),
        )],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let rows = vec![
        (
            AromaticSystemHandle::Id(AromaticSystemId(0)),
            None,
            Some(old.clone()),
        ),
        (AromaticSystemHandle::New(0), None, Some(new.clone())),
        (
            AromaticSystemHandle::Id(AromaticSystemId(0)),
            Some(old),
            Some(new.clone()),
        ),
        (AromaticSystemHandle::Id(AromaticSystemId(0)), None, None),
        (
            AromaticSystemHandle::Id(AromaticSystemId(0)),
            Some(new.clone()),
            None,
        ),
    ];
    let mut prefix = Edits::new();
    prefix.add_aromatic_system(
        vec![
            AtomHandle::Id(AtomId(5)),
            AtomHandle::Id(AtomId(3)),
            AtomHandle::Id(AtomId(4)),
        ],
        AromaticSystemForm::default(),
    );
    let mut batches = prefix.clone();
    batches.modify_aromatic_system_constraints(key, rows.clone());
    let mut singles = prefix;
    for (id, old, new) in rows {
        singles.push(Edit::ModifyAromaticSystemConstraint { id, old, new });
    }
    let mut expected = entries;
    expected.aromatic.push((
        vec![AtomId(5), AtomId(3), AtomId(4)],
        AromaticSystemForm {
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
#[case(AromaticSystemConstraintKey::ElectronCount)]
fn test_molecule_apply_aromatic_system_constraint_batches_identity(
    #[case] key: AromaticSystemConstraintKey,
) {
    let edit = Edit::ModifyAromaticSystemConstraints {
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
        Undo::RestoreAromaticSystemConstraints {
            key,
            constraints: vec![]
        }
    );
    molecule.apply_undo(undo);
    assert_eq!(molecule, Molecule::new());
}

#[rstest]
#[case::old_state(
    Some(AromaticSystemConstraintForm::undetermined(AromaticSystemConstraintKey::ElectronCount)),
    None
)]
fn test_molecule_apply_aromatic_system_constraint_batches_error(
    #[case] old: Option<AromaticSystemConstraintForm>,
    #[case] new: Option<AromaticSystemConstraintForm>,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    let key = AromaticSystemConstraintKey::ElectronCount;
    let form = AromaticSystemConstraintForm::undetermined(key);
    let mut rows: Vec<_> = (0..3)
        .map(|i| {
            (
                AromaticSystemHandle::Id(AromaticSystemId(i)),
                None,
                Some(form.clone()),
            )
        })
        .collect();
    rows[failure].1 = old;
    rows[failure].2 = new;
    if missing {
        rows[failure].0 = AromaticSystemHandle::Id(AromaticSystemId(3));
    }
    let batch = Edit::ModifyAromaticSystemConstraints { key, changes: rows };
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
            (
                vec![AtomId(8), AtomId(6), AtomId(7)],
                AromaticSystemForm::default(),
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
            kind: EntityKind::AromaticSystem,
            index: 3,
            count: 3,
        }
    } else {
        TransactionError::OldStateMismatch
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !journaled {
        for entry in &mut expected.aromatic[..failure] {
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
