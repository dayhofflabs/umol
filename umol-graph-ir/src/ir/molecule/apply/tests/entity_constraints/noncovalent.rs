use super::*;

#[rstest]
#[case(
    NoncovalentBondConstraintForm::intramolecular(false),
    NoncovalentBondConstraintForm::intramolecular(true)
)]
fn test_molecule_apply_noncovalent_bond_constraint_batches(
    #[case] old: NoncovalentBondConstraintForm,
    #[case] new: NoncovalentBondConstraintForm,
) {
    let key = old.key();
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 6],
        noncovalent: vec![([AtomId(1), AtomId(0)], NoncovalentBondForm::default())],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let rows = vec![
        (
            NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            None,
            Some(old.clone()),
        ),
        (NoncovalentBondHandle::New(0), None, Some(new.clone())),
        (
            NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            Some(old),
            Some(new.clone()),
        ),
        (NoncovalentBondHandle::Id(NoncovalentBondId(0)), None, None),
        (
            NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            Some(new.clone()),
            None,
        ),
    ];
    let mut prefix = Edits::new();
    prefix.add_noncovalent_bond(
        [AtomHandle::Id(AtomId(3)), AtomHandle::Id(AtomId(2))],
        NoncovalentBondForm::default(),
    );
    let mut batches = prefix.clone();
    batches.modify_noncovalent_bond_constraints(key, rows.clone());
    let mut singles = prefix;
    for (id, old, new) in rows {
        singles.push(Edit::ModifyNoncovalentBondConstraint { id, old, new });
    }
    let mut expected = entries;
    expected.noncovalent.push((
        [AtomId(3), AtomId(2)],
        NoncovalentBondForm {
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
#[case(NoncovalentBondConstraintKey::Intramolecular)]
fn test_molecule_apply_noncovalent_bond_constraint_batches_identity(
    #[case] key: NoncovalentBondConstraintKey,
) {
    let edit = Edit::ModifyNoncovalentBondConstraints {
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
        Undo::RestoreNoncovalentBondConstraints {
            key,
            constraints: vec![]
        }
    );
    molecule.apply_undo(undo);
    assert_eq!(molecule, Molecule::new());
}

#[rstest]
#[case::old_state(
    Some(NoncovalentBondConstraintForm::undetermined(
        NoncovalentBondConstraintKey::Intramolecular
    )),
    None
)]
fn test_molecule_apply_noncovalent_bond_constraint_batches_error(
    #[case] old: Option<NoncovalentBondConstraintForm>,
    #[case] new: Option<NoncovalentBondConstraintForm>,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    let key = NoncovalentBondConstraintKey::Intramolecular;
    let form = NoncovalentBondConstraintForm::undetermined(key);
    let mut rows: Vec<_> = (0..3)
        .map(|i| {
            (
                NoncovalentBondHandle::Id(NoncovalentBondId(i)),
                None,
                Some(form.clone()),
            )
        })
        .collect();
    rows[failure].1 = old;
    rows[failure].2 = new;
    if missing {
        rows[failure].0 = NoncovalentBondHandle::Id(NoncovalentBondId(3));
    }
    let batch = Edit::ModifyNoncovalentBondConstraints { key, changes: rows };
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 6],
        noncovalent: vec![
            ([AtomId(1), AtomId(0)], NoncovalentBondForm::default()),
            ([AtomId(3), AtomId(2)], NoncovalentBondForm::default()),
            ([AtomId(5), AtomId(4)], NoncovalentBondForm::default()),
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
            kind: EntityKind::NoncovalentBond,
            index: 3,
            count: 3,
        }
    } else {
        TransactionError::OldStateMismatch
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !journaled {
        for entry in &mut expected.noncovalent[..failure] {
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
