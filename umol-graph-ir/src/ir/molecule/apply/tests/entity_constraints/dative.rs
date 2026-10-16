use super::*;

#[rstest]
#[case(
    DativeBondConstraintForm::aromatic(false),
    DativeBondConstraintForm::aromatic(true)
)]
#[case(
    DativeBondConstraintForm::ring_membership(RingScope::All, 1),
    DativeBondConstraintForm::ring_membership(RingScope::All, 2)
)]
#[case(
    DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1),
    DativeBondConstraintForm::ring_membership(RingScope::Size(5), 2)
)]
fn test_molecule_apply_dative_bond_constraint_batches(
    #[case] old: DativeBondConstraintForm,
    #[case] new: DativeBondConstraintForm,
) {
    let key = old.key();
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 4],
        dative: vec![(
            vec![AtomId(2), AtomId(0)],
            AtomId(1),
            DativeBondForm::default(),
        )],
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let rows = vec![
        (
            DativeBondHandle::Id(DativeBondId(0)),
            None,
            Some(old.clone()),
        ),
        (DativeBondHandle::New(0), None, Some(new.clone())),
        (
            DativeBondHandle::Id(DativeBondId(0)),
            Some(old),
            Some(new.clone()),
        ),
        (DativeBondHandle::Id(DativeBondId(0)), None, None),
        (
            DativeBondHandle::Id(DativeBondId(0)),
            Some(new.clone()),
            None,
        ),
    ];
    let mut prefix = Edits::new();
    prefix.add_dative_bond(
        vec![AtomHandle::Id(AtomId(3)), AtomHandle::Id(AtomId(1))],
        AtomHandle::Id(AtomId(2)),
        DativeBondForm::default(),
    );
    let mut batches = prefix.clone();
    batches.modify_dative_bond_constraints(key, rows.clone());
    let mut singles = prefix;
    for (id, old, new) in rows {
        singles.push(Edit::ModifyDativeBondConstraint { id, old, new });
    }
    let mut expected = entries;
    expected.dative.push((
        vec![AtomId(3), AtomId(1)],
        AtomId(2),
        DativeBondForm {
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
#[case(DativeBondConstraintKey::Aromatic)]
#[case(DativeBondConstraintKey::RingMembership(RingScope::All))]
#[case(DativeBondConstraintKey::RingMembership(RingScope::Size(5)))]
fn test_molecule_apply_dative_bond_constraint_batches_identity(
    #[case] key: DativeBondConstraintKey,
) {
    let edit = Edit::ModifyDativeBondConstraints {
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
        Undo::RestoreDativeBondConstraints {
            key,
            constraints: vec![]
        }
    );
    molecule.apply_undo(undo);
    assert_eq!(molecule, Molecule::new());
}

#[rstest]
#[case::old_state(
    Some(
        DativeBondConstraintForm::undetermined(DativeBondConstraintKey::RingMembership(
            RingScope::Size(5)
        ))
    ),
    None
)]
#[case::old_key(Some(DativeBondConstraintForm::aromatic(true)), None)]
#[case::new_key(None, Some(DativeBondConstraintForm::aromatic(true)))]
#[case::ring_scope(
    None,
    Some(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1))
)]
fn test_molecule_apply_dative_bond_constraint_batches_error(
    #[case] old: Option<DativeBondConstraintForm>,
    #[case] new: Option<DativeBondConstraintForm>,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    let key = DativeBondConstraintKey::RingMembership(RingScope::Size(5));
    let form = DativeBondConstraintForm::undetermined(key);
    let mut rows: Vec<_> = (0..3)
        .map(|i| {
            (
                DativeBondHandle::Id(DativeBondId(i)),
                None,
                Some(form.clone()),
            )
        })
        .collect();
    rows[failure].1 = old;
    rows[failure].2 = new;
    if missing {
        rows[failure].0 = DativeBondHandle::Id(DativeBondId(3));
    }
    let batch = Edit::ModifyDativeBondConstraints { key, changes: rows };
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
            (
                vec![AtomId(0), AtomId(2)],
                AtomId(3),
                DativeBondForm::default(),
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
            kind: EntityKind::DativeBond,
            index: 3,
            count: 3,
        }
    } else {
        TransactionError::OldStateMismatch
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    if !journaled {
        for entry in &mut expected.dative[..failure] {
            entry.2.constraints = vec![form.clone()].into();
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
