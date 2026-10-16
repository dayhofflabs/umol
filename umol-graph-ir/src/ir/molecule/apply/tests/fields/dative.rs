use super::*;

#[rstest]
#[case::order(Edit::ModifyDativeBondOrders { changes: vec![(DativeBondHandle::Id(DativeBondId(0)), NumForm::Undetermined, NumForm::Lit(1)), (DativeBondHandle::New(0), NumForm::Undetermined, NumForm::Lit(1)), (DativeBondHandle::Id(DativeBondId(0)), NumForm::Lit(1), NumForm::LitSet(Box::new([1, 2].into())))] }, vec![Edit::ModifyDativeBondField { id: DativeBondHandle::Id(DativeBondId(0)), change: DativeBondFieldChange::Order { old: NumForm::Undetermined, new: NumForm::Lit(1) } }, Edit::ModifyDativeBondField { id: DativeBondHandle::New(0), change: DativeBondFieldChange::Order { old: NumForm::Undetermined, new: NumForm::Lit(1) } }, Edit::ModifyDativeBondField { id: DativeBondHandle::Id(DativeBondId(0)), change: DativeBondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::LitSet(Box::new([1, 2].into())) } }], DativeBondForm { order: NumForm::LitSet(Box::new([1, 2].into())), ..Default::default() }, DativeBondForm { order: NumForm::Lit(1), ..Default::default() })]
fn test_molecule_apply_dative_bond_batches(
    #[case] batch: Edit,
    #[case] singles: Vec<Edit>,
    #[case] changed: DativeBondForm,
    #[case] created: DativeBondForm,
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
    prefix.modify_atom_charges(vec![(
        AtomHandle::Id(AtomId(0)),
        NumForm::Undetermined,
        NumForm::Lit(1),
    )]);
    prefix.add_dative_bond(
        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(2))],
        AtomHandle::Id(AtomId(3)),
        DativeBondForm::default(),
    );
    let batches = Edits::from_iter(prefix.clone().into_iter().chain([batch]));
    let singles = Edits::from_iter(prefix.into_iter().chain(singles));
    let mut expected = entries;
    expected.atoms[0].charge = NumForm::Lit(1);
    expected.dative[0].2 = changed;
    expected
        .dative
        .push((vec![AtomId(0), AtomId(2)], AtomId(3), created));
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
#[case::order(Edit::ModifyDativeBondOrders { changes: (0..3).map(|id| (DativeBondHandle::Id(DativeBondId(id)), NumForm::Undetermined, NumForm::Lit(1))).collect() }, DativeBondForm { order: NumForm::Lit(1), ..Default::default() })]
fn test_molecule_apply_dative_bond_batches_error(
    #[case] mut batch: Edit,
    #[case] changed: DativeBondForm,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    match &mut batch {
        Edit::ModifyDativeBondOrders { changes } => {
            if missing {
                changes[failure].0 = DativeBondHandle::Id(DativeBondId(3));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        _ => unreachable!(),
    }
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
    let prefix = Edit::ModifyAtomCharges {
        changes: vec![(
            AtomHandle::Id(AtomId(0)),
            NumForm::Undetermined,
            NumForm::Lit(1),
        )],
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
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    molecule.apply_edit(prefix.clone(), &mut state).unwrap();
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut expected = entries;
    expected.atoms[0].charge = NumForm::Lit(1);
    if !journaled {
        for entry in &mut expected.dative[..failure] {
            entry.2 = changed.clone();
        }
    }
    assert_eq!(molecule, Molecule::from_entries(expected));
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
#[case::order(
    DativeBondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
    DativeBondForm::from_order(2),
)]
#[case::order_equivalent(
    DativeBondFieldChange::Order { old: NumForm::lit_set([1]), new: NumForm::Lit(2) },
    DativeBondForm::from_order(2),
)]
fn test_molecule_apply_edit_modify_dative_bond_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: DativeBondFieldChange,
    #[case] expected: DativeBondForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyDativeBondField {
        id: DativeBondHandle::Id(DativeBondId(0)),
        change,
    };
    let undo = if journaled {
        molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap()
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.dative[0].2 = expected;
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    let repeated = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, expected);

    if journaled {
        molecule.apply_undo(undo.unwrap());
        assert!(molecule.normalized_eq(&initial));
    }
}
