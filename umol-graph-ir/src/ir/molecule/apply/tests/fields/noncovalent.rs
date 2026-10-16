use super::*;

#[rstest]
#[case::kind(Edit::ModifyNoncovalentBondKinds { changes: vec![(NoncovalentBondHandle::Id(NoncovalentBondId(0)), NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)), (NoncovalentBondHandle::New(0), NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)), (NoncovalentBondHandle::Id(NoncovalentBondId(0)), NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic))] }, vec![Edit::ModifyNoncovalentBondField { id: NoncovalentBondHandle::Id(NoncovalentBondId(0)), change: NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindForm::Undetermined, new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond) } }, Edit::ModifyNoncovalentBondField { id: NoncovalentBondHandle::New(0), change: NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindForm::Undetermined, new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond) } }, Edit::ModifyNoncovalentBondField { id: NoncovalentBondHandle::Id(NoncovalentBondId(0)), change: NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic) } }], NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic), ..Default::default() }, NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), ..Default::default() })]
fn test_molecule_apply_noncovalent_bond_batches(
    #[case] batch: Edit,
    #[case] singles: Vec<Edit>,
    #[case] changed: NoncovalentBondForm,
    #[case] created: NoncovalentBondForm,
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
    prefix.modify_atom_charges(vec![(
        AtomHandle::Id(AtomId(0)),
        NumForm::Undetermined,
        NumForm::Lit(1),
    )]);
    prefix.add_noncovalent_bond(
        [AtomHandle::Id(AtomId(5)), AtomHandle::Id(AtomId(4))],
        NoncovalentBondForm::default(),
    );
    let batches = Edits::from_iter(prefix.clone().into_iter().chain([batch]));
    let singles = Edits::from_iter(prefix.into_iter().chain(singles));
    let mut expected = entries;
    expected.atoms[0].charge = NumForm::Lit(1);
    expected.noncovalent[0].1 = changed;
    expected.noncovalent.push(([AtomId(5), AtomId(4)], created));
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
#[case::kind(Edit::ModifyNoncovalentBondKinds { changes: (0..3).map(|id| (NoncovalentBondHandle::Id(NoncovalentBondId(id)), NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond))).collect() }, NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), ..Default::default() })]
fn test_molecule_apply_noncovalent_bond_batches_error(
    #[case] mut batch: Edit,
    #[case] changed: NoncovalentBondForm,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    match &mut batch {
        Edit::ModifyNoncovalentBondKinds { changes } => {
            if missing {
                changes[failure].0 = NoncovalentBondHandle::Id(NoncovalentBondId(3));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        _ => unreachable!(),
    }
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
    let prefix = Edit::ModifyAtomCharges {
        changes: vec![(
            AtomHandle::Id(AtomId(0)),
            NumForm::Undetermined,
            NumForm::Lit(1),
        )],
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
        for entry in &mut expected.noncovalent[..failure] {
            entry.1 = changed.clone();
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
#[case::kind(
    NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic) },
    NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
)]
fn test_molecule_apply_edit_modify_noncovalent_bond_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: NoncovalentBondFieldChange,
    #[case] expected: NoncovalentBondForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyNoncovalentBondField {
        id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
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
    removal_entries.noncovalent[0].1 = expected;
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
