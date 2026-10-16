use super::*;

#[rstest]
#[case::electrons(Edit::ModifyAromaticSystemElectrons { changes: vec![(AromaticSystemHandle::Id(AromaticSystemId(0)), ElectronCountsForm::Undetermined, ElectronCountsForm::Lit(vec![1, 2])), (AromaticSystemHandle::New(0), ElectronCountsForm::Undetermined, ElectronCountsForm::Lit(vec![1, 2])), (AromaticSystemHandle::Id(AromaticSystemId(0)), ElectronCountsForm::Lit(vec![1, 2]), ElectronCountsForm::Lit(vec![]))] }, vec![Edit::ModifyAromaticSystemField { id: AromaticSystemHandle::Id(AromaticSystemId(0)), change: AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Undetermined, new: ElectronCountsForm::Lit(vec![1, 2]) } }, Edit::ModifyAromaticSystemField { id: AromaticSystemHandle::New(0), change: AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Undetermined, new: ElectronCountsForm::Lit(vec![1, 2]) } }, Edit::ModifyAromaticSystemField { id: AromaticSystemHandle::Id(AromaticSystemId(0)), change: AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![]) } }], AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![]), ..Default::default() }, AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 2]), ..Default::default() })]
#[case::charge(Edit::ModifyAromaticSystemCharges { changes: vec![(AromaticSystemHandle::Id(AromaticSystemId(0)), NumForm::Undetermined, NumForm::Lit(1)), (AromaticSystemHandle::New(0), NumForm::Undetermined, NumForm::Lit(1)), (AromaticSystemHandle::Id(AromaticSystemId(0)), NumForm::Lit(1), NumForm::LitSet(Box::new([1, 2].into())))] }, vec![Edit::ModifyAromaticSystemField { id: AromaticSystemHandle::Id(AromaticSystemId(0)), change: AromaticSystemFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) } }, Edit::ModifyAromaticSystemField { id: AromaticSystemHandle::New(0), change: AromaticSystemFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) } }, Edit::ModifyAromaticSystemField { id: AromaticSystemHandle::Id(AromaticSystemId(0)), change: AromaticSystemFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::LitSet(Box::new([1, 2].into())) } }], AromaticSystemForm { charge: NumForm::LitSet(Box::new([1, 2].into())), ..Default::default() }, AromaticSystemForm { charge: NumForm::Lit(1), ..Default::default() })]
#[case::unpaired_electrons(Edit::ModifyAromaticSystemUnpairedElectrons { changes: vec![(AromaticSystemHandle::Id(AromaticSystemId(0)), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell()), (AromaticSystemHandle::New(0), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell()), (AromaticSystemHandle::Id(AromaticSystemId(0)), UnpairedElectronsForm::closed_shell(), UnpairedElectronsForm::from((1_u8, 2_u8)))] }, vec![Edit::ModifyAromaticSystemField { id: AromaticSystemHandle::Id(AromaticSystemId(0)), change: AromaticSystemFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm::closed_shell() } }, Edit::ModifyAromaticSystemField { id: AromaticSystemHandle::New(0), change: AromaticSystemFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm::closed_shell() } }, Edit::ModifyAromaticSystemField { id: AromaticSystemHandle::Id(AromaticSystemId(0)), change: AromaticSystemFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::closed_shell(), new: UnpairedElectronsForm::from((1_u8, 2_u8)) } }], AromaticSystemForm { unpaired_electrons: UnpairedElectronsForm::from((1_u8, 2_u8)), ..Default::default() }, AromaticSystemForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })]
fn test_molecule_apply_aromatic_system_batches(
    #[case] batch: Edit,
    #[case] singles: Vec<Edit>,
    #[case] changed: AromaticSystemForm,
    #[case] created: AromaticSystemForm,
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
    prefix.modify_atom_charges(vec![(
        AtomHandle::Id(AtomId(0)),
        NumForm::Undetermined,
        NumForm::Lit(1),
    )]);
    prefix.add_aromatic_system(
        vec![
            AtomHandle::Id(AtomId(8)),
            AtomHandle::Id(AtomId(6)),
            AtomHandle::Id(AtomId(7)),
        ],
        AromaticSystemForm::default(),
    );
    let batches = Edits::from_iter(prefix.clone().into_iter().chain([batch]));
    let singles = Edits::from_iter(prefix.into_iter().chain(singles));
    let mut expected = entries;
    expected.atoms[0].charge = NumForm::Lit(1);
    expected.aromatic[0].1 = changed;
    expected
        .aromatic
        .push((vec![AtomId(8), AtomId(6), AtomId(7)], created));
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
#[case::electrons(Edit::ModifyAromaticSystemElectrons { changes: (0..3).map(|id| (AromaticSystemHandle::Id(AromaticSystemId(id)), ElectronCountsForm::Undetermined, ElectronCountsForm::Lit(vec![1, 2]))).collect() }, AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 2]), ..Default::default() })]
#[case::charge(Edit::ModifyAromaticSystemCharges { changes: (0..3).map(|id| (AromaticSystemHandle::Id(AromaticSystemId(id)), NumForm::Undetermined, NumForm::Lit(1))).collect() }, AromaticSystemForm { charge: NumForm::Lit(1), ..Default::default() })]
#[case::unpaired_electrons(Edit::ModifyAromaticSystemUnpairedElectrons { changes: (0..3).map(|id| (AromaticSystemHandle::Id(AromaticSystemId(id)), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell())).collect() }, AromaticSystemForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })]
fn test_molecule_apply_aromatic_system_batches_error(
    #[case] mut batch: Edit,
    #[case] changed: AromaticSystemForm,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    match &mut batch {
        Edit::ModifyAromaticSystemElectrons { changes } => {
            if missing {
                changes[failure].0 = AromaticSystemHandle::Id(AromaticSystemId(3));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        Edit::ModifyAromaticSystemCharges { changes } => {
            if missing {
                changes[failure].0 = AromaticSystemHandle::Id(AromaticSystemId(3));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        Edit::ModifyAromaticSystemUnpairedElectrons { changes } => {
            if missing {
                changes[failure].0 = AromaticSystemHandle::Id(AromaticSystemId(3));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        _ => unreachable!(),
    }
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
    let prefix = Edit::ModifyAtomCharges {
        changes: vec![(
            AtomHandle::Id(AtomId(0)),
            NumForm::Undetermined,
            NumForm::Lit(1),
        )],
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
        for entry in &mut expected.aromatic[..failure] {
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
#[case::electrons(
    AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2, 2]) },
    AromaticSystemForm::from_electrons(vec![2, 2]),
)]
#[case::electron_count_length(
    AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2]) },
    AromaticSystemForm::from_electrons(vec![2]),
)]
#[case::charge(
    AromaticSystemFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) },
    AromaticSystemForm { charge: NumForm::Lit(1), ..AromaticSystemForm::from_electrons(vec![1, 2]) },
)]
#[case::unpaired_electrons(
    AromaticSystemFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
    AromaticSystemForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..AromaticSystemForm::from_electrons(vec![1, 2]) },
)]
fn test_molecule_apply_edit_modify_aromatic_system_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: AromaticSystemFieldChange,
    #[case] expected: AromaticSystemForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyAromaticSystemField {
        id: AromaticSystemHandle::Id(AromaticSystemId(0)),
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
    removal_entries.aromatic[0].1 = expected;
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
