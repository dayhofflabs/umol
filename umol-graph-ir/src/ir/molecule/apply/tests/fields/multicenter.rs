use super::*;

#[rstest]
#[case::electrons(Edit::ModifyMulticenterBondElectrons { changes: vec![(MulticenterBondHandle::Id(MulticenterBondId(0)), ElectronCountsForm::Undetermined, ElectronCountsForm::Lit(vec![1, 2])), (MulticenterBondHandle::New(0), ElectronCountsForm::Undetermined, ElectronCountsForm::Lit(vec![1, 2])), (MulticenterBondHandle::Id(MulticenterBondId(0)), ElectronCountsForm::Lit(vec![1, 2]), ElectronCountsForm::Lit(vec![]))] }, vec![Edit::ModifyMulticenterBondField { id: MulticenterBondHandle::Id(MulticenterBondId(0)), change: MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Undetermined, new: ElectronCountsForm::Lit(vec![1, 2]) } }, Edit::ModifyMulticenterBondField { id: MulticenterBondHandle::New(0), change: MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Undetermined, new: ElectronCountsForm::Lit(vec![1, 2]) } }, Edit::ModifyMulticenterBondField { id: MulticenterBondHandle::Id(MulticenterBondId(0)), change: MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![]) } }], MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![]), ..Default::default() }, MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 2]), ..Default::default() })]
#[case::charge(Edit::ModifyMulticenterBondCharges { changes: vec![(MulticenterBondHandle::Id(MulticenterBondId(0)), NumForm::Undetermined, NumForm::Lit(1)), (MulticenterBondHandle::New(0), NumForm::Undetermined, NumForm::Lit(1)), (MulticenterBondHandle::Id(MulticenterBondId(0)), NumForm::Lit(1), NumForm::LitSet(Box::new([1, 2].into())))] }, vec![Edit::ModifyMulticenterBondField { id: MulticenterBondHandle::Id(MulticenterBondId(0)), change: MulticenterBondFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) } }, Edit::ModifyMulticenterBondField { id: MulticenterBondHandle::New(0), change: MulticenterBondFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) } }, Edit::ModifyMulticenterBondField { id: MulticenterBondHandle::Id(MulticenterBondId(0)), change: MulticenterBondFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::LitSet(Box::new([1, 2].into())) } }], MulticenterBondForm { charge: NumForm::LitSet(Box::new([1, 2].into())), ..Default::default() }, MulticenterBondForm { charge: NumForm::Lit(1), ..Default::default() })]
#[case::unpaired_electrons(Edit::ModifyMulticenterBondUnpairedElectrons { changes: vec![(MulticenterBondHandle::Id(MulticenterBondId(0)), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell()), (MulticenterBondHandle::New(0), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell()), (MulticenterBondHandle::Id(MulticenterBondId(0)), UnpairedElectronsForm::closed_shell(), UnpairedElectronsForm::from((1_u8, 2_u8)))] }, vec![Edit::ModifyMulticenterBondField { id: MulticenterBondHandle::Id(MulticenterBondId(0)), change: MulticenterBondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm::closed_shell() } }, Edit::ModifyMulticenterBondField { id: MulticenterBondHandle::New(0), change: MulticenterBondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm::closed_shell() } }, Edit::ModifyMulticenterBondField { id: MulticenterBondHandle::Id(MulticenterBondId(0)), change: MulticenterBondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::closed_shell(), new: UnpairedElectronsForm::from((1_u8, 2_u8)) } }], MulticenterBondForm { unpaired_electrons: UnpairedElectronsForm::from((1_u8, 2_u8)), ..Default::default() }, MulticenterBondForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })]
fn test_molecule_apply_multicenter_bond_batches(
    #[case] batch: Edit,
    #[case] singles: Vec<Edit>,
    #[case] changed: MulticenterBondForm,
    #[case] created: MulticenterBondForm,
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
    prefix.modify_atom_charges(vec![(
        AtomHandle::Id(AtomId(0)),
        NumForm::Undetermined,
        NumForm::Lit(1),
    )]);
    prefix.add_multicenter_bond(
        vec![
            AtomHandle::Id(AtomId(0)),
            AtomHandle::Id(AtomId(2)),
            AtomHandle::Id(AtomId(3)),
        ],
        MulticenterBondForm::default(),
    );
    let batches = Edits::from_iter(prefix.clone().into_iter().chain([batch]));
    let singles = Edits::from_iter(prefix.into_iter().chain(singles));
    let mut expected = entries;
    expected.atoms[0].charge = NumForm::Lit(1);
    expected.multicenter[0].1 = changed;
    expected
        .multicenter
        .push((vec![AtomId(0), AtomId(2), AtomId(3)], created));
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
#[case::electrons(Edit::ModifyMulticenterBondElectrons { changes: (0..3).map(|id| (MulticenterBondHandle::Id(MulticenterBondId(id)), ElectronCountsForm::Undetermined, ElectronCountsForm::Lit(vec![1, 2]))).collect() }, MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 2]), ..Default::default() })]
#[case::charge(Edit::ModifyMulticenterBondCharges { changes: (0..3).map(|id| (MulticenterBondHandle::Id(MulticenterBondId(id)), NumForm::Undetermined, NumForm::Lit(1))).collect() }, MulticenterBondForm { charge: NumForm::Lit(1), ..Default::default() })]
#[case::unpaired_electrons(Edit::ModifyMulticenterBondUnpairedElectrons { changes: (0..3).map(|id| (MulticenterBondHandle::Id(MulticenterBondId(id)), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell())).collect() }, MulticenterBondForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })]
fn test_molecule_apply_multicenter_bond_batches_error(
    #[case] mut batch: Edit,
    #[case] changed: MulticenterBondForm,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    match &mut batch {
        Edit::ModifyMulticenterBondElectrons { changes } => {
            if missing {
                changes[failure].0 = MulticenterBondHandle::Id(MulticenterBondId(3));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        Edit::ModifyMulticenterBondCharges { changes } => {
            if missing {
                changes[failure].0 = MulticenterBondHandle::Id(MulticenterBondId(3));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        Edit::ModifyMulticenterBondUnpairedElectrons { changes } => {
            if missing {
                changes[failure].0 = MulticenterBondHandle::Id(MulticenterBondId(3));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        _ => unreachable!(),
    }
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
    let prefix = Edit::ModifyAtomCharges {
        changes: vec![(
            AtomHandle::Id(AtomId(0)),
            NumForm::Undetermined,
            NumForm::Lit(1),
        )],
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
        for entry in &mut expected.multicenter[..failure] {
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
    MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2, 2]) },
    MulticenterBondForm::from_electrons(vec![2, 2]),
)]
#[case::electron_count_length(
    MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2, 1, 2]) },
    MulticenterBondForm::from_electrons(vec![2, 1, 2]),
)]
#[case::charge(
    MulticenterBondFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(-1) },
    MulticenterBondForm { charge: NumForm::Lit(-1), ..MulticenterBondForm::from_electrons(vec![1, 2]) },
)]
#[case::unpaired_electrons(
    MulticenterBondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
    MulticenterBondForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..MulticenterBondForm::from_electrons(vec![1, 2]) },
)]
fn test_molecule_apply_edit_modify_multicenter_bond_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: MulticenterBondFieldChange,
    #[case] expected: MulticenterBondForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyMulticenterBondField {
        id: MulticenterBondHandle::Id(MulticenterBondId(0)),
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
    removal_entries.multicenter[0].1 = expected;
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
