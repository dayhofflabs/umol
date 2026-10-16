use super::*;

#[rstest]
#[case::order(Edit::ModifyBondOrders { changes: vec![(BondHandle::Id(BondId(0)), NumForm::Undetermined, NumForm::Lit(1)), (BondHandle::New(0), NumForm::Undetermined, NumForm::Lit(1)), (BondHandle::Id(BondId(0)), NumForm::Lit(1), NumForm::Lit(2))] }, vec![Edit::ModifyBondField { id: BondHandle::Id(BondId(0)), change: BondFieldChange::Order { old: NumForm::Undetermined, new: NumForm::Lit(1) } }, Edit::ModifyBondField { id: BondHandle::New(0), change: BondFieldChange::Order { old: NumForm::Undetermined, new: NumForm::Lit(1) } }, Edit::ModifyBondField { id: BondHandle::Id(BondId(0)), change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) } }], BondForm { order: NumForm::Lit(2), ..Default::default() }, BondForm { order: NumForm::Lit(1), ..Default::default() })]
#[case::charge(Edit::ModifyBondCharges { changes: vec![(BondHandle::Id(BondId(0)), NumForm::Undetermined, NumForm::Lit(1)), (BondHandle::New(0), NumForm::Undetermined, NumForm::Lit(1)), (BondHandle::Id(BondId(0)), NumForm::Lit(1), NumForm::Lit(2))] }, vec![Edit::ModifyBondField { id: BondHandle::Id(BondId(0)), change: BondFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) } }, Edit::ModifyBondField { id: BondHandle::New(0), change: BondFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) } }, Edit::ModifyBondField { id: BondHandle::Id(BondId(0)), change: BondFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::Lit(2) } }], BondForm { charge: NumForm::Lit(2), ..Default::default() }, BondForm { charge: NumForm::Lit(1), ..Default::default() })]
#[case::unpaired_electrons(Edit::ModifyBondUnpairedElectrons { changes: vec![(BondHandle::Id(BondId(0)), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell()), (BondHandle::New(0), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell()), (BondHandle::Id(BondId(0)), UnpairedElectronsForm::closed_shell(), UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) })] }, vec![Edit::ModifyBondField { id: BondHandle::Id(BondId(0)), change: BondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm::closed_shell() } }, Edit::ModifyBondField { id: BondHandle::New(0), change: BondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm::closed_shell() } }, Edit::ModifyBondField { id: BondHandle::Id(BondId(0)), change: BondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::closed_shell(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } } }], BondForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..Default::default() }, BondForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })]
fn test_molecule_apply_bond_batches(
    #[case] batch: Edit,
    #[case] singles: Vec<Edit>,
    #[case] changed: BondForm,
    #[case] created: BondForm,
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
    prefix.push(Edit::ModifyAtomField {
        id: AtomHandle::Id(AtomId(0)),
        change: AtomFieldChange::Charge {
            old: NumForm::Undetermined,
            new: NumForm::Lit(1),
        },
    });
    prefix.add_bond(
        AtomHandle::Id(AtomId(2)),
        AtomHandle::Id(AtomId(3)),
        BondForm::default(),
    );
    let batches = Edits::from_iter(prefix.clone().into_iter().chain([batch]));
    let singles = Edits::from_iter(prefix.into_iter().chain(singles));
    let mut expected = entries;
    expected.atoms[0].charge = NumForm::Lit(1);
    expected.bonds[0].2 = changed;
    expected.bonds.push((AtomId(2), AtomId(3), created));
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
#[case::order(Edit::ModifyBondOrders { changes: (0..3).map(|id| (BondHandle::Id(BondId(id)), NumForm::Undetermined, NumForm::Lit(1))).collect() }, BondForm { order: NumForm::Lit(1), ..Default::default() })]
#[case::charge(Edit::ModifyBondCharges { changes: (0..3).map(|id| (BondHandle::Id(BondId(id)), NumForm::Undetermined, NumForm::Lit(1))).collect() }, BondForm { charge: NumForm::Lit(1), ..Default::default() })]
#[case::unpaired_electrons(Edit::ModifyBondUnpairedElectrons { changes: (0..3).map(|id| (BondHandle::Id(BondId(id)), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell())).collect() }, BondForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })]
fn test_molecule_apply_bond_batches_error(
    #[case] mut batch: Edit,
    #[case] changed: BondForm,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    match &mut batch {
        Edit::ModifyBondOrders { changes } | Edit::ModifyBondCharges { changes } => {
            if missing {
                changes[failure].0 = BondHandle::Id(BondId(3));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        Edit::ModifyBondUnpairedElectrons { changes } => {
            if missing {
                changes[failure].0 = BondHandle::Id(BondId(3));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        _ => unreachable!(),
    }
    let entries = MoleculeEntries {
        atoms: vec![AtomForm::default(); 4],
        bonds: (0..3)
            .map(|i| (AtomId(i), AtomId(i + 1), BondForm::default()))
            .collect(),
        ..Default::default()
    };
    let source = Molecule::from_entries(entries.clone());
    let prefix = Edit::ModifyAtomField {
        id: AtomHandle::Id(AtomId(0)),
        change: AtomFieldChange::Charge {
            old: NumForm::Undetermined,
            new: NumForm::Lit(1),
        },
    };
    let error = if missing {
        TransactionError::HandleOutOfRange {
            kind: EntityKind::Bond,
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
        for (_, _, attributes) in &mut expected.bonds[..failure] {
            *attributes = changed.clone();
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
    BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
    BondForm::from_order(2),
)]
#[case::order_equivalent(
    BondFieldChange::Order { old: NumForm::lit_set([1]), new: NumForm::Lit(2) },
    BondForm::from_order(2),
)]
#[case::charge(
    BondFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) },
    BondForm { charge: NumForm::Lit(1), ..BondForm::from_order(1) },
)]
#[case::unpaired_electrons(
    BondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
    BondForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..BondForm::from_order(1) },
)]
fn test_molecule_apply_edit_modify_bond_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: BondFieldChange,
    #[case] expected: BondForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyBondField {
        id: BondHandle::Id(BondId(0)),
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
    removal_entries.bonds[0].2 = expected;
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
