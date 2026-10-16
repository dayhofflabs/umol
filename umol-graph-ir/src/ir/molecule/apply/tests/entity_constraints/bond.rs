use super::*;

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
