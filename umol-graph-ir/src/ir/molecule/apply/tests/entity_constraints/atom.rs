use super::*;

#[rstest]
#[case(AtomConstraintKey::Valence)]
#[case(AtomConstraintKey::DonatedPairs)]
#[case(AtomConstraintKey::AcceptedPairs)]
#[case(AtomConstraintKey::AromaticValence)]
#[case(AtomConstraintKey::MulticenterValence)]
#[case(AtomConstraintKey::TetrahedralStereo)]
#[case(AtomConstraintKey::Degree)]
#[case(AtomConstraintKey::TotalDegree)]
#[case(AtomConstraintKey::TotalValence)]
#[case(AtomConstraintKey::RingDegree)]
#[case(AtomConstraintKey::RingValence)]
#[case(AtomConstraintKey::TotalHydrogens)]
#[case(AtomConstraintKey::RingMembership(RingScope::All))]
#[case(AtomConstraintKey::RingMembership(RingScope::Size(5)))]
fn test_molecule_apply_atom_constraint_batches(#[case] key: AtomConstraintKey) {
    let form = AtomConstraintForm::undetermined(key);
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default()],
        ..Default::default()
    });
    let rows = vec![
        (AtomHandle::Id(AtomId(0)), None, Some(form.clone())),
        (AtomHandle::New(0), None, Some(form.clone())),
        (
            AtomHandle::Id(AtomId(0)),
            Some(form.clone()),
            Some(form.clone()),
        ),
        (AtomHandle::Id(AtomId(0)), None, None),
        (AtomHandle::Id(AtomId(0)), Some(form.clone()), None),
    ];
    let mut edits = Edits::new();
    edits.add_atom(AtomForm::default());
    edits.modify_atom_constraints(key, rows.clone());
    let mut singles = Edits::new();
    singles.add_atom(AtomForm::default());
    for (id, old, new) in rows {
        singles.push(Edit::ModifyAtomConstraint { id, old, new });
    }
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::default(),
            AtomForm {
                constraints: vec![form].into(),
                ..Default::default()
            },
        ],
        ..Default::default()
    });
    let actual = source.clone().tracked_apply(edits.clone()).unwrap();
    assert_eq!(actual.0, expected);
    assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
    let mut molecule = source.clone();
    Transaction::run(
        &mut molecule,
        RollbackMethod::Journal,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits)?;
            transaction.rollback();
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(molecule, source);
}

#[rstest]
#[case(Edit::ModifyAtomConstraints { key: AtomConstraintKey::RingMembership(RingScope::Size(5)), changes: vec![] })]
fn test_molecule_apply_atom_constraint_batches_identity(#[case] edit: Edit) {
    let source = Molecule::new();
    assert_eq!(
        source
            .clone()
            .apply(Edits::from_iter([edit.clone()]))
            .unwrap(),
        source
    );
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(edit, &mut state)
        .unwrap()
        .unwrap();
    molecule.apply_undo(undo);
    assert_eq!(molecule, source);
}

#[rstest]
fn test_molecule_apply_atom_constraint_batches_normalized() {
    let accepted = AtomConstraintForm::valence(NumForm::LitSet(Box::new([1].into())));
    let mut molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm {
            constraints: vec![AtomConstraintForm::valence(1)].into(),
            ..Default::default()
        }],
        ..Default::default()
    });
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(
            Edit::ModifyAtomConstraints {
                key: AtomConstraintKey::Valence,
                changes: vec![(AtomHandle::Id(AtomId(0)), Some(accepted.clone()), None)],
            },
            &mut state,
        )
        .unwrap()
        .unwrap();
    assert_eq!(
        undo,
        Undo::RestoreAtomConstraints {
            key: AtomConstraintKey::Valence,
            constraints: vec![(AtomId(0), Some(accepted.clone()))]
        }
    );
    molecule.apply_undo(undo);
    assert_eq!(
        molecule,
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm {
                constraints: vec![accepted].into(),
                ..Default::default()
            }],
            ..Default::default()
        })
    );
}

#[rstest]
#[case::old_key(Some(AtomConstraintForm::degree(1)), None)]
#[case::new_key(None, Some(AtomConstraintForm::degree(1)))]
#[case::old_state(Some(AtomConstraintForm::valence(1)), None)]
#[case::ring_scope(None, Some(AtomConstraintForm::ring_membership(RingScope::Size(6), 1)))]
fn test_molecule_apply_atom_constraint_batches_error(
    #[case] old: Option<AtomConstraintForm>,
    #[case] new: Option<AtomConstraintForm>,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] journaled: bool,
) {
    let key = if matches!(&new, Some(AtomConstraintForm::RingMembership(_))) {
        AtomConstraintKey::RingMembership(RingScope::Size(5))
    } else {
        AtomConstraintKey::Valence
    };
    let form = AtomConstraintForm::undetermined(key);
    let mut rows: Vec<_> = (0..3)
        .map(|id| (AtomHandle::Id(AtomId(id)), None, Some(form.clone())))
        .collect();
    rows[failure].1 = old;
    rows[failure].2 = new;
    let batch = Edit::ModifyAtomConstraints { key, changes: rows };
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 3],
        ..Default::default()
    });
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
    let mut atoms = vec![AtomForm::default(); 3];
    if !journaled {
        for atom in &mut atoms[..failure] {
            atom.constraints = vec![form.clone()].into();
        }
    }
    assert_eq!(
        molecule,
        Molecule::from_entries(MoleculeEntries {
            atoms,
            ..Default::default()
        })
    );
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
#[case::addition(None, None, Some(AtomConstraintForm::degree(3)))]
#[case::modification(
    Some(AtomConstraintForm::degree(2)),
    Some(AtomConstraintForm::degree(2)),
    Some(AtomConstraintForm::degree(3))
)]
#[case::removal(
    Some(AtomConstraintForm::degree(2)),
    Some(AtomConstraintForm::degree(2)),
    None
)]
#[case::absent(None, None, None)]
#[case::equivalent_old(Some(AtomConstraintForm::degree(2)), Some(AtomConstraintForm::degree(NumForm::lit_set([2]))), Some(AtomConstraintForm::degree(3)))]
fn test_molecule_apply_edit_modify_atom_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<AtomConstraintForm>,
    #[case] old: Option<AtomConstraintForm>,
    #[case] new: Option<AtomConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.atoms[0].constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyAtomConstraint {
        id: AtomHandle::Id(AtomId(0)),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreAtomConstraint {
            id: AtomId(0),
            key: AtomConstraintKey::Degree,
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
    removal_entries.atoms[0].constraints = new.into_iter().collect::<Vec<_>>().into();
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

#[rstest]
#[case::key(AtomConstraintForm::degree(2), AtomConstraintForm::valence(2))]
#[case::value(AtomConstraintForm::degree(3), AtomConstraintForm::degree(4))]
fn test_molecule_apply_edit_modify_atom_constraint_error(
    mut removal_entries: MoleculeEntries,
    #[case] old: AtomConstraintForm,
    #[case] new: AtomConstraintForm,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.atoms[0].constraints = vec![AtomConstraintForm::degree(2)].into();
    let initial = Molecule::from_entries(removal_entries);
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyAtomConstraint {
        id: AtomHandle::Id(AtomId(0)),
        old: Some(old),
        new: Some(new),
    };
    let result = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(result, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, initial);
}
