use super::*;

#[rstest]
#[case::literal(
    StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
    StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32)
)]
#[case::kind_change(
    StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
    StereoConfigurationForm::kinded(StereoKind::SquarePlanar, 0u32)
)]
#[case::raw_coset(
    StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 2u32),
    StereoConfigurationForm::Undetermined
)]
fn test_molecule_apply_stereo_atom_configurations(
    mut removal_entries: MoleculeEntries,
    #[case] first: StereoConfigurationForm,
    #[case] second: StereoConfigurationForm,
) {
    let source = Molecule::from_entries(removal_entries.clone());
    let mut prefix = Edits::new();
    let created =
        prefix.add_stereo_atom(AtomHandle::Id(AtomId(1)), vec![], StereoAtomForm::default());
    let rows = vec![
        (
            StereoAtomHandle::Id(StereoAtomId(0)),
            StereoConfigurationForm::Undetermined,
            first.clone(),
        ),
        (
            created,
            StereoConfigurationForm::Undetermined,
            first.clone(),
        ),
        (
            StereoAtomHandle::Id(StereoAtomId(0)),
            first.clone(),
            second.clone(),
        ),
    ];
    let mut batch = prefix.clone();
    batch.modify_stereo_atom_configurations(rows.clone());
    let mut singles = prefix;
    for (id, old, new) in rows {
        singles.push(Edit::ModifyStereoAtomField {
            id,
            change: StereoAtomFieldChange::Configuration { old, new },
        });
    }
    removal_entries.stereo_atoms[0].2.configuration = second;
    removal_entries.stereo_atoms.push((
        AtomId(1),
        vec![],
        StereoAtomForm {
            configuration: first,
            ..Default::default()
        },
    ));
    let expected = Molecule::from_entries(removal_entries);
    let actual = source.clone().tracked_apply(batch.clone()).unwrap();
    assert_eq!(actual.0, expected);
    assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut restored = source.clone();
        Transaction::run(
            &mut restored,
            recovery,
            |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(restored, source);
    }
}

#[rstest]
fn test_molecule_apply_stereo_atom_configurations_error(
    mut removal_entries: MoleculeEntries,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    removal_entries
        .stereo_atoms
        .push((AtomId(1), vec![], StereoAtomForm::default()));
    let source = Molecule::from_entries(removal_entries.clone());
    let new = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32);
    let mut changes: Vec<_> = (0..3)
        .map(|i| {
            (
                StereoAtomHandle::Id(StereoAtomId(i)),
                StereoConfigurationForm::Undetermined,
                new.clone(),
            )
        })
        .collect();
    let error = if missing {
        changes[failure].0 = StereoAtomHandle::Id(StereoAtomId(3));
        TransactionError::HandleOutOfRange {
            kind: EntityKind::StereoAtom,
            index: 3,
            count: 3,
        }
    } else {
        changes[failure].1 = new.clone();
        TransactionError::OldStateMismatch
    };
    let batch = Edit::ModifyStereoAtomConfigurations { changes };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    if !journaled {
        for entry in &mut removal_entries.stereo_atoms[..failure] {
            entry.2.configuration = new.clone();
        }
    }
    assert_eq!(molecule, Molecule::from_entries(removal_entries));
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut restored = source.clone();
        let prefix = Edit::ModifyAtomCharges {
            changes: vec![(
                AtomHandle::Id(AtomId(0)),
                NumForm::Undetermined,
                NumForm::Lit(1),
            )],
        };
        assert_eq!(
            restored.transact([Edits::from_iter([prefix, batch.clone()])], recovery),
            Err(error.clone().into())
        );
        assert_eq!(restored, source);
    }
}

#[rstest]
fn test_molecule_apply_stereo_atom_configurations_normalized(mut removal_entries: MoleculeEntries) {
    let kind = StereoKind::Tetrahedral;
    removal_entries.stereo_atoms[0].2.configuration = StereoConfigurationForm::kinded(kind, 0u32);
    let old = StereoConfigurationForm::kinded(kind, StereoCoset::lit_set([0]));
    let new = StereoConfigurationForm::kinded(kind, 1u32);
    let mut molecule = Molecule::from_entries(removal_entries.clone());
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(
            Edit::ModifyStereoAtomConfigurations {
                changes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    old.clone(),
                    new.clone(),
                )],
            },
            &mut state,
        )
        .unwrap()
        .unwrap();
    assert_eq!(
        undo,
        Undo::ModifyStereoAtomConfigurations {
            changes: vec![(StereoAtomId(0), new.clone(), old.clone())]
        }
    );
    removal_entries.stereo_atoms[0].2.configuration = new;
    assert_eq!(molecule, Molecule::from_entries(removal_entries.clone()));
    molecule.apply_undo(undo);
    removal_entries.stereo_atoms[0].2.configuration = old;
    assert_eq!(molecule, Molecule::from_entries(removal_entries));
}

#[rstest]
#[case::configuration(
    StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::Undetermined, new: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)) },
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
)]
#[case::coset_range(
    StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::Undetermined, new: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(2)) },
    StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32),
)]
fn test_molecule_apply_edit_modify_stereo_atom_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: StereoAtomFieldChange,
    #[case] expected: StereoAtomForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyStereoAtomField {
        id: StereoAtomHandle::Id(StereoAtomId(0)),
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
    removal_entries.stereo_atoms[0].2 = expected;
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

#[rstest]
#[case::configuration(
    StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)), new: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0)) },
    StereoBondForm::new(StereoKind::CisTrans, 0_u32),
)]
#[case::configuration_equivalent(
    StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::lit_set([1])), new: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0)) },
    StereoBondForm::new(StereoKind::CisTrans, 0_u32),
)]
#[case::coset_range(
    StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)), new: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(2)) },
    StereoBondForm::new(StereoKind::CisTrans, 2_u32),
)]
fn test_molecule_apply_edit_modify_stereo_bond_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: StereoBondFieldChange,
    #[case] expected: StereoBondForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyStereoBondField {
        id: StereoBondHandle::Id(StereoBondId(0)),
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
    removal_entries.stereo_bonds[0].2 = expected;
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
