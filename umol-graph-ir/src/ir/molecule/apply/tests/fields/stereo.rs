use super::*;

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
