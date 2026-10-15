use super::*;

#[rstest]
#[case::addition(
    None,
    None,
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic
    )))
)]
#[case::modification(
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic
    )))
)]
#[case::removal(
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    None
)]
#[case::absent(None, None, None)]
fn test_molecule_apply_edit_modify_stereo_atom_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<StereoAtomConstraintForm>,
    #[case] old: Option<StereoAtomConstraintForm>,
    #[case] new: Option<StereoAtomConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.stereo_atoms[0].2.constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyStereoAtomConstraint {
        id: StereoAtomHandle::Id(StereoAtomId(0)),
        kind: Some(StereoKind::CisTrans),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreStereoAtomConstraint {
            id: StereoAtomId(0),
            key: StereoAtomConstraintKey::Stereogenicity,
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
    removal_entries.stereo_atoms[0].2.constraints = new.into_iter().collect::<Vec<_>>().into();
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
#[case::addition(
    None,
    None,
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic
    )))
)]
#[case::modification(
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic
    )))
)]
#[case::removal(
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    None
)]
#[case::absent(None, None, None)]
fn test_molecule_apply_edit_modify_stereo_bond_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] current: Option<StereoBondConstraintForm>,
    #[case] old: Option<StereoBondConstraintForm>,
    #[case] new: Option<StereoBondConstraintForm>,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.stereo_bonds[0].2.constraints = current.into_iter().collect::<Vec<_>>().into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyStereoBondConstraint {
        id: StereoBondHandle::Id(StereoBondId(0)),
        kind: Some(StereoKind::Tetrahedral),
        old: old.clone(),
        new: new.clone(),
    };
    let expected_undo = if old.is_none() && new.is_none() {
        None
    } else {
        Some(Undo::RestoreStereoBondConstraint {
            id: StereoBondId(0),
            key: StereoBondConstraintKey::Stereogenicity,
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
    removal_entries.stereo_bonds[0].2.constraints = new.into_iter().collect::<Vec<_>>().into();
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
