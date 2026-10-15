use umol_perm::{Orientation, Permutation};

use super::*;
use crate::ir::{
    BooleanForm, FluxionalityForm, LigandPermutation, LigandSymmetryForm,
    OrientedLigandPermutation, StereoAtomForm, StereoConfigurationForm, StereoLigandPair,
    StereoLigandPosition, Topicity, TopicityForm, TopicityRelationForm,
};

#[rstest]
#[case::stereogenicity(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
    Stereogenicity::Stereogenic
)))]
#[case::topicity(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }))]
#[case::fluxionality(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_cycles(4, &[vec![0, 1]]).unwrap()), active: BooleanForm::Lit(false) }))]
#[case::symmetry(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_cycles(4, &[vec![0, 1]]).unwrap()), orientation: Orientation::Improper }, invariant: BooleanForm::Lit(true) }))]
fn test_molecule_apply_stereo_atom_constraints(
    mut removal_entries: MoleculeEntries,
    #[case] new: StereoAtomConstraintForm,
    #[values(None, Some(StereoKind::Tetrahedral), Some(StereoKind::Octahedral))] kind: Option<
        StereoKind,
    >,
    #[values(false, true)] kinded_host: bool,
) {
    if kinded_host {
        removal_entries.stereo_atoms[0].2.configuration =
            StereoConfigurationForm::kinded(StereoKind::SquarePlanar, 0u32);
    }
    let key = new.key();
    let old = StereoAtomConstraintForm::undetermined(key);
    let source = Molecule::from_entries(removal_entries.clone());
    let mut prefix = Edits::new();
    let created =
        prefix.add_stereo_atom(AtomHandle::Id(AtomId(1)), vec![], StereoAtomForm::default());
    let rows = vec![
        (
            StereoAtomHandle::Id(StereoAtomId(0)),
            None,
            Some(old.clone()),
        ),
        (created, None, Some(new.clone())),
        (
            StereoAtomHandle::Id(StereoAtomId(0)),
            Some(old),
            Some(new.clone()),
        ),
        (StereoAtomHandle::Id(StereoAtomId(0)), None, None),
        (
            StereoAtomHandle::Id(StereoAtomId(0)),
            Some(new.clone()),
            None,
        ),
    ];
    let mut batch = prefix.clone();
    batch.modify_stereo_atom_constraints(kind, key, rows.clone());
    let mut singles = prefix;
    for (id, old, new) in rows {
        singles.push(Edit::ModifyStereoAtomConstraint { id, kind, old, new });
    }
    removal_entries.stereo_atoms.push((
        AtomId(1),
        vec![],
        StereoAtomForm {
            constraints: vec![new].into(),
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
#[case::old_state(
    StereoAtomConstraintKey::Stereogenicity,
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    None
)]
#[case::old_key(StereoAtomConstraintKey::Stereogenicity, Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })), None)]
#[case::new_key(StereoAtomConstraintKey::Stereogenicity, None, Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })))]
#[case::pair(StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2))), None, Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })))]
#[case::permutation(StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::identity(4))), None, Some(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })))]
#[case::degree(StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::identity(4))), None, Some(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(6)), active: BooleanForm::Undetermined })))]
#[case::orientation(StereoAtomConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }), None, Some(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Improper }, invariant: BooleanForm::Undetermined })))]
fn test_molecule_apply_stereo_atom_constraints_error(
    mut removal_entries: MoleculeEntries,
    #[case] key: StereoAtomConstraintKey,
    #[case] old: Option<StereoAtomConstraintForm>,
    #[case] new: Option<StereoAtomConstraintForm>,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    removal_entries
        .stereo_atoms
        .push((AtomId(1), vec![], StereoAtomForm::default()));
    let source = Molecule::from_entries(removal_entries.clone());
    let form = StereoAtomConstraintForm::undetermined(key);
    let mut changes: Vec<_> = (0..3)
        .map(|i| {
            (
                StereoAtomHandle::Id(StereoAtomId(i)),
                None,
                Some(form.clone()),
            )
        })
        .collect();
    changes[failure].1 = old;
    changes[failure].2 = new;
    let error = if missing {
        changes[failure] = (StereoAtomHandle::Id(StereoAtomId(3)), None, None);
        TransactionError::HandleOutOfRange {
            kind: EntityKind::StereoAtom,
            index: 3,
            count: 3,
        }
    } else {
        TransactionError::OldStateMismatch
    };
    let batch = Edit::ModifyStereoAtomConstraints {
        kind: None,
        key,
        changes,
    };
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
            entry.2.constraints = vec![form.clone()].into();
        }
    }
    assert_eq!(molecule, Molecule::from_entries(removal_entries));
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut restored = source.clone();
        let prefix = Edit::ModifyStereoAtomConfigurations {
            changes: vec![(
                StereoAtomHandle::Id(StereoAtomId(0)),
                StereoConfigurationForm::Undetermined,
                StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32),
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
#[case::empty(vec![], vec![])]
#[case::noop(vec![(StereoAtomHandle::Id(StereoAtomId(0)), None, None)], vec![])]
#[case::normalized(
    vec![(StereoAtomHandle::Id(StereoAtomId(0)), Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic]))), None)],
    vec![(StereoAtomId(0), Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic]))))]
)]
fn test_molecule_apply_undo_stereo_atom_constraints(
    mut removal_entries: MoleculeEntries,
    #[case] changes: Vec<(
        StereoAtomHandle,
        Option<StereoAtomConstraintForm>,
        Option<StereoAtomConstraintForm>,
    )>,
    #[case] constraints: Vec<(StereoAtomId, Option<StereoAtomConstraintForm>)>,
) {
    let form = StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic,
    ));
    removal_entries.stereo_atoms[0].2.constraints = vec![form].into();
    let mut molecule = Molecule::from_entries(removal_entries.clone());
    let mut state = ApplicationState::new(&molecule);
    let key = StereoAtomConstraintKey::Stereogenicity;
    let undo = molecule
        .apply_edit_with_undo(
            Edit::ModifyStereoAtomConstraints {
                kind: None,
                key,
                changes,
            },
            &mut state,
        )
        .unwrap()
        .unwrap();
    assert_eq!(
        undo,
        Undo::RestoreStereoAtomConstraints {
            key,
            constraints: constraints.clone()
        }
    );
    molecule.apply_undo(undo);
    if let Some((_, constraint)) = constraints.into_iter().next() {
        removal_entries.stereo_atoms[0].2.constraints = constraint.into_iter().collect();
    }
    assert_eq!(molecule, Molecule::from_entries(removal_entries));
}

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
