use pretty_assertions::assert_eq;

use super::*;
use crate::ir::{MoleculeApplyError, RollbackMethod, Transaction};

#[rustfmt::skip]
#[rstest]
#[case::configuration_and_constraint(
    StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    StereoAtomUpdate {
        configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) },
        constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    },
    vec![
        Edit::ModifyStereoAtomField {
            id: StereoAtomHandle::Id(StereoAtomId(7)),
            change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32) },
        },
        Edit::ModifyStereoAtomConstraint {
            id: StereoAtomHandle::Id(StereoAtomId(7)),
            kind: Some(StereoKind::Tetrahedral),
            old: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
            new: None,
        },
    ],
)]
#[case::kind_assignment_raw_precondition(
    StereoAtomForm { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic]))), ..Default::default() },
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric))) },
    vec![
        Edit::ModifyStereoAtomField { id: StereoAtomHandle::Id(StereoAtomId(7)), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::Undetermined, new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined) } },
        Edit::ModifyStereoAtomConstraint { id: StereoAtomHandle::Id(StereoAtomId(7)), kind: Some(StereoKind::Tetrahedral), old: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic]))), new: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric))) },
    ],
)]
fn test_edits_update_stereo_atom(
    #[case] current: StereoAtomForm,
    #[case] update: StereoAtomUpdate,
    #[case] expected: Vec<Edit>,
) {
    let mut edits = Edits::new();
    edits.update_stereo_atom(
        StereoAtomHandle::Id(StereoAtomId(7)),
        &current,
        &update,
    );

    assert_eq!(edits.as_slice(), expected);

    let expected_attributes = current.update(&update);
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 5],
        bonds: (1..=4)
            .map(|index| (AtomId(0), AtomId(index), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|index| StereoLigand::new(AtomId(index), StereoLigandKind::Atom))
                .collect(),
            current.clone(),
        )],
        ..Default::default()
    });
    let mut applied_edits = Edits::new();
    applied_edits.update_stereo_atom(
        StereoAtomHandle::Id(StereoAtomId(0)),
        &current,
        &update,
    );
    let mut journaled = molecule.clone();
    Transaction::run(
        &mut journaled,
        RollbackMethod::Journal,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(applied_edits.clone())?;
            assert_eq!(
                transaction.probe()?.stereo_atom(StereoAtomId(0)).attributes(),
                &expected_attributes,
            );
            transaction.rollback();
            Ok(())
        },
    ).unwrap();
    assert_eq!(journaled, molecule);
    let mut editor = molecule.edit();
    editor
        .apply(applied_edits)
        .expect("stereo-atom update edits should apply");

    assert_eq!(
        editor.stereo_atom(StereoAtomId(0)).attributes(),
        &expected_attributes
    );
}

#[rustfmt::skip]
#[rstest]
#[case::empty(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate::default())]
#[case::relative(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() })]
#[case::absent_constraint_removal(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() })]
#[case::kindless_constraint(
    StereoAtomForm { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic]))), ..Default::default() },
    StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([]))), ..Default::default() },
)]
fn test_edits_update_stereo_atom_identity(
    #[case] current: StereoAtomForm,
    #[case] update: StereoAtomUpdate,
) {
    let mut edits = Edits::new();
    edits.update_stereo_atom(
        StereoAtomHandle::Id(StereoAtomId(0)),
        &current,
        &update,
    );

    assert_eq!(edits, Edits::new());
}

#[rustfmt::skip]
#[rstest]
#[case::configuration_only(
    StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() },
    vec![Edit::ModifyStereoAtomField { id: StereoAtomHandle::Id(StereoAtomId(7)), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32) } }],
)]
#[case::constraint_only(
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), ..Default::default() },
    vec![Edit::ModifyStereoAtomConstraint { id: StereoAtomHandle::Id(StereoAtomId(7)), kind: Some(StereoKind::Tetrahedral), old: None, new: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) }],
)]
#[case::clearing(
    StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
    vec![
        Edit::ModifyStereoAtomField { id: StereoAtomHandle::Id(StereoAtomId(7)), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), new: StereoConfigurationForm::Undetermined } },
    ],
)]
#[case::matching_kind_with_constraint(
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    vec![Edit::ModifyStereoAtomConstraint { id: StereoAtomHandle::Id(StereoAtomId(7)), kind: Some(StereoKind::Tetrahedral), old: None, new: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) }],
)]
#[case::different_kind_with_constraint(
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::SquarePlanar, coset: None }, constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    vec![
        Edit::ModifyStereoAtomField { id: StereoAtomHandle::Id(StereoAtomId(7)), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), new: StereoConfigurationForm::kinded(StereoKind::SquarePlanar, StereoCoset::Undetermined) } },
        Edit::ModifyStereoAtomConstraint { id: StereoAtomHandle::Id(StereoAtomId(7)), kind: Some(StereoKind::SquarePlanar), old: None, new: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    ],
)]
fn test_edits_update_stereo_atom_construction(
    #[case] current: StereoAtomForm,
    #[case] update: StereoAtomUpdate,
    #[case] expected: Vec<Edit>,
) {
    let mut edits = Edits::new();
    edits.update_stereo_atom(StereoAtomHandle::Id(StereoAtomId(7)), &current, &update);
    assert_eq!(edits.as_slice(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::configuration_and_constraint(
    StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    StereoBondUpdate {
        configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) },
        constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
    },
    vec![
        Edit::ModifyStereoBondField {
            id: StereoBondHandle::Id(StereoBondId(7)),
            change: StereoBondFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32) },
        },
        Edit::ModifyStereoBondConstraint {
            id: StereoBondHandle::Id(StereoBondId(7)),
            kind: Some(StereoKind::CisTrans),
            old: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
            new: None,
        },
    ],
)]
#[case::kind_assignment_raw_precondition(
    StereoBondForm { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic]))), ..Default::default() },
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric))) },
    vec![
        Edit::ModifyStereoBondField { id: StereoBondHandle::Id(StereoBondId(7)), change: StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Undetermined, new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined) } },
        Edit::ModifyStereoBondConstraint { id: StereoBondHandle::Id(StereoBondId(7)), kind: Some(StereoKind::CisTrans), old: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic]))), new: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric))) },
    ],
)]
fn test_edits_update_stereo_bond(
    #[case] current: StereoBondForm,
    #[case] update: StereoBondUpdate,
    #[case] expected: Vec<Edit>,
) {
    let mut edits = Edits::new();
    edits.update_stereo_bond(
        StereoBondHandle::Id(StereoBondId(7)),
        &current,
        &update,
    );

    assert_eq!(edits.as_slice(), expected);

    let expected_attributes = current.update(&update);
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(0),
            (2..=5)
                .map(|index| StereoLigand::new(AtomId(index), StereoLigandKind::Atom))
                .collect(),
            current.clone(),
        )],
        ..Default::default()
    });
    let mut applied_edits = Edits::new();
    applied_edits.update_stereo_bond(
        StereoBondHandle::Id(StereoBondId(0)),
        &current,
        &update,
    );
    let mut journaled = molecule.clone();
    Transaction::run(
        &mut journaled,
        RollbackMethod::Journal,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(applied_edits.clone())?;
            assert_eq!(
                transaction.probe()?.stereo_bond(StereoBondId(0)).attributes(),
                &expected_attributes,
            );
            transaction.rollback();
            Ok(())
        },
    ).unwrap();
    assert_eq!(journaled, molecule);
    let mut editor = molecule.edit();
    editor
        .apply(applied_edits)
        .expect("stereo-bond update edits should apply");

    assert_eq!(
        editor.stereo_bond(StereoBondId(0)).attributes(),
        &expected_attributes
    );
}

#[rustfmt::skip]
#[rstest]
#[case::empty(StereoBondForm::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate::default())]
#[case::relative(StereoBondForm::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() })]
#[case::absent_constraint_removal(StereoBondForm::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() })]
#[case::kindless_constraint(
    StereoBondForm { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic]))), ..Default::default() },
    StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([]))), ..Default::default() },
)]
fn test_edits_update_stereo_bond_identity(
    #[case] current: StereoBondForm,
    #[case] update: StereoBondUpdate,
) {
    let mut edits = Edits::new();
    edits.update_stereo_bond(
        StereoBondHandle::Id(StereoBondId(0)),
        &current,
        &update,
    );

    assert_eq!(edits, Edits::new());
}

#[rustfmt::skip]
#[rstest]
#[case::configuration_only(
    StereoBondForm::new(StereoKind::CisTrans, 0_u32),
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() },
    vec![Edit::ModifyStereoBondField { id: StereoBondHandle::Id(StereoBondId(7)), change: StereoBondFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32) } }],
)]
#[case::constraint_only(
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), ..Default::default() },
    vec![Edit::ModifyStereoBondConstraint { id: StereoBondHandle::Id(StereoBondId(7)), kind: Some(StereoKind::CisTrans), old: None, new: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) }],
)]
#[case::clearing(
    StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
    vec![
        Edit::ModifyStereoBondField { id: StereoBondHandle::Id(StereoBondId(7)), change: StereoBondFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), new: StereoConfigurationForm::Undetermined } },
    ],
)]
#[case::matching_kind_with_constraint(
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    vec![Edit::ModifyStereoBondConstraint { id: StereoBondHandle::Id(StereoBondId(7)), kind: Some(StereoKind::CisTrans), old: None, new: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) }],
)]
#[case::missing_current_kind_with_constraint(
    StereoBondForm::default(),
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    vec![
        Edit::ModifyStereoBondField { id: StereoBondHandle::Id(StereoBondId(7)), change: StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Undetermined, new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined) } },
        Edit::ModifyStereoBondConstraint { id: StereoBondHandle::Id(StereoBondId(7)), kind: Some(StereoKind::CisTrans), old: None, new: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    ],
)]
fn test_edits_update_stereo_bond_construction(
    #[case] current: StereoBondForm,
    #[case] update: StereoBondUpdate,
    #[case] expected: Vec<Edit>,
) {
    let mut edits = Edits::new();
    edits.update_stereo_bond(StereoBondHandle::Id(StereoBondId(7)), &current, &update);
    assert_eq!(edits.as_slice(), expected);
}
