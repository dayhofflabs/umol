use pretty_assertions::assert_eq;

use super::*;

#[rstest]
fn test_stereo_bond_form_new() {
    let stereo_bond = StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined);
    assert_eq!(
        stereo_bond.configuration,
        StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined)
    );
    assert_eq!(stereo_bond.constraints, StereoBondConstraintsForm::new())
}

#[rustfmt::skip]
#[rstest]
#[case::stereogenicity(
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
    StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
)]
fn test_stereo_bond_form_with_methods(
    #[case] bond: StereoBondForm,
    #[case] constraint: StereoBondConstraintForm,
    #[case] expected: StereoBondForm,
) {
    assert_eq!(bond.clone().with_constraint(constraint.clone()), expected);
    assert_eq!(bond.with_constraints([constraint]), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::folds_coset(
    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set([1])),
    Ok(StereoBondForm::new(StereoKind::CisTrans, 1u32)),
)]
#[case::empty_coset_litset_contradiction(
    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set(Vec::<u32>::new())),
    Err(Contradiction),
)]
fn test_stereo_bond_form_normalize(
    #[case] input: StereoBondForm,
    #[case] expected: Result<StereoBondForm, Contradiction>,
) {
    assert_eq!(input.normalize(), expected);
}

#[rstest]
#[case::symmetry(
    StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper },
        invariant: BooleanForm::Lit(true),
    }),
    Some(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper },
        invariant: BooleanForm::Lit(true),
    })),
)]
#[case::fluxionality(
    StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) }),
    Some(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) })),
)]
#[case::topicity(
    StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::lit_set([Topicity::Homotopic]) }),
    Some(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
)]
#[case::stereogenicity(
    StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic])),
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
)]
#[case::topicity_bottom(
    StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::lit_set([]) }),
    None,
)]
#[case::stereogenicity_bottom(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])), None)]
fn test_stereo_bond_form_normalize_constraints(
    #[case] raw: StereoBondConstraintForm,
    #[case] normalized: Option<StereoBondConstraintForm>,
    #[values(
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Undetermined)
    )]
    configuration: StereoConfigurationForm,
) {
    let form = StereoBondForm {
        configuration: configuration.clone(),
        constraints: raw.clone().into(),
    };
    let expected = if configuration.kind().is_none() {
        Ok(StereoBondForm::default())
    } else {
        normalized
            .map(|value| StereoBondForm {
                configuration,
                constraints: value.into(),
            })
            .ok_or(Contradiction)
    };
    assert_eq!(form.clone().normalize(), expected);
    assert_eq!(form.constraints.iter().collect::<Vec<_>>(), vec![&raw]);
}

#[rstest]
#[case::literal(StereogenicityForm::Lit(Stereogenicity::Symmetric))]
#[case::contradiction(StereogenicityForm::lit_set([]))]
fn test_stereo_bond_form_normalized_eq_constraints(#[case] payload: StereogenicityForm) {
    let form = StereoBondForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: StereoBondConstraintForm::Stereogenicity(payload.clone()).into(),
    };
    let top = StereoBondForm::default();
    assert_ne!(form, top);
    assert!(form.normalized_eq(&top));
    assert!(top.normalized_eq(&form));
    assert!(form.is_undetermined());
    assert!(!form.is_ground());
    assert_eq!(form.constraints.stereogenicity(), &payload);
}

#[rustfmt::skip]
#[rstest]
#[case::same_kind_narrows(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined), StereoBondForm::new(StereoKind::CisTrans, 1u32),
    Some(StereoBondForm::new(StereoKind::CisTrans, 1u32)))]
#[case::config_conflict(StereoBondForm::new(StereoKind::CisTrans, 0u32), StereoBondForm::new(StereoKind::CisTrans, 1u32), None)]
fn test_stereo_bond_form_meet(
    #[case] a: StereoBondForm,
    #[case] b: StereoBondForm,
    #[case] expected: Option<StereoBondForm>,
) {
    assert_eq!(a.meet(&b), expected);
}

#[rstest]
#[case::literal(StereogenicityForm::Lit(Stereogenicity::Symmetric))]
#[case::contradiction(StereogenicityForm::lit_set([]))]
fn test_stereo_bond_form_meet_constraints(
    #[case] payload: StereogenicityForm,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(0))] coset: StereoCoset,
) {
    let kindless = StereoBondForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: StereoBondConstraintForm::Stereogenicity(payload).into(),
    };
    let kinded = StereoBondForm {
        configuration: StereoConfigurationForm::Kinded(StereoKind::CisTrans, coset),
        constraints: StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        ))
        .into(),
    };
    assert_eq!(kindless.meet(&kindless), Some(StereoBondForm::default()));
    assert_eq!(kindless.meet(&kinded), Some(kinded.clone()));
    assert_eq!(kinded.meet(&kindless), Some(kinded.clone()));
    assert!(kindless.matches(&kinded));
    assert!(!kinded.matches(&kindless));
    assert!(kindless.is_compatible(&kinded));
    assert!(kinded.is_compatible(&kindless));
}

#[rstest]
#[case::empty_payload(StereoConfigurationForm::Kinded(
    StereoKind::CisTrans,
    StereoCoset::Undetermined
))]
fn test_stereo_bond_form_meet_error(#[case] configuration: StereoConfigurationForm) {
    let kindless = StereoBondForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([]))
            .into(),
    };
    let bottom = StereoBondForm {
        configuration,
        constraints: StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([]))
            .into(),
    };
    assert_eq!(kindless.meet(&bottom), None);
    assert_eq!(bottom.meet(&kindless), None);
    assert!(!kindless.matches(&bottom));
    assert!(!kindless.is_compatible(&bottom));
    assert!(!bottom.is_compatible(&kindless));
}

#[rstest]
#[case::kindless(StereoConfigurationForm::Undetermined)]
#[case::other_kind(StereoConfigurationForm::Kinded(
    StereoKind::Tetrahedral,
    StereoCoset::Undetermined
))]
fn test_stereo_bond_form_join_constraints(#[case] configuration: StereoConfigurationForm) {
    let a = StereoBondForm {
        configuration,
        constraints: StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        ))
        .into(),
    };
    let b = StereoBondForm {
        configuration: StereoConfigurationForm::Kinded(
            StereoKind::CisTrans,
            StereoCoset::Undetermined,
        ),
        constraints: a.constraints.clone(),
    };
    assert_eq!(a.join(&b), Ok(StereoBondForm::default()));
    assert_eq!(b.join(&a), Ok(StereoBondForm::default()));
}

#[rstest]
#[case::cis_trans(
    StereoBondForm::new(StereoKind::CisTrans, 0u32),
    StereoBondForm::new(StereoKind::CisTrans, 1u32)
)]
fn test_stereo_bond_form_swap(#[case] input: StereoBondForm, #[case] expected: StereoBondForm) {
    assert_eq!(input.swap(), Some(expected));
}

/// A within-endpoint swap is an action of the cis/trans parent group and exchanges the two
/// cosets.
#[rstest]
fn test_stereo_bond_form_reframe_by() {
    let before = [
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
    ];
    let after = [
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
    ];
    let permutation = Permutation::between(&before, &after).expect("a unique frame change");
    assert_eq!(
        StereoBondForm::new(StereoKind::CisTrans, 0u32).reframe_by(&permutation),
        Some(StereoBondForm::new(StereoKind::CisTrans, 1u32)),
    );
}

#[rstest]
#[case::degree(
    StereoBondForm::new(StereoKind::CisTrans, 0u32),
    Permutation::identity(3)
)]
#[case::kinded_outside_parent(
    StereoBondForm::new(StereoKind::CisTrans, 0u32),
    Permutation::from_image(&[1, 2, 0, 3]),
)]
#[case::kindless_outside_parent(
    StereoBondForm::default(),
    Permutation::from_image(&[1, 2, 0, 3]),
)]
fn test_stereo_bond_form_frame_transport_domain(
    #[case] form: StereoBondForm,
    #[case] permutation: Permutation,
) {
    assert_eq!(form.reframe_by(&permutation), None,);
}

#[rustfmt::skip]
#[rstest]
#[case::absolute(
    StereoBondForm::new(StereoKind::CisTrans, 0_u32),
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() },
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
)]
#[case::relative(
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() },
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
)]
#[case::undetermined(
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() },
    StereoBondForm::default(),
)]
#[case::explicit_open(
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Undetermined) }, ..Default::default() },
    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
)]
#[case::constraint_set(
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), ..Default::default() },
    StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
)]
#[case::constraint_remove(
    StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() },
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
)]
fn test_stereo_bond_form_update(
    #[case] bond: StereoBondForm,
    #[case] update: StereoBondUpdate,
    #[case] expected: StereoBondForm,
) {
    assert_eq!(bond.update(&update), expected);
}

#[rstest]
#[case::empty(StereoBondForm::new(StereoKind::CisTrans, 1_u32))]
fn test_stereo_bond_form_update_identity(#[case] bond: StereoBondForm) {
    assert_eq!(bond.update(&StereoBondUpdate::default()), bond);
}

#[rustfmt::skip]
#[rstest]
#[case::configuration_and_constraint(
    StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    StereoBondForm::default(),
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
)]
#[case::constraint_context(
    StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
)]
fn test_stereo_bond_form_difference_to(
    #[case] bond: StereoBondForm,
    #[case] other: StereoBondForm,
    #[case] expected: StereoBondUpdate,
) {
    assert_eq!(bond.difference_to(&other), expected);
}

#[rstest]
#[case::same(StereoBondForm::new(StereoKind::CisTrans, 1_u32))]
fn test_stereo_bond_form_difference_to_identity(#[case] bond: StereoBondForm) {
    assert_eq!(bond.difference_to(&bond), StereoBondUpdate::default());
}
