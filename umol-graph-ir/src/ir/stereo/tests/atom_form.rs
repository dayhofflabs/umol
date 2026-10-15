use pretty_assertions::assert_eq;

use super::*;

#[rstest]
fn test_stereo_atom_form_new() {
    let stereo_atom = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined);
    assert_eq!(
        stereo_atom.configuration,
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)
    );
    assert_eq!(stereo_atom.constraints, StereoAtomConstraintsForm::new());
}

#[rustfmt::skip]
#[rstest]
#[case::stereogenicity(
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
    StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
)]
fn test_stereo_atom_form_with_methods(
    #[case] atom: StereoAtomForm,
    #[case] constraint: StereoAtomConstraintForm,
    #[case] expected: StereoAtomForm,
) {
    assert_eq!(atom.clone().with_constraint(constraint.clone()), expected);
    assert_eq!(atom.with_constraints([constraint]), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::folds_coset(
    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([1])),
    Ok(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)),
)]
#[case::empty_coset_litset_contradiction(
    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set(Vec::<u32>::new())),
    Err(Contradiction),
)]
fn test_stereo_atom_form_normalize(
    #[case] input: StereoAtomForm,
    #[case] expected: Result<StereoAtomForm, Contradiction>,
) {
    assert_eq!(input.normalize(), expected);
}

#[rstest]
#[case::symmetry(
    StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper },
        invariant: BooleanForm::Lit(true),
    }),
    Some(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper },
        invariant: BooleanForm::Lit(true),
    })),
)]
#[case::fluxionality(
    StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) }),
    Some(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) })),
)]
#[case::topicity(
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::lit_set([Topicity::Homotopic]) }),
    Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
)]
#[case::stereogenicity(
    StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic])),
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
)]
#[case::topicity_bottom(
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::lit_set([]) }),
    None,
)]
#[case::stereogenicity_bottom(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])), None)]
fn test_stereo_atom_form_normalize_constraints(
    #[case] raw: StereoAtomConstraintForm,
    #[case] normalized: Option<StereoAtomConstraintForm>,
    #[values(
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)
    )]
    configuration: StereoConfigurationForm,
) {
    let form = StereoAtomForm {
        configuration: configuration.clone(),
        constraints: raw.clone().into(),
    };
    let expected = if configuration.kind().is_none() {
        Ok(StereoAtomForm::default())
    } else {
        normalized
            .map(|value| StereoAtomForm {
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
fn test_stereo_atom_form_normalized_eq_constraints(#[case] payload: StereogenicityForm) {
    let form = StereoAtomForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: StereoAtomConstraintForm::Stereogenicity(payload.clone()).into(),
    };
    let top = StereoAtomForm::default();
    assert_ne!(form, top);
    assert!(form.normalized_eq(&top));
    assert!(top.normalized_eq(&form));
    assert!(form.is_undetermined());
    assert!(!form.is_ground());
    assert_eq!(form.constraints.stereogenicity(), &payload);
}

#[rustfmt::skip]
#[rstest]
#[case::undetermined(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined), false)]
#[case::ground(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), true)]
fn test_stereo_atom_form_is_ground(#[case] atom: StereoAtomForm, #[case] expected: bool) {
    assert_eq!(atom.is_ground(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::same_kind_narrows(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
    Some(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)))]
#[case::different_kind(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
    StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Undetermined), None)]
#[case::config_conflict(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), None)]
fn test_stereo_atom_form_meet(
    #[case] a: StereoAtomForm,
    #[case] b: StereoAtomForm,
    #[case] expected: Option<StereoAtomForm>,
) {
    assert_eq!(a.meet(&b), expected);
}

#[rstest]
#[case::literal(StereogenicityForm::Lit(Stereogenicity::Symmetric))]
#[case::contradiction(StereogenicityForm::lit_set([]))]
fn test_stereo_atom_form_meet_constraints(
    #[case] payload: StereogenicityForm,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(0))] coset: StereoCoset,
) {
    let kindless = StereoAtomForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: StereoAtomConstraintForm::Stereogenicity(payload).into(),
    };
    let kinded = StereoAtomForm {
        configuration: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset),
        constraints: StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        ))
        .into(),
    };
    assert_eq!(kindless.meet(&kindless), Some(StereoAtomForm::default()));
    assert_eq!(kindless.meet(&kinded), Some(kinded.clone()));
    assert_eq!(kinded.meet(&kindless), Some(kinded.clone()));
    assert!(kindless.matches(&kinded));
    assert!(!kinded.matches(&kindless));
    assert!(kindless.is_compatible(&kinded));
    assert!(kinded.is_compatible(&kindless));
}

#[rustfmt::skip]
#[rstest]
#[case::same_coset(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32))]
#[case::distinct_cosets_widen(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::Tetrahedral, 2u32), StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([1, 2])))]
fn test_stereo_atom_form_join(#[case] a: StereoAtomForm, #[case] b: StereoAtomForm, #[case] expected: StereoAtomForm) {
    assert_eq!(a.join(&b), Ok(expected));
}

#[rstest]
#[case::empty_payload(StereoConfigurationForm::Kinded(
    StereoKind::Tetrahedral,
    StereoCoset::Undetermined
))]
fn test_stereo_atom_form_meet_error(#[case] configuration: StereoConfigurationForm) {
    let kindless = StereoAtomForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([]))
            .into(),
    };
    let bottom = StereoAtomForm {
        configuration,
        constraints: StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([]))
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
    StereoKind::SquarePlanar,
    StereoCoset::Undetermined
))]
fn test_stereo_atom_form_join_constraints(#[case] configuration: StereoConfigurationForm) {
    let a = StereoAtomForm {
        configuration,
        constraints: StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        ))
        .into(),
    };
    let b = StereoAtomForm {
        configuration: StereoConfigurationForm::Kinded(
            StereoKind::Tetrahedral,
            StereoCoset::Undetermined,
        ),
        constraints: a.constraints.clone(),
    };
    assert_eq!(a.join(&b), Ok(StereoAtomForm::default()));
    assert_eq!(b.join(&a), Ok(StereoAtomForm::default()));
}

#[rustfmt::skip]
#[rstest]
#[case::same_kind_match(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), true)]
#[case::different_kind(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::SquarePlanar, 1u32), false)]
fn test_stereo_atom_form_matches(
    #[case] pattern: StereoAtomForm,
    #[case] target: StereoAtomForm,
    #[case] expected: bool,
) {
    assert_eq!(pattern.matches(&target), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::open_coset(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined))]
#[case::ground(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32))]
fn test_stereo_atom_form_into_concrete(#[case] atom: StereoAtomForm) {
    assert_eq!(atom.clone().into_concrete(), atom);
}

#[rstest]
#[case::apply(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32), StereoKind::Tetrahedral.involution(), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32))]
fn test_stereo_atom_form_apply(
    #[case] input: StereoAtomForm,
    #[case] permutation: Permutation,
    #[case] expected: StereoAtomForm,
) {
    assert_eq!(input.apply(permutation), Some(expected));
}

#[rstest]
#[case::tetrahedral(
    StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
    StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)
)]
fn test_stereo_atom_form_swap(#[case] input: StereoAtomForm, #[case] expected: StereoAtomForm) {
    assert_eq!(input.swap(), Some(expected));
}

#[rstest]
#[case::chiral(
    StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
    StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)
)]
fn test_stereo_atom_form_mirror(#[case] input: StereoAtomForm, #[case] expected: StereoAtomForm) {
    assert_eq!(input.mirror(), Some(expected));
}

#[rstest]
fn test_stereo_atom_form_reframe_by_constraints() {
    let action = Permutation::from_image(&[1, 0, 2, 3]);
    let input = StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
            permutation: LigandPermutation(Permutation::from_image(&[0, 2, 1, 3])),
            active: BooleanForm::Lit(true),
        }),
    );
    let expected = StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32).with_constraint(
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
            permutation: LigandPermutation(Permutation::from_image(&[2, 1, 0, 3])),
            active: BooleanForm::Lit(true),
        }),
    );

    assert_eq!(input.reframe_by(&action), Some(expected),);
}

#[rstest]
#[case::identity(
    [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
    [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
    0,
)]
#[case::transposition(
    [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
    [StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
    1,
)]
#[case::even_cycle(
    [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
    [StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
    0,
)]
#[case::virtual_explicit_swap(
    [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen)],
    [StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
    1,
)]
fn test_stereo_atom_form_reframe_by(
    #[case] before: [StereoLigand; 4],
    #[case] after: [StereoLigand; 4],
    #[case] expected_coset: u32,
) {
    let permutation = Permutation::between(&before, &after).expect("case has a unique frame");
    assert_eq!(
        StereoAtomForm::new(StereoKind::Tetrahedral, 0u32).reframe_by(&permutation),
        Some(StereoAtomForm::new(StereoKind::Tetrahedral, expected_coset)),
    );
}

/// A form's constructors are permissive, so a coset index outside its kind's coset space is
/// representable until molecule publication checks it. The frame action is the first operation
/// that requires the index to be in range, so it declines there.
#[rstest]
#[case::atom_coset_out_of_range(StereoKind::Tetrahedral, 2u32)]
#[case::atom_coset_far_out_of_range(StereoKind::SquarePlanar, 40u32)]
fn test_stereo_atom_form_reframe_by_error(#[case] kind: StereoKind, #[case] coset: u32) {
    let form = StereoAtomForm::new(kind, coset);
    assert_eq!(form.clone().reframe_by(&Permutation::identity(4)), None);
}

#[rstest]
fn test_stereo_atom_form_reframe_by_roundtrip() {
    let before = [
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
    ];
    let after = [
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
    ];
    let permutation = Permutation::between(&before, &after).expect("a unique frame change");
    let atom = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);
    assert_eq!(
        atom.clone()
            .reframe_by(&permutation)
            .and_then(|reframed| reframed.reframe_by(&permutation.inverse())),
        Some(atom),
    );
}

#[rustfmt::skip]
#[rstest]
#[case::absolute(
    StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() },
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
)]
#[case::relative(
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() },
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
)]
#[case::undetermined(
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() },
    StereoAtomForm::default(),
)]
#[case::explicit_open(
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Undetermined) }, ..Default::default() },
    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
)]
#[case::constraint_set(
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), ..Default::default() },
    StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
)]
#[case::constraint_remove(
    StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() },
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
)]
fn test_stereo_atom_form_update(
    #[case] atom: StereoAtomForm,
    #[case] update: StereoAtomUpdate,
    #[case] expected: StereoAtomForm,
) {
    assert_eq!(atom.update(&update), expected);
}

#[rstest]
#[case::empty(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32))]
fn test_stereo_atom_form_update_identity(#[case] atom: StereoAtomForm) {
    assert_eq!(atom.update(&StereoAtomUpdate::default()), atom);
}

#[rustfmt::skip]
#[rstest]
#[case::configuration_and_constraint(
    StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    StereoAtomForm::default(),
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
)]
#[case::constraint_context(
    StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
)]
fn test_stereo_atom_form_difference_to(
    #[case] atom: StereoAtomForm,
    #[case] other: StereoAtomForm,
    #[case] expected: StereoAtomUpdate,
) {
    assert_eq!(atom.difference_to(&other), expected);
}

#[rstest]
#[case::same(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32))]
fn test_stereo_atom_form_difference_to_identity(#[case] atom: StereoAtomForm) {
    assert_eq!(atom.difference_to(&atom), StereoAtomUpdate::default());
}
