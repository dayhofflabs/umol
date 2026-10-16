use pretty_assertions::assert_eq;
use umol_edn::ParseError as EdnParseError;

use super::*;

#[rustfmt::skip]
#[rstest]
#[case::tetrahedral_ccw("Th0", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))))]
#[case::tetrahedral_cw("Th1", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))))]
#[case::open("Th*", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined)))]
#[case::square_planar("Sp2", StereoAtomDsl(StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(2))))]
#[case::octahedral("Oh6", StereoAtomDsl(StereoAtomForm::new(StereoKind::Octahedral, StereoCoset::Lit(6))))]
#[case::undetermined("*", StereoAtomDsl(StereoAtomForm::default()))]
#[case::no_canonicalization("Th~1", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1))))))]
fn test_parse_stereo_atom(#[case] input: &str, #[case] expected: StereoAtomDsl) {
    assert_eq!(parse_stereo_atom(input).unwrap(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::not_stereo("Th!", ParseError::Syntax)]
#[case::trailing_after_coset("Th1x", ParseError::TrailingInput("x".to_string()))]
#[case::unknown_predicate("Th1#x", ParseError::UnknownStereoPredicate("#x".to_string()))]
#[case::duplicate("Th1#g/#g=", ParseError::DuplicateStereoPredicate("#g".to_string()))]
#[case::undetermined_rejects_constraint("*#g/", ParseError::TrailingInput("#g/".to_string()))]
#[case::overlapping_ligand_permutation("Th1#p(0,1)(1,2)", ParseError::InvalidValue("cycle point 1 occurs more than once".to_string()))]
#[case::repeated_ligand_position("Th1#f(0,1,0)", ParseError::InvalidValue("cycle point 0 occurs more than once".to_string()))]
#[case::ligand_position_out_of_range("Th1#p(0,4)", ParseError::InvalidValue("cycle point 4 at cycle 0, position 1 is outside 0..4".to_string()))]
fn test_parse_stereo_atom_error(#[case] input: &str, #[case] expected: ParseError) {
    assert_eq!(parse_stereo_atom(input).unwrap_err(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::tetrahedral_ccw(StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))), "Th0")]
#[case::open(StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined)), "Th*")]
#[case::square_planar(StereoAtomDsl(StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(2))), "Sp2")]
#[case::octahedral(StereoAtomDsl(StereoAtomForm::new(StereoKind::Octahedral, StereoCoset::Lit(6))), "Oh6")]
fn test_stereo_atom_dsl_to_string(#[case] form: StereoAtomDsl, #[case] expected: &str) {
    assert_eq!(form.to_string(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::fluxionality("Th1#f(0,1,2)")]
#[case::ligand_symmetry_involution("Th1#p~")]
#[case::ligand_symmetry_not_present("Th1#p~!")]
#[case::ligand_symmetry_explicit("Th1#p(0,1,2)")]
#[case::topicity("Th1#o(0,1)=")]
#[case::topicity_negated("Th1#o(0,1)!'")]
#[case::topicity_lit_set("Th1#o(0,1){=,'}")]
#[case::topicity_not_set("Th1#o(0,1)!{=,'}")]
#[case::stereogenicity("Th1#g/")]
#[case::stereogenicity_lit_set("Th1#g{=,/}")]
#[case::multiple("Th1#f(0,1,2)#o(0,1)=#g/")]
fn test_stereo_atom_inline_render_identity(#[case] s: &str) {
    assert_eq!(parse_stereo_atom(s).unwrap().to_string(), s);
}

#[rustfmt::skip]
#[rstest]
#[case::topicity_open("Th1#o(0,1)*", "Th1")]
#[case::stereogenicity_open("Th1#g*", "Th1")]
fn test_stereo_atom_inline_render(#[case] input: &str, #[case] canonical: &str) {
    assert_eq!(parse_stereo_atom(input).unwrap().to_string(), canonical);
}

#[rustfmt::skip]
#[rstest]
#[case::fluxionality("Th1#f(0,1,2)",
    StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 2, 0, 3])), active: BooleanForm::Lit(true) }))]
#[case::ligand_symmetry("Th1#p(0,1,2)",
    StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 2, 0, 3])), orientation: Orientation::Proper },
        invariant: BooleanForm::Lit(true) }))]
#[case::ligand_symmetry_absent("Th1#p(0,1,2)!",
    StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 2, 0, 3])), orientation: Orientation::Proper },
        invariant: BooleanForm::Lit(false) }))]
#[case::fluxionality_absent("Th1#f(0,1,2)!",
    StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 2, 0, 3])), active: BooleanForm::Lit(false) }))]
#[case::topicity_negated("Th1#o(0,1)!'",
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Enantiotopic])) }))]
#[case::topicity_lit_set("Th1#o(0,1){=,'}",
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])) }))]
#[case::topicity_not_set("Th1#o(0,1)!{=,'}",
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])) }))]
#[case::topicity_open("Th1#o(0,1)*",
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }))]
#[case::stereogenicity("Th1#g/",
    StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)))]
fn test_stereo_atom_predicate(#[case] input: &str, #[case] expected: StereoAtomConstraintForm) {
    let dsl = parse_stereo_atom(input).unwrap();
    assert_eq!(dsl.0.constraints.iter().cloned().collect::<Vec<_>>(), vec![expected]);
}

#[rstest]
fn test_stereo_atom_predicate_involution() {
    let dsl = parse_stereo_atom("Th1#p~").unwrap();
    let expected = StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation {
            permutation: LigandPermutation(StereoKind::Tetrahedral.involution()),
            orientation: Orientation::Improper,
        },
        invariant: BooleanForm::Lit(true),
    });
    assert_eq!(
        dsl.0.constraints.iter().cloned().collect::<Vec<_>>(),
        vec![expected],
    );
}

#[rustfmt::skip]
#[rstest]
#[case::string("\"Th1\"", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))))]
#[case::keyword_ccw(":ccw", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))))]
#[case::keyword_cw(":cw", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))))]
#[case::string_square_planar("\"Sp2\"", StereoAtomDsl(StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(2))))]
fn test_stereo_atom_dsl_from_edn(#[case] input: &str, #[case] expected: StereoAtomDsl) {
    assert_eq!(StereoAtomDsl::from_edn(&read_string(input).unwrap()).unwrap(), expected);
    assert_eq!(StereoAtomDsl::from_edn_str(input).unwrap(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::unknown_keyword(":xyz", DeError::Custom("unknown stereo atom keyword :xyz".to_string()))]
#[case::wrong_type("5", DeError::TypeMismatch { expected: "string or stereo atom keyword", got: "int", path: Vec::new() })]
fn test_stereo_atom_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
    assert_eq!(StereoAtomDsl::from_edn(&read_string(input).unwrap()).unwrap_err(), expected);
}

#[rstest]
#[case::unknown_keyword(":xyz", EdnError::De(DeError::Custom("unknown stereo atom keyword :xyz".to_string())))]
#[case::keyword_trailing(":ccw 0", EdnError::Parse(EdnParseError::TrailingContent { offset: 5 }))]
#[case::string_trailing("\"*\" 0", EdnError::Parse(EdnParseError::TrailingContent { offset: 4 }))]
fn test_stereo_atom_dsl_from_edn_str_error(#[case] input: &str, #[case] expected: EdnError) {
    assert_eq!(StereoAtomDsl::from_edn_str(input), Err(expected));
}

#[rustfmt::skip]
#[rstest]
#[case::canonical_ccw(StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))), ":ccw")]
#[case::canonical_cw(StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))), ":cw")]
#[case::open_string(StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined)), "\"Th*\"")]
#[case::non_tetrahedral_string(StereoAtomDsl(StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(1))), "\"Sp1\"")]
fn test_stereo_atom_dsl_to_edn(#[case] form: StereoAtomDsl, #[case] expected: &str) {
    assert_eq!(form.to_edn(), read_string(expected).unwrap());
}

#[rstest]
#[case::symmetry(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }))]
#[case::fluxionality(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) }))]
#[case::topicity(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }))]
#[case::stereogenicity(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
    Stereogenicity::Stereogenic
)))]
#[case::contradiction(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])))]
fn test_stereo_atom_dsl_to_edn_kindless(#[case] constraint: StereoAtomConstraintForm) {
    let form = StereoAtomForm {
        constraints: constraint.clone().into(),
        ..Default::default()
    };
    let dsl = StereoAtomDsl::from_ir(&form, &StereoAtomDefaults::default());
    assert_eq!(dsl.to_string(), "*");
    assert_eq!(form.to_string(), "*");
    assert_eq!(dsl.to_edn(), Edn::Str("*".into()));
    assert_eq!(
        StereoAtomDsl::from_edn(&dsl.to_edn()).unwrap().0,
        StereoAtomForm::default()
    );
    assert_eq!(dsl.0, form);
    assert_eq!(form.constraints, constraint.into());
}

#[rstest]
fn test_stereo_atom_dsl_into_ir() {
    let form = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined);
    assert_eq!(
        StereoAtomDsl(form.clone()).into_ir(&StereoAtomDefaults::default()),
        form
    );
    assert_eq!(
        StereoAtomDsl::from_ir(&form, &StereoAtomDefaults::default()),
        StereoAtomDsl(form)
    );
}

#[rstest]
#[case::ccw("ccw", Some("Th0"))]
#[case::cw("cw", Some("Th1"))]
#[case::unknown("xyz", None)]
fn test_expand_stereo_atom_keyword(#[case] name: &str, #[case] expected: Option<&str>) {
    assert_eq!(expand_stereo_atom_keyword(name), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::empty("", StereoAtomUpdateDsl::default())]
#[case::undetermined("*", StereoAtomUpdateDsl(StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() }))]
#[case::absolute("Th1", StereoAtomUpdateDsl(StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() }))]
#[case::relative("Th", StereoAtomUpdateDsl(StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() }))]
#[case::explicit_open("Th*", StereoAtomUpdateDsl(StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Undetermined) }, ..Default::default() }))]
fn test_parse_stereo_atom_update(#[case] input: &str, #[case] expected: StereoAtomUpdateDsl) {
    assert_eq!(parse_stereo_atom_update(input).unwrap(), expected);
}

#[rstest]
#[case::undetermined_kind_with_constraint("*#o(0,1)=", ParseError::TrailingInput("#o(0,1)=".to_string()))]
fn test_parse_stereo_atom_update_error(#[case] input: &str, #[case] expected: ParseError) {
    assert_eq!(parse_stereo_atom_update(input).unwrap_err(), expected);
}

#[rstest]
#[case::undetermined_kind_with_constraint(
    "*#o(0,1)=",
    ParseError::TrailingInput("#o(0,1)=".to_string())
)]
fn test_stereo_atom_update_from_str_error(#[case] input: &str, #[case] expected: ParseError) {
    assert_eq!(input.parse::<StereoAtomUpdate>().unwrap_err(), expected);
}

#[rstest]
#[case::empty("")]
#[case::undetermined("*")]
#[case::absolute("Th1")]
#[case::relative("Th")]
#[case::explicit_open("Th*")]
#[case::topicity_removal("Th#o(0,1)*")]
#[case::topicity_change("Th#o(0,1)/")]
#[case::ligand_symmetry_removal("Th#p(0,1)*")]
#[case::ligand_symmetry("Th#p(0,1)")]
fn test_stereo_atom_update_dsl_display_roundtrip(#[case] input: &str) {
    let dsl = parse_stereo_atom_update(input).unwrap();
    assert_eq!(dsl.to_string(), input);
    assert_eq!(parse_stereo_atom_update(&dsl.to_string()).unwrap(), dsl);
}

#[rstest]
#[case::symmetry(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }))]
#[case::fluxionality(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) }))]
#[case::topicity(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }))]
#[case::stereogenicity(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
    Stereogenicity::Stereogenic
)))]
#[case::contradiction(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])))]
fn test_stereo_atom_update_dsl_to_edn_kindless(
    #[case] constraint: StereoAtomConstraintForm,
    #[values(
        StereoConfigurationUpdate::Unchanged,
        StereoConfigurationUpdate::Undetermined
    )]
    configuration: StereoConfigurationUpdate,
) {
    let expected = match configuration {
        StereoConfigurationUpdate::Unchanged => "",
        StereoConfigurationUpdate::Undetermined => "*",
        _ => unreachable!(),
    };
    let update = StereoAtomUpdate {
        configuration: configuration.clone(),
        constraints: constraint.clone().into(),
    };
    let dsl = StereoAtomUpdateDsl::from_ir(&update, &());
    assert_eq!(dsl.to_string(), expected);
    assert_eq!(update.to_string(), expected);
    assert_eq!(dsl.to_edn(), Edn::Str(expected.into()));
    assert_eq!(
        StereoAtomUpdateDsl::from_edn(&dsl.to_edn()).unwrap().0,
        StereoAtomUpdate {
            configuration,
            ..Default::default()
        }
    );
    assert_eq!(dsl.into_ir(&()), update);
    assert_eq!(update.constraints, constraint.into());
}
