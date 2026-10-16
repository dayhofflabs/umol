use pretty_assertions::assert_eq;
use umol_edn::ParseError as EdnParseError;

use super::*;

#[rustfmt::skip]
#[rstest]
#[case::cis_trans_z("Ct0", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0))))]
#[case::cis_trans_e("Ct1", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))))]
#[case::open("Ct*", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined)))]
#[case::undetermined("*", StereoBondDsl(StereoBondForm::default()))]
#[case::no_canonicalization("Ct~1", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1))))))]
fn test_parse_stereo_bond(#[case] input: &str, #[case] expected: StereoBondDsl) {
    assert_eq!(parse_stereo_bond(input).unwrap(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::not_stereo("Ct!", ParseError::Syntax)]
#[case::trailing_after_coset("Ct1x", ParseError::TrailingInput("x".to_string()))]
#[case::unknown_predicate("Ct1#x", ParseError::UnknownStereoPredicate("#x".to_string()))]
#[case::duplicate("Ct1#g/#g=", ParseError::DuplicateStereoPredicate("#g".to_string()))]
#[case::undetermined_rejects_constraint("*#g/", ParseError::TrailingInput("#g/".to_string()))]
fn test_parse_stereo_bond_error(#[case] input: &str, #[case] expected: ParseError) {
    assert_eq!(parse_stereo_bond(input).unwrap_err(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::cis_trans_z(StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0))), "Ct0")]
#[case::open(StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined)), "Ct*")]
fn test_stereo_bond_dsl_to_string(#[case] form: StereoBondDsl, #[case] expected: &str) {
    assert_eq!(form.to_string(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::fluxionality_involution("Ct1#f~")]
#[case::topicity("Ct1#o(0,1)=")]
#[case::stereogenicity("Ct1#g/")]
fn test_stereo_bond_inline_render_identity(#[case] s: &str) {
    assert_eq!(parse_stereo_bond(s).unwrap().to_string(), s);
}

#[rustfmt::skip]
#[rstest]
#[case::string("\"Ct1\"", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))))]
#[case::keyword_z(":z", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0))))]
#[case::keyword_e(":e", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))))]
fn test_stereo_bond_dsl_from_edn(#[case] input: &str, #[case] expected: StereoBondDsl) {
    assert_eq!(StereoBondDsl::from_edn(&read_string(input).unwrap()).unwrap(), expected);
    assert_eq!(StereoBondDsl::from_edn_str(input).unwrap(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::unknown_keyword(":xyz", DeError::Custom("unknown stereo bond keyword :xyz".to_string()))]
#[case::wrong_type("5", DeError::TypeMismatch { expected: "string or stereo bond keyword", got: "int", path: Vec::new() })]
fn test_stereo_bond_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
    assert_eq!(StereoBondDsl::from_edn(&read_string(input).unwrap()).unwrap_err(), expected);
}

#[rstest]
#[case::unknown_keyword(":xyz", EdnError::De(DeError::Custom("unknown stereo bond keyword :xyz".to_string())))]
#[case::keyword_trailing(":z 0", EdnError::Parse(EdnParseError::TrailingContent { offset: 3 }))]
#[case::string_trailing("\"*\" 0", EdnError::Parse(EdnParseError::TrailingContent { offset: 4 }))]
fn test_stereo_bond_dsl_from_edn_str_error(#[case] input: &str, #[case] expected: EdnError) {
    assert_eq!(StereoBondDsl::from_edn_str(input), Err(expected));
}

#[rustfmt::skip]
#[rstest]
#[case::canonical_z(StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0))), ":z")]
#[case::canonical_e(StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))), ":e")]
#[case::open_string(StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined)), "\"Ct*\"")]
fn test_stereo_bond_dsl_to_edn(#[case] form: StereoBondDsl, #[case] expected: &str) {
    assert_eq!(form.to_edn(), read_string(expected).unwrap());
}

#[rstest]
#[case::symmetry(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }))]
#[case::fluxionality(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) }))]
#[case::topicity(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }))]
#[case::stereogenicity(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
    Stereogenicity::Stereogenic
)))]
#[case::contradiction(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])))]
fn test_stereo_bond_dsl_to_edn_kindless(#[case] constraint: StereoBondConstraintForm) {
    let form = StereoBondForm {
        constraints: constraint.clone().into(),
        ..Default::default()
    };
    let dsl = StereoBondDsl::from_ir(&form, &StereoBondDefaults::default());
    assert_eq!(dsl.to_string(), "*");
    assert_eq!(form.to_string(), "*");
    assert_eq!(dsl.to_edn(), Edn::Str("*".into()));
    assert_eq!(
        StereoBondDsl::from_edn(&dsl.to_edn()).unwrap().0,
        StereoBondForm::default()
    );
    assert_eq!(dsl.0, form);
    assert_eq!(form.constraints, constraint.into());
}

#[rstest]
#[case::z("z", Some("Ct0"))]
#[case::e("e", Some("Ct1"))]
#[case::unknown("xyz", None)]
fn test_expand_stereo_bond_keyword(#[case] name: &str, #[case] expected: Option<&str>) {
    assert_eq!(expand_stereo_bond_keyword(name), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::empty("", StereoBondUpdateDsl::default())]
#[case::undetermined("*", StereoBondUpdateDsl(StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() }))]
#[case::absolute("Ct1", StereoBondUpdateDsl(StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() }))]
#[case::relative("Ct", StereoBondUpdateDsl(StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() }))]
#[case::explicit_open("Ct*", StereoBondUpdateDsl(StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Undetermined) }, ..Default::default() }))]
fn test_parse_stereo_bond_update(#[case] input: &str, #[case] expected: StereoBondUpdateDsl) {
    assert_eq!(parse_stereo_bond_update(input).unwrap(), expected);
}

#[rstest]
#[case::undetermined_kind_with_constraint("*#g/", ParseError::TrailingInput("#g/".to_string()))]
fn test_parse_stereo_bond_update_error(#[case] input: &str, #[case] expected: ParseError) {
    assert_eq!(parse_stereo_bond_update(input).unwrap_err(), expected);
}

#[rstest]
#[case::undetermined_kind_with_constraint(
    "*#g/",
    ParseError::TrailingInput("#g/".to_string())
)]
fn test_stereo_bond_update_from_str_error(#[case] input: &str, #[case] expected: ParseError) {
    assert_eq!(input.parse::<StereoBondUpdate>().unwrap_err(), expected);
}

#[rstest]
#[case::empty("")]
#[case::undetermined("*")]
#[case::absolute("Ct1")]
#[case::relative("Ct")]
#[case::explicit_open("Ct*")]
#[case::topicity_removal("Ct#o(0,1)*")]
#[case::stereogenicity_change("Ct#g/")]
fn test_stereo_bond_update_dsl_display_roundtrip(#[case] input: &str) {
    let dsl = parse_stereo_bond_update(input).unwrap();
    assert_eq!(dsl.to_string(), input);
    assert_eq!(parse_stereo_bond_update(&dsl.to_string()).unwrap(), dsl);
}

#[rstest]
#[case::symmetry(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }))]
#[case::fluxionality(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) }))]
#[case::topicity(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }))]
#[case::stereogenicity(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
    Stereogenicity::Stereogenic
)))]
#[case::contradiction(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])))]
fn test_stereo_bond_update_dsl_to_edn_kindless(
    #[case] constraint: StereoBondConstraintForm,
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
    let update = StereoBondUpdate {
        configuration: configuration.clone(),
        constraints: constraint.clone().into(),
    };
    let dsl = StereoBondUpdateDsl::from_ir(&update, &());
    assert_eq!(dsl.to_string(), expected);
    assert_eq!(update.to_string(), expected);
    assert_eq!(dsl.to_edn(), Edn::Str(expected.into()));
    assert_eq!(
        StereoBondUpdateDsl::from_edn(&dsl.to_edn()).unwrap().0,
        StereoBondUpdate {
            configuration,
            ..Default::default()
        }
    );
    assert_eq!(dsl.into_ir(&()), update);
    assert_eq!(update.constraints, constraint.into());
}
