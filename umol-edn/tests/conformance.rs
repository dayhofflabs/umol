//! EDN conformance tests — parser and streaming paths.
//!
//! Statement IDs (S2, S7, etc.) reference discussion/66-edn-spec-conformance-2026-04-01.md.

use std::borrow::Cow;
use std::collections::HashMap;
#[cfg(feature = "bignum")]
use std::str::FromStr;

use rstest::rstest;
use serde::Deserialize;
use umol_edn::config::{DuplicateKeyPolicy, ParseConfig, TagReaders};
use umol_edn::de::{from_str, from_str_with};
use umol_edn::edn::{Edn, Symbol};
use umol_edn::error::EdnError;
use umol_edn::{read_all_with, read_string_with, EdnMap, EdnSet};

fn cfg() -> ParseConfig {
    ParseConfig::default()
}

fn parse(input: &str) -> Result<Edn<'_>, EdnError> {
    read_string_with(input, &cfg())
}

fn stream<'a, T: Deserialize<'a>>(input: &'a str) -> Result<T, EdnError> {
    from_str_with(input, &cfg())
}

#[rstest]
#[case(" 1 ", Edn::Int(1))]
#[case("\t1\t", Edn::Int(1))]
#[case("\n1\n", Edn::Int(1))]
#[case("\r1\r", Edn::Int(1))]
#[case(",1,", Edn::Int(1))]
#[case(" \t\n\r,1", Edn::Int(1))]
fn test_whitespace(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(parse(input).unwrap(), expected);
}

#[rstest]
#[case(" 1 ", 1i64)]
#[case("\t1\t", 1)]
#[case("\n1\n", 1)]
#[case(",1,", 1)]
fn test_whitespace_streaming(#[case] input: &str, #[case] expected: i64) {
    assert_eq!(stream::<i64>(input).unwrap(), expected);
}

#[test]
fn test_formfeed_not_whitespace() {
    assert!(parse("\x0C1").is_err());
}

#[test]
fn test_comma_in_collections() {
    assert_eq!(
        parse("[1, 2, 3]").unwrap(),
        Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into())
    );
}

#[test]
fn test_comma_in_collections_streaming() {
    assert_eq!(stream::<Vec<i64>>("[1, 2, 3]").unwrap(), vec![1, 2, 3]);
}

#[test]
fn test_delimiters_no_whitespace() {
    assert_eq!(
        parse("[1[2]]").unwrap(),
        Edn::Vector(vec![Edn::Int(1), Edn::Vector(vec![Edn::Int(2)].into())].into())
    );
    assert_eq!(read_all_with("[1]2", &cfg()).unwrap().len(), 2);
}

#[test]
fn test_hash_is_not_delimiter() {
    assert_eq!(
        parse("foo#bar").unwrap(),
        Edn::Symbol(Symbol::new("foo#bar"))
    );
    assert_eq!(
        parse("[a#b]").unwrap(),
        Edn::Vector(vec![Edn::Symbol(Symbol::new("a#b"))].into())
    );
}

#[test]
fn test_nil() {
    assert_eq!(parse("nil").unwrap(), Edn::Nil);
}

#[test]
fn test_nil_streaming() {
    assert_eq!(stream::<Option<i64>>("nil").unwrap(), None);
}

#[rstest]
#[case("true", Edn::Bool(true))]
#[case("false", Edn::Bool(false))]
fn test_booleans(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(parse(input).unwrap(), expected);
}

#[rstest]
#[case("true", true)]
#[case("false", false)]
fn test_booleans_streaming(#[case] input: &str, #[case] expected: bool) {
    assert_eq!(stream::<bool>(input).unwrap(), expected);
}

#[rstest]
#[case("Nil")]
#[case("TRUE")]
#[case("False")]
fn test_case_sensitive(#[case] input: &str) {
    assert!(matches!(parse(input).unwrap(), Edn::Symbol(_)));
}

#[rstest]
#[case(r#""""#, "")]
#[case(r#""hello""#, "hello")]
#[case(r#""line\nbreak""#, "line\nbreak")]
#[case(r#""tab\there""#, "tab\there")]
#[case(r#""cr\rhere""#, "cr\rhere")]
#[case(r#""slash\\here""#, "slash\\here")]
#[case(r#""quote\"here""#, "quote\"here")]
fn test_string_escapes(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(
        parse(input).unwrap(),
        Edn::Str(Cow::Owned(expected.to_string()))
    );
}

#[rstest]
#[case(r#""hello""#, "hello")]
#[case(r#""line\nbreak""#, "line\nbreak")]
#[case(r#""tab\there""#, "tab\there")]
fn test_string_escapes_streaming(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(stream::<String>(input).unwrap(), expected);
}

#[rstest]
#[case(r#""\b""#)]
#[case(r#""\f""#)]
#[case(r#""\101""#)]
fn test_string_clojure_escapes_rejected(#[case] input: &str) {
    assert!(parse(input).is_err(), "Edn should reject {input}");
}

#[rstest]
#[case(r#""\u0041""#, "A")]
#[case(r#""\u00e9""#, "é")]
fn test_string_unicode_escape(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(
        parse(input).unwrap(),
        Edn::Str(Cow::Owned(expected.to_string()))
    );
    assert_eq!(stream::<String>(input).unwrap(), expected);
}

#[rstest]
#[case(r#""\b""#)]
#[case(r#""\f""#)]
fn test_string_clojure_escapes_rejected_streaming(#[case] input: &str) {
    assert!(
        stream::<String>(input).is_err(),
        "Edn streaming should reject {input}"
    );
}

#[rstest]
#[case(r#""\x""#)]
#[case(r#""\a""#)]
#[case(r#""\v""#)]
fn test_string_invalid_escape(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[test]
fn test_string_unterminated() {
    assert!(parse(r#""hello"#).is_err());
}

#[rstest]
#[case(r"\a", 'a')]
#[case(r"\Z", 'Z')]
#[case(r"\!", '!')]
#[case(r"\newline", '\n')]
#[case(r"\return", '\r')]
#[case(r"\space", ' ')]
#[case(r"\tab", '\t')]
#[case(r"\u0041", 'A')]
#[case(r"\u03B1", '\u{03B1}')]
fn test_characters(#[case] input: &str, #[case] expected: char) {
    assert_eq!(parse(input).unwrap(), Edn::Char(expected));
}

#[rstest]
#[case(r"\formfeed")]
#[case(r"\backspace")]
fn test_formfeed_backspace_rejected(#[case] input: &str) {
    assert!(parse(input).is_err(), "Edn should reject {input}");
}

#[rstest]
#[case(r"\ ")]
#[case(r"\abc")]
fn test_characters_error(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[test]
fn test_termination() {
    assert_eq!(
        parse(r"[\a 1]").unwrap(),
        Edn::Vector(vec![Edn::Char('a'), Edn::Int(1)].into())
    );
}

#[rstest]
#[case("foo", "foo")]
#[case("ns/name", "ns/name")]
#[case("my.ns/sym", "my.ns/sym")]
fn test_symbols(#[case] input: &str, #[case] name: &str) {
    assert_eq!(parse(input).unwrap(), Edn::Symbol(Symbol::new(name)));
}

#[test]
fn test_slash_alone() {
    assert_eq!(parse("/").unwrap(), Edn::Symbol(Symbol::new("/")));
}

#[rstest]
#[case(".", ".")]
#[case("*", "*")]
#[case("!", "!")]
#[case("_", "_")]
#[case("?", "?")]
#[case("$", "$")]
#[case("%", "%")]
#[case("&", "&")]
#[case("=", "=")]
#[case("<", "<")]
#[case(">", ">")]
fn test_start_chars(#[case] input: &str, #[case] name: &str) {
    assert_eq!(parse(input).unwrap(), Edn::Symbol(Symbol::new(name)));
}

#[test]
fn test_interior_chars() {
    assert_eq!(
        parse("foo+bar-baz#qux:quux'end").unwrap(),
        Edn::Symbol(Symbol::new("foo+bar-baz#qux:quux'end"))
    );
}

#[test]
fn test_sign_dot_first_char() {
    assert!(matches!(parse("+a").unwrap(), Edn::Symbol(_)));
    assert!(matches!(parse("-a").unwrap(), Edn::Symbol(_)));
    assert!(matches!(parse(".a").unwrap(), Edn::Symbol(_)));
    assert_eq!(parse("+1").unwrap(), Edn::Int(1));
    assert_eq!(parse("-1").unwrap(), Edn::Int(-1));
}

#[rstest]
#[case("ns/")]
#[case("/name")]
fn test_empty_prefix_or_name(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[rstest]
#[case("a/b/c")]
#[case("foo/1bar")]
fn test_symbol_slash_rejected(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[rstest]
#[case("foo/_bar")]
#[case("foo/bar")]
fn test_post_slash_valid(#[case] input: &str) {
    assert!(parse(input).is_ok());
}

#[rstest]
#[case(":foo", "foo")]
#[case(":a", "a")]
#[case(":ns/name", "ns/name")]
#[case(":my.ns/foo", "my.ns/foo")]
fn test_keywords(#[case] input: &str, #[case] name: &str) {
    assert_eq!(parse(input).unwrap(), Edn::keyword(name));
}

#[rstest]
#[case(":foo", "foo")]
#[case(":ns/name", "ns/name")]
fn test_keywords_streaming(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(stream::<String>(input).unwrap(), expected);
}

#[rstest]
#[case(":0")]
#[case(":0foo")]
fn test_digit_start_rejected(#[case] input: &str) {
    assert!(parse(input).is_err(), "Edn should reject {input}");
}

#[test]
fn test_hash_start_rejected() {
    assert!(parse(":#foo").is_err());
}

#[rstest]
#[case(":.foo")]
#[case(":+foo")]
#[case(":-foo")]
fn test_special_start_chars(#[case] input: &str) {
    assert!(parse(input).is_ok());
}

#[test]
fn test_namespace_must_be_symbol() {
    assert!(parse(":0/foo").is_err(), "Edn should reject :0/foo");
}

#[test]
fn test_post_slash_symbol_start() {
    assert!(parse(":foo/0bar").is_err());
    assert!(parse(":foo/#bar").is_err(), "Edn should reject :foo/#bar");
    assert!(parse(":foo/:bar").is_err(), "Edn should reject :foo/:bar");
    // Interior # and : fine:
    assert!(parse(":foo/bar#baz").is_ok());
    assert!(parse(":foo/bar:baz").is_ok());
    // Valid start chars after slash:
    assert!(parse(":foo/.bar").is_ok());
    assert!(parse(":foo/+bar").is_ok());
    assert!(parse(":foo/-bar").is_ok());
}

#[test]
fn test_double_colon_rejected() {
    assert!(parse("::foo").is_err(), "Edn should reject ::foo");
}

#[rstest]
#[case(":/")]
#[case(":/foo")]
fn test_invalid_slash(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[test]
fn test_bare_colon_error() {
    assert!(parse(": ").is_err());
}

#[rstest]
#[case("0", 0)]
#[case("1", 1)]
#[case("-1", -1)]
#[case("+1", 1)]
#[case("9223372036854775807", i64::MAX)]
#[case("-9223372036854775808", i64::MIN)]
fn test_integers(#[case] input: &str, #[case] expected: i64) {
    assert_eq!(parse(input).unwrap(), Edn::Int(expected));
}

#[rstest]
#[case("0", 0i64)]
#[case("1", 1)]
#[case("-1", -1)]
#[case("+5", 5)]
fn test_integers_streaming(#[case] input: &str, #[case] expected: i64) {
    assert_eq!(stream::<i64>(input).unwrap(), expected);
}

#[rstest]
#[case("007")]
#[case("00")]
fn test_leading_zeros_rejected(#[case] input: &str) {
    assert!(parse(input).is_err(), "Edn should reject {input}");
}

#[test]
fn test_negative_zero() {
    assert_eq!(parse("-0").unwrap(), Edn::Int(0));
}

#[cfg(not(feature = "bignum"))]
#[test]
fn test_overflow() {
    assert!(parse("9223372036854775808").is_err());
}

#[cfg(not(feature = "bignum"))]
#[rstest]
#[case("42N")]
#[case("0N")]
fn test_bigint_suffix_error(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[rstest]
#[case("1.0", 1.0)]
#[case("-3.14", -3.14)]
#[case("1e10", 1e10)]
#[case("1E10", 1e10)]
#[case("1e-10", 1e-10)]
#[case("1e+10", 1e10)]
#[case("1.5e3", 1.5e3)]
#[case("1.5E-3", 1.5e-3)]
#[case("+1.0", 1.0)]
fn test_floats(#[case] input: &str, #[case] expected: f64) {
    match parse(input).unwrap() {
        Edn::Float(v) => assert!((v - expected).abs() < f64::EPSILON),
        other => panic!("expected Float, got {other:?}"),
    }
}

#[rstest]
#[case("1.0", 1.0f64)]
#[case("-3.14", -3.14)]
#[case("1e10", 1e10)]
fn test_floats_streaming(#[case] input: &str, #[case] expected: f64) {
    let val: f64 = stream(input).unwrap();
    assert!((val - expected).abs() < 1e-10);
}

#[test]
fn test_negative_zero_float() {
    match parse("-0.0").unwrap() {
        Edn::Float(v) => {
            assert!(v == 0.0);
            assert!(v.is_sign_negative());
        }
        other => panic!("expected Float, got {other:?}"),
    }
}

#[rstest]
#[case("##NaN")]
#[case("##Inf")]
#[case("##-Inf")]
fn test_special_floats_rejected(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[rstest]
#[case("##NaN")]
#[case("##Inf")]
fn test_special_floats_rejected_streaming(#[case] input: &str) {
    assert!(stream::<f64>(input).is_err());
}

#[cfg(not(feature = "bignum"))]
#[rstest]
#[case("3.14M")]
fn test_bigdec_suffix_error(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[cfg(feature = "bignum")]
#[rstest]
#[case("42N", num_bigint::BigInt::from(42))]
#[case("0N", num_bigint::BigInt::from(0))]
#[case("-1N", num_bigint::BigInt::from(-1))]
#[case("+7N", num_bigint::BigInt::from(7))]
fn test_bigint_parse(#[case] input: &str, #[case] expected: num_bigint::BigInt) {
    assert_eq!(parse(input).unwrap(), Edn::BigInt(expected));
}

#[cfg(feature = "bignum")]
#[test]
fn test_bigint_parse_large() {
    let input = "99999999999999999999999999999N";
    let expected = num_bigint::BigInt::from_str("99999999999999999999999999999").unwrap();
    assert_eq!(parse(input).unwrap(), Edn::BigInt(expected));
}

#[cfg(feature = "bignum")]
#[rstest]
#[case("3.14M")]
#[case("0.0M")]
#[case("-1.5M")]
#[case("42M")]
fn test_bigdecimal_parse(#[case] input: &str) {
    assert!(matches!(parse(input).unwrap(), Edn::BigDecimal(_)));
}

#[cfg(feature = "bignum")]
#[rstest]
#[case("9223372036854775808")] // i64::MAX + 1
#[case("-9223372036854775809")] // i64::MIN - 1
fn test_bigint_overflow_promotes(#[case] input: &str) {
    assert!(matches!(parse(input).unwrap(), Edn::BigInt(_)));
}

#[cfg(feature = "bignum")]
#[test]
fn test_bigint_n_suffix_on_float_rejected() {
    assert!(parse("3.14N").is_err());
}

#[cfg(feature = "bignum")]
#[rstest]
#[case("42N")]
#[case("3.14M")]
#[case("99999999999999999999N")]
#[case("-0.001M")]
fn test_bignum_roundtrip(#[case] input: &str) {
    let parsed = parse(input).unwrap();
    let rendered = parsed.to_string();
    let reparsed = parse(&rendered).unwrap();
    assert_eq!(parsed, reparsed);
}

#[rstest]
#[case("()", Edn::List(vec![].into()))]
#[case("(1 2 3)", Edn::List(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()))]
#[case("[]", Edn::Vector(vec![].into()))]
#[case("[1 2 3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()))]
fn test_seqs(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(parse(input).unwrap(), expected);
}

#[test]
fn test_vector_streaming() {
    assert_eq!(stream::<Vec<i64>>("[1 2 3]").unwrap(), vec![1, 2, 3]);
}

#[test]
fn test_nested() {
    assert_eq!(
        parse("[[1 2] (3 4)]").unwrap(),
        Edn::Vector(
            vec![
                Edn::Vector(vec![Edn::Int(1), Edn::Int(2)].into()),
                Edn::List(vec![Edn::Int(3), Edn::Int(4)].into()),
            ]
            .into()
        )
    );
}

#[test]
fn test_map() {
    let mut expected = EdnMap::new();
    expected.insert(Edn::keyword("a"), Edn::Int(1));
    expected.insert(Edn::keyword("b"), Edn::Int(2));
    assert_eq!(parse("{:a 1 :b 2}").unwrap(), Edn::Map(expected));
}

#[test]
fn test_map_streaming() {
    let m: HashMap<String, i64> = stream("{:a 1 :b 2}").unwrap();
    assert_eq!(m.len(), 2);
    assert_eq!(m["a"], 1);
    assert_eq!(m["b"], 2);
}

#[test]
fn test_map_odd_elements_error() {
    assert!(parse("{:a 1 :b}").is_err());
}

#[test]
fn test_map_duplicate_key_error() {
    assert!(parse("{:a 1 :a 2}").is_err());
}

#[test]
fn test_map_duplicate_key_last_wins() {
    let config = ParseConfig {
        duplicate_keys: DuplicateKeyPolicy::LastWins,
        ..Default::default()
    };
    let result = read_string_with("{:a 1 :a 2}", &config).unwrap();
    let mut expected = EdnMap::new();
    expected.insert(Edn::keyword("a"), Edn::Int(2));
    assert_eq!(result, Edn::Map(expected));
}

#[test]
fn test_map_complex_keys() {
    let result = parse("{[1 2] :pair \"key\" :str}").unwrap();
    let mut expected = EdnMap::new();
    expected.insert(
        Edn::Vector(vec![Edn::Int(1), Edn::Int(2)].into()),
        Edn::keyword("pair"),
    );
    expected.insert(Edn::Str(Cow::Owned("key".to_string())), Edn::keyword("str"));
    assert_eq!(result, Edn::Map(expected));
}

#[test]
fn test_set() {
    let result = parse("#{1 2 3}").unwrap();
    if let Edn::Set(s) = result {
        assert_eq!(s.len(), 3);
        assert!(s.contains(&Edn::Int(1)));
        assert!(s.contains(&Edn::Int(2)));
        assert!(s.contains(&Edn::Int(3)));
    } else {
        panic!("expected Set");
    }
}

#[test]
fn test_set_empty() {
    assert_eq!(parse("#{}").unwrap(), Edn::Set(EdnSet::new()));
}

#[test]
fn test_qualified_tag() {
    let result = parse("#myapp/Person {:name \"Alice\"}").unwrap();
    let mut inner = EdnMap::new();
    inner.insert(
        Edn::keyword("name"),
        Edn::Str(Cow::Owned("Alice".to_string())),
    );
    assert_eq!(
        result,
        Edn::Tagged("myapp/Person".into(), Box::new(Edn::Map(inner)))
    );
}

#[test]
fn test_qualified_tag_streaming() {
    assert_eq!(
        stream::<Vec<i64>>("#my/tag [1 2 3]").unwrap(),
        vec![1, 2, 3]
    );
}

#[rstest]
#[case("#foo 1")]
#[case("#bar [1 2]")]
fn test_bare_tag_rejected(#[case] input: &str) {
    assert!(
        parse(input).is_err(),
        "Edn should reject bare tag in: {input}"
    );
}

#[test]
fn test_bare_tag_rejected_by_serde_default() {
    // The serde path honors the supplied config verbatim. The default config
    // rejects unknown tags, matching the native path.
    assert!(matches!(
        stream::<i64>("#foo 1"),
        Err(EdnError::InvalidTag { .. })
    ));
}

#[test]
fn test_bare_tag_passes_through_serde_with_permissive_config() {
    // Callers that want foreign types to use `#Variant` for enum dispatch or
    // to preserve arbitrary tagged literals must opt in explicitly.
    let mut config = ParseConfig::default();
    config.allow_unknown_tags = true;
    assert_eq!(from_str_with::<i64>("#foo 1", &config).unwrap(), 1);
}

#[test]
fn test_tag_without_value_error() {
    assert!(parse("#tag").is_err());
    assert!(parse("#myapp/Person").is_err());
}

#[rstest]
#[case(
    "#inst \"2024-01-01T00:00:00Z\"",
    "inst",
    Edn::Str(Cow::Borrowed("2024-01-01T00:00:00Z"))
)]
#[case(
    "#uuid \"f81d4fae-7dec-11d0-a765-00a0c91e6bf6\"",
    "uuid",
    Edn::Str(Cow::Borrowed("f81d4fae-7dec-11d0-a765-00a0c91e6bf6"))
)]
fn test_builtin_tags_parse_without_features(
    #[case] input: &str,
    #[case] tag: &str,
    #[case] inner: Edn<'_>,
) {
    // Built-in tags must parse as Tagged(...) even without chrono/uuid features.
    let result = parse(input).unwrap();
    assert_eq!(result, Edn::Tagged(tag.into(), Box::new(inner)));
}

#[cfg(feature = "chrono")]
#[test]
fn test_inst_accepted() {
    let val = parse("#inst \"2024-01-01T00:00:00Z\"").unwrap();
    assert_eq!(
        val,
        Edn::Tagged(
            "inst".into(),
            Box::new(Edn::Str(Cow::Borrowed("2024-01-01T00:00:00Z")))
        )
    );
}

#[cfg(feature = "uuid")]
#[test]
fn test_uuid_accepted() {
    let val = parse("#uuid \"f81d4fae-7dec-11d0-a765-00a0c91e6bf6\"").unwrap();
    assert_eq!(
        val,
        Edn::Tagged(
            "uuid".into(),
            Box::new(Edn::Str(Cow::Borrowed(
                "f81d4fae-7dec-11d0-a765-00a0c91e6bf6"
            )))
        )
    );
}

#[cfg(feature = "chrono")]
#[rstest]
#[case("#inst \"not-a-date\"")]
#[case("#inst \"2024-01-01\"")]
fn test_inst_invalid_rejected(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[cfg(feature = "chrono")]
#[test]
fn test_inst_non_string_rejected() {
    assert!(parse("#inst 123").is_err());
}

#[cfg(feature = "uuid")]
#[test]
fn test_uuid_invalid_rejected() {
    assert!(parse("#uuid \"not-a-uuid\"").is_err());
}

#[cfg(feature = "chrono")]
#[test]
fn test_inst_to_edn_roundtrip() {
    use umol_edn::read_string;
    let dt = chrono::DateTime::parse_from_rfc3339("2024-06-15T08:30:00+00:00").unwrap();
    let edn_val = umol_edn::tags::inst_to_edn(&dt);
    let printed = edn_val.to_string();
    let parsed = read_string(&printed).unwrap();
    assert_eq!(parsed, edn_val);
}

#[cfg(feature = "uuid")]
#[test]
fn test_uuid_to_edn_roundtrip() {
    use umol_edn::read_string;
    let id = uuid::Uuid::parse_str("f81d4fae-7dec-11d0-a765-00a0c91e6bf6").unwrap();
    let edn_val = umol_edn::tags::uuid_to_edn(&id);
    let printed = edn_val.to_string();
    let parsed = read_string(&printed).unwrap();
    assert_eq!(parsed, edn_val);
}

#[test]
fn test_custom_reader_dispatch() {
    fn double_reader(val: Edn) -> Result<Edn, EdnError> {
        match val {
            Edn::Int(n) => Ok(Edn::Int(n * 2)),
            _ => Err(EdnError::Custom("expected int".into())),
        }
    }
    let mut readers = TagReaders::default();
    readers.insert("double", double_reader);
    let config = ParseConfig {
        tag_readers: readers,
        ..Default::default()
    };
    // Bare tag "double" accepted because it's registered.
    let val = read_string_with("#double 5", &config).unwrap();
    assert_eq!(val, Edn::Int(10));
}

#[test]
fn test_custom_reader_error_propagation() {
    fn strict_reader(_val: Edn) -> Result<Edn, EdnError> {
        Err(EdnError::Custom("reader rejected value".into()))
    }
    let mut readers = TagReaders::default();
    readers.insert("strict", strict_reader);
    let config = ParseConfig {
        tag_readers: readers,
        ..Default::default()
    };
    assert!(read_string_with("#strict 1", &config).is_err());
}

#[rstest]
#[case("; comment\n1", Edn::Int(1))]
#[case("1 ; trailing", Edn::Int(1))]
#[case("; full line comment\n; another\n1", Edn::Int(1))]
fn test_comments(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(parse(input).unwrap(), expected);
}

#[test]
fn test_comments_streaming() {
    assert_eq!(stream::<i64>("; comment\n1").unwrap(), 1);
}

#[test]
fn test_comment_inside_vector() {
    assert_eq!(
        parse("[1 ; comment\n 2]").unwrap(),
        Edn::Vector(vec![Edn::Int(1), Edn::Int(2)].into())
    );
}

#[rstest]
#[case("[1 #_ 2 3]", vec![1, 3])]
#[case("[#_ #_ 1 2 3]", vec![3])]
#[case("[#_ #my/tag \"2024\" 1]", vec![1])]
fn test_discard(#[case] input: &str, #[case] expected: Vec<i64>) {
    let parsed = parse(input).unwrap();
    let items: Vec<i64> = match parsed {
        Edn::Vector(v) => v
            .iter()
            .map(|e| match e {
                Edn::Int(n) => *n,
                _ => panic!("expected int, got {e:?}"),
            })
            .collect(),
        other => panic!("expected vector, got {other:?}"),
    };
    assert_eq!(items, expected);
}

#[rstest]
#[case("[1 #_ 2 3]", vec![1, 3])]
#[case("[#_ #_ 1 2 3]", vec![3])]
fn test_discard_streaming(#[case] input: &str, #[case] expected: Vec<i64>) {
    assert_eq!(stream::<Vec<i64>>(input).unwrap(), expected);
}

#[test]
fn test_discard_in_map() {
    let result = parse("{:a #_ :skip 1}").unwrap();
    let mut expected = EdnMap::new();
    expected.insert(Edn::keyword("a"), Edn::Int(1));
    assert_eq!(result, Edn::Map(expected));
}

#[rstest]
#[case("#_")]
#[case("#_ 1")]
fn test_discard_error(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[test]
fn test_int_not_equal_to_float() {
    assert_ne!(Edn::Int(1), Edn::Float(1.0));
}

#[rstest]
#[case(vec![])]
#[case(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)])]
fn test_list_equals_vector(#[case] items: Vec<Edn<'static>>) {
    let list = Edn::List(items.clone().into());
    let vector = Edn::Vector(items.into());
    assert_eq!(list, vector);
}

#[test]
fn test_list_vector_hash_equal() {
    use std::hash::{DefaultHasher, Hash, Hasher};
    fn hash_of(v: &Edn<'_>) -> u64 {
        let mut h = DefaultHasher::new();
        v.hash(&mut h);
        h.finish()
    }
    let list = Edn::List(vec![Edn::Int(1)].into());
    let vector = Edn::Vector(vec![Edn::Int(1)].into());
    assert_eq!(hash_of(&list), hash_of(&vector));
}

#[test]
fn test_list_not_equal_to_vector_different_elements() {
    let list = Edn::List(vec![Edn::Int(1)].into());
    let vector = Edn::Vector(vec![Edn::Int(2)].into());
    assert_ne!(list, vector);
}

#[test]
fn test_parsed_list_equals_parsed_vector() {
    let list = parse("(1 2 3)").unwrap();
    let vector = parse("[1 2 3]").unwrap();
    assert_eq!(list, vector);
}

#[rstest]
#[case("nil")]
#[case("true")]
#[case("false")]
#[case("0")]
#[case("-1")]
#[case("9223372036854775807")]
#[case("1.0")]
#[case("-3.14")]
#[case("1.5e3")]
#[case(r#""hello""#)]
#[case(r#""line\nbreak""#)]
#[case(r"\a")]
#[case(r"\newline")]
#[case(r"\space")]
#[case(":foo")]
#[case(":ns/name")]
#[case("my-sym")]
#[case("/")]
#[case("()")]
#[case("(1 2 3)")]
#[case("[]")]
#[case("[1 2 3]")]
#[case("{}")]
#[case("{:a 1}")]
#[case("#{}")]
#[case("#{1 2 3}")]
#[case("#my/tag value")]
fn test_roundtrip(#[case] input: &str) {
    let parsed = parse(input).unwrap();
    let rendered = parsed.to_string();
    let reparsed = parse(&rendered).unwrap();
    assert_eq!(parsed, reparsed);
}

#[test]
fn test_error_trailing_content() {
    let err = parse("1 2").unwrap_err();
    assert!(matches!(err, EdnError::TrailingContent { .. }));
}

#[test]
fn test_error_trailing_content_streaming() {
    let err = stream::<i64>("1 2").unwrap_err();
    assert!(matches!(err, EdnError::TrailingContent { .. }));
}

#[test]
fn test_error_empty_input() {
    assert!(matches!(
        parse("").unwrap_err(),
        EdnError::UnexpectedEof { .. }
    ));
}

#[test]
fn test_error_whitespace_only() {
    assert!(matches!(
        parse("   ").unwrap_err(),
        EdnError::UnexpectedEof { .. }
    ));
}

#[rstest]
#[case("[1 2")]
#[case("{:a 1")]
#[case("(1 2")]
#[case("#{1 2")]
fn test_error_unclosed(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[test]
fn test_error_deep_nesting() {
    let depth = 100;
    let input = "[".repeat(depth) + &"]".repeat(depth);
    assert!(parse(&input).is_ok());
}

/// Fuzzer-found regressions: these inputs must not panic or hang.
#[rstest]
#[case("#_(#!V(\0\0\0\0\u{00ff}##")] // skip_atom infinite loop with null bytes in discard+collection
#[case("#_\"\\")] // skip_string panic with trailing backslash at EOF
fn test_error_deser_regression_vec(#[case] input: &str) {
    let _ = from_str::<Vec<i64>>(input);
}

/// Fuzzer-found panic slicing &str at non-char-boundary in \u escape parsing.
#[test]
fn test_error_deser_unicode_escape_multibyte_boundary() {
    let input = "\"\u{005c}u2\0`\u{07a0}";
    let _ = from_str::<String>(input);
}
