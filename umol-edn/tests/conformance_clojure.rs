//! Clojure dialect conformance tests — parser and streaming paths.
//!
//! Every test uses `Dialect::Clojure` (matches clj parser behavior).
//! Statement IDs (S2, S7, etc.) reference discussion/66-edn-spec-conformance-2026-04-01.md.

use std::borrow::Cow;
use std::collections::HashMap;

use rstest::rstest;
use serde::Deserialize;
use umol_edn::config::{Dialect, ParseConfig};
use umol_edn::edn::{Edn, Symbol};
use umol_edn::error::EdnError;
use umol_edn::{from_str_with, read_all_with, read_string_with, EdnMap};

fn cfg() -> ParseConfig {
    ParseConfig {
        dialect: Dialect::Clojure,
        ..Default::default()
    }
}

fn parse(input: &str) -> Result<Edn<'_>, EdnError> {
    read_string_with(input, &cfg())
}

fn stream<'a, T: Deserialize<'a>>(input: &'a str) -> Result<T, EdnError> {
    from_str_with(input, &cfg())
}

// ============================================================================
// S2: Whitespace
// ============================================================================

#[rstest]
#[case(" 1 ", 1i64)]
#[case("\t1\t", 1)]
#[case(",1,", 1)]
fn test_s2_whitespace_streaming(#[case] input: &str, #[case] expected: i64) {
    assert_eq!(stream::<i64>(input).unwrap(), expected);
}

// ============================================================================
// S3: Delimiters
// ============================================================================

#[test]
fn test_s3_delimiters_no_whitespace() {
    assert_eq!(
        parse("[1[2]]").unwrap(),
        Edn::Vector(vec![Edn::Int(1), Edn::Vector(vec![Edn::Int(2)])])
    );
    assert_eq!(read_all_with("[1]2", &cfg()).unwrap().len(), 2);
}

// ============================================================================
// S7: Strings — Clojure extends with \b, \f, \uNNNN, octal
// ============================================================================

#[rstest]
#[case(r#""a\tb""#, "a\tb")]
#[case(r#""a\rb""#, "a\rb")]
#[case(r#""a\nb""#, "a\nb")]
#[case(r#""a\\b""#, "a\\b")]
#[case(r#""a\"b""#, "a\"b")]
fn test_s7_string_spec_escapes(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(parse(input).unwrap(), Edn::Str(Cow::Owned(expected.to_string())));
}

#[test]
fn test_s7_string_backspace_formfeed() {
    assert_eq!(
        parse(r#""\b""#).unwrap(),
        Edn::Str(Cow::Owned("\u{0008}".to_string()))
    );
    assert_eq!(
        parse(r#""\f""#).unwrap(),
        Edn::Str(Cow::Owned("\u{000C}".to_string()))
    );
}

#[test]
fn test_s7_string_backspace_formfeed_streaming() {
    assert_eq!(stream::<String>(r#""\b""#).unwrap(), "\u{0008}");
    assert_eq!(stream::<String>(r#""\f""#).unwrap(), "\u{000C}");
}

#[test]
fn test_s7_string_unicode_escape() {
    assert_eq!(
        parse(r#""\u0041""#).unwrap(),
        Edn::Str(Cow::Owned("A".to_string()))
    );
    assert_eq!(
        parse(r#""\u03B1""#).unwrap(),
        Edn::Str(Cow::Owned("\u{03B1}".to_string()))
    );
}

#[test]
fn test_s7_string_octal() {
    assert_eq!(
        parse(r#""\0""#).unwrap(),
        Edn::Str(Cow::Owned("\0".to_string()))
    );
    assert_eq!(
        parse(r#""\101""#).unwrap(),
        Edn::Str(Cow::Owned("A".to_string()))
    );
    assert_eq!(
        parse(r#""\377""#).unwrap(),
        Edn::Str(Cow::Owned("\u{00FF}".to_string()))
    );
}

#[rstest]
#[case(r#""hello""#, "hello")]
#[case(r#""line\nbreak""#, "line\nbreak")]
#[case(r#""tab\there""#, "tab\there")]
fn test_s7_string_streaming(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(stream::<String>(input).unwrap(), expected);
}

// ============================================================================
// S8: Characters — Clojure extends with \formfeed, \backspace
// ============================================================================

#[rstest]
#[case(r"\a", 'a')]
#[case(r"\newline", '\n')]
#[case(r"\return", '\r')]
#[case(r"\space", ' ')]
#[case(r"\tab", '\t')]
#[case(r"\u0041", 'A')]
fn test_s8_characters(#[case] input: &str, #[case] expected: char) {
    assert_eq!(parse(input).unwrap(), Edn::Char(expected));
}

#[test]
fn test_s8_formfeed_backspace() {
    assert_eq!(parse(r"\formfeed").unwrap(), Edn::Char('\u{000C}'));
    assert_eq!(parse(r"\backspace").unwrap(), Edn::Char('\u{0008}'));
}

// ============================================================================
// S9: Symbols — Clojure accepts multiple slashes
// ============================================================================

#[rstest]
#[case("foo", "foo")]
#[case("ns/name", "ns/name")]
fn test_s9_symbols(#[case] input: &str, #[case] name: &str) {
    assert_eq!(parse(input).unwrap(), Edn::Symbol(Symbol::new(name)));
}

#[test]
fn test_s9d_multiple_slashes_accepted() {
    assert!(parse("a/b/c").is_ok());
}

#[rstest]
#[case("ns/")]
#[case("/name")]
fn test_s9d_empty_prefix_or_name(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[test]
fn test_s9e_post_slash_digit_rejected() {
    // Clojure also rejects digits after slash.
    assert!(parse("foo/1bar").is_err());
}

// ============================================================================
// S10: Keywords — Clojure accepts digit-start keywords
// ============================================================================

#[rstest]
#[case(":0", "0")]
#[case(":1", "1")]
#[case(":123abc", "123abc")]
fn test_s10_digit_keywords(#[case] input: &str, #[case] name: &str) {
    assert_eq!(parse(input).unwrap(), Edn::keyword(name));
}

#[test]
fn test_s10_namespace_digit_start() {
    assert!(parse(":0/foo").is_ok());
}

#[test]
fn test_s10_post_slash_lenient() {
    // Clojure accepts # and : after slash.
    assert!(parse(":foo/#bar").is_ok());
    assert!(parse(":foo/:bar").is_ok());
}

#[rstest]
#[case(":foo", "foo")]
#[case(":ns/name", "ns/name")]
fn test_s10_keywords_streaming(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(stream::<String>(input).unwrap(), expected);
}

// ============================================================================
// S11: Integers — Clojure accepts leading zeros
// ============================================================================

#[rstest]
#[case("0", 0)]
#[case("1", 1)]
#[case("-1", -1)]
#[case("+1", 1)]
#[case("007", 7)]
#[case("0000", 0)]
#[case("9223372036854775807", i64::MAX)]
#[case("-9223372036854775808", i64::MIN)]
fn test_s11_integers(#[case] input: &str, #[case] expected: i64) {
    assert_eq!(parse(input).unwrap(), Edn::Int(expected));
}

#[rstest]
#[case("0", 0i64)]
#[case("-1", -1)]
#[case("+5", 5)]
fn test_s11_integers_streaming(#[case] input: &str, #[case] expected: i64) {
    assert_eq!(stream::<i64>(input).unwrap(), expected);
}

// ============================================================================
// S12: Floats — Clojure supports ##NaN, ##Inf, ##-Inf
// ============================================================================

#[rstest]
#[case("1.0", 1.0)]
#[case("-3.14", -3.14)]
#[case("1e10", 1e10)]
fn test_s12_floats(#[case] input: &str, #[case] expected: f64) {
    match parse(input).unwrap() {
        Edn::Float(v) => assert!((v - expected).abs() < f64::EPSILON),
        other => panic!("expected Float, got {other:?}"),
    }
}

#[test]
fn test_s12_special_floats() {
    match parse("##NaN").unwrap() {
        Edn::Float(v) => assert!(v.is_nan()),
        other => panic!("expected NaN, got {other:?}"),
    }
    match parse("##Inf").unwrap() {
        Edn::Float(v) => assert!(v == f64::INFINITY),
        other => panic!("expected Inf, got {other:?}"),
    }
    match parse("##-Inf").unwrap() {
        Edn::Float(v) => assert!(v == f64::NEG_INFINITY),
        other => panic!("expected -Inf, got {other:?}"),
    }
}

// NOTE: ##NaN/##Inf only handled in deserialize_any, not deserialize_f64.
// Streaming f64 deserialization of special floats does not work.

#[rstest]
#[case("1.0", 1.0f64)]
#[case("-3.14", -3.14)]
fn test_s12_floats_streaming(#[case] input: &str, #[case] expected: f64) {
    let val: f64 = stream(input).unwrap();
    assert!((val - expected).abs() < 1e-10);
}

// ============================================================================
// S13-S16: Collections
// ============================================================================

#[test]
fn test_s14_vector_streaming() {
    assert_eq!(stream::<Vec<i64>>("[1 2 3]").unwrap(), vec![1, 2, 3]);
}

#[test]
fn test_s15_map_streaming() {
    let m: HashMap<String, i64> = stream("{:a 1 :b 2}").unwrap();
    assert_eq!(m.len(), 2);
    assert_eq!(m["a"], 1);
    assert_eq!(m["b"], 2);
}

#[test]
fn test_s15_struct_streaming() {
    #[derive(Deserialize, Debug, PartialEq)]
    struct Point {
        x: i64,
        y: i64,
    }
    let p: Point = stream("{:x 1 :y 2}").unwrap();
    assert_eq!(p, Point { x: 1, y: 2 });
}

// ============================================================================
// S17: Tags — Clojure accepts bare tags
// ============================================================================

#[test]
fn test_s17_qualified_tag() {
    let result = parse("#myapp/Person {:name \"Alice\"}").unwrap();
    let mut inner = EdnMap::new();
    inner.insert(Edn::keyword("name"), Edn::Str(Cow::Owned("Alice".to_string())));
    assert_eq!(result, Edn::Tagged("myapp/Person".into(), Box::new(Edn::Map(inner))));
}

#[test]
fn test_s17c_bare_tag_accepted() {
    let val = parse("#foo 1").unwrap();
    assert_eq!(val, Edn::Tagged("foo".into(), Box::new(Edn::Int(1))));
}

// NOTE: streaming tag stripping only works through deserialize_any, not typed
// methods (deserialize_seq, etc.). #my/tag [1 2 3] → Vec<i64> does not work.

// ============================================================================
// S19: Comments
// ============================================================================

#[test]
fn test_s19_comments_streaming() {
    assert_eq!(stream::<i64>("; comment\n1").unwrap(), 1);
}

// ============================================================================
// S20: Discard
// ============================================================================

#[test]
fn test_s20_discard() {
    assert_eq!(
        parse("[1 #_ 2 3]").unwrap(),
        Edn::Vector(vec![Edn::Int(1), Edn::Int(3)])
    );
}

#[test]
fn test_s20_discard_nested() {
    assert_eq!(
        parse("[#_ #_ 1 2 3]").unwrap(),
        Edn::Vector(vec![Edn::Int(3)])
    );
}

#[test]
fn test_s20_discard_tagged_literal() {
    assert_eq!(
        parse("[#_ #inst \"2024\" 1]").unwrap(),
        Edn::Vector(vec![Edn::Int(1)])
    );
}

#[test]
fn test_s20_discard_in_map() {
    let result = parse("{:a #_ :skip 1}").unwrap();
    let mut expected = EdnMap::new();
    expected.insert(Edn::keyword("a"), Edn::Int(1));
    assert_eq!(result, Edn::Map(expected));
}

#[test]
fn test_s20_discard_streaming() {
    assert_eq!(stream::<Vec<i64>>("[1 #_ 2 3]").unwrap(), vec![1, 3]);
}

// ============================================================================
// Round-trip
// ============================================================================

#[rstest]
#[case("nil")]
#[case("true")]
#[case("false")]
#[case("0")]
#[case("-1")]
#[case("1.0")]
#[case(r#""hello""#)]
#[case(r"\a")]
#[case(":foo")]
#[case("my-sym")]
#[case("/")]
#[case("[]")]
#[case("[1 2 3]")]
#[case("{}")]
#[case("{:a 1}")]
#[case("#{}")]
#[case("#tag value")]
fn test_roundtrip(#[case] input: &str) {
    let parsed = parse(input).unwrap();
    let rendered = parsed.to_string();
    let reparsed = parse(&rendered).unwrap();
    assert_eq!(parsed, reparsed);
}

// ============================================================================
// Errors
// ============================================================================

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
    assert!(matches!(parse("").unwrap_err(), EdnError::UnexpectedEof { .. }));
}

#[rstest]
#[case("[1 2")]
#[case("{:a 1")]
#[case("(1 2")]
#[case("#{1 2")]
fn test_error_unclosed(#[case] input: &str) {
    assert!(parse(input).is_err());
}
