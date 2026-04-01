//! EDN conformance suite.
//!
//! Systematic edge-case tests derived from `spec/edn-spec.md`.
//! Each section number maps to a spec section.

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};

use rstest::rstest;
use umol_edn::config::{Dialect, DuplicateKeyPolicy, ParseConfig};
use umol_edn::edn::{Edn, Symbol};
use umol_edn::{read_string, read_string_with, EdnError};

fn clojure_config() -> ParseConfig {
    ParseConfig {
        dialect: Dialect::Clojure,
        ..Default::default()
    }
}

fn edn_config() -> ParseConfig {
    ParseConfig {
        dialect: Dialect::Edn,
        ..Default::default()
    }
}

// -- Section 1: Whitespace and commas --------------------------------------

#[rstest]
#[case(" 1 ", Edn::Int(1))]
#[case("\t1\t", Edn::Int(1))]
#[case("\n1\n", Edn::Int(1))]
#[case("\r1\r", Edn::Int(1))]
#[case(",1,", Edn::Int(1))]
#[case(" \t\n\r,1", Edn::Int(1))]
fn test_conformance_whitespace(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(read_string(input).unwrap(), expected);
}

#[test]
fn test_conformance_formfeed_not_whitespace() {
    // Form feed (0x0C) is NOT whitespace -- it should not separate tokens.
    let result = read_string("\x0C1");
    assert!(result.is_err());
}

#[test]
fn test_conformance_comma_in_collections() {
    let result = read_string("[1, 2, 3]").unwrap();
    assert_eq!(result, Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]));
}

// -- Section 2: Comments ---------------------------------------------------

#[rstest]
#[case("; comment\n1", Edn::Int(1))]
#[case("1 ; trailing", Edn::Int(1))]
#[case("; full line comment\n; another\n1", Edn::Int(1))]
fn test_conformance_comments(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(read_string(input).unwrap(), expected);
}

#[test]
fn test_conformance_comment_at_eof() {
    let result = read_string("1 ; comment");
    assert_eq!(result.unwrap(), Edn::Int(1));
}

#[test]
fn test_conformance_comment_inside_vector() {
    let result = read_string("[1 ; comment\n 2]").unwrap();
    assert_eq!(result, Edn::Vector(vec![Edn::Int(1), Edn::Int(2)]));
}

// -- Section 3: Discard ----------------------------------------------------

#[test]
fn test_conformance_discard_nested() {
    // #_ #_ a b: both discarded, yields c
    let result = read_string_with("[#_ #_ 1 2 3]", &clojure_config()).unwrap();
    assert_eq!(result, Edn::Vector(vec![Edn::Int(3)]));
}

#[test]
fn test_conformance_discard_tagged_literal() {
    // #_ #inst "2024" discards the entire tagged form.
    let result = read_string_with("[#_ #inst \"2024\" 1]", &clojure_config()).unwrap();
    assert_eq!(result, Edn::Vector(vec![Edn::Int(1)]));
}

#[test]
fn test_conformance_discard_in_map() {
    let result = read_string_with("{:a #_ :skip 1}", &clojure_config()).unwrap();
    let mut expected = HashMap::new();
    expected.insert(Edn::keyword("a"), Edn::Int(1));
    assert_eq!(result, Edn::Map(expected));
}

#[test]
fn test_conformance_discard_edn_dialect() {
    // In Edn mode, #_ parses as tagged literal with tag "_".
    let result = read_string_with("#_ 1", &edn_config()).unwrap();
    assert_eq!(result, Edn::Tagged("_".into(), Box::new(Edn::Int(1))));
}

// -- Section 4: Nil and booleans -------------------------------------------

#[rstest]
#[case("nil", Edn::Nil)]
#[case("true", Edn::Bool(true))]
#[case("false", Edn::Bool(false))]
fn test_conformance_literals(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(read_string(input).unwrap(), expected);
}

#[rstest]
#[case("Nil")]
#[case("TRUE")]
#[case("False")]
fn test_conformance_case_sensitive_literals(#[case] input: &str) {
    // Case variants are symbols, not literals.
    let result = read_string(input).unwrap();
    assert!(matches!(result, Edn::Symbol(_)));
}

// -- Section 5: Integers ---------------------------------------------------

#[rstest]
#[case("0", 0)]
#[case("1", 1)]
#[case("-1", -1)]
#[case("+1", 1)]
#[case("007", 7)]
#[case("0000", 0)]
#[case("9223372036854775807", i64::MAX)]
#[case("-9223372036854775808", i64::MIN)]
fn test_conformance_integers(#[case] input: &str, #[case] expected: i64) {
    assert_eq!(read_string(input).unwrap(), Edn::Int(expected));
}

#[test]
fn test_conformance_integer_overflow() {
    let result = read_string("9223372036854775808");
    assert!(result.is_err());
}

#[test]
fn test_conformance_negative_zero_int() {
    assert_eq!(read_string("-0").unwrap(), Edn::Int(0));
}

#[test]
fn test_conformance_leading_plus_both_dialects() {
    // Leading + is umol-edn extension, allowed in both dialects.
    assert_eq!(
        read_string_with("+5", &edn_config()).unwrap(),
        Edn::Int(5)
    );
    assert_eq!(
        read_string_with("+5", &clojure_config()).unwrap(),
        Edn::Int(5)
    );
}

#[rstest]
#[case("42N")]
#[case("0N")]
fn test_conformance_bigint_suffix_error(#[case] input: &str) {
    assert!(read_string(input).is_err());
}

// -- Section 6: Floating-point ---------------------------------------------

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
fn test_conformance_floats(#[case] input: &str, #[case] expected: f64) {
    match read_string(input).unwrap() {
        Edn::Float(v) => assert!((v - expected).abs() < f64::EPSILON),
        other => panic!("expected Float, got {other:?}"),
    }
}

#[test]
fn test_conformance_negative_zero_float() {
    match read_string("-0.0").unwrap() {
        Edn::Float(v) => {
            assert!(v == 0.0);
            assert!(v.is_sign_negative());
        }
        other => panic!("expected Float, got {other:?}"),
    }
}

#[test]
fn test_conformance_special_floats_clojure() {
    match read_string_with("##NaN", &clojure_config()).unwrap() {
        Edn::Float(v) => assert!(v.is_nan()),
        other => panic!("expected NaN, got {other:?}"),
    }
    match read_string_with("##Inf", &clojure_config()).unwrap() {
        Edn::Float(v) => assert!(v == f64::INFINITY),
        other => panic!("expected Inf, got {other:?}"),
    }
    match read_string_with("##-Inf", &clojure_config()).unwrap() {
        Edn::Float(v) => assert!(v == f64::NEG_INFINITY),
        other => panic!("expected -Inf, got {other:?}"),
    }
}

#[rstest]
#[case("##NaN")]
#[case("##Inf")]
#[case("##-Inf")]
fn test_conformance_special_floats_edn_rejected(#[case] input: &str) {
    assert!(read_string_with(input, &edn_config()).is_err());
}

#[rstest]
#[case("3.14M")]
fn test_conformance_bigdec_suffix_error(#[case] input: &str) {
    assert!(read_string(input).is_err());
}

// -- Section 7: Strings ----------------------------------------------------

#[rstest]
#[case(r#""""#, "")]
#[case(r#""hello""#, "hello")]
#[case(r#""line\nbreak""#, "line\nbreak")]
#[case(r#""tab\there""#, "tab\there")]
#[case(r#""cr\rhere""#, "cr\rhere")]
#[case(r#""slash\\here""#, "slash\\here")]
#[case(r#""quote\"here""#, "quote\"here")]
#[case(r#""\u0041""#, "A")]
#[case(r#""\u03B1""#, "\u{03B1}")]
fn test_conformance_string_escapes(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(read_string(input).unwrap(), Edn::Str(Cow::Owned(expected.to_string())));
}

#[test]
fn test_conformance_string_backspace_formfeed_clojure() {
    assert_eq!(
        read_string_with(r#""\b""#, &clojure_config()).unwrap(),
        Edn::Str(Cow::Owned("\u{0008}".to_string()))
    );
    assert_eq!(
        read_string_with(r#""\f""#, &clojure_config()).unwrap(),
        Edn::Str(Cow::Owned("\u{000C}".to_string()))
    );
}

#[rstest]
#[case(r#""\b""#)]
#[case(r#""\f""#)]
fn test_conformance_string_bf_edn_rejected(#[case] input: &str) {
    assert!(read_string_with(input, &edn_config()).is_err());
}

#[test]
fn test_conformance_string_octal_clojure() {
    // \0 = NUL, \101 = 'A' (65), \377 = 0xFF
    assert_eq!(
        read_string_with(r#""\0""#, &clojure_config()).unwrap(),
        Edn::Str(Cow::Owned("\0".to_string()))
    );
    assert_eq!(
        read_string_with(r#""\101""#, &clojure_config()).unwrap(),
        Edn::Str(Cow::Owned("A".to_string()))
    );
    assert_eq!(
        read_string_with(r#""\377""#, &clojure_config()).unwrap(),
        Edn::Str(Cow::Owned("\u{00FF}".to_string()))
    );
}

#[test]
fn test_conformance_string_octal_edn_rejected() {
    assert!(read_string_with(r#""\0""#, &edn_config()).is_err());
    assert!(read_string_with(r#""\101""#, &edn_config()).is_err());
}

#[test]
fn test_conformance_string_unterminated() {
    assert!(read_string(r#""hello"#).is_err());
}

#[rstest]
#[case(r#""\x""#)]
#[case(r#""\a""#)]
#[case(r#""\v""#)]
fn test_conformance_string_invalid_escape(#[case] input: &str) {
    assert!(read_string(input).is_err());
}

// -- Section 8: Characters -------------------------------------------------

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
fn test_conformance_characters(#[case] input: &str, #[case] expected: char) {
    assert_eq!(read_string(input).unwrap(), Edn::Char(expected));
}

#[test]
fn test_conformance_char_formfeed_backspace_clojure() {
    assert_eq!(
        read_string_with(r"\formfeed", &clojure_config()).unwrap(),
        Edn::Char('\u{000C}')
    );
    assert_eq!(
        read_string_with(r"\backspace", &clojure_config()).unwrap(),
        Edn::Char('\u{0008}')
    );
}

#[rstest]
#[case(r"\formfeed")]
#[case(r"\backspace")]
fn test_conformance_char_fb_edn_rejected(#[case] input: &str) {
    // In Edn strict, these are multi-char errors (not valid named chars).
    assert!(read_string_with(input, &edn_config()).is_err());
}

#[test]
fn test_conformance_char_backslash_space_invalid() {
    assert!(read_string(r"\ ").is_err());
}

#[test]
fn test_conformance_char_multichar_error() {
    // \abc is not a valid named character.
    assert!(read_string(r"\abc").is_err());
}

#[test]
fn test_conformance_char_termination() {
    // Character must be followed by non-symbol char, whitespace, or EOF.
    let result = read_string(r"[\a 1]").unwrap();
    assert_eq!(
        result,
        Edn::Vector(vec![Edn::Char('a'), Edn::Int(1)])
    );
}

// -- Section 9: Keywords ---------------------------------------------------

#[rstest]
#[case(":foo", "foo")]
#[case(":ns/name", "ns/name")]
#[case(":a.b/c", "a.b/c")]
fn test_conformance_keywords(#[case] input: &str, #[case] name: &str) {
    assert_eq!(read_string(input).unwrap(), Edn::keyword(name));
}

#[test]
fn test_conformance_bare_colon_error() {
    // Bare `:` with no name is an error.
    assert!(read_string(": ").is_err());
}

#[rstest]
#[case(":0", "0")]
#[case(":1", "1")]
#[case(":123abc", "123abc")]
fn test_conformance_digit_keywords_clojure(#[case] input: &str, #[case] name: &str) {
    assert_eq!(
        read_string_with(input, &clojure_config()).unwrap(),
        Edn::keyword(name)
    );
}

#[rstest]
#[case(":0")]
#[case(":1")]
#[case(":123abc")]
fn test_conformance_digit_keywords_edn_rejected(#[case] input: &str) {
    assert!(read_string_with(input, &edn_config()).is_err());
}

// -- Section 10: Symbols ---------------------------------------------------

#[rstest]
#[case("foo", "foo")]
#[case("ns/name", "ns/name")]
#[case("my.ns/sym", "my.ns/sym")]
fn test_conformance_symbols(#[case] input: &str, #[case] name: &str) {
    assert_eq!(read_string(input).unwrap(), Edn::Symbol(Symbol::new(name)));
}

#[test]
fn test_conformance_slash_symbol() {
    // `/` alone is a valid symbol.
    assert_eq!(read_string("/").unwrap(), Edn::Symbol(Symbol::new("/")));
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
fn test_conformance_symbol_start_chars(#[case] input: &str, #[case] name: &str) {
    assert_eq!(read_string(input).unwrap(), Edn::Symbol(Symbol::new(name)));
}

#[test]
fn test_conformance_symbol_with_all_char_types() {
    // Symbol using digits, +, -, #, :, ' in non-start positions.
    let result = read_string("foo+bar-baz#qux:quux'end").unwrap();
    assert_eq!(
        result,
        Edn::Symbol(Symbol::new("foo+bar-baz#qux:quux'end"))
    );
}

// -- Section 11-12: Lists and vectors --------------------------------------

#[rstest]
#[case("()", Edn::List(vec![]))]
#[case("(1 2 3)", Edn::List(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]))]
#[case("[]", Edn::Vector(vec![]))]
#[case("[1 2 3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]))]
fn test_conformance_seqs(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(read_string(input).unwrap(), expected);
}

#[test]
fn test_conformance_nested_collections() {
    let result = read_string("[[1 2] (3 4)]").unwrap();
    assert_eq!(
        result,
        Edn::Vector(vec![
            Edn::Vector(vec![Edn::Int(1), Edn::Int(2)]),
            Edn::List(vec![Edn::Int(3), Edn::Int(4)]),
        ])
    );
}

// -- Section 13: Maps ------------------------------------------------------

#[test]
fn test_conformance_map_basic() {
    let mut expected = HashMap::new();
    expected.insert(Edn::keyword("a"), Edn::Int(1));
    expected.insert(Edn::keyword("b"), Edn::Int(2));
    assert_eq!(read_string("{:a 1 :b 2}").unwrap(), Edn::Map(expected));
}

#[test]
fn test_conformance_map_odd_elements_error() {
    assert!(read_string("{:a 1 :b}").is_err());
}

#[test]
fn test_conformance_map_duplicate_key_error() {
    assert!(read_string("{:a 1 :a 2}").is_err());
}

#[test]
fn test_conformance_map_duplicate_key_last_wins() {
    let config = ParseConfig {
        duplicate_keys: DuplicateKeyPolicy::LastWins,
        ..Default::default()
    };
    let result = read_string_with("{:a 1 :a 2}", &config).unwrap();
    let mut expected = HashMap::new();
    expected.insert(Edn::keyword("a"), Edn::Int(2));
    assert_eq!(result, Edn::Map(expected));
}

#[test]
fn test_conformance_map_complex_keys() {
    // Maps accept arbitrary EDN values as keys.
    let result = read_string("{[1 2] :pair \"key\" :str}").unwrap();
    let mut expected = HashMap::new();
    expected.insert(
        Edn::Vector(vec![Edn::Int(1), Edn::Int(2)]),
        Edn::keyword("pair"),
    );
    expected.insert(
        Edn::Str(Cow::Owned("key".to_string())),
        Edn::keyword("str"),
    );
    assert_eq!(result, Edn::Map(expected));
}

#[test]
fn test_conformance_map_contains_all_entries() {
    let result = read_string("{:z 1 :a 2}").unwrap();
    if let Edn::Map(m) = result {
        assert_eq!(m.len(), 2);
        assert_eq!(m.get(&Edn::keyword("z")), Some(&Edn::Int(1)));
        assert_eq!(m.get(&Edn::keyword("a")), Some(&Edn::Int(2)));
    } else {
        panic!("expected Map");
    }
}

// -- Section 14: Sets ------------------------------------------------------

#[test]
fn test_conformance_set_basic() {
    let result = read_string("#{1 2 3}").unwrap();
    let expected: HashSet<Edn<'_>> = [Edn::Int(1), Edn::Int(2), Edn::Int(3)]
        .into_iter()
        .collect();
    assert_eq!(result, Edn::Set(expected));
}

#[test]
fn test_conformance_set_empty() {
    assert_eq!(read_string("#{}").unwrap(), Edn::Set(HashSet::new()));
}

#[test]
fn test_conformance_set_contains_all_elements() {
    let result = read_string("#{3 1 2}").unwrap();
    if let Edn::Set(s) = result {
        assert_eq!(s.len(), 3);
        assert!(s.contains(&Edn::Int(1)));
        assert!(s.contains(&Edn::Int(2)));
        assert!(s.contains(&Edn::Int(3)));
    } else {
        panic!("expected Set");
    }
}

// -- Section 15: Tagged literals -------------------------------------------

#[test]
fn test_conformance_tagged_basic() {
    let result = read_string("#myapp/Person {:name \"Alice\"}").unwrap();
    let mut inner = HashMap::new();
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
fn test_conformance_tagged_unqualified() {
    let result = read_string("#inst \"2024-01-01\"").unwrap();
    assert_eq!(
        result,
        Edn::Tagged("inst".into(), Box::new(Edn::Str(Cow::Owned("2024-01-01".to_string()))))
    );
}

// -- Round-trip ------------------------------------------------------------

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
#[case("#tag value")]
fn test_conformance_roundtrip(#[case] input: &str) {
    let parsed = read_string(input).unwrap();
    let rendered = parsed.to_string();
    let reparsed = read_string(&rendered).unwrap();
    assert_eq!(parsed, reparsed);
}

// -- Error reporting -------------------------------------------------------

#[test]
fn test_conformance_trailing_content_error() {
    let err = read_string("1 2").unwrap_err();
    assert!(matches!(err, EdnError::TrailingContent { .. }));
}

#[test]
fn test_conformance_empty_input_error() {
    let err = read_string("").unwrap_err();
    assert!(matches!(err, EdnError::UnexpectedEof { .. }));
}

#[test]
fn test_conformance_whitespace_only_error() {
    let err = read_string("   ").unwrap_err();
    assert!(matches!(err, EdnError::UnexpectedEof { .. }));
}

#[test]
fn test_conformance_unclosed_vector() {
    assert!(read_string("[1 2").is_err());
}

#[test]
fn test_conformance_unclosed_map() {
    assert!(read_string("{:a 1").is_err());
}

#[test]
fn test_conformance_unclosed_list() {
    assert!(read_string("(1 2").is_err());
}

#[test]
fn test_conformance_unclosed_set() {
    assert!(read_string("#{1 2").is_err());
}

// -- Deep nesting ----------------------------------------------------------

#[test]
fn test_conformance_deep_nesting() {
    let depth = 100;
    let input = "[".repeat(depth) + &"]".repeat(depth);
    let result = read_string(&input);
    assert!(result.is_ok());
}
