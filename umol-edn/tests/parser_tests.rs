use pretty_assertions::assert_eq;
use rstest::rstest;
use std::collections::{BTreeMap, BTreeSet};

use umol_edn::{
    read_all, read_string, read_string_with, DuplicateKeyPolicy, Edn, EdnError, Keyword,
    ParseConfig, Reader,
};

// --- Primitives ---

#[rstest]
#[case("nil", Edn::Nil)]
#[case("true", Edn::Bool(true))]
#[case("false", Edn::Bool(false))]
fn test_read_string_literals(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(read_string(input).unwrap(), expected);
}

#[rstest]
#[case("0", 0)]
#[case("42", 42)]
#[case("-1", -1)]
#[case("+5", 5)]
#[case("9223372036854775807", i64::MAX)]
#[case("-9223372036854775808", i64::MIN)]
fn test_read_string_int(#[case] input: &str, #[case] expected: i64) {
    assert_eq!(read_string(input).unwrap(), Edn::Int(expected));
}

#[rstest]
#[case("1.0", 1.0)]
#[case("-3.14", -3.14)]
#[case("1e10", 1e10)]
#[case("1.5e-3", 1.5e-3)]
#[case("1E10", 1e10)]
fn test_read_string_float(#[case] input: &str, #[case] expected: f64) {
    assert_eq!(read_string(input).unwrap(), Edn::Float(expected));
}

#[rstest]
#[case("##NaN")]
fn test_read_string_nan(#[case] input: &str) {
    match read_string(input).unwrap() {
        Edn::Float(f) => assert!(f.is_nan()),
        other => panic!("expected Float(NaN), got {other:?}"),
    }
}

#[rstest]
#[case("##Inf", f64::INFINITY)]
#[case("##-Inf", f64::NEG_INFINITY)]
fn test_read_string_special_float(#[case] input: &str, #[case] expected: f64) {
    assert_eq!(read_string(input).unwrap(), Edn::Float(expected));
}

// --- Strings ---

#[rstest]
#[case(r#""""#, "")]
#[case(r#""hello""#, "hello")]
#[case(r#""hello world""#, "hello world")]
#[case(r#""line\nbreak""#, "line\nbreak")]
#[case(r#""tab\there""#, "tab\there")]
#[case(r#""quote\"here""#, "quote\"here")]
#[case(r#""back\\slash""#, "back\\slash")]
#[case(r#""cr\rhere""#, "cr\rhere")]
fn test_read_string_str(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(
        read_string(input).unwrap(),
        Edn::Str(std::borrow::Cow::Owned(expected.to_string()))
    );
}

#[rstest]
#[case(r#""\u0041""#, "A")]
#[case(r#""\u03BB""#, "\u{03BB}")]
fn test_read_string_unicode_escape(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(
        read_string(input).unwrap(),
        Edn::Str(std::borrow::Cow::Owned(expected.to_string()))
    );
}

// --- Characters ---

#[rstest]
#[case("\\a", 'a')]
#[case("\\Z", 'Z')]
#[case("\\newline", '\n')]
#[case("\\return", '\r')]
#[case("\\space", ' ')]
#[case("\\tab", '\t')]
#[case("\\u0041", 'A')]
fn test_read_string_char(#[case] input: &str, #[case] expected: char) {
    assert_eq!(read_string(input).unwrap(), Edn::Char(expected));
}

// --- Keywords ---

#[rstest]
#[case(":foo", "foo")]
#[case(":ns/name", "ns/name")]
#[case(":a.b/c", "a.b/c")]
fn test_read_string_keyword(#[case] input: &str, #[case] expected_name: &str) {
    let val = read_string(input).unwrap();
    match &val {
        Edn::Keyword(k) => assert_eq!(k.as_str(), expected_name),
        other => panic!("expected Keyword, got {other:?}"),
    }
}

#[test]
fn test_keyword_namespace() {
    let k = Keyword::new("ns/name");
    assert_eq!(k.namespace(), Some("ns"));
    assert_eq!(k.name(), "name");

    let k2 = Keyword::new("bare");
    assert_eq!(k2.namespace(), None);
    assert_eq!(k2.name(), "bare");
}

// --- Symbols ---

#[rstest]
#[case("foo", "foo")]
#[case("ns/name", "ns/name")]
#[case("my.ns/sym", "my.ns/sym")]
fn test_read_string_symbol(#[case] input: &str, #[case] expected_name: &str) {
    let val = read_string(input).unwrap();
    match &val {
        Edn::Symbol(s) => assert_eq!(s.as_str(), expected_name),
        other => panic!("expected Symbol, got {other:?}"),
    }
}

// --- Collections ---

#[rstest]
#[case("()", Edn::List(vec![]))]
#[case("(1 2 3)", Edn::List(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]))]
#[case("(nil true)", Edn::List(vec![Edn::Nil, Edn::Bool(true)]))]
fn test_read_string_list(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(read_string(input).unwrap(), expected);
}

#[rstest]
#[case("[]", Edn::Vector(vec![]))]
#[case("[1 2 3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]))]
fn test_read_string_vector(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(read_string(input).unwrap(), expected);
}

#[test]
fn test_read_string_map() {
    let val = read_string("{:a 1 :b 2}").unwrap();
    let mut expected = BTreeMap::new();
    expected.insert(Edn::Keyword(Keyword::new("a")), Edn::Int(1));
    expected.insert(Edn::Keyword(Keyword::new("b")), Edn::Int(2));
    assert_eq!(val, Edn::Map(expected));
}

#[test]
fn test_read_string_empty_map() {
    assert_eq!(read_string("{}").unwrap(), Edn::Map(BTreeMap::new()));
}

#[test]
fn test_read_string_set() {
    let val = read_string("#{1 2 3}").unwrap();
    let mut expected = BTreeSet::new();
    expected.insert(Edn::Int(1));
    expected.insert(Edn::Int(2));
    expected.insert(Edn::Int(3));
    assert_eq!(val, Edn::Set(expected));
}

#[test]
fn test_read_string_empty_set() {
    assert_eq!(read_string("#{}").unwrap(), Edn::Set(BTreeSet::new()));
}

// --- Nested ---

#[test]
fn test_read_string_nested() {
    let val = read_string("{:items [1 (2 3)] :flag true}").unwrap();
    let mut expected = BTreeMap::new();
    expected.insert(
        Edn::Keyword(Keyword::new("items")),
        Edn::Vector(vec![
            Edn::Int(1),
            Edn::List(vec![Edn::Int(2), Edn::Int(3)]),
        ]),
    );
    expected.insert(Edn::Keyword(Keyword::new("flag")), Edn::Bool(true));
    assert_eq!(val, Edn::Map(expected));
}

// --- Tagged literals ---

#[test]
fn test_read_string_tagged() {
    let val = read_string("#myapp/Person {:name \"Alice\"}").unwrap();
    match val {
        Edn::Tagged(tag, inner) => {
            assert_eq!(tag, "myapp/Person");
            assert!(inner.is_map());
        }
        other => panic!("expected Tagged, got {other:?}"),
    }
}

// --- Comments ---

#[rstest]
#[case("; comment\n42", Edn::Int(42))]
#[case("42 ; trailing", Edn::Int(42))]
fn test_read_string_comment(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(read_string(input).unwrap(), expected);
}

// --- Discard ---

#[rstest]
#[case("#_ foo 42", Edn::Int(42))]
#[case("[1 #_ 2 3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(3)]))]
fn test_read_string_discard(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(read_string(input).unwrap(), expected);
}

// --- Whitespace variants ---

#[rstest]
#[case("  42  ", Edn::Int(42))]
#[case("\t42\n", Edn::Int(42))]
#[case(",42,", Edn::Int(42))]
#[case("[1,,2,,3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]))]
fn test_read_string_whitespace(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(read_string(input).unwrap(), expected);
}

// --- Error cases ---

#[rstest]
#[case("")]
#[case("   ")]
fn test_read_string_error_empty(#[case] input: &str) {
    assert!(read_string(input).is_err());
}

#[test]
fn test_read_string_error_trailing() {
    let err = read_string("42 43").unwrap_err();
    match err {
        EdnError::TrailingContent { .. } => {}
        other => panic!("expected TrailingContent, got {other:?}"),
    }
}

#[test]
fn test_read_string_error_duplicate_key() {
    let err = read_string("{:a 1 :a 2}").unwrap_err();
    assert!(
        matches!(err, EdnError::Custom(_)),
        "expected error for duplicate key, got {err:?}"
    );
}

#[test]
fn test_read_string_duplicate_key_last_wins() {
    let config = ParseConfig {
        duplicate_keys: DuplicateKeyPolicy::LastWins,
        ..Default::default()
    };
    let val = read_string_with("{:a 1 :a 2}", &config).unwrap();
    let mut expected = BTreeMap::new();
    expected.insert(Edn::Keyword(Keyword::new("a")), Edn::Int(2));
    assert_eq!(val, Edn::Map(expected));
}

// --- Integer overflow ---

#[test]
fn test_read_string_error_integer_overflow() {
    assert!(read_string("99999999999999999999").is_err());
}

// --- read_all ---

#[test]
fn test_read_all() {
    let values = read_all("1 2 3").unwrap();
    assert_eq!(values, vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)]);
}

#[test]
fn test_read_all_empty() {
    let values = read_all("").unwrap();
    assert!(values.is_empty());
}

#[test]
fn test_read_all_mixed() {
    let values = read_all(":a [1 2] nil").unwrap();
    assert_eq!(values.len(), 3);
    assert!(values[0].is_keyword());
    assert!(values[1].is_vector());
    assert!(values[2].is_nil());
}

// --- Reader iterator ---

#[test]
fn test_reader_iterator() {
    let reader = Reader::new("1 :foo [3]");
    let values: Result<Vec<_>, _> = reader.collect();
    let values = values.unwrap();
    assert_eq!(values.len(), 3);
    assert_eq!(values[0], Edn::Int(1));
    assert!(values[1].is_keyword());
    assert!(values[2].is_vector());
}

#[test]
fn test_reader_empty() {
    let reader = Reader::new("");
    let values: Vec<_> = reader.collect();
    assert!(values.is_empty());
}

// --- Round-trip ---

#[rstest]
#[case("nil")]
#[case("true")]
#[case("false")]
#[case("42")]
#[case("-1")]
#[case("1.5")]
#[case(":keyword")]
#[case(":ns/name")]
#[case("symbol")]
#[case("()")]
#[case("(1 2 3)")]
#[case("[]")]
#[case("[1 2 3]")]
#[case("{}")]
#[case("\\a")]
#[case("\\newline")]
#[case("\\space")]
fn test_roundtrip(#[case] input: &str) {
    let val = read_string(input).unwrap();
    let formatted = val.to_string();
    let reparsed = read_string(&formatted).unwrap();
    assert_eq!(val, reparsed);
}

#[test]
fn test_roundtrip_string() {
    let val = read_string(r#""hello\nworld""#).unwrap();
    let formatted = val.to_string();
    let reparsed = read_string(&formatted).unwrap();
    assert_eq!(val, reparsed);
}

#[test]
fn test_roundtrip_map() {
    let val = read_string("{:a 1, :b 2}").unwrap();
    let formatted = val.to_string();
    let reparsed = read_string(&formatted).unwrap();
    assert_eq!(val, reparsed);
}

#[test]
fn test_roundtrip_set() {
    let val = read_string("#{1 2 3}").unwrap();
    let formatted = val.to_string();
    let reparsed = read_string(&formatted).unwrap();
    assert_eq!(val, reparsed);
}

#[test]
fn test_roundtrip_special_floats() {
    for input in &["##Inf", "##-Inf"] {
        let val = read_string(input).unwrap();
        let formatted = val.to_string();
        let reparsed = read_string(&formatted).unwrap();
        assert_eq!(val, reparsed);
    }
}

// --- Edn accessors ---

#[test]
fn test_edn_get() {
    let val = read_string("{:name \"Alice\" :age 30}").unwrap();
    assert_eq!(val.get("name").unwrap().as_str(), Some("Alice"));
    assert_eq!(val.get("age").unwrap().as_i64(), Some(30));
    assert!(val.get("missing").is_none());
}

#[test]
fn test_edn_numeric_narrowing() {
    let val = Edn::Int(42);
    assert_eq!(val.as_u8(), Some(42));
    assert_eq!(val.as_u16(), Some(42));
    assert_eq!(val.as_u32(), Some(42));
    assert_eq!(val.as_i32(), Some(42));

    let val_neg = Edn::Int(-1);
    assert_eq!(val_neg.as_u8(), None);
    assert_eq!(val_neg.as_i8(), Some(-1));
}

#[test]
fn test_edn_iter() {
    let val = read_string("[1 2 3]").unwrap();
    let items: Vec<_> = val.iter().collect();
    assert_eq!(items.len(), 3);
}

// --- Clojure-compatible escapes (non-strict mode) ---

#[test]
fn test_read_string_backspace_formfeed() {
    let val = read_string(r#""\b\f""#).unwrap();
    assert_eq!(
        val,
        Edn::Str(std::borrow::Cow::Owned("\u{0008}\u{000C}".to_string()))
    );
}

#[test]
fn test_read_string_char_formfeed_backspace() {
    assert_eq!(read_string("\\formfeed").unwrap(), Edn::Char('\u{000C}'));
    assert_eq!(read_string("\\backspace").unwrap(), Edn::Char('\u{0008}'));
}
