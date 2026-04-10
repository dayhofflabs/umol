//! EDN conformance tests — parser and streaming paths.

use std::borrow::Cow;
use std::collections::HashMap;
#[cfg(feature = "bignum")]
use std::str::FromStr;
use std::sync::Arc;

use rstest::rstest;
use serde::Deserialize;
use umol_edn::serde::{from_str, from_str_with};
use umol_edn::{
    read_all_with, read_string_with, DuplicateKeyPolicy, Edn, EdnError, EdnMap, EdnSet, EdnSymbol,
    ParseConfig, ParseError, TagFn, TagReaders,
};

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
#[case::spaces(" 1 ", Edn::Int(1))]
#[case::tabs("\t1\t", Edn::Int(1))]
#[case::newlines("\n1\n", Edn::Int(1))]
#[case::carriage_returns("\r1\r", Edn::Int(1))]
#[case::commas(",1,", Edn::Int(1))]
#[case::mixed(" \t\n\r,1", Edn::Int(1))]
fn test_whitespace(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(parse(input).unwrap(), expected);
}

#[rstest]
#[case::spaces(" 1 ", 1i64)]
#[case::tabs("\t1\t", 1)]
#[case::newlines("\n1\n", 1)]
#[case::commas(",1,", 1)]
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
        Edn::Symbol(EdnSymbol::new("foo#bar"))
    );
    assert_eq!(
        parse("[a#b]").unwrap(),
        Edn::Vector(vec![Edn::Symbol(EdnSymbol::new("a#b"))].into())
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
#[case::true_val("true", Edn::Bool(true))]
#[case::false_val("false", Edn::Bool(false))]
fn test_booleans(#[case] input: &str, #[case] expected: Edn<'_>) {
    assert_eq!(parse(input).unwrap(), expected);
}

#[rstest]
#[case::true_val("true", true)]
#[case::false_val("false", false)]
fn test_booleans_streaming(#[case] input: &str, #[case] expected: bool) {
    assert_eq!(stream::<bool>(input).unwrap(), expected);
}

#[rstest]
#[case::nil_capitalized("Nil")]
#[case::true_upper("TRUE")]
#[case::false_capitalized("False")]
fn test_case_sensitive(#[case] input: &str) {
    assert!(matches!(parse(input).unwrap(), Edn::Symbol(_)));
}

#[rstest]
#[case::empty(r#""""#, "")]
#[case::simple(r#""hello""#, "hello")]
#[case::newline(r#""line\nbreak""#, "line\nbreak")]
#[case::tab(r#""tab\there""#, "tab\there")]
#[case::carriage_return(r#""cr\rhere""#, "cr\rhere")]
#[case::backslash(r#""slash\\here""#, "slash\\here")]
#[case::quote(r#""quote\"here""#, "quote\"here")]
fn test_string_escapes(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(
        parse(input).unwrap(),
        Edn::Str(Cow::Owned(expected.to_string()))
    );
}

#[rstest]
#[case::simple(r#""hello""#, "hello")]
#[case::newline(r#""line\nbreak""#, "line\nbreak")]
#[case::tab(r#""tab\there""#, "tab\there")]
fn test_string_escapes_streaming(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(stream::<String>(input).unwrap(), expected);
}

#[rstest]
#[case::backspace(r#""\b""#)]
#[case::formfeed(r#""\f""#)]
#[case::octal(r#""\101""#)]
fn test_string_clojure_escapes_rejected(#[case] input: &str) {
    assert!(parse(input).is_err(), "Edn should reject {input}");
}

#[rstest]
#[case::ascii(r#""\u0041""#, "A")]
#[case::accent(r#""\u00e9""#, "é")]
fn test_string_unicode_escape(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(
        parse(input).unwrap(),
        Edn::Str(Cow::Owned(expected.to_string()))
    );
    assert_eq!(stream::<String>(input).unwrap(), expected);
}

#[rstest]
#[case::backspace(r#""\b""#)]
#[case::formfeed(r#""\f""#)]
fn test_string_clojure_escapes_rejected_streaming(#[case] input: &str) {
    assert!(
        stream::<String>(input).is_err(),
        "Edn streaming should reject {input}"
    );
}

#[rstest]
#[case::backslash_x(r#""\x""#)]
#[case::backslash_a(r#""\a""#)]
#[case::backslash_v(r#""\v""#)]
fn test_string_invalid_escape(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[rstest]
#[case::high_surrogate(r#""\uD800""#)]
#[case::low_surrogate(r#""\uDFFF""#)]
#[case::invalid_hex(r#""\u00GG""#)]
#[case::too_short(r#""\u41""#)]
fn test_string_unicode_escape_error(#[case] input: &str) {
    assert!(parse(input).is_err(), "Edn should reject {input}");
}

#[test]
fn test_string_unterminated() {
    assert!(parse(r#""hello"#).is_err());
}

#[rstest]
#[case::lowercase(r"\a", 'a')]
#[case::uppercase(r"\Z", 'Z')]
#[case::punctuation(r"\!", '!')]
#[case::newline(r"\newline", '\n')]
#[case::return_char(r"\return", '\r')]
#[case::space(r"\space", ' ')]
#[case::tab(r"\tab", '\t')]
#[case::unicode_ascii(r"\u0041", 'A')]
#[case::unicode_alpha(r"\u03B1", '\u{03B1}')]
fn test_characters(#[case] input: &str, #[case] expected: char) {
    assert_eq!(parse(input).unwrap(), Edn::Char(expected));
}

#[rstest]
#[case::formfeed(r"\formfeed")]
#[case::backspace(r"\backspace")]
fn test_formfeed_backspace_rejected(#[case] input: &str) {
    assert!(parse(input).is_err(), "Edn should reject {input}");
}

#[rstest]
#[case::backslash_space(r"\ ")]
#[case::multi_char(r"\abc")]
#[case::invalid_hex(r"\u00GG")]
#[case::surrogate(r"\uD800")]
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
#[case::bare("foo", "foo")]
#[case::namespaced("ns/name", "ns/name")]
#[case::dotted_ns("my.ns/sym", "my.ns/sym")]
fn test_symbols(#[case] input: &str, #[case] name: &str) {
    assert_eq!(parse(input).unwrap(), Edn::Symbol(EdnSymbol::new(name)));
}

#[test]
fn test_slash_alone() {
    assert_eq!(parse("/").unwrap(), Edn::Symbol(EdnSymbol::new("/")));
}

#[rstest]
#[case::dot(".", ".")]
#[case::star("*", "*")]
#[case::bang("!", "!")]
#[case::underscore("_", "_")]
#[case::question("?", "?")]
#[case::dollar("$", "$")]
#[case::percent("%", "%")]
#[case::ampersand("&", "&")]
#[case::equals("=", "=")]
#[case::less_than("<", "<")]
#[case::greater_than(">", ">")]
fn test_start_chars(#[case] input: &str, #[case] name: &str) {
    assert_eq!(parse(input).unwrap(), Edn::Symbol(EdnSymbol::new(name)));
}

#[test]
fn test_interior_chars() {
    assert_eq!(
        parse("foo+bar-baz#qux:quux'end").unwrap(),
        Edn::Symbol(EdnSymbol::new("foo+bar-baz#qux:quux'end"))
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
#[case::empty_name("ns/")]
#[case::empty_prefix("/name")]
fn test_empty_prefix_or_name(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[rstest]
#[case::double_slash("a/b/c")]
#[case::digit_after_slash("foo/1bar")]
#[case::hash_after_slash("foo/#bar")]
#[case::colon_after_slash("foo/:bar")]
fn test_symbol_slash_rejected(#[case] input: &str) {
    assert!(parse(input).is_err(), "Edn should reject {input}");
}

#[rstest]
#[case::underscore("foo/_bar")]
#[case::simple("foo/bar")]
fn test_post_slash_valid(#[case] input: &str) {
    assert!(parse(input).is_ok());
}

#[rstest]
#[case::simple(":foo", "foo")]
#[case::single_char(":a", "a")]
#[case::namespaced(":ns/name", "ns/name")]
#[case::dotted_ns(":my.ns/foo", "my.ns/foo")]
fn test_keywords(#[case] input: &str, #[case] name: &str) {
    assert_eq!(parse(input).unwrap(), Edn::keyword(name));
}

#[rstest]
#[case::bare(":foo", "foo")]
#[case::namespaced(":ns/name", "ns/name")]
fn test_keywords_streaming(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(stream::<String>(input).unwrap(), expected);
}

#[rstest]
#[case::single_digit(":0")]
#[case::digit_prefix(":0foo")]
fn test_digit_start_rejected(#[case] input: &str) {
    assert!(parse(input).is_err(), "Edn should reject {input}");
}

#[test]
fn test_hash_start_rejected() {
    assert!(parse(":#foo").is_err());
}

#[rstest]
#[case::dot(":.foo")]
#[case::plus(":+foo")]
#[case::minus(":-foo")]
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
#[case::bare_slash(":/")]
#[case::slash_prefix(":/foo")]
fn test_invalid_slash(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[test]
fn test_bare_colon_error() {
    assert!(parse(": ").is_err());
}

#[rstest]
#[case::zero("0", 0)]
#[case::positive("1", 1)]
#[case::negative("-1", -1)]
#[case::plus_sign("+1", 1)]
#[case::max("9223372036854775807", i64::MAX)]
#[case::min("-9223372036854775808", i64::MIN)]
fn test_integers(#[case] input: &str, #[case] expected: i64) {
    assert_eq!(parse(input).unwrap(), Edn::Int(expected));
}

#[rstest]
#[case::zero("0", 0i64)]
#[case::positive("1", 1)]
#[case::negative("-1", -1)]
#[case::plus_sign("+5", 5)]
fn test_integers_streaming(#[case] input: &str, #[case] expected: i64) {
    assert_eq!(stream::<i64>(input).unwrap(), expected);
}

#[rstest]
#[case::triple_zero("007")]
#[case::double_zero("00")]
fn test_leading_zeros_rejected(#[case] input: &str) {
    assert!(parse(input).is_err(), "Edn should reject {input}");
}

#[test]
fn test_negative_zero() {
    assert_eq!(parse("-0").unwrap(), Edn::Int(0));
}

#[test]
fn test_positive_zero() {
    assert_eq!(parse("+0").unwrap(), Edn::Int(0));
}

#[cfg(not(feature = "bignum"))]
#[test]
fn test_overflow() {
    assert!(parse("9223372036854775808").is_err());
}

#[cfg(not(feature = "bignum"))]
#[rstest]
#[case::positive("42N")]
#[case::zero("0N")]
#[case::float("3.14N")]
fn test_bigint_suffix_error(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[rstest]
#[case::one("1.0", 1.0)]
#[case::negative("-2.5", -2.5)]
#[case::lower_e("1e10", 1e10)]
#[case::upper_e("1E10", 1e10)]
#[case::neg_exponent("1e-10", 1e-10)]
#[case::pos_exponent("1e+10", 1e10)]
#[case::fractional_exp("1.5e3", 1.5e3)]
#[case::fractional_neg_exp("1.5E-3", 1.5e-3)]
#[case::plus_sign("+1.0", 1.0)]
fn test_floats(#[case] input: &str, #[case] expected: f64) {
    match parse(input).unwrap() {
        Edn::Float(v) => assert!((v - expected).abs() < f64::EPSILON),
        other => panic!("expected Float, got {other:?}"),
    }
}

#[rstest]
#[case::one("1.0", 1.0f64)]
#[case::negative("-2.5", -2.5)]
#[case::exponent("1e10", 1e10)]
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
#[case::nan("##NaN")]
#[case::inf("##Inf")]
#[case::neg_inf("##-Inf")]
#[case::huge_exponent("8E1313")]
#[case::overflow("1e999")]
#[case::neg_overflow("-1e999")]
fn test_special_floats_rejected(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[rstest]
#[case::nan("##NaN")]
#[case::inf("##Inf")]
fn test_special_floats_rejected_streaming(#[case] input: &str) {
    assert!(stream::<f64>(input).is_err());
}

#[cfg(not(feature = "bignum"))]
#[rstest]
#[case::decimal("3.14M")]
fn test_bigdec_suffix_error(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[cfg(feature = "bignum")]
#[rstest]
#[case::positive("42N", num_bigint::BigInt::from(42))]
#[case::zero("0N", num_bigint::BigInt::from(0))]
#[case::negative("-1N", num_bigint::BigInt::from(-1))]
#[case::plus_sign("+7N", num_bigint::BigInt::from(7))]
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
#[case::pi("3.14M")]
#[case::zero("0.0M")]
#[case::negative("-1.5M")]
#[case::integer("42M")]
fn test_bigdecimal_parse(#[case] input: &str) {
    assert!(matches!(parse(input).unwrap(), Edn::BigDecimal(_)));
}

#[cfg(feature = "bignum")]
#[rstest]
#[case::positive_overflow("9223372036854775808")] // i64::MAX + 1
#[case::negative_overflow("-9223372036854775809")] // i64::MIN - 1
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
#[case::bigint("42N")]
#[case::bigdecimal("3.14M")]
#[case::large_bigint("99999999999999999999N")]
#[case::small_bigdecimal("-0.001M")]
fn test_bignum_roundtrip(#[case] input: &str) {
    let parsed = parse(input).unwrap();
    let rendered = parsed.to_string();
    let reparsed = parse(&rendered).unwrap();
    assert_eq!(parsed, reparsed);
}

#[rstest]
#[case::empty_list("()", Edn::List(vec![].into()))]
#[case::list("(1 2 3)", Edn::List(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()))]
#[case::empty_vector("[]", Edn::Vector(vec![].into()))]
#[case::vector("[1 2 3]", Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()))]
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
fn test_set_duplicate_elements() {
    let result = parse("#{1 1 2}").unwrap();
    if let Edn::Set(s) = result {
        assert_eq!(s.len(), 2);
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
#[case::scalar("#foo 1")]
#[case::collection("#bar [1 2]")]
fn test_bare_tag_rejected(#[case] input: &str) {
    assert!(
        parse(input).is_err(),
        "Edn should reject bare tag in: {input}"
    );
}

#[test]
fn test_bare_tag_rejected_by_serde_default() {
    assert!(matches!(
        stream::<i64>("#foo 1"),
        Err(EdnError::Parse(ParseError::InvalidTag { .. }))
    ));
}

#[test]
fn test_bare_tag_passes_through_serde_with_permissive_config() {
    let config = ParseConfig {
        allow_unknown_tags: true,
        ..Default::default()
    };
    assert_eq!(from_str_with::<i64>("#foo 1", &config).unwrap(), 1);
}

#[test]
fn test_tag_without_value_error() {
    assert!(parse("#tag").is_err());
    assert!(parse("#myapp/Person").is_err());
}

#[rstest]
#[case::inst(
    "#inst \"2024-01-01T00:00:00Z\"",
    "inst",
    Edn::Str(Cow::Borrowed("2024-01-01T00:00:00Z"))
)]
#[case::uuid(
    "#uuid \"f81d4fae-7dec-11d0-a765-00a0c91e6bf6\"",
    "uuid",
    Edn::Str(Cow::Borrowed("f81d4fae-7dec-11d0-a765-00a0c91e6bf6"))
)]
fn test_builtin_tags_parse_without_features(
    #[case] input: &str,
    #[case] tag: &str,
    #[case] inner: Edn<'_>,
) {
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
#[case::not_a_date("#inst \"not-a-date\"")]
#[case::no_time("#inst \"2024-01-01\"")]
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
    use umol_edn::{inst_to_edn, read_string};
    let dt = chrono::DateTime::parse_from_rfc3339("2024-06-15T08:30:00+00:00").unwrap();
    let edn_val = inst_to_edn(&dt);
    let printed = edn_val.to_string();
    let parsed = read_string(&printed).unwrap();
    assert_eq!(parsed, edn_val);
}

#[cfg(feature = "uuid")]
#[test]
fn test_uuid_to_edn_roundtrip() {
    use umol_edn::{read_string, uuid_to_edn};
    let id = uuid::Uuid::parse_str("f81d4fae-7dec-11d0-a765-00a0c91e6bf6").unwrap();
    let edn_val = uuid_to_edn(&id);
    let printed = edn_val.to_string();
    let parsed = read_string(&printed).unwrap();
    assert_eq!(parsed, edn_val);
}

#[test]
fn test_custom_reader_dispatch() {
    fn double_reader(val: Edn) -> Result<Edn, ParseError> {
        match val {
            Edn::Int(n) => Ok(Edn::Int(n * 2)),
            _ => Err(ParseError::InvalidTag {
                offset: 0,
                tag: "double".into(),
            }),
        }
    }
    let mut readers = TagReaders::default();
    readers.insert("double", Arc::new(double_reader));
    let config = ParseConfig {
        tag_readers: readers,
        ..Default::default()
    };
    let val = read_string_with("#double 5", &config).unwrap();
    assert_eq!(val, Edn::Int(10));
}

#[test]
fn test_custom_reader_error_propagation() {
    fn strict_reader(_val: Edn) -> Result<Edn, ParseError> {
        Err(ParseError::InvalidTag {
            offset: 0,
            tag: "strict".into(),
        })
    }
    let mut readers = TagReaders::default();
    readers.insert("strict", Arc::new(strict_reader));
    let config = ParseConfig {
        tag_readers: readers,
        ..Default::default()
    };
    assert!(read_string_with("#strict 1", &config).is_err());
}

#[test]
fn test_tag_reader_captures_registry() {
    let registry: Arc<HashMap<&'static str, i64>> =
        Arc::new([("one", 1), ("two", 2)].into_iter().collect());
    let captured = registry.clone();
    let reader: TagFn = Arc::new(move |edn| match edn {
        Edn::Str(s) => Ok(Edn::Int(captured.get(s.as_ref()).copied().unwrap_or(-1))),
        _ => Err(ParseError::InvalidTag {
            offset: 0,
            tag: "lookup".into(),
        }),
    });
    let mut readers = TagReaders::default();
    readers.insert("lookup", reader);
    let config = ParseConfig {
        tag_readers: readers,
        ..Default::default()
    };
    let parsed = read_string_with(r#"#lookup "two""#, &config).unwrap();
    assert_eq!(parsed, Edn::Int(2));
}

#[rstest]
#[case::leading("; comment\n1", Edn::Int(1))]
#[case::trailing("1 ; trailing", Edn::Int(1))]
#[case::multiple("; full line comment\n; another\n1", Edn::Int(1))]
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
#[case::skip_middle("[1 #_ 2 3]", vec![1, 3])]
#[case::nested_discard("[#_ #_ 1 2 3]", vec![3])]
#[case::discard_tagged("[#_ #my/tag \"2024\" 1]", vec![1])]
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
#[case::skip_middle("[1 #_ 2 3]", vec![1, 3])]
#[case::nested_discard("[#_ #_ 1 2 3]", vec![3])]
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

#[test]
fn test_discard_map_key() {
    let result = parse("{#_ :skip :a 1}").unwrap();
    let mut expected = EdnMap::new();
    expected.insert(Edn::keyword("a"), Edn::Int(1));
    assert_eq!(result, Edn::Map(expected));
}

#[rstest]
#[case::bare_discard("#_")]
#[case::discard_only("#_ 1")]
fn test_discard_error(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[test]
fn test_int_not_equal_to_float() {
    assert_ne!(Edn::Int(1), Edn::Float(1.0));
}

#[rstest]
#[case::empty(vec![])]
#[case::three_ints(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)])]
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

#[test]
fn test_tagged_equality() {
    let a = Edn::Tagged("my/tag".into(), Box::new(Edn::Int(1)));
    let b = Edn::Tagged("my/tag".into(), Box::new(Edn::Int(1)));
    assert_eq!(a, b);
}

#[test]
fn test_tagged_inequality_tag() {
    let a = Edn::Tagged("my/tag".into(), Box::new(Edn::Int(1)));
    let b = Edn::Tagged("my/other".into(), Box::new(Edn::Int(1)));
    assert_ne!(a, b);
}

#[test]
fn test_tagged_inequality_value() {
    let a = Edn::Tagged("my/tag".into(), Box::new(Edn::Int(1)));
    let b = Edn::Tagged("my/tag".into(), Box::new(Edn::Int(2)));
    assert_ne!(a, b);
}

#[rstest]
#[case::nil("nil")]
#[case::true_val("true")]
#[case::false_val("false")]
#[case::zero("0")]
#[case::negative_int("-1")]
#[case::max_int("9223372036854775807")]
#[case::float("1.0")]
#[case::negative_float("-3.14")]
#[case::scientific("1.5e3")]
#[case::string(r#""hello""#)]
#[case::string_escape(r#""line\nbreak""#)]
#[case::char(r"\a")]
#[case::newline(r"\newline")]
#[case::space(r"\space")]
#[case::keyword(":foo")]
#[case::namespaced_keyword(":ns/name")]
#[case::symbol("my-sym")]
#[case::slash("/")]
#[case::empty_list("()")]
#[case::list("(1 2 3)")]
#[case::empty_vector("[]")]
#[case::vector("[1 2 3]")]
#[case::empty_map("{}")]
#[case::map("{:a 1}")]
#[case::empty_set("#{}")]
#[case::set("#{1 2 3}")]
#[case::tagged("#my/tag value")]
#[case::accent("\"é\"")]
#[case::greek("\"α\"")]
#[case::cjk("\"世界\"")]
fn test_roundtrip(#[case] input: &str) {
    let parsed = parse(input).unwrap();
    let rendered = parsed.to_string();
    let reparsed = parse(&rendered).unwrap();
    assert_eq!(parsed, reparsed);
}

#[rstest]
#[case::set_with_nil(Edn::Set(EdnSet::from_iter([Edn::Nil, Edn::Str(Cow::Owned("wbgpo?j".into()))])))]
#[case::list_with_accent(Edn::List(vec![Edn::Str(Cow::Owned("é".into()))].into()))]
fn test_roundtrip_regression(#[case] edn: Edn<'static>) {
    let rendered = edn.to_string();
    let reparsed = parse(&rendered).unwrap();
    assert_eq!(edn, reparsed);
}

#[test]
fn test_error_trailing_content() {
    let err = parse("1 2").unwrap_err();
    assert!(matches!(
        err,
        EdnError::Parse(ParseError::TrailingContent { .. })
    ));
}

#[test]
fn test_error_trailing_content_streaming() {
    let err = stream::<i64>("1 2").unwrap_err();
    assert!(matches!(
        err,
        EdnError::Parse(ParseError::TrailingContent { .. })
    ));
}

#[test]
fn test_error_empty_input() {
    assert!(matches!(
        parse("").unwrap_err(),
        EdnError::Parse(ParseError::UnexpectedEof { .. })
    ));
}

#[test]
fn test_error_whitespace_only() {
    assert!(matches!(
        parse("   ").unwrap_err(),
        EdnError::Parse(ParseError::UnexpectedEof { .. })
    ));
}

#[rstest]
#[case::vector("[1 2")]
#[case::map("{:a 1")]
#[case::list("(1 2")]
#[case::set("#{1 2")]
fn test_error_unclosed(#[case] input: &str) {
    assert!(parse(input).is_err());
}

#[test]
fn test_error_deep_nesting() {
    let depth = 100;
    let input = "[".repeat(depth) + &"]".repeat(depth);
    assert!(parse(&input).is_ok());
}

#[rstest]
#[case::malformed_discard("#_(#!V(\0\0\0\0\u{00ff}##")]
#[case::unterminated_string("#_\"\\")]
fn test_error_deser_regression_vec(#[case] input: &str) {
    let _ = from_str::<Vec<i64>>(input);
}

#[test]
fn test_error_deser_unicode_escape_multibyte_boundary() {
    let input = "\"\u{005c}u2\0`\u{07a0}";
    let _ = from_str::<String>(input);
}
