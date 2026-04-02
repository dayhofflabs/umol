//! Spec audit: verify umol-edn behavior against EDN spec statements.
//! Each test references a statement ID from discussion/66-edn-spec-conformance-2026-04-01.md.

use std::borrow::Cow;

use umol_edn::config::{Dialect, ParseConfig, TagReaders};
use umol_edn::edn::Edn;
use umol_edn::error::EdnError;
use umol_edn::{read_all_with, read_string, read_string_with};

fn clj() -> ParseConfig {
    ParseConfig {
        dialect: Dialect::Clojure,
        ..Default::default()
    }
}

fn edn() -> ParseConfig {
    ParseConfig {
        dialect: Dialect::Edn,
        ..Default::default()
    }
}

/// Helper: assert input parses successfully under both dialects.
fn ok_both(input: &str) {
    assert!(
        read_string_with(input, &clj()).is_ok(),
        "Clojure dialect rejected: {input}"
    );
    assert!(
        read_string_with(input, &edn()).is_ok(),
        "Edn dialect rejected: {input}"
    );
}

/// Helper: assert input fails under both dialects.
fn err_both(input: &str) {
    assert!(
        read_string_with(input, &clj()).is_err(),
        "Clojure dialect accepted: {input}"
    );
    assert!(
        read_string_with(input, &edn()).is_err(),
        "Edn dialect accepted: {input}"
    );
}

// -- S3: Delimiters need not be separated by whitespace ----------------------

#[test]
fn test_s3_delimiters_no_whitespace() {
    ok_both("[1[2]]");
    ok_both("{:a{:b 1}}");
    ok_both("(1(2))");
    // [1]2 is two values without whitespace between ] and 2.
    // read_all should parse both.
    assert_eq!(read_all_with("[1]2", &clj()).unwrap().len(), 2);
    assert_eq!(read_all_with("[1]2", &edn()).unwrap().len(), 2);
}

// -- S7: String escapes ------------------------------------------------------

#[test]
fn test_s7_string_escapes_both_dialects() {
    // The 5 spec escapes work in both dialects.
    ok_both(r#""a\tb""#);
    ok_both(r#""a\rb""#);
    ok_both(r#""a\nb""#);
    ok_both(r#""a\\b""#);
    ok_both(r#""a\"b""#);
}

#[test]
fn test_s7_string_clojure_escapes_accepted() {
    // Clojure extends with \b, \f, \uNNNN, octal.
    assert!(read_string_with(r#""\b""#, &clj()).is_ok());
    assert!(read_string_with(r#""\f""#, &clj()).is_ok());
    assert!(read_string_with(r#""\u0041""#, &clj()).is_ok());
    assert!(read_string_with(r#""\101""#, &clj()).is_ok()); // octal
}

#[test]
fn test_s7_string_clojure_escapes_rejected_edn() {
    // Edn dialect: only the 5 spec escapes.
    assert!(read_string_with(r#""\b""#, &edn()).is_err(), r#"Edn should reject \b in string"#);
    assert!(read_string_with(r#""\f""#, &edn()).is_err(), r#"Edn should reject \f in string"#);
    assert!(read_string_with(r#""\u0041""#, &edn()).is_err(), r#"Edn should reject \u in string"#);
    assert!(read_string_with(r#""\101""#, &edn()).is_err(), "Edn should reject octal in string");
}

// -- S8: Character literals --------------------------------------------------

#[test]
fn test_s8_char_both_dialects() {
    // Spec chars: \c, \newline, \return, \space, \tab, \uNNNN.
    ok_both(r"\a");
    ok_both(r"\newline");
    ok_both(r"\return");
    ok_both(r"\space");
    ok_both(r"\tab");
    ok_both(r"\u0041");
}

#[test]
fn test_s8_char_clojure_named_accepted() {
    assert!(read_string_with(r"\formfeed", &clj()).is_ok());
    assert!(read_string_with(r"\backspace", &clj()).is_ok());
}

#[test]
fn test_s8_char_clojure_named_rejected_edn() {
    assert!(read_string_with(r"\formfeed", &edn()).is_err(), r"Edn should reject \formfeed");
    assert!(read_string_with(r"\backspace", &edn()).is_err(), r"Edn should reject \backspace");
}

// -- S9b: Sign/dot first char, second must be non-numeric --------------------

#[test]
fn test_s9b_sign_dot_first_char() {
    // +a, -a, .a are symbols
    ok_both("+a");
    ok_both("-a");
    ok_both(".a");

    // +1, -1 are integers (not symbols)
    assert_eq!(
        read_string("+1").unwrap(),
        umol_edn::Edn::Int(1)
    );
    assert_eq!(
        read_string("-1").unwrap(),
        umol_edn::Edn::Int(-1)
    );
}

// -- S9d: Slash rules in symbols ---------------------------------------------

#[test]
fn test_s9d_slash_alone_valid() {
    ok_both("/");
}

#[test]
fn test_s9d_qualified_symbol_valid() {
    ok_both("ns/name");
    ok_both("my.ns/foo");
}

#[test]
fn test_s9d_empty_name_after_slash() {
    // ns/ — empty name part. Clojure and spec both reject.
    err_both("ns/");
}

#[test]
fn test_s9d_empty_prefix_before_slash() {
    // /name — empty prefix. Clojure and spec both reject.
    err_both("/name");
}

#[test]
fn test_s9d_multiple_slashes() {
    // a/b/c — multiple slashes. Spec says "once only".
    // Clojure accepts; Edn rejects.
    assert!(read_string_with("a/b/c", &clj()).is_ok());
    assert!(read_string_with("a/b/c", &edn()).is_err(), "Edn should reject a/b/c");
}

// -- S9e: Post-slash first char restriction ----------------------------------

#[test]
fn test_s9e_post_slash_digit_rejected() {
    // foo/1bar — digit after slash. Both Clojure and spec reject.
    err_both("foo/1bar");
}

#[test]
fn test_s9e_keyword_post_slash_digit_rejected() {
    // :foo/1bar — same rule applies to keywords.
    err_both(":foo/1bar");
}

#[test]
fn test_s9e_post_slash_valid_start() {
    ok_both("foo/_bar");
    ok_both("foo/bar");
}

// -- S10: Keywords follow symbol rules, prefixed with : ----------------------

#[test]
fn test_s10_bare_keyword_valid() {
    // Bare keywords (no namespace) are valid.
    ok_both(":foo");
    ok_both(":a");
}

#[test]
fn test_s10_qualified_keyword_valid() {
    ok_both(":my/foo");
    ok_both(":my.ns/foo");
}

#[test]
fn test_s10_keyword_namespace_is_symbol() {
    // Namespace must be a valid symbol — first char must be symbol-start.
    // :0/foo — namespace starts with digit. Clojure accepts; Edn rejects.
    assert!(read_string_with(":0/foo", &clj()).is_ok());
    assert!(read_string_with(":0/foo", &edn()).is_err(), "Edn should reject :0/foo");
}

#[test]
fn test_s10_keyword_first_char_restriction() {
    // After the :, the first char follows symbol-start rules.
    // :0 — digit start. Clojure accepts; Edn rejects.
    assert!(read_string_with(":0", &clj()).is_ok());
    assert!(read_string_with(":0", &edn()).is_err(), "Edn should reject :0");
    assert!(read_string_with(":0foo", &clj()).is_ok());
    assert!(read_string_with(":0foo", &edn()).is_err(), "Edn should reject :0foo");

    // :#foo — # is interior-only char. Both reject.
    err_both(":#foo");
}

#[test]
fn test_s10_keyword_special_start_chars() {
    // . + - are valid symbol-start chars, so also valid after :
    ok_both(":.foo");
    ok_both(":+foo");
    ok_both(":-foo");
}

#[test]
fn test_s10_keyword_post_slash_symbol_start() {
    // Post-slash char must be a valid symbol-start char.
    // Digits rejected in both dialects (clj also rejects):
    err_both(":foo/0bar");
    // # and : are interior-only — Edn rejects after slash, Clojure accepts:
    assert!(read_string_with(":foo/#bar", &clj()).is_ok());
    assert!(read_string_with(":foo/#bar", &edn()).is_err(), "Edn should reject :foo/#bar");
    assert!(read_string_with(":foo/:bar", &clj()).is_ok());
    assert!(read_string_with(":foo/:bar", &edn()).is_err(), "Edn should reject :foo/:bar");
    // # and : are fine as interior chars in the name part:
    ok_both(":foo/bar#baz");
    ok_both(":foo/bar:baz");
    // Valid start chars after slash:
    ok_both(":foo/.bar");
    ok_both(":foo/+bar");
    ok_both(":foo/-bar");
}

#[test]
fn test_s10b_double_colon_rejected() {
    // :: is auto-resolve in Clojure, not legal in EDN.
    assert!(
        read_string_with("::foo", &edn()).is_err(),
        "Edn dialect should reject ::foo"
    );
}

// -- S11b: No leading zeros --------------------------------------------------

#[test]
fn test_s11b_leading_zeros_rejected_edn() {
    // Spec: "No integer other than 0 may begin with 0."
    // Edn dialect: reject.
    assert!(read_string_with("007", &edn()).is_err(), "Edn dialect should reject 007");
    assert!(read_string_with("00", &edn()).is_err(), "Edn dialect should reject 00");
    // 0 alone is valid.
    assert!(read_string_with("0", &edn()).is_ok());
}

#[test]
fn test_s11b_leading_zeros_accepted_clojure() {
    // Clojure is lenient — accepts leading zeros.
    assert_eq!(read_string_with("007", &clj()).unwrap(), umol_edn::Edn::Int(7));
    assert_eq!(read_string_with("00", &clj()).unwrap(), umol_edn::Edn::Int(0));
}

// -- S17c: Reserved tags (bare unqualified tags) ----------------------------

#[test]
fn test_s17c_bare_tag_rejected_edn() {
    // Edn dialect rejects bare (unqualified) tags unless they are built-in.
    assert!(
        read_string_with("#foo 1", &edn()).is_err(),
        "Edn should reject bare #foo"
    );
    assert!(
        read_string_with("#bar [1 2]", &edn()).is_err(),
        "Edn should reject bare #bar"
    );
}

#[test]
fn test_s17c_bare_tag_accepted_clojure() {
    // Clojure dialect accepts bare tags — wraps as Tagged.
    let val = read_string_with("#foo 1", &clj()).unwrap();
    assert_eq!(val, Edn::Tagged("foo".into(), Box::new(Edn::Int(1))));
}

#[test]
fn test_s17c_qualified_tag_accepted_both() {
    // Qualified tags are always accepted in both dialects.
    ok_both("#my/foo 1");
}

#[test]
fn test_s17c_inst_accepted_edn() {
    // #inst is a built-in tag — accepted in Edn dialect when chrono feature is enabled.
    let val = read_string_with("#inst \"2024-01-01T00:00:00Z\"", &edn()).unwrap();
    assert_eq!(
        val,
        Edn::Tagged("inst".into(), Box::new(Edn::Str(Cow::Borrowed("2024-01-01T00:00:00Z"))))
    );
}

#[test]
fn test_s17c_uuid_accepted_edn() {
    // #uuid is a built-in tag — accepted in Edn dialect when uuid feature is enabled.
    let val = read_string_with(
        "#uuid \"f81d4fae-7dec-11d0-a765-00a0c91e6bf6\"",
        &edn(),
    )
    .unwrap();
    assert_eq!(
        val,
        Edn::Tagged(
            "uuid".into(),
            Box::new(Edn::Str(Cow::Borrowed("f81d4fae-7dec-11d0-a765-00a0c91e6bf6")))
        )
    );
}

#[test]
fn test_s17c_inst_invalid_rejected() {
    // Invalid RFC 3339 string is rejected by the #inst reader.
    assert!(read_string("#inst \"not-a-date\"").is_err());
    assert!(read_string("#inst \"2024-01-01\"").is_err()); // missing time
}

#[test]
fn test_s17c_uuid_invalid_rejected() {
    // Invalid UUID string is rejected by the #uuid reader.
    assert!(read_string("#uuid \"not-a-uuid\"").is_err());
}

#[test]
fn test_s17c_inst_non_string_rejected() {
    // #inst must be followed by a string.
    assert!(read_string("#inst 123").is_err());
}

// -- S20: Discard in both dialects ------------------------------------------

#[test]
fn test_s20_discard_both_dialects() {
    // #_ is discard in both dialects (spec S20).
    let clj_val = read_string_with("[1 #_ 2 3]", &clj()).unwrap();
    let edn_val = read_string_with("[1 #_ 2 3]", &edn()).unwrap();
    assert_eq!(clj_val, Edn::Vector(vec![Edn::Int(1), Edn::Int(3)]));
    assert_eq!(edn_val, Edn::Vector(vec![Edn::Int(1), Edn::Int(3)]));
}

// -- Tag reader dispatch ---------------------------------------------------

#[test]
fn test_custom_tag_reader_dispatch() {
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
    let val = read_string_with("#double 5", &config).unwrap();
    assert_eq!(val, Edn::Int(10));
}

#[test]
fn test_custom_tag_reader_error_propagation() {
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
