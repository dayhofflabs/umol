//! Parity test matrix: variant × wrapper × ser/de path.
//!
//! For each EDN-specific wrapper (`EdnKeyword`, `EdnSymbol`, `EdnList`,
//! `EdnHashSet`, `EdnTagged`, `EdnBigInt`, `EdnBigDecimal`, `Value`), this
//! suite verifies four cells:
//!
//! 1. Native `ToEdn`/`FromEdn` round-trip.
//! 2. Serde `to_string` / `from_str` round-trip over EDN.
//! 3. Serde `to_string` / `from_str` when the wrapper is a field on a
//!    normal struct (composition with the struct map dispatch).
//! 4. JSON fallback shape — the wrapper must degrade predictably when
//!    the serializer is not EDN-aware.
//!
//! It also covers the remaining Phase 4 spot-checks:
//! - `EdnTagged<T>` next to an enum-variant-tagged enum in the same struct.
//! - Wrapper composition with `#[serde(default)]`, `#[serde(rename)]`,
//!   and `#[serde(flatten)]`.

#![cfg(feature = "serde")]

use std::collections::HashMap;

use rstest::rstest;
use serde::{Deserialize, Serialize};

use umol_edn::{
    config::ParseConfig,
    de::{from_str, from_str_with},
    edn::Edn,
    native::{FromEdn, ToEdn},
    ser::to_string,
    EdnHashSet, EdnKeyword, EdnList, EdnSymbol, EdnTagged, Value,
};

/// Parse with `allow_unknown_tags = true`, for tests that exercise dynamic
/// tags (`EdnTagged<T>` with caller-chosen tag names, `#Variant` enum
/// dispatch, or `Value` containing unknown tagged literals).
fn from_str_permissive<'a, T: serde::Deserialize<'a>>(s: &'a str) -> T {
    let mut config = ParseConfig::default();
    config.allow_unknown_tags = true;
    from_str_with(s, &config).unwrap()
}

#[cfg(feature = "bignum")]
use umol_edn::{EdnBigDecimal, EdnBigInt};

#[test]
fn test_parity_keyword_native_roundtrip() {
    let k = EdnKeyword::new("foo");
    let edn = k.to_edn();
    let back = EdnKeyword::from_edn(&edn).unwrap();
    assert_eq!(k, back);
}

#[test]
fn test_parity_keyword_serde_edn_roundtrip() {
    let k = EdnKeyword::new("foo");
    let s = to_string(&k).unwrap();
    assert_eq!(s, ":foo");
    let back: EdnKeyword = from_str(&s).unwrap();
    assert_eq!(k, back);
}

#[test]
fn test_parity_keyword_in_struct() {
    #[derive(Serialize, Deserialize, PartialEq, Debug)]
    struct S {
        tag: EdnKeyword,
    }
    let v = S {
        tag: EdnKeyword::new("active"),
    };
    let s = to_string(&v).unwrap();
    assert_eq!(s, "{:tag :active}");
    let back: S = from_str(&s).unwrap();
    assert_eq!(v, back);
}

#[test]
fn test_parity_keyword_json_fallback_is_string() {
    let k = EdnKeyword::new("foo");
    let j = serde_json::to_string(&k).unwrap();
    assert_eq!(j, "\"foo\"");
    let back: EdnKeyword = serde_json::from_str(&j).unwrap();
    assert_eq!(back, k);
}

#[test]
fn test_parity_keyword_native_rejects_string() {
    let edn = Edn::Str(std::borrow::Cow::Borrowed("not a keyword"));
    assert!(EdnKeyword::from_edn(&edn).is_err());
}

#[test]
fn test_parity_keyword_serde_accepts_string_for_json_interop() {
    let kw: EdnKeyword = from_str(r#""bar""#).unwrap();
    assert_eq!(kw.as_str(), "bar");
}

#[test]
fn test_parity_symbol_native_roundtrip() {
    let s = EdnSymbol::new("sym");
    let edn = s.to_edn();
    let back = EdnSymbol::from_edn(&edn).unwrap();
    assert_eq!(s, back);
}

#[test]
fn test_parity_symbol_serde_edn_roundtrip() {
    let s = EdnSymbol::new("sym");
    let ser = to_string(&s).unwrap();
    assert_eq!(ser, "sym");
    let back: EdnSymbol = from_str(&ser).unwrap();
    assert_eq!(s, back);
}

#[test]
fn test_parity_symbol_json_fallback_is_string() {
    let s = EdnSymbol::new("sym");
    let j = serde_json::to_string(&s).unwrap();
    assert_eq!(j, "\"sym\"");
}

#[test]
fn test_parity_symbol_rejects_keyword() {
    let result: Result<EdnSymbol, _> = from_str(":foo");
    assert!(result.is_err());
}

#[test]
fn test_parity_list_native_roundtrip() {
    let l: EdnList<i64> = vec![1, 2, 3].into();
    let edn = l.to_edn();
    let back = EdnList::<i64>::from_edn(&edn).unwrap();
    assert_eq!(l, back);
}

#[test]
fn test_parity_list_serde_edn_roundtrip() {
    let l: EdnList<i64> = vec![1, 2, 3].into();
    let s = to_string(&l).unwrap();
    assert_eq!(s, "(1 2 3)");
    let back: EdnList<i64> = from_str(&s).unwrap();
    assert_eq!(l, back);
}

#[test]
fn test_parity_list_rejects_vector() {
    let result: Result<EdnList<i64>, _> = from_str("[1 2 3]");
    assert!(result.is_err());
}

#[test]
fn test_parity_list_in_struct() {
    #[derive(Serialize, Deserialize, PartialEq, Debug)]
    struct S {
        items: EdnList<String>,
    }
    let v = S {
        items: vec!["a".into(), "b".into()].into(),
    };
    let s = to_string(&v).unwrap();
    assert_eq!(s, "{:items (\"a\" \"b\")}");
    let back: S = from_str(&s).unwrap();
    assert_eq!(v, back);
}

#[test]
fn test_parity_list_json_fallback_is_array() {
    let l: EdnList<i64> = vec![1, 2, 3].into();
    let j = serde_json::to_string(&l).unwrap();
    assert_eq!(j, "[1,2,3]");
    let back: EdnList<i64> = serde_json::from_str(&j).unwrap();
    assert_eq!(l, back);
}

#[test]
fn test_parity_set_native_roundtrip() {
    let s: EdnHashSet<i64> = [1i64, 2, 3].into_iter().collect();
    let edn = s.to_edn();
    let back = EdnHashSet::<i64>::from_edn(&edn).unwrap();
    assert_eq!(s, back);
}

#[test]
fn test_parity_set_serde_edn_roundtrip() {
    let s: EdnHashSet<i64> = [1i64, 2, 3].into_iter().collect();
    let ser = to_string(&s).unwrap();
    assert!(ser.starts_with("#{") && ser.ends_with('}'));
    let back: EdnHashSet<i64> = from_str(&ser).unwrap();
    assert_eq!(s, back);
}

#[test]
fn test_parity_set_rejects_vector() {
    let result: Result<EdnHashSet<i64>, _> = from_str("[1 2 3]");
    assert!(result.is_err());
}

#[test]
fn test_parity_set_json_fallback_is_array() {
    let s: EdnHashSet<i64> = [7i64].into_iter().collect();
    let j = serde_json::to_string(&s).unwrap();
    assert_eq!(j, "[7]");
}

#[test]
fn test_parity_tagged_native_roundtrip() {
    let t = EdnTagged::new("score".to_string(), 17i64);
    let edn = t.to_edn();
    let back = EdnTagged::<i64>::from_edn(&edn).unwrap();
    assert_eq!(t, back);
}

#[test]
fn test_parity_tagged_serde_edn_roundtrip() {
    let t = EdnTagged::new("score".to_string(), 17i64);
    let s = to_string(&t).unwrap();
    assert_eq!(s, "#score 17");
    let back: EdnTagged<i64> = from_str_permissive(&s);
    assert_eq!(t, back);
}

#[test]
fn test_parity_tagged_in_struct() {
    #[derive(Serialize, Deserialize, PartialEq, Debug)]
    struct S {
        marker: EdnTagged<String>,
    }
    let v = S {
        marker: EdnTagged::new("inst", "2026-04-08".to_string()),
    };
    let s = to_string(&v).unwrap();
    assert_eq!(s, "{:marker #inst \"2026-04-08\"}");
    let back: S = from_str(&s).unwrap();
    assert_eq!(v, back);
}

#[test]
fn test_parity_tagged_json_fallback_is_tuple() {
    let t = EdnTagged::new("score".to_string(), 17i64);
    let j = serde_json::to_string(&t).unwrap();
    assert_eq!(j, "[\"score\",17]");
    let back: EdnTagged<i64> = serde_json::from_str(&j).unwrap();
    assert_eq!(back, t);
}

#[derive(Serialize, Deserialize, PartialEq, Debug)]
enum Event {
    Click(i64),
    Key(String),
}

#[derive(Serialize, Deserialize, PartialEq, Debug)]
struct MixedTagged {
    dynamic: EdnTagged<String>,
    event: Event,
}

#[test]
fn test_parity_tagged_coexist_with_enum_variant_tagged() {
    let v = MixedTagged {
        dynamic: EdnTagged::new("uuid", "abc-123".to_string()),
        event: Event::Click(17),
    };
    let s = to_string(&v).unwrap();
    assert_eq!(s, "{:dynamic #uuid \"abc-123\" :event #Click 17}");
    let back: MixedTagged = from_str_permissive(&s);
    assert_eq!(v, back);
}

#[test]
fn test_parity_tagged_coexist_nested_variant() {
    let v = MixedTagged {
        dynamic: EdnTagged::new("version", "v4".to_string()),
        event: Event::Key("Enter".into()),
    };
    let s = to_string(&v).unwrap();
    let back: MixedTagged = from_str_permissive(&s);
    assert_eq!(v, back);
}

#[cfg(feature = "bignum")]
#[test]
fn test_parity_bigint_native_roundtrip() {
    use std::str::FromStr;
    let n = EdnBigInt(num_bigint::BigInt::from_str("12345678901234567890").unwrap());
    let edn = n.to_edn();
    let back = EdnBigInt::from_edn(&edn).unwrap();
    assert_eq!(n, back);
}

#[cfg(feature = "bignum")]
#[test]
fn test_parity_bigint_serde_edn_roundtrip() {
    use std::str::FromStr;
    let n = EdnBigInt(num_bigint::BigInt::from_str("12345678901234567890").unwrap());
    let s = to_string(&n).unwrap();
    assert_eq!(s, "12345678901234567890N");
    let back: EdnBigInt = from_str(&s).unwrap();
    assert_eq!(n, back);
}

#[cfg(feature = "bignum")]
#[test]
fn test_parity_bigint_json_fallback_is_string() {
    use std::str::FromStr;
    let n = EdnBigInt(num_bigint::BigInt::from_str("17").unwrap());
    let j = serde_json::to_string(&n).unwrap();
    assert_eq!(j, "\"17\"");
    let back: EdnBigInt = serde_json::from_str(&j).unwrap();
    assert_eq!(back, n);
}

#[cfg(feature = "bignum")]
#[test]
fn test_parity_bigdecimal_serde_edn_roundtrip() {
    use bigdecimal::BigDecimal;
    use std::str::FromStr;
    let d = EdnBigDecimal::new(BigDecimal::from_str("3.14159").unwrap());
    let s = to_string(&d).unwrap();
    assert_eq!(s, "3.14159M");
    let back: EdnBigDecimal = from_str(&s).unwrap();
    assert_eq!(d, back);
}

#[rstest]
#[case("nil")]
#[case("true")]
#[case(":kw")]
#[case("sym")]
#[case("[1 2 3]")]
#[case("(1 2 3)")]
#[case("#{1 2}")]
#[case(r#"{:name "salt" :count 2}"#)]
#[case(r#"#inst "2026-04-08""#)]
fn test_parity_value_lossless_edn_roundtrip(#[case] input: &str) {
    let v: Value = from_str(input).unwrap();
    let s = to_string(&v).unwrap();
    let v2: Value = from_str(&s).unwrap();
    assert_eq!(v, v2);
}

#[test]
fn test_parity_value_as_field() {
    #[derive(Serialize, Deserialize, PartialEq, Debug)]
    struct S {
        meta: Value,
    }
    let v = S {
        meta: Value::parse("(1 :two [3 4])").unwrap(),
    };
    let s = to_string(&v).unwrap();
    let back: S = from_str(&s).unwrap();
    assert_eq!(v, back);
}

#[test]
fn test_parity_serde_default_on_wrapper() {
    #[derive(Serialize, Deserialize, PartialEq, Debug)]
    struct S {
        name: String,
        #[serde(default)]
        kind: Option<EdnKeyword>,
    }
    let v: S = from_str(r#"{:name "x"}"#).unwrap();
    assert_eq!(
        v,
        S {
            name: "x".into(),
            kind: None
        }
    );
    let v = S {
        name: "x".into(),
        kind: Some(EdnKeyword::new("active")),
    };
    let s = to_string(&v).unwrap();
    let back: S = from_str(&s).unwrap();
    assert_eq!(v, back);
}

#[test]
fn test_parity_serde_rename_on_wrapper() {
    #[derive(Serialize, Deserialize, PartialEq, Debug)]
    struct S {
        #[serde(rename = "type")]
        kind: EdnKeyword,
    }
    let v = S {
        kind: EdnKeyword::new("molecule"),
    };
    let s = to_string(&v).unwrap();
    assert_eq!(s, "{:type :molecule}");
    let back: S = from_str(&s).unwrap();
    assert_eq!(v, back);
}

#[test]
fn test_parity_serde_flatten_over_wrapper_field() {
    #[derive(Serialize, Deserialize, PartialEq, Debug)]
    struct Inner {
        kind: EdnKeyword,
        count: i64,
    }
    #[derive(Serialize, Deserialize, PartialEq, Debug)]
    struct Outer {
        name: String,
        #[serde(flatten)]
        inner: Inner,
    }
    let v = Outer {
        name: "salt".into(),
        inner: Inner {
            kind: EdnKeyword::new("mineral"),
            count: 2,
        },
    };
    let s = to_string(&v).unwrap();
    let back: Outer = from_str(&s).unwrap();
    assert_eq!(v, back);
}

#[test]
fn test_parity_multiple_wrappers_in_one_struct() {
    #[derive(Serialize, Deserialize, PartialEq, Debug)]
    struct S {
        name: EdnKeyword,
        aliases: EdnList<String>,
        ids: EdnHashSet<i64>,
        extra: Value,
    }
    let mut config = ParseConfig::default();
    config.allow_unknown_tags = true;
    let v = S {
        name: EdnKeyword::new("salt"),
        aliases: vec!["NaCl".into(), "halite".into()].into(),
        ids: [1i64, 2, 3].into_iter().collect(),
        extra: Value::parse_with("#custom {:k 1}", &config).unwrap(),
    };
    let s = to_string(&v).unwrap();
    let back: S = from_str_permissive(&s);
    assert_eq!(v, back);
}

#[test]
fn test_parity_hashmap_string_to_wrapper() {
    let mut m: HashMap<String, EdnKeyword> = HashMap::new();
    m.insert("mode".into(), EdnKeyword::new("strict"));
    m.insert("role".into(), EdnKeyword::new("admin"));
    let s = to_string(&m).unwrap();
    let back: HashMap<String, EdnKeyword> = from_str(&s).unwrap();
    assert_eq!(m, back);
}

#[rstest]
#[case::unit("nil", Edn::Nil)]
#[case::bool_true("true", Edn::Bool(true))]
#[case::bool_false("false", Edn::Bool(false))]
#[case::int("-7", Edn::Int(-7))]
#[case::float("2.5", Edn::Float(2.5))]
#[case::str(r#""abc""#, Edn::Str(std::borrow::Cow::Borrowed("abc")))]
fn test_parity_primitive_edn_roundtrip(#[case] input: &str, #[case] expected: Edn<'static>) {
    let v: Value = from_str(input).unwrap();
    assert_eq!(v.as_edn(), &expected);
}
