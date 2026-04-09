//! Property-based tests for Edn.

use std::borrow::Cow;
use std::hash::{DefaultHasher, Hash, Hasher};

use proptest::prelude::*;
use umol_edn::{read_string, Edn, EdnMap, EdnSeq, EdnSet, Keyword, Symbol};

// Strategies

/// Characters valid as the first character of a symbol name.
const SYMBOL_START: &[char] = &[
    'a', 'b', 'c', 'z', 'A', 'B', 'Z', '.', '*', '+', '!', '-', '_', '?', '$', '%', '&', '=',
    '<', '>',
];

/// Additional characters valid in the rest of a symbol name.
const SYMBOL_CONT: &[char] = &[
    'a', 'b', 'c', 'z', 'A', 'B', 'Z', '0', '1', '9', '.', '*', '+', '!', '-', '_', '?', '$',
    '%', '&', '=', '<', '>', '#', ':', '\'',
];

/// Returns true if first char is +/-/. and second char is a digit (ambiguous with numbers).
fn is_ambiguous_symbol(s: &str) -> bool {
    let mut chars = s.chars();
    match chars.next() {
        Some('+' | '-' | '.') => chars.next().is_some_and(|c| c.is_ascii_digit()),
        _ => false,
    }
}

fn symbol_name_strategy() -> impl Strategy<Value = String> {
    let simple = prop::sample::select(SYMBOL_START).prop_flat_map(|first| {
        prop::collection::vec(prop::sample::select(SYMBOL_CONT), 0..6).prop_map(move |rest| {
            let mut s = String::new();
            s.push(first);
            for c in &rest {
                s.push(*c);
            }
            s
        })
    })
    .prop_filter("ambiguous with number", |s| !is_ambiguous_symbol(s));

    // Namespace-qualified: ns/name — alpha first chars to avoid sign/dot ambiguity.
    let qualified = (
        prop::sample::select(&SYMBOL_START[..7]), // letters only (a,b,c,z,A,B,Z)
        prop::collection::vec(prop::sample::select(&SYMBOL_CONT[..15]), 0..4),
        prop::sample::select(&SYMBOL_START[..7]),
        prop::collection::vec(prop::sample::select(&SYMBOL_CONT[..15]), 0..4),
    )
        .prop_map(|(ns_first, ns_rest, name_first, name_rest)| {
            let mut s = String::new();
            s.push(ns_first);
            for c in &ns_rest {
                s.push(*c);
            }
            s.push('/');
            s.push(name_first);
            for c in &name_rest {
                s.push(*c);
            }
            s
        });

    prop_oneof![9 => simple, 1 => qualified]
}

fn keyword_name_strategy() -> impl Strategy<Value = String> {
    symbol_name_strategy().prop_filter("keywords cannot be / alone", |s| s != "/")
}

fn tag_strategy() -> impl Strategy<Value = String> {
    (
        "[a-zA-Z][a-zA-Z0-9]{0,5}",
        "[a-zA-Z][a-zA-Z0-9]{0,5}",
    )
        .prop_map(|(ns, name)| format!("{ns}/{name}"))
}

fn edn_leaf() -> impl Strategy<Value = Edn<'static>> {
    prop_oneof![
        Just(Edn::Nil),
        any::<bool>().prop_map(Edn::Bool),
        any::<i64>().prop_map(Edn::Int),
        // Finite floats only: NaN/Inf cannot roundtrip through display.
        (-1e15f64..1e15f64).prop_map(Edn::Float),
        // Printable ASCII chars (safe for roundtrip).
        (0x21u32..0x7Fu32)
            .prop_filter_map("valid char", |cp| char::from_u32(cp))
            .prop_map(Edn::Char),
        // Named characters.
        prop::sample::select(&['\n', '\r', ' ', '\t'][..]).prop_map(Edn::Char),
        "[ -~]{0,20}".prop_map(|s| Edn::Str(Cow::Owned(s))),
        keyword_name_strategy().prop_map(|s| Edn::Keyword(Keyword::owned(s))),
        symbol_name_strategy()
            .prop_filter("not a reserved word", |s| {
                s != "nil" && s != "true" && s != "false"
            })
            .prop_map(|s| Edn::Symbol(Symbol::owned(s))),
    ]
}

fn edn_strategy() -> impl Strategy<Value = Edn<'static>> {
    edn_leaf().prop_recursive(
        4,  // depth
        64, // max nodes
        8,  // items per collection
        |inner| {
            prop_oneof![
                prop::collection::vec(inner.clone(), 0..8)
                    .prop_map(|v| Edn::List(EdnSeq::from(v))),
                prop::collection::vec(inner.clone(), 0..8)
                    .prop_map(|v| Edn::Vector(EdnSeq::from(v))),
                prop::collection::vec((inner.clone(), inner.clone()), 0..4).prop_map(|pairs| {
                    let mut m = EdnMap::with_capacity(pairs.len());
                    for (k, v) in pairs {
                        m.insert(k, v);
                    }
                    Edn::Map(m)
                }),
                prop::collection::vec(inner.clone(), 0..6).prop_map(|items| {
                    let mut s = EdnSet::new();
                    for item in items {
                        s.insert(item);
                    }
                    Edn::Set(s)
                }),
                (tag_strategy(), inner).prop_map(|(tag, val)| {
                    Edn::Tagged(Cow::Owned(tag), Box::new(val))
                }),
            ]
        },
    )
}

// Properties

fn hash_of(v: &Edn<'_>) -> u64 {
    let mut h = DefaultHasher::new();
    v.hash(&mut h);
    h.finish()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1000))]

    #[test]
    fn test_display_parse_roundtrip(edn in edn_strategy()) {
        let displayed = edn.to_string();
        let reparsed = read_string(&displayed)
            .map(|v| v.into_owned())
            .map_err(|e| TestCaseError::fail(format!("{}: {}", e, displayed)))?;
        prop_assert_eq!(edn, reparsed);
    }

    #[test]
    fn test_eq_implies_hash(a in edn_strategy(), b in edn_strategy()) {
        if a == b {
            prop_assert_eq!(hash_of(&a), hash_of(&b));
        }
    }

    #[test]
    fn test_eq_reflexive(edn in edn_strategy()) {
        prop_assert!(edn == edn);
    }

    #[test]
    fn test_display_is_valid_edn(edn in edn_strategy()) {
        let displayed = edn.to_string();
        read_string(&displayed)
            .map_err(|e| TestCaseError::fail(format!("{}: {}", e, displayed)))?;
    }

    #[test]
    fn test_double_roundtrip(edn in edn_strategy()) {
        let d1 = edn.to_string();
        let p1 = read_string(&d1)
            .map(|v| v.into_owned())
            .map_err(|e| TestCaseError::fail(format!("{}: {}", e, d1)))?;
        let d2 = p1.to_string();
        let p2 = read_string(&d2)
            .map(|v| v.into_owned())
            .map_err(|e| TestCaseError::fail(format!("{}: {}", e, d2)))?;
        prop_assert_eq!(p1, p2);
        prop_assert_eq!(d1, d2);
    }
}
