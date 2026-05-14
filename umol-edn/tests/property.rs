//! Property-based tests for Edn.

use std::borrow::Cow;
use std::hash::{DefaultHasher, Hash, Hasher};

use proptest::prelude::*;
use umol_edn::{read_string, Edn, EdnKeyword, EdnMap, EdnSeq, EdnSet, EdnSymbol};

const SYMBOL_START: &[char] = &[
    'a', 'b', 'c', 'z', 'A', 'B', 'Z', '.', '*', '+', '!', '-', '_', '?', '$', '%', '&', '=', '<',
    '>',
];

const SYMBOL_CONT: &[char] = &[
    'a', 'b', 'c', 'z', 'A', 'B', 'Z', '0', '1', '9', '.', '*', '+', '!', '-', '_', '?', '$', '%',
    '&', '=', '<', '>', '#', ':', '\'',
];

fn is_ambiguous_symbol(s: &str) -> bool {
    let mut chars = s.chars();
    match chars.next() {
        Some('+' | '-' | '.') => chars.next().is_some_and(|c| c.is_ascii_digit()),
        _ => false,
    }
}

fn symbol_name_strategy() -> impl Strategy<Value = String> {
    let simple = prop::sample::select(SYMBOL_START)
        .prop_flat_map(|first| {
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
    ("[a-zA-Z][a-zA-Z0-9]{0,5}", "[a-zA-Z][a-zA-Z0-9]{0,5}")
        .prop_map(|(ns, name)| format!("{ns}/{name}"))
}

fn string_strategy() -> impl Strategy<Value = String> {
    prop_oneof![
        // ASCII without escape-worthy characters.
        "[ -~]{0,20}",
        // Strings containing characters that require escaping.
        prop::collection::vec(
            prop_oneof![
                Just('"'),
                Just('\\'),
                Just('\n'),
                Just('\r'),
                Just('\t'),
                (0x20u32..0x7Fu32).prop_filter_map("printable", char::from_u32),
            ],
            1..15,
        )
        .prop_map(|chars| chars.into_iter().collect::<String>()),
        // Non-ASCII Unicode strings.
        prop::collection::vec(
            prop_oneof![
                Just('\u{00E9}'),  // é
                Just('\u{03B1}'),  // α
                Just('\u{4E16}'),  // 世
                Just('\u{1F600}'), // 😀
                (0x80u32..0x800u32).prop_filter_map("valid", char::from_u32),
            ],
            1..10,
        )
        .prop_map(|chars| chars.into_iter().collect::<String>()),
    ]
}

fn edn_leaf() -> impl Strategy<Value = Edn<'static>> {
    prop_oneof![
        Just(Edn::Nil),
        any::<bool>().prop_map(Edn::Bool),
        any::<i64>().prop_map(Edn::Int),
        (-1e15f64..1e15f64).prop_map(Edn::Float),
        // Printable ASCII.
        (0x21u32..0x7Fu32)
            .prop_filter_map("valid char", char::from_u32)
            .prop_map(Edn::Char),
        // Named characters.
        prop::sample::select(&['\n', '\r', ' ', '\t'][..]).prop_map(Edn::Char),
        // Non-ASCII chars that roundtrip via \uXXXX.
        prop::sample::select(&['\u{03B1}', '\u{00E9}', '\u{4E16}'][..]).prop_map(Edn::Char),
        // Control chars rendered as \uXXXX.
        (1u32..0x20u32)
            .prop_filter_map("not named", |cp| {
                let c = char::from_u32(cp)?;
                if matches!(c, '\n' | '\r' | '\t') {
                    None
                } else {
                    Some(c)
                }
            })
            .prop_map(Edn::Char),
        string_strategy().prop_map(|s| Edn::Str(Cow::Owned(s))),
        keyword_name_strategy().prop_map(|s| Edn::Keyword(EdnKeyword::owned(s))),
        symbol_name_strategy()
            .prop_filter("not a reserved word", |s| {
                s != "nil" && s != "true" && s != "false"
            })
            .prop_map(|s| Edn::Symbol(EdnSymbol::owned(s))),
        // Lone slash symbol.
        Just(Edn::Symbol(EdnSymbol::owned("/".to_string()))),
    ]
}

fn edn_strategy() -> impl Strategy<Value = Edn<'static>> {
    edn_leaf().prop_recursive(
        4,  // depth
        64, // max nodes
        8,  // items per collection
        |inner| {
            prop_oneof![
                prop::collection::vec(inner.clone(), 0..8).prop_map(|v| Edn::List(EdnSeq::from(v))),
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
                (tag_strategy(), inner)
                    .prop_map(|(tag, val)| { Edn::Tagged(Cow::Owned(tag), Box::new(val)) }),
            ]
        },
    )
}

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
