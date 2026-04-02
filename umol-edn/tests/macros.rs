//! Tests for the edn! proc macro.

use std::borrow::Cow;

use rstest::rstest;
use umol_edn::edn;
use umol_edn::edn::{Edn, Symbol};
use umol_edn::EdnMap;

#[test]
fn test_nil() {
    assert_eq!(edn!(nil), Edn::Nil);
}

#[rstest]
#[case(edn!(true), Edn::Bool(true))]
#[case(edn!(false), Edn::Bool(false))]
fn test_booleans(#[case] actual: Edn<'static>, #[case] expected: Edn<'static>) {
    assert_eq!(actual, expected);
}

#[rstest]
#[case(edn!(0), Edn::Int(0))]
#[case(edn!(1), Edn::Int(1))]
#[case(edn!(-1), Edn::Int(-1))]
#[case(edn!(+5), Edn::Int(5))]
fn test_integers(#[case] actual: Edn<'static>, #[case] expected: Edn<'static>) {
    assert_eq!(actual, expected);
}

#[test]
fn test_float() {
    match edn!(3.14) {
        Edn::Float(v) => assert!((v - 3.14).abs() < f64::EPSILON),
        other => panic!("expected Float, got {other:?}"),
    }
}

#[test]
fn test_negative_float() {
    match edn!(-2.5) {
        Edn::Float(v) => assert!((v - (-2.5)).abs() < f64::EPSILON),
        other => panic!("expected Float, got {other:?}"),
    }
}

#[test]
fn test_string() {
    assert_eq!(
        edn!("hello"),
        Edn::Str(Cow::Owned("hello".to_string()))
    );
}

#[test]
fn test_keyword() {
    assert_eq!(edn!(:foo), Edn::keyword("foo"));
}

#[test]
fn test_keyword_namespaced() {
    assert_eq!(edn!(:ns/name), Edn::keyword("ns/name"));
}

#[test]
fn test_symbol() {
    assert_eq!(edn!(foo), Edn::Symbol(Symbol::new("foo")));
}

#[test]
fn test_symbol_namespaced() {
    assert_eq!(edn!(ns/name), Edn::Symbol(Symbol::new("ns/name")));
}

#[test]
fn test_slash_symbol() {
    assert_eq!(edn!(/), Edn::Symbol(Symbol::new("/")));
}

#[test]
fn test_vector() {
    assert_eq!(
        edn!([1 2 3]),
        Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into())
    );
}

#[test]
fn test_list() {
    assert_eq!(
        edn!((1 2 3)),
        Edn::List(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into())
    );
}

#[test]
fn test_map() {
    let result = edn!({:a 1 :b 2});
    let mut expected = EdnMap::new();
    expected.insert(Edn::keyword("a"), Edn::Int(1));
    expected.insert(Edn::keyword("b"), Edn::Int(2));
    assert_eq!(result, Edn::Map(expected));
}

#[test]
fn test_set() {
    let result = edn!(#{1 2 3});
    if let Edn::Set(s) = result {
        assert_eq!(s.len(), 3);
    } else {
        panic!("expected Set");
    }
}

#[test]
fn test_nested() {
    let result = edn!({:items [1 2 3] :meta {:ok true}});
    assert!(result.is_map());
    assert_eq!(result.get("items").unwrap().as_vector().unwrap().len(), 3);
}

#[test]
fn test_tagged() {
    let result = edn!(#myapp/Person {:name "Alice"});
    assert!(matches!(result, Edn::Tagged(tag, _) if tag == "myapp/Person"));
}

#[test]
fn test_special_floats() {
    assert!(matches!(edn!(##NaN), Edn::Float(v) if v.is_nan()));
    assert_eq!(edn!(##Inf), Edn::Float(f64::INFINITY));
}

#[test]
fn test_discard() {
    assert_eq!(
        edn!([1 #_ 2 3]),
        Edn::Vector(vec![Edn::Int(1), Edn::Int(3)].into())
    );
}

#[test]
fn test_empty_collections() {
    assert_eq!(edn!([]), Edn::Vector(vec![].into()));
    assert_eq!(edn!(()), Edn::List(vec![].into()));
    assert_eq!(edn!({}), Edn::Map(EdnMap::new()));
}
