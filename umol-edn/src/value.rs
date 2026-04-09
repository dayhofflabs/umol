//! Owned, serde-friendly EDN value.
//!
//! [`Value`] is an owned wrapper around [`Edn<'static>`] that carries the
//! full variant information of the internal EDN data model across the serde
//! boundary. Unlike `Edn<'a>`, which is borrow-parameterized to support
//! zero-copy parsing, `Value` has no lifetime and can be stored in struct
//! fields, sent across threads, and used as a dynamic-typed escape hatch
//! inside otherwise-typed deserialize paths.
//!
//! The underlying `Edn<'static>` is exposed via [`Value::as_edn`] and
//! [`Value::into_edn`]; [`From`] conversions round-trip without copying
//! further than what `Edn::into_owned` already requires.
//!
//! Serialize/Deserialize implementations live in `value_serde.rs` under
//! the `serde` feature.

use std::fmt;
use std::str::FromStr;

use crate::edn::Edn;
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};
use crate::reader::read_string;

/// Owned EDN value. A lossless mirror of [`Edn<'static>`] for use where a
/// fully owned, lifetime-free value is required.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Value(pub(crate) Edn<'static>);

impl Value {
    /// Wrap an already-owned `Edn<'static>` value.
    pub fn new(edn: Edn<'static>) -> Self {
        Self(edn)
    }

    /// Parse an EDN string into a `Value`.
    pub fn parse(input: &str) -> Result<Self, EdnError> {
        Ok(Self(read_string(input)?.into_owned()))
    }

    /// Borrow the inner `Edn<'static>`.
    pub fn as_edn(&self) -> &Edn<'static> {
        &self.0
    }

    /// Unwrap into the inner `Edn<'static>`.
    pub fn into_edn(self) -> Edn<'static> {
        self.0
    }
}

impl From<Edn<'static>> for Value {
    fn from(edn: Edn<'static>) -> Self {
        Self(edn)
    }
}

impl From<&Edn<'_>> for Value {
    fn from(edn: &Edn<'_>) -> Self {
        Self(edn.clone().into_owned())
    }
}

impl From<Value> for Edn<'static> {
    fn from(v: Value) -> Self {
        v.0
    }
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl FromStr for Value {
    type Err = EdnError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

impl<'de> FromEdn<'de> for Value {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        Ok(Self(edn.clone().into_owned()))
    }
}

impl ToEdn for Value {
    fn to_edn(&self) -> Edn<'_> {
        self.0.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_value_parse_primitives() {
        assert_eq!(Value::parse("nil").unwrap().as_edn(), &Edn::Nil);
        assert_eq!(Value::parse("true").unwrap().as_edn(), &Edn::Bool(true));
        assert_eq!(Value::parse("7").unwrap().as_edn(), &Edn::Int(7));
    }

    #[test]
    fn test_value_parse_roundtrip() {
        let input = "{:name \"salt\" :atoms [:Na :Cl]}";
        let v = Value::parse(input).unwrap();
        let s = v.to_string();
        let v2 = Value::parse(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[test]
    fn test_value_parse_tagged() {
        let v = Value::parse("#inst \"2026-04-08\"").unwrap();
        match v.as_edn() {
            Edn::Tagged(tag, inner) => {
                assert_eq!(tag.as_ref(), "inst");
                assert_eq!(inner.as_ref(), &Edn::Str(std::borrow::Cow::Borrowed("2026-04-08")));
            }
            other => panic!("expected tagged, got {other:?}"),
        }
    }

    #[test]
    fn test_value_from_edn_and_to_edn() {
        let edn = read_string("[1 2 3]").unwrap();
        let v = Value::from_edn(&edn).unwrap();
        let back = v.to_edn();
        assert_eq!(back, edn);
    }

    #[test]
    fn test_value_display_matches_edn() {
        let v = Value::parse(":foo").unwrap();
        assert_eq!(v.to_string(), ":foo");
    }

    #[test]
    fn test_value_from_str() {
        let v: Value = "(1 2 3)".parse().unwrap();
        assert!(matches!(v.as_edn(), Edn::List(_)));
    }
}
