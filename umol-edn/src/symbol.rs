//! Owned EDN symbol type for use in struct fields.
//!
//! `EdnSymbol` is the owned, lifetime-free counterpart to the borrowing
//! [`Symbol<'a>`](crate::Symbol) used in the `Edn` data model. Owned structs
//! can carry it to statically commit to symbol shape (rather than string or
//! keyword) and round-trip symbols losslessly through the serde layer.
//!
//! Under the `serde` feature, `Serialize` emits a bare `name` token via the
//! `SYMBOL_TOKEN` newtype-struct trick when routed through `EdnSerializer`;
//! non-EDN serializers see a transparent string. `Deserialize` via the EDN
//! tree path accepts only `Edn::Symbol`; strings and keywords are rejected.

use std::fmt;
use std::ops::Deref;

#[cfg(feature = "serde")]
use serde::de::{self, Deserializer, Visitor};
#[cfg(feature = "serde")]
use serde::ser::Serializer;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::edn::{Edn, Symbol};
use crate::error::DeError;
#[cfg(feature = "serde")]
use crate::serde_tokens::SYMBOL_TOKEN;
use crate::traits::{FromEdn, ToEdn};

/// An owned EDN symbol (`name` or `ns/name`). Always serializes as a symbol
/// in EDN.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdnSymbol(String);

impl EdnSymbol {
    pub fn new(s: impl Into<String>) -> Self {
        Self(s.into())
    }

    pub fn into_inner(self) -> String {
        self.0
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl Deref for EdnSymbol {
    type Target = str;
    fn deref(&self) -> &str {
        &self.0
    }
}

impl AsRef<str> for EdnSymbol {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for EdnSymbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<String> for EdnSymbol {
    fn from(s: String) -> Self {
        Self(s)
    }
}

impl From<&str> for EdnSymbol {
    fn from(s: &str) -> Self {
        Self(s.to_string())
    }
}

impl<'de> FromEdn<'de> for EdnSymbol {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Symbol(s) => Ok(EdnSymbol(s.as_str().to_string())),
            other => Err(DeError::TypeMismatch {
                expected: "symbol",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for EdnSymbol {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Symbol(Symbol::new(&self.0))
    }
}

#[cfg(feature = "serde")]
impl Serialize for EdnSymbol {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(SYMBOL_TOKEN, self.as_str())
    }
}

#[cfg(feature = "serde")]
struct SymbolVisitor;

#[cfg(feature = "serde")]
impl<'de> Visitor<'de> for SymbolVisitor {
    type Value = EdnSymbol;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("an EDN symbol")
    }

    fn visit_str<E: de::Error>(self, v: &str) -> Result<EdnSymbol, E> {
        Ok(EdnSymbol::new(v))
    }

    fn visit_string<E: de::Error>(self, v: String) -> Result<EdnSymbol, E> {
        Ok(EdnSymbol::new(v))
    }

    fn visit_borrowed_str<E: de::Error>(self, v: &'de str) -> Result<EdnSymbol, E> {
        Ok(EdnSymbol::new(v))
    }

    fn visit_newtype_struct<D: Deserializer<'de>>(
        self,
        deserializer: D,
    ) -> Result<EdnSymbol, D::Error> {
        String::deserialize(deserializer).map(EdnSymbol::new)
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for EdnSymbol {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_newtype_struct(SYMBOL_TOKEN, SymbolVisitor)
    }
}

#[cfg(test)]
mod tests {
    #[cfg(feature = "serde")]
    use rstest::rstest;

    use super::*;
    use crate::read_string;
    #[cfg(feature = "serde")]
    use crate::{from_str, to_string};

    #[test]
    fn test_edn_symbol_display() {
        assert_eq!(EdnSymbol::new("foo").to_string(), "foo");
    }

    #[test]
    fn test_edn_symbol_from_edn() {
        let edn = read_string("bar").unwrap();
        let sym = EdnSymbol::from_edn(&edn).unwrap();
        assert_eq!(sym.as_str(), "bar");
    }

    #[test]
    fn test_edn_symbol_from_edn_error() {
        let edn = read_string(r#""bar""#).unwrap();
        assert!(EdnSymbol::from_edn(&edn).is_err());
    }

    #[test]
    fn test_edn_symbol_to_edn() {
        let symbol = EdnSymbol::new("baz");
        let edn = symbol.to_edn();
        let Edn::Symbol(s) = &edn else {
            panic!("expected Edn::Symbol, got {:?}", edn);
        };
        assert_eq!(s.as_str(), "baz");
    }

    #[cfg(feature = "serde")]
    #[rstest]
    #[case("foo", "foo")]
    #[case("ns/foo", "ns/foo")]
    fn test_edn_symbol_serialize(#[case] input: &str, #[case] expected: &str) {
        assert_eq!(to_string(&EdnSymbol::new(input)).unwrap(), expected);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_symbol_serialize_json() {
        let json = serde_json::to_string(&EdnSymbol::new("foo")).unwrap();
        assert_eq!(json, r#""foo""#);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_symbol_serialize_in_vec() {
        let v = vec![EdnSymbol::new("a"), EdnSymbol::new("b")];
        assert_eq!(to_string(&v).unwrap(), "[a b]");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_symbol_serialize_in_struct() {
        #[derive(serde::Serialize)]
        struct Wrapper {
            tag: EdnSymbol,
            value: i64,
        }
        let w = Wrapper {
            tag: EdnSymbol::new("foo"),
            value: 5,
        };
        assert_eq!(to_string(&w).unwrap(), "{:tag foo :value 5}");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_symbol_deserialize() {
        let sym: EdnSymbol = from_str("foo").unwrap();
        assert_eq!(sym.as_str(), "foo");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_symbol_deserialize_from_string_error() {
        let result: Result<EdnSymbol, _> = from_str(r#""foo""#);
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_symbol_deserialize_from_keyword_error() {
        let result: Result<EdnSymbol, _> = from_str(":foo");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_symbol_roundtrip() {
        let sym = EdnSymbol::new("ns/foo");
        let edn = to_string(&sym).unwrap();
        let parsed: EdnSymbol = from_str(&edn).unwrap();
        assert_eq!(parsed, sym);
    }
}
