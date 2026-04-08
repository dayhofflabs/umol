//! Serde Serialize/Deserialize impls for [`EdnSymbol`].
//!
//! `EdnSymbol` lives in [`crate::symbol_owned`] (always available); the
//! impls in this file are only compiled with the `serde` feature. When
//! serialized via `EdnSerializer`, an `EdnSymbol` emits a bare `name`
//! token via the `SYMBOL_TOKEN` newtype-struct trick. Non-EDN serializers
//! see a transparent string. When deserialized via the EDN tree path,
//! only `Edn::Symbol` is accepted; strings and keywords are rejected.

use std::fmt;

use serde::de::{self, Deserializer, Visitor};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};

use crate::serde_tokens::SYMBOL_TOKEN;
use crate::symbol_owned::EdnSymbol;

impl Serialize for EdnSymbol {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(SYMBOL_TOKEN, self.as_str())
    }
}

struct SymbolVisitor;

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
        // Non-EDN deserializers fall through to the inner string payload.
        String::deserialize(deserializer).map(EdnSymbol::new)
    }
}

impl<'de> Deserialize<'de> for EdnSymbol {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_newtype_struct(SYMBOL_TOKEN, SymbolVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::de::from_str;
    use crate::ser::to_string;

    #[test]
    fn test_edn_symbol_serialize_edn() {
        let sym = EdnSymbol::new("foo");
        let edn = to_string(&sym).unwrap();
        assert_eq!(edn, "foo");
    }

    #[test]
    fn test_edn_symbol_serialize_namespaced_edn() {
        let sym = EdnSymbol::new("ns/foo");
        let edn = to_string(&sym).unwrap();
        assert_eq!(edn, "ns/foo");
    }

    #[test]
    fn test_edn_symbol_serialize_json() {
        let sym = EdnSymbol::new("foo");
        let json = serde_json::to_string(&sym).unwrap();
        assert_eq!(json, r#""foo""#);
    }

    #[test]
    fn test_edn_symbol_in_vec() {
        let v = vec![EdnSymbol::new("a"), EdnSymbol::new("b")];
        let edn = to_string(&v).unwrap();
        assert_eq!(edn, "[a b]");
    }

    #[test]
    fn test_edn_symbol_in_struct() {
        #[derive(serde::Serialize)]
        struct Wrapper {
            tag: EdnSymbol,
            value: i64,
        }
        let w = Wrapper {
            tag: EdnSymbol::new("foo"),
            value: 5,
        };
        let edn = to_string(&w).unwrap();
        assert_eq!(edn, "{:tag foo :value 5}");
    }

    #[test]
    fn test_edn_symbol_deserialize_from_edn_symbol() {
        let sym: EdnSymbol = from_str("foo").unwrap();
        assert_eq!(sym.as_str(), "foo");
    }

    #[test]
    fn test_edn_symbol_deserialize_from_edn_string_rejected() {
        // Strict: string syntax is not a symbol.
        let result: Result<EdnSymbol, _> = from_str(r#""foo""#);
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[test]
    fn test_edn_symbol_deserialize_from_edn_keyword_rejected() {
        let result: Result<EdnSymbol, _> = from_str(":foo");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[test]
    fn test_edn_symbol_roundtrip() {
        let sym = EdnSymbol::new("ns/foo");
        let edn = to_string(&sym).unwrap();
        let parsed: EdnSymbol = from_str(&edn).unwrap();
        assert_eq!(parsed, sym);
    }
}
