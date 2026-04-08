//! Serde Serialize/Deserialize impls for [`EdnBigDecimal`].
//!
//! When serialized via the EDN serializers, `EdnBigDecimal` emits `1.5M`
//! syntax via the `BIGDECIMAL_TOKEN` newtype-struct trick. Non-EDN
//! serializers see the value as a decimal string. On the deserialization
//! side, only `Edn::BigDecimal` is accepted.

use std::fmt;
use std::str::FromStr;

use bigdecimal::BigDecimal;
use serde::de::{Deserialize, Deserializer, Error as DeError, Visitor};
use serde::ser::{Serialize, Serializer};

use crate::edn::EdnBigDecimal;
use crate::serde_tokens::BIGDECIMAL_TOKEN;

impl Serialize for EdnBigDecimal {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(BIGDECIMAL_TOKEN, &self.as_inner().to_string())
    }
}

struct BigDecimalVisitor;

impl<'de> Visitor<'de> for BigDecimalVisitor {
    type Value = EdnBigDecimal;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("an EDN bigdecimal (decimal string)")
    }

    fn visit_str<E: DeError>(self, v: &str) -> Result<EdnBigDecimal, E> {
        BigDecimal::from_str(v)
            .map(EdnBigDecimal::new)
            .map_err(|e| E::custom(format!("invalid bigdecimal {v:?}: {e}")))
    }

    fn visit_string<E: DeError>(self, v: String) -> Result<EdnBigDecimal, E> {
        self.visit_str(&v)
    }

    fn visit_borrowed_str<E: DeError>(self, v: &'de str) -> Result<EdnBigDecimal, E> {
        self.visit_str(v)
    }

    fn visit_newtype_struct<D: Deserializer<'de>>(
        self,
        deserializer: D,
    ) -> Result<EdnBigDecimal, D::Error> {
        // Non-EDN deserializers fall through here: read the inner string.
        String::deserialize(deserializer).and_then(|s| self.visit_string(s))
    }
}

impl<'de> Deserialize<'de> for EdnBigDecimal {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_newtype_struct(BIGDECIMAL_TOKEN, BigDecimalVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::de::from_str;
    use crate::ser::to_string;

    #[test]
    fn test_edn_bigdecimal_serialize_edn() {
        let d = EdnBigDecimal::new(BigDecimal::from_str("3.14159265358979323846").unwrap());
        let edn = to_string(&d).unwrap();
        assert_eq!(edn, "3.14159265358979323846M");
    }

    #[test]
    fn test_edn_bigdecimal_serialize_negative() {
        let d = EdnBigDecimal::new(BigDecimal::from_str("-1.5").unwrap());
        let edn = to_string(&d).unwrap();
        assert_eq!(edn, "-1.5M");
    }

    #[test]
    fn test_edn_bigdecimal_serialize_json_falls_back_to_string() {
        let d = EdnBigDecimal::new(BigDecimal::from_str("1.5").unwrap());
        let json = serde_json::to_string(&d).unwrap();
        assert_eq!(json, "\"1.5\"");
    }

    #[test]
    fn test_edn_bigdecimal_deserialize_from_edn_bigdecimal() {
        let d: EdnBigDecimal = from_str("3.14M").unwrap();
        assert_eq!(*d.as_inner(), BigDecimal::from_str("3.14").unwrap());
    }

    #[test]
    fn test_edn_bigdecimal_deserialize_from_float_rejected() {
        let result: Result<EdnBigDecimal, _> = from_str("3.14");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[test]
    fn test_edn_bigdecimal_roundtrip() {
        let d = EdnBigDecimal::new(BigDecimal::from_str("999999.999999").unwrap());
        let edn = to_string(&d).unwrap();
        let parsed: EdnBigDecimal = from_str(&edn).unwrap();
        assert_eq!(parsed, d);
    }

    #[test]
    fn test_edn_bigdecimal_in_struct() {
        #[derive(serde::Serialize, serde::Deserialize, Debug, PartialEq)]
        struct Wrapper {
            ratio: EdnBigDecimal,
        }
        let w = Wrapper {
            ratio: EdnBigDecimal::new(BigDecimal::from_str("0.333333333333").unwrap()),
        };
        let edn = to_string(&w).unwrap();
        assert_eq!(edn, "{:ratio 0.333333333333M}");
        let parsed: Wrapper = from_str(&edn).unwrap();
        assert_eq!(parsed, w);
    }
}
