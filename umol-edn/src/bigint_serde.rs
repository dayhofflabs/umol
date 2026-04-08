//! Serde Serialize/Deserialize impls for [`EdnBigInt`].
//!
//! When serialized via the EDN serializers, `EdnBigInt` emits `123N` syntax
//! via the `BIGINT_TOKEN` newtype-struct trick. Non-EDN serializers see the
//! integer as a decimal string (JSON has no arbitrary-precision integer
//! literal). On the deserialization side, only `Edn::BigInt` is accepted.

use std::fmt;
use std::str::FromStr;

use num_bigint::BigInt;
use serde::de::{Deserialize, Deserializer, Error as DeError, Visitor};
use serde::ser::{Serialize, Serializer};

use crate::bigint_owned::EdnBigInt;
use crate::serde_tokens::BIGINT_TOKEN;

impl Serialize for EdnBigInt {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(BIGINT_TOKEN, &self.0.to_string())
    }
}

struct BigIntVisitor;

impl<'de> Visitor<'de> for BigIntVisitor {
    type Value = EdnBigInt;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("an EDN bigint (decimal string)")
    }

    fn visit_str<E: DeError>(self, v: &str) -> Result<EdnBigInt, E> {
        BigInt::from_str(v)
            .map(EdnBigInt)
            .map_err(|e| E::custom(format!("invalid bigint {v:?}: {e}")))
    }

    fn visit_string<E: DeError>(self, v: String) -> Result<EdnBigInt, E> {
        self.visit_str(&v)
    }

    fn visit_borrowed_str<E: DeError>(self, v: &'de str) -> Result<EdnBigInt, E> {
        self.visit_str(v)
    }

    fn visit_newtype_struct<D: Deserializer<'de>>(
        self,
        deserializer: D,
    ) -> Result<EdnBigInt, D::Error> {
        // Non-EDN deserializers fall through here: read the inner string.
        String::deserialize(deserializer).and_then(|s| self.visit_string(s))
    }
}

impl<'de> Deserialize<'de> for EdnBigInt {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_newtype_struct(BIGINT_TOKEN, BigIntVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::de::from_str;
    use crate::ser::to_string;

    #[test]
    fn test_edn_bigint_serialize_edn() {
        let n = EdnBigInt::new(BigInt::from_str("123456789012345678901234567890").unwrap());
        let edn = to_string(&n).unwrap();
        assert_eq!(edn, "123456789012345678901234567890N");
    }

    #[test]
    fn test_edn_bigint_serialize_negative() {
        let n = EdnBigInt::new(BigInt::from(-17));
        let edn = to_string(&n).unwrap();
        assert_eq!(edn, "-17N");
    }

    #[test]
    fn test_edn_bigint_serialize_json_falls_back_to_string() {
        let n = EdnBigInt::new(BigInt::from(17));
        let json = serde_json::to_string(&n).unwrap();
        assert_eq!(json, "\"17\"");
    }

    #[test]
    fn test_edn_bigint_deserialize_from_edn_bigint() {
        let n: EdnBigInt = from_str("123456789012345678901234567890N").unwrap();
        assert_eq!(
            n.0,
            BigInt::from_str("123456789012345678901234567890").unwrap()
        );
    }

    #[test]
    fn test_edn_bigint_deserialize_from_int_rejected() {
        let result: Result<EdnBigInt, _> = from_str("17");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[test]
    fn test_edn_bigint_roundtrip() {
        let n = EdnBigInt::new(BigInt::from_str("999999999999999999999").unwrap());
        let edn = to_string(&n).unwrap();
        let parsed: EdnBigInt = from_str(&edn).unwrap();
        assert_eq!(parsed, n);
    }

    #[test]
    fn test_edn_bigint_in_struct() {
        #[derive(serde::Serialize, serde::Deserialize, Debug, PartialEq)]
        struct Wrapper {
            count: EdnBigInt,
        }
        let w = Wrapper {
            count: EdnBigInt::new(BigInt::from_str("99999999999999999999").unwrap()),
        };
        let edn = to_string(&w).unwrap();
        assert_eq!(edn, "{:count 99999999999999999999N}");
        let parsed: Wrapper = from_str(&edn).unwrap();
        assert_eq!(parsed, w);
    }
}
