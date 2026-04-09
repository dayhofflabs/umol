//! Owned EDN arbitrary-precision integer wrapper for use in struct fields.
//!
//! `EdnBigInt` is a thin newtype around [`num_bigint::BigInt`] that opts into
//! EDN bigint syntax (`123N`) when serialized via the EDN serde layer, and
//! accepts only `Edn::BigInt` on the deserialization side. Non-EDN
//! serializers see the integer as a decimal string (JSON has no
//! arbitrary-precision integer literal).

use std::ops::{Deref, DerefMut};
use std::str::FromStr;

use num_bigint::BigInt;

#[cfg(feature = "serde")]
use std::fmt;

#[cfg(feature = "serde")]
use serde::de::{Deserialize, Deserializer, Error as DeError, Visitor};
#[cfg(feature = "serde")]
use serde::ser::{Serialize, Serializer};

use crate::edn::Edn;
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};
#[cfg(feature = "serde")]
use crate::serde_tokens::BIGINT_TOKEN;

/// An owned EDN arbitrary-precision integer.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct EdnBigInt(pub BigInt);

impl EdnBigInt {
    pub fn new(value: BigInt) -> Self {
        Self(value)
    }

    pub fn into_inner(self) -> BigInt {
        self.0
    }
}

impl Deref for EdnBigInt {
    type Target = BigInt;
    fn deref(&self) -> &BigInt {
        &self.0
    }
}

impl DerefMut for EdnBigInt {
    fn deref_mut(&mut self) -> &mut BigInt {
        &mut self.0
    }
}

impl From<BigInt> for EdnBigInt {
    fn from(v: BigInt) -> Self {
        Self(v)
    }
}

impl FromStr for EdnBigInt {
    type Err = <BigInt as FromStr>::Err;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        BigInt::from_str(s).map(Self)
    }
}

impl<'de> FromEdn<'de> for EdnBigInt {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::BigInt(n) => Ok(Self(n.clone())),
            other => Err(EdnError::TypeMismatch {
                expected: "bigint",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for EdnBigInt {
    fn to_edn(&self) -> Edn<'_> {
        Edn::BigInt(self.0.clone())
    }
}

#[cfg(feature = "serde")]
impl Serialize for EdnBigInt {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(BIGINT_TOKEN, &self.0.to_string())
    }
}

#[cfg(feature = "serde")]
struct BigIntVisitor;

#[cfg(feature = "serde")]
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
        String::deserialize(deserializer).and_then(|s| self.visit_string(s))
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for EdnBigInt {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_newtype_struct(BIGINT_TOKEN, BigIntVisitor)
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
    fn test_edn_bigint_from_edn() {
        let edn = read_string("123456789012345678901234567890N").unwrap();
        let n: EdnBigInt = EdnBigInt::from_edn(&edn).unwrap();
        assert_eq!(
            n.0,
            BigInt::from_str("123456789012345678901234567890").unwrap()
        );
    }

    #[test]
    fn test_edn_bigint_from_edn_error() {
        let edn = read_string("17").unwrap();
        let result: Result<EdnBigInt, _> = EdnBigInt::from_edn(&edn);
        assert!(result.is_err());
    }

    #[test]
    fn test_edn_bigint_to_edn() {
        let n = EdnBigInt::new(BigInt::from(17));
        assert_eq!(n.to_edn(), Edn::BigInt(BigInt::from(17)));
    }

    #[cfg(feature = "serde")]
    #[rstest]
    #[case("123456789012345678901234567890", "123456789012345678901234567890N")]
    #[case("-17", "-17N")]
    #[case("0", "0N")]
    fn test_edn_bigint_serialize(#[case] input: &str, #[case] expected: &str) {
        let n = EdnBigInt::new(BigInt::from_str(input).unwrap());
        assert_eq!(to_string(&n).unwrap(), expected);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_bigint_serialize_json() {
        let n = EdnBigInt::new(BigInt::from(17));
        let json = serde_json::to_string(&n).unwrap();
        assert_eq!(json, "\"17\"");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_bigint_deserialize() {
        let n: EdnBigInt = from_str("123456789012345678901234567890N").unwrap();
        assert_eq!(
            n.0,
            BigInt::from_str("123456789012345678901234567890").unwrap()
        );
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_bigint_deserialize_error() {
        let result: Result<EdnBigInt, _> = from_str("17");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_bigint_roundtrip() {
        let n = EdnBigInt::new(BigInt::from_str("999999999999999999999").unwrap());
        let edn = to_string(&n).unwrap();
        let parsed: EdnBigInt = from_str(&edn).unwrap();
        assert_eq!(parsed, n);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_bigint_serialize_in_struct() {
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
