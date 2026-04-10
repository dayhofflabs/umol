//! Owned EDN arbitrary-precision decimal wrapper.
//!
//! `EdnBigDecimal` is a newtype around [`bigdecimal::BigDecimal`] that adds
//! `Eq`. The upstream crate omits `Eq` despite `BigDecimal` having no
//! NaN-like values (reflexivity holds); its `Hash` impl normalizes trailing
//! zeros before hashing, consistent with `PartialEq`, so the `Eq` + `Hash`
//! contract is satisfied.
//!
//! Under the `serde` feature, `Serialize` emits `1.5M` syntax via the
//! `BIGDECIMAL_TOKEN` newtype-struct trick when routed through
//! `EdnSerializer`; non-EDN serializers see the value as a decimal string.
//! `Deserialize` via the EDN tree path accepts only `Edn::BigDecimal`.

use std::cmp::Ordering;
use std::fmt;
#[cfg(feature = "serde")]
use std::str::FromStr;

#[cfg(feature = "serde")]
use ::serde::de::{Deserialize, Deserializer, Error as DeError, Visitor};
#[cfg(feature = "serde")]
use ::serde::ser::{Serialize, Serializer};
use bigdecimal::BigDecimal;

#[cfg(feature = "serde")]
use crate::serde::BIGDECIMAL_TOKEN;

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct EdnBigDecimal(BigDecimal);

impl EdnBigDecimal {
    pub fn new(bd: BigDecimal) -> Self {
        Self(bd)
    }

    pub fn into_inner(self) -> BigDecimal {
        self.0
    }

    pub fn as_inner(&self) -> &BigDecimal {
        &self.0
    }
}

impl PartialOrd for EdnBigDecimal {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EdnBigDecimal {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl fmt::Display for EdnBigDecimal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(feature = "serde")]
impl Serialize for EdnBigDecimal {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(BIGDECIMAL_TOKEN, &self.as_inner().to_string())
    }
}

#[cfg(feature = "serde")]
struct BigDecimalVisitor;

#[cfg(feature = "serde")]
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
        String::deserialize(deserializer).and_then(|s| self.visit_string(s))
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for EdnBigDecimal {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_newtype_struct(BIGDECIMAL_TOKEN, BigDecimalVisitor)
    }
}

#[cfg(test)]
mod tests {
    use std::hash::{DefaultHasher, Hash, Hasher};
    use std::str::FromStr;

    #[cfg(feature = "serde")]
    use rstest::rstest;

    use super::*;
    #[cfg(feature = "serde")]
    use crate::serde::{from_str, to_string};

    /// Guards the assumption that `BigDecimal`'s `PartialEq` and `Hash` are
    /// consistent across representations of the same value. If this breaks on
    /// a `bigdecimal` upgrade, `EdnBigDecimal`'s `Eq` impl is unsound.
    #[test]
    fn test_edn_bigdecimal_eq() {
        let a = EdnBigDecimal::new(BigDecimal::from_str("1.0").unwrap());
        let b = EdnBigDecimal::new(BigDecimal::from_str("1.00").unwrap());
        let c = EdnBigDecimal::new(BigDecimal::from_str("1.000").unwrap());
        assert_eq!(a, b);
        assert_eq!(b, c);
        assert_eq!(a, a);
    }

    #[test]
    fn test_edn_bigdecimal_hash() {
        fn hash_bd(v: &EdnBigDecimal) -> u64 {
            let mut h = DefaultHasher::new();
            v.hash(&mut h);
            h.finish()
        }
        let a = EdnBigDecimal::new(BigDecimal::from_str("1.0").unwrap());
        let b = EdnBigDecimal::new(BigDecimal::from_str("1.00").unwrap());
        let c = EdnBigDecimal::new(BigDecimal::from_str("1.000").unwrap());
        assert_eq!(hash_bd(&a), hash_bd(&b));
        assert_eq!(hash_bd(&b), hash_bd(&c));
    }

    #[cfg(feature = "serde")]
    #[rstest]
    #[case::pi("3.14159265358979323846", "3.14159265358979323846M")]
    #[case::negative("-1.5", "-1.5M")]
    fn test_edn_bigdecimal_serialize(#[case] input: &str, #[case] expected: &str) {
        let d = EdnBigDecimal::new(BigDecimal::from_str(input).unwrap());
        assert_eq!(to_string(&d).unwrap(), expected);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_bigdecimal_serialize_json() {
        let d = EdnBigDecimal::new(BigDecimal::from_str("1.5").unwrap());
        let json = serde_json::to_string(&d).unwrap();
        assert_eq!(json, "\"1.5\"");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_bigdecimal_deserialize() {
        let d: EdnBigDecimal = from_str("3.14M").unwrap();
        assert_eq!(*d.as_inner(), BigDecimal::from_str("3.14").unwrap());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_bigdecimal_deserialize_error() {
        let result: Result<EdnBigDecimal, _> = from_str("3.14");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_bigdecimal_roundtrip() {
        let d = EdnBigDecimal::new(BigDecimal::from_str("999999.999999").unwrap());
        let edn = to_string(&d).unwrap();
        let parsed: EdnBigDecimal = from_str(&edn).unwrap();
        assert_eq!(parsed, d);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_bigdecimal_serialize_in_struct() {
        #[derive(::serde::Serialize, ::serde::Deserialize, Debug, PartialEq)]
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
