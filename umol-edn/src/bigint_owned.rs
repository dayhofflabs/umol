//! Owned EDN arbitrary-precision integer wrapper for use in struct fields.
//!
//! `EdnBigInt` is a thin newtype around [`num_bigint::BigInt`] that opts into
//! EDN bigint syntax (`123N`) when serialized via the EDN serde layer, and
//! accepts only `Edn::BigInt` on the deserialization side. Non-EDN
//! serializers see it as a string (JSON has no arbitrary-precision integer
//! literal).

use std::ops::{Deref, DerefMut};
use std::str::FromStr;

use num_bigint::BigInt;

use crate::edn::Edn;
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_string;

    #[test]
    fn test_edn_bigint_from_edn_bigint() {
        let edn = read_string("123456789012345678901234567890N").unwrap();
        let n: EdnBigInt = EdnBigInt::from_edn(&edn).unwrap();
        assert_eq!(
            n.0,
            BigInt::from_str("123456789012345678901234567890").unwrap()
        );
    }

    #[test]
    fn test_edn_bigint_from_edn_int_rejected() {
        let edn = read_string("17").unwrap();
        let result: Result<EdnBigInt, _> = EdnBigInt::from_edn(&edn);
        assert!(result.is_err());
    }

    #[test]
    fn test_edn_bigint_to_edn() {
        let n = EdnBigInt::new(BigInt::from(17));
        let edn = n.to_edn();
        assert_eq!(edn, Edn::BigInt(BigInt::from(17)));
    }
}
