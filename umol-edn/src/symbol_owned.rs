//! Owned EDN symbol type for use in struct fields.
//!
//! `EdnSymbol` is the owned, lifetime-free counterpart to the borrowing
//! [`Symbol<'a>`](crate::Symbol) used in the `Edn` data model. It exists
//! so that owned structs can carry a value that is statically known to be
//! a symbol (rather than a string or keyword) and so the EDN serde layer
//! can round-trip symbols losslessly.

use std::fmt;
use std::ops::Deref;

use crate::edn::{Edn, Symbol};
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};

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
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Symbol(s) => Ok(EdnSymbol(s.as_str().to_string())),
            other => Err(EdnError::TypeMismatch {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_string;

    #[test]
    fn test_edn_symbol_display() {
        assert_eq!(EdnSymbol::new("foo").to_string(), "foo");
    }

    #[test]
    fn test_edn_symbol_from_edn_symbol() {
        let edn = read_string("bar").unwrap();
        let sym = EdnSymbol::from_edn(&edn).unwrap();
        assert_eq!(sym.as_str(), "bar");
    }

    #[test]
    fn test_edn_symbol_from_edn_string_rejected() {
        let edn = read_string(r#""bar""#).unwrap();
        let result = EdnSymbol::from_edn(&edn);
        assert!(result.is_err());
    }

    #[test]
    fn test_edn_symbol_to_edn() {
        let sym = EdnSymbol::new("baz");
        let edn = sym.to_edn();
        if let Edn::Symbol(s) = &edn {
            assert_eq!(s.as_str(), "baz");
        } else {
            panic!("expected Edn::Symbol, got {:?}", edn);
        }
    }
}
