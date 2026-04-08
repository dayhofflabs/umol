//! Owned EDN keyword type for use in struct fields.
//!
//! `EdnKeyword` is the owned, lifetime-free counterpart to the borrowing
//! [`Keyword<'a>`](crate::Keyword) used in the `Edn` data model. It is
//! intended for use as a field type on owned structs that need to
//! preserve the keyword/string distinction at the type level.
//!
//! `FromEdn` requires the input to be `Edn::Keyword`; reading a string
//! into `EdnKeyword` is a type error. `ToEdn` always emits `Edn::Keyword`.

use std::fmt;
use std::ops::Deref;

use crate::edn::{Edn, Keyword};
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};

/// An owned EDN keyword (`:name`). Always serializes as a keyword in EDN.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdnKeyword(String);

impl EdnKeyword {
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

impl Deref for EdnKeyword {
    type Target = str;
    fn deref(&self) -> &str {
        &self.0
    }
}

impl AsRef<str> for EdnKeyword {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for EdnKeyword {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, ":{}", self.0)
    }
}

impl From<String> for EdnKeyword {
    fn from(s: String) -> Self {
        Self(s)
    }
}

impl From<&str> for EdnKeyword {
    fn from(s: &str) -> Self {
        Self(s.to_string())
    }
}

impl<'de> FromEdn<'de> for EdnKeyword {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Keyword(k) => Ok(EdnKeyword(k.as_str().to_string())),
            other => Err(EdnError::TypeMismatch {
                expected: "keyword",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for EdnKeyword {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Keyword(Keyword::new(&self.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_string;

    #[test]
    fn test_edn_keyword_display() {
        assert_eq!(EdnKeyword::new("foo").to_string(), ":foo");
    }

    #[test]
    fn test_edn_keyword_from_edn_keyword() {
        let edn = read_string(":bar").unwrap();
        let kw = EdnKeyword::from_edn(&edn).unwrap();
        assert_eq!(kw.as_str(), "bar");
    }

    #[test]
    fn test_edn_keyword_from_edn_string_rejected() {
        let edn = read_string(r#""bar""#).unwrap();
        let result = EdnKeyword::from_edn(&edn);
        assert!(result.is_err());
    }

    #[test]
    fn test_edn_keyword_to_edn() {
        let kw = EdnKeyword::new("baz");
        let edn = kw.to_edn();
        if let Edn::Keyword(k) = &edn {
            assert_eq!(k.as_str(), "baz");
        } else {
            panic!("expected Edn::Keyword, got {:?}", edn);
        }
    }
}
