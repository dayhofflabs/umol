//! Owned EDN keyword type for use in struct fields.
//!
//! `EdnKeyword` is the owned, lifetime-free counterpart to the borrowing
//! [`Keyword<'a>`](crate::Keyword) used in the `Edn` data model. `FromEdn`
//! requires the input to be `Edn::Keyword`; reading a string into
//! `EdnKeyword` is a type error. `ToEdn` always emits `Edn::Keyword`.
//!
//! Under the `serde` feature, `Serialize` emits `:keyword` syntax via the
//! `KEYWORD_TOKEN` newtype-struct trick when routed through `EdnSerializer`;
//! non-EDN serializers see a transparent string. `Deserialize` accepts any
//! string input.

use std::fmt;
use std::ops::Deref;

#[cfg(feature = "serde")]
use serde::de::Deserializer;
#[cfg(feature = "serde")]
use serde::ser::Serializer;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::edn::{Edn, Keyword};
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};
#[cfg(feature = "serde")]
use crate::serde_tokens::KEYWORD_TOKEN;

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

#[cfg(feature = "serde")]
impl Serialize for EdnKeyword {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(KEYWORD_TOKEN, self.as_str())
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for EdnKeyword {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        String::deserialize(deserializer).map(EdnKeyword::new)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::read_string;
    #[cfg(feature = "serde")]
    use crate::{from_str, to_string};

    #[test]
    fn test_edn_keyword_display() {
        assert_eq!(EdnKeyword::new("foo").to_string(), ":foo");
    }

    #[test]
    fn test_edn_keyword_from_edn() {
        let edn = read_string(":bar").unwrap();
        let kw = EdnKeyword::from_edn(&edn).unwrap();
        assert_eq!(kw.as_str(), "bar");
    }

    #[test]
    fn test_edn_keyword_from_edn_error() {
        let edn = read_string(r#""bar""#).unwrap();
        assert!(EdnKeyword::from_edn(&edn).is_err());
    }

    #[test]
    fn test_edn_keyword_to_edn() {
        let keyword = EdnKeyword::new("baz");
        let edn = keyword.to_edn();
        let Edn::Keyword(k) = &edn else {
            panic!("expected Edn::Keyword, got {:?}", edn);
        };
        assert_eq!(k.as_str(), "baz");
    }

    #[rstest]
    #[case(EdnKeyword::new("a"), EdnKeyword::new("a"), true)]
    #[case(EdnKeyword::new("a"), EdnKeyword::new("b"), false)]
    fn test_edn_keyword_eq(#[case] a: EdnKeyword, #[case] b: EdnKeyword, #[case] equal: bool) {
        assert_eq!(a == b, equal);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_keyword_serialize() {
        assert_eq!(to_string(&EdnKeyword::new("foo")).unwrap(), ":foo");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_keyword_serialize_json() {
        let json = serde_json::to_string(&EdnKeyword::new("foo")).unwrap();
        assert_eq!(json, r#""foo""#);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_keyword_serialize_in_vec() {
        let v = vec![EdnKeyword::new("a"), EdnKeyword::new("b")];
        assert_eq!(to_string(&v).unwrap(), "[:a :b]");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_keyword_serialize_in_struct() {
        #[derive(serde::Serialize)]
        struct Wrapper {
            tag: EdnKeyword,
            value: i64,
        }
        let w = Wrapper {
            tag: EdnKeyword::new("foo"),
            value: 5,
        };
        assert_eq!(to_string(&w).unwrap(), "{:tag :foo :value 5}");
    }

    #[cfg(feature = "serde")]
    #[rstest]
    #[case(":foo", "foo")]
    #[case(r#""bar""#, "bar")]
    fn test_edn_keyword_deserialize(#[case] input: &str, #[case] expected: &str) {
        let kw: EdnKeyword = from_str(input).unwrap();
        assert_eq!(kw.as_str(), expected);
    }
}
