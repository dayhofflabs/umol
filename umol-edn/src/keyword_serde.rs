//! General-purpose EDN keyword type for serde roundtrips.
//!
//! `EdnKeyword` preserves the keyword/string distinction through serde
//! boundaries, modeled after `serde_json::RawValue`. When serialized via
//! `EdnSerializer`, it emits `:keyword` syntax. Non-EDN serializers see
//! a transparent string.

use std::fmt;
use std::ops::Deref;

use serde::de::Deserializer;
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};

pub const KEYWORD_TOKEN: &str = "$edn::keyword";

/// A string value that serializes as an EDN keyword (`:name`) rather than
/// a quoted string (`"name"`) when used with `EdnSerializer`.
///
/// Non-EDN serializers treat it as a plain string.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdnKeyword(String);

impl EdnKeyword {
    pub fn new(s: impl Into<String>) -> Self {
        Self(s.into())
    }

    pub fn into_inner(self) -> String {
        self.0
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

impl Serialize for EdnKeyword {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(KEYWORD_TOKEN, &self.0)
    }
}

impl<'de> Deserialize<'de> for EdnKeyword {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        String::deserialize(deserializer).map(EdnKeyword)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn test_edn_keyword_display() {
        assert_eq!(EdnKeyword::new("foo").to_string(), ":foo");
    }

    #[test]
    fn test_edn_keyword_deref() {
        let kw = EdnKeyword::new("bar");
        assert_eq!(&*kw, "bar");
    }

    #[rstest]
    #[case(EdnKeyword::new("a"), EdnKeyword::new("a"), true)]
    #[case(EdnKeyword::new("a"), EdnKeyword::new("b"), false)]
    fn test_edn_keyword_eq(#[case] a: EdnKeyword, #[case] b: EdnKeyword, #[case] equal: bool) {
        assert_eq!(a == b, equal);
    }

    #[test]
    fn test_edn_keyword_serialize_edn() {
        let kw = EdnKeyword::new("foo");
        let edn = crate::to_string(&kw).unwrap();
        assert_eq!(edn, ":foo");
    }

    #[test]
    fn test_edn_keyword_serialize_json() {
        let kw = EdnKeyword::new("foo");
        let json = serde_json::to_string(&kw).unwrap();
        assert_eq!(json, r#""foo""#);
    }

    #[test]
    fn test_edn_keyword_in_vec() {
        let v = vec![EdnKeyword::new("a"), EdnKeyword::new("b")];
        let edn = crate::to_string(&v).unwrap();
        assert_eq!(edn, "[:a :b]");
    }

    #[test]
    fn test_edn_keyword_in_struct() {
        #[derive(serde::Serialize)]
        struct Wrapper {
            tag: EdnKeyword,
            value: i64,
        }
        let w = Wrapper {
            tag: EdnKeyword::new("foo"),
            value: 5,
        };
        let edn = crate::to_string(&w).unwrap();
        assert_eq!(edn, "{:tag :foo :value 5}");
    }

    #[test]
    fn test_edn_keyword_nested_in_tuple() {
        let t = (EdnKeyword::new("x"), 10i64, EdnKeyword::new("y"));
        let edn = crate::to_string(&t).unwrap();
        assert_eq!(edn, "[:x 10 :y]");
    }

    #[test]
    fn test_edn_keyword_deserialize_from_edn_keyword() {
        let kw: EdnKeyword = crate::from_str(":foo").unwrap();
        assert_eq!(&*kw, "foo");
    }

    #[test]
    fn test_edn_keyword_deserialize_from_edn_string() {
        let kw: EdnKeyword = crate::from_str(r#""bar""#).unwrap();
        assert_eq!(&*kw, "bar");
    }
}
