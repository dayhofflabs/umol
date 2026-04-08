//! Serde Serialize/Deserialize impls for [`EdnKeyword`].
//!
//! `EdnKeyword` lives in [`crate::keyword_owned`] (always available); the
//! impls in this file are only compiled with the `serde` feature. When
//! serialized via `EdnSerializer`, an `EdnKeyword` emits `:keyword`
//! syntax via the `KEYWORD_TOKEN` newtype-struct trick. Non-EDN
//! serializers see a transparent string.

use serde::de::Deserializer;
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};

use crate::keyword_owned::EdnKeyword;

pub const KEYWORD_TOKEN: &str = "$edn::keyword";

impl Serialize for EdnKeyword {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(KEYWORD_TOKEN, self.as_str())
    }
}

impl<'de> Deserialize<'de> for EdnKeyword {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        String::deserialize(deserializer).map(EdnKeyword::new)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::de::from_str;
    use crate::ser::to_string;

    #[test]
    fn test_edn_keyword_serialize_edn() {
        let kw = EdnKeyword::new("foo");
        let edn = to_string(&kw).unwrap();
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
        let edn = to_string(&v).unwrap();
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
        let edn = to_string(&w).unwrap();
        assert_eq!(edn, "{:tag :foo :value 5}");
    }

    #[test]
    fn test_edn_keyword_nested_in_tuple() {
        let t = (EdnKeyword::new("x"), 10i64, EdnKeyword::new("y"));
        let edn = to_string(&t).unwrap();
        assert_eq!(edn, "[:x 10 :y]");
    }

    #[rstest]
    #[case(EdnKeyword::new("a"), EdnKeyword::new("a"), true)]
    #[case(EdnKeyword::new("a"), EdnKeyword::new("b"), false)]
    fn test_edn_keyword_eq(#[case] a: EdnKeyword, #[case] b: EdnKeyword, #[case] equal: bool) {
        assert_eq!(a == b, equal);
    }

    #[test]
    fn test_edn_keyword_deserialize_from_edn_keyword() {
        let kw: EdnKeyword = from_str(":foo").unwrap();
        assert_eq!(kw.as_str(), "foo");
    }

    #[test]
    fn test_edn_keyword_deserialize_from_edn_string() {
        let kw: EdnKeyword = from_str(r#""bar""#).unwrap();
        assert_eq!(kw.as_str(), "bar");
    }
}
