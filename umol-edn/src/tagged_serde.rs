//! Serde Serialize/Deserialize impls for [`EdnTagged`].
//!
//! When serialized via the EDN serializers, `EdnTagged<T>` emits `#tag value`
//! syntax via the `TAGGED_TOKEN` tuple-struct trick. Non-EDN serializers see
//! a transparent 2-element tuple struct (`[tag, value]` in JSON). On the
//! deserialization side, only `Edn::Tagged` is accepted.

use std::fmt;
use std::marker::PhantomData;

use serde::de::{Deserialize, Deserializer, Error as DeError, SeqAccess, Visitor};
use serde::ser::{Serialize, SerializeTupleStruct, Serializer};

use crate::serde_tokens::TAGGED_TOKEN;
use crate::tagged_owned::EdnTagged;

impl<T: Serialize> Serialize for EdnTagged<T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut t = serializer.serialize_tuple_struct(TAGGED_TOKEN, 2)?;
        t.serialize_field(self.tag.as_str())?;
        t.serialize_field(&self.value)?;
        t.end()
    }
}

struct TaggedVisitor<T> {
    _marker: PhantomData<T>,
}

impl<T> Default for TaggedVisitor<T> {
    fn default() -> Self {
        Self {
            _marker: PhantomData,
        }
    }
}

impl<'de, T: Deserialize<'de>> Visitor<'de> for TaggedVisitor<T> {
    type Value = EdnTagged<T>;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("an EDN tagged literal")
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<EdnTagged<T>, A::Error> {
        let tag: String = seq
            .next_element()?
            .ok_or_else(|| A::Error::custom("missing tag in EDN tagged literal"))?;
        let value: T = seq
            .next_element()?
            .ok_or_else(|| A::Error::custom("missing value in EDN tagged literal"))?;
        Ok(EdnTagged { tag, value })
    }
}

impl<'de, T: Deserialize<'de>> Deserialize<'de> for EdnTagged<T> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_tuple_struct(TAGGED_TOKEN, 2, TaggedVisitor::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::ParseConfig;
    use crate::de::{from_str, from_str_with};
    use crate::ser::to_string;

    /// Deserialize under a config that preserves unknown tags — required
    /// when `EdnTagged<T>` carries a caller-chosen tag name that has no
    /// registered reader.
    fn from_str_permissive<'a, T: serde::Deserialize<'a>>(s: &'a str) -> T {
        let mut config = ParseConfig::default();
        config.allow_unknown_tags = true;
        from_str_with(s, &config).unwrap()
    }

    #[test]
    fn test_edn_tagged_serialize_edn() {
        let tagged: EdnTagged<i64> = EdnTagged::new("score", 17);
        let edn = to_string(&tagged).unwrap();
        assert_eq!(edn, "#score 17");
    }

    #[test]
    fn test_edn_tagged_serialize_string_value() {
        let tagged: EdnTagged<String> = EdnTagged::new("inst", "2026-04-08".to_string());
        let edn = to_string(&tagged).unwrap();
        assert_eq!(edn, "#inst \"2026-04-08\"");
    }

    #[test]
    fn test_edn_tagged_serialize_json_falls_back_to_tuple() {
        let tagged: EdnTagged<i64> = EdnTagged::new("score", 17);
        let json = serde_json::to_string(&tagged).unwrap();
        assert_eq!(json, "[\"score\",17]");
    }

    #[test]
    fn test_edn_tagged_deserialize_from_edn_tagged() {
        let tagged: EdnTagged<i64> = from_str_permissive("#score 17");
        assert_eq!(tagged.tag, "score");
        assert_eq!(tagged.value, 17);
    }

    #[test]
    fn test_edn_tagged_deserialize_string_value() {
        let tagged: EdnTagged<String> = from_str("#inst \"2026-04-08\"").unwrap();
        assert_eq!(tagged.tag, "inst");
        assert_eq!(tagged.value, "2026-04-08");
    }

    #[test]
    fn test_edn_tagged_deserialize_from_non_tagged_rejected() {
        let result: Result<EdnTagged<i64>, _> = from_str("17");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[test]
    fn test_edn_tagged_roundtrip() {
        let tagged: EdnTagged<i64> = EdnTagged::new("score", 17);
        let edn = to_string(&tagged).unwrap();
        let parsed: EdnTagged<i64> = from_str_permissive(&edn);
        assert_eq!(parsed, tagged);
    }

    #[test]
    fn test_edn_tagged_in_struct() {
        #[derive(serde::Serialize, serde::Deserialize, Debug, PartialEq)]
        struct Wrapper {
            stamp: EdnTagged<String>,
        }
        let w = Wrapper {
            stamp: EdnTagged::new("inst", "2026-04-08".to_string()),
        };
        let edn = to_string(&w).unwrap();
        assert_eq!(edn, "{:stamp #inst \"2026-04-08\"}");
        let parsed: Wrapper = from_str(&edn).unwrap();
        assert_eq!(parsed, w);
    }

    #[test]
    fn test_edn_tagged_nested_value() {
        let tagged: EdnTagged<Vec<i64>> = EdnTagged::new("scores", vec![1, 2, 3]);
        let edn = to_string(&tagged).unwrap();
        assert_eq!(edn, "#scores [1 2 3]");
        let parsed: EdnTagged<Vec<i64>> = from_str_permissive(&edn);
        assert_eq!(parsed, tagged);
    }
}
