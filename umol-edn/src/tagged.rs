//! Owned EDN tagged literal wrapper for use in struct fields.
//!
//! `EdnTagged<T>` carries a tag string and an inner value of type `T`. When
//! serialized via the EDN serde layer, it emits `#tag value` syntax via the
//! `TAGGED_TOKEN` tuple-struct trick. Non-EDN serializers see a transparent
//! 2-element tuple struct (`[tag, value]` in JSON). On the deserialization
//! side, only `Edn::Tagged` is accepted.

use std::borrow::Cow;
#[cfg(feature = "serde")]
use std::fmt;
#[cfg(feature = "serde")]
use std::marker::PhantomData;

#[cfg(feature = "serde")]
use ::serde::de::{Deserialize, Deserializer, Error as SerdeDeError, SeqAccess, Visitor};
#[cfg(feature = "serde")]
use ::serde::ser::{Serialize, SerializeTupleStruct, Serializer};

use crate::edn::Edn;
use crate::error::DeError;
#[cfg(feature = "serde")]
use crate::serde::TAGGED_TOKEN;
use crate::traits::{FromEdn, ToEdn};

/// An owned EDN tagged literal `#tag value`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EdnTagged<T> {
    pub tag: String,
    pub value: T,
}

impl<T> EdnTagged<T> {
    pub fn new(tag: impl Into<String>, value: T) -> Self {
        Self {
            tag: tag.into(),
            value,
        }
    }
}

impl<'de, T> FromEdn<'de> for EdnTagged<T>
where
    T: FromEdn<'de>,
{
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Tagged(tag, inner) => Ok(Self {
                tag: tag.to_string(),
                value: T::from_edn(inner)?,
            }),
            other => Err(DeError::TypeMismatch {
                expected: "tagged",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<T> ToEdn for EdnTagged<T>
where
    T: ToEdn,
{
    fn to_edn(&self) -> Edn<'static> {
        Edn::Tagged(
            Cow::Owned(self.tag.clone()),
            Box::new(self.value.to_edn()),
        )
    }
}

#[cfg(feature = "serde")]
impl<T: Serialize> Serialize for EdnTagged<T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut t = serializer.serialize_tuple_struct(TAGGED_TOKEN, 2)?;
        t.serialize_field(self.tag.as_str())?;
        t.serialize_field(&self.value)?;
        t.end()
    }
}

#[cfg(feature = "serde")]
struct TaggedVisitor<T> {
    _marker: PhantomData<T>,
}

#[cfg(feature = "serde")]
impl<T> Default for TaggedVisitor<T> {
    fn default() -> Self {
        Self {
            _marker: PhantomData,
        }
    }
}

#[cfg(feature = "serde")]
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

#[cfg(feature = "serde")]
impl<'de, T: Deserialize<'de>> Deserialize<'de> for EdnTagged<T> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_tuple_struct(TAGGED_TOKEN, 2, TaggedVisitor::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "serde")]
    use crate::serde::{from_str, from_str_with, to_string};
    use crate::{read_string, ParseConfig};

    /// Unknown-tag config needed for `EdnTagged<T>` with caller-chosen tag
    /// names that have no registered reader.
    #[cfg(feature = "serde")]
    fn from_str_permissive<'a, T: ::serde::Deserialize<'a>>(s: &'a str) -> T {
        let config = ParseConfig {
            allow_unknown_tags: true,
            ..Default::default()
        };
        from_str_with(s, &config).unwrap()
    }

    #[test]
    fn test_edn_tagged_from_edn() {
        let edn = read_string("#inst \"2026-04-08T00:00:00Z\"").unwrap();
        let tagged: EdnTagged<String> = EdnTagged::from_edn(&edn).unwrap();
        assert_eq!(tagged.tag, "inst");
        assert_eq!(tagged.value, "2026-04-08T00:00:00Z");
    }

    #[test]
    fn test_edn_tagged_from_edn_error() {
        let edn = read_string("\"2026-04-08T00:00:00Z\"").unwrap();
        let result: Result<EdnTagged<String>, _> = EdnTagged::from_edn(&edn);
        assert!(result.is_err());
    }

    #[test]
    fn test_edn_tagged_to_edn() {
        let tagged: EdnTagged<i64> = EdnTagged::new("score", 17);
        let edn = tagged.to_edn();
        let Edn::Tagged(tag, inner) = &edn else {
            panic!("expected Edn::Tagged, got {:?}", edn);
        };
        assert_eq!(tag.as_ref(), "score");
        assert_eq!(**inner, Edn::Int(17));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_tagged_serialize() {
        let tagged: EdnTagged<i64> = EdnTagged::new("score", 17);
        assert_eq!(to_string(&tagged).unwrap(), "#score 17");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_tagged_serialize_string_value() {
        let tagged: EdnTagged<String> = EdnTagged::new("inst", "2026-04-08T00:00:00Z".to_string());
        assert_eq!(
            to_string(&tagged).unwrap(),
            "#inst \"2026-04-08T00:00:00Z\""
        );
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_tagged_serialize_json() {
        let tagged: EdnTagged<i64> = EdnTagged::new("score", 17);
        let json = serde_json::to_string(&tagged).unwrap();
        assert_eq!(json, "[\"score\",17]");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_tagged_deserialize() {
        let tagged: EdnTagged<i64> = from_str_permissive("#score 17");
        assert_eq!(tagged.tag, "score");
        assert_eq!(tagged.value, 17);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_tagged_deserialize_string_value() {
        let tagged: EdnTagged<String> = from_str("#inst \"2026-04-08T00:00:00Z\"").unwrap();
        assert_eq!(tagged.tag, "inst");
        assert_eq!(tagged.value, "2026-04-08T00:00:00Z");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_tagged_deserialize_error() {
        let result: Result<EdnTagged<i64>, _> = from_str("17");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_tagged_roundtrip() {
        let tagged: EdnTagged<i64> = EdnTagged::new("score", 17);
        let edn = to_string(&tagged).unwrap();
        let parsed: EdnTagged<i64> = from_str_permissive(&edn);
        assert_eq!(parsed, tagged);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_tagged_serialize_in_struct() {
        #[derive(::serde::Serialize, ::serde::Deserialize, Debug, PartialEq)]
        struct Wrapper {
            stamp: EdnTagged<String>,
        }
        let w = Wrapper {
            stamp: EdnTagged::new("inst", "2026-04-08T00:00:00Z".to_string()),
        };
        let edn = to_string(&w).unwrap();
        assert_eq!(edn, "{:stamp #inst \"2026-04-08T00:00:00Z\"}");
        let parsed: Wrapper = from_str(&edn).unwrap();
        assert_eq!(parsed, w);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_tagged_nested_value() {
        let tagged: EdnTagged<Vec<i64>> = EdnTagged::new("scores", vec![1, 2, 3]);
        let edn = to_string(&tagged).unwrap();
        assert_eq!(edn, "#scores [1 2 3]");
        let parsed: EdnTagged<Vec<i64>> = from_str_permissive(&edn);
        assert_eq!(parsed, tagged);
    }
}
