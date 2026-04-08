//! Owned EDN tagged literal wrapper for use in struct fields.
//!
//! `EdnTagged<T>` carries a tag string and an inner value of type `T`. When
//! serialized via the EDN serde layer, it emits `#tag value` syntax. Non-EDN
//! serializers see it as a 2-element tuple struct (e.g. `[tag, value]` in
//! JSON). On the deserialization side, only `Edn::Tagged` is accepted.

use std::borrow::Cow;

use crate::edn::Edn;
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};

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
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Tagged(tag, inner) => Ok(Self {
                tag: tag.to_string(),
                value: T::from_edn(inner)?,
            }),
            other => Err(EdnError::TypeMismatch {
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
    fn to_edn(&self) -> Edn<'_> {
        Edn::Tagged(
            Cow::Owned(self.tag.clone()),
            Box::new(self.value.to_edn().into_owned()),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_string;

    #[test]
    fn test_edn_tagged_from_edn_tagged() {
        let edn = read_string("#inst \"2026-04-08\"").unwrap();
        let tagged: EdnTagged<String> = EdnTagged::from_edn(&edn).unwrap();
        assert_eq!(tagged.tag, "inst");
        assert_eq!(tagged.value, "2026-04-08");
    }

    #[test]
    fn test_edn_tagged_from_edn_non_tagged_rejected() {
        let edn = read_string("\"2026-04-08\"").unwrap();
        let result: Result<EdnTagged<String>, _> = EdnTagged::from_edn(&edn);
        assert!(result.is_err());
    }

    #[test]
    fn test_edn_tagged_to_edn() {
        let tagged: EdnTagged<i64> = EdnTagged::new("score", 17);
        let edn = tagged.to_edn();
        if let Edn::Tagged(tag, inner) = &edn {
            assert_eq!(tag.as_ref(), "score");
            assert_eq!(**inner, Edn::Int(17));
        } else {
            panic!("expected Edn::Tagged, got {:?}", edn);
        }
    }
}
