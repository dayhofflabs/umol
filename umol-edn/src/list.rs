//! Owned EDN list wrapper for use in struct fields.
//!
//! `EdnList<T>` is a thin wrapper around `Vec<T>` that opts into EDN list
//! syntax (`(...)`) rather than vector syntax (`[...]`) when serialized via
//! the EDN serde layer, and accepts only `Edn::List` on the deserialization
//! side. Non-EDN serializers see it as an ordinary sequence.
//!
//! Under the `serde` feature, `Serialize` emits `(...)` via the `LIST_TOKEN`
//! newtype-struct trick when routed through `EdnSerializer`; other serializers
//! see a transparent sequence. `Deserialize` via the EDN tree path accepts
//! only `Edn::List`; `Edn::Vector` and `Edn::Set` are rejected.

#[cfg(feature = "serde")]
use std::fmt;
#[cfg(feature = "serde")]
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};

#[cfg(feature = "serde")]
use ::serde::de::{Deserialize, DeserializeSeed, Deserializer, SeqAccess, Visitor};
#[cfg(feature = "serde")]
use ::serde::ser::{Serialize, SerializeSeq, Serializer};

use crate::collections::EdnSeq;
use crate::edn::Edn;
use crate::error::DeError;
#[cfg(feature = "serde")]
use crate::serde::LIST_TOKEN;
use crate::traits::{FromEdn, ToEdn};

/// An owned EDN list. Always serializes as `(...)` in EDN.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EdnList<T>(pub Vec<T>);

impl<T> EdnList<T> {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn into_inner(self) -> Vec<T> {
        self.0
    }
}

impl<T> Deref for EdnList<T> {
    type Target = Vec<T>;
    fn deref(&self) -> &Vec<T> {
        &self.0
    }
}

impl<T> DerefMut for EdnList<T> {
    fn deref_mut(&mut self) -> &mut Vec<T> {
        &mut self.0
    }
}

impl<T> From<Vec<T>> for EdnList<T> {
    fn from(v: Vec<T>) -> Self {
        Self(v)
    }
}

impl<T> FromIterator<T> for EdnList<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<'de, T> FromEdn<'de> for EdnList<T>
where
    T: FromEdn<'de>,
{
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::List(v) => {
                let mut out = Vec::with_capacity(v.len());
                for item in v.iter() {
                    out.push(T::from_edn(item)?);
                }
                Ok(EdnList(out))
            }
            other => Err(DeError::TypeMismatch {
                expected: "list",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<T> ToEdn for EdnList<T>
where
    T: ToEdn,
{
    fn to_edn(&self) -> Edn<'_> {
        let items: Vec<Edn<'static>> = self.0.iter().map(|v| v.to_edn().into_owned()).collect();
        Edn::List(EdnSeq::from(items))
    }
}

#[cfg(feature = "serde")]
impl<T: Serialize> Serialize for EdnList<T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(LIST_TOKEN, &ListPayload(&self.0))
    }
}

#[cfg(feature = "serde")]
struct ListPayload<'a, T>(&'a Vec<T>);

#[cfg(feature = "serde")]
impl<T: Serialize> Serialize for ListPayload<'_, T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;
        for v in self.0 {
            seq.serialize_element(v)?;
        }
        seq.end()
    }
}

#[cfg(feature = "serde")]
struct ListVisitor<T> {
    _marker: PhantomData<T>,
}

#[cfg(feature = "serde")]
impl<'de, T: Deserialize<'de>> Visitor<'de> for ListVisitor<T> {
    type Value = EdnList<T>;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("an EDN list")
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<EdnList<T>, A::Error> {
        let mut out = Vec::with_capacity(seq.size_hint().unwrap_or(0));
        while let Some(v) = seq.next_element_seed(PhantomDataSeed::<T>::default())? {
            out.push(v);
        }
        Ok(EdnList(out))
    }

    fn visit_newtype_struct<D: Deserializer<'de>>(
        self,
        deserializer: D,
    ) -> Result<EdnList<T>, D::Error> {
        Vec::<T>::deserialize(deserializer).map(EdnList)
    }
}

#[cfg(feature = "serde")]
struct PhantomDataSeed<T> {
    _marker: PhantomData<T>,
}

#[cfg(feature = "serde")]
impl<T> Default for PhantomDataSeed<T> {
    fn default() -> Self {
        Self {
            _marker: PhantomData,
        }
    }
}

#[cfg(feature = "serde")]
impl<'de, T: Deserialize<'de>> DeserializeSeed<'de> for PhantomDataSeed<T> {
    type Value = T;

    fn deserialize<D: Deserializer<'de>>(self, deserializer: D) -> Result<T, D::Error> {
        T::deserialize(deserializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T: Deserialize<'de>> Deserialize<'de> for EdnList<T> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_newtype_struct(
            LIST_TOKEN,
            ListVisitor {
                _marker: PhantomData,
            },
        )
    }
}

#[cfg(test)]
mod tests {
    #[cfg(feature = "serde")]
    use rstest::rstest;

    use super::*;
    use crate::read_string;
    #[cfg(feature = "serde")]
    use crate::serde::{from_str, to_string};

    #[test]
    fn test_edn_list_from_edn() {
        let edn = read_string("(1 2 3)").unwrap();
        let list: EdnList<i64> = EdnList::from_edn(&edn).unwrap();
        assert_eq!(list.0, vec![1, 2, 3]);
    }

    #[test]
    fn test_edn_list_from_edn_error() {
        let edn = read_string("[1 2 3]").unwrap();
        let result: Result<EdnList<i64>, _> = EdnList::from_edn(&edn);
        assert!(result.is_err());
    }

    #[test]
    fn test_edn_list_to_edn() {
        let list: EdnList<i64> = vec![1i64, 2, 3].into();
        let edn = list.to_edn();
        let Edn::List(v) = &edn else {
            panic!("expected Edn::List, got {:?}", edn);
        };
        assert_eq!(v.len(), 3);
    }

    #[cfg(feature = "serde")]
    #[rstest]
    #[case(vec![1, 2, 3], "(1 2 3)")]
    #[case(vec![], "()")]
    fn test_edn_list_serialize(#[case] input: Vec<i64>, #[case] expected: &str) {
        let list: EdnList<i64> = input.into();
        assert_eq!(to_string(&list).unwrap(), expected);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_list_serialize_json() {
        let list: EdnList<i64> = vec![1, 2, 3].into();
        let json = serde_json::to_string(&list).unwrap();
        assert_eq!(json, "[1,2,3]");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_list_deserialize() {
        let list: EdnList<i64> = from_str("(1 2 3)").unwrap();
        assert_eq!(list.0, vec![1, 2, 3]);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_list_deserialize_error() {
        let result: Result<EdnList<i64>, _> = from_str("[1 2 3]");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_list_roundtrip() {
        let list: EdnList<i64> = vec![1, 2, 3].into();
        let edn = to_string(&list).unwrap();
        let parsed: EdnList<i64> = from_str(&edn).unwrap();
        assert_eq!(parsed, list);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_list_serialize_in_struct() {
        #[derive(::serde::Serialize, ::serde::Deserialize, Debug, PartialEq)]
        struct Wrapper {
            steps: EdnList<i64>,
        }
        let w = Wrapper {
            steps: vec![1, 2, 3].into(),
        };
        let edn = to_string(&w).unwrap();
        assert_eq!(edn, "{:steps (1 2 3)}");
        let parsed: Wrapper = from_str(&edn).unwrap();
        assert_eq!(parsed, w);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_list_nested_in_vector() {
        let outer = vec![
            EdnList::<i64>::from(vec![1, 2]),
            EdnList::<i64>::from(vec![3, 4]),
        ];
        let edn = to_string(&outer).unwrap();
        assert_eq!(edn, "[(1 2) (3 4)]");
        let parsed: Vec<EdnList<i64>> = from_str(&edn).unwrap();
        assert_eq!(parsed, outer);
    }
}
