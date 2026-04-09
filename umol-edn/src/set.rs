//! Owned EDN set wrapper for use in struct fields.
//!
//! `EdnHashSet<T>` is a thin wrapper around `std::collections::HashSet<T>`
//! that opts into EDN set syntax (`#{...}`) when serialized via the EDN
//! serde layer, and accepts only `Edn::Set` on the deserialization side.
//! Non-EDN serializers see it as an ordinary sequence.
//!
//! Under the `serde` feature, `Serialize` emits `#{...}` via the `SET_TOKEN`
//! newtype-struct trick when routed through `EdnSerializer`; other serializers
//! see a transparent sequence. `Deserialize` via the EDN tree path accepts
//! only `Edn::Set`; `Edn::Vector` and `Edn::List` are rejected.

use std::collections::HashSet;
use std::hash::Hash;
use std::ops::{Deref, DerefMut};

#[cfg(feature = "serde")]
use std::fmt;
#[cfg(feature = "serde")]
use std::marker::PhantomData;

#[cfg(feature = "serde")]
use serde::de::{Deserialize, DeserializeSeed, Deserializer, SeqAccess, Visitor};
#[cfg(feature = "serde")]
use serde::ser::{Serialize, SerializeSeq, Serializer};

use crate::edn::Edn;
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};
#[cfg(feature = "serde")]
use crate::serde_tokens::SET_TOKEN;

/// An owned EDN set. Always serializes as `#{...}` in EDN.
#[derive(Clone, Debug, Default)]
pub struct EdnHashSet<T>(pub HashSet<T>);

impl<T: Eq + Hash> PartialEq for EdnHashSet<T> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<T: Eq + Hash> Eq for EdnHashSet<T> {}

impl<T> EdnHashSet<T> {
    pub fn new() -> Self {
        Self(HashSet::new())
    }

    pub fn into_inner(self) -> HashSet<T> {
        self.0
    }
}

impl<T> Deref for EdnHashSet<T> {
    type Target = HashSet<T>;
    fn deref(&self) -> &HashSet<T> {
        &self.0
    }
}

impl<T> DerefMut for EdnHashSet<T> {
    fn deref_mut(&mut self) -> &mut HashSet<T> {
        &mut self.0
    }
}

impl<T> From<HashSet<T>> for EdnHashSet<T> {
    fn from(set: HashSet<T>) -> Self {
        Self(set)
    }
}

impl<T: Eq + Hash> FromIterator<T> for EdnHashSet<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<'de, T> FromEdn<'de> for EdnHashSet<T>
where
    T: FromEdn<'de> + Eq + Hash,
{
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Set(s) => {
                let mut out = HashSet::with_capacity(s.len());
                for v in s.iter() {
                    out.insert(T::from_edn(v)?);
                }
                Ok(EdnHashSet(out))
            }
            other => Err(EdnError::TypeMismatch {
                expected: "set",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<T> ToEdn for EdnHashSet<T>
where
    T: ToEdn,
{
    fn to_edn(&self) -> Edn<'_> {
        let mut set = crate::collections::EdnSet::new();
        for v in &self.0 {
            set.insert(v.to_edn().into_owned());
        }
        Edn::Set(set)
    }
}

#[cfg(feature = "serde")]
impl<T: Serialize> Serialize for EdnHashSet<T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(SET_TOKEN, &SetPayload(&self.0))
    }
}

#[cfg(feature = "serde")]
struct SetPayload<'a, T>(&'a HashSet<T>);

#[cfg(feature = "serde")]
impl<T: Serialize> Serialize for SetPayload<'_, T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;
        for v in self.0 {
            seq.serialize_element(v)?;
        }
        seq.end()
    }
}

#[cfg(feature = "serde")]
struct SetVisitor<T> {
    _marker: PhantomData<T>,
}

#[cfg(feature = "serde")]
impl<'de, T> Visitor<'de> for SetVisitor<T>
where
    T: Deserialize<'de> + Eq + Hash,
{
    type Value = EdnHashSet<T>;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("an EDN set")
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<EdnHashSet<T>, A::Error> {
        let mut out = HashSet::with_capacity(seq.size_hint().unwrap_or(0));
        while let Some(v) = seq.next_element_seed(PhantomDataSeed::<T>::default())? {
            out.insert(v);
        }
        Ok(EdnHashSet(out))
    }

    fn visit_newtype_struct<D: Deserializer<'de>>(
        self,
        deserializer: D,
    ) -> Result<EdnHashSet<T>, D::Error> {
        HashSet::<T>::deserialize(deserializer).map(EdnHashSet)
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
impl<'de, T> Deserialize<'de> for EdnHashSet<T>
where
    T: Deserialize<'de> + Eq + Hash,
{
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_newtype_struct(
            SET_TOKEN,
            SetVisitor {
                _marker: PhantomData,
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_string;
    #[cfg(feature = "serde")]
    use crate::{from_str, to_string};

    #[test]
    fn test_edn_hash_set_from_edn() {
        let edn = read_string("#{1 2 3}").unwrap();
        let set: EdnHashSet<i64> = EdnHashSet::from_edn(&edn).unwrap();
        assert_eq!(set.len(), 3);
        assert!(set.contains(&1));
        assert!(set.contains(&2));
        assert!(set.contains(&3));
    }

    #[test]
    fn test_edn_hash_set_from_edn_error() {
        let edn = read_string("[1 2 3]").unwrap();
        let result: Result<EdnHashSet<i64>, _> = EdnHashSet::from_edn(&edn);
        assert!(result.is_err());
    }

    #[test]
    fn test_edn_hash_set_to_edn() {
        let set: EdnHashSet<i64> = [1i64, 2, 3].into_iter().collect();
        let edn = set.to_edn();
        let Edn::Set(s) = &edn else {
            panic!("expected Edn::Set, got {:?}", edn);
        };
        assert_eq!(s.len(), 3);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_hash_set_serialize() {
        let set: EdnHashSet<i64> = [1i64, 2, 3].into_iter().collect();
        let edn = to_string(&set).unwrap();
        assert!(edn.starts_with("#{"));
        assert!(edn.ends_with('}'));
        for elem in ["1", "2", "3"] {
            assert!(edn.contains(elem), "missing {elem} in {edn}");
        }
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_hash_set_serialize_empty() {
        let set: EdnHashSet<i64> = EdnHashSet::new();
        assert_eq!(to_string(&set).unwrap(), "#{}");
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_hash_set_serialize_json() {
        let set: EdnHashSet<i64> = [1i64, 2, 3].into_iter().collect();
        let json = serde_json::to_string(&set).unwrap();
        assert!(json.starts_with('['));
        assert!(json.ends_with(']'));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_hash_set_deserialize() {
        let set: EdnHashSet<i64> = from_str("#{1 2 3}").unwrap();
        assert_eq!(set.len(), 3);
        assert!(set.contains(&1));
        assert!(set.contains(&2));
        assert!(set.contains(&3));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_hash_set_deserialize_error() {
        let result: Result<EdnHashSet<i64>, _> = from_str("[1 2 3]");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_hash_set_roundtrip() {
        let set: EdnHashSet<i64> = [1i64, 2, 3].into_iter().collect();
        let edn = to_string(&set).unwrap();
        let parsed: EdnHashSet<i64> = from_str(&edn).unwrap();
        assert_eq!(parsed, set);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_hash_set_serialize_in_struct() {
        #[derive(serde::Serialize, serde::Deserialize, Debug, PartialEq)]
        struct Wrapper {
            tags: EdnHashSet<String>,
        }
        let w = Wrapper {
            tags: ["a".to_string(), "b".to_string()].into_iter().collect(),
        };
        let edn = to_string(&w).unwrap();
        let parsed: Wrapper = from_str(&edn).unwrap();
        assert_eq!(parsed, w);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_edn_hash_set_nested_in_vector() {
        let outer = vec![
            EdnHashSet::<i64>::from_iter([1i64, 2]),
            EdnHashSet::<i64>::from_iter([3i64, 4]),
        ];
        let edn = to_string(&outer).unwrap();
        assert!(edn.starts_with("[#{"));
        let parsed: Vec<EdnHashSet<i64>> = from_str(&edn).unwrap();
        assert_eq!(parsed.len(), 2);
    }
}
