//! Serde Serialize/Deserialize impls for [`EdnHashSet`].
//!
//! When serialized via the EDN serializers, `EdnHashSet<T>` emits `#{...}`
//! syntax via the `SET_TOKEN` newtype-struct trick. Non-EDN serializers
//! see a transparent sequence (the same as a plain `HashSet<T>`). On the
//! deserialization side, only `Edn::Set` is accepted; `Edn::Vector` and
//! `Edn::List` are rejected.

use std::collections::HashSet;
use std::fmt;
use std::hash::Hash;
use std::marker::PhantomData;

use serde::de::{Deserialize, DeserializeSeed, Deserializer, SeqAccess, Visitor};
use serde::ser::{Serialize, SerializeSeq, Serializer};

use crate::serde_tokens::SET_TOKEN;
use crate::set_owned::EdnHashSet;

impl<T: Serialize> Serialize for EdnHashSet<T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(SET_TOKEN, &SetPayload(&self.0))
    }
}

/// Payload wrapper that serializes the inner `HashSet<T>` as a sequence.
/// Used so the EDN serializer's `SET_TOKEN` arm receives a value that calls
/// `serialize_seq` (rather than `serialize_newtype_struct` again).
struct SetPayload<'a, T>(&'a HashSet<T>);

impl<T: Serialize> Serialize for SetPayload<'_, T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;
        for v in self.0 {
            seq.serialize_element(v)?;
        }
        seq.end()
    }
}

struct SetVisitor<T> {
    _marker: PhantomData<T>,
}

impl<'de, T> Visitor<'de> for SetVisitor<T>
where
    T: serde::Deserialize<'de> + Eq + Hash,
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
        // Non-EDN deserializers fall through here: read the inner HashSet.
        HashSet::<T>::deserialize(deserializer).map(EdnHashSet)
    }
}

struct PhantomDataSeed<T> {
    _marker: PhantomData<T>,
}

impl<T> Default for PhantomDataSeed<T> {
    fn default() -> Self {
        Self {
            _marker: PhantomData,
        }
    }
}

impl<'de, T: serde::Deserialize<'de>> DeserializeSeed<'de> for PhantomDataSeed<T> {
    type Value = T;

    fn deserialize<D: Deserializer<'de>>(self, deserializer: D) -> Result<T, D::Error> {
        T::deserialize(deserializer)
    }
}

impl<'de, T> serde::Deserialize<'de> for EdnHashSet<T>
where
    T: serde::Deserialize<'de> + Eq + Hash,
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
    use crate::de::from_str;
    use crate::ser::to_string;

    #[test]
    fn test_edn_hash_set_serialize_edn() {
        let set: EdnHashSet<i64> = [1i64, 2, 3].into_iter().collect();
        let edn = to_string(&set).unwrap();
        // Iteration order is unspecified; check structurally.
        assert!(edn.starts_with("#{"));
        assert!(edn.ends_with('}'));
        for elem in ["1", "2", "3"] {
            assert!(edn.contains(elem), "missing {elem} in {edn}");
        }
    }

    #[test]
    fn test_edn_hash_set_serialize_empty() {
        let set: EdnHashSet<i64> = EdnHashSet::new();
        let edn = to_string(&set).unwrap();
        assert_eq!(edn, "#{}");
    }

    #[test]
    fn test_edn_hash_set_serialize_json_falls_back_to_seq() {
        let set: EdnHashSet<i64> = [1i64, 2, 3].into_iter().collect();
        let json = serde_json::to_string(&set).unwrap();
        // JSON uses [...] for the inner HashSet; order is unspecified.
        assert!(json.starts_with('['));
        assert!(json.ends_with(']'));
    }

    #[test]
    fn test_edn_hash_set_deserialize_from_edn_set() {
        let set: EdnHashSet<i64> = from_str("#{1 2 3}").unwrap();
        assert_eq!(set.len(), 3);
        assert!(set.contains(&1));
        assert!(set.contains(&2));
        assert!(set.contains(&3));
    }

    #[test]
    fn test_edn_hash_set_deserialize_from_edn_vector_rejected() {
        let result: Result<EdnHashSet<i64>, _> = from_str("[1 2 3]");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[test]
    fn test_edn_hash_set_roundtrip() {
        let set: EdnHashSet<i64> = [1i64, 2, 3].into_iter().collect();
        let edn = to_string(&set).unwrap();
        let parsed: EdnHashSet<i64> = from_str(&edn).unwrap();
        assert_eq!(parsed, set);
    }

    #[test]
    fn test_edn_hash_set_in_struct() {
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
