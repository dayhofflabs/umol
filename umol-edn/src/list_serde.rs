//! Serde Serialize/Deserialize impls for [`EdnList`].
//!
//! When serialized via the EDN serializers, `EdnList<T>` emits `(...)` syntax
//! via the `LIST_TOKEN` newtype-struct trick. Non-EDN serializers see a
//! transparent sequence (the same as a plain `Vec<T>`). On the deserialization
//! side, only `Edn::List` is accepted; `Edn::Vector` and `Edn::Set` are
//! rejected.

use std::fmt;
use std::marker::PhantomData;

use serde::de::{Deserialize, DeserializeSeed, Deserializer, SeqAccess, Visitor};
use serde::ser::{Serialize, SerializeSeq, Serializer};

use crate::list_owned::EdnList;
use crate::serde_tokens::LIST_TOKEN;

impl<T: Serialize> Serialize for EdnList<T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_newtype_struct(LIST_TOKEN, &ListPayload(&self.0))
    }
}

struct ListPayload<'a, T>(&'a Vec<T>);

impl<T: Serialize> Serialize for ListPayload<'_, T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;
        for v in self.0 {
            seq.serialize_element(v)?;
        }
        seq.end()
    }
}

struct ListVisitor<T> {
    _marker: PhantomData<T>,
}

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

impl<'de, T: Deserialize<'de>> DeserializeSeed<'de> for PhantomDataSeed<T> {
    type Value = T;

    fn deserialize<D: Deserializer<'de>>(self, deserializer: D) -> Result<T, D::Error> {
        T::deserialize(deserializer)
    }
}

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
    use super::*;
    use crate::de::from_str;
    use crate::ser::to_string;

    #[test]
    fn test_edn_list_serialize_edn() {
        let list: EdnList<i64> = vec![1, 2, 3].into();
        let edn = to_string(&list).unwrap();
        assert_eq!(edn, "(1 2 3)");
    }

    #[test]
    fn test_edn_list_serialize_empty() {
        let list: EdnList<i64> = EdnList::new();
        let edn = to_string(&list).unwrap();
        assert_eq!(edn, "()");
    }

    #[test]
    fn test_edn_list_serialize_json_falls_back_to_seq() {
        let list: EdnList<i64> = vec![1, 2, 3].into();
        let json = serde_json::to_string(&list).unwrap();
        assert_eq!(json, "[1,2,3]");
    }

    #[test]
    fn test_edn_list_deserialize_from_edn_list() {
        let list: EdnList<i64> = from_str("(1 2 3)").unwrap();
        assert_eq!(list.0, vec![1, 2, 3]);
    }

    #[test]
    fn test_edn_list_deserialize_from_edn_vector_rejected() {
        let result: Result<EdnList<i64>, _> = from_str("[1 2 3]");
        assert!(result.is_err(), "expected error, got {:?}", result);
    }

    #[test]
    fn test_edn_list_roundtrip() {
        let list: EdnList<i64> = vec![1, 2, 3].into();
        let edn = to_string(&list).unwrap();
        let parsed: EdnList<i64> = from_str(&edn).unwrap();
        assert_eq!(parsed, list);
    }

    #[test]
    fn test_edn_list_in_struct() {
        #[derive(serde::Serialize, serde::Deserialize, Debug, PartialEq)]
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
