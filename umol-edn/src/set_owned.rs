//! Owned EDN set wrapper for use in struct fields.
//!
//! `EdnHashSet<T>` is a thin wrapper around `std::collections::HashSet<T>`
//! that opts into EDN set syntax (`#{...}`) when serialized via the EDN
//! serde layer, and accepts only `Edn::Set` on the deserialization side.
//! Non-EDN serializers see it as an ordinary sequence (the same as any
//! other `HashSet<T>` would be).

use std::collections::HashSet;
use std::hash::Hash;
use std::ops::{Deref, DerefMut};

use crate::edn::Edn;
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_string;

    #[test]
    fn test_edn_hash_set_from_edn_set() {
        let edn = read_string("#{1 2 3}").unwrap();
        let set: EdnHashSet<i64> = EdnHashSet::from_edn(&edn).unwrap();
        assert_eq!(set.len(), 3);
        assert!(set.contains(&1));
        assert!(set.contains(&2));
        assert!(set.contains(&3));
    }

    #[test]
    fn test_edn_hash_set_from_edn_vector_rejected() {
        let edn = read_string("[1 2 3]").unwrap();
        let result: Result<EdnHashSet<i64>, _> = EdnHashSet::from_edn(&edn);
        assert!(result.is_err());
    }

    #[test]
    fn test_edn_hash_set_to_edn() {
        let set: EdnHashSet<i64> = [1i64, 2, 3].into_iter().collect();
        let edn = set.to_edn();
        if let Edn::Set(s) = &edn {
            assert_eq!(s.len(), 3);
        } else {
            panic!("expected Edn::Set, got {:?}", edn);
        }
    }
}
