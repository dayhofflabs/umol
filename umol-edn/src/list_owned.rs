//! Owned EDN list wrapper for use in struct fields.
//!
//! `EdnList<T>` is a thin wrapper around `Vec<T>` that opts into EDN list
//! syntax (`(...)`) rather than vector syntax (`[...]`) when serialized via
//! the EDN serde layer, and accepts only `Edn::List` on the deserialization
//! side. Non-EDN serializers see it as an ordinary sequence (the same as
//! any other `Vec<T>`).

use std::ops::{Deref, DerefMut};

use crate::collections::EdnSeq;
use crate::edn::Edn;
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};

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
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::List(v) => {
                let mut out = Vec::with_capacity(v.len());
                for item in v.iter() {
                    out.push(T::from_edn(item)?);
                }
                Ok(EdnList(out))
            }
            other => Err(EdnError::TypeMismatch {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_string;

    #[test]
    fn test_edn_list_from_edn_list() {
        let edn = read_string("(1 2 3)").unwrap();
        let list: EdnList<i64> = EdnList::from_edn(&edn).unwrap();
        assert_eq!(list.0, vec![1, 2, 3]);
    }

    #[test]
    fn test_edn_list_from_edn_vector_rejected() {
        let edn = read_string("[1 2 3]").unwrap();
        let result: Result<EdnList<i64>, _> = EdnList::from_edn(&edn);
        assert!(result.is_err());
    }

    #[test]
    fn test_edn_list_to_edn() {
        let list: EdnList<i64> = vec![1i64, 2, 3].into();
        let edn = list.to_edn();
        if let Edn::List(v) = &edn {
            assert_eq!(v.len(), 3);
        } else {
            panic!("expected Edn::List, got {:?}", edn);
        }
    }
}
