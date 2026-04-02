//! EDN collection newtypes: EdnMap, EdnSet, EdnSeq.

use std::cmp::Ordering;
use std::collections::hash_map::{IntoIter as HashMapIntoIter, Iter as HashMapIter};
use std::collections::hash_set::{IntoIter as HashSetIntoIter, Iter as HashSetIter};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::ops::Deref;
use std::slice::Iter as SliceIter;
use std::vec::IntoIter as VecIntoIter;

use rustc_hash::{FxHashMap, FxHashSet};

use crate::edn::Edn;

// ---------------------------------------------------------------------------
// EdnMap
// ---------------------------------------------------------------------------

/// An unordered map of EDN values.
#[derive(Clone, Debug)]
pub struct EdnMap<'a>(FxHashMap<Edn<'a>, Edn<'a>>);

impl<'a> EdnMap<'a> {
    pub fn new() -> Self {
        Self(FxHashMap::default())
    }

    pub fn with_capacity(cap: usize) -> Self {
        Self(FxHashMap::with_capacity_and_hasher(cap, Default::default()))
    }

    pub fn insert(&mut self, key: Edn<'a>, value: Edn<'a>) -> Option<Edn<'a>> {
        self.0.insert(key, value)
    }

    pub fn remove(&mut self, key: &Edn<'a>) -> Option<Edn<'a>> {
        self.0.remove(key)
    }

    pub fn get(&self, key: &Edn<'a>) -> Option<&Edn<'a>> {
        self.0.get(key)
    }

    pub fn contains_key(&self, key: &Edn<'a>) -> bool {
        self.0.contains_key(key)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&Edn<'a>, &Edn<'a>)> {
        self.0.iter()
    }

    pub fn keys(&self) -> impl Iterator<Item = &Edn<'a>> {
        self.0.keys()
    }

    pub fn values(&self) -> impl Iterator<Item = &Edn<'a>> {
        self.0.values()
    }

    pub fn into_owned(self) -> EdnMap<'static> {
        EdnMap(
            self.0
                .into_iter()
                .map(|(k, v)| (k.into_owned(), v.into_owned()))
                .collect(),
        )
    }
}

impl Default for EdnMap<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> IntoIterator for EdnMap<'a> {
    type Item = (Edn<'a>, Edn<'a>);
    type IntoIter = HashMapIntoIter<Edn<'a>, Edn<'a>>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a, 'b> IntoIterator for &'b EdnMap<'a> {
    type Item = (&'b Edn<'a>, &'b Edn<'a>);
    type IntoIter = HashMapIter<'b, Edn<'a>, Edn<'a>>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> FromIterator<(Edn<'a>, Edn<'a>)> for EdnMap<'a> {
    fn from_iter<T: IntoIterator<Item = (Edn<'a>, Edn<'a>)>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl PartialEq for EdnMap<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for EdnMap<'_> {}

impl PartialOrd for EdnMap<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EdnMap<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        let mut a: Vec<_> = self.0.iter().collect();
        let mut b: Vec<_> = other.0.iter().collect();
        a.sort_by(|x, y| x.0.cmp(y.0).then_with(|| x.1.cmp(y.1)));
        b.sort_by(|x, y| x.0.cmp(y.0).then_with(|| x.1.cmp(y.1)));
        a.cmp(&b)
    }
}

impl Hash for EdnMap<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_usize(self.0.len());
        // Order-independent hash: XOR individual pair hashes.
        let mut combined = 0u64;
        for (k, v) in &self.0 {
            let mut pair_hasher = DefaultHasher::new();
            k.hash(&mut pair_hasher);
            v.hash(&mut pair_hasher);
            combined ^= pair_hasher.finish();
        }
        combined.hash(state);
    }
}

// ---------------------------------------------------------------------------
// EdnSet
// ---------------------------------------------------------------------------

/// An unordered set of EDN values.
#[derive(Clone, Debug)]
pub struct EdnSet<'a>(FxHashSet<Edn<'a>>);

impl<'a> EdnSet<'a> {
    pub fn new() -> Self {
        Self(FxHashSet::default())
    }

    pub fn insert(&mut self, value: Edn<'a>) -> bool {
        self.0.insert(value)
    }

    pub fn contains(&self, value: &Edn<'a>) -> bool {
        self.0.contains(value)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Edn<'a>> {
        self.0.iter()
    }

    pub fn into_owned(self) -> EdnSet<'static> {
        EdnSet(self.0.into_iter().map(Edn::into_owned).collect())
    }
}

impl Default for EdnSet<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> IntoIterator for EdnSet<'a> {
    type Item = Edn<'a>;
    type IntoIter = HashSetIntoIter<Edn<'a>>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a, 'b> IntoIterator for &'b EdnSet<'a> {
    type Item = &'b Edn<'a>;
    type IntoIter = HashSetIter<'b, Edn<'a>>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> FromIterator<Edn<'a>> for EdnSet<'a> {
    fn from_iter<T: IntoIterator<Item = Edn<'a>>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl PartialEq for EdnSet<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for EdnSet<'_> {}

impl PartialOrd for EdnSet<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EdnSet<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        let mut a: Vec<_> = self.0.iter().collect();
        let mut b: Vec<_> = other.0.iter().collect();
        a.sort();
        b.sort();
        a.cmp(&b)
    }
}

impl Hash for EdnSet<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_usize(self.0.len());
        let mut combined = 0u64;
        for v in &self.0 {
            let mut item_hasher = DefaultHasher::new();
            v.hash(&mut item_hasher);
            combined ^= item_hasher.finish();
        }
        combined.hash(state);
    }
}

// ---------------------------------------------------------------------------
// EdnSeq
// ---------------------------------------------------------------------------

/// An ordered sequence of EDN values (backing both lists and vectors).
#[derive(Clone, Debug)]
pub struct EdnSeq<'a>(Vec<Edn<'a>>);

pub(crate) type EdnSeqIntoIter<'a> = VecIntoIter<Edn<'a>>;

impl<'a> EdnSeq<'a> {
    pub fn new() -> Self {
        Self(Default::default())
    }

    pub fn push(&mut self, value: Edn<'a>) {
        self.0.push(value);
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn into_owned(self) -> EdnSeq<'static> {
        EdnSeq(self.0.into_iter().map(Edn::into_owned).collect())
    }
}

impl Default for EdnSeq<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> Deref for EdnSeq<'a> {
    type Target = [Edn<'a>];
    fn deref(&self) -> &[Edn<'a>] {
        &self.0
    }
}

impl<'a> From<Vec<Edn<'a>>> for EdnSeq<'a> {
    fn from(v: Vec<Edn<'a>>) -> Self {
        Self(v.into_iter().collect())
    }
}

impl<'a> FromIterator<Edn<'a>> for EdnSeq<'a> {
    fn from_iter<T: IntoIterator<Item = Edn<'a>>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<'a> IntoIterator for EdnSeq<'a> {
    type Item = Edn<'a>;
    type IntoIter = EdnSeqIntoIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a, 'b> IntoIterator for &'b EdnSeq<'a> {
    type Item = &'b Edn<'a>;
    type IntoIter = SliceIter<'b, Edn<'a>>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl PartialEq for EdnSeq<'_> {
    fn eq(&self, other: &Self) -> bool {
        **self == **other
    }
}

impl Eq for EdnSeq<'_> {}

impl PartialOrd for EdnSeq<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EdnSeq<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.iter().cmp(other.iter())
    }
}

impl Hash for EdnSeq<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.len().hash(state);
        for item in &self.0 {
            item.hash(state);
        }
    }
}
