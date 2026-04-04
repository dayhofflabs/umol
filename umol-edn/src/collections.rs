//! EDN collection newtypes: EdnMap, EdnSet, EdnSeq.

use std::cmp::Ordering;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::ops::Deref;
use std::slice::Iter as SliceIter;
use std::vec::IntoIter as VecIntoIter;

use hashbrown::hash_map::{IntoIter as HashMapIntoIter, Iter as HashMapIter};
use hashbrown::hash_set::{IntoIter as HashSetIntoIter, Iter as HashSetIter};
use hashbrown::{HashMap, HashSet};

use crate::edn::Edn;

// ---------------------------------------------------------------------------
// EdnKeyRef — borrowed key view for cross-lifetime lookups
// ---------------------------------------------------------------------------

/// Borrowed key representation for looking up values in [`EdnMap`] and [`EdnSet`]
/// without constructing an owned [`Edn`] or matching its lifetime parameter.
pub enum EdnKeyRef<'k> {
    Nil,
    Bool(bool),
    Int(i64),
    Float(f64),
    Char(char),
    Str(&'k str),
    Keyword(&'k str),
    Symbol(&'k str),
    List(&'k [Edn<'k>]),
    Vector(&'k [Edn<'k>]),
    Map(&'k EdnMap<'k>),
    Set(&'k EdnSet<'k>),
    Tagged(&'k str, &'k Edn<'k>),
    #[cfg(feature = "bignum")]
    BigInt(&'k num_bigint::BigInt),
    #[cfg(feature = "bignum")]
    BigDecimal(&'k crate::edn::EdnBigDecimal),
}

impl<'k> EdnKeyRef<'k> {
    pub fn keyword(s: &'k str) -> Self {
        Self::Keyword(s)
    }

    pub fn symbol(s: &'k str) -> Self {
        Self::Symbol(s)
    }

    pub fn str_(s: &'k str) -> Self {
        Self::Str(s)
    }
}

impl<'k> From<&'k Edn<'k>> for EdnKeyRef<'k> {
    fn from(edn: &'k Edn<'k>) -> Self {
        match edn {
            Edn::Nil => Self::Nil,
            Edn::Bool(b) => Self::Bool(*b),
            Edn::Int(n) => Self::Int(*n),
            Edn::Float(f) => Self::Float(*f),
            Edn::Char(c) => Self::Char(*c),
            Edn::Str(s) => Self::Str(s),
            Edn::Keyword(k) => Self::Keyword(k.as_str()),
            Edn::Symbol(s) => Self::Symbol(s.as_str()),
            Edn::List(v) => Self::List(v),
            Edn::Vector(v) => Self::Vector(v),
            Edn::Map(m) => Self::Map(m),
            Edn::Set(s) => Self::Set(s),
            Edn::Tagged(tag, inner) => Self::Tagged(tag, inner),
            #[cfg(feature = "bignum")]
            Edn::BigInt(n) => Self::BigInt(n),
            #[cfg(feature = "bignum")]
            Edn::BigDecimal(d) => Self::BigDecimal(d),
        }
    }
}

/// Discriminant byte, must match `variant_ord` in `edn.rs` exactly.
fn key_ref_variant_ord(k: &EdnKeyRef<'_>) -> u8 {
    match k {
        EdnKeyRef::Nil => 0,
        EdnKeyRef::Bool(_) => 1,
        EdnKeyRef::Int(_) => 2,
        #[cfg(feature = "bignum")]
        EdnKeyRef::BigInt(_) => 3,
        EdnKeyRef::Float(_) => 4,
        #[cfg(feature = "bignum")]
        EdnKeyRef::BigDecimal(_) => 5,
        EdnKeyRef::Char(_) => 6,
        EdnKeyRef::Str(_) => 7,
        EdnKeyRef::Keyword(_) => 8,
        EdnKeyRef::Symbol(_) => 9,
        EdnKeyRef::List(_) | EdnKeyRef::Vector(_) => 10,
        EdnKeyRef::Map(_) => 11,
        EdnKeyRef::Set(_) => 12,
        EdnKeyRef::Tagged(_, _) => 13,
    }
}

impl Hash for EdnKeyRef<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        key_ref_variant_ord(self).hash(state);
        match self {
            Self::Nil => {}
            Self::Bool(b) => b.hash(state),
            Self::Int(n) => n.hash(state),
            #[cfg(feature = "bignum")]
            Self::BigInt(n) => n.hash(state),
            Self::Float(f) => f.to_bits().hash(state),
            #[cfg(feature = "bignum")]
            Self::BigDecimal(d) => d.hash(state),
            Self::Char(c) => c.hash(state),
            Self::Str(s) => s.hash(state),
            Self::Keyword(s) => s.hash(state),
            Self::Symbol(s) => s.hash(state),
            Self::List(v) | Self::Vector(v) => {
                v.len().hash(state);
                for item in *v {
                    item.hash(state);
                }
            }
            Self::Map(m) => m.hash(state),
            Self::Set(s) => s.hash(state),
            Self::Tagged(tag, inner) => {
                tag.hash(state);
                (*inner).hash(state);
            }
        }
    }
}

impl hashbrown::Equivalent<Edn<'_>> for EdnKeyRef<'_> {
    fn equivalent(&self, other: &Edn<'_>) -> bool {
        match (self, other) {
            (Self::Nil, Edn::Nil) => true,
            (Self::Bool(a), Edn::Bool(b)) => a == b,
            (Self::Int(a), Edn::Int(b)) => a == b,
            #[cfg(feature = "bignum")]
            (Self::BigInt(a), Edn::BigInt(b)) => *a == b,
            (Self::Float(a), Edn::Float(b)) => a.to_bits() == b.to_bits(),
            #[cfg(feature = "bignum")]
            (Self::BigDecimal(a), Edn::BigDecimal(b)) => *a == b,
            (Self::Char(a), Edn::Char(b)) => a == b,
            (Self::Str(a), Edn::Str(b)) => *a == &**b,
            (Self::Keyword(a), Edn::Keyword(b)) => *a == b.as_str(),
            (Self::Symbol(a), Edn::Symbol(b)) => *a == b.as_str(),
            (Self::List(a) | Self::Vector(a), Edn::List(b) | Edn::Vector(b)) => *a == &**b,
            (Self::Map(a), Edn::Map(b)) => *a == b,
            (Self::Set(a), Edn::Set(b)) => *a == b,
            (Self::Tagged(ta, va), Edn::Tagged(tb, vb)) => *ta == &**tb && *va == &**vb,
            _ => false,
        }
    }
}

// ---------------------------------------------------------------------------
// EdnMap
// ---------------------------------------------------------------------------

/// An unordered map of EDN values.
#[derive(Clone, Debug)]
pub struct EdnMap<'a>(HashMap<Edn<'a>, Edn<'a>>);

impl<'a> EdnMap<'a> {
    pub fn new() -> Self {
        Self(HashMap::new())
    }

    pub fn with_capacity(cap: usize) -> Self {
        Self(HashMap::with_capacity(cap))
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

    pub fn get_ref(&self, key: EdnKeyRef<'_>) -> Option<&Edn<'a>> {
        self.0.get(&key)
    }

    pub fn contains_ref(&self, key: EdnKeyRef<'_>) -> bool {
        self.0.contains_key(&key)
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
        self.0.len() == other.0.len()
            && self.0.iter().all(|(k, v)| other.0.get(k) == Some(v))
    }
}

impl Eq for EdnMap<'_> {}

impl PartialOrd for EdnMap<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EdnMap<'_> {
    /// O(n log n) with allocation — sorts entries on each call. Used for
    /// deterministic comparison of unordered maps, not intended for hot paths.
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
pub struct EdnSet<'a>(HashSet<Edn<'a>>);

impl<'a> EdnSet<'a> {
    pub fn new() -> Self {
        Self(HashSet::new())
    }

    pub fn insert(&mut self, value: Edn<'a>) -> bool {
        self.0.insert(value)
    }

    pub fn contains(&self, value: &Edn<'a>) -> bool {
        self.0.contains(value)
    }

    pub fn contains_ref(&self, value: EdnKeyRef<'_>) -> bool {
        self.0.contains(&value)
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
        self.0.len() == other.0.len()
            && self.0.iter().all(|v| other.0.contains(v))
    }
}

impl Eq for EdnSet<'_> {}

impl PartialOrd for EdnSet<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EdnSet<'_> {
    /// O(n log n) with allocation — sorts elements on each call. Used for
    /// deterministic comparison of unordered sets, not intended for hot paths.
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
        Self(v)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::edn::{Keyword, Symbol};
    use hashbrown::Equivalent;
    use rstest::rstest;
    use std::borrow::Cow;
    use std::hash::{DefaultHasher, Hash, Hasher};

    fn hash_edn(v: &Edn<'_>) -> u64 {
        let mut h = DefaultHasher::new();
        v.hash(&mut h);
        h.finish()
    }

    fn hash_key_ref(v: &EdnKeyRef<'_>) -> u64 {
        let mut h = DefaultHasher::new();
        v.hash(&mut h);
        h.finish()
    }

    // -- Eq/Hash parity: EdnKeyRef must produce identical hashes and match Edn --

    #[rstest]
    #[case::nil(Edn::Nil)]
    #[case::bool_true(Edn::Bool(true))]
    #[case::bool_false(Edn::Bool(false))]
    #[case::int(Edn::Int(77))]
    #[case::int_neg(Edn::Int(-1))]
    #[case::float(Edn::Float(3.14))]
    #[case::float_zero(Edn::Float(0.0))]
    #[case::float_neg_zero(Edn::Float(-0.0))]
    #[case::char(Edn::Char('x'))]
    #[case::str(Edn::Str(Cow::Borrowed("hello")))]
    #[case::keyword(Edn::Keyword(Keyword::new("foo")))]
    #[case::symbol(Edn::Symbol(Symbol::new("bar")))]
    #[case::list(Edn::List(EdnSeq::from(vec![Edn::Int(1), Edn::Int(2)])))]
    #[case::vector(Edn::Vector(EdnSeq::from(vec![Edn::Bool(true)])))]
    #[case::tagged(Edn::Tagged(Cow::Borrowed("my/tag"), Box::new(Edn::Int(5))))]
    fn test_edn_key_ref_hash_parity(#[case] edn: Edn<'_>) {
        let key_ref = EdnKeyRef::from(&edn);
        assert_eq!(hash_edn(&edn), hash_key_ref(&key_ref));
        assert!(key_ref.equivalent(&edn));
    }

    #[test]
    fn test_edn_key_ref_hash_parity_map() {
        let mut m = EdnMap::new();
        m.insert(Edn::keyword("a"), Edn::Int(1));
        let edn = Edn::Map(m);
        let key_ref = EdnKeyRef::from(&edn);
        assert_eq!(hash_edn(&edn), hash_key_ref(&key_ref));
        assert!(key_ref.equivalent(&edn));
    }

    #[test]
    fn test_edn_key_ref_hash_parity_set() {
        let mut s = EdnSet::new();
        s.insert(Edn::Int(1));
        s.insert(Edn::Int(2));
        let edn = Edn::Set(s);
        let key_ref = EdnKeyRef::from(&edn);
        assert_eq!(hash_edn(&edn), hash_key_ref(&key_ref));
        assert!(key_ref.equivalent(&edn));
    }

    // -- Float edge cases --

    #[test]
    fn test_edn_key_ref_nan_payload() {
        let nan1 = Edn::Float(f64::from_bits(0x7FF8_0000_0000_0001));
        let nan2 = Edn::Float(f64::from_bits(0x7FF8_0000_0000_0002));
        let ref1 = EdnKeyRef::from(&nan1);
        let ref2 = EdnKeyRef::from(&nan2);
        assert!(ref1.equivalent(&nan1));
        assert!(!ref1.equivalent(&nan2));
        assert_ne!(hash_key_ref(&ref1), hash_key_ref(&ref2));
    }

    #[test]
    fn test_edn_key_ref_pos_neg_zero() {
        let pos = Edn::Float(0.0);
        let neg = Edn::Float(-0.0);
        let ref_pos = EdnKeyRef::from(&pos);
        let ref_neg = EdnKeyRef::from(&neg);
        assert!(ref_pos.equivalent(&pos));
        assert!(!ref_pos.equivalent(&neg));
        assert_ne!(hash_key_ref(&ref_pos), hash_key_ref(&ref_neg));
    }

    // -- Cross-lifetime lookup --

    #[test]
    fn test_edn_map_get_ref_cross_lifetime() {
        let mut m = EdnMap::new();
        m.insert(Edn::keyword("name"), Edn::Str(Cow::Owned("Alice".into())));
        m.insert(Edn::Int(10), Edn::Bool(true));

        {
            let key = String::from("name");
            assert_eq!(
                m.get_ref(EdnKeyRef::keyword(&key)),
                Some(&Edn::Str(Cow::Owned("Alice".into())))
            );
        }
        assert_eq!(m.get_ref(EdnKeyRef::Int(10)), Some(&Edn::Bool(true)));
        assert!(m.get_ref(EdnKeyRef::keyword("missing")).is_none());
    }

    #[test]
    fn test_edn_set_contains_ref_cross_lifetime() {
        let mut s = EdnSet::new();
        s.insert(Edn::keyword("x"));
        s.insert(Edn::Int(7));

        {
            let key = String::from("x");
            assert!(s.contains_ref(EdnKeyRef::keyword(&key)));
        }
        assert!(s.contains_ref(EdnKeyRef::Int(7)));
        assert!(!s.contains_ref(EdnKeyRef::keyword("y")));
    }

    // -- Map/set symmetry: get/contains and get_ref/contains_ref agree --

    #[rstest]
    #[case::nil(Edn::Nil)]
    #[case::bool(Edn::Bool(false))]
    #[case::int(Edn::Int(99))]
    #[case::float(Edn::Float(2.718))]
    #[case::char(Edn::Char('z'))]
    #[case::str(Edn::Str(Cow::Borrowed("test")))]
    #[case::keyword(Edn::Keyword(Keyword::new("k")))]
    #[case::symbol(Edn::Symbol(Symbol::new("s")))]
    fn test_edn_map_get_ref_agrees_with_get(#[case] key: Edn<'static>) {
        let mut m = EdnMap::new();
        m.insert(key.clone(), Edn::Int(1));

        let key_ref = EdnKeyRef::from(&key);
        assert_eq!(m.get(&key), m.get_ref(key_ref));
    }

    #[rstest]
    #[case::nil(Edn::Nil)]
    #[case::bool(Edn::Bool(true))]
    #[case::int(Edn::Int(0))]
    #[case::float(Edn::Float(1.0))]
    #[case::keyword(Edn::Keyword(Keyword::new("k")))]
    fn test_edn_set_contains_ref_agrees_with_contains(#[case] val: Edn<'static>) {
        let mut s = EdnSet::new();
        s.insert(val.clone());

        let key_ref = EdnKeyRef::from(&val);
        assert_eq!(s.contains(&val), s.contains_ref(key_ref));
    }

    // -- Convenience constructors --

    #[test]
    fn test_edn_key_ref_constructors() {
        let m = {
            let mut m = EdnMap::new();
            m.insert(Edn::keyword("a"), Edn::Int(1));
            m.insert(Edn::symbol("b"), Edn::Int(2));
            m.insert(Edn::string("c"), Edn::Int(3));
            m
        };

        assert_eq!(m.get_ref(EdnKeyRef::keyword("a")), Some(&Edn::Int(1)));
        assert_eq!(m.get_ref(EdnKeyRef::symbol("b")), Some(&Edn::Int(2)));
        assert_eq!(m.get_ref(EdnKeyRef::str_("c")), Some(&Edn::Int(3)));
    }

    // -- Non-equivalent variants return false --

    #[test]
    fn test_edn_key_ref_cross_variant_not_equivalent() {
        let int = Edn::Int(1);
        let float = Edn::Float(1.0);
        let ref_int = EdnKeyRef::from(&int);
        assert!(!ref_int.equivalent(&float));
    }

    // -- EdnMap: remove, contains_ref, keys, values, into_owned, trait impls --

    #[test]
    fn test_edn_map_remove() {
        let mut m = EdnMap::new();
        m.insert(Edn::keyword("a"), Edn::Int(1));
        assert_eq!(m.remove(&Edn::keyword("a")), Some(Edn::Int(1)));
        assert_eq!(m.remove(&Edn::keyword("a")), None);
    }

    #[test]
    fn test_edn_map_contains_ref() {
        let mut m = EdnMap::new();
        m.insert(Edn::keyword("x"), Edn::Int(1));
        assert!(m.contains_ref(EdnKeyRef::keyword("x")));
        assert!(!m.contains_ref(EdnKeyRef::keyword("y")));
    }

    #[test]
    fn test_edn_map_keys_values() {
        let mut m = EdnMap::new();
        m.insert(Edn::Int(1), Edn::Bool(true));
        assert_eq!(m.keys().count(), 1);
        assert_eq!(m.values().count(), 1);
    }

    #[test]
    fn test_edn_map_into_owned() {
        let mut m = EdnMap::new();
        m.insert(Edn::keyword("k"), Edn::Str(Cow::Borrowed("v")));
        let owned: EdnMap<'static> = m.into_owned();
        assert_eq!(owned.len(), 1);
    }

    #[test]
    fn test_edn_map_default() {
        let m: EdnMap<'_> = EdnMap::default();
        assert!(m.is_empty());
    }

    #[test]
    fn test_edn_map_ref_into_iter() {
        let mut m = EdnMap::new();
        m.insert(Edn::Int(1), Edn::Int(2));
        let count = (&m).into_iter().count();
        assert_eq!(count, 1);
    }

    #[test]
    fn test_edn_map_from_iter() {
        let m: EdnMap<'_> = vec![
            (Edn::Int(1), Edn::Bool(true)),
            (Edn::Int(2), Edn::Bool(false)),
        ]
        .into_iter()
        .collect();
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn test_edn_map_ord() {
        let mut m1 = EdnMap::new();
        m1.insert(Edn::Int(1), Edn::Int(10));
        let mut m2 = EdnMap::new();
        m2.insert(Edn::Int(2), Edn::Int(20));
        assert!(m1 < m2);
        assert_eq!(m1.partial_cmp(&m2), Some(std::cmp::Ordering::Less));
    }

    // -- EdnSet: is_empty, into_owned, Default, trait impls, Ord --

    #[test]
    fn test_edn_set_is_empty() {
        let s = EdnSet::new();
        assert!(s.is_empty());
        let mut s2 = EdnSet::new();
        s2.insert(Edn::Int(1));
        assert!(!s2.is_empty());
    }

    #[test]
    fn test_edn_set_default() {
        let s: EdnSet<'_> = EdnSet::default();
        assert!(s.is_empty());
    }

    #[test]
    fn test_edn_set_ref_into_iter() {
        let mut s = EdnSet::new();
        s.insert(Edn::Int(1));
        s.insert(Edn::Int(2));
        let count = (&s).into_iter().count();
        assert_eq!(count, 2);
    }

    #[test]
    fn test_edn_set_from_iter() {
        let s: EdnSet<'_> = vec![Edn::Int(1), Edn::Int(2)].into_iter().collect();
        assert_eq!(s.len(), 2);
    }

    #[test]
    fn test_edn_set_ord() {
        let mut s1 = EdnSet::new();
        s1.insert(Edn::Int(1));
        let mut s2 = EdnSet::new();
        s2.insert(Edn::Int(2));
        assert!(s1 < s2);
        assert_eq!(s1.partial_cmp(&s2), Some(std::cmp::Ordering::Less));
    }

    // -- EdnSeq: len, is_empty, Default, FromIterator, ref IntoIterator, PartialOrd --

    #[test]
    fn test_edn_seq_len_is_empty() {
        let s = EdnSeq::new();
        assert_eq!(s.len(), 0);
        assert!(s.is_empty());
        let mut s2 = EdnSeq::new();
        s2.push(Edn::Int(1));
        assert_eq!(s2.len(), 1);
        assert!(!s2.is_empty());
    }

    #[test]
    fn test_edn_seq_default() {
        let s: EdnSeq<'_> = EdnSeq::default();
        assert!(s.is_empty());
    }

    #[test]
    fn test_edn_seq_from_iter() {
        let s: EdnSeq<'_> = vec![Edn::Int(1), Edn::Int(2)].into_iter().collect();
        assert_eq!(s.len(), 2);
    }

    #[test]
    fn test_edn_seq_ref_into_iter() {
        let s: EdnSeq<'_> = vec![Edn::Int(1)].into();
        let count = (&s).into_iter().count();
        assert_eq!(count, 1);
    }

    #[test]
    fn test_edn_seq_partial_ord() {
        let a: EdnSeq<'_> = vec![Edn::Int(1)].into();
        let b: EdnSeq<'_> = vec![Edn::Int(2)].into();
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Less));
    }
}
