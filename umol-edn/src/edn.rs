//! Core EDN value type, Keyword, Symbol.

use std::borrow::Cow;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::str::FromStr;
use std::{fmt, iter};

use crate::collections::{EdnKeyRef, EdnMap, EdnSeq, EdnSet};
use crate::error::EdnError;
use crate::reader::read_string;

#[cfg(feature = "bignum")]
use bigdecimal::BigDecimal;
#[cfg(feature = "bignum")]
use num_bigint::BigInt;

/// A keyword value (`:name` or `:ns/name`).
#[derive(Clone, Debug, Eq)]
pub struct Keyword<'a>(Cow<'a, str>);

impl<'a> Keyword<'a> {
    pub fn new(name: &'a str) -> Self {
        Keyword(Cow::Borrowed(name))
    }

    pub fn owned(name: String) -> Self {
        Keyword(Cow::Owned(name))
    }

    pub fn namespaced(ns: &str, name: &str) -> Self {
        Keyword(Cow::Owned(format!("{ns}/{name}")))
    }

    pub fn name(&self) -> &str {
        match self.0.rfind('/') {
            Some(pos) if pos > 0 => &self.0[pos + 1..],
            _ => &self.0,
        }
    }

    pub fn namespace(&self) -> Option<&str> {
        match self.0.rfind('/') {
            Some(pos) if pos > 0 => Some(&self.0[..pos]),
            _ => None,
        }
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn into_cow(self) -> Cow<'a, str> {
        self.0
    }

    pub fn into_owned(self) -> Keyword<'static> {
        Keyword(Cow::Owned(self.0.into_owned()))
    }
}

impl PartialEq for Keyword<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl PartialOrd for Keyword<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Keyword<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl Hash for Keyword<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl fmt::Display for Keyword<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, ":{}", self.0)
    }
}

/// A symbol value (`name` or `ns/name`).
#[derive(Clone, Debug, Eq)]
pub struct Symbol<'a>(Cow<'a, str>);

impl<'a> Symbol<'a> {
    pub fn new(name: &'a str) -> Self {
        Symbol(Cow::Borrowed(name))
    }

    pub fn owned(name: String) -> Self {
        Symbol(Cow::Owned(name))
    }

    pub fn namespaced(ns: &str, name: &str) -> Self {
        Symbol(Cow::Owned(format!("{ns}/{name}")))
    }

    pub fn name(&self) -> &str {
        match self.0.rfind('/') {
            Some(pos) if pos > 0 => &self.0[pos + 1..],
            _ => &self.0,
        }
    }

    pub fn namespace(&self) -> Option<&str> {
        match self.0.rfind('/') {
            Some(pos) if pos > 0 => Some(&self.0[..pos]),
            _ => None,
        }
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn into_cow(self) -> Cow<'a, str> {
        self.0
    }

    pub fn into_owned(self) -> Symbol<'static> {
        Symbol(Cow::Owned(self.0.into_owned()))
    }
}

impl PartialEq for Symbol<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl PartialOrd for Symbol<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Symbol<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl Hash for Symbol<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl fmt::Display for Symbol<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Newtype adding `Eq` to `BigDecimal`. The upstream crate omits `Eq` despite
/// `BigDecimal` having no NaN-like values (reflexivity holds). Its `Hash` impl
/// normalizes trailing zeros before hashing, consistent with `PartialEq`, so
/// the `Eq` + `Hash` contract is satisfied.
#[cfg(feature = "bignum")]
#[derive(Clone, Debug, Hash)]
pub struct EdnBigDecimal(BigDecimal);

#[cfg(feature = "bignum")]
impl EdnBigDecimal {
    pub fn new(bd: BigDecimal) -> Self {
        Self(bd)
    }

    pub fn into_inner(self) -> BigDecimal {
        self.0
    }

    pub fn as_inner(&self) -> &BigDecimal {
        &self.0
    }
}

#[cfg(feature = "bignum")]
impl PartialEq for EdnBigDecimal {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

#[cfg(feature = "bignum")]
impl Eq for EdnBigDecimal {}

#[cfg(feature = "bignum")]
impl PartialOrd for EdnBigDecimal {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(feature = "bignum")]
impl Ord for EdnBigDecimal {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

#[cfg(feature = "bignum")]
impl fmt::Display for EdnBigDecimal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

// ---------------------------------------------------------------------------
// Edn enum
// ---------------------------------------------------------------------------

/// Discriminant for cross-variant comparison. List and Vector share the same
/// ordinal: the EDN spec defines sequence equality across container types.
fn variant_ord(v: &Edn<'_>) -> u8 {
    match v {
        Edn::Nil => 0,
        Edn::Bool(_) => 1,
        Edn::Int(_) => 2,
        #[cfg(feature = "bignum")]
        Edn::BigInt(_) => 3,
        Edn::Float(_) => 4,
        #[cfg(feature = "bignum")]
        Edn::BigDecimal(_) => 5,
        Edn::Char(_) => 6,
        Edn::Str(_) => 7,
        Edn::Keyword(_) => 8,
        Edn::Symbol(_) => 9,
        Edn::List(_) | Edn::Vector(_) => 10,
        Edn::Map(_) => 11,
        Edn::Set(_) => 12,
        Edn::Tagged(_, _) => 13,
    }
}

/// An EDN value.
#[derive(Clone, Debug, Default)]
pub enum Edn<'a> {
    #[default]
    Nil,
    Bool(bool),
    Int(i64),
    #[cfg(feature = "bignum")]
    BigInt(BigInt),
    Float(f64),
    #[cfg(feature = "bignum")]
    BigDecimal(EdnBigDecimal),
    Char(char),
    Str(Cow<'a, str>),
    Keyword(Keyword<'a>),
    Symbol(Symbol<'a>),
    List(EdnSeq<'a>),
    Vector(EdnSeq<'a>),
    Map(EdnMap<'a>),
    Set(EdnSet<'a>),
    Tagged(Cow<'a, str>, Box<Edn<'a>>),
}

impl FromStr for Edn<'static> {
    type Err = EdnError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        read_string(s).map(Edn::into_owned)
    }
}

impl<'a> Edn<'a> {
    // --- Constructors ---

    pub fn keyword(name: &'a str) -> Self {
        Edn::Keyword(Keyword::new(name))
    }

    pub fn symbol(name: &'a str) -> Self {
        Edn::Symbol(Symbol::new(name))
    }

    pub fn string(s: &'a str) -> Self {
        Edn::Str(Cow::Borrowed(s))
    }

    // --- Type checks ---

    pub fn is_nil(&self) -> bool {
        matches!(self, Edn::Nil)
    }

    pub fn is_bool(&self) -> bool {
        matches!(self, Edn::Bool(_))
    }

    pub fn is_int(&self) -> bool {
        matches!(self, Edn::Int(_))
    }

    pub fn is_float(&self) -> bool {
        matches!(self, Edn::Float(_))
    }

    #[cfg(feature = "bignum")]
    pub fn is_bigint(&self) -> bool {
        matches!(self, Edn::BigInt(_))
    }

    #[cfg(feature = "bignum")]
    pub fn is_bigdecimal(&self) -> bool {
        matches!(self, Edn::BigDecimal(_))
    }

    pub fn is_char(&self) -> bool {
        matches!(self, Edn::Char(_))
    }

    pub fn is_str(&self) -> bool {
        matches!(self, Edn::Str(_))
    }

    pub fn is_keyword(&self) -> bool {
        matches!(self, Edn::Keyword(_))
    }

    pub fn is_symbol(&self) -> bool {
        matches!(self, Edn::Symbol(_))
    }

    pub fn is_list(&self) -> bool {
        matches!(self, Edn::List(_))
    }

    pub fn is_vector(&self) -> bool {
        matches!(self, Edn::Vector(_))
    }

    pub fn is_map(&self) -> bool {
        matches!(self, Edn::Map(_))
    }

    pub fn is_set(&self) -> bool {
        matches!(self, Edn::Set(_))
    }

    pub fn is_tagged(&self) -> bool {
        matches!(self, Edn::Tagged(_, _))
    }

    // --- Narrowing accessors ---

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Edn::Bool(b) => Some(*b),
            _ => None,
        }
    }

    pub fn as_i64(&self) -> Option<i64> {
        match self {
            Edn::Int(n) => Some(*n),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Edn::Float(f) => Some(*f),
            Edn::Int(n) => Some(*n as f64),
            _ => None,
        }
    }

    #[cfg(feature = "bignum")]
    pub fn as_bigint(&self) -> Option<&BigInt> {
        match self {
            Edn::BigInt(n) => Some(n),
            _ => None,
        }
    }

    #[cfg(feature = "bignum")]
    pub fn as_bigdecimal(&self) -> Option<&BigDecimal> {
        match self {
            Edn::BigDecimal(d) => Some(d.as_inner()),
            _ => None,
        }
    }

    pub fn as_char(&self) -> Option<char> {
        match self {
            Edn::Char(c) => Some(*c),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Edn::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_keyword(&self) -> Option<&Keyword<'a>> {
        match self {
            Edn::Keyword(k) => Some(k),
            _ => None,
        }
    }

    pub fn as_symbol(&self) -> Option<&Symbol<'a>> {
        match self {
            Edn::Symbol(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_list(&self) -> Option<&[Edn<'a>]> {
        match self {
            Edn::List(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_vector(&self) -> Option<&[Edn<'a>]> {
        match self {
            Edn::Vector(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_map(&self) -> Option<&EdnMap<'a>> {
        match self {
            Edn::Map(m) => Some(m),
            _ => None,
        }
    }

    pub fn as_set(&self) -> Option<&EdnSet<'a>> {
        match self {
            Edn::Set(s) => Some(s),
            _ => None,
        }
    }

    // --- Numeric narrowing ---

    pub fn as_u8(&self) -> Option<u8> {
        self.as_i64().and_then(|n| u8::try_from(n).ok())
    }

    pub fn as_u16(&self) -> Option<u16> {
        self.as_i64().and_then(|n| u16::try_from(n).ok())
    }

    pub fn as_u32(&self) -> Option<u32> {
        self.as_i64().and_then(|n| u32::try_from(n).ok())
    }

    pub fn as_u64(&self) -> Option<u64> {
        self.as_i64().and_then(|n| u64::try_from(n).ok())
    }

    pub fn as_i8(&self) -> Option<i8> {
        self.as_i64().and_then(|n| i8::try_from(n).ok())
    }

    pub fn as_i16(&self) -> Option<i16> {
        self.as_i64().and_then(|n| i16::try_from(n).ok())
    }

    pub fn as_i32(&self) -> Option<i32> {
        self.as_i64().and_then(|n| i32::try_from(n).ok())
    }

    // --- Collection access ---

    /// Look up a keyword in a map by its string name.
    pub fn get_keyword(&self, key: &str) -> Option<&Edn<'a>> {
        match self {
            Edn::Map(m) => m.get_ref(EdnKeyRef::keyword(key)),
            _ => None,
        }
    }

    /// Iterate over elements of a vector, list, or set.
    pub fn iter(&self) -> Box<dyn Iterator<Item = &Edn<'a>> + '_> {
        match self {
            Edn::Vector(v) | Edn::List(v) => Box::new(v.iter()),
            Edn::Set(s) => Box::new(s.iter()),
            _ => Box::new(iter::empty()),
        }
    }

    // --- Ownership ---

    pub fn into_owned(self) -> Edn<'static> {
        match self {
            Edn::Nil => Edn::Nil,
            Edn::Bool(b) => Edn::Bool(b),
            Edn::Int(n) => Edn::Int(n),
            Edn::Float(f) => Edn::Float(f),
            Edn::Char(c) => Edn::Char(c),
            #[cfg(feature = "bignum")]
            Edn::BigInt(n) => Edn::BigInt(n),
            #[cfg(feature = "bignum")]
            Edn::BigDecimal(d) => Edn::BigDecimal(d),
            Edn::Str(s) => Edn::Str(Cow::Owned(s.into_owned())),
            Edn::Keyword(k) => Edn::Keyword(k.into_owned()),
            Edn::Symbol(s) => Edn::Symbol(s.into_owned()),
            Edn::List(v) => Edn::List(v.into_owned()),
            Edn::Vector(v) => Edn::Vector(v.into_owned()),
            Edn::Map(m) => Edn::Map(m.into_owned()),
            Edn::Set(s) => Edn::Set(s.into_owned()),
            Edn::Tagged(tag, inner) => {
                Edn::Tagged(Cow::Owned(tag.into_owned()), Box::new(inner.into_owned()))
            }
        }
    }
}

// Manual PartialEq: use f64::total_cmp for Float
impl PartialEq for Edn<'_> {
    fn eq(&self, other: &Self) -> bool {
        if variant_ord(self) != variant_ord(other) {
            return false;
        }
        match (self, other) {
            (Edn::Nil, Edn::Nil) => true,
            (Edn::Bool(a), Edn::Bool(b)) => a == b,
            (Edn::Int(a), Edn::Int(b)) => a == b,
            #[cfg(feature = "bignum")]
            (Edn::BigInt(a), Edn::BigInt(b)) => a == b,
            (Edn::Float(a), Edn::Float(b)) => a.to_bits() == b.to_bits(),
            #[cfg(feature = "bignum")]
            (Edn::BigDecimal(a), Edn::BigDecimal(b)) => a == b,
            (Edn::Char(a), Edn::Char(b)) => a == b,
            (Edn::Str(a), Edn::Str(b)) => a == b,
            (Edn::Keyword(a), Edn::Keyword(b)) => a == b,
            (Edn::Symbol(a), Edn::Symbol(b)) => a == b,
            (Edn::List(a) | Edn::Vector(a), Edn::List(b) | Edn::Vector(b)) => a == b,
            (Edn::Map(a), Edn::Map(b)) => a == b,
            (Edn::Set(a), Edn::Set(b)) => a == b,
            (Edn::Tagged(ta, va), Edn::Tagged(tb, vb)) => ta == tb && va == vb,
            _ => {
                debug_assert!(false, "variant_ord should prevent cross-variant pairing");
                false
            }
        }
    }
}

impl Eq for Edn<'_> {}

impl PartialOrd for Edn<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Edn<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        let d = variant_ord(self).cmp(&variant_ord(other));
        if d != Ordering::Equal {
            return d;
        }
        match (self, other) {
            (Edn::Nil, Edn::Nil) => Ordering::Equal,
            (Edn::Bool(a), Edn::Bool(b)) => a.cmp(b),
            (Edn::Int(a), Edn::Int(b)) => a.cmp(b),
            #[cfg(feature = "bignum")]
            (Edn::BigInt(a), Edn::BigInt(b)) => a.cmp(b),
            (Edn::Float(a), Edn::Float(b)) => a.total_cmp(b),
            #[cfg(feature = "bignum")]
            (Edn::BigDecimal(a), Edn::BigDecimal(b)) => a.cmp(b),
            (Edn::Char(a), Edn::Char(b)) => a.cmp(b),
            (Edn::Str(a), Edn::Str(b)) => a.cmp(b),
            (Edn::Keyword(a), Edn::Keyword(b)) => a.cmp(b),
            (Edn::Symbol(a), Edn::Symbol(b)) => a.cmp(b),
            (Edn::List(a) | Edn::Vector(a), Edn::List(b) | Edn::Vector(b)) => a.cmp(b),
            (Edn::Map(a), Edn::Map(b)) => a.cmp(b),
            (Edn::Set(a), Edn::Set(b)) => a.cmp(b),
            (Edn::Tagged(ta, va), Edn::Tagged(tb, vb)) => ta.cmp(tb).then_with(|| va.cmp(vb)),
            _ => {
                debug_assert!(false, "variant_ord should prevent cross-variant pairing");
                Ordering::Equal
            }
        }
    }
}

impl Hash for Edn<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        variant_ord(self).hash(state);
        match self {
            Edn::Nil => {}
            Edn::Bool(b) => b.hash(state),
            Edn::Int(n) => n.hash(state),
            #[cfg(feature = "bignum")]
            Edn::BigInt(n) => n.hash(state),
            Edn::Float(f) => f.to_bits().hash(state),
            #[cfg(feature = "bignum")]
            Edn::BigDecimal(d) => d.hash(state),
            Edn::Char(c) => c.hash(state),
            Edn::Str(s) => s.hash(state),
            Edn::Keyword(k) => k.hash(state),
            Edn::Symbol(s) => s.hash(state),
            Edn::List(v) | Edn::Vector(v) => v.hash(state),
            Edn::Map(m) => m.hash(state),
            Edn::Set(s) => s.hash(state),
            Edn::Tagged(tag, inner) => {
                tag.hash(state);
                inner.hash(state);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::hash::{DefaultHasher, Hash, Hasher};

    fn hash_of(v: &Edn<'_>) -> u64 {
        let mut h = DefaultHasher::new();
        v.hash(&mut h);
        h.finish()
    }

    #[test]
    fn test_edn_float_nan_payloads_distinguished() {
        // total_cmp follows IEEE 754 total ordering: distinct NaN payloads are not equal.
        // Hash uses to_bits, so Eq and Hash remain consistent.
        let nan1 = Edn::Float(f64::from_bits(0x7FF8_0000_0000_0001));
        let nan2 = Edn::Float(f64::from_bits(0x7FF8_0000_0000_0002));
        assert_ne!(nan1, nan2);
        assert_ne!(hash_of(&nan1), hash_of(&nan2));
    }

    #[test]
    fn test_edn_float_canonical_nan_eq_hash() {
        let nan1 = Edn::Float(f64::NAN);
        let nan2 = Edn::Float(f64::NAN);
        assert_eq!(nan1, nan2);
        assert_eq!(hash_of(&nan1), hash_of(&nan2));
    }

    #[test]
    fn test_edn_float_positive_negative_zero() {
        let pos = Edn::Float(0.0);
        let neg = Edn::Float(-0.0);
        assert_ne!(pos, neg);
        assert_ne!(hash_of(&pos), hash_of(&neg));
    }

    /// Guards the assumption that BigDecimal's PartialEq and Hash are
    /// consistent across representations of the same value. If this breaks
    /// on a bigdecimal upgrade, EdnBigDecimal's Eq impl is unsound.
    #[cfg(feature = "bignum")]
    #[test]
    fn test_edn_bigdecimal_eq_hash_across_scales() {
        use bigdecimal::BigDecimal;
        use std::str::FromStr;

        let a = EdnBigDecimal::new(BigDecimal::from_str("1.0").unwrap());
        let b = EdnBigDecimal::new(BigDecimal::from_str("1.00").unwrap());
        let c = EdnBigDecimal::new(BigDecimal::from_str("1.000").unwrap());

        assert_eq!(a, b);
        assert_eq!(b, c);
        assert_eq!(a, a); // reflexivity

        fn hash_bd(v: &EdnBigDecimal) -> u64 {
            let mut h = DefaultHasher::new();
            v.hash(&mut h);
            h.finish()
        }
        assert_eq!(hash_bd(&a), hash_bd(&b));
        assert_eq!(hash_bd(&b), hash_bd(&c));
    }

    // -- Default, FromStr --

    #[test]
    fn test_edn_default() {
        assert_eq!(Edn::default(), Edn::Nil);
    }

    #[test]
    fn test_edn_from_str() {
        use std::str::FromStr;
        let v: Edn<'static> = Edn::from_str("[1 2]").unwrap();
        assert!(v.is_vector());
        assert!(Edn::from_str("[invalid").is_err());
    }

    // -- Symbol methods --

    #[test]
    fn test_symbol_owned() {
        let s = Symbol::owned("foo".into());
        assert_eq!(s.as_str(), "foo");
    }

    #[test]
    fn test_symbol_name_namespace() {
        let s = Symbol::namespaced("ns", "bar");
        assert_eq!(s.name(), "bar");
        assert_eq!(s.namespace(), Some("ns"));

        let s = Symbol::new("simple");
        assert_eq!(s.name(), "simple");
        assert_eq!(s.namespace(), None);
    }

    #[test]
    fn test_symbol_into_owned() {
        let s = Symbol::new("hello");
        let owned: Symbol<'static> = s.into_owned();
        assert_eq!(owned.as_str(), "hello");
    }

    #[test]
    fn test_symbol_ord() {
        let a = Symbol::new("alpha");
        let b = Symbol::new("beta");
        assert!(a < b);
        assert_eq!(a.partial_cmp(&b), Some(Ordering::Less));
    }

    // -- Keyword into_owned, ordering --

    #[test]
    fn test_keyword_into_owned() {
        let k = Keyword::new("key");
        let owned: Keyword<'static> = k.into_owned();
        assert_eq!(owned.as_str(), "key");
    }

    #[test]
    fn test_keyword_ord() {
        let a = Keyword::new("alpha");
        let b = Keyword::new("beta");
        assert!(a < b);
        assert_eq!(a.partial_cmp(&b), Some(Ordering::Less));
    }

    // -- Type checks --

    #[test]
    fn test_edn_is_checks() {
        assert!(Edn::Nil.is_nil());
        assert!(Edn::Bool(true).is_bool());
        assert!(Edn::Int(1).is_int());
        assert!(Edn::Float(1.0).is_float());
        assert!(Edn::Char('x').is_char());
        assert!(Edn::Str(Cow::Borrowed("s")).is_str());
        assert!(Edn::Keyword(Keyword::new("k")).is_keyword());
        assert!(Edn::Symbol(Symbol::new("s")).is_symbol());
        assert!(Edn::List(vec![].into()).is_list());
        assert!(Edn::Vector(vec![].into()).is_vector());
        assert!(Edn::Map(EdnMap::new()).is_map());
        assert!(Edn::Set(EdnSet::new()).is_set());
        assert!(Edn::Tagged(Cow::Borrowed("t"), Box::new(Edn::Nil)).is_tagged());

        assert!(!Edn::Nil.is_bool());
        assert!(!Edn::Int(1).is_float());
    }

    // -- Accessors --

    #[test]
    fn test_edn_as_bool() {
        assert_eq!(Edn::Bool(true).as_bool(), Some(true));
        assert_eq!(Edn::Nil.as_bool(), None);
    }

    #[test]
    fn test_edn_as_i64() {
        assert_eq!(Edn::Int(7).as_i64(), Some(7));
        assert_eq!(Edn::Bool(true).as_i64(), None);
    }

    #[test]
    fn test_edn_as_f64() {
        assert_eq!(Edn::Float(3.14).as_f64(), Some(3.14));
        assert_eq!(Edn::Int(5).as_f64(), Some(5.0));
        assert_eq!(Edn::Bool(true).as_f64(), None);
    }

    #[test]
    fn test_edn_as_char() {
        assert_eq!(Edn::Char('z').as_char(), Some('z'));
        assert_eq!(Edn::Nil.as_char(), None);
    }

    #[test]
    fn test_edn_as_str() {
        assert_eq!(Edn::Str(Cow::Borrowed("hi")).as_str(), Some("hi"));
        assert_eq!(Edn::Nil.as_str(), None);
    }

    #[test]
    fn test_edn_as_keyword() {
        let k = Edn::Keyword(Keyword::new("k"));
        assert!(k.as_keyword().is_some());
        assert!(Edn::Nil.as_keyword().is_none());
    }

    #[test]
    fn test_edn_as_symbol() {
        let s = Edn::Symbol(Symbol::new("s"));
        assert!(s.as_symbol().is_some());
        assert!(Edn::Nil.as_symbol().is_none());
    }

    #[test]
    fn test_edn_as_list() {
        let l = Edn::List(vec![Edn::Int(1)].into());
        assert!(l.as_list().is_some());
        assert!(Edn::Nil.as_list().is_none());
    }

    #[test]
    fn test_edn_as_vector() {
        let v = Edn::Vector(vec![Edn::Int(1)].into());
        assert!(v.as_vector().is_some());
        assert!(Edn::Nil.as_vector().is_none());
    }

    #[test]
    fn test_edn_as_map() {
        let m = Edn::Map(EdnMap::new());
        assert!(m.as_map().is_some());
        assert!(Edn::Nil.as_map().is_none());
    }

    #[test]
    fn test_edn_as_set() {
        let s = Edn::Set(EdnSet::new());
        assert!(s.as_set().is_some());
        assert!(Edn::Nil.as_set().is_none());
    }

    // -- Numeric narrowing --

    #[test]
    fn test_edn_numeric_narrowing() {
        let v = Edn::Int(200);
        assert_eq!(v.as_u8(), Some(200));
        assert_eq!(v.as_u16(), Some(200));
        assert_eq!(v.as_u32(), Some(200));
        assert_eq!(v.as_u64(), Some(200));
        assert_eq!(v.as_i8(), None); // 200 out of range
        assert_eq!(v.as_i16(), Some(200));
        assert_eq!(v.as_i32(), Some(200));

        assert_eq!(Edn::Int(-1).as_u64(), None);
    }

    // -- get_keyword on non-map --

    #[test]
    fn test_edn_get_keyword_non_map() {
        assert!(Edn::Nil.get_keyword("x").is_none());
    }

    // -- iter --

    #[test]
    fn test_edn_iter_set() {
        let mut s = EdnSet::new();
        s.insert(Edn::Int(1));
        s.insert(Edn::Int(2));
        let edn = Edn::Set(s);
        assert_eq!(edn.iter().count(), 2);
    }

    #[test]
    fn test_edn_iter_non_collection() {
        assert_eq!(Edn::Int(5).iter().count(), 0);
    }

    // -- into_owned --

    #[test]
    fn test_edn_into_owned() {
        let values: Vec<Edn<'_>> = vec![
            Edn::Nil,
            Edn::Bool(true),
            Edn::Int(7),
            Edn::Float(3.14),
            Edn::Char('x'),
            Edn::Str(Cow::Borrowed("hello")),
            Edn::Keyword(Keyword::new("k")),
            Edn::Symbol(Symbol::new("s")),
            Edn::List(vec![Edn::Int(1)].into()),
            Edn::Vector(vec![Edn::Int(2)].into()),
            Edn::Map(EdnMap::new()),
            Edn::Set(EdnSet::new()),
            Edn::Tagged(Cow::Borrowed("tag"), Box::new(Edn::Nil)),
        ];
        for v in values {
            let cloned = v.clone();
            let owned: Edn<'static> = v.into_owned();
            assert_eq!(owned, cloned);
        }
    }

    // -- Ord across variants and within variants --

    #[test]
    fn test_edn_ord_cross_variant() {
        assert!(Edn::Nil < Edn::Bool(false));
        assert!(Edn::Bool(true) < Edn::Int(0));
        assert!(Edn::Int(0) < Edn::Float(0.0));
    }

    #[test]
    fn test_edn_ord_same_variant() {
        assert!(Edn::Bool(false) < Edn::Bool(true));
        assert!(Edn::Int(1) < Edn::Int(2));
        assert!(Edn::Float(1.0) < Edn::Float(2.0));
        assert!(Edn::Char('a') < Edn::Char('b'));
        assert!(Edn::Str(Cow::Borrowed("a")) < Edn::Str(Cow::Borrowed("b")));
        assert!(Edn::Keyword(Keyword::new("a")) < Edn::Keyword(Keyword::new("b")));
        assert!(Edn::Symbol(Symbol::new("a")) < Edn::Symbol(Symbol::new("b")));
    }

    #[test]
    fn test_edn_partial_ord() {
        assert_eq!(
            Edn::Int(1).partial_cmp(&Edn::Int(2)),
            Some(Ordering::Less)
        );
    }

    #[test]
    fn test_edn_ord_list_vector() {
        let l1 = Edn::List(vec![Edn::Int(1)].into());
        let l2 = Edn::List(vec![Edn::Int(2)].into());
        assert!(l1 < l2);

        let v1 = Edn::Vector(vec![Edn::Int(1)].into());
        let v2 = Edn::Vector(vec![Edn::Int(2)].into());
        assert!(v1 < v2);
    }

    #[test]
    fn test_edn_ord_tagged() {
        let t1 = Edn::Tagged(Cow::Borrowed("a"), Box::new(Edn::Int(1)));
        let t2 = Edn::Tagged(Cow::Borrowed("b"), Box::new(Edn::Int(1)));
        assert!(t1 < t2);

        let t3 = Edn::Tagged(Cow::Borrowed("a"), Box::new(Edn::Int(2)));
        assert!(t1 < t3);
    }

    // -- PartialEq across variants --

    #[test]
    fn test_edn_eq_cross_variant_false() {
        assert_ne!(Edn::Int(0), Edn::Float(0.0));
        assert_ne!(Edn::Nil, Edn::Bool(false));
        assert_ne!(
            Edn::Str(Cow::Borrowed("foo")),
            Edn::Keyword(Keyword::new("foo"))
        );
    }

    // -- PartialEq within each variant --

    #[test]
    fn test_edn_eq_tagged() {
        let a = Edn::Tagged(Cow::Borrowed("t"), Box::new(Edn::Int(1)));
        let b = Edn::Tagged(Cow::Borrowed("t"), Box::new(Edn::Int(1)));
        let c = Edn::Tagged(Cow::Borrowed("t"), Box::new(Edn::Int(2)));
        assert_eq!(a, b);
        assert_ne!(a, c);
    }
}
