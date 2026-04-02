//! Core EDN value type, Keyword, Symbol.

use std::borrow::Cow;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::{fmt, iter};

use crate::collections::{EdnMap, EdnSeq, EdnSet};

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

// ---------------------------------------------------------------------------
// Edn enum
// ---------------------------------------------------------------------------

/// Discriminant ordering for Edn variants (used in Ord impl).
fn variant_ord(v: &Edn<'_>) -> u8 {
    match v {
        Edn::Nil => 0,
        Edn::Bool(_) => 1,
        Edn::Int(_) => 2,
        Edn::Float(_) => 3,
        Edn::Char(_) => 4,
        Edn::Str(_) => 5,
        Edn::Keyword(_) => 6,
        Edn::Symbol(_) => 7,
        Edn::List(_) => 8,
        Edn::Vector(_) => 9,
        Edn::Map(_) => 10,
        Edn::Set(_) => 11,
        Edn::Tagged(_, _) => 12,
    }
}

/// An EDN value.
#[derive(Clone, Debug)]
pub enum Edn<'a> {
    Nil,
    Bool(bool),
    Int(i64),
    Float(f64),
    Char(char),
    Str(Cow<'a, str>),
    Keyword(Keyword<'a>),
    Symbol(Symbol<'a>),
    List(EdnSeq<'a>),
    Vector(EdnSeq<'a>),
    Map(EdnMap<'a>),
    Set(EdnSet<'a>),
    Tagged(String, Box<Edn<'a>>),
}

impl Default for Edn<'_> {
    fn default() -> Self {
        Edn::Nil
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
    pub fn get(&self, key: &str) -> Option<&Edn<'a>> {
        match self {
            Edn::Map(m) => m.get(&Edn::Keyword(Keyword::owned(key.to_string()))),
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
            Edn::Str(s) => Edn::Str(Cow::Owned(s.into_owned())),
            Edn::Keyword(k) => Edn::Keyword(k.into_owned()),
            Edn::Symbol(s) => Edn::Symbol(s.into_owned()),
            Edn::List(v) => Edn::List(v.into_owned()),
            Edn::Vector(v) => Edn::Vector(v.into_owned()),
            Edn::Map(m) => Edn::Map(m.into_owned()),
            Edn::Set(s) => Edn::Set(s.into_owned()),
            Edn::Tagged(tag, inner) => Edn::Tagged(tag, Box::new(inner.into_owned())),
        }
    }
}

// Manual PartialEq: use f64::total_cmp for Float
impl PartialEq for Edn<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
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
            (Edn::Float(a), Edn::Float(b)) => a.total_cmp(b),
            (Edn::Char(a), Edn::Char(b)) => a.cmp(b),
            (Edn::Str(a), Edn::Str(b)) => a.cmp(b),
            (Edn::Keyword(a), Edn::Keyword(b)) => a.cmp(b),
            (Edn::Symbol(a), Edn::Symbol(b)) => a.cmp(b),
            (Edn::List(a), Edn::List(b)) => a.cmp(b),
            (Edn::Vector(a), Edn::Vector(b)) => a.cmp(b),
            (Edn::Map(a), Edn::Map(b)) => a.cmp(b),
            (Edn::Set(a), Edn::Set(b)) => a.cmp(b),
            (Edn::Tagged(ta, va), Edn::Tagged(tb, vb)) => ta.cmp(tb).then_with(|| va.cmp(vb)),
            _ => unreachable!(),
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
            Edn::Float(f) => f.to_bits().hash(state),
            Edn::Char(c) => c.hash(state),
            Edn::Str(s) => s.hash(state),
            Edn::Keyword(k) => k.hash(state),
            Edn::Symbol(s) => s.hash(state),
            Edn::List(v) => v.hash(state),
            Edn::Vector(v) => v.hash(state),
            Edn::Map(m) => m.hash(state),
            Edn::Set(s) => s.hash(state),
            Edn::Tagged(tag, inner) => {
                tag.hash(state);
                inner.hash(state);
            }
        }
    }
}
