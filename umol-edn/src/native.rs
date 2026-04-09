//! Native conversion traits between Rust types and EDN values.
//!
//! This module defines the primary deserialization and serialization API for
//! umol-edn. It is independent of `serde::Deserialize` / `serde::Serialize`
//! and operates directly on `Edn<'de>` values, which means it can express
//! the full EDN data model — keywords, symbols, sets, tagged literals,
//! arbitrary-precision numbers — without the wrapper-newtype gymnastics
//! that the serde data model forces on EDN-specific features.
//!
//! ## Design
//!
//! The two traits are deliberately minimal:
//!
//! - [`FromEdn`] takes a borrowed `&Edn<'de>` (a parsed tree) and produces
//!   a value. It carries a lifetime parameter so implementations can borrow
//!   string and key references from the source buffer when zero-copy is
//!   wanted.
//! - [`ToEdn`] takes `&self` and produces an `Edn<'_>` borrowing from the
//!   value where possible.
//!
//! ## The `from_edn_str` escape hatch
//!
//! [`FromEdn::from_edn_str`] has a default implementation that parses to a
//! tree and then calls `from_edn`. This is the path that supports the full
//! EDN data model uniformly. For hot types where parse-time tree
//! construction is the bottleneck, the trait permits an override that fuses
//! parsing and deserialization in a single pass — the same architectural
//! shape as the legacy `serde::Deserialize`-based streaming path, but
//! per-type and EDN-native.
//!
//! Implementations that override `from_edn_str` must produce a value
//! equivalent to what `Self::from_edn(&read_string(input)?)` would have
//! produced. This is the contract that lets callers reach for the fast path
//! without changing semantics.
//!
//! ## Lifetimes
//!
//! `FromEdn<'de>` follows the same convention as `serde::Deserialize<'de>`:
//! `'de` is the lifetime of the source the implementation is allowed to
//! borrow from. A type that wants to borrow string slices declares
//! `MyType<'a>: FromEdn<'a>`. A type that always owns its strings is
//! `MyType: for<'de> FromEdn<'de>`.

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::hash::Hash;

use crate::collections::{EdnKeyRef, EdnMap, EdnSeq, EdnSet};
use crate::edn::{Edn, Keyword, Symbol};
use crate::error::EdnError;
use crate::reader::read_string;

#[cfg(feature = "bignum")]
use crate::edn::EdnBigDecimal;
#[cfg(feature = "bignum")]
use bigdecimal::BigDecimal;
#[cfg(feature = "bignum")]
use num_bigint::BigInt;

// Variant discriminator strings come from `Edn::kind()`.

/// Build `Self` from an EDN value.
///
/// This trait is the primary deserialization entry point for umol-edn types.
/// It replaces direct use of `serde::Deserialize` for EDN-native types,
/// allowing access to the full EDN data model (keywords, symbols, sets,
/// tagged literals, big numbers) without serde wrapper indirection.
///
/// # Implementing
///
/// At minimum, implementations must provide [`from_edn`]. The default
/// [`from_edn_str`] implementation will parse EDN source to a tree and
/// dispatch through `from_edn`, which is correct for all types and matches
/// MOL parser performance for typical molecule DSL inputs.
///
/// Performance-critical types may override [`from_edn_str`] to fuse parsing
/// and deserialization in a single pass. The override is required to be
/// observationally equivalent to the default — same value, same errors —
/// for any input.
///
/// # Lifetime
///
/// The `'de` parameter is the lifetime of the data the implementation may
/// borrow from. Owned types (`String` fields, no borrowing) implement
/// `FromEdn` for any `'de`, typically expressed as
/// `impl<'de> FromEdn<'de> for MyType`. Borrowing types
/// (`Cow<'a, str>` fields, `&'a str` fields) parameterize on `'de`:
/// `impl<'de> FromEdn<'de> for MyType<'de>`.
///
/// [`from_edn`]: FromEdn::from_edn
/// [`from_edn_str`]: FromEdn::from_edn_str
pub trait FromEdn<'de>: Sized {
    /// Build `Self` by walking a pre-parsed `Edn` tree.
    ///
    /// Implementations should inspect the tree and report a precise
    /// `EdnError` on shape mismatch. Borrowing string fields from the
    /// tree is permitted via the `'de` lifetime.
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError>;

    /// Build `Self` directly from EDN source text.
    ///
    /// The default implementation calls [`read_string`] to materialize a
    /// tree, then dispatches to [`from_edn`]. This is the canonical path
    /// and supports the full EDN data model uniformly.
    ///
    /// Override this method when single-pass parser-deserializer fusion is
    /// required for latency reasons. The override must be observationally
    /// equivalent to the default — same return value, same error variant —
    /// for every input.
    ///
    /// [`from_edn`]: FromEdn::from_edn
    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        let tree = read_string(input)?;
        Self::from_edn(&tree)
    }
}

/// Builder helper for `FromEdn` impls that consume keyword-keyed maps.
///
/// `EdnMapHelper` tracks which keys have been read so that strict-mode
/// [`finalize`](Self::finalize) can flag unconsumed keys, and it carries a
/// path prefix used in error messages for nested structures.
///
/// The helper assumes maps use keyword keys, which is the convention for the
/// molecule DSL and other umol-edn consumers. Non-keyword keys in the map are
/// ignored by `finalize`'s unknown-key check.
///
/// # Path tracking
///
/// `path` segments are pushed by callers wrapping nested helpers — e.g. an
/// outer impl can call `EdnMapHelper::with_path(inner_map, vec![":atoms".into()])`
/// when descending into a nested map field. Errors raised by the helper itself
/// (`MissingField`, `UnknownField`) carry that path. Errors raised by inner
/// `T::from_edn` calls do not currently propagate the path; that requires a
/// trait-level context parameter and is left for a follow-up.
pub struct EdnMapHelper<'m, 'de: 'm> {
    map: &'m EdnMap<'de>,
    path: Vec<String>,
    consumed: HashSet<String>,
}

impl<'m, 'de: 'm> EdnMapHelper<'m, 'de> {
    /// Create a helper rooted at `<root>` (empty path).
    pub fn new(map: &'m EdnMap<'de>) -> Self {
        Self {
            map,
            path: Vec::new(),
            consumed: HashSet::with_capacity(map.len()),
        }
    }

    /// Create a helper with an explicit path prefix for nested error messages.
    pub fn with_path(map: &'m EdnMap<'de>, path: Vec<String>) -> Self {
        Self {
            map,
            path,
            consumed: HashSet::with_capacity(map.len()),
        }
    }

    /// Read a required keyword-keyed field. Errors with `MissingField` if
    /// absent, or with whatever variant `T::from_edn` returns on the value.
    pub fn required<T: FromEdn<'de>>(&mut self, key: &str) -> Result<T, EdnError> {
        let value =
            self.map
                .get_ref(EdnKeyRef::keyword(key))
                .ok_or_else(|| EdnError::MissingField {
                    key: key.to_string(),
                    path: self.path.clone(),
                })?;
        self.consumed.insert(key.to_string());
        T::from_edn(value)
    }

    /// Read an optional keyword-keyed field. Returns `Ok(None)` if absent.
    pub fn optional<T: FromEdn<'de>>(&mut self, key: &str) -> Result<Option<T>, EdnError> {
        match self.map.get_ref(EdnKeyRef::keyword(key)) {
            Some(value) => {
                self.consumed.insert(key.to_string());
                T::from_edn(value).map(Some)
            }
            None => Ok(None),
        }
    }

    /// Strict-mode close: error with `UnknownField` if any keyword key in the
    /// map was not read via `required` or `optional`.
    pub fn finalize(self) -> Result<(), EdnError> {
        for (k, _) in self.map.iter() {
            if let Edn::Keyword(kw) = k {
                if !self.consumed.contains(kw.as_str()) {
                    return Err(EdnError::UnknownField {
                        key: kw.as_str().to_string(),
                        path: self.path,
                    });
                }
            }
        }
        Ok(())
    }
}

/// Convert `&self` into an EDN value.
///
/// This trait is the primary serialization entry point for umol-edn types.
/// It replaces direct use of `serde::Serialize` for EDN-native types,
/// allowing emission of EDN-specific constructs (keywords, sets, tagged
/// literals) without serde wrapper indirection.
///
/// # Borrowing
///
/// `to_edn` returns `Edn<'_>` with the lifetime bound to `&self`, so
/// implementations can produce zero-copy `Cow::Borrowed` strings for
/// fields that already live in the value. The returned tree is otherwise
/// owned and may be passed to a formatter, written to a string, or
/// returned to a caller for further processing.
///
/// # Infallibility
///
/// `to_edn` is infallible. Types whose serialized form depends on
/// validation should perform that validation at construction time so that
/// every reachable `&self` is serializable. Types that need fallible
/// serialization should add their own `try_to_edn` method until a concrete
/// need establishes a trait-level contract for it.
pub trait ToEdn {
    /// Produce an `Edn` value representing `self`.
    ///
    /// The returned value may borrow string and key references from
    /// `&self`. Implementations should prefer borrowed `Cow` variants over
    /// owned ones to keep allocation cost low for round-trip workflows.
    fn to_edn(&self) -> Edn<'_>;
}

impl<'de> FromEdn<'de> for bool {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Bool(b) => Ok(*b),
            other => Err(EdnError::TypeMismatch {
                expected: "bool",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for bool {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Bool(*self)
    }
}

/// FromEdn/ToEdn for integer types that always round-trip through `i64`.
///
/// Skips `i128`, `u64`, `u128`, `usize`, `isize` — those need either bignum
/// support or platform-aware fallbacks and are deferred until a concrete need
/// arises.
macro_rules! impl_int {
    ($($t:ty),* $(,)?) => {
        $(
            impl<'de> FromEdn<'de> for $t {
                fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
                    match edn {
                        Edn::Int(n) => <$t>::try_from(*n).map_err(|_| {
                            EdnError::OutOfRange {
                                value: n.to_string(),
                                target: stringify!($t),
                                path: Vec::new(),
                            }
                        }),
                        other => Err(EdnError::TypeMismatch {
                            expected: "int",
                            got: other.kind(),
                            path: Vec::new(),
                        }),
                    }
                }
            }

            impl ToEdn for $t {
                fn to_edn(&self) -> Edn<'_> {
                    Edn::Int(i64::from(*self))
                }
            }
        )*
    };
}

impl_int!(i8, i16, i32, i64, u8, u16, u32);

impl<'de> FromEdn<'de> for f64 {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Float(f) => Ok(*f),
            other => Err(EdnError::TypeMismatch {
                expected: "float",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for f64 {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Float(*self)
    }
}

impl<'de> FromEdn<'de> for f32 {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Float(f) => Ok(*f as f32),
            other => Err(EdnError::TypeMismatch {
                expected: "float",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for f32 {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Float(f64::from(*self))
    }
}

impl<'de> FromEdn<'de> for char {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Char(c) => Ok(*c),
            other => Err(EdnError::TypeMismatch {
                expected: "char",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for char {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Char(*self)
    }
}

impl<'de> FromEdn<'de> for String {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Str(s) => Ok(s.to_string()),
            other => Err(EdnError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for String {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Str(Cow::Borrowed(self.as_str()))
    }
}

impl ToEdn for str {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Str(Cow::Borrowed(self))
    }
}

/// `Cow<'de, str>` preserves zero-copy borrowing when the source string
/// itself was borrowed from the parsed buffer.
impl<'de> FromEdn<'de> for Cow<'de, str> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Str(s) => Ok(s.clone()),
            other => Err(EdnError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for Cow<'a, str> {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Str(Cow::Borrowed(self.as_ref()))
    }
}

impl<'de, T: FromEdn<'de>> FromEdn<'de> for Option<T> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Nil => Ok(None),
            other => T::from_edn(other).map(Some),
        }
    }
}

impl<T: ToEdn> ToEdn for Option<T> {
    fn to_edn(&self) -> Edn<'_> {
        match self {
            Some(t) => t.to_edn(),
            None => Edn::Nil,
        }
    }
}

impl<'de, T: FromEdn<'de>> FromEdn<'de> for Vec<T> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        let seq = match edn {
            Edn::Vector(s) | Edn::List(s) => s,
            other => {
                return Err(EdnError::TypeMismatch {
                    expected: "vector",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        seq.iter().map(T::from_edn).collect()
    }
}

impl<T: ToEdn> ToEdn for Vec<T> {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Vector(EdnSeq::from(
            self.iter().map(ToEdn::to_edn).collect::<Vec<_>>(),
        ))
    }
}

impl<'de, T: FromEdn<'de>, const N: usize> FromEdn<'de> for [T; N] {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        let seq = match edn {
            Edn::Vector(s) | Edn::List(s) => s,
            other => {
                return Err(EdnError::TypeMismatch {
                    expected: "vector",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        if seq.len() != N {
            return Err(EdnError::OutOfRange {
                value: seq.len().to_string(),
                target: "fixed-length array",
                path: Vec::new(),
            });
        }
        let parsed: Vec<T> = seq.iter().map(T::from_edn).collect::<Result<_, _>>()?;
        parsed
            .try_into()
            .map_err(|_| EdnError::Custom("array length conversion failed".into()))
    }
}

impl<T: ToEdn, const N: usize> ToEdn for [T; N] {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Vector(EdnSeq::from(
            self.iter().map(ToEdn::to_edn).collect::<Vec<_>>(),
        ))
    }
}

/// Tuple impls for arities 2 and 3 — covers coordinate triples and pairs.
/// Higher arities are not needed by the molecule DSL and can be added on
/// demand.
macro_rules! impl_tuple {
    ($($n:tt $t:ident),+) => {
        impl<'de, $($t: FromEdn<'de>),+> FromEdn<'de> for ($($t,)+) {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
                let seq = match edn {
                    Edn::Vector(s) | Edn::List(s) => s,
                    other => {
                        return Err(EdnError::TypeMismatch {
                            expected: "vector",
                            got: other.kind(),
                            path: Vec::new(),
                        });
                    }
                };
                const ARITY: usize = [$($n),+].len();
                if seq.len() != ARITY {
                    return Err(EdnError::OutOfRange {
                        value: seq.len().to_string(),
                        target: "tuple arity",
                        path: Vec::new(),
                    });
                }
                Ok(($($t::from_edn(&seq[$n])?,)+))
            }
        }

        impl<$($t: ToEdn),+> ToEdn for ($($t,)+) {
            fn to_edn(&self) -> Edn<'_> {
                Edn::Vector(EdnSeq::from(vec![$(self.$n.to_edn()),+]))
            }
        }
    };
}

impl_tuple!(0 A, 1 B);
impl_tuple!(0 A, 1 B, 2 C);
impl_tuple!(0 A, 1 B, 2 C, 3 D);

impl<'de, K, V> FromEdn<'de> for HashMap<K, V>
where
    K: FromEdn<'de> + Eq + Hash,
    V: FromEdn<'de>,
{
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        let map = match edn {
            Edn::Map(m) => m,
            other => {
                return Err(EdnError::TypeMismatch {
                    expected: "map",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        map.iter()
            .map(|(k, v)| Ok((K::from_edn(k)?, V::from_edn(v)?)))
            .collect()
    }
}

impl<K, V> ToEdn for HashMap<K, V>
where
    K: ToEdn,
    V: ToEdn,
{
    fn to_edn(&self) -> Edn<'_> {
        let mut m = EdnMap::with_capacity(self.len());
        for (k, v) in self {
            m.insert(k.to_edn(), v.to_edn());
        }
        Edn::Map(m)
    }
}

impl<'de, K, V> FromEdn<'de> for BTreeMap<K, V>
where
    K: FromEdn<'de> + Ord,
    V: FromEdn<'de>,
{
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        let map = match edn {
            Edn::Map(m) => m,
            other => {
                return Err(EdnError::TypeMismatch {
                    expected: "map",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        map.iter()
            .map(|(k, v)| Ok((K::from_edn(k)?, V::from_edn(v)?)))
            .collect()
    }
}

impl<K, V> ToEdn for BTreeMap<K, V>
where
    K: ToEdn,
    V: ToEdn,
{
    fn to_edn(&self) -> Edn<'_> {
        let mut m = EdnMap::with_capacity(self.len());
        for (k, v) in self {
            m.insert(k.to_edn(), v.to_edn());
        }
        Edn::Map(m)
    }
}

impl<'de, T> FromEdn<'de> for HashSet<T>
where
    T: FromEdn<'de> + Eq + Hash,
{
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        let set = match edn {
            Edn::Set(s) => s,
            other => {
                return Err(EdnError::TypeMismatch {
                    expected: "set",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        set.iter().map(T::from_edn).collect()
    }
}

impl<T: ToEdn> ToEdn for HashSet<T> {
    fn to_edn(&self) -> Edn<'_> {
        let mut s = EdnSet::new();
        for item in self {
            s.insert(item.to_edn());
        }
        Edn::Set(s)
    }
}

impl<'de, T> FromEdn<'de> for BTreeSet<T>
where
    T: FromEdn<'de> + Ord,
{
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        let set = match edn {
            Edn::Set(s) => s,
            other => {
                return Err(EdnError::TypeMismatch {
                    expected: "set",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        set.iter().map(T::from_edn).collect()
    }
}

impl<T: ToEdn> ToEdn for BTreeSet<T> {
    fn to_edn(&self) -> Edn<'_> {
        let mut s = EdnSet::new();
        for item in self {
            s.insert(item.to_edn());
        }
        Edn::Set(s)
    }
}

impl<'de> FromEdn<'de> for Edn<'de> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        Ok(edn.clone())
    }
}

impl<'a> ToEdn for Edn<'a> {
    fn to_edn(&self) -> Edn<'_> {
        self.clone()
    }
}

impl<'de> FromEdn<'de> for Keyword<'de> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Keyword(k) => Ok(k.clone()),
            other => Err(EdnError::TypeMismatch {
                expected: "keyword",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for Keyword<'a> {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Keyword(self.clone())
    }
}

impl<'de> FromEdn<'de> for Symbol<'de> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Symbol(s) => Ok(s.clone()),
            other => Err(EdnError::TypeMismatch {
                expected: "symbol",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for Symbol<'a> {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Symbol(self.clone())
    }
}

impl<'de> FromEdn<'de> for EdnMap<'de> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Map(m) => Ok(m.clone()),
            other => Err(EdnError::TypeMismatch {
                expected: "map",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for EdnMap<'a> {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Map(self.clone())
    }
}

impl<'de> FromEdn<'de> for EdnSet<'de> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Set(s) => Ok(s.clone()),
            other => Err(EdnError::TypeMismatch {
                expected: "set",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for EdnSet<'a> {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Set(self.clone())
    }
}

/// `EdnSeq` deserializes from either a `Vector` or a `List` — same convention
/// as `Vec<T>`. Round-tripping a `List` source through `EdnSeq` and back via
/// `to_edn` produces a `Vector`; consumers that need to preserve list-vs-vector
/// distinction should match on `Edn` directly.
impl<'de> FromEdn<'de> for EdnSeq<'de> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::Vector(s) | Edn::List(s) => Ok(s.clone()),
            other => Err(EdnError::TypeMismatch {
                expected: "vector",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for EdnSeq<'a> {
    fn to_edn(&self) -> Edn<'_> {
        Edn::Vector(self.clone())
    }
}

#[cfg(feature = "bignum")]
impl<'de> FromEdn<'de> for BigInt {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::BigInt(n) => Ok(n.clone()),
            Edn::Int(n) => Ok(BigInt::from(*n)),
            other => Err(EdnError::TypeMismatch {
                expected: "bigint",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

#[cfg(feature = "bignum")]
impl ToEdn for BigInt {
    fn to_edn(&self) -> Edn<'_> {
        Edn::BigInt(self.clone())
    }
}

#[cfg(feature = "bignum")]
impl<'de> FromEdn<'de> for BigDecimal {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::BigDecimal(d) => Ok(d.as_inner().clone()),
            other => Err(EdnError::TypeMismatch {
                expected: "bigdecimal",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

#[cfg(feature = "bignum")]
impl ToEdn for BigDecimal {
    fn to_edn(&self) -> Edn<'_> {
        Edn::BigDecimal(EdnBigDecimal::new(self.clone()))
    }
}

#[cfg(feature = "bignum")]
impl<'de> FromEdn<'de> for EdnBigDecimal {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        match edn {
            Edn::BigDecimal(d) => Ok(d.clone()),
            other => Err(EdnError::TypeMismatch {
                expected: "bigdecimal",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

#[cfg(feature = "bignum")]
impl ToEdn for EdnBigDecimal {
    fn to_edn(&self) -> Edn<'_> {
        Edn::BigDecimal(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case::bool_true(true)]
    #[case::bool_false(false)]
    fn test_bool_roundtrip(#[case] v: bool) {
        let e = v.to_edn();
        assert_eq!(bool::from_edn(&e).unwrap(), v);
    }

    #[rstest]
    #[case::zero(0i32)]
    #[case::pos(1234i32)]
    #[case::neg(-77i32)]
    fn test_i32_roundtrip(#[case] v: i32) {
        let e = v.to_edn();
        assert_eq!(i32::from_edn(&e).unwrap(), v);
    }

    #[rstest]
    #[case::pos(2.71828f64)]
    #[case::neg(-1.5f64)]
    #[case::zero(0.0f64)]
    fn test_f64_roundtrip(#[case] v: f64) {
        let e = v.to_edn();
        assert_eq!(f64::from_edn(&e).unwrap(), v);
    }

    #[rstest]
    #[case::ascii('z')]
    #[case::unicode('λ')]
    fn test_char_roundtrip(#[case] v: char) {
        let e = v.to_edn();
        assert_eq!(char::from_edn(&e).unwrap(), v);
    }

    #[test]
    fn test_string_roundtrip() {
        let s = "carbon".to_string();
        let e = s.to_edn();
        assert_eq!(String::from_edn(&e).unwrap(), s);
    }

    #[test]
    fn test_cow_str_borrows_from_source() {
        let edn = Edn::Str(Cow::Borrowed("oxygen"));
        let cow = <Cow<'_, str>>::from_edn(&edn).unwrap();
        assert!(matches!(cow, Cow::Borrowed("oxygen")));
    }

    #[test]
    fn test_bool_from_edn_error() {
        let edn = Edn::Int(1);
        let err = bool::from_edn(&edn).unwrap_err();
        assert!(matches!(
            err,
            EdnError::TypeMismatch {
                expected: "bool",
                got: "int",
                ..
            }
        ));
    }

    #[test]
    fn test_i32_from_edn_error_out_of_range() {
        let edn = Edn::Int(i64::MAX);
        let err = i32::from_edn(&edn).unwrap_err();
        assert!(matches!(err, EdnError::OutOfRange { target: "i32", .. }));
    }

    #[test]
    fn test_vec_roundtrip() {
        let v = vec![1i32, 2, 3];
        let e = v.to_edn();
        assert_eq!(Vec::<i32>::from_edn(&e).unwrap(), v);
    }

    #[test]
    fn test_vec_accepts_list_and_vector() {
        let vector = Edn::Vector(EdnSeq::from(vec![Edn::Int(1), Edn::Int(2)]));
        let list = Edn::List(EdnSeq::from(vec![Edn::Int(1), Edn::Int(2)]));
        assert_eq!(Vec::<i32>::from_edn(&vector).unwrap(), vec![1, 2]);
        assert_eq!(Vec::<i32>::from_edn(&list).unwrap(), vec![1, 2]);
    }

    #[rstest]
    #[case(Some(7))]
    #[case(None)]
    fn test_option_roundtrip(#[case] v: Option<i32>) {
        let e = v.to_edn();
        assert_eq!(<Option<i32>>::from_edn(&e).unwrap(), v);
    }

    #[test]
    fn test_array_fixed_length() {
        let arr = [1i32, 2, 3];
        let e = arr.to_edn();
        assert_eq!(<[i32; 3]>::from_edn(&e).unwrap(), arr);
    }

    #[test]
    fn test_array_wrong_length_errors() {
        let edn = Edn::Vector(EdnSeq::from(vec![Edn::Int(1), Edn::Int(2)]));
        let err = <[i32; 3]>::from_edn(&edn).unwrap_err();
        assert!(matches!(err, EdnError::OutOfRange { .. }));
    }

    #[test]
    fn test_tuple3_roundtrip() {
        let t = (1.0f64, 2.0f64, 3.0f64);
        let e = t.to_edn();
        assert_eq!(<(f64, f64, f64)>::from_edn(&e).unwrap(), t);
    }

    #[test]
    fn test_btreemap_roundtrip() {
        let mut m = BTreeMap::new();
        m.insert("a".to_string(), 1i32);
        m.insert("b".to_string(), 2i32);
        let e = m.to_edn();
        assert_eq!(<BTreeMap<String, i32>>::from_edn(&e).unwrap(), m);
    }

    #[test]
    fn test_btreeset_roundtrip() {
        let mut s = BTreeSet::new();
        s.insert(1i32);
        s.insert(2i32);
        let e = s.to_edn();
        assert_eq!(<BTreeSet<i32>>::from_edn(&e).unwrap(), s);
    }

    #[test]
    fn test_edn_passthrough() {
        let edn = Edn::Vector(EdnSeq::from(vec![
            Edn::Keyword(Keyword::new("atoms")),
            Edn::Int(2),
        ]));
        let parsed = <Edn<'_>>::from_edn(&edn).unwrap();
        assert_eq!(parsed, edn);
    }

    #[test]
    fn test_keyword_roundtrip() {
        let k = Keyword::new("atoms");
        let e = k.to_edn();
        assert_eq!(Keyword::from_edn(&e).unwrap(), k);
    }

    #[test]
    fn test_map_helper_required_and_optional() {
        let mut m = EdnMap::new();
        m.insert(Edn::keyword("name"), Edn::Str(Cow::Borrowed("water")));
        m.insert(Edn::keyword("count"), Edn::Int(2));
        let mut h = EdnMapHelper::new(&m);
        let name: String = h.required("name").unwrap();
        let count: i32 = h.required("count").unwrap();
        let charge: Option<i32> = h.optional("charge").unwrap();
        assert_eq!(name, "water");
        assert_eq!(count, 2);
        assert_eq!(charge, None);
        h.finalize().unwrap();
    }

    #[test]
    fn test_map_helper_missing_required() {
        let m = EdnMap::new();
        let mut h = EdnMapHelper::new(&m);
        let err = h.required::<String>("name").unwrap_err();
        assert!(matches!(err, EdnError::MissingField { .. }));
    }

    #[test]
    fn test_map_helper_finalize_unknown_key() {
        let mut m = EdnMap::new();
        m.insert(Edn::keyword("name"), Edn::Str(Cow::Borrowed("water")));
        m.insert(Edn::keyword("extra"), Edn::Int(0));
        let mut h = EdnMapHelper::new(&m);
        let _name: String = h.required("name").unwrap();
        let err = h.finalize().unwrap_err();
        match err {
            EdnError::UnknownField { key, .. } => assert_eq!(key, "extra"),
            other => panic!("expected UnknownField, got {other:?}"),
        }
    }

    #[test]
    fn test_from_edn_str_default() {
        let v: Vec<i32> = Vec::<i32>::from_edn_str("[1 2 3]").unwrap();
        assert_eq!(v, vec![1, 2, 3]);
    }
}
