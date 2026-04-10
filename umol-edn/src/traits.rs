//! Conversion traits between Rust types and [`Edn`] values.
//!
//! [`FromEdn`] builds a Rust value from a parsed [`Edn`] tree, and [`ToEdn`]
//! turns a Rust value into an [`Edn`] tree. Together they are the primary
//! deserialization and serialization API for umol-edn and support the full
//! EDN data model — keywords, symbols, sets, tagged literals, and
//! arbitrary-precision numbers.
//!
//! ## The two traits
//!
//! - [`FromEdn`] takes a borrowed `&Edn<'de>` and produces a value. It
//!   carries a `'de` lifetime parameter so implementations can borrow
//!   string and key references from the source buffer when zero-copy is
//!   wanted.
//! - [`ToEdn`] takes `&self` and produces an owned `Edn<'static>`.
//!
//! ## Parsing from a string
//!
//! [`FromEdn::from_edn_str`] has a default implementation that calls
//! [`read_string`] to build a tree and then dispatches to `from_edn`. This
//! is correct for every type. Hot types may override `from_edn_str` to fuse
//! parsing and deserialization in a single pass; the override must produce
//! the same value and the same errors as the default for every input.
//!
//! ## Lifetimes
//!
//! `'de` is the lifetime of the source the implementation is allowed to
//! borrow from. A type that borrows string slices declares
//! `MyType<'a>: FromEdn<'a>`. A type that always owns its strings is
//! `MyType: for<'de> FromEdn<'de>`.

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::hash::Hash;

#[cfg(feature = "bignum")]
use bigdecimal::BigDecimal;
#[cfg(feature = "bignum")]
use num_bigint::BigInt;

use crate::collections::{EdnMap, EdnSeq, EdnSet};
#[cfg(feature = "bignum")]
use crate::edn::EdnBigDecimal;
use crate::edn::{Edn, EdnKeyword, EdnSymbol};
use crate::error::{DeError, EdnError};
use crate::reader::read_string;

// Variant discriminator strings come from `Edn::kind()`.

/// Build `Self` from an EDN value.
///
/// # Implementing
///
/// At minimum, implementations must provide [`from_edn`]. The default
/// [`from_edn_str`] parses EDN source to a tree and dispatches through
/// `from_edn`, which is correct for all types.
///
/// Performance-critical types may override [`from_edn_str`] to fuse parsing
/// and deserialization in a single pass. The override must be observationally
/// equivalent to the default — same value, same errors — for every input.
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
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError>;

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
        Ok(Self::from_edn(&tree)?)
    }
}

/// Convert `&self` into an owned EDN value.
///
/// # Ownership
///
/// `to_edn` returns `Edn<'static>` — a fully owned tree with no
/// borrows. This allows callers to delegate through temporaries
/// (e.g. `self.to_ast().to_edn()`) without lifetime issues.
///
/// # Infallibility
///
/// `to_edn` is infallible. Types whose serialized form depends on
/// validation should perform that validation at construction time so that
/// every reachable `&self` is serializable. Types that need fallible
/// serialization should add their own `try_to_edn` method until a concrete
/// need establishes a trait-level contract for it.
pub trait ToEdn {
    /// Produce a fully owned `Edn` value representing `self`.
    fn to_edn(&self) -> Edn<'static>;
}

impl<'de> FromEdn<'de> for bool {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Bool(b) => Ok(*b),
            other => Err(DeError::TypeMismatch {
                expected: "bool",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for bool {
    fn to_edn(&self) -> Edn<'static> {
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
                fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                    match edn {
                        Edn::Int(n) => <$t>::try_from(*n).map_err(|_| {
                            DeError::OutOfRange {
                                value: n.to_string(),
                                target: stringify!($t),
                                path: Vec::new(),
                            }
                        }),
                        other => Err(DeError::TypeMismatch {
                            expected: "int",
                            got: other.kind(),
                            path: Vec::new(),
                        }),
                    }
                }
            }

            impl ToEdn for $t {
                fn to_edn(&self) -> Edn<'static> {
                    Edn::Int(i64::from(*self))
                }
            }
        )*
    };
}

impl_int!(i8, i16, i32, i64, u8, u16, u32);

impl<'de> FromEdn<'de> for f64 {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Float(f) => Ok(*f),
            other => Err(DeError::TypeMismatch {
                expected: "float",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for f64 {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Float(*self)
    }
}

impl<'de> FromEdn<'de> for f32 {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Float(f) => Ok(*f as f32),
            other => Err(DeError::TypeMismatch {
                expected: "float",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for f32 {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Float(f64::from(*self))
    }
}

impl<'de> FromEdn<'de> for char {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Char(c) => Ok(*c),
            other => Err(DeError::TypeMismatch {
                expected: "char",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for char {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Char(*self)
    }
}

impl<'de> FromEdn<'de> for String {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => Ok(s.to_string()),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for String {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.clone()))
    }
}

impl ToEdn for str {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_owned()))
    }
}

/// `Cow<'de, str>` preserves zero-copy borrowing when the source string
/// itself was borrowed from the parsed buffer.
impl<'de> FromEdn<'de> for Cow<'de, str> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => Ok(s.clone()),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for Cow<'a, str> {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.as_ref().to_owned()))
    }
}

impl<'de, T: FromEdn<'de>> FromEdn<'de> for Option<T> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Nil => Ok(None),
            other => T::from_edn(other).map(Some),
        }
    }
}

impl<T: ToEdn> ToEdn for Option<T> {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Some(t) => t.to_edn(),
            None => Edn::Nil,
        }
    }
}

impl<'de, T: FromEdn<'de>> FromEdn<'de> for Vec<T> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let seq = match edn {
            Edn::Vector(s) | Edn::List(s) => s,
            other => {
                return Err(DeError::TypeMismatch {
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
    fn to_edn(&self) -> Edn<'static> {
        Edn::Vector(EdnSeq::from(
            self.iter().map(ToEdn::to_edn).collect::<Vec<_>>(),
        ))
    }
}

impl<'de, T: FromEdn<'de>, const N: usize> FromEdn<'de> for [T; N] {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let seq = match edn {
            Edn::Vector(s) | Edn::List(s) => s,
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "vector",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        if seq.len() != N {
            return Err(DeError::OutOfRange {
                value: seq.len().to_string(),
                target: "fixed-length array",
                path: Vec::new(),
            });
        }
        let parsed: Vec<T> = seq.iter().map(T::from_edn).collect::<Result<_, _>>()?;
        parsed
            .try_into()
            .map_err(|_| DeError::Custom("array length conversion failed".into()))
    }
}

impl<T: ToEdn, const N: usize> ToEdn for [T; N] {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Vector(EdnSeq::from(
            self.iter().map(ToEdn::to_edn).collect::<Vec<_>>(),
        ))
    }
}

/// Tuple impls up to arity 8
macro_rules! impl_tuple {
    ($($n:tt $t:ident),+) => {
        impl<'de, $($t: FromEdn<'de>),+> FromEdn<'de> for ($($t,)+) {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                let seq = match edn {
                    Edn::Vector(s) | Edn::List(s) => s,
                    other => {
                        return Err(DeError::TypeMismatch {
                            expected: "vector",
                            got: other.kind(),
                            path: Vec::new(),
                        });
                    }
                };
                const ARITY: usize = [$($n),+].len();
                if seq.len() != ARITY {
                    return Err(DeError::OutOfRange {
                        value: seq.len().to_string(),
                        target: "tuple arity",
                        path: Vec::new(),
                    });
                }
                Ok(($($t::from_edn(&seq[$n])?,)+))
            }
        }

        impl<$($t: ToEdn),+> ToEdn for ($($t,)+) {
            fn to_edn(&self) -> Edn<'static> {
                Edn::Vector(EdnSeq::from(vec![$(self.$n.to_edn()),+]))
            }
        }
    };
}

impl_tuple!(0 A, 1 B);
impl_tuple!(0 A, 1 B, 2 C);
impl_tuple!(0 A, 1 B, 2 C, 3 D);
impl_tuple!(0 A, 1 B, 2 C, 3 D, 4 E, 5 F);
impl_tuple!(0 A, 1 B, 2 C, 3 D, 4 E, 5 F, 6 G);
impl_tuple!(0 A, 1 B, 2 C, 3 D, 4 E, 5 F, 6 G, 7 H);

impl<'de, K, V> FromEdn<'de> for HashMap<K, V>
where
    K: FromEdn<'de> + Eq + Hash,
    V: FromEdn<'de>,
{
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let map = match edn {
            Edn::Map(m) => m,
            other => {
                return Err(DeError::TypeMismatch {
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
    fn to_edn(&self) -> Edn<'static> {
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
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let map = match edn {
            Edn::Map(m) => m,
            other => {
                return Err(DeError::TypeMismatch {
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
    fn to_edn(&self) -> Edn<'static> {
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
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let set = match edn {
            Edn::Set(s) => s,
            other => {
                return Err(DeError::TypeMismatch {
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
    fn to_edn(&self) -> Edn<'static> {
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
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let set = match edn {
            Edn::Set(s) => s,
            other => {
                return Err(DeError::TypeMismatch {
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
    fn to_edn(&self) -> Edn<'static> {
        let mut s = EdnSet::new();
        for item in self {
            s.insert(item.to_edn());
        }
        Edn::Set(s)
    }
}

#[cfg(feature = "indexmap")]
impl<'de, K, V> FromEdn<'de> for indexmap::IndexMap<K, V>
where
    K: FromEdn<'de> + Eq + Hash,
    V: FromEdn<'de>,
{
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let map = match edn {
            Edn::Map(m) => m,
            other => {
                return Err(DeError::TypeMismatch {
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

#[cfg(feature = "indexmap")]
impl<K, V> ToEdn for indexmap::IndexMap<K, V>
where
    K: ToEdn,
    V: ToEdn,
{
    fn to_edn(&self) -> Edn<'static> {
        let mut m = EdnMap::with_capacity(self.len());
        for (k, v) in self {
            m.insert(k.to_edn(), v.to_edn());
        }
        Edn::Map(m)
    }
}

impl<'de> FromEdn<'de> for Edn<'de> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(edn.clone())
    }
}

impl<'a> ToEdn for Edn<'a> {
    fn to_edn(&self) -> Edn<'static> {
        self.clone().into_owned()
    }
}

impl<'de> FromEdn<'de> for EdnKeyword<'static> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Keyword(k) => Ok(k.clone().into_owned()),
            other => Err(DeError::TypeMismatch {
                expected: "keyword",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for EdnKeyword<'a> {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Keyword(self.clone().into_owned())
    }
}

impl<'de> FromEdn<'de> for EdnSymbol<'static> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Symbol(s) => Ok(s.clone().into_owned()),
            other => Err(DeError::TypeMismatch {
                expected: "symbol",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for EdnSymbol<'a> {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Symbol(self.clone().into_owned())
    }
}

impl<'de> FromEdn<'de> for EdnMap<'de> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Map(m) => Ok(m.clone()),
            other => Err(DeError::TypeMismatch {
                expected: "map",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for EdnMap<'a> {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Map(self.clone().into_owned())
    }
}

impl<'de> FromEdn<'de> for EdnSet<'de> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Set(s) => Ok(s.clone()),
            other => Err(DeError::TypeMismatch {
                expected: "set",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for EdnSet<'a> {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Set(self.clone().into_owned())
    }
}

/// `EdnSeq` deserializes from either a `Vector` or a `List` — same convention
/// as `Vec<T>`. Round-tripping a `List` source through `EdnSeq` and back via
/// `to_edn` produces a `Vector`; consumers that need to preserve list-vs-vector
/// distinction should match on `Edn` directly.
impl<'de> FromEdn<'de> for EdnSeq<'de> {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Vector(s) | Edn::List(s) => Ok(s.clone()),
            other => Err(DeError::TypeMismatch {
                expected: "vector",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl<'a> ToEdn for EdnSeq<'a> {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Vector(self.clone().into_owned())
    }
}

#[cfg(feature = "bignum")]
impl<'de> FromEdn<'de> for BigInt {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::BigInt(n) => Ok(n.clone()),
            Edn::Int(n) => Ok(BigInt::from(*n)),
            other => Err(DeError::TypeMismatch {
                expected: "bigint",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

#[cfg(feature = "bignum")]
impl ToEdn for BigInt {
    fn to_edn(&self) -> Edn<'static> {
        Edn::BigInt(self.clone())
    }
}

#[cfg(feature = "bignum")]
impl<'de> FromEdn<'de> for BigDecimal {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::BigDecimal(d) => Ok(d.as_inner().clone()),
            other => Err(DeError::TypeMismatch {
                expected: "bigdecimal",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

#[cfg(feature = "bignum")]
impl ToEdn for BigDecimal {
    fn to_edn(&self) -> Edn<'static> {
        Edn::BigDecimal(EdnBigDecimal::new(self.clone()))
    }
}

#[cfg(feature = "bignum")]
impl<'de> FromEdn<'de> for EdnBigDecimal {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::BigDecimal(d) => Ok(d.clone()),
            other => Err(DeError::TypeMismatch {
                expected: "bigdecimal",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

#[cfg(feature = "bignum")]
impl ToEdn for EdnBigDecimal {
    fn to_edn(&self) -> Edn<'static> {
        Edn::BigDecimal(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

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
            DeError::TypeMismatch {
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
        assert!(matches!(err, DeError::OutOfRange { target: "i32", .. }));
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
        assert!(matches!(err, DeError::OutOfRange { .. }));
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
            Edn::Keyword(EdnKeyword::new("atoms")),
            Edn::Int(2),
        ]));
        let parsed = <Edn<'_>>::from_edn(&edn).unwrap();
        assert_eq!(parsed, edn);
    }

    #[test]
    fn test_keyword_roundtrip() {
        let k = EdnKeyword::new("atoms");
        let e = k.to_edn();
        assert_eq!(EdnKeyword::from_edn(&e).unwrap(), k);
    }

    #[test]
    fn test_from_edn_str_default() {
        let v: Vec<i32> = Vec::<i32>::from_edn_str("[1 2 3]").unwrap();
        assert_eq!(v, vec![1, 2, 3]);
    }
}
