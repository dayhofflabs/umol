//! EDN (Extensible Data Notation) reader, writer, and Rust type bindings.
//!
//! # Data flow
//!
//! ```text
//! EDN text ── read_string ──▶ Edn tree ── FromEdn ──▶ Rust value
//! Rust value ── ToEdn ──────▶ Edn tree ── to_string_with ──▶ EDN text
//! ```
//!
//! [`FromEdn`] / [`ToEdn`] are the primary conversion traits and cover
//! the full EDN data model: keywords, symbols, lists, sets, tagged literals,
//! and arbitrary-precision numbers. [`FormatConfig`] drives the
//! pretty-printer.
//!
//! Hot types may override [`FromEdn::from_edn_str`] to fuse parsing and
//! construction in a single pass, bypassing the intermediate tree.
//!
//! # Choosing the right entry point
//!
//! ## Parsing (EDN text → Rust)
//!
//! | Function | Input | Output | Notes |
//! |---|---|---|---|
//! | [`read_string`] | `&str` | `Edn<'a>` | Zero-copy tree; borrows from input |
//! | [`read_string_with`] | `&str`, [`ParseConfig`] | `Edn<'a>` | Custom tag readers, dup-key policy |
//! | [`read_all`] | `&str` | `Vec<Edn<'a>>` | Multiple top-level forms |
//! | [`read_all_with`] | `&str`, [`ParseConfig`] | `Vec<Edn<'a>>` | Multiple forms + config |
//! | `T::from_edn(&edn)` | `&Edn` | `T` | Tree → typed value via [`FromEdn`] |
//! | `T::from_edn_str(s)` | `&str` | `T` | Fused parse+deserialize (overridable) |
//! | [`serde::from_str`] | `&str` | `T: Deserialize` | serde path (requires `serde` feature) |
//! | [`serde::from_value`] | `Edn` (owned) | `T: Deserialize` | serde path from tree |
//! | [`serde::from_value_ref`] | `&Edn` | `T: Deserialize` | serde path, borrows tree |
//!
//! ## Output (Rust → EDN text)
//!
//! | Function | Input | Output | Notes |
//! |---|---|---|---|
//! | `val.to_edn()` | `&T` | `Edn<'_>` | Tree via [`ToEdn`], borrows from value |
//! | `edn.to_string()` | `&Edn` | `String` | Compact one-line format |
//! | `edn.to_string_with(cfg)` | `&Edn`, [`FormatConfig`] | `String` | Pretty-printed |
//! | [`serde::to_string`] | `&T: Serialize` | `String` | serde path (requires `serde` feature) |
//! | [`serde::to_string_pretty`] | `&T: Serialize` | `String` | serde path, default pretty config |
//! | [`serde::to_string_with`] | `&T`, [`FormatConfig`] | `String` | serde path + custom format |
//! | [`serde::to_value`] | `&T: Serialize` | `EdnOwned` | serde → tree |
//!
//! # Type families
//!
//! **Core types** — always available, no feature flags:
//!
//! - [`Edn<'a>`](Edn) — parsed value tree. Borrows string data from the
//!   input buffer via `'a`. Call [`Edn::into_owned`] to get [`EdnOwned`]
//!   (`= Edn<'static>`) when the source outlives the tree.
//! - [`EdnKeyword<'a>`](EdnKeyword), [`EdnSymbol<'a>`](EdnSymbol) — EDN `:keyword`
//!   and `symbol` values.
//! - [`EdnMap<'a>`](EdnMap), [`EdnSet<'a>`](EdnSet),
//!   [`EdnSeq<'a>`](EdnSeq) — map, set, and ordered-sequence containers.
//! - [`EdnKeyRef`] — borrow-friendly lookup key for maps and sets without
//!   allocating an `Edn` node.
//! - [`FromEdn`], [`ToEdn`] — conversion traits.
//! - [`ParseConfig`], [`FormatConfig`] — reader and writer configuration.
//! - [`Reader`] — streaming pull parser.
//! - [`EdnMapHelper`] — keyword-keyed map reader with required/optional
//!   field tracking and unknown-key detection.
//! - [`EdnStreamDeserializer`] — byte-level streaming parser for fused
//!   `FromEdn::from_edn_str` overrides.
//!
//! **[`serde`] module** — wrapper types and serde-feature functions:
//!
//! Wrapper types ([`serde::EdnKeyword`], [`serde::EdnSymbol`],
//! [`serde::EdnList`], [`serde::EdnHashSet`], [`serde::EdnTagged`],
//! [`serde::DynEdn`]) carry EDN-only constructs through any serialization
//! format. Through the EDN serializer they preserve full fidelity; through
//! JSON or other formats they degrade to the closest equivalent.
//!
//! With the `serde` feature enabled, the module also re-exports
//! [`serde::from_str`], [`serde::to_string`], [`serde::EdnDeserializer`],
//! [`serde::EdnSerializer`], and related functions.
//!
//! # Lifetimes
//!
//! `Edn<'a>` borrows string data from the input buffer. This avoids
//! allocation for the common parse-inspect-discard pattern:
//!
//! ```rust
//! # use umol_edn::{read_string, Edn};
//! let input = r#"{:name "water" :atoms ["H" "H" "O"]}"#;
//! let edn = read_string(input).unwrap();
//! // `edn` borrows from `input` — no heap copies for string values.
//! ```
//!
//! When the tree must outlive the input, call [`Edn::into_owned`] to get
//! an [`EdnOwned`] (`Edn<'static>`).
//!
//! For [`FromEdn`] implementations: types that own all their data implement
//! `FromEdn<'de>` for any `'de`. Types that borrow strings tie `'de` to
//! their own lifetime. See the [`FromEdn`] trait docs for details.
//!
//! # Features
//!
//! | Feature | Adds |
//! |---|---|
//! | `serde` | [`serde::from_str`], [`serde::to_string`], [`serde::EdnDeserializer`], [`serde::EdnSerializer`] |
//! | `bignum` | `Edn::BigInt`, `Edn::BigDecimal`, [`serde::EdnBigInt`], [`serde::EdnBigDecimal`] |
//! | `chrono` | [`inst_to_edn`] — tag reader for `#inst` → `chrono::DateTime` |
//! | `uuid` | [`uuid_to_edn`] — tag reader for `#uuid` → `uuid::Uuid` |
//! | `macros` | `#[derive(FromEdn, ToEdn)]` and `edn!` literal macro |
//!
//! # ParseConfig / FormatConfig
//!
//! [`ParseConfig`] controls the reader: duplicate-key policy, tag reader
//! registry, unknown-tag behavior. [`FormatConfig`] controls output:
//! indentation, line width, map key sorting. Both have sensible defaults
//! via `Default::default()`.
//!
//! ```rust
//! # use umol_edn::{read_string_with, ParseConfig, DuplicateKeyPolicy};
//! let mut cfg = ParseConfig::default();
//! cfg.duplicate_keys = DuplicateKeyPolicy::LastWins;
//! let edn = read_string_with(r#"{:a 1 :a 2}"#, &cfg).unwrap();
//! ```

#[cfg(feature = "bignum")]
pub(crate) mod bigdecimal;
#[cfg(feature = "bignum")]
pub(crate) mod bigint;
pub(crate) mod collections;
pub(crate) mod config;
mod display;
pub(crate) mod dyn_edn;
pub(crate) mod edn;
pub(crate) mod error;
pub(crate) mod formatter;
pub(crate) mod list;
pub(crate) mod parser;
pub(crate) mod reader;
pub(crate) mod set;
pub(crate) mod streaming;
pub(crate) mod tagged;
pub(crate) mod tags;
pub(crate) mod traits;

#[cfg(feature = "serde")]
pub(crate) mod de;
#[cfg(feature = "serde")]
pub(crate) mod ser;

pub use collections::{EdnKeyRef, EdnMap, EdnMapHelper, EdnSeq, EdnSet};
pub use config::{DuplicateKeyPolicy, ParseConfig, TagFn, TagReaders};
pub use edn::{Edn, EdnKeyword, EdnSymbol};
pub use error::{DeError, EdnError, ParseError, SerError};
pub use formatter::FormatConfig;
pub use reader::{read_all, read_all_with, read_string, read_string_with, Reader};
pub use streaming::EdnStreamDeserializer;
#[cfg(feature = "chrono")]
pub use tags::inst_to_edn;
#[cfg(feature = "uuid")]
pub use tags::uuid_to_edn;
pub use traits::{FromEdn, ToEdn};
#[cfg(feature = "macros")]
pub use umol_edn_macros::{edn, FromEdn, ToEdn};

pub mod serde;
