//! EDN (Extensible Data Notation) reader, writer, and Rust type bindings.
//!
//! # Overview
//!
//! The reader turns EDN text into an [`Edn`] tree. [`FromEdn`] / [`ToEdn`]
//! move between that tree and Rust types and cover the full EDN data model
//! (keywords, symbols, lists, sets, tagged literals, arbitrary-precision
//! numbers). [`FormatConfig`] drives the pretty-printer.
//!
//! ```text
//! EDN text  ──read_string──▶  Edn tree  ──FromEdn──▶  Rust value
//! Rust value ──ToEdn────────▶  Edn tree  ──to_string_with(&FormatConfig)──▶ EDN text
//! ```
//!
//! Hot types may override [`FromEdn::from_edn_str`] to fuse parsing and
//! construction in a single pass when the intermediate tree would be
//! wasted.
//!
//! # Typed wrappers for EDN-only constructs
//!
//! When a value holds an EDN-only construct (keyword, symbol, list, set,
//! tagged literal, bignum) but needs to pass through another format (JSON,
//! YAML, a serde-powered pipeline), use the wrapper type so the construct
//! degrades predictably: [`EdnKeyword`], [`EdnSymbol`], [`EdnList`],
//! [`EdnHashSet`], [`EdnTagged`], [`EdnBigInt`], [`EdnBigDecimal`]. The
//! dynamic [`Value`] is lossless and the preferred choice for
//! schema-agnostic data.
//!
//! # Optional features
//!
//! - `serde` — adds [`de::from_str`] and [`ser::to_string`] for types
//!   implementing `serde::Deserialize` / `serde::Serialize`.
//! - `bignum` — enables [`EdnBigInt`] and [`EdnBigDecimal`].
//! - `macros` — re-exports `#[derive(FromEdn, ToEdn)]` and the `edn!` macro.

#[cfg(feature = "bignum")]
pub mod bigdecimal;
#[cfg(feature = "bignum")]
pub mod bigint;
pub mod collections;
pub mod config;
mod display;
pub mod edn;
pub mod error;
pub mod formatter;
pub mod keyword;
pub mod list;
pub mod parser;
pub mod reader;
pub mod set;
pub mod streaming;
pub mod symbol;
pub mod tagged;
pub mod tags;
pub mod traits;
pub mod value;

#[cfg(feature = "serde")]
pub mod de;
#[cfg(feature = "serde")]
pub mod ser;
#[cfg(feature = "serde")]
pub(crate) mod serde_tokens;

#[cfg(feature = "bignum")]
pub use bigdecimal::EdnBigDecimal;
#[cfg(feature = "bignum")]
pub use bigint::EdnBigInt;
pub use collections::{EdnKeyRef, EdnMap, EdnMapHelper, EdnSeq, EdnSet};
pub use config::{DuplicateKeyPolicy, ParseConfig, TagFn, TagReaders};
#[cfg(feature = "serde")]
pub use de::{from_str, from_str_with, from_value, EdnDeserializer};
pub use edn::{Edn, Keyword, Symbol};
pub use error::EdnError;
pub use formatter::FormatConfig;
#[cfg(feature = "serde")]
pub use formatter::{to_string_pretty, to_string_with};
pub use keyword::EdnKeyword;
pub use list::EdnList;
pub use reader::{read_all, read_all_with, read_string, read_string_with, Reader};
#[cfg(feature = "serde")]
pub use ser::{to_string, to_value, EdnSerializer};
pub use set::EdnHashSet;
pub use symbol::EdnSymbol;
pub use tagged::EdnTagged;
pub use traits::{FromEdn, ToEdn};
#[cfg(feature = "macros")]
pub use umol_edn_macros::{edn, FromEdn, ToEdn};
pub use value::Value;
