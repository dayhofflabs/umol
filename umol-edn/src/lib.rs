//! EDN (Extensible Data Notation) parser, formatter, and serde integration.

pub mod collections;
pub mod config;
mod display;
pub mod edn;
pub mod error;
pub mod formatter;
pub mod native;
pub mod parser;
pub mod reader;
pub mod tags;

#[cfg(feature = "serde")]
pub mod de;
#[cfg(feature = "serde")]
pub mod keyword_serde;
#[cfg(feature = "serde")]
pub mod ser;
#[cfg(feature = "serde")]
pub(crate) mod streaming;

#[cfg(feature = "bignum")]
pub use bigdecimal::BigDecimal;
pub use collections::{EdnKeyRef, EdnMap, EdnSeq, EdnSet};
pub use config::{DuplicateKeyPolicy, ParseConfig, TagFn, TagReaders};
#[cfg(feature = "serde")]
pub use de::{from_str, from_str_with, from_value, EdnDeserializer, StreamDeserializer};
#[cfg(feature = "bignum")]
pub use edn::EdnBigDecimal;
pub use edn::{Edn, Keyword, Symbol};
pub use error::EdnError;
pub use formatter::EdnFormatter;
pub use native::{FromEdn, ToEdn};
#[cfg(feature = "serde")]
pub use formatter::{to_string_pretty, to_string_with};
#[cfg(feature = "bignum")]
pub use num_bigint::BigInt;
pub use reader::{read_all, read_all_with, read_string, read_string_with, Reader};
#[cfg(feature = "serde")]
pub use keyword_serde::EdnKeyword;
#[cfg(feature = "serde")]
pub use ser::to_string;

#[cfg(not(feature = "macros"))]
/// Construct an `Edn<'static>` value from an EDN string literal.
#[macro_export]
macro_rules! edn {
    ($s:literal) => {
        $crate::read_string($s).expect("invalid EDN in edn! macro")
    };
}

#[cfg(feature = "macros")]
pub use umol_edn_macros::edn;
