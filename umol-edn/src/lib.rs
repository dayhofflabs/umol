//! EDN (Extensible Data Notation) parser, formatter, and serde integration.

pub mod config;
mod display;
pub mod edn;
pub mod error;
pub mod formatter;
pub mod parser;
pub mod reader;
pub mod tags;

#[cfg(feature = "serde")]
pub mod de;
#[cfg(feature = "serde")]
pub mod ser;
#[cfg(feature = "serde")]
pub(crate) mod streaming;

pub use config::{AutoResolve, Dialect, DuplicateKeyPolicy, ParseConfig, TagFn, TagReaders};
#[cfg(feature = "serde")]
pub use de::{from_str, from_str_with, from_value, EdnDeserializer, StreamDeserializer};
pub use edn::{Edn, EdnMap, EdnSet, Keyword, Symbol};
pub use error::EdnError;
pub use formatter::EdnFormatter;
#[cfg(feature = "serde")]
pub use formatter::{to_string_pretty, to_string_pretty_with};
pub use reader::{read_all, read_all_with, read_string, read_string_with, Reader};
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
