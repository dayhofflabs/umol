//! EDN (Extensible Data Notation) parser, formatter, and serde integration.

pub mod config;
mod display;
pub mod edn;
pub mod error;
pub mod formatter;
pub mod parser;
pub mod reader;

#[cfg(feature = "serde")]
pub mod de;
#[cfg(feature = "serde")]
pub mod ser;
#[cfg(feature = "serde")]
pub(crate) mod streaming;

pub use config::{Dialect, DuplicateKeyPolicy, ParseConfig};
#[cfg(feature = "serde")]
pub use de::{from_str, from_value, EdnDeserializer, StreamDeserializer};
pub use edn::{Edn, EdnMap, EdnSet, Keyword, Symbol};
pub use error::EdnError;
pub use formatter::EdnFormatter;
#[cfg(feature = "serde")]
pub use formatter::{to_string_pretty, to_string_pretty_with};
pub use reader::{read_all, read_all_with, read_string, read_string_with, Reader};
#[cfg(feature = "serde")]
pub use ser::to_string;
