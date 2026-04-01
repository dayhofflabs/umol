//! EDN (Extensible Data Notation) parser, formatter, and serde integration.

pub mod config;
mod display;
pub mod edn;
pub mod error;
pub mod parser;
pub mod reader;

#[cfg(feature = "serde")]
pub mod de;
#[cfg(feature = "serde")]
pub mod ser;

pub use config::{Dialect, DuplicateKeyPolicy, ParseConfig};
#[cfg(feature = "serde")]
pub use de::{from_str, EdnDeserializer, StreamDeserializer};
pub use edn::{Edn, Keyword, Symbol};
pub use error::EdnError;
pub use reader::{read_all, read_all_with, read_string, read_string_with, Reader};
#[cfg(feature = "serde")]
pub use ser::to_string;
