//! EDN (Extensible Data Notation) parser, formatter, and serde integration.

pub mod config;
pub mod edn;
pub mod error;
pub mod parser;
mod display;
pub mod reader;

#[cfg(feature = "serde")]
pub mod de;
#[cfg(feature = "serde")]
pub mod ser;

pub use edn::{Edn, Keyword, Symbol};
pub use config::{Dialect, DuplicateKeyPolicy, ParseConfig};
pub use error::EdnError;
pub use reader::{read_all, read_all_with, read_string, read_string_with, Reader};

#[cfg(feature = "serde")]
pub use de::from_str;
#[cfg(feature = "serde")]
pub use ser::to_string;
#[cfg(feature = "serde")]
pub use de::EdnDeserializer;
