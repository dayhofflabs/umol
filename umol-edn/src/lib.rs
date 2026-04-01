//! EDN (Extensible Data Notation) parser, formatter, and serde integration.

pub mod edn;
pub mod error;
pub mod parser;
mod display;
pub mod reader;

pub use edn::{Edn, Keyword, Symbol};
pub use error::{EdnError, Span};
pub use parser::{DuplicateKeyPolicy, ParseConfig};
pub use reader::{read_all, read_all_with, read_string, read_string_with, Reader};
