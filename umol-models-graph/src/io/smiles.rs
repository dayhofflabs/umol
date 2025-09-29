//! SMILES format parser, linter, and writer.

pub mod linter;
pub mod parser;

pub use parser::{parse_smiles, ParseError};
