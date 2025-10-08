//! SMILES format parser, linter, and writer.

pub mod config;
pub mod checker;
pub mod diagnostics;
pub mod api;
// pub mod linter;
pub mod parser;

pub use parser::{parse_smiles, ParseError};
pub use checker::check_smiles;
