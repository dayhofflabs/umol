//! SMILES format parser, linter, and writer.

pub mod config;
pub mod error;
// pub mod linter;
pub mod parser;

pub use config::SmilesIoConfig;
pub use error::ParseError;
// pub use linter::{lint_smiles, lint_smiles_with};
pub use parser::{parse_smiles, parse_smiles_with};
