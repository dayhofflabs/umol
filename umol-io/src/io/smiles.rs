//! SMILES format parser, linter, and writer.

pub mod config;
pub mod error;
pub mod parser;

pub use config::SmilesIoConfig;
pub use error::ParseError;
pub use parser::{
    parse_extended_smiles, parse_extended_smiles_bytes, parse_extended_smiles_bytes_with,
    parse_extended_smiles_with, parse_smiles_bytes_to_ast, parse_smiles_to_ast,
};
