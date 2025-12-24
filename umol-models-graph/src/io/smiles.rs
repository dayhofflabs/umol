//! SMILES format parser, linter, and writer.

pub mod config;
pub mod linter;
pub mod parser;

// TODO: move linter to crate top level
pub use linter::{lint_smiles, lint_smiles_with};
pub use parser::{parse_smiles, parse_smiles_to_sir, parse_smiles_with, ParseError};
