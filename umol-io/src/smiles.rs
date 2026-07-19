//! SMILES format parser, linter, and writer.

pub mod config;
pub mod error;
mod molecule;
mod parser;
mod reaction;

pub use config::SmilesIoConfig;
pub use error::ParseError;
pub use molecule::Smiles;
pub use parser::{
    parse_extended_reaction_smiles, parse_extended_reaction_smiles_bytes,
    parse_extended_reaction_smiles_bytes_with, parse_extended_reaction_smiles_with,
    parse_extended_smiles, parse_extended_smiles_bytes, parse_extended_smiles_bytes_with,
    parse_extended_smiles_with,
};
pub use reaction::ReactionSmiles;
