//! CTFile format parsers (CTAB, MOL, SDF).
//!
//! This module provides parsing and writing for MDL's Connection Table (CTFile) formats.

pub mod config;
pub mod error;
pub mod parser;

pub use config::CtabParseFlags;
pub use error::ParseError;
pub use parser::{
    parse_extended_mol, parse_extended_mol_bytes, parse_extended_mol_bytes_with,
    parse_extended_mol_with, parse_extended_sdf, parse_extended_sdf_bytes,
    parse_extended_sdf_bytes_with, parse_extended_sdf_with, parse_mol_bytes_to_ir, parse_mol_to_ir,
    parse_sdf, parse_sdf_bytes, parse_sdf_bytes_with, parse_sdf_with,
};
