//! SMILES format parser and writer.

pub mod linter;
pub mod parser;
// pub mod parser_old;
// removed legacy iterators/state
#[cfg(test)]
pub mod test_support;

pub use parser::parse_smiles;
pub use parser::M6Error as ParseError;
