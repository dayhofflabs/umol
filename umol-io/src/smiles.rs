//! SMILES format parser, linter, and writer.

use std::str::FromStr;

use crate::table_ir::{Molecule, SourceFormat};

pub mod config;
pub mod error;
pub mod parser;

pub use config::SmilesIoConfig;
pub use error::ParseError;
pub use parser::{
    parse_extended_smiles, parse_extended_smiles_bytes, parse_extended_smiles_bytes_with,
    parse_extended_smiles_with,
};

/// Parsed semantic value of a molecular SMILES representation.
///
/// The contained TableIR is private so it cannot drift away from the SMILES
/// representation established by parsing. This value does not preserve the
/// spelling of the source text.
#[derive(Clone, Debug, PartialEq)]
pub struct Smiles {
    table_ir: Molecule,
}

impl Smiles {
    /// Parse SMILES text with the OpenSMILES configuration.
    pub fn parse(input: &str) -> Result<Self, ParseError> {
        Self::parse_bytes(input.as_bytes())
    }

    /// Parse SMILES bytes with the OpenSMILES configuration.
    pub fn parse_bytes(input: &[u8]) -> Result<Self, ParseError> {
        Self::parse_bytes_with(input, &SmilesIoConfig::opensmiles())
    }

    /// Parse SMILES text with an explicit IO configuration.
    pub fn parse_with(input: &str, config: &SmilesIoConfig) -> Result<Self, ParseError> {
        Self::parse_bytes_with(input.as_bytes(), config)
    }

    /// Parse SMILES bytes with an explicit IO configuration.
    pub fn parse_bytes_with(input: &[u8], config: &SmilesIoConfig) -> Result<Self, ParseError> {
        parser::parse_smiles_bytes_to_table_ir_with(input, config).map(Self::from_parsed)
    }

    /// Borrow the neutral TableIR boundary value.
    pub fn as_table_ir(&self) -> &Molecule {
        &self.table_ir
    }

    /// Consume the SMILES value and return its neutral TableIR boundary value.
    pub fn into_table_ir(self) -> Molecule {
        self.table_ir
    }

    fn from_parsed(mut table_ir: Molecule) -> Self {
        table_ir.source_format = SourceFormat::SMILES;
        Self { table_ir }
    }
}

impl FromStr for Smiles {
    type Err = ParseError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        Self::parse(input)
    }
}

#[cfg(test)]
mod tests {
    use indexmap::IndexMap;
    use rstest::rstest;
    use umol_chem::element::Element;

    use super::*;
    use crate::table_ir::{Atom, Span};

    #[rstest]
    #[case::carbon(
        "C",
        Molecule {
            atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))],
            bonds: Vec::new(),
            positions: None,
            multicenter_bonds: Vec::new(),
            configuration_scope: None,
            chirality_frame: None,
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::SMILES,
        }
    )]
    #[case::wildcard(
        "*",
        Molecule {
            atoms: vec![Atom::wildcard_with_span(Span::bytes(0, 1))],
            bonds: Vec::new(),
            positions: None,
            multicenter_bonds: Vec::new(),
            configuration_scope: None,
            chirality_frame: None,
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::SMILES,
        }
    )]
    #[case::empty(
        "",
        Molecule {
            source_format: SourceFormat::SMILES,
            ..Molecule::empty()
        }
    )]
    fn test_smiles_parse(#[case] input: &str, #[case] expected: Molecule) {
        let smiles = Smiles::parse(input).unwrap();
        assert_eq!(smiles.as_table_ir(), &expected);
    }

    #[rstest]
    #[case::leading_whitespace(" C", ParseError::LeadingWhitespace)]
    #[case::invalid_element("Q", ParseError::InvalidElement { pos: 0 })]
    fn test_smiles_parse_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(Smiles::parse(input), Err(expected));
    }
}
