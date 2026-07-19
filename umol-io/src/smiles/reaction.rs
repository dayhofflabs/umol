use std::str::FromStr;

use super::config::SmilesIoConfig;
use super::error::ParseError;
use super::parser::parse_reaction;
use crate::table_ir::{Reaction, SourceFormat};

/// Parsed semantic value of a reaction SMILES representation.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionSmiles {
    table_ir: Reaction,
}

impl ReactionSmiles {
    /// Parse reaction SMILES text with the OpenSMILES configuration.
    pub fn parse(input: &str) -> Result<Self, ParseError> {
        Self::parse_bytes(input.as_bytes())
    }

    /// Parse reaction SMILES bytes with the OpenSMILES configuration.
    pub fn parse_bytes(input: &[u8]) -> Result<Self, ParseError> {
        Self::parse_bytes_with(input, &SmilesIoConfig::opensmiles())
    }

    /// Parse reaction SMILES text with an explicit IO configuration.
    pub fn parse_with(input: &str, config: &SmilesIoConfig) -> Result<Self, ParseError> {
        Self::parse_bytes_with(input.as_bytes(), config)
    }

    /// Parse reaction SMILES bytes with an explicit IO configuration.
    pub fn parse_bytes_with(input: &[u8], config: &SmilesIoConfig) -> Result<Self, ParseError> {
        parse_reaction(input, config).map(Self::from_parsed)
    }

    /// Borrow the neutral TableIR boundary value.
    pub fn as_table_ir(&self) -> &Reaction {
        &self.table_ir
    }

    /// Consume the reaction SMILES value and return its neutral TableIR boundary value.
    pub fn into_table_ir(self) -> Reaction {
        self.table_ir
    }

    fn from_parsed(mut table_ir: Reaction) -> Self {
        table_ir.source_format = SourceFormat::SMILES;
        Self { table_ir }
    }
}

impl FromStr for ReactionSmiles {
    type Err = ParseError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        Self::parse(input)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element;

    use super::*;
    use crate::table_ir::{Atom, Molecule, Span};

    #[rstest]
    #[case::empty(
        ">>",
        Reaction {
            source_format: SourceFormat::SMILES,
            ..Reaction::empty()
        }
    )]
    #[case::simple(
        "C>>C",
        Reaction {
            reactants: Molecule {
                atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))],
                source_format: SourceFormat::SMILES,
                ..Molecule::empty()
            },
            products: Molecule {
                atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))],
                source_format: SourceFormat::SMILES,
                ..Molecule::empty()
            },
            source_format: SourceFormat::SMILES,
            ..Reaction::empty()
        }
    )]
    fn test_reaction_smiles_parse(#[case] input: &str, #[case] expected: Reaction) {
        let reaction = ReactionSmiles::parse(input).unwrap();
        assert_eq!(reaction.as_table_ir(), &expected);
    }

    #[rstest]
    #[case::molecule_only("C", ParseError::MissingReactionArrow { pos: 1 })]
    #[case::leading_whitespace(" C>>C", ParseError::LeadingWhitespace)]
    fn test_reaction_smiles_parse_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(ReactionSmiles::parse(input), Err(expected));
    }
}
