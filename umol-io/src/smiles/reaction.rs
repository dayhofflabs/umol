use std::str::FromStr;

use super::config::SmilesIoConfig;
use super::error::{ParseError, ReactionSmilesRenderError};
use super::parser::parse_reaction;
use super::render::append_molecule;
use crate::table_ir::{Reaction, SourceFormat};

/// Semantic value of a reaction SMILES representation.
///
/// Owns a TableIR without changing or validating it. Rendering checks the properties it needs.
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
        let mut table_ir = parse_reaction(input, config)?;
        table_ir.source_format = SourceFormat::SMILES;
        Ok(Self::from_table_ir(table_ir))
    }

    /// Render with the OpenSMILES configuration.
    ///
    /// # Errors
    ///
    /// Fails for invalid references, unsupported fields, or unrepresentable stereo.
    pub fn render(&self) -> Result<String, ReactionSmilesRenderError> {
        self.render_with(&SmilesIoConfig::opensmiles())
    }

    /// Render reactants, agents, and products, in that order, separated by `>`.
    ///
    /// Atom.class supplies the written map labels. The derived atom_mapping index is not
    /// consumed. Traversal and marker assignment run once for each section.
    ///
    /// # Errors
    ///
    /// Reports molecular failures with section context. Unsupported reaction metadata fails
    /// rather than being omitted.
    ///
    /// # Semantic properties
    ///
    /// Preserves supported molecular semantics, section order, and every atom-class label,
    /// including repeated labels and labels on agents. Output is deterministic for the
    /// table and configuration; source spelling and ring labels need not be preserved.
    pub fn render_with(
        &self,
        config: &SmilesIoConfig,
    ) -> Result<String, ReactionSmilesRenderError> {
        let table = &self.table_ir;
        if !table.comments.is_empty() {
            return Err(ReactionSmilesRenderError::UnsupportedReaction { field: "comments" });
        }
        if !table.properties.is_empty() {
            return Err(ReactionSmilesRenderError::UnsupportedReaction {
                field: "properties",
            });
        }
        let mut output = String::with_capacity(
            table.reactants.atoms.len() + table.agents.atoms.len() + table.products.atoms.len() + 2,
        );
        append_molecule(&mut output, &table.reactants, config)
            .map_err(ReactionSmilesRenderError::Reactants)?;
        output.push('>');
        append_molecule(&mut output, &table.agents, config)
            .map_err(ReactionSmilesRenderError::Agents)?;
        output.push('>');
        append_molecule(&mut output, &table.products, config)
            .map_err(ReactionSmilesRenderError::Products)?;
        Ok(output)
    }

    /// Take ownership of a TableIR for reaction SMILES rendering.
    pub fn from_table_ir(table_ir: Reaction) -> Self {
        Self { table_ir }
    }

    /// Consume the reaction SMILES value and return its neutral TableIR boundary value.
    pub fn into_table_ir(self) -> Reaction {
        self.table_ir
    }

    /// Borrow the neutral TableIR boundary value.
    pub fn as_table_ir(&self) -> &Reaction {
        &self.table_ir
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

    use super::super::error::SmilesRenderError;
    use super::*;
    use crate::table_ir::{Atom, Bond, BondOrder, Molecule, Span};

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

    #[rstest]
    #[case::empty(">>")]
    #[case::reactants("C>>")]
    #[case::agents(">O>")]
    #[case::products(">>N")]
    #[case::components("C.O>N.Cl>F.Br")]
    #[case::labels("[C:7].[C:7]>[O:19]>[C:7]")]
    fn test_reaction_smiles_render(#[case] input: &str) {
        let reaction = ReactionSmiles::parse(input).unwrap();
        let original = reaction.clone();
        assert_eq!(reaction.render(), Ok(input.to_owned()));
        assert_eq!(reaction, original);
    }

    #[rstest]
    #[case::reactants(0, ReactionSmilesRenderError::Reactants(SmilesRenderError::AtomIndexOutOfBounds { atom: 2 }))]
    #[case::agents(1, ReactionSmilesRenderError::Agents(SmilesRenderError::AtomIndexOutOfBounds { atom: 2 }))]
    #[case::products(2, ReactionSmilesRenderError::Products(SmilesRenderError::AtomIndexOutOfBounds { atom: 2 }))]
    fn test_reaction_smiles_render_error(
        #[case] section: usize,
        #[case] expected: ReactionSmilesRenderError,
    ) {
        let mut table = Reaction::empty();
        let sections = [&mut table.reactants, &mut table.agents, &mut table.products];
        *sections.into_iter().nth(section).unwrap() = Molecule {
            atoms: vec![Atom::aliphatic_atom(Element::C)],
            bonds: vec![Bond::new(0, 2, BondOrder::Single)],
            ..Molecule::empty()
        };
        let reaction = ReactionSmiles::from_table_ir(table.clone());
        assert_eq!(reaction.as_table_ir(), &table);
        assert_eq!(reaction.render(), Err(expected));
        assert_eq!(reaction.into_table_ir(), table);
    }

    #[rstest]
    #[case::comments(Reaction { comments: vec!["note".to_owned()], ..Reaction::empty() }, "comments")]
    #[case::properties(Reaction { properties: [("name".to_owned(), "reaction".to_owned())].into(), ..Reaction::empty() }, "properties")]
    fn test_reaction_smiles_render_metadata(#[case] table: Reaction, #[case] field: &'static str) {
        assert_eq!(
            ReactionSmiles::from_table_ir(table).render(),
            Err(ReactionSmilesRenderError::UnsupportedReaction { field })
        );
    }

    #[rstest]
    fn test_reaction_smiles_render_with() {
        let reaction = ReactionSmiles::parse_with("C>N->B>O", &SmilesIoConfig::lenient()).unwrap();
        assert_eq!(
            reaction.render_with(&SmilesIoConfig::lenient()),
            Ok("C>N->B>O".to_owned())
        );
        assert_eq!(
            reaction.render(),
            Err(ReactionSmilesRenderError::Agents(
                SmilesRenderError::UnsupportedBond {
                    bond: 0,
                    field: "donation"
                }
            ))
        );
    }

    #[rstest]
    fn test_reaction_smiles_from_table_ir() {
        let mut table = ReactionSmiles::parse("[C:7].[C:7]>[O:19]>[C:7]")
            .unwrap()
            .into_table_ir();
        table.atom_mapping = [(31, (vec![900], vec![800]))].into();
        let reaction = ReactionSmiles::from_table_ir(table.clone());
        assert_eq!(reaction.as_table_ir(), &table);
        assert_eq!(reaction.render(), Ok("[C:7].[C:7]>[O:19]>[C:7]".to_owned()));
        assert_eq!(reaction.into_table_ir(), table);
    }
}
