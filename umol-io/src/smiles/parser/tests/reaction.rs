use std::collections::BTreeMap;

use map_macro::btree_map;
use pretty_assertions::assert_eq;
use rstest::*;
use umol_chem::element::Element;

use super::super::*;
use crate::table_ir::{
    Atom, Bond, BondDonation, BondOrder, ExtendedMolecule, ExtendedReaction, Molecule, Reaction,
    SourceFormat, Span,
};

#[rstest]
#[case::empty(b">>", Molecule::empty(), Molecule::empty(), Molecule::empty())]
#[case::reactants_only(
    b"C>>",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
    Molecule::empty(),
)]
#[case::products_only(
    b">>C",
    Molecule::empty(),
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
#[case::simple(
    b"C>>C",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
#[case::wildcards(
    b"*>*>*",
    Molecule { atoms: vec![Atom::wildcard_with_span(Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::wildcard_with_span(Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::wildcard_with_span(Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
)]
#[case::two_reactants(
    b"C.C>>C",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1)), Atom::aliphatic_atom_with_span(Element::C, Span::bytes(2, 3))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
#[case::two_products(
    b"C>>C.C",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1)), Atom::aliphatic_atom_with_span(Element::C, Span::bytes(2, 3))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
#[case::with_agent(
    b"C>CC>C",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1)), Atom::aliphatic_atom_with_span(Element::C, Span::bytes(1, 2))], bonds: vec![Bond { span: Some(Span::bytes(1, 2)), ..Bond::new(0, 1, BondOrder::Single) }], source_format: SourceFormat::SMILES, ..Molecule::empty() },
)]
#[case::agents_empty(
    b"C>>CC",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1)), Atom::aliphatic_atom_with_span(Element::C, Span::bytes(1, 2))], bonds: vec![Bond { span: Some(Span::bytes(1, 2)), ..Bond::new(0, 1, BondOrder::Single) }], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
fn test_parse_reaction(
    #[case] input: &[u8],
    #[case] reactants: Molecule,
    #[case] products: Molecule,
    #[case] agents: Molecule,
) {
    assert_eq!(
        parse_reaction(input, &SmilesIoConfig::opensmiles()),
        Ok(Reaction {
            reactants,
            products,
            agents,
            source_format: SourceFormat::SMILES,
            ..Reaction::empty()
        })
    );
}

#[rstest]
#[case::empty(b"", ParseError::MissingReactionArrow { pos: 0 })]
#[case::molecule_only(b"C", ParseError::MissingReactionArrow { pos: 1 })]
#[case::leading_whitespace(b" C>>C", ParseError::LeadingWhitespace)]
#[case::molecule_only_trailing_dot(b"C.", ParseError::TrailingDot { pos: 1 })]
#[case::molecule_only_one_gt(b"C.>", ParseError::TrailingDot { pos: 1 })]
#[case::agents_trailing_dot_without_products(b"C>C.", ParseError::TrailingDot { pos: 3 })]
#[case::one_gt(b"C>", ParseError::MissingReactionArrow { pos: 2 })]
#[case::one_lt(b"C<", ParseError::InvalidToken { pos: 1 })]
#[case::two_lt(b"C<<C", ParseError::InvalidToken { pos: 1 })]
#[case::gt_lt(b"C><O", ParseError::InvalidToken { pos: 2 })]
#[case::dative_bond(b"N->B>>N.B", ParseError::TrailingBond { pos: 1 })]
#[case::consecutive_dots_reactants(b"C..C>>C", ParseError::ConsecutiveDots { pos: 1 })]
#[case::consecutive_dots_products(b"C>>C..C", ParseError::ConsecutiveDots { pos: 4 })]
#[case::consecutive_dots_agents(b"C>C..C>C", ParseError::ConsecutiveDots { pos: 3 })]
#[case::leading_dot_reactants(b".C>>C", ParseError::LeadingDot { pos: 0 })]
#[case::leading_dot_products(b"C>>.C", ParseError::LeadingDot { pos: 3 })]
#[case::trailing_dot_reactants(b"C.>>C", ParseError::TrailingDot { pos: 1 })]
#[case::leading_dot_agents(b"C>.C>C", ParseError::LeadingDot { pos: 2 })]
#[case::trailing_dot_products(b"C>>C.", ParseError::TrailingDot { pos: 4 })]
#[case::trailing_dot_agents(b"C>C.>C", ParseError::TrailingDot { pos: 3 })]
fn test_parse_reaction_error(#[case] input: &[u8], #[case] expected: ParseError) {
    assert_eq!(
        parse_reaction(input, &SmilesIoConfig::opensmiles()),
        Err(expected)
    );
}

#[rstest]
#[case::dative_bond(
    b"N->B>>N.B",
    Reaction {
        reactants: Molecule {
            atoms: vec![
                Atom::aliphatic_atom_with_span(Element::N, Span::bytes(0, 1)),
                Atom::aliphatic_atom_with_span(Element::B, Span::bytes(3, 4)),
            ],
            bonds: vec![Bond {
                span: Some(Span::bytes(1, 2)),
                ..Bond::new_dative(0, 1, BondOrder::Single, BondDonation::Donating)
            }],
            source_format: SourceFormat::SMILES,
            ..Molecule::empty()
        },
        products: Molecule {
            atoms: vec![
                Atom::aliphatic_atom_with_span(Element::N, Span::bytes(0, 1)),
                Atom::aliphatic_atom_with_span(Element::B, Span::bytes(2, 3)),
            ],
            source_format: SourceFormat::SMILES,
            ..Molecule::empty()
        },
        agents: Molecule::empty(),
        source_format: SourceFormat::SMILES,
        ..Reaction::empty()
    }
)]
fn test_parse_reaction_lenient(#[case] input: &[u8], #[case] expected: Reaction) {
    assert_eq!(
        parse_reaction(input, &SmilesIoConfig::lenient()),
        Ok(expected)
    );
}

#[rstest]
#[case::mapped(b"[C:1][C:2]>>[C:1][C:2]", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![1], vec![1])))]
#[case::unmapped_reactant(b"[C:1]C>>[C:1][O:2]", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![], vec![1])))]
#[case::unmapped_product(b"[C:1][O:2]>>[C:1]C", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![1], vec![])))]
#[case::partial_mapping(b"[C:1]C>>C[O:2]", btree_map!(1 => (vec![0], vec![]), 2 => (vec![], vec![1])))]
#[case::duplicate_mapping(b"[C:1](=[O:2])[OH:2]>>[C:1](=[O:2])[O:2]C", btree_map!(1 => (vec![0], vec![0]),
    2 => (vec![1, 2], vec![1, 2])))]
fn test_parse_reaction_atom_mapping(
    #[case] input: &[u8],
    #[case] expected: BTreeMap<u32, (Vec<u32>, Vec<u32>)>,
) {
    assert_eq!(
        parse_reaction(input, &SmilesIoConfig::opensmiles()).map(|reaction| reaction.atom_mapping),
        Ok(expected)
    );
}

#[rstest]
#[case::empty(b">>", Molecule::empty(), Molecule::empty(), Molecule::empty())]
#[case::reactants_only(
    b"C>>",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
    Molecule::empty(),
)]
#[case::products_only(
    b">>C",
    Molecule::empty(),
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
#[case::simple(
    b"C>>C",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
#[case::wildcard(
    b"*>>*",
    Molecule { atoms: vec![Atom::wildcard_with_span(Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::wildcard_with_span(Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
#[case::two_reactants(
    b"C.C>>C",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1)), Atom::aliphatic_atom_with_span(Element::C, Span::bytes(2, 3))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
#[case::two_products(
    b"C>>C.C",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1)), Atom::aliphatic_atom_with_span(Element::C, Span::bytes(2, 3))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
#[case::with_agent(
    b"C>CC>C",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1)), Atom::aliphatic_atom_with_span(Element::C, Span::bytes(1, 2))], bonds: vec![Bond { span: Some(Span::bytes(1, 2)), ..Bond::new(0, 1, BondOrder::Single) }], source_format: SourceFormat::SMILES, ..Molecule::empty() },
)]
#[case::agents_empty(
    b"C>>CC",
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule { atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1)), Atom::aliphatic_atom_with_span(Element::C, Span::bytes(1, 2))], bonds: vec![Bond { span: Some(Span::bytes(1, 2)), ..Bond::new(0, 1, BondOrder::Single) }], source_format: SourceFormat::SMILES, ..Molecule::empty() },
    Molecule::empty(),
)]
fn test_parse_extended_reaction_smiles_bytes(
    #[case] input: &[u8],
    #[case] reactants: Molecule,
    #[case] products: Molecule,
    #[case] agents: Molecule,
) {
    assert_eq!(
        parse_extended_reaction_smiles_bytes(input),
        Ok(ExtendedReaction {
            reactants: ExtendedMolecule::from(reactants),
            products: ExtendedMolecule::from(products),
            agents: ExtendedMolecule::from(agents),
            source_format: SourceFormat::SMILES,
            ..ExtendedReaction::empty()
        })
    );
}

#[rstest]
#[case::empty(b"", ParseError::MissingReactionArrow { pos: 0 })]
#[case::molecule_only(b"C", ParseError::MissingReactionArrow { pos: 1 })]
#[case::leading_whitespace(b" C>>C", ParseError::LeadingWhitespace)]
#[case::molecule_only_trailing_dot(b"C.", ParseError::TrailingDot { pos: 1 })]
#[case::molecule_only_one_gt(b"C.>", ParseError::TrailingDot { pos: 1 })]
#[case::agents_trailing_dot_without_products(b"C>C.", ParseError::TrailingDot { pos: 1 })]
#[case::one_gt(b"C>", ParseError::MissingReactionArrow { pos: 2 })]
#[case::one_lt(b"C<", ParseError::InvalidToken { pos: 1 })]
#[case::two_lt(b"C<<C", ParseError::InvalidToken { pos: 1 })]
#[case::gt_lt(b"C><O", ParseError::InvalidToken { pos: 2 })]
#[case::dative_bond(b"N->B>>N.B", ParseError::TrailingBond { pos: 1 })]
#[case::consecutive_dots_reactants(b"C..C>>C", ParseError::ConsecutiveDots { pos: 1 })]
#[case::consecutive_dots_products(b"C>>C..C", ParseError::ConsecutiveDots { pos: 4 })]
#[case::consecutive_dots_agents(b"C>C..C>C", ParseError::ConsecutiveDots { pos: 3 })]
#[case::leading_dot_reactants(b".C>>C", ParseError::LeadingDot { pos: 0 })]
#[case::leading_dot_products(b"C>>.C", ParseError::LeadingDot { pos: 0 })]
#[case::trailing_dot_reactants(b"C.>>C", ParseError::TrailingDot { pos: 1 })]
#[case::leading_dot_agents(b"C>.C>C", ParseError::LeadingDot { pos: 0 })]
#[case::trailing_dot_products(b"C>>C.", ParseError::TrailingDot { pos: 1 })]
#[case::trailing_dot_agents(b"C>C.>C", ParseError::TrailingDot { pos: 3 })]
fn test_parse_extended_reaction_smiles_bytes_error(
    #[case] input: &[u8],
    #[case] expected: ParseError,
) {
    assert_eq!(parse_extended_reaction_smiles_bytes(input), Err(expected));
}

#[rstest]
#[case::dative_bond(
    b"N->B>>N.B",
    ExtendedReaction {
        reactants: ExtendedMolecule::from(Molecule {
            atoms: vec![
                Atom::aliphatic_atom_with_span(Element::N, Span::bytes(0, 1)),
                Atom::aliphatic_atom_with_span(Element::B, Span::bytes(3, 4)),
            ],
            bonds: vec![Bond {
                span: Some(Span::bytes(1, 2)),
                ..Bond::new_dative(0, 1, BondOrder::Single, BondDonation::Donating)
            }],
            source_format: SourceFormat::SMILES,
            ..Molecule::empty()
        }),
        products: ExtendedMolecule::from(Molecule {
            atoms: vec![
                Atom::aliphatic_atom_with_span(Element::N, Span::bytes(0, 1)),
                Atom::aliphatic_atom_with_span(Element::B, Span::bytes(2, 3)),
            ],
            source_format: SourceFormat::SMILES,
            ..Molecule::empty()
        }),
        agents: ExtendedMolecule::empty(),
        source_format: SourceFormat::SMILES,
        ..ExtendedReaction::empty()
    }
)]
fn test_parse_extended_reaction_smiles_bytes_with(
    #[case] input: &[u8],
    #[case] expected: ExtendedReaction,
) {
    assert_eq!(
        parse_extended_reaction_smiles_bytes_with(input, &SmilesIoConfig::lenient()),
        Ok(expected)
    );
}

#[rstest]
#[case::mapped(b"[C:1][C:2]>>[C:1][C:2]", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![1], vec![1])))]
#[case::unmapped_reactant(b"[C:1]C>>[C:1][O:2]", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![], vec![1])))]
#[case::unmapped_product(b"[C:1][O:2]>>[C:1]C", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![1], vec![])))]
#[case::partial_mapping(b"[C:1]C>>C[O:2]", btree_map!(1 => (vec![0], vec![]), 2 => (vec![], vec![1])))]
#[case::duplicate_mapping(b"[C:1](=[O:2])[OH:2]>>[C:1](=[O:2])[O:2]C", btree_map!(1 => (vec![0], vec![0]),
    2 => (vec![1, 2], vec![1, 2])))]
fn test_parse_extended_reaction_smiles_bytes_atom_mapping(
    #[case] input: &[u8],
    #[case] expected: BTreeMap<u32, (Vec<u32>, Vec<u32>)>,
) {
    assert_eq!(
        parse_extended_reaction_smiles_bytes(input).map(|reaction| reaction.atom_mapping),
        Ok(expected)
    );
}

#[rstest]
fn test_parse_reaction_cx() {
    assert_eq!(
        parse_reaction(b"C>CC>C |$r;a0;a1;p$|", &SmilesIoConfig::chemaxon()).map(|reaction| {
            (
                reaction
                    .reactants
                    .atoms
                    .into_iter()
                    .map(|atom| atom.label)
                    .collect::<Vec<_>>(),
                reaction
                    .products
                    .atoms
                    .into_iter()
                    .map(|atom| atom.label)
                    .collect::<Vec<_>>(),
                reaction
                    .agents
                    .atoms
                    .into_iter()
                    .map(|atom| atom.label)
                    .collect::<Vec<_>>(),
            )
        }),
        Ok((
            vec![Some(String::from("r"))],
            vec![Some(String::from("p"))],
            vec![Some(String::from("a0")), Some(String::from("a1"))],
        ))
    );
}

#[rstest]
fn test_parse_reaction_cx_error() {
    assert_eq!(
        parse_reaction(b"C>>C |$a;b;c$|", &SmilesIoConfig::chemaxon()),
        Err(ParseError::AtomIndexOutOfBounds { atom_idx: 2 })
    );
}

#[rstest]
fn test_parse_extended_reaction_smiles_bytes_with_cx() {
    assert_eq!(
        parse_extended_reaction_smiles_bytes_with(b"C.C>>C |f:0.1|", &SmilesIoConfig::chemaxon())
            .map(|reaction| { reaction.reactants.cx_data.and_then(|data| data.components) }),
        Ok(Some(vec![vec![0, 1]]))
    );
}
