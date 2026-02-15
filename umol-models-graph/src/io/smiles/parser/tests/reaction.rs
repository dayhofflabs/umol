use std::collections::BTreeMap;

use bstr::ByteSlice;
use map_macro::btree_map;
use pretty_assertions::assert_eq;
use rstest::*;

use super::super::*;
use super::utils::{build_extended_from_graph, build_from_graph};
use crate::table_ir::{ExtendedReaction, Reaction};

fn reaction_from_sides(
    reactants: &'static str,
    products: &'static str,
    agents: &'static str,
) -> Reaction {
    let reactants = build_from_graph(reactants);
    let products = build_from_graph(products);
    let agents = build_from_graph(agents);
    Reaction::from_molecules(reactants, products, agents)
}

fn extended_reaction_from_sides(
    reactants: &'static str,
    products: &'static str,
    agents: &'static str,
) -> ExtendedReaction {
    let reactants = build_extended_from_graph(reactants);
    let products = build_extended_from_graph(products);
    let agents = build_extended_from_graph(agents);
    ExtendedReaction::from_extended_molecules(reactants, products, agents)
}

#[rstest]
#[case::empty(b">>", reaction_from_sides("|", "|", "|"))]
#[case::reactants_only(b"C>>", reaction_from_sides("C@0 |", "|", "|"))]
#[case::products_only(b">>C", reaction_from_sides("|", "C@0 |", "|"))]
#[case::simple(b"C>>C", reaction_from_sides("C@0 |", "C@0 |", "|"))]
#[case::two_reactants(b"C.C>>C", reaction_from_sides("C@0 C@1 |", "C@0 |", "|"))]
#[case::two_products(b"C>>C.C", reaction_from_sides("C@0 |", "C@0 C@1 |", "|"))]
#[case::with_agent(b"C>CC>C", reaction_from_sides("C@0 |", "C@0 |", "C@0 C@1 | 0-1@2"))]
#[case::agents_empty(b"C>>CC", reaction_from_sides("C@0 |", "C@0 C@1 | 0-1@2", "|"))]
fn parse_reaction_smiles(#[case] input: &[u8], #[case] expected: Reaction) {
    let res = parse_reaction_smiles_bytes(input);
    assert!(
        res.is_ok(),
        "{:?} should succeed: {:?}",
        input.to_str_lossy(),
        res
    );
    let rxn = res.unwrap();
    assert_eq!(rxn.reactants.atoms.len(), expected.reactants.atoms.len());
    assert_eq!(rxn.products.atoms.len(), expected.products.atoms.len());
    assert_eq!(rxn.agents.atoms.len(), expected.agents.atoms.len());
}

#[rstest]
#[case::empty(b"", ParseError::MissingReactionArrow { pos: 0 })]
#[case::molecule_only(b"C", ParseError::MissingReactionArrow { pos: 1 })]
#[case::leading_whitespace(b" C>>C", ParseError::LeadingWhitespace)]
#[case::molecule_only_trailing_dot(b"C.", ParseError::TrailingDot { pos: 1 })]
#[case::molecule_only_one_gt(b"C.>", ParseError::TrailingDot { pos: 1 })]
#[case::no_products_trailing_dot(b"C>C.", ParseError::TrailingDot { pos: 3 })]
#[case::one_gt(b"C>", ParseError::MissingReactionArrow { pos: 2 })]
#[case::one_lt(b"C<", ParseError::InvalidToken { pos: 1 })]
#[case::two_lt(b"C<<C", ParseError::InvalidToken { pos: 1 })]
#[case::gt_lt(b"C><O", ParseError::InvalidToken { pos: 2 })]
#[case::dative_bond(b"N->B>>N.B", ParseError::TrailingBond { pos: 1 })]
#[case::consecutive_dots_reactants(b"C..C>>C", ParseError::ConsecutiveDots { pos: 1 })]
#[case::consecutive_dots_products(b"C>>C..C", ParseError::ConsecutiveDots { pos: 4 })]
#[case::consecutive_dots_agents(b"C>C..C>C", ParseError::ConsecutiveDots { pos: 3 })]
#[case::leading_dot_reactants(b".C>>C", ParseError::LeadingDot { pos: 0 })]
#[case::leading_dot_products(b".C>>C", ParseError::LeadingDot { pos: 0 })]
#[case::trailing_dot_reactants(b"C.>>C", ParseError::TrailingDot { pos: 1 })]
#[case::leading_dot_agents(b"C>.C>C", ParseError::LeadingDot { pos: 2 })]
#[case::trailing_dot_products(b"C>>C.", ParseError::TrailingDot { pos: 4 })]
#[case::trailing_dot_agents(b"C>C.>C", ParseError::TrailingDot { pos: 3 })]
#[case::wildcard(b"*>>*", ParseError::InvalidElement { pos: 0 })]
fn parse_reaction_smiles_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_reaction_smiles_bytes(input);
    assert!(
        res.is_err(),
        "{:?} should have failed, got: {:?}",
        input.to_str_lossy(),
        res.unwrap()
    );
    assert_eq!(res.unwrap_err(), expected);
}

#[rstest]
#[case::dative_bond(b"N->B>>N.B", reaction_from_sides("N@0 B@3 | 0-1->@1", "N@0 B@1 |", "|"))]
fn parse_reaction_smiles_bytes_lenient(#[case] input: &[u8], #[case] expected: Reaction) {
    let res = parse_reaction_smiles_bytes_with(
        input,
        &SmilesIoConfig::with_parse_flags(SmilesParseFlags::BASIC_MAX & SmilesParseFlags::LENIENT),
    );
    assert!(res.is_ok(), "{:?} should succeed", input.to_str_lossy());
    let rxn = res.unwrap();
    assert_eq!(rxn.reactants.atoms.len(), expected.reactants.atoms.len());
    assert_eq!(rxn.products.atoms.len(), expected.products.atoms.len());
    assert_eq!(rxn.agents.atoms.len(), expected.agents.atoms.len());
}

#[rstest]
#[case::mapped(b"[C:1][C:2]>>[C:1][C:2]", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![1], vec![1])))]
#[case::unmapped_reactant(b"[C:1]C>>[C:1][O:2]", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![], vec![1])))]
#[case::unmapped_product(b"[C:1][O:2]>>[C:1]C", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![1], vec![])))]
#[case::partial_mapping(b"[C:1]C>>C[O:2]", btree_map!(1 => (vec![0], vec![]), 2 => (vec![], vec![1])))]
#[case::duplicate_mapping(b"[C:1](=[O:2])[OH:2]>>[C:1](=[O:2])[O:2]C", btree_map!(1 => (vec![0], vec![0]),
    2 => (vec![1, 2], vec![1, 2])))]
fn parse_reaction_smiles_atom_mapping(
    #[case] input: &[u8],
    #[case] expected: BTreeMap<u32, (Vec<u32>, Vec<u32>)>,
) {
    let res = parse_reaction_smiles_bytes(input);
    assert!(res.is_ok(), "{:?} should succeed", input.to_str_lossy());
    let rxn = res.unwrap();
    assert_eq!(rxn.atom_mapping, expected);
}

#[rstest]
#[case::empty(b">>", extended_reaction_from_sides("|", "|", "|"))]
#[case::reactants_only(b"C>>", extended_reaction_from_sides("C@0 |", "|", "|"))]
#[case::products_only(b">>C", extended_reaction_from_sides("|", "C@0 |", "|"))]
#[case::simple(b"C>>C", extended_reaction_from_sides("C@0 |", "C@0 |", "|"))]
#[case::wildcard(b"*>>*", extended_reaction_from_sides("*@0 |", "*@0 |", "|"))]
#[case::two_reactants(b"C.C>>C", extended_reaction_from_sides("C@0 C@1 |", "C@0 |", "|"))]
#[case::two_products(b"C>>C.C", extended_reaction_from_sides("C@0 |", "C@0 C@1 |", "|"))]
#[case::with_agent(b"C>CC>C", extended_reaction_from_sides("C@0 |", "C@0 |", "C@0 C@1 | 0-1@2"))]
#[case::agents_empty(b"C>>CC", extended_reaction_from_sides("C@0 |", "C@0 C@1 | 0-1@2", "|"))]
fn parse_extended_reaction_smiles(#[case] input: &[u8], #[case] expected: ExtendedReaction) {
    let res = parse_extended_reaction_smiles_bytes(input);
    assert!(
        res.is_ok(),
        "{:?} should succeed: {:?}",
        input.to_str_lossy(),
        res
    );
    let rxn = res.unwrap();
    assert_eq!(rxn.reactants.atoms.len(), expected.reactants.atoms.len());
    assert_eq!(rxn.products.atoms.len(), expected.products.atoms.len());
    assert_eq!(rxn.agents.atoms.len(), expected.agents.atoms.len());
}

#[rstest]
#[case::empty(b"", ParseError::MissingReactionArrow { pos: 0 })]
#[case::molecule_only(b"C", ParseError::MissingReactionArrow { pos: 1 })]
#[case::leading_whitespace(b" C>>C", ParseError::LeadingWhitespace)]
#[case::molecule_only_trailing_dot(b"C.", ParseError::TrailingDot { pos: 1 })]
#[case::molecule_only_one_gt(b"C.>", ParseError::TrailingDot { pos: 1 })]
#[case::no_products_trailing_dot(b"C>C.", ParseError::TrailingDot { pos: 3 })]
#[case::one_gt(b"C>", ParseError::MissingReactionArrow { pos: 2 })]
#[case::one_lt(b"C<", ParseError::InvalidToken { pos: 1 })]
#[case::two_lt(b"C<<C", ParseError::InvalidToken { pos: 1 })]
#[case::gt_lt(b"C><O", ParseError::InvalidToken { pos: 2 })]
#[case::dative_bond(b"N->B>>N.B", ParseError::TrailingBond { pos: 1 })]
#[case::consecutive_dots_reactants(b"C..C>>C", ParseError::ConsecutiveDots { pos: 1 })]
#[case::consecutive_dots_products(b"C>>C..C", ParseError::ConsecutiveDots { pos: 4 })]
#[case::consecutive_dots_agents(b"C>C..C>C", ParseError::ConsecutiveDots { pos: 3 })]
#[case::leading_dot_reactants(b".C>>C", ParseError::LeadingDot { pos: 0 })]
#[case::leading_dot_products(b".C>>C", ParseError::LeadingDot { pos: 0 })]
#[case::trailing_dot_reactants(b"C.>>C", ParseError::TrailingDot { pos: 1 })]
#[case::leading_dot_agents(b"C>.C>C", ParseError::LeadingDot { pos: 2 })]
#[case::trailing_dot_products(b"C>>C.", ParseError::TrailingDot { pos: 4 })]
#[case::trailing_dot_agents(b"C>C.>C", ParseError::TrailingDot { pos: 3 })]
fn parse_extended_reaction_smiles_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_reaction_smiles_bytes(input);
    assert!(
        res.is_err(),
        "{:?} should have failed, got: {:?}",
        input.to_str_lossy(),
        res.unwrap()
    );
    assert_eq!(res.unwrap_err(), expected);
}

#[rstest]
#[case::dative_bond(b"N->B>>N.B", extended_reaction_from_sides("N@0 B@3 | 0-1->@1", "N@0 B@1 |", "|"))]
fn parse_extended_reaction_smiles_bytes_lenient(
    #[case] input: &[u8],
    #[case] expected: ExtendedReaction,
) {
    let res = parse_extended_reaction_smiles_bytes_with(
        input,
        &SmilesIoConfig::with_parse_flags(SmilesParseFlags::BASIC_MAX & SmilesParseFlags::LENIENT),
    );
    assert!(res.is_ok(), "{:?} should succeed", input.to_str_lossy());
    let rxn = res.unwrap();
    assert_eq!(rxn.reactants.atoms.len(), expected.reactants.atoms.len());
    assert_eq!(rxn.products.atoms.len(), expected.products.atoms.len());
    assert_eq!(rxn.agents.atoms.len(), expected.agents.atoms.len());
}

#[rstest]
#[case::mapped(b"[C:1][C:2]>>[C:1][C:2]", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![1], vec![1])))]
#[case::unmapped_reactant(b"[C:1]C>>[C:1][O:2]", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![], vec![1])))]
#[case::unmapped_product(b"[C:1][O:2]>>[C:1]C", btree_map!(1 => (vec![0], vec![0]), 2 => (vec![1], vec![])))]
#[case::partial_mapping(b"[C:1]C>>C[O:2]", btree_map!(1 => (vec![0], vec![]), 2 => (vec![], vec![1])))]
#[case::duplicate_mapping(b"[C:1](=[O:2])[OH:2]>>[C:1](=[O:2])[O:2]C", btree_map!(1 => (vec![0], vec![0]),
    2 => (vec![1, 2], vec![1, 2])))]
fn parse_extended_reaction_smiles_atom_mapping(
    #[case] input: &[u8],
    #[case] expected: BTreeMap<u32, (Vec<u32>, Vec<u32>)>,
) {
    let res = parse_extended_reaction_smiles_bytes(input);
    assert!(res.is_ok(), "{:?} should succeed", input.to_str_lossy());
    let rxn = res.unwrap();
    assert_eq!(rxn.atom_mapping, expected);
}

#[test]
fn parse_reaction_smiles_trailing_cx_labels_global_indices() {
    let rxn = parse_reaction_smiles_bytes_with(
        b"C>CC>C |$r;a0;a1;p$|",
        &SmilesIoConfig::with_parse_flags(SmilesParseFlags::BASIC_MAX),
    )
    .unwrap();
    assert_eq!(rxn.reactants.atoms[0].label.as_deref(), Some("r"));
    assert_eq!(rxn.agents.atoms[0].label.as_deref(), Some("a0"));
    assert_eq!(rxn.agents.atoms[1].label.as_deref(), Some("a1"));
    assert_eq!(rxn.products.atoms[0].label.as_deref(), Some("p"));
}

#[test]
fn parse_reaction_smiles_trailing_cx_invalid_index() {
    let err = parse_reaction_smiles_bytes_with(
        b"C>>C |$a;b;c$|",
        &SmilesIoConfig::with_parse_flags(SmilesParseFlags::BASIC_MAX),
    )
    .unwrap_err();
    assert_eq!(err, ParseError::AtomIndexOutOfBounds { atom_idx: 2 });
}

#[test]
fn parse_extended_reaction_smiles_trailing_cx_fragment_groups() {
    let rxn = parse_extended_reaction_smiles_bytes_with(
        b"C.C>>C |f:0.1|",
        &SmilesIoConfig::with_parse_flags(SmilesParseFlags::BASIC_MAX),
    )
    .unwrap();
    let components = rxn
        .reactants
        .cx_data
        .as_ref()
        .and_then(|d| d.components.clone())
        .unwrap();
    assert_eq!(components, vec![vec![0, 1]]);
}
