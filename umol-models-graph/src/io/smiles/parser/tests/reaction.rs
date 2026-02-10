use std::collections::BTreeMap;

use bstr::ByteSlice;
use map_macro::btree_map;
use pretty_assertions::assert_eq;
use rstest::*;

use super::super::*;
use super::utils::build_from_graph;
use crate::table_ir::Reaction;

fn reaction_from_sides(
    reactants: impl IntoIterator<Item = &'static str>,
    products: impl IntoIterator<Item = &'static str>,
    agents: impl IntoIterator<Item = &'static str>,
) -> Reaction {
    let reactants: Vec<_> = reactants.into_iter().map(|s| build_from_graph(s)).collect();
    let products: Vec<_> = products.into_iter().map(|s| build_from_graph(s)).collect();
    let agents: Vec<_> = agents.into_iter().map(|s| build_from_graph(s)).collect();
    Reaction::from_molecules(reactants, products, agents)
}

#[rstest]
#[case::empty(b">>", reaction_from_sides([], [], []))]
#[case::reactants_only(b"C>>", reaction_from_sides(["C@0 |"], [], []))]
#[case::products_only(b">>C", reaction_from_sides([], ["C@0 |"], []))]
#[case::simple(b"C>>C", reaction_from_sides(["C@0 |"], ["C@0 |"], []))]
#[case::two_reactants(b"C.C>>C", reaction_from_sides(["C@0 |", "C@0 |"], ["C@0 |"], []))]
#[case::two_products(b"C>>C.C", reaction_from_sides(["C@0 |"], ["C@0 |", "C@0 |"], []))]
#[case::with_agent(b"C>CC>C", reaction_from_sides(["C@0 |"], ["C@0 |"], ["C@0 C@1 | 0-1@2"]))]
#[case::agents_empty(b"C>>CC", reaction_from_sides(["C@0 |"], ["C@0 C@1 | 0-1@2"], []))]
fn parse_reaction_smiles(#[case] input: &[u8], #[case] expected: Reaction) {
    let res = parse_reaction_smiles_bytes(input);
    assert!(
        res.is_ok(),
        "{:?} should succeed: {:?}",
        input.to_str_lossy(),
        res
    );
    let rxn = res.unwrap();
    assert_eq!(rxn.reactants.len(), expected.reactants.len());
    assert_eq!(rxn.products.len(), expected.products.len());
    assert_eq!(rxn.agents.len(), expected.agents.len());
    for (a, b) in rxn.reactants.iter().zip(expected.reactants.iter()) {
        assert_eq!(a.atoms.len(), b.atoms.len());
    }
    for (a, b) in rxn.products.iter().zip(expected.products.iter()) {
        assert_eq!(a.atoms.len(), b.atoms.len());
    }
    for (a, b) in rxn.agents.iter().zip(expected.agents.iter()) {
        assert_eq!(a.atoms.len(), b.atoms.len());
    }
}

#[rstest]
#[case::empty(b"", ParseError::MissingReactionArrow { pos: 0 })]
#[case::molecule_only(b"C", ParseError::MissingReactionArrow { pos: 1 })]
#[case::leading_whitespace(b" C>>C", ParseError::LeadingWhitespace)]
#[case::molecule_only_trailing_dot(b"C.", ParseError::MissingReactionArrow { pos: 1 })]
#[case::molecule_only_one_gt(b"C.>", ParseError::TrailingDot { pos: 1 })]
#[case::no_products_trailing_dot(b"C>C.", ParseError::MissingReactionArrow { pos: 3 })]
#[case::one_gt(b"C>", ParseError::MissingReactionArrow { pos: 1 })]
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
#[case::dative_bond(b"N->B>>N.B", reaction_from_sides(["N@0 B@3 | 0-1->@1"], ["N@0 |", "B@0 |"], []))]
fn parse_reaction_smiles_bytes_lenient(#[case] input: &[u8], #[case] expected: Reaction) {
    let res = parse_reaction_smiles_bytes_with(
        input,
        &SmilesIoConfig::with_parse_flags(SmilesParseFlags::BASIC_MAX & SmilesParseFlags::LENIENT),
    );
    assert!(res.is_ok(), "{:?} should succeed", input.to_str_lossy());
    let rxn = res.unwrap();
    assert_eq!(rxn.reactants.len(), expected.reactants.len());
    assert_eq!(rxn.products.len(), expected.products.len());
    assert_eq!(rxn.agents.len(), expected.agents.len());
    for (a, b) in rxn.reactants.iter().zip(expected.reactants.iter()) {
        assert_eq!(a.atoms.len(), b.atoms.len());
    }
    for (a, b) in rxn.products.iter().zip(expected.products.iter()) {
        assert_eq!(a.atoms.len(), b.atoms.len());
    }
    for (a, b) in rxn.agents.iter().zip(expected.agents.iter()) {
        assert_eq!(a.atoms.len(), b.atoms.len());
    }
}

#[rstest]
#[case::mapped(b"[C:1][C:2]>>[C:1][C:2]", btree_map!(1 => (vec![(0, 0)], vec![(0, 0)]), 2 => (vec![(0, 1)], vec![(0, 1)])))]
#[case::unmapped_reactant(b"[C:1]C>>[C:1][O:2]", btree_map!(1 => (vec![(0, 0)], vec![(0, 0)]), 2 => (vec![], vec![(0, 1)])))]
#[case::unmapped_product(b"[C:1][O:2]>>[C:1]C", btree_map!(1 => (vec![(0, 0)], vec![(0, 0)]), 2 => (vec![(0, 1)], vec![])))]
#[case::partial_mapping(b"[C:1]C>>C[O:2]", btree_map!(1 => (vec![(0, 0)], vec![]), 2 => (vec![], vec![(0, 1)])))]
#[case::duplicate_mapping(b"[C:1](=[O:2])[OH:2]>>[C:1](=[O:2])[O:2]C", btree_map!(1 => (vec![(0, 0)], vec![(0, 0)]),
    2 => (vec![(0, 1), (0, 2)], vec![(0, 1), (0, 2)])))]
fn parse_reaction_smiles_atom_mapping(
    #[case] input: &[u8],
    #[case] expected: BTreeMap<u32, (Vec<(usize, usize)>, Vec<(usize, usize)>)>,
) {
    let res = parse_reaction_smiles_bytes(input);
    assert!(res.is_ok(), "{:?} should succeed", input.to_str_lossy());
    let rxn = res.unwrap();
    assert_eq!(rxn.atom_mapping, expected);
}
