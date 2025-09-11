use pretty_assertions::assert_eq;
use rstest::*;
use slog::o;
use umol::logging::setup_logger;
use umol::with_logger;
use umol_data::Element;

use crate::io::ir::{Atom, Bond, BondOrder, Chirality};
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::{chain_accept, chain, grammar, tree, tree_accept};
use crate::io::smiles::state::{Bond as SBond, ParseState};

#[fixture]
fn parse_state() -> ParseState {
    let mut state = ParseState::default();
    let root = setup_logger(slog::Level::Debug);
    state.log = Some(with_logger!(root, "io::smiles::parser"));
    state
}

#[rstest]
#[case("C", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.index = Some(0); atom})]
fn test_chain_accept_atom(mut parse_state: ParseState, #[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let parser = chain_accept::AtomParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let result = result.unwrap();
    assert_eq!(result, expected);
    assert_eq!(parse_state.current_atom, 0);
    assert_eq!(parse_state.next_atom, 1);
    assert_eq!(parse_state.next_bond, 0);
    assert!(parse_state.pending_bond.is_none());
    assert!(parse_state.error.is_none());
}

#[rstest]
#[case("=", Bond::from_order(BondOrder::Double))]
fn test_chain_accept_bond(mut parse_state: ParseState, #[case] input: &str, #[case] expected: Bond) {
    let lexer = Lexer::new(input);
    parse_state.current_atom = 0;
    parse_state.next_atom = 1;
    let parser = chain_accept::BondParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let result = result.unwrap();
    assert_eq!(result, expected);
    assert_eq!(parse_state.current_atom, 0);
    assert_eq!(parse_state.next_atom, 1);
    assert_eq!(parse_state.next_bond, 0);
    assert_eq!(parse_state.pending_bond, Some(SBond { order: BondOrder::Double, dir: None }));
}

#[rstest]
#[case("CC", 2, 1)]
#[case("CCC", 3, 2)]
#[case("C=C", 2, 1)]
#[case("C#C", 2, 1)]
#[case("C$C", 2, 1)]
#[case("C=CC", 3, 2)]
#[case("CC=C", 3, 2)]
#[case("C=C=C", 3, 2)]
#[case("CCCCCCCC", 8, 7)]
fn test_chain_accept(mut parse_state: ParseState, #[case] input: &str, #[case] atoms: usize, #[case] bonds: usize) {
    let lexer = Lexer::new(input);
    let parser = chain_accept::ChainParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(parse_state.next_atom, atoms);
    assert_eq!(parse_state.next_bond, bonds);
    assert!(parse_state.pending_bond.is_none());
    assert!(parse_state.error.is_none());
}

#[rstest]
#[case("=C", "leading bond")]
#[case("C=", "trailing bond")]
#[case("C==C", "consecutive bonds")]
fn test_chain_accept_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = chain_accept::ChainParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case("C")]
fn test_chain_atom(mut parse_state: ParseState, #[case] input: &str) {
    let lexer = Lexer::new(input);
    let parser = chain::AtomBondParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(parse_state.next_atom, 1);
    assert_eq!(parse_state.next_bond, 0);
    assert!(parse_state.pending_bond.is_none());
    assert!(parse_state.error.is_none());
    // Finalize and verify molecule IR
    parse_state.finalize_current_molecule();
    let mols = parse_state.take_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), 1);
    assert_eq!(mols[0].bonds.len(), 0);
}

#[rstest]
#[case("=")]
fn test_chain_bond(mut parse_state: ParseState, #[case] input: &str) {
    let lexer = Lexer::new(input);
    parse_state.current_atom = 0;
    parse_state.next_atom = 1;
    let parser = chain::BondParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(parse_state.pending_bond, Some(SBond { order: BondOrder::Double, dir: None }));
}

#[rstest]
#[case("CC", 2, 1)]
#[case("CCC", 3, 2)]
#[case("C=C", 2, 1)]
#[case("C#C", 2, 1)]
#[case("C$C", 2, 1)]
#[case("C=CC", 3, 2)]
#[case("CC=C", 3, 2)]
#[case("C=C=C", 3, 2)]
#[case("CCCCCCCC", 8, 7)]
fn test_chain(mut parse_state: ParseState, #[case] input: &str, #[case] atoms: usize, #[case] bonds: usize) {
    let lexer = Lexer::new(input);
    let parser = chain::ChainParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(parse_state.next_atom, atoms);
    assert_eq!(parse_state.next_bond, bonds);
    assert!(parse_state.pending_bond.is_none());
    assert!(parse_state.error.is_none());
    parse_state.finalize_current_molecule();
    let mols = parse_state.take_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), atoms);
    assert_eq!(mols[0].bonds.len(), bonds);
}

#[rstest]
#[case("=C", "leading bond")]
#[case("C=", "trailing bond")]
#[case("C==C", "consecutive bonds")]
fn test_chain_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = chain::ChainParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case("C", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.index = Some(0); atom})]
fn test_tree_accept_atom(mut parse_state: ParseState, #[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let parser = tree_accept::AtomParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let result = result.unwrap();
    assert_eq!(result, expected);
    assert_eq!(parse_state.current_atom, 0);
    assert_eq!(parse_state.next_atom, 1);
    assert_eq!(parse_state.next_bond, 0);
    assert!(parse_state.pending_bond.is_none());
    assert!(parse_state.error.is_none());
}

#[rstest]
#[case("=", Bond::from_order(BondOrder::Double))]
fn test_tree_accept_bond(mut parse_state: ParseState, #[case] input: &str, #[case] expected: Bond) {
    let lexer = Lexer::new(input);
    parse_state.current_atom = 0;
    parse_state.next_atom = 1;
    let parser = tree_accept::BondParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let result = result.unwrap();
    assert_eq!(result, expected);
    assert_eq!(parse_state.current_atom, 0);
    assert_eq!(parse_state.next_atom, 1);
    assert_eq!(parse_state.next_bond, 0);
    assert_eq!(parse_state.pending_bond, Some(SBond { order: BondOrder::Double, dir: None }));
}

#[rstest]
#[case("CC", 2, 1)]
#[case("C=C", 2, 1)]
fn test_tree_accept_chain(mut parse_state: ParseState, #[case] input: &str, #[case] atoms: usize, #[case] bonds: usize) {
    let lexer = Lexer::new(input);
    let parser = tree_accept::ChainParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(parse_state.next_atom, atoms);
    assert_eq!(parse_state.next_bond, bonds);
}

#[rstest]
#[case("CC(C)C", 4, 3)]
#[case("CC(C)", 3, 2)]
#[case("C(C)C", 3, 2)]
#[case("CC(=C)C", 4, 3)]
#[case("CC(C)(C)CC", 6, 5)]
#[case("CC(C(C)C)C", 6, 5)]
#[case("CC(CC)C", 5, 4)]
fn test_tree_accept(mut parse_state: ParseState, #[case] input: &str, #[case] atoms: usize, #[case] bonds: usize) {
    let lexer = Lexer::new(input);
    let parser = tree_accept::TreeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(parse_state.next_atom, atoms);
    assert_eq!(parse_state.next_bond, bonds);
    assert!(parse_state.pending_bond.is_none());
    assert!(parse_state.error.is_none());
}

#[rstest]
#[case("(C)CC", "leading branch with no parent")]
#[case("CC(C", "unclosed branch")]
#[case("C(C-)C", "trailing bond in branch")]
#[case("CC)C", "unmatched close paren")]
#[case("C()C", "empty branch")]
fn test_tree_accept_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = tree_accept::TreeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("c", Atom::from_aromatic_atom(Element::C))]
#[case("[C]", Atom::from_aliphatic_atom(Element::C))]
fn test_atom(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let mut state = ParseState::default();
    let parser = grammar::AtomParser::new();
    let result = parser.parse(&mut errors, &mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("[C]", Atom::from_aliphatic_atom(Element::C))]
#[case("[c]", Atom::from_aromatic_atom(Element::C))]
#[case("[13C]", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.isotope = Some(13); atom})]
#[case("[C@]", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.chirality = Some(Chirality::Clockwise); atom})]
#[case("[CH]", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.hydrogen_count = Some(1); atom})]
#[case("[C+]", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.charge = Some(1); atom})]
#[case("[C-]", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.charge = Some(-1); atom})]
#[case("[C++]", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.charge = Some(2); atom})]
#[case("[C--]", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.charge = Some(-2); atom})]
#[case("[C:1]", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.class = Some(1); atom})]
#[case("[13C@H+:1]", {let mut atom = Atom::from_aliphatic_atom(Element::C); atom.isotope = Some(13);
       atom.chirality = Some(Chirality::Clockwise); atom.hydrogen_count = Some(1); atom.charge = Some(1);
       atom.class = Some(1); atom})]
fn test_bracket_atom(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let mut state = ParseState::default();
    let parser = grammar::BracketAtomParser::new();
    let result = parser.parse(&mut errors, &mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("c", Atom::from_aromatic_atom(Element::C))]
#[case("*", Atom::default())]
fn test_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let mut state = ParseState::default();
    let parser = grammar::SymbolParser::new();
    let result = parser.parse(&mut errors, &mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("Ac", Atom::from_aliphatic_atom(Element::Ac))]
#[case("Ag", Atom::from_aliphatic_atom(Element::Ag))]
fn test_aliphatic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let mut state = ParseState::default();
    let parser = grammar::AliphaticSymbolParser::new();
    let result = parser.parse(&mut errors, &mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("p", Atom::from_aromatic_atom(Element::P))]
#[case("se", Atom::from_aromatic_atom(Element::Se))]
fn test_aromatic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let mut state = ParseState::default();
    let parser = grammar::AromaticSymbolParser::new();
    let result = parser.parse(&mut errors, &mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("*", Atom::default())]
fn test_unknown_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let mut state = ParseState::default();
    let parser = grammar::UnknownSymbolParser::new();
    let result = parser.parse(&mut errors, &mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("B", Atom::from_aliphatic_atom(Element::B))]
fn test_aliphatic_organic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let mut state = ParseState::default();
    let parser = grammar::AliphaticOrganicSymbolParser::new();
    let result = parser.parse(&mut errors, &mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("n", Atom::from_aromatic_atom(Element::N))]
#[case("o", Atom::from_aromatic_atom(Element::O))]
fn test_aromatic_organic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);

    let mut errors = Vec::new();
    let mut state = ParseState::default();
    let parser = grammar::AromaticOrganicSymbolParser::new();
    let result = parser.parse(&mut errors, &mut state, lexer).unwrap();
    assert_eq!(result, expected);
}
