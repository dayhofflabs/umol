use pretty_assertions::assert_eq;
use rstest::*;
use slog::o;
use umol::logging::setup_logger;
use umol::with_logger;
use umol_data::Element;

use crate::io::ir::{Atom, BondOrder, Chirality};
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::{branched, unbranched};
use crate::io::smiles::state::{BondSpec, ParseState};

#[fixture]
fn parse_state() -> ParseState {
    let mut state = ParseState::default();
    let root = setup_logger(slog::Level::Debug);
    state.log = Some(with_logger!(root, "io::smiles::parser"));
    state
}

#[rstest]
#[case("C")]
fn test_unbranched_atom(mut parse_state: ParseState, #[case] input: &str) {
    let lexer = Lexer::new(input);
    let parser = unbranched::AtomBondParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(parse_state.next_atom_idx, 1);
    assert_eq!(parse_state.next_bond_idx, 0);
    assert!(parse_state.staged_bond.is_none());
    assert!(parse_state.first_err.is_none());
    assert!(parse_state.drain_molecules().is_empty());
}

#[rstest]
#[case("=")]
fn test_unbranched_bond(mut parse_state: ParseState, #[case] input: &str) {
    let lexer = Lexer::new(input);
    parse_state.last_atom_idx = 0;
    parse_state.next_atom_idx = 1;
    let parser = unbranched::BondParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(
        parse_state.staged_bond,
        Some(BondSpec {
            order: BondOrder::Double,
            dir: None
        })
    );
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
fn test_unbranched(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
) {
    let lexer = Lexer::new(input);
    let parser = unbranched::ChainParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(parse_state.next_atom_idx, atoms);
    assert_eq!(parse_state.next_bond_idx, bonds);
    assert!(parse_state.staged_bond.is_none());
    assert!(parse_state.first_err.is_none());
    parse_state.finish_molecule();
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), atoms);
    assert_eq!(mols[0].bonds.len(), bonds);
}

#[rstest]
#[case("C/C", vec![Some(crate::io::ir::BondDir::Up)])]
#[case("C\\C", vec![Some(crate::io::ir::BondDir::Down)])]
#[case("C/C\\C", vec![Some(crate::io::ir::BondDir::Up), Some(crate::io::ir::BondDir::Down)])]
fn test_unbranched_bond_dirs(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] dirs: Vec<Option<crate::io::ir::BondDir>>,
) {
    let lexer = Lexer::new(input);
    let parser = unbranched::ChainParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    parse_state.finish_molecule();
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let got: Vec<Option<crate::io::ir::BondDir>> =
        mols[0].bonds.iter().map(|b| b.direction).collect();
    assert_eq!(got, dirs);
}

#[rstest]
#[case("CC.CC", 4, 2)]
fn test_unbranched_components(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
) {
    let lexer = Lexer::new(input);
    let parser = unbranched::ChainParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    parse_state.finish_molecule();
    let mols = parse_state.drain_molecules();
    let total_atoms: usize = mols.iter().map(|m| m.atoms.len()).sum();
    let total_bonds: usize = mols.iter().map(|m| m.bonds.len()).sum();
    assert_eq!(total_atoms, atoms);
    assert_eq!(total_bonds, bonds);
}

#[rstest]
#[case("[13C]C", 2, 1)]
#[case("C[O-]", 2, 1)]
#[case("[NH4+]", 1, 0)]
fn test_unbranched_bracket_atoms(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
) {
    let lexer = Lexer::new(input);
    let parser = unbranched::ChainParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    parse_state.finish_molecule();
    let mols = parse_state.drain_molecules();
    let total_atoms: usize = mols.iter().map(|m| m.atoms.len()).sum();
    let total_bonds: usize = mols.iter().map(|m| m.bonds.len()).sum();
    assert_eq!(total_atoms, atoms);
    assert_eq!(total_bonds, bonds);
}

#[rstest]
#[case("=C", "leading bond")]
#[case("C=", "trailing bond")]
#[case("C==C", "consecutive bonds")]
#[case(".CC", "leading dot")]
#[case("CC.", "trailing dot")]
#[case("C..C", "consecutive dots")]
#[case("[C", "unclosed bracket")]
#[case("C]", "unexpected close bracket")]
#[case("[C+", "unclosed bracket with charge")]
#[case("[13]", "missing symbol after isotope")]
fn test_unbranched_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = unbranched::ChainParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case("C", 1, 0)]
#[case("CC", 2, 1)]
#[case("C=C", 2, 1)]
#[case("CC(C)C", 4, 3)]
#[case("CC(C)", 3, 2)]
#[case("C(C)C", 3, 2)]
#[case("CC(=C)C", 4, 3)]
#[case("CC(C)(C)CC", 6, 5)]
#[case("CC(C(C)C)C", 6, 5)]
#[case("CC(CC)C", 5, 4)]
fn test_branched(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
) {
    let lexer = Lexer::new(input);
    let parser = branched::TreeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(parse_state.next_atom_idx, atoms);
    assert_eq!(parse_state.next_bond_idx, bonds);
    assert!(parse_state.staged_bond.is_none());
    assert!(parse_state.first_err.is_none());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), atoms);
    assert_eq!(mols[0].bonds.len(), bonds);
}

#[rstest]
#[case("C/C", vec![Some(crate::io::ir::BondDir::Up)])]
#[case("C\\C", vec![Some(crate::io::ir::BondDir::Down)])]
#[case("C(/C)C", vec![Some(crate::io::ir::BondDir::Up), None])]
fn test_branched_bond_dirs(#[case] input: &str, #[case] dirs: Vec<Option<crate::io::ir::BondDir>>) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = branched::TreeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let got: Vec<Option<crate::io::ir::BondDir>> =
        mols[0].bonds.iter().map(|b| b.direction).collect();
    assert_eq!(got, dirs);
}

#[rstest]
#[case("[13C]C", 2, 1)]
#[case("C[O-]", 2, 1)]
#[case("[NH4+]", 1, 0)]
#[case("[C@H](F)(Cl)Br", 4, 3)]
fn test_branched_bracket_atoms(#[case] input: &str, #[case] atoms: usize, #[case] bonds: usize) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = branched::TreeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), atoms);
    assert_eq!(mols[0].bonds.len(), bonds);
}

#[rstest]
#[case("C.CC", vec![(1, 0), (2, 1)])]
#[case("CC(CC.CC)C", vec![(5, 4), (2, 1)])]
#[case("CC(CC)C.CC", vec![(5, 4), (2, 1)])]
fn test_branched_components(#[case] input: &str, #[case] expected: Vec<(usize, usize)>) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = branched::TreeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert!(parse_state.staged_bond.is_none());
    assert!(parse_state.first_err.is_none());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), expected.len());
    assert_eq!(mols[0].atoms.len(), expected[0].0);
    assert_eq!(mols[0].bonds.len(), expected[0].1);
    assert_eq!(mols[1].atoms.len(), expected[1].0);
    assert_eq!(mols[1].bonds.len(), expected[1].1);
}

#[rstest]
#[case("=C", "leading bond")]
#[case("C=", "trailing bond")]
#[case("C==C", "consecutive bonds")]
#[case(".CC", "leading dot")]
#[case("CC.", "trailing dot")]
#[case("C..C", "consecutive dots")]
#[case("(C)CC", "leading branch with no parent")]
#[case("CC(C", "unclosed branch")]
#[case("C(C-)C", "trailing bond in branch")]
#[case("CC)C", "unmatched close paren")]
#[case("C()C", "empty branch")]
#[case("[C", "unclosed bracket")]
#[case("C]", "unexpected close bracket")]
#[case("[C+", "unclosed bracket with charge")]
#[case("[13]", "missing symbol after isotope")]
fn test_branched_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = branched::TreeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("c", Atom::from_aromatic_atom(Element::C))]
#[case("[C]", Atom::from_aliphatic_atom(Element::C))]
fn test_atom(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = branched::AtomParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
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
    let mut state = ParseState::default();
    let parser = branched::BracketAtomParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("[C", "unclosed bracket")]
#[case("C]", "unexpected close bracket")]
#[case("[C+", "unclosed bracket with charge")]
#[case("[13]", "missing symbol after isotope")]
fn test_bracket_atom_invalid(#[case] input: &str, #[case] _desc: &str) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = branched::BracketAtomParser::new();
    let result = parser.parse(&mut state, lexer);
    assert!(result.is_err());
}

#[rstest]
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("c", Atom::from_aromatic_atom(Element::C))]
#[case("*", Atom::default())]
fn test_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = branched::SymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("Ac", Atom::from_aliphatic_atom(Element::Ac))]
#[case("Ag", Atom::from_aliphatic_atom(Element::Ag))]
fn test_aliphatic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = branched::AliphaticSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("p", Atom::from_aromatic_atom(Element::P))]
#[case("se", Atom::from_aromatic_atom(Element::Se))]
fn test_aromatic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = branched::AromaticSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("*", Atom::default())]
fn test_unknown_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = branched::UnknownSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("B", Atom::from_aliphatic_atom(Element::B))]
fn test_aliphatic_organic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = branched::AliphaticOrganicSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("n", Atom::from_aromatic_atom(Element::N))]
#[case("o", Atom::from_aromatic_atom(Element::O))]
fn test_aromatic_organic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = branched::AromaticOrganicSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}
