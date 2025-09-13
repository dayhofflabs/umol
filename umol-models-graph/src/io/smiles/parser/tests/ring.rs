use rstest::*;

use super::fixtures::parse_state;
use crate::io::ir::{BondOrder, BondStereo, BondSymbol};
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::grammar::MoleculeParser;
use crate::io::smiles::state::ParseState;

#[rstest]
#[case("C1CC1", 3, 3)]
#[case("C12CCC1C2", 5, 6)]
#[case("C%12CC%12", 3, 3)]
fn test_ring(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), atoms);
    assert_eq!(mols[0].bonds.len(), bonds);
}

#[test]
fn test_ring_percent() {
    let mut parse_state = ParseState::default();
    let input = "C%12CC\\%12";
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert!(mols[0].bonds.iter().any(|b| b.direction.is_some()));
}

#[rstest]
#[case("[13C]1CC1", 3, 3)]
#[case("C1[C-]C1", 3, 3)]
fn test_ring_with_bracket(#[case] input: &str, #[case] atoms: usize, #[case] bonds: usize) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), atoms);
    assert_eq!(mols[0].bonds.len(), bonds);
}

#[rstest]
#[case("C1CC/C=C/C1", Some(BondStereo::Trans))]
#[case("C1CC/C=C\\C1", Some(BondStereo::Cis))]
fn test_ring_stereo(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] expected: Option<BondStereo>,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let dbl = mols[0]
        .bonds
        .iter()
        .find(|b| matches!(b.symbol, BondSymbol::Bond(BondOrder::Double)))
        .and_then(|b| b.stereo);
    assert_eq!(dbl, expected);
}

#[rstest]
#[case("F/C=C\\1CC1", Some(BondStereo::Either))]
#[case("F/C=C/1CC1", Some(BondStereo::Either))]
fn test_ring_directed_closure(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] expected: Option<BondStereo>,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let dbl = mols[0]
        .bonds
        .iter()
        .find(|b| matches!(b.symbol, BondSymbol::Bond(BondOrder::Double)))
        .and_then(|b| b.stereo);
    assert_eq!(dbl, expected);
}

#[rstest]
#[case("c1ccccc1", 6)]
#[case("c1:c:c:c:c:c1", 6)]
fn test_ring_aromatic(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] aromatic_bonds: usize,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let aromatic_count = mols[0]
        .bonds
        .iter()
        .filter(|b| b.symbol == BondSymbol::Bond(BondOrder::Aromatic))
        .count();
    assert_eq!(aromatic_count, aromatic_bonds);
}

#[test]
fn test_ring_components_merge() {
    let mut parse_state = ParseState::default();
    let input = "C1.C12.C2";
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), 3);
    assert_eq!(mols[0].bonds.len(), 2);
}

#[rstest]
#[case("C1CC", "unclosed ring index")] // no closing '1'
#[case("C11", "self-loop ring closure")] // creates a bond from atom to itself
#[case("C12C21", "two-member ring not allowed")] // forms 2-cycle
#[case("C/1CC\\1", "conflicting directions on ring bond")] // slash up at open, backslash down at close
fn test_ring_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}
