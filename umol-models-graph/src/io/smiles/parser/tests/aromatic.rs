use rstest::*;

use super::fixtures::parse_state;
use crate::io::ir::{BondOrder, BondSymbol};
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::grammar::MoleculeParser;
use crate::io::smiles::state::ParseState;

#[rstest]
#[case("cc", 2, 1, Some(BondOrder::Aromatic))]
#[case("c:c", 2, 1, Some(BondOrder::Aromatic))]
#[case("cC", 2, 1, Some(BondOrder::Single))]
fn test_aromatic_core(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
    #[case] expected_order: Option<BondOrder>,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), atoms);
    assert_eq!(mols[0].bonds.len(), bonds);
    if let Some(o) = expected_order {
        assert_eq!(mols[0].bonds[0].symbol, BondSymbol::Bond(o));
    }
}

#[test]
fn test_aromatic_ring() {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new("c1ccccc1");
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), 6);
    assert_eq!(mols[0].bonds.len(), 6);
    assert_eq!(
        mols[0].bonds[0].symbol,
        BondSymbol::Bond(BondOrder::Aromatic)
    );
}

#[rstest]
#[case("c1ccc(cc1)N", 7)]
fn test_aromatic_mixed(#[case] input: &str, #[case] bonds: usize) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].bonds.len(), bonds);
    assert!(mols[0]
        .bonds
        .iter()
        .any(|b| b.symbol == BondSymbol::Bond(BondOrder::Aromatic)));
}

#[rstest]
#[case("C:C")]
fn test_aromatic_explicit_colon(#[case] input: &str) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "Parser should accept, semantic policy TBD");
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(
        mols[0].bonds[0].symbol,
        BondSymbol::Bond(BondOrder::Aromatic)
    );
}

#[rstest]
#[case("C*C", 3, 2)]
fn test_unknown_symbol(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), atoms);
    assert_eq!(mols[0].bonds.len(), bonds);
}
