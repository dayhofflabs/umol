use rstest::*;

use super::fixtures::parse_state;
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::grammar::MoleculeParser;
use crate::io::smiles::state::ParseState;

#[rstest]
#[case("CC", 1)]
#[case("C=C", 1)]
#[case("C#C", 1)]
#[case("C$C", 1)]
#[case("c:c", 1)]
fn test_bond(mut parse_state: ParseState, #[case] input: &str, #[case] bond_count: usize) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].bonds.len(), bond_count);
}

#[rstest]
#[case("C/C", vec![Some(crate::io::ir::BondDir::Up)])]
#[case("C\\C", vec![Some(crate::io::ir::BondDir::Down)])]
#[case("C/C\\C", vec![Some(crate::io::ir::BondDir::Up), Some(crate::io::ir::BondDir::Down)])]
fn test_bond_dirs(#[case] input: &str, #[case] dirs: Vec<Option<crate::io::ir::BondDir>>) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let got: Vec<Option<crate::io::ir::BondDir>> =
        mols[0].bonds.iter().map(|b| b.direction).collect();
    assert_eq!(got, dirs);
}

#[rstest]
#[case("C1CC1", 3)]
#[case("c1ccccc1", 6)]
fn test_bond_ring(mut parse_state: ParseState, #[case] input: &str, #[case] bonds: usize) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].bonds.len(), bonds);
}

#[rstest]
#[case("cc", 1)]
#[case("c:c", 1)]
fn test_bond_aromatic(mut parse_state: ParseState, #[case] input: &str, #[case] bonds: usize) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].bonds.len(), bonds);
}

#[rstest]
#[case("=C", "leading bond")]
#[case("C=", "trailing bond")]
#[case("C==C", "consecutive bonds")]
#[case("C//C", "consecutive direction markers")]
#[case("C/\\C", "conflicting adjacent directions")]
fn test_bond_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}
