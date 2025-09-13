use rstest::*;

use super::fixtures::parse_state;
use crate::io::ir::{BondOrder, BondStereo, BondSymbol};
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::grammar::MoleculeParser;
use crate::io::smiles::state::ParseState;
#[rstest]
#[case("C", 1, 0)]
#[case("CC", 2, 1)]
#[case("CCC", 3, 2)]
#[case("C=C", 2, 1)]
fn test_chain(
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

#[rstest]
#[case("F/C=C/F", Some(BondStereo::Trans))]
#[case("F/C=C\\F", Some(BondStereo::Cis))]
#[case("F\\C=C\\F", Some(BondStereo::Trans))]
#[case("F\\C=C/F", Some(BondStereo::Cis))]
#[case("FC=C/F", Some(BondStereo::Either))]
#[case("F/C=CF", Some(BondStereo::Either))]
#[case("FC=CF", None)]
fn test_chain_ez(
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
#[case("C.C", vec![(1, 0), (1, 0)])]
fn test_chain_components(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] expected: Vec<(usize, usize)>,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), expected.len());
    for (i, (ea, eb)) in expected.iter().enumerate() {
        assert_eq!(mols[i].atoms.len(), *ea);
        assert_eq!(mols[i].bonds.len(), *eb);
    }
}

#[rstest]
#[case("C")]
#[case("C ")]
#[case("C\t")]
#[case("C\n")]
#[case("C\r")]
fn test_chain_terminator(mut parse_state: ParseState, #[case] input: &str) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), 1);
    assert_eq!(mols[0].bonds.len(), 0);
}

#[rstest]
#[case("")]
#[case(" ")]
#[case("\t")]
#[case("\n")]
#[case("\r")]
#[case("   ")]
fn test_chain_terminator_only(mut parse_state: ParseState, #[case] input: &str) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "terminator-only inputs should succeed");
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 0);
}

#[rstest]
#[case(".CC", "leading dot")]
#[case("CC.", "trailing dot")]
#[case("C..C", "consecutive dots")]
fn test_chain_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}
