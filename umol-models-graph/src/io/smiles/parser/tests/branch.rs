use rstest::*;

use super::fixtures::parse_state;
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::grammar::MoleculeParser;
use crate::io::smiles::state::ParseState;

#[rstest]
#[case("CC(C)C", 4, 3)]
#[case("CC(C)", 3, 2)]
#[case("C(C)C", 3, 2)]
#[case("CC(=C)C", 4, 3)]
#[case("CC(C)(C)CC", 6, 5)]
#[case("CC(C(C)C)C", 6, 5)]
#[case("CC(CC)C", 5, 4)]
fn test_branch(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
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
#[case("C(.C)C", vec![(2, 1), (1, 0)])]
#[case("C(C.C)C", vec![(3, 2), (1, 0)])]
fn test_branch_components(
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
    assert_eq!(mols[0].atoms.len(), expected[0].0);
    assert_eq!(mols[0].bonds.len(), expected[0].1);
    assert_eq!(mols[1].atoms.len(), expected[1].0);
    assert_eq!(mols[1].bonds.len(), expected[1].1);
}

#[rstest]
#[case("(C)CC", "leading branch with no parent")]
#[case("CC(C", "unclosed branch")]
#[case("C(C-)C", "trailing bond in branch")]
#[case("CC)C", "unmatched close paren")]
#[case("C()C", "empty branch")]
fn test_branch_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}
