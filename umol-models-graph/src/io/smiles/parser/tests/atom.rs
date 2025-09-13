use rstest::*;

use super::fixtures::parse_state;
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::grammar::MoleculeParser;
use crate::io::smiles::state::ParseState;

#[rstest]
#[case("B", 1, 0)]
#[case("C", 1, 0)]
#[case("N", 1, 0)]
#[case("O", 1, 0)]
#[case("S", 1, 0)]
#[case("P", 1, 0)]
#[case("F", 1, 0)]
#[case("Cl", 1, 0)]
#[case("Br", 1, 0)]
#[case("I", 1, 0)]
#[case("b", 1, 0)]
#[case("c", 1, 0)]
#[case("n", 1, 0)]
#[case("o", 1, 0)]
#[case("s", 1, 0)]
#[case("p", 1, 0)]
fn test_organic_atoms(
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
    let total_atoms: usize = mols.iter().map(|m| m.atoms.len()).sum();
    let total_bonds: usize = mols.iter().map(|m| m.bonds.len()).sum();
    assert_eq!(total_atoms, atoms);
    assert_eq!(total_bonds, bonds);
}

#[rstest]
#[case("[*]", 1, 0)]
fn test_unknown_atom(
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
    let total_atoms: usize = mols.iter().map(|m| m.atoms.len()).sum();
    let total_bonds: usize = mols.iter().map(|m| m.bonds.len()).sum();
    assert_eq!(total_atoms, atoms);
    assert_eq!(total_bonds, bonds);
}

#[rstest]
#[case("[13C]C", 2, 1)]
#[case("C[O-]", 2, 1)]
#[case("[NH4+]", 1, 0)]
#[case("[C@H]", 1, 0)]
#[case("[C+]", 1, 0)]
#[case("[C-]", 1, 0)]
#[case("[C++]", 1, 0)]
#[case("[C--]", 1, 0)]
#[case("[C:12]", 1, 0)]
fn test_atom_bracket(
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
    let total_atoms: usize = mols.iter().map(|m| m.atoms.len()).sum();
    let total_bonds: usize = mols.iter().map(|m| m.bonds.len()).sum();
    assert_eq!(total_atoms, atoms);
    assert_eq!(total_bonds, bonds);
}

#[rstest]
#[case("X", "unknown element symbol")] 
#[case("f", "invalid aromatic organic symbol")]
#[case("[C", "unclosed bracket")]
#[case("C]", "unexpected close bracket")]
#[case("[C+", "unclosed bracket with charge")]
#[case("[13]", "missing symbol after isotope")]
fn test_atom_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}
