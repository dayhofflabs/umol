use pretty_assertions::assert_eq;
use rstest::*;
use umol_data::Element;

use super::*;
use crate::io::ir::Atom;
use crate::io::smiles::lexer::Lexer;

#[rstest]
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("B", Atom::from_aliphatic_atom(Element::B))]
fn test_aliphatic_organic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let parser = grammar::AliphaticOrganicSymbolParser::new();
    let result = parser.parse(&mut errors, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("n", Atom::from_aromatic_atom(Element::N))]
#[case("o", Atom::from_aromatic_atom(Element::O))]
fn test_aromatic_organic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let parser = grammar::AromaticOrganicSymbolParser::new();
    let result = parser.parse(&mut errors, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("Ac", Atom::from_aliphatic_atom(Element::Ac))]
#[case("Ag", Atom::from_aliphatic_atom(Element::Ag))]
fn test_aliphatic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let parser = grammar::AliphaticSymbolParser::new();
    let result = parser.parse(&mut errors, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("p", Atom::from_aromatic_atom(Element::P))]
#[case("se", Atom::from_aromatic_atom(Element::Se))]
fn test_aromatic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let parser = grammar::AromaticSymbolParser::new();
    let result = parser.parse(&mut errors, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("*", Atom::default())]
fn test_unknown_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let parser = grammar::UnknownSymbolParser::new();
    let result = parser.parse(&mut errors, lexer).unwrap();
    assert_eq!(result, expected);
}