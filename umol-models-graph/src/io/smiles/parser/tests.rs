use pretty_assertions::assert_eq;
use rstest::*;
use umol_data::Element;

use crate::io::ir::{Atom, Chirality};
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::grammar::{
    AliphaticOrganicSymbolParser, AliphaticSymbolParser, AromaticOrganicSymbolParser,
    AromaticSymbolParser, AtomParser, BracketAtomParser, SymbolParser, UnknownSymbolParser,
};
use crate::io::smiles::state::ParseState;

#[rstest]
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("c", Atom::from_aromatic_atom(Element::C))]
#[case("[C]", Atom::from_aliphatic_atom(Element::C))]
fn test_atom(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let mut state = ParseState::default();
    let parser = AtomParser::new();
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
    let parser = BracketAtomParser::new();
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
    let parser = SymbolParser::new();
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
    let parser = AliphaticSymbolParser::new();
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
    let parser = AromaticSymbolParser::new();
    let result = parser.parse(&mut errors, &mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("*", Atom::default())]
fn test_unknown_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut errors = Vec::new();
    let mut state = ParseState::default();
    let parser = UnknownSymbolParser::new();
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
    let parser = AliphaticOrganicSymbolParser::new();
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
    let parser = AromaticOrganicSymbolParser::new();
    let result = parser.parse(&mut errors, &mut state, lexer).unwrap();
    assert_eq!(result, expected);
}
