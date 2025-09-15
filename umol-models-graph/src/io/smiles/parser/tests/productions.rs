use pretty_assertions::assert_eq;
use rstest::*;
use umol_data::Element;

use super::fixtures::parse_state;
use crate::io::ir::{Atom, BondDir, BondOrder, Chirality};
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::grammar::{
    AliphaticOrganicSymbolParser, AliphaticSymbolParser, AromaticOrganicSymbolParser,
    AromaticSymbolParser, AtomBondParser, AtomParser, BondOrderParser, BondParser,
    BracketAtomParser, ChargeParser, ChiralityParser, ClassParser, HCountParser, IsotopeParser,
    RingIndexParser, RingSpecParser, SymbolParser, UnknownSymbolParser,
};
use crate::io::smiles::state::{BondInfo, ParseState};

#[rstest]
#[case("C", { let mut a = Atom::from_aliphatic_atom(Element::C); a.implicit_h = true; a })]
#[case("c", { let mut a = Atom::from_aromatic_atom(Element::C); a.implicit_h = true; a })]
#[case("[C]", { let mut a = Atom::from_aliphatic_atom(Element::C); a.class = Some(0); a.hydrogen_count = Some(0); a.implicit_h = false; a })]
fn test_atom(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = AtomParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("[C]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.class = Some(0); atom.hydrogen_count = Some(0); atom.implicit_h = false; atom })]
#[case("[c]", { let mut atom = Atom::from_aromatic_atom(Element::C); atom.class = Some(0); atom.hydrogen_count = Some(0); atom.implicit_h = false; atom })]
#[case("[*]", { let mut atom = Atom::default(); atom.class = Some(0); atom.hydrogen_count = Some(0); atom.implicit_h = false; atom })]
#[case("[*H2+:1]", { let mut atom = Atom::default(); atom.hydrogen_count = Some(2); atom.charge = Some(1); atom.class = Some(1); atom.implicit_h = false; atom })]
#[case("[C+:1H]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.charge = Some(1); atom.class = Some(1); atom.hydrogen_count = Some(1); atom.implicit_h = false; atom })]
#[case("[C@H+]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.chirality = Some(Chirality::Clockwise); atom.hydrogen_count = Some(1); atom.charge = Some(1); atom.class = Some(0); atom.implicit_h = false; atom })]
#[case("[13C]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.isotope = Some(13); atom.class = Some(0); atom.hydrogen_count = Some(0); atom.implicit_h = false; atom })]
#[case("[C@]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.chirality = Some(Chirality::Clockwise); atom.class = Some(0); atom.hydrogen_count = Some(0); atom.implicit_h = false; atom })]
#[case("[CH]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.hydrogen_count = Some(1); atom.class = Some(0); atom.implicit_h = false; atom })]
#[case("[C+]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.charge = Some(1); atom.class = Some(0); atom.hydrogen_count = Some(0); atom.implicit_h = false; atom })]
#[case("[C-]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.charge = Some(-1); atom.class = Some(0); atom.hydrogen_count = Some(0); atom.implicit_h = false; atom })]
#[case("[C++]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.charge = Some(2); atom.class = Some(0); atom.hydrogen_count = Some(0); atom.implicit_h = false; atom })]
#[case("[C--]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.charge = Some(-2); atom.class = Some(0); atom.hydrogen_count = Some(0); atom.implicit_h = false; atom })]
#[case("[C:1]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.class = Some(1); atom.hydrogen_count = Some(0); atom.implicit_h = false; atom })]
#[case("[13C@H+:1]", { let mut atom = Atom::from_aliphatic_atom(Element::C); atom.isotope = Some(13);
       atom.chirality = Some(Chirality::Clockwise); atom.hydrogen_count = Some(1); atom.charge = Some(1);
       atom.class = Some(1); atom.implicit_h = false; atom })]
fn test_bracket_atom(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = BracketAtomParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("C", { let mut a = Atom::from_aliphatic_atom(Element::C); a.implicit_h = false; a })]
#[case("c", { let mut a = Atom::from_aromatic_atom(Element::C); a.implicit_h = false; a })]
#[case("*", { let mut a = Atom::default(); a.implicit_h = false; a })]
fn test_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = SymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("Ac", { let mut a = Atom::from_aliphatic_atom(Element::Ac); a.implicit_h = false; a })]
#[case("Ag", { let mut a = Atom::from_aliphatic_atom(Element::Ag); a.implicit_h = false; a })]
fn test_aliphatic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = AliphaticSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("p", { let mut a = Atom::from_aromatic_atom(Element::P); a.implicit_h = false; a })]
#[case("se", { let mut a = Atom::from_aromatic_atom(Element::Se); a.implicit_h = false; a })]
fn test_aromatic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = AromaticSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("*", { let mut a = Atom::default(); a.implicit_h = false; a })]
fn test_unknown_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = UnknownSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("C", { let mut a = Atom::from_aliphatic_atom(Element::C); a.implicit_h = true; a })]
#[case("B", { let mut a = Atom::from_aliphatic_atom(Element::B); a.implicit_h = true; a })]
fn test_aliphatic_organic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = AliphaticOrganicSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("n", { let mut a = Atom::from_aromatic_atom(Element::N); a.implicit_h = true; a })]
#[case("o", { let mut a = Atom::from_aromatic_atom(Element::O); a.implicit_h = true; a })]
fn test_aromatic_organic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = AromaticOrganicSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("C")]
fn test_atombond(mut parse_state: ParseState, #[case] input: &str) {
    let lexer = Lexer::new(input);
    let parser = AtomBondParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(parse_state.next_atom_idx, 1);
    assert_eq!(parse_state.next_bond_idx, 0);
    assert!(parse_state.staged_bond.is_none());
    assert!(parse_state.first_err.is_none());
    assert!(parse_state.drain_molecules().is_empty());
}

#[rstest]
#[case("13", 13u32)]
#[case("0", 0u32)]
fn test_isotope(#[case] input: &str, #[case] expected: u32) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = IsotopeParser::new();
    let got = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(got, expected);
    assert!(state.drain_molecules().is_empty());
}

#[rstest]
#[case("@", Chirality::Clockwise)]
#[case("@@", Chirality::CounterClockwise)]
#[case("@TH1", Chirality::Clockwise)]
#[case("@TH2", Chirality::CounterClockwise)]
#[case("@AL1", Chirality::Allenal { arr: 1 })]
#[case("@SP1", Chirality::SquarePlanar { arr: 1 })]
#[case("@TB10", Chirality::TrigonalBipyramidal { arr: 10 })]
#[case("@OH10", Chirality::Octahedral { arr: 10 })]
fn test_chirality(#[case] input: &str, #[case] expected: Chirality) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = ChiralityParser::new();
    let got = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(got, expected);
}

#[rstest]
#[case("H", 1u32)]
#[case("H2", 2u32)]
fn test_hcount(#[case] input: &str, #[case] expected: u32) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = HCountParser::new();
    let got = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(got, expected);
}

#[rstest]
#[case("+", 1i32)]
#[case("-", -1i32)]
#[case("++", 2i32)]
#[case("--", -2i32)]
#[case("+3", 3i32)]
#[case("-4", -4i32)]
#[case("+12", 12i32)]
#[case("-12", -12i32)]
fn test_charge(#[case] input: &str, #[case] expected: i32) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = ChargeParser::new();
    let got = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(got, expected);
}

#[rstest]
#[case(":0", 0u32)]
#[case(":12", 12u32)]
#[case(":120", 120u32)]
fn test_class(#[case] input: &str, #[case] expected: u32) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = ClassParser::new();
    let got = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(got, expected);
}

#[rstest]
#[case("X", "invalid element symbol")]
#[case("f", "invalid aromatic organic symbol")]
fn test_atom_invalid(#[case] input: &str, #[case] desc: &str) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = ClassParser::new();
    let result = parser.parse(&mut state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case("[C", "unclosed bracket")]
#[case("C]", "unexpected close bracket")]
#[case("[C+", "unclosed bracket with charge")]
#[case("[13]", "missing symbol after isotope")]
fn test_bracket_atom_invalid(#[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = BracketAtomParser::new();
    let result = parser.parse(&mut state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case("@", "not an HCount token")]
#[case(":", "not an HCount token")]
fn test_hcount_invalid(#[case] input: &str, #[case] desc: &str) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = HCountParser::new();
    let result = parser.parse(&mut state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case(":", "missing number after class colon")]
fn test_class_invalid(#[case] input: &str, #[case] desc: &str) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = ClassParser::new();
    let result = parser.parse(&mut state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case("-", BondOrder::Single, None)]
#[case("=", BondOrder::Double, None)]
#[case("#", BondOrder::Triple, None)]
#[case("$", BondOrder::Quadruple, None)]
#[case(":", BondOrder::Aromatic, None)]
#[case("/", BondOrder::Single, Some(BondDir::Up))]
#[case("\\", BondOrder::Single, Some(BondDir::Down))]
fn test_bond(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] expected_order: BondOrder,
    #[case] expected_dir: Option<BondDir>,
) {
    let lexer = Lexer::new(input);
    parse_state.next_atom_idx = 1;
    let parser = BondParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(
        parse_state.staged_bond,
        Some(BondInfo {
            order: expected_order,
            dir: expected_dir
        })
    );
    assert!(parse_state.drain_molecules().is_empty());
}

#[rstest]
#[case("-", BondOrder::Single)]
#[case("=", BondOrder::Double)]
#[case("#", BondOrder::Triple)]
#[case("$", BondOrder::Quadruple)]
fn test_bondorder(#[case] input: &str, #[case] expected: BondOrder) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = BondOrderParser::new();
    let got = parser.parse(&mut parse_state, lexer).unwrap();
    assert_eq!(got, expected);
    assert!(parse_state.drain_molecules().is_empty());
}

#[rstest]
#[case("!", "unsupported bond symbol")]
#[case("??", "unsupported token sequence")]
fn test_bond_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = BondParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(parse_state.staged_bond.is_none());
    assert!(parse_state.drain_molecules().is_empty());
}

#[rstest]
#[case("-", BondInfo { order: BondOrder::Single, dir: None })]
#[case("=", BondInfo { order: BondOrder::Double, dir: None })]
#[case(":", BondInfo { order: BondOrder::Aromatic, dir: None })]
#[case("/", BondInfo { order: BondOrder::Single, dir: Some(crate::io::ir::BondDir::Up) })]
#[case("\\", BondInfo { order: BondOrder::Single, dir: Some(crate::io::ir::BondDir::Down) })]
fn test_ringspec(#[case] input: &str, #[case] expected: BondInfo) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = RingSpecParser::new();
    let got = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(got, expected);
}

#[rstest]
#[case("1", 1u32)]
#[case("9", 9u32)]
#[case("%12", 12u32)]
fn test_ringindex(#[case] input: &str, #[case] expected: u32) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = RingIndexParser::new();
    let got = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(got, expected);
}

#[rstest]
#[case("%0", "percent ring index must be two digits and non-zero-leading")]
fn test_ringindex_invalid(#[case] input: &str, #[case] desc: &str) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = RingIndexParser::new();
    let result = parser.parse(&mut state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}
