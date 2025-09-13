use pretty_assertions::assert_eq;
use rstest::*;
use slog::o;
use umol::logging::setup_logger;
use umol::with_logger;
use umol_data::Element;

use crate::io::ctab::bond::BondStereo;
use crate::io::ir::{Atom, BondOrder, BondSymbol, Chirality};
use crate::io::smiles::lexer::Lexer;
use crate::io::smiles::parser::grammar::{
    AliphaticOrganicSymbolParser, AliphaticSymbolParser, AromaticOrganicSymbolParser,
    AromaticSymbolParser, AtomBondParser, AtomParser, BondOrderParser, BondParser,
    BracketAtomParser, ChargeParser, ChiralityParser, ClassParser, HCountParser,
    IsotopeParser, MoleculeParser, RingIndexParser, RingSpecParser, SymbolParser,
    UnknownSymbolParser,
};
use crate::io::smiles::state::{BondInfo, ParseState};

#[fixture]
fn parse_state() -> ParseState {
    let mut state = ParseState::default();
    let root = setup_logger(slog::Level::Debug);
    state.log = Some(with_logger!(root, "io::smiles::parser"));
    state
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
#[case("-", BondOrder::Single, None)]
#[case("=", BondOrder::Double, None)]
#[case("#", BondOrder::Triple, None)]
#[case("$", BondOrder::Quadruple, None)]
#[case(":", BondOrder::Aromatic, None)]
#[case("/", BondOrder::Single, Some(crate::io::ir::BondDir::Up))]
#[case("\\", BondOrder::Single, Some(crate::io::ir::BondDir::Down))]
fn test_bond(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] expected_order: BondOrder,
    #[case] expected_dir: Option<crate::io::ir::BondDir>,
) {
    let lexer = Lexer::new(input);
    // Allow staging a bond order without an existing chain
    parse_state.next_atom_idx = 1;
    let parser = BondParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    assert_eq!(
        parse_state.staged_bond,
        Some(BondInfo {
            order: expected_order,
            dir: expected_dir,
        })
    );
    assert!(parse_state.drain_molecules().is_empty());
}

#[rstest]
#[case("!", "unsupported bond glyph")]
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

// Additional microparsers in grammar order
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
fn test_class(#[case] input: &str, #[case] expected: u32) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = ClassParser::new();
    let got = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(got, expected);
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

// Invalids for selected microparsers
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
#[case("%0", "percent ring index must be two digits and non-zero-leading")]
fn test_ringindex_invalid(#[case] input: &str, #[case] desc: &str) {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = RingIndexParser::new();
    let result = parser.parse(&mut state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}
#[rstest]
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("c", Atom::from_aromatic_atom(Element::C))]
#[case("[C]", Atom::from_aliphatic_atom(Element::C))]
fn test_atom(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = AtomParser::new();
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
    let parser = BracketAtomParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
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
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("c", Atom::from_aromatic_atom(Element::C))]
#[case("*", Atom::default())]
fn test_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = SymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("Ac", Atom::from_aliphatic_atom(Element::Ac))]
#[case("Ag", Atom::from_aliphatic_atom(Element::Ag))]
fn test_aliphatic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = AliphaticSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("p", Atom::from_aromatic_atom(Element::P))]
#[case("se", Atom::from_aromatic_atom(Element::Se))]
fn test_aromatic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = AromaticSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("*", Atom::default())]
fn test_unknown_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = UnknownSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("C", Atom::from_aliphatic_atom(Element::C))]
#[case("B", Atom::from_aliphatic_atom(Element::B))]
fn test_aliphatic_organic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = AliphaticOrganicSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

#[rstest]
#[case("n", Atom::from_aromatic_atom(Element::N))]
#[case("o", Atom::from_aromatic_atom(Element::O))]
fn test_aromatic_organic_symbol(#[case] input: &str, #[case] expected: Atom) {
    let lexer = Lexer::new(input);
    let mut state = ParseState::default();
    let parser = AromaticOrganicSymbolParser::new();
    let result = parser.parse(&mut state, lexer).unwrap();
    assert_eq!(result, expected);
}

// MoleculeParser tests — single-feature (valid)
#[rstest]
#[case("C1CC1", 3, 3)]
#[case("C12CCC1C2", 5, 6)]
#[case("C%12CC%12", 3, 3)]
#[case("C$C", 2, 1)]
fn test_molecule_rings_and_special_bonds_valid(
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
#[case("C/1CC1", 1, 0)]
#[case("C1CC\\1", 0, 1)]
#[case("c:1cccc1", 0, 0)]
fn test_molecule_ring_spec_valid_linear(
    #[case] input: &str,
    #[case] up_cnt: usize,
    #[case] down_cnt: usize,
) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let dirs: Vec<_> = mols[0]
        .bonds
        .iter()
        .filter_map(|b| b.direction)
        .collect();
    let got_up = dirs.iter().filter(|d| **d == crate::io::ir::BondDir::Up).count();
    let got_down = dirs.iter().filter(|d| **d == crate::io::ir::BondDir::Down).count();
    if input != "c:1cccc1" {
        assert_eq!(got_up, up_cnt);
        assert_eq!(got_down, down_cnt);
    } else {
        // At least one aromatic bond present (colon ring spec)
        assert!(mols[0]
            .bonds
            .iter()
            .any(|b| b.symbol == BondSymbol::Bond(BondOrder::Aromatic)));
    }
}

#[rstest]
#[case("C/C", vec![Some(crate::io::ir::BondDir::Up)])]
#[case("C\\C", vec![Some(crate::io::ir::BondDir::Down)])]
#[case("C/C\\C", vec![Some(crate::io::ir::BondDir::Up), Some(crate::io::ir::BondDir::Down)])]
fn test_molecule_bond_dirs_valid_linear(#[case] input: &str, #[case] dirs: Vec<Option<crate::io::ir::BondDir>>) {
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
#[case("F/C=C/F", Some(BondStereo::Trans))]
#[case("F/C=C\\F", Some(BondStereo::Cis))]
#[case("F\\C=C\\F", Some(BondStereo::Trans))]
#[case("F\\C=C/F", Some(BondStereo::Cis))]
#[case("FC=C/F", Some(BondStereo::Either))]
#[case("F/C=CF", Some(BondStereo::Either))]
#[case("FC=CF", None)]
fn test_molecule_ez_valid_linear(#[case] input: &str, #[case] expected: Option<BondStereo>) {
    let mut parse_state = ParseState::default();
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

// Aromatic core tests
#[rstest]
#[case("cc", 2, 1, Some(BondOrder::Aromatic))]
#[case("c:c", 2, 1, Some(BondOrder::Aromatic))]
#[case("cC", 2, 1, Some(BondOrder::Single))]
fn test_molecule_aromatic_core_valid(
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
fn test_molecule_aromatic_ring_valid() {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new("c1ccccc1");
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), 6);
    assert_eq!(mols[0].bonds.len(), 6);
    // At least first bond aromatic by default (others may be as well)
    assert_eq!(
        mols[0].bonds[0].symbol,
        BondSymbol::Bond(BondOrder::Aromatic)
    );
}

#[rstest]
#[case("c1ccc(cc1)N", 7)]
fn test_molecule_aromatic_mixed_valid(#[case] input: &str, #[case] bonds: usize) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].bonds.len(), bonds);
    // Aromatic ring bonds should be aromatic
    assert!(mols[0]
        .bonds
        .iter()
        .any(|b| b.symbol == BondSymbol::Bond(BondOrder::Aromatic)));
}

#[rstest]
#[case("C:C")] // explicit colon between aliphatic atoms
fn test_molecule_aromatic_explicit_colon_valid(#[case] input: &str) {
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

// dropped as redundant with branched component tests

#[rstest]
#[case("[13C]C", 2, 1)]
#[case("C[O-]", 2, 1)]
#[case("[NH4+]", 1, 0)]
fn test_molecule_bracket_atoms_valid_linear(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
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
// MoleculeParser tests — single-feature (invalid)
fn test_molecule_linear_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
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
// MoleculeParser tests — core branching cases (valid)
fn test_molecule_branched_valid(
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

#[test]
fn test_molecule_branch_dot_components_valid() {
    // OpenBranchDot: C(.C)C -> two molecules: main CC and isolated C
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new("C(.C)C");
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 2);
    let sizes: Vec<(usize, usize)> = mols.iter().map(|m| (m.atoms.len(), m.bonds.len())).collect();
    assert!(sizes.contains(&(2, 1)));
    assert!(sizes.contains(&(1, 0)));
}

#[test]
fn test_molecule_branch_dot_tail_components_valid() {
    // BranchDotTail: C(C.C)C -> main 3 atoms/2 bonds + isolated C
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new("C(C.C)C");
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 2);
    let sizes: Vec<(usize, usize)> = mols.iter().map(|m| (m.atoms.len(), m.bonds.len())).collect();
    assert!(sizes.contains(&(3, 2)));
    assert!(sizes.contains(&(1, 0)));
}

#[test]
fn test_molecule_unknown_symbol_valid() {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new("C*C");
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms.len(), 3);
    assert_eq!(mols[0].bonds.len(), 2);
}

#[rstest]
#[case("C/C", vec![Some(crate::io::ir::BondDir::Up)])]
#[case("C\\C", vec![Some(crate::io::ir::BondDir::Down)])]
#[case("C(/C)C", vec![Some(crate::io::ir::BondDir::Up), None])]
fn test_molecule_bond_dirs_valid_branched(#[case] input: &str, #[case] dirs: Vec<Option<crate::io::ir::BondDir>>) {
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
#[case("[13C]C", 2, 1)]
#[case("C[O-]", 2, 1)]
#[case("[NH4+]", 1, 0)]
#[case("[C@H](F)(Cl)Br", 4, 3)]
fn test_molecule_bracket_atoms_valid_branched(#[case] input: &str, #[case] atoms: usize, #[case] bonds: usize) {
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
#[case("C.CC", vec![(1, 0), (2, 1)])]
#[case("CC(CC.CC)C", vec![(5, 4), (2, 1)])]
#[case("CC(CC)C.CC", vec![(5, 4), (2, 1)])]
fn test_molecule_components_valid_branched(#[case] input: &str, #[case] expected: Vec<(usize, usize)>) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
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
#[case("[C@H](F)(Cl)Br", Chirality::Clockwise)]
#[case("[C@@H](F)(Cl)Br", Chirality::CounterClockwise)]
#[case("[C@TH1H](F)(Cl)Br", Chirality::Clockwise)]
#[case("[C@TH2H](F)(Cl)Br", Chirality::CounterClockwise)]
fn test_molecule_tetra_chirality_valid(#[case] input: &str, #[case] expected: Chirality) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    // The first atom is the chiral carbon in these cases
    assert_eq!(mols[0].atoms[0].chirality, Some(expected));
}

#[rstest]
#[case("[C@]([H])(F)(Cl)Br", Chirality::Clockwise)]
#[case("[C@@]([H])(F)(Cl)Br", Chirality::CounterClockwise)]
#[case("[C@TH1]([H])(F)(Cl)Br", Chirality::Clockwise)]
#[case("[C@TH2]([H])(F)(Cl)Br", Chirality::CounterClockwise)]
fn test_molecule_tetra_chirality_explicit_h_valid(#[case] input: &str, #[case] expected: Chirality) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms[0].chirality, Some(expected));
}

#[test]
// MoleculeParser tests — single-feature tetrahedral edge case
fn test_molecule_tetra_chirality_insufficient_neighbors_edge() {
    let mut parse_state = ParseState::default();
    let input = "[C@](F)(Cl)C"; // three explicit neighbors and no bracket H
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(
        result.is_ok(),
        "{} should have succeeded syntactically",
        input
    );
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms[0].chirality, Some(Chirality::Unknown));
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
// MoleculeParser tests — whole-molecule invalids
fn test_molecule_branched_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case("[C:1]C", 2, 1)]
#[case("C[C:12]", 2, 1)]
fn test_molecule_bracket_class_valid_linear(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
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
#[case("F/C=C/F", Some(BondStereo::Trans))]
#[case("F/C=C\\F", Some(BondStereo::Cis))]
#[case("F\\C=C\\F", Some(BondStereo::Trans))]
#[case("F\\C=C/F", Some(BondStereo::Cis))]
fn test_molecule_ez_valid_branched(#[case] input: &str, #[case] expected: Option<BondStereo>) {
    let mut parse_state = ParseState::default();
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
#[case("C%0C", "invalid percent ring index")] // lexer should reject
#[case("C%09C", "invalid percent ring index leading zero")] // lexer should reject
fn test_molecule_ring_invalid_linear(
    mut parse_state: ParseState,
    #[case] input: &str,
    #[case] desc: &str,
) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case("C1C(C)CC1", 5, 5)]
#[case("C1CC(C)C1", 5, 5)]
fn test_molecule_rings_valid_branched(#[case] input: &str, #[case] atoms: usize, #[case] bonds: usize) {
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
#[case("[13C]1CC1", 3, 3)]
#[case("C1[C-]C1", 3, 3)]
// MoleculeParser tests — two-feature combinations (valid)
fn test_molecule_pair_bracket_with_rings_one_molecule_valid(
    #[case] input: &str,
    #[case] atoms: usize,
    #[case] bonds: usize,
) {
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
fn test_molecule_components_ring_merge(mut parse_state: ParseState) {
    // C1.C12.C2 should form a single 3-atom chain equivalent to CCC
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

#[test]
fn test_molecule_pair_components_and_bracket_valid() {
    let mut parse_state = ParseState::default();
    let input = "[NH4+].C[O-]";
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 2);
    assert_eq!(mols[0].atoms.len(), 1);
    assert_eq!(mols[0].bonds.len(), 0);
    assert_eq!(mols[1].atoms.len(), 2);
    assert_eq!(mols[1].bonds.len(), 1);
}

#[rstest]
#[case("[F-]/C=C\\[NH3+]", Some(BondStereo::Cis))]
#[case("F/C=C1CC\\1", Some(BondStereo::Either))]
#[case("F\\C=C1CC/1", Some(BondStereo::Either))]
fn test_molecule_pair_bracket_ez_rings_valid(
    #[case] input: &str,
    #[case] expected: Option<BondStereo>,
) {
    let mut parse_state = ParseState::default();
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

#[test]
// MoleculeParser tests — special edge cases
fn test_molecule_tetra_multi_centers_edge() {
    let mut parse_state = ParseState::default();
    let input = "Cl[C@H](F)C[C@@H](Cl)Br";
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    // Expect both marked
    assert_eq!(mols[0].atoms[1].chirality, Some(Chirality::Clockwise));
    assert_eq!(
        mols[0].atoms[4].chirality,
        Some(Chirality::CounterClockwise)
    );
}

// Non-tetrahedral chirality: SP/TB/OH

#[rstest]
#[case("[Pt@SP1](Cl)(Br)(I)F", true)]
#[case("[Pt@SP1](Cl)(Br)F", false)]
#[case("[Pt@SP1](Cl)(Br)(I)(F)N", false)]
fn test_molecule_square_planar_edge(#[case] input: &str, #[case] valid: bool) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let ch = mols[0].atoms[0].chirality;
    if valid {
        assert!(matches!(ch, Some(Chirality::SquarePlanar { .. })));
    } else {
        assert_eq!(ch, Some(Chirality::Unknown));
    }
}

#[rstest]
#[case("[P@TB1](F)(Cl)(Br)(I)N", true)]
#[case("[P@TB1](F)(Cl)(Br)(I)", false)]
#[case("[P@TB1](F)(Cl)(Br)(I)(N)O", false)]
fn test_molecule_trigonal_bipyramidal_edge(#[case] input: &str, #[case] valid: bool) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let ch = mols[0].atoms[0].chirality;
    if valid {
        assert!(matches!(ch, Some(Chirality::TrigonalBipyramidal { .. })));
    } else {
        assert_eq!(ch, Some(Chirality::Unknown));
    }
}

#[rstest]
#[case("[S@OH1](F)(Cl)(Br)(I)(N)(O)", true)]
#[case("[S@OH1](F)(Cl)(Br)(I)N", false)]
#[case("[S@OH1](F)(Cl)(Br)(I)NOF", false)]
fn test_molecule_octahedral_edge(#[case] input: &str, #[case] valid: bool) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let ch = mols[0].atoms[0].chirality;
    if valid {
        assert!(matches!(ch, Some(Chirality::Octahedral { .. })));
    } else {
        assert_eq!(ch, Some(Chirality::Unknown));
    }
}

// Allene (@AL) axial stereochemistry validation (structural viability only)
#[rstest]
// Valid: central allene carbon with two cumulated doubles; each terminal has a substituent
#[case("[C@AL1](=C([H]))=C([H])F", true)]
// Invalid: center does not have two double bonds
#[case("[C@AL1](=C)C([H])F", false)]
// Invalid: one terminal has no substituent beyond center
#[case("[C@AL1](=C)=C", false)]
fn test_molecule_allene_edge(#[case] input: &str, #[case] valid: bool) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let ch = mols[0].atoms[0].chirality;
    if valid {
        assert!(matches!(ch, Some(Chirality::Allenal { .. })));
    } else {
        assert_eq!(ch, Some(Chirality::Unknown));
    }
}

#[rstest]
// Alias: @ -> @AL1 when two incident double bonds present
#[case("[C@](=C([H]))=C([H])F", Some(1u32))]
// Alias: @@ -> @AL2 when two incident double bonds present
#[case("[C@@](=C([H]))=C([H])F", Some(2u32))]
// No alias: @ w/o allene axis (should stay tetra or downgrade per tetra rules)
#[case("[C@H](F)(Cl)Br", None)]
fn test_molecule_allene_alias_edge(#[case] input: &str, #[case] expect_arr: Option<u32>) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let ch = &mols[0].atoms[0].chirality;
    match expect_arr {
        Some(arr) => match ch {
            Some(Chirality::Allenal { arr: a }) => assert_eq!(*a, arr),
            _ => panic!("expected Allenal alias"),
        },
        None => {
            // should not be Allenal
            assert!(!matches!(ch, Some(Chirality::Allenal { .. })));
        }
    }
}

#[test]
fn test_molecule_ez_overlapping_markers_invalid() {
    let mut parse_state = ParseState::default();
    let input = "C/C=C//C=C/C";
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{} should be rejected", input);
}
