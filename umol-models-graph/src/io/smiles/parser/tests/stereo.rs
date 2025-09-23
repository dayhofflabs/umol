use rstest::*;

use crate::io::smiles::lexer_old::Lexer;
use crate::io::ir::Chirality;
use crate::io::smiles::parser::grammar::MoleculeParser;
use crate::io::smiles::state::ParseState;
use super::fixtures::parse_state;

#[rstest]
#[case("[C@H](F)(Cl)Br", Chirality::Clockwise)]
#[case("[C@@H](F)(Cl)Br", Chirality::CounterClockwise)]
#[case("[C@TH1H](F)(Cl)Br", Chirality::Clockwise)]
#[case("[C@TH2H](F)(Cl)Br", Chirality::CounterClockwise)]
fn test_stereo_tetra(#[case] input: &str, #[case] expected: Chirality) {
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
fn test_stereo_tetra_insufficient_neighbors() {
    let mut parse_state = ParseState::default();
    let input = "[C@](F)(Cl)C";
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok(), "{} should have succeeded syntactically", input);
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atoms[0].chirality, Some(Chirality::Unknown));
}

#[rstest]
#[case("[Pt@SP1](Cl)(Br)(I)F", true)]
#[case("[Pt@SP1](Cl)(Br)F", false)]
#[case("[Pt@SP1](Cl)(Br)(I)(F)N", false)]
fn test_stereo_sp(#[case] input: &str, #[case] valid: bool) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let ch = mols[0].atoms[0].chirality;
    if valid { assert!(matches!(ch, Some(Chirality::SquarePlanar { .. }))); }
    else { assert_eq!(ch, Some(Chirality::Unknown)); }
}

#[rstest]
#[case("[P@TB1](F)(Cl)(Br)(I)N", true)]
#[case("[P@TB1](F)(Cl)(Br)(I)", false)]
#[case("[P@TB1](F)(Cl)(Br)(I)(N)O", false)]
fn test_stereo_tb(#[case] input: &str, #[case] valid: bool) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let ch = mols[0].atoms[0].chirality;
    if valid { assert!(matches!(ch, Some(Chirality::TrigonalBipyramidal { .. }))); }
    else { assert_eq!(ch, Some(Chirality::Unknown)); }
}

#[rstest]
#[case("[S@OH1](F)(Cl)(Br)(I)(N)(O)", true)]
#[case("[S@OH1](F)(Cl)(Br)(I)N", false)]
#[case("[S@OH1](F)(Cl)(Br)(I)NOF", false)]
fn test_stereo_oh(#[case] input: &str, #[case] valid: bool) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let ch = mols[0].atoms[0].chirality;
    if valid { assert!(matches!(ch, Some(Chirality::Octahedral { .. }))); }
    else { assert_eq!(ch, Some(Chirality::Unknown)); }
}

#[rstest]
#[case("[C@AL1](=C([H]))=C([H])F", true)]
#[case("[C@AL1](=C)C([H])F", false)]
#[case("[C@AL1](=C)=C", false)]
fn test_stereo_allene(#[case] input: &str, #[case] valid: bool) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let ch = mols[0].atoms[0].chirality;
    if valid { assert!(matches!(ch, Some(Chirality::Allenal { .. }))); }
    else { assert_eq!(ch, Some(Chirality::Unknown)); }
}

#[rstest]
#[case("[C@](=C([H]))=C([H])F", Some(1u32))]
#[case("[C@@](=C([H]))=C([H])F", Some(2u32))]
#[case("[C@H](F)(Cl)Br", None)]
fn test_stereo_allene_alias(#[case] input: &str, #[case] expect_arr: Option<u32>) {
    let mut parse_state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_ok());
    let mols = parse_state.drain_molecules();
    assert_eq!(mols.len(), 1);
    let ch = &mols[0].atoms[0].chirality;
    match expect_arr {
        Some(arr) => match ch { Some(Chirality::Allenal { arr: a }) => assert_eq!(*a, arr), _ => panic!("expected Allenal alias"), },
        None => assert!(!matches!(ch, Some(Chirality::Allenal { .. }))),
    }
}

#[rstest]
#[case("C/C=C//C=C/C", "overlapping E/Z markers should be rejected")]
fn test_stereo_ez_overlapping_invalid(mut parse_state: ParseState, #[case] input: &str, #[case] desc: &str) {
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    let result = parser.parse(&mut parse_state, lexer);
    assert!(result.is_err(), "{}", desc);
}


