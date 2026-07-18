//! Tests for CTAB block parsing

use bstr::ByteSlice;
use nom::{Finish, Parser};
use pretty_assertions::assert_eq;
use rstest::*;
use umol_chem::element::Element;

use crate::ctfile::config::CtabParseFlags;
use crate::ctfile::error::ParseError;
use crate::ctfile::parser::{ctab_block, extended_ctab_block};
use crate::table_ir::{AtomList, AtomSymbol, BondOrder};

#[rustfmt::skip]
#[rstest]
#[case::diatomic(b"  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  3  0  0  0  0\nM  END\n",
    2, 1, Element::N, BondOrder::Triple)]
#[case::properties(b"  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  3  0  0  0  0\nM  CHG  1   2  -1\nM  END\n",
    2, 1, Element::N, BondOrder::Triple)]
#[case::short_lines(b"  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0\n    1.5400    0.0000    0.0000 C   0  0\n  1  2  1\nM  END\n",
    2, 1, Element::C, BondOrder::Single)]
#[case::no_terminal_newline(b"  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1\nM  END",
    2, 1, Element::C, BondOrder::Single)]
fn test_ctab_block(
    #[case] input: &[u8],
    #[case] expected_atoms: usize,
    #[case] expected_bonds: usize,
    #[case] expected_element0: Element,
    #[case] expected_order0: BondOrder,
) {
    let result = ctab_block(0, CtabParseFlags::BASIC).parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_ok(),
        "{:?} should parse successfully, got error: {:?}",
        input_str,
        result
    );
    let (remaining, (molecule, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} should consume all input, remaining: {:?}",
        input_str,
        remaining
    );
    assert_eq!(
        molecule.atom_count(),
        expected_atoms,
        "{:?}: atom count {} != expected {}",
        input_str,
        molecule.atom_count(),
        expected_atoms
    );
    assert_eq!(
        molecule.bond_count(),
        expected_bonds,
        "{:?}: bonds count {} != expected {}",
        input_str,
        molecule.bond_count(),
        expected_bonds
    );
    assert_eq!(
        molecule.atoms[0].element, Some(expected_element0),
        "{:?}: element0: {:?} != expected {:?}",
        input_str, molecule.atoms[0], expected_element0
    );
    assert_eq!(
        molecule.bonds[0].order, expected_order0,
        "{:?} bond order 0: {:?} != expected {:?}",
        input_str, molecule.bonds[0], expected_order0
    ) 
}

#[rstest]
#[case::als_query_property(b"  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  ALS   1  2 F Cl  Br\nM  END\n",
    ParseError::InvalidAtomLine { line: 1, col: 0 })]
#[case::legacy_atom_list(b"  2  1  1  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1 F    3   9   7   8  \nM  END\n",
    ParseError::UnsupportedLegacyAtomList { line: 4 })]
#[case::insufficient_atoms(b"  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n",
    ParseError::InvalidAtomLine { line: 2, col: 0 })]
#[case::insufficient_bonds(b"  2  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n",
    ParseError::InvalidBondLine { line: 4, col: 0 })]
#[case::missing_m_end(b"  2  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n",
    ParseError::MissingMEndTag { line: 3 })]
#[case::unicode_whitespace(b"24602\n\xA0 -OEChem-11060703412D\n\n  3  2 \xA00\xA0 0 \xA00 \xA00 \xA00 \xA00 \xA00  0999 V2000\n    2.5369   -0.1550    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0739    0.1550    0.0000 D   1  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.1550    0.0000 T   1  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\nM  ISO   2   2   2   3   3\nM  END\n",
    ParseError::InvalidCountsLine { line: 0 })]
fn test_ctab_block_invalid(#[case] input: &[u8], #[case] expected_error: ParseError) {
    let result = ctab_block(0, CtabParseFlags::BASIC).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    let error = result.finish().unwrap_err();
    assert_eq!(
        error, expected_error,
        "{:?}: error {:?} != {:?}",
        input_str, error, expected_error,
    );
}
#[rstest]
#[case::atom_list(b"  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  ALS   1  2 F Cl  Br\nM  END\n",
    2, 1, AtomSymbol::AtomList(AtomList { elements: vec![Element::Cl, Element::Br], exclusion: false }), BondOrder::Single)]
fn test_extended_ctab_block(
    #[case] input: &[u8],
    #[case] expected_atoms: usize,
    #[case] expected_bonds: usize,
    #[case] expected_symbol0: AtomSymbol,
    #[case] expected_order0: BondOrder,
) {
    let result = extended_ctab_block(0, CtabParseFlags::EXTENDED).parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_ok(),
        "{:?} should parse successfully: {:?}",
        input_str,
        result
    );
    let (remaining, (molecule, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} should consume all input, remaining: {:?}",
        input_str,
        remaining
    );
    assert_eq!(
        molecule.atom_count(),
        expected_atoms,
        "{:?}: atom count {} != expected {}",
        input_str,
        molecule.atom_count(),
        expected_atoms
    );
    assert_eq!(
        molecule.bond_count(),
        expected_bonds,
        "{:?}: bonds count {} != expected {}",
        input_str,
        molecule.bond_count(),
        expected_bonds
    );
    assert_eq!(
        molecule.atoms[0].symbol, expected_symbol0,
        "{:?}: atom symbol 0: {:?} != expected {:?}",
        input_str, molecule.atoms[0], expected_symbol0
    );
    assert_eq!(
        molecule.bonds[0].order, expected_order0,
        "{:?} bond order 0: {:?} != expected {:?}",
        input_str, molecule.bonds[0], expected_order0
    )
}

#[rstest]
#[case::legacy_atom_list(b"  2  1  1  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1 F    3   9   7   8  \nM  END\n",
    ParseError::UnsupportedLegacyAtomList { line: 4 })]
fn test_extended_ctab_block_invalid(#[case] input: &[u8], #[case] expected_error: ParseError) {
    let result = extended_ctab_block(0, CtabParseFlags::EXTENDED).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    let error = result.finish().unwrap_err();
    assert_eq!(
        error, expected_error,
        "{:?}: error {:?} != {:?}",
        input_str, error, expected_error,
    );
}
