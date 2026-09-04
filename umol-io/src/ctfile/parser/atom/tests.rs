#![allow(clippy::too_many_arguments)]

use bstr::ByteSlice;
use float_cmp::*;
use pretty_assertions::assert_eq;
use rstest::*;
use winnow::error::ErrMode;
use winnow::Parser;

use super::*;
use crate::ctfile::config::CtabParseFlags;
use crate::table_ir::atom::Chirality;
use crate::table_ir::{
    AtomExactChange, AtomInversionRetention, AtomStereoCare, RGroup, WildcardAtom,
};

#[rustfmt::skip]
#[rstest]
#[case::len_71_trailing_whitespace(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0  \n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_69(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0\n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_62(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  \n",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_61(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0 \n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_60(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0\n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_49(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0 \n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_48(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0\n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_40(b"    1.0000    2.0000    3.0000 C  -2  3 \n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_39(b"    1.0000    2.0000    3.0000 C  -2  3\n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_37(b"    1.0000    2.0000    3.0000 C  -2 \n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), None, None, None, None)]
#[case::len_36(b"    1.0000    2.0000    3.0000 C  -2\n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), None, None, None, None)]
#[case::len_35(b"    1.0000    2.0000    3.0000 C   \n",
    1.0000, 2.0000, 3.0000, Element::C, None, None, None, None, None)]
#[case::len_34(b"    1.0000    2.0000    3.0000 C  \n",
    1.0000, 2.0000, 3.0000, Element::C, None, None, None, None, None)]
#[case::len_32(b"    1.0000    2.0000    3.0000 C\n",
    1.0000, 2.0000, 3.0000, Element::C, None, None, None, None, None)]
fn test_atom_block(
    #[case] input: &[u8],
    #[case] x: f64,
    #[case] y: f64,
    #[case] z: f64,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] chirality: Option<Chirality>,
    #[case] hydrogen_count: Option<u8>,
    #[case] valence: Option<u8>,
) {
    let mut remaining = input;
    let result = atom_block(&mut remaining, 1, 0, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", input_str, result.clone().unwrap_err());
    let (atoms, positions, _) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} should have consumed all input, remaining: {:?}", input_str, remaining);
    assert_eq!(atoms.len(), 1);
    let atom = &atoms[0];
    let positions = positions.expect("positions should be present");
    let position = &positions[0];
    assert!(approx_eq!(f64, position.x, x), "{:?} x: {:?} != {:?}", input_str, position.x, x);
    assert!(approx_eq!(f64, position.y, y), "{:?} y: {:?} != {:?}", input_str, position.y, y);
    assert!(approx_eq!(f64, position.z, z), "{:?} z: {:?} != {:?}", input_str, position.z, z);
    assert_eq!(atom.element, Some(element), "{:?} element", input_str);
    assert_eq!(atom.isotope_mass, isotope_mass, "{:?} isotope_mass", input_str);
    assert_eq!(atom.charge, charge, "{:?} charge", input_str);
    assert_eq!(atom.chirality, chirality, "{:?} chirality", input_str);
    assert_eq!(atom.implicit_hydrogens, hydrogen_count, "{:?} hydrogen_count", input_str);
    assert_eq!(atom.valence, valence, "{:?} valence", input_str);
}

#[rustfmt::skip]
#[rstest]
#[case::len_69_trailing_data(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0X\n", 69)]
#[case::len_62_trailing_data(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0 1\n", 61)]
#[case::len_61_trailing_data(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  01\n", 60)]
#[case::len_60_trailing_data(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0X\n", 60)]
#[case::len_50_trailing_data(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0 4\n", 49)]
#[case::len_49_trailing_data(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  04\n", 48)]
#[case::len_48_trailing_data(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0X\n", 48)]
#[case::len_39_trailing_data(b"    1.0000    2.0000    3.0000 C  -2  3X\n", 39)]
#[case::len_38_trailing_data(b"    1.2345    2.3456    3.4567 C  -2 3\n", 37)]
#[case::len_37_trailing_data(b"    1.2345    2.3456    3.4567 C  -23\n", 36)]
#[case::len_36_trailing_data(b"    1.0000    2.0000    3.0000 C  -2X\n", 36)]
#[case::len_35_trailing_data(b"    1.2345    2.3456    3.4567 C  1\n", 34)]
#[case::len_34_trailing_data(b"    1.0000    2.0000    3.0000 C  X\n", 34)]
#[case::len_31_too_short(b"    1.0000    2.0000    3.0000 \n", 31)]
fn test_atom_block_error(#[case] input: &[u8], #[case] col: u32) {
    let mut remaining = input;
    assert_eq!(
        atom_block(&mut remaining, 1, 0, CtabParseFlags::BASIC),
        Err(ErrMode::Cut(ParseError::InvalidAtomLine { line: 0, col }))
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_71_trailing_whitespace(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0  \n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_69(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0\n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_62(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  \n",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_61(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0 \n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_60(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0\n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_49(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0 \n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_48(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0\n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_40(b"    1.0000    2.0000    3.0000 C  -2  3 \n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_39(b"    1.0000    2.0000    3.0000 C  -2  3\n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_37(b"    1.0000    2.0000    3.0000 C  -2 \n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), None, None, None, None)]
#[case::len_36(b"    1.0000    2.0000    3.0000 C  -2\n",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), None, None, None, None)]
#[case::len_35(b"    1.0000    2.0000    3.0000 C   \n",
    1.0000, 2.0000, 3.0000, Element::C, None, None, None, None, None)]
#[case::len_34(b"    1.0000    2.0000    3.0000 C  \n",
    1.0000, 2.0000, 3.0000, Element::C, None, None, None, None, None)]
#[case::len_32(b"    1.0000    2.0000    3.0000 C\n",
    1.0000, 2.0000, 3.0000, Element::C, None, None, None, None, None)]
fn test_extended_atom_block(
    #[case] input: &[u8],
    #[case] x: f64,
    #[case] y: f64,
    #[case] z: f64,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] chirality: Option<Chirality>,
    #[case] hydrogen_count: Option<u8>,
    #[case] valence: Option<u8>,
) {
    let mut remaining = input;
    let result = extended_atom_block(&mut remaining, 1, 0, CtabParseFlags::EXTENDED);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", input_str, result.clone().unwrap_err());
    let (atoms, positions, _) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} should have consumed all input, remaining: {:?}", input_str, remaining);
    assert_eq!(atoms.len(), 1);
    let atom = &atoms[0];
    let positions = positions.expect("positions should be present");
    let position = &positions[0];
    assert!(approx_eq!(f64, position.x, x), "{:?} x: {:?} != {:?}", input_str, position.x, x);
    assert!(approx_eq!(f64, position.y, y), "{:?} y: {:?} != {:?}", input_str, position.y, y);
    assert!(approx_eq!(f64, position.z, z), "{:?} z: {:?} != {:?}", input_str, position.z, z);
    assert_eq!(atom.symbol, AtomSymbol::Element(element), "{:?} symbol", input_str);
    assert_eq!(atom.isotope_mass, isotope_mass, "{:?} isotope_mass", input_str);
    assert_eq!(atom.charge, charge, "{:?} charge", input_str);
    assert_eq!(atom.chirality, chirality, "{:?} chirality", input_str);
    assert_eq!(atom.implicit_hydrogens, hydrogen_count, "{:?} hydrogen_count", input_str);
    assert_eq!(atom.valence, valence, "{:?} valence", input_str);
}

#[rustfmt::skip]
#[rstest]
#[case::len_70_malformed(b"   -1.9225   -0.6187    0.0000 R20  0  0  0  0  0  0  0  0  0  0  0  0\n", 69)]
#[case::len_69_trailing_data(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0X\n", 69)]
#[case::len_60_trailing_data(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0X\n", 60)]
#[case::len_48_trailing_data(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0X\n", 48)]
#[case::len_44_trailing_data(b"    1.2345    2.3456    3.4567 C   0  3  0 1", 43)]
#[case::len_43_trailing_data(b"    1.2345    2.3456    3.4567 C   0  3  01", 42)]
#[case::len_41_trailing_data(b"    1.2345    2.3456    3.4567 C   0  3 1", 40)]
#[case::len_40_trailing_data(b"    1.2345    2.3456    3.4567 C   0  31", 39)]
#[case::len_39_trailing_data(b"    1.0000    2.0000    3.0000 C  -2  3X\n", 39)]
#[case::len_36_trailing_data(b"    1.0000    2.0000    3.0000 C  -2X\n", 36)]
#[case::len_34_trailing_data(b"    1.0000    2.0000    3.0000 C  X\n", 34)]
#[case::len_31_too_short(b"    1.0000    2.0000    3.0000 \n", 31)]
fn test_extended_atom_block_error(#[case] input: &[u8], #[case] col: u32) {
    let mut remaining = input;
    assert_eq!(
        extended_atom_block(&mut remaining, 1, 0, CtabParseFlags::EXTENDED),
        Err(ErrMode::Cut(ParseError::InvalidAtomLine { line: 0, col }))
    );
}

#[rstest]
#[case::basic(false)]
#[case::extended(true)]
fn test_atom_block_eof_error(#[case] extended: bool) {
    let mut input: &[u8] = b"";
    let result = if extended {
        extended_atom_block(&mut input, 1, 4, CtabParseFlags::EXTENDED).map(|_| ())
    } else {
        atom_block(&mut input, 1, 4, CtabParseFlags::BASIC).map(|_| ())
    };
    assert_eq!(
        result,
        Err(ErrMode::Cut(ParseError::UnexpectedEof {
            line: 4,
            block: "atom",
        }))
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_69(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_66(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_63(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_60(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_57(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_54(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_51(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_48(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_45(b"    1.2345    2.3456    3.4567 C   0  3  1  0",
    1.2345, 2.3456, 3.4567, Element::C, None, Some(1), Some(Chirality::Clockwise), None, None)]
#[case::len_42(b"    1.2345    2.3456    3.4567 C   0  3  1",
    1.2345, 2.3456, 3.4567, Element::C, None, Some(1), Some(Chirality::Clockwise), None, None)]
#[case::len_39(b"    1.0000    2.0000    3.0000 C  -2  3",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_36(b"    1.0000    2.0000    3.0000 C  -2",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), None, None, None, None)]
#[case::len_34(b"    1.0000    2.0000    3.0000 C  ",
    1.0000, 2.0000, 3.0000, Element::C, None, None, None, None, None)]
#[case::len_33(b"    1.0000    2.0000    3.0000 C ",
    1.0000, 2.0000, 3.0000, Element::C, None, None, None, None, None)]
#[case::len_32(b"    1.0000    2.0000    3.0000 C",
    1.0000, 2.0000, 3.0000, Element::C, None, None, None, None, None)]
#[case::mass_diff_lower_bound(b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4  0  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(9), Some(1), None, None, Some(4))]
#[case::mass_diff_upper_bound(b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4  0  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(16), Some(1), None, None, Some(4))]
#[case::mass_diff_out_of_range_low(b"    1.2345    2.3456    3.4567 C  -4  3  0  0  0  4  0  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, None, Some(1), None, None, Some(4))]
#[case::mass_diff_out_of_range_high(b"    1.2345    2.3456    3.4567 C   5  3  0  0  0  4  0  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, None, Some(1), None, None, Some(4))]
#[case::charge_out_of_range_high(b"    1.2345    2.3456    3.4567 C  -2  8  0  0  0  4  0  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), None, None, None, Some(4))]
#[case::blank_mass_diff(b"    1.2345    2.3456    3.4567 C      3  0  0  0  4",
    1.2345, 2.3456, 3.4567, Element::C, None, Some(1), None, None, Some(4))]
#[case::blank_charge(b"    1.2345    2.3456    3.4567 C  -2     0  0  0  4",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), None, None, None, Some(4))]
#[case::blank_valence(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0   ", 
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, None)]
#[case::blank_stereo_parity(b"    1.2345    2.3456    3.4567 C   0  3   0  ",
    1.2345, 2.3456, 3.4567, Element::C, None, Some(1), None, None, None)]
#[case::blank_block_1(b"    1.2345    2.3456    3.4567 C  -2  3  0     0  4  0  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::blank_block_2(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4     0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::blank_block_3(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0      ",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::gaps_with_spaces_and_zeros(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4           1 0    ",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::zero_coordinates(b"    0.0000    0.0000    0.0000 C  -2  3  0  0  0  4",
    0.0000, 0.0000, 0.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::hydrogen_count(b"    0.0000    0.0000    0.0000 C   0  0  0  1  0  0  0  0  0  0  0  0",
    0.0000, 0.0000, 0.0000, Element::C, None, None, None, Some(0), None)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::H, Some(2), Some(1), None, None, Some(1))]
#[case::invalid_unused(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0XXX  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
fn test_atom_input(
    #[case] input: &[u8],
    #[case] x: f64,
    #[case] y: f64,
    #[case] z: f64,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] chirality: Option<Chirality>,
    #[case] hydrogen_count: Option<u8>,
    #[case] valence: Option<u8>,
) {
    let result = atom_input(CtabParseFlags::BASIC).parse(Input::new(input));
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", input_str, result.clone().unwrap_err());
    let (atom, position) = result.unwrap();
    assert!(
        approx_eq!(f64, position.x, x),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x,
    );
    assert!(
        approx_eq!(f64, position.y, y),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y,
    );
    assert!(
        approx_eq!(f64, position.z, z),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z,
    );
    assert_eq!(
        atom.element, Some(element),
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope_mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass
    );
    assert_eq!(
        atom.charge, charge,
        "{:?} has returned charge {:?}, expected {:?}",
        input_str, atom.charge, charge
    );
    assert_eq!(
        atom.chirality, chirality,
        "{:?} has returned chirality {:?}, expected {:?}",
        input_str, atom.chirality, chirality
    );
    assert_eq!(
        atom.implicit_hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.implicit_hydrogens, hydrogen_count
    );
    assert_eq!(
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence
    );
}

#[rustfmt::skip]
#[rstest]
#[case::too_short(b"    1.2345    2.3456    3.4567", 30)]
#[case::non_numeric_coordinate(b"    1.234a    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0", 0)]
#[case::malformed_coordinate(b"    0.1   0.0    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0", 10)]
#[case::atom_list(b"    0.7145    2.0625    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0", 31)]
#[case::pseudoatom(b"   -1.8857    2.4750    0.0000 Psd 0  0  0  0  0  0  0  0  0  0  0  0", 31)]
#[case::stereo_parity_out_of_range(b"    1.2345    2.3456    3.4567 C   0  3  4", 39)]
#[case::non_numeric_valence(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  a", 48)]
#[case::valence_out_of_range(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0 16", 48)]
#[case::non_zero_hydrogen_count_stereo_care(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4", 45)]
#[case::non_numeric_atom_map_number(b"    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  a  0  0", 60)]
#[case::len69_invalid_extended(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0XXX", 63)]
#[case::len69_non_zero_extended(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  1", 63)]
#[case::len48_invalid_extended(b"    1.2345    2.3456    3.4567 C   0  3  0  0XXX", 45)]
fn test_atom_input_error(
    #[case] input: &[u8],
    #[case] column: u32,
) {
    let error = atom_input(CtabParseFlags::BASIC)
        .parse(Input::new(input))
        .unwrap_err()
        .into_inner();
    assert_eq!(error, InputError { column });
}

#[rustfmt::skip]
#[rstest]
#[case::len_69(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_51(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_42(b"    1.0000    2.0000    3.0000 C  -2  3  1",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), Some(Chirality::Clockwise), None, None)]
#[case::len_39(b"    1.0000    2.0000    3.0000 C  -2  3",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), Some(1), None, None, None)]
#[case::len_36(b"    1.0000    2.0000    3.0000 C  -2",
    1.0000, 2.0000, 3.0000, Element::C, Some(10), None, None, None, None)]
#[case::len_34(b"    1.0000    2.0000    3.0000 C  ",
    1.0000, 2.0000, 3.0000, Element::C, None, None, None, None, None)]
#[case::mass_diff_lower_bound(b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4  0  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(9), Some(1), None, None, Some(4))]
#[case::mass_diff_upper_bound(b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4  0  0  0  0  0  0",
    1.2345, 2.3456, 3.4567, Element::C, Some(16), Some(1), None, None, Some(4))]
fn test_atom_input_strict(
    #[case] input: &[u8],
    #[case] x: f64,
    #[case] y: f64,
    #[case] z: f64,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] chirality: Option<Chirality>,
    #[case] hydrogen_count: Option<u8>,
    #[case] valence: Option<u8>,
) {
    let result = atom_input(CtabParseFlags::BASIC).parse(Input::new(input));
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", input_str, result.clone().unwrap_err());
    let (atom, position) = result.unwrap();

    assert!(
        approx_eq!(f64, position.x, x),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x,
    );
    assert!(
        approx_eq!(f64, position.y, y),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y,
    );
    assert!(
        approx_eq!(f64, position.z, z),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z,
    );
    assert_eq!(
        atom.element, Some(element),
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope_mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass
    );
    assert_eq!(
        atom.charge, charge,
        "{:?} has returned charge {:?}, expected {:?}",
        input_str, atom.charge, charge
    );
    assert_eq!(
        atom.chirality, chirality,
        "{:?} has returned chirality {:?}, expected {:?}",
        input_str, atom.chirality, chirality
    );
    assert_eq!(
        atom.implicit_hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.implicit_hydrogens, hydrogen_count
    );
    assert_eq!(
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_69_invalid_unused(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0XXX  0  0  0  0", 51)]
#[case::len_51_invalid_unused(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4", 42)]
#[case::invalid_extended(b"    1.2345    2.3456    3.4567 C  -2  3  0  0XXX  4  0  0  0  0  0  0", 45)]
#[case::invalid_unused(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0XXX", 51)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0", 31)]
#[case::non_zero_atom_map_number(b"    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  1  0  0", 60)]
fn test_atom_input_strict_error(
    #[case] input: &[u8],
    #[case] column: u32,
) {
    let error = atom_input(CtabParseFlags::BASIC & CtabParseFlags::STRICT)
        .parse(Input::new(input))
        .unwrap_err()
        .into_inner();
    assert_eq!(error, InputError { column });
}

#[rustfmt::skip]
#[rstest]
#[case::len_69(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0",
    Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::len_51(b"    1.0000    2.0000    3.0000 C  -2  3  1  0  0  4",
    Element::C, Some(10), Some(1), Some(Chirality::Clockwise), None, Some(4))]
#[case::len_42(b"    1.0000    2.0000    3.0000 C  -2  3  1",
    Element::C, Some(10), Some(1), Some(Chirality::Clockwise), None, None)]
#[case::len_39(b"    1.0000    2.0000    3.0000 C  -2  3",
    Element::C, Some(10), Some(1), None, None, None )]
#[case::len_36(b"    1.0000    2.0000    3.0000 C  -2",
    Element::C, Some(10), None, None, None, None)]
#[case::len_34(b"    1.0000    2.0000    3.0000 C  ",
    Element::C, None, None, None, None, None)]
#[case::zero_coordinates(b"    0.0000    0.0000    0.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0",
    Element::C, Some(10), Some(1), None, None, Some(4))]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0",
    Element::H, Some(2), Some(1), None, None, Some(1))]
#[case::invalid_coordinates(b"   invalid   invalid   invalid C  -2  3  0  0  0  4  0  0  0  0  0  0",
    Element::C, Some(10), Some(1), None, None, Some(4))]
fn test_atom_input_ignore_positions(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] chirality: Option<Chirality>,
    #[case] hydrogen_count: Option<u8>,
    #[case] valence: Option<u8>,
) {
    let result = atom_input(CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS)
        .parse(Input::new(input));
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", input_str, result.clone().unwrap_err());
    let (atom, _position) = result.unwrap();
    assert_eq!(
        atom.element, Some(element),
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope_mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass
    );
    assert_eq!(
        atom.charge, charge,
        "{:?} has returned charge {:?}, expected {:?}",
        input_str, atom.charge, charge
    );
    assert_eq!(
        atom.chirality, chirality,
        "{:?} has returned chirality {:?}, expected {:?}",
        input_str, atom.chirality, chirality
    );
    assert_eq!(
        atom.implicit_hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.implicit_hydrogens, hydrogen_count
    );
    assert_eq!(
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len69_invalid_coordinates(b"   invalid   invalid   invalid C  -2  3  0  0  0  4  0  0  0  0  0  0", 0)]
#[case::len60_invalid_coordinates(b"   invalid   invalid   invalid C  -2  3  0  0  0  4  0  0  0", 0)]
#[case::len42_invalid_coordinates(b"   invalid   invalid   invalid C  -2  3  0", 0)]
#[case::len39_invalid_coordinates(b"   invalid   invalid   invalid C  -2  3", 0)]
#[case::len36_invalid_coordinates(b"   invalid   invalid   invalid C  -2", 0)]
#[case::len34_invalid_coordinates(b"   invalid   invalid   invalid C  ", 0)]
fn test_atom_input_ignore_positions_strict_error(
    #[case] input: &[u8],
    #[case] column: u32,
) {
    let error = atom_input(
        CtabParseFlags::BASIC & CtabParseFlags::STRICT | CtabParseFlags::IGNORE_POSITIONS,
    )
    .parse(Input::new(input))
    .unwrap_err()
    .into_inner();
    assert_eq!(error, InputError { column });
}

#[rustfmt::skip]
#[rstest]
#[case::len_32(b"    1.0000    2.0000    3.0000 C",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), None, None, None, None, None, None, None, None, None)]
#[case::len_34(b"    1.0000    2.0000    3.0000 C  ",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), None, None, None, None, None, None, None, None, None)]
#[case::len_36(b"    1.0000    2.0000    3.0000 C  -2",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), Some(10), None, None, None, None, None, None, None, None)]
#[case::len_39(b"    1.0000    2.0000    3.0000 C  -2  3",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), Some(10), Some(1), None, None, None, None, None, None, None)]
#[case::len_48(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), Some(10), Some(1), None, None, None, None, None, None, None)]
#[case::len_60(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), Some(10), Some(1), None, None, Some(4), None, None, None, None)]
#[case::len_69(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), Some(10), Some(1), None, None, Some(4), None, None, None, None)]
#[case::stereo_parity_and_query_fields(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4  0  0  0  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), Some(10), Some(1), Some(Chirality::Clockwise), Some(1), Some(4), Some(AtomStereoCare::Care), None, None, None)]
#[case::atom_list(b"    1.0000    2.0000    3.0000 L   0  0  0  0  0  4  0  0  0  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }), None, None, None, None, Some(4), None, None, None, None)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0",
       1.2345, 2.3456, 3.4567, AtomSymbol::NamedIsotope(NamedIsotope::D), Some(2), Some(1), None, None, Some(1), None, None, None, None)]
#[case::invalid_unused(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0XXX  0  0  0",
       1.2345, 2.3456, 3.4567, AtomSymbol::Element(Element::C), Some(10), Some(1), None, None, Some(4), None, None, None, None)]
fn test_extended_atom_input(
    #[case] input: &[u8],
    #[case] x: f64,
    #[case] y: f64,
    #[case] z: f64,
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] chirality: Option<Chirality>,
    #[case] hydrogen_count: Option<u8>,
    #[case] valence: Option<u8>,
    #[case] stereo_care: Option<AtomStereoCare>,
    #[case] atom_map_num: Option<u32>,
    #[case] inversion_retention: Option<AtomInversionRetention>,
    #[case] exact_change: Option<AtomExactChange>,
) {
    let result = extended_atom_input(CtabParseFlags::EXTENDED).parse(Input::new(input));
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (atom, position) = result.unwrap();
    assert!(
        approx_eq!(f64, position.x, x),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x,
    );
    assert!(
        approx_eq!(f64, position.y, y),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y,
    );
    assert!(
        approx_eq!(f64, position.z, z),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z,
    );
    assert_eq!(
        atom.symbol, symbol,
        "{:?} has returned symbol {:?}, expected {:?}",
        input_str, atom.symbol, symbol,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{:?} has returned charge {:?}, expected {:?}",
        input_str, atom.charge, charge,
    );
    assert_eq!(
        atom.chirality, chirality,
        "{:?} has returned chirality {:?}, expected {:?}",
        input_str, atom.chirality, chirality,
    );
    assert_eq!(
        atom.implicit_hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.implicit_hydrogens, hydrogen_count,
    );
    assert_eq!(
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{:?} has returned stereo_care {:?}, expected {:?}",
        input_str, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.class, atom_map_num,
        "{:?} has returned atom_map_num {:?}, expected {:?}",
        input_str, atom.class, atom_map_num,
    );
    assert_eq!(
        atom.inversion_retention, inversion_retention,
        "{:?} has returned inversion_retention {:?}, expected {:?}",
        input_str, atom.inversion_retention, inversion_retention,
    );
    assert_eq!(
        atom.exact_change, exact_change,
        "{:?} has returned exact_change {:?}, expected {:?}",
        input_str, atom.exact_change, exact_change,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::invalid_stereo_parity(b"    1.0000    2.0000    3.0000 C  -2  3  4", 39)]
#[case::negative_hydrogen_count(b"    1.4289    0.8250    0.0000 C   0  0  0 -1  0  0  0  0  0  0  0  0", 42)]
#[case::non_numeric_hydrogen_count(b"    1.0000    2.0000    3.0000 C  -2  3  0  a", 42)]
#[case::non_numeric_valence(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  a", 48)]
#[case::stereo_care_out_of_range(b"    1.0000    2.0000    3.0000 C   0  0  0  0  2", 45)]
#[case::exact_change_out_of_range(b"    1.0000    2.0000    3.0000 C   0  0  0  0  0  0  0  0  0  0  0  2", 66)]
#[case::chemaxon_wildcard(b"    1.2345    2.3456    3.4567 AH 0  0  0  0  0  0  0  0  0  0  0  0", 31)]
#[case::pseudoatom(b"    1.2345    2.3456    3.4567 Ala 0  0  0  0  0  0  0  0  0  0  0  0", 31)]
#[case::len_30_too_short(b"    1.2345    2.3456    3.4567", 30)]
fn test_extended_atom_input_error(
    #[case] input: &[u8],
    #[case] column: u32,
) {
    let error = extended_atom_input(CtabParseFlags::EXTENDED)
        .parse(Input::new(input))
        .unwrap_err()
        .into_inner();
    assert_eq!(error, InputError { column });
}

#[rustfmt::skip]
#[rstest]
#[case::len_69(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), Some(10), Some(1), None, None, Some(4), None, None, None)]
#[case::stereo_parity_and_query_fields(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4  0  0  0  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), Some(10), Some(1), Some(Chirality::Clockwise), Some(1), Some(4), Some(AtomStereoCare::Care), None, None)]
#[case::reaction_fields(b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0           1  2  1",
       0.0000, 0.0000, 0.0000, AtomSymbol::Element(Element::C), None, None, None, None, None, None, Some(AtomInversionRetention::Retained), Some(AtomExactChange::Match))]
#[case::atom_list(b"    1.0000    2.0000    3.0000 L   0  0  0  0  0  4  0  0  0  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }), None, None, None, None, Some(4), None, None, None)]
fn test_extended_atom_input_strict(
    #[case] input: &[u8],
    #[case] x: f64,
    #[case] y: f64,
    #[case] z: f64,
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] chirality: Option<Chirality>,
    #[case] hydrogen_count: Option<u8>,
    #[case] valence: Option<u8>,
    #[case] stereo_care: Option<AtomStereoCare>,
    #[case] inversion_retention: Option<AtomInversionRetention>,
    #[case] exact_change: Option<AtomExactChange>,
) {
    let result = extended_atom_input(CtabParseFlags::STRICT).parse(Input::new(input));
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (atom, position) = result.unwrap();
    assert!(
        approx_eq!(f64, position.x, x),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x,
    );
    assert!(
        approx_eq!(f64, position.y, y),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y,
    );
    assert!(
        approx_eq!(f64, position.z, z),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z,
    );
    assert_eq!(
        atom.symbol, symbol,
        "{:?} has returned symbol {:?}, expected {:?}",
        input_str, atom.symbol, symbol,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{:?} has returned charge {:?}, expected {:?}",
        input_str, atom.charge, charge,
    );
    assert_eq!(
        atom.chirality, chirality,
        "{:?} has returned chirality {:?}, expected {:?}",
        input_str, atom.chirality, chirality,
    );
    assert_eq!(
        atom.implicit_hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.implicit_hydrogens, hydrogen_count,
    );
    assert_eq!(
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{:?} has returned stereo_care {:?}, expected {:?}",
        input_str, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.inversion_retention, inversion_retention,
        "{:?} has returned inversion_retention {:?}, expected {:?}",
        input_str, atom.inversion_retention, inversion_retention,
    );
    assert_eq!(
        atom.exact_change, exact_change,
        "{:?} has returned exact_change {:?}, expected {:?}",
        input_str, atom.exact_change, exact_change,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::invalid_unused(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0XXX  0  0  0", 51)]
fn test_extended_atom_input_strict_error(
    #[case] input: &[u8],
    #[case] column: u32,
) {
    let error = extended_atom_input(CtabParseFlags::STRICT)
        .parse(Input::new(input))
        .unwrap_err()
        .into_inner();
    assert_eq!(error, InputError { column });
}

#[rustfmt::skip]
#[rstest]
#[case::len_69(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), Some(10), Some(1), None, None, Some(4), None, None, None)]
#[case::stereo_parity_and_query_fields(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4  0  0  0  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::Element(Element::C), Some(10), Some(1), Some(Chirality::Clockwise), Some(1), Some(4), Some(AtomStereoCare::Care), None, None)]
#[case::atom_list(b"    1.0000    2.0000    3.0000 L   0  0  0  0  0  4  0  0  0  0  0  0",
       1.0000, 2.0000, 3.0000, AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }), None, None, None, None, Some(4), None, None, None)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0",
       1.2345, 2.3456, 3.4567, AtomSymbol::NamedIsotope(NamedIsotope::D), Some(2), Some(1), None, None, Some(1), None, None, None)]
#[case::chemaxon_wildcard(b"    1.2345    2.3456    3.4567 AH  0  0  0  0  0  0  0  0  0  0  0  0",
       1.2345, 2.3456, 3.4567, AtomSymbol::WildcardAtom(WildcardAtom::HeavyOrH), None, None, None, None, None, None, None, None)]
#[case::invalid_unused(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0XXX  0  0  0",
       1.2345, 2.3456, 3.4567, AtomSymbol::Element(Element::C), Some(10), Some(1), None, None, Some(4), None, None, None)]
fn test_extended_atom_input_lenient(
    #[case] input: &[u8],
    #[case] x: f64,
    #[case] y: f64,
    #[case] z: f64,
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] chirality: Option<Chirality>,
    #[case] hydrogen_count: Option<u8>,
    #[case] valence: Option<u8>,
    #[case] stereo_care: Option<AtomStereoCare>,
    #[case] inversion_retention: Option<AtomInversionRetention>,
    #[case] exact_change: Option<AtomExactChange>,
) {
    let result = extended_atom_input(CtabParseFlags::LENIENT).parse(Input::new(input));
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (atom, position) = result.unwrap();
    assert!(
        approx_eq!(f64, position.x, x),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x,
    );
    assert!(
        approx_eq!(f64, position.y, y),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y,
    );
    assert!(
        approx_eq!(f64, position.z, z),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z,
    );
    assert_eq!(
        atom.symbol, symbol,
        "{:?} has returned symbol {:?}, expected {:?}",
        input_str, atom.symbol, symbol,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{:?} has returned charge {:?}, expected {:?}",
        input_str, atom.charge, charge,
    );
    assert_eq!(
        atom.chirality, chirality,
        "{:?} has returned chirality {:?}, expected {:?}",
        input_str, atom.chirality, chirality,
    );
    assert_eq!(
        atom.implicit_hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.implicit_hydrogens, hydrogen_count,
    );
    assert_eq!(
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{:?} has returned stereo_care {:?}, expected {:?}",
        input_str, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.inversion_retention, inversion_retention,
        "{:?} has returned inversion_retention {:?}, expected {:?}",
        input_str, atom.inversion_retention, inversion_retention,
    );
    assert_eq!(
        atom.exact_change, exact_change,
        "{:?} has returned exact_change {:?}, expected {:?}",
        input_str, atom.exact_change, exact_change,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::invalid_stereo_parity(b"    1.0000    2.0000    3.0000 C  -2  3  4", 39)]
#[case::non_numeric_hydrogen_count(b"    1.0000    2.0000    3.0000 C  -2  3  0  a", 42)]
#[case::non_numeric_valence(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  a", 48)]
fn test_extended_atom_input_lenient_error(
    #[case] input: &[u8],
    #[case] column: u32,
) {
    let error = extended_atom_input(CtabParseFlags::LENIENT)
        .parse(Input::new(input))
        .unwrap_err()
        .into_inner();
    assert_eq!(error, InputError { column });
}

#[rustfmt::skip]
#[rstest]
#[case::len_69(b"    1.2345    2.3456    3.4567 Ala 0  0  0  0  0  0  0  0  0  0  0  0",
       1.2345, 2.3456, 3.4567, AtomSymbol::Pseudoatom(String::from("Ala")), None, None, None, None, None, None, None, None)]
fn test_extended_atom_input_pseudoatoms(
    #[case] input: &[u8],
    #[case] x: f64,
    #[case] y: f64,
    #[case] z: f64,
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] chirality: Option<Chirality>,
    #[case] hydrogen_count: Option<u8>,
    #[case] valence: Option<u8>,
    #[case] stereo_care: Option<AtomStereoCare>,
    #[case] inversion_retention: Option<AtomInversionRetention>,
    #[case] exact_change: Option<AtomExactChange>,
) {
    let flags = CtabParseFlags::EXTENDED | CtabParseFlags::PSEUDOATOMS;
    let result = extended_atom_input(flags).parse(Input::new(input));
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (atom, position) = result.unwrap();
    assert!(
        approx_eq!(f64, position.x, x),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x,
    );
    assert!(
        approx_eq!(f64, position.y, y),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y,
    );
    assert!(
        approx_eq!(f64, position.z, z),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z,
    );
    assert_eq!(
        atom.symbol, symbol,
        "{:?} has returned symbol {:?}, expected {:?}",
        input_str, atom.symbol, symbol,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{:?} has returned charge {:?}, expected {:?}",
        input_str, atom.charge, charge,
    );
    assert_eq!(
        atom.chirality, chirality,
        "{:?} has returned chirality {:?}, expected {:?}",
        input_str, atom.chirality, chirality,
    );
    assert_eq!(
        atom.implicit_hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.implicit_hydrogens, hydrogen_count,
    );
    assert_eq!(
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{:?} has returned stereo_care {:?}, expected {:?}",
        input_str, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.inversion_retention, inversion_retention,
        "{:?} has returned inversion_retention {:?}, expected {:?}",
        input_str, atom.inversion_retention, inversion_retention,
    );
    assert_eq!(
        atom.exact_change, exact_change,
        "{:?} has returned exact_change {:?}, expected {:?}",
        input_str, atom.exact_change, exact_change,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_69(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), None, None, None, None, None)]
#[case::zero_coordinates(b"    0.0000    0.0000    0.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), None, None, None, None, None)]
#[case::atom_list(b"    1.2345    2.3456    3.4567 L   0  0  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }), None, None, Some(4), None, None, None, None, None)]
#[case::chemaxon_wildcard(b"    1.2345    2.3456    3.4567 AH  0  0  0  0  0  0  0  0  0  0  0  0",
       AtomSymbol::WildcardAtom(WildcardAtom::HeavyOrH), None, None, None, None, None, None, None, None)]
#[case::invalid_coordinates(b"   invalid   invalid   invalid C  -2  3  0  0  0  4  0  0  0  0  0  0",
        AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), None, None, None, None, None)]
fn test_extended_atom_input_ignore_positions(
    #[case] input: &[u8],
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
    #[case] chirality: Option<Chirality>,
    #[case] hydrogen_count: Option<u8>,
    #[case] stereo_care: Option<AtomStereoCare>,
    #[case] inversion_retention: Option<AtomInversionRetention>,
    #[case] exact_change: Option<AtomExactChange>,
) {
    let flags = CtabParseFlags::LENIENT | CtabParseFlags::IGNORE_POSITIONS;
    let result = extended_atom_input(flags).parse(Input::new(input));
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (atom, position) = result.unwrap();
    assert_eq!(
        atom.symbol, symbol,
        "{:?} has returned symbol {:?}, expected {:?}",
        input_str, atom.symbol, symbol,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{:?} has returned charge {:?}, expected {:?}",
        input_str, atom.charge, charge,
    );
    assert_eq!(
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence,
    );
    assert_eq!(
        atom.chirality, chirality,
        "{:?} has returned chirality {:?}, expected {:?}",
        input_str, atom.chirality, chirality,
    );
    assert_eq!(
        atom.implicit_hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.implicit_hydrogens, hydrogen_count,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{:?} has returned stereo_care {:?}, expected {:?}",
        input_str, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.inversion_retention, inversion_retention,
        "{:?} has returned inversion_retention {:?}, expected {:?}",
        input_str, atom.inversion_retention, inversion_retention,
    );
    assert_eq!(
        atom.exact_change, exact_change,
        "{:?} has returned exact_change {:?}, expected {:?}",
        input_str, atom.exact_change, exact_change,
    );
    assert_eq!(
        position.x, 0.0,
        "{:?} has returned x {:?}, expected 0.0",
        input_str, position.x,
    );
    assert_eq!(
        position.y, 0.0,
        "{:?} has returned y {:?}, expected 0.0",
        input_str, position.y,
    );
    assert_eq!(
        position.z, 0.0,
        "{:?} has returned z {:?}, expected 0.0",
        input_str, position.z,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::invalid_coordinates(b"   invalid   invalid   invalid C  -2  3  0  0  0  4  0  0  0  0  0  0", 0)]
fn test_extended_atom_input_ignore_positions_strict_error(
    #[case] input: &[u8],
    #[case] column: u32,
) {
    let error = extended_atom_input(
        CtabParseFlags::EXTENDED & CtabParseFlags::STRICT | CtabParseFlags::IGNORE_POSITIONS,
    )
    .parse(Input::new(input))
    .unwrap_err()
    .into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::element_h(b"H  ", AtomSymbol::Element(Element::H))]
#[case::element_c(b"C  ", AtomSymbol::Element(Element::C))]
#[case::element_c_lowercase(b"c  ", AtomSymbol::Element(Element::C))]
#[case::element_pos2(b" C ", AtomSymbol::Element(Element::C))]
#[case::element_pos3(b"  C", AtomSymbol::Element(Element::C))]
#[case::element_cu(b"Cu ", AtomSymbol::Element(Element::Cu))]
#[case::element_cu_pos2(b" Cu", AtomSymbol::Element(Element::Cu))]
#[case::element_cu_lowercase(b"cu ", AtomSymbol::Element(Element::Cu))]
#[case::element_cu_uppercase(b"CU ", AtomSymbol::Element(Element::Cu))]
#[case::named_isotope_d(b"D  ", AtomSymbol::NamedIsotope(NamedIsotope::D))]
#[case::named_isotope_d_lowercase(b"d  ", AtomSymbol::NamedIsotope(NamedIsotope::D))]
#[case::named_isotope_t(b"T  ", AtomSymbol::NamedIsotope(NamedIsotope::T))]
#[case::element_h_one_character(b"H", AtomSymbol::Element(Element::H))]
#[case::element_h_two_characters(b"H ", AtomSymbol::Element(Element::H))]
#[case::element_hg_two_characters(b"Hg", AtomSymbol::Element(Element::Hg))]
fn test_parse_atom_symbol(#[case] input: &[u8], #[case] expected: AtomSymbol) {
    let result = parse_atom_symbol(input, true, 0);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let symbol = result.unwrap();
    assert_eq!(
        symbol, expected,
        "{:?} has returned symbol {:?}, expected {:?}",
        input_str, symbol, expected
    );
}

#[rstest]
#[case::empty(b"")]
#[case::blank(b"   ")]
#[case::element_invalid_1(b"Xx ")]
#[case::element_invalid_2(b"LQ ")]
#[case::wildcard_atom_a(b"A  ")]
#[case::chemaxon_wildcard_atom(b"QH ")]
#[case::atom_list(b"L  ")]
#[case::lone_pair(b"LP ")]
#[case::rgroup(b"R1 ")]
#[case::pseudoatom_ala(b"Ala")]
#[case::pseudoatom_unicode(b"\xCE\xB1 ")]
fn test_parse_atom_symbol_error(#[case] input: &[u8]) {
    assert_eq!(
        parse_atom_symbol(input, true, 0),
        Err(ErrMode::Backtrack(InputError { column: 0 }))
    );
}

#[rstest]
#[case::named_isotope_d(b"D  ")]
#[case::named_isotope_d_lowercase(b"d  ")]
#[case::named_isotope_t(b"T  ")]
fn test_parse_atom_symbol_strict_error(#[case] input: &[u8]) {
    assert_eq!(
        parse_atom_symbol(input, false, 0),
        Err(ErrMode::Backtrack(InputError { column: 0 }))
    );
}

#[rstest]
#[case::element_h(b"H  ", AtomSymbol::Element(Element::H))]
#[case::element_c(b"C  ", AtomSymbol::Element(Element::C))]
#[case::element_c_lowercase(b"c  ", AtomSymbol::Element(Element::C))]
#[case::element_pos2(b" C ", AtomSymbol::Element(Element::C))]
#[case::element_pos3(b"  C", AtomSymbol::Element(Element::C))]
#[case::element_cu(b"Cu ", AtomSymbol::Element(Element::Cu))]
#[case::element_cu_pos2(b" Cu", AtomSymbol::Element(Element::Cu))]
#[case::element_cu_lowercase(b"cu ", AtomSymbol::Element(Element::Cu))]
#[case::element_cu_uppercase(b"CU ", AtomSymbol::Element(Element::Cu))]
#[case::named_isotope_d(b"D  ", AtomSymbol::NamedIsotope(NamedIsotope::D))]
#[case::named_isotope_d_lowercase(b"d  ", AtomSymbol::NamedIsotope(NamedIsotope::D))]
#[case::named_isotope_t(b"T  ", AtomSymbol::NamedIsotope(NamedIsotope::T))]
#[case::element_h_one_character(b"H", AtomSymbol::Element(Element::H))]
#[case::element_h_two_characters(b"H ", AtomSymbol::Element(Element::H))]
#[case::element_hg_two_characters(b"Hg", AtomSymbol::Element(Element::Hg))]
#[case::wildcard_atom_a(b"A  ", AtomSymbol::WildcardAtom(WildcardAtom::Heavy))]
#[case::wildcard_atom_q(b"Q  ", AtomSymbol::WildcardAtom(WildcardAtom::Heteroatom))]
#[case::wildcard_atom_star(b"*  ", AtomSymbol::WildcardAtom(WildcardAtom::Any))]
#[case::atom_list(b"L  ", AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }))]
#[case::lone_pair(b"LP ", AtomSymbol::LonePair)]
#[case::rgroup(b"R  ", AtomSymbol::RGroup(RGroup::new(None)))]
#[case::rgroup_unlabeled(b"R# ", AtomSymbol::RGroup(RGroup::new(None)))]
#[case::rgroup_r1(b"R1 ", AtomSymbol::RGroup(RGroup::new(Some(1))))]
#[case::rgroup_r3(b"R3 ", AtomSymbol::RGroup(RGroup::new(Some(3))))]
fn test_parse_extended_atom_symbol(#[case] input: &[u8], #[case] expected: AtomSymbol) {
    let result = parse_extended_atom_symbol(input, true, true, false, true, true, false, 0);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let symbol = result.unwrap();
    assert_eq!(
        symbol, expected,
        "{:?} has returned symbol {:?}, expected {:?}",
        input_str, symbol, expected
    );
}

#[rstest]
#[case::empty(b"")]
#[case::blank(b"   ")]
#[case::element_invalid_1(b"Xx ")]
#[case::element_invalid_2(b"LQ ")]
#[case::chemaxon_wildcard_atom(b"QH ")]
#[case::pseudoatom_ala(b"Ala")]
#[case::pseudoatom_unicode(b"\xCE\xB1 ")]
fn test_extended_atom_symbol_error(#[case] input: &[u8]) {
    assert_eq!(
        parse_extended_atom_symbol(input, true, true, false, true, true, false, 0),
        Err(ErrMode::Backtrack(InputError { column: 0 }))
    );
}

#[rstest]
#[case::wildcard_atom_a(b"A  ", AtomSymbol::WildcardAtom(WildcardAtom::Heavy))]
#[case::wildcard_atom_q(b"Q  ", AtomSymbol::WildcardAtom(WildcardAtom::Heteroatom))]
#[case::wildcard_atom_star(b"*  ", AtomSymbol::WildcardAtom(WildcardAtom::Any))]
#[case::atom_list(b"L  ", AtomSymbol::AtomList(AtomList::empty()))]
#[case::lone_pair(b"LP ", AtomSymbol::LonePair)]
#[case::rgroup(b"R  ", AtomSymbol::RGroup(RGroup::new(None)))]
#[case::rgroup_unlabeled(b"R# ", AtomSymbol::RGroup(RGroup::new(None)))]
#[case::rgroup_r1(b"R1 ", AtomSymbol::RGroup(RGroup::new(Some(1))))]
#[case::rgroup_r3(b"R3 ", AtomSymbol::RGroup(RGroup::new(Some(3))))]
fn test_parse_extended_atom_symbol_strict(#[case] input: &[u8], #[case] expected: AtomSymbol) {
    let result = parse_extended_atom_symbol(input, false, true, false, true, true, false, 0);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let symbol = result.unwrap();
    assert_eq!(
        symbol, expected,
        "{:?} has returned symbol {:?}, expected {:?}",
        input_str, symbol, expected
    );
}

#[rstest]
#[case::named_isotope_d(b"D  ")]
#[case::named_isotope_d_lowercase(b"d  ")]
#[case::named_isotope_t(b"T  ")]
#[case::pseudoatom_ala(b"Ala")]
#[case::pseudoatom_unicode(b"\xCE\xB1 ")]
fn test_parse_extended_atom_symbol_strict_error(#[case] input: &[u8]) {
    assert_eq!(
        parse_extended_atom_symbol(input, false, true, false, true, true, false, 0),
        Err(ErrMode::Backtrack(InputError { column: 0 }))
    );
}

#[rustfmt::skip]
#[rstest]
#[case::chemaxon_wildcard_ah(b"AH ", AtomSymbol::WildcardAtom(WildcardAtom::HeavyOrH))]
#[case::chemaxon_wildcard_qh(b"QH ", AtomSymbol::WildcardAtom(WildcardAtom::HeteroatomOrH))]
#[case::chemaxon_wildcard_xh(b"XH ", AtomSymbol::WildcardAtom(WildcardAtom::HalogenOrH))]
#[case::chemaxon_wildcard_mh(b"MH ", AtomSymbol::WildcardAtom(WildcardAtom::MetalOrH))]
fn test_extended_atom_symbol_lenient(
    #[case] input: &[u8],
    #[case] expected: AtomSymbol,
) {
    let result = parse_extended_atom_symbol(input, true, true, true, true, true, false, 0);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let symbol = result.unwrap();
    assert_eq!(symbol, expected, "{:?} has returned symbol {:?}, expected {:?}", input_str, symbol, expected);
}

#[rstest]
#[case::empty(b"")]
#[case::blank(b"   ")]
#[case::element_invalid_1(b"Xx ")]
#[case::element_invalid_2(b"LQ ")]
#[case::pseudoatom_ala(b"Ala")]
#[case::pseudoatom_unicode(b"\xCE\xB1 ")]
fn test_extended_atom_symbol_lenient_error(#[case] input: &[u8]) {
    assert_eq!(
        parse_extended_atom_symbol(input, true, true, true, true, true, false, 0),
        Err(ErrMode::Backtrack(InputError { column: 0 }))
    );
}

#[rstest]
#[case::pseudoatom_ala(b"Ala", AtomSymbol::Pseudoatom("Ala".to_string()))]
fn test_extended_atom_symbol_pseudoatoms(#[case] input: &[u8], #[case] expected: AtomSymbol) {
    let result = parse_extended_atom_symbol(input, false, false, false, false, false, true, 0);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let symbol = result.unwrap();
    assert_eq!(
        symbol, expected,
        "{:?} has returned symbol {:?}, expected {:?}",
        input_str, symbol, expected
    );
}

#[rstest]
#[case::reserved_named_isotope(b"D  ")]
#[case::reserved_wildcard_atom(b"A  ")]
#[case::reserved_chemaxon_wildcard_atom(b"QH ")]
#[case::reserved_atom_list(b"L  ")]
#[case::reserved_lone_pair(b"LP ")]
#[case::reserved_rgroup(b"R  ")]
#[case::reserved_rgroup_unlabeled(b"R# ")]
#[case::reserved_rgroup_r1(b"R1 ")]
#[case::pseudoatom_unicode(b"\xCE\xB1 ")]
fn test_extended_atom_symbol_pseudoatoms_error(#[case] input: &[u8]) {
    assert_eq!(
        parse_extended_atom_symbol(input, false, false, false, false, false, true, 0),
        Err(ErrMode::Backtrack(InputError { column: 0 }))
    );
}
