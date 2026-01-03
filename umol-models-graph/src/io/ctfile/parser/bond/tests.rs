use bstr::ByteSlice;
use nom::error::ErrorKind as NomErrorKind;
use nom::{Err, Parser};
use pretty_assertions::assert_eq;
use rstest::*;

use super::*;
use crate::io::ctfile::config::CtabParseFlags;
use crate::table_ir::{BondDirection, BondOrder, BondReactingCenter, BondStereo, BondTopology};

#[rustfmt::skip]
#[rstest]
#[case::len_21_single(b"  2  5  1  0  0  0  0", 1, 4, BondOrder::Single, None, None)]
#[case::len_21_single_wedge(b"  1  2  1  1  0  0  0", 0, 1, BondOrder::Single, None, Some(BondDirection::Up))]
#[case::len_21_double(b"  2  5  2  0  0  0  0", 1, 4, BondOrder::Double, None, None)]
#[case::len_21_double_cis(b"  2  5  2  1  0  0  0", 1, 4, BondOrder::Double, Some(BondStereo::Cis), None)]
#[case::len_21_double_trans(b"  2  5  2  6  0  0  0", 1, 4, BondOrder::Double, Some(BondStereo::Trans), None)]
#[case::len_21_double_either(b"  2  5  2  4  0  0  0", 1, 4, BondOrder::Double, Some(BondStereo::Either), None)]
#[case::len_21_triple(b"  2  5  3  0  0  0  0", 1, 4, BondOrder::Triple, None, None)]
#[case::len_21_triple_ignored_stereo(b"  2  5  3  1  0  0  0", 1, 4, BondOrder::Triple, None, None)]
#[case::len_21_triple_empty_fields(b"  2  5  3  1         ", 1, 4, BondOrder::Triple, None, None)]
#[case::len_21_aromatic(b"  2  5  4  0  0  0  0", 1, 4, BondOrder::Aromatic, None, None)]
#[case::len_18_single(b"  1  2  1  0  0  0", 0, 1, BondOrder::Single, None, None)]
#[case::len_12_single_wedge(b"  1  3  1  1", 0, 2, BondOrder::Single, None, Some(BondDirection::Up))]
#[case::len_12_double(b"  2  5  2  0", 1, 4, BondOrder::Double, None, None)]
#[case::len_12_double_cis(b"  1  3  2  1", 0, 2, BondOrder::Double, Some(BondStereo::Cis), None)]
#[case::len_12_single_dash(b"  2  4  1  6", 1, 3, BondOrder::Single, None, Some(BondDirection::Down))]
#[case::len_12_triple_empty_fields(b"  2  5  3   ", 1, 4, BondOrder::Triple, None, None)]
#[case::non_strict_padding(b"  1  2  1  1  0  0XXX", 0, 1, BondOrder::Single, None, Some(BondDirection::Up))]
fn test_bond_input12(
    #[case] input: &[u8],
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondOrder,
    #[case] stereo: Option<BondStereo>,
    #[case] dir: Option<BondDirection>,
) {
    let result = bond_input12(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(a1, atom1, "{:?} has returned atom1 {:?}, expected {:?}", input_str, a1, atom1);
    assert_eq!(a2, atom2, "{:?} has returned atom2 {:?}, expected {:?}", input_str, a2, atom2);
    assert_eq!(bond.order, bond_type, "{:?} has returned bond type {:?}, expected {:?}", input_str, bond.order, bond_type);
    assert_eq!(bond.stereo, stereo, "{:?} has returned stereo {:?}, expected {:?}", input_str, bond.stereo, stereo);
    assert_eq!(bond.direction, dir, "{:?} has returned dir {:?}, expected {:?}", input_str, bond.direction, dir);
}

#[rustfmt::skip]
#[rstest]
#[case::len_21_extended_range_quadruple(b"  2  5  9  0  0  0  0", NomErrorKind::MapRes)]
#[case::len_21_extended_range_zero(b"  2  5  0  0  0  0  0", NomErrorKind::MapRes)]
#[case::len_12_line_too_short(b"  1  2  1 1", NomErrorKind::Eof)]
#[case::len_21_non_numeric_atom(b"  A  2  1  1  0  0  0", NomErrorKind::Digit)]
#[case::len_21_non_numeric_type(b"  1  2  A  1  0  0  0", NomErrorKind::Digit)]
#[case::len_12_non_numeric_stereo(b"  1  2  1  A", NomErrorKind::Digit)]
fn test_bond_input12_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = bond_input12(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_strict_padding(b"  1  2  1  1  0  0XXX", NomErrorKind::Verify)]
fn test_bond_input12_strict_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = bond_input12(input, CtabParseFlags::BASIC & CtabParseFlags::STRICT);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_21_quadruple(b"  2  5  9  0  0  0  0", 1, 4, BondOrder::Quadruple, None, None)]
#[case::len_21_zero(b"  2  5  0  0  0  0  0", 1, 4, BondOrder::Zero, None, None)]
fn test_bond_input12_lenient(
    #[case] input: &[u8],
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondOrder,
    #[case] stereo: Option<BondStereo>,
    #[case] dir: Option<BondDirection>,
) {
    let result = bond_input12(input, CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} has non-empty remaining",
        input_str
    );
    assert_eq!(
        a1, atom1,
        "{:?} has returned atom1 {:?}, expected {:?}",
        input_str, a1, atom1
    );
    assert_eq!(
        a2, atom2,
        "{:?} has returned atom2 {:?}, expected {:?}",
        input_str, a2, atom2
    );
    assert_eq!(
        bond.order, bond_type,
        "{:?} has returned bond type {:?}, expected {:?}",
        input_str, bond.order, bond_type
    );
    assert_eq!(
        bond.stereo, stereo,
        "{:?} has returned stereo {:?}, expected {:?}",
        input_str, bond.stereo, stereo
    );
    assert_eq!(
        bond.direction, dir,
        "{:?} has returned dir {:?}, expected {:?}",
        input_str, bond.direction, dir
    );
}

#[rstest]
#[case::single(b"  1  2  1", 0, 1, BondOrder::Single)]
#[case::double(b"  2  5  2", 1, 4, BondOrder::Double)]
#[case::triple(b"  2  5  3", 1, 4, BondOrder::Triple)]
#[case::aromatic(b"  2  5  4", 1, 4, BondOrder::Aromatic)]
fn test_bond_input9(
    #[case] input: &[u8],
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondOrder,
) {
    let input_str = input.to_str_lossy();
    let result = bond_input9(input, CtabParseFlags::BASIC);
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} has non-empty remaining",
        input_str
    );
    assert_eq!(
        a1, atom1,
        "{:?} has returned atom1 {:?}, expected {:?}",
        input_str, a1, atom1
    );
    assert_eq!(
        a2, atom2,
        "{:?} has returned atom2 {:?}, expected {:?}",
        input_str, a2, atom2
    );
    assert_eq!(
        bond.order, bond_type,
        "{:?} has returned bond type {:?}, expected {:?}",
        input_str, bond.order, bond_type
    );
    assert_eq!(
        bond.stereo, None,
        "{:?} has returned stereo {:?}, expected {:?}",
        input_str, bond.stereo, None as Option<BondStereo>
    );
}

#[rstest]
#[case::line_too_short(b"  1  2", NomErrorKind::MapRes)]
#[case::non_numeric_atom_2(b"  1  A  1", NomErrorKind::Digit)]
#[case::extended_range_quadruple(b"  1  2  9", NomErrorKind::MapRes)]
#[case::extended_range_zero(b"  1  2  0", NomErrorKind::MapRes)]
fn test_bond_input9_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let input_str = input.to_str_lossy();
    let result = bond_input9(input, CtabParseFlags::BASIC);
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::quadruple(b"  1  2  9", 0, 1, BondOrder::Quadruple)]
#[case::zero(b"  2  5  0", 1, 4, BondOrder::Zero)]
fn test_bond_input9_extended(
    #[case] input: &[u8],
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondOrder,
) {
    let input_str = input.to_str_lossy();
    let result = bond_input9(input, CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT);
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} has non-empty remaining",
        input_str
    );
    assert_eq!(
        a1, atom1,
        "{:?} has returned atom1 {:?}, expected {:?}",
        input_str, a1, atom1
    );
    assert_eq!(
        a2, atom2,
        "{:?} has returned atom2 {:?}, expected {:?}",
        input_str, a2, atom2
    );
    assert_eq!(
        bond.order, bond_type,
        "{:?} has returned bond type {:?}, expected {:?}",
        input_str, bond.order, bond_type
    );
    assert_eq!(
        bond.stereo, None,
        "{:?} has returned stereo {:?}, expected {:?}",
        input_str, bond.stereo, None as Option<BondStereo>
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_9(b"  1  2  1", 0, 1, BondOrder::Single, None, None)]
#[case::len_10_padded(b"  1  2  1 ", 0, 1, BondOrder::Single, None, None)]
#[case::len_12(b"  1  3  2  1", 0, 2, BondOrder::Double, Some(BondStereo::Cis), None)]
#[case::len_13_padded(b"  1  3  1  6 ", 0, 2, BondOrder::Single, None, Some(BondDirection::Down))]
#[case::len_18(b"  1  2  1  0  0  0", 0, 1, BondOrder::Single, None, None)]
#[case::len_to_21(b"  2  5  2  1  0  0  0", 1, 4, BondOrder::Double, Some(BondStereo::Cis), None)]
fn test_bond_input(
    #[case] input: &[u8],
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondOrder,
    #[case] stereo: Option<BondStereo>,
    #[case] dir: Option<BondDirection>,
) {
    let input_str = input.to_str_lossy();
    let mut parser = bond_input(CtabParseFlags::BASIC);
    let result = parser.parse(input);
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(a1, atom1, "{:?} has returned atom1 {:?}, expected {:?}", input_str, a1, atom1);
    assert_eq!(a2, atom2, "{:?} has returned atom2 {:?}, expected {:?}", input_str, a2, atom2);
    assert_eq!(bond.order, bond_type, "{:?} has returned bond type {:?}, expected {:?}", input_str, bond.order, bond_type);
    assert_eq!(bond.stereo, stereo, "{:?} has returned stereo {:?}, expected {:?}", input_str, bond.stereo, stereo);
    assert_eq!(bond.direction, dir, "{:?} has returned dir {:?}, expected {:?}", input_str, bond.direction, dir);
}

#[rstest]
#[case::line_too_short(b"  1  2 ", NomErrorKind::Eof)]
#[case::extended_range_quadruple(b"  1  2  9", NomErrorKind::MapRes)]
#[case::extended_range_zero(b"  1  2  0", NomErrorKind::MapRes)]
fn test_bond_input_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let input_str = input.to_str_lossy();
    let mut parser = bond_input(CtabParseFlags::BASIC);
    let result = parser.parse(input);
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::non_strict_padding(b"  1  2  1  1  0  0XXX", NomErrorKind::Verify)]
fn test_bond_input_strict_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let input_str = input.to_str_lossy();
    let mut parser = bond_input(CtabParseFlags::STRICT);
    let result = parser.parse(input);
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_9_zero_bond(b"  1  2  0", 0, 1, BondOrder::Zero)]
#[case::len_9_quadruple_bond(b"  1  2  9", 0, 1, BondOrder::Quadruple)]
fn test_bond_input_lenient(
    #[case] input: &[u8],
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondOrder,
) {
    let input_str = input.to_str_lossy();
    let result = bond_input(CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT).parse(input);
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} should consume all input", input_str);
    assert_eq!(a1, atom1, "{:?} has returned atom1 {:?}, expected {:?}", input_str, a1, atom1);
    assert_eq!(a2, atom2, "{:?} has returned atom2 {:?}, expected {:?}", input_str, a2, atom2);
    assert_eq!(bond.order, bond_type, "{:?} has returned bond type {:?}, expected {:?}", input_str, bond.order, bond_type);
}

#[rustfmt::skip]
#[rstest]
#[case::line_too_short(b"  1  2 ", NomErrorKind::Eof)]
fn test_bond_input_lenient_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let input_str = input.to_str_lossy();
    let result = bond_input(CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT).parse(input);
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::len_11(b"  1  2  1 1")]
fn test_bond_input_partial_fields(#[case] input: &[u8]) {
    let input_str = input.to_str_lossy();
    let mut parser = bond_input(CtabParseFlags::BASIC);
    let result = parser.parse(input);
    assert!(result.is_err(), "{:?} should have failed", input_str);
}

#[rstest]
#[case::len_9_padded(b"  1  2  1\n")]
#[case::len_12_padded(b"  1  3  1  1  ")]
fn test_bond_input_whitespace_padded(#[case] input: &[u8]) {
    let input_str = input.to_str_lossy();
    let mut parser = bond_input(CtabParseFlags::BASIC);
    let trimmed_input = input.trim_ascii_end();
    let result = parser.parse(trimmed_input);
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (_a1, _a2, _bond)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} has non-empty remaining",
        input_str
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_9(b"  1  2  1", 0, 1, BondOrder::Single, None, None, None, None)]
#[case::len_9_zero_bond(b"  1  2  0", 0, 1, BondOrder::Zero, None, None, None, None)]
#[case::len_9_quadruple_bond(b"  1  2  9", 0, 1, BondOrder::Quadruple, None, None, None, None)]
#[case::len_12_double_either(b"  1  2  2  3", 0, 1, BondOrder::Double, Some(BondStereo::Either), None, None, None)]
#[case::len_13_single_dash_padded(b"  1  2  1  6  ", 0, 1, BondOrder::Single, None, Some(BondDirection::Down), None, None)]
#[case::len_18_any_bond_ring(b"  1  2  8  0     1", 0, 1, BondOrder::Any, None, None, Some(BondTopology::Ring), None)]
#[case::len_18(b"  1  2  1  0  0  0", 0, 1, BondOrder::Single, None, Some(BondDirection::NotStereo), Some(BondTopology::Either), None)]
#[case::len_21_full_chain_center(b"  1  2  1  0     2  1", 0, 1, BondOrder::Single, None,
       Some(BondDirection::NotStereo), Some(BondTopology::Chain), Some(BondReactingCenter::CENTER))]
#[case::len_21_full_not_center(b"  1  2  1  0     2 -1", 0, 1, BondOrder::Single, None,
       Some(BondDirection::NotStereo), Some(BondTopology::Chain), Some(BondReactingCenter::NOT_CENTER))]
#[case::len_21_only_mandatory_fields(b"  1  2  1            ", 0, 1, BondOrder::Single, None,
       Some(BondDirection::NotStereo), None, None)]
#[case::non_strict_padding(b"  1  2  8  0XXX  1", 0, 1, BondOrder::Any, None, None, Some(BondTopology::Ring), None)]
fn test_extended_bond_input(
    #[case] input: &[u8],
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondOrder,
    #[case] stereo: Option<BondStereo>,
    #[case] dir: Option<BondDirection>,
    #[case] topology: Option<BondTopology>,
    #[case] reacting_center: Option<BondReactingCenter>,
) {
    let input_str = input.to_str_lossy();
    let mut parser = extended_bond_input(CtabParseFlags::LENIENT);
    let result = parser.parse(input);
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(a1, atom1, "{:?} has returned atom1 {:?}, expected {:?}", input_str, a1, atom1);
    assert_eq!(a2, atom2, "{:?} has returned atom2 {:?}, expected {:?}", input_str, a2, atom2);
    assert_eq!(bond.order, bond_type, "{:?} has returned bond type {:?}, expected {:?}", input_str, bond.order, bond_type);
    assert_eq!(bond.stereo, stereo, "{:?} has returned stereo {:?}, expected {:?}", input_str, bond.stereo, stereo);
    assert_eq!(bond.direction, dir, "{:?} has returned dir {:?}, expected {:?}", input_str, bond.direction, dir);
    assert_eq!(bond.topology, topology, "{:?} has returned topology {:?}, expected {:?}", input_str, bond.topology, topology);
    assert_eq!(bond.reacting_center, reacting_center, "{:?} has returned reacting_center {:?}, expected {:?}", input_str, bond.reacting_center, reacting_center);
}

#[rstest]
#[case::line_too_short(b"  1  2", NomErrorKind::MapRes)]
#[case::bond_type_above_range(b"  1  2  9", NomErrorKind::MapRes)]
#[case::bond_type_below_range(b"  1  2  0", NomErrorKind::MapRes)]
#[case::invalid_topology(b"  2  5  2  0  0  4  0", NomErrorKind::MapRes)]
#[case::non_strict_padding(b"  1  2  8  0XXX  1", NomErrorKind::Verify)]
fn test_extended_bond_input_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let input_str = input.to_str_lossy();
    let mut parser = extended_bond_input(CtabParseFlags::STRICT);
    let result = parser.parse(input);
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}
