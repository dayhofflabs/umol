use super::*;
use crate::io::ctab::bond::{BondDir, BondReactingCenter, BondStereo, BondTopology, BondType};
use nom::combinator::all_consuming;
use nom::{error, Err, Parser};
use pretty_assertions::assert_eq;
use rstest::*;

#[rustfmt::skip]
#[rstest]
#[case(b"  1  2  1  1  0  0  0", "len 21 single wedge", 0, 1, BondType::Single, None, Some(BondDir::Wedge))]
#[case(b"  2  5  1  0  0  0  0", "len 21 single unknown", 1, 4, BondType::Single, None, None)]
#[case(b"  2  5  2  1  0  0  0", "len 21 double cis", 1, 4, BondType::Double, Some(BondStereo::Cis), None)]
#[case(b"  2  5  2  6  0  0  0", "len 21 double trans", 1, 4, BondType::Double, Some(BondStereo::Trans), None)]
#[case(b"  2  5  2  4  0  0  0", "len 21 double either", 1, 4, BondType::Double, Some(BondStereo::Either), None)]
#[case(b"  2  5  2  0  0  0  0", "len 21 double none", 1, 4, BondType::Double, None, None)]
#[case(b"  2  5  3  0  0  0  0", "len 21 triple", 1, 4, BondType::Triple, None, None)]
#[case(b"  2  5  4  0  0  0  0", "len 21 aromatic", 1, 4, BondType::Aromatic, None, None)]
#[case(b"  2  5  3  1  0  0  0", "len 21 triple ignored stereo", 1, 4, BondType::Triple, None, None)]
#[case(b"  2  5  3  1         ", "len 21 triple empty fields", 1, 4, BondType::Triple, None, None)]
#[case(b"  1  2  1  0  0  0", "len 18 single", 0, 1, BondType::Single, None, None)]
#[case(b"  1  3  1  1", "len 12 single wedge", 0, 2, BondType::Single, None, Some(BondDir::Wedge))]
#[case(b"  1  3  2  1", "len 12 double cis", 0, 2, BondType::Double, Some(BondStereo::Cis), None)]
#[case(b"  2  5  2  0", "len 12 double none", 1, 4, BondType::Double, None, None)]
#[case(b"  2  4  1  6", "len 12 single dash", 1, 3, BondType::Single, None, Some(BondDir::Dash))]
#[case(b"  2  5  3   ", "len 12 triple empty fields", 1, 4, BondType::Triple, None, None)]
fn test_bond_input_standard12(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondType,
    #[case] stereo: Option<BondStereo>,
    #[case] dir: Option<BondDir>,
) {
    let result = bond_input_standard12(input);
    assert!(result.is_ok(), "{} should have succeeded", desc,);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc,);
    assert_eq!(
        a1, atom1,
        "{} has returned atom1 {:?}, expected {:?}",
        desc, a1, atom1
    );
    assert_eq!(
        a2, atom2,
        "{} has returned atom2 {:?}, expected {:?}",
        desc, a2, atom2
    );
    assert_eq!(
        bond.bond_type, bond_type,
        "{} has returned bond type {:?}, expected {:?}",
        desc, bond.bond_type, bond_type,
    );
    assert_eq!(
        bond.stereo, stereo,
        "{} has returned stereo {:?}, expected {:?}",
        desc, bond.stereo, stereo
    );
    assert_eq!(
        bond.dir, dir,
        "{} has returned dir {:?}, expected {:?}",
        desc, bond.dir, dir
    );
}

#[rstest]
#[case(b"  1  2  1  A", "len 12 non-numeric stereo", error::ErrorKind::Digit)]
#[case(b"  1  2  1 1", "len 12 line too short", error::ErrorKind::Eof)]
#[case(b"  A  2  1  1  0  0  0", "len 21 non-numeric atom", error::ErrorKind::Digit)]
#[case(b"  1  2  A  1  0  0  0", "len 21 non-numeric type", error::ErrorKind::Digit)]
fn test_bond_input_standard12_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = bond_input_standard12(input);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{} should have failed with error kind {:?}, got {:?}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case(b"  1  2  1", "single", 0, 1, BondType::Single)]
#[case(b"  2  5  2", "double", 1, 4, BondType::Double)]
#[case(b"  2  5  3", "triple", 1, 4, BondType::Triple)]
#[case(b"  2  5  4", "aromatic", 1, 4, BondType::Aromatic)]
fn test_bond_input_standard9(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondType,
) {
    let result = bond_input_standard9(input);
    assert!(result.is_ok(), "{} should have succeeded", desc,);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc,);
    assert_eq!(
        a1, atom1,
        "{} has returned atom1 {:?}, expected {:?}",
        desc, a1, atom1
    );
    assert_eq!(
        a2, atom2,
        "{} has returned atom2 {:?}, expected {:?}",
        desc, a2, atom2
    );
    assert_eq!(
        bond.bond_type, bond_type,
        "{} has returned bond type {:?}, expected {:?}",
        desc, bond.bond_type, bond_type,
    );
    assert_eq!(
        bond.stereo, None,
        "{} has returned stereo {:?}, expected {:?}",
        desc, bond.stereo, None as Option<BondStereo>
    );
}

#[rstest]
#[case(b"  1  2", "Line too short", error::ErrorKind::MapRes)]
#[case(b"  1  A  1", "Non-numeric atom 2", error::ErrorKind::Digit)]
fn test_bond_input_standard9_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = bond_input_standard9(input);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{} should have failed with error kind {:?}, got {:?}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"  1  2  1", "len 9", 0, 1, BondType::Single, None, None)]
#[case(b"  1  2  1 ", "len 10 padded", 0, 1, BondType::Single, None, None)]
#[case(b"  1  3  2  1", "len 12", 0, 2, BondType::Double, Some(BondStereo::Cis), None)]
#[case(b"  1  3  1  6 ", "len 13 padded", 0, 2, BondType::Single, None, Some(BondDir::Dash))]
#[case(b"  1  2  1  0  0  0", "len 18", 0, 1, BondType::Single, None, None)]
#[case(b"  2  5  2  1  0  0  0", "len to 21", 1, 4, BondType::Double, Some(BondStereo::Cis), None)]
fn test_bond_input_standard(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondType,
    #[case] stereo: Option<BondStereo>,
    #[case] dir: Option<BondDir>,
) {
    let mut parser = bond_input_standard();
    let result = parser.parse(input);
    assert!(result.is_ok(), "{} should have succeeded", desc,);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc,);
    assert_eq!(
        a1, atom1,
        "{} has returned atom1 {:?}, expected {:?}",
        desc, a1, atom1
    );
    assert_eq!(
        a2, atom2,
        "{} has returned atom2 {:?}, expected {:?}",
        desc, a2, atom2
    );
    assert_eq!(
        bond.bond_type, bond_type,
        "{} has returned bond type {:?}, expected {:?}",
        desc, bond.bond_type, bond_type,
    );
    assert_eq!(
        bond.stereo, stereo,
        "{} has returned stereo {:?}, expected {:?}",
        desc, bond.stereo, stereo
    );
    assert_eq!(
        bond.dir, dir,
        "{} has returned dir {:?}, expected {:?}",
        desc, bond.dir, dir
    );
}

#[rstest]
#[case(b"  1  2 ", "Line too short", error::ErrorKind::Eof)]
#[case(b"  1  2  9", "Out of range type", error::ErrorKind::MapRes)]
fn test_bond_input_standard_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = bond_input_standard();
    let result = parser.parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{} should have failed with error kind {:?}, got {:?}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case(b"  1  2  1 1", "len 11")]
fn test_bond_input_standard_partial_fields(#[case] input: &[u8], #[case] desc: &str) {
    let mut parser = bond_input_standard();
    let result = parser.parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
}

#[rstest]
#[case(b"  1  2  1\n", "len 9 padded")]
#[case(b"  1  3  1  1  ", "len 12 padded")]
fn test_bond_input_standard_whitespace_padded(#[case] input: &[u8], #[case] desc: &str) {
    let mut parser = bond_input_standard();
    let trimmed_input = input.trim_ascii_end();
    let result = parser.parse(trimmed_input);
    assert!(result.is_ok(), "{} should have succeeded", desc,);
    let (remaining, (_a1, _a2, _bond)) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc,);
}

#[rustfmt::skip]
#[rstest]
#[case(b"  1  2  1", "len 9", 0, 1, BondType::Single, None, None, None, None)]
#[case(b"  1  2  2  3", "len 12 double either", 0, 1, BondType::Double, Some(BondStereo::Either), None, None, None)]
#[case(b"  1  2  1  6  ", "len 13 single dash padded", 0, 1, BondType::Single, None, Some(BondDir::Dash), None, None)]
#[case(b"  1  2  8  0     1", "len 18 any bond, ring", 0, 1, BondType::Any, None, None, Some(BondTopology::Ring), None)]
#[case(b"  1  2  1  0     2  1", "len 21 full, chain, center", 0, 1, BondType::Single, None,
       Some(BondDir::Either), Some(BondTopology::Chain), Some(BondReactingCenter::CENTER))]
#[case(b"  1  2  1  0     2 -1", "len 21 full, not center", 0, 1, BondType::Single, None,
       Some(BondDir::Either), Some(BondTopology::Chain), Some(BondReactingCenter::NOT_CENTER))]
#[case(b"  1  2  1            ", "len 21 only mandatory fields", 0, 1, BondType::Single, None,
       Some(BondDir::Either), Some(BondTopology::Either), Some(BondReactingCenter::UNMARKED))]
fn test_bond_input(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondType,
    #[case] stereo: Option<BondStereo>,
    #[case] dir: Option<BondDir>,
    #[case] topology: Option<BondTopology>,
    #[case] reacting_center: Option<BondReactingCenter>,
) {
    let mut parser = bond_input();
    let result = parser.parse(input);
    assert!(result.is_ok(), "{} should have succeeded", desc,);
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc,);
    assert_eq!(
        a1, atom1,
        "{} has returned atom1 {:?}, expected {:?}",
        desc, a1, atom1
    );
    assert_eq!(
        a2, atom2,
        "{} has returned atom2 {:?}, expected {:?}",
        desc, a2, atom2
    );
    assert_eq!(
        bond.bond_type, bond_type,
        "{} has returned bond type {:?}, expected {:?}",
        desc, bond.bond_type, bond_type,
    );
    assert_eq!(
        bond.stereo, stereo,
        "{} has returned stereo {:?}, expected {:?}",
        desc, bond.stereo, stereo
    );
    assert_eq!(
        bond.dir, dir,
        "{} has returned dir {:?}, expected {:?}",
        desc, bond.dir, dir
    );
    assert_eq!(
        bond.topology, topology,
        "{} has returned topology {:?}, expected {:?}",
        desc, bond.topology, topology
    );
    assert_eq!(
        bond.reacting_center, reacting_center,
        "{} has returned reacting_center {:?}, expected {:?}",
        desc, bond.reacting_center, reacting_center,
    );
}

#[rstest]
#[case(b"  1  2", "Line too short", error::ErrorKind::MapRes)]
#[case(b"  1  2  9", "Out of range type", error::ErrorKind::MapRes)]
#[case(b"  2  5  2  0  0  4  0", "Invalid topology", error::ErrorKind::MapRes)]
#[case(
    b"  1  2  1  0         a",
    "trailing non-whitespace",
    error::ErrorKind::Eof
)]
fn test_bond_input_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = all_consuming(bond_input());
    let result = parser.parse(input);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{} should have failed with error kind {:?}, got {:?}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}
