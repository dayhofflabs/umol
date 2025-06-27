use super::*;
use crate::bond::{BondDir, BondStereo, BondType};
use nom::{error::ErrorKind, Err, Parser};
use rstest::rstest;

#[rstest]
#[case(
    b"  1  2  1  1  0  0  0",
    "single wedge",
    0,
    1,
    BondType::Single,
    None,
    Some(BondDir::Wedge)
)]
#[case(
    b"  2  5  1  0  0  0  0",
    "single unknown",
    1,
    4,
    BondType::Single,
    None,
    None
)]
#[case(
    b"  2  5  2  1  0  0  0",
    "double cis",
    1,
    4,
    BondType::Double,
    Some(BondStereo::Cis),
    None
)]
#[case(
    b"  2  5  2  6  0  0  0",
    "double trans",
    1,
    4,
    BondType::Double,
    Some(BondStereo::Trans),
    None
)]
#[case(
    b"  2  5  2  4  0  0  0",
    "double either",
    1,
    4,
    BondType::Double,
    Some(BondStereo::Either),
    None
)]
#[case(
    b"  2  5  2  0  0  0  0",
    "double none",
    1,
    4,
    BondType::Double,
    None,
    None
)]
#[case(b"  2  5  3  0  0  0  0", "triple", 1, 4, BondType::Triple, None, None)]
#[case(
    b"  2  5  4  0  0  0  0",
    "aromatic",
    1,
    4,
    BondType::Aromatic,
    None,
    None
)]
#[case(
    b"  2  5  3  1  0  0  0",
    "triple ignored stereo",
    1,
    4,
    BondType::Triple,
    None,
    None
)]
#[case(
    b"  2  5  3  1         ",
    "triple empty fields",
    1,
    4,
    BondType::Triple,
    None,
    None
)]
fn test_bond_input_standard21(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] atom1: usize,
    #[case] atom2: usize,
    #[case] bond_type: BondType,
    #[case] stereo: Option<BondStereo>,
    #[case] dir: Option<BondDir>,
) {
    let result = bond_input_standard21(input);
    assert!(
        result.is_ok(),
        "Parser failed for case '{}': {:?}",
        desc,
        result
    );
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}'",
        desc
    );
    assert_eq!(a1, atom1, "Mismatched atom1 for '{}'", desc);
    assert_eq!(a2, atom2, "Mismatched atom2 for '{}'", desc);
    assert_eq!(
        bond.bond_type, bond_type,
        "Mismatched bond type for '{}'",
        desc
    );
    assert_eq!(bond.stereo, stereo, "Mismatched stereo for '{}'", desc);
    assert_eq!(bond.dir, dir, "Mismatched dir for '{}'", desc);
}

#[rstest]
#[case(b"  A  2  1  1  0  0  0", "Non-numeric atom", ErrorKind::Digit)]
#[case(b"  1  2  A  1  0  0  0", "Non-numeric type", ErrorKind::Digit)]
#[case(b"  1  2  1  1  0  0 ", "Line too short", ErrorKind::Eof)]
fn test_bond_input_standard21_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = bond_input_standard21(input);
    assert!(result.is_err(), "Parser should have failed for '{}'", desc);
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    } else {
        panic!("Expected nom::Err::Error for '{}', got {:?}", desc, result);
    }
}

#[rstest]
#[case(
    b"  1  3  1  1",
    "single wedge",
    0,
    2,
    BondType::Single,
    None,
    Some(BondDir::Wedge)
)]
#[case(
    b"  1  3  2  1",
    "double cis",
    0,
    2,
    BondType::Double,
    Some(BondStereo::Cis),
    None
)]
#[case(b"  2  5  2  0", "double none", 1, 4, BondType::Double, None, None)]
#[case(
    b"  2  4  1  6",
    "single dash",
    1,
    3,
    BondType::Single,
    None,
    Some(BondDir::Dash)
)]
#[case(
    b"  2  5  3   ",
    "triple empty fields",
    1,
    4,
    BondType::Triple,
    None,
    None
)]
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
    assert!(
        result.is_ok(),
        "Parser failed for case '{}': {:?}",
        desc,
        result
    );
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}'",
        desc
    );
    assert_eq!(a1, atom1, "Mismatched atom1 for '{}'", desc);
    assert_eq!(a2, atom2, "Mismatched atom2 for '{}'", desc);
    assert_eq!(
        bond.bond_type, bond_type,
        "Mismatched bond type for '{}'",
        desc
    );
    assert_eq!(bond.stereo, stereo, "Mismatched stereo for '{}'", desc);
    assert_eq!(bond.dir, dir, "Mismatched dir for '{}'", desc);
}

#[rstest]
#[case(b"  1  2  1  A", "non-numeric stereo", ErrorKind::Digit)]
#[case(b"  1  2  1 1", "Line too short", ErrorKind::Eof)]
fn test_bond_input_standard12_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = bond_input_standard12(input);
    assert!(result.is_err(), "Parser should have failed for '{}'", desc);
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    } else {
        panic!("Expected nom::Err::Error for '{}', got {:?}", desc, result);
    }
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
    assert!(
        result.is_ok(),
        "Parser failed for case '{}': {:?}",
        desc,
        result
    );
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}'",
        desc
    );
    assert_eq!(a1, atom1, "Mismatched atom1 for '{}'", desc);
    assert_eq!(a2, atom2, "Mismatched atom2 for '{}'", desc);
    assert_eq!(
        bond.bond_type, bond_type,
        "Mismatched bond type for '{}'",
        desc
    );
    assert_eq!(bond.stereo, None, "Stereo should be None for '{}'", desc);
}

#[rstest]
#[case(b"  1  2", "Line too short", ErrorKind::MapRes)]
#[case(b"  1  A  1", "Non-numeric atom 2", ErrorKind::Digit)]
fn test_bond_input_standard9_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = bond_input_standard9(input);
    assert!(result.is_err(), "Parser should have failed for '{}'", desc);
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    } else {
        panic!("Expected nom::Err::Error for '{}', got {:?}", desc, result);
    }
}
#[rstest]
#[case(b"  1  2  1", "len 9", 0, 1, BondType::Single, None, None)]
#[case(b"  1  2  1 ", "len 10 padded", 0, 1, BondType::Single, None, None)]
#[case(
    b"  1  3  2  1",
    "len 12",
    0,
    2,
    BondType::Double,
    Some(BondStereo::Cis),
    None
)]
#[case(
    b"  1  3  1  6 ",
    "len 13 padded",
    0,
    2,
    BondType::Single,
    None,
    Some(BondDir::Dash)
)]
#[case(
    b"  2  5  2  1  0  0  0",
    "len to 21",
    1,
    4,
    BondType::Double,
    Some(BondStereo::Cis),
    None
)]
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
    assert!(
        result.is_ok(),
        "Parser failed for case '{}': {:?}",
        desc,
        result
    );
    let (remaining, (a1, a2, bond)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}'",
        desc
    );
    assert_eq!(a1, atom1, "Mismatched atom1 for '{}'", desc);
    assert_eq!(a2, atom2, "Mismatched atom2 for '{}'", desc);
    assert_eq!(
        bond.bond_type, bond_type,
        "Mismatched bond type for '{}'",
        desc
    );
    assert_eq!(bond.stereo, stereo, "Mismatched stereo for '{}'", desc);
    assert_eq!(bond.dir, dir, "Mismatched dir for '{}'", desc);
}

#[rstest]
#[case(b"  1  2 ", "Line too short", ErrorKind::Eof)]
#[case(b"  1  2  9", "Out of range type", ErrorKind::MapRes)]
fn test_bond_input_standard_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let mut parser = bond_input_standard();
    let result = parser.parse(input);
    assert!(result.is_err(), "Parser should have failed for '{}'", desc);
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    } else {
        panic!("Expected nom::Err::Error for '{}', got {:?}", desc, result);
    }
}

#[test]
fn test_bond_input_standard_partial_fields() {
    let input = b"  1  2  1 1"; // len 11
    let mut parser = bond_input_standard();
    let result = parser.parse(input);
    assert!(
        result.is_err(),
        "Parser should have failed for partial field"
    );
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code,
            ErrorKind::Eof,
            "Mismatched error kind for partial field"
        );
    } else {
        panic!(
            "Expected nom::Err::Error for partial field, got {:?}",
            result
        );
    }
}

#[rstest]
#[case(b"  1  2  1\n", "Padded 9")]
#[case(b"  1  3  1  1  ", "Padded 12")]
fn test_bond_input_standard_whitespace_padded(#[case] input: &[u8], #[case] desc: &str) {
    let mut parser = bond_input_standard();
    let trimmed_input = input.trim_ascii_end();
    let result = parser.parse(trimmed_input);
    assert!(
        result.is_ok(),
        "Parser failed for case '{}': {:?}",
        desc,
        result
    );
    let (remaining, (_a1, _a2, _bond)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining should be empty for case '{}'",
        desc
    );
}
