use super::*;
use nom::{error::ErrorKind, Err};
use rstest::rstest;

#[rstest]
#[case(b"  1   1  -1", vec![ChargeEntry { atom_index: 0, charge: -1 }])]
#[case(b"  2   1  -1   4   1", vec![ChargeEntry { atom_index: 0, charge: -1 }, ChargeEntry { atom_index: 3, charge: 1 }])]
#[case(b"  8   1   1   2   2   3   3   4   4   5   5   6   6   7   7   8   8",
   vec![
        ChargeEntry { atom_index: 0, charge: 1 },
        ChargeEntry { atom_index: 1, charge: 2 },
        ChargeEntry { atom_index: 2, charge: 3 },
        ChargeEntry { atom_index: 3, charge: 4 },
        ChargeEntry { atom_index: 4, charge: 5 },
        ChargeEntry { atom_index: 5, charge: 6 },
        ChargeEntry { atom_index: 6, charge: 7 },
        ChargeEntry { atom_index: 7, charge: 8 },
    ])]
#[case(b"  1  25  15", vec![ChargeEntry { atom_index: 24, charge: 15 }])]
fn test_charge_entries(#[case] input: &[u8], #[case] expected: Vec<ChargeEntry>) {
    let (remaining, result) = charge_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"M  CHG  1   1  -1  a", "trailing chars", ErrorKind::Digit)]
#[case(b"M  CHG  2   1  -1", "count does not match item list", ErrorKind::Digit)]
#[case(
    b"M  CHG  1   1  -1   4   1",
    "item list longer than count",
    ErrorKind::Digit
)]
#[case(b"M  CHG  0", "count is zero", ErrorKind::Digit)]
#[case(b"M  CHG  1   0 -10", "atom index is zero", ErrorKind::Digit)]
#[case(b"M  XXX  1   1  -1", "invalid property tag", ErrorKind::Digit)]
#[case(b"X  CHG  1   1  -1", "invalid prefix", ErrorKind::Digit)]
fn test_charge_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = charge_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1   2", vec![RadicalEntry { atom_index: 0, radical_type: 2 }])]
#[case(b"  2   1   1   4   3", vec![RadicalEntry { atom_index: 0, radical_type: 1 }, RadicalEntry { atom_index: 3, radical_type: 3 }])]
fn test_radical_entries(#[case] input: &[u8], #[case] expected: Vec<RadicalEntry>) {
    let (remaining, result) = radical_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"M  RAD  1   1   4", "value out of range", ErrorKind::Digit)]
#[case(b"M  RAD  1   1  -1", "value out of range", ErrorKind::Digit)]
fn test_radical_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = radical_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"  1   1  13", vec![IsotopeEntry { atom_index: 0, mass: 13 }])]
#[case(b"  2  12   2  15  14", vec![IsotopeEntry { atom_index: 11, mass: 2 }, IsotopeEntry { atom_index: 14, mass: 14 }])]
fn test_isotope_entries(#[case] input: &[u8], #[case] expected: Vec<IsotopeEntry>) {
    let (remaining, result) = isotope_entries().parse(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case(b"M  ISO  1   1  13", "value out of range", ErrorKind::Digit)]
#[case(b"M  ISO  1   1  -1", "value out of range", ErrorKind::Digit)]
fn test_isotope_entries_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = isotope_entries().parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "Mismatched error kind for {}, expected {:?}, got {}",
        desc,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}

#[rstest]
#[case(b"M  CHG  1   1  -1", PropertyEntries::ChargeEntries(vec![ChargeEntry { atom_index: 0, charge: -1 }]))]
#[case(b"M  RAD  1   1   2", PropertyEntries::RadicalEntries(vec![RadicalEntry { atom_index: 0, radical_type: 2 }]))]
#[case(b"M  ISO  1   1  13", PropertyEntries::IsotopeEntries(vec![IsotopeEntry { atom_index: 0, mass: 13 }]))]
fn test_property_input(#[case] input: &[u8], #[case] expected: PropertyEntries) {
    let (remaining, result) = property_input_standard(input).unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}
