use super::*;
use crate::atom::AtomStereoParity;
use float_cmp::approx_eq;
use nom::{error::ErrorKind, Err};
use rstest::rstest;

#[rstest]
#[case(b"H  ", "H", AtomSymbol::Element(Element::H))]
#[case(b"C  ", "C", AtomSymbol::Element(Element::C))]
#[case(b"Cu ", "Cu", AtomSymbol::Element(Element::Cu))]
#[case(b"D  ", "D", AtomSymbol::NamedIsotope(NamedIsotope::D))]
fn test_atom_symbol_standard(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected: AtomSymbol,
) {
    let result = atom_symbol_standard().parse(input);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, symbol) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(symbol, expected, "{} should have succeeded", desc);
}

#[rstest]
#[case(b"A  ", "unspecified atom", ErrorKind::MapRes)]
#[case(b"L  ", "atom list", ErrorKind::MapRes)]
#[case(b"LP ", "lone pair", ErrorKind::MapRes)]
#[case(b"R1 ", "R group", ErrorKind::MapRes)]
fn test_atom_symbol_standard_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = atom_symbol_standard().parse(input);
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
#[case(b"A  ", "A", AtomSymbol::Unspecified('A'))]
#[case(b"Q  ", "Q", AtomSymbol::Unspecified('Q'))]
#[case(b"*  ", "*", AtomSymbol::Unspecified('*'))]
#[case(b"L  ", "L", AtomSymbol::AtomList(AtomList { elements: vec![] }))]
#[case(b"LP ", "LP", AtomSymbol::LonePair)]
#[case(b"R1 ", "R1", AtomSymbol::RGroup(0))]
#[case(b"R3 ", "R3", AtomSymbol::RGroup(2))]
#[case(b"H  ", "H", AtomSymbol::Element(Element::H))]
#[case(b"C  ", "C", AtomSymbol::Element(Element::C))]
#[case(b" C ", "C", AtomSymbol::Element(Element::C))]
#[case(b"  C", "C", AtomSymbol::Element(Element::C))]
#[case(b"Cu ", "Cu", AtomSymbol::Element(Element::Cu))]
#[case(b"cu ", "cu", AtomSymbol::Element(Element::Cu))]
#[case(b"CU ", "CU", AtomSymbol::Element(Element::Cu))]
#[case(b"D  ", "D", AtomSymbol::NamedIsotope(NamedIsotope::D))]
#[case(b"d  ", "d", AtomSymbol::NamedIsotope(NamedIsotope::D))]
#[case(b"T  ", "T", AtomSymbol::NamedIsotope(NamedIsotope::T))]
fn test_atom_symbol(#[case] input: &[u8], #[case] desc: &str, #[case] expected: AtomSymbol) {
    let result = atom_symbol().parse(input);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, symbol) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(symbol, expected);
}

#[rstest]
#[case(b"   ", "empty field", ErrorKind::Alpha)]
#[case(b"H", "too short", ErrorKind::Eof)]
#[case(b"R  ", "R group index missing", ErrorKind::MapRes)]
#[case(b"R0 ", "R group index must be between 1 and 31", ErrorKind::MapRes)]
#[case(b"R32", "R group index must be between 1 and 31", ErrorKind::MapRes)]
#[case(b"Xx ", "Unknown atom symbol", ErrorKind::MapRes)]
#[case(b"LQ ", "Unknown atom symbol", ErrorKind::Eof)]
fn test_atom_symbol_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = atom_symbol().parse(input);
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
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0",
    "standard valid",
    Element::C,
    Some(10),
    1,
    Some(4),
    None
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4  0  0  0  0  0  0",
    "mass diff lower bound",
    Element::C,
    Some(9),
    1,
    Some(4),
    None
)]
#[case(
    b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4  0  0  0  0  0  0",
    "mass diff upper bound",
    Element::C,
    Some(16),
    1,
    Some(4),
    None
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -4  3  0  0  0  4  0  0  0  0  0  0",
    "mass diff out-of-range low",
    Element::C,
    None,
    1,
    Some(4),
    None
)]
#[case(
    b"    1.2345    2.3456    3.4567 C   5  3  0  0  0  4  0  0  0  0  0  0",
    "mass diff out-of-range high",
    Element::C,
    None,
    1,
    Some(4),
    None
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  8  0  0  0  4  0  0  0  0  0  0",
    "charge out-of-range high",
    Element::C,
    Some(10),
    0,
    Some(4),
    None
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  0  0  0  0  4  0  0  0  1  0  0",
    "atom map num non-zero",
    Element::C,
    Some(10),
    0,
    Some(4),
    Some(1)
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3  0     0  4  0  0  0  0  0  0",
    "blank block 1",
    Element::C,
    Some(10),
    1,
    Some(4),
    None
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4     0  0  0  0  0",
    "blank block 2",
    Element::C,
    Some(10),
    1,
    Some(4),
    None
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0      ",
    "blank block 3",
    Element::C,
    Some(10),
    1,
    Some(4),
    None
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4           1 0    ",
    "gaps with spaces and zeros",
    Element::C,
    Some(10),
    1,
    Some(4),
    Some(1)
)]
fn test_atom_input_standard69(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_element: Element,
    #[case] expected_isotope_mass: Option<u32>,
    #[case] expected_charge: i8,
    #[case] expected_valence: Option<u8>,
    #[case] expected_atom_map_num: Option<u32>,
) {
    let result = atom_input_standard69(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        desc,
        result
    );
    let (remaining, (atom, pos)) = result.unwrap();
    assert!(remaining.is_empty(), "Non-empty input for case '{}'", desc);
    assert_eq!(
        atom.element, expected_element,
        "Mismatched element for '{}'",
        desc
    );
    assert_eq!(
        atom.isotope_mass, expected_isotope_mass,
        "Mismatched isotope mass for '{}'",
        desc
    );
    assert_eq!(
        atom.charge, expected_charge,
        "Mismatched charge for '{}'",
        desc
    );
    assert_eq!(
        atom.valence, expected_valence,
        "Mismatched valence for '{}'",
        desc
    );
    assert_eq!(
        atom.atom_map_num, expected_atom_map_num,
        "Mismatched atom map num for '{}'",
        desc
    );
}

#[rstest]
#[case(
    b"    1.234a    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "non-numeric coordinate",
    ErrorKind::Eof
)]
#[case(
    b"    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  a  0  0",
    "non-numeric atom map number",
    ErrorKind::Digit
)]
#[case(
    b"    1.2345    2.3456    3.4567 L   0  0  0  0  0  0  0  0  0  0  0  0",
    "non-standard atom symbol",
    ErrorKind::MapRes
)]
fn test_atom_input_standard69_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = atom_input_standard69(input);
    assert!(result.is_err(), "{} should have failed", desc);
    let err = result.unwrap_err();
    if let Err::Error(e) = err {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    }
}

#[rstest]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4",
    "standard valid",
    Element::C,
    Some(10),
    1,
    Some(4)
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4",
    "mass diff lower bound",
    Element::C,
    Some(9),
    1,
    Some(4)
)]
#[case(
    b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4",
    "mass diff upper bound",
    Element::C,
    Some(16),
    1,
    Some(4)
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -4  3  0  0  0  4",
    "mass diff out-of-range low",
    Element::C,
    None,
    1,
    Some(4)
)]
#[case(
    b"    1.2345    2.3456    3.4567 C   5  3  0  0  0  4",
    "mass diff out-of-range high",
    Element::C,
    None,
    1,
    Some(4)
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  8  0  0  0  4",
    "charge out-of-range high",
    Element::C,
    Some(10),
    0,
    Some(4)
)]
fn test_atom_input_standard51(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_element: Element,
    #[case] expected_isotope_mass: Option<u32>,
    #[case] expected_charge: i8,
    #[case] expected_valence: Option<u8>,
) {
    let result = atom_input_standard51(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        desc,
        result
    );
    let (remaining, (atom, pos)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}': '{}'",
        desc,
        String::from_utf8_lossy(remaining)
    );
    assert_eq!(
        atom.element, expected_element,
        "Mismatched element for '{}'",
        desc
    );
    assert!(
        approx_eq!(f64, pos.x, 1.2345),
        "Mismatched x for '{}'",
        desc
    );
    assert!(
        approx_eq!(f64, pos.y, 2.3456),
        "Mismatched y for '{}'",
        desc
    );
    assert!(
        approx_eq!(f64, pos.z, 3.4567),
        "Mismatched z for '{}'",
        desc
    );
    assert_eq!(
        atom.isotope_mass, expected_isotope_mass,
        "Mismatched isotope mass for '{}'",
        desc
    );
    assert_eq!(
        atom.charge, expected_charge,
        "Mismatched charge for '{}'",
        desc
    );
    assert_eq!(
        atom.valence, expected_valence,
        "Mismatched valence for '{}'",
        desc
    );
}

#[rstest]
#[case(
    b"    1.234a    2.3456    3.4567 C  -2  3  0  0  0  4",
    "non-numeric coordinate",
    ErrorKind::Eof
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  a",
    "non-numeric valence",
    ErrorKind::Digit
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0 16",
    "out-of-range valence",
    ErrorKind::Verify
)]
#[case(
    b"    1.2345    2.3456    3.4567 L  -2  3  0  0  0  4",
    "invalid atom symbol",
    ErrorKind::MapRes
)]
fn test_atom_input_standard51_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = atom_input_standard51(input);
    assert!(result.is_err(), "{} should have failed", desc);
    let err = result.unwrap_err();
    assert!(
        matches!(err, Err::Error(_)),
        "Error should be a nom::Err::Error for '{}'",
        desc
    );
    if let Err::Error(e) = err {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    }
}

#[rstest]
#[case(
    b"    1.2345    2.3456    3.4567 C      3  0  0  0  4",
    "blank mass diff",
    None,
    1,
    Some(4)
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2     0  0  0  4",
    "blank charge",
    Some(10),
    0,
    Some(4)
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0   ",
    "blank valence",
    Some(10),
    1,
    None
)]
fn test_atom_input_standard51_empty_fields(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_isotope_mass: Option<u32>,
    #[case] expected_charge: i8,
    #[case] expected_valence: Option<u8>,
) {
    let result = atom_input_standard51(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}': '{}'",
        desc,
        String::from_utf8_lossy(remaining)
    );
    assert_eq!(
        atom.isotope_mass, expected_isotope_mass,
        "Mismatched isotope mass for '{}'",
        desc
    );
    assert_eq!(
        atom.charge, expected_charge,
        "Mismatched charge for '{}'",
        desc
    );
    assert_eq!(
        atom.valence, expected_valence,
        "Mismatched valence for '{}'",
        desc
    );
}

#[rstest]
#[case(
    b"    1.2345    2.3456    3.4567 C   0  3  1",
    "standard valid",
    Element::C,
    None,
    1,
    Some(AtomStereoParity::Odd)
)]
#[case(
    b"    1.2345    2.3456    3.4567 C   0  3   0  ",
    "blank stereo parity",
    Element::C,
    None,
    1,
    None
)]
fn test_atom_input_standard42(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_element: Element,
    #[case] expected_isotope_mass: Option<u32>,
    #[case] expected_charge: i8,
    #[case] expected_stereo_parity: Option<AtomStereoParity>,
) {
    let result = atom_input_standard42(input);
    assert!(
        result.is_ok(),
        "Parser failed for case '{}': {:?}",
        desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}': '{}'",
        desc,
        String::from_utf8_lossy(remaining)
    );
    assert_eq!(
        atom.element, expected_element,
        "Mismatched element for '{}'",
        desc
    );
    assert_eq!(
        atom.isotope_mass, expected_isotope_mass,
        "Mismatched isotope mass for '{}'",
        desc
    );
    assert_eq!(
        atom.charge, expected_charge,
        "Mismatched charge for '{}'",
        desc
    );
    assert_eq!(
        atom.stereo_parity, expected_stereo_parity,
        "Mismatched stereo parity for '{}'",
        desc
    );
}

#[rstest]
#[case(
    b"    1.2345    2.3456    3.4567 C   0  3  4",
    "stereo parity out of range",
    ErrorKind::Verify
)]
#[case(
    b"    1.234a    2.3456    3.4567 C   0  3  0",
    "non-numeric coordinate",
    ErrorKind::Eof
)]
#[case(
    b"    1.2345    2.3456    3.4567 L   0  3  0",
    "invalid atom symbol",
    ErrorKind::MapRes
)]
#[case(
    b"    1.2345    2.3456    3.4567 C   0  3  a",
    "non-numeric stereo parity",
    ErrorKind::Digit
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3     a",
    "non-numeric data in ignored block",
    ErrorKind::Verify
)]
fn test_atom_input_standard42_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = atom_input_standard42(input);
    assert!(result.is_err(), "{} should have failed", desc);
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    }
}

#[rstest]
#[case(
    b"    1.2345    2.3456    3.4567 C   0  3   ",
    "blank stereo parity",
    Element::C,
    None,
    1,
    None
)]
#[case(
    b"    1.2345    2.3456    3.4567 C      3  0",
    "blank mass diff",
    Element::C,
    None,
    1,
    None
)]
fn test_atom_input_standard42_empty_fields(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_element: Element,
    #[case] expected_isotope_mass: Option<u32>,
    #[case] expected_charge: i8,
    #[case] expected_stereo_parity: Option<AtomStereoParity>,
) {
    let result = atom_input_standard42(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}'",
        desc
    );
    assert_eq!(
        atom.element, expected_element,
        "Mismatched element for '{}'",
        desc
    );
    assert_eq!(
        atom.isotope_mass, expected_isotope_mass,
        "Mismatched isotope mass for '{}'",
        desc
    );
    assert_eq!(
        atom.charge, expected_charge,
        "Mismatched charge for '{}'",
        desc
    );
    assert_eq!(
        atom.stereo_parity, expected_stereo_parity,
        "Mismatched stereo parity for '{}'",
        desc
    );
}

#[rstest]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3  0",
    "valid numeric data in ignored gap"
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3         ",
    "whitespace in ignored gap"
)]
fn test_atom_input_standard42_ignored_gap(#[case] input: &[u8], #[case] desc: &str) {
    let result = atom_input_standard42(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining input should be empty for case '{}', but was: '{}'",
        desc,
        String::from_utf8_lossy(remaining)
    );
    assert_eq!(atom.charge, 1);
}

#[rstest]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  3",
    "standard valid",
    Element::C,
    Some(10),
    1
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -4  3",
    "mass diff out-of-range low",
    Element::C,
    None,
    1
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  8",
    "charge out-of-range high",
    Element::C,
    Some(10),
    0
)]
fn test_atom_input_standard39(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_element: Element,
    #[case] expected_isotope_mass: Option<u32>,
    #[case] expected_charge: i8,
) {
    let result = atom_input_standard39(input);
    assert!(
        result.is_ok(),
        "Parser failed for case '{}': {:?}",
        desc,
        result
    );
    let (remaining, (atom, pos)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}': '{}'",
        desc,
        String::from_utf8_lossy(remaining)
    );
    assert_eq!(
        atom.element, expected_element,
        "Mismatched element for '{}'",
        desc
    );
    assert_eq!(
        atom.isotope_mass, expected_isotope_mass,
        "Mismatched isotope mass for '{}'",
        desc
    );
    assert_eq!(
        atom.charge, expected_charge,
        "Mismatched charge for '{}'",
        desc
    );
    assert_eq!(atom.valence, None, "Valence is not None for '{}'", desc);
    assert!(
        approx_eq!(f64, pos.x, 1.2345),
        "Mismatched x for '{}'",
        desc
    );
}

#[rstest]
#[case(
    b"    1.234a    2.3456    3.4567 C  -2  3",
    "non-numeric coordinate",
    ErrorKind::Eof
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -a  3",
    "non-numeric mass diff",
    ErrorKind::Digit
)]
#[case(
    b"    1.2345    2.3456    3.4567 L  -2  3",
    "invalid atom symbol",
    ErrorKind::MapRes
)]

fn test_atom_input_standard39_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = atom_input_standard39(input);
    assert!(result.is_err(), "Parser should have failed for '{}'", desc);
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    } else {
        panic!(
            "Expected a nom::Err::Error for '{}', got {:?}",
            desc, result
        );
    }
}

#[rstest]
#[case(b"    1.2345    2.3456    3.4567 C      3", "blank mass diff", None, 1)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2   ",
    "blank charge",
    Some(10),
    0
)]
fn test_atom_input_standard39_empty_fields(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_isotope_mass: Option<u32>,
    #[case] expected_charge: i8,
) {
    let result = atom_input_standard39(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        desc,
        result
    );
    let (_, (atom, _)) = result.unwrap();
    assert_eq!(
        atom.isotope_mass, expected_isotope_mass,
        "Mismatched isotope mass for '{}'",
        desc
    );
    assert_eq!(
        atom.charge, expected_charge,
        "Mismatched charge for '{}'",
        desc
    );
}

#[rstest]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3\n", "trailing newline")]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3   ", "trailing spaces")]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3\t\t", "trailing tabs")]
fn test_atom_input_standard39_whitespace_padded(#[case] input: &[u8], #[case] desc: &str) {
    let result = all_consuming(terminated(atom_input_standard39, multispace0)).parse(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}': '{}'",
        desc,
        String::from_utf8_lossy(remaining)
    );
    assert_eq!(atom.charge, 1);
}

#[rstest]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2",
    "standard valid",
    Element::C,
    Some(10)
)]
#[case(
    b"    1.2345    2.3456    3.4567 C  -4",
    "mass diff out-of-range low",
    Element::C,
    None
)]
fn test_atom_input_standard36(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_element: Element,
    #[case] expected_isotope_mass: Option<u32>,
) {
    let result = atom_input_standard36(input);
    assert!(
        result.is_ok(),
        "Parser failed for case '{}': {:?}",
        desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}': '{}'",
        desc,
        String::from_utf8_lossy(remaining)
    );
    assert_eq!(
        atom.element, expected_element,
        "Mismatched element for '{}'",
        desc
    );
    assert_eq!(
        atom.isotope_mass, expected_isotope_mass,
        "Mismatched isotope mass for '{}'",
        desc
    );
    assert_eq!(atom.charge, 0, "Charge should be 0 for '{}'", desc);
    assert_eq!(atom.valence, None, "Valence should be None for '{}'", desc);
}

#[rstest]
#[case(
    b"    1.234a    2.3456    3.4567 C  -2",
    "non-numeric coordinate",
    ErrorKind::Eof
)]
#[case(
    b"    1.2345    2.3456    3.4567 L  -2",
    "invalid atom symbol",
    ErrorKind::MapRes
)]
fn test_atom_input_standard36_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = atom_input_standard36(input);
    assert!(result.is_err(), "{} should have failed", desc);
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    } else {
        panic!(
            "Expected a nom::Err::Error for '{}', got {:?}",
            desc, result
        );
    }
}

#[rstest]
#[case(b"    1.2345    2.3456    3.4567 C    ", "blank mass diff")]
fn test_atom_input_standard36_empty_fields(#[case] input: &[u8], #[case] desc: &str) {
    let result = atom_input_standard36(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        desc,
        result
    );
    let (_, (atom, _)) = result.unwrap();
    assert_eq!(
        atom.isotope_mass, None,
        "Mismatched isotope mass for '{}'",
        desc
    );
}

#[rstest]
#[case(
    b"    1.2345    2.3456    3.4567 C  -2  \n",
    "trailing whitespace and newline"
)]
fn test_atom_input_standard36_whitespace_padded(#[case] input: &[u8], #[case] _desc: &str) {
    let result = all_consuming(terminated(atom_input_standard36, multispace0)).parse(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        _desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining input should be empty for case '{}'",
        _desc
    );
    assert_eq!(atom.isotope_mass, Some(10));
}

#[rstest]
#[case(b"    1.2345    2.3456    3.4567 C  ", "standard valid", Element::C)]
fn test_atom_input_standard34(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_element: Element,
) {
    let result = atom_input_standard34(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}': '{}'",
        desc,
        String::from_utf8_lossy(remaining)
    );
    assert_eq!(
        atom.element, expected_element,
        "Mismatched element for '{}'",
        desc
    );
    assert_eq!(
        atom.isotope_mass, None,
        "Isotope mass should be None for '{}'",
        desc
    );
    assert_eq!(atom.charge, 0, "Charge should be 0 for '{}'", desc);
    assert_eq!(atom.valence, None, "Valence should be None for '{}'", desc);
}

#[rstest]
#[case(
    b"    1.234a    2.3456    3.4567 C  ",
    "non-numeric coordinate",
    ErrorKind::Eof
)]
#[case(
    b"    1.2345    2.3456    3.4567 L  ",
    "invalid atom symbol",
    ErrorKind::MapRes
)]
fn test_atom_input_standard34_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let result = atom_input_standard34(input);
    assert!(result.is_err(), "{} should have failed", desc);
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    } else {
        panic!(
            "Expected a nom::Err::Error for '{}', got {:?}",
            desc, result
        );
    }
}

#[rstest]
#[case(
    b"    1.2345    2.3456    3.4567 C    \n",
    "trailing whitespace and newline"
)]
fn test_atom_input_standard34_whitespace_padded(#[case] input: &[u8], #[case] desc: &str) {
    let result = all_consuming(terminated(atom_input_standard34, multispace0)).parse(input);
    assert!(
        result.is_ok(),
        "{} should have succeeded: {:?}",
        desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(remaining.is_empty(), "Remaining non-empty for {}", desc);
    assert_eq!(atom.element, Element::C);
}

#[rstest]
#[case(
    b"    1.0000    2.0000    3.0000 C  ",
    "len 34",
    Element::C,
    None,
    0,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C   ",
    "len 35 padded",
    Element::C,
    None,
    0,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2",
    "len 36",
    Element::C,
    Some(10),
    0,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  ",
    "len 38 padded",
    Element::C,
    Some(10),
    0,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3",
    "len 39",
    Element::C,
    Some(10),
    1,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0",
    "len 42 with numeric data in ignored block",
    Element::C,
    Some(10),
    1,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4",
    "len 51",
    Element::C,
    Some(10),
    1,
    Some(4)
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4 ",
    "len 52 padded",
    Element::C,
    Some(10),
    1,
    Some(4)
)]
#[case(
    b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "len 69 zeros",
    Element::C,
    None,
    0,
    None
)]
fn test_atom_input_standard(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_element: Element,
    #[case] expected_isotope_mass: Option<u32>,
    #[case] expected_charge: i8,
    #[case] expected_valence: Option<u8>,
) {
    let mut parser = atom_input_standard();
    let result = parser.parse(input);
    assert!(
        result.is_ok(),
        "Parser failed for case '{}': {:?}",
        desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}': '{}'",
        desc,
        String::from_utf8_lossy(remaining)
    );
    assert_eq!(
        atom.element, expected_element,
        "Mismatched element for '{}'",
        desc
    );
    assert_eq!(
        atom.isotope_mass, expected_isotope_mass,
        "Mismatched isotope mass for '{}'",
        desc
    );
    assert_eq!(
        atom.charge, expected_charge,
        "Mismatched charge for '{}'",
        desc
    );
    assert_eq!(
        atom.valence, expected_valence,
        "Mismatched valence for '{}'",
        desc
    );
}

#[rstest]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  a",
    "non-numeric data in ignored block",
    ErrorKind::Verify
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  a",
    "non-numeric valence",
    ErrorKind::Digit
)]
#[case(
    b"    1.0000    2.0000    3.0000 L  -2  3  0  0  0  4",
    "invalid element",
    ErrorKind::MapRes
)]
#[case(b"", "empty", ErrorKind::Eof)]
#[case(b"    1.0000    2.0000    3.0000 ", "too short", ErrorKind::Eof)]
#[case(
    b"    1.0000    2.0000    3.0000 C  a",
    "len 35 trailing non-numeric data",
    ErrorKind::Eof
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2 a",
    "len 38 trailing non-numeric data",
    ErrorKind::Eof
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2 a",
    "len 38 trailing non-numeric data",
    ErrorKind::Eof
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3 a",
    "len 41 trailing non-numeric data",
    ErrorKind::Eof
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0 a",
    "len 44 trailing non-numeric data",
    ErrorKind::Eof
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0 a",
    "len 47 trailing non-numeric data",
    ErrorKind::Eof
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0 a",
    "len 50 trailing non-numeric data",
    ErrorKind::Eof
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  0 a",
    "len 53 trailing non-numeric data",
    ErrorKind::Eof
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  0  0  0  0  0  0  0  a",
    "len 72 trailing non-numeric data",
    ErrorKind::Eof
)]
fn test_atom_input_standard_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let mut parser = atom_input_standard();
    let result = parser.parse(input);
    assert!(result.is_err(), "Parser should have failed for '{}'", desc);
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    } else {
        panic!(
            "Expected a nom::Err::Error for '{}', got {:?}",
            desc, result
        );
    }
}

#[test]
fn test_atom_input_standard_partial_fields() {
    let input = b"    1.0000    2.0000    3.0000 C  -2 3"; // len 38
    let mut parser = atom_input_standard();
    let result = parser.parse(input);
    assert!(
        result.is_err(),
        "Parser should have failed for partial field"
    );
    if let Err(Err::Error(e)) = result {
        // The charge field is 3 chars, we provided 2, fixed_width_opt should see a partial non-whitespace field and fail with Eof
        assert_eq!(
            e.code,
            ErrorKind::Eof,
            "Mismatched error kind for partial field"
        );
    } else {
        panic!(
            "Expected a nom::Err::Error for partial field, got {:?}",
            result
        );
    }
}

#[rstest]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4   \t", "len 55")]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0           ",
    "len 80"
)]
fn test_atom_input_standard_whitespace_padded(#[case] input: &[u8], #[case] desc: &str) {
    let mut parser = atom_input_standard();
    let result = parser.parse(input);
    assert!(
        result.is_ok(),
        "Parser failed for whitespace padded input: {:?}",
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(remaining.is_empty(), "Non-empty input for case '{}'", desc);
    assert_eq!(atom.charge, 1);
    assert_eq!(atom.valence, Some(4));
}

#[rstest]
#[case(
    b"    1.0000    2.0000    3.0000 C  ",
    "len 34",
    AtomSymbol::Element(Element::C),
    None,
    0,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C   ",
    "len 35 padded",
    AtomSymbol::Element(Element::C),
    None,
    0,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2",
    "len 36",
    AtomSymbol::Element(Element::C),
    Some(10),
    0,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  ",
    "len 38 padded",
    AtomSymbol::Element(Element::C),
    Some(10),
    0,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3",
    "len 39",
    AtomSymbol::Element(Element::C),
    Some(10),
    1,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0",
    "len 42 with numeric data in ignored block",
    AtomSymbol::Element(Element::C),
    Some(10),
    1,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4",
    "len 51",
    AtomSymbol::Element(Element::C),
    Some(10),
    1,
    Some(4)
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4 ",
    "len 52 padded",
    AtomSymbol::Element(Element::C),
    Some(10),
    1,
    Some(4)
)]
#[case(
    b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
    "len 69 zeros",
    AtomSymbol::Element(Element::C),
    None,
    0,
    None
)]
#[case(
    b"    1.0000    2.0000    3.0000 L   0  0  0  0  0  4",
    "atom list",
    AtomSymbol::AtomList(AtomList { elements: vec![]}),
    None,
    0,
    Some(4)
)]
fn test_atom_input(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_symbol: AtomSymbol,
    #[case] expected_isotope_mass: Option<u32>,
    #[case] expected_charge: i8,
    #[case] expected_valence: Option<u8>,
) {
    let mut parser = atom_input();
    let result = parser.parse(input);
    assert!(
        result.is_ok(),
        "Parser failed for case '{}': {:?}",
        desc,
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "Remaining non-empty for case '{}': '{}'",
        desc,
        String::from_utf8_lossy(remaining)
    );
    assert_eq!(
        atom.symbol, expected_symbol,
        "Mismatched element for '{}'",
        desc
    );
    assert_eq!(
        atom.isotope_mass, expected_isotope_mass,
        "Mismatched isotope mass for '{}'",
        desc
    );
    assert_eq!(
        atom.charge, expected_charge,
        "Mismatched charge for '{}'",
        desc
    );
    assert_eq!(
        atom.valence, expected_valence,
        "Mismatched valence for '{}'",
        desc
    );
}

#[rstest]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  a",
    "non-numeric data in ignored block",
    ErrorKind::Digit
)]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  a",
    "non-numeric valence",
    ErrorKind::Digit
)]

fn test_atom_input_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: ErrorKind,
) {
    let mut parser = atom_input();
    let result = parser.parse(input);
    assert!(result.is_err(), "Parser should have failed for '{}'", desc);
    if let Err(Err::Error(e)) = result {
        assert_eq!(
            e.code, expected_kind,
            "Mismatched error kind for '{}'",
            desc
        );
    } else {
        panic!(
            "Expected a nom::Err::Error for '{}', got {:?}",
            desc, result
        );
    }
}

#[test]
fn test_atom_input_partial_fields() {
    let input = b"    1.0000    2.0000    3.0000 C  -2 3"; // len 38
    let mut parser = atom_input();
    let result = parser.parse(input);
    assert!(
        result.is_err(),
        "Parser should have failed for partial field"
    );
    if let Err(Err::Error(e)) = result {
        // The charge field is 3 chars, we provided 2, fixed_width_opt should see a partial non-whitespace field and fail with Eof
        assert_eq!(
            e.code,
            ErrorKind::Eof,
            "Mismatched error kind for partial field"
        );
    } else {
        panic!(
            "Expected a nom::Err::Error for partial field, got {:?}",
            result
        );
    }
}

#[rstest]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4   \t", "len 55")]
#[case(
    b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0           ",
    "len 80"
)]
fn test_atom_input_whitespace_padded(#[case] input: &[u8], #[case] desc: &str) {
    let mut parser = atom_input();
    let result = parser.parse(input);
    assert!(
        result.is_ok(),
        "Parser failed for whitespace padded input: {:?}",
        result
    );
    let (remaining, (atom, _)) = result.unwrap();
    assert!(remaining.is_empty(), "Non-empty input for case '{}'", desc);
    assert_eq!(atom.charge, 1);
    assert_eq!(atom.valence, Some(4));
}
