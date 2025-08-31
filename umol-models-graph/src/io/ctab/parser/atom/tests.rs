use super::*;
use crate::io::config::ParseFlags;
use crate::io::ctab::atom::{
    AtomExactChange, AtomInversionRetention, AtomStereoCare, AtomStereoParity,
};
use crate::io::ctab::query::QueryAtom;
use float_cmp::approx_eq;
use nom::{error, Err};
use pretty_assertions::assert_eq;
use rstest::*;

#[rstest]
#[case(b"H  ", "H", AtomSymbol::Element(Element::H))]
#[case(b"C  ", "C", AtomSymbol::Element(Element::C))]
#[case(b"Cu ", "Cu", AtomSymbol::Element(Element::Cu))]
#[case(b"D  ", "D", AtomSymbol::NamedIsotope(NamedIsotope::D))]
fn test_atom_symbol(#[case] input: &[u8], #[case] desc: &str, #[case] expected: AtomSymbol) {
    let result = atom_symbol(ParseFlags::LENIENT).parse(input);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, symbol) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(symbol, expected, "{} should have succeeded", desc);
}

#[rstest]
#[case(b"A  ", "query atom", error::ErrorKind::MapRes)]
#[case(b"L  ", "atom list", error::ErrorKind::MapRes)]
#[case(b"LP ", "lone pair", error::ErrorKind::MapRes)]
#[case(b"R1 ", "R group", error::ErrorKind::MapRes)]
fn test_atom_symbol_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atom_symbol(ParseFlags::STRICT).parse(input);
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
#[case(b"A  ", "A", AtomSymbol::Query(QueryAtom::Heavy))]
#[case(b"Q  ", "Q", AtomSymbol::Query(QueryAtom::Heteroatom))]
#[case(b"QH ", "QH", AtomSymbol::Query(QueryAtom::HeteroatomOrH))]
#[case(b"*  ", "*", AtomSymbol::Query(QueryAtom::Any))]
#[case(b"L  ", "L", AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }))]
#[case(b"LP ", "LP", AtomSymbol::LonePair)]
#[case(b"R# ", "R#", AtomSymbol::RGroup(RGroup::new(None)))]
#[case(b"R1 ", "R1", AtomSymbol::RGroup(RGroup::new(Some(1))))]
#[case(b"R3 ", "R3", AtomSymbol::RGroup(RGroup::new(Some(3))))]
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
fn test_atomlike_symbol(#[case] input: &[u8], #[case] desc: &str, #[case] expected: AtomSymbol) {
    let result = atomlike_symbol(ParseFlags::LENIENT).parse(input);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, symbol) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        symbol, expected,
        "{} has returned symbol {:?}, expected {:?}",
        desc, symbol, expected
    );
}

#[rstest]
#[case(b"   ", "empty field", error::ErrorKind::MapRes)]
#[case(b"H", "too short", error::ErrorKind::Eof)]
#[case(b"Xx ", "Unknown atom symbol", error::ErrorKind::MapRes)]
#[case(b"LQ ", "Unknown atom symbol", error::ErrorKind::MapRes)]
#[case(b"A  ", "query atom", error::ErrorKind::MapRes)]
#[case(b"L  ", "atom list", error::ErrorKind::MapRes)]
#[case(b"R1 ", "R group", error::ErrorKind::MapRes)]
#[case(b"LP ", "lone pair", error::ErrorKind::MapRes)]
fn test_atomlike_symbol_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atomlike_symbol(ParseFlags::STRICT).parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
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
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0", "basic valid",
       Element::C, Some(10), 1, Some(4), None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4  0  0  0  0  0  0", "mass diff lower bound",
       Element::C, Some(9), 1, Some(4), None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4  0  0  0  0  0  0", "mass diff upper bound",
       Element::C, Some(16), 1, Some(4), None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -4  3  0  0  0  4  0  0  0  0  0  0", "mass diff out-of-range low",
       Element::C, None, 1, Some(4), None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C   5  3  0  0  0  4  0  0  0  0  0  0", "mass diff out-of-range high",
       Element::C, None, 1, Some(4), None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  8  0  0  0  4  0  0  0  0  0  0", "charge out-of-range high",
       Element::C, Some(10), 0, Some(4), None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  0  0  0  0  4  0  0  0  1  0  0", "atom map num non-zero",
       Element::C, Some(10), 0, Some(4), Some(1), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0     0  4  0  0  0  0  0  0", "blank block 1",
       Element::C, Some(10), 1, Some(4), None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4     0  0  0  0  0", "blank block 2",
       Element::C, Some(10), 1, Some(4), None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0      ", "blank block 3",
       Element::C, Some(10), 1, Some(4), None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4           1 0    ", "gaps with spaces and zeros",
       Element::C, Some(10), 1, Some(4), Some(1), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0", "named isotope",
       Element::H, Some(2), 1, Some(1), None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4  0  0  0  0  0  0", "non-strict padding",
       Element::C, Some(10), 1, Some(4), None, 1.2345, 2.3456, 3.4567)]
fn test_atom_input69(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: i8,
    #[case] valence: Option<u8>,
    #[case] atom_map_num: Option<u32>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input69(input, ParseFlags::LENIENT);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.element, element,
        "{} has returned element {:?}, expected {:?}",
        desc, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, charge,
    );
    assert_eq!(
        atom.valence, valence,
        "{} has returned valence {:?}, expected {:?}",
        desc, atom.valence, valence,
    );
    assert_eq!(
        atom.atom_map_num, atom_map_num,
        "{} has returned atom map num {:?}, expected {:?}",
        desc, atom.atom_map_num, atom_map_num,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().x, x_position),
        "{} has returned x {:?}, expected {:?}",
        desc,
        atom.position.unwrap().x,
        x_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().y, y_position),
        "{} has returned y {:?}, expected {:?}",
        desc,
        atom.position.unwrap().y,
        y_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().z, z_position),
        "{} has returned z {:?}, expected {:?}",
        desc,
        atom.position.unwrap().z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.234a    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0", "non-numeric coordinate", error::ErrorKind::Eof)]
#[case(b"    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  a  0  0", "non-numeric atom map number", error::ErrorKind::Digit)]
#[case(b"    1.2345    2.3456    3.4567 L   0  0  0  0  0  0  0  0  0  0  0  0", "atom list", error::ErrorKind::MapRes)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4  0  0  0  0  0  0", "non-strict padding", error::ErrorKind::Verify)]
fn test_atom_input69_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atom_input69(input, ParseFlags::STRICT);
    assert!(result.is_err(), "{} should have failed", desc);
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
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4", "basic valid", Element::C, Some(10), 1, Some(4), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4", "mass diff lower bound", Element::C, Some(9), 1, Some(4), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4", "mass diff upper bound", Element::C, Some(16), 1, Some(4), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -4  3  0  0  0  4", "mass diff out-of-range low", Element::C, None, 1, Some(4), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C   5  3  0  0  0  4", "mass diff out-of-range high", Element::C, None, 1, Some(4), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  8  0  0  0  4", "charge out-of-range high", Element::C, Some(10), 0, Some(4), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1", "named isotope", Element::H, Some(2), 1, Some(1), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4", "non-strict padding", Element::C, Some(10), 1, Some(4), 1.2345, 2.3456, 3.4567)]
fn test_atom_input51(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: i8,
    #[case] valence: Option<u8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input51(input, ParseFlags::LENIENT);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.element, element,
        "{} has returned element {:?}, expected {:?}",
        desc, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, charge,
    );
    assert_eq!(
        atom.valence, valence,
        "{} has returned valence {:?}, expected {:?}",
        desc, atom.valence, valence,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().x, x_position),
        "{} has returned x {:?}, expected {:?}",
        desc,
        atom.position.unwrap().x,
        x_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().y, y_position),
        "{} has returned y {:?}, expected {:?}",
        desc,
        atom.position.unwrap().y,
        y_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().z, z_position),
        "{} has returned z {:?}, expected {:?}",
        desc,
        atom.position.unwrap().z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.234a    2.3456    3.4567 C  -2  3  0  0  0  4", "non-numeric coordinate", error::ErrorKind::Eof)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  a", "non-numeric valence", error::ErrorKind::Digit)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0 16", "out-of-range valence", error::ErrorKind::Verify)]
#[case(b"    1.2345    2.3456    3.4567 L  -2  3  0  0  0  4", "invalid atom symbol", error::ErrorKind::MapRes)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4", "non-strict padding", error::ErrorKind::Verify)]
fn test_atom_input51_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atom_input51(input, ParseFlags::STRICT);
    assert!(result.is_err(), "{} should have failed", desc);
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
#[case(b"    1.2345    2.3456    3.4567 C      3  0  0  0  4", "blank mass diff", None, 1, Some(4), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2     0  0  0  4", "blank charge", Some(10), 0, Some(4), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0   ", "blank valence", Some(10), 1, None, 1.2345, 2.3456, 3.4567)]
fn test_atom_input51_empty_fields(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: i8,
    #[case] valence: Option<u8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input51(input, ParseFlags::LENIENT);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, charge,
    );
    assert_eq!(
        atom.valence, valence,
        "{} has returned valence {:?}, expected {:?}",
        desc, atom.valence, valence,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().x, x_position),
        "{} has returned x {:?}, expected {:?}",
        desc,
        atom.position.unwrap().x,
        x_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().y, y_position),
        "{} has returned y {:?}, expected {:?}",
        desc,
        atom.position.unwrap().y,
        y_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().z, z_position),
        "{} has returned z {:?}, expected {:?}",
        desc,
        atom.position.unwrap().z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.2345    2.3456    3.4567 C   0  3  1", "basic valid", Element::C, None, 1, Some(AtomStereoParity::Odd), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C   0  3   0  ", "blank stereo parity", Element::C, None, 1, None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 D  -2  3  0", "named isotope", Element::H, Some(2), 1, None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0", "non-strict padding", Element::C, Some(10), 1, None, 1.2345, 2.3456, 3.4567)]
fn test_atom_input42(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: i8,
    #[case] stereo_parity: Option<AtomStereoParity>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input42(input, ParseFlags::LENIENT);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.element, element,
        "{} has returned element {:?}, expected {:?}",
        desc, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, charge,
    );
    assert_eq!(
        atom.stereo_parity, stereo_parity,
        "{} has returned stereo parity {:?}, expected {:?}",
        desc, atom.stereo_parity, stereo_parity,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().x, x_position),
        "{} has returned x {:?}, expected {:?}",
        desc,
        atom.position.unwrap().x,
        x_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().y, y_position),
        "{} has returned y {:?}, expected {:?}",
        desc,
        atom.position.unwrap().y,
        y_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().z, z_position),
        "{} has returned z {:?}, expected {:?}",
        desc,
        atom.position.unwrap().z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.2345    2.3456    3.4567 C   0  3  4", "stereo parity out of range", error::ErrorKind::Verify)]
#[case(b"    1.234a    2.3456    3.4567 C   0  3  0", "non-numeric coordinate", error::ErrorKind::Eof)]
#[case(b"    1.2345    2.3456    3.4567 L   0  3  0", "invalid atom symbol", error::ErrorKind::MapRes)]
#[case(b"    1.2345    2.3456    3.4567 C   0  3  a", "non-numeric stereo parity", error::ErrorKind::Digit)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3     a", "non-numeric data in ignored block", error::ErrorKind::Verify)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0", "non-strict padding", error::ErrorKind::Verify)]
fn test_atom_input42_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atom_input42(input, ParseFlags::STRICT);
    assert!(result.is_err(), "{} should have failed", desc);
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
#[case(b"    1.2345    2.3456    3.4567 C   0  3   ", "blank stereo parity", Element::C, None, 1, None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C      3  0", "blank mass diff", Element::C, None, 1, None, 1.2345, 2.3456, 3.4567)]
fn test_atom_input42_empty_fields(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: i8,
    #[case] stereo_parity: Option<AtomStereoParity>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input42(input, ParseFlags::LENIENT);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.element, element,
        "{} has returned element {:?}, expected {:?}",
        desc, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, charge,
    );
    assert_eq!(
        atom.stereo_parity, stereo_parity,
        "{} has returned stereo parity {:?}, expected {:?}",
        desc, atom.stereo_parity, stereo_parity,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().x, x_position),
        "{} has returned x {:?}, expected {:?}",
        desc,
        atom.position.unwrap().x,
        x_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().y, y_position),
        "{} has returned y {:?}, expected {:?}",
        desc,
        atom.position.unwrap().y,
        y_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().z, z_position),
        "{} has returned z {:?}, expected {:?}",
        desc,
        atom.position.unwrap().z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.2345    2.3456    3.4567 C  -2  8", "charge out-of-range high", Element::C, Some(10), 0, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 D  -2  3", "named isotope", Element::H, Some(2), 1, 1.2345, 2.3456, 3.4567)]
fn test_atom_input39(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: i8,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input39(input, ParseFlags::LENIENT);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    let position = atom.position.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.element, element,
        "{} has returned element {:?}, expected {:?}",
        desc, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, charge,
    );
    assert_eq!(
        atom.valence, None,
        "{} has returned valence {:?}, expected {:?}",
        desc, atom.valence, None as Option<u8>
    );
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{} has returned x {:?}, expected {:?}",
        desc,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{} has returned y {:?}, expected {:?}",
        desc,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{} has returned z {:?}, expected {:?}",
        desc,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.234a    2.3456    3.4567 C  -2  3", "non-numeric coordinate", error::ErrorKind::Eof)]
#[case(b"    1.2345    2.3456    3.4567 C  -a  3", "non-numeric mass diff", error::ErrorKind::Digit)]
#[case(b"    1.2345    2.3456    3.4567 L  -2  3", "invalid atom symbol", error::ErrorKind::MapRes)]       
fn test_atom_input39_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atom_input39(input, ParseFlags::STRICT);
    assert!(result.is_err(), "{} should have failed", desc);
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
#[case(b"    1.2345    2.3456    3.4567 C      3", "blank mass diff", None, 1)]
#[case(b"    1.2345    2.3456    3.4567 C  -2   ", "blank charge", Some(10), 0)]
fn test_atom_input39_empty_fields(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: i8,
) {
    let result = atom_input39(input, ParseFlags::LENIENT);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, charge,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.2345    2.3456    3.4567 C  -2", "basic valid", Element::C, Some(10), 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -4", "mass diff out-of-range low", Element::C, None, 1.2345, 2.3456, 3.4567)]
fn test_atom_input36(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input36(input, ParseFlags::LENIENT);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.element, element,
        "{} has returned element {:?}, expected {:?}",
        desc, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, 0,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, 0
    );
    assert_eq!(
        atom.valence, None,
        "{} has returned valence {:?}, expected {:?}",
        desc, atom.valence, None as Option<u8>
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().x, x_position),
        "{} has returned x {:?}, expected {:?}",
        desc,
        atom.position.unwrap().x,
        x_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().y, y_position),
        "{} has returned y {:?}, expected {:?}",
        desc,
        atom.position.unwrap().y,
        y_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().z, z_position),
        "{} has returned z {:?}, expected {:?}",
        desc,
        atom.position.unwrap().z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.234a    2.3456    3.4567 C  -2", "non-numeric coordinate", error::ErrorKind::Eof)]
#[case(b"    1.2345    2.3456    3.4567 L  -2", "invalid atom symbol", error::ErrorKind::MapRes)]
fn test_atom_input36_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atom_input36(input, ParseFlags::STRICT);
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
#[case(b"    1.2345    2.3456    3.4567 C    ", "blank mass diff", None)]
fn test_atom_input36_empty_fields(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] isotope_mass: Option<u32>,
) {
    let result = atom_input36(input, ParseFlags::LENIENT);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.2345    2.3456    3.4567 C  ", "basic valid", Element::C, None, 1.2345, 2.3456, 3.4567)]
fn test_atom_input34(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input34(input, ParseFlags::LENIENT);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.element, element,
        "{} has returned element {:?}, expected {:?}",
        desc, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, 0,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, 0
    );
    assert_eq!(
        atom.valence, None,
        "{} has returned valence {:?}, expected {:?}",
        desc, atom.valence, None as Option<u8>
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().x, x_position),
        "{} has returned x {:?}, expected {:?}",
        desc,
        atom.position.unwrap().x,
        x_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().y, y_position),
        "{} has returned y {:?}, expected {:?}",
        desc,
        atom.position.unwrap().y,
        y_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().z, z_position),
        "{} has returned z {:?}, expected {:?}",
        desc,
        atom.position.unwrap().z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.234a    2.3456    3.4567 C  ", "non-numeric coordinate", error::ErrorKind::Eof)]
#[case(b"    1.2345    2.3456    3.4567 L  ", "invalid atom symbol", error::ErrorKind::MapRes)]
fn test_atom_input34_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atom_input34(input, ParseFlags::STRICT);
    assert!(result.is_err(), "{} should have failed", desc);
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
#[case(b"    1.0000    2.0000    3.0000 C  ", "len 34", Element::C, None, 0, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C   ", "len 35 padded", Element::C, None, 0, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C  -2", "len 36", Element::C, Some(10), 0, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3", "len 39", Element::C, Some(10), 1, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  1", "len 42 with stereo parity", Element::C, Some(10), 1, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4", "len 51 with query fields", Element::C, Some(10), 1, Some(4), None, 1.0000, 2.0000, 3.0000)]
#[case(b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0           1  2  1", "len 69 with reaction fields", Element::C, None, 0, None, Some(1), 0.0000, 0.0000, 0.0000)]
fn test_atom_input(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: i8,
    #[case] valence: Option<u8>,
    #[case] atom_map_num: Option<u32>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input(ParseFlags::LENIENT).parse(input);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} should consume all input", desc);
    
    assert_eq!(
        atom.element, element,
        "{} has returned element {:?}, expected {:?}",
        desc, atom.element, element
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope_mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass
    );
    assert_eq!(
        atom.charge, charge,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, charge
    );
    assert_eq!(
        atom.valence, valence,
        "{} has returned valence {:?}, expected {:?}",
        desc, atom.valence, valence
    );
    assert_eq!(
        atom.atom_map_num, atom_map_num,
        "{} has returned atom_map_num {:?}, expected {:?}",
        desc, atom.atom_map_num, atom_map_num
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().x, x_position),
        "{} has returned x {:?}, expected {:?}",
        desc,
        atom.position.unwrap().x,
        x_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().y, y_position),
        "{} has returned y {:?}, expected {:?}",
        desc,
        atom.position.unwrap().y,
        y_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().z, z_position),
        "{} has returned z {:?}, expected {:?}",
        desc,
        atom.position.unwrap().z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.234a    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0", "len 69, non-numeric coordinate", error::ErrorKind::Eof)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  a", "len 51, non-numeric valence", error::ErrorKind::Digit)]
#[case(b"    1.2345    2.3456    3.4567 L   0  3  0", "len 42, invalid atom symbol", error::ErrorKind::MapRes)]       
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4", "len 51, non-strict padding", error::ErrorKind::Verify)]
fn test_atom_input_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let result = atom_input(ParseFlags::STRICT).parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
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
#[case(b"    1.0000    2.0000    3.0000 C  ", "len 34", AtomSymbol::Element(Element::C),
       None, 0, None, None, None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C   ", "len 35 padded", AtomSymbol::Element(Element::C),
       None, 0, None, None, None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C  -2", "len 36", AtomSymbol::Element(Element::C),
       Some(10), 0, None, None, None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C  -2  ", "len 38 padded", AtomSymbol::Element(Element::C),
       Some(10), 0, None, None, None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3", "len 39", AtomSymbol::Element(Element::C),
       Some(10), 1, None, None, None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  1", "len 42 with stereo parity", AtomSymbol::Element(Element::C),
       Some(10), 1, None, Some(AtomStereoParity::Odd), None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4", "len 51 with query fields", AtomSymbol::Element(Element::C),
       Some(10), 1, Some(4), Some(AtomStereoParity::Odd), Some(1), Some(AtomStereoCare::Care), None, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4 ", "len 52 padded", AtomSymbol::Element(Element::C),
       Some(10), 1, Some(4), Some(AtomStereoParity::Either), None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0           1  2  1", "len 69 with reaction fields", AtomSymbol::Element(Element::C),
       None, 0, None, Some(AtomStereoParity::Either), None, None, Some(1), Some(AtomInversionRetention::Retained), Some(AtomExactChange::Match), 0.0000, 0.0000, 0.0000)]
#[case(b"    1.0000    2.0000    3.0000 L   0  0  0  0  0  4", "atom list", AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }),
       None, 0, Some(4), Some(AtomStereoParity::Either), None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0", "named isotope", AtomSymbol::NamedIsotope(NamedIsotope::D),
       Some(2), 1, Some(1), Some(AtomStereoParity::Either), None, None, Some(0), None, None, 1.2345, 2.3456, 3.4567)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0XXX  0  0  0", "non-strict padding", AtomSymbol::Element(Element::C),
       Some(10), 1, Some(4), Some(AtomStereoParity::Either), None, None, Some(0), None, None, 1.2345, 2.3456, 3.4567)]
fn test_atomlike_input(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: i8,
    #[case] valence: Option<u8>,
    #[case] stereo_parity: Option<AtomStereoParity>,
    #[case] hydrogen_count: Option<u8>,
    #[case] stereo_care: Option<AtomStereoCare>,
    #[case] atom_map_num: Option<u32>,
    #[case] inversion_retention: Option<AtomInversionRetention>,
    #[case] exact_change: Option<AtomExactChange>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let mut parser = atomlike_input(ParseFlags::LENIENT);
    let result = parser.parse(input);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.symbol, symbol,
        "{} has returned symbol {:?}, expected {:?}",
        desc, atom.symbol, symbol,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{} has returned isotope mass {:?}, expected {:?}",
        desc, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, charge,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, charge,
    );
    assert_eq!(
        atom.valence, valence,
        "{} has returned valence {:?}, expected {:?}",
        desc, atom.valence, valence,
    );
    assert_eq!(
        atom.stereo_parity, stereo_parity,
        "{} has returned stereo_parity {:?}, expected {:?}",
        desc, atom.stereo_parity, stereo_parity,
    );
    assert_eq!(
        atom.hydrogen_count, hydrogen_count,
        "{} has returned hydrogen_count {:?}, expected {:?}",
        desc, atom.hydrogen_count, hydrogen_count,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{} has returned stereo_care {:?}, expected {:?}",
        desc, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.atom_map_num, atom_map_num,
        "{} has returned atom_map_num {:?}, expected {:?}",
        desc, atom.atom_map_num, atom_map_num,
    );
    assert_eq!(
        atom.inversion_retention, inversion_retention,
        "{} has returned inversion_retention {:?}, expected {:?}",
        desc, atom.inversion_retention, inversion_retention,
    );
    assert_eq!(
        atom.exact_change, exact_change,
        "{} has returned exact_change {:?}, expected {:?}",
        desc, atom.exact_change, exact_change,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().x, x_position),
        "{} has returned x {:?}, expected {:?}",
        desc,
        atom.position.unwrap().x,
        x_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().y, y_position),
        "{} has returned y {:?}, expected {:?}",
        desc,
        atom.position.unwrap().y,
        y_position,
    );
    assert!(
        approx_eq!(f64, atom.position.unwrap().z, z_position),
        "{} has returned z {:?}, expected {:?}",
        desc,
        atom.position.unwrap().z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  4", "invalid stereo parity", error::ErrorKind::MapRes)]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  0  a", "non-numeric hydrogen count", error::ErrorKind::Digit)]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  a", "non-numeric valence", error::ErrorKind::Digit)]
#[case(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0XXX  0  0  0", "non-strict padding", error::ErrorKind::Verify)]
fn test_atomlike_input_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = atomlike_input(ParseFlags::STRICT);
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
#[case(b"    1.0000    2.0000    3.0000 C  -2 3", "len 38")]
fn test_atomlike_input_partial_fields(#[case] input: &[u8], #[case] desc: &str) {
    let mut parser = atomlike_input(ParseFlags::LENIENT);
    let result = parser.parse(input);
    assert!(result.is_err(), "{} should have failed", desc);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == error::ErrorKind::Eof),
        "{} should have failed with error kind {:?}, got {:?}",
        desc,
        error::ErrorKind::Eof,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rustfmt::skip]
#[rstest]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4   \t", "len 55")]
#[case(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0           ", "len 80")]
fn test_atomlike_input_whitespace_padded(#[case] input: &[u8], #[case] desc: &str) {
    let mut parser = atomlike_input(ParseFlags::LENIENT);
    let result = parser.parse(input);
    assert!(result.is_ok(), "{} should have succeeded", desc);
    let (remaining, atom) = result.unwrap();
    assert!(remaining.is_empty(), "{} has non-empty remaining", desc);
    assert_eq!(
        atom.charge, 1,
        "{} has returned charge {:?}, expected {:?}",
        desc, atom.charge, 1
    );
    assert_eq!(
        atom.valence,
        Some(4),
        "{} has returned valence {:?}, expected {:?}",
        desc,
        atom.valence,
        Some(4)
    );
}
