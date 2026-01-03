use bstr::ByteSlice;
use float_cmp::*;
use nom::error::ErrorKind as NomErrorKind;
use nom::Err;
use pretty_assertions::assert_eq;
use rstest::*;

use super::*;
use crate::io::ctfile::config::CtabParseFlags;
use crate::table_ir::{
    AtomExactChange, AtomInversionRetention, AtomStereoCare, AtomStereoParity, RGroup, WildcardAtom,
};

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
fn test_atom_symbol(#[case] input: &[u8], #[case] expected: AtomSymbol) {
    let result = atom_symbol(CtabParseFlags::BASIC).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, symbol) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} has non-empty remaining",
        input_str
    );
    assert_eq!(symbol, expected);
}

#[rstest]
#[case::empty(b"", NomErrorKind::Eof)]
#[case::blank(b"   ", NomErrorKind::Eof)]
#[case::element_invalid_1(b"Xx ", NomErrorKind::MapRes)]
#[case::element_invalid_2(b"LQ ", NomErrorKind::MapRes)]
#[case::wildcard_atom_a(b"A  ", NomErrorKind::MapRes)]
#[case::chemaxon_wildcard_atom(b"QH ", NomErrorKind::MapRes)]
#[case::atom_list(b"L  ", NomErrorKind::MapRes)]
#[case::lone_pair(b"LP ", NomErrorKind::MapRes)]
#[case::rgroup(b"R1 ", NomErrorKind::MapRes)]
#[case::pseudoatom_al(b"Ala", NomErrorKind::MapRes)]
#[case::pseudoatom_unicode(b"\xCE\xB1 ", NomErrorKind::MapRes)]
fn test_atom_symbol_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = atom_symbol(CtabParseFlags::BASIC).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::named_isotope_d(b"D  ", NomErrorKind::MapRes)]
#[case::named_isotope_d_lowercase(b"d  ", NomErrorKind::MapRes)]
#[case::named_isotope_t(b"T  ", NomErrorKind::MapRes)]
fn test_atom_symbol_minimal_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = atom_symbol(CtabParseFlags::MINIMAL).parse(input);
    assert!(result.is_err(), "{:?} should have failed", input);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
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
fn test_extended_atom_symbol(#[case] input: &[u8], #[case] expected: AtomSymbol) {
    let result = extended_atom_symbol(CtabParseFlags::EXTENDED).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, symbol) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} has non-empty remaining",
        input_str
    );
    assert_eq!(
        symbol, expected,
        "{:?} has returned symbol {:?}, expected {:?}",
        input, symbol, expected
    );
}

#[rstest]
#[case::empty(b"", NomErrorKind::Eof)]
#[case::blank(b"   ", NomErrorKind::Eof)]
#[case::element_invalid_1(b"Xx ", NomErrorKind::MapRes)]
#[case::element_invalid_2(b"LQ ", NomErrorKind::MapRes)]
#[case::chemaxon_wildcard_atom(b"QH ", NomErrorKind::MapRes)]
#[case::pseudoatom_ala(b"Ala", NomErrorKind::MapRes)]
#[case::pseudoatom_unicode(b"\xCE\xB1 ", NomErrorKind::MapRes)]
fn test_extended_atom_symbol_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let result = extended_atom_symbol(CtabParseFlags::EXTENDED).parse(input);
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

#[rstest]
#[case::named_isotope_d(b"D  ", NomErrorKind::MapRes)]
#[case::named_isotope_d_lowercase(b"d  ", NomErrorKind::MapRes)]
#[case::named_isotope_t(b"T  ", NomErrorKind::MapRes)]
#[case::wildcard_atom_a(b"A  ", NomErrorKind::MapRes)]
#[case::wildcard_atom_q(b"Q  ", NomErrorKind::MapRes)]
#[case::wildcard_atom_star(b"*  ", NomErrorKind::MapRes)]
#[case::atom_list(b"L  ", NomErrorKind::MapRes)]
#[case::lone_pair(b"LP ", NomErrorKind::MapRes)]
#[case::rgroup(b"R  ", NomErrorKind::MapRes)]
#[case::rgroup_unlabeled(b"R# ", NomErrorKind::MapRes)]
#[case::rgroup_r1(b"R1 ", NomErrorKind::MapRes)]
#[case::rgroup_r3(b"R3 ", NomErrorKind::MapRes)]
#[case::pseudoatom_ala(b"Ala", NomErrorKind::MapRes)]
#[case::pseudoatom_unicode(b"\xCE\xB1 ", NomErrorKind::MapRes)]
fn test_extended_atom_symbol_minimal_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = extended_atom_symbol(CtabParseFlags::MINIMAL).parse(input);
    assert!(result.is_err(), "{:?} should have failed", input);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
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
    let result = extended_atom_symbol(CtabParseFlags::LENIENT).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, symbol) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(symbol, expected, "{:?} has returned symbol {:?}, expected {:?}", input, symbol, expected);
}

#[rstest]
#[case::empty(b"", NomErrorKind::Eof)]
#[case::blank(b"   ", NomErrorKind::Eof)]
#[case::element_invalid_1(b"Xx ", NomErrorKind::MapRes)]
#[case::element_invalid_2(b"LQ ", NomErrorKind::MapRes)]
#[case::pseudoatom_ala(b"Ala", NomErrorKind::MapRes)]
#[case::pseudoatom_unicode(b"\xCE\xB1 ", NomErrorKind::MapRes)]
fn test_extended_atom_symbol_lenient_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = extended_atom_symbol(CtabParseFlags::LENIENT).parse(input);
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

#[rstest]
#[case::pseudoatom_ala(b"Ala", AtomSymbol::Pseudoatom("Ala".to_string()))]
fn test_extended_atom_symbol_pseudoatoms(#[case] input: &[u8], #[case] expected: AtomSymbol) {
    let result = extended_atom_symbol(CtabParseFlags::PSEUDOATOMS).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, symbol) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} has non-empty remaining",
        input_str
    );
    assert_eq!(
        symbol, expected,
        "{:?} has returned symbol {:?}, expected {:?}",
        input, symbol, expected
    );
}

#[rstest]
#[case::reserved_named_isotope(b"D  ", NomErrorKind::MapRes)]
#[case::reserved_wildcard_atom(b"A  ", NomErrorKind::MapRes)]
#[case::reserved_chemaxon_wildcard_atom(b"QH ", NomErrorKind::MapRes)]
#[case::reserved_atom_list(b"L  ", NomErrorKind::MapRes)]
#[case::reserved_lone_pair(b"LP ", NomErrorKind::MapRes)]
#[case::reserved_rgroup(b"R  ", NomErrorKind::MapRes)]
#[case::reserved_rgroup_unlabeled(b"R# ", NomErrorKind::MapRes)]
#[case::reserved_rgroup_r1(b"R1 ", NomErrorKind::MapRes)]
#[case::pseudoatom_unicode(b"\xCE\xB1 ", NomErrorKind::MapRes)]
fn test_extended_atom_symbol_pseudoatoms_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = extended_atom_symbol(CtabParseFlags::PSEUDOATOMS).parse(input);
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
#[case::len_69(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_lower_bound(b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4  0  0  0  0  0  0",
       Element::C, Some(9), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_upper_bound(b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4  0  0  0  0  0  0",
       Element::C, Some(16), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_out_of_range_low(b"    1.2345    2.3456    3.4567 C  -4  3  0  0  0  4  0  0  0  0  0  0",
       Element::C, None, Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_out_of_range_high(b"    1.2345    2.3456    3.4567 C   5  3  0  0  0  4  0  0  0  0  0  0",
       Element::C, None, Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::charge_out_of_range_high(b"    1.2345    2.3456    3.4567 C  -2  8  0  0  0  4  0  0  0  0  0  0",
       Element::C, Some(10), None, Some(4), 1.2345, 2.3456, 3.4567)]
#[case::blank_block_1(b"    1.2345    2.3456    3.4567 C  -2  3  0     0  4  0  0  0  0  0  0",
       Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::blank_block_2(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4     0  0  0  0  0",
       Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::blank_block_3(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0      ",
       Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::gaps_with_spaces_and_zeros(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4           1 0    ",
       Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0",
       Element::H, Some(2), Some(1), Some(1), 1.2345, 2.3456, 3.4567)]
#[case::non_strict_padding(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4  0  0  0  0  0  0",
       Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
fn test_atom_input69(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input69(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_numeric_coordinate(b"    1.234a    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0", NomErrorKind::Eof)]
#[case::non_numeric_atom_map_number(b"    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  a  0  0", NomErrorKind::Digit)]
#[case::atom_list(b"    1.2345    2.3456    3.4567 L   0  0  0  0  0  0  0  0  0  0  0  0", NomErrorKind::MapRes)]
fn test_atom_input69_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input69(input, CtabParseFlags::BASIC);
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
#[case::len_69(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_lower_bound(b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4  0  0  0  0  0  0",
       Element::C, Some(9), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_upper_bound(b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4  0  0  0  0  0  0",
       Element::C, Some(16), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0",
       Element::H, Some(2), Some(1), Some(1), 1.2345, 2.3456, 3.4567)]
fn test_atom_input69_strict(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input69(input, CtabParseFlags::BASIC & CtabParseFlags::STRICT);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_strict_padding(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4  0  0  0  0  0  0", NomErrorKind::Verify)]
fn test_atom_input69_strict_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input69(input, CtabParseFlags::BASIC &  CtabParseFlags::STRICT);
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
#[case::len_69(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0", Element::C, Some(10), Some(1), Some(4))]
#[case::zero_coordinates(b"    0.0000    0.0000    0.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0", Element::C, Some(10), Some(1), Some(4))]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0", Element::H, Some(2), Some(1), Some(1))]
fn test_atom_input69_ignore_positions(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
) {
    let result = atom_input69(input, CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
#[case::invalid_coordinates(b"    invalid    invalid    invalid C  -2  3  0  0  0  4  0  0  0  0  0  0", NomErrorKind::Tag)]
fn test_atom_input69_ignore_positions_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input69(
        input,
        CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS,
    );
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
#[case::len_51(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4", Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_lower_bound(b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4", Element::C, Some(9), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_upper_bound(b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4", Element::C, Some(16), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_out_of_range_low(b"    1.2345    2.3456    3.4567 C  -4  3  0  0  0  4", Element::C, None, Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_out_of_range_high(b"    1.2345    2.3456    3.4567 C   5  3  0  0  0  4", Element::C, None, Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::charge_out_of_range_high(b"    1.2345    2.3456    3.4567 C  -2  8  0  0  0  4",  Element::C, Some(10), None, Some(4), 1.2345, 2.3456, 3.4567)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1", Element::H, Some(2), Some(1), Some(1), 1.2345, 2.3456, 3.4567)]
#[case::non_strict_padding(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4", Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
fn test_atom_input51(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input51(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_numeric_coordinate(b"    1.234a    2.3456    3.4567 C  -2  3  0  0  0  4", NomErrorKind::Eof)]
#[case::non_numeric_valence(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  a", NomErrorKind::Digit)]
#[case::out_of_range_valence(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0 16", NomErrorKind::Verify)]
#[case::invalid_atom_symbol(b"    1.2345    2.3456    3.4567 L  -2  3  0  0  0  4", NomErrorKind::MapRes)]
fn test_atom_input51_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input51(input, CtabParseFlags::BASIC);
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
#[case::len_51(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4", Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_lower_bound(b"    1.2345    2.3456    3.4567 C  -3  3  0  0  0  4", Element::C, Some(9), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_upper_bound(b"    1.2345    2.3456    3.4567 C   4  3  0  0  0  4", Element::C, Some(16), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1", Element::H, Some(2), Some(1), Some(1), 1.2345, 2.3456, 3.4567)]
fn test_atom_input51_strict(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input51(input, CtabParseFlags::BASIC & CtabParseFlags::STRICT);
    let input_str = input.to_str_lossy();   
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_strict_padding(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4", NomErrorKind::Verify)]
fn test_atom_input51_strict_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input51(input, CtabParseFlags::BASIC & CtabParseFlags::STRICT);
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
#[case::len_51(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4", Element::C, Some(10), Some(1), Some(4))]
#[case::zero_coordinates(b"    0.0000    0.0000    0.0000 C  -2  3  0  0  0  4", Element::C, Some(10), Some(1), Some(4))]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1", Element::H, Some(2), Some(1), Some(1))]
fn test_atom_input51_ignore_positions(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
) {
    let result = atom_input51(
        input,
        CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS,
    );
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} has non-empty remaining",
        input_str
    );
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
#[case::invalid_coordinates(b"    invalid    invalid    invalid C  -2  3  0  0  0  4", NomErrorKind::Tag)]
fn test_atom_input51_ignore_positions_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input51(input, CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS);
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
#[case::blank_mass_diff(b"    1.2345    2.3456    3.4567 C      3  0  0  0  4", None, Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
#[case::blank_charge(b"    1.2345    2.3456    3.4567 C  -2     0  0  0  4",  Some(10), None, Some(4), 1.2345, 2.3456, 3.4567)]
#[case::blank_valence(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0   ",  Some(10), Some(1), None, 1.2345, 2.3456, 3.4567)]
fn test_atom_input51_empty_fields(
    #[case] input: &[u8],
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input51(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_42(b"    1.2345    2.3456    3.4567 C   0  3  1", Element::C, None, Some(1), 1.2345, 2.3456, 3.4567)]
#[case::blank_stereo_parity(b"    1.2345    2.3456    3.4567 C   0  3   0  ", Element::C, None, Some(1), 1.2345, 2.3456, 3.4567)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0", Element::H, Some(2), Some(1), 1.2345, 2.3456, 3.4567)]
#[case::non_strict_padding(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0", Element::C, Some(10), Some(1), 1.2345, 2.3456, 3.4567)]
fn test_atom_input42(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input42(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::stereo_parity_out_of_range(b"    1.2345    2.3456    3.4567 C   0  3  4", NomErrorKind::Verify)]
#[case::non_numeric_coordinate(b"    1.234a    2.3456    3.4567 C   0  3  0", NomErrorKind::Eof)]
#[case::invalid_atom_symbol(b"    1.2345    2.3456    3.4567 L   0  3  0", NomErrorKind::MapRes)]
#[case::non_numeric_stereo_parity(b"    1.2345    2.3456    3.4567 C   0  3  a", NomErrorKind::Digit)]
fn test_atom_input42_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input42(input, CtabParseFlags::BASIC);
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
#[case::len_42(b"    1.2345    2.3456    3.4567 C   0  3  1", Element::C, None, Some(1), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff(b"    1.2345    2.3456    3.4567 C  -2  3  0", Element::C, Some(10), Some(1), 1.2345, 2.3456, 3.4567)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0", Element::H, Some(2), Some(1), 1.2345, 2.3456, 3.4567)]
fn test_atom_input42_strict(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input42(input, CtabParseFlags::BASIC & CtabParseFlags::STRICT);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_strict_padding(b"    1.2345    2.3456    3.4567 C  -2  3     a", NomErrorKind::Verify)]
fn test_atom_input42_strict_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input42(input, CtabParseFlags::BASIC & CtabParseFlags::STRICT);
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
#[case::len_42(b"    1.2345    2.3456    3.4567 C  -2  3  0", Element::C, Some(10), Some(1))]
#[case::zero_coordinates(b"    0.0000    0.0000    0.0000 C  -2  3  0", Element::C, Some(10), Some(1))]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0", Element::H, Some(2), Some(1))]
fn test_atom_input42_ignore_positions(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
) {
    let result = atom_input42(input, CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS);
    let input_str = input.to_str_lossy();       
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
#[case::invalid_coordinates(b"    invalid    invalid    invalid C  -2  3  0", NomErrorKind::Tag)]
fn test_atom_input42_ignore_positions_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input69(input, CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS);
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
#[case::blank_stereo_parity(b"    1.2345    2.3456    3.4567 C   0  3   ", Element::C, None, Some(1), 1.2345, 2.3456, 3.4567)]
#[case::blank_mass_diff(b"    1.2345    2.3456    3.4567 C      3  0", Element::C, None, Some(1), 1.2345, 2.3456, 3.4567)]
fn test_atom_input42_empty_fields(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input42(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_39(b"    1.2345    2.3456    3.4567 C  -2  3", Element::C, Some(10), Some(1), 1.2345, 2.3456, 3.4567)]
#[case::charge_out_of_range_high(b"    1.2345    2.3456    3.4567 C  -2  8",  Element::C, Some(10), None, 1.2345, 2.3456, 3.4567)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3",  Element::H, Some(2), Some(1), 1.2345, 2.3456, 3.4567)]
fn test_atom_input39(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input39(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
        atom.valence, None,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, None as Option<u8>
    );
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_numeric_coordinate(b"    1.234a    2.3456    3.4567 C  -2  3", NomErrorKind::Eof)]
#[case::non_numeric_mass_diff(b"    1.2345    2.3456    3.4567 C  -a  3", NomErrorKind::Digit)]
#[case::invalid_atom_symbol(b"    1.2345    2.3456    3.4567 L  -2  3", NomErrorKind::MapRes)]
fn test_atom_input39_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input39(input, CtabParseFlags::BASIC);
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
#[case::len_39(b"    1.2345    2.3456    3.4567 C  -2  3", Element::C, Some(10), Some(1))]
#[case::zero_coordinates(b"    0.0000    0.0000    0.0000 C  -2  3", Element::C, Some(10), Some(1))]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3", Element::H, Some(2), Some(1))]
fn test_atom_input39_ignore_positions(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
) {
    let result = atom_input39(input, CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
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
#[case::invalid_coordinates(b"    invalid    invalid    invalid C  -2  3", NomErrorKind::Tag)]
fn test_atom_input39_ignore_positions_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input39(input, CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS);
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
#[case::blank_mass_diff(b"    1.2345    2.3456    3.4567 C      3", None, Some(1))]
#[case::blank_charge(b"    1.2345    2.3456    3.4567 C  -2   ", Some(10), None)]
fn test_atom_input39_empty_fields(
    #[case] input: &[u8],
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
) {
    let result = atom_input39(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, _position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
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
}

#[rustfmt::skip]
#[rstest]
#[case::len_36(b"    1.2345    2.3456    3.4567 C  -2", Element::C, Some(10), 1.2345, 2.3456, 3.4567)]
#[case::mass_diff_out_of_range_low(b"    1.2345    2.3456    3.4567 C  -4", Element::C, None, 1.2345, 2.3456, 3.4567)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2", Element::H, Some(2), 1.2345, 2.3456, 3.4567)]
fn test_atom_input36(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input36(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, None,
        "{:?} has returned charge {:?}, expected {:?}",
        input_str, atom.charge, None as Option<i8>
    );
    assert_eq!(
        atom.valence, None,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, None as Option<u8>
    );
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_numeric_coordinate(b"    1.234a    2.3456    3.4567 C  -2", NomErrorKind::Eof)]
#[case::invalid_atom_symbol(b"    1.2345    2.3456    3.4567 L  -2", NomErrorKind::MapRes)]
fn test_atom_input36_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input36(input, CtabParseFlags::BASIC);
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
#[case::len_36(b"    1.2345    2.3456    3.4567 C  -2", Element::C, Some(10))]
#[case::zero_coordinates(b"    0.0000    0.0000    0.0000 C  -2", Element::C, Some(10))]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2", Element::H, Some(2))]
fn test_atom_input36_ignore_positions(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
) {
    let result = atom_input36(input,  CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass,
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
#[case::invalid_coordinates(b"    invalid    invalid    invalid C  -2", NomErrorKind::Tag)]
fn test_atom_input36_ignore_positions_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input39(input, CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS);
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

#[rstest]
#[case::blank_mass_diff(b"    1.2345    2.3456    3.4567 C    ", None)]
fn test_atom_input36_empty_fields(#[case] input: &[u8], #[case] isotope_mass: Option<u32>) {
    let result = atom_input36(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, _position)) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} has non-empty remaining",
        input_str
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_34(b"    1.2345    2.3456    3.4567 C  ", Element::C, None, 1.2345, 2.3456, 3.4567)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  ", Element::H, Some(2), 1.2345, 2.3456, 3.4567)]
fn test_atom_input34(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input34(input, CtabParseFlags::BASIC);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass,
    );
    assert_eq!(
        atom.charge, None,
        "{:?} has returned charge {:?}, expected {:?}",
        input_str, atom.charge, None as Option<i8>
    );
    assert_eq!(
        atom.valence, None,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, None as Option<u8>
    );
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_numeric_coordinate(b"    1.234a    2.3456    3.4567 C  ", NomErrorKind::Eof)]
#[case::invalid_atom_symbol(b"    1.2345    2.3456    3.4567 L  ", NomErrorKind::MapRes)]
fn test_atom_input34_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input34(input, CtabParseFlags::BASIC);
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
#[case::len_34(b"    1.2345    2.3456    3.4567 C  ", Element::C, None)]
#[case::zero_coordinates(b"    0.0000    0.0000    0.0000 C  ", Element::C, None)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  ", Element::H, Some(2))]
fn test_atom_input34_ignore_positions(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
) {
    let result = atom_input34(input, CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
    assert_eq!(
        atom.element, element,
        "{:?} has returned element {:?}, expected {:?}",
        input_str, atom.element, element,
    );
    assert_eq!(
        atom.isotope_mass, isotope_mass,
        "{:?} has returned isotope mass {:?}, expected {:?}",
        input_str, atom.isotope_mass, isotope_mass,
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
#[case::invalid_coordinates(b"    invalid    invalid    invalid C  ", NomErrorKind::Tag)]
fn test_atom_input34_ignore_positions_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input51(input, CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS);
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
#[case::len_34(b"    1.0000    2.0000    3.0000 C  ", Element::C, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::len_35_padded(b"    1.0000    2.0000    3.0000 C   ", Element::C, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::len_36(b"    1.0000    2.0000    3.0000 C  -2", Element::C, Some(10), None, None, 1.0000, 2.0000, 3.0000)]
#[case::len_39(b"    1.0000    2.0000    3.0000 C  -2  3", Element::C, Some(10), Some(1), None, 1.0000, 2.0000, 3.0000)]
#[case::len_42_with_stereo_parity(b"    1.0000    2.0000    3.0000 C  -2  3  1", Element::C, Some(10), Some(1), None, 1.0000, 2.0000, 3.0000)]
#[case::len_51_with_query_fields(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4", Element::C, Some(10), Some(1), Some(4), 1.0000, 2.0000, 3.0000)]
fn test_atom_input(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input(CtabParseFlags::BASIC).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} should consume all input", input_str);

    assert_eq!(
        atom.element, element,
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
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence
    );
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_numeric_coordinate(b"    1.234a    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0", NomErrorKind::Eof)]
#[case::atom_list(b"    0.7145    2.0625    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0", NomErrorKind::MapRes)]
#[case::pseudoatom(b"   -1.8857    2.4750    0.0000 Psd 0  0  0  0  0  0  0  0  0  0  0  0", NomErrorKind::MapRes)]
#[case::non_numeric_valence(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  a", NomErrorKind::Digit)]
#[case::incorrect_yz_coordinates(b"    0.1   0.0    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0", NomErrorKind::Eof)]
fn test_atom_input_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input(CtabParseFlags::BASIC).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert! (
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_34(b"    1.0000    2.0000    3.0000 C  ", Element::C, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::len_36(b"    1.0000    2.0000    3.0000 C  -2", Element::C, Some(10), None, None, 1.0000, 2.0000, 3.0000)]
#[case::len_39(b"    1.0000    2.0000    3.0000 C  -2  3", Element::C, Some(10), Some(1), None, 1.0000, 2.0000, 3.0000)]
#[case::len_42(b"    1.0000    2.0000    3.0000 C  -2  3  1", Element::C, Some(10), Some(1), None, 1.0000, 2.0000, 3.0000)]
#[case::len_51(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4", Element::C, Some(10), Some(1), Some(4), 1.0000, 2.0000, 3.0000)]
#[case::len_69(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0", Element::C, Some(10), Some(1), Some(4), 1.2345, 2.3456, 3.4567)]
fn test_atom_input_strict(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
    #[case] x_position: f64,
    #[case] y_position: f64,
    #[case] z_position: f64,
) {
    let result = atom_input(CtabParseFlags::BASIC & CtabParseFlags::STRICT).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} should consume all input", input_str);

    assert_eq!(
        atom.element, element,
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
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence
    );
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_strict_padding(b"    1.2345    2.3456    3.4567 C  -2  3  0XXX  0  4", NomErrorKind::Verify)]
fn test_atom_input_strict_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input(CtabParseFlags::BASIC & CtabParseFlags::STRICT).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert! (
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_34(b"    1.0000    2.0000    3.0000 C  ", Element::C, None, None, None)]
#[case::len_36(b"    1.0000    2.0000    3.0000 C  -2", Element::C, Some(10), None, None)]
#[case::len_39(b"    1.0000    2.0000    3.0000 C  -2  3", Element::C, Some(10), Some(1), None)]
#[case::len_42(b"    1.0000    2.0000    3.0000 C  -2  3  1", Element::C, Some(10), Some(1), None)]
#[case::len_51(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4", Element::C, Some(10), Some(1), Some(4))]
#[case::len_69(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0", Element::C, Some(10), Some(1), Some(4))]
fn test_atom_input_ignore_positions(
    #[case] input: &[u8],
    #[case] element: Element,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
) {
    let result = atom_input(CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} should consume all input", input_str);

    assert_eq!(
        atom.element, element,
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
        atom.valence, valence,
        "{:?} has returned valence {:?}, expected {:?}",
        input_str, atom.valence, valence
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
#[case::invalid_coordinates(b"    invalid    invalid    invalid C  -2  3  0  0  0  4  0  0  0  0  0  0", NomErrorKind::Tag)]
fn test_atom_input_ignore_positions_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = atom_input69(
        input,
        CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS,
    );
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
#[case::len_69(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), Some(AtomStereoParity::Either), None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::stereo_parity_and_query_fields(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4  0  0  0  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), Some(AtomStereoParity::Odd), Some(1), Some(AtomStereoCare::Care), None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::atom_list(b"    1.0000    2.0000    3.0000 L   0  0  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }), None, None, Some(4), Some(AtomStereoParity::Either), None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0",
       AtomSymbol::NamedIsotope(NamedIsotope::D), Some(2), Some(1), Some(1), Some(AtomStereoParity::Either), None, None, None, None, None, 1.2345, 2.3456, 3.4567)]
#[case::non_strict_padding(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0XXX  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), Some(AtomStereoParity::Either), None, None, None, None, None, 1.2345, 2.3456, 3.4567)]
fn test_extended_atom_input(
    #[case] input: &[u8],
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
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
    let result = extended_atom_input(CtabParseFlags::EXTENDED).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
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
        atom.stereo_parity, stereo_parity,
        "{:?} has returned stereo_parity {:?}, expected {:?}",
        input_str, atom.stereo_parity, stereo_parity,
    );
    assert_eq!(
        atom.hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.hydrogens, hydrogen_count,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{:?} has returned stereo_care {:?}, expected {:?}",
        input_str, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.atom_map_num, atom_map_num,
        "{:?} has returned atom_map_num {:?}, expected {:?}",
        input_str, atom.atom_map_num, atom_map_num,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
                input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::invalid_stereo_parity(b"    1.0000    2.0000    3.0000 C  -2  3  4", NomErrorKind::Verify)]
#[case::non_numeric_hydrogen_count(b"    1.0000    2.0000    3.0000 C  -2  3  0  a", NomErrorKind::Digit)]
#[case::non_numeric_valence(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  a", NomErrorKind::Digit)]
#[case::chemaxon_wildcard(b"    1.2345    2.3456    3.4567 AH 0  0  0  0  0  0  0  0  0  0  0  0", NomErrorKind::MapRes)]
#[case::pseudoatom(b"    1.2345    2.3456    3.4567 Ala 0  0  0  0  0  0  0  0  0  0  0  0", NomErrorKind::MapRes)]
fn test_extended_atom_input_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = extended_atom_input(CtabParseFlags::EXTENDED).parse(input);
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
#[case::len_69(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), Some(AtomStereoParity::Either), None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::stereo_parity_and_query_fields(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4  0  0  0  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), Some(AtomStereoParity::Odd), Some(1), Some(AtomStereoCare::Care), None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::reaction_fields(b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0           1  2  1",
       AtomSymbol::Element(Element::C), None, None, None, Some(AtomStereoParity::Either), None, None, Some(1), Some(AtomInversionRetention::Retained), Some(AtomExactChange::Match), 0.0000, 0.0000, 0.0000)]
#[case::atom_list(b"    1.0000    2.0000    3.0000 L   0  0  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }), None, None, Some(4), Some(AtomStereoParity::Either), None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
fn test_extended_atom_input_strict(
    #[case] input: &[u8],
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
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
    let result = extended_atom_input(CtabParseFlags::STRICT).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
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
        atom.stereo_parity, stereo_parity,
        "{:?} has returned stereo_parity {:?}, expected {:?}",
        input_str, atom.stereo_parity, stereo_parity,
    );
    assert_eq!(
        atom.hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.hydrogens, hydrogen_count,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{:?} has returned stereo_care {:?}, expected {:?}",
        input_str, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.atom_map_num, atom_map_num,
        "{:?} has returned atom_map_num {:?}, expected {:?}",
        input_str, atom.atom_map_num, atom_map_num,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::non_strict_padding(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0XXX  0  0  0", NomErrorKind::Verify)]
fn test_extended_atom_input_strict_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = extended_atom_input(CtabParseFlags::STRICT).parse(input);
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
#[case::len_69(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), Some(AtomStereoParity::Either), None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::stereo_parity_and_query_fields(b"    1.0000    2.0000    3.0000 C  -2  3  1  2  1  4  0  0  0  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), Some(AtomStereoParity::Odd), Some(1), Some(AtomStereoCare::Care), None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::atom_list(b"    1.0000    2.0000    3.0000 L   0  0  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }), None, None, Some(4), Some(AtomStereoParity::Either), None, None, None, None, None, 1.0000, 2.0000, 3.0000)]
#[case::named_isotope(b"    1.2345    2.3456    3.4567 D  -2  3  0  0  0  1  0  0  0  0  0  0",
       AtomSymbol::NamedIsotope(NamedIsotope::D), Some(2), Some(1), Some(1), Some(AtomStereoParity::Either), None, None, None, None, None, 1.2345, 2.3456, 3.4567)]
#[case::chemaxon_wildcard(b"    1.2345    2.3456    3.4567 AH  0  0  0  0  0  0  0  0  0  0  0  0",
       AtomSymbol::WildcardAtom(WildcardAtom::HeavyOrH), None, None, None, Some(AtomStereoParity::Either), None, None, None, None, None, 1.2345, 2.3456, 3.4567)]
#[case::non_strict_padding(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0XXX  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), Some(AtomStereoParity::Either), None, None, None, None, None, 1.2345, 2.3456, 3.4567)]
fn test_extended_atom_input_lenient(
    #[case] input: &[u8],
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
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
    let result = extended_atom_input(CtabParseFlags::LENIENT).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
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
        atom.stereo_parity, stereo_parity,
        "{:?} has returned stereo_parity {:?}, expected {:?}",
        input_str, atom.stereo_parity, stereo_parity,
    );
    assert_eq!(
        atom.hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.hydrogens, hydrogen_count,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{:?} has returned stereo_care {:?}, expected {:?}",
        input_str, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.atom_map_num, atom_map_num,
        "{:?} has returned atom_map_num {:?}, expected {:?}",
        input_str, atom.atom_map_num, atom_map_num,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::invalid_stereo_parity(b"    1.0000    2.0000    3.0000 C  -2  3  4", NomErrorKind::Verify)]
#[case::non_numeric_hydrogen_count(b"    1.0000    2.0000    3.0000 C  -2  3  0  a", NomErrorKind::Digit)]
#[case::non_numeric_valence(b"    1.0000    2.0000    3.0000 C  -2  3  0  0  0  a", NomErrorKind::Digit)]
fn test_extended_atom_input_lenient_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = extended_atom_input(CtabParseFlags::LENIENT).parse(input);
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
#[case::len_69(b"    1.2345    2.3456    3.4567 Ala 0  0  0  0  0  0  0  0  0  0  0  0",
       AtomSymbol::Pseudoatom(String::from("Ala")), None, None, None, Some(AtomStereoParity::Either), None, None, None, None, None, 1.2345, 2.3456, 3.4567)]
fn test_extended_atom_input_pseudoatoms(
    #[case] input: &[u8],
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
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
    let flags = CtabParseFlags::EXTENDED | CtabParseFlags::PSEUDOATOMS;
    let result = extended_atom_input(flags).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
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
        atom.stereo_parity, stereo_parity,
        "{:?} has returned stereo_parity {:?}, expected {:?}",
        input_str, atom.stereo_parity, stereo_parity,
    );
    assert_eq!(
        atom.hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.hydrogens, hydrogen_count,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{:?} has returned stereo_care {:?}, expected {:?}",
        input_str, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.atom_map_num, atom_map_num,
        "{:?} has returned atom_map_num {:?}, expected {:?}",
        input_str, atom.atom_map_num, atom_map_num,
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
    assert!(
        approx_eq!(f64, position.x, x_position),
        "{:?} has returned x {:?}, expected {:?}",
        input_str,
        position.x,
        x_position,
    );
    assert!(
        approx_eq!(f64, position.y, y_position),
        "{:?} has returned y {:?}, expected {:?}",
        input_str,
        position.y,
        y_position,
    );
    assert!(
        approx_eq!(f64, position.z, z_position),
        "{:?} has returned z {:?}, expected {:?}",
        input_str,
        position.z,
        z_position,
    );
}

#[rustfmt::skip]
#[rstest]
#[case::len_69(b"    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), Some(AtomStereoParity::Either), None, None, None, None, None)]
#[case::zero_coordinates(b"    0.0000    0.0000    0.0000 C  -2  3  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::Element(Element::C), Some(10), Some(1), Some(4), Some(AtomStereoParity::Either), None, None, None, None, None)]
#[case::atom_list(b"    1.2345    2.3456    3.4567 L   0  0  0  0  0  4  0  0  0  0  0  0",
       AtomSymbol::AtomList(AtomList { elements: vec![], exclusion: false }), None, None, Some(4), Some(AtomStereoParity::Either), None, None, None, None, None)]
#[case::chemaxon_wildcard(b"    1.2345    2.3456    3.4567 AH  0  0  0  0  0  0  0  0  0  0  0  0",
       AtomSymbol::WildcardAtom(WildcardAtom::HeavyOrH), None, None, None, Some(AtomStereoParity::Either), None, None, None, None, None)]
fn test_extended_atom_input_ignore_positions(
    #[case] input: &[u8],
    #[case] symbol: AtomSymbol,
    #[case] isotope_mass: Option<u32>,
    #[case] charge: Option<i8>,
    #[case] valence: Option<u8>,
    #[case] stereo_parity: Option<AtomStereoParity>,
    #[case] hydrogen_count: Option<u8>,
    #[case] stereo_care: Option<AtomStereoCare>,
    #[case] atom_map_num: Option<u32>,
    #[case] inversion_retention: Option<AtomInversionRetention>,
    #[case] exact_change: Option<AtomExactChange>,
) {
    let flags = CtabParseFlags::LENIENT | CtabParseFlags::IGNORE_POSITIONS;
    let result = extended_atom_input(flags).parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, (atom, position)) = result.unwrap();
    assert!(remaining.is_empty(), "{:?} has non-empty remaining", input_str);
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
        atom.stereo_parity, stereo_parity,
        "{:?} has returned stereo_parity {:?}, expected {:?}",
        input_str, atom.stereo_parity, stereo_parity,
    );
    assert_eq!(
        atom.hydrogens, hydrogen_count,
        "{:?} has returned hydrogen_count {:?}, expected {:?}",
        input_str, atom.hydrogens, hydrogen_count,
    );
    assert_eq!(
        atom.stereo_care, stereo_care,
        "{:?} has returned stereo_care {:?}, expected {:?}",
        input_str, atom.stereo_care, stereo_care,
    );
    assert_eq!(
        atom.atom_map_num, atom_map_num,
        "{:?} has returned atom_map_num {:?}, expected {:?}",
        input_str, atom.atom_map_num, atom_map_num,
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
#[case::invalid_coordinates(b"    invalid    invalid    invalid C  -2  3  0  0  0  4  0  0  0  0  0  0", NomErrorKind::Tag)]
fn test_extended_atom_input_ignore_positions_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let result = extended_atom_input(CtabParseFlags::BASIC | CtabParseFlags::IGNORE_POSITIONS).parse(input);
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
