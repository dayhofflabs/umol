//! Convert numerical codes used in MOL files to enums
//!
//! All functions return table_ir types.

use umol_data::{Element, Isotope};

use crate::io::ctfile::error::ParseError;
use crate::table_ir::{
    AtomExactChange, AtomInversionRetention, AtomStereoCare, AtomSymbol, AttachmentPointType,
    BondDirection, BondOrder, BondReactingCenter, BondStereo, BondTopology, Chirality,
    RingBondCount, SubstitutionCount, UnsaturatedAtom,
};

/// Convert atom mass difference code (atom block)
/// 'dd' field: mass difference (-3..=4), None if 0 or value outside of this range
pub(super) fn convert_atom_mass_diff_code(code: i8) -> Option<i8> {
    match code {
        -3..=-1 | 1..=4 => Some(code),
        _ => None,
    }
}

/// Convert atom symbol and mass difference to element and isotope mass
/// 'ss' field: atom symbol, 'dd' field: mass difference
/// Processes elements and named isotopes.
/// Returns error for extended atom symbols (L, A, Q, *, LP, R#)
pub(super) fn convert_atom_symbol_mass_diff(
    symbol: AtomSymbol,
    mass_diff: Option<i8>,
) -> (Element, Option<u32>) {
    let (element, isotope_mass) = match symbol {
        AtomSymbol::Element(e) => {
            let isotope_mass = match mass_diff {
                None | Some(0) => None,
                Some(diff) => Some((e.reference_mass_number() as i8 + diff) as u32),
            };
            (e, isotope_mass)
        }
        AtomSymbol::NamedIsotope(i) => (i.element(), Some(i.mass_number())),
        _ => unreachable!("atom_symbol() should only return Element or NamedIsotope"),
    };
    (element, isotope_mass)
}

/// Convert atom charge code (includes doublet radical).
/// 'ccc' field: 0 = uncharged, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
/// 0 if outside of range.
/// Returns (charge, unpaired_e_count): for code 4, returns (None, Some(1)) for doublet radical
pub(super) fn convert_atom_charge_code(code: u8) -> (Option<i8>, Option<u8>) {
    match code {
        1..=3 | 5..=7 => (Some(4 - code as i8), None),
        4 => (None, Some(1)), // Doublet radical: 1 unpaired electron
        _ => (None, None),
    }
}

/// Validate atom isotope mass number (property block)
/// 'mmm' field: isotope mass number
pub(super) fn convert_atom_isotope_mass_number(
    element: Element,
    mass_number: u32,
    extended_isotopes: bool,
) -> Result<Option<u32>, ParseError> {
    if mass_number == 0 {
        return Ok(None);
    }
    if extended_isotopes || Isotope::is_catalogued(element, mass_number) {
        Ok(Some(mass_number))
    } else {
        Err(ParseError::InvalidIsotopeMass {
            mass: mass_number,
            element,
        })
    }
}

/// Convert atom stereo parity code to Chirality.
/// 'sss' field: 0 = not stereo, 1 = odd (clockwise), 2 = even (counter-clockwise), 3 = either
pub(super) fn convert_atom_stereo_parity_code(
    code: u8,
) -> Result<Option<Chirality>, ParseError> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(Chirality::Clockwise)),
        2 => Ok(Some(Chirality::CounterClockwise)),
        3 => Ok(Some(Chirality::Unspecified)),
        _ => Err(ParseError::InvalidStereoParity(code)),
    }
}

/// Convert atom hydrogen count code (extension: 0 in non-query atoms).
/// 'hhh' field: 0 = non-query atom, 1 = H0, 2 = H1, 3 = H2, 4 = H3, 5 = H4
/// Extended range: 6 = H5, 7 = H6, ..., 13 = H12
pub(super) fn convert_atom_hydrogen_count_code(
    code: u8,
    extended_range: bool,
) -> Result<Option<u8>, ParseError> {
    match code {
        0 => Ok(None),
        1..=5 => Ok(Some(code - 1)),
        6..=13 if extended_range => Ok(Some(code - 1)),
        _ => Err(ParseError::InvalidValenceCode(code)),
    }
}

/// Convert atom stereo care box code.
/// 'bbb' field: 0 = ignore stereo, 1 = stereo must match
/// Note: AtomStereoCare has no default value - 0 always means None
pub(super) fn convert_atom_stereo_care_code(
    code: u8,
) -> Result<Option<AtomStereoCare>, ParseError> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomStereoCare::Care)),
        _ => Err(ParseError::InvalidCode {
            field: "stereo care",
            value: code as i32,
        }),
    }
}

/// Convert atom valence code (default, explicit, explicit zero).
/// 'vvv' field: 0 = default, 1..=14 = explicit, 15 = explicit 0
/// Returns error for invalid valence codes.
pub(super) fn convert_atom_valence_code(code: u8) -> Result<Option<u8>, ParseError> {
    match code {
        0 => Ok(None),             // default/unspecified valence
        v @ 1..=14 => Ok(Some(v)), // explicit valences
        15 => Ok(Some(0)),         // explicit zero valence
        _ => Err(ParseError::InvalidValenceCode(code)),
    }
}

/// Convert atom inversion flag code.
/// 'nnn' field: 0 = not applicable, 1 = inverted, 2 = retained
/// Note: AtomInversionRetention has no default value - 0 always means None
pub(super) fn convert_atom_inversion_flag_code(
    code: u8,
) -> Result<Option<AtomInversionRetention>, ParseError> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomInversionRetention::Inverted)),
        2 => Ok(Some(AtomInversionRetention::Retained)),
        _ => Err(ParseError::InvalidCode {
            field: "inversion flag",
            value: code as i32,
        }),
    }
}

/// Convert atom exact change flag code.
/// 'eee' field: 0 = change allowed, 1 = exact change required
/// Note: AtomExactChange has no default value - 0 always means None
pub(super) fn convert_atom_exact_change_flag_code(
    code: u8,
) -> Result<Option<AtomExactChange>, ParseError> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomExactChange::Match)),
        _ => Err(ParseError::InvalidCode {
            field: "exact change flag",
            value: code as i32,
        }),
    }
}

/// Convert bond type code (basic molecules only)
/// 'ttt' field - bond type (1=Single, 2=Double, 3=Triple, 4=Aromatic)
pub(super) fn convert_bond_type_code(
    code: u8,
    extended_range: bool,
) -> Result<BondOrder, ParseError> {
    match code {
        1 => Ok(BondOrder::Single),
        2 => Ok(BondOrder::Double),
        3 => Ok(BondOrder::Triple),
        4 => Ok(BondOrder::Aromatic),
        _ => {
            if extended_range && (code == 0 || (9..=11).contains(&code)) {
                match code {
                    0 => Ok(BondOrder::Zero),
                    9 => Ok(BondOrder::Quadruple),
                    10 => Ok(BondOrder::Quintuple),
                    11 => Ok(BondOrder::Sextuple),
                    _ => unreachable!(),
                }
            } else {
                Err(ParseError::InvalidCode {
                    field: "bond type",
                    value: code as i32,
                })
            }
        }
    }
}

/// Convert bond type code
/// 'ttt' field - bond type (1=Single, 2=Double, 3=Triple, 4=Aromatic,
/// allow_queries: 5=SingleOrDouble, 6=SingleOrAromatic, 7=DoubleOrAromatic, 8=Any,
/// extended_range: 9=Quadruple, 10=Quintuple, 11=Sextuple, 0=Zero)
pub(super) fn convert_extended_bond_type_code(
    code: u8,
    extended_range: bool,
    allow_wildcards: bool,
) -> Result<BondOrder, ParseError> {
    match code {
        1 => Ok(BondOrder::Single),
        2 => Ok(BondOrder::Double),
        3 => Ok(BondOrder::Triple),
        4 => Ok(BondOrder::Aromatic),
        _ => {
            if allow_wildcards && (5..=8).contains(&code) {
                match code {
                    5 => Ok(BondOrder::SingleOrDouble),
                    6 => Ok(BondOrder::SingleOrAromatic),
                    7 => Ok(BondOrder::DoubleOrAromatic),
                    8 => Ok(BondOrder::Any),
                    _ => unreachable!(),
                }
            } else if extended_range && (code == 0 || (9..=11).contains(&code)) {
                match code {
                    0 => Ok(BondOrder::Zero),
                    9 => Ok(BondOrder::Quadruple),
                    10 => Ok(BondOrder::Quintuple),
                    11 => Ok(BondOrder::Sextuple),
                    _ => unreachable!(),
                }
            } else {
                Err(ParseError::InvalidCode {
                    field: "bond type",
                    value: code as i32,
                })
            }
        }
    }
}

/// Convert bond stereo/direction code
/// 'sss' field - can mean stereo for double bonds or direction for single bonds
/// Stereo (double bond): (0=Unknown, 1=Cis, 3|4=Either, 6=Trans)
/// Direction (single bond): (0=Not stereo, 1=Up, 3|4=Either, 6=Down)
/// NOTE: CTFile docs do not define cis/trans, 3=Double Either, 4=Single Either
pub(super) fn convert_bond_stereo_direction_code(
    code: u8,
) -> Result<(Option<BondStereo>, Option<BondDirection>), ParseError> {
    match code {
        0 => Ok((None, None)),
        1 => Ok((Some(BondStereo::Cis), Some(BondDirection::Up))),
        3 | 4 => Ok((Some(BondStereo::Either), Some(BondDirection::Either))),
        6 => Ok((Some(BondStereo::Trans), Some(BondDirection::Down))),
        _ => Err(ParseError::InvalidCode {
            field: "bond stereo/direction",
            value: code as i32,
        }),
    }
}

/// Convert bond topology code
/// 'rrr' field - bond topology (0=Either, 1=Ring, 2=Chain)
pub(super) fn convert_bond_topology_code(code: u8) -> Result<Option<BondTopology>, ParseError> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(BondTopology::Ring)),
        2 => Ok(Some(BondTopology::Chain)),
        _ => Err(ParseError::InvalidCode {
            field: "bond topology",
            value: code as i32,
        }),
    }
}

/// Convert bond reacting center code
/// 'ccc' field - bond reacting center (0=Not reacting, 1=Reacting, -1=Not a center,
/// 2=No change, 4=Bond made/broken, 8=Bond order changes)
/// If extended range, allow 3, 6, 7, 10, 11, 14, 15 (meaning no change)
pub(super) fn convert_bond_reacting_center_code(
    code: i8,
    extended_range: bool,
) -> Result<Option<BondReactingCenter>, ParseError> {
    if code == 0 {
        return Ok(None);
    }
    if code == -1 {
        return Ok(Some(BondReactingCenter::NOT_CENTER));
    }
    if !(-1..=15).contains(&code) {
        return Err(ParseError::InvalidCode {
            field: "reacting center",
            value: code as i32,
        });
    }

    // Positive codes can be partially combined:
    //   1 = a center
    //   2 = no change (cannot be combined with other flags unless extended range is set)
    //   4 = bond made/broken
    //   8 = bond order changes
    // Allowed combinations:
    // 12 = 4 + 8 (both made/broken and changes)
    // 5 = 4 + 1 (made/broken and center)
    // 9 = 8 + 1 (changes and center)
    // 13 = 12 + 1 (both made/broken and changes and center)
    let mut flags = BondReactingCenter::empty();

    if code & 2 != 0 {
        // 2 = no change
        if code != 2 && !extended_range {
            return Err(ParseError::InvalidCode {
                field: "reacting center",
                value: code as i32,
            });
        }
        flags |= BondReactingCenter::NO_CHANGE;
    } else {
        if code & 1 != 0 {
            // 1 = a center
            flags |= BondReactingCenter::CENTER;
        }
        if code & 4 != 0 {
            // 4 = bond made/broken
            flags |= BondReactingCenter::MADE_BROKEN;
        }
        if code & 8 != 0 {
            // 8 = bond order changes
            flags |= BondReactingCenter::ORDER_CHANGED;
        }
    }

    Ok(Some(flags))
}

/// Convert radical type code (property block)
/// radical type (0=no radical, 1=singlet (:), 2=doublet (. or ^), 3=triplet (^^))
/// Returns unpaired electron count: 0 for singlet, 1 for doublet, 2 for triplet
pub(super) fn convert_radical_type_code(code: u8) -> Result<Option<u8>, ParseError> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(0)), // Singlet: 0 unpaired electrons
        2 => Ok(Some(1)), // Doublet: 1 unpaired electron
        3 => Ok(Some(2)), // Triplet: 2 unpaired electrons
        _ => Err(ParseError::InvalidCode {
            field: "radical type",
            value: code as i32,
        }),
    }
}

// Convert ring bond count code (property block)
// 'vvv' field: ring bond count (-2 = as drawn (r*), -1 = no ring bonds (r0), 0 = off, 2 = r2, 3 = r3, 4 = r4+)
pub(super) fn convert_ring_bond_count_code(code: i8) -> Result<Option<RingBondCount>, ParseError> {
    match code {
        -2 => Ok(Some(RingBondCount::AsDrawn)),
        -1 => Ok(Some(RingBondCount::NoRingBonds)),
        0 => Ok(None),
        2 => Ok(Some(RingBondCount::R2)),
        3 => Ok(Some(RingBondCount::R3)),
        4 => Ok(Some(RingBondCount::R4Plus)),
        _ => Err(ParseError::InvalidCode {
            field: "ring bond count",
            value: code as i32,
        }),
    }
}

// Convert substitution count code (property block)
// 'vvv' field: substitution count (-2 = as drawn (s*), -1 = no substitution (s0), 0 = off, 1-5 = s1-s5,
// 6 = s6+),
// Extended range: 6-10 = s6-s10
pub(super) fn convert_substitution_count_code(
    code: i8,
    extended_range: bool,
) -> Result<Option<SubstitutionCount>, ParseError> {
    match code {
        -2 => Ok(Some(SubstitutionCount::AsDrawn)),
        -1 => Ok(Some(SubstitutionCount::NoSubstitution)),
        0 => Ok(None),
        1 => Ok(Some(SubstitutionCount::S1)),
        2 => Ok(Some(SubstitutionCount::S2)),
        3 => Ok(Some(SubstitutionCount::S3)),
        4 => Ok(Some(SubstitutionCount::S4)),
        5 => Ok(Some(SubstitutionCount::S5)),
        _ => {
            if extended_range {
                match code {
                    6 => Ok(Some(SubstitutionCount::S6)),
                    7 => Ok(Some(SubstitutionCount::S7)),
                    8 => Ok(Some(SubstitutionCount::S8)),
                    9 => Ok(Some(SubstitutionCount::S9)),
                    10 => Ok(Some(SubstitutionCount::S10)),
                    _ => Err(ParseError::InvalidCode {
                        field: "substitution count",
                        value: code as i32,
                    }),
                }
            } else if code == 6 {
                Ok(Some(SubstitutionCount::S6Plus))
            } else {
                Err(ParseError::InvalidCode {
                    field: "substitution count",
                    value: code as i32,
                })
            }
        }
    }
}

/// Convert unsaturated atom code (property block)
/// 'vvv' field: unsaturated (0=off, 1=on)
pub(super) fn convert_unsaturated_atom_code(
    code: u8,
) -> Result<Option<UnsaturatedAtom>, ParseError> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(UnsaturatedAtom)),
        _ => Err(ParseError::InvalidCode {
            field: "unsaturated atom",
            value: code as i32,
        }),
    }
}

/// Convert attachment point code (property block)
/// 'vvv' field - attachment point (0=none, 1=first, 2=second, 3=both)
pub(super) fn convert_attachment_point_code(
    code: u8,
) -> Result<Option<AttachmentPointType>, ParseError> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AttachmentPointType::First)),
        2 => Ok(Some(AttachmentPointType::Second)),
        3 => Ok(Some(AttachmentPointType::Both)),
        _ => Err(ParseError::InvalidCode {
            field: "attachment point",
            value: code as i32,
        }),
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_data::NamedIsotope;

    use super::*;

    #[rstest]
    #[case::zero(0, None)]
    #[case::negative(-3, Some(-3))]
    #[case::positive(4, Some(4))]
    #[case::out_of_range_high(5, None)]
    #[case(5, None)]
    fn test_convert_atom_mass_diff_code(#[case] code: i8, #[case] expected: Option<i8>) {
        assert_eq!(convert_atom_mass_diff_code(code), expected);
    }

    #[rstest]
    #[case::element(AtomSymbol::Element(Element::C), None, (Element::C, None))]
    #[case::element_mass_diff(AtomSymbol::Element(Element::C), Some(1), (Element::C, Some(13)))]
    #[case::named_isotope(AtomSymbol::NamedIsotope(NamedIsotope::D), None, (Element::H, Some(2)))]
    #[case::named_isotope_mass_diff(AtomSymbol::NamedIsotope(NamedIsotope::D), Some(3), (Element::H, Some(2)))]
    fn test_convert_atom_symbol_mass_diff(
        #[case] symbol: AtomSymbol,
        #[case] mass_diff: Option<i8>,
        #[case] expected: (Element, Option<u32>),
    ) {
        assert_eq!(convert_atom_symbol_mass_diff(symbol, mass_diff), expected);
    }

    #[rstest]
    #[case::zero(0, None, None)]
    #[case::one(1, Some(3), None)]
    #[case::doublet_radical(4, None, Some(1))] // Doublet radical: 1 unpaired electron
    #[case::minus_one(5, Some(-1), None)]
    #[case::out_of_range_high(8, None, None)]
    fn test_convert_atom_charge_code(
        #[case] code: u8,
        #[case] expected_charge: Option<i8>,
        #[case] expected_unpaired_e: Option<u8>,
    ) {
        assert_eq!(
            convert_atom_charge_code(code),
            (expected_charge, expected_unpaired_e)
        );
    }

    #[rstest]
    #[case::no_isotope(Element::C, 0, None)]
    #[case::default_isotope(Element::C, 12, Some(12))]
    #[case::isotope(Element::C, 13, Some(13))]
    fn test_convert_atom_isotope_mass_number(
        #[case] element: Element,
        #[case] mass_number: u32,
        #[case] expected: Option<u32>,
    ) {
        assert_eq!(
            convert_atom_isotope_mass_number(element, mass_number, false).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::not_catalogued(Element::C, 40)]
    fn test_convert_atom_isotope_mass_number_invalid(
        #[case] element: Element,
        #[case] mass_number: u32,
    ) {
        assert!(
            convert_atom_isotope_mass_number(element, mass_number, false).is_err(),
            "{}{} is not catalogued",
            mass_number,
            element
        );
    }

    #[rstest]
    #[case::no_isotope(Element::C, 0, None)]
    #[case::default_isotope(Element::C, 12, Some(12))]
    #[case::isotope(Element::C, 40, Some(40))]
    fn test_convert_atom_isotope_mass_number_extended(
        #[case] element: Element,
        #[case] mass_number: u32,
        #[case] expected: Option<u32>,
    ) {
        assert_eq!(
            convert_atom_isotope_mass_number(element, mass_number, true).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::none(0, None)]
    #[case::clockwise(1, Some(Chirality::Clockwise))]
    #[case::counter_clockwise(2, Some(Chirality::CounterClockwise))]
    #[case::unspecified(3, Some(Chirality::Unspecified))]
    fn test_convert_atom_stereo_parity_code(#[case] code: u8, #[case] expected: Option<Chirality>) {
        assert_eq!(convert_atom_stereo_parity_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_high(4)]
    fn test_convert_atom_stereo_parity_code_invalid(#[case] code: u8) {
        assert!(convert_atom_stereo_parity_code(code).is_err());
    }

    #[rstest]
    #[case::none(0, None)]
    #[case::zero(1, Some(0))]
    #[case::one(2, Some(1))]
    #[case::two(3, Some(2))]
    #[case::three(4, Some(3))]
    #[case::four(5, Some(4))]
    fn test_convert_atom_hydrogen_count_code(#[case] code: u8, #[case] expected: Option<u8>) {
        assert_eq!(
            convert_atom_hydrogen_count_code(code, false).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::extended_range_high(6)]
    #[case::out_of_range_high(14)]
    fn test_convert_atom_hydrogen_count_code_invalid(#[case] code: u8) {
        assert!(
            convert_atom_hydrogen_count_code(code, false).is_err(),
            "{} should have failed",
            code
        );
    }

    #[rstest]
    #[case::five(6, Some(5))]
    #[case::six(7, Some(6))]
    #[case::seven(8, Some(7))]
    #[case::eight(9, Some(8))]
    #[case::nine(10, Some(9))]
    #[case::twelve(13, Some(12))]
    fn test_convert_atom_hydrogen_count_code_extended(
        #[case] code: u8,
        #[case] expected: Option<u8>,
    ) {
        assert_eq!(
            convert_atom_hydrogen_count_code(code, true).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::none(0, None)]
    #[case::care(1, Some(AtomStereoCare::Care))]
    fn test_convert_atom_stereo_care_code(
        #[case] code: u8,
        #[case] expected: Option<AtomStereoCare>,
    ) {
        assert_eq!(convert_atom_stereo_care_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_high(2)]
    fn test_convert_atom_stereo_care_code_invalid(#[case] code: u8) {
        assert!(
            convert_atom_stereo_care_code(code).is_err(),
            "{} should have failed",
            code
        );
    }

    #[rstest]
    #[case::none(0, None)]
    #[case::one(1, Some(1))]
    #[case::fourteen(14, Some(14))]
    #[case::zero(15, Some(0))]
    fn test_convert_atom_valence_code(#[case] code: u8, #[case] expected: Option<u8>) {
        assert_eq!(convert_atom_valence_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_high(16)]
    fn test_convert_atom_valence_code_invalid(#[case] code: u8) {
        assert!(
            convert_atom_valence_code(code).is_err(),
            "{} should have failed",
            code
        );
    }

    #[rstest]
    #[case::none(0, None)]
    #[case::inverted(1, Some(AtomInversionRetention::Inverted))]
    #[case::retained(2, Some(AtomInversionRetention::Retained))]
    fn test_convert_atom_inversion_flag_code(
        #[case] code: u8,
        #[case] expected: Option<AtomInversionRetention>,
    ) {
        assert_eq!(convert_atom_inversion_flag_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_high(3)]
    fn test_convert_atom_inversion_flag_code_invalid(#[case] code: u8) {
        assert!(
            convert_atom_inversion_flag_code(code).is_err(),
            "{} should have failed",
            code
        );
    }

    #[rstest]
    #[case::none(0, None)]
    #[case::match_flag(1, Some(AtomExactChange::Match))]
    fn test_convert_atom_exact_change_flag_code(
        #[case] code: u8,
        #[case] expected: Option<AtomExactChange>,
    ) {
        assert_eq!(convert_atom_exact_change_flag_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_high(2)]
    fn test_convert_atom_exact_change_flag_code_invalid(#[case] code: u8) {
        assert!(
            convert_atom_exact_change_flag_code(code).is_err(),
            "{} should have failed",
            code
        );
    }

    #[rstest]
    #[case::single(1, BondOrder::Single)]
    #[case::double(2, BondOrder::Double)]
    #[case::triple(3, BondOrder::Triple)]
    #[case::aromatic(4, BondOrder::Aromatic)]
    fn test_convert_bond_type_code(#[case] code: u8, #[case] expected: BondOrder) {
        assert_eq!(convert_bond_type_code(code, false).unwrap(), expected);
    }

    #[rstest]
    #[case::zero(0)]
    #[case::query(5)]
    #[case::extended_range_high(9)]
    fn test_convert_bond_type_code_invalid(#[case] code: u8) {
        assert!(
            convert_bond_type_code(code, false).is_err(),
            "{} should have failed",
            code
        );
    }

    #[rstest]
    #[case::zero(0, BondOrder::Zero)]
    #[case::quadruple(9, BondOrder::Quadruple)]
    fn test_convert_bond_type_code_extended(#[case] code: u8, #[case] expected: BondOrder) {
        assert_eq!(convert_bond_type_code(code, true).unwrap(), expected);
    }

    #[rstest]
    #[case::single(1, BondOrder::Single)]
    #[case::double(2, BondOrder::Double)]
    #[case::triple(3, BondOrder::Triple)]
    #[case::aromatic(4, BondOrder::Aromatic)]
    fn test_convert_extended_bond_type_code(#[case] code: u8, #[case] expected: BondOrder) {
        assert_eq!(
            convert_extended_bond_type_code(code, false, false).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::zero(0)]
    #[case::quadruple(9)]
    #[case::out_of_range_high(12)]
    fn test_convert_extended_bond_type_code_invalid(#[case] code: u8) {
        assert!(
            convert_extended_bond_type_code(code, false, false).is_err(),
            "{} should have failed",
            code
        );
    }

    #[rstest]
    #[case::zero(0, BondOrder::Zero)]
    #[case::single(1, BondOrder::Single)]
    #[case::quadruple(9, BondOrder::Quadruple)]
    #[case::quintuple(10, BondOrder::Quintuple)]
    #[case::sextuple(11, BondOrder::Sextuple)]
    fn test_convert_extended_bond_type_code_extended(
        #[case] code: u8,
        #[case] expected: BondOrder,
    ) {
        assert_eq!(
            convert_extended_bond_type_code(code, true, false).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::single_or_double(5, BondOrder::SingleOrDouble)]
    #[case::single_or_aromatic(6, BondOrder::SingleOrAromatic)]
    #[case::double_or_aromatic(7, BondOrder::DoubleOrAromatic)]
    #[case::any(8, BondOrder::Any)]
    fn test_convert_extended_bond_type_code_query(#[case] code: u8, #[case] expected: BondOrder) {
        assert_eq!(
            convert_extended_bond_type_code(code, false, true).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::not_stereo(0, (None, None))]
    #[case::cis(1, (Some(BondStereo::Cis), Some(BondDirection::Up)))]
    #[case::either(3, (Some(BondStereo::Either), Some(BondDirection::Either)))]
    #[case::unknown(4, (Some(BondStereo::Either), Some(BondDirection::Either)))]
    #[case::trans(6, (Some(BondStereo::Trans), Some(BondDirection::Down)))]
    fn test_convert_bond_stereo_direction_code(
        #[case] code: u8,
        #[case] expected: (Option<BondStereo>, Option<BondDirection>),
    ) {
        assert_eq!(convert_bond_stereo_direction_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::unused_2(2)]
    #[case::unused_5(5)]
    #[case::out_of_range_high(7)]
    fn test_convert_bond_stereo_direction_code_invalid(#[case] code: u8) {
        assert!(convert_bond_stereo_direction_code(code).is_err());
    }

    #[rstest]
    #[case::none(0, None)]
    #[case::ring(1, Some(BondTopology::Ring))]
    #[case::chain(2, Some(BondTopology::Chain))]
    fn test_convert_bond_topology_code(#[case] code: u8, #[case] expected: Option<BondTopology>) {
        assert_eq!(convert_bond_topology_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_high(3)]
    fn test_convert_bond_topology_code_invalid(#[case] code: u8) {
        assert!(convert_bond_topology_code(code).is_err());
    }

    #[rstest]
    #[case::unmarked(0, None)]
    #[case::not_center(-1, Some(BondReactingCenter::NOT_CENTER))]
    #[case::center(1, Some(BondReactingCenter::CENTER))]
    #[case::no_change(2, Some(BondReactingCenter::NO_CHANGE))]
    #[case::made_broken(4, Some(BondReactingCenter::MADE_BROKEN))]
    #[case::order_changed(8, Some(BondReactingCenter::ORDER_CHANGED))]
    #[case::center_and_made_broken(5, Some(BondReactingCenter::CENTER | BondReactingCenter::MADE_BROKEN))]
    #[case::center_and_order_changed(9, Some(BondReactingCenter::CENTER | BondReactingCenter::ORDER_CHANGED))]
    #[case::made_broken_and_order_changed(12, Some(BondReactingCenter::MADE_BROKEN | BondReactingCenter::ORDER_CHANGED))]
    #[case::center_and_made_broken_and_order_changed(13, Some(
        BondReactingCenter::CENTER
            | BondReactingCenter::MADE_BROKEN
            | BondReactingCenter::ORDER_CHANGED
    ))]
    fn test_convert_bond_reacting_center_code(
        #[case] code: i8,
        #[case] expected: Option<BondReactingCenter>,
    ) {
        assert_eq!(
            convert_bond_reacting_center_code(code, false).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::invalid_combination(3)]
    #[case::out_of_range_high(16)]
    #[case::out_of_range_low(-2)]
    fn test_convert_bond_reacting_center_code_invalid(#[case] code: i8) {
        assert!(convert_bond_reacting_center_code(code, false).is_err());
    }

    // Extended range allows codes 3,6,7,10,11,14,15 which all map to just NO_CHANGE
    #[rstest]
    #[case::no_change_3(3, Some(BondReactingCenter::NO_CHANGE))]
    #[case::no_change_6(6, Some(BondReactingCenter::NO_CHANGE))]
    #[case::no_change_7(7, Some(BondReactingCenter::NO_CHANGE))]
    #[case::no_change_10(10, Some(BondReactingCenter::NO_CHANGE))]
    #[case::no_change_11(11, Some(BondReactingCenter::NO_CHANGE))]
    #[case::no_change(14, Some(BondReactingCenter::NO_CHANGE))]
    #[case::no_change_15(15, Some(BondReactingCenter::NO_CHANGE))]
    fn test_convert_bond_reacting_center_code_extended(
        #[case] code: i8,
        #[case] expected: Option<BondReactingCenter>,
    ) {
        assert_eq!(
            convert_bond_reacting_center_code(code, true).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::none(0, None)]
    #[case::singlet(1, Some(0))] // Singlet: 0 unpaired electrons
    #[case::doublet(2, Some(1))] // Doublet: 1 unpaired electron
    #[case::triplet(3, Some(2))] // Triplet: 2 unpaired electrons
    fn test_convert_radical_type_code(#[case] code: u8, #[case] expected: Option<u8>) {
        assert_eq!(convert_radical_type_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_high(4)]
    fn test_convert_radical_type_code_invalid(#[case] code: u8) {
        assert!(convert_radical_type_code(code).is_err());
    }

    #[rstest]
    #[case::as_drawn(-2, Some(RingBondCount::AsDrawn))]
    #[case::no_ring_bonds(-1, Some(RingBondCount::NoRingBonds))]
    #[case::r2(2, Some(RingBondCount::R2))]
    #[case::r3(3, Some(RingBondCount::R3))]
    #[case::r4_plus(4, Some(RingBondCount::R4Plus))]
    fn test_convert_ring_bond_count_code(
        #[case] code: i8,
        #[case] expected: Option<RingBondCount>,
    ) {
        assert_eq!(convert_ring_bond_count_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_high(5)]
    fn test_convert_ring_bond_count_code_invalid(#[case] code: i8) {
        assert!(convert_ring_bond_count_code(code).is_err());
    }

    #[rstest]
    #[case::as_drawn(-2, Some(SubstitutionCount::AsDrawn))]
    #[case::no_substitution(-1, Some(SubstitutionCount::NoSubstitution))]
    #[case::s1(1, Some(SubstitutionCount::S1))]
    #[case::s6_plus(6, Some(SubstitutionCount::S6Plus))]
    fn test_convert_substitution_count_code(
        #[case] code: i8,
        #[case] expected: Option<SubstitutionCount>,
    ) {
        assert_eq!(
            convert_substitution_count_code(code, false).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::out_of_range_high(7)]
    fn test_convert_substitution_count_code_invalid(#[case] code: i8) {
        assert!(convert_substitution_count_code(code, false).is_err());
    }

    #[rstest]
    #[case::s6(6, Some(SubstitutionCount::S6))]
    fn test_convert_substitution_count_code_extended(
        #[case] code: i8,
        #[case] expected: Option<SubstitutionCount>,
    ) {
        assert_eq!(
            convert_substitution_count_code(code, true).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::none(0, None)]
    #[case::unsaturated(1, Some(UnsaturatedAtom))]
    fn test_convert_unsaturated_atom_code(
        #[case] code: u8,
        #[case] expected: Option<UnsaturatedAtom>,
    ) {
        assert_eq!(convert_unsaturated_atom_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_high(2)]
    fn test_convert_unsaturated_atom_code_invalid(#[case] code: u8) {
        assert!(convert_unsaturated_atom_code(code).is_err());
    }

    #[rstest]
    #[case::none(0, None)]
    #[case::first(1, Some(AttachmentPointType::First))]
    #[case::second(2, Some(AttachmentPointType::Second))]
    #[case::both(3, Some(AttachmentPointType::Both))]
    fn test_convert_attachment_point_code(
        #[case] code: u8,
        #[case] expected: Option<AttachmentPointType>,
    ) {
        assert_eq!(convert_attachment_point_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case::out_of_range_high(4)]
    fn test_convert_attachment_point_code_invalid(#[case] code: u8) {
        assert!(convert_attachment_point_code(code).is_err());
    }
}
