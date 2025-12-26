//! Convert numerical codes used in MOL files to enums
//!
//! All functions return simple_ir types.

use umol_data::{Element, Isotope};

use crate::io::ctfile::error::SemanticError;
use crate::simple_ir::{
    AtomExactChange, AtomInversionRetention, AtomStereoCare, AtomStereoParity, AtomSymbol,
    AttachmentPointType, BondDir, BondOrder, BondReactingCenter, BondStereo, BondTopology,
    AtomRadical, RingBondCount, SubstitutionCount, UnsaturatedAtom,
};

type Result<T> = std::result::Result<T, SemanticError>;

/// Convert atom mass difference code (atom block)
/// 'dd' field: mass difference (-3..=4), None if 0 or value outside of this range
pub(crate) fn convert_atom_mass_diff_code(code: i8) -> Option<i8> {
    match code {
        -3..=-1 | 1..=4 => Some(code),
        _ => None,
    }
}

/// Convert atom symbol and mass difference to element and isotope mass
/// 'ss' field: atom symbol, 'dd' field: mass difference
/// Processes elements and named isotopes.
/// Returns error for extended atom symbols (L, A, Q, *, LP, R#)
pub(crate) fn convert_atom_symbol_mass_diff(
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
pub(crate) fn convert_atom_charge_code(code: u8) -> (i8, Option<AtomRadical>) {
    match code {
        1..=3 | 5..=7 => (4 - code as i8, None),
        4 => (0, Some(AtomRadical::Doublet)),
        _ => (0, None),
    }
}

/// Validate atom isotope mass number (property block)
/// 'mmm' field: isotope mass number
pub(crate) fn convert_atom_isotope_mass_number(
    element: Element,
    mass_number: u32,
    extended_isotopes: bool,
) -> Result<Option<u32>> {
    if mass_number == 0 {
        return Ok(None);
    }
    if extended_isotopes || Isotope::is_catalogued(element, mass_number) {
        Ok(Some(mass_number))
    } else {
        Err(SemanticError::Generic(format!(
            "Invalid isotope mass number {} for element {}",
            mass_number, element
        )))
    }
}

/// Convert atom stereo parity code (not stereo, odd, even, either or unmarked).
// 'sss' field: 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked
/// use_defaults: if true, include default values (code 0), if false, return None for defaults
pub(crate) fn convert_atom_stereo_parity_code(
    code: u8,
    use_defaults: bool,
) -> Result<Option<AtomStereoParity>> {
    match code {
        0 if use_defaults => Ok(Some(AtomStereoParity::Either)),
        0 => Ok(None),
        1 => Ok(Some(AtomStereoParity::Odd)),
        2 => Ok(Some(AtomStereoParity::Even)),
        3 => Ok(Some(AtomStereoParity::Either)),
        _ => Err(SemanticError::InvalidStereoParity(code)),
    }
}

/// Convert atom hydrogen count code (extension: 0 in non-query atoms).
/// 'hhh' field: 0 = non-query atom, 1 = H0, 2 = H1, 3 = H2, 4 = H3, 5 = H4
/// Extended range: 6 = H5, 7 = H6, ..., 13 = H12
pub(crate) fn convert_atom_hydrogen_count_code(
    code: u8,
    extended_range: bool,
) -> Result<Option<u8>> {
    match code {
        0 => Ok(None),
        1..=5 => Ok(Some(code - 1)),
        6..=13 if extended_range => Ok(Some(code - 1)),
        _ => Err(SemanticError::InvalidValenceCode(code)),
    }
}

/// Convert atom stereo care box code.
/// 'bbb' field: 0 = ignore stereo, 1 = stereo must match
/// Note: AtomStereoCare has no default value - 0 always means None
pub(crate) fn convert_atom_stereo_care_code(code: u8) -> Result<Option<AtomStereoCare>> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomStereoCare::Care)),
        _ => Err(SemanticError::Generic(format!("Invalid stereo care code '{}'", code))),
    }
}

/// Convert atom valence code (default, explicit, explicit zero).
/// 'vvv' field: 0 = default, 1..=14 = explicit, 15 = explicit 0
/// Returns error for invalid valence codes.
pub(crate) fn convert_atom_valence_code(code: u8) -> Result<Option<u8>> {
    match code {
        0 => Ok(None),             // default/unspecified valence
        v @ 1..=14 => Ok(Some(v)), // explicit valences
        15 => Ok(Some(0)),         // explicit zero valence
        _ => Err(SemanticError::InvalidValenceCode(code)),
    }
}

/// Convert atom inversion flag code.
/// 'nnn' field: 0 = not applicable, 1 = inverted, 2 = retained
/// Note: AtomInversionRetention has no default value - 0 always means None
pub(crate) fn convert_atom_inversion_flag_code(code: u8) -> Result<Option<AtomInversionRetention>> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomInversionRetention::Inverted)),
        2 => Ok(Some(AtomInversionRetention::Retained)),
        _ => Err(SemanticError::Generic(format!("Invalid inversion flag code '{}'", code))),
    }
}

/// Convert atom exact change flag code.
/// 'eee' field: 0 = change allowed, 1 = exact change required
/// Note: AtomExactChange has no default value - 0 always means None
pub(crate) fn convert_atom_exact_change_flag_code(code: u8) -> Result<Option<AtomExactChange>> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomExactChange::Match)),
        _ => Err(SemanticError::Generic(format!("Invalid exact change flag code '{}'", code))),
    }
}

/// Convert bond type code (basic molecules only)
/// 'ttt' field - bond type (1=Single, 2=Double, 3=Triple, 4=Aromatic)
pub(crate) fn convert_bond_type_code(code: u8) -> Result<BondOrder> {
    match code {
        1 => Ok(BondOrder::Single),
        2 => Ok(BondOrder::Double),
        3 => Ok(BondOrder::Triple),
        4 => Ok(BondOrder::Aromatic),
        _ => Err(SemanticError::Generic(format!("Invalid bond type code '{}'", code))),
    }
}

/// Convert bond type code
/// 'ttt' field - bond type (1=Single, 2=Double, 3=Triple, 4=Aromatic,
/// allow_queries: 5=SingleOrDouble, 6=SingleOrAromatic, 7=DoubleOrAromatic, 8=Any,
/// extended_range: 9=Quadruple, 10=Quintuple, 11=Sextuple, 0=Zero)
pub(crate) fn convert_extended_bond_type_code(
    code: u8,
    extended_range: bool,
    allow_queries: bool,
) -> Result<BondOrder> {
    match code {
        1 => Ok(BondOrder::Single),
        2 => Ok(BondOrder::Double),
        3 => Ok(BondOrder::Triple),
        4 => Ok(BondOrder::Aromatic),
        _ => {
            if allow_queries && (5..=8).contains(&code) {
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
                Err(SemanticError::Generic(format!("Invalid bond type code '{}'", code)))
            }
        }
    }
}

/// Convert bond stereo/direction code
/// 'sss' field - can mean stereo for double bonds or direction for single bonds
/// Stereo: (0=Not stereo, 1=Up, 3=Either, 4=Unknown, 6=Down)
/// Direction: (1=Up, 6=Down)
/// use_defaults: if true, include default values (code 0), if false, return None for defaults
pub(crate) fn convert_bond_stereo_dir_code(
    code: u8,
    use_defaults: bool,
) -> Result<(Option<BondStereo>, Option<BondDir>)> {
    match code {
        0 if use_defaults => Ok((Some(BondStereo::default()), Some(BondDir::default()))),
        0 => Ok((None, None)),
        1 => Ok((Some(BondStereo::Cis), Some(BondDir::Up))),
        3 | 4 => Ok((Some(BondStereo::Either), Some(BondDir::Either))),
        6 => Ok((Some(BondStereo::Trans), Some(BondDir::Down))),
        _ => Err(SemanticError::Generic(format!(
            "Invalid bond stereo/direction code '{}'",
            code
        ))),
    }
}

/// Convert bond topology code
/// 'rrr' field - bond topology (0=Either, 1=Ring, 2=Chain)
/// use_defaults: if true, include default values (code 0), if false, return None for defaults
pub(crate) fn convert_bond_topology_code(
    code: u8,
    use_defaults: bool,
) -> Result<Option<BondTopology>> {
    match code {
        0 if use_defaults => Ok(Some(BondTopology::Either)),
        0 => Ok(None),
        1 => Ok(Some(BondTopology::Ring)),
        2 => Ok(Some(BondTopology::Chain)),
        _ => Err(SemanticError::Generic(format!("Invalid bond topology code '{}'", code))),
    }
}

/// Convert bond reacting center code
/// 'ccc' field - bond reacting center (0=Not reacting, 1=Reacting, -1=Not a center,
/// 2=No change, 4=Bond made/broken, 8=Bond order changes)
/// use_defaults: if true, include default values (code 0), if false, return None for defaults
/// If extended range, allow 3, 6, 7, 10, 11, 14, 15 (meaning no change)
pub(crate) fn convert_bond_reacting_center_code(
    code: i8,
    use_defaults: bool,
    extended_range: bool,
) -> Result<Option<BondReactingCenter>> {
    if code == 0 {
        if use_defaults {
            return Ok(Some(BondReactingCenter::UNMARKED));
        } else {
            return Ok(None);
        }
    }
    if code == -1 {
        return Ok(Some(BondReactingCenter::NOT_CENTER));
    }
    if !(-1..=15).contains(&code) {
        return Err(SemanticError::Generic(format!(
            "Invalid reacting center code '{}'",
            code
        )));
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
            return Err(SemanticError::Generic(format!(
                "Invalid reacting center code combination '{}'",
                code
            )));
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
pub(crate) fn convert_radical_type_code(code: u8) -> Result<Option<AtomRadical>> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomRadical::Singlet)),
        2 => Ok(Some(AtomRadical::Doublet)),
        3 => Ok(Some(AtomRadical::Triplet)),
        _ => Err(SemanticError::Generic(format!("Invalid radical type code '{}'", code))),
    }
}

// Convert ring bond count code (property block)
// 'vvv' field: ring bond count (-2 = as drawn (r*), -1 = no ring bonds (r0), 0 = off, 2 = r2, 3 = r3, 4 = r4+)
pub(crate) fn convert_ring_bond_count_code(code: i8) -> Result<Option<RingBondCount>> {
    match code {
        -2 => Ok(Some(RingBondCount::AsDrawn)),
        -1 => Ok(Some(RingBondCount::NoRingBonds)),
        0 => Ok(None),
        2 => Ok(Some(RingBondCount::R2)),
        3 => Ok(Some(RingBondCount::R3)),
        4 => Ok(Some(RingBondCount::R4Plus)),
        _ => Err(SemanticError::Generic(format!("Invalid ring bond count code '{}'", code))),
    }
}

// Convert substitution count code (property block)
// 'vvv' field: substitution count (-2 = as drawn (s*), -1 = no substitution (s0), 0 = off, 1-5 = s1-s5,
// 6 = s6+),
// Extended range: 6-10 = s6-s10
pub(crate) fn convert_substitution_count_code(
    code: i8,
    extended_range: bool,
) -> Result<Option<SubstitutionCount>> {
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
                    _ => Err(SemanticError::Generic(format!(
                        "Invalid substitution count code '{}'",
                        code
                    ))),
                }
            } else if code == 6 {
                Ok(Some(SubstitutionCount::S6Plus))
            } else {
                Err(SemanticError::Generic(format!(
                    "Invalid substitution count code '{}'",
                    code
                )))
            }
        }
    }
}

/// Convert unsaturated atom code (property block)
/// 'vvv' field: unsaturated (0=off, 1=on)
pub(crate) fn convert_unsaturated_atom_code(code: u8) -> Result<Option<UnsaturatedAtom>> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(UnsaturatedAtom)),
        _ => Err(SemanticError::Generic(format!(
            "Invalid unsaturated atom code '{}'",
            code
        ))),
    }
}

/// Convert attachment point code (property block)
/// 'vvv' field - attachment point (0=none, 1=first, 2=second, 3=both)
pub(crate) fn convert_attachment_point_code(code: u8) -> Result<Option<AttachmentPointType>> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AttachmentPointType::First)),
        2 => Ok(Some(AttachmentPointType::Second)),
        3 => Ok(Some(AttachmentPointType::Both)),
        _ => Err(SemanticError::Generic(format!(
            "Invalid attachment point code '{}'",
            code
        ))),
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_data::NamedIsotope;

    use super::*;

    #[rstest]
    #[case(0, None)]
    #[case(-3, Some(-3))]
    #[case(4, Some(4))]
    #[case(5, None)]
    fn test_convert_atom_mass_diff_code(#[case] code: i8, #[case] expected: Option<i8>) {
        assert_eq!(convert_atom_mass_diff_code(code), expected);
    }

    #[rstest]
    #[case(AtomSymbol::Element(Element::C), None, (Element::C, None))]
    #[case(AtomSymbol::Element(Element::C), Some(1), (Element::C, Some(13)))]
    #[case(AtomSymbol::NamedIsotope(NamedIsotope::D), None, (Element::H, Some(2)))]
    #[case(AtomSymbol::NamedIsotope(NamedIsotope::D), Some(3), (Element::H, Some(2)))]
    fn test_convert_atom_symbol_mass_diff(
        #[case] symbol: AtomSymbol,
        #[case] mass_diff: Option<i8>,
        #[case] expected: (Element, Option<u32>),
    ) {
        assert_eq!(convert_atom_symbol_mass_diff(symbol, mass_diff), expected);
    }

    #[rstest]
    #[case(0, 0, None)]
    #[case(1, 3, None)]
    #[case(4, 0, Some(AtomRadical::Doublet))]
    #[case(5, -1, None)]
    #[case(8, 0, None)]
    fn test_convert_atom_charge_code(
        #[case] code: u8,
        #[case] expected_charge: i8,
        #[case] expected_radical: Option<AtomRadical>,
    ) {
        assert_eq!(
            convert_atom_charge_code(code),
            (expected_charge, expected_radical)
        );
    }

    #[rstest]
    #[case(Element::C, 0, None)]
    #[case(Element::C, 12, Some(12))]
    #[case(Element::C, 13, Some(13))]
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
    #[case(Element::C, 0, None)]
    #[case(Element::C, 12, Some(12))]
    #[case(Element::C, 40, Some(40))]
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
    #[case(Element::C, 40, "40C is not catalogued")]
    fn test_convert_atom_isotope_mass_number_invalid(
        #[case] element: Element,
        #[case] mass_number: u32,
        #[case] desc: &str,
    ) {
        assert!(
            convert_atom_isotope_mass_number(element, mass_number, false).is_err(),
            "{} should have failed",
            desc
        );
    }

    #[rstest]
    #[case(0, Some(AtomStereoParity::Either), None)]
    #[case(1, Some(AtomStereoParity::Odd), Some(AtomStereoParity::Odd))]
    #[case(2, Some(AtomStereoParity::Even), Some(AtomStereoParity::Even))]
    #[case(3, Some(AtomStereoParity::Either), Some(AtomStereoParity::Either))]
    fn test_convert_atom_stereo_parity_code(
        #[case] code: u8,
        #[case] expected_default: Option<AtomStereoParity>,
        #[case] expected_basic: Option<AtomStereoParity>,
    ) {
        assert_eq!(
            convert_atom_stereo_parity_code(code, true).unwrap(),
            expected_default
        );
        assert_eq!(
            convert_atom_stereo_parity_code(code, false).unwrap(),
            expected_basic
        );
    }

    #[rstest]
    #[case(4, "too high")]
    fn test_convert_atom_stereo_parity_code_invalid(#[case] code: u8, #[case] desc: &str) {
        assert!(
            convert_atom_stereo_parity_code(code, true).is_err(),
            "{} should have failed",
            desc
        );
    }

    #[rstest]
    #[case(0, None)]
    #[case(1, Some(0))]
    #[case(2, Some(1))]
    #[case(3, Some(2))]
    #[case(4, Some(3))]
    #[case(5, Some(4))]
    fn test_convert_atom_hydrogen_count_code(#[case] code: u8, #[case] expected: Option<u8>) {
        assert_eq!(
            convert_atom_hydrogen_count_code(code, false).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case(6, "extended range")]
    #[case(14, "out of range")]
    fn test_convert_atom_hydrogen_count_code_invalid(#[case] code: u8, #[case] desc: &str) {
        assert!(
            convert_atom_hydrogen_count_code(code, false).is_err(),
            "{} should have failed",
            desc
        );
    }

    #[rstest]
    #[case(6, Some(5))]
    #[case(7, Some(6))]
    #[case(8, Some(7))]
    #[case(9, Some(8))]
    #[case(10, Some(9))]
    #[case(13, Some(12))]
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
    #[case(0, None)]
    #[case(1, Some(AtomStereoCare::Care))]
    fn test_convert_atom_stereo_care_code(
        #[case] code: u8,
        #[case] expected: Option<AtomStereoCare>,
    ) {
        assert_eq!(convert_atom_stereo_care_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(2, "invalid")]
    fn test_convert_atom_stereo_care_code_invalid(#[case] code: u8, #[case] desc: &str) {
        assert!(
            convert_atom_stereo_care_code(code).is_err(),
            "{} should have failed",
            desc
        );
    }

    #[rstest]
    #[case(0, None)]
    #[case(1, Some(1))]
    #[case(14, Some(14))]
    #[case(15, Some(0))]
    fn test_convert_atom_valence_code(#[case] code: u8, #[case] expected: Option<u8>) {
        assert_eq!(convert_atom_valence_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(16, "too high")]
    fn test_convert_atom_valence_code_invalid(#[case] code: u8, #[case] desc: &str) {
        assert!(
            convert_atom_valence_code(code).is_err(),
            "{} should have failed",
            desc
        );
    }

    #[rstest]
    #[case(0, None)]
    #[case(1, Some(AtomInversionRetention::Inverted))]
    #[case(2, Some(AtomInversionRetention::Retained))]
    fn test_convert_atom_inversion_flag_code(
        #[case] code: u8,
        #[case] expected: Option<AtomInversionRetention>,
    ) {
        assert_eq!(convert_atom_inversion_flag_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(3, "invalid")]
    fn test_convert_atom_inversion_flag_code_invalid(#[case] code: u8, #[case] desc: &str) {
        assert!(
            convert_atom_inversion_flag_code(code).is_err(),
            "{} should have failed",
            desc
        );
    }

    #[rstest]
    #[case(0, None)]
    #[case(1, Some(AtomExactChange::Match))]
    fn test_convert_atom_exact_change_flag_code(
        #[case] code: u8,
        #[case] expected: Option<AtomExactChange>,
    ) {
        assert_eq!(convert_atom_exact_change_flag_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(2, "invalid")]
    fn test_convert_atom_exact_change_flag_code_invalid(#[case] code: u8, #[case] desc: &str) {
        assert!(
            convert_atom_exact_change_flag_code(code).is_err(),
            "{} should have failed",
            desc
        );
    }

    #[rstest]
    #[case(1, BondOrder::Single)]
    #[case(2, BondOrder::Double)]
    #[case(3, BondOrder::Triple)]
    #[case(4, BondOrder::Aromatic)]
    fn test_convert_bond_type_code(#[case] code: u8, #[case] expected: BondOrder) {
        assert_eq!(convert_bond_type_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(0u8, "zero-order bond")]
    #[case(5u8, "query bond type")]
    #[case(9u8, "extended bond type")]
    fn test_convert_bond_type_code_invalid(#[case] code: u8, #[case] desc: &str) {
        assert!(
            convert_bond_type_code(code).is_err(),
            "{} should have failed",
            desc
        );
    }

    #[rstest]
    #[case(1u8, BondOrder::Single)]
    #[case(2u8, BondOrder::Double)]
    #[case(3u8, BondOrder::Triple)]
    #[case(4u8, BondOrder::Aromatic)]
    fn test_convert_extended_bond_type_code(#[case] code: u8, #[case] expected: BondOrder) {
        assert_eq!(
            convert_extended_bond_type_code(code, false, false).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case(0u8, "zero-order bond")]
    #[case(9u8, "extended bond type")]
    #[case(12u8, "bond order outside range")]
    fn test_convert_extended_bond_type_code_invalid(#[case] code: u8, #[case] desc: &str) {
        assert!(
            convert_extended_bond_type_code(code, false, false).is_err(),
            "{} should have failed",
            desc
        );
    }

    #[rstest]
    #[case(0u8, BondOrder::Zero)]
    #[case(1u8, BondOrder::Single)]
    #[case(9u8, BondOrder::Quadruple)]
    #[case(10u8, BondOrder::Quintuple)]
    #[case(11u8, BondOrder::Sextuple)]
    fn test_convert_extended_bond_type_code_extended(#[case] code: u8, #[case] expected: BondOrder) {
        assert_eq!(
            convert_extended_bond_type_code(code, true, false).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case(5u8, BondOrder::SingleOrDouble)]
    #[case(6u8, BondOrder::SingleOrAromatic)]
    #[case(7u8, BondOrder::DoubleOrAromatic)]
    #[case(8u8, BondOrder::Any)]
    fn test_convert_extended_bond_type_code_queries(#[case] code: u8, #[case] expected: BondOrder) {
        assert_eq!(
            convert_extended_bond_type_code(code, false, true).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case(0, (Some(BondStereo::Either), Some(BondDir::Either)), (None, None))]
    #[case(1, (Some(BondStereo::Cis), Some(BondDir::Up)), (Some(BondStereo::Cis), Some(BondDir::Up)))]
    #[case(3, (Some(BondStereo::Either), Some(BondDir::Either)), (Some(BondStereo::Either), Some(BondDir::Either)))]
    #[case(4, (Some(BondStereo::Either), Some(BondDir::Either)), (Some(BondStereo::Either), Some(BondDir::Either)))]
    #[case(6, (Some(BondStereo::Trans), Some(BondDir::Down)), (Some(BondStereo::Trans), Some(BondDir::Down)))]
    fn test_convert_bond_stereo_dir_code(
        #[case] code: u8,
        #[case] expected_default: (Option<BondStereo>, Option<BondDir>),
        #[case] expected_basic: (Option<BondStereo>, Option<BondDir>),
    ) {
        assert_eq!(
            convert_bond_stereo_dir_code(code, true).unwrap(),
            expected_default
        );
        assert_eq!(
            convert_bond_stereo_dir_code(code, false).unwrap(),
            expected_basic
        );
    }

    #[rstest]
    #[case(2)]
    #[case(5)]
    #[case(7)]
    fn test_convert_bond_stereo_dir_code_invalid(#[case] code: u8) {
        assert!(convert_bond_stereo_dir_code(code, true).is_err());
    }

    #[rstest]
    #[case(0, Some(BondTopology::Either), None)]
    #[case(1, Some(BondTopology::Ring), Some(BondTopology::Ring))]
    #[case(2, Some(BondTopology::Chain), Some(BondTopology::Chain))]
    fn test_convert_bond_topology_code(
        #[case] code: u8,
        #[case] expected_default: Option<BondTopology>,
        #[case] expected_basic: Option<BondTopology>,
    ) {
        assert_eq!(
            convert_bond_topology_code(code, true).unwrap(),
            expected_default
        );
        assert_eq!(
            convert_bond_topology_code(code, false).unwrap(),
            expected_basic
        );
    }

    #[rstest]
    #[case(3)]
    fn test_convert_bond_topology_code_invalid(#[case] code: u8) {
        assert!(convert_bond_topology_code(code, true).is_err());
    }

    #[rstest]
    #[case(0, Some(BondReactingCenter::UNMARKED))]
    #[case(-1, Some(BondReactingCenter::NOT_CENTER))]
    #[case(1, Some(BondReactingCenter::CENTER))]
    #[case(2, Some(BondReactingCenter::NO_CHANGE))]
    #[case(4, Some(BondReactingCenter::MADE_BROKEN))]
    #[case(8, Some(BondReactingCenter::ORDER_CHANGED))]
    #[case(5, Some(BondReactingCenter::CENTER | BondReactingCenter::MADE_BROKEN))]
    #[case(9, Some(BondReactingCenter::CENTER | BondReactingCenter::ORDER_CHANGED))]
    #[case(12, Some(BondReactingCenter::MADE_BROKEN | BondReactingCenter::ORDER_CHANGED))]
    #[case(13, Some(
        BondReactingCenter::CENTER
            | BondReactingCenter::MADE_BROKEN
            | BondReactingCenter::ORDER_CHANGED
    ))]
    fn test_convert_bond_reacting_center_code(
        #[case] code: i8,
        #[case] expected: Option<BondReactingCenter>,
    ) {
        assert_eq!(
            convert_bond_reacting_center_code(code, true, false).unwrap(),
            expected
        );
        let expected_default = if code == 0 { None } else { expected };
        assert_eq!(
            convert_bond_reacting_center_code(code, false, false).unwrap(),
            expected_default
        );
    }

    #[rstest]
    #[case(3, Some(BondReactingCenter::NO_CHANGE))]
    #[case(15, Some(BondReactingCenter::NO_CHANGE))]
    fn test_convert_bond_reacting_center_code_extended(
        #[case] code: i8,
        #[case] expected: Option<BondReactingCenter>,
    ) {
        assert_eq!(
            convert_bond_reacting_center_code(code, true, true).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case(3)]
    #[case(16)]
    #[case(-2)]
    fn test_convert_bond_reacting_center_code_invalid(#[case] code: i8) {
        assert!(convert_bond_reacting_center_code(code, true, false).is_err());
    }

    #[rstest]
    #[case(0, None)]
    #[case(1, Some(AtomRadical::Singlet))]
    #[case(2, Some(AtomRadical::Doublet))]
    #[case(3, Some(AtomRadical::Triplet))]
    fn test_convert_radical_type_code(#[case] code: u8, #[case] expected: Option<AtomRadical>) {
        assert_eq!(convert_radical_type_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(4)]
    #[case(16)]
    fn test_convert_radical_type_code_invalid(#[case] code: u8) {
        assert!(convert_radical_type_code(code).is_err());
    }

    #[rstest]
    #[case(-2, Some(RingBondCount::AsDrawn))]
    #[case(-1, Some(RingBondCount::NoRingBonds))]
    #[case(2, Some(RingBondCount::R2))]
    #[case(3, Some(RingBondCount::R3))]
    #[case(4, Some(RingBondCount::R4Plus))]
    fn test_convert_ring_bond_count_code(
        #[case] code: i8,
        #[case] expected: Option<RingBondCount>,
    ) {
        assert_eq!(convert_ring_bond_count_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(5)]
    #[case(16)]
    fn test_convert_ring_bond_count_code_invalid(#[case] code: i8) {
        assert!(convert_ring_bond_count_code(code).is_err());
    }

    #[rstest]
    #[case(-2, Some(SubstitutionCount::AsDrawn))]
    #[case(-1, Some(SubstitutionCount::NoSubstitution))]
    #[case(1, Some(SubstitutionCount::S1))]
    #[case(6, Some(SubstitutionCount::S6Plus))]
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
    #[case(7)]
    #[case(16)]
    fn test_convert_substitution_count_code_invalid(#[case] code: i8) {
        assert!(convert_substitution_count_code(code, false).is_err());
    }

    #[rstest]
    #[case(6, Some(SubstitutionCount::S6))]
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
    #[case(0, None)]
    #[case(1, Some(UnsaturatedAtom))]
    fn test_convert_unsaturated_atom_code(
        #[case] code: u8,
        #[case] expected: Option<UnsaturatedAtom>,
    ) {
        assert_eq!(convert_unsaturated_atom_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(2)]
    fn test_convert_unsaturated_atom_code_invalid(#[case] code: u8) {
        assert!(convert_unsaturated_atom_code(code).is_err());
    }

    #[rstest]
    #[case(0, None)]
    #[case(1, Some(AttachmentPointType::First))]
    #[case(2, Some(AttachmentPointType::Second))]
    #[case(3, Some(AttachmentPointType::Both))]
    fn test_convert_attachment_point_code(
        #[case] code: u8,
        #[case] expected: Option<AttachmentPointType>,
    ) {
        assert_eq!(convert_attachment_point_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(4)]
    #[case(16)]
    fn test_convert_attachment_point_code_invalid(#[case] code: u8) {
        assert!(convert_attachment_point_code(code).is_err());
    }
}
