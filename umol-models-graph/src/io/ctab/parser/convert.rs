//! Convert numerical codes used in MOL files to enums

use crate::io::ctab::atom::{
    AtomExactChange, AtomInversionRetention, AtomStereoCare, AtomStereoParity, AtomSymbol,
};
use crate::io::ctab::bond::{BondDir, BondReactingCenter, BondStereo, BondTopology, BondType};
use umol::error::{ParseError, Result};
use umol_data::Element;

/// Convert atom mass difference code (atom block)
/// 'dd' field: mass difference (-3..=4), None if 0 or value outside of this range
pub(crate) fn convert_atom_mass_diff_code(code: i8) -> Option<i8> {
    match code {
        -3..=-1 | 1..=4 => Some(code),
        _ => None,
    }
}

/// Convert standard atom symbol and mass difference to element and isotope mass
/// 'ss' field: atom symbol, 'dd' field: mass difference
/// Processes elements and named isotopes.
/// Returns error for non-standard atom symbols (L, A, Q, *, LP, R#)
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
        _ => unreachable!("atom_symbol_standard() should only return Element or NamedIsotope"),
    };
    (element, isotope_mass)
}

/// Convert atom charge code (includes doublet radical).
/// 'ccc' field: 0 = uncharged, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
/// 0 if outside of range.
pub(crate) fn convert_atom_charge_code(code: u8) -> (i8, Option<u8>) {
    match code {
        1..=3 | 5..=7 => (4 - code as i8, None),
        4 => (0, Some(2)),
        _ => (0, None),
    }
}

/// Convert atom stereo parity code (not stereo, odd, even, either or unmarked).
// 'sss' field: 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked
pub(crate) fn convert_atom_stereo_parity_code(code: u8) -> Result<Option<AtomStereoParity>> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomStereoParity::Odd)),
        2 => Ok(Some(AtomStereoParity::Even)),
        3 => Ok(Some(AtomStereoParity::Either)),
        _ => Err(ParseError::Invalid(format!("Invalid stereo parity code '{}'", code)).into()),
    }
}

/// Convert atom hydrogen count code (non-standard: 0 in non-query atoms).
/// 'hhh' field: 0 = non-query atom, 1 = H0, 2 = H1, 3 = H2, 4 = H3, 5 = H4
pub(crate) fn convert_atom_hydrogen_count_code(code: u8) -> Result<Option<u8>> {
    match code {
        0 => Ok(None),
        1..=5 => Ok(Some(code - 1)),
        _ => Err(ParseError::Invalid(format!("Invalid hydrogen count code '{}'", code)).into()),
    }
}

/// Convert atom stereo care box code.
/// 'bbb' field: 0 = ignore stereo, 1 = stereo must match
pub(crate) fn convert_atom_stereo_care_code(code: u8) -> Result<Option<AtomStereoCare>> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomStereoCare::Care)),
        _ => Err(ParseError::Invalid(format!("Invalid stereo care code '{}'", code)).into()),
    }
}

/// Convert atom valence code (default, explicit, explicit zero).
/// 'vvv' field: 0 = default, 1..=14 = explicit, 15 = explicit 0
/// Returns error for invalid valence codes.
pub(crate) fn convert_atom_valence_code(code: u8) -> Result<Option<u8>> {
    match code {
        0 => Ok(None),                   // default/unspecified valence
        v @ 1..=14 => Ok(Some(v as u8)), // explicit valences
        15 => Ok(Some(0)),               // explicit zero valence
        _ => Err(ParseError::Invalid(format!("Invalid valence code '{}'", code)).into()),
    }
}

/// Convert atom inversion flag code.
/// 'nnn' field: 0 = not applicable, 1 = inverted, 2 = retained
pub(crate) fn convert_atom_inversion_flag_code(code: u8) -> Result<Option<AtomInversionRetention>> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomInversionRetention::Inverted)),
        2 => Ok(Some(AtomInversionRetention::Retained)),
        _ => Err(ParseError::Invalid(format!("Invalid inversion flag code '{}'", code)).into()),
    }
}

/// Convert atom exact change flag code.
/// 'eee' field: 0 = change allowed, 1 = exact change required
pub(crate) fn convert_atom_exact_change_flag_code(code: u8) -> Result<Option<AtomExactChange>> {
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomExactChange::Match)),
        _ => Err(ParseError::Invalid(format!("Invalid exact change flag code '{}'", code)).into()),
    }
}

/// Convert bond type code (standard molecules only)
/// 'ttt' field - bond type (1=Single, 2=Double, 3=Triple, 4=Aromatic)
pub(crate) fn convert_bond_type_code_standard(code: u8) -> Result<BondType> {
    match code {
        1 => Ok(BondType::Single),
        2 => Ok(BondType::Double),
        3 => Ok(BondType::Triple),
        4 => Ok(BondType::Aromatic),
        _ => Err(ParseError::Invalid(format!("Invalid standard bond type code '{}'", code)).into()),
    }
}

/// Convert bond type code
/// 'ttt' field - bond type (1=Single, 2=Double, 3=Triple, 4=Aromatic, 5=SingleOrDouble,
/// 6=SingleOrAromatic, 7=DoubleOrAromatic, 8=Any)
pub(crate) fn convert_bond_type_code(code: u8) -> Result<BondType> {
    match code {
        1 => Ok(BondType::Single),
        2 => Ok(BondType::Double),
        3 => Ok(BondType::Triple),
        4 => Ok(BondType::Aromatic),
        5 => Ok(BondType::SingleOrDouble),
        6 => Ok(BondType::SingleOrAromatic),
        7 => Ok(BondType::DoubleOrAromatic),
        8 => Ok(BondType::Any),
        _ => Err(ParseError::Invalid(format!("Invalid bond type code '{}'", code)).into()),
    }
}

/// Convert bond stereo/direction code
/// 'sss' field - can mean stereo for double bonds or direction for single bonds
/// Stereo: (0=Not stereo, 1=Up, 3=Either, 4=Unknown, 6=Down)
/// Direction: (1=Up, 6=Down)
pub(crate) fn convert_bond_stereo_dir_code(
    code: u8,
) -> Result<(Option<BondStereo>, Option<BondDir>)> {
    match code {
        0 => Ok((None, None)),
        1 => Ok((Some(BondStereo::Cis), Some(BondDir::Wedge))),
        3 | 4 => Ok((Some(BondStereo::Either), Some(BondDir::Either))),
        6 => Ok((Some(BondStereo::Trans), Some(BondDir::Dash))),
        _ => Err(
            ParseError::Invalid(format!("Invalid bond stereo/direction code '{}'", code)).into(),
        ),
    }
}

/// Convert bond topology code
/// 'rrr' field - bond topology (0=Either, 1=Ring, 2=Chain)
pub(crate) fn convert_bond_topology_code(code: u8) -> Result<Option<BondTopology>> {
    match code {
        0 => Ok(Some(BondTopology::Either)),
        1 => Ok(Some(BondTopology::Ring)),
        2 => Ok(Some(BondTopology::Chain)),
        _ => Err(ParseError::Invalid(format!("Invalid bond topology code '{}'", code)).into()),
    }
}

/// Convert bond reacting center code
/// 'ccc' field - bond reacting center (0=Not reacting, 1=Reacting, -1=Not a center,
/// 2=No change, 4=Bond made/broken, 8=Bond order changes)
pub(crate) fn convert_bond_reacting_center_code(code: i8) -> Result<Option<BondReactingCenter>> {
    if code == 0 {
        return Ok(Some(BondReactingCenter::UNMARKED));
    }
    if code == -1 {
        return Ok(Some(BondReactingCenter::NOT_CENTER));
    }
    if code < -1 || code > 15 {
        return Err(
            ParseError::Invalid(format!("Invalid reacting center code \'{}\'", code)).into(),
        );
    }

    // Positive codes can be partially combined:
    //   1 = a center
    //   2 = no change (cannot be combined with other flags)
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
        if code != 2 {
            return Err(ParseError::Invalid(format!(
                "Invalid reacting center code combination \'{}\'",
                code
            ))
            .into());
        }
        flags |= BondReactingCenter::NO_CHANGE;
    }

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

    Ok(Some(flags))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use umol_data::NamedIsotope;

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
    #[case(4, 0, Some(2))]
    #[case(5, -1, None)]
    #[case(8, 0, None)]
    fn test_convert_atom_charge_code(
        #[case] code: u8,
        #[case] expected_charge: i8,
        #[case] expected_radical: Option<u8>,
    ) {
        assert_eq!(
            convert_atom_charge_code(code),
            (expected_charge, expected_radical)
        );
    }

    #[rstest]
    #[case(0, None)]
    #[case(1, Some(AtomStereoParity::Odd))]
    #[case(2, Some(AtomStereoParity::Even))]
    #[case(3, Some(AtomStereoParity::Either))]
    fn test_convert_atom_stereo_parity_code(
        #[case] code: u8,
        #[case] expected: Option<AtomStereoParity>,
    ) {
        assert_eq!(convert_atom_stereo_parity_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(4, "too high")]
    fn test_convert_atom_stereo_parity_code_invalid(#[case] code: u8, #[case] desc: &str) {
        assert!(
            convert_atom_stereo_parity_code(code).is_err(),
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
        assert_eq!(convert_atom_hydrogen_count_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(6, "too high")]
    fn test_convert_atom_hydrogen_count_code_invalid(#[case] code: u8, #[case] desc: &str) {
        assert!(
            convert_atom_hydrogen_count_code(code).is_err(),
            "{} should have failed",
            desc
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

    #[test]
    fn test_convert_bond_type_code() {
        assert_eq!(convert_bond_type_code(1).unwrap(), BondType::Single);
        assert_eq!(convert_bond_type_code(2).unwrap(), BondType::Double);
        assert_eq!(convert_bond_type_code(3).unwrap(), BondType::Triple);
        assert_eq!(convert_bond_type_code(4).unwrap(), BondType::Aromatic);
        assert_eq!(convert_bond_type_code(5).unwrap(), BondType::SingleOrDouble);
        assert_eq!(
            convert_bond_type_code(6).unwrap(),
            BondType::SingleOrAromatic
        );
        assert_eq!(
            convert_bond_type_code(7).unwrap(),
            BondType::DoubleOrAromatic
        );
        assert_eq!(convert_bond_type_code(8).unwrap(), BondType::Any);
        assert!(convert_bond_type_code(9).is_err());
    }

    #[rstest]
    #[case(1, BondType::Single)]
    #[case(2, BondType::Double)]
    #[case(3, BondType::Triple)]
    #[case(4, BondType::Aromatic)]
    fn test_convert_bond_type_code_standard(#[case] code: u8, #[case] expected: BondType) {
        assert_eq!(convert_bond_type_code_standard(code).unwrap(), expected);
    }

    #[rstest]
    #[case(0)]
    #[case(5)]
    #[case(8)]
    fn test_convert_bond_type_code_standard_invalid(#[case] code: u8) {
        assert!(convert_bond_type_code_standard(code).is_err());
    }

    #[rstest]
    #[case(0, (None, None))]
    #[case(1, (Some(BondStereo::Cis), Some(BondDir::Wedge)))]
    #[case(3, (Some(BondStereo::Either), Some(BondDir::Either)))]
    #[case(4, (Some(BondStereo::Either), Some(BondDir::Either)))]
    #[case(6, (Some(BondStereo::Trans), Some(BondDir::Dash)))]
    fn test_convert_bond_stereo_dir_code(
        #[case] code: u8,
        #[case] expected: (Option<BondStereo>, Option<BondDir>),
    ) {
        assert_eq!(convert_bond_stereo_dir_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(2)]
    #[case(5)]
    #[case(7)]
    fn test_convert_bond_stereo_dir_code_invalid(#[case] code: u8) {
        assert!(convert_bond_stereo_dir_code(code).is_err());
    }

    #[rstest]
    #[case(0, BondTopology::Either)]
    #[case(1, BondTopology::Ring)]
    #[case(2, BondTopology::Chain)]
    fn test_convert_bond_topology_code(#[case] code: u8, #[case] topology: BondTopology) {
        assert_eq!(convert_bond_topology_code(code).unwrap(), Some(topology));
    }

    #[rstest]
    #[case(3)]
    fn test_convert_bond_topology_code_invalid(#[case] code: u8) {
        assert!(convert_bond_topology_code(code).is_err());
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
        assert_eq!(convert_bond_reacting_center_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(3)]
    #[case(16)]
    #[case(-2)]
    fn test_convert_bond_reacting_center_code_invalid(#[case] code: i8) {
        assert!(convert_bond_reacting_center_code(code).is_err());
    }
}
