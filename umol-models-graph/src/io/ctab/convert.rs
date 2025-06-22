//! Convert numerical codes used in MOL files to enums

use crate::atom::AtomStereoParity;
use crate::bond::{BondDir, BondReactingCenter, BondStereo, BondTopology, BondType};
use umol::error::{ParseError, Result};

/// Convert atom mass difference code (atom block)
/// 'dd' field: mass difference (-3..=4), 0 if value outside of this range
pub(crate) fn convert_atom_mass_diff_code(code: i8) -> Result<Option<i8>> {
    match code {
        0 => Ok(None),
        -3..=-1 | 1..=4 => Ok(Some(code)),
        _ => Err(ParseError::Invalid(format!("Invalid mass difference code '{}'", code)).into()),
    }
}

/// Convert atom valence code (default, explicit, explicit zero).
// 'vvv' field: 0 = default, 1..=14 = explicit, 15 = explicit 0
pub(crate) fn convert_atom_valence_code(code: u8) -> Result<Option<u8>> {
    match code {
        0 => Ok(None),                   // default/unspecified valence
        v @ 1..=14 => Ok(Some(v as u8)), // explicit valences
        15 => Ok(Some(0)),               // explicit zero valence
        _ => Err(ParseError::Invalid(format!("Invalid valence code '{}'", code)).into()),
    }
}

/// Convert atom charge code
/// 'ccc' field: 0 = uncharged, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
pub(crate) fn convert_atom_charge_code(code: u8) -> Result<i8> {
    match code {
        0 => Ok(0),
        1..=3 | 5..=7 => Ok(4 - code as i8),
        4 => Ok(0), // Code 4 is doublet radical, not charge
        _ => Err(ParseError::Invalid(format!("Invalid charge code '{}'", code)).into()),
    }
}

/// Convert atom radical code
/// 'ccc' field: 4 = doublet radical
pub(crate) fn convert_atom_radical_code(code: u8) -> Result<Option<u8>> {
    match code {
        4 => Ok(Some(2)), // Code 4 is doublet radical
        0..=3 | 5..=7 => Ok(None),
        _ => Err(ParseError::Invalid(format!(
            "Invalid code '{}' passed to parse_radical_code",
            code
        ))
        .into()),
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

/// Convert bond stereo code
/// 'sss' field - bond stereo (0=Not stereo, 1=Up, 3=Either, 4=Unknown, 6=Down)
/// Used for double bond Cis/Trans by convention
pub(crate) fn convert_bond_stereo_code(code: u8) -> Result<Option<BondStereo>> {
    match code {
        0 => Ok(None),                     // Not stereo
        1 => Ok(Some(BondStereo::Cis)),    // Up (used for Cis/Trans by convention)
        3 | 4 => Ok(Some(BondStereo::Either)), // Either or Unknown
        6 => Ok(Some(BondStereo::Trans)),  // Down (used for Cis/Trans by convention)
        _ => Err(ParseError::Invalid(format!("Invalid bond stereo code '{}'", code)).into()),
    }
}

/// Convert bond direction code
/// 'sss' field - bond direction (0=Not stereo, 1=Up, 6=Down)
/// Used for single bond Wedge/Dash
pub(crate) fn convert_bond_dir_code(code: u8) -> Result<Option<BondDir>> {
    match code {
        0 => Ok(None),                 // Not stereo
        1 => Ok(Some(BondDir::Wedge)), // Up
        6 => Ok(Some(BondDir::Dash)),  // Down
        // Note: Codes 3 (Either) and 4 (Unknown) exist but aren't typically used for wedge/dash dir
        _ => Err(ParseError::Invalid(format!("Invalid bond direction code '{}'", code)).into()),
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

    #[rstest]
    #[case(0, None)]
    #[case(-3, Some(-3))]
    #[case(4, Some(4))]
    fn test_atom_convert_mass_diff_code(#[case] code: i8, #[case] expected: Option<i8>) {
        assert_eq!(convert_atom_mass_diff_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(5, "too high")]
    fn test_atom_convert_mass_diff_code_invalid(#[case] code: i8, #[case] expected: &str) {
        assert!(
            convert_atom_mass_diff_code(code).is_err(),
            "{} should have failed",
            expected
        );
    }

    #[rstest]
    #[case(0, None)]
    #[case(1, Some(1))]
    #[case(14, Some(14))]
    #[case(15, Some(0))]
    fn test_atom_convert_valence_code(#[case] code: u8, #[case] expected: Option<u8>) {
        assert_eq!(convert_atom_valence_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(16, "too high")]
    fn test_atom_convert_valence_code_invalid(#[case] code: u8, #[case] expected: &str) {
        assert!(
            convert_atom_valence_code(code).is_err(),
            "{} should have failed",
            expected
        );
    }

    #[rstest]
    #[case(0, 0)]
    #[case(1, 3)]
    #[case(4, 0)]
    #[case(5, -1)]
    fn test_atom_convert_charge_code(#[case] code: u8, #[case] expected: i8) {
        assert_eq!(convert_atom_charge_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(8, "too high")]
    fn test_atom_convert_charge_code_invalid(#[case] code: u8, #[case] expected: &str) {
        assert!(
            convert_atom_charge_code(code).is_err(),
            "{} should have failed",
            expected
        );
    }

    #[rstest]
    #[case(4, Some(2))]
    #[case(1, None)]
    fn test_atom_convert_radical_code(#[case] code: u8, #[case] expected: Option<u8>) {
        assert_eq!(convert_atom_radical_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(8, "too high")]
    fn test_atom_convert_radical_code_invalid(#[case] code: u8, #[case] expected: &str) {
        assert!(
            convert_atom_radical_code(code).is_err(),
            "{} should have failed",
            expected
        );
    }

    #[rstest]
    #[case(0, None)]
    #[case(1, Some(AtomStereoParity::Odd))]
    #[case(2, Some(AtomStereoParity::Even))]
    #[case(3, Some(AtomStereoParity::Either))]
    fn test_atom_convert_stereo_parity_code(
        #[case] code: u8,
        #[case] expected: Option<AtomStereoParity>,
    ) {
        assert_eq!(convert_atom_stereo_parity_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(4, "too high")]
    fn test_atom_convert_stereo_parity_code_invalid(#[case] code: u8, #[case] expected: &str) {
        assert!(
            convert_atom_stereo_parity_code(code).is_err(),
            "{} should have failed",
            expected
        );
    }

    #[rstest]
    #[case(0, None)]
    #[case(1, Some(0))]
    #[case(2, Some(1))]
    #[case(3, Some(2))]
    #[case(4, Some(3))]
    #[case(5, Some(4))]
    fn test_atom_convert_hydrogen_count_code(#[case] code: u8, #[case] expected: Option<u8>) {
        assert_eq!(convert_atom_hydrogen_count_code(code).unwrap(), expected);
    }

    #[rstest]
    #[case(6, "too high")]
    fn test_atom_convert_hydrogen_count_code_invalid(#[case] code: u8, #[case] expected: &str) {
        assert!(
            convert_atom_hydrogen_count_code(code).is_err(),
            "{} should have failed",
            expected
        );
    }

    #[test]
    fn test_bond_convert_type_code() {
        assert_eq!(convert_bond_type_code(1).unwrap(), BondType::Single);
        assert_eq!(convert_bond_type_code(2).unwrap(), BondType::Double);
        assert_eq!(convert_bond_type_code(3).unwrap(), BondType::Triple);
        assert_eq!(convert_bond_type_code(4).unwrap(), BondType::Aromatic);
        assert_eq!(convert_bond_type_code(5).unwrap(), BondType::SingleOrDouble);
        assert_eq!(convert_bond_type_code(6).unwrap(), BondType::SingleOrAromatic);
        assert_eq!(convert_bond_type_code(7).unwrap(), BondType::DoubleOrAromatic);
        assert_eq!(convert_bond_type_code(8).unwrap(), BondType::Any);
        assert!(convert_bond_type_code(9).is_err());
    }

    #[test]
    fn test_bond_convert_stereo_code() {
        assert_eq!(convert_bond_stereo_code(0).unwrap(), None);
        assert_eq!(convert_bond_stereo_code(1).unwrap(), Some(BondStereo::Cis));
        assert_eq!(convert_bond_stereo_code(3).unwrap(), Some(BondStereo::Either));
        assert_eq!(convert_bond_stereo_code(6).unwrap(), Some(BondStereo::Trans));
        assert!(convert_bond_stereo_code(2).is_err());
    }

    #[test]
    fn test_bond_convert_dir_code() {
        assert_eq!(convert_bond_dir_code(0).unwrap(), None);
        assert_eq!(convert_bond_dir_code(1).unwrap(), Some(BondDir::Wedge));
        assert_eq!(convert_bond_dir_code(6).unwrap(), Some(BondDir::Dash));
        assert!(convert_bond_dir_code(3).is_err());
    }

    #[test]
    fn test_bond_convert_topology_code() {
        assert_eq!(
            convert_bond_topology_code(0).unwrap(),
            Some(BondTopology::Either)
        );
        assert_eq!(
            convert_bond_topology_code(1).unwrap(),
            Some(BondTopology::Ring)
        );
        assert_eq!(
            convert_bond_topology_code(2).unwrap(),
            Some(BondTopology::Chain)
        );
        assert!(convert_bond_topology_code(3).is_err());
    }

    #[test]
    fn test_bond_convert_reacting_center_code() {
        // Test individual flags
        assert_eq!(
            convert_bond_reacting_center_code(0).unwrap(),
            Some(BondReactingCenter::UNMARKED)
        );
        assert_eq!(
            convert_bond_reacting_center_code(-1).unwrap(),
            Some(BondReactingCenter::NOT_CENTER)
        );
        assert_eq!(
            convert_bond_reacting_center_code(1).unwrap(),
            Some(BondReactingCenter::CENTER)
        );
        assert_eq!(
            convert_bond_reacting_center_code(2).unwrap(),
            Some(BondReactingCenter::NO_CHANGE)
        );
        assert_eq!(
            convert_bond_reacting_center_code(4).unwrap(),
            Some(BondReactingCenter::MADE_BROKEN)
        );
        assert_eq!(
            convert_bond_reacting_center_code(8).unwrap(),
            Some(BondReactingCenter::ORDER_CHANGED)
        );

        // Test combinations from the spec
        // 5 = (4 + 1) -> MADE_BROKEN | CENTER
        assert_eq!(
            convert_bond_reacting_center_code(5).unwrap(),
            Some(BondReactingCenter::CENTER | BondReactingCenter::MADE_BROKEN)
        );
        // 9 = (8 + 1) -> ORDER_CHANGED | CENTER
        assert_eq!(
            convert_bond_reacting_center_code(9).unwrap(),
            Some(BondReactingCenter::CENTER | BondReactingCenter::ORDER_CHANGED)
        );
        // 12 = (4 + 8) -> MADE_BROKEN | ORDER_CHANGED
        assert_eq!(
            convert_bond_reacting_center_code(12).unwrap(),
            Some(BondReactingCenter::MADE_BROKEN | BondReactingCenter::ORDER_CHANGED)
        );
        // 13 = (4 + 8 + 1) -> MADE_BROKEN | ORDER_CHANGED | CENTER
        assert_eq!(
            convert_bond_reacting_center_code(13).unwrap(),
            Some(
                BondReactingCenter::CENTER
                    | BondReactingCenter::MADE_BROKEN
                    | BondReactingCenter::ORDER_CHANGED
            )
        );

        // Test invalid codes
        assert!(convert_bond_reacting_center_code(3).is_err());
        assert!(convert_bond_reacting_center_code(16).is_err());
        assert!(convert_bond_reacting_center_code(-2).is_err());
    }
}
