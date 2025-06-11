//! Bond parsing functions for MOL files.

use crate::{BondDir, BondReactingCenter, BondStereo, BondTopology, BondType};
use umol::error::FormatError;
use umol::Result;

/// Parse bond type code (single, double, triple, aromatic, single or double, single or aromatic, double or aromatic, any).
pub(crate) fn parse_bond_type_code(code: u8) -> Result<BondType> {
    // 'ttt' field - bond type (1=Single, 2=Double, etc.)
    match code {
        1 => Ok(BondType::Single),
        2 => Ok(BondType::Double),
        3 => Ok(BondType::Triple),
        4 => Ok(BondType::Aromatic),
        5 => Ok(BondType::SingleOrDouble),
        6 => Ok(BondType::SingleOrAromatic),
        7 => Ok(BondType::DoubleOrAromatic),
        8 => Ok(BondType::Any),
        _ => {
            Err(FormatError::InvalidMolFormat(format!("Invalid bond type code '{}'", code)).into())
        }
    }
}

/// Parse bond stereo code (not stereo, cis, either, trans).
pub(crate) fn parse_bond_stereo_code(code: u8) -> Result<Option<BondStereo>> {
    // 'sss' field - bond stereo (0=Not stereo, 1=Up, 3=Either, 6=Down)
    // Used for double bond Cis/Trans by convention
    match code {
        0 => Ok(None),                     // Not stereo
        1 => Ok(Some(BondStereo::Cis)),    // Up (used for Cis/Trans by convention)
        3 => Ok(Some(BondStereo::Either)), // Either
        6 => Ok(Some(BondStereo::Trans)),  // Down (used for Cis/Trans by convention)
        _ => Err(
            FormatError::InvalidMolFormat(format!("Invalid bond stereo code '{}'", code)).into(),
        ),
    }
}

/// Parse bond direction code (not stereo, wedge, dash).
pub(crate) fn parse_bond_dir_code(code: u8) -> Result<Option<BondDir>> {
    // 'sss' field - bond stereo (0=Not stereo, 1=Up, 6=Down)
    // Used for single bond Wedge/Dash
    match code {
        0 => Ok(None),                 // Not stereo
        1 => Ok(Some(BondDir::Wedge)), // Up
        6 => Ok(Some(BondDir::Dash)),  // Down
        // Note: Codes 3 (Either) and 4 (Unknown) exist but aren't typically used for wedge/dash dir
        _ => Err(
            FormatError::InvalidMolFormat(format!("Invalid bond direction code '{}'", code)).into(),
        ),
    }
}

/// Parse bond topology code (chain, ring, either).
pub(crate) fn parse_bond_topology_code(code: u8) -> Result<Option<BondTopology>> {
    // 'rrr' field - bond topology (0=Either, 1=Ring, 2=Chain)
    match code {
        0 => Ok(Some(BondTopology::Either)),
        1 => Ok(Some(BondTopology::Ring)),
        2 => Ok(Some(BondTopology::Chain)),
        _ => Err(
            FormatError::InvalidMolFormat(format!("Invalid bond topology code '{}'", code)).into(),
        ),
    }
}

/// Parse reacting center code (not reacting, reacting).
pub(crate) fn parse_reacting_center_code(code: i8) -> Result<Option<BondReactingCenter>> {
    if code == 0 {
        return Ok(Some(BondReactingCenter::UNMARKED));
    }
    if code == -1 {
        return Ok(Some(BondReactingCenter::NOT_CENTER));
    }
    if code < -1 || code > 15 {
        return Err(FormatError::InvalidMolFormat(format!(
            "Invalid reacting center code \'{}\'",
            code
        ))
        .into());
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
            return Err(FormatError::InvalidMolFormat(format!(
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

    #[test]
    fn test_parse_bond_type_code() {
        assert_eq!(parse_bond_type_code(1).unwrap(), BondType::Single);
        assert_eq!(parse_bond_type_code(2).unwrap(), BondType::Double);
        assert_eq!(parse_bond_type_code(3).unwrap(), BondType::Triple);
        assert_eq!(parse_bond_type_code(4).unwrap(), BondType::Aromatic);
        assert_eq!(parse_bond_type_code(5).unwrap(), BondType::SingleOrDouble);
        assert_eq!(parse_bond_type_code(6).unwrap(), BondType::SingleOrAromatic);
        assert_eq!(parse_bond_type_code(7).unwrap(), BondType::DoubleOrAromatic);
        assert_eq!(parse_bond_type_code(8).unwrap(), BondType::Any);
        assert!(parse_bond_type_code(9).is_err());
    }

    #[test]
    fn test_parse_bond_stereo_code() {
        assert_eq!(parse_bond_stereo_code(0).unwrap(), None);
        assert_eq!(parse_bond_stereo_code(1).unwrap(), Some(BondStereo::Cis));
        assert_eq!(parse_bond_stereo_code(3).unwrap(), Some(BondStereo::Either));
        assert_eq!(parse_bond_stereo_code(6).unwrap(), Some(BondStereo::Trans));
        assert!(parse_bond_stereo_code(2).is_err());
    }

    #[test]
    fn test_parse_bond_dir_code() {
        assert_eq!(parse_bond_dir_code(0).unwrap(), None);
        assert_eq!(parse_bond_dir_code(1).unwrap(), Some(BondDir::Wedge));
        assert_eq!(parse_bond_dir_code(6).unwrap(), Some(BondDir::Dash));
        assert!(parse_bond_dir_code(3).is_err());
    }

    #[test]
    fn test_parse_bond_topology_code() {
        assert_eq!(
            parse_bond_topology_code(0).unwrap(),
            Some(BondTopology::Either)
        );
        assert_eq!(
            parse_bond_topology_code(1).unwrap(),
            Some(BondTopology::Ring)
        );
        assert_eq!(
            parse_bond_topology_code(2).unwrap(),
            Some(BondTopology::Chain)
        );
        assert!(parse_bond_topology_code(3).is_err());
    }

    #[test]
    fn test_parse_reacting_center_code() {
        // Test individual flags
        assert_eq!(
            parse_reacting_center_code(0).unwrap(),
            Some(BondReactingCenter::UNMARKED)
        );
        assert_eq!(
            parse_reacting_center_code(-1).unwrap(),
            Some(BondReactingCenter::NOT_CENTER)
        );
        assert_eq!(
            parse_reacting_center_code(1).unwrap(),
            Some(BondReactingCenter::CENTER)
        );
        assert_eq!(
            parse_reacting_center_code(2).unwrap(),
            Some(BondReactingCenter::NO_CHANGE)
        );
        assert_eq!(
            parse_reacting_center_code(4).unwrap(),
            Some(BondReactingCenter::MADE_BROKEN)
        );
        assert_eq!(
            parse_reacting_center_code(8).unwrap(),
            Some(BondReactingCenter::ORDER_CHANGED)
        );

        // Test combinations from the spec
        // 5 = (4 + 1) -> MADE_BROKEN | CENTER
        assert_eq!(
            parse_reacting_center_code(5).unwrap(),
            Some(BondReactingCenter::CENTER | BondReactingCenter::MADE_BROKEN)
        );
        // 9 = (8 + 1) -> ORDER_CHANGED | CENTER
        assert_eq!(
            parse_reacting_center_code(9).unwrap(),
            Some(BondReactingCenter::CENTER | BondReactingCenter::ORDER_CHANGED)
        );
        // 12 = (4 + 8) -> MADE_BROKEN | ORDER_CHANGED
        assert_eq!(
            parse_reacting_center_code(12).unwrap(),
            Some(BondReactingCenter::MADE_BROKEN | BondReactingCenter::ORDER_CHANGED)
        );
        // 13 = (4 + 8 + 1) -> MADE_BROKEN | ORDER_CHANGED | CENTER
        assert_eq!(
            parse_reacting_center_code(13).unwrap(),
            Some(
                BondReactingCenter::CENTER
                    | BondReactingCenter::MADE_BROKEN
                    | BondReactingCenter::ORDER_CHANGED
            )
        );

        // Test invalid codes
        assert!(parse_reacting_center_code(3).is_err());
        assert!(parse_reacting_center_code(16).is_err());
        assert!(parse_reacting_center_code(-2).is_err());
    }
}
