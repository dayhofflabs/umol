//! Bond parsing functions for MOL files.

use crate::{BondDir, BondStereo, BondType};
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
}
