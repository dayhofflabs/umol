//! Atom parsing functions for MOL files.

use crate::AtomStereoParity;
use fixed_width::{from_bytes_with_fields, FieldSet};
use umol::error::FormatError;
use umol::{Error, Result};
use umol_data::Element;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum AtomSymbol {
    Element(Element),
    AtomList,
    Unspecified(char),
    LonePair,
    RGroup(u8),
}

/// Parse atom, bond, or sgroup index, adjusting for 1-based indexing.
pub(crate) fn parse_index(s: &[u8]) -> Result<usize> {
    let index1 = from_bytes_with_fields::<usize>(
        s,
        FieldSet::Seq(vec![FieldSet::new_field(0..s.len()).name("index")]),
    )
    .map_err(|_| {
        Error::from(FormatError::InvalidMolFormat(format!(
            "Invalid index entry: {}",
            String::from_utf8_lossy(s)
        )))
    })?;

    let index = index1.checked_sub(1).ok_or_else(|| {
        Error::from(FormatError::InvalidMolFormat(format!(
            "Invalid index '{}'",
            index1
        )))
    })?;

    Ok(index)
}

/// Parse atom symbol (element, atom list, unspecified, lone pair, R group).
pub(crate) fn parse_atom_symbol(s: &[u8]) -> Result<AtomSymbol> {
    let symbol = from_bytes_with_fields::<String>(
        s,
        FieldSet::Seq(vec![FieldSet::new_field(0..s.len()).name("symbol")]),
    );

    match symbol {
        Ok(symbol) => match symbol.as_str() {
            "L" => Ok(AtomSymbol::AtomList),
            "A" => Ok(AtomSymbol::Unspecified('A')),
            "Q" => Ok(AtomSymbol::Unspecified('Q')),
            "*" => Ok(AtomSymbol::Unspecified('*')),
            "LP" => Ok(AtomSymbol::LonePair),
            element if Element::is_element(element) => Element::from_symbol(element)
                .map(AtomSymbol::Element)
                .ok_or_else(|| {
                    FormatError::InvalidMolFormat(format!("Invalid element symbol: '{}'", element))
                        .into()
                }),
            rgroup if rgroup.bytes().nth(0) == Some(b'R') => {
                let idx = parse_index(rgroup[1..].as_bytes())?;
                Ok(AtomSymbol::RGroup(idx as u8))
            }
            _ => Err(FormatError::InvalidMolFormat(format!(
                "Invalid atom symbol: '{}'",
                String::from_utf8_lossy(s)
            ))
            .into()),
        },
        Err(_) => Err(FormatError::InvalidMolFormat(format!(
            "Invalid atom symbol: '{}'",
            String::from_utf8_lossy(s)
        ))
        .into()),
    }
}

/// Parse valence code (default, explicit, explicit zero).
pub(crate) fn parse_valence_code(code: u8) -> Result<Option<u8>> {
    // 'vvv' field - 0 = default, 1-14 = explicit, 15 = explicit 0
    match code {
        0 => Ok(None),                   // Code 0 means default/unspecified valence
        v @ 1..=14 => Ok(Some(v as u8)), // Codes 1-14 are explicit valences
        15 => Ok(Some(0)),               // Code 15 means explicit zero valence
        _ => Err(FormatError::InvalidMolFormat(format!("Invalid valence code '{}'", code)).into()),
    }
}

/// Parse charge code (uncharged, +3, +2, +1, doublet radical, -1, -2, -3).
pub(crate) fn parse_charge_code(code: u8) -> Result<i8> {
    // 'ccc' field - 0 = uncharged, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
    match code {
        0 => Ok(0),
        1 => Ok(3),
        2 => Ok(2),
        3 => Ok(1),
        4 => Ok(0), // Code 4 is doublet radical, not charge. Set charge = 0 here, see also parse_radical_code
        5 => Ok(-1),
        6 => Ok(-2),
        7 => Ok(-3),
        _ => Err(FormatError::InvalidMolFormat(format!("Invalid charge code '{}'", code)).into()),
    }
}

/// Parse radical code (doublet radical, other charge codes are not radicals by default).
pub(crate) fn parse_radical_code(code: u8) -> Result<Option<u8>> {
    // Check 'ccc' field specifically for radical code 4
    match code {
        4 => Ok(Some(2)),                      // Code 4 is doublet radical
        0 | 1 | 2 | 3 | 5 | 6 | 7 => Ok(None), // Other valid charge codes are not radicals by default
        _ => Err(FormatError::InvalidMolFormat(format!(
            "Invalid code '{}' passed to parse_radical_code",
            code
        ))
        .into()),
    }
}

/// Parse stereo parity code (not stereo, odd, even, either or unmarked).
pub(crate) fn parse_stereo_parity_code(code: u8) -> Result<Option<AtomStereoParity>> {
    // 'sss' field - 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked
    match code {
        0 => Ok(None),
        1 => Ok(Some(AtomStereoParity::Odd)),
        2 => Ok(Some(AtomStereoParity::Even)),
        3 => Ok(Some(AtomStereoParity::Either)), // Treat 'either or unmarked' as Either
        _ => Err(
            FormatError::InvalidMolFormat(format!("Invalid stereo parity code '{}'", code)).into(),
        ),
    }
}

/// Parse hydrogen count code (non-standard: 0 in non-query atoms).
pub(crate) fn parse_hydrogen_count_code(code: u8) -> Result<Option<u8>> {
    match code {
        0 => Ok(None), // non-standard: 0 in non-query atoms
        1 => Ok(Some(0)),
        2 => Ok(Some(1)),
        3 => Ok(Some(2)),
        4 => Ok(Some(3)),
        5 => Ok(Some(4)),
        _ => Err(
            FormatError::InvalidMolFormat(format!("Invalid hydrogen count code '{}'", code)).into(),
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_index() {
        assert_eq!(parse_index(b"  1").unwrap(), 0);
        assert_eq!(parse_index(b" 1 ").unwrap(), 0);
        assert_eq!(parse_index(b"1  ").unwrap(), 0);
        assert_eq!(parse_index(b"1 ").unwrap(), 0);
        assert_eq!(parse_index(b"1").unwrap(), 0);
        assert_eq!(parse_index(b"  2").unwrap(), 1);
        assert_eq!(parse_index(b"  3").unwrap(), 2);
        assert_eq!(parse_index(b"001").unwrap(), 0);
        assert!(parse_index(b"0").is_err());
        assert!(parse_index(b"a").is_err());
    }

    #[test]
    fn test_parse_atom_symbol() {
        assert_eq!(
            parse_atom_symbol(b"H").unwrap(),
            AtomSymbol::Element(Element::H)
        );
        assert_eq!(
            parse_atom_symbol(b"h").unwrap(),
            AtomSymbol::Element(Element::H)
        );
        assert_eq!(
            parse_atom_symbol(b"Cu").unwrap(),
            AtomSymbol::Element(Element::Cu)
        );
        assert_eq!(
            parse_atom_symbol(b"CU").unwrap(),
            AtomSymbol::Element(Element::Cu)
        );
        assert_eq!(
            parse_atom_symbol(b"cu").unwrap(),
            AtomSymbol::Element(Element::Cu)
        );
        assert_eq!(parse_atom_symbol(b"L").unwrap(), AtomSymbol::AtomList);
        assert_eq!(
            parse_atom_symbol(b"A").unwrap(),
            AtomSymbol::Unspecified('A')
        );
        assert_eq!(
            parse_atom_symbol(b"Q").unwrap(),
            AtomSymbol::Unspecified('Q')
        );
        assert_eq!(
            parse_atom_symbol(b"*").unwrap(),
            AtomSymbol::Unspecified('*')
        );
        assert_eq!(parse_atom_symbol(b"LP").unwrap(), AtomSymbol::LonePair);
        assert_eq!(parse_atom_symbol(b"R1").unwrap(), AtomSymbol::RGroup(0));
        assert_eq!(parse_atom_symbol(b"R10").unwrap(), AtomSymbol::RGroup(9));
        assert!(parse_atom_symbol(b"X").is_err());
        assert!(parse_atom_symbol(b"R0").is_err());
    }

    #[test]
    fn test_parse_valence_code() {
        assert_eq!(parse_valence_code(0).unwrap(), None);
        assert_eq!(parse_valence_code(1).unwrap(), Some(1));
        assert_eq!(parse_valence_code(14).unwrap(), Some(14));
        assert_eq!(parse_valence_code(15).unwrap(), Some(0));
        assert!(parse_valence_code(16).is_err());
    }

    #[test]
    fn test_parse_charge_code() {
        assert_eq!(parse_charge_code(0).unwrap(), 0);
        assert_eq!(parse_charge_code(1).unwrap(), 3);
        assert_eq!(parse_charge_code(4).unwrap(), 0);
        assert_eq!(parse_charge_code(5).unwrap(), -1);
        assert!(parse_charge_code(8).is_err());
    }

    #[test]
    fn test_parse_radical_code() {
        assert_eq!(parse_radical_code(4).unwrap(), Some(2));
        assert_eq!(parse_radical_code(1).unwrap(), None);
        assert!(parse_radical_code(8).is_err());
    }

    #[test]
    fn test_parse_stereo_parity_code() {
        assert_eq!(parse_stereo_parity_code(0).unwrap(), None);
        assert_eq!(
            parse_stereo_parity_code(1).unwrap(),
            Some(AtomStereoParity::Odd)
        );
        assert_eq!(
            parse_stereo_parity_code(2).unwrap(),
            Some(AtomStereoParity::Even)
        );
        assert_eq!(
            parse_stereo_parity_code(3).unwrap(),
            Some(AtomStereoParity::Either)
        );
    }

    #[test]
    fn test_parse_hydrogen_count_code() {
        assert_eq!(parse_hydrogen_count_code(0).unwrap(), None);
        assert_eq!(parse_hydrogen_count_code(1).unwrap(), Some(0));
        assert!(parse_hydrogen_count_code(6).is_err());
    }
}
