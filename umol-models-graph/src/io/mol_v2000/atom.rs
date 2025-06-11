//! Atom parsing functions for MOL files.

use crate::AtomStereoParity;
use fixed_width::{from_bytes_with_fields, FieldSet};
use umol::error::FormatError;
use umol::{Error, Result};
use umol_data::{Element, NamedIsotope};

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

/// Process atom symbol
pub(crate) fn process_atom_symbol(symbol: &str) -> Result<AtomSymbol> {
    match symbol {
        "L" => Ok(AtomSymbol::AtomList),
        "A" => Ok(AtomSymbol::Unspecified('A')),
        "Q" => Ok(AtomSymbol::Unspecified('Q')),
        "*" => Ok(AtomSymbol::Unspecified('*')),
        "LP" => Ok(AtomSymbol::LonePair),
        isotope if NamedIsotope::is_named_isotope(isotope) => NamedIsotope::from_symbol(isotope)
            .map(|isotope| AtomSymbol::Element(isotope.element()))
            .ok_or_else(|| {
                FormatError::InvalidMolFormat(format!(
                    "Invalid named isotope symbol: '{}'",
                    isotope
                ))
                .into()
            }),
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
            symbol
        ))
        .into()),
    }
}

pub(crate) fn process_isotope_symbol(symbol: &str) -> Result<Option<u32>> {
    if NamedIsotope::is_named_isotope(symbol) {
        NamedIsotope::from_symbol(symbol)
            .map(|isotope| Some(isotope.mass_number()))
            .ok_or_else(|| {
                FormatError::InvalidMolFormat(format!("Invalid named isotope symbol: '{}'", symbol))
                    .into()
            })
    } else {
        Ok(None)
    }
}

/// Process mass difference code in atom block (-3, ..., 4).
pub(crate) fn process_mass_diff_code(code: i8) -> Result<Option<i8>> {
    match code {
        0 => Ok(None),
        -3..=-1 | 1..=4 => Ok(Some(code)),
        _ => {
            Err(FormatError::InvalidMolFormat(format!("Invalid mass diff code '{}'", code)).into())
        }
    }
}

/// Process valence code (default, explicit, explicit zero).
pub(crate) fn process_valence_code(code: u8) -> Result<Option<u8>> {
    // 'vvv' field - 0 = default, 1-14 = explicit, 15 = explicit 0
    match code {
        0 => Ok(None),                   // default/unspecified valence
        v @ 1..=14 => Ok(Some(v as u8)), // explicit valences
        15 => Ok(Some(0)),               // explicit zero valence
        _ => Err(FormatError::InvalidMolFormat(format!("Invalid valence code '{}'", code)).into()),
    }
}

/// Process charge code (uncharged, +3, +2, +1, doublet radical, -1, -2, -3).
pub(crate) fn process_charge_code(code: u8) -> Result<i8> {
    // 'ccc' field - 0 = uncharged, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
    match code {
        0 => Ok(0),
        1..=3 | 5..=7 => Ok(4 - code as i8),
        4 => Ok(0), // Code 4 is doublet radical, not charge. Set charge = 0 here, see also parse_radical_code
        _ => Err(FormatError::InvalidMolFormat(format!("Invalid charge code '{}'", code)).into()),
    }
}

/// Process radical code (doublet radical, other charge codes are not radicals by default).
pub(crate) fn process_radical_code(code: u8) -> Result<Option<u8>> {
    // Check 'ccc' field specifically for radical code 4
    match code {
        4 => Ok(Some(2)),          // Code 4 is doublet radical
        0..=3 | 5..=7 => Ok(None), // Other valid charge codes are not radicals by default
        _ => Err(FormatError::InvalidMolFormat(format!(
            "Invalid code '{}' passed to parse_radical_code",
            code
        ))
        .into()),
    }
}

/// Process stereo parity code (not stereo, odd, even, either or unmarked).
pub(crate) fn process_stereo_parity_code(code: u8) -> Result<Option<AtomStereoParity>> {
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

/// Process hydrogen count code (non-standard: 0 in non-query atoms).
pub(crate) fn process_hydrogen_count_code(code: u8) -> Result<Option<u8>> {
    match code {
        0 => Ok(None), // non-standard: 0 in non-query atoms
        1..=5 => Ok(Some(code - 1)),
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
    fn test_process_atom_symbol() {
        assert_eq!(
            process_atom_symbol("H").unwrap(),
            AtomSymbol::Element(Element::H)
        );
        assert_eq!(
            process_atom_symbol("h").unwrap(),
            AtomSymbol::Element(Element::H)
        );
        assert_eq!(
            process_atom_symbol("Cu").unwrap(),
            AtomSymbol::Element(Element::Cu)
        );
        assert_eq!(
            process_atom_symbol("CU").unwrap(),
            AtomSymbol::Element(Element::Cu)
        );
        assert_eq!(
            process_atom_symbol("cu").unwrap(),
            AtomSymbol::Element(Element::Cu)
        );
        assert_eq!(process_atom_symbol("L").unwrap(), AtomSymbol::AtomList);
        assert_eq!(
            process_atom_symbol("A").unwrap(),
            AtomSymbol::Unspecified('A')
        );
        assert_eq!(
            process_atom_symbol("Q").unwrap(),
            AtomSymbol::Unspecified('Q')
        );

        assert_eq!(
            process_atom_symbol("D").unwrap(),
            AtomSymbol::Element(Element::H)
        );
        assert_eq!(
            process_atom_symbol("T").unwrap(),
            AtomSymbol::Element(Element::H)
        );
        assert_eq!(
            process_atom_symbol("*").unwrap(),
            AtomSymbol::Unspecified('*')
        );
        assert_eq!(process_atom_symbol("LP").unwrap(), AtomSymbol::LonePair);
        assert_eq!(process_atom_symbol("R1").unwrap(), AtomSymbol::RGroup(0));
        assert_eq!(process_atom_symbol("R10").unwrap(), AtomSymbol::RGroup(9));
        assert!(process_atom_symbol("X").is_err());
        assert!(process_atom_symbol("R0").is_err());
    }

    #[test]
    fn test_process_isotope_symbol() {
        assert_eq!(process_isotope_symbol("D").unwrap(), Some(2));
        assert_eq!(process_isotope_symbol("d").unwrap(), Some(2));
        assert_eq!(process_isotope_symbol("H").unwrap(), None);
        assert!(process_isotope_symbol("X").is_err());
    }

    #[test]
    fn test_process_mass_diff_code() {
        assert_eq!(process_mass_diff_code(0).unwrap(), None);
        assert_eq!(process_mass_diff_code(-3).unwrap(), Some(-3));
        assert_eq!(process_mass_diff_code(4).unwrap(), Some(4));
        assert!(process_mass_diff_code(5).is_err());
    }

    #[test]
    fn test_process_valence_code() {
        assert_eq!(process_valence_code(0).unwrap(), None);
        assert_eq!(process_valence_code(1).unwrap(), Some(1));
        assert_eq!(process_valence_code(14).unwrap(), Some(14));
        assert_eq!(process_valence_code(15).unwrap(), Some(0));
        assert!(process_valence_code(16).is_err());
    }

    #[test]
    fn test_process_charge_code() {
        assert_eq!(process_charge_code(0).unwrap(), 0);
        assert_eq!(process_charge_code(1).unwrap(), 3);
        assert_eq!(process_charge_code(4).unwrap(), 0);
        assert_eq!(process_charge_code(5).unwrap(), -1);
        assert!(process_charge_code(8).is_err());
    }

    #[test]
    fn test_process_radical_code() {
        assert_eq!(process_radical_code(4).unwrap(), Some(2));
        assert_eq!(process_radical_code(1).unwrap(), None);
        assert!(process_radical_code(8).is_err());
    }

    #[test]
    fn test_process_stereo_parity_code() {
        assert_eq!(process_stereo_parity_code(0).unwrap(), None);
        assert_eq!(
            process_stereo_parity_code(1).unwrap(),
            Some(AtomStereoParity::Odd)
        );
        assert_eq!(
            process_stereo_parity_code(2).unwrap(),
            Some(AtomStereoParity::Even)
        );
        assert_eq!(
            process_stereo_parity_code(3).unwrap(),
            Some(AtomStereoParity::Either)
        );
    }

    #[test]
    fn test_process_hydrogen_count_code() {
        assert_eq!(process_hydrogen_count_code(0).unwrap(), None);
        assert_eq!(process_hydrogen_count_code(1).unwrap(), Some(0));
        assert!(process_hydrogen_count_code(6).is_err());
    }
}
