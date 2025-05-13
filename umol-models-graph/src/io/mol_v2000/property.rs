//! Property parsing functions for MOL files.

use super::atom::parse_index;
use crate::{Atom, Bond, SGroup, SGroupType};
use fixed_width::{from_bytes_with_fields, FieldSet, Reader};
use serde::de::DeserializeOwned;
use umol::error::FormatError;
use umol::{Error, Result};

/// Type for M line parsers for atom properties
pub(crate) type MAtomParserFn = fn(&mut Vec<Atom>, &[u8]) -> Result<()>;
/// Type for M line SGroup property parsers
pub(crate) type MSGroupParserFn = fn(
    sgroups: &mut Vec<SGroup>,
    value: &[u8],
    atoms: &[Atom],
    bonds: &[(usize, usize, Bond)],
) -> Result<()>;

/// Format: `M  XXXnn8 aaa vvv ...`
/// `count` is the number of fields per entry (e.g., 2 for `aaa vvv`).
pub(crate) fn parse_m_pairs<'de, T: DeserializeOwned>(input: &'de [u8]) -> Result<Vec<(usize, T)>> {
    let count: usize = from_bytes_with_fields(
        &input[6..9],
        FieldSet::Seq(vec![FieldSet::new_field(0..3).name("count")]),
    )
    .map_err(|_| {
        Error::from(FormatError::InvalidMolFormat(
            "Invalid count in M line".to_string(),
        ))
    })?;

    if count > 8 {
        return Err(Error::from(FormatError::InvalidMolFormat(
            "More than 8 pairs in M line".to_string(),
        )));
    }

    let pair_bytes = Reader::from_bytes(&input[9..])
        .width(8)
        .byte_reader()
        .map(|line| {
            line.map_err(|e| {
                Error::from(FormatError::InvalidMolFormat(format!(
                    "Invalid pair format in M line: {:?}",
                    e
                )))
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let mut pairs = Vec::with_capacity(count);
    for pair in pair_bytes {
        let (index_str, value) = from_bytes_with_fields::<(String, T)>(
            &pair,
            FieldSet::Seq(vec![
                FieldSet::new_field(1..4).name("index"),
                FieldSet::new_field(5..8).name("value"),
            ]),
        )
        .map_err(|e| {
            Error::from(FormatError::InvalidMolFormat(format!(
                "Invalid pair format in M line: {:?}",
                e
            )))
        })?;

        let index = parse_index(index_str.as_bytes())?;
        pairs.push((index, value));
    }

    if pairs.len() != count {
        return Err(Error::from(FormatError::InvalidMolFormat(format!(
            "Invalid number of pairs in M line: {} (expected {})",
            pairs.len(),
            count
        ))));
    }

    Ok(pairs)
}

/// Charge property parser
/// Format: `M  CHGnn8 aaa vvv ...`
/// TODO: Rewrite to return editing commands
pub(crate) fn parse_m_chg(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
    let pairs = parse_m_pairs::<i8>(input)?;

    for (idx, charge) in pairs {
        let atom = &mut atoms[idx];
        if charge >= -15 && charge <= 15 {
            atom.charge = charge as i8;
        } else {
            return Err(Error::from(FormatError::InvalidMolFormat(format!(
                "Charge value {} for atom {} out of range in M  CHG line.",
                charge, idx
            ))));
        }
    }
    Ok(())
}

/// Isotope property parser
/// Format: `M  ISOnn8 aaa vvv ...`
/// TODO: Rewrite to return editing commands
pub(crate) fn parse_m_iso(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
    let pairs = parse_m_pairs::<i32>(input)?;

    for (idx, mass_diff) in pairs {
        let atom = &mut atoms[idx];
        if mass_diff == 0 {
            atom.isotope_mass = None;
        } else if mass_diff >= -18 && mass_diff <= 12 {
            atom.isotope_mass = Some(
                (atom.element.reference_atomic_mass() as i32 + mass_diff)
                    .try_into()
                    .unwrap(),
            );
        } else {
            return Err(Error::from(FormatError::InvalidMolFormat(format!(
                "Mass difference {} for atom {} out of range in M  ISO line.",
                mass_diff, idx
            ))));
        }
    }
    Ok(())
}

/// Radical property parser
/// Format: `M  RADnn8 aaa vvv ...`
/// TODO: Rewrite to return editing commands
pub(crate) fn parse_m_rad(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
    let pairs = parse_m_pairs::<i64>(input)?;

    for (idx, rad) in pairs {
        let atom = &mut atoms[idx];
        match rad {
            1 => atom.radical = Some(1), // Singlet
            2 => atom.radical = Some(2), // Doublet
            3 => atom.radical = Some(3), // Triplet
            0 => atom.radical = None,    // Explicitly non-radical
            _ => {
                return Err(Error::from(FormatError::InvalidMolFormat(format!(
                    "Invalid radical value {} for atom {} in M RAD line, ignoring.",
                    rad, idx
                ))));
            }
        }
    }
    Ok(())
}

/// SGroup type parser
/// Format: M  STYnn8 sss ttt ...
pub(crate) fn parse_m_sty(
    sgroups: &mut Vec<SGroup>,
    input: &[u8],
    _atoms: &[Atom],
    _bonds: &[(usize, usize, Bond)],
) -> Result<()> {
    let pairs = parse_m_pairs::<String>(input)?;

    for (index, group_type) in pairs {
        let sgroup = SGroup {
            id: index,
            group_type: match group_type.as_str() {
                "GEN" => SGroupType::Generic,
                "MUL" => SGroupType::MultipleGroup,
                "SRU" => SGroupType::SRU,
                "SUP" => SGroupType::Superatom,
                "DAT" => SGroupType::Data,
                _ => SGroupType::Unknown(group_type),
            },
            label: None,
            subscript: None,
            atom_indices: Vec::new(),
            bond_indices: Vec::new(),
        };

        if index < sgroups.len() {
            sgroups[index] = sgroup;
        } else if index == sgroups.len() {
            sgroups.push(sgroup);
        } else {
            return Err(Error::from(FormatError::InvalidMolFormat(format!(
                "SGroup ID {} out of range in M STY line (expected up to {})",
                index,
                sgroups.len()
            ))));
        }
    }
    Ok(())
}

/// SGroup atom list parser
/// Format: `M  SAL sssn15 aaa ...`
pub(crate) fn parse_m_sal(
    sgroups: &mut Vec<SGroup>,
    input: &[u8],
    _atoms: &[Atom],
    _bonds: &[(usize, usize, Bond)],
) -> Result<()> {
    let (group, count): (usize, usize) = from_bytes_with_fields(
        &input[7..13],
        FieldSet::Seq(vec![
            FieldSet::new_field(0..3).name("index"),
            FieldSet::new_field(3..6).name("count"),
        ]),
    )
    .map_err(|_| {
        Error::from(FormatError::InvalidMolFormat(
            "Invalid index and count in M SAL line".to_string(),
        ))
    })?;

    if group >= sgroups.len() {
        return Err(FormatError::InvalidMolFormat(format!(
            "SGroup ID {} out of range in M SAL line",
            group
        ))
        .into());
    }

    let indices_bytes = Reader::from_bytes(&input[13..])
        .width(4)
        .byte_reader()
        .map(|line| {
            line.map_err(|e| {
                Error::from(FormatError::InvalidMolFormat(format!(
                    "Invalid index format in M SAL line: {:?}",
                    e
                )))
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let mut indices = Vec::with_capacity(count);
    for index in indices_bytes {
        let index_str = from_bytes_with_fields::<String>(
            &index,
            FieldSet::Seq(vec![FieldSet::new_field(1..4).name("index")]),
        )
        .map_err(|e| {
            Error::from(FormatError::InvalidMolFormat(format!(
                "Invalid index format in M SAL line: {:?}",
                e
            )))
        })?;
        let index = parse_index(index_str.as_bytes())?;
        indices.push(index);
    }

    let sgroup = &mut sgroups[group];
    sgroup.atom_indices.extend(indices);

    Ok(())
}

/// SGroup bond list parser
/// Format: `M  SBL sssn15 bbb ...`
pub(crate) fn parse_m_sbl(
    sgroups: &mut Vec<SGroup>,
    input: &[u8],
    _atoms: &[Atom],
    _bonds: &[(usize, usize, Bond)],
) -> Result<()> {
    let (group, count): (usize, usize) = from_bytes_with_fields(
        &input[7..13],
        FieldSet::Seq(vec![
            FieldSet::new_field(0..3).name("index"),
            FieldSet::new_field(3..6).name("count"),
        ]),
    )
    .map_err(|_| {
        Error::from(FormatError::InvalidMolFormat(
            "Invalid index and count in M SBL line".to_string(),
        ))
    })?;

    let indices_bytes = Reader::from_bytes(&input[13..])
        .width(4)
        .byte_reader()
        .map(|line| {
            line.map_err(|e| {
                Error::from(FormatError::InvalidMolFormat(format!(
                    "Invalid index format in M SBL line: {:?}",
                    e
                )))
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let mut indices = Vec::with_capacity(count);
    for index in indices_bytes {
        let index_str = from_bytes_with_fields::<String>(
            &index,
            FieldSet::Seq(vec![FieldSet::new_field(1..4).name("index")]),
        )
        .map_err(|e| {
            Error::from(FormatError::InvalidMolFormat(format!(
                "Invalid index format in M SBL line: {:?}",
                e
            )))
        })?;
        let index = parse_index(index_str.as_bytes())?;
        indices.push(index);
    }

    let sgroup = &mut sgroups[group];
    sgroup.bond_indices.extend(indices);

    Ok(())
}

/// SGroup label parser (M SLB)
/// Format: `M  SLB sss vvv ...`
pub(crate) fn parse_m_slb(
    sgroups: &mut Vec<SGroup>,
    input: &[u8],
    _atoms: &[Atom],
    _bonds: &[(usize, usize, Bond)],
) -> Result<()> {
    let pairs = parse_m_pairs::<String>(input)?;
    for (index, label) in pairs {
        let sgroup = &mut sgroups[index];
        sgroup.label = Some(label);
    }

    Ok(())
}

/// Atom alias parser
/// 
pub(crate) fn parse_a_prop(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
    let (index, alias) = from_bytes_with_fields::<(usize, String)>(
        &input[3..],
        FieldSet::Seq(vec![
            FieldSet::new_field(0..3).name("index"),
            FieldSet::new_field(4..80).name("alias"),
        ]),
    )
    .map_err(|_| {
        Error::from(FormatError::InvalidMolFormat(
            "Invalid index and alias in A line".to_string(),
        ))
    })?;

    let atom = &mut atoms[index];
    atom.properties.insert("molFileAlias".to_string(), alias);
    Ok(())
}

// /// Atom value parser
// pub(crate) fn parse_v_prop(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
//     let (index, value) = from_bytes_with_fields::<(usize, String)>(
//         &input[3..],
//         FieldSet::Seq(vec![
//             FieldSet::new_field(0..3).name("index"),
//             FieldSet::new_field(4..80).name("value"),
//         ]),
//     )
//     .map_err(|_| {
//         Error::from(FormatError::InvalidMolFormat(
//             "Invalid index and value in V line".to_string(),
//         ))
//     })?;

//     let atom = &mut atoms[index];
//     atom.properties.insert("molFileValue".to_string(), value);
//     Ok(())
// }

// /// Group abbreviation parser
// pub(crate) fn parse_g_prop(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
//     let (index, abbr) = from_bytes_with_fields::<(usize, String)>(
//         &input[3..],
//         FieldSet::Seq(vec![
//             FieldSet::new_field(0..3).name("index"),
//             FieldSet::new_field(4..80).name("abbreviation"),
//         ]),
//     )
//     Ok(())
// }

#[cfg(test)]
mod tests {
    use umol_data::Element;

    use super::*;

    #[test]
    fn test_parse_pairs() {
        let line = b"M  CHG  2   1   1   2  -1";
        let pairs = parse_m_pairs::<i64>(line).unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0, 0);
        assert_eq!(pairs[0].1, 1);
        assert_eq!(pairs[1].0, 1);
        assert_eq!(pairs[1].1, -1);
    }

    #[test]
    fn test_parse_chg() {
        let line = b"M  CHG  2   1   1   2  -1";
        let mut atoms = vec![Atom::new(Element::O), Atom::new(Element::C)];
        parse_m_chg(&mut atoms, line).unwrap();
        assert_eq!(atoms[0].charge, 1);
        assert_eq!(atoms[1].charge, -1);
    }

    #[test]
    fn test_parse_iso() {
        let line = b"M  ISO  1   1   2";
        let mut atoms = vec![Atom::new(Element::O), Atom::new(Element::C)];
        parse_m_iso(&mut atoms, line).unwrap();
        assert_eq!(atoms[0].isotope_mass, Some(18));
        assert_eq!(atoms[1].isotope_mass, None);
    }

    #[test]
    fn test_parse_rad() {
        let line = b"M  RAD  2   1   1   2   0";
        let mut atoms = vec![Atom::new(Element::H), Atom::new(Element::C)];
        parse_m_rad(&mut atoms, line).unwrap();
        assert_eq!(atoms[0].radical, Some(1));
        assert_eq!(atoms[1].radical, None);
    }

    #[test]
    fn test_parse_sty() {
        let line = b"M  STY  2   1 SUP   2 DAT";
        let mut sgroups = vec![];
        parse_m_sty(&mut sgroups, line, &[], &[]).unwrap();
        assert_eq!(sgroups[0].group_type, SGroupType::Superatom);
        assert_eq!(sgroups[1].group_type, SGroupType::Data);
    }

    // #[test]
    // fn test_parse_sal() {
    //     let line = b"M   SAL  2   1   1   2";
    //     let mut sgroups = vec![];
    //     parse_m_sal(&mut sgroups, line, &[], &[]).unwrap();
    //     assert_eq!(sgroups[0].atom_indices, vec![0, 1]);
    // }
}
