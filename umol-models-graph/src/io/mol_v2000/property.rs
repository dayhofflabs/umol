//! Property parsing functions for MOL files.

//! Implementation status for the Property Block MOL v2000 file
//!   https://en.wikipedia.org/wiki/Chemical_table_file
//!
//! | Property   | Symbol | Implementation | Class   | RDKit* | ChemAxon** | CDK*** | Notes                 |
//! |------------|--------|----------------|---------|--------|------------|--------|-----------------------|
//! | Atom Alias | A      | x              | ISIS    | x      | x          | x      |                       |
//! | Atom Value | V      | x              | ISIS    | x      | x          | x      |                       |
//! ! Group Abbr | G      | -              | ISIS    | -      | -          | -      | Outdated, use `M SUP` |
//! | Charge     | CHG    | x              | Generic | x      | x          | x      |                       |
//! | Radical    | RAD    | x              | Generic | x      | x          | x      |                       |
//! | Isotope    | ISO    | x              | Generic | x      | x          | x      |                       |
//! | Ring Bonds | RBC    | -              | Query   | x      | x          | x      |                       |
//! | Subs Count | SUB    | -              | Query   | x      | x          | x      |                       |
//! | Unsat Atom | UNS    | -              | Query   | x      | x          | x      |                       |
//! | Link Atom  | LIN    | -              | Query   | x      | x          | x      |                       |
//! | Atom List  | ALS    | -              | Query   | x      | x          | x      |                       |
//! | Att Point  | APO    | -              | RGroup  | x      | x          | x      |                       |
//! | Att Order  | AAL    | -              | RGroup  | -      | -          | x      |                       |
//! | Lab Loc    | RGP    | -              | RGroup  | x      | x          | x      |                       |
//! | Logic      | LOG    | -              | RGroup  | -      | x          | x      |                       |
//! | Sgrp Type  | STY    | x              | SGroup  | x      | x          | x      |                       |
//! | Sgrp Subt  | SST    | -              | SGroup  | x      | x          | x      |                       |
//! | Sgrp Label | SLB    | x              | SGroup  | x      | -          | x      |                       |
//! | Sgrp Conn  | SCN    | -              | SGroup  | x      | x          | x      |                       |
//! | Sgrp Expan | SDS    | -              | SGroup  | x      | x          | x      |                       |
//! | Sgrp Atoms | SAL    | x              | SGroup  | x      | x          | x      |                       |
//! | Sgrp Bonds | SBL    | x              | SGroup  | x      | x          | x      |                       |
//! | Sgrp Parnt | SPA    | -              | SGroup  | x      | x          | x      |                       |
//! | Sgrp Subs  | SMT    | -              | SGroup  | x      | x          | x      |                       |
//! | Sgrp Corr  | CRS    | -              | SGroup  | x      | -          | x      |                       |
//! | Sgrp Disp  | SDI    | -              | SGroup  | x      | x          | x      |                       |
//! | Sup Bd Vec | SBV    | -              | SGroup  | x      | -          | x      |                       |
//! | Data Flds  | SDT    | -              | SGroup  | x      | x          | x      |                       |
//! | Data Disp  | SDD    | -              | SGroup  | x      | x          | x      |                       |
//! | Data Sgrp  | SCD    | -              | SGroup  | x      | x          | x      | Continued data line   |
//! | Data Sgrp  | SED    | -              | SGroup  | x      | x          | x      | End of data line      |
//! | Sgrp Hier  | SPL    | -              | SGroup  | x      | x          | x      | Parent list           |
//! | Sgrp Comp# | SNC    | -              | SGroup  | x      | x          | x      |                       |
//! | 3D Feat    | $3D    | -              | 3D      | -      | -          | x      |                       |
//! | Phantom    | PXA    | -              | ISIS    | x      | -          | -      |                       |
//! | Sup Att Pt | SAP    | -              | ISIS    | x      | x          | -      |                       |
//! | Sup Class  | SCL    | -              | ISIS    | x      | -          | -      |                       |
//! | Regno      | REG    | -              | ISIS    | -      | -          | -      |                       |
//! | Sgrp Brkt  | SBT    | -              | ISIS    | x      | x          | x      |                       |
//! | 0-Order Bd | ZBO    | -              | Bd Ext  | x      | -          | -      | DOI:10.1021/ci200488k |
//! | Virt Hs    | ZCH    | -              | Bd Ext  | x      | -          | -      | DOI:10.1021/ci200488k |                     |
//! | Marvin SM  | MRV    | -              | Marvin  | x      | x          | -      |                       |
//! | Atom Label | ZZC    | -              | ADC     | -      | -          | x      |                       |
//! | Skip       | SKIP   | -              | Generic | ?      | -          | x      |                       |
//! | End        | END    | x              | Generic | x      | x          | x      |                       |
//! 
//! * RDKit: https://www.rdkit.org/docs/GettingStartedInPython.html#writing-molecules
//! ** ChemAxon: https://docs.chemaxon.com/display/docs/formats_mdl-molfiles-rgfiles-sdfiles-rxnfiles-rdfiles-formats.md
//! *** CDK: https://cdk.github.io/cdk/latest/docs/api/org/openscience/cdk/io/MDLV2000Reader.html

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
pub(crate) fn parse_m_chg(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
    let pairs = parse_m_pairs::<i8>(input)?;

    for (idx, charge) in pairs {
        let atom = &mut atoms[idx];
        println!("DEBUG CHG PROP: INDEX: {}, VALUE: {}", idx, charge);
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
pub(crate) fn parse_m_iso(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
    let pairs = parse_m_pairs::<i32>(input)?;

    for (idx, mass_diff) in pairs {
        let atom = &mut atoms[idx];
        println!("DEBUG ISO PROP: INDEX: {}, VALUE: {}", idx, mass_diff);
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
pub(crate) fn parse_m_rad(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
    let pairs = parse_m_pairs::<i64>(input)?;

    for (idx, rad) in pairs {
        let atom = &mut atoms[idx];
        println!("DEBUG RAD PROP: INDEX: {}, VALUE: {}", idx, rad);
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
                _ => SGroupType::Unknown(group_type.clone()),
            },
            label: None,
            subscript: None,
            atom_indices: Vec::new(),
            bond_indices: Vec::new(),
        };

        println!("DEBUG STY PROP: INDEX: {}, VALUE: {}", index, group_type);
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
    let (group_str, count): (String, usize) = from_bytes_with_fields(
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

    let group = parse_index(group_str.as_bytes())?;
    if group > sgroups.len() {
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

    println!("DEBUG SAL PROP: GROUP: {}, INDICES: {:?}", group, indices);
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
    let (group_str, count): (String, usize) = from_bytes_with_fields(
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

    let group = parse_index(group_str.as_bytes())?;

    if group > sgroups.len() {
        return Err(FormatError::InvalidMolFormat(format!(
            "SGroup ID {} out of range in M SBL line",
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

    println!("DEBUG SBL PROP: GROUP: {}, INDICES: {:?}", group, indices);

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
        println!("DEBUG SLB PROP: INDEX: {}, VALUE: {}", index, label);
        sgroup.label = Some(label);
    }
    Ok(())
}

/// Atom alias parser
/// Format: `A  aaa\nx...` (constructed using `combine_next_n` from 2 lines)
pub(crate) fn parse_a_prop(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
    let (index_str, alias) = from_bytes_with_fields::<(String, String)>(
        &input[3..],
        FieldSet::Seq(vec![
            FieldSet::new_field(0..3).name("index"),
            FieldSet::new_field(4..input.len() - 3).name("alias"),
        ]),
    )
    .map_err(|_| {
        Error::from(FormatError::InvalidMolFormat(
            "Invalid index and alias in A line".to_string(),
        ))
    })?;

    let index = parse_index(index_str.as_bytes())?;
    if index >= atoms.len() {
        return Err(FormatError::InvalidMolFormat(format!(
            "Atom index {} out of range in A line",
            index
        ))
        .into());
    }

    println!("DEBUG A PROP: INDEX: {}, VALUE: {}", index, alias);

    let atom = &mut atoms[index];
    atom.properties.insert("molFileAlias".to_string(), alias);
    Ok(())
}

/// Atom value parser
pub(crate) fn parse_v_prop(atoms: &mut Vec<Atom>, input: &[u8]) -> Result<()> {
    let (index_str, value) = from_bytes_with_fields::<(String, String)>(
        &input[3..],
        FieldSet::Seq(vec![
            FieldSet::new_field(0..3).name("index"),
            FieldSet::new_field(4..input.len() - 3).name("value"),
        ]),
    )
    .map_err(|_| {
        Error::from(FormatError::InvalidMolFormat(
            "Invalid index and value in V line".to_string(),
        ))
    })?;

    let index = parse_index(index_str.as_bytes())?;
    if index >= atoms.len() {
        return Err(FormatError::InvalidMolFormat(format!(
            "Atom index {} out of range in V line",
            index
        ))
        .into());
    }

    println!("DEBUG V PROP: INDEX: {}, VALUE: {}", index, value);

    let atom = &mut atoms[index];
    atom.properties.insert("molFileValue".to_string(), value);

    Ok(())
}

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
        let line = b"M  STY  3   1 SUP   2 SUP   3 SUP";
        let mut sgroups = vec![];
        parse_m_sty(&mut sgroups, line, &[], &[]).unwrap();
        assert_eq!(sgroups[0].group_type, SGroupType::Superatom);
        assert_eq!(sgroups[1].group_type, SGroupType::Superatom);
        assert_eq!(sgroups[2].group_type, SGroupType::Superatom);
    }

    #[test]
    fn test_parse_sal() {
        let line = b"M  SAL   1  1   1";
        let mut sgroups = vec![SGroup {
            id: 0,
            group_type: SGroupType::Superatom,
            label: None,
            subscript: None,
            atom_indices: Vec::new(),
            bond_indices: Vec::new(),
        }];
        parse_m_sal(&mut sgroups, line, &[], &[]).unwrap();
        assert_eq!(sgroups[0].atom_indices, vec![0]);
    }

    #[test]
    fn test_parse_sbl() {
        let line = b"M  SBL   1  1   2";
        let mut sgroups = vec![SGroup {
            id: 0,
            group_type: SGroupType::Superatom,
            label: None,
            subscript: None,
            atom_indices: Vec::new(),
            bond_indices: Vec::new(),
        }];
        parse_m_sbl(&mut sgroups, line, &[], &[]).unwrap();
        assert_eq!(sgroups[0].bond_indices, vec![1]);
    }

    #[test]
    fn test_parse_slb() {
        let line = b"M  SLB  1   1   1";
        let mut sgroups = vec![SGroup {
            id: 0,
            group_type: SGroupType::Superatom,
            label: None,
            subscript: None,
            atom_indices: Vec::new(),
            bond_indices: Vec::new(),
        }];
        parse_m_slb(&mut sgroups, line, &[], &[]).unwrap();
        assert!(matches!(&sgroups[0].label, Some(label) if label == "1"));
    }

    #[test]
    fn test_parse_a_prop() {
        let line = b"A    1 CF3";
        let mut atoms = vec![
            Atom::new(Element::C),
            Atom::new(Element::F),
            Atom::new(Element::F),
            Atom::new(Element::F),
        ];
        parse_a_prop(&mut atoms, line).unwrap();
        assert_eq!(atoms[0].properties["molFileAlias"], "CF3");
    }

    #[test]
    fn test_parse_v_prop() {
        let line = b"V    1 *";
        let mut atoms = vec![Atom::new(Element::C)];
        parse_v_prop(&mut atoms, line).unwrap();
        assert_eq!(atoms[0].properties["molFileValue"], "*");
    }

}
