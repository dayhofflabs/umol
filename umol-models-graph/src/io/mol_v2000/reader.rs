//! MOL v2000 file reader.

use super::atom::{
    parse_atom_symbol, parse_charge_code, parse_hydrogen_count_code, parse_radical_code,
    parse_stereo_parity_code, parse_valence_code, AtomSymbol,
};
use super::bond::{parse_bond_dir_code, parse_bond_stereo_code, parse_bond_type_code};
use super::conformer::is_3d;
use super::property::{
    parse_a_prop, parse_m_chg, parse_m_iso, parse_m_rad, parse_m_sal, parse_m_sbl, parse_m_slb,
    parse_m_smt, parse_m_sty, MAtomParserFn, MSGroupParserFn,
};
use crate::io::utils::{detect_line_break, CombineNextN};
use crate::{Atom, Bond, BondType, Conformer, Molecule, Point3D, SGroup};
use fixed_width::Reader;
use fixed_width::{FieldSet, FixedWidth, LineBreak};
use serde::Deserialize;
use std::collections::HashMap;
use std::io::BufRead;
use std::result::Result as StdResult;
use umol::error::FormatError;
use umol::{Error, Result};

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
pub(crate) struct CountsLine {
    atoms: i32,      // 'aaa' - number of atoms (max 255)
    bonds: i32,      // 'bbb' - number of bonds (max 255)
    atom_lists: i32, // 'lll' - number of atom lists (max 30)
    // fff is obsolete, skipping
    chiral_flag: i32,   // 'ccc' - chiral flag (0=not chiral, 1=chiral)
    stext_entries: i32, // 'sss' - number of stext entries
    // xxx is obsolete, skipping
    // rrr is obsolete, skipping
    // ppp is obsolete, skipping
    // iii is obsolete, skipping
    properties_lines: i32, // 'mmm' - number of additional properties lines
    version: String,       // 'vvvvvv' - version stamp (V2000)
}

impl FixedWidth for CountsLine {
    fn fields() -> FieldSet {
        FieldSet::Seq(vec![
            FieldSet::new_field(0..3).name("atoms"),
            FieldSet::new_field(3..6).name("bonds"),
            FieldSet::new_field(6..9).name("atom_lists"),
            // Skip obsolete field fff (9..12)
            FieldSet::new_field(12..15).name("chiral_flag"),
            FieldSet::new_field(15..18).name("stext_entries"),
            // Skip obsolete field xxx (18..21)
            // Skip obsolete field rrr (21..24)
            // Skip obsolete field ppp (24..27)
            // Skip obsolete field iii (27..30)
            FieldSet::new_field(30..33).name("properties_lines"),
            FieldSet::new_field(33..39).name("version"),
        ])
    }
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
pub(crate) struct AtomLine {
    x: f64,
    y: f64,
    z: f64,
    symbol: String,     // 'aaa' field - atom symbol (see AtomSymbol enum)
    mass_diff: i8, // 'dd' field - mass difference (-3, -2, -1, 0, 1, 2, 3, 4), 0 if value outside of this range
    charge: u8, // 'ccc' field - 0 = uncharged, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
    stereo_parity: u8, // 'sss' field - 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked
    hydrogen_count: u8, // 'hhh' field - 1 = H0, 2 = H1, 3 = H2, 4 = H3, 5 = H4
    stereo_care: u8, // 'bbb' field - 0 = ignore stereo, 1 = stereo in query must match
    valence: u8, // 'vvv' field - 0 = default, 1-14 = explicit, 15 = explicit 0
    // Skipping obsolete or unused fields (HHH, rrr, iii)
    atom_mapping: u8, // 'mmm' field - 1 = number of atoms
    inversion: u8,    // 'nnn' field - 0 = property not applied, 2 = inverted, 3 = retained
    exact_change: u8, // 'eee' field - 0 = property not applied, 1 = charge in query must match
}

impl FixedWidth for AtomLine {
    fn fields() -> FieldSet {
        FieldSet::Seq(vec![
            FieldSet::new_field(0..10).name("x"),
            FieldSet::new_field(10..20).name("y"),
            FieldSet::new_field(20..30).name("z"),
            // Position 30 is a space
            FieldSet::new_field(31..34).name("symbol"),
            FieldSet::new_field(34..36).name("mass_diff"),
            FieldSet::new_field(36..39).name("charge"),
            FieldSet::new_field(39..42).name("stereo_parity"),
            FieldSet::new_field(42..45).name("hydrogen_count"),
            FieldSet::new_field(45..48).name("stereo_care"),
            FieldSet::new_field(48..51).name("valence"),
            // Skipping fields HHH, rrr, iii (51..60)
            FieldSet::new_field(60..63).name("atom_mapping"),
            FieldSet::new_field(63..66).name("inversion"),
            FieldSet::new_field(66..69).name("exact_change"),
        ])
    }
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
pub(crate) struct BondLine {
    atom1: usize,    // '111' - atom 1 number
    atom2: usize,    // '222' - atom 2 number
    bond_type: u8,   // 'ttt' - bond type (1=Single, 2=Double, etc.)
    bond_stereo: u8, // 'sss' - bond stereo (0=Not stereo, 1=Up, etc.)
    // xxx is not used, skipping
    bond_topology: u8,   // 'rrr' - bond topology (0=Either, 1=Ring, 2=Chain)
    reacting_center: u8, // 'ccc' - reacting center status
}

impl FixedWidth for BondLine {
    fn fields() -> FieldSet {
        FieldSet::Seq(vec![
            FieldSet::new_field(0..3).name("atom1"),
            FieldSet::new_field(3..6).name("atom2"),
            FieldSet::new_field(6..9).name("bond_type"),
            FieldSet::new_field(9..12).name("bond_stereo"),
            // Skip unused field xxx (12..15)
            FieldSet::new_field(15..18).name("bond_topology"),
            FieldSet::new_field(18..21).name("reacting_center"),
        ])
    }
}

// Parse V2000 MOL file.
pub fn read_mol_v2000(mut reader: impl BufRead) -> Result<Molecule> {
    let mut molecule = Molecule::new();
    let line_break: LineBreak = detect_line_break(&mut reader)
        .map_err(|e| Error::from(FormatError::IoError(e)))?
        .try_into()
        .map_err(|e| {
            Error::from(FormatError::InvalidMolFormat(format!(
                "Unsupported line break: {}",
                e
            )))
        })?;

    // Parse header block: molecule name, info line, comment line
    let mut name_line = String::new();
    reader
        .read_line(&mut name_line)
        .map_err(|e| Error::from(FormatError::IoError(e)))?;
    if !name_line.trim().is_empty() {
        molecule.set_prop("mol_name".to_string(), name_line.trim_end().to_string());
        // Trim newline
    }
    let mut _line_buffer = String::new();
    // Read and discard info line and comment line
    for _ in 0..2 {
        reader
            .read_line(&mut _line_buffer)
            .map_err(|e| Error::from(FormatError::IoError(e)))?;
        _line_buffer.clear();
    }

    println!("DEBUG: BEFORE COUNTS LINE");
    println!(
        "DEBUG: BUFFER: {:?}",
        String::from_utf8_lossy(reader.fill_buf().unwrap())
    );

    // Counts line
    const COUNTS_LINE_WIDTH: usize = 39;
    let mut counts_buffer = String::new();
    let length = reader
        .read_line(&mut counts_buffer)
        .map_err(|e| {
            Err(FormatError::InvalidMolFormat(format!(
                "Failed to read counts line: {}",
                e
            )).into())
        })?;

    if length < COUNTS_LINE_WIDTH {
        return Err(FormatError::InvalidMolFormat(format!(
                    "Counts line too short: found {}, expected {}",
            length, COUNTS_LINE_WIDTH
        )));
    }

    let counts_data = Reader::from_string(counts_buffer)
        .width(COUNTS_LINE_WIDTH)
        .linebreak(line_break.clone())
        .byte_reader()
        .filter_map(StdResult::ok)
        .next()
        .unwrap();

    println!(
        "DEBUG: COUNTS DATA: {:?}",
        String::from_utf8_lossy(&counts_data)
    );
    println!(
        "DEBUG: BUFFER: {:?}",
        String::from_utf8_lossy(reader.fill_buf().unwrap())
    );

    let counts_data = fixed_width::from_bytes::<CountsLine>(&counts_data).map_err(|e| {
        FormatError::InvalidMolFormat(format!("Failed to parse counts line: {}", e))
    })?;
    if counts_data.version != "V2000" {
        return Err(
            FormatError::InvalidMolFormat("Only V2000 format supported".to_string()).into(),
        );
    }

    println!("DEBUG: AFTER COUNTS LINE");
    println!(
        "DEBUG: BUFFER: {:?}",
        String::from_utf8_lossy(reader.fill_buf().unwrap())
    );

    let num_atoms: usize = counts_data.atoms as usize;
    let num_bonds: usize = counts_data.bonds as usize;

    println!("DEBUG: NUM ATOMS: {:?}", num_atoms);
    println!("DEBUG: NUM BONDS: {:?}", num_bonds);
    println!("DEBUG: BEFORE ATOM BLOCK");
    println!(
        "DEBUG: BUFFER: {:?}",
        String::from_utf8_lossy(reader.fill_buf().unwrap())
    );

    // Atom block
    const ATOM_LINE_WIDTH: usize = 69;
    let mut atom_buffer = vec![0; num_atoms * (ATOM_LINE_WIDTH + line_break.byte_width())];
    reader.read_exact(&mut atom_buffer).map_err(|e| {
        Error::from(FormatError::InvalidMolFormat(format!(
            "Failed to read atom block: {}",
            e
        )))
    })?;

    println!(
        "DEBUG: ATOM BUFFER: {:?}, CAPACITY: {:?}",
        String::from_utf8_lossy(&atom_buffer),
        atom_buffer.capacity()
    );

    let mut atom_reader = Reader::from_bytes(atom_buffer)
        .width(ATOM_LINE_WIDTH)
        .linebreak(line_break.clone());

    let (mut atoms, positions) = atom_reader.byte_reader().try_fold(
        (Vec::new(), Vec::new()),
        |(mut atoms, mut positions), res| -> Result<_> {
            let bytes = res.map_err(|e| {
                Error::from(FormatError::InvalidMolFormat(format!(
                    "Failed to parse atom line: {}",
                    e
                )))
            })?;
            let atom_data = fixed_width::from_bytes::<AtomLine>(&bytes).map_err(|e| {
                Error::from(FormatError::InvalidMolFormat(format!(
                    "Failed to parse atom line: {}",
                    e
                )))
            })?;
            let atom_symbol = parse_atom_symbol(&atom_data.symbol)?;
            if let AtomSymbol::Element(element) = atom_symbol {
                let mut atom = Atom::new(element);
                atom.charge = parse_charge_code(atom_data.charge)?;
                atom.radical = parse_radical_code(atom_data.charge)?;
                atom.isotope_mass = if atom_data.mass_diff == 0 {
                    None
                } else {
                    Some(
                        (atom.element.reference_atomic_mass() as i64 + atom_data.mass_diff as i64)
                            .try_into()
                            .unwrap(),
                    )
                };
                atom.stereo_parity = parse_stereo_parity_code(atom_data.stereo_parity)?;
                atom.valence = parse_valence_code(atom_data.valence)?;
                atom.hydrogen_count = parse_hydrogen_count_code(atom_data.hydrogen_count)?;
                atom.atom_map_num = if atom_data.atom_mapping == 0 {
                    None
                } else {
                    Some(atom_data.atom_mapping as u32)
                };
                let point = Point3D::new(atom_data.x, atom_data.y, atom_data.z);
                positions.push(point);
                atoms.push(atom);
            } else {
                // TODO: Handle atom lists, lone pairs, R groups, etc.
                Err(Error::from(FormatError::InvalidMolFormat(format!(
                    "Invalid atom symbol: '{}'",
                    atom_data.symbol
                ))))?
            }
            Ok((atoms, positions))
        },
    )?;

    println!("DEBUG: AFTER ATOM BLOCK");
    println!("DEBUG: ATOMS: {:?}", atoms);
    println!("DEBUG: POSITIONS: {:?}", positions);
    println!("DEBUG: BEFORE BOND BLOCK");
    println!(
        "DEBUG: BUFFER: {:?}",
        String::from_utf8_lossy(reader.fill_buf().unwrap())
    );

    // Bond block
    const BOND_LINE_WIDTH: usize = 21;

    let mut bond_buffer = vec![0; num_bonds * (BOND_LINE_WIDTH + line_break.byte_width())];
    reader.read_exact(&mut bond_buffer).map_err(|e| {
        Error::from(FormatError::InvalidMolFormat(format!(
            "Failed to read bond block: {}",
            e
        )))
    })?;

    println!(
        "DEBUG: BOND BUFFER: {:?}, CAPACITY: {:?}",
        String::from_utf8_lossy(&bond_buffer),
        bond_buffer.capacity()
    );

    let mut bond_reader = Reader::from_bytes(bond_buffer)
        .width(BOND_LINE_WIDTH)
        .linebreak(line_break.clone());

    let bonds = bond_reader.byte_reader().try_fold(
        Vec::with_capacity(num_bonds),
        |mut bonds, res| -> Result<_> {
            let bytes = res.map_err(|e| {
                Error::from(FormatError::InvalidMolFormat(format!(
                    "Failed to parse bond line: {}",
                    e
                )))
            })?;
            let bond_data = fixed_width::from_bytes::<BondLine>(&bytes).map_err(|e| {
                Error::from(FormatError::InvalidMolFormat(format!(
                    "Failed to parse bond line: {}",
                    e
                )))
            })?;

            let idx1 = bond_data.atom1.checked_sub(1).ok_or_else(|| {
                Error::from(FormatError::InvalidMolFormat("Atom index 1 out of bounds".to_string()))
            })?;
            let idx2 = bond_data.atom2.checked_sub(1).ok_or_else(|| {
                Error::from(FormatError::InvalidMolFormat("Atom index 2 out of bounds".to_string()))
            })?;

            let bond_type = parse_bond_type_code(bond_data.bond_type)?;
            let mut bond = Bond::new(bond_type);

            match bond_type {
                BondType::Single => {
                    bond.dir = parse_bond_dir_code(bond_data.bond_stereo)?;
                }
                BondType::Double => {
                    bond.stereo = parse_bond_stereo_code(bond_data.bond_stereo)?;
                }
                _ => {}
            }
            // TODO: Process bond_data.bond_topology, bond_data.reacting_center if needed

            bonds.push((idx1, idx2, bond));
            Ok(bonds)
        },
    )?;

    println!("DEBUG: AFTER BOND BLOCK");
    println!("DEBUG: BONDS: {:?}", bonds);
    println!("DEBUG: AFTER BOND BLOCK");

    // Atom M property parsers
    let m_atom_parsers: HashMap<&'static [u8], MAtomParserFn> = [
        (b"CHG", parse_m_chg as MAtomParserFn),
        (b"ISO", parse_m_iso as MAtomParserFn),
        (b"RAD", parse_m_rad as MAtomParserFn),
        // TODO: Add other atom property parsers here
    ]
    .iter()
    .cloned()
    .collect();

    // SGroup M property parsers
    let m_sgroup_parsers: HashMap<&'static [u8], MSGroupParserFn> = [
        (b"STY", parse_m_sty as MSGroupParserFn),
        (b"SAL", parse_m_sal as MSGroupParserFn),
        (b"SBL", parse_m_sbl as MSGroupParserFn),
        (b"SLB", parse_m_slb as MSGroupParserFn),
        (b"SMT", parse_m_smt as MSGroupParserFn),
        // TODO: Add other SGroup property parsers here
    ]
    .iter()
    .cloned()
    .collect();

    // SGroup definitions
    let mut sgroups: Vec<SGroup> = Vec::new();

    println!("DEBUG: BEFORE PROPERTIES BLOCK");
    println!(
        "DEBUG: BUFFER: {:?}",
        String::from_utf8_lossy(reader.fill_buf().unwrap())
    );

    // Properties block
    let mut properties_reader = Reader::from_reader(reader).linebreak(line_break.clone());

    let terminated = properties_reader
        .byte_reader()
        .combine_next_n(
            |line| {
                if let Ok(line) = line {
                    if line.starts_with(b"A  ") || line.starts_with(b"G  ") {
                        Some(1u16)
                    } else {
                        None
                    }
                } else {
                    None
                }
            },
            " ",
        )
        .try_fold(false, |mut terminated, res| -> Result<bool> {
            let line = res.map_err(|e| {
                Error::from(FormatError::InvalidMolFormat(format!(
                    "Failed to read properties line: {}",
                    e
                )))
            })?;

            if line.starts_with(b"M  END") {
                // --- ADD Properties Block DEBUG ---
                println!("DEBUG Properties Block - Found M END");
                // --- END Properties Block DEBUG ---
                terminated = true;
            } else if terminated {
                return Err(Error::from(FormatError::InvalidMolFormat(
                    "Data found after M END".to_string(),
                )));
            } else if line.starts_with(b"M  ") {
                let key = &line[3..6];
                let value = &line[6..];
                if let Some(parser_fn) = m_atom_parsers.get(key) {
                    parser_fn(&mut atoms, value)?;
                } else if let Some(parser_fn) = m_sgroup_parsers.get(key) {
                    parser_fn(&mut sgroups, value, &atoms, &bonds)?;
                } else {
                    return Err(Error::from(FormatError::InvalidMolFormat(
                        "Malformed M line".to_string(),
                    )));
                }
            // } else if line.starts_with(b"V  ") {
            //     // TODO: Handle V property
            // } else if line.starts_with(b"A  ") {
            //     parse_a_prop(&mut atoms, &line[3..])?;
            // } else if line.starts_with(b"G  ") {
            //     // TODO: Handle G property
            // } else if !line.iter().all(|&b| b.is_ascii_whitespace()) {
            //     return Err(Error::from(FormatError::InvalidMolFormat(format!(
            //         "Invalid line in properties block: {:?}",
            //         line
            //     ))));
            // }
            Ok(terminated)
        })?;

    println!("DEBUG: AFTER PROPERTIES BLOCK");
    println!("DEBUG: TERMINATED: {:?}", terminated);

    if !terminated {
        return Err(Error::from(FormatError::InvalidMolFormat(
            "M END not found".to_string(),
        )));
    }

    // Add atoms to molecule
    for atom in atoms.into_iter() {
        molecule.add_atom(atom);
    }

    // Add bonds to molecule
    for (idx1, idx2, bond) in bonds {
        molecule.add_bond(idx1, idx2, bond);
    }

    // Add SGroups to molecule
    molecule.sgroups = sgroups;

    // Add conformer to molecule
    if num_atoms > 0 {
        let has_3d = is_3d(&positions);
        let mut conformer = Conformer::new(num_atoms, has_3d);
        for (idx, pos) in positions.into_iter().enumerate() {
            conformer.set_position(idx, pos);
        }
        molecule.add_conformer(conformer)?;
    }

    Ok(molecule)
}
