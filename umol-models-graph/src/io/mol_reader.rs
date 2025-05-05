//! MOL file reader.

use crate::atom::{Atom, AtomStereoParity};
use crate::bond::{Bond, BondDir, BondStereo, BondType};
use crate::conformer::{Conformer, Point3D};
use crate::molecule::Molecule;
use fixed_width::{FieldSet, FixedWidth, LineBreak, Reader};
use serde::Deserialize;
use std::collections::HashMap;
use std::io::BufRead;
use std::result::Result as StdResult;
use umol::error::{DataError, FormatError};

use umol::{Error, Result};
use umol_data::Element;

#[derive(Debug)]
#[allow(dead_code)]
enum AtomSymbol {
    Element(Element),
    AtomList,
    Unspecified(char),
    LonePair,
    RGroup(u8),
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct CountsLine {
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
struct AtomLine {
    x: f64,
    y: f64,
    z: f64,
    symbol: String,     // 'aaa' field - atom symbol (see AtomSymbol enum)
    mass_diff: i8, // 'dd' field - mass difference (-3, -2, -1, 0, 1, 2, 3, 4), 0 if value outside of this range
    charge: i8, // 'ccc' field - 0 = uncharged, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
    stereo_parity: i8, // 'sss' field - 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked
    hydrogen_count: i8, // 'hhh' field - 1 = H0, 2 = H1, 3 = H2, 4 = H3, 5 = H4
    stereo_care: i8, // 'bbb' field - 0 = ignore stereo, 1 = stereo in query must match
    valence: i8, // 'vvv' field - 0 = default, 1-14 = explicit, 15 = explicit 0
    // Skipping obsolete or unused fields (HHH, rrr, iii)
    atom_mapping: i8, // 'mmm' field - 1 = number of atoms
    inversion: i8,    // 'nnn' field - 0 = property not applied, 2 = inverted, 3 = retained
    exact_change: i8, // 'eee' field - 0 = property not applied, 1 = charge in query must match
}

impl FixedWidth for AtomLine {
    fn fields() -> FieldSet {
        FieldSet::Seq(vec![
            FieldSet::new_field(0..10).name("x"),
            FieldSet::new_field(10..20).name("y"),
            FieldSet::new_field(20..30).name("z"),
            FieldSet::new_field(31..34).name("symbol"),
            FieldSet::new_field(34..36).name("mass_diff"),
            FieldSet::new_field(36..39).name("charge"),
            FieldSet::new_field(39..42).name("stereo_parity"),
            FieldSet::new_field(42..45).name("hydrogen_count"),
            FieldSet::new_field(45..48).name("stereo_care"),
            FieldSet::new_field(48..51).name("valence"),
            // Skipping fields at positions 51-57 (HHH, rrr, iii)
            FieldSet::new_field(57..60).name("atom_mapping"),
            FieldSet::new_field(60..63).name("inversion"),
            FieldSet::new_field(63..66).name("exact_change"),
        ])
    }
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct BondLine {
    first_atom: i32,  // '111' - first atom number
    second_atom: i32, // '222' - second atom number
    bond_type: i32,   // 'ttt' - bond type (1=Single, 2=Double, etc.)
    bond_stereo: i32, // 'sss' - bond stereo (0=Not stereo, 1=Up, etc.)
    // xxx is not used, skipping
    bond_topology: i32,   // 'rrr' - bond topology (0=Either, 1=Ring, 2=Chain)
    reacting_center: i32, // 'ccc' - reacting center status
}

impl FixedWidth for BondLine {
    fn fields() -> FieldSet {
        FieldSet::Seq(vec![
            FieldSet::new_field(0..3).name("first_atom"),
            FieldSet::new_field(3..6).name("second_atom"),
            FieldSet::new_field(6..9).name("bond_type"),
            FieldSet::new_field(9..12).name("bond_stereo"),
            // Skip unused field xxx (12..15)
            FieldSet::new_field(15..18).name("bond_topology"),
            FieldSet::new_field(18..21).name("reacting_center"),
        ])
    }
}

fn detect_line_break(reader: &mut impl BufRead) -> Result<LineBreak> {
    let buffer = match reader.fill_buf() {
        Ok(buffer) if !buffer.is_empty() => buffer,
        _ => return Err(FormatError::InvalidMolFormat("Empty file".to_string()).into()),
    };
    for i in 0..buffer.len() - 1 {
        if buffer[i] == b'\r' && buffer[i + 1] == b'\n' {
            return Ok(LineBreak::CRLF);
        } else if buffer[i] == b'\n' {
            return Ok(LineBreak::Newline);
        }
    }
    Err(FormatError::InvalidMolFormat("Incomplete file".to_string()).into())
}

// TODO: The mix of match-based and subsequent processing is not ideal.
fn parse_atom_symbol(input: &str) -> Result<AtomSymbol> {
    let trimmed = input.trim();
    match trimmed {
        "L" => return Ok(AtomSymbol::AtomList),
        "A" => return Ok(AtomSymbol::Unspecified('A')),
        "Q" => return Ok(AtomSymbol::Unspecified('Q')),
        "*" => return Ok(AtomSymbol::Unspecified('*')),
        "LP" => return Ok(AtomSymbol::LonePair),
        _ => {}
    }
    if trimmed.starts_with('R') && trimmed.len() > 1 {
        if let Ok(number) = trimmed[1..].parse::<u8>() {
            return Ok(AtomSymbol::RGroup(number));
        }
    }
    Element::from_symbol(trimmed)
        .ok_or_else(|| {
            FormatError::InvalidMolFormat(format!("Unknown atom symbol: '{}'", trimmed)).into()
        })
        .map(AtomSymbol::Element)
}

fn parse_valence_code(code: i8) -> Result<Option<u8>> {
    // 'vvv' field - 0 = default, 1-14 = explicit, 15 = explicit 0
    match code {
        0 => Ok(None),                   // Code 0 means default/unspecified valence
        v @ 1..=14 => Ok(Some(v as u8)), // Codes 1-14 are explicit valences
        15 => Ok(Some(0)),               // Code 15 means explicit zero valence
        _ => Err(FormatError::InvalidMolFormat(format!("Invalid valence code '{}'", code)).into()),
    }
}

fn parse_charge_code(code: i8) -> Result<i8> {
    // 'ccc' field - 0 = uncharged, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
    match code {
        0 => Ok(0),
        1 => Ok(3),
        2 => Ok(2),
        3 => Ok(1),
        4 => Ok(0), // NOTE: Code 4 is doublet radical, not charge. Treat as 0 charge for now.
        5 => Ok(-1),
        6 => Ok(-2),
        7 => Ok(-3),
        _ => Err(FormatError::InvalidMolFormat(format!("Invalid charge code '{}'", code)).into()),
    }
}

fn parse_stereo_parity_code(code: i8) -> Result<Option<AtomStereoParity>> {
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

fn parse_hydrogen_count_code(code: i8) -> Result<Option<u8>> {
    match code {
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

fn parse_bond_type_code(code: i32) -> Result<BondType> {
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

fn parse_bond_stereo_code(code: i32) -> Result<Option<BondStereo>> {
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

fn parse_bond_dir_code(code: i32) -> Result<Option<BondDir>> {
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

/// Parses the common `M  XXX nn8 aaa vvv ...` structure from M lines.
/// Expects `parts` from `line.split_whitespace().collect()`.
/// `entry_size` is the number of fields per entry (e.g., 2 for `aaa vvv`).
fn parse_m_pairs(parts: &[&str], entry_size: usize) -> Result<Vec<(usize, i64)>> {
    if parts.len() < 3 {
        return Err(Error::from(FormatError::InvalidMolFormat(
            "M line too short".to_string(),
        )));
    }
    let count: usize = parts[2].parse().map_err(|_| {
        Error::from(FormatError::InvalidMolFormat(
            "Invalid count in M line".to_string(),
        ))
    })?;

    // Expected length: M, TYPE, count, then count * entry_size fields
    if parts.len() < 3 + count * entry_size {
        return Err(Error::from(FormatError::InvalidMolFormat(format!(
            "Insufficient data pairs in M line '{}' (expected {}, found {})",
            parts[1], // Property type
            count,
            (parts.len() - 3).checked_div(entry_size).unwrap_or(0)
        ))));
    }

    let mut pairs = Vec::with_capacity(count);
    for i in 0..count {
        let start_index = 3 + i * entry_size;
        let mol_idx_str = parts[start_index];
        let val_str = parts[start_index + 1]; // Assumes entry_size >= 2

        let mol_idx: usize = mol_idx_str.trim().parse().map_err(|_| {
            Error::from(FormatError::InvalidMolFormat(format!(
                "Invalid mol index '{}' in M line '{}'",
                mol_idx_str, parts[1]
            )))
        })?;
        let value: i64 = val_str.trim().parse().map_err(|_| {
            Error::from(FormatError::InvalidMolFormat(format!(
                "Invalid value '{}' in M line '{}'",
                val_str, parts[1]
            )))
        })?;
        pairs.push((mol_idx, value));
    }
    Ok(pairs)
}

/// Type for M line parsers
type MParserFn = fn(&mut Vec<Atom>, &[&str]) -> Result<()>;

/// Charge property parser
fn parse_m_chg(atoms: &mut Vec<Atom>, parts: &[&str]) -> Result<()> {
    let pairs = parse_m_pairs(parts, 2)?;
    for (mol_idx, charge_val) in pairs {
        let atom_vec_idx = mol_idx.checked_sub(1).ok_or_else(|| {
            Error::from(FormatError::InvalidMolFormat(format!(
                "M CHG line references invalid zero atom index",
            )))
        })?;
        if atom_vec_idx >= atoms.len() {
            return Err(Error::from(FormatError::InvalidMolFormat(format!(
                "M CHG line references out-of-bounds atom index {}",
                mol_idx
            ))));
        }
        let atom = &mut atoms[atom_vec_idx];

        if charge_val >= i8::MIN as i64 && charge_val <= i8::MAX as i64 {
            atom.formal_charge = charge_val as i8;
        } else {
            // TODO: Use logging instead of eprintln
            eprintln!(
                "Warning: Charge value {} for atom {} out of range (i8) in M CHG line, ignoring.",
                charge_val, mol_idx
            );
        }
    }
    Ok(())
}

/// Isotope property parser
fn parse_m_iso(atoms: &mut Vec<Atom>, parts: &[&str]) -> Result<()> {
    let pairs = parse_m_pairs(parts, 2)?;
    for (mol_idx, isotope_mass_val) in pairs {
        let atom_vec_idx = mol_idx.checked_sub(1).ok_or_else(|| {
            Error::from(FormatError::InvalidMolFormat(format!(
                "M ISO line references invalid zero atom index",
            )))
        })?;
        if atom_vec_idx >= atoms.len() {
            return Err(Error::from(FormatError::InvalidMolFormat(format!(
                "M ISO line references out-of-bounds atom index {}",
                mol_idx
            ))));
        }
        let atom = &mut atoms[atom_vec_idx];

        if isotope_mass_val < 0 {
            return Err(Error::from(FormatError::InvalidMolFormat(format!(
                "Invalid negative isotope mass {} for atom {} in M ISO line",
                isotope_mass_val, mol_idx
            ))));
        }

        // Calculate difference from standard mass
        let standard_mass = atom.element.atomic_mass(); // Assumes atomic_mass returns f64
                                                        // Use round() for comparison, but store as integer diff
        let diff_f64 = (isotope_mass_val as f64 - standard_mass).round();

        if diff_f64 >= i8::MIN as f64 && diff_f64 <= i8::MAX as f64 {
            let diff = diff_f64 as i8;
            if diff != 0 {
                atom.mass_difference = Some(diff);
            } else {
                atom.mass_difference = None; // Explicit 0 difference means default isotope
            }
        } else {
            // TODO: Use logging instead of eprintln
            eprintln!(
                "Warning: Calculated isotope mass difference ({}) for atom {} out of range (i8), ignoring M ISO.",
                diff_f64, mol_idx
            );
            atom.mass_difference = None; // Reset if out of range
        }
    }
    Ok(())
}

/// Radical property parser
fn parse_m_rad(atoms: &mut Vec<Atom>, parts: &[&str]) -> Result<()> {
    let pairs = parse_m_pairs(parts, 2)?;
    for (mol_idx, radical_val) in pairs {
        let atom_vec_idx = mol_idx.checked_sub(1).ok_or_else(|| {
            Error::from(FormatError::InvalidMolFormat(format!(
                "M RAD line references invalid zero atom index",
            )))
        })?;
        if atom_vec_idx >= atoms.len() {
            return Err(Error::from(FormatError::InvalidMolFormat(format!(
                "M RAD line references out-of-bounds atom index {}",
                mol_idx
            ))));
        }
        let atom = &mut atoms[atom_vec_idx];

        match radical_val {
            1 => atom.radical = Some(1), // Singlet
            2 => atom.radical = Some(2), // Doublet
            3 => atom.radical = Some(3), // Triplet
            0 => atom.radical = None,    // Explicitly non-radical
            _ => {
                // TODO: Use logging instead of eprintln
                eprintln!(
                    "Warning: Invalid radical value {} for atom {} in M RAD line, ignoring.",
                    radical_val, mol_idx
                );
                atom.radical = None; // Treat invalid as non-radical
            }
        }
    }
    Ok(())
}

fn is_3d(positions: &[Point3D]) -> bool {
    positions
        .iter()
        .any(|pos| pos.x * pos.x + pos.y * pos.y + pos.z * pos.z > f64::EPSILON)
}

// Parse V2000 MOL file.
pub fn read_mol_v2000(mut reader: impl BufRead) -> Result<Molecule> {
    let mut molecule = Molecule::new();
    let line_break = detect_line_break(&mut reader)?;

    // Parse header block: molecule name, info line, comment line
    let mut name_line = String::new();
    reader
        .read_line(&mut name_line)
        .map_err(|e| Error::from(FormatError::IoError(e)))?;
    if !name_line.trim().is_empty() {
        molecule.set_prop("mol_name".to_string(), name_line.to_string());
    }
    let mut _line_buffer = String::new();
    for _ in 0..1 {
        reader
            .read_line(&mut _line_buffer)
            .map_err(|e| Error::from(FormatError::IoError(e)))?;
        _line_buffer.clear();
    }

    // Counts line
    const COUNTS_LINE_WIDTH: usize = 39;
    let counts_data = Reader::from_reader(&mut reader)
        .width(COUNTS_LINE_WIDTH)
        .linebreak(line_break.clone())
        .byte_reader()
        .filter_map(StdResult::ok)
        .next()
        .ok_or(Error::from(FormatError::InvalidMolFormat(
            "Counts line not found".to_string(),
        )))?;
    let counts_data = fixed_width::from_bytes::<CountsLine>(&counts_data).map_err(|e| {
        FormatError::InvalidMolFormat(format!("Failed to parse counts line: {}", e))
    })?;
    if counts_data.version != "V2000" {
        return Err(
            FormatError::InvalidMolFormat("Only V2000 format supported".to_string()).into(),
        );
    }
    let num_atoms: usize = counts_data.atoms as usize;
    let num_bonds: usize = counts_data.bonds as usize;

    // Atom block
    const ATOM_LINE_WIDTH: usize = 66;
    let mut atom_reader = Reader::from_reader(&mut reader)
        .width(ATOM_LINE_WIDTH)
        .linebreak(line_break);
    let (mut atoms, positions) = atom_reader.byte_reader().take(num_atoms).try_fold(
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
                atom.formal_charge = parse_charge_code(atom_data.charge)?;
                atom.mass_difference = if atom_data.mass_diff == 0 {
                    None
                } else {
                    Some(atom_data.mass_diff as i8)
                };
                atom.stereo_parity = parse_stereo_parity_code(atom_data.stereo_parity)?;
                atom.valence = parse_valence_code(atom_data.valence)?;
                atom.explicit_hydrogens = parse_hydrogen_count_code(atom_data.hydrogen_count)?;
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

    // Bond block
    const BOND_LINE_WIDTH: usize = 21;
    let mut bond_reader = Reader::from_reader(&mut reader).width(BOND_LINE_WIDTH);
    let bonds = bond_reader.byte_reader().take(num_bonds).try_fold(
        Vec::with_capacity(num_bonds),
        |mut bonds, bytes_result| -> Result<_> {
            let bytes = bytes_result.map_err(|e| {
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

            let idx1 = bond_data.first_atom as usize;
            let idx2 = bond_data.second_atom as usize;

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

    // Define M property parsers before using them
    let m_parsers: HashMap<&'static str, MParserFn> = [
        ("CHG", parse_m_chg as MParserFn),
        ("ISO", parse_m_iso as MParserFn),
        ("RAD", parse_m_rad as MParserFn),
        // TODO: Add other parsers here
    ]
    .iter()
    .cloned()
    .collect();

    // Properties block - Parse M lines using try_for_each
    let terminated = reader
        .lines()
        .map(|line| line.map_err(|e| Error::from(FormatError::IoError(e))))
        .try_fold(false, |mut terminated, line| -> Result<bool> {
            let line = line?;
            if line.starts_with("M  END") {
                terminated = true;
            } else if line.starts_with("M  ") {
                if terminated {
                    return Err(Error::from(FormatError::InvalidMolFormat(
                        "M line found after M END".to_string(),
                    )));
                }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let Some(parser_fn) = m_parsers.get(parts[1]) {
                        parser_fn(&mut atoms, &parts)?;
                    } else {
                        // Optional: Log unknown M property type
                        // eprintln!("Warning: Unknown M property type '{}'", parts[1]);
                    }
                }
                // Silently ignore malformed M lines (parts.len() < 2)
            }
            // Silently ignore non-"M  " lines within the properties block
            Ok(terminated)
        })?;

    if !terminated {
        return Err(Error::from(FormatError::InvalidMolFormat(
            "M END not found".to_string(),
        )));
    }

    // Add atoms (potentially modified by properties) to the molecule
    for (idx, atom) in atoms.into_iter().enumerate() {
        molecule.add_atom(idx + 1, atom);
    }

    // Add bonds to the molecule
    for (idx1, idx2, bond) in bonds {
        molecule.add_bond(idx1, idx2, bond)?;
    }

    // Add conformer
    if num_atoms > 0 {
        let has_3d = is_3d(&positions);
        let mut conformer = Conformer::new(num_atoms, has_3d);
        for (ext_idx, pos) in positions.into_iter().enumerate() {
            let mol_idx = ext_idx + 1;
            let graph_idx = *molecule
                .external_indices
                .get(&mol_idx)
                .ok_or_else(|| Error::from(DataError::MissingAtomIndex(mol_idx)))?;
            conformer.set_position(graph_idx, pos);
        }
        molecule.add_conformer(conformer)?;
    }

    Ok(molecule)
}
