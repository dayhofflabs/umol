//! CTFile format parsers (CTAB, MOL, SDF).
//!
//! This module provides the main entry points for parsing Connection Table formats.

use std::borrow::Cow;

use bstr::{join, ByteSlice};
use indexmap::IndexMap;
use nom::bytes::complete::{is_not, tag};
use nom::character::complete::{line_ending, multispace0};
use nom::combinator::{all_consuming, map_parser, opt};
use nom::sequence::terminated;
use nom::{Err as NomErr, Parser};

use self::accumulator::MoleculeProperties;
use self::atom::{atom_block, extended_atom_block};
pub use self::atom::{atom_input, extended_atom_input}; // NOTE: Re-exported for benchmarks
use self::bond::{bond_block, extended_bond_block};
pub use self::bond::{bond_input, extended_bond_input}; // NOTE: Re-exported for benchmarks
pub use self::counts::counts_input;
use self::counts::Counts;
use self::header::header_block;
pub use self::legacy_atom_list::legacy_atom_list_input; // NOTE: Re-exported for benchmarks
use self::properties::{atom_alias_input, PropertyEntries};
pub use self::properties::{extended_property_input, property_input}; // NOTE: Re-exported for benchmarks
use self::sdf_data::{sdf_data_block, sdf_delimiter};
use super::config::{CtabParseFlags, CtfileIoConfig};
use super::error::ParseError;
use crate::io::utils::normalize_whitespace;
use crate::position::Point3D;
use crate::table_ir::bond::Bond;
use crate::table_ir::source::SourceFormat;
use crate::table_ir::{Atom, AtomSymbol, ExtendedAtom, ExtendedBond, ExtendedMolecule, Molecule};

mod accumulator;
mod atom;
mod bond;
mod context;
mod convert;
mod counts;
mod header;
mod legacy_atom_list;
mod properties;
mod rgroup;
mod sdf_data;
mod sgroup;
mod utils;

/// Parse CTAB block (basic parser, optimized for performance, basic molecules only)
///
/// This parser is optimized for basic molecules without query features.
/// It will fail if the CTAB contains query atoms, query bonds, or other advanced features.
pub fn ctab_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Molecule, Error = ParseError> + use<'inp> {
    debug_assert!(
        CtabParseFlags::BASIC_MAX.contains(flags),
        "flags must be a subset of BASIC_MAX"
    );
    let legacy_atom_lists = flags.contains(CtabParseFlags::LEGACY_ATOM_LISTS);

    move |input: &'inp [u8]| {
        let (remaining, (counts, line_offset)) = counts_block(line_offset, flags).parse(input)?;
        let atom_count = counts.atom_count;
        let bond_count = counts.bond_count;
        let atom_list_count = counts.atom_list_count;

        if !legacy_atom_lists && atom_list_count > 0 {
            return Err(NomErr::Error(ParseError::InvalidPropertyLine {
                line: line_offset + atom_count + bond_count,
                col: 0,
            }));
        }

        let (remaining, (atoms, positions, line_offset)) =
            atom_block(atom_count, line_offset, flags).parse(remaining)?;

        let (remaining, (bonds, line_offset)) =
            bond_block(bond_count, line_offset, flags).parse(remaining)?;

        let (remaining, (legacy_properties, line_offset)) =
            legacy_atom_list_block(atom_list_count, line_offset, flags).parse(remaining)?;

        let (remaining, (properties, line_offset)) =
            properties_block(line_offset, flags).parse(remaining)?;

        let properties = if !legacy_properties.is_empty() {
            properties.into_iter().chain(legacy_properties).collect()
        } else {
            properties
        };

        let (remaining, _) = end_block(line_offset, flags).parse(remaining)?;

        let molecule = build_molecule(atoms, bonds, positions, properties, flags)
            .map_err(|e| NomErr::Error(e))?;
        Ok((remaining, molecule))
    }
}

/// Parse CTAB block (general parser, handles all features including queries)
///
/// This parser handles all CTAB features including query atoms, query bonds,
/// R-groups, S-groups, and all property lines.
pub fn extended_ctab_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = ExtendedMolecule, Error = ParseError> + use<'inp> {
    let legacy_atom_lists = flags.contains(CtabParseFlags::LEGACY_ATOM_LISTS);

    move |input: &'inp [u8]| {
        let (remaining, (counts, line_offset)) = counts_block(line_offset, flags).parse(input)?;
        let atom_count = counts.atom_count;
        let bond_count = counts.bond_count;
        let atom_list_count = counts.atom_list_count;

        if !legacy_atom_lists && atom_list_count > 0 {
            return Err(NomErr::Error(ParseError::InvalidPropertyLine {
                line: line_offset + atom_count + bond_count,
                col: 0,
            }));
        }

        let (remaining, (atoms, positions, line_offset)) =
            extended_atom_block(atom_count, line_offset, flags).parse(remaining)?;

        let (remaining, (bonds, line_offset)) =
            extended_bond_block(bond_count, line_offset, flags).parse(remaining)?;

        let (remaining, (legacy_properties, line_offset)) = if legacy_atom_lists {
            legacy_atom_list_block(atom_list_count, line_offset, flags).parse(remaining)?
        } else {
            (remaining, (Vec::new(), line_offset + atom_list_count))
        };

        let (remaining, (properties, line_offset)) =
            extended_properties_block(line_offset, flags).parse(remaining)?;

        let properties = if !legacy_properties.is_empty() {
            properties.into_iter().chain(legacy_properties).collect()
        } else {
            properties
        };

        let (remaining, _) = end_block(line_offset, flags).parse(remaining)?;

        let extended = build_extended_molecule(
            atoms,
            bonds,
            positions,
            properties,
            flags,
        )
        .map_err(|e| NomErr::Error(e))?;
        Ok((remaining, extended))
    }
}

/// Parse counts block
fn counts_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Counts, u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let (remaining, counts) =
            map_parser(terminated(is_not("\r\n"), line_ending), counts_input(flags))
                .parse(input)
                .map_err(|e| NomErr::Error(ParseError::counts_from_nom(e, line_offset)))?;

        Ok((remaining, (counts, line_offset + 1)))
    }
}

// Parse legacy atom list block
fn legacy_atom_list_block<'inp>(
    atom_list_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<PropertyEntries>, u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let mut properties = Vec::with_capacity(atom_list_count as usize);
        let mut lines_iter = input.lines_with_terminator();
        let mut offset = 0;

        for i in 0..atom_list_count {
            let line = lines_iter.next().ok_or_else(|| {
                NomErr::Error(ParseError::UnexpectedEof {
                    line: line_offset + i,
                    block: "legacy atom list",
                })
            })?;

            let (_, property) = all_consuming(legacy_atom_list_input(flags))
                .parse(line)
                .map_err(|e| {
                    NomErr::Error(ParseError::legacy_atom_list_from_nom(
                        e,
                        line_offset + i,
                        line,
                    ))
                })?;
            properties.push(property);
            offset += line.len();
        }

        let remaining = &input[offset..];
        Ok((remaining, (properties, line_offset + atom_list_count)))
    }
}

/// Parse properties block (basic properties only)
fn properties_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<PropertyEntries>, u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let mut properties = Vec::new();
        let mut lines_iter = input.lines_with_terminator();
        let mut offset = 0;
        let mut line_count = 0;

        while let Some(line) = lines_iter.next() {
            if line.starts_with(b"M  END") {
                // Include M  END line in remaining for termination
                break;
            }

            // Handle atom alias (two-line property)
            if line.starts_with(b"A  ") {
                if let Some(next_line) = lines_iter.next() {
                    let combined_line = join(b"\n", [line, next_line]);
                    let line_bytes = line.len() + next_line.len();

                    if let Ok((_, property)) =
                        all_consuming(atom_alias_input()).parse(&combined_line)
                    {
                        properties.push(property);
                        offset += line_bytes;
                        line_count += 2;
                    } else {
                        break; // Backtrack
                    };
                } else {
                    break; // Incomplete atom alias
                }
            } else {
                match all_consuming(property_input(flags)).parse(line) {
                    Ok((_, property)) => {
                        properties.push(property);
                        offset += line.len();
                        line_count += 1;
                    }
                    Err(_) => {
                        return Err(NomErr::Error(ParseError::InvalidPropertyLine {
                            line: line_offset + line_count,
                            col: 0,
                        }));
                    }
                }
            }
        }

        let remaining = &input[offset..];
        Ok((remaining, (properties, line_offset + line_count)))
    }
}

/// Parse extended properties block
fn extended_properties_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<PropertyEntries>, u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let mut properties = Vec::new();
        let mut lines_iter = input.lines_with_terminator();
        let mut offset = 0;
        let mut line_count = 0;

        while let Some(line) = lines_iter.next() {
            if line.starts_with(b"M  END") {
                break; // Include M  END line in remaining for termination
            }

            // Handle atom alias (two-line property)
            if line.starts_with(b"A  ") {
                if let Some(next_line) = lines_iter.next() {
                    let combined_line = join(b"\n", [line, next_line]);
                    let line_bytes = line.len() + next_line.len();

                    if let Ok((_, property)) =
                        all_consuming(atom_alias_input()).parse(&combined_line)
                    {
                        properties.push(property);
                        offset += line_bytes;
                        line_count += 2;
                    } else {
                        break; // Backtrack
                    };
                } else {
                    break; // Incomplete atom alias
                }
            } else if let Ok((_, property)) =
                all_consuming(extended_property_input(flags)).parse(line)
            {
                properties.push(property);
                offset += line.len();
                line_count += 1;
            } else {
                break; // Backtrack
            }
        }

        let remaining = &input[offset..];
        Ok((remaining, (properties, line_offset + line_count)))
    }
}

/// Parse end block (M END line)
/// Accepts "M  END" followed by optional whitespace, or nothing.
fn end_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = ((), u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let no_v2000_end_tags = flags.contains(CtabParseFlags::NO_V2000_END_TAGS);
        if no_v2000_end_tags {
            let (remaining, found) = opt(terminated(tag("M  END"), multispace0))
                .parse(input)
                .map_err(|e| NomErr::Error(ParseError::m_end_from_nom(e, line_offset + 1)))?;
            let line_count = if found.is_some() { 1 } else { 0 };
            return Ok((remaining, ((), line_offset + line_count)));
        } else {
            let (remaining, _) = terminated(tag("M  END"), multispace0)
                .parse(input)
                .map_err(|e| NomErr::Error(ParseError::m_end_from_nom(e, line_offset + 1)))?;
            return Ok((remaining, ((), line_offset + 1)));
        }
    }
}

/// Build Molecule
fn build_molecule(
    atoms: Vec<Atom>,
    bonds: Vec<(usize, usize, Bond)>,
    positions: Option<Vec<Point3D>>,
    properties: Vec<PropertyEntries>,
    flags: CtabParseFlags,
) -> Result<Molecule, ParseError> {
    let bonds: Vec<Bond> = bonds
        .into_iter()
        .map(|(idx1, idx2, mut bond)| {
            bond.set_atoms(idx1 as u32, idx2 as u32);
            bond
        })
        .collect();

    let mut molecule = Molecule {
        atoms,
        bonds,
        rings: Vec::new(),
        positions,
        comments: Vec::new(),
        properties: IndexMap::new(),
        source_format: SourceFormat::MOL,
    };

    let mut acc = MoleculeProperties::new();
    for entry in properties {
        acc.add_entry(entry, flags)?;
    }
    acc.update_molecule(&mut molecule, flags)?;

    Ok(molecule)
}

/// Build extended molecule
fn build_extended_molecule(
    atoms: Vec<ExtendedAtom>,
    bonds: Vec<(usize, usize, ExtendedBond)>,
    positions: Option<Vec<Point3D>>,
    properties: Vec<PropertyEntries>,
    flags: CtabParseFlags,
) -> Result<ExtendedMolecule, ParseError> {
    let bonds: Vec<ExtendedBond> = bonds
        .into_iter()
        .map(|(idx1, idx2, mut bond)| {
            bond.set_atoms(idx1 as u32, idx2 as u32);
            bond
        })
        .collect();

    let mut molecule = ExtendedMolecule {
        atoms,
        bonds,
        rings: Vec::new(),
        positions,
        fragments: Vec::new(),
        links: Vec::new(),
        electrons: None,
        comments: Vec::new(),
        properties: IndexMap::new(),
        ctfile_data: None,
        source_format: SourceFormat::MOL,
    };

    let mut acc = MoleculeProperties::new();
    for entry in properties {
        acc.add_entry(entry, flags)?;
    }

    acc.update_extended_molecule(&mut molecule, flags)?;

    Ok(molecule)
}

/// Check if extended molecule contains extended features
///
/// Returns true if the molecule contains features not supported in basic MOL format:
/// - Extended atom symbols (atom lists, R-groups, etc.)
/// - Extended bond types (SingleOrDouble, Any, Zero, etc.)
/// - S-groups, R-groups
pub fn has_extended_features(molecule: &ExtendedMolecule) -> bool {
    for atom in &molecule.atoms {
        match &atom.symbol {
            AtomSymbol::WildcardAtom(_)
            | AtomSymbol::AtomList(_)
            | AtomSymbol::RGroup(_)
            | AtomSymbol::LonePair
            | AtomSymbol::Pseudoatom(_) => return true,
            AtomSymbol::Element(_) | AtomSymbol::NamedIsotope(_) => {}
        }
    }

    for bond in &molecule.bonds {
        if bond.order.is_query() || bond.order.is_extended() {
            return true;
        }
    }

    if !molecule.sgroups().is_empty() {
        return true;
    }

    if !molecule.rgroups().is_empty() {
        return true;
    }

    false
}

/// Parse MOL bytes into a Molecule with options (optimized, basic molecules only)
pub fn parse_mol_bytes_with(input: &[u8], config: &CtfileIoConfig) -> Result<Molecule, ParseError> {
    let flags = config.parse_flags;

    let data: Cow<'_, [u8]> = if flags.contains(CtabParseFlags::UNICODE) {
        normalize_whitespace(input)
    } else {
        Cow::Borrowed(input)
    };

    let (remaining, comments) = header_block()
        .parse(&*data)
        .map_err(|e| ParseError::header_from_nom(e, 0))?;

    let (_, mut molecule) = ctab_block(3, flags).parse(remaining)?;

    molecule.comments = comments;

    Ok(molecule)
}

/// Parse MOL bytes into a Molecule
///
/// Parses MOL file with basic flags. Returns error if molecule contains
/// extended atom/bond types (use parse_extended_mol_bytes for those).
pub fn parse_mol_bytes(input: &[u8]) -> Result<Molecule, ParseError> {
    let config = CtfileIoConfig::basic();
    parse_mol_bytes_with(input, &config)
}

/// Parse MOL string into a Molecule with options (optimized, basic molecules only)
pub fn parse_mol_with(input: &str, config: &CtfileIoConfig) -> Result<Molecule, ParseError> {
    parse_mol_bytes_with(input.as_bytes(), config)
}

/// Parse MOL string into a Molecule (optimized, basic molecules only)
///
/// Optimized parsing function for basic molecules.
/// Fails if the MOL file contains extended features.
pub fn parse_mol(input: &str) -> Result<Molecule, ParseError> {
    parse_mol_bytes(input.as_bytes())
}

/// Parse MOL bytes into an ExtendedMolecule with options
pub fn parse_extended_mol_bytes_with(
    input: &[u8],
    config: &CtfileIoConfig,
) -> Result<ExtendedMolecule, ParseError> {
    let flags = config.parse_flags;

    let data: Cow<'_, [u8]> = if flags.contains(CtabParseFlags::UNICODE) {
        normalize_whitespace(input)
    } else {
        Cow::Borrowed(input)
    };

    let (remaining, comments) = header_block()
        .parse(&*data)
        .map_err(|e| ParseError::header_from_nom(e, 0))?;

    let (_, mut molecule) = extended_ctab_block(3, flags).parse(remaining)?;

    molecule.comments = comments;

    Ok(molecule)
}

/// Parse MOL bytes into an ExtendedMolecule
///
/// Generic parsing function that handles both basic and extended molecules.
pub fn parse_extended_mol_bytes(input: &[u8]) -> Result<ExtendedMolecule, ParseError> {
    let config = CtfileIoConfig::extended();
    parse_extended_mol_bytes_with(input, &config)
}

/// Parse MOL string into an ExtendedMolecule with options
pub fn parse_extended_mol_with(
    input: &str,
    config: &CtfileIoConfig,
) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_mol_bytes_with(input.as_bytes(), config)
}

/// Parse MOL string into an ExtendedMolecule
pub fn parse_extended_mol(input: &str) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_mol_bytes(input.as_bytes())
}

/// Parse single SDF compound into Molecule (basic, optimized)
fn parse_sdf_molecule<'inp>(
    input: &'inp [u8],
    line_offset: u32,
    config: &CtfileIoConfig,
) -> Result<(&'inp [u8], Molecule, u32), ParseError> {
    let flags = config.parse_flags;

    let (remaining, comments) = header_block()
        .parse(input)
        .map_err(|e| ParseError::header_from_nom(e, line_offset))?;

    let (remaining, mut molecule) = ctab_block(line_offset + 3, flags).parse(remaining)?;

    let consumed_for_mol = input.len() - remaining.len();
    let mol_lines = input[..consumed_for_mol].lines_with_terminator().count() as u32;

    let (remaining, data_fields) = sdf_data_block()
        .parse(remaining)
        .map_err(|e| ParseError::sdf_data_from_nom(e, line_offset + mol_lines))?;

    let (remaining, _) = opt(sdf_delimiter())
        .parse(remaining)
        .map_err(|e| ParseError::delimiter_from_nom(e, line_offset + mol_lines))?;

    molecule.comments = comments;
    molecule.properties = data_fields;

    let consumed = input.len() - remaining.len();
    let lines_consumed = input[..consumed].lines_with_terminator().count() as u32;

    Ok((remaining, molecule, lines_consumed))
}

/// Parse single SDF compound into ExtendedMolecule
fn parse_sdf_extended_molecule<'inp>(
    input: &'inp [u8],
    line_offset: u32,
    config: &CtfileIoConfig,
) -> Result<(&'inp [u8], ExtendedMolecule, u32), ParseError> {
    let flags = config.parse_flags;

    let (remaining, comments) = header_block()
        .parse(input)
        .map_err(|e| ParseError::header_from_nom(e, line_offset))?;

    let (remaining, mut molecule) = extended_ctab_block(line_offset + 3, flags).parse(remaining)?;

    let consumed_for_mol = input.len() - remaining.len();
    let mol_lines = input[..consumed_for_mol].lines_with_terminator().count() as u32;

    let (remaining, data_fields) = sdf_data_block()
        .parse(remaining)
        .map_err(|e| ParseError::sdf_data_from_nom(e, line_offset + mol_lines))?;

    let (remaining, _) = opt(sdf_delimiter())
        .parse(remaining)
        .map_err(|e| ParseError::delimiter_from_nom(e, line_offset + mol_lines))?;

    molecule.comments = comments;
    molecule.properties = data_fields;

    let consumed = input.len() - remaining.len();
    let lines_consumed = input[..consumed].lines_with_terminator().count() as u32;

    Ok((remaining, molecule, lines_consumed))
}

/// Parse SDF bytes into Vec<Molecule> with config (basic, optimized)
pub fn parse_sdf_bytes_with(
    input: &[u8],
    config: &CtfileIoConfig,
) -> Result<Vec<Molecule>, ParseError> {
    let flags = config.parse_flags;

    let data: Cow<'_, [u8]> = if flags.contains(CtabParseFlags::UNICODE) {
        normalize_whitespace(input)
    } else {
        Cow::Borrowed(input)
    };

    let mut molecules = Vec::new();
    let mut line_offset = 0u32;
    let mut remaining: &[u8] = &data;

    while !remaining.trim_ascii().is_empty() {
        let (rem, molecule, lines_consumed) = parse_sdf_molecule(remaining, line_offset, config)?;
        molecules.push(molecule);
        line_offset += lines_consumed;
        remaining = rem;
    }

    Ok(molecules)
}

/// Parse SDF bytes into Vec<Molecule>
pub fn parse_sdf_bytes(input: &[u8]) -> Result<Vec<Molecule>, ParseError> {
    let config = CtfileIoConfig::basic();
    parse_sdf_bytes_with(input, &config)
}

/// Parse SDF string into Vec<Molecule> with config
pub fn parse_sdf_with(input: &str, config: &CtfileIoConfig) -> Result<Vec<Molecule>, ParseError> {
    parse_sdf_bytes_with(input.as_bytes(), config)
}

/// Parse SDF string into Vec<Molecule>
pub fn parse_sdf(input: &str) -> Result<Vec<Molecule>, ParseError> {
    parse_sdf_bytes(input.as_bytes())
}

/// Parse SDF bytes into Vec<ExtendedMolecule> with config
pub fn parse_extended_sdf_bytes_with(
    input: &[u8],
    config: &CtfileIoConfig,
) -> Result<Vec<ExtendedMolecule>, ParseError> {
    let flags = config.parse_flags;

    let data: Cow<'_, [u8]> = if flags.contains(CtabParseFlags::UNICODE) {
        normalize_whitespace(input)
    } else {
        Cow::Borrowed(input)
    };

    let mut molecules = Vec::new();
    let mut line_offset = 0u32;
    let mut remaining: &[u8] = &data;

    while !remaining.trim_ascii().is_empty() {
        let (rem, molecule, lines_consumed) =
            parse_sdf_extended_molecule(remaining, line_offset, config)?;
        molecules.push(molecule);
        line_offset += lines_consumed;
        remaining = rem;
    }

    Ok(molecules)
}

/// Parse SDF bytes into Vec<ExtendedMolecule>
pub fn parse_extended_sdf_bytes(input: &[u8]) -> Result<Vec<ExtendedMolecule>, ParseError> {
    let config = CtfileIoConfig::extended();
    parse_extended_sdf_bytes_with(input, &config)
}

/// Parse SDF string into Vec<ExtendedMolecule> with config
pub fn parse_extended_sdf_with(
    input: &str,
    config: &CtfileIoConfig,
) -> Result<Vec<ExtendedMolecule>, ParseError> {
    parse_extended_sdf_bytes_with(input.as_bytes(), config)
}

/// Parse SDF string into Vec<ExtendedMolecule>
pub fn parse_extended_sdf(input: &str) -> Result<Vec<ExtendedMolecule>, ParseError> {
    parse_extended_sdf_bytes(input.as_bytes())
}

/// Iterator for lazy streaming SDF parsing into Molecule
pub struct SdfIter<'inp> {
    data: Cow<'inp, [u8]>,
    offset: usize,
    line_offset: u32,
    config: CtfileIoConfig,
}

impl<'inp> SdfIter<'inp> {
    fn new(input: &'inp [u8], config: CtfileIoConfig) -> Self {
        let data = if config.parse_flags.contains(CtabParseFlags::UNICODE) {
            normalize_whitespace(input)
        } else {
            Cow::Borrowed(input)
        };
        Self {
            data,
            offset: 0,
            line_offset: 0,
            config,
        }
    }

    fn remaining(&self) -> &[u8] {
        &self.data[self.offset..]
    }
}

impl<'inp> Iterator for SdfIter<'inp> {
    type Item = Result<Molecule, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let remaining = self.remaining();
        if remaining.trim_ascii().is_empty() {
            return None;
        }

        match parse_sdf_molecule(remaining, self.line_offset, &self.config) {
            Ok((rem, molecule, lines_consumed)) => {
                self.offset = self.data.len() - rem.len();
                self.line_offset += lines_consumed;
                Some(Ok(molecule))
            }
            Err(e) => {
                self.offset = self.data.len();
                Some(Err(e))
            }
        }
    }
}

/// Iterator for lazy streaming SDF parsing into ExtendedMolecule
pub struct ExtendedSdfIter<'inp> {
    data: Cow<'inp, [u8]>,
    offset: usize,
    line_offset: u32,
    config: CtfileIoConfig,
}

impl<'inp> ExtendedSdfIter<'inp> {
    fn new(input: &'inp [u8], config: CtfileIoConfig) -> Self {
        let data = if config.parse_flags.contains(CtabParseFlags::UNICODE) {
            normalize_whitespace(input)
        } else {
            Cow::Borrowed(input)
        };
        Self {
            data,
            offset: 0,
            line_offset: 0,
            config,
        }
    }

    fn remaining(&self) -> &[u8] {
        &self.data[self.offset..]
    }
}

impl<'inp> Iterator for ExtendedSdfIter<'inp> {
    type Item = Result<ExtendedMolecule, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let remaining = self.remaining();
        if remaining.trim_ascii().is_empty() {
            return None;
        }

        match parse_sdf_extended_molecule(remaining, self.line_offset, &self.config) {
            Ok((rem, molecule, lines_consumed)) => {
                self.offset = self.data.len() - rem.len();
                self.line_offset += lines_consumed;
                Some(Ok(molecule))
            }
            Err(e) => {
                self.offset = self.data.len();
                Some(Err(e))
            }
        }
    }
}

/// Create a lazy streaming iterator for SDF bytes into Molecule (basic)
pub fn parse_sdf_iter(input: &[u8]) -> SdfIter<'_> {
    SdfIter::new(input, CtfileIoConfig::basic_max())
}

/// Create a lazy streaming iterator for SDF bytes into Molecule with config
pub fn parse_sdf_iter_with(input: &[u8], config: CtfileIoConfig) -> SdfIter<'_> {
    SdfIter::new(input, config)
}

/// Create a lazy streaming iterator for SDF bytes into ExtendedMolecule
pub fn parse_extended_sdf_iter(input: &[u8]) -> ExtendedSdfIter<'_> {
    ExtendedSdfIter::new(input, CtfileIoConfig::extended())
}

/// Create a lazy streaming iterator for SDF bytes into ExtendedMolecule with config
pub fn parse_extended_sdf_iter_with(input: &[u8], config: CtfileIoConfig) -> ExtendedSdfIter<'_> {
    ExtendedSdfIter::new(input, config)
}

#[cfg(test)]
mod tests;
