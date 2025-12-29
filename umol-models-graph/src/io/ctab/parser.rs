//! CTab format parser.
//!
//! This module provides the main entry points for parsing Connection Table (CTAB) format,
//! which is the core molecular structure representation used in MOL, SDF, and other formats.

use bstr::{join, ByteSlice};
use nom::bytes::complete::{is_not, tag};
use nom::character::complete::{line_ending, multispace0};
use nom::combinator::{all_consuming, map_parser, opt};
use nom::sequence::terminated;
use nom::{Err as NomErr, Parser};

use self::accumulator::MoleculeProperties;
pub use self::atom::{atom_input, extended_atom_input}; // NOTE: Re-exported for benchmarks
pub use self::bond::{bond_input, extended_bond_input}; // NOTE: Re-exported for benchmarks
pub use self::counts::counts_input;
use self::counts::Counts;
pub use self::legacy_atom_list::legacy_atom_list_input; // NOTE: Re-exported for benchmarks
use self::properties::{atom_alias_input, PropertyEntries};
pub use self::properties::{extended_property_input, property_input}; // NOTE: Re-exported for benchmarks
use super::config::CtabParseFlags;
use crate::io::ctfile::error::ParseError;
use crate::position::{all_zero, Point3D};
use crate::table_ir::bond::Bond;
use crate::table_ir::source::SourceFormat;
use crate::table_ir::{Atom, ExtendedAtom, ExtendedBond, ExtendedMolecule, Molecule};

mod accumulator;
mod atom;
mod bond;
mod context;
mod convert;
mod counts;
mod legacy_atom_list;
mod properties;
mod rgroup;
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

        let extended = build_extended_molecule(atoms, bonds, positions, properties, flags)
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

/// Parse atom block (basic atoms only)
fn atom_block<'inp>(
    atom_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<Atom>, Option<Vec<Point3D>>, u32), Error = ParseError>
       + use<'inp> {
    move |input: &'inp [u8]| {
        let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
        let mut atoms = Vec::with_capacity(atom_count as usize);
        let mut positions = Vec::with_capacity(if ignore_positions {
            0
        } else {
            atom_count as usize
        });
        let mut lines_iter = input.lines_with_terminator();
        let mut offset = 0;

        for i in 0..atom_count {
            let line = lines_iter.next().ok_or_else(|| {
                NomErr::Error(ParseError::UnexpectedEof {
                    line: line_offset + i,
                    block: "atom",
                })
            })?;

            let (_, (atom, pos)) = all_consuming(atom_input(flags))
                .parse(line)
                .map_err(|e| NomErr::Error(ParseError::atom_from_nom(e, line_offset + i, line)))?;
            atoms.push(atom);
            if !ignore_positions {
                positions.push(pos);
            }
            offset += line.len();
        }

        let remaining = &input[offset..];
        if ignore_positions || (atom_count > 1 && all_zero(&positions)) {
            Ok((remaining, (atoms, None, line_offset + atom_count)))
        } else {
            Ok((
                remaining,
                (atoms, Some(positions), line_offset + atom_count),
            ))
        }
    }
}

/// Parse extended atom block
fn extended_atom_block<'inp>(
    atom_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<
    &'inp [u8],
    Output = (Vec<ExtendedAtom>, Option<Vec<Point3D>>, u32),
    Error = ParseError,
> + use<'inp> {
    move |input: &'inp [u8]| {
        let ignore_positions = flags.contains(CtabParseFlags::IGNORE_POSITIONS);
        let mut atoms = Vec::with_capacity(atom_count as usize);
        let mut positions = Vec::with_capacity(if ignore_positions {
            0
        } else {
            atom_count as usize
        });
        let mut lines_iter = input.lines_with_terminator();
        let mut offset = 0;

        for i in 0..atom_count {
            let line = lines_iter.next().ok_or_else(|| {
                NomErr::Error(ParseError::UnexpectedEof {
                    line: line_offset + i,
                    block: "atom",
                })
            })?;

            let (_, (atom, pos)) = all_consuming(extended_atom_input(flags))
                .parse(line)
                .map_err(|e| NomErr::Error(ParseError::atom_from_nom(e, line_offset + i, line)))?;
            atoms.push(atom);
            if !ignore_positions {
                positions.push(pos);
            }
            offset += line.len();
        }

        let remaining = &input[offset..];
        if ignore_positions || (atom_count > 1 && all_zero(&positions)) {
            Ok((remaining, (atoms, None, line_offset + atom_count)))
        } else {
            Ok((
                remaining,
                (atoms, Some(positions), line_offset + atom_count),
            ))
        }
    }
}

/// Parse bond block (basic bonds only)
fn bond_block<'inp>(
    bond_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<(usize, usize, Bond)>, u32), Error = ParseError> + use<'inp>
{
    move |input: &'inp [u8]| {
        let mut bonds = Vec::with_capacity(bond_count as usize);
        let mut lines_iter = input.lines_with_terminator();
        let mut offset = 0;

        for i in 0..bond_count {
            let line = lines_iter.next().ok_or_else(|| {
                NomErr::Error(ParseError::UnexpectedEof {
                    line: line_offset + i,
                    block: "bond",
                })
            })?;

            let (_, (atom1, atom2, bond)) = all_consuming(bond_input(flags))
                .parse(line)
                .map_err(|e| NomErr::Error(ParseError::bond_from_nom(e, line_offset + i, line)))?;
            bonds.push((atom1, atom2, bond));
            offset += line.len();
        }

        let remaining = &input[offset..];
        Ok((remaining, (bonds, line_offset + bond_count)))
    }
}

/// Parse extended bond block
fn extended_bond_block<'inp>(
    bond_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<(usize, usize, ExtendedBond)>, u32), Error = ParseError>
       + use<'inp> {
    move |input: &'inp [u8]| {
        let mut bonds = Vec::with_capacity(bond_count as usize);
        let mut lines_iter = input.lines_with_terminator();
        let mut offset = 0;

        for i in 0..bond_count {
            let line = lines_iter.next().ok_or_else(|| {
                NomErr::Error(ParseError::UnexpectedEof {
                    line: line_offset + i,
                    block: "bond",
                })
            })?;

            let (_, (atom1, atom2, bond)) = all_consuming(extended_bond_input(flags))
                .parse(line)
                .map_err(|e| {
                NomErr::Error(ParseError::bond_from_nom(e, line_offset + i, line))
            })?;
            bonds.push((atom1, atom2, bond));
            offset += line.len();
        }

        let remaining = &input[offset..];
        Ok((remaining, (bonds, line_offset + bond_count)))
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
        source_format: SourceFormat::MOL,
        positions,
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
        source_format: SourceFormat::MOL,
        positions,
        fragments: Vec::new(),
        links: Vec::new(),
        electrons: None,
        properties: Vec::new(),
        comments: Vec::new(),
        ctfile_data: None,
    };

    let mut acc = MoleculeProperties::new();
    for entry in properties {
        acc.add_entry(entry, flags)?;
    }

    acc.update_extended_molecule(&mut molecule, flags)?;

    Ok(molecule)
}

#[cfg(test)]
mod tests;
