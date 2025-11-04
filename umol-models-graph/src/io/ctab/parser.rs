//! CTab format parser.
//!
//! This module provides the main entry points for parsing Connection Table (CTAB) format,
//! which is the core molecular structure representation used in MOL, SDF, and other formats.

use bstr::{join, ByteSlice};
use nom::bytes::complete::{is_not, tag};
use nom::character::complete::{line_ending, multispace0};
use nom::combinator::{all_consuming, map_parser, opt, value};
use nom::error::ErrorKind;
use nom::sequence::terminated;
use nom::{error, Err, Parser};
use umol::Result;

use self::accumulator::MoleculeProperties;
pub use self::atom::{atom_input, atomlike_input};
pub use self::bond::{bond_input, bondlike_input};
pub use self::counts::{counts_input, Counts};
pub use self::properties::{
    atom_alias_input, basic_property_input, legacy_atom_list_input, property_input, PropertyEntries,
};
use super::atom::{Atom, AtomLike};
use super::bond::{Bond, BondLike};
use super::config::CtabParseFlags;
use super::molecule::{Molecule, MoleculeLike};

mod accumulator;
mod atom;
mod bond;
mod context;
mod convert;
mod counts;
mod properties;
mod sgroup;
mod utils;

/// Parse CTAB block (basic parser, optimized for performance, basic molecules only)
///
/// Parses from the counts line through M END, returning a Molecule.
/// This parser is optimized for basic molecules without query features.
/// It will fail if the CTAB contains query atoms, query bonds, or other advanced features.
pub fn basic_ctab_block<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Molecule, Error = error::Error<&'inp [u8]>> + use<'inp, 'fl> {
    move |input: &'inp [u8]| {
        let (remaining, counts) = counts_block(flags).parse(input)?;
        let atom_count = counts.atom_count as usize;
        let bond_count = counts.bond_count as usize;
        let (remaining, atoms) = atom_block(atom_count, flags).parse(remaining)?;
        let (remaining, bonds) = bond_block(bond_count, flags).parse(remaining)?;
        let (remaining, properties) = basic_properties_block(flags).parse(remaining)?;
        let (remaining, _) = end_block().parse(remaining)?;
        let molecule = build_molecule(atoms, bonds, properties)
            .map_err(|_| Err::Error(error::Error::new(remaining, ErrorKind::MapRes)))?;
        Ok((remaining, molecule))
    }
}

/// Parse CTAB block (general parser, handles all features including queries)
///
/// Parses from the counts line through M END, returning a complete Molecule.
/// This parser handles all CTAB features including query atoms, query bonds,
/// R-groups, S-groups, and all property lines.
pub fn ctab_block<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = MoleculeLike, Error = error::Error<&'inp [u8]>> + use<'inp, 'fl>
{
    move |input: &'inp [u8]| {
        let legacy_features = flags.contains(CtabParseFlags::LEGACY_FEATURES);
        let (remaining, counts) = counts_block(flags).parse(input)?;
        let atom_count = counts.atom_count as usize;
        let bond_count = counts.bond_count as usize;
        let (remaining, atoms) = atomlike_block(atom_count, flags).parse(remaining)?;
        let (remaining, bonds) = bondlike_block(bond_count, flags).parse(remaining)?;
        let (remaining, legacy_properties) = if legacy_features {
            legacy_atom_list_block(flags).parse(remaining)?
        } else {
            (remaining, Vec::new())
        };
        let (remaining, properties) = properties_block(flags).parse(remaining)?;
        let (remaining, _) = end_block().parse(remaining)?;
        let properties = properties.into_iter().chain(legacy_properties).collect();
        let molecule = build_moleculelike(atoms, bonds, properties)
            .map_err(|_| Err::Error(error::Error::new(remaining, ErrorKind::MapRes)))?;
        Ok((remaining, molecule))
    }
}

/// Parse counts block
fn counts_block<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Counts, Error = error::Error<&'inp [u8]>> + use<'inp, 'fl> {
    map_parser(terminated(is_not("\r\n"), line_ending), counts_input(flags))
}

/// Parse atom block (basic atoms only)
fn atom_block<'inp, 'fl>(
    atom_count: usize,
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Vec<Atom>, Error = error::Error<&'inp [u8]>> + use<'inp, 'fl>
{
    move |input: &'inp [u8]| {
        let mut atoms = Vec::with_capacity(atom_count);
        let mut lines_iter = input.lines_with_terminator();
        let mut consumed = 0;

        for _ in 0..atom_count {
            let line = lines_iter
                .next()
                .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

            let (_, atom) = all_consuming(atom_input(flags)).parse(line)?;
            atoms.push(atom);
            consumed += line.len();
        }

        let remaining = &input[consumed..];
        Ok((remaining, atoms))
    }
}

/// Parse atom-like block
fn atomlike_block<'inp, 'fl>(
    atom_count: usize,
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Vec<AtomLike>, Error = error::Error<&'inp [u8]>> + use<'inp, 'fl>
{
    move |input: &'inp [u8]| {
        let mut atoms = Vec::with_capacity(atom_count);
        let mut lines_iter = input.lines_with_terminator();
        let mut consumed = 0;

        for _ in 0..atom_count {
            let line = lines_iter
                .next()
                .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

            let (_, atom) = all_consuming(atomlike_input(flags)).parse(line)?;
            atoms.push(atom);
            consumed += line.len();
        }

        let remaining = &input[consumed..];
        Ok((remaining, atoms))
    }
}

/// Parse bond block (basic bonds only)
fn bond_block<'inp, 'fl>(
    bond_count: usize,
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Vec<(usize, usize, Bond)>, Error = error::Error<&'inp [u8]>>
       + use<'inp, 'fl> {
    move |input: &'inp [u8]| {
        let mut bonds = Vec::with_capacity(bond_count);
        let mut lines_iter = input.lines_with_terminator();
        let mut consumed = 0;

        for _ in 0..bond_count {
            let line = lines_iter
                .next()
                .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

            let (_, (atom1, atom2, bond)) = all_consuming(bond_input(flags)).parse(line)?;
            bonds.push((atom1, atom2, bond));
            consumed += line.len();
        }

        let remaining = &input[consumed..];
        Ok((remaining, bonds))
    }
}

/// Parse bond-like block
fn bondlike_block<'inp, 'fl>(
    bond_count: usize,
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Vec<(usize, usize, BondLike)>, Error = error::Error<&'inp [u8]>>
       + use<'inp, 'fl> {
    move |input: &'inp [u8]| {
        let mut bonds = Vec::with_capacity(bond_count);
        let mut lines_iter = input.lines_with_terminator();
        let mut consumed = 0;

        for _bond_idx in 0..bond_count {
            let line = lines_iter
                .next()
                .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

            let (_, (atom1, atom2, bond)) = all_consuming(bondlike_input(flags)).parse(line)?;
            bonds.push((atom1, atom2, bond));
            consumed += line.len();
        }

        let remaining = &input[consumed..];
        Ok((remaining, bonds))
    }
}

// Parse legacy atom list block
fn legacy_atom_list_block<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Vec<PropertyEntries>, Error = error::Error<&'inp [u8]>>
       + use<'inp, 'fl> {
    move |input: &'inp [u8]| {
        let mut properties = Vec::new();
        let lines_iter = input.lines_with_terminator();
        let mut consumed = 0;

        for line in lines_iter {
            if let Ok((_, property)) = all_consuming(legacy_atom_list_input(flags)).parse(line) {
                properties.push(property);
                consumed += line.len();
            } else {
                break; // Backtrack
            }
        }

        let remaining = &input[consumed..];
        Ok((remaining, properties))
    }
}

/// Parse properties block (basic properties only)
fn basic_properties_block<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Vec<PropertyEntries>, Error = error::Error<&'inp [u8]>>
       + use<'inp, 'fl> {
    move |input: &'inp [u8]| {
        let mut properties = Vec::new();
        let mut lines_iter = input.lines_with_terminator();
        let mut consumed = 0;

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
                        consumed += line_bytes;
                    } else {
                        break; // Backtrack
                    };
                } else {
                    break; // Incomplete atom alias
                }
            } else {
                match all_consuming(basic_property_input(flags)).parse(line) {
                    Ok((_, property)) => {
                        properties.push(property);
                        consumed += line.len();
                    }
                    Err(_) => {
                        return Err(Err::Error(error::Error::new(
                            &input[consumed..],
                            ErrorKind::Tag,
                        )));
                    }
                }
            }
        }

        let remaining = &input[consumed..];
        Ok((remaining, properties))
    }
}

/// Parse properties block
fn properties_block<'inp, 'fl>(
    flags: &'fl CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Vec<PropertyEntries>, Error = error::Error<&'inp [u8]>>
       + use<'inp, 'fl> {
    move |input: &'inp [u8]| {
        let mut properties = Vec::new();
        let mut lines_iter = input.lines_with_terminator();
        let mut consumed = 0;

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
                        consumed += line_bytes;
                    } else {
                        break; // Backtrack
                    };
                } else {
                    break; // Incomplete atom alias
                }
            } else if let Ok((_, property)) = all_consuming(property_input(flags)).parse(line) {
                properties.push(property);
                consumed += line.len();
            } else {
                break; // Backtrack
            }
        }

        let remaining = &input[consumed..];
        Ok((remaining, properties))
    }
}

/// Parse end block (M END line)
fn end_block<'inp>() -> impl Parser<&'inp [u8], Output = (), Error = error::Error<&'inp [u8]>> {
    value((), opt(terminated(tag("M  END"), multispace0)))
}

/// Build molecule
fn build_molecule(
    atoms: Vec<Atom>,
    bonds: Vec<(usize, usize, Bond)>,
    properties: Vec<PropertyEntries>,
) -> Result<Molecule> {
    let mut molecule = Molecule::new();

    for atom in atoms {
        molecule.add_atom(atom);
    }

    for (idx1, idx2, bond) in bonds {
        molecule.add_bond(idx1, idx2, bond);
    }

    let mut acc = MoleculeProperties::new();
    for entry in properties {
        acc.add_entry(entry, CtabParseFlags::BASIC)?;
    }
    acc.update_molecule(&mut molecule, CtabParseFlags::BASIC)?;

    Ok(molecule)
}

/// Build molecule-like structure
fn build_moleculelike(
    atoms: Vec<AtomLike>,
    bonds: Vec<(usize, usize, BondLike)>,
    properties: Vec<PropertyEntries>,
) -> Result<MoleculeLike> {
    let mut molecule = MoleculeLike::new();

    for atom in atoms {
        molecule.add_atom(atom);
    }

    for (idx1, idx2, bond) in bonds {
        molecule.add_bond(idx1, idx2, bond);
    }

    let mut acc = MoleculeProperties::new();
    for entry in properties {
        acc.add_entry(entry, CtabParseFlags::LENIENT)?;
    }

    acc.update_moleculelike(&mut molecule, CtabParseFlags::LENIENT)?;

    Ok(molecule)
}

#[cfg(test)]
mod tests;
