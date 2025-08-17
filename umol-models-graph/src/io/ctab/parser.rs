//! CTab format parser.
//!
//! This module provides the main entry points for parsing Connection Table (CTAB) format,
//! which is the core molecular structure representation used in MOL, SDF, and other formats.

mod accumulator;
mod atom;
mod bond;
mod context;
mod convert;
mod counts;
mod properties;
mod sgroup;
mod utils;

pub use self::atom::{atom_input, atom_input_standard};
pub use self::bond::{bond_input, bond_input_standard};
pub use self::counts::counts_input;
pub use self::properties::{legacy_atom_list_input, property_input, property_input_standard};

use bstr::{join, ByteSlice};
use nom::character::complete::line_ending;
use nom::sequence::terminated;
use nom::{error, Err, Parser};

use self::accumulator::MoleculeProperties;
use self::properties::PropertyEntries;
use self::utils::remaining_input;
use super::atom::{Atom, AtomStandard};
use super::bond::{Bond, BondStandard};
use super::molecule::{Molecule, MoleculeStandard};

/// Parse CTAB block (general parser, handles all features including queries)
///
/// Parses from the counts line through M END, returning a complete Molecule.
/// This parser handles all CTAB features including query atoms, query bonds,
/// R-groups, S-groups, and all property lines.
pub fn ctab_block<'a>() -> impl Parser<&'a [u8], Output = Molecule, Error = error::Error<&'a [u8]>>
{
    move |input: &'a [u8]| {
        let (remaining, counts) = terminated(counts_input(), line_ending).parse(input)?;
        let atom_count = counts.atoms() as usize;
        let bond_count = counts.bonds() as usize;
        let (remaining, atoms) = atom_block(atom_count).parse(remaining)?;
        let (remaining, bonds) = bond_block(bond_count).parse(remaining)?;
        let (remaining, legacy_properties) = legacy_atom_list_block().parse(remaining)?;
        let (remaining, properties) = properties_block().parse(remaining)?;
        let properties = properties.into_iter().chain(legacy_properties).collect();
        let molecule = build_molecule(atoms, bonds, properties);
        Ok((remaining, molecule))
    }
}

/// Parse CTAB block (standard parser, optimized for performance, standard molecules only)
///
/// Parses from the counts line through M END, returning a MoleculeStandard.
/// This parser is optimized for standard molecules without query features.
/// It will fail if the CTAB contains query atoms, query bonds, or other non-standard features.
pub fn ctab_block_standard<'a>(
) -> impl Parser<&'a [u8], Output = MoleculeStandard, Error = error::Error<&'a [u8]>> + 'a {
    move |input: &'a [u8]| {
        let (remaining, counts) = terminated(counts_input(), line_ending).parse(input)?;
        let atom_count = counts.atoms() as usize;
        let bond_count = counts.bonds() as usize;
        let (remaining, atoms) = atom_block_standard(atom_count).parse(remaining)?;
        let (remaining, bonds) = bond_block_standard(bond_count).parse(remaining)?;
        let (remaining, properties) = properties_block_standard().parse(remaining)?;
        let molecule = build_molecule_standard(atoms, bonds, properties);
        Ok((remaining, molecule))
    }
}

/// Parse atom block
fn atom_block<'a>(
    atom_count: usize,
) -> impl Parser<&'a [u8], Output = Vec<Atom>, Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let mut atoms = Vec::with_capacity(atom_count);
        let mut lines_iter = input.lines();

        for _ in 0..atom_count {
            let line = lines_iter
                .next()
                .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

            let (_, atom) = atom_input().parse(line)?;
            atoms.push(atom);
        }

        // Calculate remaining input safely
        let remaining = remaining_input(input, lines_iter.as_bytes().as_ptr());

        Ok((remaining, atoms))
    }
}

/// Parse atom block (standard parser)
fn atom_block_standard<'a>(
    atom_count: usize,
) -> impl Parser<&'a [u8], Output = Vec<AtomStandard>, Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let mut atoms = Vec::with_capacity(atom_count);
        let mut lines_iter = input.lines();

        for _ in 0..atom_count {
            let line = lines_iter
                .next()
                .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

            let (_, atom) = atom_input_standard().parse(line)?;
            atoms.push(atom);
        }

        // Calculate remaining input safely
        let remaining = remaining_input(input, lines_iter.as_bytes().as_ptr());

        Ok((remaining, atoms))
    }
}

/// Parse bond block
fn bond_block<'a>(
    bond_count: usize,
) -> impl Parser<&'a [u8], Output = Vec<(usize, usize, Bond)>, Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let mut bonds = Vec::with_capacity(bond_count);
        let mut lines_iter = input.lines();

        for _ in 0..bond_count {
            let line = lines_iter
                .next()
                .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

            // Parse the bond from the line
            let (_, (atom1, atom2, bond)) = bond_input().parse(line)?;
            bonds.push((atom1, atom2, bond));
        }

        // Calculate remaining input safely
        let remaining = remaining_input(input, lines_iter.as_bytes().as_ptr());

        Ok((remaining, bonds))
    }
}

/// Parse bond block (standard parser)
fn bond_block_standard<'a>(
    bond_count: usize,
) -> impl Parser<&'a [u8], Output = Vec<(usize, usize, BondStandard)>, Error = error::Error<&'a [u8]>>
{
    move |input: &'a [u8]| {
        let mut bonds = Vec::with_capacity(bond_count);
        let mut lines_iter = input.lines();

        for _ in 0..bond_count {
            let line = lines_iter
                .next()
                .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

            let (_, (atom1, atom2, bond)) = bond_input_standard().parse(line)?;
            bonds.push((atom1, atom2, bond));
        }

        // Calculate remaining input safely
        let remaining = remaining_input(input, lines_iter.as_bytes().as_ptr());

        Ok((remaining, bonds))
    }
}

/// Parse legacy atom list block
fn legacy_atom_list_block<'a>(
) -> impl Parser<&'a [u8], Output = Vec<PropertyEntries>, Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let mut properties = Vec::new();
        let mut lines_iter = input.lines();

        while let Some(line) = lines_iter.next() {
            match legacy_atom_list_input().parse(line) {
                Ok((_, property)) => {
                    properties.push(property);
                }
                Err(_) => {
                    // Backtrack safely
                    let remaining = remaining_input(input, line.as_ptr());
                    return Ok((remaining, properties));
                }
            }
        }

        // Calculate remaining input safely
        let remaining_bytes = lines_iter.as_bytes();
        let input_start = input.as_ptr() as usize;
        let remaining_start = remaining_bytes.as_ptr() as usize;

        // Ensure we don't underflow when calculating offset
        if remaining_start >= input_start {
            let len = remaining_start - input_start;
            if len <= input.len() {
                let remaining = &input[len..];
                Ok((remaining, properties))
            } else {
                // Remaining pointer is beyond input, return empty slice
                Ok((&input[input.len()..], properties))
            }
        } else {
            // Remaining pointer is before input start, something went wrong
            // Return empty slice to avoid underflow
            Ok((&input[input.len()..], properties))
        }
    }
}

/// Parse properties block
fn properties_block<'a>(
) -> impl Parser<&'a [u8], Output = Vec<PropertyEntries>, Error = error::Error<&'a [u8]>> + 'a {
    move |input: &'a [u8]| {
        let mut properties = Vec::new();
        let mut lines_iter = input.lines();

        while let Some(line) = lines_iter.next() {
            if line.starts_with(b"M  END") {
                break;
            }

            // Handle atom alias (two-line property)
            let (input_line, combined_lines): (&[u8], Vec<u8>);
            if line.starts_with(b"A  ") {
                if let Some(next_line) = lines_iter.next() {
                    combined_lines = join(b"\n", &[line, next_line]);
                    input_line = &combined_lines;
                } else {
                    let remaining = remaining_input(input, line.as_ptr());
                    return Ok((remaining, properties));
                }
            } else {
                input_line = line;
            }

            match property_input().parse(input_line) {
                Ok((_, property)) => {
                    properties.push(property);
                }
                Err(_) => {
                    let remaining = remaining_input(input, line.as_ptr());
                    return Ok((remaining, properties));
                }
            };
        }

        let remaining = remaining_input(input, lines_iter.as_bytes().as_ptr());
        Ok((remaining, properties))
    }
}

/// Parse properties block (standard parser)
fn properties_block_standard<'a>(
) -> impl Parser<&'a [u8], Output = Vec<PropertyEntries>, Error = error::Error<&'a [u8]>> + 'a {
    move |input: &'a [u8]| {
        let mut properties = Vec::new();
        let mut lines_iter = input.lines();

        while let Some(line) = lines_iter.next() {
            if line.starts_with(b"M  END") {
                break;
            }

            // Handle atom alias (two-line property)
            let (input_line, combined_lines): (&[u8], Vec<u8>);
            if line.starts_with(b"A  ") {
                if let Some(next_line) = lines_iter.next() {
                    combined_lines = join(b"\n", &[line, next_line]);
                    input_line = &combined_lines;
                } else {
                    let remaining = remaining_input(input, line.as_ptr());
                    return Ok((remaining, properties));
                }
            } else {
                input_line = line;
            }

            match property_input_standard().parse(input_line) {
                Ok((_, property)) => {
                    properties.push(property);
                }
                Err(_) => {
                    let remaining = remaining_input(input, line.as_ptr());
                    return Ok((remaining, properties));
                }
            };
        }

        let remaining = remaining_input(input, lines_iter.as_bytes().as_ptr());
        Ok((remaining, properties))
    }
}

/// Build molecule from parsed components
fn build_molecule(
    atoms: Vec<Atom>,
    bonds: Vec<(usize, usize, Bond)>,
    properties: Vec<PropertyEntries>,
) -> Molecule {
    let mut molecule = Molecule::new();

    for atom in atoms {
        molecule.add_atom(atom);
    }

    for (idx1, idx2, bond) in bonds {
        molecule.add_bond(idx1, idx2, bond);
    }

    let mut acc = MoleculeProperties::new();
    for entry in properties {
        if let Err(e) = acc.add_entry(entry) {
            eprintln!("Warning: Failed to add property entry: {}", e);
        }
    }

    if let Err(e) = acc.apply(&mut molecule) {
        eprintln!("Warning: Failed to apply properties: {}", e);
    }

    molecule
}

/// Build standard molecule from parsed standard components
fn build_molecule_standard(
    atoms: Vec<AtomStandard>,
    bonds: Vec<(usize, usize, BondStandard)>,
    properties: Vec<PropertyEntries>,
) -> MoleculeStandard {
    let mut molecule = MoleculeStandard::new();

    for atom in atoms {
        molecule.add_atom(atom);
    }

    for (idx1, idx2, bond) in bonds {
        molecule.add_bond(idx1, idx2, bond);
    }

    let mut acc = MoleculeProperties::new();
    for entry in properties {
        if let Err(e) = acc.add_entry(entry) {
            eprintln!("Warning: Failed to add property entry: {}", e);
        }
    }
    if let Err(e) = acc.apply_standard(&mut molecule) {
        eprintln!("Warning: Failed to apply properties: {}", e);
    }

    molecule
}

#[cfg(test)]
mod tests;
