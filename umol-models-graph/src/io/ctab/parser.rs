//! CTab format parser.
//!
//! This module provides the main entry points for parsing Connection Table (CTAB) format,
//! which is the core molecular structure representation used in MOL, SDF, and other formats.

use bstr::{join, ByteSlice};
use nom::bytes::complete::tag;
use nom::character::complete::line_ending;
use nom::combinator::{opt, value};
use nom::sequence::terminated;
use nom::{error, Err, Parser};

// mod accumulator;
mod atom;
mod bond;
// mod context;
mod convert;
mod counts;
// mod properties;
// mod sgroup;
mod utils;

// use self::accumulator::MoleculeProperties;
// pub use self::atom::{atom_input, atomlike_input};
// pub use self::bond::{bond_input, bondlike_input};
// pub use self::counts::{counts_input, Counts};
// pub use self::properties::{
//     atom_alias_entry, basic_property_input, legacy_atom_list_input, property_input, PropertyEntries,
// };

// use super::atom::{Atom, AtomLike};
// use super::bond::{Bond, BondLike};
// use super::molecule::{Molecule, MoleculeLike};

// /// Parse CTAB block (standard parser, optimized for performance, standard molecules only)
// ///
// /// Parses from the counts line through M END, returning a MoleculeStandard.
// /// This parser is optimized for standard molecules without query features.
// /// It will fail if the CTAB contains query atoms, query bonds, or other non-standard features.
// pub fn basic_ctab_block<'a>(
// ) -> impl Parser<&'a [u8], Output = Molecule, Error = error::Error<&'a [u8]>> + 'a {
//     move |input: &'a [u8]| {
//         let (remaining, counts) = terminated(counts_input(), line_ending).parse(input)?;
//         let atom_count = counts.atoms() as usize;
//         let bond_count = counts.bonds() as usize;
//         let (remaining, atoms) = atom_block(atom_count).parse(remaining)?;
//         let (remaining, bonds) = bond_block(bond_count).parse(remaining)?;
//         let (remaining, properties) = basic_properties_block().parse(remaining)?;
//         let (remaining, _) = end_block().parse(remaining)?;

//         let molecule = build_molecule(atoms, bonds, properties);
//         Ok((remaining, molecule))
//     }
// }

// /// Parse CTAB block (general parser, handles all features including queries)
// ///
// /// Parses from the counts line through M END, returning a complete Molecule.
// /// This parser handles all CTAB features including query atoms, query bonds,
// /// R-groups, S-groups, and all property lines.
// pub fn ctab_block<'a>(
// ) -> impl Parser<&'a [u8], Output = MoleculeLike, Error = error::Error<&'a [u8]>> {
//     move |input: &'a [u8]| {
//         let (remaining, counts) = terminated(counts_input(), line_ending).parse(input)?;
//         let atom_count = counts.atoms() as usize;
//         let bond_count = counts.bonds() as usize;
//         let (remaining, atoms) = atomlike_block(atom_count).parse(remaining)?;
//         let (remaining, bonds) = bondlike_block(bond_count).parse(remaining)?;
//         let (remaining, legacy_properties) = legacy_atom_list_block().parse(remaining)?;
//         let (remaining, properties) = properties_block().parse(remaining)?;
//         let (remaining, _) = end_block().parse(remaining)?;

//         let properties = properties.into_iter().chain(legacy_properties).collect();
//         let molecule = build_moleculelike(atoms, bonds, properties);
//         Ok((remaining, molecule))
//     }
// }

// /// Parse atom block (basic atoms only)
// fn atom_block<'a>(
//     atom_count: usize,
//     allow_unicode: bool,
//     allow_named_isotopes: bool,
//     strict_padding: bool,
// ) -> impl Parser<&'a [u8], Output = Vec<Atom>, Error = error::Error<&'a [u8]>> {
//     move |input: &'a [u8]| {
//         let mut atoms = Vec::with_capacity(atom_count);
//         let mut lines_iter = input.lines_with_terminator();
//         let mut consumed = 0;

//         for _atom_idx in 0..atom_count {
//             let line = lines_iter
//                 .next()
//                 .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

//             let (_, atom) = atom_input(true, true, true).parse(line)?;
//             atoms.push(atom);
//             consumed += line.len();
//         }

//         let remaining = &input[consumed..];
//         Ok((remaining, atoms))
//     }
// }

// /// Parse atom-like block
// fn atomlike_block<'a>(
//     atom_count: usize,
//     allow_unicode: bool,
//     allow_named_isotopes: bool,
//     allow_rgroups: bool,
//     allow_queries: bool,
//     allow_extended_queries: bool,
//     allow_electrons: bool,
//     strict_padding: bool,
// ) -> impl Parser<&'a [u8], Output = Vec<AtomLike>, Error = error::Error<&'a [u8]>> {
//     move |input: &'a [u8]| {
//         let mut atoms = Vec::with_capacity(atom_count);
//         let mut lines_iter = input.lines_with_terminator();
//         let mut consumed = 0;

//         for _ in 0..atom_count {
//             let line = lines_iter
//                 .next()
//                 .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

//             let (_, atom) = atomlike_input(
//                 allow_unicode,
//                 allow_named_isotopes,
//                 allow_rgroups,
//                 allow_queries,
//                 allow_extended_queries,
//                 allow_electrons,
//                 strict_padding,
//             )
//             .parse(line)?;
//             atoms.push(atom);
//             consumed += line.len();
//         }

//         let remaining = &input[consumed..];
//         Ok((remaining, atoms))
//     }
// }

// /// Parse bond block (basic bonds only)
// fn bond_block<'a>(
//     bond_count: usize,
// ) -> impl Parser<&'a [u8], Output = Vec<(usize, usize, Bond)>, Error = error::Error<&'a [u8]>> {
//     move |input: &'a [u8]| {
//         let mut bonds = Vec::with_capacity(bond_count);
//         let mut lines_iter = input.lines_with_terminator();
//         let mut consumed = 0;

//         for _bond_idx in 0..bond_count {
//             let line = lines_iter
//                 .next()
//                 .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

//             let content = line.trim_end_with(|c| c == '\r' || c == '\n');
//             let (_, (atom1, atom2, bond)) = bond_input().parse(content)?;
//             bonds.push((atom1, atom2, bond));
//             consumed += line.len();
//         }

//         let remaining = &input[consumed..];
//         Ok((remaining, bonds))
//     }
// }

// /// Parse bond-like block
// fn bondlike_block<'a>(
//     bond_count: usize,
// ) -> impl Parser<&'a [u8], Output = Vec<(usize, usize, BondLike)>, Error = error::Error<&'a [u8]>> {
//     move |input: &'a [u8]| {
//         let mut bonds = Vec::with_capacity(bond_count);
//         let mut lines_iter = input.lines_with_terminator();
//         let mut consumed = 0;

//         for _bond_idx in 0..bond_count {
//             let line = lines_iter
//                 .next()
//                 .ok_or_else(|| Err::Error(error::Error::new(input, error::ErrorKind::Eof)))?;

//             let content = line.trim_end_with(|c| c == '\r' || c == '\n');
//             let (_, (atom1, atom2, bond)) = bondlike_input().parse(content)?;
//             bonds.push((atom1, atom2, bond));
//             consumed += line.len();
//         }

//         let remaining = &input[consumed..];
//         Ok((remaining, bonds))
//     }
// }

// // Parse legacy atom list block
// fn legacy_atom_list_block<'a>(
// ) -> impl Parser<&'a [u8], Output = Vec<PropertyEntries>, Error = error::Error<&'a [u8]>> {
//     move |input: &'a [u8]| {
//         let mut properties = Vec::new();
//         let mut lines_iter = input.lines_with_terminator();
//         let mut consumed = 0;

//         while let Some(line) = lines_iter.next() {
//             let content = line.trim_end_with(|c| c == '\r' || c == '\n');

//             match legacy_atom_list_input().parse(content) {
//                 Ok((_, property)) => {
//                     properties.push(property);
//                     consumed += line.len();
//                 }
//                 Err(_) => break, // Backtrack
//             }
//         }

//         let remaining = &input[consumed..];
//         Ok((remaining, properties))
//     }
// }

// /// Parse properties block (basic properties only)
// fn basic_properties_block<'a>(
// ) -> impl Parser<&'a [u8], Output = Vec<PropertyEntries>, Error = error::Error<&'a [u8]>> + 'a {
//     move |input: &'a [u8]| {
//         let mut properties = Vec::new();
//         let mut lines_iter = input.lines_with_terminator();
//         let mut consumed = 0;

//         while let Some(line) = lines_iter.next() {
//             let content = line.trim_end_with(|c| c == '\r' || c == '\n');

//             if content.starts_with(b"M  END") {
//                 // Include M  END line in remaining for termination
//                 break;
//             }

//             // Handle atom alias (two-line property)
//             if content.starts_with(b"A  ") {
//                 if let Some(next_line) = lines_iter.next() {
//                     let next_content = next_line.trim_end_with(|c| c == '\r' || c == '\n');
//                     let input_for_parser = join(b"\n", &[content, next_content]);
//                     let line_bytes = line.len() + next_line.len();

//                     // TODO: Replace with simple atom alias parser
//                     match basic_property_input().parse(&input_for_parser) {
//                         Ok((_, property)) => {
//                             properties.push(property);
//                             consumed += line_bytes;
//                         }
//                         Err(_) => break, // Backtrack
//                     };
//                 } else {
//                     break; // Incomplete atom alias
//                 }
//             } else {
//                 match basic_property_input().parse(content) {
//                     Ok((_, property)) => {
//                         properties.push(property);
//                         consumed += line.len();
//                     }
//                     Err(_) => break, // Backtrack
//                 }
//             }
//         }

//         let remaining = &input[consumed..];
//         Ok((remaining, properties))
//     }
// }

// /// Parse properties block
// fn properties_block<'a>(
// ) -> impl Parser<&'a [u8], Output = Vec<PropertyEntries>, Error = error::Error<&'a [u8]>> + 'a {
//     move |input: &'a [u8]| {
//         let mut properties = Vec::new();
//         let mut lines_iter = input.lines_with_terminator();
//         let mut consumed = 0;

//         while let Some(line) = lines_iter.next() {
//             let content = line.trim_end_with(|c| c == '\r' || c == '\n');

//             if content.starts_with(b"M  END") {
//                 break; // Include M  END line in remaining for termination
//             }

//             // Handle atom alias (two-line property)
//             if content.starts_with(b"A  ") {
//                 if let Some(next_line) = lines_iter.next() {
//                     let next_content = next_line.trim_end_with(|c| c == '\r' || c == '\n');
//                     let input_for_parser = join(b"\n", &[content, next_content]);
//                     let line_bytes = line.len() + next_line.len();

//                     // TODO: Replace with simple atom alias parser
//                     match property_input().parse(&input_for_parser) {
//                         Ok((_, property)) => {
//                             properties.push(property);
//                             consumed += line_bytes;
//                         }
//                         Err(_) => break, // Backtrack
//                     };
//                 } else {
//                     break; // Incomplete atom alias
//                 }
//             } else {
//                 match property_input().parse(content) {
//                     Ok((_, property)) => {
//                         properties.push(property);
//                         consumed += line.len();
//                     }
//                     Err(_) => break, // Backtrack
//                 }
//             }
//         }

//         let remaining = &input[consumed..];
//         Ok((remaining, properties))
//     }
// }

// /// Parse end block (M END line)
// fn end_block<'a>() -> impl Parser<&'a [u8], Output = (), Error = error::Error<&'a [u8]>> {
//     value((), opt(terminated(tag("M  END"), opt(line_ending))))
// }

// /// Build molecule
// fn build_molecule(
//     atoms: Vec<Atom>,
//     bonds: Vec<(usize, usize, Bond)>,
//     properties: Vec<PropertyEntries>,
// ) -> Molecule {
//     let mut molecule = Molecule::new();

//     for atom in atoms {
//         molecule.add_atom(atom);
//     }

//     for (idx1, idx2, bond) in bonds {
//         molecule.add_bond(idx1, idx2, bond);
//     }

//     let mut acc = MoleculeProperties::new();
//     for entry in properties {
//         if let Err(e) = acc.add_entry(entry) {
//             eprintln!("Warning: Failed to add property entry: {}", e);
//         }
//     }
//     if let Err(e) = acc.update_molecule(&mut molecule) {
//         eprintln!("Warning: Failed to update molecule: {}", e);
//     }

//     molecule
// }

// /// Build molecule-like structure
// fn build_moleculelike(
//     atoms: Vec<AtomLike>,
//     bonds: Vec<(usize, usize, BondLike)>,
//     properties: Vec<PropertyEntries>,
// ) -> MoleculeLike {
//     let mut molecule = MoleculeLike::new();

//     for atom in atoms {
//         molecule.add_atom(atom);
//     }

//     for (idx1, idx2, bond) in bonds {
//         molecule.add_bond(idx1, idx2, bond);
//     }

//     let mut acc = MoleculeProperties::new();
//     for entry in properties {
//         if let Err(e) = acc.add_entry(entry) {
//             eprintln!("Warning: Failed to add property entry: {}", e);
//         }
//     }

//     if let Err(e) = acc.update_moleculelike(&mut molecule) {
//         eprintln!("Warning: Failed to update molecule-like structure: {}", e);
//     }

//     molecule
// }

// #[cfg(test)]
// mod tests;