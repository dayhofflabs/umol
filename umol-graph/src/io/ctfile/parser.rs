//! CTFile format parsers (CTAB, MOL, SDF).
//!
//! This module provides the main entry points for parsing Connection Table (MDL) formats.

use std::borrow::Cow;

use indexmap::IndexMap;
use nom::character::complete::multispace0;
use nom::sequence::terminated;
use nom::{Err, Parser};

use self::accumulator::PropertyAccumulator;
use self::atom::{atom_block, extended_atom_block};
pub use self::atom::{atom_input, extended_atom_input}; // NOTE: Re-exported for benchmarks
use self::bond::{bond_block, extended_bond_block};
pub use self::bond::{bond_input, extended_bond_input}; // NOTE: Re-exported for benchmarks
use self::counts::counts_block;
pub use self::counts::counts_input; // NOTE: Re-exported for benchmarks
use self::header::header_block;
use self::legacy_atom_list::legacy_atom_list_block;
pub use self::legacy_atom_list::legacy_atom_list_input; // NOTE: Re-exported for benchmarks
use self::properties::{extended_properties_block, properties_block, PropertyEntries};
pub use self::properties::{extended_property_input, property_input}; // NOTE: Re-exported for benchmarks
use self::sdf_data::sdf_data_block;
use super::config::{CtabParseFlags, CtfileIoConfig};
use super::error::{CtfileError, ParseError};
use crate::io::utils::normalize_whitespace;
use crate::ops::config::ChemistryModel;
use crate::ops::resolver::Resolver;
use crate::ops::solution::Solution;
use crate::position::Point3D;
use crate::table_ir::bond::Bond;
use crate::table_ir::source::SourceFormat;
use crate::table_ir::{Atom, AtomSymbol, ExtendedAtom, ExtendedBond, ExtendedMolecule, Molecule};
use umol_ast::ast::{MoleculeAst, TryIntoAst};

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
/// This parser is optimized for basic molecules. For extended molecules, use extended_ctab_block.
pub fn ctab_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Molecule, u32), Error = ParseError> + use<'inp> {
    debug_assert!(
        CtabParseFlags::BASIC_MAX.contains(flags),
        "flags must be a subset of BASIC_MAX"
    );
    let legacy_atom_lists = flags.contains(CtabParseFlags::LEGACY_ATOM_LISTS);

    move |input: &'inp [u8]| {
        let (remaining, (counts, molecule_properties, line_offset)) =
            counts_block(line_offset, flags).parse(input)?;
        let atom_count = counts.atom_count;
        let bond_count = counts.bond_count;
        let atom_list_count = counts.atom_list_count;

        if !legacy_atom_lists && atom_list_count > 0 {
            return Err(Err::Error(ParseError::UnsupportedLegacyAtomList {
                line: line_offset + atom_count + bond_count,
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

        let properties = if !legacy_properties.is_empty() || !molecule_properties.is_empty() {
            properties
                .into_iter()
                .chain(molecule_properties)
                .chain(legacy_properties)
                .collect()
        } else {
            properties
        };

        let molecule =
            build_molecule(atoms, bonds, positions, properties, flags).map_err(Err::Error)?;
        Ok((remaining, (molecule, line_offset)))
    }
}

/// Parse CTAB block (general parser, handles all features including queries)
///
/// This parser handles extended molecules. For basic molecules, use ctab_block.
pub fn extended_ctab_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (ExtendedMolecule, u32), Error = ParseError> + use<'inp> {
    let legacy_atom_lists = flags.contains(CtabParseFlags::LEGACY_ATOM_LISTS);

    move |input: &'inp [u8]| {
        let (remaining, (counts, molecule_properties, line_offset)) =
            counts_block(line_offset, flags).parse(input)?;
        let atom_count = counts.atom_count;
        let bond_count = counts.bond_count;
        let atom_list_count = counts.atom_list_count;

        if !legacy_atom_lists && atom_list_count > 0 {
            return Err(Err::Error(ParseError::UnsupportedLegacyAtomList {
                line: line_offset + atom_count + bond_count,
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

        let properties = if !legacy_properties.is_empty() || !molecule_properties.is_empty() {
            properties
                .into_iter()
                .chain(legacy_properties)
                .chain(molecule_properties)
                .collect()
        } else {
            properties
        };

        let extended = build_extended_molecule(atoms, bonds, positions, properties, flags)
            .map_err(Err::Error)?;
        Ok((remaining, (extended, line_offset)))
    }
}

/// Build Molecule
fn build_molecule(
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    positions: Option<Vec<Point3D>>,
    properties: Vec<PropertyEntries>,
    flags: CtabParseFlags,
) -> Result<Molecule, ParseError> {
    let mut molecule = Molecule {
        atoms,
        bonds,
        rings: Vec::new(),
        positions,
        multicenter_bonds: Vec::new(),
        comments: Vec::new(),
        properties: IndexMap::new(),
        stereo_interpretation: None,
        source_format: SourceFormat::MOL,
    };

    let mut acc = PropertyAccumulator::new();
    for entry in properties {
        acc.add_entry(entry, flags)?;
    }
    acc.update_molecule(&mut molecule, flags)?;

    Ok(molecule)
}

/// Build extended molecule
fn build_extended_molecule(
    atoms: Vec<ExtendedAtom>,
    bonds: Vec<ExtendedBond>,
    positions: Option<Vec<Point3D>>,
    properties: Vec<PropertyEntries>,
    flags: CtabParseFlags,
) -> Result<ExtendedMolecule, ParseError> {
    let mut molecule = ExtendedMolecule {
        atoms,
        bonds,
        rings: Vec::new(),
        positions,
        multicenter_bonds: Vec::new(),
        stereo_interpretation: None,
        comments: Vec::new(),
        properties: IndexMap::new(),
        ctfile_data: None,
        cx_data: None,
        source_format: SourceFormat::MOL,
    };

    let mut acc = PropertyAccumulator::new();
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

/// Parse MOL to a resolved [`MoleculeAst`] using default IO config and
/// [`ChemistryModel::default`].
pub fn parse_mol(input: &str) -> Result<MoleculeAst, CtfileError> {
    parse_mol_bytes(input.as_bytes())
}

/// Parse MOL bytes to a resolved [`MoleculeAst`] using default IO config and
/// [`ChemistryModel::default`].
pub fn parse_mol_bytes(input: &[u8]) -> Result<MoleculeAst, CtfileError> {
    parse_mol_bytes_with(input, &CtfileIoConfig::basic(), &ChemistryModel::default())
}

/// Parse MOL to a resolved [`MoleculeAst`] with explicit IO config and
/// chemistry model.
pub fn parse_mol_with(
    input: &str,
    io_config: &CtfileIoConfig,
    model: &ChemistryModel,
) -> Result<MoleculeAst, CtfileError> {
    parse_mol_bytes_with(input.as_bytes(), io_config, model)
}

/// Parse MOL bytes to a resolved [`MoleculeAst`] with explicit IO config and
/// chemistry model.
pub fn parse_mol_bytes_with(
    input: &[u8],
    io_config: &CtfileIoConfig,
    model: &ChemistryModel,
) -> Result<MoleculeAst, CtfileError> {
    let table_mol = parse_mol_bytes_to_table_ir_with(input, io_config)?;
    let mut ast: MoleculeAst = (&table_mol)
        .try_into_ast(&())
        .expect("table_ir → MoleculeAst lift is currently infallible");
    match Resolver::new(model).resolve(&mut ast)? {
        Solution::Determined(()) => Ok(ast),
        Solution::Underdetermined(()) => Err(CtfileError::ResolveUnderdetermined),
        Solution::Contradictory(c) => Err(CtfileError::ResolveContradictory(c)),
    }
}

/// Parse MOL to [`MoleculeAst`] without running the solver.
pub fn parse_mol_to_ast(input: &str) -> Result<MoleculeAst, CtfileError> {
    parse_mol_bytes_to_ast(input.as_bytes())
}

/// Parse MOL bytes to [`MoleculeAst`] without running the solver.
pub fn parse_mol_bytes_to_ast(input: &[u8]) -> Result<MoleculeAst, CtfileError> {
    let table_mol = parse_mol_bytes_to_table_ir(input)?;
    let ast: MoleculeAst = (&table_mol)
        .try_into_ast(&())
        .expect("table_ir → MoleculeAst lift is currently infallible");
    Ok(ast)
}

/// Parse MOL bytes to `table_ir::Molecule` with options (optimized, basic molecules only)
pub fn parse_mol_bytes_to_table_ir_with(
    input: &[u8],
    config: &CtfileIoConfig,
) -> Result<Molecule, ParseError> {
    let flags = config.parse_flags;

    let data: Cow<'_, [u8]> = if flags.contains(CtabParseFlags::UNICODE) {
        normalize_whitespace(input)
    } else {
        Cow::Borrowed(input)
    };

    let (remaining, (comments, line_offset)) = header_block(0).parse(&*data)?;

    let (_, (mut molecule, _line_offset)) =
        terminated(ctab_block(line_offset, flags), multispace0).parse(remaining)?;

    molecule.comments = comments;

    Ok(molecule)
}

/// Parse MOL bytes to `table_ir::Molecule` with basic flags.
pub fn parse_mol_bytes_to_table_ir(input: &[u8]) -> Result<Molecule, ParseError> {
    let config = CtfileIoConfig::basic();
    parse_mol_bytes_to_table_ir_with(input, &config)
}

/// Parse MOL string to `table_ir::Molecule` with options.
pub fn parse_mol_to_table_ir_with(
    input: &str,
    config: &CtfileIoConfig,
) -> Result<Molecule, ParseError> {
    parse_mol_bytes_to_table_ir_with(input.as_bytes(), config)
}

/// Parse MOL string to `table_ir::Molecule` (optimized, basic molecules only).
pub fn parse_mol_to_table_ir(input: &str) -> Result<Molecule, ParseError> {
    parse_mol_bytes_to_table_ir(input.as_bytes())
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

    let (remaining, (comments, line_offset)) = header_block(0).parse(&*data)?;

    let (_, (mut molecule, _line_offset)) =
        terminated(extended_ctab_block(line_offset, flags), multispace0).parse(remaining)?;
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
) -> Result<(&'inp [u8], (Molecule, u32)), ParseError> {
    let flags = config.parse_flags;

    let (remaining, (comments, line_offset)) = header_block(line_offset).parse(input)?;

    let (remaining, (mut molecule, line_offset)) =
        ctab_block(line_offset, flags).parse(remaining)?;

    let (remaining, (mut sdf_data, line_offset)) = sdf_data_block(line_offset).parse(remaining)?;

    molecule.comments = comments;
    molecule.properties.append(&mut sdf_data);

    Ok((remaining, (molecule, line_offset)))
}

/// Parse single SDF compound into ExtendedMolecule
fn parse_sdf_extended_molecule<'inp>(
    input: &'inp [u8],
    line_offset: u32,
    config: &CtfileIoConfig,
) -> Result<(&'inp [u8], (ExtendedMolecule, u32)), ParseError> {
    let flags = config.parse_flags;

    let (remaining, (comments, line_offset)) = header_block(line_offset).parse(input)?;

    let (remaining, (mut molecule, line_offset)) =
        extended_ctab_block(line_offset, flags).parse(remaining)?;

    let (remaining, (mut sdf_data, line_offset)) = sdf_data_block(line_offset).parse(remaining)?;

    molecule.comments = comments;
    molecule.properties.append(&mut sdf_data);

    Ok((remaining, (molecule, line_offset)))
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
        let (new_remaining, (molecule, new_line_offset)) =
            parse_sdf_molecule(remaining, line_offset, config)?;
        molecules.push(molecule);
        remaining = new_remaining;
        line_offset = new_line_offset;
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
        let (new_remaining, (molecule, new_line_offset)) =
            parse_sdf_extended_molecule(remaining, line_offset, config)?;
        molecules.push(molecule);
        remaining = new_remaining;
        line_offset = new_line_offset;
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

#[cfg(test)]
mod tests;
