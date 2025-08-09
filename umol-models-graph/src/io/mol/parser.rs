//! MOL file parser

use nom::bytes::complete::tag;
use nom::character::complete::line_ending;
use nom::combinator::opt;
use nom::multi::{count, many0};
use nom::sequence::terminated;
use nom::{error, IResult, Parser};

use super::super::ctab::atom::{Atom, AtomStandard, AtomSymbol};
use super::super::ctab::bond::{Bond, BondStandard};
use super::super::ctab::molecule::{Header, Molecule, MoleculeStandard, ParsedMol};

use super::super::ctab::parser::accumulator::MoleculeProperties;
use super::super::ctab::parser::atom::{atom_input, atom_input_standard};
use super::super::ctab::parser::bond::{bond_input, bond_input_standard};
use super::super::ctab::parser::counts::counts_input;
use super::super::ctab::parser::header::header;
use super::super::ctab::parser::properties::{
    legacy_atom_list_input, property_input, property_input_standard, PropertyEntries,
};

/// Parse MOL block (general parser, handles all features)
pub fn mol_block(input: &[u8]) -> IResult<&[u8], ParsedMol, error::Error<&[u8]>> {
    let (remaining, header) = header().parse(input)?;
    let (remaining, counts) = terminated(counts_input(), line_ending).parse(remaining)?;
    let atom_count = counts.atoms() as usize;
    let bond_count = counts.bonds() as usize;
    let (remaining, atoms) = atom_block(atom_count).parse(remaining)?;
    let (remaining, bonds) = bond_block(bond_count).parse(remaining)?;
    let (remaining, legacy_properties) = legacy_atom_list_block(remaining)?;
    let (remaining, properties) = properties_block(remaining)?;
    let properties = properties.into_iter().chain(legacy_properties).collect();
    let (molecule, is_query) = build_molecule(header, atoms, bonds, properties);
    Ok((remaining, ParsedMol::new(molecule, is_query)))
}

/// Parse MOL block (standard parser, optimized for performance, standard molecules only)
pub fn mol_block_standard(input: &[u8]) -> IResult<&[u8], MoleculeStandard, error::Error<&[u8]>> {
    let (remaining, header) = header().parse(input)?;
    let (remaining, counts) = terminated(counts_input(), line_ending).parse(remaining)?;
    let atom_count = counts.atoms() as usize;
    let bond_count = counts.bonds() as usize;
    let (remaining, atoms) = atom_block_standard(atom_count).parse(remaining)?;
    let (remaining, bonds) = bond_block_standard(bond_count).parse(remaining)?;
    let (remaining, properties) = properties_block_standard(remaining)?;
    let molecule = build_molecule_standard(header, atoms, bonds, properties);
    Ok((remaining, molecule))
}

/// Parse atom block
fn atom_block<'a>(
    atom_count: usize,
) -> impl Parser<&'a [u8], Output = Vec<Atom>, Error = error::Error<&'a [u8]>> {
    count(
        |input| {
            let (input, atom) = terminated(atom_input(), line_ending).parse(input)?;
            Ok((input, atom))
        },
        atom_count,
    )
}

/// Parse atom block (standard parser)
fn atom_block_standard<'a>(
    atom_count: usize,
) -> impl Parser<&'a [u8], Output = Vec<AtomStandard>, Error = error::Error<&'a [u8]>> {
    count(
        |input| {
            let (input, atom) = terminated(atom_input_standard(), line_ending).parse(input)?;
            Ok((input, atom))
        },
        atom_count,
    )
}

/// Parse bond block
fn bond_block<'a>(
    bond_count: usize,
) -> impl Parser<&'a [u8], Output = Vec<(usize, usize, Bond)>, Error = error::Error<&'a [u8]>> {
    count(
        |input| {
            let (input, (atom1, atom2, bond)) =
                terminated(bond_input(), line_ending).parse(input)?;
            Ok((input, (atom1, atom2, bond)))
        },
        bond_count,
    )
}

/// Parse bond block (standard parser)
fn bond_block_standard<'a>(
    bond_count: usize,
) -> impl Parser<&'a [u8], Output = Vec<(usize, usize, BondStandard)>, Error = error::Error<&'a [u8]>>
{
    count(
        |input| {
            let (input, (atom1, atom2, bond)) =
                terminated(bond_input_standard(), line_ending).parse(input)?;
            Ok((input, (atom1, atom2, bond)))
        },
        bond_count,
    )
}

/// Parse legacy atom list block
fn legacy_atom_list_block(
    input: &[u8],
) -> IResult<&[u8], Vec<PropertyEntries>, error::Error<&[u8]>> {
    let (input, legacy_properties) =
        many0(terminated(legacy_atom_list_input(), line_ending)).parse(input)?;
    Ok((input, legacy_properties))
}

/// Parse properties block
fn properties_block(input: &[u8]) -> IResult<&[u8], Vec<PropertyEntries>, error::Error<&[u8]>> {
    let (input, properties) = many0(terminated(property_input(), line_ending)).parse(input)?;
    let (input, _) = opt(terminated(tag("M  END"), opt(line_ending))).parse(input)?;
    Ok((input, properties))
}

/// Parse properties block (standard parser)
fn properties_block_standard(
    input: &[u8],
) -> IResult<&[u8], Vec<PropertyEntries>, error::Error<&[u8]>> {
    let (input, properties) =
        many0(terminated(property_input_standard(), line_ending)).parse(input)?;
    let (input, _) = opt(terminated(tag("M  END"), opt(line_ending))).parse(input)?;
    Ok((input, properties))
}

/// Detect if molecule contains query features
/// TODO: Update with remaining properties
fn detect_query_features(
    atoms: &[Atom],
    _bonds: &[(usize, usize, Bond)],
    properties: &[PropertyEntries],
) -> bool {
    // Check atoms for query symbols
    for atom in atoms {
        match &atom.symbol {
            AtomSymbol::Element(_) | AtomSymbol::NamedIsotope(_) => {
                // Standard atoms
            }
            AtomSymbol::AtomList(_)
            | AtomSymbol::Query(_)
            | AtomSymbol::LonePair
            | AtomSymbol::RGroup(_) => {
                return true; // Query atom found
            }
        }

        // Check for query-specific atom properties
        if atom.attachment_point.is_some()
            || atom.attachment_order.is_some()
            || atom.ring_bond_count.is_some()
            || atom.substitution_count.is_some()
            || atom.unsaturated.is_some()
            || atom.link_atom.is_some()
        {
            return true;
        }
    }

    // Check properties for query-specific entries
    for property in properties {
        // TODO: Implement query property detection in PropertyEntries
        match property {
            PropertyEntries::AtomListEntry(_)
            | PropertyEntries::AttachmentPointEntries(_)
            | PropertyEntries::AtomAttachmentOrderEntry(_)
            | PropertyEntries::RingBondCountEntries(_)
            | PropertyEntries::SubstitutionCountEntries(_)
            | PropertyEntries::UnsaturatedAtomEntries(_)
            | PropertyEntries::LinkAtomEntries(_) => {
                return true; // Query property found
            }
            _ => {
                // Standard properties
            }
        }
    }

    // TODO: Check bonds for query features when implemented

    false
}

/// Build molecule from parsed components
fn build_molecule(
    header: Header,
    atoms: Vec<Atom>,
    bonds: Vec<(usize, usize, Bond)>,
    properties: Vec<PropertyEntries>,
) -> (Molecule, bool) {
    // Detect query features
    let is_query = detect_query_features(&atoms, &bonds, &properties);

    // Create molecule
    let mut molecule = Molecule::new();
    molecule.header = header;

    // Add atoms to molecule
    for atom in atoms {
        molecule.add_atom(atom);
    }

    // Add bonds to molecule
    for (idx1, idx2, bond) in bonds {
        molecule.add_bond(idx1, idx2, bond);
    }

    // Create a new MoleculeProperties accumulator
    let mut acc = MoleculeProperties::new();
    for entry in properties {
        if let Err(e) = acc.add_entry(entry) {
            eprintln!("Warning: Failed to add property entry: {}", e);
        }
    }

    // Apply the accumulated properties to the molecule
    if let Err(e) = acc.apply(&mut molecule) {
        eprintln!("Warning: Failed to apply properties: {}", e);
    }

    (molecule, is_query)
}

/// Build standard molecule from parsed standard components
fn build_molecule_standard(
    header: Header,
    atoms: Vec<AtomStandard>,
    bonds: Vec<(usize, usize, BondStandard)>,
    properties: Vec<PropertyEntries>,
) -> MoleculeStandard {
    // Create molecule
    let mut molecule = MoleculeStandard::new();
    molecule.header = header;

    // Add atoms to molecule
    for atom in atoms {
        molecule.add_atom(atom);
    }

    // Add bonds to molecule
    for (idx1, idx2, bond) in bonds {
        molecule.add_bond(idx1, idx2, bond);
    }

    // Create a new MoleculeProperties accumulator
    let mut acc = MoleculeProperties::new();
    for entry in properties {
        if let Err(e) = acc.add_entry(entry) {
            eprintln!("Warning: Failed to add property entry: {}", e);
        }
    }

    // Apply the accumulated properties to the molecule
    if let Err(e) = acc.apply_standard(&mut molecule) {
        eprintln!("Warning: Failed to apply properties: {}", e);
    }

    molecule
}
