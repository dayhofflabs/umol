//! MOL file parser

use nom::bytes::complete::tag;
use nom::character::complete::line_ending;
use nom::combinator::opt;
use nom::multi::{count, many0};
use nom::sequence::terminated;
use nom::{error, IResult, Parser};

use super::super::ctab::atom::{Atom, AtomStandard, AtomSymbol};
use super::super::ctab::bond::{Bond, BondStandard};
use super::super::ctab::conformer::{Conformer, Point3D};
use super::super::ctab::molecule::{AtomIndex, Header, Molecule, MoleculeStandard, ParsedMol};

use super::super::ctab::parser::apply::Apply;
use super::super::ctab::parser::atom::{atom_input, atom_input_standard};
use super::super::ctab::parser::bond::{bond_input, bond_input_standard};
use super::super::ctab::parser::counts::counts_input;
use super::super::ctab::parser::header::header;
use super::super::ctab::parser::properties::{property_input, property_input_standard, PropertyEntries};

/// Parse MOL block (general parser, handles all features)
pub fn mol_block<'a>(input: &'a [u8]) -> IResult<&'a [u8], ParsedMol, error::Error<&'a [u8]>> {
    let (remaining, header) = header().parse(input)?;
    let (remaining, counts) = terminated(counts_input(), line_ending).parse(remaining)?;
    let atom_count = counts.atoms() as usize;
    let bond_count = counts.bonds() as usize;
    let (remaining, atoms_and_coords) = atom_block(atom_count).parse(remaining)?;
    let (remaining, bonds) = bond_block(bond_count).parse(remaining)?;
    let (remaining, properties) = properties_block(remaining)?;
    let (molecule, is_query) = build_molecule(header, atoms_and_coords, bonds, properties);
    Ok((remaining, ParsedMol::new(molecule, is_query)))
}

/// Parse MOL block (standard parser, optimized for performance, standard molecules only)
pub fn mol_block_standard<'a>(input: &'a [u8]) -> IResult<&'a [u8], MoleculeStandard, error::Error<&'a [u8]>> {
    let (remaining, header) = header().parse(input)?;
    let (remaining, counts) = terminated(counts_input(), line_ending).parse(remaining)?;
    let atom_count = counts.atoms() as usize;
    let bond_count = counts.bonds() as usize;
    let (remaining, atoms_and_coords) = atom_block_standard(atom_count).parse(remaining)?;
    let (remaining, bonds) = bond_block_standard(bond_count).parse(remaining)?;
    let (remaining, properties) = properties_block_standard(remaining)?;
    let molecule = build_molecule_standard(header, atoms_and_coords, bonds, properties);
    Ok((remaining, molecule))
}

/// Parse atom block
fn atom_block<'a>(
    atom_count: usize,
) -> impl Parser<&'a [u8], Output = Vec<(Atom, Point3D)>, Error = error::Error<&'a [u8]>> {
    count(
        |input| {
            let (input, (atom, coord)) = terminated(atom_input(), line_ending).parse(input)?;
            Ok((input, (atom, coord)))
        },
        atom_count,
    )
}

/// Parse atom block (standard parser)
fn atom_block_standard<'a>(
    atom_count: usize,
) -> impl Parser<&'a [u8], Output = Vec<(AtomStandard, Point3D)>, Error = error::Error<&'a [u8]>> {
    count(
        |input| {
            let (input, (atom, coord)) = terminated(atom_input_standard(), line_ending).parse(input)?;
            Ok((input, (atom, coord)))
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
) -> impl Parser<&'a [u8], Output = Vec<(usize, usize, BondStandard)>, Error = error::Error<&'a [u8]>> {
    count(
        |input| {
            let (input, (atom1, atom2, bond)) =
                terminated(bond_input_standard(), line_ending).parse(input)?;
            Ok((input, (atom1, atom2, bond)))
        },
        bond_count,
    )
}

/// Parse properties block (until M  END or end of input)
fn properties_block<'a>(
    input: &'a [u8],
) -> IResult<&'a [u8], Vec<PropertyEntries>, error::Error<&'a [u8]>> {
    let (input, properties) = many0(terminated(property_input, line_ending)).parse(input)?;
    let (input, _) = opt(terminated(tag("M  END"), opt(line_ending))).parse(input)?;
    Ok((input, properties))
}

/// Parse properties block (standard parser)
fn properties_block_standard<'a>(
    input: &'a [u8],
) -> IResult<&'a [u8], Vec<PropertyEntries>, error::Error<&'a [u8]>> {
    let mut properties = Vec::new();
    let mut remaining = input;
    
    loop {
        // Check if we've reached M  END or end of input
        if remaining.is_empty() || remaining.starts_with(b"M  END") {
            break;
        }
        
        // Parse a property line
        match terminated(property_input_standard, line_ending).parse(remaining) {
            Ok((new_remaining, property)) => {
                properties.push(property);
                remaining = new_remaining;
            }
            Err(_) => {
                // If we can't parse as a property, we might have hit M  END or EOF
                break;
            }
        }
    }
    
    // Consume M  END if present
    let (remaining, _) = opt(terminated(tag("M  END"), opt(line_ending))).parse(remaining)?;
    
    Ok((remaining, properties))
}

/// Detect if molecule contains query features
fn detect_query_features(
    atoms: &[(Atom, Point3D)],
    _bonds: &[(usize, usize, Bond)],
    properties: &[PropertyEntries],
) -> bool {
    // Check atoms for query symbols
    for (atom, _) in atoms {
        match &atom.symbol {
            AtomSymbol::Element(_) | AtomSymbol::NamedIsotope(_) => {
                // Standard atoms
            }
            AtomSymbol::AtomList(_)
            | AtomSymbol::Unspecified(_)
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
    atoms_and_coords: Vec<(Atom, Point3D)>,
    bonds: Vec<(usize, usize, Bond)>,
    properties: Vec<PropertyEntries>,
) -> (Molecule, bool) {
    // Detect query features
    let is_query = detect_query_features(&atoms_and_coords, &bonds, &properties);

    // Create molecule
    let mut molecule = Molecule::new();
    molecule.header = header;

    // Extract coordinates for conformer
    let coordinates: Vec<Point3D> = atoms_and_coords.iter().map(|(_, coord)| *coord).collect();

    // Add atoms to molecule
    for (atom, _) in atoms_and_coords {
        molecule.add_atom(atom);
    }

    // Add bonds to molecule
    for (idx1, idx2, bond) in bonds {
        molecule.add_bond(idx1, idx2, bond);
    }

    // Apply properties
    for property in properties {
        if let Err(e) = property.apply(&mut molecule) {
            // For now, ignore property application errors and continue
            // In the future, we might want to collect warnings
            eprintln!("Warning: Failed to apply property: {}", e);
        }
    }

    // Add conformer if we have coordinates
    if !coordinates.is_empty() {
        let conformer = Conformer::from_positions(coordinates);
        if let Err(e) = molecule.add_conformer(conformer) {
            // This shouldn't happen if we parsed correctly, but handle gracefully
            eprintln!("Warning: Failed to add conformer: {}", e);
        }
    }

    (molecule, is_query)
}

/// Build standard molecule from parsed standard components
fn build_molecule_standard(
    header: Header,
    atoms_and_coords: Vec<(AtomStandard, Point3D)>,
    bonds: Vec<(usize, usize, BondStandard)>,
    properties: Vec<PropertyEntries>,
) -> MoleculeStandard {
    // Create molecule
    let mut molecule = MoleculeStandard::new();
    molecule.header = header;

    // Extract coordinates for conformer
    let coordinates: Vec<Point3D> = atoms_and_coords.iter().map(|(_, coord)| *coord).collect();

    // Add atoms to molecule
    for (atom, _) in atoms_and_coords {
        molecule.add_atom(atom);
    }

    // Add bonds to molecule (convert BondStandard to Bond)
    for (idx1, idx2, bond_standard) in bonds {
        // Convert BondStandard to Bond by creating a new Bond with the same type
        let bond = Bond::new(bond_standard.bond_type);
        molecule.graph.add_edge(AtomIndex::new(idx1), AtomIndex::new(idx2), bond);
    }

    // Apply properties - for now, ignore properties in standard parser
    // Standard molecules should only have basic properties that are handled in atom parsing
    if !properties.is_empty() {
        eprintln!("Warning: Properties found in standard parser - some may be ignored");
    }

    // Add conformer if we have coordinates
    if !coordinates.is_empty() {
        let conformer = Conformer::from_positions(coordinates);
        molecule.conformers.push(conformer);
    }

    molecule
}
