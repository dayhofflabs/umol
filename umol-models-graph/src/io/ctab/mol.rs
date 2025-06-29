//! MOL file parser

use nom::bytes::complete::tag;
use nom::character::complete::line_ending;
use nom::combinator::opt;
use nom::multi::{count, many0};
use nom::sequence::terminated;
use nom::{error, IResult, Parser};

use crate::atom::{Atom, AtomSymbol};
use crate::bond::Bond;
use crate::conformer::{Conformer, Point3D};
use crate::molecule::{Header, Molecule, ParsedMol};

use super::apply::Apply;
use super::atom::atom_input;
use super::bond::bond_input;
use super::counts::counts_input;
use super::header::header;
use super::properties::{property_input, PropertyEntries};

/// Parse MOL block
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

/// Parse properties block (until M  END or end of input)
fn properties_block<'a>(
    input: &'a [u8],
) -> IResult<&'a [u8], Vec<PropertyEntries>, error::Error<&'a [u8]>> {
    let (input, properties) = many0(terminated(property_input, line_ending)).parse(input)?;
    let (input, _) = opt(terminated(tag("M  END"), opt(line_ending))).parse(input)?;
    Ok((input, properties))
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
