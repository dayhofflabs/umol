//! MOL file parser

use bstr::ByteSlice;
use nom::bytes::complete::tag;
use nom::character::complete::{line_ending, multispace0};
use nom::combinator::opt;
use nom::multi::many0;
use nom::sequence::{preceded, terminated};
use nom::{error, Err, IResult, Parser};

use super::super::ctab::atom::{Atom, AtomStandard, AtomSymbol};
use super::super::ctab::bond::BondType;
use super::super::ctab::bond::{Bond, BondStandard};
use super::super::ctab::molecule::{Molecule, MoleculeStandard, ParsedMol};

use super::super::ctab::parser::accumulator::MoleculeProperties;
use super::super::ctab::parser::atom::{atom_input, atom_input_standard};
use super::super::ctab::parser::bond::{bond_input, bond_input_standard};
use super::super::ctab::parser::counts::counts_input;
use super::super::ctab::parser::properties::{
    legacy_atom_list_input, property_input, property_input_standard, PropertyEntries,
};

/// Parse MOL block (general parser, handles all features)
pub fn mol_block(input: &[u8]) -> IResult<&[u8], ParsedMol, error::Error<&[u8]>> {
    let (remaining, counts) = terminated(counts_input(), line_ending).parse(input)?;
    let atom_count = counts.atoms() as usize;
    let bond_count = counts.bonds() as usize;
    let (remaining, atoms) = atom_block(atom_count).parse(remaining)?;
    let (remaining, bonds) = bond_block(bond_count).parse(remaining)?;
    let (remaining, legacy_properties) = legacy_atom_list_block(remaining)?;
    let (remaining, properties) = properties_block(remaining)?;
    let (remaining, _) = footer(remaining)?;
    let properties = properties.into_iter().chain(legacy_properties).collect();
    let (molecule, is_query) = build_molecule(atoms, bonds, properties);
    Ok((remaining, ParsedMol::new(molecule, is_query)))
}

/// Parse MOL block (standard parser, optimized for performance, standard molecules only)
pub fn mol_block_standard(input: &[u8]) -> IResult<&[u8], MoleculeStandard, error::Error<&[u8]>> {
    let (remaining, counts) = terminated(counts_input(), line_ending).parse(input)?;
    let atom_count = counts.atoms() as usize;
    let bond_count = counts.bonds() as usize;
    let (remaining, atoms) = atom_block_standard(atom_count).parse(remaining)?;
    let (remaining, bonds) = bond_block_standard(bond_count).parse(remaining)?;
    let (remaining, properties) = properties_block_standard(remaining)?;
    let (remaining, _) = footer(remaining)?;
    let molecule = build_molecule_standard(atoms, bonds, properties);
    Ok((remaining, molecule))
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

        let len = lines_iter.as_bytes().as_ptr() as usize - input.as_ptr() as usize;
        let remaining = &input[len..];
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

            // Parse the atom from the line
            let (_, atom) = atom_input_standard().parse(line)?;
            atoms.push(atom);
        }

        let len = lines_iter.as_bytes().as_ptr() as usize - input.as_ptr() as usize;
        let remaining = &input[len..];
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

            let (_, (atom1, atom2, bond)) = bond_input().parse(line)?;
            bonds.push((atom1, atom2, bond));
        }

        let len = lines_iter.as_bytes().as_ptr() as usize - input.as_ptr() as usize;
        let remaining = &input[len..];
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

        let len = lines_iter.as_bytes().as_ptr() as usize - input.as_ptr() as usize;
        let remaining = &input[len..];
        Ok((remaining, bonds))
    }
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

/// Parse optional footer content (whitespace and SDF record delimiter)
fn footer(input: &[u8]) -> IResult<&[u8], (), error::Error<&[u8]>> {
    let (input, _) = multispace0.parse(input)?;
    let (input, _) = opt(preceded(tag("$$$$"), multispace0)).parse(input)?;
    Ok((input, ()))
}

/// Detect if molecule contains non-standard/query features
///
/// Returns true if the molecule contains any features that are not supported
/// in the standard MOL format, including:
/// - Query atom symbols (atom lists, R-groups, etc.)
/// - Query bond types (SingleOrDouble, Any, Zero, etc.)
/// - Query-specific properties (ring bond count, topology, etc.)
/// - Extended properties not in property_input_standard
fn detect_nonstandard_features(
    atoms: &[Atom],
    bonds: &[(usize, usize, Bond)],
    properties: &[PropertyEntries],
) -> bool {
    // Check atoms for non-standard features
    for atom in atoms {
        match &atom.symbol {
            AtomSymbol::AtomList(_)
            | AtomSymbol::Query(_)
            | AtomSymbol::LonePair
            | AtomSymbol::RGroup(_) => {
                return true;
            }
            AtomSymbol::Element(_) | AtomSymbol::NamedIsotope(_) => {}
        }

        // Check for non-standard atom properties
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

    // Check bonds for non-standard features
    for (_from, _to, bond) in bonds {
        match bond.bond_type {
            BondType::SingleOrDouble
            | BondType::SingleOrAromatic
            | BondType::DoubleOrAromatic
            | BondType::Any
            | BondType::Zero => {
                return true;
            }
            BondType::Single | BondType::Double | BondType::Triple | BondType::Aromatic => {}
        }

        // Check for non-standard bond properties
        if bond.topology.is_some() || bond.reacting_center.is_some() {
            return true;
        }
    }
    // Check properties for non-standard entries
    for property in properties {
        match property {
            PropertyEntries::AtomListEntry(_)
            | PropertyEntries::AttachmentPointEntries(_)
            | PropertyEntries::AtomAttachmentOrderEntry(_)
            | PropertyEntries::RingBondCountEntries(_)
            | PropertyEntries::SubstitutionCountEntries(_)
            | PropertyEntries::UnsaturatedAtomEntries(_)
            | PropertyEntries::LinkAtomEntries(_)
            | PropertyEntries::RGroupLabelEntries(_)
            | PropertyEntries::RGroupLogicEntry(_)
            | PropertyEntries::SGroupConnectivityEntries(_)
            | PropertyEntries::SGroupExpansionEntries(_)
            | PropertyEntries::SGroupParentAtomEntry(_)
            | PropertyEntries::SGroupCorrespondenceEntry(_)
            | PropertyEntries::SGroupDisplayInfoEntry(_)
            | PropertyEntries::SGroupConnectingBondEntry(_)
            | PropertyEntries::SGroupDataDescriptionEntry(_)
            | PropertyEntries::SGroupDataDisplayEntry(_)
            | PropertyEntries::SGroupHierarchyEntries(_)
            | PropertyEntries::SGroupComponentEntries(_)
            | PropertyEntries::SGroupDataEntry(_) => return true,
            _ => {}
        }
    }
    false
}

/// Build molecule from parsed components
fn build_molecule(
    atoms: Vec<Atom>,
    bonds: Vec<(usize, usize, Bond)>,
    properties: Vec<PropertyEntries>,
) -> (Molecule, bool) {
    let is_query = detect_nonstandard_features(&atoms, &bonds, &properties);
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

    (molecule, is_query)
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

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_footer_empty() {
//         let result = footer(b"");
//         assert!(result.is_ok());
//         let (remaining, _) = result.unwrap();
//         assert!(remaining.is_empty());
//     }

//     #[test]
//     fn test_footer_whitespace() {
//         let result = footer(b"\n\n  \t\n");
//         assert!(result.is_ok());
//         let (remaining, _) = result.unwrap();
//         assert!(remaining.is_empty());
//     }

//     #[test]
//     fn test_footer_sdf_delimiter() {
//         let result = footer(b"\n$$$$\n");
//         assert!(result.is_ok());
//         let (remaining, _) = result.unwrap();
//         assert!(remaining.is_empty());
//     }

//     #[test]
//     fn test_footer_sdf_delimiter_with_whitespace() {
//         let result = footer(b"\n\n$$$$\n  \n");
//         assert!(result.is_ok());
//         let (remaining, _) = result.unwrap();
//         assert!(remaining.is_empty());
//     }

//     #[test]
//     fn test_mol_with_sdf_footer() {
//         let mol_content = b"test\nRDKit\ncomment\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n\n$$$$\n";
//         let result = mol_block(mol_content);
//         assert!(
//             result.is_ok(),
//             "MOL with SDF footer should parse successfully"
//         );
//         let (remaining, _) = result.unwrap();
//         assert!(
//             remaining.is_empty(),
//             "All content should be consumed including SDF footer"
//         );
//     }
// }
