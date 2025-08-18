//! MOL file parser

use bstr::{join, ByteSlice};
use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::multispace0;
use nom::combinator::{all_consuming, complete, map, opt, value};
use nom::sequence::terminated;
use nom::{error, Parser};
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use serde::{Deserialize, Serialize};

use crate::io::ctab::atom::{Atom, AtomStandard};
use crate::io::ctab::bond::{Bond, BondStandard};
use crate::io::ctab::molecule::{Molecule, MoleculeStandard};
use crate::io::ctab::parser::{
    atom_alias_entry, atom_input, atom_input_standard, bond_input, bond_input_standard,
    counts_input, ctab_block, ctab_block_standard, legacy_atom_list_input, property_input,
    property_input_standard, Counts, PropertyEntries,
};
use crate::io::mol::parser::header::{header_input, header, Header};

pub mod header;
use umol::error::DataError;
use umol::Result;

/// Complete MOL file structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MolFile {
    pub header: Header,
    pub molecule: Molecule,
}

impl MolFile {
    /// Create a new MOL file
    pub fn new(header: Header, molecule: Molecule) -> Self {
        Self { header, molecule }
    }
}

/// Complete MOL file structure for standard molecules
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MolFileStandard {
    pub header: Header,
    pub molecule: MoleculeStandard,
}

impl MolFileStandard {
    /// Create a new MOL file
    pub fn new(header: Header, molecule: MoleculeStandard) -> Self {
        Self { header, molecule }
    }
}

/// Parse optional footer content (whitespace and SDF record delimiter)
fn footer<'a>() -> impl Parser<&'a [u8], Output = (), Error = error::Error<&'a [u8]>> {
    value((), (multispace0, opt((tag("$$$$"), multispace0))))
}

/// Parse complete MOL file (header + CTAB block)
fn mol_file<'a>() -> impl Parser<&'a [u8], Output = MolFile, Error = error::Error<&'a [u8]>> {
    map(
        terminated((header::header(), ctab_block()), footer()),
        |(header, molecule)| MolFile::new(header, molecule),
    )
}

/// Parse complete MOL file (header + CTAB block) for standard molecules
fn mol_file_standard<'a>(
) -> impl Parser<&'a [u8], Output = MolFileStandard, Error = error::Error<&'a [u8]>> {
    map(
        terminated((header::header(), ctab_block_standard()), footer()),
        |(header, molecule)| MolFileStandard::new(header, molecule),
    )
}

/// Check if molecule contains non-standardquery features
///
/// Returns true if the molecule contains any features that are not supported
/// in the standard MOL format, including:
/// - Query atom symbols (atom lists, R-groups, etc.)
/// - Query bond types (SingleOrDouble, Any, Zero, etc.)
/// - Query-specific properties (ring bond count, topology, etc.)
pub fn has_nonstandard_features(molecule: &Molecule) -> bool {
    // Check atoms for non-standard features
    for node_idx in molecule.graph.node_indices() {
        if let Some(atom) = molecule.graph.node_weight(node_idx) {
            // Non-standard atom features
            if !atom.symbol.is_standard()
                || atom.attachment_point.is_some()
                || atom.attachment_order.is_some()
                || atom.ring_bond_count.is_some()
                || atom.substitution_count.is_some()
                || atom.unsaturated.is_some()
                || atom.link_atom.is_some()
            {
                return true;
            }
        }
    }

    // Non-standard bond features
    for edge_ref in molecule.graph.edge_references() {
        if let Some(bond) = molecule.graph.edge_weight(edge_ref.id()) {
            if !bond.bond_type.is_standard()
                || bond.topology.map_or(false, |t| !t.is_default())
                || bond.reacting_center.map_or(false, |r| !r.is_default())
            {
                return true;
            }
        }
    }

    // S-groups (non-standard feature)
    if !molecule.sgroups.is_empty() {
        return true;
    }

    false
}

/// Parse a MOL string into a Molecule
///
/// This is the main parsing function that handles both standard and query molecules.
pub fn parse_mol_str(input: &str) -> Result<Molecule> {
    parse_mol(input.as_bytes())
}

/// Parse a MOL string into a MoleculeStandard (optimized, standard molecules only)
///
/// This is the high-performance parsing function for standard molecules.
/// It will fail if the MOL file contains query features.
pub fn parse_mol_standard_str(input: &str) -> Result<MoleculeStandard> {
    parse_mol_standard(input.as_bytes())
}

/// Parse MOL bytes into a Molecule
///
/// This is the primary parsing function that handles both standard and query molecules.
pub fn parse_mol(input: &[u8]) -> Result<Molecule> {
    parse_mol_file(input).map(|mol_file| mol_file.molecule)
}

/// Parse MOL bytes into a MoleculeStandard (optimized, standard molecules only)
///
/// This is the high-performance parsing function for standard molecules.
/// It will fail if the MOL file contains query features.
pub fn parse_mol_standard(input: &[u8]) -> Result<MoleculeStandard> {
    all_consuming(terminated(
        map(
            complete(terminated((header(), ctab_block_standard()), footer())),
            |(_, molecule)| molecule,
        ),
        multispace0,
    ))
    .parse(input)
    .map(|(_, molecule)| molecule)
    .map_err(|e| {
        DataError::InvalidMolFormat(format!("Standard MOL parsing failed: {:?}", e)).into()
    })
}

/// Parse MOL bytes into a MolFile (includes header information)
///
/// Use this when you need access to the MOL file header (name, program info, comment)
/// in addition to the molecular structure.
pub fn parse_mol_file(input: &[u8]) -> Result<MolFile> {
    all_consuming(complete(terminated(mol_file(), multispace0)))
        .parse(input)
        .map(|(_, mol_file)| mol_file)
        .map_err(|e| {
            DataError::InvalidMolFormat(format!("MOL file parsing failed: {:?}", e)).into()
        })
}

/// Parse MOL bytes into a MolFileStandard (includes header information)
///
/// This is the high-performance parsing function for standard molecules.
/// It will fail if the MOL file contains query features.
pub fn parse_mol_file_standard(input: &[u8]) -> Result<MolFileStandard> {
    all_consuming(complete(terminated(mol_file_standard(), multispace0)))
        .parse(input)
        .map(|(_, mol_file)| mol_file)
        .map_err(|e| {
            DataError::InvalidMolFormat(format!("MOL file parsing failed: {:?}", e)).into()
        })
}

/// Represents the parsed content of a single logical line or multi-line block from a MOL file.
#[derive(Debug, Clone, PartialEq)]
pub enum ParsedMolLine {
    Header(String),
    Counts(Counts),
    Atom(Atom),
    Bond(Bond),
    Property(PropertyEntries),
    LegacyAtomList,
    End,
    Empty,
    Unknown(String),
}

/// Represents the parsed content of a single logical line or multi-line block from a MOL file.
#[derive(Debug, Clone, PartialEq)]
pub enum ParsedMolStandardLine {
    Header(String),
    Counts(Counts),
    Atom(AtomStandard),
    Bond(BondStandard),
    Property(PropertyEntries),
    End,
    Empty,
    Unknown(String),
}

/// Progressively parse a MOL file, line by line, for diagnostic purposes.
///
/// This function processes one logical line at a time, handling multi-line blocks
/// like atom aliases, and returns a vector of `ParsedMolLine` enums representing
/// the content of each line. It uses the general parser that supports query features.
pub fn parse_mol_progressive(input: &[u8]) -> Vec<ParsedMolLine> {
    let mut results = Vec::new();
    let mut lines = input.lines();

    for _ in 0..3 {
        if let Some(line) = lines.next() {
            results.push(
                header_input()
                    .parse(line)
                    .map(|(_, h)| ParsedMolLine::Header(h))
                    .unwrap_or(ParsedMolLine::Unknown(line.to_str_lossy().to_string())),
            );
        }
    }

    let (atom_count, bond_count) = if let Some(line) = lines.next() {
        match counts_input().parse(line) {
            Ok((_, counts)) => {
                let (ac, bc) = (counts.atoms() as usize, counts.bonds() as usize);
                results.push(ParsedMolLine::Counts(counts));
                (ac, bc)
            }
            Err(_) => {
                results.push(ParsedMolLine::Unknown(line.to_str_lossy().to_string()));
                (0, 0)
            }
        }
    } else {
        (0, 0)
    };

    for _ in 0..atom_count {
        if let Some(line) = lines.next() {
            results.push(
                atom_input()
                    .parse(line)
                    .map(|(_, atom)| ParsedMolLine::Atom(atom))
                    .unwrap_or(ParsedMolLine::Unknown(line.to_str_lossy().to_string())),
            );
        }
    }

    for _ in 0..bond_count {
        if let Some(line) = lines.next() {
            results.push(
                bond_input()
                    .parse(line)
                    .map(|(_, (_, _, bond))| ParsedMolLine::Bond(bond))
                    .unwrap_or(ParsedMolLine::Unknown(line.to_str_lossy().to_string())),
            );
        }
    }

    let mut lines = lines.peekable();
    while let Some(line) = lines.next() {
        if line.trim().is_empty() {
            results.push(ParsedMolLine::Empty);
            continue;
        }

        if line.starts_with(b"A  ") {
            if let Some(next_line) = lines.next() {
                let combined = join(b"\n", &[&line[3..], next_line]);
                results.push(
                    atom_alias_entry()
                        .parse(&combined)
                        .map(|(_, entry)| {
                            ParsedMolLine::Property(PropertyEntries::AtomAliasEntry(entry))
                        })
                        .unwrap_or(ParsedMolLine::Unknown(combined.to_str_lossy().to_string())),
                );
            } else {
                results.push(ParsedMolLine::Unknown(line.to_str_lossy().to_string()));
            }
        } else if line.starts_with(b"M  END") {
            results.push(ParsedMolLine::End);
        } else {
            // Here we must try legacy list first as it has no prefix
            let mut property_parser = alt((
                map(legacy_atom_list_input(), |p| p),
                map(property_input(), |p| p),
            ));
            results.push(
                property_parser
                    .parse(line)
                    .map(|(_, prop)| ParsedMolLine::Property(prop))
                    .unwrap_or(ParsedMolLine::Unknown(line.to_str_lossy().to_string())),
            );
        }
    }

    results
}

/// Progressively parse a standard MOL file, line by line, for diagnostic purposes.
///
/// This function processes one logical line at a time and returns a vector of
/// `ParsedMolStandardLine` enums. It uses the optimized standard parser and will not
/// correctly handle query features.
pub fn parse_mol_progressive_standard(input: &[u8]) -> Vec<ParsedMolStandardLine> {
    let mut results = Vec::new();
    let mut lines = input.lines();

    for _ in 0..3 {
        if let Some(line) = lines.next() {
            results.push(
                header_input()
                    .parse(line)
                    .map(|(_, h)| ParsedMolStandardLine::Header(h))
                    .unwrap_or(ParsedMolStandardLine::Unknown(
                        line.to_str_lossy().to_string(),
                    )),
            );
        }
    }

    let (atom_count, bond_count) = if let Some(line) = lines.next() {
        match counts_input().parse(line) {
            Ok((_, counts)) => {
                let (ac, bc) = (counts.atoms() as usize, counts.bonds() as usize);
                results.push(ParsedMolStandardLine::Counts(counts));
                (ac, bc)
            }
            Err(_) => {
                results.push(ParsedMolStandardLine::Unknown(
                    line.to_str_lossy().to_string(),
                ));
                (0, 0)
            }
        }
    } else {
        (0, 0)
    };

    for _ in 0..atom_count {
        if let Some(line) = lines.next() {
            results.push(
                atom_input_standard()
                    .parse(line)
                    .map(|(_, atom)| ParsedMolStandardLine::Atom(atom))
                    .unwrap_or(ParsedMolStandardLine::Unknown(
                        line.to_str_lossy().to_string(),
                    )),
            );
        }
    }

    for _ in 0..bond_count {
        if let Some(line) = lines.next() {
            results.push(
                bond_input_standard()
                    .parse(line)
                    .map(|(_, (_, _, bond))| ParsedMolStandardLine::Bond(bond))
                    .unwrap_or(ParsedMolStandardLine::Unknown(
                        line.to_str_lossy().to_string(),
                    )),
            );
        }
    }

    let mut lines = lines.peekable();
    while let Some(line) = lines.next() {
        if line.trim().is_empty() {
            results.push(ParsedMolStandardLine::Empty);
            continue;
        }

        if line.starts_with(b"A  ") {
            if let Some(next_line) = lines.next() {
                let combined = join(b"\n", &[&line[3..], next_line]);
                results.push(
                    atom_alias_entry()
                        .parse(&combined)
                        .map(|(_, entry)| {
                            ParsedMolStandardLine::Property(PropertyEntries::AtomAliasEntry(entry))
                        })
                        .unwrap_or(ParsedMolStandardLine::Unknown(
                            combined.to_str_lossy().to_string(),
                        )),
                );
            } else {
                results.push(ParsedMolStandardLine::Unknown(
                    line.to_str_lossy().to_string(),
                ));
            }
        } else if line.starts_with(b"M  END") {
            results.push(ParsedMolStandardLine::End);
        } else {
            results.push(
                property_input_standard()
                    .parse(line)
                    .map(|(_, prop)| ParsedMolStandardLine::Property(prop))
                    .unwrap_or(ParsedMolStandardLine::Unknown(
                        line.to_str_lossy().to_string(),
                    )),
            );
        }
    }

    results
}

#[cfg(test)]
mod tests;
