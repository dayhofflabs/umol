//! MOL file parser

use nom::combinator::{complete, map};
use nom::{error, Parser};
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use serde::{Deserialize, Serialize};

use crate::io::config::ParseFlags;
use crate::io::ctab::molecule::{Molecule, MoleculeLike};
use crate::io::ctab::parser::{basic_ctab_block, ctab_block};
use crate::io::mol::parser::header::Header;

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

/// Complete MOL file structure (includes header information)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MolFileLike {
    pub header: Header,
    pub molecule: MoleculeLike,
}

impl MolFileLike {
    /// Create a new MOL file
    pub fn new(header: Header, molecule: MoleculeLike) -> Self {
        Self { header, molecule }
    }
}



/// Parse complete MOL file (header + CTAB block)
pub(crate) fn mol_file<'a>(
) -> impl Parser<&'a [u8], Output = MolFile, Error = error::Error<&'a [u8]>> {
    map(
        (
            header::header(),
            basic_ctab_block(ParseFlags::BASIC | ParseFlags::DEBUG),
        ),
        |(header, molecule)| MolFile::new(header, molecule),
    )
}

/// Parse complete MOL file (header + CTAB block)
pub(crate) fn mol_file_moleculelike<'a>(
) -> impl Parser<&'a [u8], Output = MolFileLike, Error = error::Error<&'a [u8]>> {
    map(
        (header::header(), ctab_block(ParseFlags::LENIENT | ParseFlags::DEBUG)),
        |(header, molecule)| MolFileLike::new(header, molecule),
    )
}

/// Check if molecule contains advanced features
///
/// Returns -ue if the molecule contains any features that are not supported
/// in the basic MOL format, including:
/// - Query atom symbols (atom lists, R-groups, etc.)
/// - Query bond types (SingleOrDouble, Any, Zero, etc.)
/// - Query-specific properties (ring bond count, topology, etc.)
pub fn has_advanced_features(molecule: &MoleculeLike) -> bool {
    // Check atoms for advanced features
    for node_idx in molecule.graph.node_indices() {
        if let Some(atomlike) = molecule.graph.node_weight(node_idx) {
            // Advanced atom features
            if atomlike.symbol.is_atomlike()
                || atomlike.attachment_point.is_some()
                || atomlike.attachment_order.is_some()
                || atomlike.ring_bond_count.is_some()
                || atomlike.substitution_count.
                is_some()
                || atomlike.unsaturated.is_some()
                || atomlike.link_atom.is_some()
            {
                return true
            }
        }
    }

    // Advanced bond features
    for edge_ref in molecule.graph.edge_references() {
        if let Some(bond) = molecule.graph.edge_weight(edge_ref.id()) {
            if bond.bond_type.is_bondlike()
                || bond.topology.map_or(false, |t| !t.is_default())
                || bond.reacting_center.map_or(false, |r| !r.is_default())
            {
                return true;
            }
        }
    }

    // S-groups (advanced features)
    if !molecule.sgroups.is_empty() {
        return true;
    }

    false
}

/// Parse a MOL string into a Molecule
///
/// This is the main parsing function that handles both basic and query molecules.
pub fn parse_mol_moleculelike_str(input: &str) -> Result<MoleculeLike> {
    parse_mol_moleculelike(input.as_bytes())
}

/// Parse a MOL string into a Molecule (optimized, basic molecules only)
///
/// This is the high-performance parsing function for basic molecules.
/// It will fail if the MOL file contains query features.
pub fn parse_mol_str(input: &str) -> Result<Molecule> {
    parse_mol(input.as_bytes())
}

/// Parse MOL bytes into a Molecule
///
/// This is the primary parsing function that handles both basic and query molecules.
/// Stops parsing at M  END and ignores any remaining input (e.g., SDF properties).
pub fn parse_mol_moleculelike(input: &[u8]) -> Result<MoleculeLike> {
    complete(mol_file_moleculelike())
        .parse(input)
        .map(|(_, mol_file)| mol_file.molecule)
        .map_err(|e| {
            DataError::InvalidMolFormat(format!("MOL file parsing failed: {:?}", e)).into()
        })
}

/// Parse MOL bytes into a Molecule (optimized, basic molecules only)
///
/// This is the high-performance parsing function for basic molecules.
/// It will fail if the MOL file contains query features.
/// Stops parsing at M  END and ignores any remaining input (e.g., SDF properties).
pub fn parse_mol(input: &[u8]) -> Result<Molecule> {
    complete(mol_file())
        .parse(input)
        .map(|(_, mol_file)| mol_file.molecule)
        .map_err(|e| DataError::InvalidMolFormat(format!("MOL parsing failed: {:?}", e)).into())
}

#[cfg(test)]
mod tests;
