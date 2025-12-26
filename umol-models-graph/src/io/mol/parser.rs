//! MOL file parser

use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use serde::{Deserialize, Serialize};

use crate::io::ctab::config::{CtabParseFlags, MolIoConfig};
use crate::io::ctab::molecule::{Molecule, MoleculeLike};
use crate::io::ctab::parser::{basic_ctab_block, ctab_block};
use crate::io::ctfile::error::ParseError;
use crate::io::mol::parser::header::Header;
use crate::io::utils::normalize_whitespace;

pub mod header;
use nom::Parser;

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

/// Check if molecule contains advanced features
///
/// Returns true if the molecule contains any features that are not supported
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
                || atomlike.substitution_count.is_some()
                || atomlike.unsaturated.is_some()
                || atomlike.link_atom.is_some()
            {
                return true;
            }
        }
    }

    // Advanced bond features
    for edge_ref in molecule.graph.edge_references() {
        if let Some(bond) = molecule.graph.edge_weight(edge_ref.id()) {
            if bond.bond_type.is_bondlike()
                || bond.topology.is_some_and(|t| !t.is_default())
                || bond.reacting_center.is_some_and(|r| !r.is_default())
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
pub fn parse_mol_moleculelike_str(input: &str) -> std::result::Result<MoleculeLike, ParseError> {
    parse_mol_moleculelike(input.as_bytes())
}

/// Parse a MOL string into a Molecule (optimized, basic molecules only)
///
/// This is the high-performance parsing function for basic molecules.
/// It will fail if the MOL file contains query features.
pub fn parse_mol_str(input: &str) -> std::result::Result<Molecule, ParseError> {
    parse_mol(input.as_bytes())
}

/// Parse MOL bytes into a Molecule
///
/// This is the primary parsing function that handles both basic and query molecules.
/// Stops parsing at M  END and ignores trailing input.
pub fn parse_mol_moleculelike(input: &[u8]) -> std::result::Result<MoleculeLike, ParseError> {
    let config = MolIoConfig::lenient();
    let flags = config.parse_flags;

    let bytes = if flags.contains(CtabParseFlags::UNICODE) {
        normalize_whitespace(input).into_owned()
    } else {
        input.to_vec()
    };

    let (remaining, _header) = header::header()
        .parse(&bytes)
        .map_err(|e| ParseError::from_nom(e, 0, &bytes))?;

    let (_, molecule) = ctab_block(remaining, &flags, 3)?;

    Ok(molecule)
}

/// Parse MOL bytes into a Molecule (optimized, basic molecules only)
///
/// This is the optimized parsing function for basic molecules.
/// It will fail if the MOL file contains query features.
/// Stops parsing at M  END and ignores trailing input.
pub fn parse_mol(input: &[u8]) -> std::result::Result<Molecule, ParseError> {
    let config = MolIoConfig::basic();
    let flags = config.parse_flags;

    let (remaining, _header) = header::header()
        .parse(input)
        .map_err(|e| ParseError::from_nom(e, 0, input))?;

    let (_, molecule) = basic_ctab_block(remaining, &flags, 3)?;

    Ok(molecule)
}

#[cfg(test)]
mod tests;
