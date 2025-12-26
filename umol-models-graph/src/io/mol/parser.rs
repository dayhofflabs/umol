//! MOL file parser

use serde::{Deserialize, Serialize};

use crate::io::ctab::config::{CtabParseFlags, MolIoConfig};
use crate::io::ctab::parser::{ctab_block, extended_ctab_block};
use crate::io::ctfile::error::ParseError;
use crate::io::mol::parser::header::Header;
use crate::io::utils::normalize_whitespace;
use crate::simple_ir::{ExtendedMolecule, Molecule};

pub mod header;
use nom::Parser;

/// Complete MOL file structure (basic molecules)
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

/// Extended MOL file structure (includes query features, S-groups, R-groups)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtendedMolFile {
    pub header: Header,
    pub molecule: ExtendedMolecule,
}

impl ExtendedMolFile {
    /// Create a new extended MOL file
    pub fn new(header: Header, molecule: ExtendedMolecule) -> Self {
        Self { header, molecule }
    }
}

/// Check if extended molecule contains advanced features
///
/// Returns true if the molecule contains any features that are not supported
/// in the basic MOL format, including:
/// - Query atom symbols (atom lists, R-groups, etc.)
/// - Query bond types (SingleOrDouble, Any, Zero, etc.)
/// - S-groups, R-groups
pub fn has_advanced_features(molecule: &ExtendedMolecule) -> bool {
    use crate::simple_ir::AtomSymbol;

    // Check atoms for query features
    for atom in &molecule.atoms {
        match &atom.symbol {
            AtomSymbol::Query(_)
            | AtomSymbol::AtomList(_)
            | AtomSymbol::RGroup(_)
            | AtomSymbol::LonePair
            | AtomSymbol::Variable(_)
            | AtomSymbol::Pseudoatom(_)
            | AtomSymbol::Unknown => return true,
            AtomSymbol::Element(_) | AtomSymbol::NamedIsotope(_) => {}
        }
    }

    // Check bonds for query features
    for bond in &molecule.bonds {
        if bond.order.is_query() || bond.order.is_extended() {
            return true;
        }
    }

    // S-groups (advanced features)
    if !molecule.sgroups.is_empty() {
        return true;
    }

    // R-groups
    if !molecule.rgroups.is_empty() {
        return true;
    }

    false
}

/// Parse a MOL string into an ExtendedMolecule
///
/// This is the main parsing function that handles both basic and query molecules.
pub fn parse_extended_mol_str(input: &str) -> std::result::Result<ExtendedMolecule, ParseError> {
    parse_extended_mol(input.as_bytes())
}

/// Parse a MOL string into a Molecule (optimized, basic molecules only)
///
/// This is the high-performance parsing function for basic molecules.
/// It will fail if the MOL file contains query features.
pub fn parse_mol_str(input: &str) -> std::result::Result<Molecule, ParseError> {
    parse_mol(input.as_bytes())
}

/// Parse MOL bytes into an ExtendedMolecule
///
/// This is the primary parsing function that handles both basic and query molecules.
/// Stops parsing at M  END and ignores trailing input.
pub fn parse_extended_mol(input: &[u8]) -> std::result::Result<ExtendedMolecule, ParseError> {
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

    let (_, molecule) = extended_ctab_block(remaining, &flags, 3)?;

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

    let (_, molecule) = ctab_block(remaining, &flags, 3)?;

    Ok(molecule)
}

#[cfg(test)]
mod tests;
