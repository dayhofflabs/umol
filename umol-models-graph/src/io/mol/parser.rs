//! MOL file parser

use nom::combinator::map;
use nom::error::Error as NomError;

use crate::io::ctab::config::{CtabParseFlags, MolIoConfig};
use crate::io::ctab::parser::{ctab_block, extended_ctab_block};
use crate::io::ctfile::error::ParseError;
use crate::io::mol::parser::header::Header;
use crate::io::utils::normalize_whitespace;
use crate::table_ir::{ExtendedMolecule, Molecule};

pub mod header;
use nom::Parser;

/// Complete MOL file structure (basic molecules)
#[derive(Debug, Clone)]
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

/// Parse complete MOL file (header + CTAB block)
pub(crate) fn mol_file<'inp>(
    flags: CtabParseFlags,
    line_num: u32,
) -> impl Parser<&'inp [u8], Output = MolFile, Error = ParseError> + use<'inp> {
    map(
        (header::header(), ctab_block(flags, line_num)),
        |(header, molecule)| MolFile::new(header, molecule),
    )
}

/// Extended MOL file structure (includes query features, S-groups, R-groups)
#[derive(Debug, Clone)]
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

/// Check if extended molecule actually contains extended features
///
/// Return true if the extended molecule contains any features that
/// are not supported in the basic MOL format, including:
/// - Extended atom symbols (atom lists, R-groups, etc.)
/// - Extended bond types (SingleOrDouble, Any, Zero, etc.)
/// - S-groups, R-groups
pub fn has_extended_features(molecule: &ExtendedMolecule) -> bool {
    use crate::table_ir::AtomSymbol;

    // Check atoms for query features
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

    // Check bonds for query features
    for bond in &molecule.bonds {
        if bond.order.is_query() || bond.order.is_extended() {
            return true;
        }
    }

    // S-groups (advanced features)
    if !molecule.sgroups().is_empty() {
        return true;
    }

    // R-groups
    if !molecule.rgroups().is_empty() {
        return true;
    }

    false
}

/// Parse MOL bytes into a Molecule with options (optimized, basic molecules only)
pub fn parse_mol_bytes_with(input: &[u8], config: &MolIoConfig) -> Result<Molecule, ParseError> {
    let flags = config.parse_flags;

    // let (remaining, _header) = header::header()
    //     .parse(input)
    //     .map_err(|e| ParseError::header_from_nom(e, 0))?;

    // let (_, molecule) = ctab_block(remaining, &flags, 3)?;

    Ok(molecule)
}

/// Parse MOL bytes into a Molecule (optimized, basic molecules only)
///
/// Optimized parsing function for basic molecules.
/// Fails if the MOL file contains extended features.
/// Stops parsing at M  END and ignores trailing input.
pub fn parse_mol_bytes(input: &[u8]) -> Result<Molecule, ParseError> {
    let config = MolIoConfig::basic();
    parse_mol_bytes_with(input, &config)
}

/// Parse MOL string into a Molecule with options (optimized, basic molecules only)
pub fn parse_mol_with(input: &str, config: &MolIoConfig) -> Result<Molecule, ParseError> {
    parse_mol_bytes_with(input.as_bytes(), config)
}

/// Parse MOL string into a Molecule (optimized, basic molecules only)
///
/// Optimized parsing function for basic molecules.
/// Fails if the MOL file contains extended features.
pub fn parse_mol(input: &str) -> Result<Molecule, ParseError> {
    parse_mol_bytes(input.as_bytes())
}

/// Parse MOL bytes into an ExtendedMolecule
///
/// Generic parsing function that handles both basic and extended molecules.
/// Stops parsing at M  END and ignores trailing input.
pub fn parse_extended_mol_bytes_with(
    input: &[u8],
    config: &MolIoConfig,
) -> Result<ExtendedMolecule, ParseError> {
    let config = MolIoConfig::extended();
    let flags = config.parse_flags;

    Ok(molecule)
}

/// Parse MOL bytes into an ExtendedMolecule
///
/// Generic parsing function that handles both basic and extended molecules.
/// Stops parsing at M  END and ignores trailing input.
pub fn parse_extended_mol_bytes(input: &[u8]) -> Result<ExtendedMolecule, ParseError> {
    let config = MolIoConfig::extended();
    parse_extended_mol_bytes_with(input, &config)
}

/// Parse MOL string into an ExtendedMolecule with options
///
/// Generic parsing function that handles both basic and extended molecules.
pub fn parse_extended_mol_with(
    input: &str,
    config: &MolIoConfig,
) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_mol_bytes_with(input.as_bytes(), config)
}

/// Parse MOL string into an ExtendedMolecule
///
/// Generic parsing function that handles both basic and extended molecules.
pub fn parse_extended_mol(input: &str) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_mol_bytes(input.as_bytes())
}

#[cfg(test)]
mod tests;
