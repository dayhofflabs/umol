//! MOL file parser

use nom::bytes::complete::tag;
use nom::character::complete::multispace0;
use nom::combinator::{map, opt, value};
use nom::sequence::terminated;
use nom::{error, Err, Parser};
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use serde::{Deserialize, Serialize};

use crate::io::ctab::atom::AtomSymbol;
use crate::io::ctab::bond::BondType;
use crate::io::ctab::molecule::{Molecule, MoleculeStandard};
use crate::io::ctab::parser::{ctab_block, ctab_block_standard};

pub mod header;
use header::{header, Header};
use umol::error::DataError;
use umol::Result;

/// Complete MOL file structure (Header + Molecule)
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

/// Parse MOL block (general parser, handles all features)
///
/// This function parses just the CTAB portion (no header), starting from the counts line.
/// Use `parse_mol_file()` for complete MOL files with headers.
#[allow(dead_code)]
fn mol_block<'a>() -> impl Parser<&'a [u8], Output = Molecule, Error = error::Error<&'a [u8]>>
{
    move |input: &'a [u8]| {
        let (remaining, molecule) = ctab_block().parse(input)?;
        let (remaining, _) = footer().parse(remaining)?;
        Ok((remaining, molecule))
    }
}

/// Parse MOL block (standard parser, optimized for performance, standard molecules only)
///
/// This function parses just the CTAB portion (no header), starting from the counts line.
/// Use `parse_mol_file()` for complete MOL files with headers.
fn mol_block_standard<'a>(
) -> impl Parser<&'a [u8], Output = MoleculeStandard, Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let (remaining, molecule) = ctab_block_standard().parse(input)?;
        let (remaining, _) = footer().parse(remaining)?;
        Ok((remaining, molecule))
    }
}

/// Check if molecule contains non-standard/query features
///
/// Returns true if the molecule contains any features that are not supported
/// in the standard MOL format, including:
/// - Query atom symbols (atom lists, R-groups, etc.)
/// - Query bond types (SingleOrDouble, Any, Zero, etc.)
/// - Query-specific properties (ring bond count, topology, etc.)
pub fn has_query_features(molecule: &Molecule) -> bool {
    // Check atoms for non-standard features
    for node_idx in molecule.graph.node_indices() {
        if let Some(atom) = molecule.graph.node_weight(node_idx) {
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
    }

    // Check bonds for non-standard features
    for edge_ref in molecule.graph.edge_references() {
        if let Some(bond) = molecule.graph.edge_weight(edge_ref.id()) {
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

            // Check for non-default bond properties (actual query features)
            if let Some(topology) = &bond.topology {
                if !topology.is_default() {
                    return true;
                }
            }
            if let Some(reacting_center) = &bond.reacting_center {
                if !reacting_center.is_default() {
                    return true;
                }
            }
        }
    }

    // Check for S-groups (non-standard feature)
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
    match mol_file().parse(input) {
        Ok((remaining, mol_file)) => {
            // Check if there's unexpected remaining data
            if !remaining.is_empty() && !remaining.iter().all(|&b| b.is_ascii_whitespace()) {
                return Err(DataError::InvalidMolFormat(format!(
                    "Unexpected data after MOL block: {} bytes remaining",
                    remaining.len()
                ))
                .into());
            }
            Ok(mol_file.molecule)
        }
        Err(Err::Error(e)) | Err(Err::Failure(e)) => {
            Err(DataError::InvalidMolFormat(format!("MOL parsing failed: {:?}", e)).into())
        }
        Err(Err::Incomplete(_)) => {
            Err(DataError::InvalidMolFormat("Incomplete MOL data".to_string()).into())
        }
    }
}

/// Parse MOL bytes into a MoleculeStandard (optimized, standard molecules only)
///
/// This is the high-performance parsing function for standard molecules.
/// It will fail if the MOL file contains query features.
pub fn parse_mol_standard(input: &[u8]) -> Result<MoleculeStandard> {
    // Parse header first to skip it, then use the standard CTAB parser
    match header().parse(input) {
        Ok((remaining, _)) => {
            // Parse the CTAB portion with the standard parser
            match mol_block_standard().parse(remaining) {
                Ok((remaining, molecule)) => {
                    // Check if there's unexpected remaining data
                    if !remaining.is_empty() && !remaining.iter().all(|&b| b.is_ascii_whitespace())
                    {
                        return Err(DataError::InvalidMolFormat(format!(
                            "Unexpected data after MOL block: {} bytes remaining",
                            remaining.len()
                        ))
                        .into());
                    }
                    Ok(molecule)
                }
                Err(Err::Error(e)) | Err(Err::Failure(e)) => Err(DataError::InvalidMolFormat(
                    format!("Standard MOL parsing failed: {:?}", e),
                )
                .into()),
                Err(Err::Incomplete(_)) => {
                    Err(DataError::InvalidMolFormat("Incomplete MOL data".to_string()).into())
                }
            }
        }
        Err(Err::Error(e)) | Err(Err::Failure(e)) => {
            Err(DataError::InvalidMolFormat(format!("Standard MOL parsing failed: {:?}", e)).into())
        }
        Err(Err::Incomplete(_)) => {
            Err(DataError::InvalidMolFormat("Incomplete MOL data".to_string()).into())
        }
    }
}

/// Parse MOL bytes into a MolFile (includes header information)
///
/// Use this when you need access to the MOL file header (name, program info, comment)
/// in addition to the molecular structure.
pub fn parse_mol_file(input: &[u8]) -> Result<MolFile> {
    match mol_file().parse(input) {
        Ok((remaining, mol_file)) => {
            // Check if there's unexpected remaining data
            if !remaining.is_empty() && !remaining.iter().all(|&b| b.is_ascii_whitespace()) {
                return Err(DataError::InvalidMolFormat(format!(
                    "Unexpected data after MOL block: {} bytes remaining",
                    remaining.len()
                ))
                .into());
            }
            Ok(mol_file)
        }
        Err(Err::Error(e)) | Err(Err::Failure(e)) => {
            Err(DataError::InvalidMolFormat(format!("MOL file parsing failed: {:?}", e)).into())
        }
        Err(Err::Incomplete(_)) => {
            Err(DataError::InvalidMolFormat("Incomplete MOL data".to_string()).into())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[test]
    fn test_footer_empty() {
        let result = footer().parse(b"");
        assert!(result.is_ok());
        let (remaining, _) = result.unwrap();
        assert!(remaining.is_empty());
    }

    #[test]
    fn test_footer_whitespace() {
        let result = footer().parse(b"\n\n  \t\n");
        assert!(result.is_ok());
        let (remaining, _) = result.unwrap();
        assert!(remaining.is_empty());
    }

    #[test]
    fn test_footer_sdf_delimiter() {
        let result = footer().parse(b"\n$$$$\n");
        assert!(result.is_ok());
        let (remaining, _) = result.unwrap();
        assert!(remaining.is_empty());
    }

    #[test]
    fn test_footer_sdf_delimiter_with_whitespace() {
        let result = footer().parse(b"\n\n$$$$\n  \n");
        assert!(result.is_ok());
        let (remaining, _) = result.unwrap();
        assert!(remaining.is_empty());
    }

    #[test]
    fn test_parse_mol_file() {
        // Test parsing complete MOL file with header extraction
        let mol_content = b"Ethane\nRDKit\nTest molecule\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n";
        let result = parse_mol_file(mol_content);
        assert!(
            result.is_ok(),
            "Complete MOL file should parse successfully"
        );

        let mol_file = result.unwrap();
        // Check header
        assert_eq!(mol_file.header.name, "Ethane");
        assert_eq!(mol_file.header.program_info, "RDKit");
        assert_eq!(mol_file.header.comment, "Test molecule");

        // Check molecule
        assert_eq!(mol_file.molecule.graph.node_count(), 2);
        assert_eq!(mol_file.molecule.graph.edge_count(), 1);
    }

    #[test]
    fn test_has_query_features() {
        // Test standard molecule
        let standard_mol = b"Methane\nRDKit          3D\nGenerated by RDKit\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0\nM  END\n";
        let molecule = parse_mol(standard_mol).unwrap();
        assert!(!has_query_features(&molecule), "Standard molecule should not have query features");
        
        // Test simple ethane molecule
        let ethane_mol = b"Ethane\nRDKit\nTest molecule\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n";
        let molecule2 = parse_mol(ethane_mol).unwrap();
        assert!(!has_query_features(&molecule2), "Simple ethane should not have query features");
    }

    #[rstest]
    #[case(b"Methane\nRDKit          3D\nGenerated by RDKit\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0\nM  END\n",
       "valid")]
    fn test_parse_mol(#[case] mol_str: &[u8], #[case] desc: &str) {
        let result = parse_mol(mol_str);
        assert!(
            result.is_ok(),
            "{} should have succeeded: {:?}",
            desc,
            result.err()
        );

        let molecule = result.unwrap();
        assert!(
            !has_query_features(&molecule),
            "{} should not have query features",
            desc
        );

        assert_eq!(molecule.atom_count(), 1, "{} should have 1 atom", desc);
        assert_eq!(molecule.bond_count(), 0, "{} should have 0 bonds", desc);
    }

    #[rstest]
    #[case(b"Ethane\nRDKit          3D\nGenerated by RDKit\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n",
      "standard ethane")]
    fn test_parse_mol_standard(#[case] mol_str: &[u8], #[case] desc: &str) {
        let result = parse_mol_standard(mol_str);
        assert!(
            result.is_ok(),
            "{} should have succeeded: {:?}",
            desc,
            result.err()
        );

        let molecule = result.unwrap();
        assert_eq!(molecule.atom_count(), 2, "{} should have 2 atoms", desc);
        assert_eq!(molecule.bond_count(), 1, "{} should have 1 bond", desc);
    }

    #[test]
    fn test_parse_mol_standard_query() {
        // MOL with query atom 'A' (any atom)
        let query_mol = b"Query\nRDKit          3D\nWith query atom\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

        let result = parse_mol_standard(query_mol);
        assert!(
            result.is_err(),
            "Standard parser should fail on query atoms"
        );

        let error = result.unwrap_err();
        let error_string = format!("{}", error);
        assert!(
            error_string.contains("MOL parsing failed")
                || error_string.contains("Standard MOL parsing failed"),
            "Error should mention parsing failure: {}",
            error_string
        );
    }
}
