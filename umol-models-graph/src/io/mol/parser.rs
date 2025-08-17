//! MOL file parser

use nom::bytes::complete::tag;
use nom::character::complete::multispace0;
use nom::combinator::{all_consuming, complete, map, opt, value};
use nom::sequence::terminated;
use nom::{error, Parser};
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use serde::{Deserialize, Serialize};

use crate::io::ctab::molecule::{Molecule, MoleculeStandard};
use crate::io::ctab::parser::{ctab_block, ctab_block_standard};

pub mod header;
use header::{header, Header};
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
fn mol_file_standard<'a>() -> impl Parser<&'a [u8], Output = MolFileStandard, Error = error::Error<&'a [u8]>> {
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
        assert!(
            !has_nonstandard_features(&molecule),
            "Standard molecule should not have query features"
        );

        // Test simple ethane molecule
        let ethane_mol = b"Ethane\nRDKit\nTest molecule\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n";
        let molecule2 = parse_mol(ethane_mol).unwrap();
        assert!(
            !has_nonstandard_features(&molecule2),
            "Simple ethane should not have query features"
        );
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
            !has_nonstandard_features(&molecule),
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

    #[rstest]
    #[case(b"Ethane\nRDKit\nTest molecule\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n",
      "standard ethane")]
    fn test_parse_mol_file_standard(#[case] mol_str: &[u8], #[case] desc: &str) {
        let result = parse_mol_file_standard(mol_str);
        assert!(result.is_ok(), "{} should have succeeded: {:?}", desc, result.err());
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
