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
use crate::io::mol::parser::header::header_input;

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

fn mol_input<'a>() -> impl Parser<&'a [u8], Output = ParsedMolLine, Error = error::Error<&'a [u8]>>
{
    alt((
        map(counts_input(), ParsedMolLine::Counts),
        map(atom_input(), ParsedMolLine::Atom),
        map(bond_input(), |(_, _, bond)| ParsedMolLine::Bond(bond)),
        map(legacy_atom_list_input(), ParsedMolLine::Property),
        map(tag("M  END"), |_| ParsedMolLine::End),
        map(property_input(), ParsedMolLine::Property),
    ))
}

fn mol_input_standard<'a>(
) -> impl Parser<&'a [u8], Output = ParsedMolStandardLine, Error = error::Error<&'a [u8]>> {
    alt((
        map(counts_input(), ParsedMolStandardLine::Counts),
        map(atom_input_standard(), ParsedMolStandardLine::Atom),
        map(bond_input_standard(), |(_, _, bond)| {
            ParsedMolStandardLine::Bond(bond)
        }),
        map(tag("M  END"), |_| ParsedMolStandardLine::End),
        map(property_input_standard(), ParsedMolStandardLine::Property),
    ))
}

/// Progressively parse a MOL file, line by line, for diagnostic purposes.
///
/// This function processes one logical line at a time, handling multi-line blocks
/// like atom aliases, and returns a vector of `ParsedMolLine` enums representing
/// the content of each line. It uses the general parser that supports query features.
pub fn parse_mol_progressive(input: &[u8]) -> Vec<ParsedMolLine> {
    let mut results = Vec::new();
    let mut lines = input.lines().peekable();

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

    while let Some(line) = lines.next() {
        if line.trim().is_empty() {
            results.push(ParsedMolLine::Empty);
            continue;
        }

        if line.starts_with(b"A  ") {
            if let Some(next_line) = lines.next() {
                let combined_lines = join(b"\n", &[&line[3..], next_line]);
                results.push(
                    atom_alias_entry()
                        .parse(&combined_lines)
                        .map(|(_, alias_entry)| {
                            ParsedMolLine::Property(PropertyEntries::AtomAliasEntry(alias_entry))
                        })
                        .unwrap_or(ParsedMolLine::Unknown(line.to_str_lossy().to_string())),
                );
                continue;
            }
        }
        results.push(
            mol_input()
                .parse(line)
                .map(|(_, parsed)| parsed)
                .unwrap_or(ParsedMolLine::Unknown(line.to_str_lossy().to_string())),
        );
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
    let mut lines = input.lines().peekable();

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

    while let Some(line) = lines.next() {
        if line.trim().is_empty() {
            results.push(ParsedMolStandardLine::Empty);
            continue;
        }

        if line.starts_with(b"A  ") {
            if let Some(next_line) = lines.next() {
                let combined_lines = join(b"\n", &[&line[..3], next_line]);
                results.push(
                    atom_alias_entry()
                        .parse(&combined_lines)
                        .map(|(_, alias_entry)| {
                            ParsedMolStandardLine::Property(PropertyEntries::AtomAliasEntry(alias_entry))
                        })
                        .unwrap_or(ParsedMolStandardLine::Unknown(line.to_str_lossy().to_string())),
                );
                continue;
            }
        }

        results.push(
            mol_input_standard()
                .parse(line)
                .map(|(_, parsed)| parsed)
                .unwrap_or(ParsedMolStandardLine::Unknown(
                    line.to_str_lossy().to_string(),
                )),
        );
    }
    results
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
        assert!(
            result.is_ok(),
            "{} should have succeeded: {:?}",
            desc,
            result.err()
        );
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

    #[test]
    fn test_progressive_parser_standard() {
        let mol = b"MyMol\n\nComment\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n";
        let result = parse_mol_progressive_standard(mol);
        assert!(matches!(result[0], ParsedMolStandardLine::Header(_)));
        assert!(matches!(result[1], ParsedMolStandardLine::Header(_)));
        assert!(matches!(result[2], ParsedMolStandardLine::Header(_)));
        assert!(matches!(result[3], ParsedMolStandardLine::Counts(_)));
        assert!(matches!(result[4], ParsedMolStandardLine::Atom(_)));
        assert!(matches!(result[5], ParsedMolStandardLine::End));
    }

    #[test]
    fn test_progressive_parser_alias() {
        let mol_with_alias = b"MyMol\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nA    1\nCF3\nM  END\n";
        let result = parse_mol_progressive(mol_with_alias);
        println!("{:?}", result);
        assert!(
            matches!(
                result[5],
                ParsedMolLine::Property(PropertyEntries::AtomAliasEntry(_))
            ),
            "Failed to parse two-line atom alias. Got: {:?}",
            result[5]
        );
    }

    #[test]
    fn test_progressive_parser_legacy_list() {
        let mol_with_legacy = b"MyMol\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n  1 F    2   6   8\nM  END\n";
        let result = parse_mol_progressive(mol_with_legacy);
        assert!(
            matches!(
                result[5],
                ParsedMolLine::Property(PropertyEntries::AtomListEntry(_))
            ),
            "Failed to parse legacy atom list. Got: {:?}",
            result[5]
        );
    }

    #[rstest]
    #[case("src/io/mol/data/glycine-short-lines.mol", "glycine, short lines")]
    fn test_progressive_parser_standard_from_path(#[case] path: &str, #[case] desc: &str) {
        let input = std::fs::read(path).unwrap();
        let result = parse_mol_progressive_standard(&input);
        println!("{:?}", result);
        assert_eq!(result.len(), 10, "{} should have 10 lines", desc);
    }

}
