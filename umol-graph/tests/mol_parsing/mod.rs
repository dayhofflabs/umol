//! MOL file parsing conformance tests.
//!
//! This module runs parsing tests against a collection of MOL files from various sources.
//! Files are organized by parse result category:
//! - molecule: passes all 4 parsers
//! - molecule_lenient: needs lenient flags for basic parser
//! - extended_molecule: needs extended parser
//! - extended_molecule_lenient: needs lenient flags for extended parser
//! - invalid: fails all parsers

use std::fs::read;
use std::path::{Component, Path, PathBuf};

use insta::{assert_yaml_snapshot, Settings};
use rstest::*;
use serde::Serialize;
use umol_graph::io::ctfile::config::CtfileIoConfig;
use umol_graph::io::ctfile::error::ParseError;
use umol_graph::io::ctfile::parser::{parse_extended_mol_bytes_with, parse_mol_bytes_with};
use umol_graph::table_ir::{ExtendedMolecule, Molecule};

/// Category based on which parsers succeed
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
enum Category {
    Molecule,
    MoleculeLenient,
    ExtendedMolecule,
    ExtendedMoleculeLenient,
    Invalid,
    Bug,
}

impl Category {
    fn from_dir_name(name: &str) -> Option<Self> {
        match name {
            "molecule" => Some(Category::Molecule),
            "molecule_lenient" => Some(Category::MoleculeLenient),
            "extended_molecule" => Some(Category::ExtendedMolecule),
            "extended_molecule_lenient" => Some(Category::ExtendedMoleculeLenient),
            "invalid" => Some(Category::Invalid),
            "bug" => Some(Category::Bug),
            _ => None,
        }
    }

    fn from_parse_results(mol: bool, mol_len: bool, ext: bool, ext_len: bool) -> Self {
        match (mol, mol_len, ext, ext_len) {
            (true, true, true, true) => Category::Molecule,
            (false, true, false, true) => Category::MoleculeLenient,
            (false, false, true, true) => Category::ExtendedMolecule,
            (false, false, false, true) => Category::ExtendedMoleculeLenient,
            (false, false, false, false) => Category::Invalid,
            _ => Category::Bug,
        }
    }
}

/// Summary of molecule statistics for snapshot comparison
#[derive(Serialize)]
struct MoleculeSummary {
    sum_formula: String,
    atom_count: usize,
    bond_count: usize,
    property_count: usize,
}

/// Extended summary including extended features
#[derive(Serialize)]
struct ExtendedMoleculeSummary {
    sum_formula: String,
    atom_count: usize,
    bond_count: usize,
    extended_atoms: usize,
    extended_bonds: usize,
    property_count: usize,
    rgroup_count: usize,
    sgroup_count: usize,
}

/// Error summary for failed parses
#[derive(Serialize)]
struct ParseErrorSummary {
    error_type: String,
    message: String,
}

/// Result of a single parser invocation
#[derive(Serialize)]
struct ParseResult {
    success: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    summary: Option<MoleculeSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    extended_summary: Option<ExtendedMoleculeSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<ParseErrorSummary>,
}

/// Combined results from all 4 parsers for a single file
#[derive(Serialize)]
struct FileParseResults {
    expected_category: Category,
    category: Category,
    molecule: ParseResult,
    molecule_lenient: ParseResult,
    extended_molecule: ParseResult,
    extended_molecule_lenient: ParseResult,
}

impl From<&Molecule> for MoleculeSummary {
    fn from(mol: &Molecule) -> Self {
        Self {
            sum_formula: mol.sum_formula(),
            atom_count: mol.atom_count(),
            bond_count: mol.bond_count(),
            property_count: mol.property_count(),
        }
    }
}

impl From<&ExtendedMolecule> for ExtendedMoleculeSummary {
    fn from(mol: &ExtendedMolecule) -> Self {
        Self {
            sum_formula: mol.sum_formula(),
            atom_count: mol.atom_count(),
            bond_count: mol.bond_count(),
            extended_atoms: mol.extended_atom_count(),
            extended_bonds: mol.extended_bond_count(),
            property_count: mol.property_count(),
            rgroup_count: mol.rgroup_count(),
            sgroup_count: mol.sgroup_count(),
        }
    }
}

fn error_type_name(e: &ParseError) -> String {
    // Extract variant name from Debug representation
    let debug = format!("{:?}", e);
    // Take first word (variant name) before any { or (
    debug
        .split(['{', '('])
        .next()
        .unwrap_or(&debug)
        .trim()
        .to_string()
}

fn parse_with_basic_flags(bytes: &[u8]) -> ParseResult {
    let config = CtfileIoConfig::basic();
    match parse_mol_bytes_with(bytes, &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: Some(MoleculeSummary::from(&mol)),
            extended_summary: None,
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            extended_summary: None,
            error: Some(ParseErrorSummary {
                error_type: error_type_name(&e),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_with_lenient_flags(bytes: &[u8]) -> ParseResult {
    let config = CtfileIoConfig::basic_lenient();
    match parse_mol_bytes_with(bytes, &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: Some(MoleculeSummary::from(&mol)),
            extended_summary: None,
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            extended_summary: None,
            error: Some(ParseErrorSummary {
                error_type: error_type_name(&e),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_extended_with_extended_flags(bytes: &[u8]) -> ParseResult {
    let config = CtfileIoConfig::extended();
    match parse_extended_mol_bytes_with(bytes, &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: None,
            extended_summary: Some(ExtendedMoleculeSummary::from(&mol)),
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            extended_summary: None,
            error: Some(ParseErrorSummary {
                error_type: error_type_name(&e),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_extended_with_lenient_flags(bytes: &[u8]) -> ParseResult {
    let config = CtfileIoConfig::extended_lenient();
    match parse_extended_mol_bytes_with(bytes, &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: None,
            extended_summary: Some(ExtendedMoleculeSummary::from(&mol)),
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            extended_summary: None,
            error: Some(ParseErrorSummary {
                error_type: error_type_name(&e),
                message: e.to_string(),
            }),
        },
    }
}

fn extract_expected_category(path: &Path) -> Category {
    // Path structure: .../data/<category>/<source>/<file>.mol
    // We need to find the category directory (child of "data")
    let components: Vec<_> = path.components().collect();
    for (i, comp) in components.iter().enumerate() {
        if let Component::Normal(name) = comp {
            if name.to_str() == Some("data") && i + 1 < components.len() {
                if let Component::Normal(category_name) = &components[i + 1] {
                    if let Some(cat) = Category::from_dir_name(category_name.to_str().unwrap_or(""))
                    {
                        return cat;
                    }
                }
            }
        }
    }
    panic!(
        "Could not determine expected category from path: {:?}",
        path
    );
}

fn parse_file(path: &Path) -> FileParseResults {
    let bytes = read(path).expect("Failed to read file");
    let expected_category = extract_expected_category(path);

    let molecule = parse_with_basic_flags(&bytes);
    let molecule_lenient = parse_with_lenient_flags(&bytes);
    let extended_molecule = parse_extended_with_extended_flags(&bytes);
    let extended_molecule_lenient = parse_extended_with_lenient_flags(&bytes);

    let category = Category::from_parse_results(
        molecule.success,
        molecule_lenient.success,
        extended_molecule.success,
        extended_molecule_lenient.success,
    );

    FileParseResults {
        expected_category,
        category,
        molecule,
        molecule_lenient,
        extended_molecule,
        extended_molecule_lenient,
    }
}

fn run_conformance_test(file_path: &PathBuf) {
    let source = file_path
        .parent()
        .unwrap()
        .file_name()
        .unwrap()
        .to_str()
        .unwrap();
    let filename = file_path.file_name().unwrap().to_str().unwrap();

    let results = parse_file(file_path);

    // Validate category matches expectation
    assert_eq!(
        results.expected_category, results.category,
        "Category mismatch for {:?}: expected {:?}, got {:?}",
        file_path, results.expected_category, results.category
    );

    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("mol_parsing");

    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(base.join("snapshots"));
    settings.set_snapshot_suffix(format!("{}_{}", source, filename));
    settings.bind(|| {
        assert_yaml_snapshot!(results);
    });
}

#[rstest]
fn test_conformance(#[files("tests/mol_parsing/data/**/*.mol")] file_path: PathBuf) {
    run_conformance_test(&file_path);
}
