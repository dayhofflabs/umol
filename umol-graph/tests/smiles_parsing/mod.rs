//! SMILES parsing conformance tests.
//!
//! This module runs parsing tests against individual SMILES files from various sources.
//! Each .smiles file contains a single SMILES string (with optional comment header).
//!
//! Files are organized by parse result category:
//! - basic_opensmiles: passes strict OpenSMILES parser (no wildcards)
//! - opensmiles: passes extended OpenSMILES parser (with wildcards)
//! - basic_chemaxon: requires CHEMAXON_EXTENSIONS (basic parser)
//! - chemaxon: requires CHEMAXON_EXTENSIONS (extended parser)
//! - chemaxon_invalid: SMILES part parses, but CX block is invalid/unhandled
//! - invalid: fails all parsers

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::LazyLock;

use insta::{assert_yaml_snapshot, Settings};
use regex::Regex;
use rstest::*;
use serde::Serialize;
use umol_graph::io::smiles::config::SmilesIoConfig;
use umol_graph::io::smiles::error::ParseError;
use umol_graph::io::smiles::{
    parse_extended_smiles_bytes_with, parse_smiles, parse_smiles_bytes_with,
};
use umol_graph::table_ir::{ExtendedMolecule, Molecule};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
enum Category {
    BasicOpensmiles,
    Opensmiles,
    BasicChemaxon,
    Chemaxon,
    Invalid,
    ChemaxonInvalid,
    Bug,
}

impl Category {
    fn from_dir_name(name: &str) -> Option<Self> {
        match name {
            "basic_opensmiles" => Some(Category::BasicOpensmiles),
            "opensmiles" => Some(Category::Opensmiles),
            "basic_chemaxon" => Some(Category::BasicChemaxon),
            "chemaxon" => Some(Category::Chemaxon),
            "invalid" => Some(Category::Invalid),
            "chemaxon_invalid" => Some(Category::ChemaxonInvalid),
            "bug" => Some(Category::Bug),
            _ => None,
        }
    }

    fn from_parse_result(
        has_cx: bool,
        basic_ok: bool,
        opensmiles_ok: bool,
        basic_chemaxon_ok: bool,
        chemaxon_ok: bool,
    ) -> Self {
        match (
            has_cx,
            basic_ok,
            opensmiles_ok,
            basic_chemaxon_ok,
            chemaxon_ok,
        ) {
            (false, true, true, true, true) => Category::BasicOpensmiles,
            (true, true, true, true, true) => Category::BasicChemaxon,
            (false, false, true, false, true) => Category::Opensmiles,
            (true, true, true, false, true) | (true, false, true, false, true) => {
                Category::Chemaxon
            }
            (false, false, false, false, false) | (true, false, false, false, false) => {
                Category::Invalid
            }
            (true, true, true, false, false) | (true, false, true, false, false) => {
                Category::ChemaxonInvalid
            }

            // Anything else is either a hierarchy violation or a CX-vs-SMILES inconsistency.
            _ => Category::Bug,
        }
    }
}

#[derive(Serialize)]
struct MoleculeSummary {
    sum_formula: String,
    atom_count: usize,
    bond_count: usize,
}

impl From<&Molecule> for MoleculeSummary {
    fn from(mol: &Molecule) -> Self {
        Self {
            sum_formula: mol.sum_formula(),
            atom_count: mol.atom_count(),
            bond_count: mol.bond_count(),
        }
    }
}

impl From<&ExtendedMolecule> for MoleculeSummary {
    fn from(mol: &ExtendedMolecule) -> Self {
        Self {
            sum_formula: mol.sum_formula(),
            atom_count: mol.atom_count(),
            bond_count: mol.bond_count(),
        }
    }
}

#[derive(Serialize)]
struct ParseErrorSummary {
    error_type: String,
    message: String,
}

#[derive(Serialize)]
struct ParseResult {
    success: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    summary: Option<MoleculeSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<ParseErrorSummary>,
}

#[derive(Serialize)]
struct FileParseResults {
    expected_category: Category,
    category: Category,
    basic_opensmiles: ParseResult,
    #[serde(skip_serializing_if = "Option::is_none")]
    opensmiles: Option<ParseResult>,
}

fn error_type_name(e: &ParseError) -> String {
    let debug = format!("{:?}", e);
    debug
        .split(|c| c == '{' || c == '(')
        .next()
        .unwrap_or(&debug)
        .trim()
        .to_string()
}

fn extract_expected_category(path: &Path) -> Category {
    let components: Vec<_> = path.components().collect();
    for (i, comp) in components.iter().enumerate() {
        if let std::path::Component::Normal(name) = comp {
            if name.to_str() == Some("data") && i + 1 < components.len() {
                if let std::path::Component::Normal(category_name) = &components[i + 1] {
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

fn read_smiles_from_file(path: &Path) -> String {
    let content = fs::read_to_string(path).expect("Failed to read file");
    content
        .lines()
        .find(|line| !line.trim().is_empty() && !line.trim().starts_with('#'))
        .map(|s| s.trim().to_string())
        .unwrap_or_default()
}

fn has_cx_extensions(smiles: &str) -> bool {
    static CX_ANNOTATIONS_RE: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"^\S+\s+\|.*\|").expect("CX annotations regex"));
    CX_ANNOTATIONS_RE.is_match(smiles)
}

fn parse_with_strict(smiles: &str) -> ParseResult {
    match parse_smiles(smiles) {
        Ok(mol) => ParseResult {
            success: true,
            summary: Some(MoleculeSummary::from(&mol)),
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            error: Some(ParseErrorSummary {
                error_type: error_type_name(&e),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_with_opensmiles(smiles: &str) -> ParseResult {
    let config = SmilesIoConfig::opensmiles();
    match parse_extended_smiles_bytes_with(smiles.as_bytes(), &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: Some(MoleculeSummary::from(&mol)),
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            error: Some(ParseErrorSummary {
                error_type: error_type_name(&e),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_with_basic_chemaxon(smiles: &str) -> ParseResult {
    let config = SmilesIoConfig::basic_chemaxon();
    match parse_smiles_bytes_with(smiles.as_bytes(), &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: Some(MoleculeSummary::from(&mol)),
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            error: Some(ParseErrorSummary {
                error_type: error_type_name(&e),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_with_chemaxon(smiles: &str) -> ParseResult {
    let config = SmilesIoConfig::chemaxon();
    match parse_extended_smiles_bytes_with(smiles.as_bytes(), &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: Some(MoleculeSummary::from(&mol)),
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            error: Some(ParseErrorSummary {
                error_type: error_type_name(&e),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_file(path: &Path) -> FileParseResults {
    let expected_category = extract_expected_category(path);
    let smiles = read_smiles_from_file(path);
    let has_cx = has_cx_extensions(&smiles);

    // Run all parsers to verify parser hierarchy
    let basic_opensmiles = parse_with_strict(&smiles);
    let opensmiles_result = parse_with_opensmiles(&smiles);
    let basic_chemaxon_result = parse_with_basic_chemaxon(&smiles);
    let chemaxon_result = parse_with_chemaxon(&smiles);

    let category = Category::from_parse_result(
        has_cx,
        basic_opensmiles.success,
        opensmiles_result.success,
        basic_chemaxon_result.success,
        chemaxon_result.success,
    );

    // Only include extended result in output if it differs from basic
    let opensmiles = if basic_opensmiles.success {
        None
    } else {
        Some(opensmiles_result)
    };

    FileParseResults {
        expected_category,
        category,
        basic_opensmiles,
        opensmiles,
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

    assert_eq!(
        results.expected_category, results.category,
        "Category mismatch for {:?}: expected {:?}, got {:?}",
        file_path, results.expected_category, results.category
    );

    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("smiles_parsing");

    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(base.join("snapshots"));
    settings.set_snapshot_suffix(format!("{}_{}", source, filename));
    settings.bind(|| {
        assert_yaml_snapshot!(results);
    });
}

#[rstest]
fn test_conformance(#[files("tests/smiles_parsing/data/**/*.smiles")] file_path: PathBuf) {
    run_conformance_test(&file_path);
}
