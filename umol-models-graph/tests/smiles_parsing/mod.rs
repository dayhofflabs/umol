//! SMILES parsing conformance tests.
//!
//! This module runs parsing tests against individual SMILES files from various sources.
//! Each .smiles file contains a single SMILES string (with optional comment header).
//!
//! Files are organized by parse result category:
//! - opensmiles_strict: passes strict OpenSMILES parser
//! - invalid: fails parser

use std::fs;
use std::path::{Path, PathBuf};

use insta::{assert_yaml_snapshot, Settings};
use rstest::rstest;
use serde::Serialize;
use umol_models_graph::io::smiles::error::ParseError;
use umol_models_graph::io::smiles::parse_smiles;
use umol_models_graph::table_ir::Molecule;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
enum Category {
    OpensmilesStrict,
    Invalid,
    #[allow(dead_code)] // Placeholder for future configs
    Bug,
}

impl Category {
    fn from_dir_name(name: &str) -> Option<Self> {
        match name {
            "opensmiles_strict" => Some(Category::OpensmilesStrict),
            "invalid" => Some(Category::Invalid),
            _ => None,
        }
    }

    fn from_parse_result(strict_ok: bool) -> Self {
        if strict_ok {
            Category::OpensmilesStrict
        } else {
            Category::Invalid
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
    opensmiles_strict: ParseResult,
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

fn parse_with_strict(smiles: &str) -> ParseResult {
    match parse_smiles(smiles.as_bytes()) {
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

    let opensmiles_strict = parse_with_strict(&smiles);
    let category = Category::from_parse_result(opensmiles_strict.success);

    FileParseResults {
        expected_category,
        category,
        opensmiles_strict,
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
