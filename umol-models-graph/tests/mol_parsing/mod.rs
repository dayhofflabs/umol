use std::path::{Path, PathBuf};

use insta::{assert_yaml_snapshot, Settings};
use rstest::*;
use serde::Serialize;
use umol_models_graph::io::mol::parser::{parse_extended_mol_bytes, parse_mol_bytes};
use umol_models_graph::table_ir::{ExtendedMolecule, Molecule};

#[derive(Serialize)]
struct MoleculeSummary {
    atom_count: usize,
    bond_count: usize,
    sum_formula: String,
}

impl From<&ExtendedMolecule> for MoleculeSummary {
    fn from(molecule: &ExtendedMolecule) -> Self {
        Self {
            atom_count: molecule.atom_count(),
            bond_count: molecule.bond_count(),
            sum_formula: molecule.sum_formula(),
        }
    }
}

impl From<&Molecule> for MoleculeSummary {
    fn from(molecule: &Molecule) -> Self {
        Self {
            atom_count: molecule.atom_count(),
            bond_count: molecule.bond_count(),
            sum_formula: molecule.sum_formula(),
        }
    }
}

fn store_summary_snapshot(summary: MoleculeSummary) {
    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(Path::new("snapshots"));
    settings.bind(|| {
        assert_yaml_snapshot!(summary);
    });
}

fn test_parse_extended_mol(path: &Path, expected_success: bool) {
    let mol_bytes = std::fs::read(path).unwrap();
    let result = parse_extended_mol_bytes(&mol_bytes);

    if expected_success {
        if let Err(e) = &result {
            eprintln!("Parse error for {}: {:?}", path.display(), e);
        }
        assert!(
            result.is_ok(),
            "Expected parsing to succeed, but it failed for file: {} with error: {:?}",
            path.display(),
            result.as_ref().err()
        );

        let summary = MoleculeSummary::from(&result.unwrap());
        store_summary_snapshot(summary);
    } else {
        assert!(
            result.is_err(),
            "Expected parsing to fail, but it succeeded for file: {}",
            path.display()
        );
    }
}

fn test_parse_mol(path: &Path, expected_success: bool) {
    let mol_bytes = std::fs::read(path).unwrap();
    let result = parse_mol_bytes(&mol_bytes);

    if expected_success {
        if let Err(ref e) = result {
            eprintln!("Parse error for {}: {:?}", path.display(), e);
        }
        assert!(
            result.is_ok(),
            "Expected parsing to succeed, but it failed for file: {} with error: {:?}",
            path.display(),
            result.as_ref().err()
        );

        let summary = MoleculeSummary::from(&result.unwrap());
        store_summary_snapshot(summary);
    } else {
        assert!(
            result.is_err(),
            "Expected parsing to fail, but it succeeded for file: {}",
            path.display()
        );
    }
}

// Helper functions for test organization
fn extract_source_and_filename(path: &Path) -> (String, String) {
    let source = path
        .parent()
        .unwrap()
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();
    let filename = path.file_name().unwrap().to_str().unwrap().to_string();
    (source, filename)
}

fn run_parse_mol_test(file_path: &Path, expected_success: bool) {
    let (source, filename) = extract_source_and_filename(file_path);
    let mut settings = insta::Settings::clone_current();
    settings.set_snapshot_suffix(format!("parse_mol_{}_{}", source, filename));
    let _guard = settings.bind_to_scope();
    test_parse_mol(file_path, expected_success);
}

fn run_parse_extended_mol_test(file_path: &Path, expected_success: bool) {
    let (source, filename) = extract_source_and_filename(file_path);
    let mut settings = insta::Settings::clone_current();
    settings.set_snapshot_suffix(format!("parse_extended_mol_{}_{}", source, filename));
    let _guard = settings.bind_to_scope();
    test_parse_extended_mol(file_path, expected_success);
}

// AUTO-GENERATED TESTS - managed by build.rs
// To add a new source:
// 1. Add files to tests/mol_parsing/data_raw/newsource/
// 2. Run: cargo run --bin mol_classifier -- --sort
// 3. Tests will be automatically generated based on the organized structure

include!(concat!(env!("OUT_DIR"), "/generated_tests.rs"));
