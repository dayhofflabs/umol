use std::path::{Path, PathBuf};

use insta::{assert_yaml_snapshot, Settings};
use rstest::*;
use serde::Serialize;
use umol_models_graph::io::ctab::molecule::{Molecule, MoleculeStandard};
use umol_models_graph::io::mol::parser::{parse_mol_file, parse_mol_file_standard};

#[allow(dead_code)]
enum TestMode {
    Summary,
    Full,
    // SGroups, etc. would go here
}

#[derive(Serialize)]
struct MoleculeSummary {
    atom_count: usize,
    bond_count: usize,
    sum_formula: String,
    graph6: String,
}

impl From<&Molecule> for MoleculeSummary {
    fn from(molecule: &Molecule) -> Self {
        Self {
            atom_count: molecule.atom_count(),
            bond_count: molecule.bond_count(),
            sum_formula: molecule.sum_formula(),
            graph6: molecule.graph6(),
        }
    }
}

impl From<&MoleculeStandard> for MoleculeSummary {
    fn from(molecule: &MoleculeStandard) -> Self {
        Self {
            atom_count: molecule.atom_count(),
            bond_count: molecule.bond_count(),
            sum_formula: molecule.sum_formula(),
            graph6: molecule.graph6(),
        }
    }
}

#[derive(Serialize)]
struct FullSnapshot<'a> {
    category: &'a str,
    filename: &'a str,
    molecule: serde_yaml::Value,
}

fn get_category(path: &Path) -> &str {
    path.parent()
        .unwrap()
        .parent()
        .unwrap()
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
}

fn store_summary_snapshot(summary: MoleculeSummary) {
    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(Path::new("snapshots").join("summary"));
    settings.bind(|| {
        assert_yaml_snapshot!(summary);
    });
}

fn store_full_snapshot(path: &Path, category: &str, value: serde_yaml::Value) {
    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(Path::new("snapshots").join("full"));
    settings.bind(|| {
        let snapshot = FullSnapshot {
            category,
            filename: path.file_name().unwrap().to_str().unwrap(),
            molecule: value,
        };
        assert_yaml_snapshot!(snapshot);
    });
}

fn run_test(path: &Path, mode: TestMode) {
    let path_str = path.to_str().unwrap();
    let category = get_category(path);
    let expected_success = path_str.contains("/valid/") || path_str.contains("data/");

    let mol_bytes = std::fs::read(path).unwrap();
    let result = parse_mol_file(&mol_bytes);
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
    } else {
        assert!(
            result.is_err(),
            "Expected parsing to fail, but it succeeded for file: {}",
            path.display()
        );
        return;
    }

    match mode {
        TestMode::Summary => {
            let summary = MoleculeSummary::from(&result.unwrap().molecule);
            store_summary_snapshot(summary);
        }
        TestMode::Full => {
            let serialized = serde_yaml::to_value(&result.unwrap().molecule).unwrap();
            store_full_snapshot(path, category, serialized);
        }
    }
}

fn run_test_standard(path: &Path, mode: TestMode) {
    let path_str = path.to_str().unwrap();
    let category = get_category(path);
    let expected_success = path_str.contains("/valid/") || path_str.contains("data/");

    let mol_bytes = std::fs::read(path).unwrap();

    let result = parse_mol_file_standard(&mol_bytes);
    if expected_success {
        assert!(
            result.is_ok(),
            "Expected parsing to succeed, but it failed for file: {}",
            path.display()
        );
    } else {
        assert!(
            result.is_err(),
            "Expected parsing to fail, but it succeeded for file: {}",
            path.display()
        );
        return;
    }
    match mode {
        TestMode::Summary => {
            let summary = MoleculeSummary::from(&result.unwrap().molecule);
            store_summary_snapshot(summary);
        }
        TestMode::Full => {
            let serialized = serde_yaml::to_value(&result.unwrap().molecule).unwrap();
            store_full_snapshot(path, category, serialized);
        }
    }
}

// Makes snapshot files unique for each input file
macro_rules! set_snapshot_suffix {
    ($($expr:expr),*) => {
        let mut settings = insta::Settings::clone_current();
        settings.set_snapshot_suffix(format!($($expr,)*));
        let _guard = settings.bind_to_scope();
    }
}

#[rstest]
fn test_summary_basic_standard(#[files("tests/mol_parsing/data/basic/*.mol")] file_path: PathBuf) {
    set_snapshot_suffix!("{}", file_path.file_name().unwrap().to_str().unwrap());
    run_test_standard(&file_path, TestMode::Summary);
}

#[rstest]
fn test_summary_basic(#[files("tests/mol_parsing/data/basic/*.mol")] file_path: PathBuf) {
    set_snapshot_suffix!("{}", file_path.file_name().unwrap().to_str().unwrap());
    run_test(&file_path, TestMode::Summary);
}

#[rstest]
fn test_summary_properties(#[files("tests/mol_parsing/data/properties/*.mol")] file_path: PathBuf) {
    set_snapshot_suffix!("{}", file_path.file_name().unwrap().to_str().unwrap());
    run_test(&file_path, TestMode::Summary);
}

#[rstest]
fn test_summary_query(#[files("tests/mol_parsing/data/query/*.mol")] file_path: PathBuf) {
    set_snapshot_suffix!("{}", file_path.file_name().unwrap().to_str().unwrap());
    run_test(&file_path, TestMode::Summary);
}

#[rstest]
fn test_summary_rgroups(#[files("tests/mol_parsing/data/rgroups/*.mol")] file_path: PathBuf) {
    set_snapshot_suffix!("{}", file_path.file_name().unwrap().to_str().unwrap());
    run_test(&file_path, TestMode::Summary);
}

#[rstest]
fn test_summary_sgroups(#[files("tests/mol_parsing/data/sgroups/*.mol")] file_path: PathBuf) {
    set_snapshot_suffix!("{}", file_path.file_name().unwrap().to_str().unwrap());
    run_test(&file_path, TestMode::Summary);
}
