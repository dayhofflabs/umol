use std::path::{Path, PathBuf};

use insta::{assert_yaml_snapshot, Settings};
use rstest::*;
use serde::Serialize;
use umol_models_graph::io::ctab::molecule::{Molecule, MoleculeLike};
use umol_models_graph::io::mol::parser::{parse_mol, parse_mol_moleculelike};

#[derive(Serialize)]
struct MoleculeSummary {
    atom_count: usize,
    bond_count: usize,
    sum_formula: String,
    graph6: String,
}

impl From<&MoleculeLike> for MoleculeSummary {
    fn from(molecule: &MoleculeLike) -> Self {
        Self {
            atom_count: molecule.atom_count(),
            bond_count: molecule.bond_count(),
            sum_formula: molecule.sum_formula(),
            graph6: molecule.graph6(),
        }
    }
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

fn store_summary_snapshot(summary: MoleculeSummary) {
    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(Path::new("snapshots"));
    settings.bind(|| {
        assert_yaml_snapshot!(summary);
    });
}

fn test_parse_mol_moleculelike(path: &Path, expected_success: bool) {
    let mol_bytes = std::fs::read(path).unwrap();
    let result = parse_mol_moleculelike(&mol_bytes);

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
    let result = parse_mol(&mol_bytes);

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

fn run_parse_mol_moleculelike_test(file_path: &Path, expected_success: bool) {
    let (source, filename) = extract_source_and_filename(file_path);
    let mut settings = insta::Settings::clone_current();
    settings.set_snapshot_suffix(format!("parse_mol_moleculelike_{}_{}", source, filename));
    let _guard = settings.bind_to_scope();
    test_parse_mol_moleculelike(file_path, expected_success);
}

// parse_mol tests: should succeed on molecule/*, fail on invalid/*

// Pure basic sources (100% basic) - should succeed
#[rstest]
fn test_parse_mol_molecule_jmol(
    #[files("tests/mol_parsing/data/molecule/jmol/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_molecule_nist(
    #[files("tests/mol_parsing/data/molecule/nist/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

// High-volume basic sources - should succeed
#[rstest]
fn test_parse_mol_molecule_rdkit(
    #[files("tests/mol_parsing/data/molecule/rdkit/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_molecule_openbabel(
    #[files("tests/mol_parsing/data/molecule/openbabel/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_molecule_cdk(
    #[files("tests/mol_parsing/data/molecule/cdk/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_molecule_chemspider(
    #[files("tests/mol_parsing/data/molecule/chemspider/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_molecule_indigo(
    #[files("tests/mol_parsing/data/molecule/indigo/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_molecule_ketcher(
    #[files("tests/mol_parsing/data/molecule/ketcher/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_molecule_marvin(
    #[files("tests/mol_parsing/data/molecule/marvin/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_molecule_reaxys(
    #[files("tests/mol_parsing/data/molecule/reaxys/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_molecule_scifinder(
    #[files("tests/mol_parsing/data/molecule/scifinder/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_molecule_chebi(
    #[files("tests/mol_parsing/data/molecule/chebi/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, true);
}

// Invalid files (should fail)
#[rstest]
fn test_parse_mol_invalid_chembl(
    #[files("tests/mol_parsing/data/invalid/chembl/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_invalid_indigo(
    #[files("tests/mol_parsing/data/invalid/indigo/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_invalid_ketcher(
    #[files("tests/mol_parsing/data/invalid/ketcher/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_invalid_rdkit(
    #[files("tests/mol_parsing/data/invalid/rdkit/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_test(&file_path, false);
}

// parse_mol_moleculelike tests: should succeed on molecule/* + moleculelike/*, fail on invalid/*

// Basic files (should also work with extended parser) - should succeed
#[rstest]
fn test_parse_mol_moleculelike_molecule_jmol(
    #[files("tests/mol_parsing/data/molecule/jmol/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_nist(
    #[files("tests/mol_parsing/data/molecule/nist/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_rdkit(
    #[files("tests/mol_parsing/data/molecule/rdkit/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_openbabel(
    #[files("tests/mol_parsing/data/molecule/openbabel/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_cdk(
    #[files("tests/mol_parsing/data/molecule/cdk/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_chemspider(
    #[files("tests/mol_parsing/data/molecule/chemspider/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_indigo(
    #[files("tests/mol_parsing/data/molecule/indigo/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_ketcher(
    #[files("tests/mol_parsing/data/molecule/ketcher/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_marvin(
    #[files("tests/mol_parsing/data/molecule/marvin/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_reaxys(
    #[files("tests/mol_parsing/data/molecule/reaxys/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_scifinder(
    #[files("tests/mol_parsing/data/molecule/scifinder/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_molecule_chebi(
    #[files("tests/mol_parsing/data/molecule/chebi/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

// Extended feature files (only work with extended parser) - should succeed
#[rstest]
fn test_parse_mol_moleculelike_moleculelike_chebi(
    #[files("tests/mol_parsing/data/moleculelike/chebi/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_moleculelike_rhea(
    #[files("tests/mol_parsing/data/moleculelike/rhea/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_moleculelike_rdkit(
    #[files("tests/mol_parsing/data/moleculelike/rdkit/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_moleculelike_cdk(
    #[files("tests/mol_parsing/data/moleculelike/cdk/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_moleculelike_indigo(
    #[files("tests/mol_parsing/data/moleculelike/indigo/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_moleculelike_ketcher(
    #[files("tests/mol_parsing/data/moleculelike/ketcher/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_moleculelike_marvin(
    #[files("tests/mol_parsing/data/moleculelike/marvin/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_moleculelike_scifinder(
    #[files("tests/mol_parsing/data/moleculelike/scifinder/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

#[rstest]
fn test_parse_mol_moleculelike_moleculelike_openbabel(
    #[files("tests/mol_parsing/data/moleculelike/openbabel/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, true);
}

// Invalid files (should fail with extended parser too)
#[rstest]
fn test_parse_mol_moleculelike_invalid_chembl(
    #[files("tests/mol_parsing/data/invalid/chembl/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_moleculelike_invalid_indigo(
    #[files("tests/mol_parsing/data/invalid/indigo/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_moleculelike_invalid_ketcher(
    #[files("tests/mol_parsing/data/invalid/ketcher/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_moleculelike_invalid_rdkit(
    #[files("tests/mol_parsing/data/invalid/rdkit/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_moleculelike_invalid_marvin(
    #[files("tests/mol_parsing/data/invalid/marvin/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_moleculelike_invalid_cdk(
    #[files("tests/mol_parsing/data/invalid/cdk/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_moleculelike_invalid_chemspider(
    #[files("tests/mol_parsing/data/invalid/chemspider/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_moleculelike_invalid_rhea(
    #[files("tests/mol_parsing/data/invalid/rhea/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, false);
}

#[rstest]
fn test_parse_mol_moleculelike_invalid_openbabel(
    #[files("tests/mol_parsing/data/invalid/openbabel/*.mol")] file_path: PathBuf,
) {
    run_parse_mol_moleculelike_test(&file_path, false);
}
