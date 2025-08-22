use std::path::Path;

use insta::{assert_yaml_snapshot, Settings};
use serde::Serialize;

#[allow(dead_code)]
enum TestMode {
    Summary,
    Full,
}

#[derive(Serialize)]
struct SdfSummary<'a> {
    category: &'a str,
    filename: &'a str,
    compound_count: usize,
    compounds: Vec<CompoundSummary>,
}

#[derive(Serialize)]
struct CompoundSummary {
    sum_formula: String,
    atom_count: usize,
    bond_count: usize,
    graph6: String,
    data_fields: Vec<(String, String)>,
}

#[derive(Serialize)]
struct SdfFullSnapshot<'a> {
    category: &'a str,
    filename: &'a str,
    compounds: Vec<CompoundFull>,
}

#[derive(Serialize)]
struct CompoundFull {
    mol_data: serde_yaml::Value,
    data_fields: Vec<(String, String)>,
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

fn run_test_sdf(path: &Path, mode: TestMode) {
    let path_str = path.to_str().unwrap();
    let category = get_category(path);
    let sdf_bytes = std::fs::read(path).unwrap();
    
    // Parse SDF file
    let sdf_file = match umol_models_graph::io::sdf::parse_sdf(&sdf_bytes) {
        Ok(sdf) => sdf,
        Err(e) => panic!("Failed to parse SDF file {}: {}", path_str, e),
    };

    match mode {
        TestMode::Summary => {
            let mut compound_summaries = Vec::new();
            
            for compound in &sdf_file.compounds {
                let molecule = &compound.mol_file.molecule;
                
                compound_summaries.push(CompoundSummary {
                    sum_formula: molecule.sum_formula(),
                    atom_count: molecule.atom_count(),
                    bond_count: molecule.bond_count(),
                    graph6: molecule.graph6(),
                    data_fields: compound.data_fields.clone(),
                });
            }
            
            let summary = SdfSummary {
                category,
                filename: path.file_name().unwrap().to_str().unwrap(),
                compound_count: sdf_file.compounds.len(),
                compounds: compound_summaries,
            };
            
            store_summary_snapshot(summary);
        },
        TestMode::Full => {
            let mut compound_fulls = Vec::new();
            
            for compound in &sdf_file.compounds {
                let molecule = &compound.mol_file.molecule;
                
                compound_fulls.push(CompoundFull {
                    mol_data: serde_yaml::to_value(&molecule).expect("Failed to serialize molecule"),
                    data_fields: compound.data_fields.clone(),
                });
            }
            
            let full_snapshot = SdfFullSnapshot {
                category,
                filename: path.file_name().unwrap().to_str().unwrap(),
                compounds: compound_fulls,
            };
            
            store_full_snapshot(path, category, full_snapshot);
        }
    }
}

fn store_summary_snapshot(summary: SdfSummary) {
    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(Path::new("snapshots").join("summary"));
    settings.bind(|| {
        assert_yaml_snapshot!(summary);
    });
}

fn store_full_snapshot(_path: &Path, _category: &str, snapshot: SdfFullSnapshot) {
    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(Path::new("snapshots").join("full"));
    settings.bind(|| {
        assert_yaml_snapshot!(snapshot);
    });
}

fn run_test_invalid_sdf(path: &std::path::Path) {
    let path_str = path.to_str().unwrap();
    let sdf_bytes = std::fs::read(path).unwrap();
    
    // These files should fail to parse
    let result = umol_models_graph::io::sdf::parse_sdf(&sdf_bytes);
    assert!(result.is_err(), "Invalid SDF file {} should fail to parse", path_str);
    
    // Create snapshot of the error for regression testing
    let error_msg = format!("{}", result.unwrap_err());
    insta::with_settings!({
        description => format!("Error for invalid SDF file: {}", path.file_name().unwrap().to_str().unwrap()),
    }, {
        insta::assert_snapshot!(format!("invalid_sdf_error_{}", path.file_stem().unwrap().to_str().unwrap()), error_msg);
    });
}

#[test]
fn test_summaries() {
    insta::glob!("data/*.sdf", |path| {
        run_test_sdf(path, TestMode::Summary);
    });
}

#[test]
fn test_full() {
    insta::glob!("data/*.sdf", |path| {
        run_test_sdf(path, TestMode::Full);
    });
}

#[test] 
fn test_invalid_files() {
    // Test that known invalid SDF files fail for the right reasons
    insta::glob!("data/invalid/*.sdf", |path| {
        run_test_invalid_sdf(path);
    });
}