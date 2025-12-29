use std::fs;
use std::path::{Path, PathBuf};

use indexmap::IndexMap;
use insta::{assert_yaml_snapshot, Settings};
use rstest::*;
use serde::Serialize;
use umol_models_graph::io::sdf::parser::parse_sdf_bytes;

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
    data_fields: IndexMap<String, String>,
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

fn run_test_sdf(path: &Path) {
    let path_str = path.to_str().unwrap();
    let category = get_category(path);
    let sdf_bytes = fs::read(path).unwrap();

    let sdf_file = match parse_sdf_bytes(&sdf_bytes) {
        Ok(sdf) => sdf,
        Err(e) => panic!("Failed to parse SDF file {}: {}", path_str, e),
    };

    let mut compound_summaries = Vec::new();

    for compound in &sdf_file.compounds {
        let molecule = &compound.mol_file.molecule;

        compound_summaries.push(CompoundSummary {
            sum_formula: molecule.sum_formula(),
            atom_count: molecule.atom_count(),
            bond_count: molecule.bond_count(),
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
}

fn store_summary_snapshot(summary: SdfSummary) {
    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(Path::new("snapshots").join("summary"));
    settings.bind(|| {
        assert_yaml_snapshot!(summary);
    });
}


#[allow(dead_code)]
fn run_test_invalid_sdf(path: &Path) {
    let path_str = path.to_str().unwrap();
    let sdf_bytes = fs::read(path).unwrap();

    let result = parse_sdf_bytes(&sdf_bytes);
    assert!(
        result.is_err(),
        "Invalid SDF file {} should fail to parse",
        path_str
    );

    let error_msg = format!("{}", result.unwrap_err());
    insta::with_settings!({
        description => format!("Error for invalid SDF file: {}", path.file_name().unwrap().to_str().unwrap()),
    }, {
        insta::assert_snapshot!(format!("invalid_sdf_error_{}", path.file_stem().unwrap().to_str().unwrap()), error_msg);
    });
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
fn test_summary(#[files("tests/sdf_parsing/data/*.sdf")] file_path: PathBuf) {
    set_snapshot_suffix!("{}", file_path.file_name().unwrap().to_str().unwrap());
    run_test_sdf(&file_path);
}
