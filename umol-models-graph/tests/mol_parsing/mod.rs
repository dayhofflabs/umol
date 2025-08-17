use insta::assert_yaml_snapshot;
use serde::Serialize;
use std::fs::read;
use std::path::Path;

use umol_models_graph::io::mol::parser::parse_mol_file;

#[derive(Serialize)]
struct Snapshot<'a> {
    category: &'a str,
    filename: &'a str,
    expected_success: bool,
    success: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    molecule: Option<serde_yaml::Value>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<String>,
}

fn run_test(path: &Path) {
    let mol_bytes = read(path).unwrap();
    let result = parse_mol_file(&mol_bytes);

    let path_str = path.to_str().unwrap();
    let expected_success = path_str.contains("/valid/");
    let category = path
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .file_name()
        .unwrap()
        .to_str()
        .unwrap();

    let snapshot = match result {
        Ok(mol_file) => Snapshot {
            category,
            filename: path.file_name().unwrap().to_str().unwrap(),
            expected_success,
            success: true,
            molecule: Some(serde_yaml::to_value(&mol_file.molecule).unwrap()),
            error: None,
        },
        Err(e) => Snapshot {
            category,
            filename: path.file_name().unwrap().to_str().unwrap(),
            expected_success,
            success: false,
            molecule: None,
            error: Some(e.to_string()),
        },
    };

    assert_yaml_snapshot!(snapshot);
}

#[test]
fn test_basic() {
    insta::glob!("data/basic/valid/*.mol", run_test);
}

#[test]
fn test_properties() {
    insta::glob!("data/with_properties/valid/*.mol", run_test);
}

#[test]
fn test_query() {
    insta::glob!("data/query/valid/*.mol", run_test);
}

#[test]
fn test_rgroups() {
    insta::glob!("data/with_rgroups/valid/*.mol", run_test);
}

#[test]
fn test_sgroups() {
    insta::glob!("data/with_sgroups/valid/*.mol", run_test);
}

#[test]
fn test_edge_cases() {
    insta::glob!("data/uncategorized/*.mol", run_test);
}
