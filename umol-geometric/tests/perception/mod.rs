use std::fs;
use std::path::{Path, PathBuf};

use insta::{assert_yaml_snapshot, Settings};
use rstest::rstest;
use serde::Serialize;
use umol_geometric::io::xyz::parse_xyz;
use umol_geometric::molecule::Molecule;
use umol_msym::Thresholds;

#[derive(Serialize)]
struct PerceptionResult {
    point_group: String,
    order: usize,
    equivalence_sets: Vec<Vec<String>>,
}

fn perceive(mol: &Molecule) -> PerceptionResult {
    let sym = mol.perceive_symmetry(Thresholds::default()).unwrap();

    let equivalence_sets: Vec<Vec<String>> = sym
        .equivalence_sets()
        .iter()
        .map(|set| {
            set.iter()
                .map(|&i| format!("{}({})", sym.element(i).symbol(), i))
                .collect()
        })
        .collect();

    PerceptionResult {
        point_group: sym.point_group().to_string(),
        order: sym.point_group().order(),
        equivalence_sets,
    }
}

fn run_perception_test(file_path: &Path) {
    let content = fs::read_to_string(file_path).unwrap();
    let (mol, _comment) = parse_xyz(&content).unwrap();
    let result = perceive(&mol);

    let filename = file_path.file_stem().unwrap().to_str().unwrap();

    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("perception");

    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(base.join("snapshots"));
    settings.set_snapshot_suffix(filename.to_string());
    settings.bind(|| {
        assert_yaml_snapshot!(result);
    });
}

#[rstest]
fn test_perception(#[files("tests/perception/data/**/*.xyz")] file_path: PathBuf) {
    run_perception_test(&file_path);
}
