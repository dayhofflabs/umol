use std::fs;
use std::path::{Path, PathBuf};

use insta::{assert_yaml_snapshot, Settings};
use rstest::rstest;
use serde::Serialize;
use umol_geometric::io::xyz::parse_xyz;
use umol_geometric::molecule::Molecule;
use umol_msym::Thresholds;

#[derive(Serialize)]
#[serde(untagged)]
enum PerceptionResult {
    Ok {
        point_group: String,
        order: usize,
        equivalence_sets: Vec<Vec<String>>,
    },
    Err {
        error_code: i32,
        error: String,
    },
}

fn perceive(mol: &Molecule) -> PerceptionResult {
    let sym = match mol.symmetrize(Thresholds::default()) {
        Ok(sym) => sym,
        Err(e) => {
            return PerceptionResult::Err {
                error_code: e.code,
                error: e.to_string(),
            };
        }
    };

    let equivalence_sets: Vec<Vec<String>> = sym
        .equivalence_sets()
        .iter()
        .map(|set| {
            set.iter()
                .map(|&i| format!("{}({})", sym.element(i).symbol(), i))
                .collect()
        })
        .collect();

    PerceptionResult::Ok {
        point_group: sym.point_group().to_string(),
        order: sym.point_group().order(),
        equivalence_sets,
    }
}

fn run_perception_test(file_path: &Path) {
    let content = fs::read_to_string(file_path).unwrap();
    let (mol, _comment) = parse_xyz(&content).unwrap();
    let result = perceive(&mol);

    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("perception")
        .join("data");
    let rel = file_path.strip_prefix(&data_dir).unwrap();
    let suffix = rel
        .with_extension("")
        .to_str()
        .unwrap()
        .replace('/', "_");

    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("perception");

    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(base.join("snapshots"));
    settings.set_snapshot_suffix(suffix);
    settings.bind(|| {
        assert_yaml_snapshot!(result);
    });
}

// Bump to force recompile when data files change (rstest evaluates globs at compile time).
const _REFRESH: u32 = 11;

#[rstest]
fn test_perception(#[files("tests/perception/data/**/*.xyz")] file_path: PathBuf) {
    run_perception_test(&file_path);
}
