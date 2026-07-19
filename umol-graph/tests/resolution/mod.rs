//! Resolution conformance suite.
//!
//! Runs each `.edn` test input through the atom-typing and counts resolver
//! configurations, producing insta EDN snapshots.

use std::borrow::Cow;
use std::fs;
use std::path::{Component, Path, PathBuf};

use insta::{assert_snapshot, Settings};
use rstest::*;
use umol_ast::ast::{FromAst, IntoAst, MoleculeAst};
use umol_ast::dsl::{MoleculeDefaults, MoleculeDsl, MoleculeOverrides};
use umol_edn::{FormatConfig, FromEdn, ToEdn};
use umol_graph::ops::model::{AromaticityModel, ChemistryModel, StereoModel, ValenceModel};
use umol_graph::ops::resolve::Resolver;
use umol_graph::ops::valence::ValenceTable;
use umol_utils::solution::Solution;

#[derive(FromEdn)]
struct TestInput {
    input: MoleculeDsl,
    #[edn(default)]
    config_overrides: MoleculeOverrides,
}

#[derive(ToEdn)]
struct TestResults {
    atom_typing: ResolveResult,
    counts: ResolveResult,
}

#[derive(ToEdn)]
struct ResolveResult {
    success: bool,
    output: Option<MoleculeDsl>,
    error: Option<String>,
}

fn raise(input: &MoleculeDsl, defaults: &MoleculeDefaults) -> MoleculeAst {
    input.clone().into_ast(defaults)
}

fn lower(ast: &MoleculeAst) -> MoleculeDsl {
    let cfg = MoleculeDefaults::zeroed();
    MoleculeDsl::from_ast(ast, &cfg)
}

fn resolve_test(
    input: &MoleculeDsl,
    chemistry: &ChemistryModel,
    defaults: &MoleculeDefaults,
) -> ResolveResult {
    let mut ast = raise(input, defaults);
    match Resolver::new(chemistry).resolve(&mut ast) {
        Ok(Solution::Determined(())) => ResolveResult {
            success: true,
            output: Some(lower(&ast)),
            error: None,
        },
        Ok(Solution::Underdetermined(())) => ResolveResult {
            success: false,
            output: None,
            error: Some("resolution underdetermined".to_string()),
        },
        Ok(Solution::Contradictory(c)) => ResolveResult {
            success: false,
            output: None,
            error: Some(format!("resolution contradictory: {}", c)),
        },
        Err(e) => ResolveResult {
            success: false,
            output: None,
            error: Some(e.to_string()),
        },
    }
}

fn atom_typing_chemistry() -> ChemistryModel {
    ChemistryModel::default()
}

fn counts_chemistry() -> ChemistryModel {
    ChemistryModel {
        valence: ValenceModel::Counts {
            table: Cow::Borrowed(ValenceTable::default_table()),
        },
        aromaticity: AromaticityModel::daylight(),
        stereo: StereoModel::default(),
    }
}

fn extract_category(path: &Path) -> String {
    let components: Vec<_> = path.components().collect();
    for (i, comp) in components.iter().enumerate() {
        if let Component::Normal(name) = comp {
            if name.to_str() == Some("data") && i + 1 < components.len() {
                if let Component::Normal(dir_name) = &components[i + 1] {
                    return dir_name.to_str().unwrap_or("unknown").to_string();
                }
            }
        }
    }
    "unknown".to_string()
}

fn run_conformance_test(file_path: &Path) {
    let content = fs::read_to_string(file_path).expect("failed to read test file");
    let test_input = TestInput::from_edn_str(&content).expect("failed to parse EDN input");
    let defaults = MoleculeDefaults::default().with_overrides(test_input.config_overrides);

    let atom_typing = resolve_test(&test_input.input, &atom_typing_chemistry(), &defaults);
    let counts = resolve_test(&test_input.input, &counts_chemistry(), &defaults);

    let results = TestResults {
        atom_typing,
        counts,
    };

    assert!(
        results.atom_typing.success,
        "atom-typing resolution did not succeed: {:?}",
        results.atom_typing.error
    );
    assert!(
        results.counts.success,
        "counts resolution did not succeed: {:?}",
        results.counts.error
    );

    let source_dir = extract_category(file_path);
    let filename = file_path.file_stem().unwrap().to_str().unwrap();

    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("resolution");

    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(base.join("snapshots"));
    settings.set_snapshot_suffix(format!("{}_{}", source_dir, filename));
    settings.bind(|| {
        assert_snapshot!(results.to_edn().to_string_with(&FormatConfig::default()));
    });
}

#[rstest]
fn test_conformance(#[files("tests/resolution/data/**/*.edn")] file_path: PathBuf) {
    run_conformance_test(&file_path);
}
