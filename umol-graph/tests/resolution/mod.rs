//! Resolution conformance suite.
//!
//! Runs each `.edn` test input through multiple resolver configurations,
//! producing insta EDN snapshots.

use std::fs;
use std::path::{Component, Path, PathBuf};

use insta::{assert_snapshot, Settings};
use rstest::*;
use umol_edn::{FormatConfig, FromEdn, ToEdn};
use umol_graph::ast::config::{
    ImplicitHydrogenMode, MoleculeAstConfig, MoleculeAstConfigOverrides,
};
use umol_graph::dsl::molecule::MoleculeAstWrapper;
use umol_graph::solver::propagate::{AromaticityConfig, Solution, Solver, ValenceStrategy};
use umol_graph::solver::valence::ValenceTable;

#[derive(FromEdn)]
struct TestInput {
    input: MoleculeAstWrapper,
    #[edn(default)]
    config_overrides: MoleculeAstConfigOverrides,
}

#[derive(ToEdn)]
struct TestResults {
    atom_typing: ResolveResult,
    counts: ResolveResult,
}

#[derive(ToEdn)]
struct ResolveResult {
    success: bool,
    output: Option<MoleculeAstWrapper>,
    error: Option<String>,
}

fn resolve_test(
    input: &MoleculeAstWrapper,
    solver: &Solver,
    config: &MoleculeAstConfig,
) -> ResolveResult {
    let mut ast = input.ast().clone();
    ast.coerce(config);
    match solver.resolve(&mut ast) {
        Ok(Solution::Determined(()) | Solution::Underdetermined(())) => {
            ast.release(config);
            ResolveResult {
                success: true,
                output: Some(MoleculeAstWrapper::new(ast, input.metadata().clone())),
                error: None,
            }
        }
        Ok(Solution::Contradictory) => ResolveResult {
            success: false,
            output: None,
            error: Some("contradictory".to_string()),
        },
        Err(e) => ResolveResult {
            success: false,
            output: None,
            error: Some(e.to_string()),
        },
    }
}

fn atom_typing_solver() -> Solver {
    Solver::default()
}

fn counts_solver(config: &MoleculeAstConfig) -> Solver {
    Solver {
        valence: ValenceStrategy::Counts {
            table: ValenceTable::default_table().clone(),
            allow_implicit_hydrogens: !matches!(
                config.atom.implicit_h_mode,
                ImplicitHydrogenMode::Zero
            ),
        },
        aromaticity: AromaticityConfig::daylight(),
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
    let config = MoleculeAstConfig::open().with_overrides(test_input.config_overrides);

    let atom_typing = resolve_test(
        &test_input.input,
        &atom_typing_solver(),
        &config,
    );
    let counts = resolve_test(
        &test_input.input,
        &counts_solver(&config),
        &config,
    );

    let results = TestResults {
        atom_typing,
        counts,
    };

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
