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
use umol_graph::dsl::molecule::MoleculeDsl;
use umol_graph::ops::aromaticity::AromaticityTheory;
use umol_graph::ops::chemistry::Chemistry;
use umol_graph::ops::propagate::ValenceTheory;
use umol_graph::ops::resolve::Resolver;
use umol_graph::ops::solution::Solution;
use umol_graph::ops::valence::ValenceTable;

#[derive(FromEdn)]
struct TestInput {
    input: MoleculeDsl,
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
    output: Option<MoleculeDsl>,
    error: Option<String>,
}

fn resolve_test(
    input: &MoleculeDsl,
    chemistry: &Chemistry,
    config: &MoleculeAstConfig,
) -> ResolveResult {
    let (mut ast, metadata) = input.lower_parts().unwrap();
    ast.coerce(config);
    match Resolver::new(chemistry).resolve(ast) {
        Ok(Solution::Determined(mut ast)) => {
            ast.release(&MoleculeAstConfig::zeroed());
            ResolveResult {
                success: true,
                output: Some(MoleculeDsl::new(ast, metadata)),
                error: None,
            }
        }
        Ok(Solution::Underdetermined(_)) => ResolveResult {
            success: false,
            output: None,
            error: Some("resolution underdetermined".to_string()),
        },
        Ok(Solution::Contradictory) => ResolveResult {
            success: false,
            output: None,
            error: Some("resolution contradictory".to_string()),
        },
        Err(e) => ResolveResult {
            success: false,
            output: None,
            error: Some(e.to_string()),
        },
    }
}

fn atom_typing_chemistry() -> Chemistry {
    Chemistry::default()
}

fn counts_chemistry(config: &MoleculeAstConfig) -> Chemistry {
    Chemistry {
        valence: ValenceTheory::Counts {
            table: ValenceTable::default_table().clone(),
            allow_implicit_hydrogens: !matches!(
                config.atom.implicit_h_mode,
                ImplicitHydrogenMode::Zero
            ),
        },
        aromaticity: AromaticityTheory::daylight(),
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
        &atom_typing_chemistry(),
        &config,
    );
    let counts = resolve_test(
        &test_input.input,
        &counts_chemistry(&config),
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
