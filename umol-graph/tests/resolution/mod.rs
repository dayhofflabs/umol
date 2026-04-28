//! Resolution conformance suite.
//!
//! Runs each `.edn` test input through the atom-typing and counts resolver
//! configurations, producing insta EDN snapshots.

use std::fs;
use std::path::{Component, Path, PathBuf};

use insta::{assert_snapshot, Settings};
use rstest::*;
use umol_ast::ast::{FromAst, IntoAst, MoleculeAst};
use umol_ast::dsl::{
    ImplicitHydrogensDefault, MoleculeDefaults, MoleculeDsl, MoleculeOverrides,
};
use umol_edn::{FormatConfig, FromEdn, ToEdn};
use umol_graph::ops::config::{AromaticityModel, ChemistryModel, ValenceModel};
use umol_graph::ops::resolver::Resolver;
use umol_graph::ops::solution::Solution;
use umol_graph::ops::valence::ValenceTable;

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

fn lower(input: &MoleculeDsl, defaults: &MoleculeDefaults) -> Result<MoleculeAst, String> {
    input
        .clone()
        .into_ast(defaults)
        .map_err(|e| format!("lowering: {}", e))
}

fn raise(ast: &MoleculeAst) -> Result<MoleculeDsl, String> {
    let zeroed = MoleculeDefaults::zeroed();
    MoleculeDsl::from_ast(ast, &zeroed).map_err(|e| format!("raising: {}", e))
}

fn resolve_test(
    input: &MoleculeDsl,
    chemistry: &ChemistryModel,
    defaults: &MoleculeDefaults,
) -> ResolveResult {
    let mut ast = match lower(input, defaults) {
        Ok(ast) => ast,
        Err(error) => {
            return ResolveResult {
                success: false,
                output: None,
                error: Some(error),
            }
        }
    };
    match Resolver::new(chemistry).resolve(&mut ast) {
        Ok(Solution::Determined(())) => match raise(&ast) {
            Ok(out) => ResolveResult {
                success: true,
                output: Some(out),
                error: None,
            },
            Err(error) => ResolveResult {
                success: false,
                output: None,
                error: Some(error),
            },
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

fn counts_chemistry(defaults: &MoleculeDefaults) -> ChemistryModel {
    ChemistryModel {
        valence: ValenceModel::Counts {
            table: ValenceTable::default_table().clone(),
            allow_implicit_hydrogens: !matches!(
                defaults.atom.implicit_hydrogens,
                ImplicitHydrogensDefault::Zero,
            ),
        },
        aromaticity: AromaticityModel::daylight(),
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
    let defaults = MoleculeDefaults::verbatim().with_overrides(test_input.config_overrides);

    let atom_typing = resolve_test(&test_input.input, &atom_typing_chemistry(), &defaults);
    let counts = resolve_test(&test_input.input, &counts_chemistry(&defaults), &defaults);

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
