//! Resolution conformance suite.
//!
//! Runs each `.edn` test input through multiple resolver configurations,
//! producing insta EDN snapshots.

use std::fs;
use std::path::{Component, Path, PathBuf};

use insta::{assert_snapshot, Settings};
use rstest::*;
use serde::{Deserialize, Serialize};
use umol_graph::dsl::ast::{FromAst, ToAst};
use umol_graph::dsl::config::{
    ImplicitHydrogenMode, MoleculeDslConfig, MoleculeDslConfigOverrides,
};
use umol_edn::{from_str as edn_from_str, to_string_pretty as edn_to_string_pretty};
use umol_graph::dsl::molecule::MoleculeAst;
use umol_graph::graph_ir::config_data::ValenceTable;
use umol_graph::graph_ir::molecule_builder::{MoleculeBuilder, ResolutionContext};
use umol_graph::graph_ir::{resolve_molecule_with, ResolveConfig, ValenceStrategy};

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct TestInput {
    input: MoleculeAst,
    #[serde(default)]
    config_overrides: MoleculeDslConfigOverrides,
    #[serde(default)]
    context: ResolutionContext,
}

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct TestResults {
    atom_typing: ResolveResult,
    counts: ResolveResult,
}

#[derive(Serialize)]
struct ResolveResult {
    success: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    output: Option<MoleculeAst>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<String>,
}

fn resolve_with_config(
    input: &MoleculeAst,
    dsl_config: &MoleculeDslConfig,
    context: &ResolutionContext,
    resolve_config: &ResolveConfig,
) -> ResolveResult {
    let builder = MoleculeBuilder::from_ast(input.clone(), dsl_config);
    match builder {
        Err(e) => ResolveResult {
            success: false,
            output: None,
            error: Some(e.to_string()),
        },
        Ok(mut builder) => {
            builder.set_resolution_context(context.clone());
            let result = resolve_molecule_with(&mut builder, resolve_config)
                .and_then(|()| builder.build(resolve_config));
            match result {
                Ok(mol) => {
                    let mut ast = mol.to_ast(&MoleculeDslConfig::zeroed());
                    ast.aromatic_systems.sort_by(|a, b| a.atoms.cmp(&b.atoms));
                    ResolveResult {
                        success: true,
                        output: Some(ast),
                        error: None,
                    }
                }
                Err(e) => ResolveResult {
                    success: false,
                    output: None,
                    error: Some(e.to_string()),
                },
            }
        }
    }
}

fn counts_config(dsl_config: &MoleculeDslConfig) -> ResolveConfig {
    let mut config = ResolveConfig::default();
    config.valence.strategy = ValenceStrategy::Counts {
        table: ValenceTable::default_table().clone(),
        allow_implicit_hydrogens: !matches!(
            dsl_config.atom.implicit_h_mode,
            ImplicitHydrogenMode::Zero
        ),
    };
    config
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
    let test_input: TestInput = edn_from_str(&content).expect("failed to parse EDN input");
    let config = MoleculeDslConfig::open().with_overrides(test_input.config_overrides);

    let atom_typing = resolve_with_config(
        &test_input.input,
        &config,
        &test_input.context,
        &ResolveConfig::default(),
    );
    let counts = resolve_with_config(
        &test_input.input,
        &config,
        &test_input.context,
        &counts_config(&config),
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
        assert_snapshot!(edn_to_string_pretty(&results).unwrap());
    });
}

#[rstest]
fn test_conformance(#[files("tests/resolution/data/**/*.edn")] file_path: PathBuf) {
    run_conformance_test(&file_path);
}
