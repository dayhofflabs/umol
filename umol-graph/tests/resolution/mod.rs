//! Resolution conformance suite.
//!
//! Runs each `.edn` test input through the atom-typing and counts resolver
//! configurations, producing insta EDN snapshots.

use std::borrow::Cow;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Component, Path, PathBuf};

use insta::{assert_snapshot, Settings};
use rstest::*;
use umol_edn::{FormatConfig, FromEdn, ToEdn};
use umol_graph::ops::model::{
    AromaticityModel, ChemistryModel, StereoModel, ValenceModel, ValenceTieBreak,
};
use umol_graph::ops::resolve::Resolver;
use umol_graph::ops::valence::{ResolveReport, ValenceTable};
use umol_graph::ops::validate::ConnectivityModel;
use umol_graph_ir::dsl::{AtomDefaults, AtomDsl, MoleculeDefaults, MoleculeDsl, MoleculeOverrides};
use umol_graph_ir::ir::{FromIr, IntoIr, Molecule};
use umol_utils::solution::Solution;

#[derive(FromEdn)]
struct TestInput {
    input: MoleculeDsl,
    #[edn(default)]
    config_overrides: MoleculeOverrides,
}

#[derive(ToEdn)]
struct TestResults {
    atom_typing_strict: ResolveResult,
    atom_typing_most_saturated: ResolveResult,
    counts_strict: ResolveResult,
    counts_most_saturated: ResolveResult,
}

#[derive(ToEdn)]
struct ResolveResult {
    success: bool,
    output: Option<MoleculeDsl>,
    tie_breaks: Vec<u32>,
    unresolved: Option<BTreeMap<u32, Vec<String>>>,
    error: Option<String>,
}

fn tie_break_ids(report: &ResolveReport) -> Vec<u32> {
    report
        .tie_breaks
        .iter()
        .map(|atom| atom.index() as u32)
        .collect()
}

fn raise(input: &MoleculeDsl, defaults: &MoleculeDefaults) -> Molecule {
    input.clone().into_ir(defaults)
}

fn lower(molecule: &Molecule) -> MoleculeDsl {
    let cfg = MoleculeDefaults::ground();
    MoleculeDsl::from_ir(molecule, &cfg)
}

fn resolve_test(
    input: &MoleculeDsl,
    chemistry: &ChemistryModel,
    defaults: &MoleculeDefaults,
) -> ResolveResult {
    let mut molecule = raise(input, defaults);
    match Resolver::new(chemistry).resolve(&mut molecule) {
        Ok(Solution::Determined(report)) => ResolveResult {
            success: true,
            output: Some(lower(&molecule)),
            tie_breaks: tie_break_ids(&report),
            unresolved: None,
            error: None,
        },
        Ok(Solution::Underdetermined(report)) => ResolveResult {
            success: false,
            output: None,
            tie_breaks: tie_break_ids(&report),
            unresolved: Some(
                report
                    .unresolved
                    .iter()
                    .map(|(atom, forms)| {
                        (
                            atom.index() as u32,
                            forms
                                .iter()
                                .map(|form| {
                                    AtomDsl::from_ir(form, &AtomDefaults::ground()).to_string()
                                })
                                .collect(),
                        )
                    })
                    .collect(),
            ),
            error: None,
        },
        Ok(Solution::Contradictory(c)) => ResolveResult {
            success: false,
            output: None,
            tie_breaks: Vec::new(),
            unresolved: None,
            error: Some(format!("resolution contradictory: {}", c)),
        },
        Err(e) => ResolveResult {
            success: false,
            output: None,
            tie_breaks: Vec::new(),
            unresolved: None,
            error: Some(e.to_string()),
        },
    }
}

fn atom_typing_chemistry(tie_break: ValenceTieBreak) -> ChemistryModel {
    ChemistryModel {
        valence: ValenceModel {
            tie_break,
            ..ValenceModel::default()
        },
        ..ChemistryModel::default()
    }
}

fn counts_chemistry(tie_break: ValenceTieBreak) -> ChemistryModel {
    ChemistryModel {
        valence: ValenceModel {
            tie_break,
            ..ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table()))
        },
        connectivity: ConnectivityModel::default(),
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

    let results = TestResults {
        atom_typing_strict: resolve_test(
            &test_input.input,
            &atom_typing_chemistry(ValenceTieBreak::Strict),
            &defaults,
        ),
        atom_typing_most_saturated: resolve_test(
            &test_input.input,
            &atom_typing_chemistry(ValenceTieBreak::MostSaturated),
            &defaults,
        ),
        counts_strict: resolve_test(
            &test_input.input,
            &counts_chemistry(ValenceTieBreak::Strict),
            &defaults,
        ),
        counts_most_saturated: resolve_test(
            &test_input.input,
            &counts_chemistry(ValenceTieBreak::MostSaturated),
            &defaults,
        ),
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

#[rstest]
#[ignore]
fn generate_f420_input() {
    let smiles = "C[C@H](OP(=O)(O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2nc(=O)nc(=O)c-2cc2ccc(O)cc21)C(=O)N[C@@H](CCC(=O)N[C@@H](CCC(=O)O)C(=O)O)C(=O)O";
    let table = umol_io::smiles::Smiles::parse(smiles)
        .unwrap()
        .into_table_ir();
    let molecule: Molecule = (&table).try_into_ir(&()).unwrap();
    let dsl = MoleculeDsl::from_ir(&molecule, &MoleculeDefaults::default());
    println!(
        "{{:input {}}}",
        dsl.to_edn().to_string_with(&FormatConfig::default())
    );
}
