//! Resolution conformance suite.
//!
//! Runs each `.toml` test input through multiple resolver configurations,
//! producing insta YAML snapshots with per-atom AtomTypeSpec notation.

use std::path::{Component, Path, PathBuf};

use insta::{assert_yaml_snapshot, Settings};
use rstest::rstest;
use serde::{Deserialize, Serialize};
use umol_models_graph::graph_ir::{
    resolve_molecule_with, AtomTypeQuery, Molecule, ResolutionError, ResolveConfig,
    TopologyNodeRef, TopologyProjection, ValenceStrategy,
};
use umol_models_graph::table_ir::{
    Atom as TableAtom, Bond as TableBond, BondDonation, BondOrder, Molecule as TableMolecule,
    UnpairedElectrons,
};

#[derive(Deserialize)]
struct TestInput {
    atoms: String,
    bonds: Option<String>,
    dative: Option<String>,
    #[serde(default = "default_true")]
    implicit_h: bool,
}

fn default_true() -> bool {
    true
}

#[derive(Serialize)]
struct FileResolveResults {
    category: String,
    atom_typing: ResolveResult,
    counts: ResolveResult,
}

#[derive(Serialize)]
struct ResolveResult {
    success: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    summary: Option<ResolveSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<ErrorSummary>,
}

#[derive(Serialize)]
struct ResolveSummary {
    atom_count: usize,
    bond_count: usize,
    topology: String,
    atoms: Vec<String>,
    bonds: Vec<String>,
}

#[derive(Serialize)]
struct ErrorSummary {
    error_type: String,
    message: String,
}

fn parse_atom_token(token: &str, implicit_h: bool) -> TableAtom {
    let query: AtomTypeQuery = token
        .parse()
        .unwrap_or_else(|e| panic!("invalid atom query '{}': {}", token, e));

    let mut atom = TableAtom::from_element(query.element);
    atom.charge = query.charge;
    atom.hydrogens = if implicit_h {
        Some(query.hydrogens.unwrap_or_else(|| {
            panic!(
                "missing hydrogen token in '{}' while implicit_h=true; specify H0 explicitly",
                token
            )
        }))
    } else {
        Some(query.hydrogens.unwrap_or(0))
    };
    atom.lone_pairs = query.lone_pairs;

    if let Some(unpaired) = query.unpaired_electrons {
        atom.unpaired_electrons = Some(UnpairedElectrons::new(unpaired, query.multiplicity));
    }

    atom
}

fn split_bond_indices(left: &str, right: &str, token: &str) -> (u32, u32) {
    let a = left
        .parse()
        .unwrap_or_else(|_| panic!("invalid bond index in '{}'", token));
    let b = right
        .parse()
        .unwrap_or_else(|_| panic!("invalid bond index in '{}'", token));
    (a, b)
}

fn parse_bond_token(token: &str) -> TableBond {
    for (sep, order) in [
        ('-', BondOrder::Single),
        ('=', BondOrder::Double),
        ('#', BondOrder::Triple),
        (':', BondOrder::Aromatic),
    ] {
        if let Some((l, r)) = token.split_once(sep) {
            let (a, b) = split_bond_indices(l, r, token);
            return TableBond::new(a, b, order);
        }
    }
    panic!("invalid bond token: '{}'", token);
}

fn parse_dative_token(token: &str) -> TableBond {
    let (l, r) = token
        .split_once("->")
        .unwrap_or_else(|| panic!("dative bond must use '->' notation: '{}'", token));
    let (a, b) = split_bond_indices(l, r, token);
    TableBond::new_dative(a, b, BondOrder::Single, BondDonation::Donating)
}

fn build_table_molecule(input: &TestInput) -> TableMolecule {
    let mut mol = TableMolecule::empty();

    for token in input.atoms.split_whitespace() {
        mol.atoms.push(parse_atom_token(token, input.implicit_h));
    }

    if let Some(ref bonds) = input.bonds {
        let trimmed = bonds.trim();
        if !trimmed.is_empty() {
            for token in trimmed.split_whitespace() {
                mol.bonds.push(parse_bond_token(token));
            }
        }
    }

    if let Some(ref dative) = input.dative {
        let trimmed = dative.trim();
        if !trimmed.is_empty() {
            for token in trimmed.split_whitespace() {
                mol.bonds.push(parse_dative_token(token));
            }
        }
    }

    mol
}

fn bond_order_symbol(order: u8) -> String {
    match order {
        0 => ".".to_string(),
        1 => "-".to_string(),
        2 => "=".to_string(),
        3 => "#".to_string(),
        4 => "$".to_string(),
        other => format!("~{}", other),
    }
}

fn summarize(mol: &Molecule) -> ResolveSummary {
    let tg = mol.topology_graph(TopologyProjection::ordinary());
    let (topology, order) = tg.to_graph6_canonical().unwrap();

    let mut canon_pos = vec![0usize; mol.atom_count()];
    for (pos, nidx) in order.iter().enumerate() {
        if let Some(TopologyNodeRef::Atom(ai)) = tg.node_ref(*nidx) {
            canon_pos[ai.index()] = pos;
        }
    }

    let atoms: Vec<String> = order
        .iter()
        .filter_map(|nidx| {
            if let Some(TopologyNodeRef::Atom(ai)) = tg.node_ref(*nidx) {
                Some(mol.atom(ai).unwrap().to_spec().to_string())
            } else {
                None
            }
        })
        .collect();

    let mut bonds: Vec<(usize, usize, u8)> = mol
        .bond_indices()
        .map(|idx| {
            let (a, b) = mol.bond_atom_indices(idx).unwrap();
            let ca = canon_pos[a.index()];
            let cb = canon_pos[b.index()];
            let order = mol.bond(idx).unwrap().order();
            (ca.min(cb), ca.max(cb), order)
        })
        .collect();
    bonds.sort();
    let bonds: Vec<String> = bonds
        .iter()
        .map(|&(_, _, o)| bond_order_symbol(o))
        .collect();

    ResolveSummary {
        atom_count: mol.atom_count(),
        bond_count: mol.bond_count(),
        topology,
        atoms,
        bonds,
    }
}

fn error_summary(e: &ResolutionError) -> ErrorSummary {
    let debug = format!("{:?}", e);
    let error_type = debug
        .split(|c| c == '{' || c == '(')
        .next()
        .unwrap_or(&debug)
        .trim()
        .to_string();
    ErrorSummary {
        error_type,
        message: e.to_string(),
    }
}

fn resolve_with_config(table_mol: &TableMolecule, config: &ResolveConfig) -> ResolveResult {
    match resolve_molecule_with(table_mol, config) {
        Ok(mol) => ResolveResult {
            success: true,
            summary: Some(summarize(&mol)),
            error: None,
        },
        Err(e) => ResolveResult {
            success: false,
            summary: None,
            error: Some(error_summary(&e)),
        },
    }
}

fn atom_typing_config() -> ResolveConfig {
    let mut config = ResolveConfig::default();
    config.valence.strategy = ValenceStrategy::AtomTyping;
    config
}

fn counts_config(enable_implicit_hydrogens: bool) -> ResolveConfig {
    let mut config = ResolveConfig::default();
    config.valence.strategy = ValenceStrategy::Counts;
    config.valence.enable_implicit_hydrogens = enable_implicit_hydrogens;
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

fn resolve_file(path: &Path) -> FileResolveResults {
    let content = std::fs::read_to_string(path).expect("failed to read test file");
    let input: TestInput = toml::from_str(&content).expect("failed to parse TOML input");
    let table_mol = build_table_molecule(&input);
    let category = extract_category(path);

    let atom_typing = resolve_with_config(&table_mol, &atom_typing_config());
    let counts = resolve_with_config(&table_mol, &counts_config(input.implicit_h));

    FileResolveResults {
        category,
        atom_typing,
        counts,
    }
}

fn run_conformance_test(file_path: &PathBuf) {
    let source_dir = extract_category(file_path);
    let filename = file_path.file_stem().unwrap().to_str().unwrap();

    let results = resolve_file(file_path);

    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("resolution");

    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(base.join("snapshots"));
    settings.set_snapshot_suffix(format!("{}_{}", source_dir, filename));
    settings.bind(|| {
        assert_yaml_snapshot!(results);
    });
}

#[rstest]
// Keep recursive glob so category folders and recently added files are auto-discovered at compile time (refresh marker v2).
fn test_conformance(#[files("tests/resolution/data/**/*.toml")] file_path: PathBuf) {
    run_conformance_test(&file_path);
}
