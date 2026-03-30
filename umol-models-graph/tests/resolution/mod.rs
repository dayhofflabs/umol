//! Resolution conformance suite.
//!
//! Runs each `.toml` test input through multiple resolver configurations,
//! producing insta YAML snapshots with per-atom AtomTypeSpec notation.

use std::fs;
use std::path::{Component, Path, PathBuf};

use insta::{assert_yaml_snapshot, Settings};
use rstest::rstest;
use serde::{Deserialize, Serialize};
use umol_data::SpinMultiplicity;
use umol_models_graph::atom::ImplicitHydrogens;
use umol_models_graph::graph_ir::config_data::ValenceTable;
use umol_models_graph::graph_ir::rings::{RingRelation, RingSet};
use umol_models_graph::graph_ir::{
    resolve_molecule_with, AromaticConstraint, AtomTypeQuery, HydrogenConstraint, Molecule,
    ResolutionError, ResolveConfig, TopologyNodeRef, TopologyProjection, ValenceStrategy,
};
use umol_models_graph::table_ir::{
    Atom as TableAtom, Bond as TableBond, BondDonation, BondOrder, Molecule as TableMolecule,
};

#[derive(Deserialize)]
#[serde(rename_all = "lowercase")]
enum ImplicitHydrogenMode {
    Zero,
    Normal,
    Provided,
}

fn default_implicit_h_mode() -> ImplicitHydrogenMode {
    ImplicitHydrogenMode::Provided
}

#[derive(Deserialize)]
#[serde(rename_all = "lowercase")]
enum ChargeMode {
    Zero,
    Provided,
}

fn default_charge_mode() -> ChargeMode {
    ChargeMode::Zero
}

#[derive(Deserialize)]
#[serde(rename_all = "lowercase")]
enum AromaticMode {
    None,
    Any,
    Provided,
}

fn default_aromatic_mode() -> AromaticMode {
    AromaticMode::Provided
}

#[derive(Deserialize)]
struct TestInput {
    atoms: String,
    bonds: Option<String>,
    dative: Option<String>,
    #[serde(default = "default_implicit_h_mode")]
    implicit_h: ImplicitHydrogenMode,
    #[serde(default = "default_charge_mode")]
    charge: ChargeMode,
    #[serde(default = "default_aromatic_mode")]
    aromatic: AromaticMode,
}

#[derive(Serialize)]
struct AromaticSystemSummary {
    atoms: Vec<usize>,
    electrons: u8,
    charge: i8,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    ring_graph: Vec<String>,
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
    #[serde(skip_serializing_if = "Vec::is_empty")]
    aromatic_systems: Vec<AromaticSystemSummary>,
}

#[derive(Serialize)]
struct ErrorSummary {
    error_type: String,
    message: String,
}

fn parse_atom_token(
    token: &str,
    implicit_h: &ImplicitHydrogenMode,
    charge_mode: &ChargeMode,
    aromatic_mode: &AromaticMode,
) -> TableAtom {
    let query: AtomTypeQuery = token
        .parse()
        .unwrap_or_else(|e| panic!("invalid atom query '{}': {}", token, e));

    let mut atom = TableAtom::from_element(query.element);
    atom.charge = match query.charge {
        Some(charge) => Some(charge),
        None => match charge_mode {
            ChargeMode::Zero => Some(0),
            ChargeMode::Provided => None,
        },
    };
    atom.implicit_hydrogens = match query.implicit_hydrogens {
        Some(HydrogenConstraint::Hydrogens(h)) => Some(ImplicitHydrogens::Hydrogens(h)),
        Some(HydrogenConstraint::Normal) => Some(ImplicitHydrogens::Normal),
        Some(HydrogenConstraint::Any) => None,
        None => match implicit_h {
            ImplicitHydrogenMode::Zero => Some(ImplicitHydrogens::Hydrogens(0)),
            ImplicitHydrogenMode::Normal => Some(ImplicitHydrogens::Normal),
            ImplicitHydrogenMode::Provided => None,
        },
    };
    atom.lone_pairs = query.lone_pairs;
    atom.unpaired_electrons = query.unpaired_electrons;
    atom.multiplicity = query.multiplicity;

    atom.aromatic = match query.aromatic_valence {
        Some(AromaticConstraint::Any | AromaticConstraint::Valence(_)) => Some(true),
        Some(AromaticConstraint::None) => Some(false),
        None => match aromatic_mode {
            AromaticMode::None => Some(false),
            AromaticMode::Any => Some(true),
            AromaticMode::Provided => None,
        },
    };

    atom
}

fn parse_atom_tokens(input: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let mut i = 0usize;

    while i < input.len() {
        // Skip leading whitespace.
        let c = input[i..].chars().next().expect("index in bounds");
        if c.is_ascii_whitespace() {
            i += 1;
            continue;
        }

        // Parse one ?{...} atom token.
        if !input[i..].starts_with("?{") {
            panic!("invalid atoms token stream near: '{}'", &input[i..]);
        }
        let start = i;
        i += 2;
        while i < input.len() {
            let c = input[i..].chars().next().expect("index in bounds");
            i += c.len_utf8();
            if c == '}' {
                break;
            }
        }
        if !input[start..i].ends_with('}') {
            panic!("unterminated atom token in atoms string");
        }
        tokens.push(&input[start..i]);

        continue;
    }

    tokens
}

fn parse_bond_token(token: &str) -> TableBond {
    for (sep, order) in [
        ('-', BondOrder::Single),
        ('=', BondOrder::Double),
        ('#', BondOrder::Triple),
        (':', BondOrder::Aromatic),
    ] {
        if let Some((l, r)) = token.split_once(sep) {
            let (idx_str, charge, mult) = parse_bond_suffix(r, token);
            let (a, b) = split_bond_indices(l, idx_str, token);
            let mut bond = TableBond::new(a, b, order);
            bond.charge = charge;
            bond.multiplicity = mult;
            return bond;
        }
    }
    panic!("invalid bond token: '{}'", token);
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

fn parse_bond_suffix<'a>(
    right: &'a str,
    token: &str,
) -> (&'a str, Option<i8>, Option<SpinMultiplicity>) {
    let (rest, mult) = if let Some((before_x, m_str)) = right.rsplit_once('x') {
        let m: u8 = m_str
            .parse()
            .unwrap_or_else(|_| panic!("invalid multiplicity in bond token '{}'", token));
        let mult = SpinMultiplicity::from_multiplicity(m).unwrap_or_else(|| {
            panic!("invalid multiplicity value {} in bond token '{}'", m, token)
        });
        (before_x, Some(mult))
    } else {
        (right, None)
    };

    if let Some(pos) = rest.find('+').or_else(|| {
        rest.char_indices()
            .skip(1)
            .find(|(_, c)| *c == '-')
            .map(|(i, _)| i)
    }) {
        let idx_str = &rest[..pos];
        let charge: i8 = rest[pos..]
            .parse()
            .unwrap_or_else(|_| panic!("invalid charge in bond token '{}'", token));
        (idx_str, Some(charge), mult)
    } else {
        (rest, None, mult)
    }
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

    for token in parse_atom_tokens(&input.atoms) {
        mol.atoms.push(parse_atom_token(
            token,
            &input.implicit_h,
            &input.charge,
            &input.aromatic,
        ));
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

fn bond_order_symbol(order: u8, charge: i8, multiplicity: SpinMultiplicity) -> String {
    let base = match order {
        0 => ".".to_string(),
        1 => "-".to_string(),
        2 => "=".to_string(),
        3 => "#".to_string(),
        4 => "$".to_string(),
        other => format!("~{}", other),
    };
    let charge_str = match charge {
        0 => String::new(),
        c if c > 0 => format!("+{}", c),
        c => format!("{}", c),
    };
    let mult_str = if multiplicity == SpinMultiplicity::Singlet {
        String::new()
    } else {
        format!("x{}", multiplicity.multiplicity())
    };
    format!("{}{}{}", base, charge_str, mult_str)
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
                Some(mol.atom(ai).unwrap().to_string())
            } else {
                None
            }
        })
        .collect();

    let mut bonds: Vec<(usize, usize, u8, i8, SpinMultiplicity)> = mol
        .bond_indices()
        .map(|idx| {
            let (a, b) = mol.bond_atom_indices(idx).unwrap();
            let ca = canon_pos[a.index()];
            let cb = canon_pos[b.index()];
            let bond = mol.bond(idx).unwrap();
            (
                ca.min(cb),
                ca.max(cb),
                bond.order(),
                bond.charge(),
                bond.multiplicity(),
            )
        })
        .collect();
    bonds.sort();
    let bonds: Vec<String> = bonds
        .iter()
        .map(|&(_, _, o, c, m)| bond_order_symbol(o, c, m))
        .collect();

    let mut aromatic_systems: Vec<AromaticSystemSummary> = mol
        .aromatic_systems()
        .map(|system| {
            let system_atoms: Vec<_> = system.atoms().collect();
            let mut atoms: Vec<usize> = system_atoms
                .iter()
                .map(|ai| canon_pos[ai.index()])
                .collect();
            atoms.sort_unstable();
            let display_ring_set = RingSet::induced_from_molecule_atoms(mol, &system_atoms);
            let mut ring_graph: Vec<String> = display_ring_set
                .ring_graph()
                .edges()
                .iter()
                .map(|edge| {
                    format!(
                        "{}-{}:{}",
                        edge.source.index(),
                        edge.target.index(),
                        ring_relation_code(edge.relation)
                    )
                })
                .collect();
            ring_graph.sort_unstable();

            AromaticSystemSummary {
                atoms,
                electrons: system.electron_count(),
                charge: system.charge(),
                ring_graph,
            }
        })
        .collect();
    aromatic_systems.sort_by(|a, b| {
        a.atoms
            .cmp(&b.atoms)
            .then(a.electrons.cmp(&b.electrons))
            .then(a.charge.cmp(&b.charge))
    });

    ResolveSummary {
        atom_count: mol.atom_count(),
        bond_count: mol.bond_count(),
        topology,
        atoms,
        bonds,
        aromatic_systems,
    }
}

fn ring_relation_code(relation: RingRelation) -> &'static str {
    match relation {
        RingRelation::Fused => "f",
        RingRelation::Spiro => "s",
        RingRelation::Bridged => "b",
        RingRelation::MultiSpiro => "m",
        RingRelation::Noncontiguous => "n",
        RingRelation::Disjoint | RingRelation::Identical => "d",
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
    ResolveConfig::default()
}

fn counts_config(mode: &ImplicitHydrogenMode) -> ResolveConfig {
    let mut config = ResolveConfig::default();
    config.valence.strategy = ValenceStrategy::Counts {
        table: ValenceTable::default_table().clone(),
        allow_implicit_hydrogens: !matches!(mode, ImplicitHydrogenMode::Zero),
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

fn resolve_file(path: &Path) -> FileResolveResults {
    let content = fs::read_to_string(path).expect("failed to read test file");
    let input: TestInput = toml::from_str(&content).expect("failed to parse TOML input");
    let table_mol = build_table_molecule(&input);
    let category = extract_category(path);

    let atom_typing = resolve_with_config(&table_mol, &atom_typing_config());
    let counts = resolve_with_config(&table_mol, &counts_config(&input.implicit_h));

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
// Keep recursive glob for fast auto-discovery of new files at compile time (refresh marker v9).
fn test_conformance(#[files("tests/resolution/data/**/*.toml")] file_path: PathBuf) {
    run_conformance_test(&file_path);
}
