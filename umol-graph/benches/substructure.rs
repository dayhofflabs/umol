//! End-to-end substructure matching over the OpenSMILES corpus
//! (`umol-io/tests/smiles_parsing/data/opensmiles`, ~9k molecules).
//!
//! Sweeps both `MoleculeAst::substructure_matches` strategies against all six
//! subgraph-isomorphism algorithms, over three patterns of increasing cost. The
//! corpus, patterns, and matching semantics (element-only atoms, any-bonds,
//! all-embeddings enumeration) mirror `scripts/rdkit_substructure_baseline.py`
//! so the timings are directly comparable to the actual-RDKit baseline.

use std::fs::read_to_string;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_ast::ast::SubstructureMatchAlgorithm::{GraphAndOverlays, Incidence};
use umol_ast::ast::{
    AtomAst, AtomId, BondAst, MoleculeAst, MoleculeEntries, SubstructureMatchAlgorithm,
    SubstructureMatchConfig, ValueAst,
};
use umol_chem::element::Element;
use umol_graph::ingest::ingest_smiles;
use umol_graph_core::SubgraphIsomorphismAlgorithm::{
    ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
};
use umol_graph_core::{
    RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm, ARCMATCH_DEFAULT_PATH_LENGTH,
};
use walkdir::WalkDir;

const SUBISO: [SubgraphIsomorphismAlgorithm; 6] = [
    Vf2,
    Ullmann,
    Ri,
    ArcMatch {
        path_length: ARCMATCH_DEFAULT_PATH_LENGTH,
    },
    Vf2Rdkit,
    RayKirsch,
];
const STRATEGIES: [SubstructureMatchAlgorithm; 2] = [GraphAndOverlays, Incidence];

fn subiso_name(algorithm: SubgraphIsomorphismAlgorithm) -> &'static str {
    match algorithm {
        Vf2 => "vf2",
        Ullmann => "ullmann",
        Ri => "ri",
        ArcMatch { .. } => "arcmatch",
        Vf2Rdkit => "vf2rdkit",
        RayKirsch => "raykirsch",
    }
}

fn strategy_name(strategy: SubstructureMatchAlgorithm) -> &'static str {
    match strategy {
        GraphAndOverlays => "graph_overlays",
        Incidence => "incidence",
    }
}

fn load_corpus() -> Vec<MoleculeAst> {
    let dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../umol-io/tests/smiles_parsing/data/opensmiles"
    );
    let mut molecules = Vec::new();
    for entry in WalkDir::new(dir)
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| e.path().extension().is_some_and(|x| x == "smiles"))
    {
        let Ok(content) = read_to_string(entry.path()) else {
            continue;
        };
        if let Some(smiles) = content.lines().nth(1) {
            if !smiles.is_empty() {
                if let Ok(molecule) = ingest_smiles(smiles) {
                    molecules.push(molecule);
                }
            }
        }
    }
    molecules
}

fn carbon() -> AtomAst {
    AtomAst::from_element(Element::C)
}

fn any_bond() -> BondAst {
    BondAst::new(ValueAst::Undetermined)
}

fn pattern(atoms: Vec<AtomAst>, bonds: Vec<(u32, u32, BondAst)>) -> MoleculeAst {
    let bond_list = bonds
        .into_iter()
        .map(|(s, t, b)| (AtomId(s), AtomId(t), b))
        .collect();
    MoleculeAst::from_entries(MoleculeEntries {
        atoms,
        bonds: bond_list,
        ..Default::default()
    })
}

/// `C(C)C(C)N` — 5 atoms, all any-bonds.
fn pattern_branched() -> MoleculeAst {
    pattern(
        vec![
            carbon(),
            carbon(),
            carbon(),
            carbon(),
            AtomAst::from_element(Element::N),
        ],
        vec![
            (0, 1, any_bond()),
            (0, 2, any_bond()),
            (2, 3, any_bond()),
            (2, 4, any_bond()),
        ],
    )
}

/// `c1ccccc1O` — 6-ring + hydroxyl, all any-bonds.
fn pattern_phenol() -> MoleculeAst {
    pattern(
        vec![
            carbon(),
            carbon(),
            carbon(),
            carbon(),
            carbon(),
            carbon(),
            AtomAst::from_element(Element::O),
        ],
        vec![
            (0, 1, any_bond()),
            (1, 2, any_bond()),
            (2, 3, any_bond()),
            (3, 4, any_bond()),
            (4, 5, any_bond()),
            (5, 0, any_bond()),
            (5, 6, any_bond()),
        ],
    )
}

/// Fused 5-6 bicyclic carbon skeleton (ring A 0-1-5-6-7-8, ring B 1-2-3-4-5,
/// fused edge 1-5), all any-bonds.
fn pattern_bicyclic() -> MoleculeAst {
    pattern(
        (0..9).map(|_| carbon()).collect(),
        vec![
            (0, 1, any_bond()),
            (1, 2, any_bond()),
            (2, 3, any_bond()),
            (3, 4, any_bond()),
            (4, 5, any_bond()),
            (1, 5, any_bond()),
            (5, 6, any_bond()),
            (6, 7, any_bond()),
            (7, 8, any_bond()),
            (8, 0, any_bond()),
        ],
    )
}

fn substructure_benchmark(c: &mut Criterion) {
    let corpus = load_corpus();
    let patterns: [(&str, MoleculeAst); 3] = [
        ("branched", pattern_branched()),
        ("phenol", pattern_phenol()),
        ("bicyclic", pattern_bicyclic()),
    ];

    let mut group = c.benchmark_group("substructure");
    for (pattern_name, pat) in &patterns {
        for strategy in STRATEGIES {
            for algorithm in SUBISO {
                let id = format!(
                    "{pattern_name}/{}/{}",
                    strategy_name(strategy),
                    subiso_name(algorithm)
                );
                group.bench_function(id, |b| {
                    b.iter(|| {
                        for target in &corpus {
                            black_box(pat.substructure_matches(
                                target,
                                SubstructureMatchConfig {
                                    match_algorithm: strategy,
                                    subgraph_isomorphism_algorithm: algorithm,
                                    relevant_cycle_algorithm:
                                        RelevantCycleEnumerationAlgorithm::Vismara,
                                },
                            ));
                        }
                    });
                });
            }
        }
    }
    group.finish();
}

criterion_group!(benches, substructure_benchmark);
criterion_main!(benches);
