//! End-to-end substructure matching over the OpenSMILES corpus
//! (`umol-io/tests/smiles_parsing/data/opensmiles`, ~9k molecules).
//!
//! Groups:
//! - `substructure`: one-shot `Molecule::substructure_matches` across both strategies
//!   and all four subgraph-isomorphism algorithms on three patterns (`branched`, `phenol`,
//!   `bicyclic`).
//! - `substructure/prepared`: prepared-pattern (`matches`) and prepared-pattern +
//!   prepared-host (`matches_in_host`) reuse over the same corpus, plus one-time
//!   `PreparedSubstructureHost` corpus build cost.
//! - `substructure/screening`: first-match (`has_*`) and count (`count_*`) screening over
//!   the corpus (including a 0-hit tetracyclic steroid core).
//! - `substructure/symmetry`: symmetry-reduced matching (`none`, `pattern`, `host`,
//!   `both`) on symmetric ring queries against the OpenSMILES corpus and coronene
//!   (prepared hosts cache lazy `TargetSymmetry` once across iterations).
//! - `substructure/overlay_counts`: small aromatic-ring and multicenter cases measuring
//!   nonuniform, uniform, and undetermined electron counts.

use std::fs::read_to_string;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_chem::element::Element;
use umol_graph::ingest::ingest_smiles;
use umol_graph_core::SubgraphIsomorphismAlgorithm::{BitsetRi, Glasgow, RayKirsch, Vf2Rdkit};
use umol_graph_core::{
    AutomorphismAlgorithm, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::SubstructureMatchAlgorithm::{GraphAndOverlays, Incidence};
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomForm, AtomId, BondForm, ElectronCountsForm, Molecule, MoleculeEntries,
    MulticenterBondForm, NumForm, PreparedSubstructureHost, PreparedSubstructurePattern,
    SubstructureMatchAlgorithm, SubstructureMatchConfig,
};
use walkdir::WalkDir;

const SUBISO: [SubgraphIsomorphismAlgorithm; 4] = [Vf2Rdkit, RayKirsch, BitsetRi, Glasgow];
const STRATEGIES: [SubstructureMatchAlgorithm; 2] = [GraphAndOverlays, Incidence];
const SYMMETRY_MODES: [(&str, bool, bool); 4] = [
    ("none", false, false),
    ("pattern", true, false),
    ("host", false, true),
    ("both", true, true),
];

fn match_config(
    strategy: SubstructureMatchAlgorithm,
    algorithm: SubgraphIsomorphismAlgorithm,
    reduce_pattern_symmetry: bool,
    reduce_host_symmetry: bool,
) -> SubstructureMatchConfig {
    SubstructureMatchConfig {
        match_algorithm: strategy,
        subgraph_isomorphism_algorithm: algorithm,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
        reduce_pattern_symmetry,
        reduce_host_symmetry,
    }
}

fn subiso_name(algorithm: SubgraphIsomorphismAlgorithm) -> &'static str {
    match algorithm {
        Vf2Rdkit => "vf2rdkit",
        RayKirsch => "raykirsch",
        BitsetRi => "bitset_ri",
        Glasgow => "glasgow",
    }
}

fn strategy_name(strategy: SubstructureMatchAlgorithm) -> &'static str {
    match strategy {
        GraphAndOverlays => "graph_overlays",
        Incidence => "incidence",
    }
}

fn load_corpus() -> Vec<Molecule> {
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

fn carbon() -> AtomForm {
    AtomForm::from_element(Element::C)
}

fn any_bond() -> BondForm {
    BondForm::new(NumForm::Undetermined)
}

fn pattern(atoms: Vec<AtomForm>, bonds: Vec<(u32, u32, BondForm)>) -> Molecule {
    let bond_list = bonds
        .into_iter()
        .map(|(s, t, b)| (AtomId(s), AtomId(t), b))
        .collect();
    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds: bond_list,
        ..Default::default()
    })
}

/// `C(C)C(C)N` — 5 atoms, all any-bonds.
fn pattern_branched() -> Molecule {
    pattern(
        vec![
            carbon(),
            carbon(),
            carbon(),
            carbon(),
            AtomForm::from_element(Element::N),
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
fn pattern_phenol() -> Molecule {
    pattern(
        vec![
            carbon(),
            carbon(),
            carbon(),
            carbon(),
            carbon(),
            carbon(),
            AtomForm::from_element(Element::O),
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
fn pattern_bicyclic() -> Molecule {
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

/// 17-carbon tetracyclic steroid (gonane) skeleton, all any-bonds (0 hits on OpenSMILES).
fn pattern_steroid() -> Molecule {
    pattern(
        (0..17).map(|_| carbon()).collect(),
        vec![
            (0, 1, any_bond()),
            (1, 2, any_bond()),
            (2, 3, any_bond()),
            (3, 4, any_bond()),
            (4, 5, any_bond()),
            (5, 0, any_bond()),
            (4, 6, any_bond()),
            (6, 7, any_bond()),
            (7, 8, any_bond()),
            (8, 9, any_bond()),
            (9, 5, any_bond()),
            (8, 10, any_bond()),
            (10, 11, any_bond()),
            (11, 12, any_bond()),
            (12, 13, any_bond()),
            (13, 9, any_bond()),
            (12, 14, any_bond()),
            (14, 15, any_bond()),
            (15, 16, any_bond()),
            (16, 13, any_bond()),
        ],
    )
}

/// 6-carbon ring, all any-bonds (12-fold $D_6$ pattern symmetry).
fn pattern_ring6() -> Molecule {
    pattern(
        vec![carbon(); 6],
        (0..6).map(|i| (i, (i + 1) % 6, any_bond())).collect(),
    )
}

/// 24-carbon coronene skeleton ($D_6$ host symmetry: central 6-ring fused to 6 outer 6-rings).
fn host_coronene() -> Molecule {
    let mut bonds = Vec::with_capacity(30);
    for i in 0..6u32 {
        bonds.push((i, (i + 1) % 6, BondForm::from_order(1)));
    }
    for i in 0..6u32 {
        let spoke = 6 + i;
        bonds.push((i, spoke, BondForm::from_order(1)));
    }
    for i in 0..6u32 {
        let spoke_a = 6 + i;
        let spoke_b = 6 + (i + 1) % 6;
        let rim_a = 12 + 2 * i;
        let rim_b = 12 + 2 * i + 1;
        bonds.push((spoke_a, rim_a, BondForm::from_order(1)));
        bonds.push((rim_a, rim_b, BondForm::from_order(1)));
        bonds.push((rim_b, spoke_b, BondForm::from_order(1)));
    }
    pattern(vec![carbon(); 24], bonds)
}

fn substructure_benchmark(c: &mut Criterion) {
    let corpus = load_corpus();
    let patterns: [(&str, Molecule); 3] = [
        ("branched", pattern_branched()),
        ("phenol", pattern_phenol()),
        ("bicyclic", pattern_bicyclic()),
    ];

    let mut group = c.benchmark_group("substructure");
    for (pattern_name, pat) in &patterns {
        for strategy in STRATEGIES {
            for algorithm in SUBISO {
                let config = match_config(strategy, algorithm, false, false);
                let id = format!(
                    "{pattern_name}/{}/{}",
                    strategy_name(strategy),
                    subiso_name(algorithm)
                );
                group.bench_function(id, |b| {
                    b.iter(|| {
                        for target in &corpus {
                            black_box(pat.substructure_matches(target, config));
                        }
                    });
                });
            }
        }
    }
    group.finish();
}

fn prepared_substructure_benchmark(c: &mut Criterion) {
    let corpus = load_corpus();
    let prepared_corpus: Vec<PreparedSubstructureHost> =
        corpus.iter().map(PreparedSubstructureHost::new).collect();
    let patterns: [(&str, Molecule); 3] = [
        ("branched", pattern_branched()),
        ("phenol", pattern_phenol()),
        ("bicyclic", pattern_bicyclic()),
    ];

    let mut group = c.benchmark_group("substructure/prepared");
    group.bench_function("build_hosts/corpus", |b| {
        b.iter(|| {
            let hosts: Vec<PreparedSubstructureHost> =
                corpus.iter().map(PreparedSubstructureHost::new).collect();
            black_box(hosts)
        });
    });

    for (pattern_name, pat) in &patterns {
        for strategy in STRATEGIES {
            for algorithm in SUBISO {
                let prepared_pat = PreparedSubstructurePattern::new(
                    pat,
                    match_config(strategy, algorithm, false, false),
                );
                let suffix = format!(
                    "{pattern_name}/{}/{}",
                    strategy_name(strategy),
                    subiso_name(algorithm)
                );
                group.bench_function(format!("pattern/{suffix}"), |b| {
                    b.iter(|| {
                        for target in &corpus {
                            black_box(prepared_pat.matches(target));
                        }
                    });
                });
                group.bench_function(format!("both/{suffix}"), |b| {
                    b.iter(|| {
                        for host in &prepared_corpus {
                            black_box(prepared_pat.matches_in_host(host));
                        }
                    });
                });
            }
        }
    }
    group.finish();
}

fn screening_substructure_benchmark(c: &mut Criterion) {
    let corpus = load_corpus();
    let prepared_corpus: Vec<PreparedSubstructureHost> =
        corpus.iter().map(PreparedSubstructureHost::new).collect();
    let patterns: [(&str, Molecule); 4] = [
        ("branched", pattern_branched()),
        ("phenol", pattern_phenol()),
        ("bicyclic", pattern_bicyclic()),
        ("steroid", pattern_steroid()),
    ];

    let mut group = c.benchmark_group("substructure/screening");
    for (pattern_name, pat) in &patterns {
        for algorithm in SUBISO {
            let config = match_config(GraphAndOverlays, algorithm, false, false);
            let prepared_pat = PreparedSubstructurePattern::new(pat, config);
            let alg = subiso_name(algorithm);
            group.bench_function(format!("direct_has/{pattern_name}/{alg}"), |b| {
                b.iter(|| {
                    let mut hits = 0usize;
                    for target in &corpus {
                        hits += usize::from(pat.has_substructure_match(target, config));
                    }
                    black_box(hits)
                });
            });
            group.bench_function(format!("prepared_pattern_has/{pattern_name}/{alg}"), |b| {
                b.iter(|| {
                    let mut hits = 0usize;
                    for target in &corpus {
                        hits += usize::from(prepared_pat.has_match(target));
                    }
                    black_box(hits)
                });
            });
            group.bench_function(format!("prepared_both_has/{pattern_name}/{alg}"), |b| {
                b.iter(|| {
                    let mut hits = 0usize;
                    for host in &prepared_corpus {
                        hits += usize::from(prepared_pat.has_match_in_host(host));
                    }
                    black_box(hits)
                });
            });
            group.bench_function(format!("prepared_both_count/{pattern_name}/{alg}"), |b| {
                b.iter(|| {
                    let mut total = 0u64;
                    for host in &prepared_corpus {
                        total += prepared_pat.count_matches_in_host(host);
                    }
                    black_box(total)
                });
            });
        }
    }
    group.finish();
}

fn symmetry_substructure_benchmark(c: &mut Criterion) {
    let corpus = load_corpus();
    let prepared_corpus: Vec<PreparedSubstructureHost> =
        corpus.iter().map(PreparedSubstructureHost::new).collect();
    let ring6 = pattern_ring6();
    let coronene = host_coronene();
    let prepared_coronene = PreparedSubstructureHost::new(&coronene);

    let mut group = c.benchmark_group("substructure/symmetry");
    for (mode, reduce_pattern, reduce_host) in SYMMETRY_MODES {
        let prepared_ring6 = PreparedSubstructurePattern::new(
            &ring6,
            match_config(GraphAndOverlays, BitsetRi, reduce_pattern, reduce_host),
        );
        group.bench_function(format!("ring6_corpus/bitset_ri/{mode}"), |b| {
            b.iter(|| {
                let mut total = 0u64;
                for host in &prepared_corpus {
                    total += prepared_ring6.count_matches_in_host(host);
                }
                black_box(total)
            });
        });
        for algorithm in SUBISO {
            let prepared = PreparedSubstructurePattern::new(
                &ring6,
                match_config(GraphAndOverlays, algorithm, reduce_pattern, reduce_host),
            );
            group.bench_function(
                format!("ring6_coronene/{}/{mode}", subiso_name(algorithm)),
                |b| {
                    b.iter(|| black_box(prepared.matches_in_host(black_box(&prepared_coronene))));
                },
            );
        }
    }
    group.finish();
}

fn overlay_substructure_benchmark(c: &mut Criterion) {
    let aromatic = MoleculeEntries {
        atoms: vec![carbon(); 6],
        bonds: (0..6)
            .map(|i| (AtomId(i), AtomId((i + 1) % 6), any_bond()))
            .collect(),
        aromatic: vec![(
            (0..6).map(AtomId).collect(),
            AromaticSystemForm::from_electrons(vec![2, 1, 0, 1, 0, 1]),
        )],
        ..Default::default()
    };
    let multicenter = MoleculeEntries {
        atoms: vec![carbon(); 3],
        multicenter: vec![(
            (0..3).map(AtomId).collect(),
            MulticenterBondForm::from_electrons(vec![2, 0, 1]),
        )],
        ..Default::default()
    };
    let mut cases = Vec::new();
    for (name, source) in [("aromatic_ring", aromatic), ("multicenter", multicenter)] {
        for (suffix, electrons) in [
            (
                "uniform",
                ElectronCountsForm::Lit(vec![1; source.atoms.len()]),
            ),
            ("undetermined", ElectronCountsForm::Undetermined),
        ] {
            let mut entries = source.clone();
            for (_, attributes) in &mut entries.aromatic {
                attributes.electrons = electrons.clone();
            }
            for (_, attributes) in &mut entries.multicenter {
                attributes.electrons = electrons.clone();
            }
            cases.push((format!("{name}_{suffix}"), Molecule::from_entries(entries)));
        }
        cases.push((name.to_owned(), Molecule::from_entries(source)));
    }
    let mut group = c.benchmark_group("substructure/overlay_counts");
    for (name, molecule) in cases {
        for strategy in STRATEGIES {
            let config = match_config(strategy, BitsetRi, false, false);
            group.bench_function(format!("{name}/{}", strategy_name(strategy)), |b| {
                b.iter(|| black_box(molecule.substructure_matches(black_box(&molecule), config)));
            });
        }
    }
    group.finish();
}

criterion_group!(
    benches,
    substructure_benchmark,
    prepared_substructure_benchmark,
    screening_substructure_benchmark,
    symmetry_substructure_benchmark,
    overlay_substructure_benchmark
);
criterion_main!(benches);
