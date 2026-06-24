//! Throughput benchmark for the implemented fingerprint featurizers over the
//! basic_opensmiles corpus (`umol-io/tests/smiles_parsing/data/basic_opensmiles`,
//! ~9k molecules), parsed and resolved to ground.
//!
//! Two groups. `fingerprint` runs the circular featurizers — WL, ECFP4, Morgan
//! (radius 2) — over the whole corpus; they are linear in graph size and fast.
//! `fingerprint_structural` runs the heavier methods — the unhashed substructure
//! screen (a nauty canonical form per enumerated subgraph) and the RDKit
//! pattern-fingerprint replica (13 subgraph-isomorphism matches per molecule) —
//! over a bounded prefix of the corpus with reduced sampling, since their
//! per-molecule cost is far higher. Adjust the consts to widen either sweep.

use std::fs::read_to_string;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_ast::ast::MoleculeAst;
use umol_graph::fingerprint::{
    EcfpFeaturizer, MorganFeaturizer, PatternFingerprinter, SubstructureFeaturizer, WlFeaturizer,
};
use umol_graph::parse::parse_smiles;
use umol_graph_core::{RefinementRounds, RefinementXxh3Scheme};
use walkdir::WalkDir;

/// Circular fingerprint radius (ECFP4 / Morgan radius 2).
const CIRCULAR_RADIUS: u32 = 2;
/// WL refinement rounds.
const WL_ROUNDS: u32 = 3;
/// Bond bound for the substructure screen; its connected-subgraph enumeration grows
/// combinatorially, so a modest diameter keeps the sweep tractable.
const SUBSTRUCTURE_MAX_BONDS: u32 = 4;
/// Corpus prefix size for the heavier structural methods.
const STRUCTURAL_PREFIX: usize = 300;
/// Reduced sample count for the structural group (each iteration is expensive).
const STRUCTURAL_SAMPLES: usize = 10;

fn load_corpus() -> Vec<MoleculeAst> {
    let dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../umol-io/tests/smiles_parsing/data/basic_opensmiles"
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
                if let Ok(molecule) = parse_smiles(smiles) {
                    molecules.push(molecule);
                }
            }
        }
    }
    molecules
}

fn circular_benchmark(c: &mut Criterion) {
    let corpus = load_corpus();
    let size = corpus.len();
    let mut group = c.benchmark_group("fingerprint");

    let wl = WlFeaturizer {
        rounds: RefinementRounds::Fixed(WL_ROUNDS),
        scheme: RefinementXxh3Scheme::albatross(),
    };
    group.bench_function(BenchmarkId::new("wl", size), |b| {
        b.iter(|| {
            for molecule in &corpus {
                black_box(wl.featurize(molecule));
            }
        });
    });

    let ecfp = EcfpFeaturizer::new(CIRCULAR_RADIUS);
    group.bench_function(BenchmarkId::new("ecfp", size), |b| {
        b.iter(|| {
            for molecule in &corpus {
                black_box(ecfp.featurize(molecule));
            }
        });
    });

    let morgan = MorganFeaturizer::new(CIRCULAR_RADIUS);
    group.bench_function(BenchmarkId::new("morgan", size), |b| {
        b.iter(|| {
            for molecule in &corpus {
                black_box(morgan.featurize(molecule));
            }
        });
    });

    group.finish();
}

fn structural_benchmark(c: &mut Criterion) {
    let corpus = load_corpus();
    let corpus = &corpus[..STRUCTURAL_PREFIX.min(corpus.len())];
    let size = corpus.len();
    let mut group = c.benchmark_group("fingerprint_structural");
    group.sample_size(STRUCTURAL_SAMPLES);

    let substructure = SubstructureFeaturizer::new(SUBSTRUCTURE_MAX_BONDS);
    group.bench_function(BenchmarkId::new("substructure", size), |b| {
        b.iter(|| {
            for molecule in corpus {
                black_box(substructure.featurize(molecule).unwrap());
            }
        });
    });

    let pattern = PatternFingerprinter::new();
    group.bench_function(BenchmarkId::new("pattern", size), |b| {
        b.iter(|| {
            for molecule in corpus {
                black_box(pattern.fingerprint(molecule).unwrap());
            }
        });
    });

    group.finish();
}

criterion_group!(benches, circular_benchmark, structural_benchmark);
criterion_main!(benches);
