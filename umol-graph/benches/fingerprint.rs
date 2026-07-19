//! Throughput benchmark for the implemented fingerprint featurizers over the
//! OpenSMILES corpus (`umol-io/tests/smiles_parsing/data/opensmiles`,
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
use umol_ast::ast::{
    AtomDelta, AtomId, BondDelta, BondId, Delta, Deltas, MoleculeAst, ReactionAst,
};
use umol_graph::fingerprint::{
    featurize_reaction, EcfpFeaturizer, Featurizer, MorganFeaturizer, PatternFingerprinter,
    ReactionCombinator, SubstructureFeaturizer, WlFeaturizer,
};
use umol_graph::hash::RefinementXxh3Scheme;
use umol_graph::ingest::ingest_smiles;
use umol_graph_core::RefinementRounds;
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

fn ethanol_deoxygenation(molecule: &MoleculeAst) -> ReactionAst {
    ReactionAst::new(
        molecule.clone(),
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Remove {
                id: AtomId(2),
                ast: molecule.atom(AtomId(2)).ast.clone(),
            }),
            Delta::Bond(BondDelta::Remove {
                id: BondId(1),
                atoms: [AtomId(1), AtomId(2)],
                ast: molecule.bond(BondId(1)).ast.clone(),
            }),
        ]),
    )
}

fn fixture_benchmark(c: &mut Criterion) {
    let molecule = ingest_smiles("CCO").unwrap();
    let reaction = ethanol_deoxygenation(&molecule);
    let wl = WlFeaturizer {
        rounds: RefinementRounds::Fixed(WL_ROUNDS),
        scheme: RefinementXxh3Scheme::albatross(),
    };
    let ecfp = EcfpFeaturizer::new(CIRCULAR_RADIUS);
    let morgan = MorganFeaturizer::new(CIRCULAR_RADIUS);
    let pattern = PatternFingerprinter::new();
    let substructure = SubstructureFeaturizer::new(2);
    let reaction_featurizer = Featurizer::Morgan(MorganFeaturizer::new(1));
    let mut group = c.benchmark_group("fingerprint_fixture");

    group.bench_function("wl", |b| {
        b.iter(|| black_box(wl.featurize(black_box(&molecule))))
    });
    group.bench_function("ecfp", |b| {
        b.iter(|| black_box(ecfp.featurize(black_box(&molecule))))
    });
    group.bench_function("morgan", |b| {
        b.iter(|| black_box(morgan.featurize(black_box(&molecule))))
    });
    group.bench_function("pattern", |b| {
        b.iter(|| black_box(pattern.fingerprint(black_box(&molecule)).unwrap()))
    });
    group.bench_function("structural", |b| {
        b.iter(|| black_box(substructure.featurize(black_box(&molecule)).unwrap()))
    });
    group.bench_function("reaction_difference", |b| {
        b.iter(|| {
            black_box(
                featurize_reaction(
                    black_box(&reaction),
                    &reaction_featurizer,
                    ReactionCombinator::Difference,
                )
                .unwrap(),
            )
        })
    });
    group.bench_function("reaction_disjoint_union", |b| {
        b.iter(|| {
            black_box(
                featurize_reaction(
                    black_box(&reaction),
                    &reaction_featurizer,
                    ReactionCombinator::DisjointUnion,
                )
                .unwrap(),
            )
        })
    });

    group.finish();
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

criterion_group!(
    benches,
    fixture_benchmark,
    circular_benchmark,
    structural_benchmark
);
criterion_main!(benches);
