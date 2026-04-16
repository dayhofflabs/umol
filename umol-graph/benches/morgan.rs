//! Morgan fingerprint benchmark: direct MoleculeAst access vs MorganTarget view.
//!
//! Loads the conformance corpus (~9k SMILES), lifts each parsed table_ir
//! molecule into a MoleculeAst (no resolution), then benchmarks ECFP4
//! (radius 2) fingerprint computation.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use umol_graph::ast::molecule::MoleculeAst;
use umol_graph::ast::morgan::{
    morgan_direct, morgan_view, morgan_view_opt, MorganTarget, MorganTargetOpt,
};
use umol_graph::io::smiles::parse_smiles;

fn load_smiles() -> Vec<String> {
    let data_dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/smiles_parsing/data/basic_opensmiles"
    );
    let mut smiles_list = Vec::new();
    for entry in walkdir::WalkDir::new(data_dir)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.path()
                .extension()
                .map_or(false, |ext| ext == "smiles")
        })
    {
        let content = std::fs::read_to_string(entry.path()).unwrap();
        if let Some(s) = content.lines().nth(1) {
            if !s.is_empty() {
                smiles_list.push(s.to_string());
            }
        }
    }
    smiles_list
}

fn load_corpus(smiles_list: &[String]) -> Vec<MoleculeAst> {
    let mut asts = Vec::new();
    for s in smiles_list {
        if let Ok(table_mol) = parse_smiles(s) {
            asts.push(MoleculeAst::from_table_molecule(&table_mol));
        }
    }
    asts
}

fn smiles_parsing_benchmark(c: &mut Criterion) {
    let smiles_list = load_smiles();
    let corpus_size = smiles_list.len();

    c.bench_function(
        &format!("smiles_parse/{corpus_size}"),
        |b| {
            b.iter(|| {
                for s in &smiles_list {
                    black_box(parse_smiles(s));
                }
            });
        },
    );
}

fn morgan_benchmark(c: &mut Criterion) {
    let smiles_list = load_smiles();
    let asts = load_corpus(&smiles_list);
    let corpus_size = asts.len();

    let mut group = c.benchmark_group("morgan_ecfp4");

    group.bench_function(
        BenchmarkId::new("direct", corpus_size),
        |b| {
            b.iter(|| {
                for ast in &asts {
                    black_box(morgan_direct(ast, 2));
                }
            });
        },
    );

    // Pre-build views (amortized cost)
    let targets: Vec<MorganTarget> = asts.iter().map(MorganTarget::new).collect();

    group.bench_function(
        BenchmarkId::new("view", corpus_size),
        |b| {
            b.iter(|| {
                for target in &targets {
                    black_box(morgan_view(target, 2));
                }
            });
        },
    );

    group.bench_function(
        BenchmarkId::new("view_with_build", corpus_size),
        |b| {
            b.iter(|| {
                for ast in &asts {
                    let target = MorganTarget::new(ast);
                    black_box(morgan_view(&target, 2));
                }
            });
        },
    );

    let targets_opt: Vec<MorganTargetOpt> = asts.iter().map(MorganTargetOpt::new).collect();

    group.bench_function(
        BenchmarkId::new("view_opt", corpus_size),
        |b| {
            b.iter(|| {
                for target in &targets_opt {
                    black_box(morgan_view_opt(target, 2));
                }
            });
        },
    );

    group.bench_function(
        BenchmarkId::new("view_opt_with_build", corpus_size),
        |b| {
            b.iter(|| {
                for ast in &asts {
                    let target = MorganTargetOpt::new(ast);
                    black_box(morgan_view_opt(&target, 2));
                }
            });
        },
    );

    group.finish();
}

fn morgan_ecfp6_benchmark(c: &mut Criterion) {
    let smiles_list = load_smiles();
    let asts = load_corpus(&smiles_list);
    let corpus_size = asts.len();

    let mut group = c.benchmark_group("morgan_ecfp6");

    let targets_opt: Vec<MorganTargetOpt> = asts.iter().map(MorganTargetOpt::new).collect();

    group.bench_function(
        BenchmarkId::new("view_opt", corpus_size),
        |b| {
            b.iter(|| {
                for target in &targets_opt {
                    black_box(morgan_view_opt(target, 3));
                }
            });
        },
    );

    group.finish();
}

criterion_group!(benches, smiles_parsing_benchmark, morgan_benchmark, morgan_ecfp6_benchmark);
criterion_main!(benches);
