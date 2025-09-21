//! Benchmarks for SMILES linting

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_models_graph::io::smiles::linter::lint_smiles;

fn opensmiles_linting(c: &mut Criterion) {
    let inputs = [
        ("short", "C1=CC=CC=C1"),
        ("with_ws", " C - C . 1 [CH3] C ( = O ) N "),
        ("rings", "c1ccccc1.c1ccncc1.%12C1CCCCC1"),
        (
            "long",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        ),
    ];

    let mut group = c.benchmark_group("opensmiles_linting");
    for (name, s) in inputs.iter() {
        group.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let _ = lint_smiles(input);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, opensmiles_linting);
criterion_main!(benches);


