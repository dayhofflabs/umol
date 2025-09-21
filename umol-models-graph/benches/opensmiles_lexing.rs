//! Benchmarks for SMILES lexing-only and segmentation

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_models_graph::io::smiles::lexer::Lexer;
use umol_models_graph::io::smiles::iterators::Segments;

fn opensmiles_lexing(c: &mut Criterion) {
    let inputs = [
        ("short", "C1=CC=CC=C1"),
        ("with_brackets", "C[CH3]C(=O)N"),
        ("rings", "c1ccccc1.c1ccncc1"),
        (
            "long",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        ),
    ];

    // Logos lexing only
    let mut group_lex = c.benchmark_group("opensmiles_lexing/lexer_iter");
    for (name, s) in inputs.iter() {
        group_lex.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                // Iterate all tokens and discard; count to prevent optimization
                let mut n = 0usize;
                for tok in Lexer::new(input) { let _ = tok; n += 1; }
                std::hint::black_box(n);
            })
        });
    }
    group_lex.finish();

    // Segmentation (lexer + grouping) for comparison
    let mut group_seg = c.benchmark_group("opensmiles_lexing/segments");
    for (name, s) in inputs.iter() {
        group_seg.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let mut n = 0usize;
                for _ in Segments::new(input) { n += 1; }
                std::hint::black_box(n);
            })
        });
    }
    group_seg.finish();
}

criterion_group!(benches, opensmiles_lexing);
criterion_main!(benches);


