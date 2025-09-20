//! Benchmarks for SMILES parsing

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_models_graph::io::smiles::lexer::Lexer;
use umol_models_graph::io::smiles::linter::lint_smiles_parse;
use umol_models_graph::io::smiles::parser::grammar::MoleculeParser;
use umol_models_graph::io::smiles::state::ParseState;

fn smiles_parsing(c: &mut Criterion) {
    let inputs = [
        ("short", "C1=CC=CC=C1"),
        ("with_brackets", "C[CH3]C(=O)N"),
        ("rings", "c1ccccc1.c1ccncc1"),
        (
            "long",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        ),
    ];

    // Parse only
    let mut group_parse = c.benchmark_group("smiles_parsing/parse_only");
    for (name, s) in inputs.iter() {
        group_parse.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let mut state = ParseState::default();
                let parser = MoleculeParser::new();
                let lexer = Lexer::new(input);
                let _ = parser.parse(&mut state, lexer);
            })
        });
    }
    group_parse.finish();

    // Lint + parse (lint_smiles_parse)
    let mut group_full = c.benchmark_group("smiles_parsing/lint_plus_parse");
    for (name, s) in inputs.iter() {
        group_full.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let _ = lint_smiles_parse(input);
            })
        });
    }
    group_full.finish();
}

criterion_group!(benches, smiles_parsing);
criterion_main!(benches);


