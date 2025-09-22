//! Benchmarks for SMILES parsing

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_models_graph::io::smiles::lexer::Lexer;
use umol_models_graph::io::smiles::linter::{lint_smiles, lint_smiles_parse, lint_smiles_parse_fast};
use umol_models_graph::io::smiles::state::{ParseState, ParserMode};
use umol_models_graph::io::smiles::parser::grammar::MoleculeParser;
// (removed duplicate ParseState import)

fn opensmiles_parsing(c: &mut Criterion) {
    let inputs = [
        ("short", "C1=CC=CC=C1"),
        ("with_brackets", "C[CH3]C(=O)N"),
        ("rings", "c1ccccc1.c1ccncc1"),
        (
            "long",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        ),
    ];

    // Lex-only baseline
    let mut group_lex = c.benchmark_group("opensmiles_parsing/lex_only");
    for (name, s) in inputs.iter() {
        group_lex.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let mut n = 0usize;
                for tok in Lexer::new(input) { let _ = tok; n += 1; }
                std::hint::black_box(n);
            })
        });
    }
    group_lex.finish();

    // Parse only
    let mut group_parse = c.benchmark_group("opensmiles_parsing/parse_only");
    for (name, s) in inputs.iter() {
        group_parse.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let mut state = ParseState::default();
                let parser = MoleculeParser::new();
                let lexer = Lexer::new(input);
                let _ = parser.parse(&mut state, lexer);
            })
        });
    }
    group_parse.finish();

    // Lint + parse (lint_smiles_parse)
    let mut group_full = c.benchmark_group("opensmiles_parsing/lint_plus_parse");
    for (name, s) in inputs.iter() {
        group_full.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let _ = lint_smiles_parse(black_box(input));
            })
        });
    }
    group_full.finish();

    // Lint + parser fast mode (no IR)
    let mut group_fast = c.benchmark_group("opensmiles_parsing/lint_plus_parse_fast");
    for (name, s) in inputs.iter() {
        group_fast.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let _ = lint_smiles_parse_fast(black_box(input));
            })
        });
    }
    group_fast.finish();

    // Linter only
    let mut group_lint = c.benchmark_group("opensmiles_parsing/lint_only");
    for (name, s) in inputs.iter() {
        group_lint.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let _ = lint_smiles(black_box(input));
            })
        });
    }
    group_lint.finish();

    // Parser minimal (no IR, no diags, increment counter in every action)
    let mut group_min = c.benchmark_group("opensmiles_parsing/parse_minimal");
    for (name, s) in inputs.iter() {
        group_min.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let mut state = ParseState::with_mode(ParserMode::Minimal);
                let parser = MoleculeParser::new();
                let lexer = Lexer::new(input);
                let _ = parser.parse(&mut state, lexer);
                std::hint::black_box(state.action_count);
            })
        });
    }
    group_min.finish();
}

criterion_group!(benches, opensmiles_parsing);
criterion_main!(benches);
