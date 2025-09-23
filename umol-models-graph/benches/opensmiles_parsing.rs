//! Benchmarks for SMILES parsing

use std::hint::black_box;

use bstr::{ByteSlice, ByteVec};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::prelude::*;
use rand::seq::SliceRandom;
use rand_chacha::ChaCha12Rng;
use umol_models_graph::io::smiles::lexer::Lexer;
use umol_models_graph::io::smiles::lexer_old as legacy_lex;
use umol_models_graph::io::smiles::parser::grammar::MoleculeParser;
use umol_models_graph::io::smiles::state::{ParseState, ParserMode};

fn opensmiles_parsing(c: &mut Criterion) {
    // Chain-only corpus, bare atoms (organic-only mix omitting bare H)

    let mut rng = ChaCha12Rng::seed_from_u64(20250922);
    let mut mix_20 = Vec::from_slice(b"CNOFPSICNOFPSICNOFPS");
    mix_20.shuffle(&mut rng);
    let mut mix_50 = Vec::from_slice(b"CNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSIC");
    mix_50.shuffle(&mut rng);
    let mut mix_100 = Vec::from_slice(b"CNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICN");
    mix_100.shuffle(&mut rng);

    let inputs = [
        ("chain_empty", &b""[..]),
        ("chain_c_1", &b"C"[..]),
        ("chain_c_5", &b"CCCCC"[..]),
        ("chain_c_10", &b"CCCCCCCCCC"[..]),
        ("chain_c_50", &b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("chain_c_100", &b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("chain_c_1000", &b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("chain_mix_20", &mix_20[..]),
        ("chain_mix_50", &mix_50[..]),
        ("chain_mix_100", &mix_100[..]),
    ];

    // Lex-only baseline
    let mut group_lex = c.benchmark_group("opensmiles_parsing/lex_only");
    for (name, s) in inputs.iter() {
        group_lex.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let mut n = 0usize;
                for tok in Lexer::new(input.as_bytes()) {
                    let _ = tok;
                    n += 1;
                }
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
                let lexer = Lexer::new(input.as_bytes());
                let _ = parser.parse(&mut state, lexer);
            })
        });
    }
    group_parse.finish();

    // Parser minimal (no IR, no diags, increment counter in every action)
    let mut group_min = c.benchmark_group("opensmiles_parsing/parse_minimal");
    for (name, s) in inputs.iter() {
        group_min.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let mut state = ParseState::with_mode(ParserMode::Minimal);
                let parser = MoleculeParser::new();
                let lexer = Lexer::new(input.as_bytes());
                let _ = parser.parse(&mut state, lexer);
                std::hint::black_box(state.action_count);
            })
        });
    }
    group_min.finish();
}

criterion_group!(benches, opensmiles_parsing);
criterion_main!(benches);
