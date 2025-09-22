//! Benchmarks for SMILES parsing

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_models_graph::io::smiles::lexer::Lexer;
use umol_models_graph::io::smiles::linter::{lint_smiles, lint_smiles_parse, lint_smiles_parse_fast};
use umol_models_graph::io::smiles::state::{ParseState, ParserMode};
use umol_models_graph::io::smiles::parser::grammar::MoleculeParser;
// (removed duplicate ParseState import)

fn opensmiles_parsing(c: &mut Criterion) {
    // Chain-only corpus, bare atoms (organic-only mix omitting bare H)
    const MIX7: &str = "CNOFPSI"; // simple deterministic mix
    const CHAIN_MIX_20: &str = "CNOFPSICNOFPSICNOFPS"; // 7+7+6 = 20 atoms
    const CHAIN_MIX_50: &str = "CNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSIC"; // 7*7 + 1 = 50
    const CHAIN_MIX_100: &str = "CNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICN"; // 7*14 + 2 = 100

    let inputs = [
        ("chain_empty", ""),
        ("chain_c_1", "C"),
        ("chain_c_5", "CCCCC"),
        ("chain_c_10", "CCCCCCCCCC"),
        ("chain_c_50", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        ("chain_c_100", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        ("chain_c_1000", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        ("chain_mix_20", CHAIN_MIX_20),
        ("chain_mix_50", CHAIN_MIX_50),
        ("chain_mix_100", CHAIN_MIX_100),
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
