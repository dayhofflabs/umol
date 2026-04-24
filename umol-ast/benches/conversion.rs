//! Conversion benchmarks for `MoleculeDsl` ↔ `MoleculeAst`.
//!
//! Measures the `FromAst` and `IntoAst` paths on `MoleculeDsl` — separate
//! from `FromEdn`/`ToEdn`, which bypass these traits via `MoleculeInput`.
//! Used as a regression net for `Metadata` / DSL-newtype refactors.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, Criterion, Throughput};
use umol_ast::ast::{FromAst, IntoAst};
use umol_ast::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_edn::FromEdn;

#[path = "fixtures.rs"]
mod fixtures;
use fixtures::{
    MOL_BENZENE, MOL_INDOLE, MOL_LARGE_ALL_IDS, MOL_LARGE_NO_IDS, MOL_LARGE_PARTIAL_IDS, MOL_SMALL,
    MOL_WITH_CONSTRAINTS,
};

fn cases() -> [(&'static str, &'static str); 7] {
    [
        ("small", MOL_SMALL),
        ("benzene", MOL_BENZENE),
        ("indole", MOL_INDOLE),
        ("with_constraints", MOL_WITH_CONSTRAINTS),
        ("large_no_ids", MOL_LARGE_NO_IDS.as_str()),
        ("large_all_ids", MOL_LARGE_ALL_IDS.as_str()),
        ("large_partial_ids", MOL_LARGE_PARTIAL_IDS.as_str()),
    ]
}

fn bench_molecule_from_ast(c: &mut Criterion) {
    let cfg = MoleculeDefaults::zeroed();
    let mut g = c.benchmark_group("molecule_from_ast");
    for (label, source) in cases() {
        let ast = MoleculeDsl::from_edn_str(source).unwrap().into_parts().0;
        g.throughput(Throughput::Bytes(source.len() as u64));
        g.bench_function(label, |b| {
            b.iter(|| MoleculeDsl::from_ast(black_box(&ast), &cfg).unwrap())
        });
    }
    g.finish();
}

fn bench_molecule_into_ast(c: &mut Criterion) {
    let cfg = MoleculeDefaults::zeroed();
    let mut g = c.benchmark_group("molecule_into_ast");
    for (label, source) in cases() {
        let dsl = MoleculeDsl::from_edn_str(source).unwrap();
        g.throughput(Throughput::Bytes(source.len() as u64));
        g.bench_function(label, |b| {
            b.iter_batched(
                || dsl.clone(),
                |dsl| dsl.into_ast(&cfg).unwrap(),
                BatchSize::SmallInput,
            )
        });
    }
    g.finish();
}

criterion_group!(conversion, bench_molecule_from_ast, bench_molecule_into_ast);
criterion_main!(conversion);
