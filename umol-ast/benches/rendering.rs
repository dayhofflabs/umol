//! Rendering benchmarks for the DSL layer.
//!
//! Covers `MoleculeDsl::to_edn` (AST + Metadata → `Edn` tree) over the same
//! fixture set as the parsing benches. The render path exercises the
//! per-entity id lookups (`atom_id`, `bond_id`, ...) on every ref — the
//! dominant cost on id-heavy molecules.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use umol_ast::dsl::molecule::MoleculeDsl;
use umol_edn::{FromEdn, ToEdn};

#[path = "fixtures.rs"]
mod fixtures;
use fixtures::{
    MOL_BENZENE, MOL_INDOLE, MOL_LARGE_ALL_IDS, MOL_LARGE_NO_IDS, MOL_LARGE_PARTIAL_IDS, MOL_SMALL,
    MOL_WITH_CONSTRAINTS,
};

fn bench_case(c: &mut Criterion, label: &str, source: &str) {
    let dsl: MoleculeDsl = MoleculeDsl::from_edn_str(source).unwrap();
    let mut g = c.benchmark_group("molecule_render");
    g.throughput(Throughput::Bytes(source.len() as u64));
    g.bench_function(label, |b| b.iter(|| black_box(&dsl).to_edn()));
    g.finish();
}

fn bench_molecule_render(c: &mut Criterion) {
    bench_case(c, "small", MOL_SMALL);
    bench_case(c, "benzene", MOL_BENZENE);
    bench_case(c, "indole", MOL_INDOLE);
    bench_case(c, "with_constraints", MOL_WITH_CONSTRAINTS);
    bench_case(c, "large_no_ids", MOL_LARGE_NO_IDS.as_str());
    bench_case(c, "large_all_ids", MOL_LARGE_ALL_IDS.as_str());
    bench_case(c, "large_partial_ids", MOL_LARGE_PARTIAL_IDS.as_str());
}

criterion_group!(render, bench_molecule_render);
criterion_main!(render);
