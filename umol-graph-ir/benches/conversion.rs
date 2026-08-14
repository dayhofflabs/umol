//! Conversion benchmarks for `MoleculeDsl` ↔ `Molecule`.
//!
//! Measures the `FromIr` and `IntoIr` paths on `MoleculeDsl` — separate
//! from `FromEdn`/`ToEdn`, which bypass these traits via `MoleculeInput`.
//! Used as a regression net for `MoleculeMetadata` / DSL-newtype refactors.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, Criterion, Throughput};
use umol_edn::FromEdn;
use umol_graph_ir::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_graph_ir::ir::{FromIr, IntoIr};

#[path = "fixtures.rs"]
mod fixtures;
use fixtures::{
    MOL_BENZENE, MOL_DIBORANE, MOL_INDOLE, MOL_LARGE_ALL_IDS, MOL_LARGE_NO_IDS,
    MOL_LARGE_PARTIAL_IDS, MOL_SMALL, MOL_WITH_CONSTRAINTS,
};

fn cases() -> [(&'static str, &'static str); 8] {
    [
        ("small", MOL_SMALL),
        ("benzene", MOL_BENZENE),
        ("indole", MOL_INDOLE),
        ("diborane", MOL_DIBORANE),
        ("with_constraints", MOL_WITH_CONSTRAINTS),
        ("large_no_ids", MOL_LARGE_NO_IDS.as_str()),
        ("large_all_ids", MOL_LARGE_ALL_IDS.as_str()),
        ("large_partial_ids", MOL_LARGE_PARTIAL_IDS.as_str()),
    ]
}

fn bench_molecule_from_ir(c: &mut Criterion) {
    let cfg = MoleculeDefaults::ground();
    let mut g = c.benchmark_group("molecule_from_ir");
    for (label, source) in cases() {
        let molecule = MoleculeDsl::from_edn_str(source).unwrap().into_parts().0;
        g.throughput(Throughput::Bytes(source.len() as u64));
        g.bench_function(label, |b| {
            b.iter(|| MoleculeDsl::from_ir(black_box(&molecule), &cfg))
        });
    }
    g.finish();
}

fn bench_molecule_into_ir(c: &mut Criterion) {
    let cfg = MoleculeDefaults::zeroed();
    let mut g = c.benchmark_group("molecule_into_ir");
    for (label, source) in cases() {
        let dsl = MoleculeDsl::from_edn_str(source).unwrap();
        g.throughput(Throughput::Bytes(source.len() as u64));
        g.bench_function(label, |b| {
            b.iter_batched(
                || dsl.clone(),
                |dsl| dsl.into_ir(&cfg),
                BatchSize::SmallInput,
            )
        });
    }
    g.finish();
}

criterion_group!(conversion, bench_molecule_from_ir, bench_molecule_into_ir);
criterion_main!(conversion);
