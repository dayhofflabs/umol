//! End-to-end SMILES ingest benchmarks — parse, raise, and resolve — over
//! representative single molecules: an unbranched alkane (localized valence
//! only), benzene and pyridine (joint aromatic selection), bare methane
//! (plural admission collapsed by the tie-break), and the bare-nitrogen
//! fused heteroaromatics quinoline and purine (assignment search over a
//! flexible fused component).

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_graph::ingest::ingest_smiles;

fn bench_ingest_smiles(c: &mut Criterion) {
    let mut g = c.benchmark_group("resolve_ingest_smiles");

    for (name, smiles) in [
        ("methane", "C"),
        ("octane", "CCCCCCCC"),
        ("benzene", "c1ccccc1"),
        ("pyridine", "c1ccncc1"),
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("quinoline", "c1ccc2ccccc2n1"),
        ("purine", "c1ncc2ncnc2n1"),
    ] {
        g.bench_function(name, |b| {
            b.iter(|| ingest_smiles(black_box(smiles)).unwrap())
        });
    }

    g.finish();
}

criterion_group!(resolve, bench_ingest_smiles);
criterion_main!(resolve);
