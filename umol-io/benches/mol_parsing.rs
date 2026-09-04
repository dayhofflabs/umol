//! Benchmarks for top-level MOL and SDF parsing.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_io::ctfile::parser::{
    parse_extended_mol_bytes, parse_mol_bytes_to_table_ir, parse_sdf_bytes,
};

const CAFFEINE_MOL: &[u8] =
    include_bytes!("../tests/mol_parsing/data/molecule/chemdoodle/caffeine.mol");
const COPOLYMER_MOL: &[u8] =
    include_bytes!("../tests/mol_parsing/data/extended_molecule/rdkit/Sgroups_Copolymer_01.mol");
const COMPONENTS_SDF: &[u8] =
    include_bytes!("../tests/sdf_parsing/data/molecule/wwpdb/components-pub-10.sdf");

fn mol_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("mol_parsing");

    group.bench_with_input(
        BenchmarkId::new("basic_mol", "caffeine"),
        &CAFFEINE_MOL,
        |b, input| b.iter(|| parse_mol_bytes_to_table_ir(black_box(input)).unwrap()),
    );
    group.bench_with_input(
        BenchmarkId::new("extended_mol", "copolymer_sgroup"),
        &COPOLYMER_MOL,
        |b, input| b.iter(|| parse_extended_mol_bytes(black_box(input)).unwrap()),
    );
    group.bench_with_input(
        BenchmarkId::new("sdf", "components_10"),
        &COMPONENTS_SDF,
        |b, input| b.iter(|| parse_sdf_bytes(black_box(input)).unwrap()),
    );

    group.finish();
}

criterion_group!(benches, mol_parsing);
criterion_main!(benches);
