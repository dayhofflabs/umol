//! Benchmarks comparing the three molecule DSL parsing paths:
//! canonical fused (`parse_molecule_dsl`), serde streaming
//! (`parse_molecule_dsl_serde`), and native tree
//! (`parse_molecule_dsl_tree`). The fused path is the canonical entry
//! point per the Phase 3 architectural decision; the other two are
//! retained as regression detectors.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_graph::dsl::molecule::{
    parse_molecule_dsl, parse_molecule_dsl_serde, parse_molecule_dsl_tree,
};

const EMPTY: &str = r#"{:atoms [] :bonds []}"#;

const WATER: &str = r#"{:atoms ["O" "H" "H"] :bonds [[0 1 :single] [0 2 :single]]}"#;

const TAGGED_ETHANOL: &str = r#"{:atoms [[:C1 "C #h3"] [:C2 "C #h2"] [:O "O #h1"]]
                                   :bonds [[:C1 :C2 :single] [:C2 :O :single]]}"#;

const ALIASED_BENZENE: &str = r#"{:atoms [:ch :ch :ch :ch :ch :ch]
                                    :bonds [[0 1 :single] [1 2 :single] [2 3 :single]
                                            [3 4 :single] [4 5 :single] [5 0 :single]]
                                    :aliases [:ch "C #h1 #v2 #a1"]}"#;

const C20_CHAIN: &str = r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"
                                     "C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
                             :bonds [[0 1 :single] [1 2 :single] [2 3 :single] [3 4 :single]
                                     [4 5 :single] [5 6 :single] [6 7 :single] [7 8 :single]
                                     [8 9 :single] [9 10 :single] [10 11 :single] [11 12 :single]
                                     [12 13 :single] [13 14 :single] [14 15 :single] [15 16 :single]
                                     [16 17 :single] [17 18 :single] [18 19 :single]]}"#;

fn bench_molecule_dsl(c: &mut Criterion) {
    let cases: &[(&str, &str)] = &[
        ("empty", EMPTY),
        ("water", WATER),
        ("tagged_ethanol", TAGGED_ETHANOL),
        ("aliased_benzene", ALIASED_BENZENE),
        ("c20_chain", C20_CHAIN),
    ];

    let mut group = c.benchmark_group("molecule_dsl_fused");
    for (name, input) in cases {
        group.bench_function(*name, |b| {
            b.iter(|| parse_molecule_dsl(black_box(input)).unwrap())
        });
    }
    group.finish();

    let mut group = c.benchmark_group("molecule_dsl_serde");
    for (name, input) in cases {
        group.bench_function(*name, |b| {
            b.iter(|| parse_molecule_dsl_serde(black_box(input)).unwrap())
        });
    }
    group.finish();

    let mut group = c.benchmark_group("molecule_dsl_tree");
    for (name, input) in cases {
        group.bench_function(*name, |b| {
            b.iter(|| parse_molecule_dsl_tree(black_box(input)).unwrap())
        });
    }
    group.finish();
}

criterion_group!(benches, bench_molecule_dsl);
criterion_main!(benches);
