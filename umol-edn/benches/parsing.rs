use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use serde::Deserialize;
use umol_edn::{from_str, read_all, read_string, EdnFormatter};

const MOLECULE_SMALL: &str = r#"{:atoms [C O] :bonds [["0" "1" :single]]}"#;

const MOLECULE_LARGE: &str = r#"
{:atoms [C C C C C C C C C C]
 :bonds [[0 1 :single] [1 2 :single] [2 3 :single] [3 4 :single]
         [4 5 :single] [5 0 :single] [3 6 :single] [6 7 :single]
         [7 8 :single] [8 9 :single] [9 4 :single]]
 :config-overrides {:atom {:implicit-h-mode :normal :charge-mode :zero :aromatic-valence-mode :aromatic}}
 :context {:atom_candidates {}
           :atom_aromatic_hints {0 true 1 true 2 true 3 true 4 true 5 true 6 true 7 true 8 true 9 true}
           :bond_aromatic_hints {0 true 1 true 2 true 3 true 4 true 5 true 6 true 7 true 8 true 9 true 10 true}
           :atom_normal_implicit_hydrogens [0 1 2 3 4 5 6 7 8 9]}}
"#;

fn keyword_heavy() -> String {
    let mut s = String::from("{");
    for i in 0..200 {
        if i > 0 {
            s.push(' ');
        }
        s.push_str(&format!(":key-{i} {i}"));
    }
    s.push('}');
    s
}

fn deeply_nested(depth: usize) -> String {
    let open: String = "[".repeat(depth);
    let close: String = "]".repeat(depth);
    format!("{open}1{close}")
}

fn many_values(count: usize) -> String {
    (0..count)
        .map(|i| format!(":val-{i}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn bench_parse_atoms(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse_atoms");
    group.bench_function("nil", |b| b.iter(|| read_string(black_box("nil"))));
    group.bench_function("int", |b| b.iter(|| read_string(black_box("12345"))));
    group.bench_function("float", |b| b.iter(|| read_string(black_box("3.14159e10"))));
    group.bench_function("string_short", |b| {
        b.iter(|| read_string(black_box(r#""hello""#)))
    });
    group.bench_function("string_escaped", |b| {
        b.iter(|| read_string(black_box(r#""line\nwith\ttabs\u0041""#)))
    });
    group.bench_function("keyword", |b| {
        b.iter(|| read_string(black_box(":ns/some-name")))
    });
    group.bench_function("symbol", |b| {
        b.iter(|| read_string(black_box("my.ns/some-symbol")))
    });
    group.finish();
}

fn bench_parse_collections(c: &mut Criterion) {
    let keyword_input = keyword_heavy();
    let nested_50 = deeply_nested(50);
    let nested_200 = deeply_nested(200);

    let mut group = c.benchmark_group("parse_collections");
    group.bench_function("molecule_small", |b| {
        b.iter(|| read_string(black_box(MOLECULE_SMALL)))
    });
    group.bench_function("molecule_large", |b| {
        b.iter(|| read_string(black_box(MOLECULE_LARGE)))
    });
    group.bench_function("keyword_map_200", |b| {
        b.iter(|| read_string(black_box(&keyword_input)))
    });
    group.bench_function("nested_50", |b| {
        b.iter(|| read_string(black_box(&nested_50)))
    });
    group.bench_function("nested_200", |b| {
        b.iter(|| read_string(black_box(&nested_200)))
    });
    group.finish();
}

fn bench_read_all(c: &mut Criterion) {
    let stream_100 = many_values(100);
    let stream_1000 = many_values(1000);

    let mut group = c.benchmark_group("read_all");
    group.bench_function("100_values", |b| {
        b.iter(|| read_all(black_box(&stream_100)))
    });
    group.bench_function("1000_values", |b| {
        b.iter(|| read_all(black_box(&stream_1000)))
    });
    group.finish();
}

fn bench_roundtrip(c: &mut Criterion) {
    let mut group = c.benchmark_group("roundtrip");
    group.bench_function("molecule_large", |b| {
        b.iter(|| {
            let edn = read_string(black_box(MOLECULE_LARGE)).unwrap();
            let s = edn.to_string();
            let _ = read_string(black_box(&s)).unwrap();
        })
    });
    group.finish();
}

fn bench_display(c: &mut Criterion) {
    let edn_small = read_string(MOLECULE_SMALL).unwrap();
    let edn_large = read_string(MOLECULE_LARGE).unwrap();
    let fmt = EdnFormatter::default();

    let mut group = c.benchmark_group("display");
    group.bench_function("to_string_small", |b| {
        b.iter(|| black_box(&edn_small).to_string())
    });
    group.bench_function("to_string_large", |b| {
        b.iter(|| black_box(&edn_large).to_string())
    });
    group.bench_function("formatter_small", |b| {
        b.iter(|| black_box(&edn_small).to_string_with(&fmt))
    });
    group.bench_function("formatter_large", |b| {
        b.iter(|| black_box(&edn_large).to_string_with(&fmt))
    });
    group.finish();
}

#[allow(dead_code)]
#[derive(Deserialize)]
struct MoleculeProxy {
    atoms: Vec<String>,
    bonds: Vec<(String, String, String)>,
}

// TODO: Remove JSON reference benchmarks (and serde_json dependency) once EDN parsing stabilizes
const MOLECULE_SMALL_JSON: &str = r#"{"atoms":["C","O"],"bonds":[["0","1","single"]]}"#;

fn bench_serde(c: &mut Criterion) {
    let mut group = c.benchmark_group("serde");

    // EDN: parse to Edn value tree.
    group.bench_function("edn_parse_to_edn", |b| {
        b.iter(|| read_string(black_box(MOLECULE_SMALL)).unwrap())
    });

    // EDN: streaming from_str (bypasses Edn tree).
    group.bench_function("edn_from_str_struct", |b| {
        b.iter(|| from_str::<MoleculeProxy>(black_box(MOLECULE_SMALL)).unwrap())
    });

    // EDN: from pre-parsed Edn (Edn→struct cost only).
    let edn = read_string(MOLECULE_SMALL).unwrap();
    group.bench_function("edn_to_struct", |b| {
        b.iter(|| {
            let val = black_box(&edn).clone();
            umol_edn::from_value::<MoleculeProxy>(val).unwrap()
        })
    });

    // JSON reference: streaming from_str.
    group.bench_function("json_from_str_struct", |b| {
        b.iter(|| serde_json::from_str::<MoleculeProxy>(black_box(MOLECULE_SMALL_JSON)).unwrap())
    });

    // JSON reference: parse to Value.
    group.bench_function("json_parse_to_value", |b| {
        b.iter(|| {
            serde_json::from_str::<serde_json::Value>(black_box(MOLECULE_SMALL_JSON)).unwrap()
        })
    });

    // JSON reference: Value→struct (two-step).
    let json_val: serde_json::Value = serde_json::from_str(MOLECULE_SMALL_JSON).unwrap();
    group.bench_function("json_value_to_struct", |b| {
        b.iter(|| {
            let val = black_box(&json_val).clone();
            serde_json::from_value::<MoleculeProxy>(val).unwrap()
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_parse_atoms,
    bench_parse_collections,
    bench_read_all,
    bench_roundtrip,
    bench_display,
    bench_serde,
);
criterion_main!(benches);
