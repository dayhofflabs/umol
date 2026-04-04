use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use serde::{Deserialize, Serialize};
use umol_edn::{from_str, read_all, read_string, Edn, EdnFormatter, EdnKeyRef, EdnMap};

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

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

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
    let nested_100 = deeply_nested(100);

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
    group.bench_function("nested_100", |b| {
        b.iter(|| read_string(black_box(&nested_100)))
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

// ---------------------------------------------------------------------------
// Display / formatting
// ---------------------------------------------------------------------------

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

    let float_heavy = Edn::Vector(
        (0..100)
            .map(|i| Edn::Float(i as f64 * 0.1 + 0.001))
            .collect::<Vec<_>>()
            .into(),
    );
    group.bench_function("float_heavy_100", |b| {
        b.iter(|| black_box(&float_heavy).to_string())
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

// ---------------------------------------------------------------------------
// Serde deserialization
// ---------------------------------------------------------------------------

#[allow(dead_code)]
#[derive(Deserialize, Serialize)]
struct MoleculeProxy {
    atoms: Vec<String>,
    bonds: Vec<(String, String, String)>,
}

fn bench_deserialize(c: &mut Criterion) {
    let mut group = c.benchmark_group("deserialize");

    group.bench_function("streaming_struct", |b| {
        b.iter(|| from_str::<MoleculeProxy>(black_box(MOLECULE_SMALL)).unwrap())
    });

    group.bench_function("tree_then_struct", |b| {
        b.iter(|| {
            let edn = read_string(black_box(MOLECULE_SMALL)).unwrap();
            umol_edn::from_value::<MoleculeProxy>(edn).unwrap()
        })
    });

    let edn = read_string(MOLECULE_SMALL).unwrap();
    group.bench_function("from_value_only", |b| {
        b.iter(|| {
            let val = black_box(&edn).clone();
            umol_edn::from_value::<MoleculeProxy>(val).unwrap()
        })
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Serde serialization
// ---------------------------------------------------------------------------

fn bench_serialize(c: &mut Criterion) {
    let proxy = from_str::<MoleculeProxy>(MOLECULE_SMALL).unwrap();

    let mut group = c.benchmark_group("serialize");
    group.bench_function("to_string_small", |b| {
        b.iter(|| umol_edn::to_string(black_box(&proxy)).unwrap())
    });

    let large_proxy = MoleculeProxy {
        atoms: (0..100).map(|i| format!("C{i}")).collect(),
        bonds: (0..99)
            .map(|i| (i.to_string(), (i + 1).to_string(), "single".into()))
            .collect(),
    };
    group.bench_function("to_string_large", |b| {
        b.iter(|| umol_edn::to_string(black_box(&large_proxy)).unwrap())
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Map lookup
// ---------------------------------------------------------------------------

fn make_keyword_map(n: usize) -> (EdnMap<'static>, Vec<String>) {
    let mut m = EdnMap::with_capacity(n);
    let mut keys = Vec::with_capacity(n);
    for i in 0..n {
        let k = format!("key-{i}");
        m.insert(Edn::keyword(&k).into_owned(), Edn::Int(i as i64));
        keys.push(k);
    }
    (m, keys)
}

fn make_mixed_map(n: usize) -> (EdnMap<'static>, Vec<Edn<'static>>) {
    let mut m = EdnMap::with_capacity(n);
    let mut keys = Vec::with_capacity(n);
    for i in 0..n {
        let key = match i % 3 {
            0 => Edn::keyword(&format!("k-{i}")).into_owned(),
            1 => Edn::Int(i as i64),
            _ => Edn::Str(format!("s-{i}").into()),
        };
        m.insert(key.clone(), Edn::Int(i as i64));
        keys.push(key);
    }
    (m, keys)
}

fn bench_map_lookup(c: &mut Criterion) {
    let mut group = c.benchmark_group("map_lookup");

    for &size in &[8, 32, 128, 512, 1024] {
        let (kw_map, kw_keys) = make_keyword_map(size);
        let (mixed_map, mixed_keys) = make_mixed_map(size);

        group.bench_function(&format!("get_keyword_{size}"), |b| {
            b.iter(|| {
                for k in &kw_keys {
                    black_box(kw_map.get_ref(EdnKeyRef::keyword(k)));
                }
            })
        });

        group.bench_function(&format!("get_ref_mixed_{size}"), |b| {
            b.iter(|| {
                for k in &mixed_keys {
                    black_box(mixed_map.get_ref(EdnKeyRef::from(k)));
                }
            })
        });

        group.bench_function(&format!("get_exact_mixed_{size}"), |b| {
            b.iter(|| {
                for k in &mixed_keys {
                    black_box(mixed_map.get(k));
                }
            })
        });

        let miss_keys: Vec<String> = (0..size).map(|i| format!("miss-{i}")).collect();
        group.bench_function(&format!("get_ref_miss_{size}"), |b| {
            b.iter(|| {
                for k in &miss_keys {
                    black_box(kw_map.get_ref(EdnKeyRef::keyword(k)));
                }
            })
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------

criterion_group!(
    benches,
    bench_parse_atoms,
    bench_parse_collections,
    bench_read_all,
    bench_display,
    bench_roundtrip,
    bench_deserialize,
    bench_serialize,
    bench_map_lookup,
);
criterion_main!(benches);
