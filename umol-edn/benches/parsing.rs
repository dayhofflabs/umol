use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use serde::{Deserialize, Serialize};
use umol_edn::serde::{
    from_str, from_str_with, from_value, to_string, DynEdn, EdnList, EdnSet, EdnStreamDeserializer,
    EdnTagged,
};
use umol_edn::{
    read_all, read_string, Edn, EdnKeyRef, EdnKeyword, EdnMap, EdnSymbol, FormatConfig, ParseConfig,
};

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

fn comment_heavy() -> String {
    let mut s = String::from("[");
    for i in 0..50 {
        s.push_str(&format!("; comment line {i}\n"));
        s.push_str(&format!("#_ :discarded-{i} "));
        s.push_str(&format!("{i} "));
    }
    s.push(']');
    s
}

fn set_50() -> String {
    let mut s = String::from("#{");
    for i in 0..50 {
        if i > 0 {
            s.push(' ');
        }
        s.push_str(&format!(":item-{i}"));
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
    group.bench_function("string_unicode", |b| {
        b.iter(|| read_string(black_box("\"héllo wörld αβγ 世界\"")))
    });
    group.finish();
}

fn bench_parse_collections(c: &mut Criterion) {
    let keyword_input = keyword_heavy();
    let nested_50 = deeply_nested(50);
    let nested_100 = deeply_nested(100);
    let comments = comment_heavy();
    let set = set_50();

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
    group.bench_function("comment_heavy", |b| {
        b.iter(|| read_string(black_box(&comments)))
    });
    group.bench_function("set_50", |b| {
        b.iter(|| read_string(black_box(&set)))
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

fn bench_display(c: &mut Criterion) {
    let edn_small = read_string(MOLECULE_SMALL).unwrap();
    let edn_large = read_string(MOLECULE_LARGE).unwrap();
    let fmt = FormatConfig::default();

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

    let escape_heavy = Edn::Str(
        "line1\nline2\nline3\ttab\there\r\nquote\"end\\back\nmore\tescapes\n".into(),
    );
    group.bench_function("string_escape_heavy", |b| {
        b.iter(|| black_box(&escape_heavy).to_string())
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

#[allow(dead_code)]
#[derive(Deserialize, Serialize)]
struct Molecule {
    atoms: Vec<String>,
    bonds: Vec<(String, String, String)>,
}

const MOLECULE_SMALL_JSON: &str = r#"{"atoms":["C","O"],"bonds":[["0","1","single"]]}"#;

fn bench_deserialize(c: &mut Criterion) {
    let mut group = c.benchmark_group("deserialize");

    group.bench_function("streaming_struct", |b| {
        b.iter(|| from_str::<Molecule>(black_box(MOLECULE_SMALL)).unwrap())
    });

    group.bench_function("tree_then_struct", |b| {
        b.iter(|| {
            let edn = read_string(black_box(MOLECULE_SMALL)).unwrap();
            from_value::<Molecule>(edn).unwrap()
        })
    });

    let edn = read_string(MOLECULE_SMALL).unwrap();
    group.bench_function("from_value_only", |b| {
        b.iter(|| {
            let val = black_box(&edn).clone();
            from_value::<Molecule>(val).unwrap()
        })
    });

    group.bench_function("json_from_str_struct", |b| {
        b.iter(|| serde_json::from_str::<Molecule>(black_box(MOLECULE_SMALL_JSON)).unwrap())
    });

    group.bench_function("json_parse_to_value", |b| {
        b.iter(|| {
            serde_json::from_str::<serde_json::Value>(black_box(MOLECULE_SMALL_JSON)).unwrap()
        })
    });

    let json_val: serde_json::Value = serde_json::from_str(MOLECULE_SMALL_JSON).unwrap();
    group.bench_function("json_value_to_struct", |b| {
        b.iter(|| {
            let val = black_box(&json_val).clone();
            serde_json::from_value::<Molecule>(val).unwrap()
        })
    });

    group.finish();
}

fn bench_serialize(c: &mut Criterion) {
    let proxy = from_str::<Molecule>(MOLECULE_SMALL).unwrap();

    let mut group = c.benchmark_group("serialize");
    group.bench_function("to_string_small", |b| {
        b.iter(|| to_string(black_box(&proxy)).unwrap())
    });

    let large_proxy = Molecule {
        atoms: (0..100).map(|i| format!("C{i}")).collect(),
        bonds: (0..99)
            .map(|i| (i.to_string(), (i + 1).to_string(), "single".into()))
            .collect(),
    };
    group.bench_function("to_string_large", |b| {
        b.iter(|| to_string(black_box(&large_proxy)).unwrap())
    });

    group.finish();
}

fn many_molecules(count: usize) -> String {
    let mut s = String::with_capacity(count * (MOLECULE_SMALL.len() + 1));
    for _ in 0..count {
        s.push_str(MOLECULE_SMALL);
        s.push('\n');
    }
    s
}

fn bench_stream_throughput(c: &mut Criterion) {
    let stream_1k = many_molecules(1_000);
    let stream_10k = many_molecules(10_000);

    let mut group = c.benchmark_group("stream_throughput");

    for (label, stream) in [("1k", &stream_1k), ("10k", &stream_10k)] {
        group.throughput(Throughput::Bytes(stream.len() as u64));

        group.bench_function(format!("direct_stream_{label}"), |b| {
            b.iter(|| {
                let iter = EdnStreamDeserializer::<Molecule>::new(black_box(stream));
                let mut count = 0usize;
                for r in iter {
                    black_box(r.unwrap());
                    count += 1;
                }
                count
            })
        });

        group.bench_function(format!("read_all_then_from_value_{label}"), |b| {
            b.iter(|| {
                let all = read_all(black_box(stream)).unwrap();
                let mut count = 0usize;
                for v in all {
                    black_box(from_value::<Molecule>(v).unwrap());
                    count += 1;
                }
                count
            })
        });
    }

    group.finish();
}

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

        group.bench_function(format!("get_keyword_{size}"), |b| {
            b.iter(|| {
                for k in &kw_keys {
                    black_box(kw_map.get_ref(EdnKeyRef::keyword(k)));
                }
            })
        });

        group.bench_function(format!("get_ref_mixed_{size}"), |b| {
            b.iter(|| {
                for k in &mixed_keys {
                    black_box(mixed_map.get_ref(EdnKeyRef::from(k)));
                }
            })
        });

        group.bench_function(format!("get_exact_mixed_{size}"), |b| {
            b.iter(|| {
                for k in &mixed_keys {
                    black_box(mixed_map.get(k));
                }
            })
        });

        let miss_keys: Vec<String> = (0..size).map(|i| format!("miss-{i}")).collect();
        group.bench_function(format!("get_ref_miss_{size}"), |b| {
            b.iter(|| {
                for k in &miss_keys {
                    black_box(kw_map.get_ref(EdnKeyRef::keyword(k)));
                }
            })
        });
    }

    group.finish();
}

const VALUE_MIXED: &str = r#"{:name :salt
                              :sym chem/NaCl
                              :list (1 2 3 4 5)
                              :set #{:a :b :c :d}
                              :tagged #score 99
                              :atoms [C C C C C]
                              :nested {:k1 1 :k2 [2 3] :k3 #{:x :y}}}"#;

fn permissive_config() -> ParseConfig {
    ParseConfig {
        allow_unknown_tags: true,
        ..Default::default()
    }
}

fn bench_value_native(c: &mut Criterion) {
    let mut group = c.benchmark_group("value_native");

    let cfg = permissive_config();
    group.bench_function("parse_mixed", |b| {
        b.iter(|| DynEdn::parse_with(black_box(VALUE_MIXED), &cfg).unwrap())
    });

    let value = DynEdn::parse_with(VALUE_MIXED, &cfg).unwrap();
    group.bench_function("display_mixed", |b| {
        b.iter(|| black_box(&value).to_string())
    });

    group.bench_function("clone_mixed", |b| b.iter(|| black_box(&value).clone()));

    group.finish();
}

fn bench_value_serde(c: &mut Criterion) {
    let mut group = c.benchmark_group("value_serde");

    let cfg = permissive_config();

    group.bench_function("from_str_mixed", |b| {
        b.iter(|| from_str_with::<DynEdn>(black_box(VALUE_MIXED), &cfg).unwrap())
    });

    let value: DynEdn = from_str_with(VALUE_MIXED, &cfg).unwrap();

    group.bench_function("to_string_mixed", |b| {
        b.iter(|| to_string(black_box(&value)).unwrap())
    });

    group.bench_function("roundtrip_mixed", |b| {
        b.iter(|| {
            let v: DynEdn = from_str_with(black_box(VALUE_MIXED), &cfg).unwrap();
            to_string(&v).unwrap()
        })
    });

    group.bench_function("json_to_string_mixed", |b| {
        b.iter(|| serde_json::to_string(black_box(&value)).unwrap())
    });

    group.finish();
}

#[derive(Serialize, Deserialize)]
struct WrapperHeavy {
    name: EdnKeyword<'static>,
    ns: EdnSymbol<'static>,
    aliases: EdnList<String>,
    ids: EdnSet<i64>,
    marker: EdnTagged<String>,
}

fn wrapper_heavy_fixture() -> WrapperHeavy {
    WrapperHeavy {
        name: EdnKeyword::new("salt"),
        ns: EdnSymbol::new("chem/NaCl"),
        aliases: vec![
            "NaCl".into(),
            "halite".into(),
            "rock-salt".into(),
            "table-salt".into(),
        ]
        .into(),
        ids: (0..16).collect(),
        marker: EdnTagged::new("inst", "2026-04-08T00:00:00Z".to_string()),
    }
}

fn bench_wrappers_serde(c: &mut Criterion) {
    let mut group = c.benchmark_group("wrappers_serde");

    let fixture = wrapper_heavy_fixture();
    let serialized = to_string(&fixture).unwrap();

    group.bench_function("serialize_wrappers", |b| {
        b.iter(|| to_string(black_box(&fixture)).unwrap())
    });

    let cfg = permissive_config();

    group.bench_function("deserialize_wrappers", |b| {
        b.iter(|| from_str_with::<WrapperHeavy>(black_box(&serialized), &cfg).unwrap())
    });

    group.bench_function("roundtrip_wrappers", |b| {
        b.iter(|| {
            let s = to_string(black_box(&fixture)).unwrap();
            from_str_with::<WrapperHeavy>(&s, &cfg).unwrap()
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_parse_atoms,
    bench_parse_collections,
    bench_read_all,
    bench_display,
    bench_roundtrip,
    bench_deserialize,
    bench_serialize,
    bench_stream_throughput,
    bench_map_lookup,
    bench_value_native,
    bench_value_serde,
    bench_wrappers_serde,
);
criterion_main!(benches);
