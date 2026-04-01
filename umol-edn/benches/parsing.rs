use criterion::{black_box, criterion_group, criterion_main, Criterion};
use umol_edn::{read_all, read_string, EdnFormatter};

const MOLECULE_SMALL: &str = r#"{:atoms [C O] :bonds [[:0 :1 :single]]}"#;

const MOLECULE_LARGE: &str = r#"
{:atoms [C C C C C C C C C C]
 :bonds [[:0 :1 :single] [:1 :2 :single] [:2 :3 :single] [:3 :4 :single]
         [:4 :5 :single] [:5 :0 :single] [:3 :6 :single] [:6 :7 :single]
         [:7 :8 :single] [:8 :9 :single] [:9 :4 :single]]
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

criterion_group!(
    benches,
    bench_parse_atoms,
    bench_parse_collections,
    bench_read_all,
    bench_roundtrip,
);
criterion_main!(benches);
