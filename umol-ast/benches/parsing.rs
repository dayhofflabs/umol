//! Parsing benchmarks for the DSL layer.
//!
//! Covers the two entry paths for every parseable DSL type:
//!
//! - Tree: `read_string` (umol-edn) → `<T as FromEdn>::from_edn`.
//! - Streaming: `<T as FromEdn>::from_edn_str`.
//!
//! Groups mirror the type hierarchy: entity-string DSLs, tree-shaped
//! constraint DSLs, and the full molecule map DSL at a few sizes.

use std::hint::black_box;

use criterion::measurement::WallTime;
use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use umol_ast::dsl::aromatic::AromaticSystemDsl;
use umol_ast::dsl::atom::AtomDsl;
use umol_ast::dsl::bond::BondDsl;
use umol_ast::dsl::constraint::{ConstraintDsl, ConstraintsDsl};
use umol_ast::dsl::dative::DativeBondDsl;
use umol_ast::dsl::molecule::MoleculeDsl;
use umol_ast::dsl::multicenter::MulticenterBondDsl;
use umol_ast::dsl::noncovalent::NoncovalentBondDsl;
use umol_edn::{read_string, FromEdn};

fn bench_pair<T>(group: &mut criterion::BenchmarkGroup<'_, WallTime>, label: &str, source: &'static str)
where
    T: for<'de> FromEdn<'de>,
{
    group.throughput(Throughput::Bytes(source.len() as u64));
    group.bench_function(format!("{label}/tree"), |b| {
        b.iter(|| {
            let edn = read_string(black_box(source)).unwrap();
            T::from_edn(&edn).unwrap()
        })
    });
    group.bench_function(format!("{label}/stream"), |b| {
        b.iter(|| T::from_edn_str(black_box(source)).unwrap())
    });
}

// -- Entity DSLs --------------------

fn bench_atom_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("atom_dsl");
    bench_pair::<AtomDsl>(&mut g, "plain", r##""C""##);
    bench_pair::<AtomDsl>(&mut g, "charge", r##""N#c+""##);
    bench_pair::<AtomDsl>(&mut g, "full", r##""C#c+1#R+#v4""##);
    g.finish();
}

fn bench_bond_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("bond_dsl");
    bench_pair::<BondDsl>(&mut g, "single", r##""1""##);
    bench_pair::<BondDsl>(&mut g, "aromatic", r##""1#a""##);
    bench_pair::<BondDsl>(&mut g, "full", r##""2#c+1#u1#s2#R+""##);
    g.finish();
}

fn bench_dative_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("dative_dsl");
    bench_pair::<DativeBondDsl>(&mut g, "empty", r##""""##);
    bench_pair::<DativeBondDsl>(&mut g, "ring", r##""#R#r6""##);
    g.finish();
}

fn bench_aromatic_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("aromatic_dsl");
    bench_pair::<AromaticSystemDsl>(&mut g, "empty", r##""""##);
    bench_pair::<AromaticSystemDsl>(&mut g, "full", r##""#c0#u0#s1#e6""##);
    g.finish();
}

fn bench_multicenter_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("multicenter_dsl");
    bench_pair::<MulticenterBondDsl>(&mut g, "empty", r##""""##);
    bench_pair::<MulticenterBondDsl>(&mut g, "full", r##""#c0#u0#s1#e2""##);
    g.finish();
}

fn bench_noncovalent_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("noncovalent_dsl");
    bench_pair::<NoncovalentBondDsl>(&mut g, "literal", r##""Hbd""##);
    bench_pair::<NoncovalentBondDsl>(&mut g, "set", r##""{Hbd,Ion}""##);
    bench_pair::<NoncovalentBondDsl>(&mut g, "bind", r##""(?k :: {Hbd,Xbd})""##);
    g.finish();
}

// -- Constraint DSLs --------------------
//
// `ConstraintDsl` and `ConstraintsDsl` don't implement the `FromEdn` trait
// with tree → (ast, metadata) resolution the way `MoleculeDsl` does; they're
// still round-trippable via `from_edn`/streaming for the DSL form itself.
// Benchmarked directly via their `FromEdn` impls.

fn bench_constraint_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("constraint_dsl");
    bench_pair::<ConstraintDsl>(&mut g, "atom_leaf", r##"{:atom [0 {:valence 4}]}"##);
    bench_pair::<ConstraintDsl>(&mut g, "bond_leaf_flag", r##"{:bond [0 :aromatic]}"##);
    bench_pair::<ConstraintDsl>(
        &mut g,
        "molecule_connected",
        r##"{:connected [0 1 2 3 4 5]}"##,
    );
    bench_pair::<ConstraintDsl>(
        &mut g,
        "molecule_sub_pattern",
        r##"{:sub-pattern {:anchor {:atoms [[0 0]]} :pattern {:atoms ["N"] :bonds []}}}"##,
    );
    bench_pair::<ConstraintDsl>(
        &mut g,
        "nested_combinators",
        r##"{:and [{:or [{:atom [0 {:valence 3}]} {:atom [0 {:valence 4}]}]} {:not {:connected [0 1]}}]}"##,
    );
    g.finish();
}

fn bench_constraints_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("constraints_dsl");
    bench_pair::<ConstraintsDsl>(&mut g, "empty", r##"[]"##);
    bench_pair::<ConstraintsDsl>(
        &mut g,
        "mixed_small",
        r##"[{:atom [0 {:valence 4}]} {:connected [0 1]} {:not {:bond [0 :aromatic]}}]"##,
    );
    g.finish();
}

// -- MoleculeDsl --------------------

const MOL_SMALL: &str = r##"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"##;

const MOL_BENZENE: &str = r##"{:atoms ["C" "C" "C" "C" "C" "C"]
 :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
 :aromatic [{:atoms [0 1 2 3 4 5] :type "#e6"}]}"##;

const MOL_INDOLE: &str = r##"{:atoms [[:n "N"] [:c2 "C"] [:c3 "C"] [:c3a "C"] [:c4 "C"] [:c5 "C"] [:c6 "C"] [:c7 "C"] [:c7a "C"]]
 :bonds [[:n :c2 "1"] [:c2 :c3 "1"] [:c3 :c3a "1"] [:c3a :c4 "1"] [:c4 :c5 "1"] [:c5 :c6 "1"] [:c6 :c7 "1"] [:c7 :c7a "1"] [:c7a :n "1"] [:c3a :c7a "1"]]
 :aromatic [{:atoms [:n :c2 :c3 :c3a :c4 :c5 :c6 :c7 :c7a] :type "#e10"}]}"##;

const MOL_WITH_CONSTRAINTS: &str = r##"{:atoms [[:c1 "C"] [:c2 "C"] [:o "O"]]
 :bonds [{:id :b1 :a :c1 :b :c2 :type "1"} {:id :b2 :a :c2 :b :o :type "1"}]
 :constraints [{:connected [:c1 :c2 :o]}
               {:bond-order-sum {:bonds [:b1 :b2] :sum 2}}
               {:not {:atom [:c1 {:valence 3}]}}]}"##;

fn bench_molecule_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("molecule_dsl");
    bench_pair::<MoleculeDsl>(&mut g, "small", MOL_SMALL);
    bench_pair::<MoleculeDsl>(&mut g, "benzene", MOL_BENZENE);
    bench_pair::<MoleculeDsl>(&mut g, "indole", MOL_INDOLE);
    bench_pair::<MoleculeDsl>(&mut g, "with_constraints", MOL_WITH_CONSTRAINTS);
    g.finish();
}

criterion_group!(
    entity,
    bench_atom_dsl,
    bench_bond_dsl,
    bench_dative_dsl,
    bench_aromatic_dsl,
    bench_multicenter_dsl,
    bench_noncovalent_dsl,
);
criterion_group!(constraint, bench_constraint_dsl, bench_constraints_dsl);
criterion_group!(molecule, bench_molecule_dsl);
criterion_main!(entity, constraint, molecule);
