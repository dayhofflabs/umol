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
use umol_ast::dsl::{
    AromaticSystemDsl, AtomDsl, BondDsl, ConstraintDsl, ConstraintsDsl, DativeBondDsl, MoleculeDsl,
    MulticenterBondDsl, NoncovalentBondDsl,
};
use umol_edn::{read_string, FromEdn};

#[path = "fixtures.rs"]
mod fixtures;
use fixtures::{
    MOL_BENZENE, MOL_DIBORANE, MOL_INDOLE, MOL_LARGE_ALL_IDS, MOL_LARGE_NO_IDS,
    MOL_LARGE_PARTIAL_IDS, MOL_SMALL, MOL_WITH_CONSTRAINTS,
};

fn bench_pair<T>(
    group: &mut criterion::BenchmarkGroup<'_, WallTime>,
    label: &str,
    source: &'static str,
) where
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
    bench_pair::<AtomDsl>(&mut g, "full", r##""C#c+1#h3#u0#s1#R+#v4""##);
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
    bench_pair::<DativeBondDsl>(&mut g, "undetermined", r##""*""##);
    bench_pair::<DativeBondDsl>(&mut g, "ring", r##""1#R(6)""##);
    g.finish();
}

fn bench_aromatic_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("aromatic_dsl");
    bench_pair::<AromaticSystemDsl>(&mut g, "undetermined", r##""*""##);
    bench_pair::<AromaticSystemDsl>(&mut g, "full", r##""[1,1,1]#c0#u0#s1#e6""##);
    g.finish();
}

fn bench_multicenter_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("multicenter_dsl");
    bench_pair::<MulticenterBondDsl>(&mut g, "undetermined", r##""*""##);
    bench_pair::<MulticenterBondDsl>(&mut g, "full", r##""[1,1,1]#c0#u0#s1#e2""##);
    g.finish();
}

fn bench_noncovalent_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("noncovalent_dsl");
    bench_pair::<NoncovalentBondDsl>(&mut g, "literal", r##""Hbd""##);
    bench_pair::<NoncovalentBondDsl>(&mut g, "undetermined", r##""*""##);
    bench_pair::<NoncovalentBondDsl>(&mut g, "intramolecular", r##""Hbd#I""##);
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
    bench_pair::<ConstraintDsl>(
        &mut g,
        "bond_leaf_flag",
        r##"{:bond [0 {:aromatic true}]}"##,
    );
    bench_pair::<ConstraintDsl>(
        &mut g,
        "molecule_connected",
        r##"{:connected {:atoms [0 1 2 3 4 5]}}"##,
    );
    bench_pair::<ConstraintDsl>(
        &mut g,
        "molecule_sub_pattern",
        r##"{:sub-pattern {:anchor {:atoms [[0 0]]} :pattern {:atoms ["N"] :bonds []}}}"##,
    );
    bench_pair::<ConstraintDsl>(
        &mut g,
        "nested_combinators",
        r##"{:and [{:or [{:atom [0 {:valence 3}]} {:atom [0 {:valence 4}]}]} {:not {:connected {:atoms [0 1]}}}]}"##,
    );
    g.finish();
}

fn bench_constraints_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("constraints_dsl");
    bench_pair::<ConstraintsDsl>(&mut g, "empty", r##"[]"##);
    bench_pair::<ConstraintsDsl>(
        &mut g,
        "mixed_small",
        r##"[{:atom [0 {:valence 4}]} {:connected {:atoms [0 1]}} {:not {:bond [0 {:aromatic true}]}}]"##,
    );
    g.finish();
}

// -- MoleculeDsl --------------------

fn bench_molecule_dsl(c: &mut Criterion) {
    let mut g = c.benchmark_group("molecule_dsl");
    bench_pair::<MoleculeDsl>(&mut g, "small", MOL_SMALL);
    bench_pair::<MoleculeDsl>(&mut g, "benzene", MOL_BENZENE);
    bench_pair::<MoleculeDsl>(&mut g, "indole", MOL_INDOLE);
    bench_pair::<MoleculeDsl>(&mut g, "diborane", MOL_DIBORANE);
    bench_pair::<MoleculeDsl>(&mut g, "with_constraints", MOL_WITH_CONSTRAINTS);
    bench_pair::<MoleculeDsl>(&mut g, "large_no_ids", MOL_LARGE_NO_IDS.as_str());
    bench_pair::<MoleculeDsl>(&mut g, "large_all_ids", MOL_LARGE_ALL_IDS.as_str());
    bench_pair::<MoleculeDsl>(&mut g, "large_partial_ids", MOL_LARGE_PARTIAL_IDS.as_str());
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
