//! End-to-end SMILES ingest benchmarks — parse, raise, and resolve — over
//! representative single molecules: an unbranched alkane (localized valence
//! only), benzene and pyridine (joint aromatic selection), bare methane
//! (plural admission collapsed by the tie-break), and the bare-nitrogen
//! fused heteroaromatics quinoline and purine (assignment search over a
//! flexible fused component).
//!
//! The `smiles_roundtrip/resolve` group isolates resolution with the SMILES counts
//! model. Parsing, raising, resolver construction, and input cloning are outside timing.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use umol_graph::ingest::ingest_smiles;
use umol_graph::ops::model::{ChemistryModel, ValenceModel};
use umol_graph::ops::resolve::{IsotopePolicy, ResolveConfig, Resolver};
use umol_graph_ir::ir::{Molecule, TryIntoIr};
use umol_io::smiles::Smiles;
use umol_utils::solution::Solution;

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

fn bench_resolve(c: &mut Criterion) {
    let model = ChemistryModel {
        valence: ValenceModel::smiles(),
        ..ChemistryModel::default()
    };
    let resolver = Resolver::with_config(
        &model,
        ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        },
    );
    let mut group = c.benchmark_group("smiles_roundtrip/resolve");
    for (name, input) in [
        (
            "chain_64",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        ),
        ("branched", "CC(C)(F)C(Cl)CC"),
        ("components", "CC.O.[Na+]"),
        ("aromatic_fused", "c1ccc2ccccc2c1"),
        ("aromatic_lone_pair", "[nH]1cccc1"),
        ("charged", "[NH4+]"),
        ("radical", "[CH3]"),
        ("tetra_four", "[C@](F)(Cl)(Br)I"),
        ("tetra_ring", "C[C@H]1CCCCO1"),
        ("tetra_explicit_h", "C[C@]1([H])CCCCO1"),
        ("tetra_lone_pair", "C[S@](=O)CC"),
        ("alkene_four", "F/C(Cl)=C(/Br)I"),
        ("shared_chain", "C/C=C/C=C/C"),
        ("partial_triene", "C/C=C/C=CC=C/C"),
        ("shared_branch", "F/C=C(/C=C/F)C=C"),
        ("shared_cycle", "C1/C=C/C=C/CCC1"),
    ] {
        let table = Smiles::parse(input).unwrap().into_table_ir();
        let molecule: Molecule = (&table).try_into_ir(&()).unwrap();
        let mut checked = molecule.clone();
        assert!(
            matches!(
                resolver.resolve(&mut checked).unwrap(),
                Solution::Determined(_)
            ),
            "{name}"
        );
        group.bench_function(name, |b| {
            b.iter_batched_ref(
                || molecule.clone(),
                |molecule| resolver.resolve(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_valence_project(c: &mut Criterion) {
    for (label, valence) in [
        ("counts", ValenceModel::smiles()),
        ("atom_typing", ValenceModel::default()),
    ] {
        let model = ChemistryModel {
            valence,
            ..Default::default()
        };
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                isotope: IsotopePolicy::Natural,
                ..Default::default()
            },
        );
        let mut group = c.benchmark_group(format!("smiles_roundtrip/valence_project/{label}"));
        for (name, input) in [
            ("methane", "C"),
            (
                "chain_64",
                "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            ),
            ("methyl", "[CH3]"),
            ("ammonium", "[NH4+]"),
        ] {
            let mut source = ingest_smiles(input).unwrap();
            assert!(matches!(
                resolver.resolve(&mut source).unwrap(),
                Solution::Determined(_)
            ));
            let mut checked = source.clone();
            assert!(matches!(
                resolver
                    .valence
                    .project(&mut checked, resolver.tie_break)
                    .unwrap(),
                Solution::Determined(_)
            ));
            assert!(matches!(
                resolver.resolve(&mut checked).unwrap(),
                Solution::Determined(_)
            ));
            assert_eq!(checked, source);
            group.bench_function(name, |b| {
                b.iter_batched_ref(
                    || source.clone(),
                    |molecule| {
                        resolver
                            .valence
                            .project(black_box(molecule), resolver.tie_break)
                            .unwrap()
                    },
                    BatchSize::SmallInput,
                );
            });
            group.bench_function(format!("{name}_and_resolve"), |b| {
                b.iter_batched_ref(
                    || source.clone(),
                    |molecule| {
                        let projected = resolver
                            .valence
                            .project(black_box(molecule), resolver.tie_break)
                            .unwrap();
                        let resolved = resolver.resolve(molecule).unwrap();
                        black_box((projected, resolved))
                    },
                    BatchSize::SmallInput,
                );
            });
        }
        group.finish();
    }
}

criterion_group!(
    resolve,
    bench_ingest_smiles,
    bench_resolve,
    bench_valence_project
);
criterion_main!(resolve);
