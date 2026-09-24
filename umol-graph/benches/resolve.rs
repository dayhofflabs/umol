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

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_graph::ingest::ingest_smiles;
use umol_graph::ops::model::{ChemistryModel, ValenceModel, ValenceTieBreak};
use umol_graph::ops::resolve::{IsotopePolicy, ProjectFlags, ResolveConfig, Resolver};
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

fn bench_aromaticity_project(c: &mut Criterion) {
    let model = ChemistryModel::default();
    let resolver = Resolver::new(&model);
    let mut group = c.benchmark_group("smiles_roundtrip/aromaticity_project");
    for (name, input) in [
        ("benzene", "c1ccccc1"),
        ("pyrrole", "[nH]1cccc1"),
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("biphenyl", "c1ccccc1-c1ccccc1"),
    ] {
        let source = ingest_smiles(input).unwrap();
        let mut checked = source.clone();
        assert_eq!(
            resolver.aromaticity.project(&mut checked),
            Ok(Solution::Determined(()))
        );
        assert!(!checked.has_aromatic_systems());
        group.bench_function(name, |b| {
            b.iter_batched_ref(
                || source.clone(),
                |molecule| resolver.aromaticity.project(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_stereo_project(c: &mut Criterion) {
    let model = ChemistryModel::default();
    let resolver = Resolver::new(&model);
    let mut group = c.benchmark_group("smiles_roundtrip/stereo_project");
    for (name, input) in [
        ("tetrahedral", "N[C@H](F)Cl"),
        ("explicit_hydrogen", "[H][C@](F)(Cl)Br"),
        ("lone_pair", "C[S@](=O)CC"),
        ("cis_trans", "C/C=C/C"),
        ("mixed", "N[C@H](F)/C=C/C"),
    ] {
        let source = ingest_smiles(input).unwrap();
        assert!(source.has_stereo_atoms() || source.has_stereo_bonds());
        let mut checked = source.clone();
        assert_eq!(
            resolver.stereo.project(&mut checked),
            Ok(Solution::Determined(()))
        );
        assert!(!checked.has_stereo_atoms() && !checked.has_stereo_bonds());
        group.bench_function(name, |b| {
            b.iter_batched_ref(
                || source.clone(),
                |molecule| resolver.stereo.project(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_project(c: &mut Criterion) {
    let mut group = c.benchmark_group("smiles_roundtrip/project");
    for (source_name, valence) in [
        ("typing", ValenceModel::default()),
        ("counts", ValenceModel::smiles()),
    ] {
        for (policy_name, policy) in [
            ("strict", ValenceTieBreak::Strict),
            ("saturated", ValenceTieBreak::MostSaturated),
        ] {
            let mut model = ChemistryModel {
                valence: valence.clone(),
                ..Default::default()
            };
            model.valence.tie_break = policy;
            let resolver = Resolver::with_config(
                &model,
                ResolveConfig {
                    isotope: IsotopePolicy::Natural,
                    ..Default::default()
                },
            );
            for (name, input) in [
                ("octane", "CCCCCCCC"),
                (
                    "chain_64",
                    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
                ),
                ("benzene", "c1ccccc1"),
                ("naphthalene", "c1ccc2ccccc2c1"),
                ("stereo", "N[C@H](F)/C=C/C"),
                ("combined", "[13CH3][C@H](F)/C=C/c1ccccc1"),
                ("charged", "[NH4+]"),
                ("radical", "[CH3]"),
            ] {
                let source = ingest_smiles(input).unwrap();
                let mut checked = source.clone();
                assert_eq!(
                    resolver.project(&mut checked, ProjectFlags::all()),
                    Ok(Solution::Determined(()))
                );
                assert!(
                    !checked.has_stereo_atoms()
                        && !checked.has_stereo_bonds()
                        && !checked.has_aromatic_systems()
                );
                for (actual, original) in checked.atoms().iter().zip(source.atoms().iter()) {
                    assert_eq!(actual.implicit_hydrogens(), original.implicit_hydrogens());
                }
                group.bench_function(
                    BenchmarkId::new(format!("{source_name}_{policy_name}"), name),
                    |b| {
                        b.iter_batched_ref(
                            || source.clone(),
                            |molecule| {
                                resolver
                                    .project(black_box(molecule), ProjectFlags::all())
                                    .unwrap()
                            },
                            BatchSize::SmallInput,
                        );
                    },
                );
            }
        }
    }
    group.finish();
}

fn bench_resolve_project_ownership(c: &mut Criterion) {
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
    let dense = ["[13CH3][C@H](F)/C=C/c1ccccc1"; 8].join(".");
    let mut group = c.benchmark_group("resolve_project_ownership");
    for (name, smiles) in [("sparse8", "CCCCCCCC"), ("dense88", dense.as_str())] {
        let table = Smiles::parse(smiles).unwrap().into_table_ir();
        let raw: Molecule = (&table).try_into_ir(&()).unwrap();
        let projected = ingest_smiles(smiles).unwrap();
        let mut checked = raw.clone();
        assert!(matches!(
            resolver.resolve(&mut checked).unwrap(),
            Solution::Determined(_)
        ));
        let mut checked = projected.clone();
        assert_eq!(
            resolver.project(&mut checked, ProjectFlags::all()),
            Ok(Solution::Determined(()))
        );
        for (ownership, unique) in [("unique", true), ("shared", false)] {
            group.bench_function(
                BenchmarkId::new(format!("{name}/{ownership}"), "resolve"),
                |b| {
                    b.iter_batched_ref(
                        || {
                            if unique {
                                let table = Smiles::parse(smiles).unwrap().into_table_ir();
                                (&table).try_into_ir(&()).unwrap()
                            } else {
                                raw.clone()
                            }
                        },
                        |molecule| black_box(resolver.resolve(molecule).unwrap()),
                        BatchSize::SmallInput,
                    );
                },
            );
            group.bench_function(
                BenchmarkId::new(format!("{name}/{ownership}"), "project"),
                |b| {
                    b.iter_batched_ref(
                        || {
                            if unique {
                                ingest_smiles(smiles).unwrap()
                            } else {
                                projected.clone()
                            }
                        },
                        |molecule| {
                            black_box(resolver.project(molecule, ProjectFlags::all()).unwrap())
                        },
                        BatchSize::SmallInput,
                    );
                },
            );
        }
    }
    group.finish();
}

criterion_group!(
    resolve,
    bench_ingest_smiles,
    bench_resolve,
    bench_aromaticity_project,
    bench_stereo_project,
    bench_project,
    bench_resolve_project_ownership
);
criterion_main!(resolve);
