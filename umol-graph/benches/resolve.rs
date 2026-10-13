//! End-to-end SMILES ingest benchmarks — parse, raise, and resolve — over
//! representative single molecules: an unbranched alkane (localized valence
//! only), benzene and pyridine (joint aromatic selection), bare methane
//! (plural admission collapsed by the tie-break), and the bare-nitrogen
//! fused heteroaromatics quinoline and purine (assignment search over a
//! flexible fused component).
//!
//! The `smiles_roundtrip/resolve` group isolates resolution with the SMILES counts
//! model. Parsing, raising, resolver construction, and input cloning are outside timing.
//! The `valence_readings` group measures counts/registry admission, registry
//! classification of resolved molecules, and complete resolution with atom typing.
//! `aromatic_selection/hueckel` stresses joint selection with 128–4096 admitted
//! assignments. `aromaticity_resolve/hueckel` isolates standalone resolution of
//! projected aromatic molecules. Setup and input cloning are outside timing.

use std::collections::BTreeSet;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph::ingest::ingest_smiles;
use umol_graph::ops::model::{AromaticityModel, ChemistryModel, ValenceModel, ValenceTieBreak};
use umol_graph::ops::resolve::{
    AromaticityResolver, IsotopePolicy, ProjectFlags, ResolveConfig, Resolver,
};
use umol_graph::ops::valence::{AtomTypeRegistry, AtomTypingValence, CountsValence, ValenceTable};
use umol_graph_ir::ir::{AsLit, ElementForm, Molecule, Reframe, TryIntoIr, UnpairedElectronsForm};
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
                resolver.resolve_into(&mut checked).unwrap(),
                Solution::Determined(_)
            ),
            "{name}"
        );
        group.bench_function(name, |b| {
            b.iter_batched_ref(
                || molecule.clone(),
                |molecule| resolver.resolve_into(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_valence(c: &mut Criterion) {
    let typing = AtomTypingValence::new(AtomTypeRegistry::default_registry());
    let counts = CountsValence::new(ValenceTable::smiles_table());
    let mut model = ChemistryModel::default();
    model.valence.tie_break = ValenceTieBreak::MostSaturated;
    let resolver = Resolver::with_config(
        &model,
        ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        },
    );
    let mut group = c.benchmark_group("valence_readings");
    for (name, input) in [
        ("octane", "CCCCCCCC"),
        (
            "chain_64",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        ),
        ("benzene", "c1ccccc1"),
        ("pyridine", "c1ccncc1"),
        ("charged", "[NH4+]"),
        ("radical", "[CH3]"),
    ] {
        let table = Smiles::parse(input).unwrap().into_table_ir();
        let molecule: Molecule = (&table).try_into_ir(&()).unwrap();
        let ground = ingest_smiles(input).unwrap();
        assert!(matches!(typing.admit(&molecule), Solution::Determined(_)));
        assert!(matches!(counts.admit(&molecule), Solution::Determined(_)));
        for id in ground.atoms().ids() {
            assert_eq!(
                typing.classify_molecule_atom(&ground, id),
                Solution::Determined(())
            );
        }
        let mut checked = molecule.clone();
        assert!(matches!(
            resolver.resolve_into(&mut checked).unwrap(),
            Solution::Determined(_)
        ));
        group.bench_function(BenchmarkId::new("counts_admit", name), |b| {
            b.iter(|| counts.admit(black_box(&molecule)));
        });
        group.bench_function(BenchmarkId::new("typing_admit", name), |b| {
            b.iter(|| typing.admit(black_box(&molecule)));
        });
        group.bench_function(BenchmarkId::new("typing_classify", name), |b| {
            b.iter(|| {
                let molecule = black_box(&ground);
                for id in molecule.atoms().ids() {
                    black_box(typing.classify_molecule_atom(molecule, id));
                }
            });
        });
        group.bench_function(BenchmarkId::new("typing_resolve", name), |b| {
            b.iter_batched_ref(
                || molecule.clone(),
                |molecule| resolver.resolve_into(black_box(molecule)).unwrap(),
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
            resolver.aromaticity.project_into(&mut checked),
            Ok(Solution::Determined(()))
        );
        assert!(!checked.has_aromatic_systems());
        group.bench_function(name, |b| {
            b.iter_batched_ref(
                || source.clone(),
                |molecule| {
                    resolver
                        .aromaticity
                        .project_into(black_box(molecule))
                        .unwrap()
                },
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_aromatic_selection(c: &mut Criterion) {
    let mut model = ChemistryModel::default();
    model.valence.tie_break = ValenceTieBreak::MostSaturated;
    let resolver = Resolver::with_config(
        &model,
        ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        },
    );
    let mut group = c.benchmark_group("aromatic_selection/hueckel");
    for (name, smiles, nitrogen_count, expected_combinations) in [
        ("heptaazanaphthalene", "c1ccc2ccccc2c1", 7, 128),
        ("octaazanaphthalene", "c1ccc2ccccc2c1", 8, 256),
        ("decaazaanthracene", "c1ccc2cc3ccccc3cc2c1", 10, 1024),
        ("dodecaazatetracene", "c1ccc2cc3cc4ccccc4cc3cc2c1", 12, 4096),
    ] {
        let table = Smiles::parse(smiles).unwrap().into_table_ir();
        let mut molecule: Molecule = (&table).try_into_ir(&()).unwrap();
        let nitrogens: Vec<_> = molecule
            .atoms()
            .iter()
            .filter(|atom| atom.neighbors().count() == 2)
            .map(|atom| atom.id())
            .take(nitrogen_count)
            .collect();
        assert_eq!(nitrogens.len(), nitrogen_count, "{name}");
        for atom in molecule.atoms().ids().collect::<Vec<_>>() {
            molecule.atom_mut(atom).attributes_mut().unpaired_electrons =
                UnpairedElectronsForm::closed_shell();
        }
        for atom in nitrogens {
            molecule.atom_mut(atom).attributes_mut().element = ElementForm::Lit(Element::N);
        }
        let Solution::Determined(state) = resolver.valence.admit(&molecule).unwrap() else {
            panic!("{name}: valence admission must be determined");
        };
        let mut flexible_atoms = 0;
        let mut combinations = 1;
        for (_, forms) in state.completions.iter() {
            let contributions: Option<BTreeSet<_>> = forms
                .iter()
                .map(|form| form.constraints.aromatic_valence().as_lit())
                .collect();
            let contributions =
                contributions.expect("fixture completions have literal contributions");
            if contributions.len() > 1 {
                flexible_atoms += 1;
                combinations *= forms.len();
            }
        }
        assert_eq!(flexible_atoms, nitrogen_count, "{name}");
        assert_eq!(combinations, expected_combinations, "{name}");
        assert!(
            matches!(
                resolver.resolve(molecule.clone()).unwrap(),
                Solution::Determined(_)
            ),
            "{name}"
        );
        group.bench_function(BenchmarkId::new(name, combinations), |b| {
            b.iter_batched(
                || molecule.clone(),
                |molecule| resolver.resolve(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_aromaticity_resolve(c: &mut Criterion) {
    let resolver = AromaticityResolver::new(&AromaticityModel::daylight());
    let mut group = c.benchmark_group("aromaticity_resolve/hueckel");
    for (name, smiles) in [
        ("benzene", "c1ccccc1"),
        ("pyridine", "c1ccncc1"),
        ("pyrrole", "[nH]1cccc1"),
        ("furan", "o1cccc1"),
        ("thiophene", "s1cccc1"),
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("quinoline", "c1ccc2ccccc2n1"),
        ("purine", "c1ncc2ncnc2n1"),
        ("biphenyl", "c1ccccc1-c1ccccc1"),
        ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
    ] {
        let expected = ingest_smiles(smiles).unwrap();
        assert!(expected.has_aromatic_systems(), "{name}");
        let Solution::Determined(molecule) = resolver.project(expected.clone()).unwrap() else {
            panic!("{name}: aromaticity projection must be determined");
        };
        assert!(!molecule.has_aromatic_systems(), "{name}");
        let Solution::Determined(actual) = resolver.resolve(molecule.clone()).unwrap() else {
            panic!("{name}: standalone aromaticity resolution must be determined");
        };
        assert!(actual.framed_eq(&expected), "{name}");
        group.bench_function(name, |b| {
            b.iter_batched(
                || molecule.clone(),
                |molecule| resolver.resolve(black_box(molecule)).unwrap(),
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
            resolver.stereo.project_into(&mut checked),
            Ok(Solution::Determined(()))
        );
        assert!(!checked.has_stereo_atoms() && !checked.has_stereo_bonds());
        group.bench_function(name, |b| {
            b.iter_batched_ref(
                || source.clone(),
                |molecule| resolver.stereo.project_into(black_box(molecule)).unwrap(),
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
                    resolver.project_into(&mut checked, ProjectFlags::all()),
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
                                    .project_into(black_box(molecule), ProjectFlags::all())
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
            resolver.resolve_into(&mut checked).unwrap(),
            Solution::Determined(_)
        ));
        let mut checked = projected.clone();
        assert_eq!(
            resolver.project_into(&mut checked, ProjectFlags::all()),
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
                        |molecule| black_box(resolver.resolve_into(molecule).unwrap()),
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
                            black_box(
                                resolver
                                    .project_into(molecule, ProjectFlags::all())
                                    .unwrap(),
                            )
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
    bench_valence,
    bench_aromaticity_project,
    bench_aromatic_selection,
    bench_aromaticity_resolve,
    bench_stereo_project,
    bench_project,
    bench_resolve_project_ownership
);
criterion_main!(resolve);
