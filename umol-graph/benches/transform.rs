//! Molecular transformations through consuming and in-place operations.
//! Parsing, construction, and cloning are outside timing.
//! Inputs share storage with their fixtures unless labeled unique; shared mutation includes
//! storage detachment.

#[path = "transform/stereo.rs"]
mod stereo;

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph::ops::model::AromaticityModel;
use umol_graph::ops::transform::{
    Aromatizer, HydrogenFolder, HydrogenUnfolder, KekulizeConfig, Kekulizer, Transformer,
};
use umol_graph::ops::validate::ValenceInvariantsValidator;
use umol_graph_ir::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_graph_ir::ir::{AtomForm, AtomId, BondForm, IntoIr, Molecule, MoleculeEntries};
use umol_graph_ir::mol_dsl_concrete;
use umol_utils::solution::Solution;

fn bench_transform(c: &mut Criterion) {
    let aromatizer = Aromatizer::new(&AromaticityModel::daylight());
    let defaults = MoleculeDefaults::concrete();
    let mut group = c.benchmark_group("transform");
    for (name, aromatic_source, kekule_source) in [
        (
            "benzene",
            include_str!("../tests/kekulization/data/benzene_aromatic_input.edn"),
            include_str!("../tests/kekulization/data/benzene_kekulized_expected.edn"),
        ),
        (
            "pyridine",
            include_str!("../tests/kekulization/data/pyridine_aromatic_input.edn"),
            include_str!("../tests/kekulization/data/pyridine_kekulized_expected.edn"),
        ),
        (
            "pyrrole",
            include_str!("../tests/kekulization/data/pyrrole_aromatic_input.edn"),
            include_str!("../tests/kekulization/data/pyrrole_kekulized_expected.edn"),
        ),
        (
            "furan",
            include_str!("../tests/kekulization/data/furan_aromatic_input.edn"),
            include_str!("../tests/kekulization/data/furan_kekulized_expected.edn"),
        ),
        (
            "thiophene",
            include_str!("../tests/kekulization/data/thiophene_aromatic_input.edn"),
            include_str!("../tests/kekulization/data/thiophene_kekulized_expected.edn"),
        ),
        (
            "azulene",
            include_str!("../tests/kekulization/data/azulene_aromatic_input.edn"),
            include_str!("../tests/kekulization/data/azulene_kekulized_expected.edn"),
        ),
        (
            "indole_prescribed_donor",
            include_str!("../tests/kekulization/data/indole_prescribed_donor_aromatic_input.edn"),
            include_str!(
                "../tests/kekulization/data/indole_prescribed_donor_kekulized_expected.edn"
            ),
        ),
    ] {
        let aromatic = aromatic_source
            .parse::<MoleculeDsl>()
            .unwrap()
            .into_ir(&defaults);
        let kekule = kekule_source
            .parse::<MoleculeDsl>()
            .unwrap()
            .into_ir(&defaults);
        assert!(aromatic.has_aromatic_systems(), "{name}");
        assert!(!kekule.has_aromatic_systems(), "{name}");
        assert_eq!(
            ValenceInvariantsValidator.validate(&aromatic),
            Ok(Solution::Determined(())),
            "{name} aromatic input"
        );
        assert_eq!(
            ValenceInvariantsValidator.validate(&kekule),
            Ok(Solution::Determined(())),
            "{name} localized output"
        );
        let kekulizer = Kekulizer::new(KekulizeConfig::default(), aromatic.atoms().ids().collect());
        assert_eq!(
            kekulizer.transform(aromatic.clone()).unwrap(),
            kekule,
            "{name}"
        );
        let mut localized = aromatic.clone();
        kekulizer.transform_into(&mut localized).unwrap();
        assert_eq!(localized, kekule, "{name}");

        let perceived = aromatizer.transform(kekule.clone()).unwrap();
        assert_eq!(
            ValenceInvariantsValidator.validate(&perceived),
            Ok(Solution::Determined(())),
            "{name} aromatized output"
        );
        assert_eq!(
            perceived
                .aromatic_systems()
                .iter()
                .map(|system| (system.atom_ids().collect::<Vec<_>>(), system.attributes()))
                .collect::<Vec<_>>(),
            aromatic
                .aromatic_systems()
                .iter()
                .map(|system| (system.atom_ids().collect::<Vec<_>>(), system.attributes()))
                .collect::<Vec<_>>(),
            "{name}"
        );
        let mut in_place = kekule.clone();
        aromatizer.transform_into(&mut in_place).unwrap();
        assert_eq!(in_place, perceived, "{name}");

        group.bench_function(BenchmarkId::new("kekulize", name), |b| {
            b.iter_batched(
                || aromatic.clone(),
                |molecule| kekulizer.transform(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("kekulize_into", name), |b| {
            b.iter_batched_ref(
                || aromatic.clone(),
                |molecule| kekulizer.transform_into(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("aromatize", name), |b| {
            b.iter_batched(
                || kekule.clone(),
                |molecule| aromatizer.transform(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("aromatize_into", name), |b| {
            b.iter_batched_ref(
                || kekule.clone(),
                |molecule| aromatizer.transform_into(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_fold_hydrogens(c: &mut Criterion) {
    let mut group = c.benchmark_group("fold_hydrogens");
    for (name, input, expected) in [
        (
            "methane",
            mol_dsl_concrete!(r#"{:atoms ["C#h3" "H"] :bonds [[0 1 "1"]]}"#),
            mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#),
        ),
        (
            "methanol",
            mol_dsl_concrete!(
                r#"{
                :atoms ["H" "C" "H" "H" "O" "H"]
                :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [1 4 "1"] [4 5 "1"]]
            }"#
            ),
            mol_dsl_concrete!(r#"{:atoms ["C#h3" "O#h"] :bonds [[0 1 "1"]]}"#),
        ),
        (
            "tetrahedral",
            mol_dsl_concrete!(
                r#"{
                :atoms ["H" "C#T0" "F" "Cl" "Br"]
                :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [1 4 "1"]]
                :stereo-atoms [{:site 1 :ligands [0 2 3 4] :attrs "Th0"}]
            }"#
            ),
            mol_dsl_concrete!(
                r#"{
                :atoms ["C#h#T1" "F" "Cl" "Br"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
                :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th0"}]
            }"#
            ),
        ),
        (
            "unchanged",
            mol_dsl_concrete!(r#"{:atoms ["C#h3" "O#h"] :bonds [[0 1 "1"]]}"#),
            mol_dsl_concrete!(r#"{:atoms ["C#h3" "O#h"] :bonds [[0 1 "1"]]}"#),
        ),
    ] {
        assert_eq!(
            HydrogenFolder.transform(input.clone()).unwrap(),
            expected,
            "{name}"
        );
        let mut in_place = input.clone();
        HydrogenFolder.transform_into(&mut in_place).unwrap();
        assert_eq!(in_place, expected, "{name}");
        group.bench_function(BenchmarkId::new("transform", name), |b| {
            b.iter_batched(
                || input.clone(),
                |molecule| HydrogenFolder.transform(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_into", name), |b| {
            b.iter_batched_ref(
                || input.clone(),
                |molecule| HydrogenFolder.transform_into(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        if name == "methane" {
            group.bench_function(BenchmarkId::new("transform_unique", name), |b| {
                b.iter_batched(
                    || mol_dsl_concrete!(r#"{:atoms ["C#h3" "H"] :bonds [[0 1 "1"]]}"#),
                    |molecule| HydrogenFolder.transform(black_box(molecule)).unwrap(),
                    BatchSize::SmallInput,
                );
            });
            group.bench_function(BenchmarkId::new("transform_into_unique", name), |b| {
                b.iter_batched_ref(
                    || mol_dsl_concrete!(r#"{:atoms ["C#h3" "H"] :bonds [[0 1 "1"]]}"#),
                    |molecule| HydrogenFolder.transform_into(black_box(molecule)).unwrap(),
                    BatchSize::SmallInput,
                );
            });
        }
    }
    group.finish();
}

fn bench_unfold_hydrogens(c: &mut Criterion) {
    let mut group = c.benchmark_group("unfold_hydrogens");
    for (name, input, expected) in [
        (
            "methane",
            mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#),
            mol_dsl_concrete!(
                r#"{
                :atoms ["C" "H" "H" "H" "H"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            }"#
            ),
        ),
        (
            "methanol",
            mol_dsl_concrete!(r#"{:atoms ["C#h3" "O#h"] :bonds [[0 1 "1"]]}"#),
            mol_dsl_concrete!(
                r#"{
                :atoms ["C" "O" "H" "H" "H" "H"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"]]
            }"#
            ),
        ),
        (
            "tetrahedral",
            mol_dsl_concrete!(
                r#"{
                :atoms ["C#h#T1" "F" "Cl" "Br"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
                :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th0"}]
            }"#
            ),
            mol_dsl_concrete!(
                r#"{
                :atoms ["C#T1" "F" "Cl" "Br" "H"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [4 1 2 3] :attrs "Th0"}]
            }"#
            ),
        ),
        (
            "unchanged",
            mol_dsl_concrete!(r#"{:atoms ["C" "H"] :bonds [[0 1 "1"]]}"#),
            mol_dsl_concrete!(r#"{:atoms ["C" "H"] :bonds [[0 1 "1"]]}"#),
        ),
    ] {
        assert_eq!(
            HydrogenUnfolder.transform(input.clone()).unwrap(),
            expected,
            "{name}"
        );
        let mut in_place = input.clone();
        HydrogenUnfolder.transform_into(&mut in_place).unwrap();
        assert_eq!(in_place, expected, "{name}");
        group.bench_function(BenchmarkId::new("transform", name), |b| {
            b.iter_batched(
                || input.clone(),
                |molecule| HydrogenUnfolder.transform(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_into", name), |b| {
            b.iter_batched_ref(
                || input.clone(),
                |molecule| {
                    HydrogenUnfolder
                        .transform_into(black_box(molecule))
                        .unwrap()
                },
                BatchSize::SmallInput,
            );
        });
        if name == "methane" {
            group.bench_function(BenchmarkId::new("transform_unique", name), |b| {
                b.iter_batched(
                    || mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#),
                    |molecule| HydrogenUnfolder.transform(black_box(molecule)).unwrap(),
                    BatchSize::SmallInput,
                );
            });
            group.bench_function(BenchmarkId::new("transform_into_unique", name), |b| {
                b.iter_batched_ref(
                    || mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#),
                    |molecule| {
                        HydrogenUnfolder
                            .transform_into(black_box(molecule))
                            .unwrap()
                    },
                    BatchSize::SmallInput,
                );
            });
        }
    }
    group.finish();
}

fn bench_hydrogens_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("hydrogens_scaling");
    for size in [16usize, 64, 256, 1024] {
        let mut implicit = MoleculeEntries::default();
        let mut explicit = MoleculeEntries::default();
        for i in 0..size {
            let count = if i == 0 || i + 1 == size { 3 } else { 2 };
            implicit.atoms.push(
                AtomForm::from_element(Element::C)
                    .with_implicit_hydrogens(count)
                    .into_concrete(),
            );
            explicit
                .atoms
                .push(AtomForm::from_element(Element::C).into_concrete());
            if i != 0 {
                let bond = (
                    AtomId::from(i - 1),
                    AtomId::from(i),
                    BondForm::from_order(1).into_concrete(),
                );
                implicit.bonds.push(bond.clone());
                explicit.bonds.push(bond);
            }
        }
        for i in 0..size {
            let count = if i == 0 || i + 1 == size { 3 } else { 2 };
            for _ in 0..count {
                let hydrogen = AtomId::from(explicit.atoms.len());
                explicit
                    .atoms
                    .push(AtomForm::from_element(Element::H).into_concrete());
                explicit.bonds.push((
                    AtomId::from(i),
                    hydrogen,
                    BondForm::from_order(1).into_concrete(),
                ));
            }
        }
        let implicit = Molecule::from_entries(implicit);
        let explicit = Molecule::from_entries(explicit);
        assert_eq!(
            HydrogenFolder.transform(explicit.clone()).unwrap(),
            implicit
        );
        assert_eq!(
            HydrogenUnfolder.transform(implicit.clone()).unwrap(),
            explicit
        );
        let mut folded = explicit.clone();
        HydrogenFolder.transform_into(&mut folded).unwrap();
        assert_eq!(folded, implicit);
        let mut unfolded = implicit.clone();
        HydrogenUnfolder.transform_into(&mut unfolded).unwrap();
        assert_eq!(unfolded, explicit);

        group.bench_function(BenchmarkId::new("fold", size), |b| {
            b.iter_batched(
                || explicit.clone(),
                |molecule| HydrogenFolder.transform(black_box(molecule)).unwrap(),
                BatchSize::LargeInput,
            );
        });
        group.bench_function(BenchmarkId::new("fold_into", size), |b| {
            b.iter_batched_ref(
                || explicit.clone(),
                |molecule| HydrogenFolder.transform_into(black_box(molecule)).unwrap(),
                BatchSize::LargeInput,
            );
        });
        group.bench_function(BenchmarkId::new("unfold", size), |b| {
            b.iter_batched(
                || implicit.clone(),
                |molecule| HydrogenUnfolder.transform(black_box(molecule)).unwrap(),
                BatchSize::LargeInput,
            );
        });
        group.bench_function(BenchmarkId::new("unfold_into", size), |b| {
            b.iter_batched_ref(
                || implicit.clone(),
                |molecule| {
                    HydrogenUnfolder
                        .transform_into(black_box(molecule))
                        .unwrap()
                },
                BatchSize::LargeInput,
            );
        });
    }
    group.finish();
}

criterion_group!(
    transform,
    bench_transform,
    bench_fold_hydrogens,
    bench_unfold_hydrogens,
    bench_hydrogens_scaling,
    stereo::bench_stereo
);
criterion_main!(transform);
