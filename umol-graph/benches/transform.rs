//! Kekulization and aromatization of paired fixtures, through consuming and
//! in-place operations. Parsing, construction, and cloning are outside timing.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_graph::ops::model::AromaticityModel;
use umol_graph::ops::transform::{Aromatizer, KekulizeConfig, Kekulizer, Transformer};
use umol_graph::ops::validate::ValenceInvariantsValidator;
use umol_graph_ir::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_graph_ir::ir::IntoIr;
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

criterion_group!(transform, bench_transform);
criterion_main!(transform);
