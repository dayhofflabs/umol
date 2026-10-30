//! Molecular transformations through consuming and in-place operations.
//! Parsing, construction, and cloning are outside timing.
//! Inputs share storage with their fixtures unless labeled unique; shared mutation includes
//! storage detachment.

#[path = "transform/stereo.rs"]
mod stereo;

use std::collections::BTreeSet;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph::ops::model::{AromaticityModel, StereoModel};
use umol_graph::ops::stereo::{StereoSelection, StereoSite};
use umol_graph::ops::transform::{
    Aromatizer, EnumerateStereoConfig, HydrogenFolder, HydrogenUnfolder, KekulizeConfig, Kekulizer,
    RemoveStereoConfig, RemoveStereoMode, StereoAdder, StereoEnumerator, StereoRemover,
    Transformer,
};
use umol_graph::ops::validate::ValenceInvariantsValidator;
use umol_graph_ir::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomForm, AtomId, BondForm, BondId, IntoIr, Molecule, MoleculeEntries,
    StereoAtomForm, StereoAtomId, StereoBondId, StereoCoset, StereoKind, StereoLigand,
    StereoLigandKind, Stereogenicity, TetrahedralStereoForm,
};
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

fn bench_stereo_adder(c: &mut Criterion) {
    let transformer = StereoAdder::new(&StereoModel::default());
    let retained = mol_dsl_concrete!(
        r#"{:atoms ["C#T0" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th1"}]}"#
    );
    let mut cases = vec![
        (
            "tetrahedral",
            mol_dsl_concrete!(
                r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#
            ),
            mol_dsl_concrete!(
                r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#
            ),
        ),
        (
            "cis_trans",
            mol_dsl_concrete!(
                r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#
            ),
            mol_dsl_concrete!(
                r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct*"}]}"#
            ),
        ),
        (
            "source_evidence",
            mol_dsl_concrete!(
                r#"{:atoms ["C#T1" "C#h" "F" "Cl" "Br" "C#h" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "2#C0"] [5 6 "1"]]}"#
            ),
            mol_dsl_concrete!(
                r#"{:atoms ["C#T1" "C#h" "F" "Cl" "Br" "C#h" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "2#C0"] [5 6 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
            :stereo-bonds [{:site 4 :ligands [0 [:h 1] 6 [:h 5]] :attrs "Ct0"}]}"#
            ),
        ),
        ("retained", retained.clone(), retained),
    ];
    let mut entries = MoleculeEntries::default();
    let mut stereo_atoms = Vec::new();
    for component in 0..16 {
        let base = component * 5;
        entries.atoms.extend(
            [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
                .map(|element| AtomForm::from_element(element).into_concrete()),
        );
        entries.bonds.extend((1..=4).map(|ligand| {
            (
                AtomId(base),
                AtomId(base + ligand),
                BondForm::from_order(1).into_concrete(),
            )
        }));
        stereo_atoms.push((
            AtomId(base),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(base + ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        ));
    }
    let input = Molecule::from_entries(entries.clone());
    entries.stereo_atoms = stereo_atoms;
    cases.push(("sites_16", input, Molecule::from_entries(entries)));

    let mut group = c.benchmark_group("stereo_adder");
    for (name, input, expected) in cases {
        let entries = MoleculeEntries {
            atoms: input
                .atoms()
                .iter()
                .map(|atom| atom.attributes().clone())
                .collect(),
            bonds: input
                .bonds()
                .iter()
                .map(|bond| {
                    let [a, b] = bond.atom_ids();
                    (a, b, bond.attributes().clone())
                })
                .collect(),
            stereo_atoms: input
                .stereo_atoms()
                .iter()
                .map(|stereo| {
                    (
                        stereo.site_id(),
                        stereo.ligand_ids().to_vec(),
                        stereo.attributes().clone(),
                    )
                })
                .collect(),
            stereo_bonds: input
                .stereo_bonds()
                .iter()
                .map(|stereo| {
                    (
                        stereo.site_id(),
                        stereo.ligand_ids().to_vec(),
                        stereo.attributes().clone(),
                    )
                })
                .collect(),
            ..Default::default()
        };
        assert_eq!(Molecule::from_entries(entries.clone()), input, "{name}");
        assert_eq!(
            transformer.transform(input.clone()),
            Ok(expected.clone()),
            "{name}"
        );
        let mut in_place = input.clone();
        assert_eq!(transformer.transform_into(&mut in_place), Ok(()), "{name}");
        assert_eq!(in_place, expected, "{name}");
        assert_eq!(
            transformer.transform(Molecule::from_entries(entries.clone())),
            Ok(expected.clone()),
            "{name}"
        );
        let mut unique = Molecule::from_entries(entries.clone());
        assert_eq!(transformer.transform_into(&mut unique), Ok(()), "{name}");
        assert_eq!(unique, expected, "{name}");
        assert_eq!(
            expected.clone().edit().finish(),
            Ok(expected.clone()),
            "{name}"
        );
        eprintln!("stereo_adder_input\t{name}\tatoms={}\tbonds={}\tretained_atoms={}\tretained_bonds={}\tadded_atoms={}\tadded_bonds={}\talgorithm=Nauty\tpara_stereo=false\tminimum_ring_size=8",
            input.atoms().count(), input.bonds().count(), input.stereo_atoms().count(), input.stereo_bonds().count(),
            expected.stereo_atoms().count() - input.stereo_atoms().count(), expected.stereo_bonds().count() - input.stereo_bonds().count());

        group.bench_function(BenchmarkId::new("transform_shared", name), |b| {
            b.iter_batched(
                || input.clone(),
                |molecule| transformer.transform(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_into_shared", name), |b| {
            b.iter_batched_ref(
                || input.clone(),
                |molecule| transformer.transform_into(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_unique", name), |b| {
            b.iter_batched(
                || Molecule::from_entries(entries.clone()),
                |molecule| transformer.transform(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_into_unique", name), |b| {
            b.iter_batched_ref(
                || Molecule::from_entries(entries.clone()),
                |molecule| transformer.transform_into(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_stereo_remover(c: &mut Criterion) {
    let mixed = mol_dsl_concrete!(
        r#"{:atoms ["C#T1#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#C0#R!"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#
    );
    let mut cases = vec![
        (
            "mixed_all",
            mixed.clone(),
            RemoveStereoConfig::default(),
            mol_dsl_concrete!(
                r#"{:atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]}"#
            ),
        ),
        (
            "mixed_undetermined",
            mixed,
            RemoveStereoConfig {
                mode: RemoveStereoMode::Undetermined,
                ..Default::default()
            },
            mol_dsl_concrete!(
                r#"{:atoms ["C#T1#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]
                :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]}"#
            ),
        ),
    ];
    let mut entries = MoleculeEntries::default();
    let mut without_stereo = MoleculeEntries::default();
    let mut without_first = MoleculeEntries::default();
    for component in 0..16 {
        let base = component * 5;
        let plain = [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| AtomForm::from_element(element).into_concrete());
        let mut asserted = plain.clone();
        asserted[0] = asserted[0]
            .clone()
            .with_constraint(AtomConstraintForm::tetrahedral_stereo(
                TetrahedralStereoForm::Stereo(StereoCoset::Undetermined),
            ));
        entries.atoms.extend(asserted.clone());
        without_stereo.atoms.extend(plain.clone());
        without_first
            .atoms
            .extend(if component == 0 { plain } else { asserted });
        let bonds: Vec<_> = (1..=4)
            .map(|ligand| {
                (
                    AtomId(base),
                    AtomId(base + ligand),
                    BondForm::from_order(1).into_concrete(),
                )
            })
            .collect();
        entries.bonds.extend(bonds.clone());
        without_stereo.bonds.extend(bonds.clone());
        without_first.bonds.extend(bonds);
        let stereo = (
            AtomId(base),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(base + ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        );
        entries.stereo_atoms.push(stereo.clone());
        if component > 0 {
            without_first.stereo_atoms.push(stereo);
        }
    }
    let input = Molecule::from_entries(entries);
    let stripped = Molecule::from_entries(without_stereo);
    cases.push((
        "all_16",
        input.clone(),
        RemoveStereoConfig::default(),
        stripped.clone(),
    ));
    cases.push((
        "selected_1",
        input.clone(),
        RemoveStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
                ..Default::default()
            },
            ..Default::default()
        },
        Molecule::from_entries(without_first),
    ));
    cases.push((
        "stereogenic_16",
        input,
        RemoveStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
                ..Default::default()
            },
            ..Default::default()
        },
        stripped,
    ));

    let mut group = c.benchmark_group("stereo_remover");
    for (name, input, config, expected) in cases {
        let transformer = StereoRemover::with_config(&StereoModel::default(), config.clone());
        let entries = MoleculeEntries {
            atoms: input
                .atoms()
                .iter()
                .map(|atom| atom.attributes().clone())
                .collect(),
            bonds: input
                .bonds()
                .iter()
                .map(|bond| {
                    let [a, b] = bond.atom_ids();
                    (a, b, bond.attributes().clone())
                })
                .collect(),
            stereo_atoms: input
                .stereo_atoms()
                .iter()
                .map(|stereo| {
                    (
                        stereo.site_id(),
                        stereo.ligand_ids().to_vec(),
                        stereo.attributes().clone(),
                    )
                })
                .collect(),
            stereo_bonds: input
                .stereo_bonds()
                .iter()
                .map(|stereo| {
                    (
                        stereo.site_id(),
                        stereo.ligand_ids().to_vec(),
                        stereo.attributes().clone(),
                    )
                })
                .collect(),
            ..Default::default()
        };
        assert_eq!(Molecule::from_entries(entries.clone()), input, "{name}");
        assert_eq!(
            transformer.transform(input.clone()),
            Ok(expected.clone()),
            "{name}"
        );
        let mut in_place = input.clone();
        assert_eq!(transformer.transform_into(&mut in_place), Ok(()), "{name}");
        assert_eq!(in_place, expected, "{name}");
        assert_eq!(
            transformer.transform(Molecule::from_entries(entries.clone())),
            Ok(expected.clone()),
            "{name}"
        );
        let mut unique = Molecule::from_entries(entries.clone());
        assert_eq!(transformer.transform_into(&mut unique), Ok(()), "{name}");
        assert_eq!(unique, expected, "{name}");
        assert_eq!(
            expected.clone().edit().finish(),
            Ok(expected.clone()),
            "{name}"
        );
        eprintln!("stereo_remover_input\t{name}\tatoms={}\tbonds={}\tstereo_atoms={}\tstereo_bonds={}\tremoved_stereo_atoms={}\tremoved_stereo_bonds={}\tmode={:?}\tselected_sites={:?}\tclasses={:?}\talgorithm=Nauty\tpara_stereo=false",
            input.atoms().count(), input.bonds().count(), input.stereo_atoms().count(), input.stereo_bonds().count(),
            input.stereo_atoms().count() - expected.stereo_atoms().count(), input.stereo_bonds().count() - expected.stereo_bonds().count(),
            config.mode, config.selection.sites.as_ref().map(BTreeSet::len), config.selection.stereogenicities);
        group.bench_function(BenchmarkId::new("transform_shared", name), |b| {
            b.iter_batched(
                || input.clone(),
                |molecule| transformer.transform(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_into_shared", name), |b| {
            b.iter_batched_ref(
                || input.clone(),
                |molecule| transformer.transform_into(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_unique", name), |b| {
            b.iter_batched(
                || Molecule::from_entries(entries.clone()),
                |molecule| transformer.transform(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_into_unique", name), |b| {
            b.iter_batched_ref(
                || Molecule::from_entries(entries.clone()),
                |molecule| transformer.transform_into(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_stereo_enumerator(c: &mut Criterion) {
    let mixed = mol_dsl_concrete!(
        r#"{:atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th*"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#
    );
    let mixed_expected: Vec<_> = [(0, 0), (0, 1), (1, 0), (1, 1)]
        .into_iter()
        .map(|(atom, bond)| {
            let mut result = mixed.clone();
            result
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = (StereoKind::Tetrahedral, atom).into();
            result
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration = (StereoKind::CisTrans, bond).into();
            result
        })
        .collect();
    let evidenced = mol_dsl_concrete!(
        r#"{:atoms ["C#D4#T0" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!#C1"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th*"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#
    );
    let evidenced_expected = mol_dsl_concrete!(
        r#"{:atoms ["C#D4#T0" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!#C1"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th1"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct1"}]}"#
    );
    let symmetric = mol_dsl_concrete!(
        r#"{:atoms ["C" "H" "H" "H" "H" "C" "H" "H" "H" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "1"] [5 7 "1"] [5 8 "1"] [5 9 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"} {:site 5 :ligands [6 7 8 9] :attrs "Th*"}]}"#
    );
    let symmetric_expected: Vec<_> = [(0, 0), (0, 1), (1, 0), (1, 1)]
        .into_iter()
        .map(|(first, second)| {
            let mut result = symmetric.clone();
            result
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = (StereoKind::Tetrahedral, first).into();
            result
                .stereo_atom_mut(StereoAtomId(1))
                .attributes_mut()
                .configuration = (StereoKind::Tetrahedral, second).into();
            result
        })
        .collect();
    let mut cases = vec![
        (
            "mixed",
            mixed,
            EnumerateStereoConfig::default(),
            vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))],
            mixed_expected,
            true,
        ),
        (
            "evidenced",
            evidenced,
            EnumerateStereoConfig::default(),
            vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))],
            vec![evidenced_expected],
            true,
        ),
        (
            "symmetric_2",
            symmetric,
            EnumerateStereoConfig::default(),
            vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(5))],
            symmetric_expected,
            true,
        ),
    ];
    for count in [8u32, 128] {
        let mut entries = MoleculeEntries::default();
        for component in 0..count {
            let base = component * 5;
            entries.atoms.extend(
                [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
                    .map(|element| AtomForm::from_element(element).into_concrete()),
            );
            entries.bonds.extend((1..5).map(|ligand| {
                (
                    AtomId(base),
                    AtomId(base + ligand),
                    BondForm::from_order(1).into_concrete(),
                )
            }));
            entries.stereo_atoms.push((
                AtomId(base),
                (1..5)
                    .map(|ligand| StereoLigand::new(AtomId(base + ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            ));
        }
        let input = Molecule::from_entries(entries);
        let order: Vec<_> = (0..count)
            .map(|index| StereoSite::Atom(AtomId(index * 5)))
            .collect();
        let full = count == 8;
        let expected: Vec<_> = (0u32..if full { 256 } else { 4 })
            .map(|bits| {
                let mut result = input.clone();
                for index in 0..count {
                    let shift = count - 1 - index;
                    let value = if shift < u32::BITS {
                        (bits >> shift) & 1
                    } else {
                        0
                    };
                    result
                        .stereo_atom_mut(StereoAtomId(index))
                        .attributes_mut()
                        .configuration = (StereoKind::Tetrahedral, value).into();
                }
                result
            })
            .collect();
        if full {
            cases.push((
                "all_8",
                input.clone(),
                EnumerateStereoConfig::default(),
                order.clone(),
                expected.clone(),
                true,
            ));
            cases.push((
                "classified_8",
                input,
                EnumerateStereoConfig {
                    selection: StereoSelection {
                        stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
                        ..Default::default()
                    },
                    ..Default::default()
                },
                order,
                expected,
                true,
            ));
        } else {
            let selected_expected: Vec<_> = [0, 1]
                .into_iter()
                .map(|value| {
                    let mut result = input.clone();
                    result
                        .stereo_atom_mut(StereoAtomId(0))
                        .attributes_mut()
                        .configuration = (StereoKind::Tetrahedral, value).into();
                    result
                })
                .collect();
            cases.push((
                "selected_1_of_128",
                input.clone(),
                EnumerateStereoConfig {
                    selection: StereoSelection {
                        sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
                        ..Default::default()
                    },
                    ..Default::default()
                },
                order.clone(),
                selected_expected,
                true,
            ));
            cases.push((
                "prefix_128",
                input,
                EnumerateStereoConfig::default(),
                order,
                expected,
                false,
            ));
        }
    }
    let mut group = c.benchmark_group("stereo_enumerator");
    for (name, input, config, order, expected, full) in cases {
        let engine = StereoEnumerator::with_config(&StereoModel::default(), config.clone(), order);
        let entries = MoleculeEntries {
            atoms: input
                .atoms()
                .iter()
                .map(|atom| atom.attributes().clone())
                .collect(),
            bonds: input
                .bonds()
                .iter()
                .map(|bond| {
                    let [a, b] = bond.atom_ids();
                    (a, b, bond.attributes().clone())
                })
                .collect(),
            stereo_atoms: input
                .stereo_atoms()
                .iter()
                .map(|stereo| {
                    (
                        stereo.site_id(),
                        stereo.ligand_ids().to_vec(),
                        stereo.attributes().clone(),
                    )
                })
                .collect(),
            stereo_bonds: input
                .stereo_bonds()
                .iter()
                .map(|stereo| {
                    (
                        stereo.site_id(),
                        stereo.ligand_ids().to_vec(),
                        stereo.attributes().clone(),
                    )
                })
                .collect(),
            ..Default::default()
        };
        assert_eq!(Molecule::from_entries(entries.clone()), input, "{name}");
        assert_eq!(
            engine.transform(input.clone()),
            Ok(expected[0].clone()),
            "{name}"
        );
        let mut in_place = input.clone();
        assert_eq!(engine.transform_into(&mut in_place), Ok(()), "{name}");
        assert_eq!(in_place, expected[0], "{name}");
        assert_eq!(
            engine.transform(Molecule::from_entries(entries.clone())),
            Ok(expected[0].clone()),
            "{name}"
        );
        let mut unique = Molecule::from_entries(entries.clone());
        assert_eq!(engine.transform_into(&mut unique), Ok(()), "{name}");
        assert_eq!(unique, expected[0], "{name}");
        assert_eq!(
            engine.transform_iter(&input).take(4).collect::<Vec<_>>(),
            expected[..expected.len().min(4)],
            "{name}"
        );
        if full {
            assert_eq!(
                engine.transform_iter(&input).collect::<Vec<_>>(),
                expected,
                "{name}"
            );
        }
        for result in &expected {
            assert_eq!(result.clone().edit().finish(), Ok(result.clone()), "{name}");
        }
        let selected = config.selection.sites.as_ref().map_or(
            input.stereo_atoms().count() + input.stereo_bonds().count(),
            BTreeSet::len,
        );
        eprintln!("stereo_enumerator_input\t{name}\tatoms={}\tbonds={}\tstereo_atoms={}\tstereo_bonds={}\tselected_sites={}\tdomain_product=2^{}\tfull_results={:?}\tprefix_results={}\tclasses={:?}\talgorithm=Nauty\tpara_stereo=false",
            input.atoms().count(), input.bonds().count(), input.stereo_atoms().count(), input.stereo_bonds().count(), selected,
            selected, full.then_some(expected.len()), expected.len().min(4), config.selection.stereogenicities);
        group.bench_function(BenchmarkId::new("transform_shared", name), |b| {
            b.iter_batched(
                || input.clone(),
                |molecule| engine.transform(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_into_shared", name), |b| {
            b.iter_batched_ref(
                || input.clone(),
                |molecule| engine.transform_into(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_unique", name), |b| {
            b.iter_batched(
                || Molecule::from_entries(entries.clone()),
                |molecule| engine.transform(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("transform_into_unique", name), |b| {
            b.iter_batched_ref(
                || Molecule::from_entries(entries.clone()),
                |molecule| engine.transform_into(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("iter_setup", name), |b| {
            b.iter(|| engine.transform_iter(black_box(&input)));
        });
        group.bench_function(BenchmarkId::new("iter_first", name), |b| {
            b.iter_batched(
                || engine.transform_iter(&input),
                |mut run| run.next().unwrap(),
                BatchSize::SmallInput,
            );
        });
        group.bench_function(BenchmarkId::new("iter_prefix_4", name), |b| {
            b.iter_batched(
                || engine.transform_iter(&input),
                |run| run.take(4).collect::<Vec<_>>(),
                BatchSize::SmallInput,
            );
        });
        if full {
            group.bench_function(BenchmarkId::new("iter_all", name), |b| {
                b.iter_batched(
                    || engine.transform_iter(&input),
                    |run| run.collect::<Vec<_>>(),
                    BatchSize::SmallInput,
                );
            });
        }
    }
    group.finish();
}

criterion_group!(
    transform,
    bench_transform,
    bench_fold_hydrogens,
    bench_unfold_hydrogens,
    bench_hydrogens_scaling,
    bench_stereo_adder,
    bench_stereo_remover,
    bench_stereo_enumerator,
    stereo::bench_stereo
);
criterion_main!(transform);
