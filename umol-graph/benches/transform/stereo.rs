//! Read-only perception and symmetry costs over explicit stereo domains.

#[path = "../../src/ops/transform/tests/stereo/fixtures.rs"]
mod fixtures;

use std::collections::{BTreeMap, BTreeSet};
use std::hint::black_box;

use criterion::{BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph::ops::model::StereoModel;
use umol_graph::ops::stereo::{StereoDerivation, StereoPerceiver, StereoSelection, StereoSite};
use umol_graph::ops::validate::StereoConformanceValidator;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomForm, AtomId, BondConstraintForm, BondForm, BondId, CisTransStereoForm,
    ElementForm, Entity, GraphSymmetryConfig, LigandPermutation, Molecule, MoleculeEntries,
    NumForm, OrientedLigandPermutation, RingConfig, RingModel, StereoAtomConstraintForm,
    StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId, StereoConfigurationForm,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind, StereoLigandPosition, StereoTerm,
    StereogenicSiteKind, Stereogenicity, StereogenicityForm,
};
use umol_graph_ir::mol_dsl_concrete;
use umol_perm::{Orientation, Permutation};

pub(super) fn bench_stereo(c: &mut Criterion) {
    let mut cases = vec![
        (
            "fixture_tetrahedral".to_owned(),
            fixtures::tetrahedral(),
            None,
            None,
        ),
        (
            "fixture_prostereogenic".to_owned(),
            fixtures::prostereogenic(),
            None,
            Some(vec![Stereogenicity::Prostereogenic]),
        ),
        (
            "fixture_nonstereogenic".to_owned(),
            fixtures::nonstereogenic(),
            None,
            Some(vec![Stereogenicity::Nonstereogenic]),
        ),
        (
            "fixture_implicit_hydrogen".to_owned(),
            fixtures::implicit_hydrogen(),
            None,
            None,
        ),
        (
            "fixture_lone_pair".to_owned(),
            fixtures::lone_pair(),
            None,
            None,
        ),
        (
            "fixture_hydrogen_and_lone_pair".to_owned(),
            fixtures::hydrogen_and_lone_pair(),
            None,
            None,
        ),
        (
            "fixture_cis_trans".to_owned(),
            fixtures::cis_trans(),
            None,
            None,
        ),
        ("fixture_imine".to_owned(), fixtures::imine(), None, None),
        (
            "fixture_remote_stereo".to_owned(),
            fixtures::remote_stereo(),
            None,
            None,
        ),
        (
            "fixture_small_ring".to_owned(),
            fixtures::small_ring(),
            Some(6),
            None,
        ),
        (
            "fixture_square_planar".to_owned(),
            fixtures::square_planar(),
            None,
            None,
        ),
        (
            "fixture_evidenced_stereo".to_owned(),
            fixtures::evidenced_stereo(),
            None,
            None,
        ),
        (
            "fixture_partial_stereo".to_owned(),
            fixtures::partial_stereo(),
            None,
            None,
        ),
    ];
    let mut tetrahedral_literal = fixtures::tetrahedral();
    tetrahedral_literal
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32);
    cases.push((
        "fixture_tetrahedral_literal".into(),
        tetrahedral_literal,
        None,
        Some(vec![Stereogenicity::Stereogenic]),
    ));
    let mut cis_trans_literal = fixtures::cis_trans();
    cis_trans_literal
        .stereo_bond_mut(StereoBondId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::CisTrans, 0u32);
    cases.push((
        "fixture_cis_trans_literal".into(),
        cis_trans_literal,
        None,
        Some(vec![Stereogenicity::Stereogenic]),
    ));

    for carbons in [1usize, 8, 32, 128] {
        let mut entries = MoleculeEntries {
            atoms: [Element::C, Element::F, Element::Cl, Element::Br]
                .into_iter()
                .map(|element| AtomForm::from_element(element).into_concrete())
                .collect(),
            bonds: (1..4)
                .map(|ligand| {
                    (
                        AtomId(0),
                        AtomId(ligand),
                        BondForm::from_order(1).into_concrete(),
                    )
                })
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            )],
            ..Default::default()
        };
        for carbon in 0..carbons {
            entries.atoms.push(
                AtomForm::from_element(Element::C)
                    .with_implicit_hydrogens(if carbon + 1 == carbons { 3 } else { 2 })
                    .into_concrete(),
            );
            entries.bonds.push((
                if carbon == 0 {
                    AtomId(0)
                } else {
                    AtomId::from(3 + carbon)
                },
                AtomId::from(4 + carbon),
                BondForm::from_order(1).into_concrete(),
            ));
        }
        cases.push((
            format!("atoms_{}", entries.atoms.len()),
            Molecule::from_entries(entries),
            None,
            None,
        ));
    }

    cases.extend([
        (
            "bromochloroethanol".to_owned(),
            mol_dsl_concrete!(
                r#"{
        :atoms ["C" "C#h3" "Br" "Cl" "O#h"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]
    }"#
            ),
            None,
            Some(vec![Stereogenicity::Stereogenic]),
        ),
        (
            "dimethylcyclohexane_1_2".to_owned(),
            mol_dsl_concrete!(
                r#"{
        :atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2" "C#h3" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 2 7 [:h 1]] :attrs "Th*"}]
    }"#
            ),
            None,
            None,
        ),
        (
            "trimethylcyclohexane_1_2_3".to_owned(),
            mol_dsl_concrete!(
                r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h3" "C#h3" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [1 7 "1"] [2 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 2 7 [:h 1]] :attrs "Th*"}
                       {:site 2 :ligands [1 3 8 [:h 2]] :attrs "Th*"}]
    }"#
            ),
            None,
            None,
        ),
        (
            "tetramethylcyclobutane".to_owned(),
            mol_dsl_concrete!(
                r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h3" "C#h3" "C#h3" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]
                [0 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 2 5 [:h 1]] :attrs "Th*"}
                       {:site 2 :ligands [1 3 6 [:h 2]] :attrs "Th*"}
                       {:site 3 :ligands [0 2 7 [:h 3]] :attrs "Th*"}]
    }"#
            ),
            None,
            None,
        ),
        (
            "decalin".to_owned(),
            mol_dsl_concrete!(
                r#"{
        :atoms ["C#h" "C#h2" "C#h2" "C#h2" "C#h2" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 0 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 9 [:h 0]] :attrs "Th*"}
                       {:site 5 :ligands [0 4 6 [:h 5]] :attrs "Th*"}]
    }"#
            ),
            None,
            None,
        ),
    ]);

    for size in [4usize, 6, 8, 12] {
        let entries = MoleculeEntries {
            atoms: (0..size)
                .map(|atom| {
                    AtomForm::from_element(Element::C)
                        .with_implicit_hydrogens(if atom < 2 { 1 } else { 2 })
                        .into_concrete()
                })
                .collect(),
            bonds: (0..size)
                .map(|bond| {
                    let form = if bond == 0 {
                        BondForm::from_order(2).with_constraint(BondConstraintForm::CisTransStereo(
                            CisTransStereoForm::Stereo,
                        ))
                    } else {
                        BondForm::from_order(1)
                    };
                    (
                        AtomId::from(bond),
                        AtomId::from((bond + 1) % size),
                        form.into_concrete(),
                    )
                })
                .collect(),
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId::from(size - 1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            )],
            ..Default::default()
        };
        cases.push((
            format!("ring_{size}"),
            Molecule::from_entries(entries),
            Some(size),
            None,
        ));
    }

    let model = StereoModel::default();
    let perception = StereoPerceiver::new(&model);
    let config = GraphSymmetryConfig::default();
    eprintln!("stereo_config\talgorithm=Nauty\tinclude_constraints=false\tmax_classification_trials=4096\tminimum_ring_size={}", model.stereo_bond_min_ring_size);
    eprintln!("stereo_columns\tcase\tatoms\tbonds\tstereo_atoms\tstereo_bonds\tring_size\tproper_atom_orbits_size_count\tconfigurations\tperceived_atoms\tperceived_bonds\tinconsistencies\tclassifications");
    let mut recognition_cases = Vec::new();
    for (name, molecule, ring_size, _) in &cases {
        if matches!(
            name.as_str(),
            "fixture_square_planar" | "fixture_partial_stereo" | "fixture_hydrogen_and_lone_pair"
        ) {
            continue;
        }
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: molecule
                .atoms()
                .iter()
                .map(|atom| {
                    let mut form = atom.attributes().clone();
                    form.constraints = Default::default();
                    form
                })
                .collect(),
            bonds: molecule
                .bonds()
                .iter()
                .map(|bond| {
                    let [a, b] = bond.atom_ids();
                    let mut form = bond.attributes().clone();
                    form.constraints = Default::default();
                    (a, b, form)
                })
                .collect(),
            ..Default::default()
        });
        let selection = StereoSelection::default();
        let expected = StereoDerivation {
            stereo_atoms: molecule
                .stereo_atoms()
                .iter()
                .map(|stereo| {
                    (
                        stereo.site_id(),
                        stereo.ligand_ids().to_vec(),
                        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
                    )
                })
                .collect(),
            stereo_bonds: if ring_size
                .is_some_and(|size| size < model.stereo_bond_min_ring_size as usize)
            {
                vec![]
            } else {
                molecule
                    .stereo_bonds()
                    .iter()
                    .map(|stereo| {
                        (
                            stereo.site_id(),
                            stereo.ligand_ids().to_vec(),
                            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
                        )
                    })
                    .collect()
            },
            ..Default::default()
        };
        assert_eq!(
            perception.recognize(
                &source,
                &selection,
                AutomorphismAlgorithm::TrailIr,
                RingConfig::default(),
                4096
            ),
            Ok(expected.clone()),
            "{name}"
        );
        recognition_cases.push((name.clone(), source, selection, expected));
    }
    let mut group = c.benchmark_group("stereo");
    for (name, molecule, ring_size, expected_classifications) in cases {
        let original = molecule.clone();
        assert_eq!(molecule.clone().edit().finish(), Ok(molecule.clone()));
        if let Some(size) = ring_size {
            assert_eq!(
                molecule
                    .rings(RingModel::default(), RingConfig::default())
                    .bond(BondId(0))
                    .smallest_ring_size(),
                Some(size)
            );
        }
        let derivation = perception
            .derive(&molecule, RingConfig::default())
            .into_determined()
            .expect("stereo benchmark fixture has a determined derivation");
        if name == "fixture_evidenced_stereo" {
            assert_eq!(
                derivation.stereo_atoms,
                vec![(
                    AtomId(0),
                    (1..=4)
                        .map(|ligand| { StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom) })
                        .collect(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined)
                )]
            );
            assert_eq!(
                derivation.stereo_bonds,
                vec![(
                    BondId(4),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined)
                )]
            );
            assert_eq!(derivation.inconsistencies, vec![]);
        }
        let graph = molecule.graph_symmetry(&config).into_graph_symmetry();
        let classifications: Vec<_> = molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| stereo.with_symmetry(&graph).stereogenicity())
            .chain(
                molecule
                    .stereo_bonds()
                    .iter()
                    .map(|stereo| stereo.with_symmetry(&graph).stereogenicity()),
            )
            .collect();
        if let Some(expected) = expected_classifications {
            assert_eq!(classifications, expected, "{name}");
        }
        let mut orbit_histogram = BTreeMap::<usize, usize>::new();
        for atom in molecule.atoms().ids() {
            let orbit = molecule.proper_orbit_of(&graph, atom);
            if orbit[0] == atom {
                *orbit_histogram.entry(orbit.len()).or_default() += 1;
            }
        }
        let configurations: Vec<_> = molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| format!("{:?}", stereo.attributes().configuration))
            .chain(
                molecule
                    .stereo_bonds()
                    .iter()
                    .map(|stereo| format!("{:?}", stereo.attributes().configuration)),
            )
            .collect();
        eprintln!("stereo_input\t{name}\t{}\t{}\t{}\t{}\t{ring_size:?}\t{orbit_histogram:?}\t{configurations:?}\t{}\t{}\t{:?}\t{classifications:?}",
            molecule.atoms().count(), molecule.bonds().count(), molecule.stereo_atoms().count(), molecule.stereo_bonds().count(),
            derivation.stereo_atoms.len(), derivation.stereo_bonds.len(), derivation.inconsistencies);
        assert_eq!(molecule, original);
        drop(original);

        group.bench_function(BenchmarkId::new("unchanged", &name), |b| {
            b.iter(|| black_box(black_box(&molecule)));
        });
        group.bench_function(BenchmarkId::new("perception", &name), |b| {
            b.iter(|| black_box(perception.derive(black_box(&molecule), RingConfig::default())));
        });
        group.bench_function(BenchmarkId::new("graph_symmetry", &name), |b| {
            b.iter(|| black_box(black_box(&molecule).graph_symmetry(black_box(&config))));
        });
        group.bench_function(BenchmarkId::new("stereogenicity", &name), |b| {
            b.iter(|| {
                for stereo in molecule.stereo_atoms().iter() {
                    black_box(stereo.with_symmetry(black_box(&graph)).stereogenicity());
                }
                for stereo in molecule.stereo_bonds().iter() {
                    black_box(stereo.with_symmetry(black_box(&graph)).stereogenicity());
                }
            });
        });
    }
    let meso = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th1"}]
    }"#
    );
    let mut chiral = meso.clone();
    chiral
        .stereo_atom_mut(StereoAtomId(1))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32);
    let mut raw_term = meso.clone();
    raw_term
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(
        StereoKind::Tetrahedral,
        StereoCoset::term(StereoTerm::Lit(0)),
    );
    let mut raw_singleton = meso.clone();
    raw_singleton.atom_mut(AtomId(0)).attributes_mut().element = ElementForm::lit_set([Element::C]);
    raw_singleton
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration =
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::lit_set([0]));
    let mut constrained = meso.clone();
    constrained.atom_mut(AtomId(2)).attributes_mut().constraints =
        AtomConstraintForm::Valence(NumForm::lit_set([1])).into();
    constrained
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints = StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([
        Stereogenicity::Stereogenic,
    ]))
    .into();
    let prostereogenic = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 5] :attrs "Th0"}
                       {:site 1 :ligands [0 2 6 7] :attrs "Th0"}]
    }"#
    );
    let trichloropentane_meso = mol_dsl_concrete!(
        r#"{
        :atoms ["C#h3" "C" "C" "C" "C#h3" "Cl" "Cl" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]
                [1 8 "1"] [2 9 "1"] [3 10 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 5 8] :attrs "Th0"}
                       {:site 2 :ligands [1 3 6 9] :attrs "Th0"}
                       {:site 3 :ligands [2 4 7 10] :attrs "Th0"}]
    }"#
    );
    let mut trichloropentane_chiral = trichloropentane_meso.clone();
    trichloropentane_chiral
        .stereo_atom_mut(StereoAtomId(2))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32);
    let ligand_swap = OrientedLigandPermutation {
        permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
        orientation: Orientation::Proper,
    };
    for (name, molecule) in [
        ("readings_meso", meso),
        ("readings_chiral", chiral),
        ("readings_prostereogenic", prostereogenic),
        ("readings_trichloropentane_meso", trichloropentane_meso),
        ("readings_trichloropentane_chiral", trichloropentane_chiral),
        ("readings_raw_term", raw_term),
        ("readings_raw_singleton", raw_singleton),
        ("readings_constraints", constrained),
    ] {
        assert_eq!(molecule.clone().edit().finish(), Ok(molecule.clone()));
        let graph = molecule.graph_symmetry(&config).into_graph_symmetry();
        let proper: Vec<_> = molecule
            .atoms()
            .ids()
            .map(|id| molecule.proper_orbit_of(&graph, id))
            .collect();
        let star: Vec<_> = molecule
            .atoms()
            .ids()
            .map(|id| molecule.star_orbit_of(&graph, id))
            .collect();
        eprintln!("stereo_readings_input\t{name}\tcounts={:?}\tchiral={}\tproper={proper:?}\tstar={star:?}",
            molecule.entity_counts(), molecule.is_chiral(&graph));
        group.bench_function(BenchmarkId::new("graph_symmetry", name), |b| {
            b.iter(|| black_box(black_box(&molecule).graph_symmetry(black_box(&config))));
        });
        group.bench_function(BenchmarkId::new("stereogenicity", name), |b| {
            b.iter(|| {
                let molecule = black_box(&molecule);
                let graph = black_box(&graph);
                for id in molecule.stereo_atoms().ids() {
                    black_box(
                        molecule
                            .stereo_atom(id)
                            .with_symmetry(graph)
                            .stereogenicity(),
                    );
                }
                for id in molecule.stereo_bonds().ids() {
                    black_box(
                        molecule
                            .stereo_bond(id)
                            .with_symmetry(graph)
                            .stereogenicity(),
                    );
                }
            });
        });
        group.bench_function(BenchmarkId::new("build_molecule_readings", name), |b| {
            b.iter(|| {
                let molecule = black_box(&molecule);
                let graph = molecule
                    .graph_symmetry(black_box(&config))
                    .into_graph_symmetry();
                black_box(molecule.is_chiral(&graph));
                black_box(molecule.same_proper_orbit(&graph, AtomId(0), AtomId(1)));
                black_box(molecule.same_star_orbit(&graph, AtomId(0), AtomId(1)));
                for id in molecule.atoms().ids() {
                    black_box(molecule.proper_orbit_of(&graph, id));
                    black_box(molecule.star_orbit_of(&graph, id));
                }
            });
        });
        group.bench_function(BenchmarkId::new("build_stereo_readings", name), |b| {
            b.iter(|| {
                let molecule = black_box(&molecule);
                let graph = molecule
                    .graph_symmetry(black_box(&config))
                    .into_graph_symmetry();
                for stereo in molecule.stereo_atoms().iter() {
                    let stereo = stereo.with_symmetry(&graph);
                    black_box(stereo.stereogenicity());
                    for first in 0..stereo.kind().degree() {
                        for second in first + 1..stereo.kind().degree() {
                            black_box(stereo.topicity(
                                StereoLigandPosition::from(first),
                                StereoLigandPosition::from(second),
                            ));
                        }
                    }
                    black_box(stereo.ligand_symmetry(ligand_swap));
                }
                for stereo in molecule.stereo_bonds().iter() {
                    let stereo = stereo.with_symmetry(&graph);
                    black_box(stereo.stereogenicity());
                    for first in 0..stereo.kind().degree() {
                        for second in first + 1..stereo.kind().degree() {
                            black_box(stereo.topicity(
                                StereoLigandPosition::from(first),
                                StereoLigandPosition::from(second),
                            ));
                        }
                    }
                    black_box(stereo.ligand_symmetry(ligand_swap));
                }
            });
        });
    }
    group.finish();

    let mut group = c.benchmark_group("stereo_recognition");
    for (name, source, selection, expected) in recognition_cases {
        let original = source.clone();
        let sites: BTreeSet<_> = expected
            .stereo_atoms
            .iter()
            .map(|(site, _, _)| StereoSite::Atom(*site))
            .chain(
                expected
                    .stereo_bonds
                    .iter()
                    .map(|(site, _, _)| StereoSite::Bond(*site)),
            )
            .collect();
        let mut editor = source.clone().edit();
        for (site, ligands, form) in &expected.stereo_atoms {
            editor.add_stereo_atom(*site, ligands, form.clone());
        }
        for (site, ligands, form) in &expected.stereo_bonds {
            editor.add_stereo_bond(*site, ligands, form.clone());
        }
        let prospective = editor.finish().unwrap();
        let graph = prospective.graph_symmetry(&config).into_graph_symmetry();
        let mut classifications = prospective
            .stereogenic_site_kinds(&graph)
            .into_iter()
            .map(|(entity, kind)| {
                let site = match entity {
                    Entity::StereoAtom(id) => {
                        StereoSite::Atom(prospective.stereo_atom(id).site_id())
                    }
                    Entity::StereoBond(id) => {
                        StereoSite::Bond(prospective.stereo_bond(id).site_id())
                    }
                    _ => unreachable!(),
                };
                (site, kind)
            })
            .collect::<Vec<_>>();
        classifications.sort_unstable_by_key(|&(site, _)| site);
        let requested: Vec<_> = sites.iter().copied().collect();
        assert_eq!(
            perception.classify(&prospective, &requested, &graph),
            Ok(classifications.clone())
        );
        let mut filtered_selection = selection.clone();
        filtered_selection.stereogenic_site_kinds =
            Some(BTreeSet::from([StereogenicSiteKind::Constitutional]));
        let mut filtered = expected.clone();
        filtered.stereo_atoms.retain(|(site, _, _)| {
            classifications[classifications
                .binary_search_by_key(&StereoSite::Atom(*site), |&(site, _)| site)
                .unwrap()]
            .1 == StereogenicSiteKind::Constitutional
        });
        filtered.stereo_bonds.retain(|(site, _, _)| {
            classifications[classifications
                .binary_search_by_key(&StereoSite::Bond(*site), |&(site, _)| site)
                .unwrap()]
            .1 == StereogenicSiteKind::Constitutional
        });
        assert_eq!(
            perception.recognize(
                &source,
                &filtered_selection,
                AutomorphismAlgorithm::TrailIr,
                RingConfig::default(),
                4096
            ),
            Ok(filtered)
        );
        assert_eq!(source, original);
        eprintln!("stereo_recognition_input\t{name}\tatoms={}\tbonds={}\tselected_sites={:?}\tcandidates={}\tclassifications={classifications:?}\tgraph_calls={}\tlocal_projections={}",
            source.atoms().count(), source.bonds().count(), selection.sites.as_ref().map(BTreeSet::len), sites.len(), usize::from(!sites.is_empty()), sites.len());

        group.bench_function(BenchmarkId::new("local_candidates", &name), |b| {
            b.iter(|| {
                black_box(
                    perception
                        .recognize(
                            black_box(&source),
                            black_box(&selection),
                            AutomorphismAlgorithm::TrailIr,
                            RingConfig::default(),
                            4096,
                        )
                        .unwrap(),
                )
            });
        });
        group.bench_function(BenchmarkId::new("filtered_candidates", &name), |b| {
            b.iter(|| {
                black_box(
                    perception
                        .recognize(
                            black_box(&source),
                            black_box(&filtered_selection),
                            AutomorphismAlgorithm::TrailIr,
                            RingConfig::default(),
                            4096,
                        )
                        .unwrap(),
                )
            });
        });
        group.bench_function(BenchmarkId::new("graph_symmetry", &name), |b| {
            b.iter(|| black_box(black_box(&prospective).graph_symmetry(black_box(&config))));
        });
        group.bench_function(BenchmarkId::new("local_classification", &name), |b| {
            b.iter(|| {
                let molecule = black_box(&prospective);
                let graph = black_box(&graph);
                for stereo in molecule.stereo_atoms().iter() {
                    black_box(
                        molecule
                            .stereo_atom(stereo.id())
                            .with_symmetry(graph)
                            .stereogenicity(),
                    );
                }
                for stereo in molecule.stereo_bonds().iter() {
                    black_box(
                        molecule
                            .stereo_bond(stereo.id())
                            .with_symmetry(graph)
                            .stereogenicity(),
                    );
                }
            });
        });
        let requested: Vec<_> = sites.iter().copied().collect();
        group.bench_function(BenchmarkId::new("classify", &name), |b| {
            b.iter(|| {
                black_box(
                    perception
                        .classify(
                            black_box(&prospective),
                            &requested,
                            &prospective.graph_symmetry(&config).into_graph_symmetry(),
                        )
                        .unwrap(),
                )
            });
        });
    }
    let validator = StereoConformanceValidator::new(&model);
    for (name, molecule) in [
        (
            "chlorofluoropropane",
            mol_dsl_concrete!(
                r#"{:atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
                :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 3 4 5] :attrs "Th0"}
                               {:site 1 :ligands [0 2 6 7] :attrs "Th0#g'"}]}"#
            ),
        ),
        (
            "trichloropentane_meso",
            mol_dsl_concrete!(
                r#"{:atoms ["C#h3" "C" "C" "C" "C#h3" "Cl" "Cl" "Cl" "H" "H" "H"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]
                        [1 8 "1"] [2 9 "1"] [3 10 "1"]]
                :stereo-atoms [{:site 1 :ligands [0 2 5 8] :attrs "Th0"}
                               {:site 2 :ligands [1 3 6 9] :attrs "Th0"}
                               {:site 3 :ligands [2 4 7 10] :attrs "Th0"}]}"#
            ),
        ),
        (
            "inositol",
            mol_dsl_concrete!(
                r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"
                        "O#h" "O#h" "O#h" "O#h" "O#h" "O#h"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                        [0 6 "1"] [1 7 "1"] [2 8 "1"] [3 9 "1"] [4 10 "1"] [5 11 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                               {:site 1 :ligands [0 2 7 [:h 1]] :attrs "Th*"}
                               {:site 2 :ligands [1 3 8 [:h 2]] :attrs "Th*"}
                               {:site 3 :ligands [2 4 9 [:h 3]] :attrs "Th*"}
                               {:site 4 :ligands [3 5 10 [:h 4]] :attrs "Th*"}
                               {:site 5 :ligands [0 4 11 [:h 5]] :attrs "Th*"}]}"#
            ),
        ),
        (
            "bis_fluoroethyl_difluorohexane",
            mol_dsl_concrete!(
                r#"{:atoms ["C#h3" "C#h3" "C#h3" "C#h3" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"
                        "F" "F" "F" "F"]
                :bonds [[8 9 "1"] [8 4 "1"] [8 5 "1"] [9 6 "1"] [9 7 "1"] [4 0 "1"] [5 1 "1"]
                        [6 2 "1"] [7 3 "1"] [4 10 "1"] [5 11 "1"] [6 12 "1"] [7 13 "1"]]
                :stereo-atoms [{:site 4 :ligands [0 8 10 [:h 4]] :attrs "Th*"}
                               {:site 5 :ligands [1 8 11 [:h 5]] :attrs "Th*"}
                               {:site 6 :ligands [2 9 12 [:h 6]] :attrs "Th*"}
                               {:site 7 :ligands [3 9 13 [:h 7]] :attrs "Th*"}
                               {:site 8 :ligands [4 5 9 [:h 8]] :attrs "Th*"}
                               {:site 9 :ligands [6 7 8 [:h 9]] :attrs "Th*"}]}"#
            ),
        ),
        (
            "tris_hydroxyethyl_methane",
            mol_dsl_concrete!(
                r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h3" "C#h3" "C#h3" "O#h" "O#h" "O#h"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [2 5 "1"] [3 6 "1"]
                        [1 7 "1"] [2 8 "1"] [3 9 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th*"}
                               {:site 1 :ligands [0 4 7 [:h 1]] :attrs "Th*"}
                               {:site 2 :ligands [0 5 8 [:h 2]] :attrs "Th*"}
                               {:site 3 :ligands [0 6 9 [:h 3]] :attrs "Th*"}]}"#
            ),
        ),
    ] {
        let sites: BTreeSet<_> = molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| StereoSite::Atom(stereo.site_id()))
            .collect();
        let requested: Vec<_> = sites.iter().copied().collect();
        group.bench_function(BenchmarkId::new("classify", name), |b| {
            b.iter(|| {
                black_box(
                    perception
                        .classify(
                            black_box(&molecule),
                            &requested,
                            &molecule.graph_symmetry(&config).into_graph_symmetry(),
                        )
                        .unwrap(),
                )
            });
        });
        group.bench_function(BenchmarkId::new("validate", name), |b| {
            b.iter(|| black_box(validator.validate(black_box(&molecule))));
        });
    }
    group.finish();
}
