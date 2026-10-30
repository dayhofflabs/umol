//! Read-only perception and symmetry costs over explicit stereo domains.

#[path = "../../src/ops/transform/tests/stereo/fixtures.rs"]
mod fixtures;

use std::collections::{BTreeMap, BTreeSet};
use std::hint::black_box;

use criterion::{BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph::ops::model::StereoModel;
use umol_graph::ops::stereo::{StereoDerivation, StereoPerception, StereoSelection, StereoSite};
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondConstraintForm, BondForm, BondId, CisTransStereoForm,
    ConstitutionColoring, GraphSymmetryConfig, Molecule, MoleculeEntries, RingConfig, RingModel,
    StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId, StereoConfigurationForm,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind, Stereogenicity,
};

pub(super) fn bench_stereo(c: &mut Criterion) {
    let mut cases = vec![
        (
            "fixture_tetrahedral".to_owned(),
            fixtures::tetrahedral(),
            None,
            None,
        ),
        (
            "fixture_prochiral".to_owned(),
            fixtures::prochiral(),
            None,
            Some(vec![Stereogenicity::Prochiral]),
        ),
        (
            "fixture_symmetric".to_owned(),
            fixtures::symmetric(),
            None,
            Some(vec![Stereogenicity::Symmetric]),
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

    for sites in [1usize, 2, 4, 8, 16] {
        let mut entries = MoleculeEntries::default();
        for component in 0..16 {
            let base = component * 5;
            entries.atoms.extend(
                [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
                    .into_iter()
                    .map(|element| AtomForm::from_element(element).into_concrete()),
            );
            entries.bonds.extend((1..=4).map(|ligand| {
                (
                    AtomId::from(base),
                    AtomId::from(base + ligand),
                    BondForm::from_order(1).into_concrete(),
                )
            }));
            if component < sites {
                entries.stereo_atoms.push((
                    AtomId::from(base),
                    (1..=4)
                        .map(|ligand| {
                            StereoLigand::new(AtomId::from(base + ligand), StereoLigandKind::Atom)
                        })
                        .collect(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
                ));
            }
        }
        cases.push((
            format!("sites_{sites}"),
            Molecule::from_entries(entries),
            None,
            None,
        ));
    }

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
                            CisTransStereoForm::Stereo(StereoCoset::Undetermined),
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
    let perception = StereoPerception::new(&model);
    let config = GraphSymmetryConfig {
        coloring: ConstitutionColoring::full(),
        iterate_to_fixpoint: model.para_stereo,
        max_iterations: 16,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    eprintln!("stereo_config\talgorithm=Nauty\tcoloring=full\tpara_stereo={}\tmax_iterations=16\tminimum_ring_size={}",
        model.para_stereo, model.stereo_bond_minimum_ring_size);
    eprintln!("stereo_columns\tcase\tatoms\tbonds\tstereo_atoms\tstereo_bonds\tring_size\tproper_atom_orbits_size_count\tconfigurations\tperceived_atoms\tperceived_bonds\tskipped_bonds\tinconsistencies\tclassifications");
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
        let sites: BTreeSet<_> = molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| StereoSite::Atom(stereo.site_id()))
            .chain(
                molecule
                    .stereo_bonds()
                    .iter()
                    .map(|stereo| StereoSite::Bond(stereo.site_id())),
            )
            .collect();
        let selection = StereoSelection {
            sites: name.starts_with("sites_").then_some(sites),
            ..Default::default()
        };
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
                .is_some_and(|size| size < model.stereo_bond_minimum_ring_size as usize)
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
                AutomorphismAlgorithm::Nauty,
                RingConfig::default()
            ),
            Ok(expected.clone()),
            "{name}"
        );
        recognition_cases.push((name.clone(), source, selection, expected));
    }
    let mut group = c.benchmark_group("stereo_baseline");
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
        let derivation = perception.derive(&molecule);
        if name == "fixture_evidenced_stereo" {
            assert_eq!(
                derivation.stereo_atoms,
                vec![(
                    AtomId(0),
                    (1..=4)
                        .map(|ligand| { StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom) })
                        .collect(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)
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
                    StereoBondForm::new(StereoKind::CisTrans, 0u32)
                )]
            );
            assert_eq!(derivation.inconsistencies, vec![]);
        }
        let graph = molecule.graph_symmetry(&config);
        let symmetries: Vec<_> = molecule
            .stereo_atoms()
            .ids()
            .map(|id| molecule.stereo_atom_symmetry(&graph, id))
            .chain(
                molecule
                    .stereo_bonds()
                    .ids()
                    .map(|id| molecule.stereo_bond_symmetry(&graph, id)),
            )
            .collect();
        let classifications: Vec<_> = symmetries
            .iter()
            .map(|symmetry| symmetry.stereogenicity())
            .collect();
        if let Some(expected) = expected_classifications {
            assert_eq!(classifications, expected, "{name}");
        }
        let mut orbit_histogram = BTreeMap::<usize, usize>::new();
        for atom in molecule.atoms().ids() {
            let orbit = graph.proper_orbit_of(atom);
            if orbit[0] == atom {
                *orbit_histogram.entry(orbit.len()).or_default() += 1;
            }
        }
        let configurations: Vec<_> = symmetries
            .iter()
            .map(|symmetry| format!("{:?}:{:?}", symmetry.kind(), symmetry.coset()))
            .collect();
        eprintln!("stereo_input\t{name}\t{}\t{}\t{}\t{}\t{ring_size:?}\t{orbit_histogram:?}\t{configurations:?}\t{}\t{}\t{:?}\t{:?}\t{classifications:?}",
            molecule.atoms().count(), molecule.bonds().count(), molecule.stereo_atoms().count(), molecule.stereo_bonds().count(),
            derivation.stereo_atoms.len(), derivation.stereo_bonds.len(), derivation.skipped_stereo_bonds, derivation.inconsistencies);
        assert_eq!(molecule, original);
        drop(original);

        group.bench_function(BenchmarkId::new("unchanged", &name), |b| {
            b.iter(|| black_box(black_box(&molecule)));
        });
        group.bench_function(BenchmarkId::new("perception", &name), |b| {
            b.iter(|| black_box(perception.derive(black_box(&molecule))));
        });
        group.bench_function(BenchmarkId::new("graph_symmetry", &name), |b| {
            b.iter(|| black_box(black_box(&molecule).graph_symmetry(black_box(&config))));
        });
        group.bench_function(BenchmarkId::new("local_symmetry", &name), |b| {
            b.iter(|| {
                let molecule = black_box(&molecule);
                let graph = black_box(&graph);
                for id in molecule.stereo_atoms().ids() {
                    black_box(molecule.stereo_atom_symmetry(graph, id));
                }
                for id in molecule.stereo_bonds().ids() {
                    black_box(molecule.stereo_bond_symmetry(graph, id));
                }
            });
        });
        group.bench_function(BenchmarkId::new("classification", &name), |b| {
            b.iter(|| {
                for symmetry in black_box(&symmetries) {
                    black_box(symmetry.stereogenicity());
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
        let graph = prospective.graph_symmetry(&config);
        let classifications: BTreeMap<_, _> = prospective
            .stereo_atoms()
            .iter()
            .map(|stereo| {
                (
                    StereoSite::Atom(stereo.site_id()),
                    prospective
                        .stereo_atom_symmetry(&graph, stereo.id())
                        .stereogenicity(),
                )
            })
            .chain(prospective.stereo_bonds().iter().map(|stereo| {
                (
                    StereoSite::Bond(stereo.site_id()),
                    prospective
                        .stereo_bond_symmetry(&graph, stereo.id())
                        .stereogenicity(),
                )
            }))
            .collect();
        assert_eq!(
            perception.classify(&prospective, &sites, AutomorphismAlgorithm::Nauty),
            Ok(classifications.clone())
        );
        let mut filtered_selection = selection.clone();
        filtered_selection.stereogenicities = Some(BTreeSet::from([Stereogenicity::Stereogenic]));
        let mut filtered = expected.clone();
        filtered.stereo_atoms.retain(|(site, _, _)| {
            classifications[&StereoSite::Atom(*site)] == Stereogenicity::Stereogenic
        });
        filtered.stereo_bonds.retain(|(site, _, _)| {
            classifications[&StereoSite::Bond(*site)] == Stereogenicity::Stereogenic
        });
        assert_eq!(
            perception.recognize(
                &source,
                &filtered_selection,
                AutomorphismAlgorithm::Nauty,
                RingConfig::default()
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
                            AutomorphismAlgorithm::Nauty,
                            RingConfig::default(),
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
                            AutomorphismAlgorithm::Nauty,
                            RingConfig::default(),
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
                            .stereo_atom_symmetry(graph, stereo.id())
                            .stereogenicity(),
                    );
                }
                for stereo in molecule.stereo_bonds().iter() {
                    black_box(
                        molecule
                            .stereo_bond_symmetry(graph, stereo.id())
                            .stereogenicity(),
                    );
                }
            });
        });
        group.bench_function(BenchmarkId::new("classify", &name), |b| {
            b.iter(|| {
                black_box(
                    perception
                        .classify(
                            black_box(&prospective),
                            black_box(&sites),
                            AutomorphismAlgorithm::Nauty,
                        )
                        .unwrap(),
                )
            });
        });
    }
    group.finish();
}
