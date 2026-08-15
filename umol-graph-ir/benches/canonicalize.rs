//! Aggregate canonicalization benchmarks.
//!
//! Criterion ids include the measured graph's node and edge counts.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph_core::{AutomorphismAlgorithm, Correspondence};
use umol_graph_ir::ir::canonicalize::{
    constitution_partition_descriptors, constraint_blocks, initial_class_keys,
    partition_descriptors, rank_initial_classes, structure_partition, AutomorphismAdapter,
    OrderedPartition,
};
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomForm, AtomId, BondForm, BondId,
    Canonicalize, CanonicalizeContext, CanonicalizeLevel, Constraint, DativeBondForm, DativeBondId,
    IncidenceLevel, Molecule, MoleculeConstraint, MoleculeCorrespondence, MoleculeEntries,
    MulticenterBondForm, MulticenterBondId, NoncovalentBondForm, NoncovalentBondId,
    NoncovalentBondKind, NumForm, StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
};

const ALGORITHM: AutomorphismAlgorithm = AutomorphismAlgorithm::Nauty;

struct CorpusCase {
    name: &'static str,
    molecule: Molecule,
}

fn atom(element: Element) -> AtomForm {
    AtomForm::from_element(element)
}

fn bond(first: u32, second: u32, order: u8) -> (AtomId, AtomId, BondForm) {
    (AtomId(first), AtomId(second), BondForm::from_order(order))
}

fn ligand(atom: u32) -> StereoLigand {
    StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)
}

fn implicit_hydrogen(site: u32) -> StereoLigand {
    StereoLigand::new(AtomId(site), StereoLigandKind::ImplicitHydrogen)
}

fn carbon_graph(atom_count: usize, edges: &[(u32, u32)]) -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![atom(Element::C); atom_count],
        bonds: edges
            .iter()
            .map(|&(first, second)| bond(first, second, 1))
            .collect(),
        ..Default::default()
    })
}

fn ordinary_naphthalene() -> Molecule {
    carbon_graph(
        10,
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 4),
        ],
    )
}

fn disconnected_rings() -> Molecule {
    carbon_graph(
        12,
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 6),
        ],
    )
}

fn overlay_heavy() -> Molecule {
    let mut atoms = [
        Element::C,
        Element::C,
        Element::C,
        Element::C,
        Element::N,
        Element::O,
        Element::F,
        Element::Cl,
    ]
    .into_iter()
    .map(atom)
    .collect::<Vec<_>>();
    atoms[0].constraints = AtomConstraintForm::valence(4).into();

    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 2),
            bond(2, 3, 1),
            bond(3, 0, 1),
            bond(1, 4, 1),
            bond(1, 5, 1),
            bond(2, 6, 1),
            bond(2, 7, 1),
        ],
        dative: vec![(
            vec![AtomId(4), AtomId(5)],
            AtomId(3),
            DativeBondForm::from_order(1),
        )],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
            AromaticSystemForm::default(),
        )],
        multicenter: vec![(
            vec![AtomId(0), AtomId(4), AtomId(5)],
            MulticenterBondForm::default(),
        )],
        noncovalent: vec![(
            AtomId(6),
            AtomId(7),
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(1),
            vec![ligand(0), ligand(2), ligand(4), ligand(5)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        stereo_bonds: vec![(
            BondId(1),
            vec![ligand(0), ligand(4), ligand(3), ligand(6)],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
        )],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]),
            sum: NumForm::Lit(0),
        })
        .into(),
    })
}

fn tetrahedral_stereo() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(atom)
            .collect(),
        bonds: vec![bond(0, 1, 1), bond(0, 2, 1), bond(0, 3, 1), bond(0, 4, 1)],
        stereo_atoms: vec![(
            AtomId(0),
            vec![ligand(1), ligand(2), ligand(3), ligand(4)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    })
}

fn meso_dichlorobutane() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::Cl,
            Element::Cl,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        bonds: vec![
            bond(0, 1, 1),
            bond(0, 2, 1),
            bond(1, 3, 1),
            bond(0, 4, 1),
            bond(1, 5, 1),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![ligand(1), ligand(2), ligand(4), implicit_hydrogen(0)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(1),
                vec![ligand(0), ligand(3), ligand(5), implicit_hydrogen(1)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
}

fn para_stereo_trichloropentane() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::Cl,
            Element::Cl,
            Element::Cl,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 1),
            bond(2, 3, 1),
            bond(3, 4, 1),
            bond(1, 5, 1),
            bond(2, 6, 1),
            bond(3, 7, 1),
        ],
        stereo_atoms: vec![
            (
                AtomId(1),
                vec![ligand(0), ligand(2), ligand(5), implicit_hydrogen(1)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(2),
                vec![ligand(1), ligand(3), ligand(6), implicit_hydrogen(2)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(3),
                vec![ligand(2), ligand(4), ligand(7), implicit_hydrogen(3)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
}

fn para_stereo_cascade() -> Molecule {
    let outer_ligands = vec![ligand(10), ligand(11), ligand(12), ligand(13)];

    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::F,
            Element::Cl,
            Element::Br,
            Element::I,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [2, 3, 4, 5].map(ligand).into(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(1),
                [6, 8, 7, 9].map(ligand).into(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(2),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(3),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::CisTrans, 0u32),
            ),
            (
                AtomId(4),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Axial, 0u32),
            ),
            (
                AtomId(5),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
            ),
            (
                AtomId(6),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(7),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::CisTrans, 0u32),
            ),
            (
                AtomId(8),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Axial, 0u32),
            ),
            (
                AtomId(9),
                outer_ligands,
                StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
            ),
        ],
        ..Default::default()
    })
}

fn corpus() -> [CorpusCase; 7] {
    [
        CorpusCase {
            name: "ordinary_naphthalene",
            molecule: ordinary_naphthalene(),
        },
        CorpusCase {
            name: "disconnected_rings",
            molecule: disconnected_rings(),
        },
        CorpusCase {
            name: "overlay_heavy",
            molecule: overlay_heavy(),
        },
        CorpusCase {
            name: "tetrahedral_stereo",
            molecule: tetrahedral_stereo(),
        },
        CorpusCase {
            name: "meso_dichlorobutane",
            molecule: meso_dichlorobutane(),
        },
        CorpusCase {
            name: "para_stereo_trichloropentane",
            molecule: para_stereo_trichloropentane(),
        },
        CorpusCase {
            name: "para_stereo_cascade",
            molecule: para_stereo_cascade(),
        },
    ]
}

const LEVELS: [IncidenceLevel; 3] = [
    IncidenceLevel::Topology,
    IncidenceLevel::Constitution,
    IncidenceLevel::Full,
];

const OPERATIONS: [(&str, CanonicalizeLevel, bool); 5] = [
    ("topology", CanonicalizeLevel::Topology, false),
    ("constitution", CanonicalizeLevel::Constitution, false),
    ("structure", CanonicalizeLevel::Structure, false),
    ("para_stereo_structure", CanonicalizeLevel::Structure, true),
    ("full", CanonicalizeLevel::Full, true),
];

fn level_name(level: IncidenceLevel) -> &'static str {
    match level {
        IncidenceLevel::Topology => "topology",
        IncidenceLevel::Constitution => "constitution",
        IncidenceLevel::Full => "full",
    }
}

fn graph_size(nodes: usize, edges: usize) -> String {
    format!("n{nodes}_e{edges}")
}

fn incidence_level(level: CanonicalizeLevel) -> IncidenceLevel {
    match level {
        CanonicalizeLevel::Topology => IncidenceLevel::Topology,
        CanonicalizeLevel::Constitution => IncidenceLevel::Constitution,
        CanonicalizeLevel::Structure | CanonicalizeLevel::Full => IncidenceLevel::Full,
    }
}

fn reverse_correspondence(molecule: &Molecule) -> MoleculeCorrespondence {
    fn reverse<Id>(count: usize) -> Correspondence<Id>
    where
        Id: Copy + Ord + From<usize>,
    {
        let images = (0..count).rev().map(Id::from).collect::<Vec<_>>();
        Correspondence::from_images(&images, count)
    }

    MoleculeCorrespondence::new(
        reverse::<AtomId>(molecule.atoms().count()),
        reverse::<BondId>(molecule.bonds().count()),
        reverse::<DativeBondId>(molecule.dative_bonds().count()),
        reverse::<AromaticSystemId>(molecule.aromatic_systems().count()),
        reverse::<MulticenterBondId>(molecule.multicenter_bonds().count()),
        reverse::<NoncovalentBondId>(molecule.noncovalent_bonds().count()),
        reverse::<StereoAtomId>(molecule.stereo_atoms().count()),
        reverse::<StereoBondId>(molecule.stereo_bonds().count()),
    )
}

fn bench_incidence_construction(c: &mut Criterion) {
    let corpus = corpus();

    for level in LEVELS {
        let mut group = c.benchmark_group(format!(
            "canonicalize/incidence_construction/{}",
            level_name(level)
        ));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(level);
            let size = graph_size(
                incidence.graph().node_count(),
                incidence.graph().edge_count(),
            );
            group.bench_function(BenchmarkId::new(case.name, size), |b| {
                b.iter(|| black_box(&case.molecule).incidence_graph(level))
            });
        }
        group.finish();
    }
}

fn bench_initial_class_construction(c: &mut Criterion) {
    let corpus = corpus();

    for level in LEVELS {
        let mut group = c.benchmark_group(format!(
            "canonicalize/initial_class_construction/{}",
            level_name(level)
        ));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(level);
            let size = graph_size(
                incidence.graph().node_count(),
                incidence.graph().edge_count(),
            );
            group.bench_function(BenchmarkId::new(case.name, size), |b| {
                b.iter(|| {
                    let (entity_keys, incidence_keys) =
                        initial_class_keys(black_box(&case.molecule), black_box(&incidence))
                            .expect("benchmark corpus normalizes");
                    rank_initial_classes(&entity_keys, &incidence_keys)
                })
            });
        }
        group.finish();
    }
}

fn bench_adapter_construction(c: &mut Criterion) {
    let corpus = corpus();

    for level in LEVELS {
        let mut group = c.benchmark_group(format!(
            "canonicalize/adapter_construction/{}",
            level_name(level)
        ));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(level);
            let (entity_keys, incidence_keys) = initial_class_keys(&case.molecule, &incidence)
                .expect("benchmark corpus normalizes");
            let classes = rank_initial_classes(&entity_keys, &incidence_keys);
            let adapter = AutomorphismAdapter::new(&incidence, &classes);
            let size = graph_size(adapter.graph().node_count(), adapter.graph().edge_count());
            group.bench_function(BenchmarkId::new(case.name, size), |b| {
                b.iter(|| AutomorphismAdapter::new(black_box(&incidence), black_box(&classes)))
            });
        }
        group.finish();
    }
}

fn bench_adapter_labeling(c: &mut Criterion) {
    let corpus = corpus();

    for level in LEVELS {
        let mut group = c.benchmark_group(format!(
            "canonicalize/adapter_labeling/{}",
            level_name(level)
        ));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(level);
            let (entity_keys, incidence_keys) = initial_class_keys(&case.molecule, &incidence)
                .expect("benchmark corpus normalizes");
            let classes = rank_initial_classes(&entity_keys, &incidence_keys);
            let adapter = AutomorphismAdapter::new(&incidence, &classes);
            let size = graph_size(adapter.graph().node_count(), adapter.graph().edge_count());
            group.bench_function(BenchmarkId::new(case.name, size), |b| {
                b.iter(|| black_box(adapter.automorphisms(ALGORITHM)))
            });
        }
        group.finish();
    }
}

fn bench_refinement(c: &mut Criterion) {
    let corpus = corpus();

    for (operation, level, para_stereo) in OPERATIONS {
        let mut group = c.benchmark_group(format!("canonicalize/refinement/{operation}"));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(incidence_level(level));
            let (entity_keys, incidence_keys) = initial_class_keys(&case.molecule, &incidence)
                .expect("benchmark corpus normalizes");
            let classes = rank_initial_classes(&entity_keys, &incidence_keys);
            let adapter = AutomorphismAdapter::new(&incidence, &classes);
            let size = graph_size(adapter.graph().node_count(), adapter.graph().edge_count());
            group.bench_function(BenchmarkId::new(case.name, size), |b| match level {
                CanonicalizeLevel::Topology => {
                    let descriptors =
                        partition_descriptors(&adapter, &entity_keys, &incidence_keys);
                    b.iter(|| {
                        OrderedPartition::from_descriptors(black_box(&descriptors))
                            .refine(black_box(adapter.graph()))
                    })
                }
                CanonicalizeLevel::Constitution => {
                    let descriptors =
                        constitution_partition_descriptors(&adapter, &entity_keys, &incidence);
                    b.iter(|| {
                        OrderedPartition::from_descriptors(black_box(&descriptors))
                            .refine(black_box(adapter.graph()))
                    })
                }
                CanonicalizeLevel::Structure | CanonicalizeLevel::Full => b.iter(|| {
                    structure_partition(
                        black_box(&case.molecule),
                        black_box(&incidence),
                        black_box(&adapter),
                        black_box(&entity_keys),
                        para_stereo,
                    )
                    .expect("benchmark corpus refines")
                }),
            });
        }
        group.finish();
    }
}

fn bench_constraint_key_construction(c: &mut Criterion) {
    let corpus = corpus();
    let mut group = c.benchmark_group("canonicalize/constraint_key_construction");

    for case in &corpus {
        let counts = molecule_counts(&case.molecule)
            .into_iter()
            .map(|count| count.to_string())
            .collect::<Vec<_>>()
            .join("_");
        group.bench_function(BenchmarkId::new(case.name, counts), |b| {
            b.iter(|| constraint_blocks(black_box(&case.molecule)))
        });
    }

    group.finish();
}

fn bench_remapping(c: &mut Criterion) {
    let corpus = corpus();
    let mut group = c.benchmark_group("canonicalize/remapping");

    for case in &corpus {
        let correspondence = reverse_correspondence(&case.molecule);
        let counts = molecule_counts(&case.molecule)
            .into_iter()
            .map(|count| count.to_string())
            .collect::<Vec<_>>()
            .join("_");
        group.bench_function(BenchmarkId::new(case.name, counts), |b| {
            b.iter(|| black_box(&case.molecule).remap(black_box(&correspondence)))
        });
    }

    group.finish();
}

fn bench_canonicalize(c: &mut Criterion) {
    let corpus = corpus();

    for (operation, level, para_stereo) in OPERATIONS {
        let context = CanonicalizeContext {
            para_stereo,
            automorphism_algorithm: ALGORITHM,
        };
        let mut group = c.benchmark_group(format!("canonicalize/operation/{operation}"));
        for case in &corpus {
            let counts = molecule_counts(&case.molecule)
                .into_iter()
                .map(|count| count.to_string())
                .collect::<Vec<_>>()
                .join("_");
            group.bench_function(BenchmarkId::new(case.name, counts), |b| {
                b.iter_batched(
                    || case.molecule.clone(),
                    |molecule| {
                        black_box(
                            molecule
                                .canonicalize_by(level, &context)
                                .expect("benchmark corpus canonicalizes"),
                        )
                    },
                    BatchSize::SmallInput,
                )
            });
        }
        group.finish();
    }
}

fn molecule_counts(molecule: &Molecule) -> [usize; 8] {
    [
        molecule.atoms().count(),
        molecule.bonds().count(),
        molecule.dative_bonds().count(),
        molecule.aromatic_systems().count(),
        molecule.multicenter_bonds().count(),
        molecule.noncovalent_bonds().count(),
        molecule.stereo_atoms().count(),
        molecule.stereo_bonds().count(),
    ]
}

criterion_group!(
    canonicalize,
    bench_incidence_construction,
    bench_initial_class_construction,
    bench_adapter_construction,
    bench_adapter_labeling,
    bench_refinement,
    bench_constraint_key_construction,
    bench_remapping,
    bench_canonicalize,
);
criterion_main!(canonicalize);
