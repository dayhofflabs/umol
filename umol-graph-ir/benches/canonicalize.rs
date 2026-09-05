//! Aggregate canonicalization and participant-frame benchmarks.
//!
//! Criterion ids include the measured graph's node and edge counts.

use std::fmt::Debug;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_edn::FromEdn;
use umol_graph_core::{AutomorphismAlgorithm, EdgeId, GraphRemapping, NodeId, Remapping};
use umol_graph_ir::dsl::MoleculeDsl;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomForm, AtomId, BondForm, BondId,
    Canonicalize, CanonicalizeContext, Constraint, DativeBondForm, DativeBondId, IncidenceLevel,
    Molecule, MoleculeConstraint, MoleculeEntries, MoleculeRemapping, MulticenterBondForm,
    MulticenterBondId, NoncovalentBondForm, NoncovalentBondId, NoncovalentBondKind, NumForm,
    Reframe, StereoAtomConstraintForm, StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind, StereoLigandPair, Topicity,
    TopicityForm, TopicityRelationForm,
};

const ALGORITHM: AutomorphismAlgorithm = AutomorphismAlgorithm::Nauty;

struct CorpusCase {
    name: &'static str,
    molecule: Molecule,
}

// Extended C/H propane network, seed flask 0.
const FEATURE_FREE_CONNECTED: &str = r#"
{:atoms ["H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"]
 :bonds [[0 8 "1#c0#u0#s"]
         [1 8 "1#c0#u0#s"]
         [2 8 "1#c0#u0#s"]
         [3 9 "1#c0#u0#s"]
         [4 9 "1#c0#u0#s"]
         [5 9 "1#c0#u0#s"]
         [6 10 "1#c0#u0#s"]
         [7 10 "1#c0#u0#s"]
         [8 10 "1#c0#u0#s"]
         [9 10 "1#c0#u0#s"]]}
"#;

// Extended C/H propane network, product flask 72.
const FEATURE_FREE_DISCONNECTED: &str = r#"
{:atoms ["H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u#s2"
         "C#i=#c0#h0#n0#u#s2"]
 :bonds [[0 8 "1#c0#u0#s"]
         [1 8 "1#c0#u0#s"]
         [2 8 "1#c0#u0#s"]
         [3 8 "1#c0#u0#s"]
         [4 9 "1#c0#u0#s"]
         [5 9 "1#c0#u0#s"]
         [6 10 "1#c0#u0#s"]
         [7 10 "1#c0#u0#s"]
         [9 10 "1#c0#u0#s"]]}
"#;

// Extended C/H ethane network, product flask 99.
const SYMMETRY_HEAVY_RADICALS: &str = r#"
{:atoms ["H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "C#i=#c0#h0#n0#u#s2"
         "C#i=#c0#h0#n0#u#s2"]
 :bonds [[6 7 "3#c0#u0#s"]]}
"#;

const SCALING_CASES: [(&str, &str, [usize; 8]); 3] = [
    (
        "feature_free_connected",
        FEATURE_FREE_CONNECTED,
        [11, 10, 0, 0, 0, 0, 0, 0],
    ),
    (
        "feature_free_disconnected",
        FEATURE_FREE_DISCONNECTED,
        [11, 9, 0, 0, 0, 0, 0, 0],
    ),
    (
        "symmetry_heavy_radicals",
        SYMMETRY_HEAVY_RADICALS,
        [8, 1, 0, 0, 0, 0, 0, 0],
    ),
];

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

fn carbon_path(atom_count: usize) -> Molecule {
    carbon_graph(
        atom_count,
        &(0..atom_count as u32 - 1)
            .map(|first| (first, first + 1))
            .collect::<Vec<_>>(),
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
            vec![AtomId(5), AtomId(4)],
            AtomId(3),
            DativeBondForm::from_order(1),
        )],
        aromatic: vec![(
            vec![AtomId(3), AtomId(2), AtomId(1), AtomId(0)],
            AromaticSystemForm::default(),
        )],
        multicenter: vec![(
            vec![AtomId(5), AtomId(4), AtomId(0)],
            MulticenterBondForm::default(),
        )],
        noncovalent: vec![(
            [AtomId(7), AtomId(6)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(1),
            vec![ligand(5), ligand(4), ligand(2), ligand(0)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        stereo_bonds: vec![(
            BondId(1),
            vec![ligand(6), ligand(3), ligand(4), ligand(0)],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
        )],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]),
            sum: NumForm::Lit(0),
        })
        .into(),
    })
}

fn large_aromatic() -> Molecule {
    const ATOM_COUNT: usize = 128;

    Molecule::from_entries(MoleculeEntries {
        atoms: vec![atom(Element::C); ATOM_COUNT],
        aromatic: vec![(
            (0..ATOM_COUNT).rev().map(AtomId::from).collect::<Vec<_>>(),
            AromaticSystemForm::from_electrons((0..ATOM_COUNT).map(|value| value as i64).collect()),
        )],
        ..Default::default()
    })
}

fn reframing_corpus() -> Vec<CorpusCase> {
    vec![
        CorpusCase {
            name: "empty",
            molecule: Molecule::default(),
        },
        CorpusCase {
            name: "overlay_heavy",
            molecule: overlay_heavy(),
        },
        CorpusCase {
            name: "large_aromatic",
            molecule: large_aromatic(),
        },
    ]
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

fn cis_trans_stereo_bond() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
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
        bonds: vec![
            bond(0, 1, 2),
            bond(0, 2, 1),
            bond(0, 3, 1),
            bond(1, 4, 1),
            bond(1, 5, 1),
        ],
        stereo_bonds: vec![(
            BondId(0),
            vec![ligand(3), ligand(2), ligand(5), ligand(4)],
            StereoBondForm::new(StereoKind::CisTrans, 1u32),
        )],
        ..Default::default()
    })
}

fn mixed_atom_and_bond_stereo() -> Molecule {
    let stereo_atom = tetrahedral_stereo();
    let stereo_bond = cis_trans_stereo_bond();
    Molecule::combine_all([&stereo_atom, &stereo_bond]).0
}

fn frame_relative_stereo_constraint() -> Molecule {
    let constraint = StereoAtomConstraintForm::Topicity(TopicityForm {
        pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
    });
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(atom)
            .collect(),
        bonds: vec![bond(0, 1, 1), bond(0, 2, 1), bond(0, 3, 1), bond(0, 4, 1)],
        stereo_atoms: vec![(
            AtomId(0),
            vec![ligand(1), ligand(2), ligand(3), ligand(4)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                .with_constraint(constraint.clone()),
        )],
        constraints: Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, constraint)
            .into(),
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
        bonds: [
            (0, [2, 3, 4, 5]),
            (1, [6, 8, 7, 9]),
            (2, [10, 11, 12, 13]),
            (3, [10, 11, 12, 13]),
            (4, [10, 11, 12, 13]),
            (5, [10, 11, 12, 13]),
            (6, [10, 11, 12, 13]),
            (7, [10, 11, 12, 13]),
            (8, [10, 11, 12, 13]),
            (9, [10, 11, 12, 13]),
        ]
        .into_iter()
        .flat_map(|(site, ligands)| ligands.into_iter().map(move |ligand| bond(site, ligand, 1)))
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
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
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
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
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

fn corpus() -> Vec<CorpusCase> {
    let mut corpus = vec![
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
            name: "cis_trans_stereo_bond",
            molecule: cis_trans_stereo_bond(),
        },
        CorpusCase {
            name: "mixed_atom_and_bond_stereo",
            molecule: mixed_atom_and_bond_stereo(),
        },
        CorpusCase {
            name: "frame_relative_stereo_constraint",
            molecule: frame_relative_stereo_constraint(),
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
    ];
    corpus.extend(retained_scaling_corpus());
    corpus
}

fn retained_scaling_corpus() -> Vec<CorpusCase> {
    SCALING_CASES
        .map(|(name, source, expected_counts)| {
            let molecule = MoleculeDsl::from_edn_str(source)
                .unwrap_or_else(|error| panic!("benchmark case {name} must parse: {error}"))
                .into_parts()
                .0;
            assert_eq!(
                molecule_counts(&molecule),
                expected_counts,
                "benchmark case {name} changed shape"
            );
            assert!(
                !molecule.has_constraints()
                    && molecule
                        .atoms()
                        .iter()
                        .all(|atom| atom.attributes.constraints.is_empty())
                    && molecule
                        .bonds()
                        .iter()
                        .all(|bond| bond.attributes.constraints.is_empty()),
                "benchmark case {name} must remain feature-free"
            );
            CorpusCase { name, molecule }
        })
        .into()
}

const LEVELS: [IncidenceLevel; 3] = [
    IncidenceLevel::Topology,
    IncidenceLevel::Constitution,
    IncidenceLevel::Full,
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

fn reverse_remapping(molecule: &Molecule) -> MoleculeRemapping {
    fn reverse<Id>(count: usize) -> Remapping<Id>
    where
        Id: Copy + Debug + Into<usize> + From<usize>,
    {
        let images = (0..count).rev().map(Id::from).collect::<Vec<_>>();
        Remapping::new(images).unwrap()
    }

    MoleculeRemapping::new(
        GraphRemapping::new(
            reverse::<NodeId>(molecule.atoms().count()),
            reverse::<EdgeId>(molecule.bonds().count()),
        ),
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

fn bench_remapping(c: &mut Criterion) {
    let corpus = corpus();
    let mut group = c.benchmark_group("canonicalize/remapping");

    for case in &corpus {
        let correspondence = reverse_remapping(&case.molecule);
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

    for (name, context) in [
        (
            "without_para_stereo",
            CanonicalizeContext {
                para_stereo: false,
                automorphism_algorithm: ALGORITHM,
            },
        ),
        (
            "with_para_stereo",
            CanonicalizeContext {
                para_stereo: true,
                automorphism_algorithm: ALGORITHM,
            },
        ),
    ] {
        let mut group = c.benchmark_group(format!("canonicalize/complete/{name}"));
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
                                .canonicalize(&context)
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

fn bench_retained_scaling_cases(c: &mut Criterion) {
    let context = CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: ALGORITHM,
    };
    let cases = retained_scaling_corpus()
        .into_iter()
        .map(|case| {
            let renumbered = case.molecule.remap(&reverse_remapping(&case.molecule));
            let additional_atom = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::O)],
                ..Default::default()
            });
            let (unequal, _) = Molecule::combine_all([&case.molecule, &additional_atom]);
            assert!(case.molecule.canonical_eq(&renumbered, &context));
            assert!(!case.molecule.canonical_eq(&unequal, &context));
            (case, renumbered, unequal)
        })
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("canonicalize/scaling/canonicalize");
    for (case, _, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .canonicalize(&context)
                            .expect("retained scaling case canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/scaling/canonicalize_with_remapping");
    for (case, _, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .canonicalize_with_remapping(&context)
                            .expect("retained scaling case canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/scaling/canonical_hash");
    for (case, _, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .canonical_hash(&context)
                            .expect("retained scaling case canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/scaling/canonical_eq_equal");
    for (case, renumbered, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(&case.molecule).canonical_eq(black_box(renumbered), &context))
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/scaling/canonical_eq_unequal");
    for (case, _, unequal) in &cases {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(&case.molecule).canonical_eq(black_box(unequal), &context))
        });
    }
    group.finish();
}

fn bench_topology_path_scaling(c: &mut Criterion) {
    let context = CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: ALGORITHM,
    };
    let cases = [8, 16, 32, 64, 77, 128].map(|atom_count| (atom_count, carbon_path(atom_count)));
    let mut group = c.benchmark_group("canonicalize/scaling/topology_path");

    for (atom_count, molecule) in &cases {
        group.bench_function(BenchmarkId::from_parameter(atom_count), |b| {
            b.iter_batched(
                || molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .canonicalize(&context)
                            .expect("topology path canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

fn bench_reframe(c: &mut Criterion) {
    let corpus = reframing_corpus();

    let mut group = c.benchmark_group("reframe/representative_action");
    for case in &corpus {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(&case.molecule).representative_action())
        });
    }
    group.finish();

    let mut group = c.benchmark_group("reframe/reframe_with_action");
    for case in &corpus {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .reframe_with_action()
                            .expect("benchmark corpus reframes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("reframe/reframe");
    for case in &corpus {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| black_box(molecule.reframe().expect("benchmark corpus reframes")),
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let representatives = corpus
        .iter()
        .map(|case| {
            case.molecule
                .clone()
                .reframe()
                .expect("benchmark corpus reframes")
        })
        .collect::<Vec<_>>();
    let mut group = c.benchmark_group("reframe/framed_eq");
    for (case, representative) in corpus.iter().zip(&representatives) {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(&case.molecule).framed_eq(black_box(representative)))
        });
    }
    group.finish();
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
    bench_remapping,
    bench_canonicalize,
    bench_retained_scaling_cases,
    bench_topology_path_scaling,
    bench_reframe,
);
criterion_main!(canonicalize);
