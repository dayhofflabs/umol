//! Baselines for molecule canonicalization carriers.
//!
//! The raw-topology path and the incidence paths do not yet implement the same
//! semantics: the former labels only the stored atom/bond graph, while the
//! latter includes the entity kinds selected by `IncidenceNodeSelection`.
//! Criterion ids include the measured graph's node and edge counts. The exact
//! compact-versus-incidence comparison belongs to the later canonicalization
//! carrier benchmark once a compact overlay-aware labeling path exists.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomForm, AtomId, BondForm, BondId, ConstitutionColoring, DativeBondForm,
    Entity, IncidenceGraph, IncidenceNodeSelection, Molecule, MoleculeColoring, MoleculeEntries,
    MulticenterBondForm, NoncovalentBondForm, NoncovalentBondKind, StereoAtomForm, StereoBondForm,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
};

const ALGORITHM: AutomorphismAlgorithm = AutomorphismAlgorithm::Nauty;

struct CorpusCase {
    name: &'static str,
    molecule: Molecule,
}

#[derive(Clone, Copy)]
enum Level {
    Topology,
    Constitution,
    Full,
}

impl Level {
    const ALL: [Self; 3] = [Self::Topology, Self::Constitution, Self::Full];

    fn name(self) -> &'static str {
        match self {
            Self::Topology => "topology",
            Self::Constitution => "constitution",
            Self::Full => "full",
        }
    }

    fn selection(self) -> IncidenceNodeSelection {
        match self {
            Self::Topology => IncidenceNodeSelection::topological(),
            Self::Constitution => IncidenceNodeSelection::constitution(),
            Self::Full => IncidenceNodeSelection::full(),
        }
    }
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
    Molecule::from_entries(MoleculeEntries {
        atoms: [
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
        .collect(),
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
        ..Default::default()
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

fn corpus() -> [CorpusCase; 6] {
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
    ]
}

fn graph_size(nodes: usize, edges: usize) -> String {
    format!("n{nodes}_e{edges}")
}

fn incidence_colors(molecule: &Molecule, incidence: &IncidenceGraph) -> Vec<u64> {
    let coloring = ConstitutionColoring::full();
    incidence
        .graph()
        .node_ids()
        .map(|node| coloring.color(molecule, incidence.entity(node)))
        .collect()
}

fn bench_raw_topology_labeling(c: &mut Criterion) {
    let corpus = corpus();
    let coloring = ConstitutionColoring::full();
    let mut group = c.benchmark_group("canonicalization/raw_topology_labeling");

    for case in &corpus {
        let graph = case.molecule.raw_graph();
        let colors: Vec<u64> = graph
            .node_ids()
            .map(|node| coloring.color(&case.molecule, Entity::Atom(AtomId::from(node))))
            .collect();
        let size = graph_size(graph.node_count(), graph.edge_count());
        group.bench_function(BenchmarkId::new(case.name, size), |b| {
            b.iter(|| black_box(graph.automorphisms(|node| colors[node.index()], ALGORITHM)))
        });
    }

    group.finish();
}

fn bench_incidence_construction(c: &mut Criterion) {
    let corpus = corpus();

    for level in Level::ALL {
        let selection = level.selection();
        let mut group = c.benchmark_group(format!(
            "canonicalization/incidence_construction/{}",
            level.name()
        ));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(selection);
            let size = graph_size(
                incidence.graph().node_count(),
                incidence.graph().edge_count(),
            );
            group.bench_function(BenchmarkId::new(case.name, size), |b| {
                b.iter(|| black_box(&case.molecule).incidence_graph(selection))
            });
        }
        group.finish();
    }
}

fn bench_incidence_labeling(c: &mut Criterion) {
    let corpus = corpus();

    for level in Level::ALL {
        let selection = level.selection();
        let mut group = c.benchmark_group(format!(
            "canonicalization/incidence_labeling/{}",
            level.name()
        ));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(selection);
            let graph = incidence.graph();
            let colors = incidence_colors(&case.molecule, &incidence);
            let size = graph_size(graph.node_count(), graph.edge_count());
            group.bench_function(BenchmarkId::new(case.name, size), |b| {
                b.iter(|| black_box(graph.automorphisms(|node| colors[node.index()], ALGORITHM)))
            });
        }
        group.finish();
    }
}

criterion_group!(
    canonicalization,
    bench_raw_topology_labeling,
    bench_incidence_construction,
    bench_incidence_labeling,
);
criterion_main!(canonicalization);
