//! Molecule canonicalization carrier benchmarks.
//!
//! The raw-topology path labels only atoms in the stored graph and is therefore
//! a performance reference, not an exact semantic alternative to the incidence
//! paths. The exact adapter retains single-role localized endpoints as direct
//! edges and subdivides role-bearing or duplicate incidences. Criterion ids
//! include the measured graph's node and edge counts.

use std::collections::BTreeMap;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_graph_core::{AutomorphismAlgorithm, Correspondence, Graph, NodeId};
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemId, AtomForm, AtomId, BondForm, BondId, ConstitutionColoring,
    DativeBondForm, DativeBondId, Entity, Incidence, IncidenceGraph, IncidenceLevel, Molecule,
    MoleculeColoring, MoleculeCorrespondence, MoleculeEntries, MulticenterBondForm,
    MulticenterBondId, NoncovalentBondForm, NoncovalentBondId, NoncovalentBondKind, StereoAtomForm,
    StereoAtomId, StereoBondForm, StereoBondId, StereoCoset, StereoKind, StereoLigand,
    StereoLigandKind,
};

const ALGORITHM: AutomorphismAlgorithm = AutomorphismAlgorithm::Nauty;

mod canonicalization_cases;

use canonicalization_cases::{corpus, level_name, LEVELS};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum AdapterColor {
    Entity(u64),
    Incidence(Incidence),
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

fn build_adapter(molecule: &Molecule, incidence: &IncidenceGraph) -> (Graph, Vec<AdapterColor>) {
    let source = incidence.graph();
    let entity_colors = incidence_colors(molecule, incidence);
    let mut colors = entity_colors
        .into_iter()
        .map(AdapterColor::Entity)
        .collect::<Vec<_>>();
    let direct_pair_counts = source
        .edge_ids()
        .filter(|&edge| {
            matches!(
                incidence.incidence(edge),
                Incidence::BondEndpoint | Incidence::NoncovalentEndpoint
            )
        })
        .fold(BTreeMap::<[NodeId; 2], usize>::new(), |mut counts, edge| {
            *counts.entry(source.edge_endpoints(edge)).or_default() += 1;
            counts
        });
    let mut edges = Vec::new();
    for edge in source.edge_ids() {
        let endpoints = source.edge_endpoints(edge);
        let direct = matches!(
            incidence.incidence(edge),
            Incidence::BondEndpoint | Incidence::NoncovalentEndpoint
        ) && direct_pair_counts[&endpoints] == 1;
        if direct {
            edges.push([endpoints[0].0, endpoints[1].0]);
            continue;
        }

        let occurrence = colors.len() as u32;
        colors.push(AdapterColor::Incidence(incidence.incidence(edge).clone()));
        edges.push([endpoints[0].0, occurrence]);
        edges.push([occurrence, endpoints[1].0]);
    }

    (Graph::new(colors.len(), &edges), colors)
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

    for level in LEVELS {
        let mut group = c.benchmark_group(format!(
            "canonicalization/incidence_construction/{}",
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

fn bench_incidence_labeling(c: &mut Criterion) {
    let corpus = corpus();

    for level in LEVELS {
        let mut group = c.benchmark_group(format!(
            "canonicalization/incidence_labeling/{}",
            level_name(level)
        ));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(level);
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

fn bench_adapter_construction(c: &mut Criterion) {
    let corpus = corpus();

    for level in LEVELS {
        let mut group = c.benchmark_group(format!(
            "canonicalization/adapter_construction/{}",
            level_name(level)
        ));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(level);
            let (adapter, _) = build_adapter(&case.molecule, &incidence);
            let size = graph_size(adapter.node_count(), adapter.edge_count());
            group.bench_function(BenchmarkId::new(case.name, size), |b| {
                b.iter(|| build_adapter(black_box(&case.molecule), black_box(&incidence)))
            });
        }
        group.finish();
    }
}

fn bench_adapter_labeling(c: &mut Criterion) {
    let corpus = corpus();

    for level in LEVELS {
        let mut group = c.benchmark_group(format!(
            "canonicalization/adapter_labeling/{}",
            level_name(level)
        ));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(level);
            let (graph, colors) = build_adapter(&case.molecule, &incidence);
            let size = graph_size(graph.node_count(), graph.edge_count());
            group.bench_function(BenchmarkId::new(case.name, size), |b| {
                b.iter(|| black_box(graph.automorphisms(|node| &colors[node.index()], ALGORITHM)))
            });
        }
        group.finish();
    }
}

fn bench_remapping(c: &mut Criterion) {
    let corpus = corpus();
    let mut group = c.benchmark_group("canonicalization/remapping");

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
    canonicalization,
    bench_raw_topology_labeling,
    bench_incidence_construction,
    bench_incidence_labeling,
    bench_adapter_construction,
    bench_adapter_labeling,
    bench_remapping,
);
criterion_main!(canonicalization);
