//! Cycle enumeration, minimum cycle basis, unique ring families, and shortest-cycle benchmarks.

#[allow(dead_code)]
mod fixtures;

use std::ops::ControlFlow;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_graph_core::{
    EdgeId, MinimumCycleBasisAlgorithm, NodeId, RelevantCycleEnumerationAlgorithm,
    ShortestCycleAlgorithm, SimpleCycleEnumerationAlgorithm, UniqueRingFamilyAlgorithm,
};

fn relevant_cycle_enumeration(c: &mut Criterion) {
    let mut group = c.benchmark_group("relevant_cycles");
    for (name, graph) in fixtures::cycle_corpus() {
        group.bench_with_input(BenchmarkId::new("vismara", name), &graph, |b, graph| {
            b.iter(|| {
                graph.enumerate_relevant_cycles(
                    usize::MAX,
                    RelevantCycleEnumerationAlgorithm::Vismara,
                )
            });
        });
    }
    group.finish();
}

fn simple_cycle_enumeration(c: &mut Criterion) {
    let mut bounded = c.benchmark_group("simple_cycles_bounded_8");
    for (name, graph) in fixtures::cycle_corpus() {
        bounded.bench_with_input(BenchmarkId::new("read_tarjan", name), &graph, |b, graph| {
            b.iter(|| {
                graph.enumerate_simple_cycles(8, SimpleCycleEnumerationAlgorithm::ReadTarjan)
            });
        });
    }
    bounded.finish();

    let cases = [
        ("path_6", fixtures::path(6)),
        ("hexagon", fixtures::cycle(6)),
        ("naphthalene", fixtures::naphthalene()),
        ("prismane", fixtures::prismane()),
        ("cubane", fixtures::cubane()),
        ("adamantane", fixtures::adamantane()),
    ];
    let mut unbounded = c.benchmark_group("simple_cycles_unbounded");
    for (name, graph) in cases {
        unbounded.bench_with_input(BenchmarkId::new("read_tarjan", name), &graph, |b, graph| {
            b.iter(|| {
                graph.enumerate_simple_cycles(
                    usize::MAX,
                    SimpleCycleEnumerationAlgorithm::ReadTarjan,
                )
            });
        });
    }
    unbounded.finish();
}

fn minimum_cycle_basis(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimum_cycle_basis");
    for (name, graph) in fixtures::cycle_corpus() {
        group.bench_with_input(BenchmarkId::new("horton", name), &graph, |b, graph| {
            b.iter(|| graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton));
        });
    }
    group.finish();
}

fn unique_ring_families(c: &mut Criterion) {
    let mut decomposition = c.benchmark_group("unique_ring_families/decomposition");
    for (name, graph) in fixtures::cycle_corpus() {
        decomposition.bench_with_input(BenchmarkId::new("kolodzik", name), &graph, |b, graph| {
            b.iter(|| graph.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik));
        });
    }
    decomposition.finish();

    let families = fixtures::cycle_corpus()
        .into_iter()
        .map(|(name, graph)| {
            (
                name,
                graph.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik),
            )
        })
        .collect::<Vec<_>>();
    let mut emission = c.benchmark_group("unique_ring_families/lazy_emission");
    for (name, families) in &families {
        emission.bench_with_input(
            BenchmarkId::new("kolodzik", name),
            families,
            |b, families| {
                b.iter(|| {
                    let mut count = 0;
                    for id in families.ids() {
                        let _: ControlFlow<()> = families.visit_relevant_cycles(id, |_| {
                            count += 1;
                            ControlFlow::Continue(())
                        });
                    }
                    count
                });
            },
        );
    }
    emission.finish();
}

fn shortest_cycle(c: &mut Criterion) {
    let graphs = [
        ("hexagon", fixtures::cycle(6)),
        ("naphthalene", fixtures::naphthalene()),
        ("cubane", fixtures::cubane()),
        ("dodecahedron", fixtures::dodecahedron()),
        ("c60", fixtures::fullerene_c60()),
    ];

    let mut group = c.benchmark_group("shortest_cycle");
    for (name, g) in &graphs {
        let eid = EdgeId(0);
        group.bench_function(format!("{name}/edge"), |b| {
            b.iter(|| g.shortest_cycle_through_edge(eid, ShortestCycleAlgorithm::Bfs));
        });
        let nid = NodeId(0);
        group.bench_function(format!("{name}/node"), |b| {
            b.iter(|| g.shortest_cycle_through_node(nid, ShortestCycleAlgorithm::Bfs));
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    relevant_cycle_enumeration,
    simple_cycle_enumeration,
    minimum_cycle_basis,
    unique_ring_families,
    shortest_cycle,
);
criterion_main!(benches);
