//! Connectivity, traversal, and maximum independent set benchmarks.

#[allow(dead_code)]
mod fixtures;

use std::hint::black_box;
use std::ops::ControlFlow;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_graph_core::{
    BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm, Graph,
    MaximumIndependentSetAlgorithm, NeighborhoodAlgorithm, NodeId,
};

fn connected_components(c: &mut Criterion) {
    let graphs = [
        ("hexagon", fixtures::cycle(6)),
        ("cubane", fixtures::cubane()),
        ("dodecahedron", fixtures::dodecahedron()),
        ("c60", fixtures::fullerene_c60()),
    ];

    let mut group = c.benchmark_group("connected_components");
    for (name, g) in &graphs {
        group.bench_function(*name, |b| {
            b.iter(|| g.enumerate_connected_components(ConnectedComponentsAlgorithm::Bfs));
        });
    }
    group.finish();
}

fn traversal(c: &mut Criterion) {
    let tree_edges: Vec<_> = (1..255).map(|node| [(node - 1) / 2, node]).collect();
    let disconnected_edges: Vec<_> = (0..4)
        .flat_map(|component| {
            (0..6).map(move |node| [6 * component + node, 6 * component + (node + 1) % 6])
        })
        .collect();
    let graphs = [
        ("path_64", fixtures::path(64)),
        ("path_1024", fixtures::path(1024)),
        ("binary_tree_255", Graph::new(255, &tree_edges)),
        ("cycle_64", fixtures::cycle(64)),
        ("four_hexagons", Graph::new(24, &disconnected_edges)),
        ("path_and_isolated", Graph::new(4, &[[0, 1], [1, 2]])),
    ];
    let mut components = c.benchmark_group("traversal/components");
    for (name, graph) in &graphs {
        components.bench_function(*name, |b| {
            b.iter(|| {
                black_box(graph).enumerate_connected_components(ConnectedComponentsAlgorithm::Bfs)
            });
        });
    }
    components.finish();
    let mut neighborhood = c.benchmark_group("traversal/neighborhood");
    for (name, graph) in &graphs {
        neighborhood.bench_function(*name, |b| {
            b.iter(|| {
                black_box(graph).neighborhood(NodeId(0), usize::MAX, NeighborhoodAlgorithm::Bfs)
            });
        });
    }
    neighborhood.finish();
    let mut depth_first = c.benchmark_group("traversal/depth_first");
    for (name, graph) in &graphs {
        depth_first.bench_function(*name, |b| {
            b.iter(|| {
                black_box(graph).visit_depth_first(graph.node_ids(), |event| {
                    black_box(event);
                    ControlFlow::<()>::Continue(())
                })
            });
        });
    }
    depth_first.finish();
    let mut breadth_first = c.benchmark_group("traversal/breadth_first");
    for (name, graph) in &graphs {
        breadth_first.bench_function(*name, |b| {
            b.iter(|| {
                black_box(graph).visit_breadth_first(graph.node_ids(), None, |event| {
                    black_box(event);
                    ControlFlow::<()>::Continue(())
                })
            });
        });
    }
    breadth_first.finish();
    let mut depth_first_events = c.benchmark_group("traversal/depth_first_events");
    for (name, graph) in &graphs {
        depth_first_events.bench_function(*name, |b| {
            b.iter(|| black_box(graph).enumerate_depth_first_events(graph.node_ids()));
        });
    }
    depth_first_events.finish();
    let mut breadth_first_events = c.benchmark_group("traversal/breadth_first_events");
    for (name, graph) in &graphs {
        breadth_first_events.bench_function(*name, |b| {
            b.iter(|| black_box(graph).enumerate_breadth_first_events(graph.node_ids(), None));
        });
    }
    breadth_first_events.finish();
    let mut component_visitor = c.benchmark_group("traversal/component_visitor");
    for (name, graph) in &graphs {
        component_visitor.bench_function(*name, |b| {
            b.iter(|| {
                black_box(graph).visit_connected_components(
                    ConnectedComponentsAlgorithm::Bfs,
                    |component| {
                        black_box(component);
                        ControlFlow::<()>::Continue(())
                    },
                )
            });
        });
    }
    component_visitor.finish();
}

fn biconnected_components(c: &mut Criterion) {
    let graphs = [
        ("hexagon", fixtures::cycle(6)),
        ("naphthalene", fixtures::naphthalene()),
        ("cubane", fixtures::cubane()),
        ("dodecahedron", fixtures::dodecahedron()),
        ("c60", fixtures::fullerene_c60()),
    ];

    let mut group = c.benchmark_group("biconnected_components");
    for (name, g) in &graphs {
        group.bench_function(*name, |b| {
            b.iter(|| g.enumerate_biconnected_components(BiconnectedComponentsAlgorithm::Tarjan));
        });
    }
    group.finish();
}

fn maximum_independent_set(c: &mut Criterion) {
    let graphs = [
        ("path_6", fixtures::path(6)),
        ("hexagon", fixtures::cycle(6)),
        ("cubane", fixtures::cubane()),
        ("petersen", fixtures::petersen()),
        ("dodecahedron", fixtures::dodecahedron()),
    ];

    let mut group = c.benchmark_group("maximum_independent_set");
    for (name, g) in &graphs {
        group.bench_function(*name, |b| {
            b.iter(|| g.maximum_independent_set(MaximumIndependentSetAlgorithm::BranchAndBound));
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    connected_components,
    traversal,
    biconnected_components,
    maximum_independent_set,
);
criterion_main!(benches);
