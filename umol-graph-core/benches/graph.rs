//! Graph addition, removal, and pushout benchmarks on paths.
//!
//! Addition and removal exclude input construction and final destruction. Addition
//! covers unique and shared storage; removal starts from shared storage. Pushout
//! includes destruction of its result.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_graph_core::{Correspondence, Graph, GraphCorrespondence, NodeId};

fn graph(node_count: usize) -> Graph {
    let edges: Vec<_> = (1..node_count as u32)
        .map(|node| [node - 1, node])
        .collect();
    Graph::new(node_count, &edges)
}

fn extend_nodes(c: &mut Criterion) {
    let mut group = c.benchmark_group("graph_extend_nodes");
    for size in [64usize, 1024] {
        let original = graph(size);
        for count in [1usize, 16] {
            for shared in [false, true] {
                let ownership = if shared { "shared" } else { "unique" };
                let fixture = format!("path_{size}/batch_{count}/{ownership}");
                group.bench_function(BenchmarkId::new("add_node", &fixture), |b| {
                    b.iter_batched_ref(
                        || {
                            if shared {
                                original.clone()
                            } else {
                                graph(size)
                            }
                        },
                        |graph| {
                            for _ in 0..black_box(count) {
                                black_box(graph.add_node());
                            }
                            black_box(graph);
                        },
                        BatchSize::SmallInput,
                    );
                });
                group.bench_function(BenchmarkId::new("extend_nodes", &fixture), |b| {
                    b.iter_batched_ref(
                        || {
                            if shared {
                                original.clone()
                            } else {
                                graph(size)
                            }
                        },
                        |graph| {
                            let _ = black_box(graph.extend_nodes(black_box(count)));
                            black_box(graph);
                        },
                        BatchSize::SmallInput,
                    );
                });
            }
        }
    }
    group.finish();
}

fn extend_edges(c: &mut Criterion) {
    let mut group = c.benchmark_group("graph_extend_edges");
    for size in [64usize, 1024] {
        let original = graph(size);
        for count in [1usize, 16] {
            let edges: Vec<_> = (0..count)
                .map(|index| [NodeId::from(index), NodeId::from(size - index - 1)])
                .collect();
            for shared in [false, true] {
                let ownership = if shared { "shared" } else { "unique" };
                let fixture = format!("path_{size}/batch_{count}/{ownership}");
                group.bench_function(BenchmarkId::new("add_edge", &fixture), |b| {
                    b.iter_batched_ref(
                        || {
                            if shared {
                                original.clone()
                            } else {
                                graph(size)
                            }
                        },
                        |graph| {
                            for &[a, b] in black_box(&edges) {
                                black_box(graph.add_edge(a, b));
                            }
                            black_box(graph);
                        },
                        BatchSize::SmallInput,
                    );
                });
                group.bench_function(BenchmarkId::new("extend_edges", &fixture), |b| {
                    b.iter_batched_ref(
                        || {
                            if shared {
                                original.clone()
                            } else {
                                graph(size)
                            }
                        },
                        |graph| {
                            let _ = black_box(graph.extend_edges(black_box(&edges)));
                            black_box(graph);
                        },
                        BatchSize::SmallInput,
                    );
                });
            }
        }
    }
    group.finish();
}

fn extend(c: &mut Criterion) {
    let mut group = c.benchmark_group("graph_extend");
    for size in [64usize, 1024] {
        let original = graph(size);
        for count in [1usize, 16] {
            let edges: Vec<_> = (0..count)
                .map(|index| [NodeId::from(size + index - 1), NodeId::from(size + index)])
                .collect();
            for shared in [false, true] {
                let ownership = if shared { "shared" } else { "unique" };
                let fixture = format!("path_{size}/batch_{count}/{ownership}");
                group.bench_function(BenchmarkId::new("add_node/add_edge", &fixture), |b| {
                    b.iter_batched_ref(
                        || {
                            if shared {
                                original.clone()
                            } else {
                                graph(size)
                            }
                        },
                        |graph| {
                            for _ in 0..black_box(count) {
                                black_box(graph.add_node());
                            }
                            for &[a, b] in black_box(&edges) {
                                black_box(graph.add_edge(a, b));
                            }
                            black_box(graph);
                        },
                        BatchSize::SmallInput,
                    );
                });
                group.bench_function(BenchmarkId::new("extend", &fixture), |b| {
                    b.iter_batched_ref(
                        || {
                            if shared {
                                original.clone()
                            } else {
                                graph(size)
                            }
                        },
                        |graph| {
                            let _ = black_box(graph.extend(black_box(count), black_box(&edges)));
                            black_box(graph);
                        },
                        BatchSize::SmallInput,
                    );
                });
            }
        }
    }
    group.finish();
}

fn remove_cascading(c: &mut Criterion) {
    let mut group = c.benchmark_group("graph_remove_cascading");
    for size in [8, 64] {
        let graph = graph(size);
        let removed = [NodeId((size / 2) as u32)];
        group.bench_function(BenchmarkId::new("remove_cascading/path", size), |b| {
            b.iter_batched(
                || graph.clone(),
                |mut graph| {
                    graph.remove_cascading(black_box(&removed), &[]);
                    black_box(graph)
                },
                BatchSize::SmallInput,
            )
        });
        group.bench_function(
            BenchmarkId::new("tracked_remove_cascading/path", size),
            |b| {
                b.iter_batched(
                    || graph.clone(),
                    |mut graph| {
                        let compaction = graph.tracked_remove_cascading(black_box(&removed), &[]);
                        black_box((graph, compaction))
                    },
                    BatchSize::SmallInput,
                )
            },
        );
    }
    group.finish();
}

fn pushout(c: &mut Criterion) {
    let mut group = c.benchmark_group("graph_pushout");
    for size in [8, 64] {
        let graph = graph(size);
        let overlap = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId((size - 1) as u32), NodeId(0))], size, size).unwrap(),
            Correspondence::new(vec![], size - 1, size - 1).unwrap(),
        );
        group.bench_function(BenchmarkId::new("pushout/path_pair", size), |b| {
            b.iter(|| black_box(&graph).pushout(black_box(&graph), black_box(&overlap)))
        });
        group.bench_function(BenchmarkId::new("tracked_pushout/path_pair", size), |b| {
            b.iter(|| black_box(&graph).tracked_pushout(black_box(&graph), black_box(&overlap)))
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    extend_nodes,
    extend_edges,
    extend,
    remove_cascading,
    pushout,
);
criterion_main!(benches);
