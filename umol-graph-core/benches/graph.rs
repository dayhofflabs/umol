//! Graph and correspondence operations benchmarks.
//!
//! Addition and removal exclude input construction and final destruction. Addition
//! covers unique and shared storage; removal starts from shared storage. Pushout
//! includes destruction of its result.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_graph_core::{Compaction, Correspondence, EdgeId, Graph, GraphCorrespondence, NodeId};

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
                            black_box(graph.extend_nodes(black_box(count)));
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
                            black_box(graph.extend_edges(black_box(&edges)));
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
                            black_box(graph.extend(black_box(count), black_box(&edges)));
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

fn remove(c: &mut Criterion) {
    let mut group = c.benchmark_group("graph_remove");
    for size in [8, 64] {
        let graph = graph(size);
        let middle = (size / 2) as u32;
        group.bench_function(BenchmarkId::new("remove/path", size), |b| {
            b.iter_batched(
                || {
                    (
                        graph.clone(),
                        vec![NodeId(middle)],
                        vec![EdgeId(middle - 1), EdgeId(middle)],
                    )
                },
                |(mut graph, nodes, edges)| {
                    graph.remove(black_box(nodes), black_box(edges));
                    black_box(graph)
                },
                BatchSize::SmallInput,
            )
        });
        group.bench_function(BenchmarkId::new("tracked_remove/path", size), |b| {
            b.iter_batched(
                || {
                    (
                        graph.clone(),
                        vec![NodeId(middle)],
                        vec![EdgeId(middle - 1), EdgeId(middle)],
                    )
                },
                |(mut graph, nodes, edges)| {
                    let compaction = graph.tracked_remove(black_box(nodes), black_box(edges));
                    black_box((graph, compaction))
                },
                BatchSize::SmallInput,
            )
        });
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
                || (graph.clone(), removed.to_vec()),
                |(mut graph, removed)| {
                    graph.remove_cascading(black_box(removed), vec![]);
                    black_box(graph)
                },
                BatchSize::SmallInput,
            )
        });
        group.bench_function(
            BenchmarkId::new("tracked_remove_cascading/path", size),
            |b| {
                b.iter_batched(
                    || (graph.clone(), removed.to_vec()),
                    |(mut graph, removed)| {
                        let compaction = graph.tracked_remove_cascading(black_box(removed), vec![]);
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

fn induced_edges(c: &mut Criterion) {
    let mut group = c.benchmark_group("induced_edges");
    for (name, size, selected, expected) in [
        ("ring_6", 6, 6, 6),
        ("ring_64", 64, 64, 64),
        ("ring_1024", 1024, 1024, 1024),
        ("ring_4096_subset_6", 4096, 6, 5),
    ] {
        let edges: Vec<[u32; 2]> = (0..size).map(|node| [node, (node + 1) % size]).collect();
        let graph = Graph::new(size as usize, &edges);
        let nodes: Vec<NodeId> = (0..selected).rev().map(NodeId).collect();
        assert_eq!(graph.induced_edges(&nodes).count(), expected);
        group.bench_function(name, |b| {
            b.iter(|| black_box(graph.induced_edges(black_box(&nodes)).count()))
        });
    }
    group.finish();
}

fn correspondence_updates(c: &mut Criterion) {
    let mut group = c.benchmark_group("correspondence_updates");
    for size in [256usize, 4096, 65536] {
        let correspondence = Correspondence::<NodeId>::identity(size);
        let compaction = Compaction::new(size, vec![NodeId::from(size / 2)]).unwrap();
        let step = Correspondence::from(&compaction);
        group.bench_function(BenchmarkId::new("compact_right", size), |b| {
            b.iter_batched(
                || correspondence.clone(),
                |mut value| {
                    value.compact_right(black_box(&compaction)).unwrap();
                    black_box(value)
                },
                BatchSize::SmallInput,
            )
        });
        group.bench_function(BenchmarkId::new("compose", size), |b| {
            b.iter(|| black_box(correspondence.compose(black_box(&step)).unwrap()))
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    extend_nodes,
    extend_edges,
    extend,
    remove,
    remove_cascading,
    pushout,
    induced_edges,
    correspondence_updates,
);
criterion_main!(benches);
