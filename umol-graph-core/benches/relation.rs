//! Synthetic relation-storage scaling benchmarks.
//!
//! Sparse rows have distinct references. Repeated rows share reference zero, with half of each
//! factor repeating zero and half repeating a row-specific reference. Coincidence hits select
//! the last row; misses keep the shared anchor and change one participant. Incidence measures
//! returning the indexed slice, not traversing its contents.
//! Construction excludes input cloning and output destruction. Permutation excludes cloning
//! and destruction of the working set; each iteration starts in the original frame.

use std::array;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_graph_core::{
    EdgeId, FixedFixedBirelationSet, FixedRelationSet, FixedVarBirelationSet, NodeId,
    ParticipantPosition, RelationId, VarRelationSet, VarVarBirelationSet,
};

fn fixed(c: &mut Criterion) {
    let mut group = c.benchmark_group("relation/fixed");
    for count in [1usize, 64, 1024] {
        for repeated in [false, true] {
            let pattern = if repeated { "repeated" } else { "sparse" };
            let fixture = format!("rows={count}/{pattern}");
            let entries: Vec<([NodeId; 4], usize)> = (0..count)
                .map(|row| {
                    let first = array::from_fn(|position| {
                        NodeId(if repeated {
                            if position < 2 {
                                0
                            } else {
                                row as u32 + 1
                            }
                        } else {
                            (row * 4 + position) as u32
                        })
                    });
                    (first, row)
                })
                .collect();
            let relations = FixedRelationSet::new(entries.clone());
            let id = RelationId((count - 1) as u32);
            let query_1: Vec<_> = relations.participants(id).iter().rev().copied().collect();
            let order_1: Vec<_> = (0..4)
                .rev()
                .map(|position| ParticipantPosition(position as u32))
                .collect();
            let node = relations.participants(id)[0];
            let mut missing = query_1.clone();
            missing[0] = NodeId(u32::MAX);
            let expected_incidence: Vec<_> = if repeated {
                (0..count).map(|row| RelationId(row as u32)).collect()
            } else {
                vec![id]
            };
            assert_eq!(relations.incident(node), expected_incidence);
            assert_eq!(relations.coincident(node, &query_1), Some(id));
            assert_eq!(relations.coincident(node, &missing), None);
            group.bench_function(BenchmarkId::new("new", &fixture), |b| {
                b.iter_batched(
                    || entries.clone(),
                    |entries| black_box(FixedRelationSet::new(black_box(entries))),
                    BatchSize::LargeInput,
                )
            });
            group.bench_function(BenchmarkId::new("incident", &fixture), |b| {
                b.iter(|| black_box(relations.incident(black_box(node))))
            });
            group.bench_function(BenchmarkId::new("coincident_hit", &fixture), |b| {
                b.iter(|| black_box(relations.coincident(black_box(node), black_box(&query_1))))
            });
            group.bench_function(BenchmarkId::new("coincident_miss", &fixture), |b| {
                b.iter(|| black_box(relations.coincident(black_box(node), black_box(&missing))))
            });
            group.bench_function(BenchmarkId::new("permute_with", &fixture), |b| {
                b.iter_batched_ref(
                    || relations.clone(),
                    |relations| {
                        relations.permute_with(black_box(id), black_box(&order_1));
                        black_box(relations.participants(id));
                    },
                    BatchSize::LargeInput,
                )
            });
        }
    }
    group.finish();
}

fn var(c: &mut Criterion) {
    let mut group = c.benchmark_group("relation/var");
    for count in [1usize, 64, 1024] {
        for width in [4usize, 32] {
            for repeated in [false, true] {
                let pattern = if repeated { "repeated" } else { "sparse" };
                let fixture = format!("rows={count}/width={width}/{pattern}");
                let entries: Vec<(Vec<NodeId>, usize)> = (0..count)
                    .map(|row| {
                        let first = (0..width)
                            .map(|position| {
                                NodeId(if repeated {
                                    if position < width / 2 {
                                        0
                                    } else {
                                        row as u32 + 1
                                    }
                                } else {
                                    (row * width + position) as u32
                                })
                            })
                            .collect();
                        (first, row)
                    })
                    .collect();
                let relations = VarRelationSet::new(entries.clone());
                let id = RelationId((count - 1) as u32);
                let query_1: Vec<_> = relations.participants(id).iter().rev().copied().collect();
                let order_1: Vec<_> = (0..width)
                    .rev()
                    .map(|position| ParticipantPosition(position as u32))
                    .collect();
                let node = relations.participants(id)[0];
                let mut missing = query_1.clone();
                missing[0] = NodeId(u32::MAX);
                let expected_incidence: Vec<_> = if repeated {
                    (0..count).map(|row| RelationId(row as u32)).collect()
                } else {
                    vec![id]
                };
                assert_eq!(relations.incident(node), expected_incidence);
                assert_eq!(relations.coincident(node, &query_1), Some(id));
                assert_eq!(relations.coincident(node, &missing), None);
                group.bench_function(BenchmarkId::new("new", &fixture), |b| {
                    b.iter_batched(
                        || entries.clone(),
                        |entries| black_box(VarRelationSet::new(black_box(entries))),
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("incident", &fixture), |b| {
                    b.iter(|| black_box(relations.incident(black_box(node))))
                });
                group.bench_function(BenchmarkId::new("coincident_hit", &fixture), |b| {
                    b.iter(|| black_box(relations.coincident(black_box(node), black_box(&query_1))))
                });
                group.bench_function(BenchmarkId::new("coincident_miss", &fixture), |b| {
                    b.iter(|| black_box(relations.coincident(black_box(node), black_box(&missing))))
                });
                group.bench_function(BenchmarkId::new("permute_with", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.permute_with(black_box(id), black_box(&order_1));
                            black_box(relations.participants(id));
                        },
                        BatchSize::LargeInput,
                    )
                });
            }
        }
    }
    group.finish();
}

fn fixed_fixed(c: &mut Criterion) {
    let mut group = c.benchmark_group("relation/fixed_fixed");
    for count in [1usize, 64, 1024] {
        for repeated in [false, true] {
            let pattern = if repeated { "repeated" } else { "sparse" };
            let fixture = format!("rows={count}/{pattern}");
            let entries: Vec<([NodeId; 4], [EdgeId; 4], usize)> = (0..count)
                .map(|row| {
                    let first = array::from_fn(|position| {
                        NodeId(if repeated {
                            if position < 2 {
                                0
                            } else {
                                row as u32 + 1
                            }
                        } else {
                            (row * 4 + position) as u32
                        })
                    });
                    let second = array::from_fn(|position| {
                        EdgeId(if repeated {
                            if position < 2 {
                                0
                            } else {
                                row as u32 + 1
                            }
                        } else {
                            (row * 4 + position) as u32
                        })
                    });
                    (first, second, row)
                })
                .collect();
            let relations = FixedFixedBirelationSet::new(entries.clone());
            let id = RelationId((count - 1) as u32);
            let query_1: Vec<_> = relations.participants_1(id).iter().rev().copied().collect();
            let order_1: Vec<_> = (0..4)
                .rev()
                .map(|position| ParticipantPosition(position as u32))
                .collect();
            let query_2: Vec<_> = relations.participants_2(id).iter().rev().copied().collect();
            let order_2: Vec<_> = (0..4)
                .rev()
                .map(|position| ParticipantPosition(position as u32))
                .collect();
            let node = relations.participants_1(id)[0];
            let edge = relations.participants_2(id)[0];
            let mut missing = query_1.clone();
            missing[0] = NodeId(u32::MAX);
            let expected_incidence: Vec<_> = if repeated {
                (0..count).map(|row| RelationId(row as u32)).collect()
            } else {
                vec![id]
            };
            assert_eq!(relations.incident(node), expected_incidence);
            assert_eq!(relations.coincident(node, &query_1, &query_2), Some(id));
            assert_eq!(relations.coincident(node, &missing, &query_2), None);
            assert_eq!(relations.incident_edge(edge), expected_incidence);
            assert_eq!(
                relations.coincident_edge(edge, &query_1, &query_2),
                Some(id)
            );
            assert_eq!(relations.coincident_edge(edge, &missing, &query_2), None);
            group.bench_function(BenchmarkId::new("new", &fixture), |b| {
                b.iter_batched(
                    || entries.clone(),
                    |entries| black_box(FixedFixedBirelationSet::new(black_box(entries))),
                    BatchSize::LargeInput,
                )
            });
            group.bench_function(BenchmarkId::new("incident", &fixture), |b| {
                b.iter(|| black_box(relations.incident(black_box(node))))
            });
            group.bench_function(BenchmarkId::new("coincident_hit", &fixture), |b| {
                b.iter(|| {
                    black_box(relations.coincident(
                        black_box(node),
                        black_box(&query_1),
                        black_box(&query_2),
                    ))
                })
            });
            group.bench_function(BenchmarkId::new("coincident_miss", &fixture), |b| {
                b.iter(|| {
                    black_box(relations.coincident(
                        black_box(node),
                        black_box(&missing),
                        black_box(&query_2),
                    ))
                })
            });
            group.bench_function(BenchmarkId::new("incident_edge", &fixture), |b| {
                b.iter(|| black_box(relations.incident_edge(black_box(edge))))
            });
            group.bench_function(BenchmarkId::new("coincident_edge_hit", &fixture), |b| {
                b.iter(|| {
                    black_box(relations.coincident_edge(
                        black_box(edge),
                        black_box(&query_1),
                        black_box(&query_2),
                    ))
                })
            });
            group.bench_function(BenchmarkId::new("coincident_edge_miss", &fixture), |b| {
                b.iter(|| {
                    black_box(relations.coincident_edge(
                        black_box(edge),
                        black_box(&missing),
                        black_box(&query_2),
                    ))
                })
            });
            group.bench_function(BenchmarkId::new("permute_1_with", &fixture), |b| {
                b.iter_batched_ref(
                    || relations.clone(),
                    |relations| {
                        relations.permute_1_with(black_box(id), black_box(&order_1));
                        black_box(relations.participants_1(id));
                    },
                    BatchSize::LargeInput,
                )
            });
            group.bench_function(BenchmarkId::new("permute_2_with", &fixture), |b| {
                b.iter_batched_ref(
                    || relations.clone(),
                    |relations| {
                        relations.permute_2_with(black_box(id), black_box(&order_2));
                        black_box(relations.participants_2(id));
                    },
                    BatchSize::LargeInput,
                )
            });
        }
    }
    group.finish();
}

fn fixed_var(c: &mut Criterion) {
    let mut group = c.benchmark_group("relation/fixed_var");
    for count in [1usize, 64, 1024] {
        for width in [4usize, 32] {
            for repeated in [false, true] {
                let pattern = if repeated { "repeated" } else { "sparse" };
                let fixture = format!("rows={count}/width={width}/{pattern}");
                let entries: Vec<([EdgeId; 2], Vec<NodeId>, usize)> = (0..count)
                    .map(|row| {
                        let first = array::from_fn(|position| {
                            EdgeId(if repeated {
                                if position < 1 {
                                    0
                                } else {
                                    row as u32 + 1
                                }
                            } else {
                                (row * 2 + position) as u32
                            })
                        });
                        let second = (0..width)
                            .map(|position| {
                                NodeId(if repeated {
                                    if position < width / 2 {
                                        0
                                    } else {
                                        row as u32 + 1
                                    }
                                } else {
                                    (row * width + position) as u32
                                })
                            })
                            .collect();
                        (first, second, row)
                    })
                    .collect();
                let relations = FixedVarBirelationSet::new(entries.clone());
                let id = RelationId((count - 1) as u32);
                let query_1: Vec<_> = relations.participants_1(id).iter().rev().copied().collect();
                let order_1: Vec<_> = (0..2)
                    .rev()
                    .map(|position| ParticipantPosition(position as u32))
                    .collect();
                let query_2: Vec<_> = relations.participants_2(id).iter().rev().copied().collect();
                let order_2: Vec<_> = (0..width)
                    .rev()
                    .map(|position| ParticipantPosition(position as u32))
                    .collect();
                let edge = relations.participants_1(id)[0];
                let node = relations.participants_2(id)[0];
                let mut missing = query_1.clone();
                missing[0] = EdgeId(u32::MAX);
                let expected_incidence: Vec<_> = if repeated {
                    (0..count).map(|row| RelationId(row as u32)).collect()
                } else {
                    vec![id]
                };
                assert_eq!(relations.incident(node), expected_incidence);
                assert_eq!(relations.coincident(node, &query_1, &query_2), Some(id));
                assert_eq!(relations.coincident(node, &missing, &query_2), None);
                assert_eq!(relations.incident_edge(edge), expected_incidence);
                assert_eq!(
                    relations.coincident_edge(edge, &query_1, &query_2),
                    Some(id)
                );
                assert_eq!(relations.coincident_edge(edge, &missing, &query_2), None);
                group.bench_function(BenchmarkId::new("new", &fixture), |b| {
                    b.iter_batched(
                        || entries.clone(),
                        |entries| black_box(FixedVarBirelationSet::new(black_box(entries))),
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("incident", &fixture), |b| {
                    b.iter(|| black_box(relations.incident(black_box(node))))
                });
                group.bench_function(BenchmarkId::new("coincident_hit", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident(
                            black_box(node),
                            black_box(&query_1),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("coincident_miss", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident(
                            black_box(node),
                            black_box(&missing),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("incident_edge", &fixture), |b| {
                    b.iter(|| black_box(relations.incident_edge(black_box(edge))))
                });
                group.bench_function(BenchmarkId::new("coincident_edge_hit", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_edge(
                            black_box(edge),
                            black_box(&query_1),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("coincident_edge_miss", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_edge(
                            black_box(edge),
                            black_box(&missing),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("permute_1_with", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.permute_1_with(black_box(id), black_box(&order_1));
                            black_box(relations.participants_1(id));
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("permute_2_with", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.permute_2_with(black_box(id), black_box(&order_2));
                            black_box(relations.participants_2(id));
                        },
                        BatchSize::LargeInput,
                    )
                });
            }
        }
    }
    group.finish();
}

fn var_var(c: &mut Criterion) {
    let mut group = c.benchmark_group("relation/var_var");
    for count in [1usize, 64, 1024] {
        for width in [4usize, 32] {
            for repeated in [false, true] {
                let pattern = if repeated { "repeated" } else { "sparse" };
                let fixture = format!("rows={count}/width={width}/{pattern}");
                let entries: Vec<(Vec<NodeId>, Vec<EdgeId>, usize)> = (0..count)
                    .map(|row| {
                        let first = (0..width)
                            .map(|position| {
                                NodeId(if repeated {
                                    if position < width / 2 {
                                        0
                                    } else {
                                        row as u32 + 1
                                    }
                                } else {
                                    (row * width + position) as u32
                                })
                            })
                            .collect();
                        let second = (0..width)
                            .map(|position| {
                                EdgeId(if repeated {
                                    if position < width / 2 {
                                        0
                                    } else {
                                        row as u32 + 1
                                    }
                                } else {
                                    (row * width + position) as u32
                                })
                            })
                            .collect();
                        (first, second, row)
                    })
                    .collect();
                let relations = VarVarBirelationSet::new(entries.clone());
                let id = RelationId((count - 1) as u32);
                let query_1: Vec<_> = relations.participants_1(id).iter().rev().copied().collect();
                let order_1: Vec<_> = (0..width)
                    .rev()
                    .map(|position| ParticipantPosition(position as u32))
                    .collect();
                let query_2: Vec<_> = relations.participants_2(id).iter().rev().copied().collect();
                let order_2: Vec<_> = (0..width)
                    .rev()
                    .map(|position| ParticipantPosition(position as u32))
                    .collect();
                let node = relations.participants_1(id)[0];
                let edge = relations.participants_2(id)[0];
                let mut missing = query_1.clone();
                missing[0] = NodeId(u32::MAX);
                let expected_incidence: Vec<_> = if repeated {
                    (0..count).map(|row| RelationId(row as u32)).collect()
                } else {
                    vec![id]
                };
                assert_eq!(relations.incident(node), expected_incidence);
                assert_eq!(relations.coincident(node, &query_1, &query_2), Some(id));
                assert_eq!(relations.coincident(node, &missing, &query_2), None);
                assert_eq!(relations.incident_edge(edge), expected_incidence);
                assert_eq!(
                    relations.coincident_edge(edge, &query_1, &query_2),
                    Some(id)
                );
                assert_eq!(relations.coincident_edge(edge, &missing, &query_2), None);
                group.bench_function(BenchmarkId::new("new", &fixture), |b| {
                    b.iter_batched(
                        || entries.clone(),
                        |entries| black_box(VarVarBirelationSet::new(black_box(entries))),
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("incident", &fixture), |b| {
                    b.iter(|| black_box(relations.incident(black_box(node))))
                });
                group.bench_function(BenchmarkId::new("coincident_hit", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident(
                            black_box(node),
                            black_box(&query_1),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("coincident_miss", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident(
                            black_box(node),
                            black_box(&missing),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("incident_edge", &fixture), |b| {
                    b.iter(|| black_box(relations.incident_edge(black_box(edge))))
                });
                group.bench_function(BenchmarkId::new("coincident_edge_hit", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_edge(
                            black_box(edge),
                            black_box(&query_1),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("coincident_edge_miss", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_edge(
                            black_box(edge),
                            black_box(&missing),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("permute_1_with", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.permute_1_with(black_box(id), black_box(&order_1));
                            black_box(relations.participants_1(id));
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("permute_2_with", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.permute_2_with(black_box(id), black_box(&order_2));
                            black_box(relations.participants_2(id));
                        },
                        BatchSize::LargeInput,
                    )
                });
            }
        }
    }
    group.finish();
}

criterion_group!(benches, fixed, var, fixed_fixed, fixed_var, var_var);
criterion_main!(benches);
