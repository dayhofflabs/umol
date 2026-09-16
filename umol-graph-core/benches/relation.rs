//! Synthetic relation-storage scaling benchmarks.
//!
//! Sparse rows have distinct references. Repeated rows share reference zero, with half of each
//! factor repeating zero and half repeating a row-specific reference. Coincidence hits select
//! the last row; misses keep the shared anchor and change one participant. Incidence measures
//! returning the indexed slice, not traversing its contents.
//! Construction excludes input cloning and output destruction. Permutation excludes cloning
//! and destruction of the working set; each iteration starts in the original frame. Replacement
//! likewise excludes setup cloning and final destruction, but includes rebuilding the whole index.
//! Whole-factor replacement supplies new references (retaining zero in repeated fixtures); local
//! replacement changes the first participant. Both exercise first, middle, and last rows.
//! Variable-row replacement measures empty, shorter, equal-length, and longer final rows.
//! Local insertion/removal starts from a fresh set on every iteration. All variable mutations
//! include offset adjustment and any movement of later participants as well as index rebuilding.

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
            for (position, row) in [("first", 0), ("middle", count / 2), ("last", count - 1)] {
                let id = RelationId(row as u32);
                let fixture = format!("{fixture}/row={position}");
                let participant = NodeId((count * 4) as u32);
                let participants = array::from_fn(|index| {
                    if repeated && index < 2 {
                        NodeId(0)
                    } else {
                        NodeId((count * 4 + index) as u32)
                    }
                });
                let mut replaced = relations.clone();
                replaced.replace_participants(id, participants);
                let mut expected = entries.clone();
                expected[row].0 = participants;
                for node in entries[row].0.into_iter().chain(participants) {
                    let incidence: Vec<_> = expected
                        .iter()
                        .enumerate()
                        .filter(|(_, (parts, _))| parts.contains(&node))
                        .map(|(i, _)| RelationId(i as u32))
                        .collect();
                    assert_eq!(replaced.incident(node), incidence);
                }
                assert_eq!(replaced.into_entries(), expected);
                let mut replaced = relations.clone();
                replaced.replace_participant(id, ParticipantPosition(0), participant);
                let mut expected = entries.clone();
                expected[row].0[0] = participant;
                for node in [entries[row].0[0], participant] {
                    let incidence: Vec<_> = expected
                        .iter()
                        .enumerate()
                        .filter(|(_, (parts, _))| parts.contains(&node))
                        .map(|(i, _)| RelationId(i as u32))
                        .collect();
                    assert_eq!(replaced.incident(node), incidence);
                }
                assert_eq!(replaced.into_entries(), expected);
                group.bench_function(BenchmarkId::new("replace_participants", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.replace_participants(black_box(id), black_box(participants));
                            black_box(relations);
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("replace_participant", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.replace_participant(
                                black_box(id),
                                black_box(ParticipantPosition(0)),
                                black_box(participant),
                            );
                            black_box(relations);
                        },
                        BatchSize::LargeInput,
                    )
                });
            }
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
                for (position, row) in [("first", 0), ("middle", count / 2), ("last", count - 1)] {
                    let id = RelationId(row as u32);
                    let fixture = format!("{fixture}/row={position}");
                    for length in [0, width / 2, width, width * 2] {
                        let participants: Vec<_> = (0..length)
                            .map(|index| {
                                if repeated && index < length / 2 {
                                    NodeId(0)
                                } else {
                                    NodeId((count * width + index) as u32)
                                }
                            })
                            .collect();
                        let mut expected = entries.clone();
                        expected[row].0 = participants.clone();
                        let mut replaced = relations.clone();
                        replaced.replace_participants(id, &participants);
                        for node in entries[row].0.iter().chain(&participants) {
                            let incidence: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (parts, _))| parts.contains(node))
                                .map(|(i, _)| RelationId::from(i))
                                .collect();
                            assert_eq!(replaced.incident(*node), incidence);
                        }
                        assert_eq!(replaced.into_entries(), expected);
                        group.bench_function(
                            BenchmarkId::new(
                                "replace_participants",
                                format!("{fixture}/length={length}"),
                            ),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participants(
                                            black_box(id),
                                            black_box(&participants),
                                        );
                                        black_box(relations);
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    let participant = NodeId((count * width) as u32);
                    let local_position = ParticipantPosition(0);
                    let mut replaced = relations.clone();
                    replaced.replace_participant(id, local_position, participant);
                    let mut inserted = relations.clone();
                    inserted.insert_participant(id, local_position, participant);
                    let mut removed = relations.clone();
                    removed.remove_participant(id, local_position);
                    let mut replace_expected = entries.clone();
                    replace_expected[row].0[0] = participant;
                    let mut insert_expected = entries.clone();
                    insert_expected[row].0.insert(0, participant);
                    let mut remove_expected = entries.clone();
                    remove_expected[row].0.remove(0);
                    for (changed, expected) in [
                        (replaced, replace_expected),
                        (inserted, insert_expected),
                        (removed, remove_expected),
                    ] {
                        for node in entries[row].0.iter().copied().chain([participant]) {
                            let incidence: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (parts, _))| parts.contains(&node))
                                .map(|(i, _)| RelationId::from(i))
                                .collect();
                            assert_eq!(changed.incident(node), incidence);
                        }
                        assert_eq!(changed.into_entries(), expected);
                    }
                    group.bench_function(BenchmarkId::new("replace_participant", &fixture), |b| {
                        b.iter_batched_ref(
                            || relations.clone(),
                            |relations| {
                                relations.replace_participant(
                                    black_box(id),
                                    black_box(local_position),
                                    black_box(participant),
                                );
                                black_box(relations);
                            },
                            BatchSize::LargeInput,
                        )
                    });
                    group.bench_function(BenchmarkId::new("insert_participant", &fixture), |b| {
                        b.iter_batched_ref(
                            || relations.clone(),
                            |relations| {
                                relations.insert_participant(
                                    black_box(id),
                                    black_box(local_position),
                                    black_box(participant),
                                );
                                black_box(relations);
                            },
                            BatchSize::LargeInput,
                        )
                    });
                    group.bench_function(BenchmarkId::new("remove_participant", &fixture), |b| {
                        b.iter_batched_ref(
                            || relations.clone(),
                            |relations| {
                                relations
                                    .remove_participant(black_box(id), black_box(local_position));
                                black_box(relations);
                            },
                            BatchSize::LargeInput,
                        )
                    });
                }
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
            for (position, row) in [("first", 0), ("middle", count / 2), ("last", count - 1)] {
                let id = RelationId::from(row);
                let fixture = format!("{fixture}/row={position}");
                let first: [NodeId; 4] = array::from_fn(|index| {
                    NodeId(if repeated && index < 2 {
                        0
                    } else {
                        (count * 4 + index) as u32
                    })
                });
                let second: [EdgeId; 4] = array::from_fn(|index| {
                    EdgeId(if repeated && index < 2 {
                        0
                    } else {
                        (count * 4 + index) as u32
                    })
                });
                let node = NodeId((count * 4) as u32);
                let edge = EdgeId((count * 4) as u32);
                let mut both = relations.clone();
                both.replace_participants(id, first, second);
                let mut factor_1 = relations.clone();
                factor_1.replace_participants_1(id, first);
                let mut factor_2 = relations.clone();
                factor_2.replace_participants_2(id, second);
                let mut local_1 = relations.clone();
                local_1.replace_participant_1(id, ParticipantPosition(0), node);
                let mut local_2 = relations.clone();
                local_2.replace_participant_2(id, ParticipantPosition(0), edge);
                let mut local_first = entries[row].0;
                local_first[0] = node;
                let mut local_second = entries[row].1;
                local_second[0] = edge;
                for (changed, a, b) in [
                    (both, first, second),
                    (factor_1, first, entries[row].1),
                    (factor_2, entries[row].0, second),
                    (local_1, local_first, entries[row].1),
                    (local_2, entries[row].0, local_second),
                ] {
                    let mut expected = entries.clone();
                    expected[row].0 = a;
                    expected[row].1 = b;
                    for node in entries[row].0.into_iter().chain(a) {
                        let incidence: Vec<_> = expected
                            .iter()
                            .enumerate()
                            .filter(|(_, (parts, _, _))| parts.contains(&node))
                            .map(|(i, _)| RelationId::from(i))
                            .collect();
                        assert_eq!(changed.incident(node), incidence);
                    }
                    for edge in entries[row].1.into_iter().chain(b) {
                        let incidence: Vec<_> = expected
                            .iter()
                            .enumerate()
                            .filter(|(_, (_, parts, _))| parts.contains(&edge))
                            .map(|(i, _)| RelationId::from(i))
                            .collect();
                        assert_eq!(changed.incident_edge(edge), incidence);
                    }
                    assert_eq!(changed.into_entries(), expected);
                }
                group.bench_function(BenchmarkId::new("replace_participants", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.replace_participants(
                                black_box(id),
                                black_box(first),
                                black_box(second),
                            );
                            black_box(relations);
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("replace_participants_1", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.replace_participants_1(black_box(id), black_box(first));
                            black_box(relations);
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("replace_participants_2", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.replace_participants_2(black_box(id), black_box(second));
                            black_box(relations);
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("replace_participant_1", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.replace_participant_1(
                                black_box(id),
                                black_box(ParticipantPosition(0)),
                                black_box(node),
                            );
                            black_box(relations);
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("replace_participant_2", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.replace_participant_2(
                                black_box(id),
                                black_box(ParticipantPosition(0)),
                                black_box(edge),
                            );
                            black_box(relations);
                        },
                        BatchSize::LargeInput,
                    )
                });
            }
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
