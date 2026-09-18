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
//! Whole-relation addition/removal uses fresh sets, excluding setup cloning and final destruction.
//! It includes column changes, payload removal, incidence rebuilding, and compaction construction.

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
            assert_eq!(relations.incident_to_node(node), expected_incidence);
            assert_eq!(relations.coincident_to_node(node, &query_1), Some(id));
            assert_eq!(relations.coincident_to_node(node, &missing), None);
            group.bench_function(BenchmarkId::new("new", &fixture), |b| {
                b.iter_batched(
                    || entries.clone(),
                    |entries| black_box(FixedRelationSet::new(black_box(entries))),
                    BatchSize::LargeInput,
                )
            });
            group.bench_function(BenchmarkId::new("incident_to_node", &fixture), |b| {
                b.iter(|| black_box(relations.incident_to_node(black_box(node))))
            });
            group.bench_function(BenchmarkId::new("coincident_to_node_hit", &fixture), |b| {
                b.iter(|| {
                    black_box(relations.coincident_to_node(black_box(node), black_box(&query_1)))
                })
            });
            group.bench_function(BenchmarkId::new("coincident_to_node_miss", &fixture), |b| {
                b.iter(|| {
                    black_box(relations.coincident_to_node(black_box(node), black_box(&missing)))
                })
            });
            group.bench_function(BenchmarkId::new("permute_participants", &fixture), |b| {
                b.iter_batched_ref(
                    || relations.clone(),
                    |relations| {
                        relations.permute_participants(black_box(id), black_box(&order_1));
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
                    assert_eq!(replaced.incident_to_node(node), incidence);
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
                    assert_eq!(replaced.incident_to_node(node), incidence);
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

fn fixed_rows(c: &mut Criterion) {
    let mut group = c.benchmark_group("relation/fixed_rows");
    for count in [0usize, 1, 64, 1024] {
        for repeated in [false, true] {
            let pattern = if repeated { "repeated" } else { "sparse" };
            let fixture = format!("rows={count}/{pattern}");
            let entries: Vec<([NodeId; 4], usize)> = (0..count)
                .map(|row| {
                    (
                        array::from_fn(|position| {
                            NodeId(if repeated {
                                if position < 2 {
                                    0
                                } else {
                                    row as u32 + 1
                                }
                            } else {
                                (row * 4 + position) as u32
                            })
                        }),
                        row,
                    )
                })
                .collect();
            let relations = FixedRelationSet::new(entries.clone());
            let participants = array::from_fn(|position| {
                NodeId(if repeated {
                    if position < 2 {
                        0
                    } else {
                        count as u32 + 1
                    }
                } else {
                    (count * 4 + position) as u32
                })
            });
            let mut added = relations.clone();
            assert_eq!(added.add(participants, count), RelationId::from(count));
            let mut expected = entries.clone();
            expected.push((participants, count));
            assert_eq!(added.into_entries(), expected);
            group.bench_function(BenchmarkId::new("add", &fixture), |b| {
                b.iter_batched_ref(
                    || relations.clone(),
                    |relations| black_box(relations.add(black_box(participants), black_box(count))),
                    BatchSize::LargeInput,
                )
            });
            let mut removals = vec![("empty", vec![])];
            if count > 0 {
                removals.extend([
                    ("first", vec![RelationId(0)]),
                    ("middle", vec![RelationId::from(count / 2)]),
                    ("last", vec![RelationId::from(count - 1)]),
                    (
                        "alternating",
                        (0..count).step_by(2).map(RelationId::from).collect(),
                    ),
                    ("all", (0..count).map(RelationId::from).collect()),
                ]);
            }
            for (selection, ids) in removals {
                let fixture = format!("{fixture}/selection={selection}");
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| !ids.contains(&RelationId::from(*i)))
                    .map(|(_, entry)| *entry)
                    .collect();
                let mut removed = relations.clone();
                let compaction = removed.tracked_remove(&ids);
                assert_eq!(compaction.source_count(), count);
                assert_eq!(compaction.result_count(), expected.len());
                for node in entries.iter().flat_map(|(row, _)| row) {
                    let incidence: Vec<_> = expected
                        .iter()
                        .enumerate()
                        .filter(|(_, (row, _))| row.contains(node))
                        .map(|(i, _)| RelationId::from(i))
                        .collect();
                    assert_eq!(removed.incident_to_node(*node), incidence);
                }
                assert_eq!(removed.into_entries(), expected);
                group.bench_function(BenchmarkId::new("remove", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.remove(black_box(&ids));
                            black_box(relations.count());
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("tracked_remove", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| black_box(relations.tracked_remove(black_box(&ids))),
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
                assert_eq!(relations.incident_to_node(node), expected_incidence);
                assert_eq!(relations.coincident_to_node(node, &query_1), Some(id));
                assert_eq!(relations.coincident_to_node(node, &missing), None);
                group.bench_function(BenchmarkId::new("new", &fixture), |b| {
                    b.iter_batched(
                        || entries.clone(),
                        |entries| black_box(VarRelationSet::new(black_box(entries))),
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("incident_to_node", &fixture), |b| {
                    b.iter(|| black_box(relations.incident_to_node(black_box(node))))
                });
                group.bench_function(BenchmarkId::new("coincident_to_node_hit", &fixture), |b| {
                    b.iter(|| {
                        black_box(
                            relations.coincident_to_node(black_box(node), black_box(&query_1)),
                        )
                    })
                });
                group.bench_function(BenchmarkId::new("coincident_to_node_miss", &fixture), |b| {
                    b.iter(|| {
                        black_box(
                            relations.coincident_to_node(black_box(node), black_box(&missing)),
                        )
                    })
                });
                group.bench_function(BenchmarkId::new("permute_participants", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.permute_participants(black_box(id), black_box(&order_1));
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
                            assert_eq!(replaced.incident_to_node(*node), incidence);
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
                            assert_eq!(changed.incident_to_node(node), incidence);
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

fn var_rows(c: &mut Criterion) {
    let mut group = c.benchmark_group("relation/var_rows");
    for count in [0usize, 1, 64, 1024] {
        for width in [0usize, 4, 32] {
            for repeated in [false, true] {
                let pattern = if repeated { "repeated" } else { "sparse" };
                let fixture = format!("rows={count}/width={width}/{pattern}");
                let entries: Vec<(Vec<NodeId>, usize)> = (0..count)
                    .map(|row| {
                        (
                            (0..width)
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
                                .collect(),
                            row,
                        )
                    })
                    .collect();
                let relations = VarRelationSet::new(entries.clone());
                let participants: Vec<_> = (0..width)
                    .map(|position| {
                        NodeId(if repeated {
                            if position < width / 2 {
                                0
                            } else {
                                count as u32 + 1
                            }
                        } else {
                            (count * width + position) as u32
                        })
                    })
                    .collect();
                let mut added = relations.clone();
                assert_eq!(added.add(&participants, count), RelationId::from(count));
                let mut expected = entries.clone();
                expected.push((participants.clone(), count));
                assert_eq!(added.into_entries(), expected);
                group.bench_function(BenchmarkId::new("add", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            black_box(relations.add(black_box(&participants), black_box(count)))
                        },
                        BatchSize::LargeInput,
                    )
                });
                let mut removals = vec![("empty", vec![])];
                if count > 0 {
                    removals.extend([
                        ("first", vec![RelationId(0)]),
                        ("middle", vec![RelationId::from(count / 2)]),
                        ("last", vec![RelationId::from(count - 1)]),
                        (
                            "alternating",
                            (0..count).step_by(2).map(RelationId::from).collect(),
                        ),
                        ("all", (0..count).map(RelationId::from).collect()),
                    ]);
                }
                for (selection, ids) in removals {
                    let fixture = format!("{fixture}/selection={selection}");
                    let expected: Vec<_> = entries
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| !ids.contains(&RelationId::from(*i)))
                        .map(|(_, entry)| entry.clone())
                        .collect();
                    let mut removed = relations.clone();
                    let compaction = removed.tracked_remove(&ids);
                    assert_eq!(compaction.source_count(), count);
                    assert_eq!(compaction.result_count(), expected.len());
                    for node in entries
                        .first()
                        .into_iter()
                        .chain(entries.last())
                        .flat_map(|(row, _)| row)
                    {
                        let incidence: Vec<_> = expected
                            .iter()
                            .enumerate()
                            .filter(|(_, (row, _))| row.contains(node))
                            .map(|(i, _)| RelationId::from(i))
                            .collect();
                        assert_eq!(removed.incident_to_node(*node), incidence);
                    }
                    assert_eq!(removed.into_entries(), expected);
                    group.bench_function(BenchmarkId::new("remove", &fixture), |b| {
                        b.iter_batched_ref(
                            || relations.clone(),
                            |relations| {
                                relations.remove(black_box(&ids));
                                black_box(relations.count());
                            },
                            BatchSize::LargeInput,
                        )
                    });
                    group.bench_function(BenchmarkId::new("tracked_remove", &fixture), |b| {
                        b.iter_batched_ref(
                            || relations.clone(),
                            |relations| black_box(relations.tracked_remove(black_box(&ids))),
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
            assert_eq!(relations.incident_to_node(node), expected_incidence);
            assert_eq!(
                relations.coincident_to_node(node, &query_1, &query_2),
                Some(id)
            );
            assert_eq!(relations.coincident_to_node(node, &missing, &query_2), None);
            assert_eq!(relations.incident_to_edge(edge), expected_incidence);
            assert_eq!(
                relations.coincident_to_edge(edge, &query_1, &query_2),
                Some(id)
            );
            assert_eq!(relations.coincident_to_edge(edge, &missing, &query_2), None);
            group.bench_function(BenchmarkId::new("new", &fixture), |b| {
                b.iter_batched(
                    || entries.clone(),
                    |entries| black_box(FixedFixedBirelationSet::new(black_box(entries))),
                    BatchSize::LargeInput,
                )
            });
            group.bench_function(BenchmarkId::new("incident_to_node", &fixture), |b| {
                b.iter(|| black_box(relations.incident_to_node(black_box(node))))
            });
            group.bench_function(BenchmarkId::new("coincident_to_node_hit", &fixture), |b| {
                b.iter(|| {
                    black_box(relations.coincident_to_node(
                        black_box(node),
                        black_box(&query_1),
                        black_box(&query_2),
                    ))
                })
            });
            group.bench_function(BenchmarkId::new("coincident_to_node_miss", &fixture), |b| {
                b.iter(|| {
                    black_box(relations.coincident_to_node(
                        black_box(node),
                        black_box(&missing),
                        black_box(&query_2),
                    ))
                })
            });
            group.bench_function(BenchmarkId::new("incident_to_edge", &fixture), |b| {
                b.iter(|| black_box(relations.incident_to_edge(black_box(edge))))
            });
            group.bench_function(BenchmarkId::new("coincident_to_edge_hit", &fixture), |b| {
                b.iter(|| {
                    black_box(relations.coincident_to_edge(
                        black_box(edge),
                        black_box(&query_1),
                        black_box(&query_2),
                    ))
                })
            });
            group.bench_function(BenchmarkId::new("coincident_to_edge_miss", &fixture), |b| {
                b.iter(|| {
                    black_box(relations.coincident_to_edge(
                        black_box(edge),
                        black_box(&missing),
                        black_box(&query_2),
                    ))
                })
            });
            group.bench_function(BenchmarkId::new("permute_participants_1", &fixture), |b| {
                b.iter_batched_ref(
                    || relations.clone(),
                    |relations| {
                        relations.permute_participants_1(black_box(id), black_box(&order_1));
                        black_box(relations.participants_1(id));
                    },
                    BatchSize::LargeInput,
                )
            });
            group.bench_function(BenchmarkId::new("permute_participants_2", &fixture), |b| {
                b.iter_batched_ref(
                    || relations.clone(),
                    |relations| {
                        relations.permute_participants_2(black_box(id), black_box(&order_2));
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
                        assert_eq!(changed.incident_to_node(node), incidence);
                    }
                    for edge in entries[row].1.into_iter().chain(b) {
                        let incidence: Vec<_> = expected
                            .iter()
                            .enumerate()
                            .filter(|(_, (_, parts, _))| parts.contains(&edge))
                            .map(|(i, _)| RelationId::from(i))
                            .collect();
                        assert_eq!(changed.incident_to_edge(edge), incidence);
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

fn fixed_fixed_rows(c: &mut Criterion) {
    let mut group = c.benchmark_group("relation/fixed_fixed_rows");
    for count in [0usize, 1, 64, 1024] {
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
                    let second = first.map(|node| EdgeId(node.0));
                    (first, second, row)
                })
                .collect();
            let relations = FixedFixedBirelationSet::new(entries.clone());
            let participants = array::from_fn(|position| {
                NodeId(if repeated {
                    if position < 2 {
                        0
                    } else {
                        count as u32 + 1
                    }
                } else {
                    (count * 4 + position) as u32
                })
            });
            let participants_2 = participants.map(|node| EdgeId(node.0));
            let mut added = relations.clone();
            assert_eq!(
                added.add(participants, participants_2, count),
                RelationId::from(count)
            );
            let mut expected = entries.clone();
            expected.push((participants, participants_2, count));
            assert_eq!(added.into_entries(), expected);
            group.bench_function(BenchmarkId::new("add", &fixture), |b| {
                b.iter_batched_ref(
                    || relations.clone(),
                    |relations| {
                        black_box(relations.add(
                            black_box(participants),
                            black_box(participants_2),
                            black_box(count),
                        ))
                    },
                    BatchSize::LargeInput,
                )
            });
            let mut removals = vec![("empty", vec![])];
            if count > 0 {
                removals.extend([
                    ("first", vec![RelationId(0)]),
                    ("middle", vec![RelationId::from(count / 2)]),
                    ("last", vec![RelationId::from(count - 1)]),
                    (
                        "alternating",
                        (0..count).step_by(2).map(RelationId::from).collect(),
                    ),
                    ("all", (0..count).map(RelationId::from).collect()),
                ]);
            }
            for (selection, ids) in removals {
                let fixture = format!("{fixture}/selection={selection}");
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| !ids.contains(&RelationId::from(*i)))
                    .map(|(_, entry)| *entry)
                    .collect();
                let mut removed = relations.clone();
                let compaction = removed.tracked_remove(&ids);
                assert_eq!(compaction.source_count(), count);
                assert_eq!(compaction.result_count(), expected.len());
                for node in entries.iter().flat_map(|(row, _, _)| row) {
                    let incidence: Vec<_> = expected
                        .iter()
                        .enumerate()
                        .filter(|(_, (row, _, _))| row.contains(node))
                        .map(|(i, _)| RelationId::from(i))
                        .collect();
                    assert_eq!(removed.incident_to_node(*node), incidence);
                    assert_eq!(removed.incident_to_edge(EdgeId(node.0)), incidence);
                }
                assert_eq!(removed.into_entries(), expected);
                group.bench_function(BenchmarkId::new("remove", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.remove(black_box(&ids));
                            black_box(relations.count());
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("tracked_remove", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| black_box(relations.tracked_remove(black_box(&ids))),
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
                assert_eq!(relations.incident_to_node(node), expected_incidence);
                assert_eq!(
                    relations.coincident_to_node(node, &query_1, &query_2),
                    Some(id)
                );
                assert_eq!(relations.coincident_to_node(node, &missing, &query_2), None);
                assert_eq!(relations.incident_to_edge(edge), expected_incidence);
                assert_eq!(
                    relations.coincident_to_edge(edge, &query_1, &query_2),
                    Some(id)
                );
                assert_eq!(relations.coincident_to_edge(edge, &missing, &query_2), None);
                group.bench_function(BenchmarkId::new("new", &fixture), |b| {
                    b.iter_batched(
                        || entries.clone(),
                        |entries| black_box(FixedVarBirelationSet::new(black_box(entries))),
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("incident_to_node", &fixture), |b| {
                    b.iter(|| black_box(relations.incident_to_node(black_box(node))))
                });
                group.bench_function(BenchmarkId::new("coincident_to_node_hit", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_to_node(
                            black_box(node),
                            black_box(&query_1),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("coincident_to_node_miss", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_to_node(
                            black_box(node),
                            black_box(&missing),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("incident_to_edge", &fixture), |b| {
                    b.iter(|| black_box(relations.incident_to_edge(black_box(edge))))
                });
                group.bench_function(BenchmarkId::new("coincident_to_edge_hit", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_to_edge(
                            black_box(edge),
                            black_box(&query_1),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("coincident_to_edge_miss", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_to_edge(
                            black_box(edge),
                            black_box(&missing),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("permute_participants_1", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.permute_participants_1(black_box(id), black_box(&order_1));
                            black_box(relations.participants_1(id));
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("permute_participants_2", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.permute_participants_2(black_box(id), black_box(&order_2));
                            black_box(relations.participants_2(id));
                        },
                        BatchSize::LargeInput,
                    )
                });

                for (row_name, row) in [("first", 0), ("middle", count / 2), ("last", count - 1)] {
                    let id = RelationId::from(row);
                    let fixture = format!("{fixture}/row={row_name}");
                    let first = [EdgeId((count * 2) as u32), EdgeId((count * 2 + 1) as u32)];
                    let participant_1 = first[0];
                    let participant_2 = NodeId((count * width) as u32);
                    let position = ParticipantPosition(0);
                    for replacement_len in [0, width / 2, width, width * 2] {
                        let fixture = format!("{fixture}/replacement={replacement_len}");
                        let second: Vec<_> = (0..replacement_len)
                            .map(|position| {
                                NodeId(if repeated && position < replacement_len / 2 {
                                    0
                                } else {
                                    (count * width + position) as u32
                                })
                            })
                            .collect();
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.replace_participants(id, first, &second);
                        expected[row].0 = first;
                        expected[row].1 = second.clone();

                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("replace_participants", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participants(
                                            black_box(id),
                                            black_box(first),
                                            black_box(&second),
                                        );
                                        black_box(relations.participants_1(id));
                                        black_box(relations.participants_2(id));
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    for replacement_len in [0, width / 2, width, width * 2] {
                        let fixture = format!("{fixture}/replacement={replacement_len}");
                        let second: Vec<_> = (0..replacement_len)
                            .map(|position| {
                                NodeId(if repeated && position < replacement_len / 2 {
                                    0
                                } else {
                                    (count * width + position) as u32
                                })
                            })
                            .collect();
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.replace_participants_2(id, &second);
                        expected[row].1 = second.clone();

                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("replace_participants_2", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participants_2(
                                            black_box(id),
                                            black_box(&second),
                                        );
                                        black_box(relations.participants_1(id));
                                        black_box(relations.participants_2(id));
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.replace_participants_1(id, first);
                        expected[row].0 = first;

                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("replace_participants_1", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participants_1(
                                            black_box(id),
                                            black_box(first),
                                        );
                                        black_box(relations.participants_1(id));
                                        black_box(relations.participants_2(id));
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.replace_participant_1(id, position, participant_1);
                        expected[row].0[0] = participant_1;

                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("replace_participant_1", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participant_1(
                                            black_box(id),
                                            black_box(position),
                                            black_box(participant_1),
                                        );
                                        black_box(relations.participants_1(id));
                                        black_box(relations.participants_2(id));
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.replace_participant_2(id, position, participant_2);
                        expected[row].1[0] = participant_2;

                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("replace_participant_2", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participant_2(
                                            black_box(id),
                                            black_box(position),
                                            black_box(participant_2),
                                        );
                                        black_box(relations.participants_1(id));
                                        black_box(relations.participants_2(id));
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.insert_participant_2(id, position, participant_2);
                        expected[row].1.insert(0, participant_2);

                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("insert_participant_2", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.insert_participant_2(
                                            black_box(id),
                                            black_box(position),
                                            black_box(participant_2),
                                        );
                                        black_box(relations.participants_1(id));
                                        black_box(relations.participants_2(id));
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.remove_participant_2(id, position);
                        expected[row].1.remove(0);

                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("remove_participant_2", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.remove_participant_2(
                                            black_box(id),
                                            black_box(position),
                                        );
                                        black_box(relations.participants_1(id));
                                        black_box(relations.participants_2(id));
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                }
            }
        }
    }
    group.finish();
}

fn fixed_var_rows(c: &mut Criterion) {
    let mut group = c.benchmark_group("relation/fixed_var_rows");
    for count in [0usize, 1, 64, 1024] {
        for width in [0usize, 4, 32] {
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
                let first: [EdgeId; 2] = array::from_fn(|position| {
                    EdgeId(if repeated {
                        if position < 1 {
                            0
                        } else {
                            count as u32 + 1
                        }
                    } else {
                        (count * 2 + position) as u32
                    })
                });
                let participants: Vec<_> = (0..width)
                    .map(|position| {
                        NodeId(if repeated {
                            if position < width / 2 {
                                0
                            } else {
                                count as u32 + 1
                            }
                        } else {
                            (count * width + position) as u32
                        })
                    })
                    .collect();
                let mut added = relations.clone();
                assert_eq!(
                    added.add(first, &participants, count),
                    RelationId::from(count)
                );
                let mut expected = entries.clone();
                expected.push((first, participants.clone(), count));
                assert_eq!(added.into_entries(), expected);
                group.bench_function(BenchmarkId::new("add", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            black_box(relations.add(
                                black_box(first),
                                black_box(&participants),
                                black_box(count),
                            ))
                        },
                        BatchSize::LargeInput,
                    )
                });
                let mut removals = vec![("empty", vec![])];
                if count > 0 {
                    removals.extend([
                        ("first", vec![RelationId(0)]),
                        ("middle", vec![RelationId::from(count / 2)]),
                        ("last", vec![RelationId::from(count - 1)]),
                        (
                            "alternating",
                            (0..count).step_by(2).map(RelationId::from).collect(),
                        ),
                        ("all", (0..count).map(RelationId::from).collect()),
                    ]);
                }
                for (selection, ids) in removals {
                    let fixture = format!("{fixture}/selection={selection}");
                    let expected: Vec<_> = entries
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| !ids.contains(&RelationId::from(*i)))
                        .map(|(_, entry)| entry.clone())
                        .collect();
                    let mut removed = relations.clone();
                    let compaction = removed.tracked_remove(&ids);
                    assert_eq!(compaction.source_count(), count);
                    assert_eq!(compaction.result_count(), expected.len());
                    for node in entries
                        .first()
                        .into_iter()
                        .chain(entries.last())
                        .flat_map(|(_, row, _)| row)
                    {
                        let incidence: Vec<_> = expected
                            .iter()
                            .enumerate()
                            .filter(|(_, (_, row, _))| row.contains(node))
                            .map(|(i, _)| RelationId::from(i))
                            .collect();
                        assert_eq!(removed.incident_to_node(*node), incidence);
                    }
                    for edge in entries
                        .first()
                        .into_iter()
                        .chain(entries.last())
                        .flat_map(|(row, _, _)| row)
                    {
                        let incidence: Vec<_> = expected
                            .iter()
                            .enumerate()
                            .filter(|(_, (row, _, _))| row.contains(edge))
                            .map(|(i, _)| RelationId::from(i))
                            .collect();
                        assert_eq!(removed.incident_to_edge(*edge), incidence);
                    }
                    assert_eq!(removed.into_entries(), expected);
                    group.bench_function(BenchmarkId::new("remove", &fixture), |b| {
                        b.iter_batched_ref(
                            || relations.clone(),
                            |relations| {
                                relations.remove(black_box(&ids));
                                black_box(relations.count());
                            },
                            BatchSize::LargeInput,
                        )
                    });
                    group.bench_function(BenchmarkId::new("tracked_remove", &fixture), |b| {
                        b.iter_batched_ref(
                            || relations.clone(),
                            |relations| black_box(relations.tracked_remove(black_box(&ids))),
                            BatchSize::LargeInput,
                        )
                    });
                }
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
                assert_eq!(relations.incident_to_node(node), expected_incidence);
                assert_eq!(
                    relations.coincident_to_node(node, &query_1, &query_2),
                    Some(id)
                );
                assert_eq!(relations.coincident_to_node(node, &missing, &query_2), None);
                assert_eq!(relations.incident_to_edge(edge), expected_incidence);
                assert_eq!(
                    relations.coincident_to_edge(edge, &query_1, &query_2),
                    Some(id)
                );
                assert_eq!(relations.coincident_to_edge(edge, &missing, &query_2), None);
                group.bench_function(BenchmarkId::new("new", &fixture), |b| {
                    b.iter_batched(
                        || entries.clone(),
                        |entries| black_box(VarVarBirelationSet::new(black_box(entries))),
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("incident_to_node", &fixture), |b| {
                    b.iter(|| black_box(relations.incident_to_node(black_box(node))))
                });
                group.bench_function(BenchmarkId::new("coincident_to_node_hit", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_to_node(
                            black_box(node),
                            black_box(&query_1),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("coincident_to_node_miss", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_to_node(
                            black_box(node),
                            black_box(&missing),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("incident_to_edge", &fixture), |b| {
                    b.iter(|| black_box(relations.incident_to_edge(black_box(edge))))
                });
                group.bench_function(BenchmarkId::new("coincident_to_edge_hit", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_to_edge(
                            black_box(edge),
                            black_box(&query_1),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("coincident_to_edge_miss", &fixture), |b| {
                    b.iter(|| {
                        black_box(relations.coincident_to_edge(
                            black_box(edge),
                            black_box(&missing),
                            black_box(&query_2),
                        ))
                    })
                });
                group.bench_function(BenchmarkId::new("permute_participants_1", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.permute_participants_1(black_box(id), black_box(&order_1));
                            black_box(relations.participants_1(id));
                        },
                        BatchSize::LargeInput,
                    )
                });
                group.bench_function(BenchmarkId::new("permute_participants_2", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            relations.permute_participants_2(black_box(id), black_box(&order_2));
                            black_box(relations.participants_2(id));
                        },
                        BatchSize::LargeInput,
                    )
                });

                for (row_name, row) in [("first", 0), ("middle", count / 2), ("last", count - 1)] {
                    let id = RelationId::from(row);
                    let fixture = format!("{fixture}/row={row_name}");
                    let node = NodeId((count * width) as u32);
                    let edge = EdgeId((count * width) as u32);
                    let position = ParticipantPosition(0);
                    for (length_1, length_2) in [
                        (0, 0),
                        (width / 2, width / 2),
                        (width, width),
                        (width * 2, width * 2),
                        (width / 2, width * 2),
                        (width * 2, width / 2),
                    ] {
                        let fixture = format!("{fixture}/replacement={length_1},{length_2}");
                        let first: Vec<_> = (0..length_1)
                            .map(|position| {
                                NodeId(if repeated && position < length_1 / 2 {
                                    0
                                } else {
                                    (count * width + position) as u32
                                })
                            })
                            .collect();
                        let second: Vec<_> = (0..length_2)
                            .map(|position| {
                                EdgeId(if repeated && position < length_2 / 2 {
                                    0
                                } else {
                                    (count * width + position) as u32
                                })
                            })
                            .collect();
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.replace_participants(id, &first, &second);
                        expected[row].0 = first.clone();
                        expected[row].1 = second.clone();
                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("replace_participants", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participants(
                                            black_box(id),
                                            black_box(&first),
                                            black_box(&second),
                                        );
                                        black_box(relations);
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    for length_1 in [0, width / 2, width, width * 2] {
                        let fixture = format!("{fixture}/replacement={length_1}");
                        let first: Vec<_> = (0..length_1)
                            .map(|position| {
                                NodeId(if repeated && position < length_1 / 2 {
                                    0
                                } else {
                                    (count * width + position) as u32
                                })
                            })
                            .collect();
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.replace_participants_1(id, &first);
                        expected[row].0 = first.clone();
                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("replace_participants_1", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participants_1(
                                            black_box(id),
                                            black_box(&first),
                                        );
                                        black_box(relations);
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    for length_2 in [0, width / 2, width, width * 2] {
                        let fixture = format!("{fixture}/replacement={length_2}");
                        let second: Vec<_> = (0..length_2)
                            .map(|position| {
                                EdgeId(if repeated && position < length_2 / 2 {
                                    0
                                } else {
                                    (count * width + position) as u32
                                })
                            })
                            .collect();
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.replace_participants_2(id, &second);
                        expected[row].1 = second.clone();
                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("replace_participants_2", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participants_2(
                                            black_box(id),
                                            black_box(&second),
                                        );
                                        black_box(relations);
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.replace_participant_1(id, position, node);
                        expected[row].0[0] = node;
                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("replace_participant_1", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participant_1(
                                            black_box(id),
                                            black_box(position),
                                            black_box(node),
                                        );
                                        black_box(relations);
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.replace_participant_2(id, position, edge);
                        expected[row].1[0] = edge;
                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("replace_participant_2", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.replace_participant_2(
                                            black_box(id),
                                            black_box(position),
                                            black_box(edge),
                                        );
                                        black_box(relations);
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.insert_participant_1(id, position, node);
                        expected[row].0.insert(0, node);
                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("insert_participant_1", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.insert_participant_1(
                                            black_box(id),
                                            black_box(position),
                                            black_box(node),
                                        );
                                        black_box(relations);
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.insert_participant_2(id, position, edge);
                        expected[row].1.insert(0, edge);
                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("insert_participant_2", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.insert_participant_2(
                                            black_box(id),
                                            black_box(position),
                                            black_box(edge),
                                        );
                                        black_box(relations);
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.remove_participant_1(id, position);
                        expected[row].0.remove(0);
                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("remove_participant_1", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.remove_participant_1(
                                            black_box(id),
                                            black_box(position),
                                        );
                                        black_box(relations);
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                    {
                        let mut changed = relations.clone();
                        let mut expected = entries.clone();
                        changed.remove_participant_2(id, position);
                        expected[row].1.remove(0);
                        assert_eq!(changed.clone().into_entries(), expected);
                        for node in entries[row].0.iter().chain(&expected[row].0) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (first, _, _))| first.contains(node))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_node(*node), incident);
                        }
                        for edge in entries[row].1.iter().chain(&expected[row].1) {
                            let incident: Vec<_> = expected
                                .iter()
                                .enumerate()
                                .filter(|(_, (_, second, _))| second.contains(edge))
                                .map(|(index, _)| RelationId::from(index))
                                .collect();
                            assert_eq!(changed.incident_to_edge(*edge), incident);
                        }
                        group.bench_function(
                            BenchmarkId::new("remove_participant_2", &fixture),
                            |b| {
                                b.iter_batched_ref(
                                    || relations.clone(),
                                    |relations| {
                                        relations.remove_participant_2(
                                            black_box(id),
                                            black_box(position),
                                        );
                                        black_box(relations);
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                }
            }
        }
    }
    group.finish();
}

fn var_var_rows(c: &mut Criterion) {
    let mut group = c.benchmark_group("relation/var_var_rows");
    for count in [0usize, 1, 64, 1024] {
        for (width_1, width_2) in [
            (0usize, 0usize),
            (0, 4),
            (0, 32),
            (4, 0),
            (4, 4),
            (4, 32),
            (32, 0),
            (32, 4),
            (32, 32),
        ] {
            for repeated in [false, true] {
                let pattern = if repeated { "repeated" } else { "sparse" };
                let fixture = format!("rows={count}/widths={width_1},{width_2}/{pattern}");
                let entries: Vec<(Vec<NodeId>, Vec<EdgeId>, usize)> = (0..count)
                    .map(|row| {
                        let first = (0..width_1)
                            .map(|position| {
                                NodeId(if repeated {
                                    if position < width_1 / 2 {
                                        0
                                    } else {
                                        row as u32 + 1
                                    }
                                } else {
                                    (row * width_1 + position) as u32
                                })
                            })
                            .collect();
                        let second = (0..width_2)
                            .map(|position| {
                                EdgeId(if repeated {
                                    if position < width_2 / 2 {
                                        0
                                    } else {
                                        row as u32 + 1
                                    }
                                } else {
                                    (row * width_2 + position) as u32
                                })
                            })
                            .collect();
                        (first, second, row)
                    })
                    .collect();
                let relations = VarVarBirelationSet::new(entries.clone());
                let first: Vec<_> = (0..width_1)
                    .map(|position| {
                        NodeId(if repeated {
                            if position < width_1 / 2 {
                                0
                            } else {
                                count as u32 + 1
                            }
                        } else {
                            (count * width_1 + position) as u32
                        })
                    })
                    .collect();
                let second: Vec<_> = (0..width_2)
                    .map(|position| {
                        EdgeId(if repeated {
                            if position < width_2 / 2 {
                                0
                            } else {
                                count as u32 + 1
                            }
                        } else {
                            (count * width_2 + position) as u32
                        })
                    })
                    .collect();
                let mut added = relations.clone();
                assert_eq!(added.add(&first, &second, count), RelationId::from(count));
                let mut expected = entries.clone();
                expected.push((first.clone(), second.clone(), count));
                assert_eq!(added.into_entries(), expected);
                group.bench_function(BenchmarkId::new("add", &fixture), |b| {
                    b.iter_batched_ref(
                        || relations.clone(),
                        |relations| {
                            black_box(relations.add(
                                black_box(&first),
                                black_box(&second),
                                black_box(count),
                            ))
                        },
                        BatchSize::LargeInput,
                    )
                });
                let mut removals = vec![("empty", vec![])];
                if count > 0 {
                    removals.extend([
                        ("first", vec![RelationId(0)]),
                        ("middle", vec![RelationId::from(count / 2)]),
                        ("last", vec![RelationId::from(count - 1)]),
                        (
                            "alternating",
                            (0..count).step_by(2).map(RelationId::from).collect(),
                        ),
                        ("all", (0..count).map(RelationId::from).collect()),
                    ]);
                }
                for (selection, ids) in removals {
                    let fixture = format!("{fixture}/selection={selection}");
                    let expected: Vec<_> = entries
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| !ids.contains(&RelationId::from(*i)))
                        .map(|(_, entry)| entry.clone())
                        .collect();
                    let mut removed = relations.clone();
                    let compaction = removed.tracked_remove(&ids);
                    assert_eq!(compaction.source_count(), count);
                    assert_eq!(compaction.result_count(), expected.len());
                    for node in entries
                        .first()
                        .into_iter()
                        .chain(entries.last())
                        .flat_map(|(row, _, _)| row)
                    {
                        let incidence: Vec<_> = expected
                            .iter()
                            .enumerate()
                            .filter(|(_, (row, _, _))| row.contains(node))
                            .map(|(i, _)| RelationId::from(i))
                            .collect();
                        assert_eq!(removed.incident_to_node(*node), incidence);
                    }
                    for edge in entries
                        .first()
                        .into_iter()
                        .chain(entries.last())
                        .flat_map(|(_, row, _)| row)
                    {
                        let incidence: Vec<_> = expected
                            .iter()
                            .enumerate()
                            .filter(|(_, (_, row, _))| row.contains(edge))
                            .map(|(i, _)| RelationId::from(i))
                            .collect();
                        assert_eq!(removed.incident_to_edge(*edge), incidence);
                    }
                    assert_eq!(removed.into_entries(), expected);
                    group.bench_function(BenchmarkId::new("remove", &fixture), |b| {
                        b.iter_batched_ref(
                            || relations.clone(),
                            |relations| {
                                relations.remove(black_box(&ids));
                                black_box(relations.count());
                            },
                            BatchSize::LargeInput,
                        )
                    });
                    group.bench_function(BenchmarkId::new("tracked_remove", &fixture), |b| {
                        b.iter_batched_ref(
                            || relations.clone(),
                            |relations| black_box(relations.tracked_remove(black_box(&ids))),
                            BatchSize::LargeInput,
                        )
                    });
                }
            }
        }
    }
    group.finish();
}

criterion_group!(
    benches,
    fixed,
    fixed_rows,
    var,
    var_rows,
    fixed_fixed,
    fixed_fixed_rows,
    fixed_var,
    fixed_var_rows,
    var_var,
    var_var_rows
);
criterion_main!(benches);
