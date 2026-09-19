//! Restoration benchmarks over synthetic storage fixtures.
//!
//! Row reconstruction consumes entries, scatters saved rows by original id, fills survivor
//! slots in order, and rebuilds incidence. Participant reconstruction expands surviving
//! references and rebuilds incidence. Graph reconstruction expands endpoints and rebuilds CSR.
//! Saved rows/edges are supplied in reverse order. Original entries remain independent oracles;
//! endpoint and incidence checks run outside timing, including in Criterion's test mode.
//!
//! Setup cloning, fixture capture, and final result destruction are excluded. Entry extraction,
//! temporary allocations, old storage disposal, and constructor/index work are timed. Graph
//! setup creates either a unique CSR or a clone sharing the retained fixture's CSR. These valid-
//! input baselines omit restoration's panic guards and identity fast paths.
//! Native Graph::restore uses the same ownership setup, includes panic guards and CSR rebuilding,
//! and exercises its identity fast path. Payloads are usize; variable entries are materialized
//! as per-row Vecs by into_entries.
//!
//! Reconstruction can panic on manipulated inputs and is not a contract-compliant substitute
//! for restoration. The comparison measures the cost of guarded full reconstruction to inform
//! whether storage-specific optimization is worthwhile. Compare replacement algorithms under
//! the same no-panic contract. The graph `isolates` case removes a trailing degree-zero node;
//! it does not measure a single bonded-atom removal or establish a typical application workload.
//!
//! Fixed- and variable-relation row restoration also measure retained capacity: setup removes rows from a
//! fresh original set, while the ordinary setup clones the compacted set without its spare
//! capacity. Removal stays outside timing for both native and reconstruction paths. Native
//! restoration reuses payload and offset/fixed columns. Variable restoration repacks a flat
//! participant buffer without per-survivor vectors; saved participants are copied into that
//! buffer through the old buffer. All non-identity native relation operations rebuild incidence.

use std::array;
use std::collections::{BTreeMap, BTreeSet};
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_graph_core::{
    Compaction, EdgeId, FixedFixedBirelationSet, FixedRelationSet, FixedVarBirelationSet, Graph,
    GraphCompaction, Neighbor, NodeId, ParticipantRefs, RelationId, RelationParticipant,
    VarRelationSet, VarVarBirelationSet,
};

fn removal_cases(count: usize) -> Vec<(&'static str, Vec<usize>)> {
    if count == 0 {
        return vec![("identity", vec![])];
    }
    vec![
        ("identity", vec![]),
        ("first", vec![0]),
        ("middle", vec![count / 2]),
        ("last", vec![count - 1]),
        ("sparse", (1..count).step_by(16).collect()),
        ("interleaved", (1..count).step_by(2).collect()),
        ("all", (0..count).collect()),
    ]
}

fn rebuild_rows<E>(
    survivors: Vec<E>,
    compaction: &Compaction<RelationId>,
    removed: Vec<(RelationId, E)>,
) -> Vec<E> {
    let mut slots: Vec<_> = (0..compaction.source_count()).map(|_| None).collect();
    for (id, entry) in removed {
        slots[id.index()] = Some(entry);
    }
    let mut survivors = survivors.into_iter();
    slots
        .into_iter()
        .map(|entry| entry.unwrap_or_else(|| survivors.next().expect("fixture covers every slot")))
        .collect()
}

fn assert_incidence<'a>(
    reference_count: usize,
    rows: impl Iterator<Item = Vec<ParticipantRefs>>,
    node_incidence: impl Fn(NodeId) -> &'a [RelationId],
    edge_incidence: impl Fn(EdgeId) -> &'a [RelationId],
) {
    let mut nodes: BTreeMap<NodeId, Vec<RelationId>> = BTreeMap::new();
    let mut edges: BTreeMap<EdgeId, Vec<RelationId>> = BTreeMap::new();
    for (row, participants) in rows.enumerate() {
        for refs in participants {
            if let Some(node) = refs.node {
                nodes.entry(node).or_default().push(RelationId::from(row));
            }
            if let Some(edge) = refs.edge {
                edges.entry(edge).or_default().push(RelationId::from(row));
            }
        }
    }
    for incidence in nodes.values_mut().chain(edges.values_mut()) {
        incidence.dedup();
    }
    for node in (0..=reference_count).map(NodeId::from) {
        assert_eq!(
            node_incidence(node),
            nodes.get(&node).map_or(&[][..], Vec::as_slice)
        );
    }
    for edge in (0..=reference_count).map(EdgeId::from) {
        assert_eq!(
            edge_incidence(edge),
            edges.get(&edge).map_or(&[][..], Vec::as_slice)
        );
    }
}

fn rebuild_graph(
    graph: &Graph,
    compaction: &GraphCompaction,
    removed: &[(EdgeId, [NodeId; 2])],
) -> Graph {
    let mut endpoints = vec![[0; 2]; compaction.edges().source_count()];
    for edge in graph.edge_ids() {
        let original = compaction.uncompact_edge(edge);
        endpoints[original.index()] = graph
            .edge_endpoints(edge)
            .map(|node| compaction.uncompact_node(node).0);
    }
    for &(edge, nodes) in removed {
        endpoints[edge.index()] = nodes.map(|node| node.0);
    }
    Graph::new(compaction.nodes().source_count(), &endpoints)
}

fn graph(c: &mut Criterion) {
    let mut group = c.benchmark_group("restore/graph");
    for count in [0usize, 1, 64, 1024] {
        for shared_references in [false, true] {
            let pattern = if shared_references {
                "shared"
            } else {
                "sparse"
            };
            let node_count = if count == 0 { 0 } else { count + 2 };
            let endpoints: Vec<_> = (0..count as u32)
                .map(|row| {
                    if shared_references {
                        if row % 4 == 0 {
                            [0, 0]
                        } else {
                            [0, 1]
                        }
                    } else {
                        match row % 4 {
                            0 => [row, row],
                            1 => [row - 1, row],
                            2 => [row - 2, row - 1],
                            _ => [row, row + 1],
                        }
                    }
                })
                .collect();
            let original = Graph::new(node_count, &endpoints);
            let mut cases = removal_cases(node_count);
            if count > 0 {
                cases.push(("isolates", vec![node_count - 1]));
                cases.push(("edges_only", vec![]));
            }
            for (removal, positions) in cases {
                let nodes: Vec<_> = positions.into_iter().map(NodeId::from).collect();
                let edges: Vec<_> = if removal == "isolates" {
                    vec![]
                } else {
                    removal_cases(count)
                        .into_iter()
                        .find(|(name, _)| {
                            *name
                                == if removal == "edges_only" {
                                    "all"
                                } else {
                                    removal
                                }
                        })
                        .unwrap()
                        .1
                        .into_iter()
                        .map(EdgeId::from)
                        .collect()
                };
                let mut compacted = original.clone();
                let compaction = compacted.tracked_remove_cascading(&nodes, &edges);
                let removed: Vec<_> = endpoints
                    .iter()
                    .enumerate()
                    .rev()
                    .filter(|(i, pair)| {
                        edges.contains(&EdgeId::from(*i))
                            || pair.iter().any(|&n| nodes.contains(&NodeId(n)))
                    })
                    .map(|(i, pair)| (EdgeId::from(i), pair.map(NodeId)))
                    .collect();
                let rebuilt = rebuild_graph(&compacted, &compaction, &removed);
                let mut restored = compacted.clone();
                restored.restore(&compaction, &removed);
                for checked in [&rebuilt, &restored] {
                    assert_eq!(checked.node_count(), node_count);
                    assert_eq!(checked.edge_count(), endpoints.len());
                    assert_eq!(
                        checked
                            .edge_ids()
                            .map(|id| checked.edge_endpoints(id))
                            .collect::<Vec<_>>(),
                        endpoints
                            .iter()
                            .map(|pair| pair.map(NodeId))
                            .collect::<Vec<_>>()
                    );
                    let mut expected_neighbors = vec![vec![]; node_count];
                    for (edge, &[a, b]) in endpoints.iter().enumerate() {
                        expected_neighbors[a as usize].push(Neighbor {
                            node: NodeId(b),
                            edge: EdgeId::from(edge),
                        });
                        expected_neighbors[b as usize].push(Neighbor {
                            node: NodeId(a),
                            edge: EdgeId::from(edge),
                        });
                    }
                    for (node, mut expected) in expected_neighbors.into_iter().enumerate() {
                        expected.sort_unstable_by_key(|neighbor| (neighbor.node, neighbor.edge));
                        let mut actual = checked.neighbors(NodeId::from(node)).to_vec();
                        actual.sort_unstable_by_key(|neighbor| (neighbor.node, neighbor.edge));
                        assert_eq!(actual, expected);
                    }
                }
                let compacted_endpoints: Vec<_> = compacted
                    .edge_ids()
                    .map(|edge| compacted.edge_endpoints(edge).map(|node| node.0))
                    .collect();
                for shared_csr in [false, true] {
                    let ownership = if shared_csr { "shared" } else { "unique" };
                    let fixture = format!("edges={count}/{pattern}/{removal}/csr={ownership}");
                    group.bench_function(BenchmarkId::new("rebuild", &fixture), |b| {
                        b.iter_batched(
                            || {
                                if shared_csr {
                                    compacted.clone()
                                } else {
                                    Graph::new(compacted.node_count(), &compacted_endpoints)
                                }
                            },
                            |mut graph| {
                                graph = rebuild_graph(
                                    black_box(&graph),
                                    black_box(&compaction),
                                    black_box(&removed),
                                );
                                black_box(graph)
                            },
                            BatchSize::LargeInput,
                        )
                    });
                    group.bench_function(BenchmarkId::new("restore", &fixture), |b| {
                        b.iter_batched(
                            || {
                                if shared_csr {
                                    compacted.clone()
                                } else {
                                    Graph::new(compacted.node_count(), &compacted_endpoints)
                                }
                            },
                            |mut graph| {
                                graph.restore(black_box(&compaction), black_box(&removed));
                                black_box(graph)
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

fn fixed(c: &mut Criterion) {
    let mut group = c.benchmark_group("restore/fixed");
    for count in [0usize, 1, 64, 1024] {
        let width = 4usize;
        for shared in [false, true] {
            let pattern = if shared { "shared" } else { "sparse" };
            let fixture = format!("rows={count}/width={width}/{pattern}");
            let entries: Vec<_> = (0..count)
                .map(|row| {
                    let first = array::from_fn::<_, 4, _>(|position| {
                        NodeId::from(
                            2 * if shared {
                                if position < 2 {
                                    0
                                } else {
                                    row + 1
                                }
                            } else {
                                row * 4 + position
                            },
                        )
                    });
                    (first, row)
                })
                .collect();
            let original = FixedRelationSet::new(entries.clone());
            let check = |rebuilt: &FixedRelationSet<NodeId, usize, 4>, expected: &[_]| {
                assert_eq!(rebuilt.clone().into_entries(), expected);
                assert_incidence(
                    2 * count * width,
                    expected
                        .iter()
                        .map(|(first, _)| first.iter().map(|p| p.refs()).collect()),
                    |node| rebuilt.incident_to_node(node),
                    |edge| rebuilt.incident_to_edge(edge),
                );
            };
            for (removal, positions) in removal_cases(count) {
                let ids: Vec<_> = positions.into_iter().map(RelationId::from).collect();
                let removed: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .rev()
                    .filter(|(row, _)| ids.contains(&RelationId::from(*row)))
                    .map(|(row, entry)| (RelationId::from(row), *entry))
                    .collect();
                let mut compacted = original.clone();
                let compaction = compacted.tracked_remove(&ids);
                let rebuilt = FixedRelationSet::new(rebuild_rows(
                    compacted.clone().into_entries(),
                    &compaction,
                    removed.clone(),
                ));
                check(&rebuilt, &entries);
                let saved: Vec<_> = removed
                    .iter()
                    .map(|(id, (row, data))| (*id, *row, *data))
                    .collect();
                for retained in [false, true] {
                    let setup = || {
                        if retained {
                            let mut relations = original.clone();
                            relations.remove(&ids);
                            relations
                        } else {
                            compacted.clone()
                        }
                    };
                    let mut restored = setup();
                    restored.restore(&compaction, saved.clone());
                    check(&restored, &entries);
                    let rebuild_name = if retained {
                        "rows_rebuild_retained"
                    } else {
                        "rows_rebuild"
                    };
                    let restore_name = if retained {
                        "rows_restore_retained"
                    } else {
                        "rows_restore"
                    };
                    group.bench_function(
                        BenchmarkId::new(rebuild_name, format!("{fixture}/{removal}")),
                        |b| {
                            b.iter_batched(
                                || (setup(), removed.clone()),
                                |(relations, removed)| {
                                    black_box(FixedRelationSet::new(rebuild_rows(
                                        black_box(relations).into_entries(),
                                        black_box(&compaction),
                                        black_box(removed),
                                    )))
                                },
                                BatchSize::LargeInput,
                            )
                        },
                    );
                    group.bench_function(
                        BenchmarkId::new(restore_name, format!("{fixture}/{removal}")),
                        |b| {
                            b.iter_batched(
                                || (setup(), saved.clone()),
                                |(mut relations, removed)| {
                                    relations.restore(black_box(&compaction), black_box(removed));
                                    black_box(relations)
                                },
                                BatchSize::LargeInput,
                            )
                        },
                    );
                }
            }
            let reference_count = 2 * count * width;
            for (removal, positions) in removal_cases(reference_count) {
                let removed: BTreeSet<_> = positions.iter().copied().collect();
                let compaction = GraphCompaction::new(
                    Compaction::new(
                        reference_count,
                        positions.iter().copied().map(NodeId::from).collect(),
                    )
                    .unwrap(),
                    Compaction::new(
                        reference_count,
                        positions.iter().copied().map(EdgeId::from).collect(),
                    )
                    .unwrap(),
                );
                let expected: Vec<_> = entries
                    .iter()
                    .filter(|entry| {
                        entry.0.iter().map(|p| p.refs()).all(|refs| {
                            refs.node.is_none_or(|id| !removed.contains(&id.index()))
                                && refs.edge.is_none_or(|id| !removed.contains(&id.index()))
                        })
                    })
                    .cloned()
                    .collect();
                let compacted = original.compact(&compaction);
                let rebuild = |relations: FixedRelationSet<NodeId, usize, 4>| {
                    FixedRelationSet::new(
                        relations
                            .into_entries()
                            .into_iter()
                            .map(|(first, data)| {
                                let first = first.map(|p| p.uncompact(&compaction));
                                (first, data)
                            })
                            .collect(),
                    )
                };
                check(&rebuild(compacted.clone()), &expected);
                let mut restored = compacted.clone();
                restored.restore_participants(&compaction);
                check(&restored, &expected);
                group.bench_function(
                    BenchmarkId::new("participants_restore", format!("{fixture}/{removal}")),
                    |b| {
                        b.iter_batched(
                            || compacted.clone(),
                            |mut relations| {
                                relations.restore_participants(black_box(&compaction));
                                black_box(relations)
                            },
                            BatchSize::LargeInput,
                        )
                    },
                );
                group.bench_function(
                    BenchmarkId::new("participants_rebuild", format!("{fixture}/{removal}")),
                    |b| {
                        b.iter_batched(
                            || compacted.clone(),
                            |relations| black_box(rebuild(black_box(relations))),
                            BatchSize::LargeInput,
                        )
                    },
                );
            }
        }
    }
    group.finish();
}

fn var(c: &mut Criterion) {
    let mut group = c.benchmark_group("restore/var");
    for count in [0usize, 1, 64, 1024] {
        for width in [4usize, 32] {
            for shared in [false, true] {
                let pattern = if shared { "shared" } else { "sparse" };
                let fixture = format!("rows={count}/width={width}/{pattern}");
                let entries: Vec<_> = (0..count)
                    .map(|row| {
                        let first = (0..(if row % 8 == 7 { 0 } else { width - row % 3 }))
                            .map(|position| {
                                EdgeId::from(
                                    2 * if shared {
                                        if position < width / 2 {
                                            0
                                        } else {
                                            row + 1
                                        }
                                    } else {
                                        row * width + position
                                    },
                                )
                            })
                            .collect::<Vec<_>>();
                        (first, row)
                    })
                    .collect();
                let original = VarRelationSet::new(entries.clone());
                let check = |rebuilt: &VarRelationSet<EdgeId, usize>, expected: &[_]| {
                    assert_eq!(rebuilt.clone().into_entries(), expected);
                    assert_incidence(
                        2 * count * width,
                        expected
                            .iter()
                            .map(|(first, _)| first.iter().map(|p| p.refs()).collect()),
                        |node| rebuilt.incident_to_node(node),
                        |edge| rebuilt.incident_to_edge(edge),
                    );
                };
                for (removal, positions) in removal_cases(count) {
                    let ids: Vec<_> = positions.into_iter().map(RelationId::from).collect();
                    let removed: Vec<_> = entries
                        .iter()
                        .enumerate()
                        .rev()
                        .filter(|(row, _)| ids.contains(&RelationId::from(*row)))
                        .map(|(row, entry)| (RelationId::from(row), entry.clone()))
                        .collect();
                    let mut compacted = original.clone();
                    let compaction = compacted.tracked_remove(&ids);
                    let rebuilt = VarRelationSet::new(rebuild_rows(
                        compacted.clone().into_entries(),
                        &compaction,
                        removed.clone(),
                    ));
                    check(&rebuilt, &entries);
                    let saved: Vec<_> = removed
                        .iter()
                        .map(|(id, (row, data))| (*id, row.clone(), *data))
                        .collect();
                    for retained in [false, true] {
                        let setup = || {
                            if retained {
                                let mut relations = original.clone();
                                relations.remove(&ids);
                                relations
                            } else {
                                compacted.clone()
                            }
                        };
                        let mut restored = setup();
                        restored.restore(&compaction, saved.clone());
                        check(&restored, &entries);
                        let rebuild_name = if retained {
                            "rows_rebuild_retained"
                        } else {
                            "rows_rebuild"
                        };
                        let restore_name = if retained {
                            "rows_restore_retained"
                        } else {
                            "rows_restore"
                        };
                        group.bench_function(
                            BenchmarkId::new(rebuild_name, format!("{fixture}/{removal}")),
                            |b| {
                                b.iter_batched(
                                    || (setup(), removed.clone()),
                                    |(relations, removed)| {
                                        black_box(VarRelationSet::new(rebuild_rows(
                                            black_box(relations).into_entries(),
                                            black_box(&compaction),
                                            black_box(removed),
                                        )))
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                        group.bench_function(
                            BenchmarkId::new(restore_name, format!("{fixture}/{removal}")),
                            |b| {
                                b.iter_batched(
                                    || (setup(), saved.clone()),
                                    |(mut relations, removed)| {
                                        relations
                                            .restore(black_box(&compaction), black_box(removed));
                                        black_box(relations)
                                    },
                                    BatchSize::LargeInput,
                                )
                            },
                        );
                    }
                }
                let reference_count = 2 * count * width;
                for (removal, positions) in removal_cases(reference_count) {
                    let removed: BTreeSet<_> = positions.iter().copied().collect();
                    let compaction = GraphCompaction::new(
                        Compaction::new(
                            reference_count,
                            positions.iter().copied().map(NodeId::from).collect(),
                        )
                        .unwrap(),
                        Compaction::new(
                            reference_count,
                            positions.iter().copied().map(EdgeId::from).collect(),
                        )
                        .unwrap(),
                    );
                    let expected: Vec<_> = entries
                        .iter()
                        .filter(|entry| {
                            entry.0.iter().map(|p| p.refs()).all(|refs| {
                                refs.node.is_none_or(|id| !removed.contains(&id.index()))
                                    && refs.edge.is_none_or(|id| !removed.contains(&id.index()))
                            })
                        })
                        .cloned()
                        .collect();
                    let compacted = original.compact(&compaction);
                    let rebuild = |relations: VarRelationSet<EdgeId, usize>| {
                        VarRelationSet::new(
                            relations
                                .into_entries()
                                .into_iter()
                                .map(|(first, data)| {
                                    let first = first
                                        .into_iter()
                                        .map(|p| p.uncompact(&compaction))
                                        .collect();
                                    (first, data)
                                })
                                .collect(),
                        )
                    };
                    check(&rebuild(compacted.clone()), &expected);
                    let mut restored = compacted.clone();
                    restored.restore_participants(&compaction);
                    check(&restored, &expected);
                    group.bench_function(
                        BenchmarkId::new("participants_restore", format!("{fixture}/{removal}")),
                        |b| {
                            b.iter_batched(
                                || compacted.clone(),
                                |mut relations| {
                                    relations.restore_participants(black_box(&compaction));
                                    black_box(relations)
                                },
                                BatchSize::LargeInput,
                            )
                        },
                    );
                    group.bench_function(
                        BenchmarkId::new("participants_rebuild", format!("{fixture}/{removal}")),
                        |b| {
                            b.iter_batched(
                                || compacted.clone(),
                                |relations| black_box(rebuild(black_box(relations))),
                                BatchSize::LargeInput,
                            )
                        },
                    );
                }
            }
        }
    }
    group.finish();
}

fn fixed_fixed(c: &mut Criterion) {
    let mut group = c.benchmark_group("restore/fixed_fixed");
    for count in [0usize, 1, 64, 1024] {
        let width = 4usize;
        for shared in [false, true] {
            let pattern = if shared { "shared" } else { "sparse" };
            let fixture = format!("rows={count}/width={width}/{pattern}");
            let entries: Vec<_> = (0..count)
                .map(|row| {
                    let first = array::from_fn::<_, 4, _>(|position| {
                        NodeId::from(
                            2 * if shared {
                                if position < 2 {
                                    0
                                } else {
                                    row + 1
                                }
                            } else {
                                row * 4 + position
                            },
                        )
                    });
                    let second = array::from_fn::<_, 4, _>(|position| {
                        EdgeId::from(
                            2 * if shared {
                                if position < 2 {
                                    0
                                } else {
                                    row + 1
                                }
                            } else {
                                row * 4 + position
                            },
                        )
                    });
                    (first, second, row)
                })
                .collect();
            let original = FixedFixedBirelationSet::new(entries.clone());
            let check = |rebuilt: &FixedFixedBirelationSet<NodeId, 4, EdgeId, 4, usize>,
                         expected: &[_]| {
                assert_eq!(rebuilt.clone().into_entries(), expected);
                assert_incidence(
                    2 * count * width,
                    expected.iter().map(|(first, second, _)| {
                        first
                            .iter()
                            .map(|p| p.refs())
                            .chain(second.iter().map(|p| p.refs()))
                            .collect()
                    }),
                    |node| rebuilt.incident_to_node(node),
                    |edge| rebuilt.incident_to_edge(edge),
                );
            };
            for (removal, positions) in removal_cases(count) {
                let ids: Vec<_> = positions.into_iter().map(RelationId::from).collect();
                let removed: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .rev()
                    .filter(|(row, _)| ids.contains(&RelationId::from(*row)))
                    .map(|(row, entry)| (RelationId::from(row), *entry))
                    .collect();
                let mut compacted = original.clone();
                let compaction = compacted.tracked_remove(&ids);
                let rebuilt = FixedFixedBirelationSet::new(rebuild_rows(
                    compacted.clone().into_entries(),
                    &compaction,
                    removed.clone(),
                ));
                check(&rebuilt, &entries);
                group.bench_function(
                    BenchmarkId::new("rows_rebuild", format!("{fixture}/{removal}")),
                    |b| {
                        b.iter_batched(
                            || (compacted.clone(), removed.clone()),
                            |(relations, removed)| {
                                black_box(FixedFixedBirelationSet::new(rebuild_rows(
                                    black_box(relations).into_entries(),
                                    black_box(&compaction),
                                    black_box(removed),
                                )))
                            },
                            BatchSize::LargeInput,
                        )
                    },
                );
            }
            let reference_count = 2 * count * width;
            for (removal, positions) in removal_cases(reference_count) {
                let removed: BTreeSet<_> = positions.iter().copied().collect();
                let compaction = GraphCompaction::new(
                    Compaction::new(
                        reference_count,
                        positions.iter().copied().map(NodeId::from).collect(),
                    )
                    .unwrap(),
                    Compaction::new(
                        reference_count,
                        positions.iter().copied().map(EdgeId::from).collect(),
                    )
                    .unwrap(),
                );
                let expected: Vec<_> = entries
                    .iter()
                    .filter(|entry| {
                        entry
                            .0
                            .iter()
                            .map(|p| p.refs())
                            .chain(entry.1.iter().map(|p| p.refs()))
                            .all(|refs| {
                                refs.node.is_none_or(|id| !removed.contains(&id.index()))
                                    && refs.edge.is_none_or(|id| !removed.contains(&id.index()))
                            })
                    })
                    .cloned()
                    .collect();
                let compacted = original.compact(&compaction);
                let rebuild = |relations: FixedFixedBirelationSet<NodeId, 4, EdgeId, 4, usize>| {
                    FixedFixedBirelationSet::new(
                        relations
                            .into_entries()
                            .into_iter()
                            .map(|(first, second, data)| {
                                let first = first.map(|p| p.uncompact(&compaction));
                                let second = second.map(|p| p.uncompact(&compaction));
                                (first, second, data)
                            })
                            .collect(),
                    )
                };
                check(&rebuild(compacted.clone()), &expected);
                group.bench_function(
                    BenchmarkId::new("participants_rebuild", format!("{fixture}/{removal}")),
                    |b| {
                        b.iter_batched(
                            || compacted.clone(),
                            |relations| black_box(rebuild(black_box(relations))),
                            BatchSize::LargeInput,
                        )
                    },
                );
            }
        }
    }
    group.finish();
}

fn fixed_var(c: &mut Criterion) {
    let mut group = c.benchmark_group("restore/fixed_var");
    for count in [0usize, 1, 64, 1024] {
        for width in [4usize, 32] {
            for shared in [false, true] {
                let pattern = if shared { "shared" } else { "sparse" };
                let fixture = format!("rows={count}/width={width}/{pattern}");
                let entries: Vec<_> = (0..count)
                    .map(|row| {
                        let first = array::from_fn::<_, 4, _>(|position| {
                            NodeId::from(
                                2 * if shared {
                                    if position < 2 {
                                        0
                                    } else {
                                        row + 1
                                    }
                                } else {
                                    row * 4 + position
                                },
                            )
                        });
                        let second = (0..(if row % 8 == 6 {
                            0
                        } else {
                            width - (row + 1) % 3
                        }))
                            .map(|position| {
                                EdgeId::from(
                                    2 * if shared {
                                        if position < width / 2 {
                                            0
                                        } else {
                                            row + 1
                                        }
                                    } else {
                                        row * width + position
                                    },
                                )
                            })
                            .collect::<Vec<_>>();
                        (first, second, row)
                    })
                    .collect();
                let original = FixedVarBirelationSet::new(entries.clone());
                let check = |rebuilt: &FixedVarBirelationSet<NodeId, 4, EdgeId, usize>,
                             expected: &[_]| {
                    assert_eq!(rebuilt.clone().into_entries(), expected);
                    assert_incidence(
                        2 * count * width,
                        expected.iter().map(|(first, second, _)| {
                            first
                                .iter()
                                .map(|p| p.refs())
                                .chain(second.iter().map(|p| p.refs()))
                                .collect()
                        }),
                        |node| rebuilt.incident_to_node(node),
                        |edge| rebuilt.incident_to_edge(edge),
                    );
                };
                for (removal, positions) in removal_cases(count) {
                    let ids: Vec<_> = positions.into_iter().map(RelationId::from).collect();
                    let removed: Vec<_> = entries
                        .iter()
                        .enumerate()
                        .rev()
                        .filter(|(row, _)| ids.contains(&RelationId::from(*row)))
                        .map(|(row, entry)| (RelationId::from(row), entry.clone()))
                        .collect();
                    let mut compacted = original.clone();
                    let compaction = compacted.tracked_remove(&ids);
                    let rebuilt = FixedVarBirelationSet::new(rebuild_rows(
                        compacted.clone().into_entries(),
                        &compaction,
                        removed.clone(),
                    ));
                    check(&rebuilt, &entries);
                    group.bench_function(
                        BenchmarkId::new("rows_rebuild", format!("{fixture}/{removal}")),
                        |b| {
                            b.iter_batched(
                                || (compacted.clone(), removed.clone()),
                                |(relations, removed)| {
                                    black_box(FixedVarBirelationSet::new(rebuild_rows(
                                        black_box(relations).into_entries(),
                                        black_box(&compaction),
                                        black_box(removed),
                                    )))
                                },
                                BatchSize::LargeInput,
                            )
                        },
                    );
                }
                let reference_count = 2 * count * width;
                for (removal, positions) in removal_cases(reference_count) {
                    let removed: BTreeSet<_> = positions.iter().copied().collect();
                    let compaction = GraphCompaction::new(
                        Compaction::new(
                            reference_count,
                            positions.iter().copied().map(NodeId::from).collect(),
                        )
                        .unwrap(),
                        Compaction::new(
                            reference_count,
                            positions.iter().copied().map(EdgeId::from).collect(),
                        )
                        .unwrap(),
                    );
                    let expected: Vec<_> = entries
                        .iter()
                        .filter(|entry| {
                            entry
                                .0
                                .iter()
                                .map(|p| p.refs())
                                .chain(entry.1.iter().map(|p| p.refs()))
                                .all(|refs| {
                                    refs.node.is_none_or(|id| !removed.contains(&id.index()))
                                        && refs.edge.is_none_or(|id| !removed.contains(&id.index()))
                                })
                        })
                        .cloned()
                        .collect();
                    let compacted = original.compact(&compaction);
                    let rebuild = |relations: FixedVarBirelationSet<NodeId, 4, EdgeId, usize>| {
                        FixedVarBirelationSet::new(
                            relations
                                .into_entries()
                                .into_iter()
                                .map(|(first, second, data)| {
                                    let first = first.map(|p| p.uncompact(&compaction));
                                    let second = second
                                        .into_iter()
                                        .map(|p| p.uncompact(&compaction))
                                        .collect();
                                    (first, second, data)
                                })
                                .collect(),
                        )
                    };
                    check(&rebuild(compacted.clone()), &expected);
                    group.bench_function(
                        BenchmarkId::new("participants_rebuild", format!("{fixture}/{removal}")),
                        |b| {
                            b.iter_batched(
                                || compacted.clone(),
                                |relations| black_box(rebuild(black_box(relations))),
                                BatchSize::LargeInput,
                            )
                        },
                    );
                }
            }
        }
    }
    group.finish();
}

fn var_var(c: &mut Criterion) {
    let mut group = c.benchmark_group("restore/var_var");
    for count in [0usize, 1, 64, 1024] {
        for width in [4usize, 32] {
            for shared in [false, true] {
                let pattern = if shared { "shared" } else { "sparse" };
                let fixture = format!("rows={count}/width={width}/{pattern}");
                let entries: Vec<_> = (0..count)
                    .map(|row| {
                        let first = (0..(if row % 8 == 7 { 0 } else { width - row % 3 }))
                            .map(|position| {
                                NodeId::from(
                                    2 * if shared {
                                        if position < width / 2 {
                                            0
                                        } else {
                                            row + 1
                                        }
                                    } else {
                                        row * width + position
                                    },
                                )
                            })
                            .collect::<Vec<_>>();
                        let second = (0..(if row % 8 == 6 {
                            0
                        } else {
                            width - (row + 1) % 3
                        }))
                            .map(|position| {
                                EdgeId::from(
                                    2 * if shared {
                                        if position < width / 2 {
                                            0
                                        } else {
                                            row + 1
                                        }
                                    } else {
                                        row * width + position
                                    },
                                )
                            })
                            .collect::<Vec<_>>();
                        (first, second, row)
                    })
                    .collect();
                let original = VarVarBirelationSet::new(entries.clone());
                let check = |rebuilt: &VarVarBirelationSet<NodeId, EdgeId, usize>,
                             expected: &[_]| {
                    assert_eq!(rebuilt.clone().into_entries(), expected);
                    assert_incidence(
                        2 * count * width,
                        expected.iter().map(|(first, second, _)| {
                            first
                                .iter()
                                .map(|p| p.refs())
                                .chain(second.iter().map(|p| p.refs()))
                                .collect()
                        }),
                        |node| rebuilt.incident_to_node(node),
                        |edge| rebuilt.incident_to_edge(edge),
                    );
                };
                for (removal, positions) in removal_cases(count) {
                    let ids: Vec<_> = positions.into_iter().map(RelationId::from).collect();
                    let removed: Vec<_> = entries
                        .iter()
                        .enumerate()
                        .rev()
                        .filter(|(row, _)| ids.contains(&RelationId::from(*row)))
                        .map(|(row, entry)| (RelationId::from(row), entry.clone()))
                        .collect();
                    let mut compacted = original.clone();
                    let compaction = compacted.tracked_remove(&ids);
                    let rebuilt = VarVarBirelationSet::new(rebuild_rows(
                        compacted.clone().into_entries(),
                        &compaction,
                        removed.clone(),
                    ));
                    check(&rebuilt, &entries);
                    group.bench_function(
                        BenchmarkId::new("rows_rebuild", format!("{fixture}/{removal}")),
                        |b| {
                            b.iter_batched(
                                || (compacted.clone(), removed.clone()),
                                |(relations, removed)| {
                                    black_box(VarVarBirelationSet::new(rebuild_rows(
                                        black_box(relations).into_entries(),
                                        black_box(&compaction),
                                        black_box(removed),
                                    )))
                                },
                                BatchSize::LargeInput,
                            )
                        },
                    );
                }
                let reference_count = 2 * count * width;
                for (removal, positions) in removal_cases(reference_count) {
                    let removed: BTreeSet<_> = positions.iter().copied().collect();
                    let compaction = GraphCompaction::new(
                        Compaction::new(
                            reference_count,
                            positions.iter().copied().map(NodeId::from).collect(),
                        )
                        .unwrap(),
                        Compaction::new(
                            reference_count,
                            positions.iter().copied().map(EdgeId::from).collect(),
                        )
                        .unwrap(),
                    );
                    let expected: Vec<_> = entries
                        .iter()
                        .filter(|entry| {
                            entry
                                .0
                                .iter()
                                .map(|p| p.refs())
                                .chain(entry.1.iter().map(|p| p.refs()))
                                .all(|refs| {
                                    refs.node.is_none_or(|id| !removed.contains(&id.index()))
                                        && refs.edge.is_none_or(|id| !removed.contains(&id.index()))
                                })
                        })
                        .cloned()
                        .collect();
                    let compacted = original.compact(&compaction);
                    let rebuild = |relations: VarVarBirelationSet<NodeId, EdgeId, usize>| {
                        VarVarBirelationSet::new(
                            relations
                                .into_entries()
                                .into_iter()
                                .map(|(first, second, data)| {
                                    let first = first
                                        .into_iter()
                                        .map(|p| p.uncompact(&compaction))
                                        .collect();
                                    let second = second
                                        .into_iter()
                                        .map(|p| p.uncompact(&compaction))
                                        .collect();
                                    (first, second, data)
                                })
                                .collect(),
                        )
                    };
                    check(&rebuild(compacted.clone()), &expected);
                    group.bench_function(
                        BenchmarkId::new("participants_rebuild", format!("{fixture}/{removal}")),
                        |b| {
                            b.iter_batched(
                                || compacted.clone(),
                                |relations| black_box(rebuild(black_box(relations))),
                                BatchSize::LargeInput,
                            )
                        },
                    );
                }
            }
        }
    }
    group.finish();
}

criterion_group!(benches, graph, fixed, var, fixed_fixed, fixed_var, var_var);
criterion_main!(benches);
