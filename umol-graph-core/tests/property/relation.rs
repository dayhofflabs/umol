//! Construction, queries, participant mutation, exact-size iteration, and transport.
//!
//! The storage laws documented on [FixedRelationSet], [VarRelationSet],
//! [FixedFixedBirelationSet], [FixedVarBirelationSet], and [VarVarBirelationSet]
//! preserve complete input rows. Incidence follows [RelationParticipant::refs];
//! coincidence compares complete participant multisets independently in each factor.
//! A direct row scan supplies expected incidence, and occurrence counts supply multiset equality.
//! Generated participants may reference a node, an edge, both, or neither; labels distinguish
//! values with identical references. Small id/label domains exercise duplicates and shared refs.
//! Variable factors include empty sequences. Queries include independent inputs, reversed rows,
//! and changed labels with unchanged references to distinguish value equality from incidence.
//! Replacement laws from [FixedRelationSet::replace_participants] and
//! [FixedRelationSet::replace_participant] compare sequences of edits against independently edited
//! rows, scanning incidence after every step. Local/whole replacement equivalence supplements that
//! reference check; equality alone would not detect a stale index.
//! [VarRelationSet::replace_participants] and its single-participant methods use the same row
//! model with variable-length factors and positional payloads. Mixed edit sequences exercise
//! growing, shrinking, and empty rows; every step checks both reference spaces and all row
//! offsets through public readers. Insertion followed by removal also restores the prior value.
//! [FixedFixedBirelationSet::replace_participants] and its factor/local variants preserve
//! both distinguished factors against a row model. Generated references can overlap within
//! and across factors; direct union scans check that editing one factor retains references
//! in the other. Factor/local edits also agree with replacing both factors together.
//! [FixedVarBirelationSet::replace_participants] and its factor/local variants combine
//! shared-reference preservation with variable-buffer resizing. Independent rows check later-row
//! offsets after growing, shrinking, and emptying a factor. Insertion/removal roundtrips and
//! factor/local equivalence supplement the union scans and complete-value coincidence checks.
//! [VarVarBirelationSet::replace_participants] and its factor/local variants use the same
//! independent row model with two independently resized factors, including empty factors.
//! Mixed sequences check both sets of row boundaries through public readers after each edit;
//! each factor's insertion/removal roundtrip also checks restored union incidence.
//! Transport exercises the identity/composition laws of [FixedRelationSet::try_map] and its
//! peers, plus the positional preservation law of [FixedRelationSet::remap] and its peers.
//! It uses permutations of eight node and edge ids; unit cases cover partial mappings.

use std::iter;

use proptest::prelude::*;
use proptest::test_runner::{Config, TestCaseResult, TestRunner};
use rstest::rstest;
use umol_graph_core::{
    Correspondence, EdgeId, FixedFixedBirelationSet, FixedRelationSet, FixedVarBirelationSet,
    GraphCompaction, GraphCorrespondence, GraphRemapping, NodeId, ParticipantPosition,
    ParticipantRefs, RelationId, RelationParticipant, Remapping, VarRelationSet,
    VarVarBirelationSet,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct TestData(usize);

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct TestParticipant {
    node: Option<NodeId>,
    edge: Option<EdgeId>,
    label: u8,
}

impl RelationParticipant for TestParticipant {
    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        Some(Self {
            node: match self.node {
                Some(node) => Some(node.compact(compaction)?),
                None => None,
            },
            edge: match self.edge {
                Some(edge) => Some(edge.compact(compaction)?),
                None => None,
            },
            label: self.label,
        })
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
        Self {
            node: self.node.map(|node| node.uncompact(compaction)),
            edge: self.edge.map(|edge| edge.uncompact(compaction)),
            label: self.label,
        }
    }

    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self> {
        Some(Self {
            node: match self.node {
                Some(node) => Some(node.try_map(correspondence)?),
                None => None,
            },
            edge: match self.edge {
                Some(edge) => Some(edge.try_map(correspondence)?),
                None => None,
            },
            label: self.label,
        })
    }

    fn remap(self, remapping: &GraphRemapping) -> Self {
        Self {
            node: self.node.map(|node| node.remap(remapping)),
            edge: self.edge.map(|edge| edge.remap(remapping)),
            label: self.label,
        }
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: self.node,
            edge: self.edge,
        }
    }
}

fn participant_strategy() -> impl Strategy<Value = TestParticipant> {
    (prop::option::of(0u32..3), prop::option::of(0u32..3), 0u8..2).prop_map(
        |(node, edge, label)| TestParticipant {
            node: node.map(NodeId),
            edge: edge.map(EdgeId),
            label,
        },
    )
}

fn same_multiset(left: &[TestParticipant], right: &[TestParticipant]) -> bool {
    left.len() == right.len()
        && left.iter().all(|item| {
            left.iter().filter(|value| *value == item).count()
                == right.iter().filter(|value| *value == item).count()
        })
}

fn participant_mapping_strategy() -> impl Strategy<Value = (GraphCorrespondence, GraphRemapping)> {
    (
        Just((0..8).map(NodeId).collect::<Vec<_>>()).prop_shuffle(),
        Just((0..8).map(EdgeId).collect::<Vec<_>>()).prop_shuffle(),
    )
        .prop_map(|(nodes, edges)| {
            (
                GraphCorrespondence::new(
                    Correspondence::from_images(&nodes, 8),
                    Correspondence::from_images(&edges, 8),
                ),
                GraphRemapping::new(
                    Remapping::new(nodes).expect("permutation images"),
                    Remapping::new(edges).expect("permutation images"),
                ),
            )
        })
}

fn assert_relation_ids(
    mut iterator: impl ExactSizeIterator<Item = RelationId>,
    count: usize,
    prefix: usize,
) -> TestCaseResult {
    prop_assert_eq!(iterator.len(), count);
    prop_assert_eq!(iterator.size_hint(), (count, Some(count)));
    for index in 0..prefix {
        prop_assert_eq!(iterator.next(), Some(RelationId(index as u32)));
        let remaining = count - index - 1;
        prop_assert_eq!(iterator.len(), remaining);
        prop_assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    prop_assert_eq!(
        iterator.collect::<Vec<_>>(),
        (prefix..count)
            .map(|index| RelationId(index as u32))
            .collect::<Vec<_>>(),
    );
    Ok(())
}

fn assert_data_iter_mut<'a>(
    mut iterator: impl ExactSizeIterator<Item = &'a mut TestData>,
    count: usize,
    prefix: usize,
) -> TestCaseResult {
    prop_assert_eq!(iterator.len(), count);
    prop_assert_eq!(iterator.size_hint(), (count, Some(count)));
    for index in 0..prefix {
        let item = iterator.next().expect("expected generated data item");
        prop_assert_eq!(*item, TestData(index));
        item.0 += count;
        let remaining = count - index - 1;
        prop_assert_eq!(iterator.len(), remaining);
        prop_assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    for index in prefix..count {
        let item = iterator.next().expect("expected remaining data item");
        prop_assert_eq!(*item, TestData(index));
        item.0 += count;
        let remaining = count - index - 1;
        prop_assert_eq!(iterator.len(), remaining);
        prop_assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    prop_assert_eq!(iterator.next(), None);
    prop_assert_eq!(iterator.len(), 0);
    Ok(())
}

#[rstest]
fn test_fixed_relation_set_new_incidence() {
    let strategy = (
        prop::collection::vec(
            (prop::array::uniform3(participant_strategy()), any::<u8>()),
            0..8,
        ),
        prop::collection::vec(participant_strategy(), 0..5),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_fixed_relation_set_new_incidence"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(entries, query)| {
            let relations = FixedRelationSet::<TestParticipant, u8, 3>::new(entries.clone());
            prop_assert_eq!(relations.count(), entries.len());
            for (index, (participants, data)) in entries.iter().enumerate() {
                let id = RelationId(index as u32);
                prop_assert_eq!(relations.participants(id), participants.as_slice());
                prop_assert_eq!(relations.data(id), data);
            }
            for key in 0..4 {
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(_, (participants, _))| {
                        participants
                            .iter()
                            .any(|participant| participant.node == Some(NodeId(key)))
                    })
                    .map(|(index, _)| RelationId(index as u32))
                    .collect();
                prop_assert_eq!(relations.incident_to_node(NodeId(key)), &expected);
                prop_assert_eq!(
                    relations.has_incident_to_node(NodeId(key)),
                    !expected.is_empty()
                );
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(_, (participants, _))| {
                        participants
                            .iter()
                            .any(|participant| participant.edge == Some(EdgeId(key)))
                    })
                    .map(|(index, _)| RelationId(index as u32))
                    .collect();
                prop_assert_eq!(relations.incident_to_edge(EdgeId(key)), &expected);
                prop_assert_eq!(
                    relations.has_incident_to_edge(EdgeId(key)),
                    !expected.is_empty()
                );
            }
            let queries = iter::once(query).chain(entries.iter().flat_map(|(participants, _)| {
                let query: Vec<_> = participants.iter().rev().copied().collect();
                let mut changed = query.clone();
                if let Some(participant) = changed.first_mut() {
                    participant.label ^= 1;
                }
                [query, changed]
            }));
            for query in queries {
                for (index, (participants, _)) in entries.iter().enumerate() {
                    prop_assert_eq!(
                        relations.is_coincident(RelationId(index as u32), &query),
                        same_multiset(participants, &query)
                    );
                }
                for key in 0..4 {
                    let expected = entries
                        .iter()
                        .position(|(participants, _)| {
                            participants
                                .iter()
                                .any(|participant| participant.node == Some(NodeId(key)))
                                && same_multiset(participants, &query)
                        })
                        .map(|index| RelationId(index as u32));
                    prop_assert_eq!(relations.coincident_to_node(NodeId(key), &query), expected);
                    let expected = entries
                        .iter()
                        .position(|(participants, _)| {
                            participants
                                .iter()
                                .any(|participant| participant.edge == Some(EdgeId(key)))
                                && same_multiset(participants, &query)
                        })
                        .map(|index| RelationId(index as u32));
                    prop_assert_eq!(relations.coincident_to_edge(EdgeId(key), &query), expected);
                }
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

fn assert_fixed_relation_rows(
    relations: &FixedRelationSet<TestParticipant, u8, 3>,
    entries: &[([TestParticipant; 3], u8)],
    previous: &[TestParticipant],
    query: &[TestParticipant],
) -> TestCaseResult {
    prop_assert_eq!(relations.count(), entries.len());
    prop_assert_eq!(
        relations.ids().collect::<Vec<_>>(),
        (0..entries.len()).map(RelationId::from).collect::<Vec<_>>()
    );
    prop_assert_eq!(
        relations
            .iter()
            .map(|(id, row, data)| (id, *row, *data))
            .collect::<Vec<_>>(),
        entries
            .iter()
            .enumerate()
            .map(|(i, (row, data))| (RelationId::from(i), *row, *data))
            .collect::<Vec<_>>()
    );
    let mut queries = vec![previous.to_vec(), query.to_vec()];
    for (index, (row, data)) in entries.iter().enumerate() {
        prop_assert_eq!(relations.participants(RelationId::from(index)), row);
        prop_assert_eq!(relations.data(RelationId::from(index)), data);
        let reversed: Vec<_> = row.iter().rev().copied().collect();
        let mut changed = reversed.clone();
        changed[0].label ^= 1;
        queries.extend([reversed, changed]);
    }
    for query in &queries {
        for (index, (row, _)) in entries.iter().enumerate() {
            prop_assert_eq!(
                relations.is_coincident(RelationId::from(index), query),
                same_multiset(row, query)
            );
        }
    }
    for key in 0..4 {
        let nodes: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (row, _))| row.iter().any(|p| p.node == Some(NodeId(key))))
            .map(|(index, _)| RelationId::from(index))
            .collect();
        let edges: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (row, _))| row.iter().any(|p| p.edge == Some(EdgeId(key))))
            .map(|(index, _)| RelationId::from(index))
            .collect();
        prop_assert_eq!(relations.incident_to_node(NodeId(key)), &nodes);
        prop_assert_eq!(
            relations.has_incident_to_node(NodeId(key)),
            !nodes.is_empty()
        );
        prop_assert_eq!(relations.incident_to_edge(EdgeId(key)), &edges);
        prop_assert_eq!(
            relations.has_incident_to_edge(EdgeId(key)),
            !edges.is_empty()
        );
        for query in &queries {
            let node_hit = nodes
                .iter()
                .copied()
                .find(|id| same_multiset(&entries[id.index()].0, query));
            let edge_hit = edges
                .iter()
                .copied()
                .find(|id| same_multiset(&entries[id.index()].0, query));
            prop_assert_eq!(relations.coincident_to_node(NodeId(key), query), node_hit);
            prop_assert_eq!(relations.coincident_to_edge(EdgeId(key), query), edge_hit);
        }
    }
    Ok(())
}

#[rstest]
fn test_fixed_relation_set_replace_participants() {
    let strategy = (
        prop::collection::vec(
            (prop::array::uniform3(participant_strategy()), any::<u8>()),
            1..8,
        ),
        prop::collection::vec(
            (0usize..7, prop::array::uniform3(participant_strategy())),
            1..12,
        ),
        prop::collection::vec(participant_strategy(), 0..5),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_fixed_relation_set_replace_participants"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(mut entries, edits, query)| {
            let mut relations = FixedRelationSet::new(entries.clone());
            for (row, participants) in edits {
                let row = row % entries.len();
                let previous = entries[row].0;
                entries[row].0 = participants;
                relations.replace_participants(RelationId::from(row), participants);
                assert_fixed_relation_rows(&relations, &entries, &previous, &query)?;
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

#[rstest]
fn test_fixed_relation_set_replace_participant() {
    let strategy = (
        prop::collection::vec(
            (prop::array::uniform3(participant_strategy()), any::<u8>()),
            1..8,
        ),
        prop::collection::vec((0usize..7, 0usize..3, participant_strategy()), 1..12),
        prop::collection::vec(participant_strategy(), 0..5),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_fixed_relation_set_replace_participant"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(mut entries, edits, query)| {
            let mut relations = FixedRelationSet::new(entries.clone());
            let mut whole = relations.clone();
            for (row, position, participant) in edits {
                let row = row % entries.len();
                let previous = entries[row].0;
                entries[row].0[position] = participant;
                relations.replace_participant(
                    RelationId::from(row),
                    ParticipantPosition(position as u32),
                    participant,
                );
                whole.replace_participants(RelationId::from(row), entries[row].0);
                prop_assert_eq!(&relations, &whole);
                assert_fixed_relation_rows(&relations, &entries, &previous, &query)?;
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

#[rstest]
fn test_var_relation_set_new_incidence() {
    let strategy = (
        prop::collection::vec(
            (
                prop::collection::vec(participant_strategy(), 0..5),
                any::<u8>(),
            ),
            0..8,
        ),
        prop::collection::vec(participant_strategy(), 0..5),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_var_relation_set_new_incidence"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(entries, query)| {
            let relations = VarRelationSet::<TestParticipant, u8>::new(entries.clone());
            prop_assert_eq!(relations.count(), entries.len());
            for (index, (participants, data)) in entries.iter().enumerate() {
                let id = RelationId(index as u32);
                prop_assert_eq!(relations.participants(id), participants.as_slice());
                prop_assert_eq!(relations.data(id), data);
            }
            for key in 0..4 {
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(_, (participants, _))| {
                        participants
                            .iter()
                            .any(|participant| participant.node == Some(NodeId(key)))
                    })
                    .map(|(index, _)| RelationId(index as u32))
                    .collect();
                prop_assert_eq!(relations.incident_to_node(NodeId(key)), &expected);
                prop_assert_eq!(
                    relations.has_incident_to_node(NodeId(key)),
                    !expected.is_empty()
                );
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(_, (participants, _))| {
                        participants
                            .iter()
                            .any(|participant| participant.edge == Some(EdgeId(key)))
                    })
                    .map(|(index, _)| RelationId(index as u32))
                    .collect();
                prop_assert_eq!(relations.incident_to_edge(EdgeId(key)), &expected);
                prop_assert_eq!(
                    relations.has_incident_to_edge(EdgeId(key)),
                    !expected.is_empty()
                );
            }
            let queries = iter::once(query).chain(entries.iter().flat_map(|(participants, _)| {
                let query: Vec<_> = participants.iter().rev().copied().collect();
                let mut changed = query.clone();
                if let Some(participant) = changed.first_mut() {
                    participant.label ^= 1;
                }
                [query, changed]
            }));
            for query in queries {
                for (index, (participants, _)) in entries.iter().enumerate() {
                    prop_assert_eq!(
                        relations.is_coincident(RelationId(index as u32), &query),
                        same_multiset(participants, &query)
                    );
                }
                for key in 0..4 {
                    let expected = entries
                        .iter()
                        .position(|(participants, _)| {
                            participants
                                .iter()
                                .any(|participant| participant.node == Some(NodeId(key)))
                                && same_multiset(participants, &query)
                        })
                        .map(|index| RelationId(index as u32));
                    prop_assert_eq!(relations.coincident_to_node(NodeId(key), &query), expected);
                    let expected = entries
                        .iter()
                        .position(|(participants, _)| {
                            participants
                                .iter()
                                .any(|participant| participant.edge == Some(EdgeId(key)))
                                && same_multiset(participants, &query)
                        })
                        .map(|index| RelationId(index as u32));
                    prop_assert_eq!(relations.coincident_to_edge(EdgeId(key), &query), expected);
                }
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

fn assert_var_relation_rows(
    relations: &VarRelationSet<TestParticipant, Vec<u8>>,
    entries: &[(Vec<TestParticipant>, Vec<u8>)],
    previous: &[TestParticipant],
    query: &[TestParticipant],
) -> TestCaseResult {
    prop_assert_eq!(relations.count(), entries.len());
    prop_assert_eq!(
        relations.ids().collect::<Vec<_>>(),
        (0..entries.len()).map(RelationId::from).collect::<Vec<_>>()
    );
    prop_assert_eq!(
        relations
            .iter()
            .map(|(id, row, data)| (id, row.to_vec(), data.clone()))
            .collect::<Vec<_>>(),
        entries
            .iter()
            .enumerate()
            .map(|(i, (row, data))| (RelationId::from(i), row.clone(), data.clone()))
            .collect::<Vec<_>>()
    );
    let mut queries = vec![previous.to_vec(), query.to_vec()];
    for (index, (row, data)) in entries.iter().enumerate() {
        prop_assert_eq!(relations.participants(RelationId::from(index)), row);
        prop_assert_eq!(relations.data(RelationId::from(index)), data);
        let reversed: Vec<_> = row.iter().rev().copied().collect();
        let mut changed = reversed.clone();
        if let Some(participant) = changed.first_mut() {
            participant.label ^= 1;
        }
        queries.extend([reversed, changed]);
    }
    for query in &queries {
        for (index, (row, _)) in entries.iter().enumerate() {
            prop_assert_eq!(
                relations.is_coincident(RelationId::from(index), query),
                same_multiset(row, query)
            );
        }
    }
    for key in 0..4 {
        let nodes: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (row, _))| row.iter().any(|p| p.node == Some(NodeId(key))))
            .map(|(index, _)| RelationId::from(index))
            .collect();
        let edges: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (row, _))| row.iter().any(|p| p.edge == Some(EdgeId(key))))
            .map(|(index, _)| RelationId::from(index))
            .collect();
        prop_assert_eq!(relations.incident_to_node(NodeId(key)), &nodes);
        prop_assert_eq!(
            relations.has_incident_to_node(NodeId(key)),
            !nodes.is_empty()
        );
        prop_assert_eq!(relations.incident_to_edge(EdgeId(key)), &edges);
        prop_assert_eq!(
            relations.has_incident_to_edge(EdgeId(key)),
            !edges.is_empty()
        );
        for query in &queries {
            let node_hit = nodes
                .iter()
                .copied()
                .find(|id| same_multiset(&entries[id.index()].0, query));
            let edge_hit = edges
                .iter()
                .copied()
                .find(|id| same_multiset(&entries[id.index()].0, query));
            prop_assert_eq!(relations.coincident_to_node(NodeId(key), query), node_hit);
            prop_assert_eq!(relations.coincident_to_edge(EdgeId(key), query), edge_hit);
        }
    }
    Ok(())
}

#[rstest]
fn test_var_relation_set_replace_participants() {
    let strategy = (
        prop::collection::vec(
            (
                prop::collection::vec(participant_strategy(), 0..6),
                prop::collection::vec(any::<u8>(), 0..6),
            ),
            1..8,
        ),
        prop::collection::vec(
            (
                0usize..7,
                prop::collection::vec(participant_strategy(), 0..9),
            ),
            1..16,
        ),
        prop::collection::vec(participant_strategy(), 0..6),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_var_relation_set_replace_participants"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(mut entries, edits, query)| {
            let mut relations = VarRelationSet::new(entries.clone());
            for (row, participants) in edits {
                let row = row % entries.len();
                let previous = entries[row].0.clone();
                relations.replace_participants(RelationId::from(row), &participants);
                entries[row].0 = participants;
                assert_var_relation_rows(&relations, &entries, &previous, &query)?;
                let prior = relations.clone();
                relations.replace_participants(RelationId::from(row), &entries[row].0);
                prop_assert_eq!(&relations, &prior);
                assert_var_relation_rows(&relations, &entries, &previous, &query)?;
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

#[rstest]
fn test_var_relation_set_replace_participants_local_edits() {
    let strategy = (
        prop::collection::vec(
            (
                prop::collection::vec(participant_strategy(), 0..6),
                prop::collection::vec(any::<u8>(), 0..6),
            ),
            1..8,
        ),
        prop::collection::vec(
            (
                0u8..4,
                0usize..7,
                0usize..16,
                participant_strategy(),
                prop::collection::vec(participant_strategy(), 0..9),
            ),
            1..24,
        ),
        prop::collection::vec(participant_strategy(), 0..6),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_var_relation_set_replace_participants_local_edits"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(mut entries, edits, query)| {
            let mut relations = VarRelationSet::new(entries.clone());
            let mut whole = relations.clone();
            for (operation, row, position, participant, replacement) in edits {
                let row = row % entries.len();
                let id = RelationId::from(row);
                let previous = entries[row].0.clone();
                let len = previous.len();
                let operation = if len == 0 && operation != 0 {
                    2
                } else {
                    operation
                };
                match operation {
                    0 => {
                        entries[row].0 = replacement;
                        relations.replace_participants(id, &entries[row].0);
                    }
                    1 => {
                        let position = position % len;
                        entries[row].0[position] = participant;
                        relations.replace_participant(
                            id,
                            ParticipantPosition(position as u32),
                            participant,
                        );
                    }
                    2 => {
                        let position = position % (len + 1);
                        entries[row].0.insert(position, participant);
                        let prior = relations.clone();
                        relations.insert_participant(
                            id,
                            ParticipantPosition(position as u32),
                            participant,
                        );
                        let mut restored = relations.clone();
                        restored.remove_participant(id, ParticipantPosition(position as u32));
                        prop_assert_eq!(&restored, &prior);
                        let mut original_entries = entries.clone();
                        original_entries[row].0 = previous.clone();
                        assert_var_relation_rows(&restored, &original_entries, &previous, &query)?;
                    }
                    _ => {
                        let position = position % len;
                        entries[row].0.remove(position);
                        relations.remove_participant(id, ParticipantPosition(position as u32));
                    }
                }
                whole.replace_participants(id, &entries[row].0);
                prop_assert_eq!(&relations, &whole);
                assert_var_relation_rows(&relations, &entries, &previous, &query)?;
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

#[rstest]
fn test_fixed_fixed_birelation_set_new_incidence() {
    let strategy = (
        prop::collection::vec(
            (
                prop::array::uniform3(participant_strategy()),
                prop::array::uniform3(participant_strategy()),
                any::<u8>(),
            ),
            0..8,
        ),
        (
            prop::collection::vec(participant_strategy(), 0..5),
            prop::collection::vec(participant_strategy(), 0..5),
        ),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_fixed_fixed_birelation_set_new_incidence"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(entries, query)| {
            let relations =
                FixedFixedBirelationSet::<TestParticipant, 3, TestParticipant, 3, u8>::new(
                    entries.clone(),
                );
            prop_assert_eq!(relations.count(), entries.len());
            for (index, (first, second, data)) in entries.iter().enumerate() {
                let id = RelationId(index as u32);
                prop_assert_eq!(relations.participants_1(id), first.as_slice());
                prop_assert_eq!(relations.participants_2(id), second.as_slice());
                prop_assert_eq!(relations.data(id), data);
            }
            for key in 0..4 {
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(_, (first, second, _))| {
                        first
                            .iter()
                            .chain(second)
                            .any(|participant| participant.node == Some(NodeId(key)))
                    })
                    .map(|(index, _)| RelationId(index as u32))
                    .collect();
                prop_assert_eq!(relations.incident_to_node(NodeId(key)), &expected);
                prop_assert_eq!(
                    relations.has_incident_to_node(NodeId(key)),
                    !expected.is_empty()
                );
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(_, (first, second, _))| {
                        first
                            .iter()
                            .chain(second)
                            .any(|participant| participant.edge == Some(EdgeId(key)))
                    })
                    .map(|(index, _)| RelationId(index as u32))
                    .collect();
                prop_assert_eq!(relations.incident_to_edge(EdgeId(key)), &expected);
                prop_assert_eq!(
                    relations.has_incident_to_edge(EdgeId(key)),
                    !expected.is_empty()
                );
            }
            let queries = iter::once(query).chain(entries.iter().flat_map(|(first, second, _)| {
                let first: Vec<_> = first.iter().rev().copied().collect();
                let second: Vec<_> = second.iter().rev().copied().collect();
                let mut changed_first = first.clone();
                let mut changed_second = second.clone();
                if let Some(participant) = changed_first.first_mut() {
                    participant.label ^= 1;
                }
                if let Some(participant) = changed_second.first_mut() {
                    participant.label ^= 1;
                }
                [
                    (changed_first, second.clone()),
                    (first.clone(), changed_second),
                    (first, second),
                ]
            }));
            for (query_1, query_2) in queries {
                for (index, (first, second, _)) in entries.iter().enumerate() {
                    prop_assert_eq!(
                        relations.is_coincident(RelationId(index as u32), &query_1, &query_2),
                        same_multiset(first, &query_1) && same_multiset(second, &query_2)
                    );
                }
                for key in 0..4 {
                    let expected = entries
                        .iter()
                        .position(|(first, second, _)| {
                            first
                                .iter()
                                .chain(second)
                                .any(|participant| participant.node == Some(NodeId(key)))
                                && same_multiset(first, &query_1)
                                && same_multiset(second, &query_2)
                        })
                        .map(|index| RelationId(index as u32));
                    prop_assert_eq!(
                        relations.coincident_to_node(NodeId(key), &query_1, &query_2),
                        expected
                    );
                    let expected = entries
                        .iter()
                        .position(|(first, second, _)| {
                            first
                                .iter()
                                .chain(second)
                                .any(|participant| participant.edge == Some(EdgeId(key)))
                                && same_multiset(first, &query_1)
                                && same_multiset(second, &query_2)
                        })
                        .map(|index| RelationId(index as u32));
                    prop_assert_eq!(
                        relations.coincident_to_edge(EdgeId(key), &query_1, &query_2),
                        expected
                    );
                }
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

fn assert_fixed_fixed_birelation_rows(
    relations: &FixedFixedBirelationSet<TestParticipant, 2, TestParticipant, 3, Vec<u8>>,
    entries: &[([TestParticipant; 2], [TestParticipant; 3], Vec<u8>)],
    previous: &([TestParticipant; 2], [TestParticipant; 3]),
    query: &(Vec<TestParticipant>, Vec<TestParticipant>),
) -> TestCaseResult {
    prop_assert_eq!(relations.count(), entries.len());
    prop_assert_eq!(
        relations.ids().collect::<Vec<_>>(),
        (0..entries.len()).map(RelationId::from).collect::<Vec<_>>()
    );
    prop_assert_eq!(
        relations
            .iter()
            .map(|(id, a, b, data)| (id, *a, *b, data.clone()))
            .collect::<Vec<_>>(),
        entries
            .iter()
            .enumerate()
            .map(|(i, (a, b, data))| (RelationId::from(i), *a, *b, data.clone()))
            .collect::<Vec<_>>()
    );
    let mut queries = vec![(previous.0.to_vec(), previous.1.to_vec()), query.clone()];
    for (index, (a, b, data)) in entries.iter().enumerate() {
        let id = RelationId::from(index);
        prop_assert_eq!(relations.participants_1(id), a);
        prop_assert_eq!(relations.participants_2(id), b);
        prop_assert_eq!(relations.data(id), data);
        let reversed_1: Vec<_> = a.iter().rev().copied().collect();
        let reversed_2: Vec<_> = b.iter().rev().copied().collect();
        let mut changed_1 = reversed_1.clone();
        changed_1[0].label ^= 1;
        let mut changed_2 = reversed_2.clone();
        changed_2[0].label ^= 1;
        queries.extend([
            (changed_1, reversed_2.clone()),
            (reversed_1.clone(), changed_2),
            (reversed_1, reversed_2),
        ]);
    }
    for (query_1, query_2) in &queries {
        for (index, (a, b, _)) in entries.iter().enumerate() {
            prop_assert_eq!(
                relations.is_coincident(RelationId::from(index), query_1, query_2),
                same_multiset(a, query_1) && same_multiset(b, query_2)
            );
        }
    }
    for key in 0..4 {
        let nodes: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (a, b, _))| a.iter().chain(b).any(|p| p.node == Some(NodeId(key))))
            .map(|(index, _)| RelationId::from(index))
            .collect();
        let edges: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (a, b, _))| a.iter().chain(b).any(|p| p.edge == Some(EdgeId(key))))
            .map(|(index, _)| RelationId::from(index))
            .collect();
        prop_assert_eq!(relations.incident_to_node(NodeId(key)), &nodes);
        prop_assert_eq!(
            relations.has_incident_to_node(NodeId(key)),
            !nodes.is_empty()
        );
        prop_assert_eq!(relations.incident_to_edge(EdgeId(key)), &edges);
        prop_assert_eq!(
            relations.has_incident_to_edge(EdgeId(key)),
            !edges.is_empty()
        );
        for (query_1, query_2) in &queries {
            let node_hit = nodes.iter().copied().find(|id| {
                same_multiset(&entries[id.index()].0, query_1)
                    && same_multiset(&entries[id.index()].1, query_2)
            });
            let edge_hit = edges.iter().copied().find(|id| {
                same_multiset(&entries[id.index()].0, query_1)
                    && same_multiset(&entries[id.index()].1, query_2)
            });
            prop_assert_eq!(
                relations.coincident_to_node(NodeId(key), query_1, query_2),
                node_hit
            );
            prop_assert_eq!(
                relations.coincident_to_edge(EdgeId(key), query_1, query_2),
                edge_hit
            );
        }
    }
    Ok(())
}

#[rstest]
fn test_fixed_fixed_birelation_set_replace_participants() {
    let strategy = (
        prop::collection::vec(
            (
                prop::array::uniform2(participant_strategy()),
                prop::array::uniform3(participant_strategy()),
                prop::collection::vec(any::<u8>(), 0..6),
            ),
            1..8,
        ),
        prop::collection::vec(
            (
                0usize..7,
                prop::array::uniform2(participant_strategy()),
                prop::array::uniform3(participant_strategy()),
            ),
            1..16,
        ),
        (
            prop::collection::vec(participant_strategy(), 0..5),
            prop::collection::vec(participant_strategy(), 0..5),
        ),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_fixed_fixed_birelation_set_replace_participants"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(mut entries, edits, query)| {
            let mut relations = FixedFixedBirelationSet::new(entries.clone());
            for (row, first, second) in edits {
                let row = row % entries.len();
                let previous = (entries[row].0, entries[row].1);
                relations.replace_participants(RelationId::from(row), first, second);
                entries[row].0 = first;
                entries[row].1 = second;
                assert_fixed_fixed_birelation_rows(&relations, &entries, &previous, &query)?;
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

#[rstest]
fn test_fixed_fixed_birelation_set_replace_participants_local_edits() {
    let strategy = (
        prop::collection::vec(
            (
                prop::array::uniform2(participant_strategy()),
                prop::array::uniform3(participant_strategy()),
                prop::collection::vec(any::<u8>(), 0..6),
            ),
            1..8,
        ),
        prop::collection::vec(
            (
                0u8..5,
                0usize..7,
                0usize..6,
                prop::array::uniform2(participant_strategy()),
                prop::array::uniform3(participant_strategy()),
            ),
            1..24,
        ),
        (
            prop::collection::vec(participant_strategy(), 0..5),
            prop::collection::vec(participant_strategy(), 0..5),
        ),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_fixed_fixed_birelation_set_replace_participants_local_edits"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(mut entries, edits, query)| {
            let mut relations = FixedFixedBirelationSet::new(entries.clone());
            let mut whole = relations.clone();
            for (operation, row, position, first, second) in edits {
                let row = row % entries.len();
                let id = RelationId::from(row);
                let previous = (entries[row].0, entries[row].1);
                match operation {
                    0 => {
                        entries[row].0 = first;
                        entries[row].1 = second;
                        relations.replace_participants(id, first, second);
                    }
                    1 => {
                        entries[row].0 = first;
                        relations.replace_participants_1(id, first);
                    }
                    2 => {
                        entries[row].1 = second;
                        relations.replace_participants_2(id, second);
                    }
                    3 => {
                        let position = position % 2;
                        entries[row].0[position] = first[0];
                        relations.replace_participant_1(
                            id,
                            ParticipantPosition(position as u32),
                            first[0],
                        );
                    }
                    _ => {
                        let position = position % 3;
                        entries[row].1[position] = second[0];
                        relations.replace_participant_2(
                            id,
                            ParticipantPosition(position as u32),
                            second[0],
                        );
                    }
                }
                whole.replace_participants(id, entries[row].0, entries[row].1);
                prop_assert_eq!(&relations, &whole);
                assert_fixed_fixed_birelation_rows(&relations, &entries, &previous, &query)?;
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

#[rstest]
fn test_fixed_var_birelation_set_new_incidence() {
    let strategy = (
        prop::collection::vec(
            (
                prop::array::uniform3(participant_strategy()),
                prop::collection::vec(participant_strategy(), 0..5),
                any::<u8>(),
            ),
            0..8,
        ),
        (
            prop::collection::vec(participant_strategy(), 0..5),
            prop::collection::vec(participant_strategy(), 0..5),
        ),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_fixed_var_birelation_set_new_incidence"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(entries, query)| {
            let relations = FixedVarBirelationSet::<TestParticipant, 3, TestParticipant, u8>::new(
                entries.clone(),
            );
            prop_assert_eq!(relations.count(), entries.len());
            for (index, (first, second, data)) in entries.iter().enumerate() {
                let id = RelationId(index as u32);
                prop_assert_eq!(relations.participants_1(id), first.as_slice());
                prop_assert_eq!(relations.participants_2(id), second.as_slice());
                prop_assert_eq!(relations.data(id), data);
            }
            for key in 0..4 {
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(_, (first, second, _))| {
                        first
                            .iter()
                            .chain(second)
                            .any(|participant| participant.node == Some(NodeId(key)))
                    })
                    .map(|(index, _)| RelationId(index as u32))
                    .collect();
                prop_assert_eq!(relations.incident_to_node(NodeId(key)), &expected);
                prop_assert_eq!(
                    relations.has_incident_to_node(NodeId(key)),
                    !expected.is_empty()
                );
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(_, (first, second, _))| {
                        first
                            .iter()
                            .chain(second)
                            .any(|participant| participant.edge == Some(EdgeId(key)))
                    })
                    .map(|(index, _)| RelationId(index as u32))
                    .collect();
                prop_assert_eq!(relations.incident_to_edge(EdgeId(key)), &expected);
                prop_assert_eq!(
                    relations.has_incident_to_edge(EdgeId(key)),
                    !expected.is_empty()
                );
            }
            let queries = iter::once(query).chain(entries.iter().flat_map(|(first, second, _)| {
                let first: Vec<_> = first.iter().rev().copied().collect();
                let second: Vec<_> = second.iter().rev().copied().collect();
                let mut changed_first = first.clone();
                let mut changed_second = second.clone();
                if let Some(participant) = changed_first.first_mut() {
                    participant.label ^= 1;
                }
                if let Some(participant) = changed_second.first_mut() {
                    participant.label ^= 1;
                }
                [
                    (changed_first, second.clone()),
                    (first.clone(), changed_second),
                    (first, second),
                ]
            }));
            for (query_1, query_2) in queries {
                for (index, (first, second, _)) in entries.iter().enumerate() {
                    prop_assert_eq!(
                        relations.is_coincident(RelationId(index as u32), &query_1, &query_2),
                        same_multiset(first, &query_1) && same_multiset(second, &query_2)
                    );
                }
                for key in 0..4 {
                    let expected = entries
                        .iter()
                        .position(|(first, second, _)| {
                            first
                                .iter()
                                .chain(second)
                                .any(|participant| participant.node == Some(NodeId(key)))
                                && same_multiset(first, &query_1)
                                && same_multiset(second, &query_2)
                        })
                        .map(|index| RelationId(index as u32));
                    prop_assert_eq!(
                        relations.coincident_to_node(NodeId(key), &query_1, &query_2),
                        expected
                    );
                    let expected = entries
                        .iter()
                        .position(|(first, second, _)| {
                            first
                                .iter()
                                .chain(second)
                                .any(|participant| participant.edge == Some(EdgeId(key)))
                                && same_multiset(first, &query_1)
                                && same_multiset(second, &query_2)
                        })
                        .map(|index| RelationId(index as u32));
                    prop_assert_eq!(
                        relations.coincident_to_edge(EdgeId(key), &query_1, &query_2),
                        expected
                    );
                }
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

fn assert_fixed_var_birelation_rows(
    relations: &FixedVarBirelationSet<TestParticipant, 2, TestParticipant, Vec<u8>>,
    entries: &[([TestParticipant; 2], Vec<TestParticipant>, Vec<u8>)],
    previous: &([TestParticipant; 2], Vec<TestParticipant>),
    query: &(Vec<TestParticipant>, Vec<TestParticipant>),
) -> TestCaseResult {
    prop_assert_eq!(relations.count(), entries.len());
    prop_assert_eq!(
        relations.ids().collect::<Vec<_>>(),
        (0..entries.len()).map(RelationId::from).collect::<Vec<_>>()
    );
    prop_assert_eq!(
        relations
            .iter()
            .map(|(id, a, b, data)| (id, *a, b.to_vec(), data.clone()))
            .collect::<Vec<_>>(),
        entries
            .iter()
            .enumerate()
            .map(|(i, (a, b, data))| (RelationId::from(i), *a, b.clone(), data.clone()))
            .collect::<Vec<_>>()
    );
    let mut queries = vec![(previous.0.to_vec(), previous.1.to_vec()), query.clone()];
    for (index, (a, b, data)) in entries.iter().enumerate() {
        let id = RelationId::from(index);
        prop_assert_eq!(relations.participants_1(id), a);
        prop_assert_eq!(relations.participants_2(id), b);
        prop_assert_eq!(relations.data(id), data);
        let reversed_1: Vec<_> = a.iter().rev().copied().collect();
        let reversed_2: Vec<_> = b.iter().rev().copied().collect();
        let mut changed_1 = reversed_1.clone();
        changed_1[0].label ^= 1;
        let mut changed_2 = reversed_2.clone();
        if let Some(participant) = changed_2.first_mut() {
            participant.label ^= 1;
        }
        queries.extend([
            (changed_1, reversed_2.clone()),
            (reversed_1.clone(), changed_2),
            (reversed_1, reversed_2),
        ]);
    }
    for (query_1, query_2) in &queries {
        for (index, (a, b, _)) in entries.iter().enumerate() {
            prop_assert_eq!(
                relations.is_coincident(RelationId::from(index), query_1, query_2),
                same_multiset(a, query_1) && same_multiset(b, query_2)
            );
        }
    }
    for key in 0..4 {
        let nodes: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (a, b, _))| a.iter().chain(b).any(|p| p.node == Some(NodeId(key))))
            .map(|(index, _)| RelationId::from(index))
            .collect();
        let edges: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (a, b, _))| a.iter().chain(b).any(|p| p.edge == Some(EdgeId(key))))
            .map(|(index, _)| RelationId::from(index))
            .collect();
        prop_assert_eq!(relations.incident_to_node(NodeId(key)), &nodes);
        prop_assert_eq!(
            relations.has_incident_to_node(NodeId(key)),
            !nodes.is_empty()
        );
        prop_assert_eq!(relations.incident_to_edge(EdgeId(key)), &edges);
        prop_assert_eq!(
            relations.has_incident_to_edge(EdgeId(key)),
            !edges.is_empty()
        );
        for (query_1, query_2) in &queries {
            let node_hit = nodes.iter().copied().find(|id| {
                same_multiset(&entries[id.index()].0, query_1)
                    && same_multiset(&entries[id.index()].1, query_2)
            });
            let edge_hit = edges.iter().copied().find(|id| {
                same_multiset(&entries[id.index()].0, query_1)
                    && same_multiset(&entries[id.index()].1, query_2)
            });
            prop_assert_eq!(
                relations.coincident_to_node(NodeId(key), query_1, query_2),
                node_hit
            );
            prop_assert_eq!(
                relations.coincident_to_edge(EdgeId(key), query_1, query_2),
                edge_hit
            );
        }
    }
    Ok(())
}

#[rstest]
fn test_fixed_var_birelation_set_replace_participants() {
    let strategy = (
        prop::collection::vec(
            (
                prop::array::uniform2(participant_strategy()),
                prop::collection::vec(participant_strategy(), 0..7),
                prop::collection::vec(any::<u8>(), 0..6),
            ),
            1..8,
        ),
        prop::collection::vec(
            (
                0usize..7,
                prop::array::uniform2(participant_strategy()),
                prop::collection::vec(participant_strategy(), 0..7),
            ),
            1..16,
        ),
        (
            prop::collection::vec(participant_strategy(), 0..5),
            prop::collection::vec(participant_strategy(), 0..5),
        ),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_fixed_var_birelation_set_replace_participants"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(mut entries, edits, query)| {
            let mut relations = FixedVarBirelationSet::new(entries.clone());
            for (row, first, second) in edits {
                let row = row % entries.len();
                let previous = (entries[row].0, entries[row].1.clone());
                relations.replace_participants(RelationId::from(row), first, &second);
                entries[row].0 = first;
                entries[row].1 = second;
                assert_fixed_var_birelation_rows(&relations, &entries, &previous, &query)?;
                let before = relations.clone();
                relations.replace_participants(
                    RelationId::from(row),
                    entries[row].0,
                    &entries[row].1,
                );
                prop_assert_eq!(&relations, &before);
                assert_fixed_var_birelation_rows(&relations, &entries, &previous, &query)?;
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

#[rstest]
fn test_fixed_var_birelation_set_replace_participants_local_edits() {
    let strategy = (
        prop::collection::vec(
            (
                prop::array::uniform2(participant_strategy()),
                prop::collection::vec(participant_strategy(), 0..7),
                prop::collection::vec(any::<u8>(), 0..6),
            ),
            1..8,
        ),
        prop::collection::vec(
            (
                0u8..7,
                0usize..7,
                0usize..8,
                participant_strategy(),
                prop::array::uniform2(participant_strategy()),
                prop::collection::vec(participant_strategy(), 0..7),
            ),
            1..24,
        ),
        (
            prop::collection::vec(participant_strategy(), 0..5),
            prop::collection::vec(participant_strategy(), 0..5),
        ),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_fixed_var_birelation_set_replace_participants_local_edits"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(mut entries, edits, query)| {
            let mut relations = FixedVarBirelationSet::new(entries.clone());
            let mut whole = relations.clone();
            for (operation, row, position, participant, first, second) in edits {
                let row = row % entries.len();
                let id = RelationId::from(row);
                let previous = (entries[row].0, entries[row].1.clone());
                match operation {
                    0 => {
                        entries[row].0 = first;
                        relations.replace_participants(id, first, &second);
                        entries[row].1 = second;
                    }
                    1 => {
                        entries[row].0 = first;
                        relations.replace_participants_1(id, first);
                    }
                    2 => {
                        relations.replace_participants_2(id, &second);
                        entries[row].1 = second;
                    }
                    3 => {
                        let position = position % 2;
                        entries[row].0[position] = first[0];
                        relations.replace_participant_1(
                            id,
                            ParticipantPosition(position as u32),
                            first[0],
                        );
                    }
                    4 if !entries[row].1.is_empty() => {
                        let position = position % entries[row].1.len();
                        entries[row].1[position] = participant;
                        relations.replace_participant_2(
                            id,
                            ParticipantPosition(position as u32),
                            participant,
                        );
                    }
                    5 => {
                        let position = position % (entries[row].1.len() + 1);
                        let before = relations.clone();
                        entries[row].1.insert(position, participant);
                        relations.insert_participant_2(
                            id,
                            ParticipantPosition(position as u32),
                            participant,
                        );
                        let mut restored = relations.clone();
                        restored.remove_participant_2(id, ParticipantPosition(position as u32));
                        prop_assert_eq!(&restored, &before);
                        let mut before_entries = entries.clone();
                        before_entries[row].1.remove(position);
                        assert_fixed_var_birelation_rows(
                            &restored,
                            &before_entries,
                            &previous,
                            &query,
                        )?;
                    }
                    6 if !entries[row].1.is_empty() => {
                        let position = position % entries[row].1.len();
                        entries[row].1.remove(position);
                        relations.remove_participant_2(id, ParticipantPosition(position as u32));
                    }
                    _ => continue,
                }
                whole.replace_participants(id, entries[row].0, &entries[row].1);
                prop_assert_eq!(&relations, &whole);
                assert_fixed_var_birelation_rows(&relations, &entries, &previous, &query)?;
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

#[rstest]
fn test_var_var_birelation_set_new_incidence() {
    let strategy = (
        prop::collection::vec(
            (
                prop::collection::vec(participant_strategy(), 0..5),
                prop::collection::vec(participant_strategy(), 0..5),
                any::<u8>(),
            ),
            0..8,
        ),
        (
            prop::collection::vec(participant_strategy(), 0..5),
            prop::collection::vec(participant_strategy(), 0..5),
        ),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_var_var_birelation_set_new_incidence"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(entries, query)| {
            let relations =
                VarVarBirelationSet::<TestParticipant, TestParticipant, u8>::new(entries.clone());
            prop_assert_eq!(relations.count(), entries.len());
            for (index, (first, second, data)) in entries.iter().enumerate() {
                let id = RelationId(index as u32);
                prop_assert_eq!(relations.participants_1(id), first.as_slice());
                prop_assert_eq!(relations.participants_2(id), second.as_slice());
                prop_assert_eq!(relations.data(id), data);
            }
            for key in 0..4 {
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(_, (first, second, _))| {
                        first
                            .iter()
                            .chain(second)
                            .any(|participant| participant.node == Some(NodeId(key)))
                    })
                    .map(|(index, _)| RelationId(index as u32))
                    .collect();
                prop_assert_eq!(relations.incident_to_node(NodeId(key)), &expected);
                prop_assert_eq!(
                    relations.has_incident_to_node(NodeId(key)),
                    !expected.is_empty()
                );
                let expected: Vec<_> = entries
                    .iter()
                    .enumerate()
                    .filter(|(_, (first, second, _))| {
                        first
                            .iter()
                            .chain(second)
                            .any(|participant| participant.edge == Some(EdgeId(key)))
                    })
                    .map(|(index, _)| RelationId(index as u32))
                    .collect();
                prop_assert_eq!(relations.incident_to_edge(EdgeId(key)), &expected);
                prop_assert_eq!(
                    relations.has_incident_to_edge(EdgeId(key)),
                    !expected.is_empty()
                );
            }
            let queries = iter::once(query).chain(entries.iter().flat_map(|(first, second, _)| {
                let first: Vec<_> = first.iter().rev().copied().collect();
                let second: Vec<_> = second.iter().rev().copied().collect();
                let mut changed_first = first.clone();
                let mut changed_second = second.clone();
                if let Some(participant) = changed_first.first_mut() {
                    participant.label ^= 1;
                }
                if let Some(participant) = changed_second.first_mut() {
                    participant.label ^= 1;
                }
                [
                    (changed_first, second.clone()),
                    (first.clone(), changed_second),
                    (first, second),
                ]
            }));
            for (query_1, query_2) in queries {
                for (index, (first, second, _)) in entries.iter().enumerate() {
                    prop_assert_eq!(
                        relations.is_coincident(RelationId(index as u32), &query_1, &query_2),
                        same_multiset(first, &query_1) && same_multiset(second, &query_2)
                    );
                }
                for key in 0..4 {
                    let expected = entries
                        .iter()
                        .position(|(first, second, _)| {
                            first
                                .iter()
                                .chain(second)
                                .any(|participant| participant.node == Some(NodeId(key)))
                                && same_multiset(first, &query_1)
                                && same_multiset(second, &query_2)
                        })
                        .map(|index| RelationId(index as u32));
                    prop_assert_eq!(
                        relations.coincident_to_node(NodeId(key), &query_1, &query_2),
                        expected
                    );
                    let expected = entries
                        .iter()
                        .position(|(first, second, _)| {
                            first
                                .iter()
                                .chain(second)
                                .any(|participant| participant.edge == Some(EdgeId(key)))
                                && same_multiset(first, &query_1)
                                && same_multiset(second, &query_2)
                        })
                        .map(|index| RelationId(index as u32));
                    prop_assert_eq!(
                        relations.coincident_to_edge(EdgeId(key), &query_1, &query_2),
                        expected
                    );
                }
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

fn assert_var_var_birelation_rows(
    relations: &VarVarBirelationSet<TestParticipant, TestParticipant, Vec<u8>>,
    entries: &[(Vec<TestParticipant>, Vec<TestParticipant>, Vec<u8>)],
    previous: &(Vec<TestParticipant>, Vec<TestParticipant>),
    query: &(Vec<TestParticipant>, Vec<TestParticipant>),
) -> TestCaseResult {
    prop_assert_eq!(relations.count(), entries.len());
    prop_assert_eq!(
        relations.ids().collect::<Vec<_>>(),
        (0..entries.len()).map(RelationId::from).collect::<Vec<_>>()
    );
    prop_assert_eq!(
        relations
            .iter()
            .map(|(id, a, b, data)| (id, a.to_vec(), b.to_vec(), data.clone()))
            .collect::<Vec<_>>(),
        entries
            .iter()
            .enumerate()
            .map(|(i, (a, b, data))| (RelationId::from(i), a.clone(), b.clone(), data.clone()))
            .collect::<Vec<_>>()
    );
    let mut queries = vec![(previous.0.to_vec(), previous.1.to_vec()), query.clone()];
    for (index, (a, b, data)) in entries.iter().enumerate() {
        let id = RelationId::from(index);
        prop_assert_eq!(relations.participants_1(id), a);
        prop_assert_eq!(relations.participants_2(id), b);
        prop_assert_eq!(relations.data(id), data);
        let reversed_1: Vec<_> = a.iter().rev().copied().collect();
        let reversed_2: Vec<_> = b.iter().rev().copied().collect();
        let mut changed_1 = reversed_1.clone();
        if let Some(participant) = changed_1.first_mut() {
            participant.label ^= 1;
        }
        let mut changed_2 = reversed_2.clone();
        if let Some(participant) = changed_2.first_mut() {
            participant.label ^= 1;
        }
        queries.extend([
            (changed_1, reversed_2.clone()),
            (reversed_1.clone(), changed_2),
            (reversed_1, reversed_2),
        ]);
    }
    for (query_1, query_2) in &queries {
        for (index, (a, b, _)) in entries.iter().enumerate() {
            prop_assert_eq!(
                relations.is_coincident(RelationId::from(index), query_1, query_2),
                same_multiset(a, query_1) && same_multiset(b, query_2)
            );
        }
    }
    for key in 0..4 {
        let nodes: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (a, b, _))| a.iter().chain(b).any(|p| p.node == Some(NodeId(key))))
            .map(|(index, _)| RelationId::from(index))
            .collect();
        let edges: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (a, b, _))| a.iter().chain(b).any(|p| p.edge == Some(EdgeId(key))))
            .map(|(index, _)| RelationId::from(index))
            .collect();
        prop_assert_eq!(relations.incident_to_node(NodeId(key)), &nodes);
        prop_assert_eq!(
            relations.has_incident_to_node(NodeId(key)),
            !nodes.is_empty()
        );
        prop_assert_eq!(relations.incident_to_edge(EdgeId(key)), &edges);
        prop_assert_eq!(
            relations.has_incident_to_edge(EdgeId(key)),
            !edges.is_empty()
        );
        for (query_1, query_2) in &queries {
            let node_hit = nodes.iter().copied().find(|id| {
                same_multiset(&entries[id.index()].0, query_1)
                    && same_multiset(&entries[id.index()].1, query_2)
            });
            let edge_hit = edges.iter().copied().find(|id| {
                same_multiset(&entries[id.index()].0, query_1)
                    && same_multiset(&entries[id.index()].1, query_2)
            });
            prop_assert_eq!(
                relations.coincident_to_node(NodeId(key), query_1, query_2),
                node_hit
            );
            prop_assert_eq!(
                relations.coincident_to_edge(EdgeId(key), query_1, query_2),
                edge_hit
            );
        }
    }
    Ok(())
}

#[rstest]
fn test_var_var_birelation_set_replace_participants() {
    let strategy = (
        prop::collection::vec(
            (
                prop::collection::vec(participant_strategy(), 0..5),
                prop::collection::vec(participant_strategy(), 0..7),
                prop::collection::vec(any::<u8>(), 0..6),
            ),
            1..8,
        ),
        prop::collection::vec(
            (
                0usize..7,
                prop::collection::vec(participant_strategy(), 0..5),
                prop::collection::vec(participant_strategy(), 0..7),
            ),
            1..16,
        ),
        (
            prop::collection::vec(participant_strategy(), 0..5),
            prop::collection::vec(participant_strategy(), 0..5),
        ),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_var_var_birelation_set_replace_participants"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(mut entries, edits, query)| {
            let mut relations = VarVarBirelationSet::new(entries.clone());
            for (row, first, second) in edits {
                let row = row % entries.len();
                let previous = (entries[row].0.clone(), entries[row].1.clone());
                relations.replace_participants(RelationId::from(row), &first, &second);
                entries[row].0 = first;
                entries[row].1 = second;
                assert_var_var_birelation_rows(&relations, &entries, &previous, &query)?;
                let before = relations.clone();
                relations.replace_participants(
                    RelationId::from(row),
                    &entries[row].0,
                    &entries[row].1,
                );
                prop_assert_eq!(&relations, &before);
                assert_var_var_birelation_rows(&relations, &entries, &previous, &query)?;
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

#[rstest]
fn test_var_var_birelation_set_replace_participants_local_edits() {
    let strategy = (
        prop::collection::vec(
            (
                prop::collection::vec(participant_strategy(), 0..5),
                prop::collection::vec(participant_strategy(), 0..7),
                prop::collection::vec(any::<u8>(), 0..6),
            ),
            1..8,
        ),
        prop::collection::vec(
            (
                0u8..9,
                0usize..7,
                0usize..8,
                participant_strategy(),
                prop::collection::vec(participant_strategy(), 0..5),
                prop::collection::vec(participant_strategy(), 0..7),
            ),
            1..24,
        ),
        (
            prop::collection::vec(participant_strategy(), 0..5),
            prop::collection::vec(participant_strategy(), 0..5),
        ),
    );
    let config = Config {
        source_file: Some(file!()),
        test_name: Some(concat!(
            module_path!(),
            "::test_var_var_birelation_set_replace_participants_local_edits"
        )),
        ..Config::default()
    };
    TestRunner::new(config)
        .run(&strategy, |(mut entries, edits, query)| {
            let mut relations = VarVarBirelationSet::new(entries.clone());
            let mut whole = relations.clone();
            for (operation, row, position, participant, first, second) in edits {
                let row = row % entries.len();
                let id = RelationId::from(row);
                let previous = (entries[row].0.clone(), entries[row].1.clone());
                match operation {
                    0 => {
                        relations.replace_participants(id, &first, &second);
                        entries[row].0 = first;
                        entries[row].1 = second;
                    }
                    1 => {
                        relations.replace_participants_1(id, &first);
                        entries[row].0 = first;
                    }
                    2 => {
                        relations.replace_participants_2(id, &second);
                        entries[row].1 = second;
                    }
                    3 if !entries[row].0.is_empty() => {
                        let position = position % entries[row].0.len();
                        entries[row].0[position] = participant;
                        relations.replace_participant_1(
                            id,
                            ParticipantPosition(position as u32),
                            participant,
                        );
                    }
                    4 if !entries[row].1.is_empty() => {
                        let position = position % entries[row].1.len();
                        entries[row].1[position] = participant;
                        relations.replace_participant_2(
                            id,
                            ParticipantPosition(position as u32),
                            participant,
                        );
                    }
                    5 => {
                        let position = position % (entries[row].1.len() + 1);
                        let before = relations.clone();
                        entries[row].1.insert(position, participant);
                        relations.insert_participant_2(
                            id,
                            ParticipantPosition(position as u32),
                            participant,
                        );
                        let mut restored = relations.clone();
                        restored.remove_participant_2(id, ParticipantPosition(position as u32));
                        prop_assert_eq!(&restored, &before);
                        let mut before_entries = entries.clone();
                        before_entries[row].1.remove(position);
                        assert_var_var_birelation_rows(
                            &restored,
                            &before_entries,
                            &previous,
                            &query,
                        )?;
                    }
                    6 if !entries[row].1.is_empty() => {
                        let position = position % entries[row].1.len();
                        entries[row].1.remove(position);
                        relations.remove_participant_2(id, ParticipantPosition(position as u32));
                    }
                    7 => {
                        let position = position % (entries[row].0.len() + 1);
                        let before = relations.clone();
                        entries[row].0.insert(position, participant);
                        relations.insert_participant_1(
                            id,
                            ParticipantPosition(position as u32),
                            participant,
                        );
                        let mut restored = relations.clone();
                        restored.remove_participant_1(id, ParticipantPosition(position as u32));
                        prop_assert_eq!(&restored, &before);
                        let mut before_entries = entries.clone();
                        before_entries[row].0.remove(position);
                        assert_var_var_birelation_rows(
                            &restored,
                            &before_entries,
                            &previous,
                            &query,
                        )?;
                    }
                    8 if !entries[row].0.is_empty() => {
                        let position = position % entries[row].0.len();
                        entries[row].0.remove(position);
                        relations.remove_participant_1(id, ParticipantPosition(position as u32));
                    }
                    _ => continue,
                }
                whole.replace_participants(id, &entries[row].0, &entries[row].1);
                prop_assert_eq!(&relations, &whole);
                assert_var_var_birelation_rows(&relations, &entries, &previous, &query)?;
            }
            prop_assert_eq!(relations.into_entries(), entries);
            Ok(())
        })
        .unwrap();
}

proptest! {

    #[test]
    fn test_fixed_relation_set_map_composition(
        rows in prop::collection::vec((prop::array::uniform2(0u32..8), prop::array::uniform2(0u32..8), any::<[u32; 2]>()), 0..8),
        (first, remapping) in participant_mapping_strategy(),
        (second, _) in participant_mapping_strategy(),
    ) {
        let input: FixedRelationSet<NodeId, [u32; 2], 2> = FixedRelationSet::new(rows.into_iter().map(|(_edges, nodes, data)| (nodes.map(NodeId), data)).collect());
        let mapped = input.map(&first);
        prop_assert_eq!(&mapped, &input.remap(&remapping));
        let composed = first.compose(&second).unwrap();
        prop_assert_eq!(mapped.map(&second), input.map(&composed));
        let inverse = GraphCorrespondence::new(first.nodes().reverse(), first.edges().reverse());
        prop_assert_eq!(mapped.map(&inverse), input);
    }

    #[test]
    fn test_var_relation_set_map_composition(
        rows in prop::collection::vec((prop::array::uniform2(0u32..8), prop::array::uniform2(0u32..8), any::<[u32; 2]>()), 0..8),
        (first, remapping) in participant_mapping_strategy(),
        (second, _) in participant_mapping_strategy(),
    ) {
        let input: VarRelationSet<NodeId, [u32; 2]> = VarRelationSet::new(rows.into_iter().map(|(_edges, nodes, data)| (nodes.map(NodeId).to_vec(), data)).collect());
        let mapped = input.map(&first);
        prop_assert_eq!(&mapped, &input.remap(&remapping));
        let composed = first.compose(&second).unwrap();
        prop_assert_eq!(mapped.map(&second), input.map(&composed));
        let inverse = GraphCorrespondence::new(first.nodes().reverse(), first.edges().reverse());
        prop_assert_eq!(mapped.map(&inverse), input);
    }

    #[test]
    fn test_fixed_fixed_birelation_set_map_composition(
        rows in prop::collection::vec((prop::array::uniform2(0u32..8), prop::array::uniform2(0u32..8), any::<[u32; 2]>()), 0..8),
        (first, remapping) in participant_mapping_strategy(),
        (second, _) in participant_mapping_strategy(),
    ) {
        let input: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, [u32; 2]> = FixedFixedBirelationSet::new(rows.into_iter().map(|(edges, nodes, data)| (edges.map(EdgeId), nodes.map(NodeId), data)).collect());
        let mapped = input.map(&first);
        prop_assert_eq!(&mapped, &input.remap(&remapping));
        let composed = first.compose(&second).unwrap();
        prop_assert_eq!(mapped.map(&second), input.map(&composed));
        let inverse = GraphCorrespondence::new(first.nodes().reverse(), first.edges().reverse());
        prop_assert_eq!(mapped.map(&inverse), input);
    }

    #[test]
    fn test_fixed_var_birelation_set_map_composition(
        rows in prop::collection::vec((prop::array::uniform2(0u32..8), prop::array::uniform2(0u32..8), any::<[u32; 2]>()), 0..8),
        (first, remapping) in participant_mapping_strategy(),
        (second, _) in participant_mapping_strategy(),
    ) {
        let input: FixedVarBirelationSet<EdgeId, 2, NodeId, [u32; 2]> = FixedVarBirelationSet::new(rows.into_iter().map(|(edges, nodes, data)| (edges.map(EdgeId), nodes.map(NodeId).to_vec(), data)).collect());
        let mapped = input.map(&first);
        prop_assert_eq!(&mapped, &input.remap(&remapping));
        let composed = first.compose(&second).unwrap();
        prop_assert_eq!(mapped.map(&second), input.map(&composed));
        let inverse = GraphCorrespondence::new(first.nodes().reverse(), first.edges().reverse());
        prop_assert_eq!(mapped.map(&inverse), input);
    }

    #[test]
    fn test_var_var_birelation_set_map_composition(
        rows in prop::collection::vec((prop::array::uniform2(0u32..8), prop::array::uniform2(0u32..8), any::<[u32; 2]>()), 0..8),
        (first, remapping) in participant_mapping_strategy(),
        (second, _) in participant_mapping_strategy(),
    ) {
        let input: VarVarBirelationSet<EdgeId, NodeId, [u32; 2]> = VarVarBirelationSet::new(rows.into_iter().map(|(edges, nodes, data)| (edges.map(EdgeId).to_vec(), nodes.map(NodeId).to_vec(), data)).collect());
        let mapped = input.map(&first);
        prop_assert_eq!(&mapped, &input.remap(&remapping));
        let composed = first.compose(&second).unwrap();
        prop_assert_eq!(mapped.map(&second), input.map(&composed));
        let inverse = GraphCorrespondence::new(first.nodes().reverse(), first.edges().reverse());
        prop_assert_eq!(mapped.map(&inverse), input);
    }
    #[test]
    fn test_fixed_relation_set_iterators_exact_size(
        count in 0usize..16,
        prefix in any::<usize>(),
    ) {
        let entries = (0..count)
            .map(|index| {
                (
                    [NodeId((2 * index) as u32), NodeId((2 * index + 1) as u32)],
                    TestData(index),
                )
            })
            .collect::<Vec<_>>();
        let mut relations = FixedRelationSet::<NodeId, TestData, 2>::new(entries);
        let prefix = prefix.min(count);

        assert_relation_ids(relations.ids(), count, prefix)?;
        assert_data_iter_mut(
            relations.iter_mut().map(|(_, _, data)| data),
            count,
            prefix,
        )?;
        for index in 0..count {
            prop_assert_eq!(
                relations.data(RelationId(index as u32)),
                &TestData(index + count),
            );
        }
    }

    #[test]
    fn test_var_relation_set_iterators_exact_size(
        count in 0usize..16,
        prefix in any::<usize>(),
    ) {
        let entries = (0..count)
            .map(|index| {
                (
                    vec![NodeId((2 * index) as u32), NodeId((2 * index + 1) as u32)],
                    TestData(index),
                )
            })
            .collect::<Vec<_>>();
        let mut relations = VarRelationSet::<NodeId, TestData>::new(entries);
        let prefix = prefix.min(count);

        assert_relation_ids(relations.ids(), count, prefix)?;
        assert_data_iter_mut(
            relations.iter_mut().map(|(_, _, data)| data),
            count,
            prefix,
        )?;
        for index in 0..count {
            prop_assert_eq!(
                relations.data(RelationId(index as u32)),
                &TestData(index + count),
            );
        }
    }

    #[test]
    fn test_fixed_fixed_birelation_set_iterators_exact_size(
        count in 0usize..16,
        prefix in any::<usize>(),
    ) {
        let entries = (0..count)
            .map(|index| {
                (
                    [NodeId(index as u32)],
                    [EdgeId(index as u32)],
                    TestData(index),
                )
            })
            .collect::<Vec<_>>();
        let mut relations = FixedFixedBirelationSet::<NodeId, 1, EdgeId, 1, TestData, >::new(entries);
        let prefix = prefix.min(count);

        assert_relation_ids(relations.ids(), count, prefix)?;
        assert_data_iter_mut(
            relations.iter_mut().map(|(_, _, _, data)| data),
            count,
            prefix,
        )?;
        for index in 0..count {
            prop_assert_eq!(
                relations.data(RelationId(index as u32)),
                &TestData(index + count),
            );
        }
    }

    #[test]
    fn test_fixed_var_birelation_set_iterators_exact_size(
        count in 0usize..16,
        prefix in any::<usize>(),
    ) {
        let entries = (0..count)
            .map(|index| {
                (
                    [NodeId(index as u32)],
                    vec![EdgeId(index as u32)],
                    TestData(index),
                )
            })
            .collect::<Vec<_>>();
        let mut relations = FixedVarBirelationSet::<NodeId, 1, EdgeId, TestData, >::new(entries);
        let prefix = prefix.min(count);

        assert_relation_ids(relations.ids(), count, prefix)?;
        assert_data_iter_mut(
            relations.iter_mut().map(|(_, _, _, data)| data),
            count,
            prefix,
        )?;
        for index in 0..count {
            prop_assert_eq!(
                relations.data(RelationId(index as u32)),
                &TestData(index + count),
            );
        }
    }

    #[test]
    fn test_var_var_birelation_set_iterators_exact_size(
        count in 0usize..16,
        prefix in any::<usize>(),
    ) {
        let entries = (0..count)
            .map(|index| {
                (
                    vec![NodeId(index as u32)],
                    vec![EdgeId(index as u32)],
                    TestData(index),
                )
            })
            .collect::<Vec<_>>();
        let mut relations = VarVarBirelationSet::<NodeId, EdgeId, TestData, >::new(entries);
        let prefix = prefix.min(count);

        assert_relation_ids(relations.ids(), count, prefix)?;
        assert_data_iter_mut(
            relations.iter_mut().map(|(_, _, _, data)| data),
            count,
            prefix,
        )?;
        for index in 0..count {
            prop_assert_eq!(
                relations.data(RelationId(index as u32)),
                &TestData(index + count),
            );
        }
    }
}
