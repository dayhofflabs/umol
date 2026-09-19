use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::panic::{catch_unwind, AssertUnwindSafe};

use pretty_assertions::assert_eq;
use rstest::{fixture, rstest};

use crate::{
    Compaction, Correspondence, EdgeId, GraphCompaction, GraphCorrespondence, GraphRemapping,
    NodeId, ParticipantPosition, ParticipantRefs, RelationId, RelationParticipant,
    RelationPullbackCorrespondence, RelationPushoutCorrespondence, Remapping, VarRelationSet,
};

#[fixture]
fn participant_correspondence() -> GraphCorrespondence {
    GraphCorrespondence::new(
        Correspondence::new(vec![(NodeId(0), NodeId(5)), (NodeId(2), NodeId(1))], 4, 6).unwrap(),
        Correspondence::new(vec![(EdgeId(0), EdgeId(6)), (EdgeId(2), EdgeId(3))], 4, 7).unwrap(),
    )
}

fn hash<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct PositionLabels(Vec<u32>);

#[derive(Debug, PartialEq, Eq)]
struct NonCloneData(Vec<u32>);

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct References {
    node: Option<NodeId>,
    edge: Option<EdgeId>,
    label: u8,
}

impl RelationParticipant for References {
    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self> {
        Some(Self {
            node: match self.node {
                Some(id) => Some(id.try_map(correspondence)?),
                None => None,
            },
            edge: match self.edge {
                Some(id) => Some(id.try_map(correspondence)?),
                None => None,
            },
            ..self
        })
    }

    fn remap(self, remapping: &GraphRemapping) -> Self {
        Self {
            node: self.node.map(|id| id.remap(remapping)),
            edge: self.edge.map(|id| id.remap(remapping)),
            ..self
        }
    }

    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        Some(Self {
            node: match self.node {
                Some(id) => Some(id.compact(compaction)?),
                None => None,
            },
            edge: match self.edge {
                Some(id) => Some(id.compact(compaction)?),
                None => None,
            },
            ..self
        })
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
        Self {
            node: self.node.map(|id| id.uncompact(compaction)),
            edge: self.edge.map(|id| id.uncompact(compaction)),
            ..self
        }
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: self.node,
            edge: self.edge,
        }
    }
}

fn assert_exact_size<T>(mut iterator: impl ExactSizeIterator<Item = T>, expected: Vec<T>)
where
    T: Debug + PartialEq,
{
    assert_eq!(iterator.len(), expected.len());
    assert_eq!(iterator.size_hint(), (expected.len(), Some(expected.len())));
    while let Some(expected_item) = expected.get(expected.len() - iterator.len()) {
        let previous = iterator.len();
        assert_eq!(iterator.next().as_ref(), Some(expected_item));
        let remaining = iterator.len();
        assert_eq!(remaining, previous - 1);
        assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    assert_eq!(iterator.next(), None);
    assert_eq!(iterator.len(), 0);
    assert_eq!(iterator.size_hint(), (0, Some(0)));
}

#[fixture]
fn var_relation_set_mutation_entries() -> Vec<(Vec<NodeId>, PositionLabels)> {
    vec![
        (
            vec![NodeId(0), NodeId(1), NodeId(0)],
            PositionLabels(vec![7, 11, 13]),
        ),
        (
            vec![NodeId(1), NodeId(1), NodeId(2)],
            PositionLabels(vec![17, 19, 23]),
        ),
        (vec![], PositionLabels(vec![29])),
        (vec![NodeId(3), NodeId(4)], PositionLabels(vec![31, 37])),
    ]
}

fn assert_var_relation_rows(
    relations: &VarRelationSet<NodeId, PositionLabels>,
    entries: &[(Vec<NodeId>, PositionLabels)],
) {
    assert_eq!(relations.count(), entries.len());
    assert_eq!(
        relations.ids().collect::<Vec<_>>(),
        (0..entries.len()).map(RelationId::from).collect::<Vec<_>>()
    );
    for (index, (participants, data)) in entries.iter().enumerate() {
        assert_eq!(
            relations.participants(RelationId::from(index)),
            participants
        );
        assert_eq!(relations.data(RelationId::from(index)), data);
        let query: Vec<_> = participants.iter().rev().copied().collect();
        assert!(relations.is_coincident(RelationId::from(index), &query));
        for node in [0, 1, 2, 3, 4, 5, u32::MAX].map(NodeId) {
            let expected = entries
                .iter()
                .position(|(row, _)| {
                    row.contains(&node)
                        && row.len() == query.len()
                        && row.iter().all(|p| {
                            row.iter().filter(|q| *q == p).count()
                                == query.iter().filter(|q| *q == p).count()
                        })
                })
                .map(RelationId::from);
            assert_eq!(relations.coincident_to_node(node, &query), expected);
        }
    }
    for node in [0, 1, 2, 3, 4, 5, u32::MAX].map(NodeId) {
        let expected: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (row, _))| row.contains(&node))
            .map(|(i, _)| RelationId::from(i))
            .collect();
        assert_eq!(relations.incident_to_node(node), expected);
        assert_eq!(relations.has_incident_to_node(node), !expected.is_empty());
        assert_eq!(relations.incident_to_edge(EdgeId(node.0)), &[]);
        assert!(!relations.has_incident_to_edge(EdgeId(node.0)));
    }
    assert_eq!(relations.clone().into_entries(), entries);
}

#[fixture]
fn var_relation_set_compaction_input() -> VarRelationSet<NodeId, &'static str> {
    VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(2), NodeId(4)], "keep"),
        (vec![NodeId(1), NodeId(3)], "drop"),
    ])
}

#[rstest]
fn test_var_relation_set_new() {
    let rs: VarRelationSet<NodeId, &str> = VarRelationSet::new(vec![(
        vec![
            NodeId(0),
            NodeId(1),
            NodeId(2),
            NodeId(3),
            NodeId(4),
            NodeId(5),
        ],
        "benzene",
    )]);
    assert_eq!(rs.count(), 1);
    assert_eq!(rs.data(RelationId(0)), &"benzene");
    assert_eq!(
        rs.participants(RelationId(0)),
        &[
            NodeId(0),
            NodeId(1),
            NodeId(2),
            NodeId(3),
            NodeId(4),
            NodeId(5)
        ]
    );
}

#[rstest]
#[case::empty(vec![], vec![], vec![])]
#[case::repeated(
    vec![
        (vec![NodeId(2), NodeId(0), NodeId(2)], "first"),
        (vec![NodeId(0), NodeId(2), NodeId(2)], "duplicate"),
        (vec![NodeId(3), NodeId(3), NodeId(3)], "other"),
    ],
    vec![RelationId(0), RelationId(1)],
    vec![RelationId(2)],
)]
#[case::empty_factors(vec![(vec![], "empty")], vec![], vec![])]
fn test_var_relation_set_new_incidence(
    #[case] entries: Vec<(Vec<NodeId>, &str)>,
    #[case] at_two: Vec<RelationId>,
    #[case] at_three: Vec<RelationId>,
) {
    let relations = VarRelationSet::<NodeId, &str>::new(entries.clone());
    assert_eq!(relations.count(), entries.len());
    for (index, (participants, data)) in entries.iter().enumerate() {
        let id = RelationId(index as u32);
        assert_eq!(relations.participants(id), participants);
        assert_eq!(relations.data(id), data);
    }
    assert_eq!(relations.incident_to_node(NodeId(2)), at_two);
    assert_eq!(relations.incident_to_node(NodeId(3)), at_three);
    assert_eq!(relations.incident_to_node(NodeId(4)), &[]);
    assert_eq!(relations.incident_to_edge(EdgeId(2)), &[]);
    assert_eq!(
        relations.has_incident_to_node(NodeId(2)),
        !at_two.is_empty()
    );
    assert!(!relations.has_incident_to_edge(EdgeId(2)));
    assert_eq!(relations.into_entries(), entries);
}

#[rstest]
#[case::roundtrip(
    vec![
        (vec![NodeId(2), NodeId(0)], "first"),
        (vec![NodeId(4), NodeId(3), NodeId(1)], "second"),
    ],
)]
fn test_var_relation_set_into_entries(#[case] entries: Vec<(Vec<NodeId>, &str)>) {
    let rs: VarRelationSet<NodeId, &str> = VarRelationSet::new(entries.clone());
    assert_eq!(rs.into_entries(), entries);
}

#[rstest]
#[case::first(RelationId(0), true)]
#[case::out_of_range(RelationId(1), false)]
fn test_var_relation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
    let rs: VarRelationSet<NodeId, ()> =
        VarRelationSet::new(vec![(vec![NodeId(0), NodeId(1)], ())]);
    assert_eq!(rs.contains(id), expected);
}

#[rstest]
fn test_var_relation_set_relation_ids() {
    assert_exact_size(VarRelationSet::<NodeId, ()>::default().ids(), vec![]);
    let rs: VarRelationSet<NodeId, ()> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], ()),
        (vec![NodeId(1), NodeId(2)], ()),
    ]);
    assert_exact_size(rs.ids(), vec![RelationId(0), RelationId(1)]);
}

#[rstest]
fn test_var_relation_set_iter() {
    let empty = VarRelationSet::<NodeId, i32>::default();
    assert_eq!(empty.iter().collect::<Vec<_>>(), vec![]);

    let rs: VarRelationSet<NodeId, i32> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], 1),
        (vec![NodeId(2), NodeId(3), NodeId(4)], 2),
        (vec![NodeId(5)], 3),
    ]);
    assert_eq!(rs.iter().len(), 3);
    assert_eq!(
        rs.iter().collect::<Vec<_>>(),
        vec![
            (RelationId(0), [NodeId(0), NodeId(1)].as_slice(), &1),
            (
                RelationId(1),
                [NodeId(2), NodeId(3), NodeId(4)].as_slice(),
                &2
            ),
            (RelationId(2), [NodeId(5)].as_slice(), &3),
        ],
    );
}

#[rstest]
fn test_var_relation_set_iter_mut() {
    let mut empty = VarRelationSet::<NodeId, i32>::default();
    assert_eq!(empty.iter_mut().len(), 0);

    let mut rs: VarRelationSet<NodeId, i32> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], 1),
        (vec![NodeId(2), NodeId(3), NodeId(4)], 2),
        (vec![NodeId(5)], 3),
    ]);
    let arities: Vec<usize> = rs
        .iter_mut()
        .map(|(_, participants, data)| {
            *data *= 10;
            participants.len()
        })
        .collect();
    assert_eq!(arities, vec![2, 3, 1]);
    assert_eq!(rs.data(RelationId(0)), &10);
    assert_eq!(rs.data(RelationId(1)), &20);
    assert_eq!(rs.data(RelationId(2)), &30);
}

#[rstest]
fn test_var_relation_set_participants_ordered() {
    let rs: VarRelationSet<NodeId, ()> = VarRelationSet::new(vec![
        (vec![NodeId(5), NodeId(2), NodeId(0), NodeId(3)], ()),
        (vec![NodeId(4), NodeId(1)], ()),
    ]);
    assert_eq!(
        rs.participants(RelationId(0)),
        &[NodeId(5), NodeId(2), NodeId(0), NodeId(3)]
    );
    assert_eq!(rs.participants(RelationId(1)), &[NodeId(4), NodeId(1)]);
    assert_eq!(rs.incident_to_node(NodeId(0)), &[RelationId(0)]);
    assert_eq!(rs.incident_to_node(NodeId(4)), &[RelationId(1)]);
}

#[rstest]
fn test_var_relation_set_variable_arity() {
    let rs: VarRelationSet<NodeId, &str> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], "pair"),
        (vec![NodeId(2), NodeId(3), NodeId(4), NodeId(5)], "quad"),
    ]);
    assert_eq!(rs.participants(RelationId(0)), &[NodeId(0), NodeId(1)]);
    assert_eq!(
        rs.participants(RelationId(1)),
        &[NodeId(2), NodeId(3), NodeId(4), NodeId(5)]
    );
}

#[rstest]
fn test_var_relation_set_data_mut() {
    let mut rs: VarRelationSet<NodeId, i32> =
        VarRelationSet::new(vec![(vec![NodeId(0), NodeId(1), NodeId(2)], 1)]);
    *rs.data_mut(RelationId(0)) = 99;
    assert_eq!(rs.data(RelationId(0)), &99);
}

#[rstest]
fn test_var_relation_set_incidence() {
    let rs: VarRelationSet<NodeId, ()> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1), NodeId(2)], ()),
        (vec![NodeId(2), NodeId(3), NodeId(4)], ()),
    ]);
    assert_eq!(rs.incident_to_node(NodeId(0)), &[RelationId(0)]);
    assert_eq!(
        rs.incident_to_node(NodeId(2)),
        &[RelationId(0), RelationId(1)]
    );
    assert_eq!(rs.incident_to_node(NodeId(4)), &[RelationId(1)]);
    assert!(rs.has_incident_to_node(NodeId(0)));
    assert!(!rs.has_incident_to_node(NodeId(7)));
}

#[rstest]
fn test_var_relation_set_edge_incidence() {
    let rs: VarRelationSet<EdgeId, &str> = VarRelationSet::new(vec![
        (vec![EdgeId(0), EdgeId(2)], "a"),
        (vec![EdgeId(1), EdgeId(2)], "b"),
    ]);
    assert_eq!(
        rs.incident_to_edge(EdgeId(2)),
        &[RelationId(0), RelationId(1)]
    );
    assert_eq!(rs.incident_to_edge(EdgeId(0)), &[RelationId(0)]);
    assert!(rs.has_incident_to_edge(EdgeId(2)));
    assert!(!rs.has_incident_to_edge(EdgeId(5)));
    assert!(rs.incident_to_node(NodeId(0)).is_empty());
    assert!(!rs.has_incident_to_node(NodeId(0)));
}

#[rstest]
#[case::exact(vec![NodeId(0), NodeId(1), NodeId(2)], Some(RelationId(0)))]
#[case::reordered(vec![NodeId(2), NodeId(0), NodeId(1)], Some(RelationId(0)))]
#[case::second(vec![NodeId(3), NodeId(4)], Some(RelationId(1)))]
#[case::subset(vec![NodeId(0), NodeId(1)], None)]
#[case::superset(vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)], None)]
fn test_var_relation_set_coincident_to_node(
    #[case] query: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
) {
    let rs: VarRelationSet<NodeId, ()> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1), NodeId(2)], ()),
        (vec![NodeId(3), NodeId(4)], ()),
    ]);
    assert_eq!(
        query
            .first()
            .and_then(|&anchor| rs.coincident_to_node(anchor, &query)),
        expected,
    );
}

#[rstest]
#[case::reordered(NodeId(2), vec![NodeId(2), NodeId(2), NodeId(0)], Some(RelationId(0)), true)]
#[case::multiplicity(NodeId(2), vec![NodeId(0), NodeId(0), NodeId(2)], None, false)]
#[case::short(NodeId(2), vec![NodeId(2), NodeId(0)], None, false)]
#[case::long(NodeId(2), vec![NodeId(2), NodeId(0), NodeId(2), NodeId(2)], None, false)]
#[case::absent_anchor(NodeId(4), vec![NodeId(2), NodeId(2), NodeId(0)], None, true)]
fn test_var_relation_set_coincident_to_node_multiplicity(
    #[case] anchor: NodeId,
    #[case] query: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
    #[case] coincides: bool,
) {
    let relations = VarRelationSet::<NodeId, &str>::new(vec![
        (vec![NodeId(2), NodeId(0), NodeId(2)], "first"),
        (vec![NodeId(0), NodeId(2), NodeId(2)], "duplicate"),
        (vec![NodeId(3), NodeId(3), NodeId(3)], "other"),
    ]);
    assert_eq!(relations.coincident_to_node(anchor, &query), expected);
    assert_eq!(relations.is_coincident(RelationId(0), &query), coincides);
}

#[rstest]
#[case::empty_set(vec![], vec![NodeId(2), NodeId(2)], vec![(vec![NodeId(2), NodeId(2)], "added")])]
#[case::empty_row(vec![], vec![], vec![(vec![], "added")])]
#[case::adjacent_empty(vec![(vec![], "first"), (vec![], "second")], vec![], vec![(vec![], "first"), (vec![], "second"), (vec![], "added")])]
#[case::after_empty(vec![(vec![NodeId(0)], "first"), (vec![], "second")], vec![NodeId(2), NodeId(0)], vec![(vec![NodeId(0)], "first"), (vec![], "second"), (vec![NodeId(2), NodeId(0)], "added")])]
#[case::coinciding(vec![(vec![NodeId(2), NodeId(2)], "old")], vec![NodeId(2), NodeId(2)], vec![(vec![NodeId(2), NodeId(2)], "old"), (vec![NodeId(2), NodeId(2)], "added")])]
#[case::sparse(vec![(vec![NodeId(0)], "old")], vec![NodeId(u32::MAX), NodeId(0), NodeId(u32::MAX)], vec![(vec![NodeId(0)], "old"), (vec![NodeId(u32::MAX), NodeId(0), NodeId(u32::MAX)], "added")])]
#[case::empty_after_nonempty(vec![(vec![NodeId(1), NodeId(0)], "old")], vec![], vec![(vec![NodeId(1), NodeId(0)], "old"), (vec![], "added")])]
fn test_var_relation_set_add(
    #[case] entries: Vec<(Vec<NodeId>, &'static str)>,
    #[case] participants: Vec<NodeId>,
    #[case] expected: Vec<(Vec<NodeId>, &'static str)>,
) {
    let count = entries.len();
    let mut relations = VarRelationSet::new(entries);
    assert_eq!(
        relations.add(&participants, "added"),
        RelationId::from(count)
    );
    assert_eq!(relations.count(), count + 1);
    for (index, (row, data)) in expected.iter().enumerate() {
        assert_eq!(relations.participants(RelationId::from(index)), row);
        assert_eq!(relations.data(RelationId::from(index)), data);
    }
    for node in [NodeId(0), NodeId(1), NodeId(2), NodeId(u32::MAX)] {
        let incidence: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (row, _))| row.contains(&node))
            .map(|(i, _)| RelationId::from(i))
            .collect();
        assert_eq!(relations.incident_to_node(node), incidence);
        assert_eq!(relations.incident_to_edge(EdgeId(node.0)), &[]);
    }
    assert_eq!(relations.into_entries(), expected);
}

#[rstest]
fn test_var_relation_set_add_payload() {
    let mut relations = VarRelationSet::default();
    assert_eq!(
        relations.add(&[EdgeId(2), EdgeId(0)], NonCloneData(vec![7, 11])),
        RelationId(0)
    );
    assert_eq!(relations.add(&[], NonCloneData(vec![13])), RelationId(1));
    assert_eq!(
        relations.add(&[EdgeId(2), EdgeId(2)], NonCloneData(vec![17, 19])),
        RelationId(2)
    );
    assert_eq!(relations.incident_to_edge(EdgeId(0)), &[RelationId(0)]);
    assert_eq!(
        relations.incident_to_edge(EdgeId(2)),
        &[RelationId(0), RelationId(2)]
    );
    assert_eq!(relations.incident_to_node(NodeId(2)), &[]);
    assert_eq!(
        relations.into_entries(),
        vec![
            (vec![EdgeId(2), EdgeId(0)], NonCloneData(vec![7, 11])),
            (vec![], NonCloneData(vec![13])),
            (vec![EdgeId(2), EdgeId(2)], NonCloneData(vec![17, 19])),
        ]
    );
}

#[rstest]
#[case::first_empty(vec![RelationId(0)], vec![1, 2, 3, 4, 5])]
#[case::first_nonempty(vec![RelationId(1)], vec![0, 2, 3, 4, 5])]
#[case::middle_empty(vec![RelationId(2)], vec![0, 1, 3, 4, 5])]
#[case::last(vec![RelationId(5)], vec![0, 1, 2, 3, 4])]
#[case::unordered_repeated(vec![RelationId(4), RelationId(1), RelationId(4), RelationId(2)], vec![0, 3, 5])]
#[case::all_nonempty(vec![RelationId(1), RelationId(3), RelationId(5)], vec![0, 2, 4])]
#[case::all_empty(vec![RelationId(0), RelationId(2), RelationId(4)], vec![1, 3, 5])]
#[case::all(vec![RelationId(5), RelationId(3), RelationId(1), RelationId(0), RelationId(2), RelationId(4)], vec![])]
fn test_var_relation_set_tracked_remove(
    #[case] ids: Vec<RelationId>,
    #[case] survivors: Vec<usize>,
) {
    let entries = vec![
        (vec![], PositionLabels(vec![7])),
        (
            vec![NodeId(0), NodeId(1), NodeId(0)],
            PositionLabels(vec![11, 13, 17]),
        ),
        (vec![], PositionLabels(vec![19])),
        (vec![NodeId(2), NodeId(3)], PositionLabels(vec![23, 29])),
        (vec![], PositionLabels(vec![31])),
        (vec![NodeId(1)], PositionLabels(vec![37])),
    ];
    let expected: Vec<_> = survivors.iter().map(|&i| entries[i].clone()).collect();
    let mut relations = VarRelationSet::new(entries);
    let mut plain = relations.clone();
    plain.remove(&ids);
    let compaction = relations.tracked_remove(&ids);
    assert_eq!(relations, plain);
    assert_eq!(compaction.source_count(), 6);
    assert_eq!(compaction.result_count(), survivors.len());
    let removed: Vec<_> = (0..6)
        .filter(|i| !survivors.contains(i))
        .map(RelationId::from)
        .collect();
    assert_eq!(compaction.removed(), removed);
    for old in 0..=6 {
        assert_eq!(
            compaction.compact(RelationId::from(old)),
            survivors
                .iter()
                .position(|&i| i == old)
                .map(RelationId::from)
        );
    }
    assert_var_relation_rows(&relations, &expected);
    assert_eq!(relations.into_entries(), expected);
}

#[rstest]
#[case::empty(vec![])]
#[case::empty_rows(vec![(vec![], PositionLabels(vec![7])), (vec![], PositionLabels(vec![11]))])]
#[case::mixed(vec![(vec![], PositionLabels(vec![7])), (vec![NodeId(2), NodeId(0)], PositionLabels(vec![11, 13]))])]
fn test_var_relation_set_tracked_remove_identity(
    #[case] entries: Vec<(Vec<NodeId>, PositionLabels)>,
) {
    let mut relations = VarRelationSet::new(entries.clone());
    let original = relations.clone();
    assert_eq!(
        relations.tracked_remove(&[]),
        Compaction::identity(entries.len())
    );
    assert_eq!(relations, original);
    assert_var_relation_rows(&relations, &entries);
}

#[rstest]
fn test_var_relation_set_tracked_remove_payload() {
    let mut relations = VarRelationSet::default();
    relations.add(&[EdgeId(2), EdgeId(0)], NonCloneData(vec![7, 11]));
    relations.add(&[], NonCloneData(vec![13]));
    relations.add(
        &[EdgeId(2), EdgeId(2), EdgeId(3)],
        NonCloneData(vec![17, 19, 23]),
    );
    relations.add(&[], NonCloneData(vec![29]));
    let compaction = relations.tracked_remove(&[RelationId(0), RelationId(3)]);
    assert_eq!(compaction.source_count(), 4);
    assert_eq!(compaction.result_count(), 2);
    assert_eq!(compaction.removed(), &[RelationId(0), RelationId(3)]);
    assert_eq!(relations.participants(RelationId(0)), &[]);
    assert_eq!(
        relations.participants(RelationId(1)),
        &[EdgeId(2), EdgeId(2), EdgeId(3)]
    );
    assert_eq!(relations.incident_to_edge(EdgeId(0)), &[]);
    assert_eq!(relations.incident_to_edge(EdgeId(2)), &[RelationId(1)]);
    assert_eq!(relations.incident_to_edge(EdgeId(3)), &[RelationId(1)]);
    assert_eq!(
        relations.into_entries(),
        vec![
            (vec![], NonCloneData(vec![13])),
            (
                vec![EdgeId(2), EdgeId(2), EdgeId(3)],
                NonCloneData(vec![17, 19, 23])
            ),
        ]
    );
}

#[rstest]
#[case::empty(vec![], vec![RelationId(0)])]
#[case::end(vec![(vec![], "empty")], vec![RelationId(1)])]
#[case::mixed(vec![(vec![NodeId(0)], "first"), (vec![], "empty"), (vec![NodeId(2)], "last")], vec![RelationId(0), RelationId(3)])]
#[case::distant(vec![(vec![NodeId(0)], "first")], vec![RelationId(u32::MAX), RelationId(0)])]
fn test_var_relation_set_tracked_remove_error(
    #[case] entries: Vec<(Vec<NodeId>, &'static str)>,
    #[case] ids: Vec<RelationId>,
) {
    let original = VarRelationSet::new(entries);
    let mut relations = original.clone();
    let panic = catch_unwind(AssertUnwindSafe(|| relations.tracked_remove(&ids))).unwrap_err();
    let message = panic.downcast_ref::<String>().unwrap();
    assert!(message.starts_with("removed relations belong to the source set"));
    assert_eq!(relations, original);
    for node in [NodeId(0), NodeId(2)] {
        assert_eq!(
            relations.incident_to_node(node),
            original.incident_to_node(node)
        );
    }
    let panic = catch_unwind(AssertUnwindSafe(|| relations.remove(&ids))).unwrap_err();
    assert_eq!(panic.downcast_ref::<String>(), Some(message));
    assert_eq!(relations, original);
}

#[rstest]
#[case::interleaved(
    vec![(vec![], "survivor"), (vec![NodeId(3), NodeId(1), NodeId(3)], "changed")],
    5, vec![RelationId(0), RelationId(2), RelationId(4)],
    vec![(RelationId(4), vec![NodeId(4)], "last"), (RelationId(0), vec![NodeId(0), NodeId(1)], "first"), (RelationId(2), vec![], "middle")],
    vec![(vec![NodeId(0), NodeId(1)], "first"), (vec![], "survivor"), (vec![], "middle"), (vec![NodeId(3), NodeId(1), NodeId(3)], "changed"), (vec![NodeId(4)], "last")],
)]
#[case::all(
    vec![], 3, vec![RelationId(0), RelationId(1), RelationId(2)],
    vec![(RelationId(2), vec![NodeId(2), NodeId(1), NodeId(2), NodeId(0), NodeId(2)], "last"), (RelationId(0), vec![], "empty"), (RelationId(1), vec![NodeId(1)], "middle")],
    vec![(vec![], "empty"), (vec![NodeId(1)], "middle"), (vec![NodeId(2), NodeId(1), NodeId(2), NodeId(0), NodeId(2)], "last")],
)]
#[case::empty_rows(
    vec![(vec![], "survivor")], 3, vec![RelationId(0), RelationId(2)],
    vec![(RelationId(2), vec![], "last"), (RelationId(0), vec![], "first")],
    vec![(vec![], "first"), (vec![], "survivor"), (vec![], "last")],
)]
#[case::coincident(
    vec![(vec![NodeId(0), NodeId(1)], "survivor")], 2, vec![RelationId(0)],
    vec![(RelationId(0), vec![NodeId(0), NodeId(1)], "restored")],
    vec![(vec![NodeId(0), NodeId(1)], "restored"), (vec![NodeId(0), NodeId(1)], "survivor")],
)]
fn test_var_relation_set_restore(
    #[case] entries: Vec<(Vec<NodeId>, &'static str)>,
    #[case] count: usize,
    #[case] ids: Vec<RelationId>,
    #[case] removed: Vec<(RelationId, Vec<NodeId>, &'static str)>,
    #[case] expected: Vec<(Vec<NodeId>, &'static str)>,
) {
    let mut relations = VarRelationSet::new(entries);
    relations.restore(&Compaction::new(count, ids).unwrap(), removed);
    assert_eq!(relations.clone().into_entries(), expected);
    for node in 0..6 {
        let incidence: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (row, _))| row.contains(&NodeId(node)))
            .map(|(id, _)| RelationId::from(id))
            .collect();
        assert_eq!(relations.incident_to_node(NodeId(node)), incidence);
        assert_eq!(relations.incident_to_edge(EdgeId(node)), &[]);
    }
}

#[rstest]
#[case::empty(VarRelationSet::<NodeId, &str>::new(vec![]))]
#[case::empty_rows(VarRelationSet::<NodeId, &str>::new(vec![(vec![], "first"), (vec![], "second")]))]
#[case::nonempty(VarRelationSet::new(vec![(vec![NodeId(0), NodeId(2)], "first")]))]
fn test_var_relation_set_restore_identity(#[case] input: VarRelationSet<NodeId, &'static str>) {
    let mut relations = input.clone();
    relations.restore(&Compaction::identity(input.count()), vec![]);
    assert_eq!(relations, input);
}

#[rstest]
#[case::edges(vec![EdgeId(2), EdgeId(2)], vec![EdgeId(0), EdgeId(2)])]
#[case::empty_rows(vec![], vec![])]
#[case::empty_survivor(vec![], vec![EdgeId(0), EdgeId(2), EdgeId(0)])]
#[case::empty_restored(vec![EdgeId(2)], vec![])]
fn test_var_relation_set_restore_payload(
    #[case] surviving: Vec<EdgeId>,
    #[case] removed: Vec<EdgeId>,
) {
    let mut relations = VarRelationSet::default();
    relations.add(&surviving, NonCloneData(vec![7, 11]));
    relations.restore(
        &Compaction::new(2, vec![RelationId(0)]).unwrap(),
        vec![(RelationId(0), removed.clone(), NonCloneData(vec![13, 17]))],
    );
    let expected = vec![
        (removed, NonCloneData(vec![13, 17])),
        (surviving, NonCloneData(vec![7, 11])),
    ];
    for edge in 0..4 {
        let incidence: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (row, _))| row.contains(&EdgeId(edge)))
            .map(|(id, _)| RelationId::from(id))
            .collect();
        assert_eq!(relations.incident_to_edge(EdgeId(edge)), incidence);
    }
    assert_eq!(relations.into_entries(), expected);
}

#[rstest]
#[case::missing(3, vec![RelationId(0), RelationId(2)], vec![])]
#[case::duplicate(3, vec![RelationId(0), RelationId(2)], vec![RelationId(0), RelationId(0)])]
#[case::survivor_slot(3, vec![RelationId(0), RelationId(2)], vec![RelationId(1), RelationId(2)])]
#[case::out_of_range(3, vec![RelationId(0), RelationId(2)], vec![RelationId(u32::MAX)])]
#[case::too_many_survivors(1, vec![RelationId(0)], vec![RelationId(0)])]
#[case::too_few_survivors(8, vec![RelationId(0)], vec![RelationId(0)])]
#[case::identity_with_entries(1, vec![], vec![RelationId(0)])]
#[case::oversized(usize::MAX, vec![RelationId(0)], vec![])]
fn test_var_relation_set_restore_malformed(
    #[case] count: usize,
    #[case] ids: Vec<RelationId>,
    #[case] removed: Vec<RelationId>,
) {
    let mut relations = VarRelationSet::default();
    relations.add(&[NodeId(0)], NonCloneData(vec![7]));
    relations.restore(
        &Compaction::new(count, ids).unwrap(),
        removed
            .into_iter()
            .map(|id| (id, vec![NodeId(u32::MAX)], NonCloneData(vec![11])))
            .collect(),
    );
}

#[rstest]
#[case::nodes(vec![(vec![NodeId(0), NodeId(1), NodeId(0)], "first"), (vec![], "empty"), (vec![NodeId(2)], "second")], vec![(vec![NodeId(1), NodeId(3), NodeId(1)], "first"), (vec![], "empty"), (vec![NodeId(4)], "second")])]
#[case::edges(vec![(vec![EdgeId(0), EdgeId(1), EdgeId(0)], "first")], vec![(vec![EdgeId(0), EdgeId(2), EdgeId(0)], "first")])]
#[case::empty_rows(vec![(Vec::<NodeId>::new(), "first")], vec![(vec![], "first")])]
#[case::empty(Vec::<(Vec<NodeId>, &str)>::new(), vec![])]
fn test_var_relation_set_restore_participants<P: RelationParticipant + Debug>(
    #[case] entries: Vec<(Vec<P>, &'static str)>,
    #[case] expected: Vec<(Vec<P>, &'static str)>,
) {
    let compaction = GraphCompaction::new(
        Compaction::new(5, vec![NodeId(0), NodeId(2)]).unwrap(),
        Compaction::new(4, vec![EdgeId(1)]).unwrap(),
    );
    let mut relations = VarRelationSet::new(entries);
    relations.restore_participants(&compaction);
    assert_eq!(relations.clone().into_entries(), expected);
    for id in 0..6 {
        let nodes: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (row, _))| row.iter().any(|p| p.refs().node == Some(NodeId(id))))
            .map(|(id, _)| RelationId::from(id))
            .collect();
        let edges: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (row, _))| row.iter().any(|p| p.refs().edge == Some(EdgeId(id))))
            .map(|(id, _)| RelationId::from(id))
            .collect();
        assert_eq!(relations.incident_to_node(NodeId(id)), nodes);
        assert_eq!(relations.incident_to_edge(EdgeId(id)), edges);
    }
}

#[rstest]
fn test_var_relation_set_restore_participants_references() {
    let compacted = vec![
        References {
            node: Some(NodeId(0)),
            edge: Some(EdgeId(1)),
            label: 1,
        },
        References {
            node: Some(NodeId(0)),
            edge: None,
            label: 2,
        },
        References {
            node: None,
            edge: Some(EdgeId(1)),
            label: 3,
        },
        References {
            node: None,
            edge: None,
            label: 4,
        },
        References {
            node: Some(NodeId(0)),
            edge: Some(EdgeId(1)),
            label: 5,
        },
    ];
    let mut relations = VarRelationSet::default();
    relations.add(&compacted, NonCloneData(vec![7, 11]));
    let payload = relations.data(RelationId(0)) as *const NonCloneData;
    relations.restore_participants(&GraphCompaction::new(
        Compaction::new(3, vec![NodeId(0)]).unwrap(),
        Compaction::new(4, vec![EdgeId(1)]).unwrap(),
    ));
    assert_eq!(
        relations.data(RelationId(0)) as *const NonCloneData,
        payload
    );
    assert_eq!(
        relations.participants(RelationId(0)),
        &[
            References {
                node: Some(NodeId(1)),
                edge: Some(EdgeId(2)),
                label: 1
            },
            References {
                node: Some(NodeId(1)),
                edge: None,
                label: 2
            },
            References {
                node: None,
                edge: Some(EdgeId(2)),
                label: 3
            },
            References {
                node: None,
                edge: None,
                label: 4
            },
            References {
                node: Some(NodeId(1)),
                edge: Some(EdgeId(2)),
                label: 5
            },
        ]
    );
    assert_eq!(relations.data(RelationId(0)), &NonCloneData(vec![7, 11]));
    assert_eq!(relations.incident_to_node(NodeId(1)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_edge(EdgeId(2)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_node(NodeId(0)), &[]);
    assert_eq!(relations.incident_to_edge(EdgeId(1)), &[]);
}

#[rstest]
#[case::empty(VarRelationSet::<NodeId, &str>::new(vec![]))]
#[case::empty_rows(VarRelationSet::<NodeId, &str>::new(vec![(vec![], "first")]))]
#[case::nonempty(VarRelationSet::new(vec![(vec![NodeId(0), NodeId(2)], "first")]))]
fn test_var_relation_set_restore_participants_identity(
    #[case] input: VarRelationSet<NodeId, &'static str>,
) {
    let mut relations = input.clone();
    relations.restore_participants(&GraphCompaction::new(
        Compaction::identity(3),
        Compaction::identity(0),
    ));
    assert_eq!(relations, input);
}

#[rstest]
#[case::node_domain(3, 3, Some(NodeId(2)), None)]
#[case::edge_domain(3, 3, None, Some(EdgeId(2)))]
#[case::both_domains(3, 3, Some(NodeId(u32::MAX)), Some(EdgeId(u32::MAX)))]
#[case::oversized_nodes(usize::MAX, 3, Some(NodeId(0)), None)]
#[case::oversized_edges(3, usize::MAX, None, Some(EdgeId(0)))]
fn test_var_relation_set_restore_participants_malformed(
    #[case] nodes: usize,
    #[case] edges: usize,
    #[case] node: Option<NodeId>,
    #[case] edge: Option<EdgeId>,
) {
    let mut relations = VarRelationSet::default();
    relations.add(
        &[
            References {
                node: Some(NodeId(0)),
                edge: Some(EdgeId(0)),
                label: 1,
            },
            References {
                node,
                edge,
                label: 2,
            },
        ],
        NonCloneData(vec![7]),
    );
    relations.restore_participants(&GraphCompaction::new(
        Compaction::new(nodes, vec![NodeId(0)]).unwrap(),
        Compaction::new(edges, vec![EdgeId(0)]).unwrap(),
    ));
}

#[rstest]
#[case::second_relation(RelationId(1), vec![ParticipantPosition(2), ParticipantPosition(0), ParticipantPosition(1)],
    vec![NodeId(4), NodeId(2), NodeId(3)])]
#[case::last_relation(RelationId(2), vec![ParticipantPosition(1), ParticipantPosition(0)], vec![NodeId(6), NodeId(5)])]
fn test_var_relation_set_permute_participants(
    #[case] id: RelationId,
    #[case] order: Vec<ParticipantPosition>,
    #[case] expected: Vec<NodeId>,
) {
    let mut rs: VarRelationSet<NodeId, &str> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], "a"),
        (vec![NodeId(2), NodeId(3), NodeId(4)], "b"),
        (vec![NodeId(5), NodeId(6)], "c"),
    ]);
    let incidence_before: Vec<Vec<RelationId>> = (0..7)
        .map(|i| rs.incident_to_node(NodeId(i)).to_vec())
        .collect();

    rs.permute_participants(id, &order);

    assert_eq!(rs.participants(id), expected.as_slice());
    assert_eq!(rs.data(id), &["a", "b", "c"][id.index()]);
    for other in [RelationId(0), RelationId(1), RelationId(2)] {
        if other != id {
            let stored: Vec<Vec<NodeId>> = vec![
                vec![NodeId(0), NodeId(1)],
                vec![NodeId(2), NodeId(3), NodeId(4)],
                vec![NodeId(5), NodeId(6)],
            ];
            assert_eq!(rs.participants(other), stored[other.index()].as_slice());
        }
    }
    let incidence_after: Vec<Vec<RelationId>> = (0..7)
        .map(|i| rs.incident_to_node(NodeId(i)).to_vec())
        .collect();
    assert_eq!(incidence_after, incidence_before);
}

#[rstest]
fn test_var_relation_set_permute_participants_identity() {
    let input: VarRelationSet<NodeId, &str> = VarRelationSet::new(vec![
        (vec![NodeId(2), NodeId(0), NodeId(1)], "a"),
        (vec![NodeId(4), NodeId(3)], "b"),
    ]);
    let mut permuted = input.clone();
    permuted.permute_participants(
        RelationId(0),
        &[
            ParticipantPosition(0),
            ParticipantPosition(1),
            ParticipantPosition(2),
        ],
    );
    assert_eq!(permuted, input);
}

#[rstest]
#[case::order_too_short(vec![ParticipantPosition(0), ParticipantPosition(1)])]
#[case::position_out_of_range(vec![ParticipantPosition(0), ParticipantPosition(1), ParticipantPosition(9)])]
#[case::position_repeated(vec![ParticipantPosition(1), ParticipantPosition(1), ParticipantPosition(0)])]
#[should_panic(expected = "permute")]
fn test_var_relation_set_permute_participants_error(#[case] order: Vec<ParticipantPosition>) {
    let mut rs: VarRelationSet<NodeId, &str> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], "a"),
        (vec![NodeId(2), NodeId(3), NodeId(4)], "b"),
    ]);
    rs.permute_participants(RelationId(1), &order);
}

#[rstest]
#[case::grow_first(RelationId(0), vec![NodeId(5), NodeId(1), NodeId(5), NodeId(3)])]
#[case::shrink_first(RelationId(0), vec![NodeId(5)])]
#[case::clear_first(RelationId(0), vec![])]
#[case::reorder_middle(RelationId(1), vec![NodeId(2), NodeId(1), NodeId(1)])]
#[case::shared_middle(RelationId(1), vec![NodeId(1), NodeId(3), NodeId(5), NodeId(3)])]
#[case::clear_middle(RelationId(1), vec![])]
#[case::grow_empty(RelationId(2), vec![NodeId(5), NodeId(5)])]
#[case::grow_last(RelationId(3), vec![NodeId(u32::MAX), NodeId(3), NodeId(u32::MAX)])]
#[case::clear_last(RelationId(3), vec![])]
fn test_var_relation_set_replace_participants(
    var_relation_set_mutation_entries: Vec<(Vec<NodeId>, PositionLabels)>,
    #[case] id: RelationId,
    #[case] participants: Vec<NodeId>,
) {
    let mut entries = var_relation_set_mutation_entries;
    let mut relations = VarRelationSet::new(entries.clone());
    relations.replace_participants(id, &participants);
    entries[id.index()].0 = participants;
    assert_var_relation_rows(&relations, &entries);
}

#[rstest]
#[case::first(RelationId(0))]
#[case::middle(RelationId(1))]
#[case::empty(RelationId(2))]
#[case::last(RelationId(3))]
fn test_var_relation_set_replace_participants_identity(
    var_relation_set_mutation_entries: Vec<(Vec<NodeId>, PositionLabels)>,
    #[case] id: RelationId,
) {
    let entries = var_relation_set_mutation_entries;
    let input = VarRelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participants(id, &entries[id.index()].0);
    assert_eq!(relations, input);
    assert_var_relation_rows(&relations, &entries);
}

#[rstest]
fn test_var_relation_set_replace_participants_edges() {
    let mut relations = VarRelationSet::new(vec![
        (vec![EdgeId(0), EdgeId(1), EdgeId(0)], 7),
        (vec![EdgeId(1)], 11),
    ]);
    relations.replace_participants(RelationId(0), &[EdgeId(3), EdgeId(3)]);
    assert_eq!(relations.incident_to_edge(EdgeId(0)), &[]);
    assert_eq!(relations.incident_to_edge(EdgeId(1)), &[RelationId(1)]);
    assert_eq!(relations.incident_to_edge(EdgeId(3)), &[RelationId(0)]);
    assert_eq!(
        relations.coincident_to_edge(EdgeId(3), &[EdgeId(3), EdgeId(3)]),
        Some(RelationId(0))
    );
    assert_eq!(
        relations.into_entries(),
        vec![(vec![EdgeId(3), EdgeId(3)], 7), (vec![EdgeId(1)], 11)]
    );
}

#[rstest]
#[case::past_end(RelationId(4))]
#[case::maximum(RelationId(u32::MAX))]
#[should_panic]
fn test_var_relation_set_replace_participants_error(
    var_relation_set_mutation_entries: Vec<(Vec<NodeId>, PositionLabels)>,
    #[case] id: RelationId,
) {
    VarRelationSet::new(var_relation_set_mutation_entries).replace_participants(id, &[]);
}

#[rstest]
#[should_panic]
fn test_var_relation_set_replace_participants_empty() {
    VarRelationSet::<NodeId, NonCloneData>::default().replace_participants(RelationId(0), &[]);
}

#[rstest]
#[case::first(RelationId(0), ParticipantPosition(0), NodeId(5), vec![NodeId(5), NodeId(1), NodeId(0)])]
#[case::shared(RelationId(1), ParticipantPosition(0), NodeId(3), vec![NodeId(3), NodeId(1), NodeId(2)])]
#[case::middle(RelationId(1), ParticipantPosition(1), NodeId(5), vec![NodeId(1), NodeId(5), NodeId(2)])]
#[case::last(RelationId(3), ParticipantPosition(1), NodeId(u32::MAX), vec![NodeId(3), NodeId(u32::MAX)])]
fn test_var_relation_set_replace_participant(
    var_relation_set_mutation_entries: Vec<(Vec<NodeId>, PositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
    #[case] participant: NodeId,
    #[case] expected: Vec<NodeId>,
) {
    let mut entries = var_relation_set_mutation_entries;
    let mut relations = VarRelationSet::new(entries.clone());
    relations.replace_participant(id, position, participant);
    entries[id.index()].0 = expected;
    assert_var_relation_rows(&relations, &entries);
}

#[rstest]
#[case::first(ParticipantPosition(0))]
#[case::middle(ParticipantPosition(1))]
#[case::last(ParticipantPosition(2))]
fn test_var_relation_set_replace_participant_identity(
    var_relation_set_mutation_entries: Vec<(Vec<NodeId>, PositionLabels)>,
    #[case] position: ParticipantPosition,
) {
    let entries = var_relation_set_mutation_entries;
    let input = VarRelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participant(RelationId(1), position, entries[1].0[position.index()]);
    assert_eq!(relations, input);
    assert_var_relation_rows(&relations, &entries);
}

#[rstest]
fn test_var_relation_set_replace_participant_edges() {
    let mut relations =
        VarRelationSet::new(vec![(vec![EdgeId(3), EdgeId(3)], 7), (vec![EdgeId(1)], 11)]);
    relations.replace_participant(RelationId(0), ParticipantPosition(0), EdgeId(4));
    assert_eq!(relations.incident_to_edge(EdgeId(1)), &[RelationId(1)]);
    assert_eq!(relations.incident_to_edge(EdgeId(3)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_edge(EdgeId(4)), &[RelationId(0)]);
    assert_eq!(
        relations.coincident_to_edge(EdgeId(4), &[EdgeId(3), EdgeId(4)]),
        Some(RelationId(0))
    );
    assert_eq!(
        relations.into_entries(),
        vec![(vec![EdgeId(4), EdgeId(3)], 7), (vec![EdgeId(1)], 11)]
    );
}

#[rstest]
#[case::invalid_id(RelationId(4), ParticipantPosition(0))]
#[case::maximum_id(RelationId(u32::MAX), ParticipantPosition(0))]
#[case::at_length(RelationId(1), ParticipantPosition(3))]
#[case::maximum_position(RelationId(1), ParticipantPosition(u32::MAX))]
#[case::empty_row(RelationId(2), ParticipantPosition(0))]
#[should_panic]
fn test_var_relation_set_replace_participant_error(
    var_relation_set_mutation_entries: Vec<(Vec<NodeId>, PositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
) {
    VarRelationSet::new(var_relation_set_mutation_entries).replace_participant(
        id,
        position,
        NodeId(5),
    );
}

#[rstest]
#[should_panic]
fn test_var_relation_set_replace_participant_empty() {
    VarRelationSet::<NodeId, NonCloneData>::default().replace_participant(
        RelationId(0),
        ParticipantPosition(0),
        NodeId(5),
    );
}

#[rstest]
#[case::first(RelationId(0), ParticipantPosition(0), NodeId(5), vec![NodeId(5), NodeId(0), NodeId(1), NodeId(0)])]
#[case::middle(RelationId(1), ParticipantPosition(1), NodeId(3), vec![NodeId(1), NodeId(3), NodeId(1), NodeId(2)])]
#[case::append(RelationId(1), ParticipantPosition(3), NodeId(5), vec![NodeId(1), NodeId(1), NodeId(2), NodeId(5)])]
#[case::duplicate(RelationId(1), ParticipantPosition(0), NodeId(1), vec![NodeId(1), NodeId(1), NodeId(1), NodeId(2)])]
#[case::empty_row(RelationId(2), ParticipantPosition(0), NodeId(5), vec![NodeId(5)])]
#[case::last(RelationId(3), ParticipantPosition(2), NodeId(u32::MAX), vec![NodeId(3), NodeId(4), NodeId(u32::MAX)])]
fn test_var_relation_set_insert_participant(
    var_relation_set_mutation_entries: Vec<(Vec<NodeId>, PositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
    #[case] participant: NodeId,
    #[case] expected: Vec<NodeId>,
) {
    let mut entries = var_relation_set_mutation_entries;
    let mut relations = VarRelationSet::new(entries.clone());
    relations.insert_participant(id, position, participant);
    entries[id.index()].0 = expected;
    assert_var_relation_rows(&relations, &entries);
}

#[rstest]
fn test_var_relation_set_insert_participant_edges() {
    let mut relations =
        VarRelationSet::new(vec![(vec![EdgeId(4), EdgeId(3)], 7), (vec![EdgeId(1)], 11)]);
    relations.insert_participant(RelationId(1), ParticipantPosition(1), EdgeId(4));
    assert_eq!(relations.incident_to_edge(EdgeId(1)), &[RelationId(1)]);
    assert_eq!(relations.incident_to_edge(EdgeId(3)), &[RelationId(0)]);
    assert_eq!(
        relations.incident_to_edge(EdgeId(4)),
        &[RelationId(0), RelationId(1)]
    );
    assert_eq!(
        relations.coincident_to_edge(EdgeId(4), &[EdgeId(4), EdgeId(1)]),
        Some(RelationId(1))
    );
    assert_eq!(
        relations.into_entries(),
        vec![
            (vec![EdgeId(4), EdgeId(3)], 7),
            (vec![EdgeId(1), EdgeId(4)], 11)
        ]
    );
}

#[rstest]
#[case::invalid_id(RelationId(4), ParticipantPosition(0))]
#[case::maximum_id(RelationId(u32::MAX), ParticipantPosition(0))]
#[case::past_length(RelationId(1), ParticipantPosition(4))]
#[case::maximum_position(RelationId(1), ParticipantPosition(u32::MAX))]
#[case::past_empty(RelationId(2), ParticipantPosition(1))]
#[should_panic]
fn test_var_relation_set_insert_participant_error(
    var_relation_set_mutation_entries: Vec<(Vec<NodeId>, PositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
) {
    VarRelationSet::new(var_relation_set_mutation_entries).insert_participant(
        id,
        position,
        NodeId(5),
    );
}

#[rstest]
#[should_panic]
fn test_var_relation_set_insert_participant_empty() {
    VarRelationSet::<NodeId, NonCloneData>::default().insert_participant(
        RelationId(0),
        ParticipantPosition(0),
        NodeId(5),
    );
}

#[rstest]
#[case::first(RelationId(0), ParticipantPosition(0), vec![NodeId(1), NodeId(0)])]
#[case::shared(RelationId(1), ParticipantPosition(0), vec![NodeId(1), NodeId(2)])]
#[case::middle(RelationId(1), ParticipantPosition(1), vec![NodeId(1), NodeId(2)])]
#[case::unique(RelationId(1), ParticipantPosition(2), vec![NodeId(1), NodeId(1)])]
#[case::last(RelationId(3), ParticipantPosition(1), vec![NodeId(3)])]
fn test_var_relation_set_remove_participant(
    var_relation_set_mutation_entries: Vec<(Vec<NodeId>, PositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
    #[case] expected: Vec<NodeId>,
) {
    let mut entries = var_relation_set_mutation_entries;
    let mut relations = VarRelationSet::new(entries.clone());
    relations.remove_participant(id, position);
    entries[id.index()].0 = expected;
    assert_var_relation_rows(&relations, &entries);
}

#[rstest]
fn test_var_relation_set_remove_participant_edges() {
    let mut relations =
        VarRelationSet::new(vec![(vec![EdgeId(4)], 7), (vec![EdgeId(1), EdgeId(4)], 11)]);
    relations.remove_participant(RelationId(0), ParticipantPosition(0));
    assert_eq!(relations.incident_to_edge(EdgeId(1)), &[RelationId(1)]);
    assert_eq!(relations.incident_to_edge(EdgeId(4)), &[RelationId(1)]);
    assert!(relations.is_coincident(RelationId(0), &[]));
    assert_eq!(
        relations.coincident_to_edge(EdgeId(4), &[EdgeId(4), EdgeId(1)]),
        Some(RelationId(1))
    );
    assert_eq!(
        relations.into_entries(),
        vec![(vec![], 7), (vec![EdgeId(1), EdgeId(4)], 11)]
    );
}

#[rstest]
#[case::invalid_id(RelationId(4), ParticipantPosition(0))]
#[case::maximum_id(RelationId(u32::MAX), ParticipantPosition(0))]
#[case::at_length(RelationId(1), ParticipantPosition(3))]
#[case::maximum_position(RelationId(1), ParticipantPosition(u32::MAX))]
#[case::empty_row(RelationId(2), ParticipantPosition(0))]
#[should_panic]
fn test_var_relation_set_remove_participant_error(
    var_relation_set_mutation_entries: Vec<(Vec<NodeId>, PositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
) {
    VarRelationSet::new(var_relation_set_mutation_entries).remove_participant(id, position);
}

#[rstest]
#[should_panic]
fn test_var_relation_set_remove_participant_empty() {
    VarRelationSet::<NodeId, NonCloneData>::default()
        .remove_participant(RelationId(0), ParticipantPosition(0));
}

#[rstest]
#[case::rows(VarRelationSet::new(vec![(vec![NodeId(2), NodeId(0)], vec![7, 11]), (vec![NodeId(2), NodeId(0)], vec![13, 17])]),
    VarRelationSet::new(vec![(vec![NodeId(1), NodeId(5)], vec![7, 11]), (vec![NodeId(1), NodeId(5)], vec![13, 17])]))]
fn test_var_relation_set_map(
    participant_correspondence: GraphCorrespondence,
    #[case] input: VarRelationSet<NodeId, Vec<u32>>,
    #[case] expected: VarRelationSet<NodeId, Vec<u32>>,
) {
    assert_eq!(input.map(&participant_correspondence), expected);
    assert_eq!(
        input.try_map(&participant_correspondence),
        Some(expected.clone())
    );
    let reverse = GraphCorrespondence::new(
        participant_correspondence.nodes().reverse(),
        participant_correspondence.edges().reverse(),
    );
    assert_eq!(expected.map(&reverse), input);
    let composed = participant_correspondence.compose(&reverse).unwrap();
    assert_eq!(input.map(&composed), input);
    assert_eq!(
        expected.incident_to_node(NodeId(1)),
        &[RelationId(0), RelationId(1)]
    );
}

#[rstest]
#[case::empty(VarRelationSet::new(vec![]))]
#[case::rows(VarRelationSet::new(vec![(vec![NodeId(2), NodeId(0)], vec![7, 11]), (vec![NodeId(2), NodeId(0)], vec![13, 17])]))]
fn test_var_relation_set_map_identity(#[case] input: VarRelationSet<NodeId, Vec<u32>>) {
    let identity = GraphCorrespondence::new(
        Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2), NodeId(3)], 4),
        Correspondence::from_images(&[EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)], 4),
    );
    assert_eq!(input.try_map(&identity), Some(input.clone()));
    assert_eq!(input.map(&identity), input);
}

#[rstest]
#[should_panic(expected = "correspondence must cover every participant reference")]
fn test_var_relation_set_map_error(participant_correspondence: GraphCorrespondence) {
    let node = 1;

    let input: VarRelationSet<NodeId, Vec<u32>> =
        VarRelationSet::new(vec![(vec![NodeId(node), NodeId(0)], vec![7, 11])]);
    input.map(&participant_correspondence);
}

#[rstest]
#[case::missing_node(1)]
#[case::outside_node(4)]
fn test_var_relation_set_try_map_error(
    participant_correspondence: GraphCorrespondence,
    #[case] node: u32,
) {
    let input: VarRelationSet<NodeId, Vec<u32>> = VarRelationSet::new(vec![
        (vec![NodeId(2), NodeId(0)], vec![7, 11]),
        (vec![NodeId(node), NodeId(0)], vec![13, 17]),
    ]);
    assert_eq!(input.try_map(&participant_correspondence), None);
}

#[rstest]
fn test_var_relation_set_remap() {
    let rs: VarRelationSet<EdgeId, PositionLabels> = VarRelationSet::new(vec![(
        vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        PositionLabels(vec![20, 21, 22]),
    )]);
    let remapping = GraphRemapping::new(
        Remapping::empty(),
        Remapping::new(vec![EdgeId(2), EdgeId(0), EdgeId(1)]).expect("permutation images"),
    );
    let out = rs.remap(&remapping);
    assert_eq!(
        out.participants(RelationId(0)),
        &[EdgeId(2), EdgeId(0), EdgeId(1)]
    );
    assert_eq!(out.data(RelationId(0)), &PositionLabels(vec![20, 21, 22]));
}

#[rstest]
#[case::covered(vec![EdgeId(2), EdgeId(0), EdgeId(1)], true)]
#[case::uncovered_edge(vec![EdgeId(1), EdgeId(0)], false)]
fn test_var_relation_set_try_remap(#[case] edges: Vec<EdgeId>, #[case] covered: bool) {
    let rs: VarRelationSet<EdgeId, PositionLabels> = VarRelationSet::new(vec![(
        vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        PositionLabels(vec![20, 21, 22]),
    )]);
    let remapping = GraphRemapping::new(
        Remapping::empty(),
        Remapping::new(edges).expect("permutation images"),
    );
    let expected = covered.then(|| rs.remap(&remapping));
    assert_eq!(rs.try_remap(&remapping), expected);
}

#[rstest]
#[case::empty(VarRelationSet::default())]
#[case::rows(
    VarRelationSet::new(vec![(vec![NodeId(0), NodeId(2), NodeId(4)], "keep"), (vec![NodeId(1), NodeId(3)], "drop")]),
)]
fn test_var_relation_set_compact_identity(#[case] input: VarRelationSet<NodeId, &'static str>) {
    let compaction = GraphCompaction::new(Compaction::identity(5), Compaction::empty());
    assert_eq!(input.compact(&compaction), input);
    assert_eq!(
        input.tracked_compact(&compaction),
        (input.clone(), Compaction::identity(input.count())),
    );
}

#[rstest]
#[case::partial(
    vec![NodeId(1)],
    VarRelationSet::new(vec![(vec![NodeId(0), NodeId(1), NodeId(3)], "keep")]),
    vec![RelationId(1)],
)]
#[case::all(
    vec![NodeId(0), NodeId(1)],
    VarRelationSet::default(),
    vec![RelationId(0), RelationId(1)],
)]
fn test_var_relation_set_tracked_compact(
    var_relation_set_compaction_input: VarRelationSet<NodeId, &'static str>,
    #[case] removed_nodes: Vec<NodeId>,
    #[case] expected: VarRelationSet<NodeId, &'static str>,
    #[case] removed_relations: Vec<RelationId>,
) {
    let input = var_relation_set_compaction_input;
    let compaction = GraphCompaction::new(
        Compaction::new(5, removed_nodes).unwrap(),
        Compaction::empty(),
    );
    let (output, witness) = input.tracked_compact(&compaction);
    assert_eq!(input.compact(&compaction), expected);
    assert_eq!(output, expected);
    assert_eq!(
        witness,
        Compaction::new(2, removed_relations.clone()).unwrap()
    );
    let survivors = (0..2)
        .map(RelationId)
        .filter(|id| !removed_relations.contains(id))
        .collect::<Vec<_>>();
    for (idx, &old) in survivors.iter().enumerate() {
        assert_eq!(witness.compact(old), Some(RelationId::from(idx)));
    }
}

/// Paired rows retain the left frame; the combination sees both original payload frames.
#[rstest]
#[case::both_frames_kept(
    VarRelationSet::<NodeId, PositionLabels>::new(
        vec![(vec![NodeId(7), NodeId(3)], PositionLabels(vec![70, 30]))]),
    VarRelationSet::<NodeId, PositionLabels>::new(
        vec![(vec![NodeId(3), NodeId(7)], PositionLabels(vec![30, 70]))]),
    [NodeId(7), NodeId(3)],
    PositionLabels(vec![70, 30]),
    PositionLabels(vec![30, 70]),
)]
fn test_var_relation_set_tracked_pushout_coincidence_frame(
    #[case] left: VarRelationSet<NodeId, PositionLabels>,
    #[case] right: VarRelationSet<NodeId, PositionLabels>,
    #[case] expected_frame: [NodeId; 2],
    #[case] expected_left_seen: PositionLabels,
    #[case] expected_right_seen: PositionLabels,
) {
    let mut seen = None;
    let (object, correspondence) = left
        .tracked_pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident_to_node(n, q)),
            |(_, a), (_, b)| {
                seen = Some((a.clone(), b.clone()));
                Some(a.clone())
            },
        )
        .expect("combine never rejects here");
    assert_eq!(
        correspondence,
        RelationPushoutCorrespondence {
            left: Correspondence::new(vec![(RelationId(0), RelationId(0))], 1, 1).unwrap(),
            right: Correspondence::new(vec![(RelationId(0), RelationId(0))], 1, 1).unwrap(),
        }
    );

    assert_eq!(object.count(), 1, "the two entries coincide");
    assert_eq!(
        object.participants(RelationId(0)),
        expected_frame.as_slice()
    );
    assert_eq!(
        seen,
        Some((expected_left_seen, expected_right_seen)),
        "the two payloads as `combine` received them"
    );
}

#[rstest]
fn test_var_relation_set_tracked_pushout() {
    let left =
        VarRelationSet::<NodeId, i32>::new(vec![(vec![NodeId(0), NodeId(1), NodeId(2)], 10)]);
    let right = VarRelationSet::<NodeId, i32>::new(vec![
        (vec![NodeId(0), NodeId(1), NodeId(2)], 5),
        (vec![NodeId(3), NodeId(4)], 20),
    ]);
    let (object, glue) = left
        .tracked_pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident_to_node(n, q)),
            |(_, a), (_, b)| Some(a + b),
        )
        .expect("no ⊥");
    assert_eq!(
        left.pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident_to_node(n, q)),
            |(_, a), (_, b)| Some(a + b),
        ),
        Some(object.clone()),
    );
    assert_eq!(
        object,
        VarRelationSet::new(vec![
            (vec![NodeId(0), NodeId(1), NodeId(2)], 15),
            (vec![NodeId(3), NodeId(4)], 20)
        ])
    );
    assert_eq!(
        glue,
        RelationPushoutCorrespondence {
            left: Correspondence::new(vec![(RelationId(0), RelationId(0))], 1, 2,).unwrap(),
            right: Correspondence::new(
                vec![
                    (RelationId(0), RelationId(0)),
                    (RelationId(1), RelationId(1))
                ],
                2,
                2,
            )
            .unwrap(),
        },
    );
    assert_eq!(
        left.tracked_pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident_to_node(n, q)),
            |_, _| None,
        ),
        None
    );
    assert_eq!(
        left.pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident_to_node(n, q)),
            |_, _| None,
        ),
        None
    );
}

/// Repeated callback targets violate the injective correspondence required by pushout.
#[rstest]
#[should_panic(expected = "correspondence images must be unique")]
fn test_var_relation_set_tracked_pushout_repeated_coincidence() {
    let left: VarRelationSet<NodeId, PositionLabels> = VarRelationSet::new(vec![(
        vec![NodeId(0), NodeId(1)],
        PositionLabels(vec![1, 1]),
    )]);
    let right: VarRelationSet<NodeId, PositionLabels> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], PositionLabels(vec![2, 2])),
        (vec![NodeId(0), NodeId(1)], PositionLabels(vec![4, 4])),
    ]);

    left.tracked_pushout(
        &right,
        |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident_to_node(n, q)),
        |(_, a), (_, b)| {
            Some(PositionLabels(
                a.0.iter().zip(&b.0).map(|(x, y)| x + y).collect(),
            ))
        },
    );
}

#[rstest]
#[case::combined(Some(15))]
#[case::incompatible(None)]
fn test_var_relation_set_tracked_pullback(#[case] combined: Option<i32>) {
    let left: VarRelationSet<NodeId, i32> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], 10),
        (vec![NodeId(2), NodeId(3)], 20),
    ]);
    let right: VarRelationSet<NodeId, i32> = VarRelationSet::new(vec![
        (vec![NodeId(4), NodeId(5)], 30),
        (vec![NodeId(0), NodeId(1)], 5),
    ]);
    let result = left.tracked_pullback(
        &right,
        |set, parts: &[NodeId]| {
            parts
                .first()
                .and_then(|&id| set.coincident_to_node(id, parts))
        },
        |(_, a), (_, b)| combined.map(|_| a + b),
    );
    let plain = left.pullback(
        &right,
        |set, parts: &[NodeId]| {
            parts
                .first()
                .and_then(|&id| set.coincident_to_node(id, parts))
        },
        |(_, a), (_, b)| combined.map(|_| a + b),
    );
    let expected = combined.map(|value| {
        (
            VarRelationSet::new(vec![(vec![NodeId(0), NodeId(1)], value)]),
            RelationPullbackCorrespondence {
                left: Correspondence::new(vec![(RelationId(0), RelationId(0))], 1, 2).unwrap(),
                right: Correspondence::new(vec![(RelationId(0), RelationId(1))], 1, 2).unwrap(),
            },
        )
    });
    assert_eq!(plain, expected.as_ref().map(|(object, _)| object.clone()));
    assert_eq!(result, expected);
}

#[rstest]
fn test_var_relation_set_default() {
    let rs = VarRelationSet::<NodeId, ()>::default();
    assert_eq!(rs.count(), 0);
    assert!(!rs.has_incident_to_node(NodeId(0)));
}

#[rstest]
fn test_var_relation_set_hash() {
    let entries = vec![
        (vec![NodeId(2), NodeId(0)], "first"),
        (vec![NodeId(4), NodeId(3), NodeId(1)], "second"),
    ];
    let left: VarRelationSet<NodeId, &str> = VarRelationSet::new(entries.clone());
    let right: VarRelationSet<NodeId, &str> = VarRelationSet::new(entries);
    assert_eq!(left, right);
    assert_eq!(hash(&left), hash(&right));
}
