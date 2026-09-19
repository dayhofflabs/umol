use std::fmt::Debug;
use std::panic::{catch_unwind, AssertUnwindSafe};

use pretty_assertions::assert_eq;
use rstest::{fixture, rstest};

use crate::{
    Compaction, Correspondence, EdgeId, FixedFixedBirelationSet, GraphCompaction,
    GraphCorrespondence, GraphRemapping, NodeId, ParticipantPosition, ParticipantRefs, RelationId,
    RelationParticipant, RelationPullbackCorrespondence, RelationPushoutCorrespondence, Remapping,
};

#[fixture]
fn participant_correspondence() -> GraphCorrespondence {
    GraphCorrespondence::new(
        Correspondence::new(vec![(NodeId(0), NodeId(5)), (NodeId(2), NodeId(1))], 4, 6).unwrap(),
        Correspondence::new(vec![(EdgeId(0), EdgeId(6)), (EdgeId(2), EdgeId(3))], 4, 7).unwrap(),
    )
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct BiPositionLabels {
    factor_1: Vec<u32>,
    factor_2: Vec<u32>,
}

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
fn fixed_fixed_birelation_set_mutation_entries() -> Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>
{
    vec![
        (
            [NodeId(0), NodeId(1)],
            [NodeId(1), NodeId(2), NodeId(1)],
            BiPositionLabels {
                factor_1: vec![7, 11],
                factor_2: vec![13, 17, 19],
            },
        ),
        (
            [NodeId(2), NodeId(3)],
            [NodeId(3), NodeId(4), NodeId(3)],
            BiPositionLabels {
                factor_1: vec![23, 29],
                factor_2: vec![31, 37, 41],
            },
        ),
        (
            [NodeId(4), NodeId(5)],
            [NodeId(5), NodeId(0), NodeId(5)],
            BiPositionLabels {
                factor_1: vec![43, 47],
                factor_2: vec![53, 59, 61],
            },
        ),
    ]
}

fn assert_fixed_fixed_birelation_rows<const N1: usize, const N2: usize>(
    relations: &FixedFixedBirelationSet<NodeId, N1, NodeId, N2, BiPositionLabels>,
    entries: &[([NodeId; N1], [NodeId; N2], BiPositionLabels)],
) {
    assert_eq!(relations.count(), entries.len());
    assert_eq!(
        relations.ids().collect::<Vec<_>>(),
        (0..entries.len()).map(RelationId::from).collect::<Vec<_>>()
    );
    for (index, (first, second, data)) in entries.iter().enumerate() {
        let id = RelationId::from(index);
        assert_eq!(relations.participants_1(id), first);
        assert_eq!(relations.participants_2(id), second);
        assert_eq!(relations.data(id), data);
        let query_1: Vec<_> = first.iter().rev().copied().collect();
        let query_2: Vec<_> = second.iter().rev().copied().collect();
        assert!(relations.is_coincident(id, &query_1, &query_2));
        for node in [0, 1, 2, 3, 4, 5, 6, 7, u32::MAX].map(NodeId) {
            let expected = entries
                .iter()
                .position(|(a, b, _)| {
                    a.iter().chain(b).any(|p| *p == node)
                        && a.iter().all(|p| {
                            a.iter().filter(|q| *q == p).count()
                                == query_1.iter().filter(|q| *q == p).count()
                        })
                        && b.iter().all(|p| {
                            b.iter().filter(|q| *q == p).count()
                                == query_2.iter().filter(|q| *q == p).count()
                        })
                })
                .map(RelationId::from);
            assert_eq!(
                relations.coincident_to_node(node, &query_1, &query_2),
                expected
            );
        }
    }
    for node in [0, 1, 2, 3, 4, 5, 6, 7, u32::MAX].map(NodeId) {
        let expected: Vec<_> = entries
            .iter()
            .enumerate()
            .filter(|(_, (a, b, _))| a.iter().chain(b).any(|p| *p == node))
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
fn fixed_fixed_birelation_set_compaction_input(
) -> FixedFixedBirelationSet<NodeId, 1, NodeId, 2, &'static str> {
    FixedFixedBirelationSet::new(vec![
        ([NodeId(0)], [NodeId(2), NodeId(4)], "keep"),
        ([NodeId(1)], [NodeId(5), NodeId(6)], "drop"),
    ])
}

#[rstest]
fn test_fixed_fixed_birelation_set_new() {
    let rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 2, &str> =
        FixedFixedBirelationSet::new(vec![([NodeId(0)], [NodeId(2), NodeId(1)], "x")]);
    assert_eq!(rs.count(), 1);
    assert_eq!(rs.participants_1(RelationId(0)), &[NodeId(0)]);
    assert_eq!(rs.participants_2(RelationId(0)), &[NodeId(2), NodeId(1)]);
    assert_eq!(rs.data(RelationId(0)), &"x");
}

#[rstest]
#[case::empty(vec![], vec![], vec![])]
#[case::repeated(
    vec![
        ([NodeId(2), NodeId(0), NodeId(2)], [NodeId(2), NodeId(1), NodeId(2)], "first"),
        ([NodeId(0), NodeId(2), NodeId(2)], [NodeId(1), NodeId(2), NodeId(2)], "duplicate"),
        ([NodeId(3), NodeId(3), NodeId(3)], [NodeId(3), NodeId(3), NodeId(3)], "other"),
    ],
    vec![RelationId(0), RelationId(1)],
    vec![RelationId(2)],
)]
fn test_fixed_fixed_birelation_set_new_incidence(
    #[case] entries: Vec<([NodeId; 3], [NodeId; 3], &str)>,
    #[case] at_two: Vec<RelationId>,
    #[case] at_three: Vec<RelationId>,
) {
    let relations = FixedFixedBirelationSet::<NodeId, 3, NodeId, 3, &str>::new(entries.clone());
    assert_eq!(relations.count(), entries.len());
    for (index, (first, second, data)) in entries.iter().enumerate() {
        let id = RelationId(index as u32);
        assert_eq!(relations.participants_1(id), first);
        assert_eq!(relations.participants_2(id), second);
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
        ([NodeId(2)], [NodeId(4), NodeId(1)], "first"),
        ([NodeId(3)], [NodeId(5), NodeId(0)], "second"),
    ],
)]
fn test_fixed_fixed_birelation_set_into_entries(
    #[case] entries: Vec<([NodeId; 1], [NodeId; 2], &str)>,
) {
    let rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 2, &str> =
        FixedFixedBirelationSet::new(entries.clone());
    assert_eq!(rs.into_entries(), entries);
}

#[rstest]
#[case::first(RelationId(0), true)]
#[case::out_of_range(RelationId(1), false)]
fn test_fixed_fixed_birelation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
    let rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, &str> =
        FixedFixedBirelationSet::new(vec![([NodeId(0)], [NodeId(1)], "x")]);
    assert_eq!(rs.contains(id), expected);
}

#[rstest]
fn test_fixed_fixed_birelation_set_relation_ids() {
    assert_exact_size(
        FixedFixedBirelationSet::<NodeId, 1, NodeId, 1, &str>::default().ids(),
        vec![],
    );
    let rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, &str> =
        FixedFixedBirelationSet::new(vec![
            ([NodeId(0)], [NodeId(1)], "a"),
            ([NodeId(2)], [NodeId(3)], "b"),
        ]);
    assert_exact_size(rs.ids(), vec![RelationId(0), RelationId(1)]);
}

#[rstest]
fn test_fixed_fixed_birelation_set_iter() {
    let empty = FixedFixedBirelationSet::<NodeId, 1, NodeId, 1, i32>::default();
    assert_eq!(empty.iter().collect::<Vec<_>>(), vec![]);

    let rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, i32> =
        FixedFixedBirelationSet::new(vec![
            ([NodeId(0)], [NodeId(1)], 1),
            ([NodeId(2)], [NodeId(3)], 2),
        ]);
    assert_eq!(rs.iter().len(), 2);
    assert_eq!(
        rs.iter().collect::<Vec<_>>(),
        vec![
            (RelationId(0), &[NodeId(0)], &[NodeId(1)], &1),
            (RelationId(1), &[NodeId(2)], &[NodeId(3)], &2),
        ],
    );
}

#[rstest]
fn test_fixed_fixed_birelation_set_iter_mut() {
    let mut empty = FixedFixedBirelationSet::<NodeId, 1, NodeId, 1, i32>::default();
    assert_eq!(empty.iter_mut().len(), 0);

    let mut rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, i32> =
        FixedFixedBirelationSet::new(vec![
            ([NodeId(0)], [NodeId(1)], 1),
            ([NodeId(2)], [NodeId(3)], 2),
        ]);
    for (_, first, second, data) in rs.iter_mut() {
        assert_eq!(second[0].0, first[0].0 + 1);
        *data *= 10;
    }
    assert_eq!(rs.data(RelationId(0)), &10);
    assert_eq!(rs.data(RelationId(1)), &20);
}

#[rstest]
fn test_fixed_fixed_birelation_set_data_mut() {
    let mut rs: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, i32> =
        FixedFixedBirelationSet::new(vec![([NodeId(0)], [NodeId(1)], 1)]);
    *rs.data_mut(RelationId(0)) = 99;
    assert_eq!(rs.data(RelationId(0)), &99);
}

#[rstest]
fn test_fixed_fixed_birelation_set_incidence() {
    let rs: FixedFixedBirelationSet<NodeId, 1, EdgeId, 1, &str> =
        FixedFixedBirelationSet::new(vec![([NodeId(0)], [EdgeId(7)], "x")]);
    assert_eq!(rs.incident_to_node(NodeId(0)), &[RelationId(0)]);
    assert_eq!(rs.incident_to_edge(EdgeId(7)), &[RelationId(0)]);
    assert!(rs.has_incident_to_node(NodeId(0)));
    assert!(rs.has_incident_to_edge(EdgeId(7)));
    assert!(!rs.has_incident_to_node(NodeId(5)));
}

#[rstest]
#[case::exact(vec![NodeId(0), NodeId(1)], vec![NodeId(2)], Some(RelationId(0)))]
#[case::reordered_factor(vec![NodeId(1), NodeId(0)], vec![NodeId(2)], Some(RelationId(0)))]
#[case::second(vec![NodeId(3), NodeId(4)], vec![NodeId(5)], Some(RelationId(1)))]
#[case::absent(vec![NodeId(0), NodeId(1)], vec![NodeId(9)], None)]
fn test_fixed_fixed_birelation_set_coincident_to_node(
    #[case] query_1: Vec<NodeId>,
    #[case] query_2: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
) {
    let rs: FixedFixedBirelationSet<NodeId, 2, NodeId, 1, ()> = FixedFixedBirelationSet::new(vec![
        ([NodeId(0), NodeId(1)], [NodeId(2)], ()),
        ([NodeId(3), NodeId(4)], [NodeId(5)], ()),
    ]);
    assert_eq!(
        query_1
            .first()
            .and_then(|&anchor| rs.coincident_to_node(anchor, &query_1, &query_2)),
        expected,
    );
}

#[rstest]
#[case::reordered(
    NodeId(2),
    vec![NodeId(2), NodeId(2), NodeId(0)],
    vec![NodeId(2), NodeId(2), NodeId(1)],
    Some(RelationId(0)),
    true,
)]
#[case::multiplicity(
    NodeId(2),
    vec![NodeId(0), NodeId(0), NodeId(2)],
    vec![NodeId(2), NodeId(2), NodeId(1)],
    None,
    false,
)]
#[case::short(
    NodeId(2),
    vec![NodeId(2), NodeId(0)],
    vec![NodeId(2), NodeId(2), NodeId(1)],
    None,
    false,
)]
#[case::long(
    NodeId(2),
    vec![NodeId(2), NodeId(0), NodeId(2), NodeId(2)],
    vec![NodeId(2), NodeId(2), NodeId(1)],
    None,
    false,
)]
#[case::absent_anchor(
    NodeId(4),
    vec![NodeId(2), NodeId(2), NodeId(0)],
    vec![NodeId(2), NodeId(2), NodeId(1)],
    None,
    true,
)]
#[case::second_multiplicity(
    NodeId(2),
    vec![NodeId(2), NodeId(2), NodeId(0)],
    vec![NodeId(2), NodeId(1), NodeId(1)],
    None,
    false,
)]
#[case::factor_swap(
    NodeId(2),
    vec![NodeId(2), NodeId(2), NodeId(1)],
    vec![NodeId(2), NodeId(2), NodeId(0)],
    None,
    false,
)]
fn test_fixed_fixed_birelation_set_coincident_to_node_multiplicity(
    #[case] anchor: NodeId,
    #[case] query: Vec<NodeId>,
    #[case] query_2: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
    #[case] coincides: bool,
) {
    let relations = FixedFixedBirelationSet::<NodeId, 3, NodeId, 3, &str>::new(vec![
        (
            [NodeId(2), NodeId(0), NodeId(2)],
            [NodeId(2), NodeId(1), NodeId(2)],
            "first",
        ),
        (
            [NodeId(0), NodeId(2), NodeId(2)],
            [NodeId(1), NodeId(2), NodeId(2)],
            "duplicate",
        ),
        (
            [NodeId(3), NodeId(3), NodeId(3)],
            [NodeId(3), NodeId(3), NodeId(3)],
            "other",
        ),
    ]);
    assert_eq!(
        relations.coincident_to_node(anchor, &query, &query_2),
        expected
    );
    assert_eq!(
        relations.is_coincident(RelationId(0), &query, &query_2),
        coincides
    );
}

#[rstest]
#[case::empty(vec![], [NodeId(2), NodeId(2)], [NodeId(2), NodeId(3)])]
#[case::shared(vec![([NodeId(0), NodeId(2)], [NodeId(2), NodeId(1)], "old")], [NodeId(2), NodeId(2)], [NodeId(2), NodeId(3)])]
#[case::coinciding(vec![([NodeId(2), NodeId(2)], [NodeId(2), NodeId(3)], "old")], [NodeId(2), NodeId(2)], [NodeId(2), NodeId(3)])]
#[case::reordered(vec![([NodeId(0), NodeId(1)], [NodeId(2), NodeId(3)], "old")], [NodeId(1), NodeId(0)], [NodeId(3), NodeId(2)])]
#[case::sparse(vec![], [NodeId(u32::MAX)], [NodeId(u32::MAX), NodeId(0)])]
#[case::first_empty(vec![([], [NodeId(0), NodeId(1)], "old")], [], [NodeId(1), NodeId(1)])]
#[case::second_empty(vec![([NodeId(0), NodeId(1)], [], "old")], [NodeId(1), NodeId(1)], [])]
#[case::both_empty(vec![([], [], "first"), ([], [], "second")], [], [])]
fn test_fixed_fixed_birelation_set_add<const N1: usize, const N2: usize>(
    #[case] entries: Vec<([NodeId; N1], [NodeId; N2], &'static str)>,
    #[case] first: [NodeId; N1],
    #[case] second: [NodeId; N2],
) {
    let mut relations = FixedFixedBirelationSet::new(entries.clone());
    assert_eq!(
        relations.add(first, second, "added"),
        RelationId::from(entries.len())
    );
    let mut expected = entries;
    expected.push((first, second, "added"));
    assert_eq!(relations.count(), expected.len());
    for (i, (a, b, data)) in expected.iter().enumerate() {
        assert_eq!(relations.participants_1(RelationId::from(i)), a);
        assert_eq!(relations.participants_2(RelationId::from(i)), b);
        assert_eq!(relations.data(RelationId::from(i)), data);
    }
    for node in [0, 1, 2, 3, u32::MAX].map(NodeId) {
        let incidence: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (a, b, _))| a.contains(&node) || b.contains(&node))
            .map(|(i, _)| RelationId::from(i))
            .collect();
        assert_eq!(relations.incident_to_node(node), incidence);
        assert_eq!(relations.incident_to_edge(EdgeId(node.0)), &[]);
    }
    assert_eq!(relations.into_entries(), expected);
}

#[rstest]
fn test_fixed_fixed_birelation_set_add_payload() {
    let mut relations = FixedFixedBirelationSet::default();
    assert_eq!(
        relations.add(
            [NodeId(0)],
            [EdgeId(2), EdgeId(0)],
            NonCloneData(vec![7, 11])
        ),
        RelationId(0)
    );
    assert_eq!(
        relations.add(
            [NodeId(1)],
            [EdgeId(2), EdgeId(2)],
            NonCloneData(vec![13, 17])
        ),
        RelationId(1)
    );
    assert_eq!(relations.incident_to_node(NodeId(0)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_node(NodeId(1)), &[RelationId(1)]);
    assert_eq!(relations.incident_to_node(NodeId(2)), &[]);
    assert_eq!(relations.incident_to_edge(EdgeId(0)), &[RelationId(0)]);
    assert_eq!(
        relations.incident_to_edge(EdgeId(2)),
        &[RelationId(0), RelationId(1)]
    );
    assert_eq!(
        relations.into_entries(),
        vec![
            (
                [NodeId(0)],
                [EdgeId(2), EdgeId(0)],
                NonCloneData(vec![7, 11])
            ),
            (
                [NodeId(1)],
                [EdgeId(2), EdgeId(2)],
                NonCloneData(vec![13, 17])
            ),
        ]
    );
}

#[rstest]
#[case::first(vec![RelationId(0)], vec![1, 2])]
#[case::middle(vec![RelationId(1)], vec![0, 2])]
#[case::last(vec![RelationId(2)], vec![0, 1])]
#[case::unordered_repeated(vec![RelationId(2), RelationId(0), RelationId(2)], vec![1])]
#[case::all(vec![RelationId(2), RelationId(0), RelationId(1)], vec![])]
fn test_fixed_fixed_birelation_set_tracked_remove(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] ids: Vec<RelationId>,
    #[case] survivors: Vec<usize>,
) {
    let entries = fixed_fixed_birelation_set_mutation_entries;
    let expected: Vec<_> = survivors.iter().map(|&i| entries[i].clone()).collect();
    let mut relations = FixedFixedBirelationSet::new(entries);
    let mut plain = relations.clone();
    plain.remove(&ids);
    let compaction = relations.tracked_remove(&ids);
    assert_eq!(relations, plain);
    assert_eq!(compaction.source_count(), 3);
    assert_eq!(compaction.result_count(), survivors.len());
    let removed: Vec<_> = (0..3)
        .filter(|i| !survivors.contains(i))
        .map(RelationId::from)
        .collect();
    assert_eq!(compaction.removed(), removed);
    for old in 0..=3 {
        assert_eq!(
            compaction.compact(RelationId::from(old)),
            survivors
                .iter()
                .position(|&i| i == old)
                .map(RelationId::from)
        );
    }
    assert_fixed_fixed_birelation_rows(&relations, &expected);
    assert_eq!(relations.into_entries(), expected);
}

#[rstest]
#[case::empty(FixedFixedBirelationSet::<NodeId, 0, NodeId, 0, &str>::default())]
#[case::both_empty(FixedFixedBirelationSet::<NodeId, 0, NodeId, 0, &str>::new(vec![([], [], "first"), ([], [], "second")]))]
#[case::nonempty(FixedFixedBirelationSet::new(vec![([NodeId(0)], [NodeId(2), NodeId(0)], "first")]))]
fn test_fixed_fixed_birelation_set_tracked_remove_identity<const N1: usize, const N2: usize>(
    #[case] input: FixedFixedBirelationSet<NodeId, N1, NodeId, N2, &'static str>,
) {
    let mut relations = input.clone();
    assert_eq!(
        relations.tracked_remove(&[]),
        Compaction::identity(input.count())
    );
    assert_eq!(relations, input);
    for node in [NodeId(0), NodeId(2)] {
        assert_eq!(
            relations.incident_to_node(node),
            input.incident_to_node(node)
        );
    }
}

#[rstest]
#[case::first_empty([], [NodeId(0), NodeId(1)])]
#[case::second_empty([NodeId(0), NodeId(1)], [])]
#[case::both_empty([], [])]
fn test_fixed_fixed_birelation_set_tracked_remove_arity<const N1: usize, const N2: usize>(
    #[case] first: [NodeId; N1],
    #[case] second: [NodeId; N2],
) {
    let mut relations = FixedFixedBirelationSet::new(vec![
        (first, second, "first"),
        (first, second, "second"),
        (first, second, "third"),
    ]);
    let compaction = relations.tracked_remove(&[RelationId(1)]);
    assert_eq!(compaction.source_count(), 3);
    assert_eq!(compaction.result_count(), 2);
    assert_eq!(compaction.removed(), &[RelationId(1)]);
    for node in [NodeId(0), NodeId(1), NodeId(2)] {
        let expected = if first.contains(&node) || second.contains(&node) {
            vec![RelationId(0), RelationId(1)]
        } else {
            vec![]
        };
        assert_eq!(relations.incident_to_node(node), expected);
    }
    assert_eq!(
        relations.into_entries(),
        vec![(first, second, "first"), (first, second, "third")]
    );
}

#[rstest]
fn test_fixed_fixed_birelation_set_tracked_remove_payload() {
    let mut relations = FixedFixedBirelationSet::default();
    relations.add(
        [NodeId(0)],
        [EdgeId(2), EdgeId(0)],
        NonCloneData(vec![7, 11]),
    );
    relations.add(
        [NodeId(1)],
        [EdgeId(2), EdgeId(2)],
        NonCloneData(vec![13, 17]),
    );
    relations.add(
        [NodeId(2)],
        [EdgeId(3), EdgeId(2)],
        NonCloneData(vec![19, 23]),
    );
    let compaction = relations.tracked_remove(&[RelationId(0)]);
    assert_eq!(compaction.source_count(), 3);
    assert_eq!(compaction.result_count(), 2);
    assert_eq!(compaction.removed(), &[RelationId(0)]);
    assert_eq!(relations.incident_to_node(NodeId(0)), &[]);
    assert_eq!(relations.incident_to_node(NodeId(1)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_node(NodeId(2)), &[RelationId(1)]);
    assert_eq!(relations.incident_to_edge(EdgeId(0)), &[]);
    assert_eq!(
        relations.incident_to_edge(EdgeId(2)),
        &[RelationId(0), RelationId(1)]
    );
    assert_eq!(relations.incident_to_edge(EdgeId(3)), &[RelationId(1)]);
    assert_eq!(
        relations.into_entries(),
        vec![
            (
                [NodeId(1)],
                [EdgeId(2), EdgeId(2)],
                NonCloneData(vec![13, 17])
            ),
            (
                [NodeId(2)],
                [EdgeId(3), EdgeId(2)],
                NonCloneData(vec![19, 23])
            ),
        ]
    );
}

#[rstest]
#[case::empty(vec![], vec![RelationId(0)])]
#[case::end(vec![([NodeId(0)], [NodeId(2)], "first")], vec![RelationId(1)])]
#[case::mixed(vec![([NodeId(0)], [NodeId(2)], "first"), ([NodeId(2)], [NodeId(0)], "second")], vec![RelationId(0), RelationId(2)])]
#[case::distant(vec![([NodeId(0)], [NodeId(2)], "first")], vec![RelationId(u32::MAX), RelationId(0)])]
fn test_fixed_fixed_birelation_set_tracked_remove_error(
    #[case] entries: Vec<([NodeId; 1], [NodeId; 1], &'static str)>,
    #[case] ids: Vec<RelationId>,
) {
    let original = FixedFixedBirelationSet::new(entries);
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
    vec![([NodeId(1), NodeId(1)], [NodeId(3)], "changed")], 3,
    vec![RelationId(0), RelationId(2)],
    vec![(RelationId(2), [NodeId(4), NodeId(4)], [NodeId(1)], "last"), (RelationId(0), [NodeId(0), NodeId(1)], [NodeId(0)], "first")],
    vec![([NodeId(0), NodeId(1)], [NodeId(0)], "first"), ([NodeId(1), NodeId(1)], [NodeId(3)], "changed"), ([NodeId(4), NodeId(4)], [NodeId(1)], "last")],
)]
#[case::all(vec![], 2, vec![RelationId(0), RelationId(1)],
    vec![(RelationId(1), [NodeId(1)], [NodeId(2), NodeId(1)], "second"), (RelationId(0), [NodeId(0)], [NodeId(1), NodeId(0)], "first")],
    vec![([NodeId(0)], [NodeId(1), NodeId(0)], "first"), ([NodeId(1)], [NodeId(2), NodeId(1)], "second")])]
#[case::empty_first(vec![([], [NodeId(1)], "survivor")], 2, vec![RelationId(0)],
    vec![(RelationId(0), [], [NodeId(0)], "restored")], vec![([], [NodeId(0)], "restored"), ([], [NodeId(1)], "survivor")])]
#[case::empty_second(vec![([NodeId(1)], [], "survivor")], 2, vec![RelationId(1)],
    vec![(RelationId(1), [NodeId(0)], [], "restored")], vec![([NodeId(1)], [], "survivor"), ([NodeId(0)], [], "restored")])]
#[case::empty_both(vec![], 1, vec![RelationId(0)], vec![(RelationId(0), [], [], "restored")], vec![([], [], "restored")])]
#[case::coincident(vec![([NodeId(1)], [NodeId(1)], "survivor")], 2, vec![RelationId(0)],
    vec![(RelationId(0), [NodeId(1)], [NodeId(1)], "restored")], vec![([NodeId(1)], [NodeId(1)], "restored"), ([NodeId(1)], [NodeId(1)], "survivor")])]
fn test_fixed_fixed_birelation_set_restore<const N1: usize, const N2: usize>(
    #[case] entries: Vec<([NodeId; N1], [NodeId; N2], &'static str)>,
    #[case] count: usize,
    #[case] ids: Vec<RelationId>,
    #[case] removed: Vec<(RelationId, [NodeId; N1], [NodeId; N2], &'static str)>,
    #[case] expected: Vec<([NodeId; N1], [NodeId; N2], &'static str)>,
) {
    let mut relations = FixedFixedBirelationSet::new(entries);
    relations.restore(&Compaction::new(count, ids).unwrap(), removed);
    for node in 0..6 {
        let incidence: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (first, second, _))| {
                first.contains(&NodeId(node)) || second.contains(&NodeId(node))
            })
            .map(|(id, _)| RelationId::from(id))
            .collect();
        assert_eq!(relations.incident_to_node(NodeId(node)), incidence);
        assert_eq!(relations.incident_to_edge(EdgeId(node)), &[]);
    }
    assert_eq!(relations.into_entries(), expected);
}

#[rstest]
#[case::empty(FixedFixedBirelationSet::<NodeId, 1, EdgeId, 2, &str>::new(vec![]))]
#[case::empty_first(FixedFixedBirelationSet::new(vec![([], [EdgeId(0)], "first")]))]
#[case::empty_second(FixedFixedBirelationSet::new(vec![([NodeId(0)], [], "first")]))]
#[case::empty_both(FixedFixedBirelationSet::new(vec![([], [], "first")]))]
#[case::nonempty(FixedFixedBirelationSet::new(vec![([NodeId(0)], [EdgeId(2), EdgeId(2)], "first")]))]
fn test_fixed_fixed_birelation_set_restore_identity<const N1: usize, const N2: usize>(
    #[case] input: FixedFixedBirelationSet<NodeId, N1, EdgeId, N2, &'static str>,
) {
    let mut relations = input.clone();
    relations.restore(&Compaction::identity(input.count()), vec![]);
    assert_eq!(relations, input);
}

#[rstest]
#[case::mixed([NodeId(2), NodeId(2)], [EdgeId(0)], [NodeId(0), NodeId(2)], [EdgeId(2)])]
#[case::empty_first([], [EdgeId(0)], [], [EdgeId(2)])]
#[case::empty_second([NodeId(2)], [], [NodeId(0)], [])]
#[case::empty_both([], [], [], [])]
fn test_fixed_fixed_birelation_set_restore_payload<const N1: usize, const N2: usize>(
    #[case] first: [NodeId; N1],
    #[case] second: [EdgeId; N2],
    #[case] saved_first: [NodeId; N1],
    #[case] saved_second: [EdgeId; N2],
) {
    let mut relations = FixedFixedBirelationSet::default();
    relations.add(first, second, NonCloneData(vec![7, 11]));
    relations.restore(
        &Compaction::new(2, vec![RelationId(0)]).unwrap(),
        vec![(
            RelationId(0),
            saved_first,
            saved_second,
            NonCloneData(vec![13, 17]),
        )],
    );
    let expected = vec![
        (saved_first, saved_second, NonCloneData(vec![13, 17])),
        (first, second, NonCloneData(vec![7, 11])),
    ];
    for id in 0..4 {
        let nodes: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (row, _, _))| row.contains(&NodeId(id)))
            .map(|(i, _)| RelationId::from(i))
            .collect();
        let edges: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (_, row, _))| row.contains(&EdgeId(id)))
            .map(|(i, _)| RelationId::from(i))
            .collect();
        assert_eq!(relations.incident_to_node(NodeId(id)), nodes);
        assert_eq!(relations.incident_to_edge(EdgeId(id)), edges);
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
fn test_fixed_fixed_birelation_set_restore_malformed(
    #[case] count: usize,
    #[case] ids: Vec<RelationId>,
    #[case] removed: Vec<RelationId>,
) {
    let mut relations = FixedFixedBirelationSet::default();
    relations.add([NodeId(0)], [EdgeId(0)], NonCloneData(vec![7]));
    relations.restore(
        &Compaction::new(count, ids).unwrap(),
        removed
            .into_iter()
            .map(|id| {
                (
                    id,
                    [NodeId(u32::MAX)],
                    [EdgeId(u32::MAX)],
                    NonCloneData(vec![11]),
                )
            })
            .collect(),
    );
}

#[rstest]
#[case::mixed(vec![([NodeId(0), NodeId(1)], [EdgeId(1)], "first"), ([NodeId(2), NodeId(0)], [EdgeId(0)], "second")],
    vec![([NodeId(1), NodeId(3)], [EdgeId(2)], "first"), ([NodeId(4), NodeId(1)], [EdgeId(0)], "second")])]
#[case::shared_nodes(vec![([NodeId(0)], [NodeId(0), NodeId(1)], "first")], vec![([NodeId(1)], [NodeId(1), NodeId(3)], "first")])]
#[case::shared_edges(vec![([EdgeId(1)], [EdgeId(1), EdgeId(0)], "first")], vec![([EdgeId(2)], [EdgeId(2), EdgeId(0)], "first")])]
#[case::empty_first(vec![([] as [NodeId; 0], [EdgeId(1)], "first")], vec![([], [EdgeId(2)], "first")])]
#[case::empty_second(vec![([NodeId(0)], [] as [EdgeId; 0], "first")], vec![([NodeId(1)], [], "first")])]
#[case::empty_both(vec![([] as [NodeId; 0], [] as [EdgeId; 0], "first")], vec![([], [], "first")])]
#[case::empty(Vec::<([NodeId; 1], [EdgeId; 2], &str)>::new(), vec![])]
fn test_fixed_fixed_birelation_set_restore_participants<
    L1: RelationParticipant + Debug,
    const N1: usize,
    L2: RelationParticipant + Debug,
    const N2: usize,
>(
    #[case] entries: Vec<([L1; N1], [L2; N2], &'static str)>,
    #[case] expected: Vec<([L1; N1], [L2; N2], &'static str)>,
) {
    let mut relations = FixedFixedBirelationSet::new(entries);
    relations.restore_participants(&GraphCompaction::new(
        Compaction::new(5, vec![NodeId(0), NodeId(2)]).unwrap(),
        Compaction::new(4, vec![EdgeId(1)]).unwrap(),
    ));
    for id in 0..6 {
        let nodes: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (first, second, _))| {
                first
                    .iter()
                    .map(|p| p.refs())
                    .chain(second.iter().map(|p| p.refs()))
                    .any(|p| p.node == Some(NodeId(id)))
            })
            .map(|(i, _)| RelationId::from(i))
            .collect();
        let edges: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (first, second, _))| {
                first
                    .iter()
                    .map(|p| p.refs())
                    .chain(second.iter().map(|p| p.refs()))
                    .any(|p| p.edge == Some(EdgeId(id)))
            })
            .map(|(i, _)| RelationId::from(i))
            .collect();
        assert_eq!(relations.incident_to_node(NodeId(id)), nodes);
        assert_eq!(relations.incident_to_edge(EdgeId(id)), edges);
    }
    assert_eq!(relations.into_entries(), expected);
}

#[rstest]
fn test_fixed_fixed_birelation_set_restore_participants_references() {
    let mut relations = FixedFixedBirelationSet::default();
    relations.add(
        [
            References {
                node: Some(NodeId(0)),
                edge: Some(EdgeId(1)),
                label: 1,
            },
            References {
                node: None,
                edge: None,
                label: 2,
            },
        ],
        [
            References {
                node: Some(NodeId(0)),
                edge: None,
                label: 3,
            },
            References {
                node: None,
                edge: Some(EdgeId(1)),
                label: 4,
            },
            References {
                node: Some(NodeId(0)),
                edge: Some(EdgeId(1)),
                label: 5,
            },
        ],
        NonCloneData(vec![7, 11]),
    );
    relations.restore_participants(&GraphCompaction::new(
        Compaction::new(3, vec![NodeId(0)]).unwrap(),
        Compaction::new(4, vec![EdgeId(1)]).unwrap(),
    ));
    assert_eq!(relations.incident_to_node(NodeId(1)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_edge(EdgeId(2)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_node(NodeId(0)), &[]);
    assert_eq!(relations.incident_to_edge(EdgeId(1)), &[]);
    assert_eq!(
        relations.into_entries(),
        vec![(
            [
                References {
                    node: Some(NodeId(1)),
                    edge: Some(EdgeId(2)),
                    label: 1
                },
                References {
                    node: None,
                    edge: None,
                    label: 2
                }
            ],
            [
                References {
                    node: Some(NodeId(1)),
                    edge: None,
                    label: 3
                },
                References {
                    node: None,
                    edge: Some(EdgeId(2)),
                    label: 4
                },
                References {
                    node: Some(NodeId(1)),
                    edge: Some(EdgeId(2)),
                    label: 5
                }
            ],
            NonCloneData(vec![7, 11])
        )]
    );
}

#[rstest]
#[case::empty(FixedFixedBirelationSet::<NodeId, 1, EdgeId, 2, &str>::new(vec![]))]
#[case::empty_first(FixedFixedBirelationSet::new(vec![([], [EdgeId(0)], "first")]))]
#[case::empty_second(FixedFixedBirelationSet::new(vec![([NodeId(0)], [], "first")]))]
#[case::empty_both(FixedFixedBirelationSet::new(vec![([], [], "first")]))]
#[case::nonempty(FixedFixedBirelationSet::new(vec![([NodeId(0)], [EdgeId(2), EdgeId(2)], "first")]))]
fn test_fixed_fixed_birelation_set_restore_participants_identity<
    const N1: usize,
    const N2: usize,
>(
    #[case] input: FixedFixedBirelationSet<NodeId, N1, EdgeId, N2, &'static str>,
) {
    let mut relations = input.clone();
    relations.restore_participants(&GraphCompaction::new(
        Compaction::identity(3),
        Compaction::identity(3),
    ));
    assert_eq!(relations, input);
}

#[rstest]
#[case::first_node(3, 3, Some(NodeId(2)), None, false)]
#[case::second_node(3, 3, Some(NodeId(2)), None, true)]
#[case::first_edge(3, 3, None, Some(EdgeId(2)), false)]
#[case::second_edge(3, 3, None, Some(EdgeId(2)), true)]
#[case::both_domains(3, 3, Some(NodeId(u32::MAX)), Some(EdgeId(u32::MAX)), true)]
#[case::oversized_nodes(usize::MAX, 3, Some(NodeId(0)), None, false)]
#[case::oversized_edges(3, usize::MAX, None, Some(EdgeId(0)), true)]
fn test_fixed_fixed_birelation_set_restore_participants_malformed(
    #[case] nodes: usize,
    #[case] edges: usize,
    #[case] node: Option<NodeId>,
    #[case] edge: Option<EdgeId>,
    #[case] second: bool,
) {
    let valid = References {
        node: Some(NodeId(0)),
        edge: Some(EdgeId(0)),
        label: 1,
    };
    let invalid = References {
        node,
        edge,
        label: 2,
    };
    let mut relations = FixedFixedBirelationSet::default();
    relations.add(
        [valid, if second { valid } else { invalid }],
        [if second { invalid } else { valid }],
        NonCloneData(vec![7]),
    );
    relations.restore_participants(&GraphCompaction::new(
        Compaction::new(nodes, vec![NodeId(0)]).unwrap(),
        Compaction::new(edges, vec![EdgeId(0)]).unwrap(),
    ));
}

#[rstest]
fn test_fixed_fixed_birelation_set_permute_participants_1() {
    let mut rs: FixedFixedBirelationSet<NodeId, 3, EdgeId, 2, &str> =
        FixedFixedBirelationSet::new(vec![(
            [NodeId(0), NodeId(1), NodeId(2)],
            [EdgeId(7), EdgeId(8)],
            "a",
        )]);
    rs.permute_participants_1(
        RelationId(0),
        &[
            ParticipantPosition(2),
            ParticipantPosition(0),
            ParticipantPosition(1),
        ],
    );
    assert_eq!(
        rs.participants_1(RelationId(0)),
        &[NodeId(2), NodeId(0), NodeId(1)]
    );
    assert_eq!(rs.participants_2(RelationId(0)), &[EdgeId(7), EdgeId(8)]);
    assert_eq!(rs.data(RelationId(0)), &"a");
}

#[rstest]
fn test_fixed_fixed_birelation_set_permute_participants_2() {
    let mut rs: FixedFixedBirelationSet<NodeId, 3, EdgeId, 2, &str> =
        FixedFixedBirelationSet::new(vec![(
            [NodeId(0), NodeId(1), NodeId(2)],
            [EdgeId(7), EdgeId(8)],
            "a",
        )]);
    rs.permute_participants_2(
        RelationId(0),
        &[ParticipantPosition(1), ParticipantPosition(0)],
    );
    assert_eq!(
        rs.participants_1(RelationId(0)),
        &[NodeId(0), NodeId(1), NodeId(2)]
    );
    assert_eq!(rs.participants_2(RelationId(0)), &[EdgeId(8), EdgeId(7)]);
    assert_eq!(rs.data(RelationId(0)), &"a");
}

#[rstest]
#[case::reordered(RelationId(0), [NodeId(1), NodeId(0)], [NodeId(1), NodeId(1), NodeId(2)])]
#[case::repeated(RelationId(1), [NodeId(6), NodeId(6)], [NodeId(6), NodeId(6), NodeId(6)])]
#[case::shared_second(RelationId(0), [NodeId(6), NodeId(7)], [NodeId(1), NodeId(2), NodeId(1)])]
#[case::shared_first(RelationId(0), [NodeId(0), NodeId(1)], [NodeId(6), NodeId(7), NodeId(6)])]
#[case::removed_both(RelationId(0), [NodeId(6), NodeId(7)], [NodeId(2), NodeId(7), NodeId(7)])]
#[case::sparse_last(RelationId(2), [NodeId(u32::MAX), NodeId(5)], [NodeId(u32::MAX), NodeId(0), NodeId(6)])]
fn test_fixed_fixed_birelation_set_replace_participants(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] first: [NodeId; 2],
    #[case] second: [NodeId; 3],
) {
    let mut entries = fixed_fixed_birelation_set_mutation_entries;
    let mut relations = FixedFixedBirelationSet::new(entries.clone());
    relations.replace_participants(id, first, second);
    entries[id.index()].0 = first;
    entries[id.index()].1 = second;
    assert_fixed_fixed_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::both_empty([], [])]
#[case::first_empty([], [NodeId(1), NodeId(2), NodeId(1)])]
#[case::second_empty([NodeId(0), NodeId(1)], [])]
#[case::both_present([NodeId(1), NodeId(0)], [NodeId(2), NodeId(1), NodeId(1)])]
fn test_fixed_fixed_birelation_set_replace_participants_identity<
    const N1: usize,
    const N2: usize,
>(
    #[case] first: [NodeId; N1],
    #[case] second: [NodeId; N2],
) {
    let entries = vec![(
        first,
        second,
        BiPositionLabels {
            factor_1: vec![7, 11],
            factor_2: vec![13, 17, 19],
        },
    )];
    let input = FixedFixedBirelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participants(RelationId(0), first, second);
    assert_eq!(relations, input);
    assert_fixed_fixed_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_fixed_birelation_set_replace_participants_edges() {
    let mut relations = FixedFixedBirelationSet::new(vec![(
        [NodeId(0), NodeId(1)],
        [EdgeId(1), EdgeId(2), EdgeId(1)],
        7,
    )]);
    relations.replace_participants(
        RelationId(0),
        [NodeId(5), NodeId(5)],
        [EdgeId(6), EdgeId(6), EdgeId(6)],
    );
    let first = [NodeId(5), NodeId(5)];
    let second = [EdgeId(6), EdgeId(6), EdgeId(6)];
    for key in 0..7 {
        let nodes = if first.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if second.contains(&EdgeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        assert_eq!(relations.incident_to_node(NodeId(key)), nodes);
        assert_eq!(relations.incident_to_edge(EdgeId(key)), edges);
        assert_eq!(
            relations.coincident_to_node(NodeId(key), &first, &second),
            nodes.first().copied()
        );
        assert_eq!(
            relations.coincident_to_edge(EdgeId(key), &first, &second),
            edges.first().copied()
        );
    }
    assert_eq!(relations.into_entries(), vec![(first, second, 7)]);
}

#[rstest]
#[case::past_end(RelationId(3))]
#[case::maximum(RelationId(u32::MAX))]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participants_error(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] id: RelationId,
) {
    FixedFixedBirelationSet::new(fixed_fixed_birelation_set_mutation_entries).replace_participants(
        id,
        [NodeId(6), NodeId(7)],
        [NodeId(6), NodeId(7), NodeId(6)],
    );
}

#[rstest]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participants_empty() {
    FixedFixedBirelationSet::<NodeId, 2, NodeId, 3, NonCloneData>::default().replace_participants(
        RelationId(0),
        [NodeId(6), NodeId(7)],
        [NodeId(6), NodeId(7), NodeId(6)],
    );
}

#[rstest]
#[case::shared_other_factor(RelationId(0), [NodeId(6), NodeId(7)])]
#[case::reordered(RelationId(0), [NodeId(1), NodeId(0)])]
#[case::repeated(RelationId(1), [NodeId(6), NodeId(6)])]
#[case::sparse_last(RelationId(2), [NodeId(u32::MAX), NodeId(5)])]
fn test_fixed_fixed_birelation_set_replace_participants_1(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] participants: [NodeId; 2],
) {
    let mut entries = fixed_fixed_birelation_set_mutation_entries;
    let mut relations = FixedFixedBirelationSet::new(entries.clone());
    relations.replace_participants_1(id, participants);
    entries[id.index()].0 = participants;
    assert_fixed_fixed_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::both_empty([], [])]
#[case::first_empty([], [NodeId(1), NodeId(2), NodeId(1)])]
#[case::second_empty([NodeId(0), NodeId(1)], [])]
#[case::both_present([NodeId(1), NodeId(0)], [NodeId(2), NodeId(1), NodeId(1)])]
fn test_fixed_fixed_birelation_set_replace_participants_1_identity<
    const N1: usize,
    const N2: usize,
>(
    #[case] first: [NodeId; N1],
    #[case] second: [NodeId; N2],
) {
    let entries = vec![(
        first,
        second,
        BiPositionLabels {
            factor_1: vec![7, 11],
            factor_2: vec![13, 17, 19],
        },
    )];
    let input = FixedFixedBirelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participants_1(RelationId(0), first);
    assert_eq!(relations, input);
    assert_fixed_fixed_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_fixed_birelation_set_replace_participants_1_edges() {
    let mut relations = FixedFixedBirelationSet::new(vec![(
        [NodeId(0), NodeId(1)],
        [EdgeId(1), EdgeId(2), EdgeId(1)],
        7,
    )]);
    relations.replace_participants_1(RelationId(0), [NodeId(5), NodeId(5)]);
    let first = [NodeId(5), NodeId(5)];
    let second = [EdgeId(1), EdgeId(2), EdgeId(1)];
    for key in 0..7 {
        let nodes = if first.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if second.contains(&EdgeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        assert_eq!(relations.incident_to_node(NodeId(key)), nodes);
        assert_eq!(relations.incident_to_edge(EdgeId(key)), edges);
        assert_eq!(
            relations.coincident_to_node(NodeId(key), &first, &second),
            nodes.first().copied()
        );
        assert_eq!(
            relations.coincident_to_edge(EdgeId(key), &first, &second),
            edges.first().copied()
        );
    }
    assert_eq!(relations.into_entries(), vec![(first, second, 7)]);
}

#[rstest]
#[case::past_end(RelationId(3))]
#[case::maximum(RelationId(u32::MAX))]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participants_1_error(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] id: RelationId,
) {
    FixedFixedBirelationSet::new(fixed_fixed_birelation_set_mutation_entries)
        .replace_participants_1(id, [NodeId(6), NodeId(7)]);
}

#[rstest]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participants_1_empty() {
    FixedFixedBirelationSet::<NodeId, 2, NodeId, 3, NonCloneData>::default()
        .replace_participants_1(RelationId(0), [NodeId(6), NodeId(7)]);
}

#[rstest]
#[case::shared_other_factor(RelationId(0), [NodeId(6), NodeId(7), NodeId(6)])]
#[case::reordered(RelationId(0), [NodeId(1), NodeId(1), NodeId(2)])]
#[case::repeated(RelationId(1), [NodeId(6), NodeId(6), NodeId(6)])]
#[case::sparse_last(RelationId(2), [NodeId(u32::MAX), NodeId(5), NodeId(0)])]
fn test_fixed_fixed_birelation_set_replace_participants_2(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] participants: [NodeId; 3],
) {
    let mut entries = fixed_fixed_birelation_set_mutation_entries;
    let mut relations = FixedFixedBirelationSet::new(entries.clone());
    relations.replace_participants_2(id, participants);
    entries[id.index()].1 = participants;
    assert_fixed_fixed_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::both_empty([], [])]
#[case::first_empty([], [NodeId(1), NodeId(2), NodeId(1)])]
#[case::second_empty([NodeId(0), NodeId(1)], [])]
#[case::both_present([NodeId(1), NodeId(0)], [NodeId(2), NodeId(1), NodeId(1)])]
fn test_fixed_fixed_birelation_set_replace_participants_2_identity<
    const N1: usize,
    const N2: usize,
>(
    #[case] first: [NodeId; N1],
    #[case] second: [NodeId; N2],
) {
    let entries = vec![(
        first,
        second,
        BiPositionLabels {
            factor_1: vec![7, 11],
            factor_2: vec![13, 17, 19],
        },
    )];
    let input = FixedFixedBirelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participants_2(RelationId(0), second);
    assert_eq!(relations, input);
    assert_fixed_fixed_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_fixed_birelation_set_replace_participants_2_edges() {
    let mut relations = FixedFixedBirelationSet::new(vec![(
        [NodeId(0), NodeId(1)],
        [EdgeId(1), EdgeId(2), EdgeId(1)],
        7,
    )]);
    relations.replace_participants_2(RelationId(0), [EdgeId(6), EdgeId(6), EdgeId(6)]);
    let first = [NodeId(0), NodeId(1)];
    let second = [EdgeId(6), EdgeId(6), EdgeId(6)];
    for key in 0..7 {
        let nodes = if first.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if second.contains(&EdgeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        assert_eq!(relations.incident_to_node(NodeId(key)), nodes);
        assert_eq!(relations.incident_to_edge(EdgeId(key)), edges);
        assert_eq!(
            relations.coincident_to_node(NodeId(key), &first, &second),
            nodes.first().copied()
        );
        assert_eq!(
            relations.coincident_to_edge(EdgeId(key), &first, &second),
            edges.first().copied()
        );
    }
    assert_eq!(relations.into_entries(), vec![(first, second, 7)]);
}

#[rstest]
#[case::past_end(RelationId(3))]
#[case::maximum(RelationId(u32::MAX))]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participants_2_error(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] id: RelationId,
) {
    FixedFixedBirelationSet::new(fixed_fixed_birelation_set_mutation_entries)
        .replace_participants_2(id, [NodeId(6), NodeId(7), NodeId(6)]);
}

#[rstest]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participants_2_empty() {
    FixedFixedBirelationSet::<NodeId, 2, NodeId, 3, NonCloneData>::default()
        .replace_participants_2(RelationId(0), [NodeId(6), NodeId(7), NodeId(6)]);
}

#[rstest]
#[case::first(RelationId(0), ParticipantPosition(0), NodeId(6), [NodeId(6), NodeId(1)])]
#[case::shared_other_factor(RelationId(0), ParticipantPosition(1), NodeId(7), [NodeId(0), NodeId(7)])]
#[case::middle_row(RelationId(1), ParticipantPosition(0), NodeId(4), [NodeId(4), NodeId(3)])]
#[case::sparse_last(RelationId(2), ParticipantPosition(1), NodeId(u32::MAX), [NodeId(4), NodeId(u32::MAX)])]
fn test_fixed_fixed_birelation_set_replace_participant_1(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
    #[case] participant: NodeId,
    #[case] expected: [NodeId; 2],
) {
    let mut entries = fixed_fixed_birelation_set_mutation_entries;
    let mut relations = FixedFixedBirelationSet::new(entries.clone());
    relations.replace_participant_1(id, position, participant);
    entries[id.index()].0 = expected;
    assert_fixed_fixed_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::position_0(ParticipantPosition(0))]
#[case::position_1(ParticipantPosition(1))]
fn test_fixed_fixed_birelation_set_replace_participant_1_identity(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] position: ParticipantPosition,
) {
    let entries = fixed_fixed_birelation_set_mutation_entries;
    let input = FixedFixedBirelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participant_1(RelationId(0), position, entries[0].0[position.index()]);
    assert_eq!(relations, input);
    assert_fixed_fixed_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_fixed_birelation_set_replace_participant_1_edges() {
    let mut relations = FixedFixedBirelationSet::new(vec![(
        [NodeId(0), NodeId(1)],
        [EdgeId(1), EdgeId(2), EdgeId(1)],
        7,
    )]);
    relations.replace_participant_1(RelationId(0), ParticipantPosition(1), NodeId(5));
    let first = [NodeId(0), NodeId(5)];
    let second = [EdgeId(1), EdgeId(2), EdgeId(1)];
    for key in 0..7 {
        let nodes = if first.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if second.contains(&EdgeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        assert_eq!(relations.incident_to_node(NodeId(key)), nodes);
        assert_eq!(relations.incident_to_edge(EdgeId(key)), edges);
        assert_eq!(
            relations.coincident_to_node(NodeId(key), &first, &second),
            nodes.first().copied()
        );
        assert_eq!(
            relations.coincident_to_edge(EdgeId(key), &first, &second),
            edges.first().copied()
        );
    }
    assert_eq!(relations.into_entries(), vec![(first, second, 7)]);
}

#[rstest]
fn test_fixed_fixed_birelation_set_replace_participant_1_other_zero_arity() {
    let mut relations = FixedFixedBirelationSet::<NodeId, 2, NodeId, 0, _>::new(vec![(
        [NodeId(0), NodeId(1)],
        [],
        7,
    )]);
    relations.replace_participant_1(RelationId(0), ParticipantPosition(1), NodeId(5));
    assert_eq!(relations.incident_to_node(NodeId(0)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_node(NodeId(1)), &[]);
    assert_eq!(relations.incident_to_node(NodeId(5)), &[RelationId(0)]);
    assert_eq!(
        relations.into_entries(),
        vec![([NodeId(0), NodeId(5)], [], 7)]
    );
}

#[rstest]
#[case::past_end(RelationId(3), ParticipantPosition(0))]
#[case::maximum_id(RelationId(u32::MAX), ParticipantPosition(0))]
#[case::at_length(RelationId(0), ParticipantPosition(2))]
#[case::maximum_position(RelationId(0), ParticipantPosition(u32::MAX))]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participant_1_error(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
) {
    FixedFixedBirelationSet::new(fixed_fixed_birelation_set_mutation_entries)
        .replace_participant_1(id, position, NodeId(6));
}

#[rstest]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participant_1_zero_arity() {
    let mut relations =
        FixedFixedBirelationSet::<NodeId, 0, NodeId, 1, ()>::new(vec![([], [NodeId(1)], ())]);
    relations.replace_participant_1(RelationId(0), ParticipantPosition(0), NodeId(6));
}

#[rstest]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participant_1_empty() {
    FixedFixedBirelationSet::<NodeId, 2, NodeId, 3, NonCloneData>::default().replace_participant_1(
        RelationId(0),
        ParticipantPosition(0),
        NodeId(6),
    );
}

#[rstest]
#[case::first(RelationId(0), ParticipantPosition(0), NodeId(6), [NodeId(6), NodeId(2), NodeId(1)])]
#[case::shared_other_factor(RelationId(0), ParticipantPosition(2), NodeId(7), [NodeId(1), NodeId(2), NodeId(7)])]
#[case::middle_row(RelationId(1), ParticipantPosition(0), NodeId(4), [NodeId(4), NodeId(4), NodeId(3)])]
#[case::sparse_last(RelationId(2), ParticipantPosition(2), NodeId(u32::MAX), [NodeId(5), NodeId(0), NodeId(u32::MAX)])]
fn test_fixed_fixed_birelation_set_replace_participant_2(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
    #[case] participant: NodeId,
    #[case] expected: [NodeId; 3],
) {
    let mut entries = fixed_fixed_birelation_set_mutation_entries;
    let mut relations = FixedFixedBirelationSet::new(entries.clone());
    relations.replace_participant_2(id, position, participant);
    entries[id.index()].1 = expected;
    assert_fixed_fixed_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::position_0(ParticipantPosition(0))]
#[case::position_1(ParticipantPosition(1))]
#[case::position_2(ParticipantPosition(2))]
fn test_fixed_fixed_birelation_set_replace_participant_2_identity(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] position: ParticipantPosition,
) {
    let entries = fixed_fixed_birelation_set_mutation_entries;
    let input = FixedFixedBirelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participant_2(RelationId(0), position, entries[0].1[position.index()]);
    assert_eq!(relations, input);
    assert_fixed_fixed_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_fixed_birelation_set_replace_participant_2_edges() {
    let mut relations = FixedFixedBirelationSet::new(vec![(
        [NodeId(0), NodeId(1)],
        [EdgeId(1), EdgeId(2), EdgeId(1)],
        7,
    )]);
    relations.replace_participant_2(RelationId(0), ParticipantPosition(2), EdgeId(6));
    let first = [NodeId(0), NodeId(1)];
    let second = [EdgeId(1), EdgeId(2), EdgeId(6)];
    for key in 0..7 {
        let nodes = if first.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if second.contains(&EdgeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        assert_eq!(relations.incident_to_node(NodeId(key)), nodes);
        assert_eq!(relations.incident_to_edge(EdgeId(key)), edges);
        assert_eq!(
            relations.coincident_to_node(NodeId(key), &first, &second),
            nodes.first().copied()
        );
        assert_eq!(
            relations.coincident_to_edge(EdgeId(key), &first, &second),
            edges.first().copied()
        );
    }
    assert_eq!(relations.into_entries(), vec![(first, second, 7)]);
}

#[rstest]
fn test_fixed_fixed_birelation_set_replace_participant_2_other_zero_arity() {
    let mut relations = FixedFixedBirelationSet::<NodeId, 0, NodeId, 2, _>::new(vec![(
        [],
        [NodeId(0), NodeId(1)],
        7,
    )]);
    relations.replace_participant_2(RelationId(0), ParticipantPosition(1), NodeId(5));
    assert_eq!(relations.incident_to_node(NodeId(0)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_node(NodeId(1)), &[]);
    assert_eq!(relations.incident_to_node(NodeId(5)), &[RelationId(0)]);
    assert_eq!(
        relations.into_entries(),
        vec![([], [NodeId(0), NodeId(5)], 7)]
    );
}

#[rstest]
#[case::past_end(RelationId(3), ParticipantPosition(0))]
#[case::maximum_id(RelationId(u32::MAX), ParticipantPosition(0))]
#[case::at_length(RelationId(0), ParticipantPosition(3))]
#[case::maximum_position(RelationId(0), ParticipantPosition(u32::MAX))]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participant_2_error(
    fixed_fixed_birelation_set_mutation_entries: Vec<([NodeId; 2], [NodeId; 3], BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
) {
    FixedFixedBirelationSet::new(fixed_fixed_birelation_set_mutation_entries)
        .replace_participant_2(id, position, NodeId(6));
}

#[rstest]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participant_2_zero_arity() {
    let mut relations =
        FixedFixedBirelationSet::<NodeId, 1, NodeId, 0, ()>::new(vec![([NodeId(1)], [], ())]);
    relations.replace_participant_2(RelationId(0), ParticipantPosition(0), NodeId(6));
}

#[rstest]
#[should_panic]
fn test_fixed_fixed_birelation_set_replace_participant_2_empty() {
    FixedFixedBirelationSet::<NodeId, 2, NodeId, 3, NonCloneData>::default().replace_participant_2(
        RelationId(0),
        ParticipantPosition(0),
        NodeId(6),
    );
}

#[rstest]
#[case::rows(FixedFixedBirelationSet::new(vec![([EdgeId(2), EdgeId(0)], [NodeId(2), NodeId(0)], vec![7, 11]), ([EdgeId(2), EdgeId(0)], [NodeId(2), NodeId(0)], vec![13, 17])]),
    FixedFixedBirelationSet::new(vec![([EdgeId(3), EdgeId(6)], [NodeId(1), NodeId(5)], vec![7, 11]), ([EdgeId(3), EdgeId(6)], [NodeId(1), NodeId(5)], vec![13, 17])]))]
fn test_fixed_fixed_birelation_set_map(
    participant_correspondence: GraphCorrespondence,
    #[case] input: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, Vec<u32>>,
    #[case] expected: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, Vec<u32>>,
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
    assert_eq!(
        expected.incident_to_edge(EdgeId(3)),
        expected.incident_to_node(NodeId(1))
    );
}

#[rstest]
#[case::empty(FixedFixedBirelationSet::new(vec![]))]
#[case::rows(FixedFixedBirelationSet::new(vec![([EdgeId(2), EdgeId(0)], [NodeId(2), NodeId(0)], vec![7, 11]), ([EdgeId(2), EdgeId(0)], [NodeId(2), NodeId(0)], vec![13, 17])]))]
fn test_fixed_fixed_birelation_set_map_identity(
    #[case] input: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, Vec<u32>>,
) {
    let identity = GraphCorrespondence::new(
        Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2), NodeId(3)], 4),
        Correspondence::from_images(&[EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)], 4),
    );
    assert_eq!(input.try_map(&identity), Some(input.clone()));
    assert_eq!(input.map(&identity), input);
}

#[rstest]
#[should_panic(expected = "correspondence must cover every participant reference")]
fn test_fixed_fixed_birelation_set_map_error(participant_correspondence: GraphCorrespondence) {
    let node = 1;
    let edge = 2;
    let input: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, Vec<u32>> =
        FixedFixedBirelationSet::new(vec![(
            [EdgeId(edge), EdgeId(0)],
            [NodeId(node), NodeId(0)],
            vec![7, 11],
        )]);
    input.map(&participant_correspondence);
}

#[rstest]
#[case::missing_node(1, 2)]
#[case::outside_node(4, 2)]
#[case::missing_edge(2, 1)]
#[case::outside_edge(2, 4)]
fn test_fixed_fixed_birelation_set_try_map_error(
    participant_correspondence: GraphCorrespondence,
    #[case] node: u32,
    #[case] edge: u32,
) {
    let input: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, Vec<u32>> =
        FixedFixedBirelationSet::new(vec![
            ([EdgeId(2), EdgeId(0)], [NodeId(2), NodeId(0)], vec![7, 11]),
            (
                [EdgeId(edge), EdgeId(0)],
                [NodeId(node), NodeId(0)],
                vec![13, 17],
            ),
        ]);
    assert_eq!(input.try_map(&participant_correspondence), None);
}

#[rstest]
fn test_fixed_fixed_birelation_set_remap() {
    let rs: FixedFixedBirelationSet<NodeId, 2, EdgeId, 2, BiPositionLabels> =
        FixedFixedBirelationSet::new(vec![(
            [NodeId(0), NodeId(1)],
            [EdgeId(0), EdgeId(1)],
            BiPositionLabels {
                factor_1: vec![10, 11],
                factor_2: vec![20, 21],
            },
        )]);
    let remapping = GraphRemapping::new(
        Remapping::new(vec![NodeId(1), NodeId(0)]).expect("permutation images"),
        Remapping::new(vec![EdgeId(1), EdgeId(0)]).expect("permutation images"),
    );
    let out = rs.remap(&remapping);
    assert_eq!(out.participants_1(RelationId(0)), &[NodeId(1), NodeId(0)]);
    assert_eq!(out.participants_2(RelationId(0)), &[EdgeId(1), EdgeId(0)]);
    assert_eq!(
        out.data(RelationId(0)),
        &BiPositionLabels {
            factor_1: vec![10, 11],
            factor_2: vec![20, 21],
        }
    );
}

#[rstest]
#[case::covered(vec![NodeId(1), NodeId(0)], vec![EdgeId(1), EdgeId(0)], true)]
#[case::uncovered_node(vec![NodeId(0)], vec![EdgeId(1), EdgeId(0)], false)]
#[case::uncovered_edge(vec![NodeId(1), NodeId(0)], vec![EdgeId(0)], false)]
fn test_fixed_fixed_birelation_set_try_remap(
    #[case] nodes: Vec<NodeId>,
    #[case] edges: Vec<EdgeId>,
    #[case] covered: bool,
) {
    let rs: FixedFixedBirelationSet<NodeId, 2, EdgeId, 2, BiPositionLabels> =
        FixedFixedBirelationSet::new(vec![(
            [NodeId(0), NodeId(1)],
            [EdgeId(0), EdgeId(1)],
            BiPositionLabels {
                factor_1: vec![10, 11],
                factor_2: vec![20, 21],
            },
        )]);
    let remapping = GraphRemapping::new(
        Remapping::new(nodes).expect("permutation images"),
        Remapping::new(edges).expect("permutation images"),
    );
    let expected = covered.then(|| rs.remap(&remapping));
    assert_eq!(rs.try_remap(&remapping), expected);
}

#[rstest]
#[case::empty(FixedFixedBirelationSet::default())]
#[case::rows(
    FixedFixedBirelationSet::new(vec![([NodeId(0)], [NodeId(2), NodeId(4)], "keep"), ([NodeId(1)], [NodeId(5), NodeId(6)], "drop")]),
)]
fn test_fixed_fixed_birelation_set_compact_identity(
    #[case] input: FixedFixedBirelationSet<NodeId, 1, NodeId, 2, &'static str>,
) {
    let compaction = GraphCompaction::new(Compaction::identity(7), Compaction::empty());
    assert_eq!(input.compact(&compaction), input);
    assert_eq!(
        input.tracked_compact(&compaction),
        (input.clone(), Compaction::identity(input.count())),
    );
}

#[rstest]
#[case::partial(
    vec![NodeId(1)],
    FixedFixedBirelationSet::new(vec![([NodeId(0)], [NodeId(1), NodeId(3)], "keep")]),
    vec![RelationId(1)],
)]
#[case::all(
    vec![NodeId(0), NodeId(1)],
    FixedFixedBirelationSet::default(),
    vec![RelationId(0), RelationId(1)],
)]
fn test_fixed_fixed_birelation_set_tracked_compact(
    fixed_fixed_birelation_set_compaction_input: FixedFixedBirelationSet<
        NodeId,
        1,
        NodeId,
        2,
        &'static str,
    >,
    #[case] removed_nodes: Vec<NodeId>,
    #[case] expected: FixedFixedBirelationSet<NodeId, 1, NodeId, 2, &'static str>,
    #[case] removed_relations: Vec<RelationId>,
) {
    let input = fixed_fixed_birelation_set_compaction_input;
    let compaction = GraphCompaction::new(
        Compaction::new(7, removed_nodes).unwrap(),
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

#[rstest]
fn test_fixed_fixed_birelation_set_tracked_pushout() {
    let left = FixedFixedBirelationSet::<NodeId, 1, NodeId, 2, i32>::new(vec![(
        [NodeId(0)],
        [NodeId(1), NodeId(2)],
        10,
    )]);
    let right = FixedFixedBirelationSet::<NodeId, 1, NodeId, 2, i32>::new(vec![
        ([NodeId(0)], [NodeId(1), NodeId(2)], 5),
        ([NodeId(3)], [NodeId(4), NodeId(5)], 20),
    ]);
    let (object, glue) = left
        .tracked_pushout(
            &right,
            |set: &_, q1: &[NodeId], q2: &_| {
                q1.first().and_then(|&n| set.coincident_to_node(n, q1, q2))
            },
            |(_, _, a), (_, _, b)| Some(a + b),
        )
        .expect("no ⊥");
    assert_eq!(
        left.pushout(
            &right,
            |set: &_, q1: &[NodeId], q2: &_| {
                q1.first().and_then(|&n| set.coincident_to_node(n, q1, q2))
            },
            |(_, _, a), (_, _, b)| Some(a + b),
        ),
        Some(object.clone()),
    );
    assert_eq!(
        object,
        FixedFixedBirelationSet::new(vec![
            ([NodeId(0)], [NodeId(1), NodeId(2)], 15),
            ([NodeId(3)], [NodeId(4), NodeId(5)], 20)
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
            |set: &_, q1: &[NodeId], q2: &_| {
                q1.first().and_then(|&n| set.coincident_to_node(n, q1, q2))
            },
            |_, _| None,
        ),
        None
    );
    assert_eq!(
        left.pushout(
            &right,
            |set: &_, q1: &[NodeId], q2: &_| {
                q1.first().and_then(|&n| set.coincident_to_node(n, q1, q2))
            },
            |_, _| None,
        ),
        None
    );
}

#[rstest]
#[case::combined(Some(15))]
#[case::incompatible(None)]
fn test_fixed_fixed_birelation_set_tracked_pullback(#[case] combined: Option<i32>) {
    let left: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, i32> =
        FixedFixedBirelationSet::new(vec![
            ([NodeId(0)], [NodeId(1)], 10),
            ([NodeId(2)], [NodeId(3)], 20),
        ]);
    let right: FixedFixedBirelationSet<NodeId, 1, NodeId, 1, i32> =
        FixedFixedBirelationSet::new(vec![
            ([NodeId(4)], [NodeId(5)], 30),
            ([NodeId(0)], [NodeId(1)], 5),
        ]);
    let result = left.tracked_pullback(
        &right,
        |set, first: &[NodeId], second: &[NodeId]| {
            first
                .first()
                .and_then(|&id| set.coincident_to_node(id, first, second))
        },
        |(_, _, a), (_, _, b)| combined.map(|_| a + b),
    );
    let plain = left.pullback(
        &right,
        |set, first: &[NodeId], second: &[NodeId]| {
            first
                .first()
                .and_then(|&id| set.coincident_to_node(id, first, second))
        },
        |(_, _, a), (_, _, b)| combined.map(|_| a + b),
    );
    let expected = combined.map(|value| {
        (
            FixedFixedBirelationSet::new(vec![([NodeId(0)], [NodeId(1)], value)]),
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
fn test_fixed_fixed_birelation_set_default() {
    let rs = FixedFixedBirelationSet::<NodeId, 1, NodeId, 1, ()>::default();
    assert_eq!(rs.count(), 0);
    assert!(!rs.has_incident_to_node(NodeId(0)));
}
