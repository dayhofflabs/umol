use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::panic::{catch_unwind, AssertUnwindSafe};

use pretty_assertions::assert_eq;
use rstest::{fixture, rstest};

use crate::{
    Compaction, Correspondence, EdgeId, FixedVarBirelationSet, Graph, GraphCompaction,
    GraphCorrespondence, GraphRemapping, NodeId, ParticipantPosition, RelationId,
    RelationPullbackCorrespondence, RelationPushoutCorrespondence, Remapping,
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
struct BiPositionLabels {
    factor_1: Vec<u32>,
    factor_2: Vec<u32>,
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

#[derive(Debug, PartialEq, Eq)]
struct NonCloneData(Vec<u32>);

#[fixture]
fn fixed_var_birelation_set_mutation_entries() -> Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>
{
    vec![
        (
            [NodeId(0), NodeId(1)],
            vec![NodeId(1), NodeId(2), NodeId(1)],
            BiPositionLabels {
                factor_1: vec![7, 11],
                factor_2: vec![13, 17, 19],
            },
        ),
        (
            [NodeId(2), NodeId(3)],
            vec![],
            BiPositionLabels {
                factor_1: vec![23, 29],
                factor_2: vec![31, 37, 41],
            },
        ),
        (
            [NodeId(4), NodeId(5)],
            vec![NodeId(5)],
            BiPositionLabels {
                factor_1: vec![43, 47],
                factor_2: vec![53, 59, 61],
            },
        ),
    ]
}

fn assert_fixed_var_birelation_rows<const N1: usize>(
    relations: &FixedVarBirelationSet<NodeId, N1, NodeId, BiPositionLabels>,
    entries: &[([NodeId; N1], Vec<NodeId>, BiPositionLabels)],
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
                        && b.len() == query_2.len()
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
fn fixed_var_birelation_set_compaction_input(
) -> FixedVarBirelationSet<NodeId, 1, NodeId, &'static str> {
    FixedVarBirelationSet::new(vec![
        ([NodeId(0)], vec![NodeId(2), NodeId(4)], "keep"),
        ([NodeId(5)], vec![NodeId(1), NodeId(3)], "drop"),
    ])
}

#[rstest]
fn test_fixed_var_birelation_set_new() {
    let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> = FixedVarBirelationSet::new(vec![(
        [EdgeId(0)],
        vec![NodeId(1), NodeId(2), NodeId(3)],
        "ct",
    )]);
    assert_eq!(rs.count(), 1);
    assert_eq!(rs.participants_1(RelationId(0)), &[EdgeId(0)]);
    assert_eq!(
        rs.participants_2(RelationId(0)),
        &[NodeId(1), NodeId(2), NodeId(3)]
    );
    assert_eq!(rs.data(RelationId(0)), &"ct");
}

#[rstest]
#[case::empty(vec![], vec![], vec![])]
#[case::repeated(
    vec![
        ([NodeId(2), NodeId(0), NodeId(2)], vec![NodeId(2), NodeId(1), NodeId(2)], "first"),
        ([NodeId(0), NodeId(2), NodeId(2)], vec![NodeId(1), NodeId(2), NodeId(2)], "duplicate"),
        ([NodeId(3), NodeId(3), NodeId(3)], vec![NodeId(3), NodeId(3), NodeId(3)], "other"),
    ],
    vec![RelationId(0), RelationId(1)],
    vec![RelationId(2)],
)]
#[case::empty_second_factor(
    vec![([NodeId(2), NodeId(2), NodeId(2)], vec![], "empty")],
    vec![RelationId(0)],
    vec![],
)]
fn test_fixed_var_birelation_set_new_incidence(
    #[case] entries: Vec<([NodeId; 3], Vec<NodeId>, &str)>,
    #[case] at_two: Vec<RelationId>,
    #[case] at_three: Vec<RelationId>,
) {
    let relations = FixedVarBirelationSet::<NodeId, 3, NodeId, &str>::new(entries.clone());
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
        ([EdgeId(2)], vec![NodeId(3), NodeId(1)], "first"),
        ([EdgeId(4)], vec![NodeId(5), NodeId(0)], "second"),
    ],
)]
fn test_fixed_var_birelation_set_into_entries(
    #[case] entries: Vec<([EdgeId; 1], Vec<NodeId>, &str)>,
) {
    let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> =
        FixedVarBirelationSet::new(entries.clone());
    assert_eq!(rs.into_entries(), entries);
}

#[rstest]
#[case::first(RelationId(0), true)]
#[case::out_of_range(RelationId(1), false)]
fn test_fixed_var_birelation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
    let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> =
        FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![NodeId(1)], "ct")]);
    assert_eq!(rs.contains(id), expected);
}

#[rstest]
fn test_fixed_var_birelation_set_relation_ids() {
    assert_exact_size(
        FixedVarBirelationSet::<EdgeId, 1, NodeId, &str>::default().ids(),
        vec![],
    );
    let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> = FixedVarBirelationSet::new(vec![
        ([EdgeId(0)], vec![NodeId(1)], "a"),
        ([EdgeId(1)], vec![NodeId(2)], "b"),
    ]);
    assert_exact_size(rs.ids(), vec![RelationId(0), RelationId(1)]);
}

#[rstest]
fn test_fixed_var_birelation_set_iter() {
    let empty = FixedVarBirelationSet::<EdgeId, 1, NodeId, i32>::default();
    assert_eq!(empty.iter().collect::<Vec<_>>(), vec![]);

    let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, i32> = FixedVarBirelationSet::new(vec![
        ([EdgeId(0)], vec![NodeId(1), NodeId(3)], 1),
        ([EdgeId(1)], vec![NodeId(2)], 2),
    ]);
    assert_eq!(rs.iter().len(), 2);
    assert_eq!(
        rs.iter().collect::<Vec<_>>(),
        vec![
            (
                RelationId(0),
                &[EdgeId(0)],
                [NodeId(1), NodeId(3)].as_slice(),
                &1
            ),
            (RelationId(1), &[EdgeId(1)], [NodeId(2)].as_slice(), &2),
        ],
    );
}

#[rstest]
fn test_fixed_var_birelation_set_iter_mut() {
    let mut empty = FixedVarBirelationSet::<EdgeId, 1, NodeId, i32>::default();
    assert_eq!(empty.iter_mut().len(), 0);

    let mut rs: FixedVarBirelationSet<EdgeId, 1, NodeId, i32> = FixedVarBirelationSet::new(vec![
        ([EdgeId(0)], vec![NodeId(1), NodeId(3)], 1),
        ([EdgeId(1)], vec![NodeId(2)], 2),
    ]);
    let arities: Vec<usize> = rs
        .iter_mut()
        .map(|(_, _, ligands, data)| {
            *data *= 10;
            ligands.len()
        })
        .collect();
    assert_eq!(arities, vec![2, 1]);
    assert_eq!(rs.data(RelationId(0)), &10);
    assert_eq!(rs.data(RelationId(1)), &20);
}

#[rstest]
fn test_fixed_var_birelation_set_data_mut() {
    let mut rs: FixedVarBirelationSet<EdgeId, 1, NodeId, i32> =
        FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![NodeId(1)], 1)]);
    *rs.data_mut(RelationId(0)) = 99;
    assert_eq!(rs.data(RelationId(0)), &99);
}

#[rstest]
fn test_fixed_var_birelation_set_incidence() {
    let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> =
        FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![NodeId(1), NodeId(2)], "ct")]);
    assert_eq!(rs.incident_to_edge(EdgeId(0)), &[RelationId(0)]);
    assert_eq!(rs.incident_to_node(NodeId(2)), &[RelationId(0)]);
    assert!(rs.has_incident_to_edge(EdgeId(0)));
    assert!(rs.has_incident_to_node(NodeId(1)));
    assert!(rs.incident_to_node(NodeId(0)).is_empty());
    assert!(!rs.has_incident_to_node(NodeId(0)));
}

#[rstest]
#[case::exact(vec![NodeId(0)], vec![NodeId(1)], Some(RelationId(0)))]
#[case::role_swap(vec![NodeId(1)], vec![NodeId(0)], None)]
#[case::multiset_reordered(vec![NodeId(3)], vec![NodeId(5), NodeId(4), NodeId(4)], Some(RelationId(1)))]
#[case::wrong_multiplicity(vec![NodeId(3)], vec![NodeId(4), NodeId(5)], None)]
#[case::absent(vec![NodeId(0)], vec![NodeId(2)], None)]
fn test_fixed_var_birelation_set_coincident_to_node(
    #[case] query_1: Vec<NodeId>,
    #[case] query_2: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
) {
    // Each factor matches as a multiset; repeated values and factor roles remain distinct.
    let rs: FixedVarBirelationSet<NodeId, 1, NodeId, ()> = FixedVarBirelationSet::new(vec![
        ([NodeId(0)], vec![NodeId(1)], ()),
        ([NodeId(3)], vec![NodeId(4), NodeId(4), NodeId(5)], ()),
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
fn test_fixed_var_birelation_set_coincident_to_node_multiplicity(
    #[case] anchor: NodeId,
    #[case] query: Vec<NodeId>,
    #[case] query_2: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
    #[case] coincides: bool,
) {
    let relations = FixedVarBirelationSet::<NodeId, 3, NodeId, &str>::new(vec![
        (
            [NodeId(2), NodeId(0), NodeId(2)],
            vec![NodeId(2), NodeId(1), NodeId(2)],
            "first",
        ),
        (
            [NodeId(0), NodeId(2), NodeId(2)],
            vec![NodeId(1), NodeId(2), NodeId(2)],
            "duplicate",
        ),
        (
            [NodeId(3), NodeId(3), NodeId(3)],
            vec![NodeId(3), NodeId(3), NodeId(3)],
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
#[case::exact(vec![NodeId(1), NodeId(2)], Some(RelationId(0)))]
#[case::reordered(vec![NodeId(2), NodeId(1)], Some(RelationId(0)))]
#[case::absent(vec![NodeId(1), NodeId(3)], None)]
fn test_fixed_var_birelation_set_find_by_participants_edge_anchor(
    #[case] ligands: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
) {
    // Stereo-bond-like: factor1 is an `EdgeId` site, so the anchor routes through `incident_to_edge`.
    let rs: FixedVarBirelationSet<EdgeId, 1, NodeId, ()> =
        FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![NodeId(1), NodeId(2)], ())]);
    assert_eq!(
        rs.coincident_to_edge(EdgeId(0), &[EdgeId(0)], &ligands),
        expected,
    );
}

#[rstest]
#[case::empty(vec![], [NodeId(2), NodeId(2)], vec![NodeId(2), NodeId(3)])]
#[case::shared(vec![([NodeId(0), NodeId(2)], vec![NodeId(2), NodeId(1)], "old")], [NodeId(2), NodeId(2)], vec![NodeId(2), NodeId(3), NodeId(2)])]
#[case::coinciding(vec![([NodeId(2)], vec![NodeId(2), NodeId(3)], "old")], [NodeId(2)], vec![NodeId(2), NodeId(3)])]
#[case::reordered(vec![([NodeId(0), NodeId(1)], vec![NodeId(2), NodeId(3)], "old")], [NodeId(1), NodeId(0)], vec![NodeId(3), NodeId(2)])]
#[case::sparse(vec![], [NodeId(u32::MAX)], vec![NodeId(u32::MAX), NodeId(0)])]
#[case::first_empty(vec![([], vec![NodeId(0), NodeId(1)], "old")], [], vec![NodeId(1)])]
#[case::second_empty(vec![([NodeId(0)], vec![NodeId(1), NodeId(2)], "old")], [NodeId(1)], vec![])]
#[case::adjacent_empty(vec![([NodeId(0)], vec![], "first"), ([NodeId(1)], vec![], "second")], [NodeId(2)], vec![])]
#[case::after_empty(vec![([NodeId(0)], vec![], "first"), ([NodeId(1)], vec![], "second")], [NodeId(2)], vec![NodeId(3), NodeId(2)])]
#[case::both_empty(vec![([], vec![], "first")], [], vec![])]
fn test_fixed_var_birelation_set_add<const N1: usize>(
    #[case] entries: Vec<([NodeId; N1], Vec<NodeId>, &'static str)>,
    #[case] first: [NodeId; N1],
    #[case] second: Vec<NodeId>,
) {
    let mut relations = FixedVarBirelationSet::new(entries.clone());
    assert_eq!(
        relations.add(first, &second, "added"),
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
fn test_fixed_var_birelation_set_add_payload() {
    let mut relations = FixedVarBirelationSet::default();
    assert_eq!(
        relations.add(
            [EdgeId(0)],
            &[NodeId(2), NodeId(0)],
            NonCloneData(vec![7, 11])
        ),
        RelationId(0)
    );
    assert_eq!(
        relations.add([EdgeId(1)], &[], NonCloneData(vec![13])),
        RelationId(1)
    );
    assert_eq!(
        relations.add(
            [EdgeId(2)],
            &[NodeId(2), NodeId(2)],
            NonCloneData(vec![17, 19])
        ),
        RelationId(2)
    );
    assert_eq!(relations.incident_to_edge(EdgeId(0)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_edge(EdgeId(1)), &[RelationId(1)]);
    assert_eq!(relations.incident_to_edge(EdgeId(2)), &[RelationId(2)]);
    assert_eq!(relations.incident_to_node(NodeId(0)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_node(NodeId(1)), &[]);
    assert_eq!(
        relations.incident_to_node(NodeId(2)),
        &[RelationId(0), RelationId(2)]
    );
    assert_eq!(
        relations.into_entries(),
        vec![
            (
                [EdgeId(0)],
                vec![NodeId(2), NodeId(0)],
                NonCloneData(vec![7, 11])
            ),
            ([EdgeId(1)], vec![], NonCloneData(vec![13])),
            (
                [EdgeId(2)],
                vec![NodeId(2), NodeId(2)],
                NonCloneData(vec![17, 19])
            ),
        ]
    );
}

#[rstest]
#[case::first(vec![RelationId(0)], vec![1, 2])]
#[case::empty_second_factor(vec![RelationId(1)], vec![0, 2])]
#[case::last(vec![RelationId(2)], vec![0, 1])]
#[case::unordered_repeated(vec![RelationId(2), RelationId(0), RelationId(2)], vec![1])]
#[case::all(vec![RelationId(2), RelationId(0), RelationId(1)], vec![])]
fn test_fixed_var_birelation_set_tracked_remove(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] ids: Vec<RelationId>,
    #[case] survivors: Vec<usize>,
) {
    let entries = fixed_var_birelation_set_mutation_entries;
    let expected: Vec<_> = survivors.iter().map(|&i| entries[i].clone()).collect();
    let mut relations = FixedVarBirelationSet::new(entries);
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
    assert_fixed_var_birelation_rows(&relations, &expected);
    assert_eq!(relations.into_entries(), expected);
}

#[rstest]
#[case::empty(FixedVarBirelationSet::<NodeId, 0, NodeId, &str>::default())]
#[case::both_empty(FixedVarBirelationSet::<NodeId, 0, NodeId, &str>::new(vec![([], vec![], "first"), ([], vec![], "second")]))]
#[case::nonempty(FixedVarBirelationSet::new(vec![([NodeId(0)], vec![NodeId(2), NodeId(0)], "first")]))]
fn test_fixed_var_birelation_set_tracked_remove_identity<const N1: usize>(
    #[case] input: FixedVarBirelationSet<NodeId, N1, NodeId, &'static str>,
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
#[case::first_empty([], vec![NodeId(0), NodeId(1)])]
#[case::second_empty([NodeId(0), NodeId(1)], vec![])]
#[case::both_empty([], vec![])]
#[case::overlap([NodeId(0)], vec![NodeId(0), NodeId(1)])]
fn test_fixed_var_birelation_set_tracked_remove_arity<const N1: usize>(
    #[case] first: [NodeId; N1],
    #[case] second: Vec<NodeId>,
) {
    let mut relations = FixedVarBirelationSet::new(vec![
        (first, vec![], "first"),
        (first, second, "second"),
        (first, vec![], "third"),
        (first, vec![NodeId(2)], "fourth"),
    ]);
    let compaction = relations.tracked_remove(&[RelationId(1), RelationId(2)]);
    assert_eq!(compaction.source_count(), 4);
    assert_eq!(compaction.result_count(), 2);
    assert_eq!(compaction.removed(), &[RelationId(1), RelationId(2)]);
    let expected = vec![(first, vec![], "first"), (first, vec![NodeId(2)], "fourth")];
    for (i, (a, b, data)) in expected.iter().enumerate() {
        assert_eq!(relations.participants_1(RelationId::from(i)), a);
        assert_eq!(relations.participants_2(RelationId::from(i)), b);
        assert_eq!(relations.data(RelationId::from(i)), data);
    }
    for node in [NodeId(0), NodeId(1), NodeId(2)] {
        let incidence: Vec<_> = expected
            .iter()
            .enumerate()
            .filter(|(_, (a, b, _))| a.contains(&node) || b.contains(&node))
            .map(|(i, _)| RelationId::from(i))
            .collect();
        assert_eq!(relations.incident_to_node(node), incidence);
    }
    assert_eq!(relations.into_entries(), expected);
}

#[rstest]
fn test_fixed_var_birelation_set_tracked_remove_payload() {
    let mut relations = FixedVarBirelationSet::default();
    relations.add(
        [EdgeId(0)],
        &[NodeId(2), NodeId(0)],
        NonCloneData(vec![7, 11]),
    );
    relations.add([EdgeId(1)], &[], NonCloneData(vec![13]));
    relations.add(
        [EdgeId(2)],
        &[NodeId(2), NodeId(2), NodeId(3)],
        NonCloneData(vec![17, 19, 23]),
    );
    let compaction = relations.tracked_remove(&[RelationId(0)]);
    assert_eq!(compaction.source_count(), 3);
    assert_eq!(compaction.result_count(), 2);
    assert_eq!(compaction.removed(), &[RelationId(0)]);
    assert_eq!(relations.incident_to_edge(EdgeId(0)), &[]);
    assert_eq!(relations.incident_to_edge(EdgeId(1)), &[RelationId(0)]);
    assert_eq!(relations.incident_to_edge(EdgeId(2)), &[RelationId(1)]);
    assert_eq!(relations.incident_to_node(NodeId(0)), &[]);
    assert_eq!(relations.incident_to_node(NodeId(2)), &[RelationId(1)]);
    assert_eq!(relations.incident_to_node(NodeId(3)), &[RelationId(1)]);
    assert_eq!(
        relations.into_entries(),
        vec![
            ([EdgeId(1)], vec![], NonCloneData(vec![13])),
            (
                [EdgeId(2)],
                vec![NodeId(2), NodeId(2), NodeId(3)],
                NonCloneData(vec![17, 19, 23])
            ),
        ]
    );
}

#[rstest]
#[case::empty(vec![], vec![RelationId(0)])]
#[case::end(vec![([NodeId(0)], vec![], "first")], vec![RelationId(1)])]
#[case::mixed(vec![([NodeId(0)], vec![NodeId(2)], "first"), ([NodeId(2)], vec![], "second")], vec![RelationId(0), RelationId(2)])]
#[case::distant(vec![([NodeId(0)], vec![NodeId(2)], "first")], vec![RelationId(u32::MAX), RelationId(0)])]
fn test_fixed_var_birelation_set_tracked_remove_error(
    #[case] entries: Vec<([NodeId; 1], Vec<NodeId>, &'static str)>,
    #[case] ids: Vec<RelationId>,
) {
    let original = FixedVarBirelationSet::new(entries);
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
fn test_fixed_var_birelation_set_permute_participants_1() {
    let mut rs: FixedVarBirelationSet<EdgeId, 2, NodeId, &str> = FixedVarBirelationSet::new(vec![
        (
            [EdgeId(4), EdgeId(5)],
            vec![NodeId(0), NodeId(1), NodeId(2)],
            "a",
        ),
        ([EdgeId(6), EdgeId(7)], vec![NodeId(3), NodeId(4)], "b"),
    ]);
    rs.permute_participants_1(
        RelationId(1),
        &[ParticipantPosition(1), ParticipantPosition(0)],
    );
    assert_eq!(rs.participants_1(RelationId(0)), &[EdgeId(4), EdgeId(5)]);
    assert_eq!(rs.participants_1(RelationId(1)), &[EdgeId(7), EdgeId(6)]);
    assert_eq!(rs.participants_2(RelationId(1)), &[NodeId(3), NodeId(4)]);
    assert_eq!(rs.data(RelationId(1)), &"b");
}

#[rstest]
fn test_fixed_var_birelation_set_permute_participants_1_identity() {
    let input: FixedVarBirelationSet<EdgeId, 2, NodeId, &str> = FixedVarBirelationSet::new(vec![(
        [EdgeId(5), EdgeId(4)],
        vec![NodeId(2), NodeId(0), NodeId(1)],
        "a",
    )]);
    let mut permuted = input.clone();
    permuted.permute_participants_1(
        RelationId(0),
        &[ParticipantPosition(0), ParticipantPosition(1)],
    );
    permuted.permute_participants_2(
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
fn test_fixed_var_birelation_set_permute_participants_2() {
    let mut rs: FixedVarBirelationSet<EdgeId, 2, NodeId, &str> = FixedVarBirelationSet::new(vec![
        (
            [EdgeId(4), EdgeId(5)],
            vec![NodeId(0), NodeId(1), NodeId(2)],
            "a",
        ),
        ([EdgeId(6), EdgeId(7)], vec![NodeId(3), NodeId(4)], "b"),
    ]);
    let incidence_before: Vec<Vec<RelationId>> = (0..5)
        .map(|i| rs.incident_to_node(NodeId(i)).to_vec())
        .collect();

    rs.permute_participants_2(
        RelationId(0),
        &[
            ParticipantPosition(2),
            ParticipantPosition(0),
            ParticipantPosition(1),
        ],
    );

    assert_eq!(
        rs.participants_2(RelationId(0)),
        &[NodeId(2), NodeId(0), NodeId(1)]
    );
    assert_eq!(rs.participants_2(RelationId(1)), &[NodeId(3), NodeId(4)]);
    assert_eq!(rs.participants_1(RelationId(0)), &[EdgeId(4), EdgeId(5)]);
    assert_eq!(rs.data(RelationId(0)), &"a");
    let incidence_after: Vec<Vec<RelationId>> = (0..5)
        .map(|i| rs.incident_to_node(NodeId(i)).to_vec())
        .collect();
    assert_eq!(incidence_after, incidence_before);
}

#[rstest]
#[case::order_too_short(vec![ParticipantPosition(0), ParticipantPosition(1)])]
#[case::position_out_of_range(vec![ParticipantPosition(0), ParticipantPosition(1), ParticipantPosition(3)])]
#[case::position_repeated(vec![ParticipantPosition(2), ParticipantPosition(2), ParticipantPosition(0)])]
#[should_panic(expected = "permute")]
fn test_fixed_var_birelation_set_permute_participants_2_error(
    #[case] order: Vec<ParticipantPosition>,
) {
    let mut rs: FixedVarBirelationSet<EdgeId, 2, NodeId, &str> =
        FixedVarBirelationSet::new(vec![(
            [EdgeId(4), EdgeId(5)],
            vec![NodeId(0), NodeId(1), NodeId(2)],
            "a",
        )]);
    rs.permute_participants_2(RelationId(0), &order);
}

#[rstest]
#[case::reordered(RelationId(0), [NodeId(1), NodeId(0)], vec![NodeId(1), NodeId(1), NodeId(2)])]
#[case::empty_first_row(RelationId(0), [NodeId(6), NodeId(7)], vec![])]
#[case::grow_middle(RelationId(1), [NodeId(6), NodeId(6)], vec![NodeId(6), NodeId(6), NodeId(7), NodeId(6)])]
#[case::shrink_first(RelationId(0), [NodeId(6), NodeId(7)], vec![NodeId(1)])]
#[case::sparse_last(RelationId(2), [NodeId(u32::MAX), NodeId(5)], vec![NodeId(u32::MAX), NodeId(0)])]
fn test_fixed_var_birelation_set_replace_participants(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] first: [NodeId; 2],
    #[case] second: Vec<NodeId>,
) {
    let mut entries = fixed_var_birelation_set_mutation_entries;
    let mut relations = FixedVarBirelationSet::new(entries.clone());
    relations.replace_participants(id, first, &second);
    entries[id.index()].0 = first;
    entries[id.index()].1 = second;
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::both_empty([], vec![])]
#[case::first_empty([], vec![NodeId(1), NodeId(2), NodeId(1)])]
#[case::second_empty([NodeId(0), NodeId(1)], vec![])]
#[case::both_present([NodeId(1), NodeId(0)], vec![NodeId(2), NodeId(1), NodeId(1)])]
fn test_fixed_var_birelation_set_replace_participants_identity<const N1: usize>(
    #[case] first: [NodeId; N1],
    #[case] second: Vec<NodeId>,
) {
    let entries = vec![(
        first,
        second.clone(),
        BiPositionLabels {
            factor_1: vec![7, 11],
            factor_2: vec![13, 17, 19],
        },
    )];
    let input = FixedVarBirelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participants(RelationId(0), first, &second);
    assert_eq!(relations, input);
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::grow(vec![], vec![NodeId(1), NodeId(1), NodeId(2)])]
#[case::empty(vec![NodeId(1), NodeId(2)], vec![])]
fn test_fixed_var_birelation_set_replace_participants_zero_arity(
    #[case] second: Vec<NodeId>,
    #[case] replacement: Vec<NodeId>,
) {
    let mut entries = vec![
        (
            [],
            second,
            BiPositionLabels {
                factor_1: vec![],
                factor_2: vec![7, 11],
            },
        ),
        (
            [],
            vec![NodeId(3)],
            BiPositionLabels {
                factor_1: vec![],
                factor_2: vec![13],
            },
        ),
    ];
    let mut relations = FixedVarBirelationSet::<NodeId, 0, NodeId, _>::new(entries.clone());
    relations.replace_participants(RelationId(0), [], &replacement);
    entries[0].1 = replacement;
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_var_birelation_set_replace_participants_edges() {
    let mut relations =
        FixedVarBirelationSet::new(vec![([EdgeId(1)], vec![NodeId(1), NodeId(2)], 7)]);
    relations.replace_participants(
        RelationId(0),
        [EdgeId(5)],
        &[NodeId(6), NodeId(6), NodeId(6)],
    );
    let first = [EdgeId(5)];
    let second = vec![NodeId(6), NodeId(6), NodeId(6)];
    for key in 0..7 {
        let nodes = if second.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if first.contains(&EdgeId(key)) {
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
#[case::maximum_id(RelationId(u32::MAX))]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participants_error(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
) {
    FixedVarBirelationSet::new(fixed_var_birelation_set_mutation_entries).replace_participants(
        id,
        [NodeId(6), NodeId(7)],
        &[NodeId(6)],
    );
}

#[rstest]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participants_empty() {
    FixedVarBirelationSet::<NodeId, 2, NodeId, NonCloneData>::default().replace_participants(
        RelationId(0),
        [NodeId(6), NodeId(7)],
        &[NodeId(6)],
    );
}

#[rstest]
#[case::shared_second(RelationId(0), [NodeId(6), NodeId(7)])]
#[case::reordered(RelationId(0), [NodeId(1), NodeId(0)])]
#[case::empty_second(RelationId(1), [NodeId(6), NodeId(6)])]
#[case::sparse_last(RelationId(2), [NodeId(u32::MAX), NodeId(5)])]
fn test_fixed_var_birelation_set_replace_participants_1(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] first: [NodeId; 2],
) {
    let mut entries = fixed_var_birelation_set_mutation_entries;
    let mut relations = FixedVarBirelationSet::new(entries.clone());
    relations.replace_participants_1(id, first);
    entries[id.index()].0 = first;
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::both_empty([], vec![])]
#[case::first_empty([], vec![NodeId(1), NodeId(2), NodeId(1)])]
#[case::second_empty([NodeId(0), NodeId(1)], vec![])]
#[case::both_present([NodeId(1), NodeId(0)], vec![NodeId(2), NodeId(1), NodeId(1)])]
fn test_fixed_var_birelation_set_replace_participants_1_identity<const N1: usize>(
    #[case] first: [NodeId; N1],
    #[case] second: Vec<NodeId>,
) {
    let entries = vec![(
        first,
        second.clone(),
        BiPositionLabels {
            factor_1: vec![7, 11],
            factor_2: vec![13, 17, 19],
        },
    )];
    let input = FixedVarBirelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participants_1(RelationId(0), first);
    assert_eq!(relations, input);
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_var_birelation_set_replace_participants_1_edges() {
    let mut relations =
        FixedVarBirelationSet::new(vec![([EdgeId(1)], vec![NodeId(1), NodeId(2)], 7)]);
    relations.replace_participants_1(RelationId(0), [EdgeId(5)]);
    let first = [EdgeId(5)];
    let second = vec![NodeId(1), NodeId(2)];
    for key in 0..7 {
        let nodes = if second.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if first.contains(&EdgeId(key)) {
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
#[case::maximum_id(RelationId(u32::MAX))]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participants_1_error(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
) {
    FixedVarBirelationSet::new(fixed_var_birelation_set_mutation_entries)
        .replace_participants_1(id, [NodeId(6), NodeId(7)]);
}

#[rstest]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participants_1_empty() {
    FixedVarBirelationSet::<NodeId, 2, NodeId, NonCloneData>::default()
        .replace_participants_1(RelationId(0), [NodeId(6), NodeId(7)]);
}

#[rstest]
#[case::shared_first(RelationId(0), vec![NodeId(6), NodeId(7)])]
#[case::reordered(RelationId(0), vec![NodeId(1), NodeId(1), NodeId(2)])]
#[case::empty_first_row(RelationId(0), vec![])]
#[case::grow_middle(RelationId(1), vec![NodeId(6), NodeId(6), NodeId(7), NodeId(6)])]
#[case::sparse_last(RelationId(2), vec![NodeId(u32::MAX), NodeId(5), NodeId(0)])]
fn test_fixed_var_birelation_set_replace_participants_2(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] second: Vec<NodeId>,
) {
    let mut entries = fixed_var_birelation_set_mutation_entries;
    let mut relations = FixedVarBirelationSet::new(entries.clone());
    relations.replace_participants_2(id, &second);
    entries[id.index()].1 = second;
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::both_empty([], vec![])]
#[case::first_empty([], vec![NodeId(1), NodeId(2), NodeId(1)])]
#[case::second_empty([NodeId(0), NodeId(1)], vec![])]
#[case::both_present([NodeId(1), NodeId(0)], vec![NodeId(2), NodeId(1), NodeId(1)])]
fn test_fixed_var_birelation_set_replace_participants_2_identity<const N1: usize>(
    #[case] first: [NodeId; N1],
    #[case] second: Vec<NodeId>,
) {
    let entries = vec![(
        first,
        second.clone(),
        BiPositionLabels {
            factor_1: vec![7, 11],
            factor_2: vec![13, 17, 19],
        },
    )];
    let input = FixedVarBirelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participants_2(RelationId(0), &second);
    assert_eq!(relations, input);
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_var_birelation_set_replace_participants_2_edges() {
    let mut relations =
        FixedVarBirelationSet::new(vec![([EdgeId(1)], vec![NodeId(1), NodeId(2)], 7)]);
    relations.replace_participants_2(RelationId(0), &[NodeId(6), NodeId(6), NodeId(6)]);
    let first = [EdgeId(1)];
    let second = vec![NodeId(6), NodeId(6), NodeId(6)];
    for key in 0..7 {
        let nodes = if second.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if first.contains(&EdgeId(key)) {
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
fn test_fixed_var_birelation_set_replace_participants_2_zero_arity() {
    let mut relations = FixedVarBirelationSet::<NodeId, 0, NodeId, _>::new(vec![(
        [],
        vec![NodeId(1), NodeId(2)],
        BiPositionLabels {
            factor_1: vec![],
            factor_2: vec![7, 11],
        },
    )]);
    relations.replace_participants_2(RelationId(0), &[NodeId(6), NodeId(6), NodeId(6)]);
    let entries = vec![(
        [],
        vec![NodeId(6), NodeId(6), NodeId(6)],
        BiPositionLabels {
            factor_1: vec![],
            factor_2: vec![7, 11],
        },
    )];
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::past_end(RelationId(3))]
#[case::maximum_id(RelationId(u32::MAX))]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participants_2_error(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
) {
    FixedVarBirelationSet::new(fixed_var_birelation_set_mutation_entries)
        .replace_participants_2(id, &[NodeId(6)]);
}

#[rstest]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participants_2_empty() {
    FixedVarBirelationSet::<NodeId, 2, NodeId, NonCloneData>::default()
        .replace_participants_2(RelationId(0), &[NodeId(6)]);
}

#[rstest]
#[case::first(RelationId(0), ParticipantPosition(0), NodeId(6), [NodeId(6), NodeId(1)])]
#[case::shared_second(RelationId(0), ParticipantPosition(1), NodeId(7), [NodeId(0), NodeId(7)])]
#[case::empty_second(RelationId(1), ParticipantPosition(0), NodeId(4), [NodeId(4), NodeId(3)])]
#[case::sparse_last(RelationId(2), ParticipantPosition(1), NodeId(u32::MAX), [NodeId(4), NodeId(u32::MAX)])]
fn test_fixed_var_birelation_set_replace_participant_1(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
    #[case] participant: NodeId,
    #[case] expected: [NodeId; 2],
) {
    let mut entries = fixed_var_birelation_set_mutation_entries;
    let mut relations = FixedVarBirelationSet::new(entries.clone());
    relations.replace_participant_1(id, position, participant);
    entries[id.index()].0 = expected;
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::first(ParticipantPosition(0))]
#[case::last(ParticipantPosition(1))]
fn test_fixed_var_birelation_set_replace_participant_1_identity(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] position: ParticipantPosition,
) {
    let entries = fixed_var_birelation_set_mutation_entries;
    let input = FixedVarBirelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participant_1(RelationId(0), position, entries[0].0[position.index()]);
    assert_eq!(relations, input);
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_var_birelation_set_replace_participant_1_edges() {
    let mut relations =
        FixedVarBirelationSet::new(vec![([EdgeId(1)], vec![NodeId(1), NodeId(2)], 7)]);
    relations.replace_participant_1(RelationId(0), ParticipantPosition(0), EdgeId(5));
    let first = [EdgeId(5)];
    let second = vec![NodeId(1), NodeId(2)];
    for key in 0..7 {
        let nodes = if second.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if first.contains(&EdgeId(key)) {
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
#[case::past_end(RelationId(3), ParticipantPosition(0))]
#[case::maximum_id(RelationId(u32::MAX), ParticipantPosition(0))]
#[case::past_position(RelationId(0), ParticipantPosition(2))]
#[case::maximum_position(RelationId(0), ParticipantPosition(u32::MAX))]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participant_1_error(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
) {
    FixedVarBirelationSet::new(fixed_var_birelation_set_mutation_entries).replace_participant_1(
        id,
        position,
        NodeId(6),
    );
}

#[rstest]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participant_1_empty() {
    FixedVarBirelationSet::<NodeId, 2, NodeId, NonCloneData>::default().replace_participant_1(
        RelationId(0),
        ParticipantPosition(0),
        NodeId(6),
    );
}

#[rstest]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participant_1_zero_arity() {
    let mut relations = FixedVarBirelationSet::new(vec![([], vec![NodeId(1)], 7)]);
    relations.replace_participant_1(RelationId(0), ParticipantPosition(0), NodeId(2));
}

#[rstest]
#[case::shared_first(RelationId(0), ParticipantPosition(0), NodeId(6), vec![NodeId(6), NodeId(2), NodeId(1)])]
#[case::removed_reference(RelationId(0), ParticipantPosition(1), NodeId(7), vec![NodeId(1), NodeId(7), NodeId(1)])]
#[case::last_position(RelationId(0), ParticipantPosition(2), NodeId(7), vec![NodeId(1), NodeId(2), NodeId(7)])]
#[case::sparse_last_row(RelationId(2), ParticipantPosition(0), NodeId(u32::MAX), vec![NodeId(u32::MAX)])]
fn test_fixed_var_birelation_set_replace_participant_2(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
    #[case] participant: NodeId,
    #[case] expected: Vec<NodeId>,
) {
    let mut entries = fixed_var_birelation_set_mutation_entries;
    let mut relations = FixedVarBirelationSet::new(entries.clone());
    relations.replace_participant_2(id, position, participant);
    entries[id.index()].1 = expected;
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::first(ParticipantPosition(0))]
#[case::last(ParticipantPosition(2))]
fn test_fixed_var_birelation_set_replace_participant_2_identity(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] position: ParticipantPosition,
) {
    let entries = fixed_var_birelation_set_mutation_entries;
    let input = FixedVarBirelationSet::new(entries.clone());
    let mut relations = input.clone();
    relations.replace_participant_2(RelationId(0), position, entries[0].1[position.index()]);
    assert_eq!(relations, input);
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_var_birelation_set_replace_participant_2_edges() {
    let mut relations =
        FixedVarBirelationSet::new(vec![([EdgeId(1)], vec![NodeId(1), NodeId(2)], 7)]);
    relations.replace_participant_2(RelationId(0), ParticipantPosition(1), NodeId(6));
    let first = [EdgeId(1)];
    let second = vec![NodeId(1), NodeId(6)];
    for key in 0..7 {
        let nodes = if second.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if first.contains(&EdgeId(key)) {
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
fn test_fixed_var_birelation_set_replace_participant_2_zero_arity() {
    let mut relations = FixedVarBirelationSet::<NodeId, 0, NodeId, _>::new(vec![(
        [],
        vec![NodeId(1), NodeId(2)],
        BiPositionLabels {
            factor_1: vec![],
            factor_2: vec![7, 11],
        },
    )]);
    relations.replace_participant_2(RelationId(0), ParticipantPosition(1), NodeId(6));
    let entries = vec![(
        [],
        vec![NodeId(1), NodeId(6)],
        BiPositionLabels {
            factor_1: vec![],
            factor_2: vec![7, 11],
        },
    )];
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::past_end(RelationId(3), ParticipantPosition(0))]
#[case::maximum_id(RelationId(u32::MAX), ParticipantPosition(0))]
#[case::past_position(RelationId(0), ParticipantPosition(3))]
#[case::maximum_position(RelationId(0), ParticipantPosition(u32::MAX))]
#[case::empty_factor(RelationId(1), ParticipantPosition(0))]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participant_2_error(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
) {
    FixedVarBirelationSet::new(fixed_var_birelation_set_mutation_entries).replace_participant_2(
        id,
        position,
        NodeId(6),
    );
}

#[rstest]
#[should_panic]
fn test_fixed_var_birelation_set_replace_participant_2_empty() {
    FixedVarBirelationSet::<NodeId, 2, NodeId, NonCloneData>::default().replace_participant_2(
        RelationId(0),
        ParticipantPosition(0),
        NodeId(6),
    );
}

#[rstest]
#[case::prepend(RelationId(0), ParticipantPosition(0), NodeId(6), vec![NodeId(6), NodeId(1), NodeId(2), NodeId(1)])]
#[case::middle(RelationId(0), ParticipantPosition(1), NodeId(1), vec![NodeId(1), NodeId(1), NodeId(2), NodeId(1)])]
#[case::append(RelationId(0), ParticipantPosition(3), NodeId(7), vec![NodeId(1), NodeId(2), NodeId(1), NodeId(7)])]
#[case::empty_middle(RelationId(1), ParticipantPosition(0), NodeId(3), vec![NodeId(3)])]
#[case::sparse_last(RelationId(2), ParticipantPosition(1), NodeId(u32::MAX), vec![NodeId(5), NodeId(u32::MAX)])]
fn test_fixed_var_birelation_set_insert_participant_2(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
    #[case] participant: NodeId,
    #[case] expected: Vec<NodeId>,
) {
    let mut entries = fixed_var_birelation_set_mutation_entries;
    let mut relations = FixedVarBirelationSet::new(entries.clone());
    relations.insert_participant_2(id, position, participant);
    entries[id.index()].1 = expected;
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_var_birelation_set_insert_participant_2_edges() {
    let mut relations =
        FixedVarBirelationSet::new(vec![([EdgeId(1)], vec![NodeId(1), NodeId(2)], 7)]);
    relations.insert_participant_2(RelationId(0), ParticipantPosition(1), NodeId(6));
    let first = [EdgeId(1)];
    let second = vec![NodeId(1), NodeId(6), NodeId(2)];
    for key in 0..7 {
        let nodes = if second.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if first.contains(&EdgeId(key)) {
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
#[case::empty(vec![], ParticipantPosition(0), vec![NodeId(6)])]
#[case::middle(vec![NodeId(1), NodeId(2)], ParticipantPosition(1), vec![NodeId(1), NodeId(6), NodeId(2)])]
fn test_fixed_var_birelation_set_insert_participant_2_zero_arity(
    #[case] second: Vec<NodeId>,
    #[case] position: ParticipantPosition,
    #[case] expected: Vec<NodeId>,
) {
    let mut relations = FixedVarBirelationSet::<NodeId, 0, NodeId, _>::new(vec![(
        [],
        second,
        BiPositionLabels {
            factor_1: vec![],
            factor_2: vec![7, 11],
        },
    )]);
    relations.insert_participant_2(RelationId(0), position, NodeId(6));
    let entries = vec![(
        [],
        expected,
        BiPositionLabels {
            factor_1: vec![],
            factor_2: vec![7, 11],
        },
    )];
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::past_end(RelationId(3), ParticipantPosition(0))]
#[case::maximum_id(RelationId(u32::MAX), ParticipantPosition(0))]
#[case::past_position(RelationId(0), ParticipantPosition(4))]
#[case::maximum_position(RelationId(0), ParticipantPosition(u32::MAX))]
#[case::empty_factor(RelationId(1), ParticipantPosition(1))]
#[should_panic]
fn test_fixed_var_birelation_set_insert_participant_2_error(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
) {
    FixedVarBirelationSet::new(fixed_var_birelation_set_mutation_entries).insert_participant_2(
        id,
        position,
        NodeId(6),
    );
}

#[rstest]
#[should_panic]
fn test_fixed_var_birelation_set_insert_participant_2_empty() {
    FixedVarBirelationSet::<NodeId, 2, NodeId, NonCloneData>::default().insert_participant_2(
        RelationId(0),
        ParticipantPosition(0),
        NodeId(6),
    );
}

#[rstest]
#[case::shared_first(RelationId(0), ParticipantPosition(0), vec![NodeId(2), NodeId(1)])]
#[case::removed_reference(RelationId(0), ParticipantPosition(1), vec![NodeId(1), NodeId(1)])]
#[case::last_position(RelationId(0), ParticipantPosition(2), vec![NodeId(1), NodeId(2)])]
#[case::sole_participant(RelationId(2), ParticipantPosition(0), vec![])]
fn test_fixed_var_birelation_set_remove_participant_2(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
    #[case] expected: Vec<NodeId>,
) {
    let mut entries = fixed_var_birelation_set_mutation_entries;
    let mut relations = FixedVarBirelationSet::new(entries.clone());
    relations.remove_participant_2(id, position);
    entries[id.index()].1 = expected;
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
fn test_fixed_var_birelation_set_remove_participant_2_edges() {
    let mut relations =
        FixedVarBirelationSet::new(vec![([EdgeId(1)], vec![NodeId(1), NodeId(2)], 7)]);
    relations.remove_participant_2(RelationId(0), ParticipantPosition(0));
    let first = [EdgeId(1)];
    let second = vec![NodeId(2)];
    for key in 0..7 {
        let nodes = if second.contains(&NodeId(key)) {
            vec![RelationId(0)]
        } else {
            vec![]
        };
        let edges = if first.contains(&EdgeId(key)) {
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
fn test_fixed_var_birelation_set_remove_participant_2_zero_arity() {
    let mut relations = FixedVarBirelationSet::<NodeId, 0, NodeId, _>::new(vec![(
        [],
        vec![NodeId(1), NodeId(2)],
        BiPositionLabels {
            factor_1: vec![],
            factor_2: vec![7, 11],
        },
    )]);
    relations.remove_participant_2(RelationId(0), ParticipantPosition(0));
    let entries = vec![(
        [],
        vec![NodeId(2)],
        BiPositionLabels {
            factor_1: vec![],
            factor_2: vec![7, 11],
        },
    )];
    assert_fixed_var_birelation_rows(&relations, &entries);
}

#[rstest]
#[case::past_end(RelationId(3), ParticipantPosition(0))]
#[case::maximum_id(RelationId(u32::MAX), ParticipantPosition(0))]
#[case::past_position(RelationId(0), ParticipantPosition(3))]
#[case::maximum_position(RelationId(0), ParticipantPosition(u32::MAX))]
#[case::empty_factor(RelationId(1), ParticipantPosition(0))]
#[should_panic]
fn test_fixed_var_birelation_set_remove_participant_2_error(
    fixed_var_birelation_set_mutation_entries: Vec<([NodeId; 2], Vec<NodeId>, BiPositionLabels)>,
    #[case] id: RelationId,
    #[case] position: ParticipantPosition,
) {
    FixedVarBirelationSet::new(fixed_var_birelation_set_mutation_entries)
        .remove_participant_2(id, position);
}

#[rstest]
#[should_panic]
fn test_fixed_var_birelation_set_remove_participant_2_empty() {
    FixedVarBirelationSet::<NodeId, 2, NodeId, NonCloneData>::default()
        .remove_participant_2(RelationId(0), ParticipantPosition(0));
}

#[rstest]
#[case::rows(FixedVarBirelationSet::new(vec![([EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![7, 11]), ([EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![13, 17])]),
    FixedVarBirelationSet::new(vec![([EdgeId(3), EdgeId(6)], vec![NodeId(1), NodeId(5)], vec![7, 11]), ([EdgeId(3), EdgeId(6)], vec![NodeId(1), NodeId(5)], vec![13, 17])]))]
fn test_fixed_var_birelation_set_map(
    participant_correspondence: GraphCorrespondence,
    #[case] input: FixedVarBirelationSet<EdgeId, 2, NodeId, Vec<u32>>,
    #[case] expected: FixedVarBirelationSet<EdgeId, 2, NodeId, Vec<u32>>,
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
#[case::empty(FixedVarBirelationSet::new(vec![]))]
#[case::rows(FixedVarBirelationSet::new(vec![([EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![7, 11]), ([EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![13, 17])]))]
fn test_fixed_var_birelation_set_map_identity(
    #[case] input: FixedVarBirelationSet<EdgeId, 2, NodeId, Vec<u32>>,
) {
    let identity = GraphCorrespondence::new(
        Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2), NodeId(3)], 4),
        Correspondence::from_images(&[EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)], 4),
    );
    assert_eq!(input.try_map(&identity), Some(input.clone()));
    assert_eq!(input.map(&identity), input);
}

#[rstest]
#[case::forward(vec![NodeId(0), NodeId(1)], vec![NodeId(1), NodeId(2)])]
#[case::reversed(vec![NodeId(1), NodeId(0)], vec![NodeId(2), NodeId(1)])]
fn test_fixed_var_birelation_set_map_tracked_pushout(
    #[case] participants: Vec<NodeId>,
    #[case] expected_participants: Vec<NodeId>,
) {
    let left = Graph::new(2, &[[0, 1]]);
    let right = Graph::new(2, &[[0, 1]]);
    let overlap = GraphCorrespondence::new(
        Correspondence::new(vec![(NodeId(1), NodeId(0))], 2, 2).unwrap(),
        Correspondence::new(vec![], 1, 1).unwrap(),
    );
    let (_, pushout) = left.tracked_pushout(&right, &overlap);

    let relations: FixedVarBirelationSet<EdgeId, 1, NodeId, Vec<u32>> =
        FixedVarBirelationSet::new(vec![([EdgeId(0)], participants, vec![7, 11])]);
    let expected =
        FixedVarBirelationSet::new(vec![([EdgeId(1)], expected_participants, vec![7, 11])]);
    assert_eq!(relations.map(&pushout.right), expected);
}

#[rstest]
#[should_panic(expected = "correspondence must cover every participant reference")]
fn test_fixed_var_birelation_set_map_error(participant_correspondence: GraphCorrespondence) {
    let node = 1;
    let edge = 2;
    let input: FixedVarBirelationSet<EdgeId, 2, NodeId, Vec<u32>> =
        FixedVarBirelationSet::new(vec![(
            [EdgeId(edge), EdgeId(0)],
            vec![NodeId(node), NodeId(0)],
            vec![7, 11],
        )]);
    input.map(&participant_correspondence);
}

#[rstest]
#[case::missing_node(1, 2)]
#[case::outside_node(4, 2)]
#[case::missing_edge(2, 1)]
#[case::outside_edge(2, 4)]
fn test_fixed_var_birelation_set_try_map_error(
    participant_correspondence: GraphCorrespondence,
    #[case] node: u32,
    #[case] edge: u32,
) {
    let input: FixedVarBirelationSet<EdgeId, 2, NodeId, Vec<u32>> =
        FixedVarBirelationSet::new(vec![
            (
                [EdgeId(2), EdgeId(0)],
                vec![NodeId(2), NodeId(0)],
                vec![7, 11],
            ),
            (
                [EdgeId(edge), EdgeId(0)],
                vec![NodeId(node), NodeId(0)],
                vec![13, 17],
            ),
        ]);
    assert_eq!(input.try_map(&participant_correspondence), None);
}

#[rstest]
fn test_fixed_var_birelation_set_remap() {
    let rs: FixedVarBirelationSet<EdgeId, 2, NodeId, BiPositionLabels> =
        FixedVarBirelationSet::new(vec![(
            [EdgeId(0), EdgeId(1)],
            vec![NodeId(0), NodeId(1), NodeId(2)],
            BiPositionLabels {
                factor_1: vec![30, 31],
                factor_2: vec![40, 41, 42],
            },
        )]);
    let remapping = GraphRemapping::new(
        Remapping::new(vec![NodeId(2), NodeId(0), NodeId(1)]).expect("permutation images"),
        Remapping::new(vec![EdgeId(2), EdgeId(0), EdgeId(1)]).expect("permutation images"),
    );
    let out = rs.remap(&remapping);
    assert_eq!(out.participants_1(RelationId(0)), &[EdgeId(2), EdgeId(0)]);
    assert_eq!(
        out.participants_2(RelationId(0)),
        &[NodeId(2), NodeId(0), NodeId(1)]
    );
    assert_eq!(
        out.data(RelationId(0)),
        &BiPositionLabels {
            factor_1: vec![30, 31],
            factor_2: vec![40, 41, 42],
        }
    );
}

#[rstest]
#[case::covered(
    vec![NodeId(2), NodeId(0), NodeId(1)],
    vec![EdgeId(2), EdgeId(0), EdgeId(1)],
    true,
)]
#[case::uncovered_node(vec![NodeId(1), NodeId(0)], vec![EdgeId(1), EdgeId(0)], false)]
#[case::uncovered_edge(vec![NodeId(2), NodeId(0), NodeId(1)], vec![EdgeId(0)], false)]
fn test_fixed_var_birelation_set_try_remap(
    #[case] nodes: Vec<NodeId>,
    #[case] edges: Vec<EdgeId>,
    #[case] covered: bool,
) {
    let rs: FixedVarBirelationSet<EdgeId, 2, NodeId, BiPositionLabels> =
        FixedVarBirelationSet::new(vec![(
            [EdgeId(0), EdgeId(1)],
            vec![NodeId(0), NodeId(1), NodeId(2)],
            BiPositionLabels {
                factor_1: vec![30, 31],
                factor_2: vec![40, 41, 42],
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
#[case::empty(FixedVarBirelationSet::default())]
#[case::rows(
    FixedVarBirelationSet::new(vec![([NodeId(0)], vec![NodeId(2), NodeId(4)], "keep"), ([NodeId(5)], vec![NodeId(1), NodeId(3)], "drop")]),
)]
fn test_fixed_var_birelation_set_compact_identity(
    #[case] input: FixedVarBirelationSet<NodeId, 1, NodeId, &'static str>,
) {
    let compaction = GraphCompaction::new(Compaction::identity(6), Compaction::empty());
    assert_eq!(input.compact(&compaction), input);
    assert_eq!(
        input.tracked_compact(&compaction),
        (input.clone(), Compaction::identity(input.count())),
    );
}

#[rstest]
#[case::partial(
    vec![NodeId(1)],
    FixedVarBirelationSet::new(vec![([NodeId(0)], vec![NodeId(1), NodeId(3)], "keep")]),
    vec![RelationId(1)],
)]
#[case::all(
    vec![NodeId(0), NodeId(1)],
    FixedVarBirelationSet::default(),
    vec![RelationId(0), RelationId(1)],
)]
fn test_fixed_var_birelation_set_tracked_compact(
    fixed_var_birelation_set_compaction_input: FixedVarBirelationSet<
        NodeId,
        1,
        NodeId,
        &'static str,
    >,
    #[case] removed_nodes: Vec<NodeId>,
    #[case] expected: FixedVarBirelationSet<NodeId, 1, NodeId, &'static str>,
    #[case] removed_relations: Vec<RelationId>,
) {
    let input = fixed_var_birelation_set_compaction_input;
    let compaction = GraphCompaction::new(
        Compaction::new(6, removed_nodes).unwrap(),
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
fn test_fixed_var_birelation_set_tracked_pushout() {
    // coincidence requires *both* factors equal (site + members).
    let left = FixedVarBirelationSet::<NodeId, 1, NodeId, i32>::new(vec![(
        [NodeId(0)],
        vec![NodeId(1), NodeId(2)],
        10,
    )]);
    let right = FixedVarBirelationSet::<NodeId, 1, NodeId, i32>::new(vec![
        ([NodeId(0)], vec![NodeId(1), NodeId(2)], 5),
        ([NodeId(3)], vec![NodeId(4), NodeId(5)], 20),
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
        FixedVarBirelationSet::new(vec![
            ([NodeId(0)], vec![NodeId(1), NodeId(2)], 15),
            ([NodeId(3)], vec![NodeId(4), NodeId(5)], 20)
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
fn test_fixed_var_birelation_set_tracked_pullback(#[case] combined: Option<i32>) {
    let left: FixedVarBirelationSet<NodeId, 1, NodeId, i32> = FixedVarBirelationSet::new(vec![
        ([NodeId(0)], vec![NodeId(1), NodeId(2)], 10),
        ([NodeId(2)], vec![NodeId(3)], 20),
    ]);
    let right: FixedVarBirelationSet<NodeId, 1, NodeId, i32> = FixedVarBirelationSet::new(vec![
        ([NodeId(4)], vec![NodeId(5)], 30),
        ([NodeId(0)], vec![NodeId(1), NodeId(2)], 5),
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
            FixedVarBirelationSet::new(vec![([NodeId(0)], vec![NodeId(1), NodeId(2)], value)]),
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
fn test_fixed_var_birelation_set_default() {
    let rs = FixedVarBirelationSet::<EdgeId, 1, NodeId, ()>::default();
    assert_eq!(rs.count(), 0);
    assert!(!rs.has_incident_to_edge(EdgeId(0)));
}

#[rstest]
fn test_fixed_var_birelation_set_hash() {
    let entries = vec![
        ([EdgeId(2)], vec![NodeId(3), NodeId(1)], "first"),
        ([EdgeId(4)], vec![NodeId(5), NodeId(0)], "second"),
    ];
    let left: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> =
        FixedVarBirelationSet::new(entries.clone());
    let right: FixedVarBirelationSet<EdgeId, 1, NodeId, &str> = FixedVarBirelationSet::new(entries);
    assert_eq!(left, right);
    assert_eq!(hash(&left), hash(&right));
}
