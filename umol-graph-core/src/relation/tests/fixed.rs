use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};

use pretty_assertions::assert_eq;
use rstest::{fixture, rstest};

use crate::{
    Compaction, Correspondence, EdgeId, FixedRelationSet, GraphCompaction, GraphCorrespondence,
    GraphRemapping, NodeId, ParticipantPosition, RelationId, RelationPullbackCorrespondence,
    RelationPushoutCorrespondence, Remapping,
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

#[rstest]
fn test_fixed_relation_set_new() {
    let rs: FixedRelationSet<NodeId, &str, 2> = FixedRelationSet::new(vec![
        ([NodeId(0), NodeId(1)], "dative"),
        ([NodeId(1), NodeId(2)], "noncov"),
    ]);
    assert_eq!(rs.count(), 2);
    assert_eq!(rs.data(RelationId(0)), &"dative");
    assert_eq!(rs.participants(RelationId(0)), &[NodeId(0), NodeId(1)]);
    assert_eq!(rs.participants(RelationId(1)), &[NodeId(1), NodeId(2)]);
}

#[rstest]
#[case::empty(vec![], vec![], vec![])]
#[case::repeated(
    vec![
        ([NodeId(2), NodeId(0), NodeId(2)], "first"),
        ([NodeId(0), NodeId(2), NodeId(2)], "duplicate"),
        ([NodeId(3), NodeId(3), NodeId(3)], "other"),
    ],
    vec![RelationId(0), RelationId(1)],
    vec![RelationId(2)],
)]
fn test_fixed_relation_set_new_incidence(
    #[case] entries: Vec<([NodeId; 3], &str)>,
    #[case] at_two: Vec<RelationId>,
    #[case] at_three: Vec<RelationId>,
) {
    let relations = FixedRelationSet::<NodeId, &str, 3>::new(entries.clone());
    assert_eq!(relations.count(), entries.len());
    for (index, (participants, data)) in entries.iter().enumerate() {
        let id = RelationId(index as u32);
        assert_eq!(relations.participants(id), participants);
        assert_eq!(relations.data(id), data);
    }
    assert_eq!(relations.incident(NodeId(2)), at_two);
    assert_eq!(relations.incident(NodeId(3)), at_three);
    assert_eq!(relations.incident(NodeId(4)), &[]);
    assert_eq!(relations.incident_edge(EdgeId(2)), &[]);
    assert_eq!(relations.has_incident(NodeId(2)), !at_two.is_empty());
    assert!(!relations.has_incident_edge(EdgeId(2)));
    assert_eq!(relations.into_entries(), entries);
}

#[rstest]
fn test_fixed_relation_set_hash() {
    let entries = vec![
        ([NodeId(2), NodeId(0)], "first"),
        ([NodeId(3), NodeId(1)], "second"),
    ];
    let left: FixedRelationSet<NodeId, &str, 2> = FixedRelationSet::new(entries.clone());
    let right: FixedRelationSet<NodeId, &str, 2> = FixedRelationSet::new(entries);
    assert_eq!(left, right);
    assert_eq!(hash(&left), hash(&right));
}

#[rstest]
#[case::roundtrip(
    vec![([NodeId(2), NodeId(0)], "first"), ([NodeId(3), NodeId(1)], "second")],
)]
fn test_fixed_relation_set_into_entries(#[case] entries: Vec<([NodeId; 2], &str)>) {
    let rs: FixedRelationSet<NodeId, &str, 2> = FixedRelationSet::new(entries.clone());
    assert_eq!(rs.into_entries(), entries);
}

#[rstest]
fn test_fixed_relation_set_data_mut() {
    let mut rs: FixedRelationSet<NodeId, i32, 2> =
        FixedRelationSet::new(vec![([NodeId(0), NodeId(1)], 1)]);
    *rs.data_mut(RelationId(0)) = 99;
    assert_eq!(rs.data(RelationId(0)), &99);
}

#[rstest]
fn test_fixed_relation_set_iter() {
    let empty = FixedRelationSet::<NodeId, i32, 2>::default();
    assert_eq!(empty.iter().collect::<Vec<_>>(), vec![]);

    let rs: FixedRelationSet<NodeId, i32, 2> = FixedRelationSet::new(vec![
        ([NodeId(0), NodeId(1)], 1),
        ([NodeId(1), NodeId(2)], 2),
    ]);
    assert_eq!(rs.iter().len(), 2);
    assert_eq!(
        rs.iter().collect::<Vec<_>>(),
        vec![
            (RelationId(0), &[NodeId(0), NodeId(1)], &1),
            (RelationId(1), &[NodeId(1), NodeId(2)], &2),
        ],
    );
}

#[rstest]
fn test_fixed_relation_set_iter_mut() {
    let mut empty = FixedRelationSet::<NodeId, i32, 2>::default();
    assert_eq!(empty.iter_mut().len(), 0);

    let mut rs: FixedRelationSet<NodeId, i32, 2> = FixedRelationSet::new(vec![
        ([NodeId(0), NodeId(1)], 1),
        ([NodeId(1), NodeId(2)], 2),
        ([NodeId(2), NodeId(3)], 3),
    ]);
    assert_eq!(rs.iter_mut().len(), 3);
    for (id, participants, data) in rs.iter_mut() {
        assert_eq!(participants[0], NodeId(id.index() as u32));
        *data *= 10;
    }
    assert_eq!(rs.data(RelationId(0)), &10);
    assert_eq!(rs.data(RelationId(1)), &20);
    assert_eq!(rs.data(RelationId(2)), &30);
}

#[rstest]
fn test_fixed_relation_set_participants_ordered() {
    let rs: FixedRelationSet<NodeId, &str, 2> = FixedRelationSet::new(vec![
        ([NodeId(2), NodeId(0)], "a"),
        ([NodeId(3), NodeId(1)], "b"),
    ]);
    assert_eq!(rs.participants(RelationId(0)), &[NodeId(2), NodeId(0)]);
    assert_eq!(rs.participants(RelationId(1)), &[NodeId(3), NodeId(1)]);
}

#[rstest]
#[case::rotation(vec![ParticipantPosition(2), ParticipantPosition(0), ParticipantPosition(1)], [NodeId(2), NodeId(0), NodeId(1)])]
#[case::transposition(vec![ParticipantPosition(1), ParticipantPosition(0), ParticipantPosition(2)], [NodeId(1), NodeId(0), NodeId(2)])]
#[case::reversal(vec![ParticipantPosition(2), ParticipantPosition(1), ParticipantPosition(0)], [NodeId(2), NodeId(1), NodeId(0)])]
fn test_fixed_relation_set_permute_with(
    #[case] order: Vec<ParticipantPosition>,
    #[case] expected: [NodeId; 3],
) {
    let mut rs: FixedRelationSet<NodeId, &str, 3> = FixedRelationSet::new(vec![
        ([NodeId(0), NodeId(1), NodeId(2)], "a"),
        ([NodeId(3), NodeId(4), NodeId(5)], "b"),
    ]);
    let incidence_before: Vec<Vec<RelationId>> =
        (0..6).map(|i| rs.incident(NodeId(i)).to_vec()).collect();

    rs.permute_with(RelationId(0), &order);

    assert_eq!(rs.participants(RelationId(0)), &expected);
    assert_eq!(
        rs.participants(RelationId(1)),
        &[NodeId(3), NodeId(4), NodeId(5)]
    );
    assert_eq!(rs.data(RelationId(0)), &"a");
    assert_eq!(rs.data(RelationId(1)), &"b");
    let incidence_after: Vec<Vec<RelationId>> =
        (0..6).map(|i| rs.incident(NodeId(i)).to_vec()).collect();
    assert_eq!(incidence_after, incidence_before);
}

#[rstest]
#[case::ordered_factor(FixedRelationSet::<NodeId, &str, 3>::new(vec![
    ([NodeId(2), NodeId(0), NodeId(1)], "a"),
]))]
#[case::unordered_factor(FixedRelationSet::<NodeId, &str, 3>::new(vec![
    ([NodeId(2), NodeId(0), NodeId(1)], "a"),
]))]
fn test_fixed_relation_set_permute_with_identity(
    #[case] input: FixedRelationSet<NodeId, &'static str, 3>,
) {
    let mut permuted = input.clone();
    permuted.permute_with(
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
#[case::order_too_long(vec![ParticipantPosition(0), ParticipantPosition(1), ParticipantPosition(2), ParticipantPosition(0)])]
#[case::position_out_of_range(vec![ParticipantPosition(0), ParticipantPosition(1), ParticipantPosition(3)])]
#[case::position_repeated(vec![ParticipantPosition(0), ParticipantPosition(0), ParticipantPosition(1)])]
#[should_panic(expected = "permute")]
fn test_fixed_relation_set_permute_with_error(#[case] order: Vec<ParticipantPosition>) {
    let mut rs: FixedRelationSet<NodeId, &str, 3> =
        FixedRelationSet::new(vec![([NodeId(0), NodeId(1), NodeId(2)], "a")]);
    rs.permute_with(RelationId(0), &order);
}

#[rstest]
fn test_fixed_relation_set_incidence() {
    let rs: FixedRelationSet<NodeId, (), 2> = FixedRelationSet::new(vec![
        ([NodeId(0), NodeId(1)], ()),
        ([NodeId(0), NodeId(2)], ()),
        ([NodeId(2), NodeId(3)], ()),
    ]);
    assert_eq!(rs.incident(NodeId(0)), &[RelationId(0), RelationId(1)]);
    assert_eq!(rs.incident(NodeId(1)), &[RelationId(0)]);
    assert_eq!(rs.incident(NodeId(2)), &[RelationId(1), RelationId(2)]);
    assert_eq!(rs.incident(NodeId(3)), &[RelationId(2)]);
    assert!(rs.has_incident(NodeId(0)));
    assert!(!rs.has_incident(NodeId(5)));
}

#[rstest]
#[case::first(RelationId(0), true)]
#[case::last(RelationId(1), true)]
#[case::out_of_range(RelationId(2), false)]
fn test_fixed_relation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
    let rs: FixedRelationSet<NodeId, (), 2> = FixedRelationSet::new(vec![
        ([NodeId(0), NodeId(1)], ()),
        ([NodeId(1), NodeId(2)], ()),
    ]);
    assert_eq!(rs.contains(id), expected);
}

#[rstest]
fn test_fixed_relation_set_relation_ids() {
    assert_exact_size(FixedRelationSet::<NodeId, (), 2>::default().ids(), vec![]);
    let rs: FixedRelationSet<NodeId, (), 2> = FixedRelationSet::new(vec![
        ([NodeId(0), NodeId(1)], ()),
        ([NodeId(1), NodeId(2)], ()),
    ]);
    assert_exact_size(rs.ids(), vec![RelationId(0), RelationId(1)]);
}

#[fixture]
fn fixed_relation_set_compaction_input() -> FixedRelationSet<NodeId, &'static str, 2> {
    FixedRelationSet::new(vec![
        ([NodeId(0), NodeId(2)], "keep"),
        ([NodeId(1), NodeId(3)], "drop"),
    ])
}

#[rstest]
#[case::partial(
    vec![NodeId(1)],
    FixedRelationSet::new(vec![([NodeId(0), NodeId(1)], "keep")]),
    vec![RelationId(1)],
)]
#[case::all(
    vec![NodeId(0), NodeId(1)],
    FixedRelationSet::default(),
    vec![RelationId(0), RelationId(1)],
)]
fn test_fixed_relation_set_tracked_compact(
    fixed_relation_set_compaction_input: FixedRelationSet<NodeId, &'static str, 2>,
    #[case] removed_nodes: Vec<NodeId>,
    #[case] expected: FixedRelationSet<NodeId, &'static str, 2>,
    #[case] removed_relations: Vec<RelationId>,
) {
    let input = fixed_relation_set_compaction_input;
    let compaction = GraphCompaction::new(
        Compaction::new(4, removed_nodes).unwrap(),
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
#[case::empty(FixedRelationSet::default())]
#[case::rows(
    FixedRelationSet::new(vec![([NodeId(0), NodeId(2)], "keep"), ([NodeId(1), NodeId(3)], "drop")]),
)]
fn test_fixed_relation_set_compact_identity(
    #[case] input: FixedRelationSet<NodeId, &'static str, 2>,
) {
    let compaction = GraphCompaction::new(Compaction::identity(4), Compaction::empty());
    assert_eq!(input.compact(&compaction), input);
    assert_eq!(
        input.tracked_compact(&compaction),
        (input.clone(), Compaction::identity(input.count())),
    );
}

#[rstest]
#[case::rows(FixedRelationSet::new(vec![([NodeId(2), NodeId(0)], vec![7, 11]), ([NodeId(2), NodeId(0)], vec![13, 17])]),
    FixedRelationSet::new(vec![([NodeId(1), NodeId(5)], vec![7, 11]), ([NodeId(1), NodeId(5)], vec![13, 17])]))]
fn test_fixed_relation_set_map(
    participant_correspondence: GraphCorrespondence,
    #[case] input: FixedRelationSet<NodeId, Vec<u32>, 2>,
    #[case] expected: FixedRelationSet<NodeId, Vec<u32>, 2>,
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
        expected.incident(NodeId(1)),
        &[RelationId(0), RelationId(1)]
    );
}

#[rstest]
#[case::empty(FixedRelationSet::new(vec![]))]
#[case::rows(FixedRelationSet::new(vec![([NodeId(2), NodeId(0)], vec![7, 11]), ([NodeId(2), NodeId(0)], vec![13, 17])]))]
fn test_fixed_relation_set_map_identity(#[case] input: FixedRelationSet<NodeId, Vec<u32>, 2>) {
    let identity = GraphCorrespondence::new(
        Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2), NodeId(3)], 4),
        Correspondence::from_images(&[EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)], 4),
    );
    assert_eq!(input.try_map(&identity), Some(input.clone()));
    assert_eq!(input.map(&identity), input);
}

#[rstest]
#[case::missing_node(1)]
#[case::outside_node(4)]
fn test_fixed_relation_set_try_map_error(
    participant_correspondence: GraphCorrespondence,
    #[case] node: u32,
) {
    let input: FixedRelationSet<NodeId, Vec<u32>, 2> = FixedRelationSet::new(vec![
        ([NodeId(2), NodeId(0)], vec![7, 11]),
        ([NodeId(node), NodeId(0)], vec![13, 17]),
    ]);
    assert_eq!(input.try_map(&participant_correspondence), None);
}

#[rstest]
#[should_panic(expected = "correspondence must cover every participant reference")]
fn test_fixed_relation_set_map_error(participant_correspondence: GraphCorrespondence) {
    let node = 1;

    let input: FixedRelationSet<NodeId, Vec<u32>, 2> =
        FixedRelationSet::new(vec![([NodeId(node), NodeId(0)], vec![7, 11])]);
    input.map(&participant_correspondence);
}

#[rstest]
fn test_fixed_relation_set_remap() {
    let rs: FixedRelationSet<NodeId, PositionLabels, 2> =
        FixedRelationSet::new(vec![([NodeId(0), NodeId(1)], PositionLabels(vec![10, 11]))]);
    let remapping = GraphRemapping::new(
        Remapping::new(vec![NodeId(1), NodeId(0)]).expect("permutation images"),
        Remapping::empty(),
    );
    let out = rs.remap(&remapping);
    assert_eq!(out.participants(RelationId(0)), &[NodeId(1), NodeId(0)]);
    assert_eq!(out.data(RelationId(0)), &PositionLabels(vec![10, 11]));
}

#[rstest]
#[case::covered(vec![NodeId(1), NodeId(0)], true)]
#[case::uncovered_node(vec![NodeId(0)], false)]
fn test_fixed_relation_set_try_remap(#[case] nodes: Vec<NodeId>, #[case] covered: bool) {
    let rs: FixedRelationSet<NodeId, PositionLabels, 2> =
        FixedRelationSet::new(vec![([NodeId(0), NodeId(1)], PositionLabels(vec![10, 11]))]);
    let remapping = GraphRemapping::new(
        Remapping::new(nodes).expect("permutation images"),
        Remapping::empty(),
    );
    let expected = covered.then(|| rs.remap(&remapping));
    assert_eq!(rs.try_remap(&remapping), expected);
}

#[rstest]
fn test_fixed_relation_set_default() {
    let rs = FixedRelationSet::<NodeId, (), 2>::default();
    assert_eq!(rs.count(), 0);
    assert!(!rs.has_incident(NodeId(0)));
}

#[rstest]
#[case::exact(vec![NodeId(0), NodeId(1)], Some(RelationId(0)))]
#[case::reordered(vec![NodeId(1), NodeId(0)], Some(RelationId(0)))]
#[case::second(vec![NodeId(2), NodeId(3)], Some(RelationId(1)))]
#[case::absent(vec![NodeId(0), NodeId(3)], None)]
#[case::wrong_arity(vec![NodeId(0)], None)]
fn test_fixed_relation_set_coincident(
    #[case] query: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
) {
    let rs: FixedRelationSet<NodeId, (), 2> = FixedRelationSet::new(vec![
        ([NodeId(0), NodeId(1)], ()),
        ([NodeId(2), NodeId(3)], ()),
    ]);
    assert_eq!(
        query
            .first()
            .and_then(|&anchor| rs.coincident(anchor, &query)),
        expected,
    );
}

#[rstest]
#[case::reordered(NodeId(2), vec![NodeId(2), NodeId(2), NodeId(0)], Some(RelationId(0)), true)]
#[case::multiplicity(NodeId(2), vec![NodeId(0), NodeId(0), NodeId(2)], None, false)]
#[case::short(NodeId(2), vec![NodeId(2), NodeId(0)], None, false)]
#[case::long(NodeId(2), vec![NodeId(2), NodeId(0), NodeId(2), NodeId(2)], None, false)]
#[case::absent_anchor(NodeId(4), vec![NodeId(2), NodeId(2), NodeId(0)], None, true)]
fn test_fixed_relation_set_coincident_multiplicity(
    #[case] anchor: NodeId,
    #[case] query: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
    #[case] coincides: bool,
) {
    let relations = FixedRelationSet::<NodeId, &str, 3>::new(vec![
        ([NodeId(2), NodeId(0), NodeId(2)], "first"),
        ([NodeId(0), NodeId(2), NodeId(2)], "duplicate"),
        ([NodeId(3), NodeId(3), NodeId(3)], "other"),
    ]);
    assert_eq!(relations.coincident(anchor, &query), expected);
    assert_eq!(relations.is_coincident(RelationId(0), &query), coincides);
}

#[rstest]
fn test_fixed_relation_set_tracked_pushout() {
    // same-space glue: self {01}=10 {23}=20 ; right {01}=5 (coincides) {45}=30 (new); combine=sum.
    let left = FixedRelationSet::<NodeId, i32, 2>::new(vec![
        ([NodeId(0), NodeId(1)], 10),
        ([NodeId(2), NodeId(3)], 20),
    ]);
    let right = FixedRelationSet::<NodeId, i32, 2>::new(vec![
        ([NodeId(0), NodeId(1)], 5),
        ([NodeId(4), NodeId(5)], 30),
    ]);
    let (object, glue) = left
        .tracked_pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
            |(_, a), (_, b)| Some(a + b),
        )
        .expect("no ⊥");
    assert_eq!(
        left.pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
            |(_, a), (_, b)| Some(a + b),
        ),
        Some(object.clone()),
    );
    assert_eq!(
        object,
        FixedRelationSet::new(vec![
            ([NodeId(0), NodeId(1)], 15),
            ([NodeId(2), NodeId(3)], 20),
            ([NodeId(4), NodeId(5)], 30)
        ])
    );
    assert_eq!(
        glue,
        RelationPushoutCorrespondence {
            left: Correspondence::new(
                vec![
                    (RelationId(0), RelationId(0)),
                    (RelationId(1), RelationId(1))
                ],
                2,
                3,
            )
            .unwrap(),
            right: Correspondence::new(
                vec![
                    (RelationId(0), RelationId(0)),
                    (RelationId(1), RelationId(2))
                ],
                2,
                3,
            )
            .unwrap(),
        },
    );
}

#[rstest]
fn test_fixed_relation_set_tracked_pushout_error() {
    // combine returns ⊥ on the coincidence → the whole glue is inadmissible.
    let left = FixedRelationSet::<NodeId, i32, 2>::new(vec![([NodeId(0), NodeId(1)], 10)]);
    let right = FixedRelationSet::<NodeId, i32, 2>::new(vec![([NodeId(0), NodeId(1)], 5)]);
    assert_eq!(
        left.tracked_pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
            |_, _| None
        ),
        None
    );
    assert_eq!(
        left.pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
            |_, _| None
        ),
        None
    );
}

#[rstest]
#[case::combined(Some(15))]
#[case::incompatible(None)]
fn test_fixed_relation_set_tracked_pullback(#[case] combined: Option<i32>) {
    let left: FixedRelationSet<NodeId, i32, 2> = FixedRelationSet::new(vec![
        ([NodeId(0), NodeId(1)], 10),
        ([NodeId(2), NodeId(3)], 20),
    ]);
    let right: FixedRelationSet<NodeId, i32, 2> = FixedRelationSet::new(vec![
        ([NodeId(4), NodeId(5)], 30),
        ([NodeId(0), NodeId(1)], 5),
    ]);
    let result = left.tracked_pullback(
        &right,
        |set, parts: &[NodeId]| parts.first().and_then(|&id| set.coincident(id, parts)),
        |(_, a), (_, b)| combined.map(|_| a + b),
    );
    let plain = left.pullback(
        &right,
        |set, parts: &[NodeId]| parts.first().and_then(|&id| set.coincident(id, parts)),
        |(_, a), (_, b)| combined.map(|_| a + b),
    );
    let expected = combined.map(|value| {
        (
            FixedRelationSet::new(vec![([NodeId(0), NodeId(1)], value)]),
            RelationPullbackCorrespondence {
                left: Correspondence::new(vec![(RelationId(0), RelationId(0))], 1, 2).unwrap(),
                right: Correspondence::new(vec![(RelationId(0), RelationId(1))], 1, 2).unwrap(),
            },
        )
    });
    assert_eq!(plain, expected.as_ref().map(|(object, _)| object.clone()));
    assert_eq!(result, expected);
}
