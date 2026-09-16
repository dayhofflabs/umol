use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};

use pretty_assertions::assert_eq;
use rstest::{fixture, rstest};

use crate::{
    Compaction, Correspondence, EdgeId, GraphCompaction, GraphCorrespondence, GraphRemapping,
    NodeId, ParticipantPosition, RelationId, RelationPullbackCorrespondence,
    RelationPushoutCorrespondence, Remapping, VarRelationSet,
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

/// What frame a pushout coincidence comes out in, and whether the two sides reaching `combine`
/// are in the same frame. Recorded so the answer is read off a run rather than argued.
///
/// Both sides carry the same two participants in opposite order with a position-labelled
/// payload. `combine` records what it was handed and returns the left payload untouched, so the
/// test observes the inputs rather than any realignment.
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
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
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
/// Two right entries coinciding with one left entry is rejected, not merged: the right
/// coprojection must be injective, and `Correspondence` asserts that.
///
/// This is why `pushout` may read the left payload out of its own output buffer — the entry can
/// never be merged into twice. `pullback` reads the same payload from the source instead, and
/// the two agree because the case that would separate them cannot be constructed.
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
        |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
        |(_, a), (_, b)| {
            Some(PositionLabels(
                a.0.iter().zip(&b.0).map(|(x, y)| x + y).collect(),
            ))
        },
    );
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
    assert_eq!(relations.incident(NodeId(2)), at_two);
    assert_eq!(relations.incident(NodeId(3)), at_three);
    assert_eq!(relations.incident(NodeId(4)), &[]);
    assert_eq!(relations.incident_edge(EdgeId(2)), &[]);
    assert_eq!(relations.has_incident(NodeId(2)), !at_two.is_empty());
    assert!(!relations.has_incident_edge(EdgeId(2)));
    assert_eq!(relations.into_entries(), entries);
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
fn test_var_relation_set_data_mut() {
    let mut rs: VarRelationSet<NodeId, i32> =
        VarRelationSet::new(vec![(vec![NodeId(0), NodeId(1), NodeId(2)], 1)]);
    *rs.data_mut(RelationId(0)) = 99;
    assert_eq!(rs.data(RelationId(0)), &99);
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
    assert_eq!(rs.incident(NodeId(0)), &[RelationId(0)]);
    assert_eq!(rs.incident(NodeId(4)), &[RelationId(1)]);
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
#[case::second_relation(RelationId(1), vec![ParticipantPosition(2), ParticipantPosition(0), ParticipantPosition(1)],
    vec![NodeId(4), NodeId(2), NodeId(3)])]
#[case::last_relation(RelationId(2), vec![ParticipantPosition(1), ParticipantPosition(0)], vec![NodeId(6), NodeId(5)])]
fn test_var_relation_set_permute_with(
    #[case] id: RelationId,
    #[case] order: Vec<ParticipantPosition>,
    #[case] expected: Vec<NodeId>,
) {
    let mut rs: VarRelationSet<NodeId, &str> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], "a"),
        (vec![NodeId(2), NodeId(3), NodeId(4)], "b"),
        (vec![NodeId(5), NodeId(6)], "c"),
    ]);
    let incidence_before: Vec<Vec<RelationId>> =
        (0..7).map(|i| rs.incident(NodeId(i)).to_vec()).collect();

    rs.permute_with(id, &order);

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
    let incidence_after: Vec<Vec<RelationId>> =
        (0..7).map(|i| rs.incident(NodeId(i)).to_vec()).collect();
    assert_eq!(incidence_after, incidence_before);
}

#[rstest]
fn test_var_relation_set_permute_with_identity() {
    let input: VarRelationSet<NodeId, &str> = VarRelationSet::new(vec![
        (vec![NodeId(2), NodeId(0), NodeId(1)], "a"),
        (vec![NodeId(4), NodeId(3)], "b"),
    ]);
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
#[case::position_out_of_range(vec![ParticipantPosition(0), ParticipantPosition(1), ParticipantPosition(9)])]
#[case::position_repeated(vec![ParticipantPosition(1), ParticipantPosition(1), ParticipantPosition(0)])]
#[should_panic(expected = "permute")]
fn test_var_relation_set_permute_with_error(#[case] order: Vec<ParticipantPosition>) {
    let mut rs: VarRelationSet<NodeId, &str> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], "a"),
        (vec![NodeId(2), NodeId(3), NodeId(4)], "b"),
    ]);
    rs.permute_with(RelationId(1), &order);
}

#[rstest]
fn test_var_relation_set_incidence() {
    let rs: VarRelationSet<NodeId, ()> = VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(1), NodeId(2)], ()),
        (vec![NodeId(2), NodeId(3), NodeId(4)], ()),
    ]);
    assert_eq!(rs.incident(NodeId(0)), &[RelationId(0)]);
    assert_eq!(rs.incident(NodeId(2)), &[RelationId(0), RelationId(1)]);
    assert_eq!(rs.incident(NodeId(4)), &[RelationId(1)]);
    assert!(rs.has_incident(NodeId(0)));
    assert!(!rs.has_incident(NodeId(7)));
}

#[rstest]
fn test_var_relation_set_edge_incidence() {
    let rs: VarRelationSet<EdgeId, &str> = VarRelationSet::new(vec![
        (vec![EdgeId(0), EdgeId(2)], "a"),
        (vec![EdgeId(1), EdgeId(2)], "b"),
    ]);
    assert_eq!(rs.incident_edge(EdgeId(2)), &[RelationId(0), RelationId(1)]);
    assert_eq!(rs.incident_edge(EdgeId(0)), &[RelationId(0)]);
    assert!(rs.has_incident_edge(EdgeId(2)));
    assert!(!rs.has_incident_edge(EdgeId(5)));
    assert!(rs.incident(NodeId(0)).is_empty());
    assert!(!rs.has_incident(NodeId(0)));
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

#[fixture]
fn var_relation_set_compaction_input() -> VarRelationSet<NodeId, &'static str> {
    VarRelationSet::new(vec![
        (vec![NodeId(0), NodeId(2), NodeId(4)], "keep"),
        (vec![NodeId(1), NodeId(3)], "drop"),
    ])
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
        expected.incident(NodeId(1)),
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
#[should_panic(expected = "correspondence must cover every participant reference")]
fn test_var_relation_set_map_error(participant_correspondence: GraphCorrespondence) {
    let node = 1;

    let input: VarRelationSet<NodeId, Vec<u32>> =
        VarRelationSet::new(vec![(vec![NodeId(node), NodeId(0)], vec![7, 11])]);
    input.map(&participant_correspondence);
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
fn test_var_relation_set_default() {
    let rs = VarRelationSet::<NodeId, ()>::default();
    assert_eq!(rs.count(), 0);
    assert!(!rs.has_incident(NodeId(0)));
}

#[rstest]
#[case::exact(vec![NodeId(0), NodeId(1), NodeId(2)], Some(RelationId(0)))]
#[case::reordered(vec![NodeId(2), NodeId(0), NodeId(1)], Some(RelationId(0)))]
#[case::second(vec![NodeId(3), NodeId(4)], Some(RelationId(1)))]
#[case::subset(vec![NodeId(0), NodeId(1)], None)]
#[case::superset(vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)], None)]
fn test_var_relation_set_coincident(
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
fn test_var_relation_set_coincident_multiplicity(
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
    assert_eq!(relations.coincident(anchor, &query), expected);
    assert_eq!(relations.is_coincident(RelationId(0), &query), coincides);
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
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
            |_, _| None,
        ),
        None
    );
    assert_eq!(
        left.pushout(
            &right,
            |set: &_, q: &[NodeId]| q.first().and_then(|&n| set.coincident(n, q)),
            |_, _| None,
        ),
        None
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
