use std::fmt::Debug;

use pretty_assertions::assert_eq;
use rstest::{fixture, rstest};

use crate::{
    Compaction, Correspondence, EdgeId, GraphCompaction, GraphCorrespondence, GraphRemapping,
    NodeId, ParticipantPosition, RelationId, RelationPullbackCorrespondence,
    RelationPushoutCorrespondence, Remapping, VarVarBirelationSet,
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
fn var_var_birelation_set_compaction_input() -> VarVarBirelationSet<NodeId, NodeId, &'static str> {
    VarVarBirelationSet::new(vec![
        (vec![NodeId(0), NodeId(2)], vec![NodeId(4)], "keep"),
        (vec![NodeId(5)], vec![NodeId(1)], "drop"),
    ])
}

#[rstest]
fn test_var_var_birelation_set_new() {
    let rs: VarVarBirelationSet<NodeId, EdgeId, &str> =
        VarVarBirelationSet::new(vec![(vec![NodeId(0), NodeId(1)], vec![EdgeId(5)], "y")]);
    assert_eq!(rs.count(), 1);
    assert_eq!(rs.participants_1(RelationId(0)), &[NodeId(0), NodeId(1)]);
    assert_eq!(rs.participants_2(RelationId(0)), &[EdgeId(5)]);
    assert_eq!(rs.data(RelationId(0)), &"y");
}

#[rstest]
#[case::empty(vec![], vec![], vec![])]
#[case::repeated(
    vec![
        (vec![NodeId(2), NodeId(0), NodeId(2)], vec![NodeId(2), NodeId(1), NodeId(2)], "first"),
        (
            vec![NodeId(0), NodeId(2), NodeId(2)],
            vec![NodeId(1), NodeId(2), NodeId(2)],
            "duplicate",
        ),
        (vec![NodeId(3), NodeId(3), NodeId(3)], vec![NodeId(3), NodeId(3), NodeId(3)], "other"),
    ],
    vec![RelationId(0), RelationId(1)],
    vec![RelationId(2)],
)]
#[case::empty_factors(vec![(vec![], vec![], "empty")], vec![], vec![])]
fn test_var_var_birelation_set_new_incidence(
    #[case] entries: Vec<(Vec<NodeId>, Vec<NodeId>, &str)>,
    #[case] at_two: Vec<RelationId>,
    #[case] at_three: Vec<RelationId>,
) {
    let relations = VarVarBirelationSet::<NodeId, NodeId, &str>::new(entries.clone());
    assert_eq!(relations.count(), entries.len());
    for (index, (first, second, data)) in entries.iter().enumerate() {
        let id = RelationId(index as u32);
        assert_eq!(relations.participants_1(id), first);
        assert_eq!(relations.participants_2(id), second);
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
#[case::roundtrip(
    vec![
        (vec![NodeId(2), NodeId(0)], vec![EdgeId(4), EdgeId(1)], "first"),
        (vec![NodeId(5), NodeId(3)], vec![EdgeId(6), EdgeId(2)], "second"),
    ],
)]
fn test_var_var_birelation_set_into_entries(
    #[case] entries: Vec<(Vec<NodeId>, Vec<EdgeId>, &str)>,
) {
    let rs: VarVarBirelationSet<NodeId, EdgeId, &str> = VarVarBirelationSet::new(entries.clone());
    assert_eq!(rs.into_entries(), entries);
}

#[rstest]
#[case::first(RelationId(0), true)]
#[case::out_of_range(RelationId(1), false)]
fn test_var_var_birelation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
    let rs: VarVarBirelationSet<NodeId, EdgeId, &str> =
        VarVarBirelationSet::new(vec![(vec![NodeId(0)], vec![EdgeId(1)], "y")]);
    assert_eq!(rs.contains(id), expected);
}

#[rstest]
fn test_var_var_birelation_set_relation_ids() {
    assert_exact_size(
        VarVarBirelationSet::<NodeId, EdgeId, &str>::default().ids(),
        vec![],
    );
    let rs: VarVarBirelationSet<NodeId, EdgeId, &str> = VarVarBirelationSet::new(vec![
        (vec![NodeId(0)], vec![EdgeId(1)], "a"),
        (vec![NodeId(2)], vec![EdgeId(3)], "b"),
    ]);
    assert_exact_size(rs.ids(), vec![RelationId(0), RelationId(1)]);
}

#[rstest]
fn test_var_var_birelation_set_iter() {
    let empty = VarVarBirelationSet::<NodeId, EdgeId, i32>::default();
    assert_eq!(empty.iter().collect::<Vec<_>>(), vec![]);

    let rs: VarVarBirelationSet<NodeId, EdgeId, i32> = VarVarBirelationSet::new(vec![
        (vec![NodeId(0), NodeId(4)], vec![EdgeId(1)], 1),
        (vec![NodeId(2)], vec![EdgeId(3)], 2),
    ]);
    assert_eq!(rs.iter().len(), 2);
    assert_eq!(
        rs.iter().collect::<Vec<_>>(),
        vec![
            (
                RelationId(0),
                [NodeId(0), NodeId(4)].as_slice(),
                [EdgeId(1)].as_slice(),
                &1,
            ),
            (
                RelationId(1),
                [NodeId(2)].as_slice(),
                [EdgeId(3)].as_slice(),
                &2
            ),
        ],
    );
}

#[rstest]
fn test_var_var_birelation_set_iter_mut() {
    let mut empty = VarVarBirelationSet::<NodeId, EdgeId, i32>::default();
    assert_eq!(empty.iter_mut().len(), 0);

    let mut rs: VarVarBirelationSet<NodeId, EdgeId, i32> = VarVarBirelationSet::new(vec![
        (vec![NodeId(0), NodeId(4)], vec![EdgeId(1)], 1),
        (vec![NodeId(2)], vec![EdgeId(3)], 2),
    ]);
    let arities: Vec<(usize, usize)> = rs
        .iter_mut()
        .map(|(_, first, second, data)| {
            *data *= 10;
            (first.len(), second.len())
        })
        .collect();
    assert_eq!(arities, vec![(2, 1), (1, 1)]);
    assert_eq!(rs.data(RelationId(0)), &10);
    assert_eq!(rs.data(RelationId(1)), &20);
}

#[rstest]
fn test_var_var_birelation_set_data_mut() {
    let mut rs: VarVarBirelationSet<NodeId, EdgeId, i32> =
        VarVarBirelationSet::new(vec![(vec![NodeId(0)], vec![EdgeId(1)], 1)]);
    *rs.data_mut(RelationId(0)) = 99;
    assert_eq!(rs.data(RelationId(0)), &99);
}

#[rstest]
fn test_var_var_birelation_set_incidence() {
    let rs: VarVarBirelationSet<NodeId, EdgeId, &str> =
        VarVarBirelationSet::new(vec![(vec![NodeId(0), NodeId(1)], vec![EdgeId(5)], "y")]);
    assert_eq!(rs.incident(NodeId(1)), &[RelationId(0)]);
    assert_eq!(rs.incident_edge(EdgeId(5)), &[RelationId(0)]);
    assert!(rs.has_incident(NodeId(0)));
    assert!(rs.has_incident_edge(EdgeId(5)));
    assert!(!rs.has_incident_edge(EdgeId(0)));
}

#[rstest]
#[case::exact(vec![NodeId(0), NodeId(1)], vec![NodeId(2), NodeId(3)], Some(RelationId(0)))]
#[case::reordered(vec![NodeId(1), NodeId(0)], vec![NodeId(3), NodeId(2)], Some(RelationId(0)))]
#[case::role_swap(vec![NodeId(2), NodeId(3)], vec![NodeId(0), NodeId(1)], None)]
#[case::absent(vec![NodeId(0), NodeId(1)], vec![NodeId(2), NodeId(9)], None)]
fn test_var_var_birelation_set_coincident(
    #[case] query_1: Vec<NodeId>,
    #[case] query_2: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
) {
    let rs: VarVarBirelationSet<NodeId, NodeId, ()> = VarVarBirelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], vec![NodeId(2), NodeId(3)], ()),
        (vec![NodeId(4)], vec![NodeId(5)], ()),
    ]);
    assert_eq!(
        query_1
            .first()
            .and_then(|&anchor| rs.coincident(anchor, &query_1, &query_2)),
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
fn test_var_var_birelation_set_coincident_multiplicity(
    #[case] anchor: NodeId,
    #[case] query: Vec<NodeId>,
    #[case] query_2: Vec<NodeId>,
    #[case] expected: Option<RelationId>,
    #[case] coincides: bool,
) {
    let relations = VarVarBirelationSet::<NodeId, NodeId, &str>::new(vec![
        (
            vec![NodeId(2), NodeId(0), NodeId(2)],
            vec![NodeId(2), NodeId(1), NodeId(2)],
            "first",
        ),
        (
            vec![NodeId(0), NodeId(2), NodeId(2)],
            vec![NodeId(1), NodeId(2), NodeId(2)],
            "duplicate",
        ),
        (
            vec![NodeId(3), NodeId(3), NodeId(3)],
            vec![NodeId(3), NodeId(3), NodeId(3)],
            "other",
        ),
    ]);
    assert_eq!(relations.coincident(anchor, &query, &query_2), expected);
    assert_eq!(
        relations.is_coincident(RelationId(0), &query, &query_2),
        coincides
    );
}

#[rstest]
fn test_var_var_birelation_set_permute_1_with() {
    let mut rs: VarVarBirelationSet<NodeId, EdgeId, &str> = VarVarBirelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], vec![EdgeId(9)], "a"),
        (
            vec![NodeId(2), NodeId(3), NodeId(4)],
            vec![EdgeId(7), EdgeId(8)],
            "b",
        ),
    ]);
    rs.permute_1_with(
        RelationId(1),
        &[
            ParticipantPosition(2),
            ParticipantPosition(0),
            ParticipantPosition(1),
        ],
    );
    assert_eq!(rs.participants_1(RelationId(0)), &[NodeId(0), NodeId(1)]);
    assert_eq!(
        rs.participants_1(RelationId(1)),
        &[NodeId(4), NodeId(2), NodeId(3)]
    );
    assert_eq!(rs.participants_2(RelationId(1)), &[EdgeId(7), EdgeId(8)]);
    assert_eq!(rs.data(RelationId(1)), &"b");
}

#[rstest]
fn test_var_var_birelation_set_permute_2_with() {
    let mut rs: VarVarBirelationSet<NodeId, EdgeId, &str> = VarVarBirelationSet::new(vec![
        (vec![NodeId(0), NodeId(1)], vec![EdgeId(9)], "a"),
        (
            vec![NodeId(2), NodeId(3), NodeId(4)],
            vec![EdgeId(7), EdgeId(8)],
            "b",
        ),
    ]);
    rs.permute_2_with(
        RelationId(1),
        &[ParticipantPosition(1), ParticipantPosition(0)],
    );
    assert_eq!(rs.participants_2(RelationId(0)), &[EdgeId(9)]);
    assert_eq!(rs.participants_2(RelationId(1)), &[EdgeId(8), EdgeId(7)]);
    assert_eq!(
        rs.participants_1(RelationId(1)),
        &[NodeId(2), NodeId(3), NodeId(4)]
    );
    assert_eq!(rs.data(RelationId(1)), &"b");
}

#[rstest]
#[case::rows(VarVarBirelationSet::new(vec![(vec![EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![7, 11]), (vec![EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![13, 17])]),
    VarVarBirelationSet::new(vec![(vec![EdgeId(3), EdgeId(6)], vec![NodeId(1), NodeId(5)], vec![7, 11]), (vec![EdgeId(3), EdgeId(6)], vec![NodeId(1), NodeId(5)], vec![13, 17])]))]
fn test_var_var_birelation_set_map(
    participant_correspondence: GraphCorrespondence,
    #[case] input: VarVarBirelationSet<EdgeId, NodeId, Vec<u32>>,
    #[case] expected: VarVarBirelationSet<EdgeId, NodeId, Vec<u32>>,
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
    assert_eq!(
        expected.incident_edge(EdgeId(3)),
        expected.incident(NodeId(1))
    );
}

#[rstest]
#[case::empty(VarVarBirelationSet::new(vec![]))]
#[case::rows(VarVarBirelationSet::new(vec![(vec![EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![7, 11]), (vec![EdgeId(2), EdgeId(0)], vec![NodeId(2), NodeId(0)], vec![13, 17])]))]
fn test_var_var_birelation_set_map_identity(
    #[case] input: VarVarBirelationSet<EdgeId, NodeId, Vec<u32>>,
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
fn test_var_var_birelation_set_map_error(participant_correspondence: GraphCorrespondence) {
    let node = 1;
    let edge = 2;
    let input: VarVarBirelationSet<EdgeId, NodeId, Vec<u32>> = VarVarBirelationSet::new(vec![(
        vec![EdgeId(edge), EdgeId(0)],
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
fn test_var_var_birelation_set_try_map_error(
    participant_correspondence: GraphCorrespondence,
    #[case] node: u32,
    #[case] edge: u32,
) {
    let input: VarVarBirelationSet<EdgeId, NodeId, Vec<u32>> = VarVarBirelationSet::new(vec![
        (
            vec![EdgeId(2), EdgeId(0)],
            vec![NodeId(2), NodeId(0)],
            vec![7, 11],
        ),
        (
            vec![EdgeId(edge), EdgeId(0)],
            vec![NodeId(node), NodeId(0)],
            vec![13, 17],
        ),
    ]);
    assert_eq!(input.try_map(&participant_correspondence), None);
}

#[rstest]
fn test_var_var_birelation_set_remap() {
    let rs: VarVarBirelationSet<NodeId, EdgeId, BiPositionLabels> =
        VarVarBirelationSet::new(vec![(
            vec![NodeId(0), NodeId(1)],
            vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            BiPositionLabels {
                factor_1: vec![50, 51],
                factor_2: vec![60, 61, 62],
            },
        )]);
    let remapping = GraphRemapping::new(
        Remapping::new(vec![NodeId(1), NodeId(0)]).expect("permutation images"),
        Remapping::new(vec![EdgeId(2), EdgeId(0), EdgeId(1)]).expect("permutation images"),
    );
    let out = rs.remap(&remapping);
    assert_eq!(out.participants_1(RelationId(0)), &[NodeId(1), NodeId(0)]);
    assert_eq!(
        out.participants_2(RelationId(0)),
        &[EdgeId(2), EdgeId(0), EdgeId(1)]
    );
    assert_eq!(
        out.data(RelationId(0)),
        &BiPositionLabels {
            factor_1: vec![50, 51],
            factor_2: vec![60, 61, 62],
        }
    );
}

#[rstest]
#[case::covered(
    vec![NodeId(1), NodeId(0)],
    vec![EdgeId(2), EdgeId(0), EdgeId(1)],
    true,
)]
#[case::uncovered_node(vec![NodeId(0)], vec![EdgeId(2), EdgeId(0), EdgeId(1)], false)]
#[case::uncovered_edge(vec![NodeId(1), NodeId(0)], vec![EdgeId(1), EdgeId(0)], false)]
fn test_var_var_birelation_set_try_remap(
    #[case] nodes: Vec<NodeId>,
    #[case] edges: Vec<EdgeId>,
    #[case] covered: bool,
) {
    let rs: VarVarBirelationSet<NodeId, EdgeId, BiPositionLabels> =
        VarVarBirelationSet::new(vec![(
            vec![NodeId(0), NodeId(1)],
            vec![EdgeId(0), EdgeId(1), EdgeId(2)],
            BiPositionLabels {
                factor_1: vec![50, 51],
                factor_2: vec![60, 61, 62],
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
#[case::empty(VarVarBirelationSet::default())]
#[case::rows(
    VarVarBirelationSet::new(vec![(vec![NodeId(0), NodeId(2)], vec![NodeId(4)], "keep"), (vec![NodeId(5)], vec![NodeId(1)], "drop")]),
)]
fn test_var_var_birelation_set_compact_identity(
    #[case] input: VarVarBirelationSet<NodeId, NodeId, &'static str>,
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
    VarVarBirelationSet::new(vec![(vec![NodeId(0), NodeId(1)], vec![NodeId(3)], "keep")]),
    vec![RelationId(1)],
)]
#[case::all(
    vec![NodeId(0), NodeId(1)],
    VarVarBirelationSet::default(),
    vec![RelationId(0), RelationId(1)],
)]
fn test_var_var_birelation_set_tracked_compact(
    var_var_birelation_set_compaction_input: VarVarBirelationSet<NodeId, NodeId, &'static str>,
    #[case] removed_nodes: Vec<NodeId>,
    #[case] expected: VarVarBirelationSet<NodeId, NodeId, &'static str>,
    #[case] removed_relations: Vec<RelationId>,
) {
    let input = var_var_birelation_set_compaction_input;
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
fn test_var_var_birelation_set_tracked_pushout() {
    let left = VarVarBirelationSet::<NodeId, NodeId, i32>::new(vec![(
        vec![NodeId(0), NodeId(1)],
        vec![NodeId(2), NodeId(3)],
        10,
    )]);
    let right = VarVarBirelationSet::<NodeId, NodeId, i32>::new(vec![
        (vec![NodeId(0), NodeId(1)], vec![NodeId(2), NodeId(3)], 5),
        (vec![NodeId(4)], vec![NodeId(5)], 20),
    ]);
    let (object, glue) = left
        .tracked_pushout(
            &right,
            |set: &_, q1: &[NodeId], q2: &_| q1.first().and_then(|&n| set.coincident(n, q1, q2)),
            |(_, _, a), (_, _, b)| Some(a + b),
        )
        .expect("no ⊥");
    assert_eq!(
        left.pushout(
            &right,
            |set: &_, q1: &[NodeId], q2: &_| {
                q1.first().and_then(|&n| set.coincident(n, q1, q2))
            },
            |(_, _, a), (_, _, b)| Some(a + b),
        ),
        Some(object.clone()),
    );
    assert_eq!(
        object,
        VarVarBirelationSet::new(vec![
            (vec![NodeId(0), NodeId(1)], vec![NodeId(2), NodeId(3)], 15),
            (vec![NodeId(4)], vec![NodeId(5)], 20)
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
                q1.first().and_then(|&n| set.coincident(n, q1, q2))
            },
            |_, _| None,
        ),
        None
    );
    assert_eq!(
        left.pushout(
            &right,
            |set: &_, q1: &[NodeId], q2: &_| {
                q1.first().and_then(|&n| set.coincident(n, q1, q2))
            },
            |_, _| None,
        ),
        None
    );
}

#[rstest]
#[case::combined(Some(15))]
#[case::incompatible(None)]
fn test_var_var_birelation_set_tracked_pullback(#[case] combined: Option<i32>) {
    let left: VarVarBirelationSet<NodeId, NodeId, i32> = VarVarBirelationSet::new(vec![
        (vec![NodeId(0)], vec![NodeId(1)], 10),
        (vec![NodeId(2)], vec![NodeId(3)], 20),
    ]);
    let right: VarVarBirelationSet<NodeId, NodeId, i32> = VarVarBirelationSet::new(vec![
        (vec![NodeId(4)], vec![NodeId(5)], 30),
        (vec![NodeId(0)], vec![NodeId(1)], 5),
    ]);
    let result = left.tracked_pullback(
        &right,
        |set, first: &[NodeId], second: &[NodeId]| {
            first
                .first()
                .and_then(|&id| set.coincident(id, first, second))
        },
        |(_, _, a), (_, _, b)| combined.map(|_| a + b),
    );
    let plain = left.pullback(
        &right,
        |set, first: &[NodeId], second: &[NodeId]| {
            first
                .first()
                .and_then(|&id| set.coincident(id, first, second))
        },
        |(_, _, a), (_, _, b)| combined.map(|_| a + b),
    );
    let expected = combined.map(|value| {
        (
            VarVarBirelationSet::new(vec![(vec![NodeId(0)], vec![NodeId(1)], value)]),
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
fn test_var_var_birelation_set_default() {
    let rs = VarVarBirelationSet::<NodeId, EdgeId, ()>::default();
    assert_eq!(rs.count(), 0);
    assert!(!rs.has_incident(NodeId(0)));
}
