use pretty_assertions::assert_eq;
use rstest::{fixture, rstest};

use crate::{
    Compaction, Correspondence, EdgeId, GraphCompaction, GraphCorrespondence, NodeId,
    ParticipantRefs, RelationParticipant,
};

#[fixture]
fn participant_correspondence() -> GraphCorrespondence {
    GraphCorrespondence::new(
        Correspondence::new(vec![(NodeId(0), NodeId(5)), (NodeId(2), NodeId(1))], 4, 6).unwrap(),
        Correspondence::new(vec![(EdgeId(0), EdgeId(6)), (EdgeId(2), EdgeId(3))], 4, 7).unwrap(),
    )
}

#[rstest]
#[case::mapped(NodeId(2), Some(NodeId(1)))]
#[case::unmatched(NodeId(1), None)]
#[case::outside(NodeId(4), None)]
fn test_node_id_try_map(
    participant_correspondence: GraphCorrespondence,
    #[case] id: NodeId,
    #[case] expected: Option<NodeId>,
) {
    assert_eq!(id.try_map(&participant_correspondence), expected);
}

#[rstest]
#[case::before_removed(NodeId(0), Some(NodeId(0)))]
#[case::removed(NodeId(1), None)]
#[case::after_removed(NodeId(2), Some(NodeId(1)))]
fn test_node_id_compact(#[case] id: NodeId, #[case] expected: Option<NodeId>) {
    let compaction = GraphCompaction::new(
        Compaction::new(3, vec![NodeId(1)]).unwrap(),
        Compaction::empty(),
    );
    assert_eq!(id.compact(&compaction), expected);
}

#[rstest]
#[case::before_gap(NodeId(0), NodeId(0))]
#[case::after_gap(NodeId(1), NodeId(2))]
fn test_node_id_uncompact(#[case] id: NodeId, #[case] expected: NodeId) {
    let compaction = GraphCompaction::new(
        Compaction::new(3, vec![NodeId(1)]).unwrap(),
        Compaction::empty(),
    );
    assert_eq!(id.uncompact(&compaction), expected);
}

#[rstest]
fn test_node_id_refs() {
    assert_eq!(
        NodeId(3).refs(),
        ParticipantRefs {
            node: Some(NodeId(3)),
            edge: None,
        }
    );
}

#[rstest]
#[case::mapped(EdgeId(2), Some(EdgeId(3)))]
#[case::unmatched(EdgeId(1), None)]
#[case::outside(EdgeId(4), None)]
fn test_edge_id_try_map(
    participant_correspondence: GraphCorrespondence,
    #[case] id: EdgeId,
    #[case] expected: Option<EdgeId>,
) {
    assert_eq!(id.try_map(&participant_correspondence), expected);
}

#[rstest]
#[case::removed(EdgeId(0), None)]
#[case::after_removed(EdgeId(2), Some(EdgeId(1)))]
fn test_edge_id_compact(#[case] id: EdgeId, #[case] expected: Option<EdgeId>) {
    let compaction = GraphCompaction::new(
        Compaction::empty(),
        Compaction::new(3, vec![EdgeId(0)]).unwrap(),
    );
    assert_eq!(id.compact(&compaction), expected);
}

#[rstest]
#[case::before_gap(EdgeId(0), EdgeId(0))]
#[case::after_gap(EdgeId(1), EdgeId(2))]
fn test_edge_id_uncompact(#[case] id: EdgeId, #[case] expected: EdgeId) {
    let compaction = GraphCompaction::new(
        Compaction::empty(),
        Compaction::new(3, vec![EdgeId(1)]).unwrap(),
    );
    assert_eq!(id.uncompact(&compaction), expected);
}

#[rstest]
fn test_edge_id_refs() {
    assert_eq!(
        EdgeId(2).refs(),
        ParticipantRefs {
            node: None,
            edge: Some(EdgeId(2)),
        }
    );
}
