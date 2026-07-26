//! Cross-validation properties for complete common-subgraph enumeration.
//!
//! The modular-product and direct searches use different intermediate
//! representations but must return the same complete, sorted correspondence
//! vector for every labeled graph pair and embedding kind.

use proptest::prelude::*;
use umol_graph_core::{
    CommonSubgraphEnumerationAlgorithm, EdgeId, EmbeddingKind, Graph, GraphCorrespondence, NodeId,
};

#[derive(Clone, Debug)]
struct LabeledGraph {
    graph: Graph,
    node_labels: Vec<u8>,
    edge_labels: Vec<u8>,
}

fn labeled_graph(max_nodes: usize) -> impl Strategy<Value = LabeledGraph> {
    (0..=max_nodes)
        .prop_flat_map(|node_count| {
            let pairs: Vec<[u32; 2]> = (0..node_count as u32)
                .flat_map(|first| {
                    ((first + 1)..node_count as u32).map(move |second| [first, second])
                })
                .collect();
            (
                Just(node_count),
                Just(pairs.clone()),
                prop::collection::vec(0_u8..3, node_count),
                prop::collection::vec(
                    prop_oneof![Just(None), (0_u8..3).prop_map(Some)],
                    pairs.len(),
                ),
            )
        })
        .prop_map(|(node_count, pairs, node_labels, edge_labels)| {
            let (edges, edge_labels): (Vec<_>, Vec<_>) = pairs
                .into_iter()
                .zip(edge_labels)
                .filter_map(|(edge, label)| label.map(|label| (edge, label)))
                .unzip();
            LabeledGraph {
                graph: Graph::new(node_count, &edges),
                node_labels,
                edge_labels,
            }
        })
}

fn enumerate(
    left: &LabeledGraph,
    right: &LabeledGraph,
    embedding: EmbeddingKind,
    algorithm: CommonSubgraphEnumerationAlgorithm,
) -> Vec<GraphCorrespondence> {
    left.graph.enumerate_common_subgraphs(
        &right.graph,
        &mut |left_node: NodeId, right_node: NodeId| {
            left.node_labels[left_node.index()] == right.node_labels[right_node.index()]
        },
        &mut |left_edge: EdgeId, right_edge: EdgeId| {
            left.edge_labels[left_edge.index()] == right.edge_labels[right_edge.index()]
        },
        embedding,
        algorithm,
    )
}

proptest! {
    #[test]
    fn test_graph_enumerate_common_subgraphs_cross_validation(
        left in labeled_graph(4),
        right in labeled_graph(5),
    ) {
        for embedding in [EmbeddingKind::Induced, EmbeddingKind::Monomorphism] {
            let modular_product = enumerate(
                &left,
                &right,
                embedding,
                CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
            );
            let direct = enumerate(
                &left,
                &right,
                embedding,
                CommonSubgraphEnumerationAlgorithm::DirectBacktracking,
            );
            prop_assert_eq!(direct, modular_product);
        }
    }
}
