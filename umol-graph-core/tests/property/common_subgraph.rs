//! Cross-validation and delivery-agreement properties for common-subgraph
//! enumeration.
//!
//! The modular-product and direct searches use different intermediate
//! representations but must return the same complete, sorted correspondence
//! vector for every labeled graph pair and embedding kind; and for every
//! selector — complete and maximal — the visitor emissions must equal the
//! eager enumeration as normalized node-pair sets.

use std::ops::ControlFlow;

use proptest::prelude::*;
use umol_graph_core::{
    CommonSubgraphEnumerationAlgorithm, EdgeId, EmbeddingKind, Graph, GraphCorrespondence,
    MaximalCommonSubgraphAlgorithm, NodeId,
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

/// Sorted, per-emission-sorted node-pair sets collected through the complete
/// visitor, with label-equality node and edge matching.
fn visited(
    left: &LabeledGraph,
    right: &LabeledGraph,
    embedding: EmbeddingKind,
    algorithm: CommonSubgraphEnumerationAlgorithm,
) -> Vec<Vec<(NodeId, NodeId)>> {
    let mut found: Vec<Vec<(NodeId, NodeId)>> = Vec::new();
    let _: ControlFlow<()> = left.graph.visit_common_subgraphs(
        &right.graph,
        &mut |left_node: NodeId, right_node: NodeId| {
            left.node_labels[left_node.index()] == right.node_labels[right_node.index()]
        },
        &mut |left_edge: EdgeId, right_edge: EdgeId| {
            left.edge_labels[left_edge.index()] == right.edge_labels[right_edge.index()]
        },
        embedding,
        algorithm,
        |pairs| {
            let mut pairs = pairs.to_vec();
            pairs.sort_unstable();
            found.push(pairs);
            ControlFlow::Continue(())
        },
    );
    found.sort();
    found
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

    #[test]
    fn test_graph_visit_common_subgraphs_agreement(
        left in labeled_graph(4),
        right in labeled_graph(5),
    ) {
        for embedding in [EmbeddingKind::Induced, EmbeddingKind::Monomorphism] {
            for algorithm in [
                CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
                CommonSubgraphEnumerationAlgorithm::DirectBacktracking,
            ] {
                let mut enumerated: Vec<Vec<(NodeId, NodeId)>> =
                    enumerate(&left, &right, embedding, algorithm)
                        .iter()
                        .map(|subgraph| subgraph.nodes().matched_pairs().to_vec())
                        .collect();
                enumerated.sort();
                prop_assert_eq!(
                    visited(&left, &right, embedding, algorithm),
                    enumerated,
                    "{:?} {:?}: visitor disagrees with enumeration",
                    algorithm,
                    embedding
                );
            }
        }
    }

    #[test]
    fn test_graph_visit_maximal_common_subgraphs_agreement(
        left in labeled_graph(4),
        right in labeled_graph(5),
    ) {
        for embedding in [EmbeddingKind::Induced, EmbeddingKind::Monomorphism] {
            let mut found: Vec<Vec<(NodeId, NodeId)>> = Vec::new();
            let _: ControlFlow<()> = left.graph.visit_maximal_common_subgraphs(
                &right.graph,
                &mut |left_node: NodeId, right_node: NodeId| {
                    left.node_labels[left_node.index()] == right.node_labels[right_node.index()]
                },
                &mut |left_edge: EdgeId, right_edge: EdgeId| {
                    left.edge_labels[left_edge.index()] == right.edge_labels[right_edge.index()]
                },
                embedding,
                MaximalCommonSubgraphAlgorithm::BronKerbosch,
                |pairs| {
                    let mut pairs = pairs.to_vec();
                    pairs.sort_unstable();
                    found.push(pairs);
                    ControlFlow::Continue(())
                },
            );
            found.sort();
            let mut enumerated: Vec<Vec<(NodeId, NodeId)>> = left
                .graph
                .enumerate_maximal_common_subgraphs(
                    &right.graph,
                    &mut |left_node: NodeId, right_node: NodeId| {
                        left.node_labels[left_node.index()] == right.node_labels[right_node.index()]
                    },
                    &mut |left_edge: EdgeId, right_edge: EdgeId| {
                        left.edge_labels[left_edge.index()] == right.edge_labels[right_edge.index()]
                    },
                    embedding,
                    MaximalCommonSubgraphAlgorithm::BronKerbosch,
                )
                .iter()
                .map(|subgraph| subgraph.nodes().matched_pairs().to_vec())
                .collect();
            enumerated.sort();
            prop_assert_eq!(
                found,
                enumerated,
                "{:?}: visitor disagrees with enumeration",
                embedding
            );
        }
    }
}
