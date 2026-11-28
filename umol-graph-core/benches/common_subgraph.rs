//! Common-subgraph enumeration benchmarks (collection, first match, and visitation count).

#[allow(dead_code)]
mod fixtures;

use std::hint::black_box;
use std::ops::ControlFlow;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_graph_core::{CommonSubgraphEnumerationAlgorithm, EdgeId, EmbeddingKind, Graph, NodeId};

struct CommonSubgraphCase {
    name: &'static str,
    expected_output_count: usize,
    left: Graph,
    right: Graph,
    left_node_labels: Vec<u8>,
    right_node_labels: Vec<u8>,
    left_edge_labels: Vec<u8>,
    right_edge_labels: Vec<u8>,
    embedding: EmbeddingKind,
}

fn common_subgraph_enumeration(c: &mut Criterion) {
    let triangle = fixtures::cycle(3);
    let path_4 = fixtures::path(4);
    let cycle_4 = fixtures::cycle(4);
    let path_5 = fixtures::path(5);
    let cycle_6 = fixtures::cycle(6);
    let benzene = fixtures::cycle(6);
    let naphthalene = fixtures::naphthalene();

    let cases = [
        CommonSubgraphCase {
            name: "dense_compatible/triangle_triangle/induced",
            expected_output_count: 34,
            left_node_labels: vec![0; triangle.node_count()],
            right_node_labels: vec![0; triangle.node_count()],
            left_edge_labels: vec![0; triangle.edge_count()],
            right_edge_labels: vec![0; triangle.edge_count()],
            left: triangle.clone(),
            right: triangle,
            embedding: EmbeddingKind::Induced,
        },
        CommonSubgraphCase {
            name: "structural/path_4_cycle_4/induced",
            expected_output_count: 69,
            left_node_labels: vec![0; path_4.node_count()],
            right_node_labels: vec![0; cycle_4.node_count()],
            left_edge_labels: vec![0; path_4.edge_count()],
            right_edge_labels: vec![0; cycle_4.edge_count()],
            left: path_4.clone(),
            right: cycle_4.clone(),
            embedding: EmbeddingKind::Induced,
        },
        CommonSubgraphCase {
            name: "structural/path_4_cycle_4/monomorphism",
            expected_output_count: 209,
            left_node_labels: vec![0; path_4.node_count()],
            right_node_labels: vec![0; cycle_4.node_count()],
            left_edge_labels: vec![0; path_4.edge_count()],
            right_edge_labels: vec![0; cycle_4.edge_count()],
            left: path_4,
            right: cycle_4,
            embedding: EmbeddingKind::Monomorphism,
        },
        CommonSubgraphCase {
            name: "label_selective/path_5_cycle_6/induced",
            expected_output_count: 109,
            left_node_labels: vec![0, 1, 0, 1, 0],
            right_node_labels: vec![0, 1, 0, 1, 0, 1],
            left_edge_labels: fixtures::alternating_edge_labels(&path_5),
            right_edge_labels: fixtures::alternating_edge_labels(&cycle_6),
            left: path_5,
            right: cycle_6,
            embedding: EmbeddingKind::Induced,
        },
        CommonSubgraphCase {
            name: "molecular/benzene_naphthalene/induced",
            expected_output_count: 1_957,
            left_node_labels: benzene
                .node_ids()
                .map(|node| benzene.degree(node) as u8)
                .collect(),
            right_node_labels: naphthalene
                .node_ids()
                .map(|node| naphthalene.degree(node) as u8)
                .collect(),
            left_edge_labels: vec![0; benzene.edge_count()],
            right_edge_labels: vec![0; naphthalene.edge_count()],
            left: benzene,
            right: naphthalene,
            embedding: EmbeddingKind::Induced,
        },
    ];

    let mut group = c.benchmark_group("common_subgraph_enumeration");
    for case in cases {
        for (algorithm_name, algorithm) in [
            (
                "modular_product_backtracking",
                CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
            ),
            (
                "direct_backtracking",
                CommonSubgraphEnumerationAlgorithm::DirectBacktracking,
            ),
        ] {
            let mut node_match = |left: NodeId, right: NodeId| {
                case.left_node_labels[left.index()] == case.right_node_labels[right.index()]
            };
            let mut edge_match = |left: EdgeId, right: EdgeId| {
                case.left_edge_labels[left.index()] == case.right_edge_labels[right.index()]
            };
            let output = case.left.enumerate_common_subgraphs(
                &case.right,
                &mut node_match,
                &mut edge_match,
                case.embedding,
                algorithm,
            );
            assert_eq!(
                output.len(),
                case.expected_output_count,
                "unexpected output count for {} with {algorithm_name}",
                case.name,
            );

            group.bench_with_input(
                BenchmarkId::new(case.name, algorithm_name),
                &algorithm,
                |b, &algorithm| {
                    b.iter(|| {
                        let mut node_match = |left: NodeId, right: NodeId| {
                            case.left_node_labels[left.index()]
                                == case.right_node_labels[right.index()]
                        };
                        let mut edge_match = |left: EdgeId, right: EdgeId| {
                            case.left_edge_labels[left.index()]
                                == case.right_edge_labels[right.index()]
                        };
                        black_box(case.left.enumerate_common_subgraphs(
                            &case.right,
                            &mut node_match,
                            &mut edge_match,
                            case.embedding,
                            algorithm,
                        ))
                    });
                },
            );
            group.bench_with_input(
                BenchmarkId::new(case.name, format!("{algorithm_name}/first")),
                &algorithm,
                |b, &algorithm| {
                    b.iter(|| {
                        let mut node_match = |left: NodeId, right: NodeId| {
                            case.left_node_labels[left.index()]
                                == case.right_node_labels[right.index()]
                        };
                        let mut edge_match = |left: EdgeId, right: EdgeId| {
                            case.left_edge_labels[left.index()]
                                == case.right_edge_labels[right.index()]
                        };
                        case.left.visit_common_subgraphs(
                            &case.right,
                            &mut node_match,
                            &mut edge_match,
                            case.embedding,
                            algorithm,
                            |_| ControlFlow::Break(()),
                        )
                    });
                },
            );
            group.bench_with_input(
                BenchmarkId::new(case.name, format!("{algorithm_name}/count")),
                &algorithm,
                |b, &algorithm| {
                    b.iter(|| {
                        let mut node_match = |left: NodeId, right: NodeId| {
                            case.left_node_labels[left.index()]
                                == case.right_node_labels[right.index()]
                        };
                        let mut edge_match = |left: EdgeId, right: EdgeId| {
                            case.left_edge_labels[left.index()]
                                == case.right_edge_labels[right.index()]
                        };
                        let mut count = 0usize;
                        let _: ControlFlow<()> = case.left.visit_common_subgraphs(
                            &case.right,
                            &mut node_match,
                            &mut edge_match,
                            case.embedding,
                            algorithm,
                            |_| {
                                count += 1;
                                ControlFlow::Continue(())
                            },
                        );
                        count
                    });
                },
            );
        }
    }
    group.finish();
}

criterion_group!(benches, common_subgraph_enumeration);
criterion_main!(benches);
