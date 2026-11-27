//! Unit tests for the partitions and the search.

use std::ops::ControlFlow;

use rstest::rstest;

use super::*;

fn direct_graph_adapter(source: &Graph) -> AutomorphismAdapter {
    AutomorphismAdapter {
        graph: source.clone(),
        node_colors: vec![0; source.node_count()],
        edge_colors: vec![0; source.edge_count()],
    }
}

#[rstest]
#[case::empty(
    4,
    Vec::new(),
    vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
)]
#[case::one_generator(
    4,
    vec![vec![NodeId(1), NodeId(0), NodeId(2), NodeId(3)]],
    vec![NodeId(0), NodeId(0), NodeId(2), NodeId(3)],
)]
#[case::generated_subgroup(
    4,
    vec![
        vec![NodeId(1), NodeId(0), NodeId(2), NodeId(3)],
        vec![NodeId(0), NodeId(2), NodeId(1), NodeId(3)],
    ],
    vec![NodeId(0), NodeId(0), NodeId(0), NodeId(3)],
)]
fn test_orbits_from_generators(
    #[case] node_count: usize,
    #[case] generators: Vec<Vec<NodeId>>,
    #[case] expected: Vec<NodeId>,
) {
    assert_eq!(orbits_from_generators(node_count, &generators), expected);
}

#[rstest]
#[case::ordered_classes(
        vec![2_u32, 0, 2, 1, 0],
        OrderedPartition {
            cells: vec![
                vec![NodeId(1), NodeId(4)],
                vec![NodeId(3)],
                vec![NodeId(0), NodeId(2)],
            ],
        },
    )]
fn test_ordered_partition_from_descriptors(
    #[case] descriptors: Vec<u32>,
    #[case] expected: OrderedPartition,
) {
    assert_eq!(OrderedPartition::from_descriptors(&descriptors), expected,);
}

#[rstest]
#[case::empty(Vec::new(), OrderedPartition { cells: Vec::new() })]
#[case::repeated(
    vec![2, 0, 2, 1, 0],
    OrderedPartition {
        cells: vec![
            vec![NodeId(1), NodeId(4)],
            vec![NodeId(3)],
            vec![NodeId(0), NodeId(2)],
        ],
    },
)]
#[case::ordered_distinct(
    vec![0, 1, 2],
    OrderedPartition {
        cells: vec![vec![NodeId(0)], vec![NodeId(1)], vec![NodeId(2)]],
    },
)]
fn test_ordered_partition_from_color_ranks(
    #[case] colors: Vec<u32>,
    #[case] expected: OrderedPartition,
) {
    assert_eq!(OrderedPartition::from_color_ranks(colors), expected);
}

#[rstest]
#[case::empty(
    OrderedPartition { cells: Vec::new() },
    Graph::new(0, &[]),
    Vec::new(),
    OrderedPartition { cells: Vec::new() },
)]
#[case::discrete(
    OrderedPartition {
        cells: vec![vec![NodeId(0)], vec![NodeId(1)], vec![NodeId(2)]],
    },
    Graph::new(3, &[[0, 1], [1, 2]]),
    vec![0, 0],
    OrderedPartition {
        cells: vec![vec![NodeId(0)], vec![NodeId(1)], vec![NodeId(2)]],
    },
)]
#[case::unsplit(
    OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)]],
    },
    Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
    vec![0, 0, 0, 0],
    OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)]],
    },
)]
#[case::edge_color_square(
    OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)]],
    },
    Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
    vec![1, 0, 0, 0],
    OrderedPartition {
        cells: vec![
            vec![NodeId(2), NodeId(3)],
            vec![NodeId(0), NodeId(1)],
        ],
    },
)]
#[case::edge_color_cascade(
    OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4), NodeId(5)]],
    },
    Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
    vec![1, 0, 0, 0, 0, 0],
    OrderedPartition {
        cells: vec![
            vec![NodeId(3), NodeId(4)],
            vec![NodeId(2), NodeId(5)],
            vec![NodeId(0), NodeId(1)],
        ],
    },
)]
#[case::edge_color_multiplicity(
    OrderedPartition {
        cells: vec![
            vec![NodeId(0), NodeId(1)],
            vec![NodeId(2), NodeId(3), NodeId(4)],
        ],
    },
    Graph::new(5, &[[0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4]]),
    vec![0, 0, 1, 0, 1, 1],
    OrderedPartition {
        cells: vec![
            vec![NodeId(0)],
            vec![NodeId(1)],
            vec![NodeId(2)],
            vec![NodeId(3)],
            vec![NodeId(4)],
        ],
    },
)]
#[case::multiply_split(
    OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4)]],
    },
    Graph::new(5, &[[0, 1], [1, 2], [2, 3], [3, 4]]),
    vec![0, 0, 0, 0],
    OrderedPartition {
        cells: vec![
            vec![NodeId(2)],
            vec![NodeId(1), NodeId(3)],
            vec![NodeId(0), NodeId(4)],
        ],
    },
)]
#[case::disconnected(
    OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4)]],
    },
    Graph::new(5, &[[0, 1], [2, 3]]),
    vec![0, 0],
    OrderedPartition {
        cells: vec![
            vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
            vec![NodeId(4)],
        ],
    },
)]
#[case::individualized(
    OrderedPartition {
        cells: vec![
            vec![NodeId(2)],
            vec![NodeId(0), NodeId(1), NodeId(3)],
        ],
    },
    Graph::new(4, &[[0, 1], [0, 2], [0, 3]]),
    vec![0, 0, 0],
    OrderedPartition {
        cells: vec![
            vec![NodeId(2)],
            vec![NodeId(0)],
            vec![NodeId(1), NodeId(3)],
        ],
    },
)]
fn test_ordered_partition_refine(
    #[case] partition: OrderedPartition,
    #[case] graph: Graph,
    #[case] edge_colors: Vec<u32>,
    #[case] expected: OrderedPartition,
) {
    assert_eq!(partition.clone().refine(&graph, &edge_colors), expected);
    assert_eq!(
        partition
            .clone()
            .refine_by_active_cells(&graph, &edge_colors),
        expected,
    );
    assert_eq!(
        partition.refine_by_all_cells(&graph, &edge_colors),
        expected,
    );
}

#[rstest]
#[case::naphthalene(
    10,
    &[
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
        (5, 6), (6, 7), (7, 8), (8, 9), (9, 4),
    ],
    vec![0; 11],
)]
#[case::naphthalene_kekule(
    10,
    &[
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
        (5, 6), (6, 7), (7, 8), (8, 9), (9, 4),
    ],
    vec![2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1],
)]
#[case::disconnected_rings(
    12,
    &[
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
        (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 6),
    ],
    vec![0; 12],
)]
#[case::long_path(
    77,
    &(0_u32..76)
        .map(|first| (first, first + 1))
        .collect::<Vec<_>>(),
    vec![0; 76],
)]
#[case::long_colored_path(
    77,
    &(0_u32..76)
        .map(|first| (first, first + 1))
        .collect::<Vec<_>>(),
    (0_u32..76).map(|index| (index * 3 + 1) % 5).collect(),
)]
fn test_ordered_partition_refine_active_cell_equivalence(
    #[case] atom_count: usize,
    #[case] edges: &[(u32, u32)],
    #[case] edge_colors: Vec<u32>,
) {
    let graph = Graph::new(
        atom_count,
        &edges
            .iter()
            .map(|&(first, second)| [first, second])
            .collect::<Vec<_>>(),
    );
    let partition = OrderedPartition::from_descriptors(&vec![0_u32; atom_count]);

    assert_eq!(
        partition.clone().refine(&graph, &edge_colors),
        partition.clone().refine_by_all_cells(&graph, &edge_colors),
    );
    assert_eq!(
        partition
            .clone()
            .refine_by_active_cells(&graph, &edge_colors),
        partition.refine_by_all_cells(&graph, &edge_colors),
    );
}

#[rstest]
#[case::first_cell(
        OrderedPartition {
            cells: vec![vec![NodeId(0), NodeId(3)], vec![NodeId(1), NodeId(2)]],
        },
        0,
        NodeId(3),
        OrderedPartition {
            cells: vec![
                vec![NodeId(3)],
                vec![NodeId(0)],
                vec![NodeId(1), NodeId(2)],
            ],
        },
    )]
fn test_ordered_partition_individualize(
    #[case] partition: OrderedPartition,
    #[case] cell_index: usize,
    #[case] node: NodeId,
    #[case] expected: OrderedPartition,
) {
    assert_eq!(partition.individualize(cell_index, node), expected);
}

#[rstest]
#[case::singleton(
    1,
    CanonicalSearchStats {
        initial_residual_cell_sizes: vec![],
        refinement_calls: 1,
        branch_order_calls: 0,
        backend_calls: 0,
        visited_leaves: 1,
        leaf_comparisons: 0,
        orbit_pruned_branches: 0,
    },
)]
#[case::orbit_pruned(
    2,
    CanonicalSearchStats {
        initial_residual_cell_sizes: vec![2],
        refinement_calls: 2,
        branch_order_calls: 1,
        backend_calls: 1,
        visited_leaves: 1,
        leaf_comparisons: 0,
        orbit_pruned_branches: 1,
    },
)]
fn test_canonical_search_stats(#[case] node_count: usize, #[case] expected: CanonicalSearchStats) {
    let source = Graph::new(node_count, &[]);
    let adapter = direct_graph_adapter(&source);
    let leaf_candidate = |order: &[NodeId]| CanonicalCandidate {
        key: order.to_vec(),
        entity_order: order.to_vec(),
    };
    let actual = canonical_search(
        &adapter,
        adapter.node_colors(),
        AutomorphismAlgorithm::Nauty,
        &leaf_candidate,
    );

    assert_eq!(actual.stats, expected);
}

#[rstest]
#[case::empty(vec![], vec![], vec![vec![]])]
#[case::discrete(vec![0, 1, 2], vec![0, 1, 2], vec![vec![0, 1, 2]])]
#[case::unpruned(
    vec![0, 0, 0],
    vec![0, 0, 0],
    vec![vec![0, 1, 2], vec![0, 2, 1], vec![1, 0, 2], vec![1, 2, 0], vec![2, 0, 1], vec![2, 1, 0]],
)]
#[case::incompatible(vec![0, 0, 0], vec![0, 1, 2], vec![])]
fn test_visit_remapping_candidates(
    #[case] source_colors: Vec<u32>,
    #[case] target_colors: Vec<u32>,
    #[case] expected: Vec<Vec<u32>>,
) {
    let adapter = AutomorphismAdapter {
        graph: Graph::new(source_colors.len(), &[]),
        node_colors: vec![0; source_colors.len()],
        edge_colors: Vec::new(),
    };
    let mut actual = Vec::new();
    let flow: ControlFlow<()> = visit_remapping_candidates(
        &adapter,
        OrderedPartition::from_descriptors(&source_colors),
        &adapter,
        OrderedPartition::from_descriptors(&target_colors),
        AutomorphismAlgorithm::Nauty,
        &mut |images: &[NodeId]| {
            actual.push(images.iter().map(|node| node.0).collect::<Vec<_>>());
            ControlFlow::Continue(())
        },
        true,
        &|generators| generators.clear(),
    );
    actual.sort();

    assert_eq!(flow, ControlFlow::Continue(()));
    assert_eq!(actual, expected);
}

#[rstest]
#[case::first(vec![NodeId(0), NodeId(1), NodeId(2)])]
#[case::later(vec![NodeId(1), NodeId(0), NodeId(2)])]
fn test_visit_remapping_candidates_termination(#[case] accepted: Vec<NodeId>) {
    let adapter = AutomorphismAdapter {
        graph: Graph::new(3, &[]),
        node_colors: vec![0; 3],
        edge_colors: Vec::new(),
    };
    let partition = OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2)]],
    };
    let mut found = false;
    let mut calls_after_break = 0;
    let flow = visit_remapping_candidates(
        &adapter,
        partition.clone(),
        &adapter,
        partition,
        AutomorphismAlgorithm::Nauty,
        &mut |images: &[NodeId]| {
            if found {
                calls_after_break += 1;
            }
            if images == accepted {
                found = true;
                ControlFlow::Break(images.to_vec())
            } else {
                ControlFlow::Continue(())
            }
        },
        true,
        &|generators| generators.clear(),
    );

    assert_eq!(flow, ControlFlow::Break(accepted));
    assert_eq!((found, calls_after_break), (true, 0));
}

#[rstest]
#[case::empty(Vec::new(), 2, vec![0, 1], vec![0, 1], 2)]
#[case::distinct(
    vec![(20, 0), (10, 1), (30, 2)],
    0,
    vec![0, 0, 0],
    vec![1, 0, 2],
    3
)]
#[case::equal(
    vec![(20, 0), (10, 1), (20, 2)],
    0,
    vec![0, 0, 0],
    vec![1, 0, 1],
    2
)]
#[case::first_above_zero(
    vec![(20, 2), (10, 3), (20, 4)],
    2,
    vec![0, 1, 0, 0, 0],
    vec![0, 1, 3, 2, 3],
    4
)]
fn test_rank_colors(
    #[case] colors: Vec<(i32, usize)>,
    #[case] first: u32,
    #[case] mut ranks: Vec<u32>,
    #[case] expected_ranks: Vec<u32>,
    #[case] expected_next: u32,
) {
    assert_eq!(rank_colors(colors, first, &mut ranks), expected_next);
    assert_eq!(ranks, expected_ranks);
}
