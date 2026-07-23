//! Definition-level exhaustive cycle operations for small test graphs.

use std::cmp::Ordering;

use umol_graph_core::{EdgeId, Graph, NodeId};

pub(super) fn enumerate_cycles(graph: &Graph) -> Vec<Vec<EdgeId>> {
    assert!(
        graph.edge_count() < usize::BITS as usize,
        "the exhaustive oracle requires an edge subset to fit in usize"
    );

    let mut cycles = Vec::new();
    for mask in 1..(1usize << graph.edge_count()) {
        let edges: Vec<EdgeId> = graph
            .edge_ids()
            .filter(|edge| mask & (1usize << edge.index()) != 0)
            .collect();
        if is_cycle(graph, &edges) {
            cycles.push(edges);
        }
    }
    cycles.sort_by(|left, right| compare_cycles(left, right));
    cycles
}

pub(super) fn cycle_space_rank(graph: &Graph) -> usize {
    graph.edge_count() + component_count(graph) - graph.node_count()
}

pub(super) fn are_linearly_independent(cycles: &[Vec<EdgeId>], edge_count: usize) -> bool {
    cycle_vector_rank(cycles, edge_count) == cycles.len()
}

pub(super) fn enumerate_cycle_bases(graph: &Graph) -> Vec<Vec<Vec<EdgeId>>> {
    let cycles = enumerate_cycles(graph);
    let rank = cycle_space_rank(graph);
    let mut selected = Vec::with_capacity(rank);
    let mut bases = Vec::new();
    collect_cycle_bases(
        &cycles,
        graph.edge_count(),
        rank,
        0,
        &mut selected,
        &mut bases,
    );
    bases
}

pub(super) fn minimum_cycle_bases(graph: &Graph) -> Vec<Vec<Vec<EdgeId>>> {
    let mut bases = enumerate_cycle_bases(graph);
    let minimum_weight = bases
        .iter()
        .map(|basis| basis.iter().map(Vec::len).sum())
        .min()
        .expect("the graph cycle set spans its cycle space");
    bases.retain(|basis| basis.iter().map(Vec::len).sum::<usize>() == minimum_weight);
    bases
}

pub(super) fn relevant_cycles(graph: &Graph) -> Vec<Vec<EdgeId>> {
    let mut cycles: Vec<Vec<EdgeId>> = minimum_cycle_bases(graph).into_iter().flatten().collect();
    cycles.sort_by(|left, right| compare_cycles(left, right));
    cycles.dedup();
    cycles
}

pub(super) fn unique_ring_families(graph: &Graph) -> Vec<Vec<Vec<EdgeId>>> {
    let cycles = enumerate_cycles(graph);
    let relevant = relevant_cycles(graph);
    let mut visited = vec![false; relevant.len()];
    let mut families = Vec::new();

    for start in 0..relevant.len() {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut stack = vec![start];
        let mut family = Vec::new();
        while let Some(current) = stack.pop() {
            family.push(relevant[current].clone());
            for candidate in 0..relevant.len() {
                if !visited[candidate]
                    && are_urf_related(
                        &relevant[current],
                        &relevant[candidate],
                        &cycles,
                        graph.edge_count(),
                    )
                {
                    visited[candidate] = true;
                    stack.push(candidate);
                }
            }
        }
        family.sort_by(|left, right| compare_cycles(left, right));
        families.push(family);
    }

    families.sort_by(|left, right| compare_cycles(&left[0], &right[0]));
    families
}

fn cycle_vector_rank(cycles: &[Vec<EdgeId>], edge_count: usize) -> usize {
    let mut basis: Vec<Option<Vec<bool>>> = vec![None; edge_count];
    let mut rank = 0;

    for cycle in cycles {
        let mut row = vec![false; edge_count];
        for &edge in cycle {
            row[edge.index()] ^= true;
        }

        for pivot in 0..edge_count {
            if !row[pivot] {
                continue;
            }
            if let Some(basis_row) = &basis[pivot] {
                for column in pivot..edge_count {
                    row[column] ^= basis_row[column];
                }
            } else {
                basis[pivot] = Some(row);
                rank += 1;
                break;
            }
        }
    }

    rank
}

fn collect_cycle_bases(
    cycles: &[Vec<EdgeId>],
    edge_count: usize,
    rank: usize,
    start: usize,
    selected: &mut Vec<Vec<EdgeId>>,
    bases: &mut Vec<Vec<Vec<EdgeId>>>,
) {
    if selected.len() == rank {
        bases.push(selected.clone());
        return;
    }
    let remaining = rank - selected.len();
    if cycles.len() - start < remaining {
        return;
    }

    for index in start..=cycles.len() - remaining {
        selected.push(cycles[index].clone());
        if are_linearly_independent(selected, edge_count) {
            collect_cycle_bases(cycles, edge_count, rank, index + 1, selected, bases);
        }
        selected.pop();
    }
}

fn are_urf_related(
    first: &[EdgeId],
    second: &[EdgeId],
    cycles: &[Vec<EdgeId>],
    edge_count: usize,
) -> bool {
    if first.len() != second.len() || !first.iter().any(|edge| second.binary_search(edge).is_ok()) {
        return false;
    }

    let smaller: Vec<Vec<EdgeId>> = cycles
        .iter()
        .filter(|cycle| cycle.len() < first.len())
        .cloned()
        .collect();
    let smaller_rank = cycle_vector_rank(&smaller, edge_count);
    let difference: Vec<EdgeId> = (0..edge_count)
        .map(|edge| EdgeId(edge as u32))
        .filter(|edge| first.binary_search(edge).is_ok() ^ second.binary_search(edge).is_ok())
        .collect();
    let mut with_difference = smaller;
    with_difference.push(difference);
    cycle_vector_rank(&with_difference, edge_count) == smaller_rank
}

fn compare_cycles(first: &[EdgeId], second: &[EdgeId]) -> Ordering {
    first
        .len()
        .cmp(&second.len())
        .then_with(|| first.cmp(second))
}

fn is_cycle(graph: &Graph, edges: &[EdgeId]) -> bool {
    let mut degrees = vec![0usize; graph.node_count()];
    for &edge in edges {
        let [first, second] = graph.edge_endpoints(edge);
        degrees[first.index()] += 1;
        degrees[second.index()] += 1;
    }
    if degrees.iter().any(|&degree| degree != 0 && degree != 2) {
        return false;
    }

    let Some(start) = degrees
        .iter()
        .position(|&degree| degree != 0)
        .map(|node| NodeId(node as u32))
    else {
        return false;
    };
    let mut visited = vec![false; graph.node_count()];
    let mut stack = vec![start];
    visited[start.index()] = true;
    while let Some(node) = stack.pop() {
        for &edge in edges {
            let [first, second] = graph.edge_endpoints(edge);
            let neighbor = if first == node {
                Some(second)
            } else if second == node {
                Some(first)
            } else {
                None
            };
            if let Some(neighbor) = neighbor {
                if !visited[neighbor.index()] {
                    visited[neighbor.index()] = true;
                    stack.push(neighbor);
                }
            }
        }
    }

    degrees
        .iter()
        .enumerate()
        .all(|(node, &degree)| degree == 0 || visited[node])
}

fn component_count(graph: &Graph) -> usize {
    let mut count = 0;
    let mut visited = vec![false; graph.node_count()];

    for start in graph.node_ids() {
        if visited[start.index()] {
            continue;
        }
        count += 1;
        visited[start.index()] = true;
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            for neighbor in graph.neighbors(node) {
                if !visited[neighbor.node.index()] {
                    visited[neighbor.node.index()] = true;
                    stack.push(neighbor.node);
                }
            }
        }
    }

    count
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::empty(0, &[], vec![])]
    #[case::isolated(2, &[], vec![])]
    #[case::loop_edge(1, &[[0, 0]], vec![vec![EdgeId(0)]])]
    #[case::parallel_pair(
        2,
        &[[0, 1], [0, 1]],
        vec![vec![EdgeId(0), EdgeId(1)]],
    )]
    #[case::parallel_triple(
        2,
        &[[0, 1], [0, 1], [0, 1]],
        vec![
            vec![EdgeId(0), EdgeId(1)],
            vec![EdgeId(0), EdgeId(2)],
            vec![EdgeId(1), EdgeId(2)],
        ],
    )]
    #[case::triangle(
        3,
        &[[0, 1], [1, 2], [0, 2]],
        vec![vec![EdgeId(0), EdgeId(1), EdgeId(2)]],
    )]
    #[case::square_with_diagonal(
        4,
        &[[0, 1], [1, 2], [2, 3], [0, 3], [0, 2]],
        vec![
            vec![EdgeId(0), EdgeId(1), EdgeId(4)],
            vec![EdgeId(2), EdgeId(3), EdgeId(4)],
            vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
        ],
    )]
    #[case::disconnected_loops(
        2,
        &[[0, 0], [1, 1]],
        vec![vec![EdgeId(0)], vec![EdgeId(1)]],
    )]
    fn test_enumerate_cycles(
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
        #[case] expected: Vec<Vec<EdgeId>>,
    ) {
        assert_eq!(enumerate_cycles(&Graph::new(node_count, edges)), expected);
    }

    #[rstest]
    #[case::empty(0, &[], 0)]
    #[case::isolated(3, &[], 0)]
    #[case::path(3, &[[0, 1], [1, 2]], 0)]
    #[case::triangle(3, &[[0, 1], [1, 2], [0, 2]], 1)]
    #[case::square_with_diagonal(4, &[[0, 1], [1, 2], [2, 3], [0, 3], [0, 2]], 2)]
    #[case::loop_edge(1, &[[0, 0]], 1)]
    #[case::parallel_triple(2, &[[0, 1], [0, 1], [0, 1]], 2)]
    #[case::disconnected(4, &[[0, 1], [1, 2], [0, 2], [3, 3]], 2)]
    fn test_cycle_space_rank(
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
        #[case] expected: usize,
    ) {
        assert_eq!(cycle_space_rank(&Graph::new(node_count, edges)), expected);
    }

    #[rstest]
    #[case::empty(&[], 0, true)]
    #[case::single(&[vec![EdgeId(0), EdgeId(1), EdgeId(4)]], 5, true)]
    #[case::independent(
        &[
            vec![EdgeId(0), EdgeId(1), EdgeId(4)],
            vec![EdgeId(2), EdgeId(3), EdgeId(4)],
        ],
        5,
        true,
    )]
    #[case::dependent(
        &[
            vec![EdgeId(0), EdgeId(1), EdgeId(4)],
            vec![EdgeId(2), EdgeId(3), EdgeId(4)],
            vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
        ],
        5,
        false,
    )]
    #[case::duplicate(
        &[
            vec![EdgeId(0), EdgeId(1)],
            vec![EdgeId(0), EdgeId(1)],
        ],
        2,
        false,
    )]
    #[case::loops(&[vec![EdgeId(0)], vec![EdgeId(1)]], 2, true)]
    fn test_are_linearly_independent(
        #[case] cycles: &[Vec<EdgeId>],
        #[case] edge_count: usize,
        #[case] expected: bool,
    ) {
        assert_eq!(are_linearly_independent(cycles, edge_count), expected);
    }

    #[rstest]
    #[case::path(3, &[[0, 1], [1, 2]], vec![vec![]])]
    #[case::square_with_diagonal(
        4,
        &[[0, 1], [1, 2], [2, 3], [0, 3], [0, 2]],
        vec![
            vec![
                vec![EdgeId(0), EdgeId(1), EdgeId(4)],
                vec![EdgeId(2), EdgeId(3), EdgeId(4)],
            ],
            vec![
                vec![EdgeId(0), EdgeId(1), EdgeId(4)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
            ],
            vec![
                vec![EdgeId(2), EdgeId(3), EdgeId(4)],
                vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
            ],
        ],
    )]
    #[case::parallel_triple(
        2,
        &[[0, 1], [0, 1], [0, 1]],
        vec![
            vec![
                vec![EdgeId(0), EdgeId(1)],
                vec![EdgeId(0), EdgeId(2)],
            ],
            vec![
                vec![EdgeId(0), EdgeId(1)],
                vec![EdgeId(1), EdgeId(2)],
            ],
            vec![
                vec![EdgeId(0), EdgeId(2)],
                vec![EdgeId(1), EdgeId(2)],
            ],
        ],
    )]
    fn test_enumerate_cycle_bases(
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
        #[case] expected: Vec<Vec<Vec<EdgeId>>>,
    ) {
        assert_eq!(
            enumerate_cycle_bases(&Graph::new(node_count, edges)),
            expected
        );
    }

    #[rstest]
    #[case::path(3, &[[0, 1], [1, 2]], vec![vec![]])]
    #[case::square_with_diagonal(
        4,
        &[[0, 1], [1, 2], [2, 3], [0, 3], [0, 2]],
        vec![vec![
            vec![EdgeId(0), EdgeId(1), EdgeId(4)],
            vec![EdgeId(2), EdgeId(3), EdgeId(4)],
        ]],
    )]
    #[case::unequal_theta(
        6,
        &[
            [0, 1],
            [1, 4],
            [0, 2],
            [2, 4],
            [0, 3],
            [3, 5],
            [5, 4],
        ],
        vec![
            vec![
                vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
                vec![EdgeId(0), EdgeId(1), EdgeId(4), EdgeId(5), EdgeId(6)],
            ],
            vec![
                vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
                vec![EdgeId(2), EdgeId(3), EdgeId(4), EdgeId(5), EdgeId(6)],
            ],
        ],
    )]
    fn test_minimum_cycle_bases(
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
        #[case] expected: Vec<Vec<Vec<EdgeId>>>,
    ) {
        assert_eq!(
            minimum_cycle_bases(&Graph::new(node_count, edges)),
            expected
        );
    }

    #[rstest]
    #[case::path(3, &[[0, 1], [1, 2]], vec![])]
    #[case::square_with_diagonal(
        4,
        &[[0, 1], [1, 2], [2, 3], [0, 3], [0, 2]],
        vec![
            vec![EdgeId(0), EdgeId(1), EdgeId(4)],
            vec![EdgeId(2), EdgeId(3), EdgeId(4)],
        ],
    )]
    #[case::unequal_theta(
        6,
        &[
            [0, 1],
            [1, 4],
            [0, 2],
            [2, 4],
            [0, 3],
            [3, 5],
            [5, 4],
        ],
        vec![
            vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)],
            vec![EdgeId(0), EdgeId(1), EdgeId(4), EdgeId(5), EdgeId(6)],
            vec![EdgeId(2), EdgeId(3), EdgeId(4), EdgeId(5), EdgeId(6)],
        ],
    )]
    fn test_relevant_cycles(
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
        #[case] expected: Vec<Vec<EdgeId>>,
    ) {
        assert_eq!(relevant_cycles(&Graph::new(node_count, edges)), expected);
    }

    #[rstest]
    #[case::path(3, &[[0, 1], [1, 2]], vec![])]
    #[case::disconnected_loops(
        2,
        &[[0, 0], [1, 1]],
        vec![
            vec![vec![EdgeId(0)]],
            vec![vec![EdgeId(1)]],
        ],
    )]
    #[case::parallel_triple(
        2,
        &[[0, 1], [0, 1], [0, 1]],
        vec![
            vec![vec![EdgeId(0), EdgeId(1)]],
            vec![vec![EdgeId(0), EdgeId(2)]],
            vec![vec![EdgeId(1), EdgeId(2)]],
        ],
    )]
    #[case::square_with_diagonal(
        4,
        &[[0, 1], [1, 2], [2, 3], [0, 3], [0, 2]],
        vec![
            vec![vec![EdgeId(0), EdgeId(1), EdgeId(4)]],
            vec![vec![EdgeId(2), EdgeId(3), EdgeId(4)]],
        ],
    )]
    #[case::unequal_theta(
        6,
        &[
            [0, 1],
            [1, 4],
            [0, 2],
            [2, 4],
            [0, 3],
            [3, 5],
            [5, 4],
        ],
        vec![
            vec![vec![EdgeId(0), EdgeId(1), EdgeId(2), EdgeId(3)]],
            vec![
                vec![EdgeId(0), EdgeId(1), EdgeId(4), EdgeId(5), EdgeId(6)],
                vec![EdgeId(2), EdgeId(3), EdgeId(4), EdgeId(5), EdgeId(6)],
            ],
        ],
    )]
    fn test_unique_ring_families(
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
        #[case] expected: Vec<Vec<Vec<EdgeId>>>,
    ) {
        assert_eq!(
            unique_ring_families(&Graph::new(node_count, edges)),
            expected
        );
    }
}
