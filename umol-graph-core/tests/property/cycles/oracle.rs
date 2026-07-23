//! Definition-level exhaustive cycle operations for small test graphs.

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
    cycles.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));
    cycles
}

pub(super) fn cycle_space_rank(graph: &Graph) -> usize {
    graph.edge_count() + component_count(graph) - graph.node_count()
}

pub(super) fn are_linearly_independent(cycles: &[Vec<EdgeId>], edge_count: usize) -> bool {
    let mut basis: Vec<Option<Vec<bool>>> = vec![None; edge_count];

    for cycle in cycles {
        let mut row = vec![false; edge_count];
        for &edge in cycle {
            row[edge.index()] ^= true;
        }

        let mut inserted = false;
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
                inserted = true;
                break;
            }
        }
        if !inserted {
            return false;
        }
    }

    true
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
}
