//! Binary cycle-space operations.

use std::cmp::Ordering;
use std::collections::HashSet;

use bitvec::prelude::*;

use super::relevant::ShortestPathDag;
use super::{Cycle, MinimumCycleBasis};
use crate::algorithms::connected::ConnectedComponentsAlgorithm;
use crate::graph::{EdgeId, Graph};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct EdgeVector {
    bits: BitVec,
}

impl EdgeVector {
    pub(super) fn from_edges(edge_count: usize, edges: impl IntoIterator<Item = EdgeId>) -> Self {
        let mut bits = bitvec![0; edge_count];
        for edge in edges {
            assert!(
                edge.index() < edge_count,
                "edge vector contains an out-of-bounds edge"
            );
            bits.set(edge.index(), true);
        }
        Self { bits }
    }

    pub(super) fn from_cycle(edge_count: usize, cycle: &Cycle) -> Self {
        Self::from_edges(edge_count, cycle.edges().iter().copied())
    }

    pub(super) fn symmetric_difference_assign(&mut self, other: &Self) {
        assert_eq!(
            self.bits.len(),
            other.bits.len(),
            "edge vectors must have equal dimensions"
        );
        self.bits ^= &other.bits;
    }

    fn highest_edge(&self) -> Option<usize> {
        self.bits.iter_ones().next_back()
    }
}

#[derive(Clone, Debug)]
pub(super) struct CycleVectorBasis {
    rows: Vec<Option<EdgeVector>>,
    rank: usize,
}

impl CycleVectorBasis {
    pub(super) fn new(edge_count: usize) -> Self {
        Self {
            rows: vec![None; edge_count],
            rank: 0,
        }
    }

    pub(super) fn insert(&mut self, mut vector: EdgeVector) -> bool {
        self.reduce(&mut vector);
        let Some(pivot) = vector.highest_edge() else {
            return false;
        };
        self.rows[pivot] = Some(vector);
        self.rank += 1;
        true
    }

    pub(super) fn is_independent(&self, mut vector: EdgeVector) -> bool {
        self.reduce(&mut vector);
        vector.highest_edge().is_some()
    }

    pub(super) fn reduced(&self, mut vector: EdgeVector) -> EdgeVector {
        self.reduce(&mut vector);
        vector
    }

    fn reduce(&self, vector: &mut EdgeVector) {
        while let Some(pivot) = vector.highest_edge() {
            let Some(row) = &self.rows[pivot] else {
                return;
            };
            vector.symmetric_difference_assign(row);
        }
    }

    #[cfg(test)]
    pub(super) fn rank(&self) -> usize {
        self.rank
    }
}

pub(super) fn cycle_space_rank(graph: &Graph) -> usize {
    graph.edge_count()
        + graph
            .connected_components(ConnectedComponentsAlgorithm::Bfs)
            .len()
        - graph.node_count()
}

pub(super) fn minimum_cycle_basis_horton(source: &Graph) -> MinimumCycleBasis {
    let mut loops = Vec::new();
    let mut loopless_edges = Vec::new();
    let mut edge_sources = Vec::new();
    let mut endpoint_pairs = HashSet::new();
    let mut has_parallel_edges = false;

    for edge in source.edge_ids() {
        let [first, second] = source.edge_endpoints(edge);
        if first == second {
            loops.push(Cycle::normalized(source, vec![first], vec![edge]));
            continue;
        }
        has_parallel_edges |= !endpoint_pairs.insert([first, second]);
        loopless_edges.push([first.0, second.0]);
        edge_sources.push(edge);
    }

    let mut cycles = if loops.is_empty() && !has_parallel_edges {
        minimum_cycle_basis_simple(source)
    } else {
        let loopless = Graph::new(source.node_count(), &loopless_edges);
        if has_parallel_edges {
            let subdivision = loopless.subdivide_edges();
            minimum_cycle_basis_simple(subdivision.graph())
                .into_iter()
                .map(|cycle| cycle.map_subdivision(source, &subdivision, &edge_sources))
                .collect()
        } else {
            minimum_cycle_basis_simple(&loopless)
                .into_iter()
                .map(|cycle| cycle.map_edges(source, &edge_sources))
                .collect()
        }
    };
    cycles.extend(loops);
    cycles.sort_by(compare_cycles);

    assert_eq!(
        cycles.len(),
        cycle_space_rank(source),
        "Horton candidates must span the source cycle space"
    );
    let total_length = cycles.iter().map(Cycle::length).sum();
    MinimumCycleBasis {
        cycles,
        total_length,
    }
}

fn minimum_cycle_basis_simple(graph: &Graph) -> Vec<Cycle> {
    let target_rank = cycle_space_rank(graph);
    if target_rank == 0 {
        return Vec::new();
    }

    let mut candidates = HashSet::new();
    for root in graph.node_ids() {
        let paths = ShortestPathDag::new(graph, root);
        for edge in graph.edge_ids() {
            let [first, second] = graph.edge_endpoints(edge);
            let (Some(first_path), Some(second_path)) =
                (paths.path_to(first), paths.path_to(second))
            else {
                continue;
            };
            if let Some(cycle) = first_path.cycle_with(graph, &second_path, edge) {
                candidates.insert(cycle);
            }
        }
    }

    let mut candidates = candidates.into_iter().collect::<Vec<_>>();
    candidates.sort_by(compare_cycles);
    let mut vectors = CycleVectorBasis::new(graph.edge_count());
    let mut cycles = Vec::with_capacity(target_rank);
    for candidate in candidates {
        if vectors.insert(EdgeVector::from_cycle(graph.edge_count(), &candidate)) {
            cycles.push(candidate);
            if cycles.len() == target_rank {
                break;
            }
        }
    }
    assert_eq!(
        cycles.len(),
        target_rank,
        "Horton candidates must span the cycle space"
    );
    cycles
}

fn compare_cycles(first: &Cycle, second: &Cycle) -> Ordering {
    first
        .length()
        .cmp(&second.length())
        .then_with(|| first.nodes().cmp(second.nodes()))
        .then_with(|| first.edges().cmp(second.edges()))
}

#[cfg(test)]
mod tests {
    use bitvec::prelude::{bitvec, Lsb0};
    use proptest::prelude::*;
    use rstest::rstest;

    use super::{cycle_space_rank, CycleVectorBasis, EdgeVector};
    use crate::algorithms::connected::ConnectedComponentsAlgorithm;
    use crate::algorithms::cycles::SimpleCycleEnumerationAlgorithm;
    use crate::graph::Graph;

    #[rstest]
    #[case::overlap(
        EdgeVector { bits: bitvec![1, 1, 0, 0] },
        EdgeVector { bits: bitvec![0, 1, 1, 0] },
        EdgeVector { bits: bitvec![1, 0, 1, 0] },
    )]
    #[case::identity(
        EdgeVector { bits: bitvec![1, 0, 1] },
        EdgeVector { bits: bitvec![0, 0, 0] },
        EdgeVector { bits: bitvec![1, 0, 1] },
    )]
    #[case::self_inverse(
        EdgeVector { bits: bitvec![1, 1, 1] },
        EdgeVector { bits: bitvec![1, 1, 1] },
        EdgeVector { bits: bitvec![0, 0, 0] },
    )]
    fn test_edge_vector_symmetric_difference_assign(
        #[case] mut vector: EdgeVector,
        #[case] other: EdgeVector,
        #[case] expected: EdgeVector,
    ) {
        vector.symmetric_difference_assign(&other);
        assert_eq!(vector, expected);
    }

    #[rstest]
    #[case::independent(
        4,
        vec![
            EdgeVector { bits: bitvec![1, 1, 0, 0] },
            EdgeVector { bits: bitvec![0, 1, 1, 0] },
            EdgeVector { bits: bitvec![0, 0, 1, 1] },
        ],
        vec![true, true, true],
        3,
    )]
    #[case::dependent(
        3,
        vec![
            EdgeVector { bits: bitvec![1, 1, 0] },
            EdgeVector { bits: bitvec![0, 1, 1] },
            EdgeVector { bits: bitvec![1, 0, 1] },
        ],
        vec![true, true, false],
        2,
    )]
    #[case::duplicate(
        2,
        vec![
            EdgeVector { bits: bitvec![1, 1] },
            EdgeVector { bits: bitvec![1, 1] },
        ],
        vec![true, false],
        1,
    )]
    fn test_cycle_vector_basis_insert(
        #[case] edge_count: usize,
        #[case] vectors: Vec<EdgeVector>,
        #[case] expected: Vec<bool>,
        #[case] expected_rank: usize,
    ) {
        let mut basis = CycleVectorBasis::new(edge_count);
        let actual = vectors
            .into_iter()
            .map(|vector| basis.insert(vector))
            .collect::<Vec<_>>();
        assert_eq!(actual, expected);
        assert_eq!(basis.rank(), expected_rank);
    }

    #[rstest]
    #[case::independent(
        4,
        vec![
            EdgeVector { bits: bitvec![1, 1, 0, 0] },
            EdgeVector { bits: bitvec![0, 1, 1, 0] },
        ],
        EdgeVector { bits: bitvec![0, 0, 1, 1] },
        true,
    )]
    #[case::dependent(
        3,
        vec![
            EdgeVector { bits: bitvec![1, 1, 0] },
            EdgeVector { bits: bitvec![0, 1, 1] },
        ],
        EdgeVector { bits: bitvec![1, 0, 1] },
        false,
    )]
    fn test_cycle_vector_basis_is_independent(
        #[case] edge_count: usize,
        #[case] basis_vectors: Vec<EdgeVector>,
        #[case] candidate: EdgeVector,
        #[case] expected: bool,
    ) {
        let mut basis = CycleVectorBasis::new(edge_count);
        for vector in basis_vectors {
            assert!(basis.insert(vector));
        }
        assert_eq!(basis.is_independent(candidate), expected);
    }

    #[rstest]
    #[case::empty(Graph::new(0, &[]), 0)]
    #[case::isolated(Graph::new(3, &[]), 0)]
    #[case::tree(Graph::new(4, &[[0, 1], [1, 2], [1, 3]]), 0)]
    #[case::triangle(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), 1)]
    #[case::two_components(
        Graph::new(6, &[[0, 1], [1, 2], [0, 2], [3, 4], [4, 5], [3, 5]]),
        2,
    )]
    #[case::loop_and_digon(Graph::new(2, &[[0, 0], [0, 1], [0, 1]]), 2)]
    fn test_cycle_space_rank(#[case] graph: Graph, #[case] expected: usize) {
        assert_eq!(cycle_space_rank(&graph), expected);
    }

    proptest! {
        #[test]
        fn test_cycle_vector_basis_rank(
            (node_count, endpoints) in (
                0usize..=5,
                prop::collection::vec((0..5u32, 0..5u32), 0..=7),
            ),
        ) {
            let edges = endpoints
                .into_iter()
                .filter(|&(first, second)| {
                    first < node_count as u32 && second < node_count as u32
                })
                .map(|(first, second)| [first, second])
                .collect::<Vec<_>>();
            let graph = Graph::new(
                node_count,
                &edges,
            );
            let mut basis = CycleVectorBasis::new(graph.edge_count());
            for cycle in graph.enumerate_simple_cycles(
                usize::MAX,
                SimpleCycleEnumerationAlgorithm::ReadTarjan,
            ) {
                basis.insert(EdgeVector::from_cycle(graph.edge_count(), &cycle));
            }
            let expected = graph.edge_count()
                + graph
                    .connected_components(ConnectedComponentsAlgorithm::Bfs)
                    .len()
                - graph.node_count();

            prop_assert_eq!(basis.rank(), expected);
            prop_assert_eq!(cycle_space_rank(&graph), expected);
        }
    }
}
