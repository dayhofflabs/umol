//! Binary cycle-space operations.

use bitvec::prelude::*;

use super::Cycle;
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
        while let Some(pivot) = vector.highest_edge() {
            let Some(row) = &self.rows[pivot] else {
                self.rows[pivot] = Some(vector);
                self.rank += 1;
                return true;
            };
            vector.symmetric_difference_assign(row);
        }
        false
    }

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
