//! Complete and maximal common-subgraph enumeration.
//!
//! Complete enumeration supports ordered clique backtracking over the modular
//! product and direct backtracking over partial node mappings. Maximal
//! enumeration uses Bron--Kerbosch with pivoting. See
//! [Bron and Kerbosch (1973)](https://doi.org/10.1145/362342.362367).

use std::ops::ControlFlow;

use bitvec::prelude::*;

use crate::correspondence::{Correspondence, GraphCorrespondence};
use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CommonSubgraphEnumerationAlgorithm {
    /// Enumerate every clique of the modular product.
    ModularProductBacktracking,
    /// Enumerate partial injective node mappings directly.
    DirectBacktracking,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaximalCommonSubgraphAlgorithm {
    BronKerbosch,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmbeddingKind {
    Induced,
    Monomorphism,
}

impl Graph {
    /// Visits every common subgraph of `self` and `other` — the empty one
    /// included — as a left↔right node-pair slice until traversal completes or
    /// the visitor returns [`ControlFlow::Break`]. The slice borrows search
    /// state, so it is only valid for the duration of the call. Traversal is
    /// deterministic for a fixed graph representation, but its order is not a
    /// canonical ordering contract.
    pub fn visit_common_subgraphs<B, N, E, F>(
        &self,
        other: &Graph,
        node_match: &mut N,
        edge_match: &mut E,
        embedding: EmbeddingKind,
        alg: CommonSubgraphEnumerationAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
        F: FnMut(&[(NodeId, NodeId)]) -> ControlFlow<B>,
    {
        match alg {
            CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking => {
                let (pairs, neighbors) =
                    self.modular_product(other, node_match, edge_match, embedding);
                all_cliques(
                    &neighbors,
                    &pairs,
                    &mut Vec::new(),
                    bitvec![1; pairs.len()],
                    &mut Vec::new(),
                    &mut visitor,
                )
            }
            CommonSubgraphEnumerationAlgorithm::DirectBacktracking => {
                direct_backtracking(self, other, node_match, edge_match, embedding, &mut visitor)
            }
        }
    }

    /// Collects every common subgraph as a [`GraphCorrespondence`] by
    /// collecting [`Graph::visit_common_subgraphs`], sorted by node pairs.
    pub fn enumerate_common_subgraphs<N, E>(
        &self,
        other: &Graph,
        node_match: &mut N,
        edge_match: &mut E,
        embedding: EmbeddingKind,
        alg: CommonSubgraphEnumerationAlgorithm,
    ) -> Vec<GraphCorrespondence>
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
    {
        let mut subgraphs = Vec::new();
        let _: ControlFlow<()> =
            self.visit_common_subgraphs(other, node_match, edge_match, embedding, alg, |pairs| {
                let nodes =
                    Correspondence::new(pairs.to_vec(), self.node_count(), other.node_count())
                        .expect("common-subgraph node pairs form a valid correspondence");
                subgraphs.push(
                    GraphCorrespondence::induce(self, other, nodes)
                        .expect("a common-subgraph pairing induces a unique graph correspondence"),
                );
                ControlFlow::Continue(())
            });
        subgraphs.sort_by(|x, y| x.nodes().matched_pairs().cmp(y.nodes().matched_pairs()));
        subgraphs.dedup();
        subgraphs
    }

    /// Visits every maximal common subgraph of `self` and `other` as a
    /// left↔right node-pair slice until traversal completes or the visitor
    /// returns [`ControlFlow::Break`]. The slice borrows search state, so it
    /// is only valid for the duration of the call. Traversal is deterministic
    /// for a fixed graph representation, but its order is not a canonical
    /// ordering contract.
    pub fn visit_maximal_common_subgraphs<B, N, E, F>(
        &self,
        other: &Graph,
        node_match: &mut N,
        edge_match: &mut E,
        embedding: EmbeddingKind,
        alg: MaximalCommonSubgraphAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
        F: FnMut(&[(NodeId, NodeId)]) -> ControlFlow<B>,
    {
        let (pairs, neighbors) = self.modular_product(other, node_match, edge_match, embedding);
        let count = pairs.len();
        match alg {
            MaximalCommonSubgraphAlgorithm::BronKerbosch => bron_kerbosch(
                &neighbors,
                &pairs,
                &mut Vec::new(),
                bitvec![1; count],
                bitvec![0; count],
                &mut Vec::new(),
                &mut visitor,
            ),
        }
    }

    /// Collects every maximal common subgraph as a [`GraphCorrespondence`] by
    /// collecting [`Graph::visit_maximal_common_subgraphs`], sorted by node
    /// pairs.
    pub fn enumerate_maximal_common_subgraphs<N, E>(
        &self,
        other: &Graph,
        node_match: &mut N,
        edge_match: &mut E,
        embedding: EmbeddingKind,
        alg: MaximalCommonSubgraphAlgorithm,
    ) -> Vec<GraphCorrespondence>
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
    {
        let mut subgraphs = Vec::new();
        let _: ControlFlow<()> = self.visit_maximal_common_subgraphs(
            other,
            node_match,
            edge_match,
            embedding,
            alg,
            |pairs| {
                let nodes =
                    Correspondence::new(pairs.to_vec(), self.node_count(), other.node_count())
                        .expect("common-subgraph node pairs form a valid correspondence");
                subgraphs.push(
                    GraphCorrespondence::induce(self, other, nodes)
                        .expect("a common-subgraph pairing induces a unique graph correspondence"),
                );
                ControlFlow::Continue(())
            },
        );
        subgraphs.sort_by(|x, y| x.nodes().matched_pairs().cmp(y.nodes().matched_pairs()));
        subgraphs.dedup();
        subgraphs
    }

    pub(super) fn modular_product<N, E>(
        &self,
        other: &Graph,
        node_match: &mut N,
        edge_match: &mut E,
        embedding: EmbeddingKind,
    ) -> (Vec<(NodeId, NodeId)>, Vec<BitVec>)
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
    {
        let mut pairs: Vec<(NodeId, NodeId)> = Vec::new();
        for a in self.node_ids() {
            for b in other.node_ids() {
                if node_match(a, b) {
                    pairs.push((a, b));
                }
            }
        }
        let count = pairs.len();

        let mut neighbors = vec![bitvec![0; count]; count];
        for i in 0..count {
            let (a1, b1) = pairs[i];
            for j in (i + 1)..count {
                let (a2, b2) = pairs[j];
                if a1 == a2 || b1 == b2 {
                    continue;
                }
                let agree = match (self.find_edge(a1, a2), other.find_edge(b1, b2)) {
                    (Some(ea), Some(eb)) => edge_match(ea, eb),
                    (None, None) => true,
                    _ => matches!(embedding, EmbeddingKind::Monomorphism),
                };
                if agree {
                    neighbors[i].set(j, true);
                    neighbors[j].set(i, true);
                }
            }
        }
        (pairs, neighbors)
    }
}

fn direct_backtracking<B, N, E>(
    left: &Graph,
    right: &Graph,
    node_match: &mut N,
    edge_match: &mut E,
    embedding: EmbeddingKind,
    emit: &mut impl FnMut(&[(NodeId, NodeId)]) -> ControlFlow<B>,
) -> ControlFlow<B>
where
    N: FnMut(NodeId, NodeId) -> bool,
    E: FnMut(EdgeId, EdgeId) -> bool,
{
    let candidates = left
        .node_ids()
        .map(|left_node| {
            right
                .node_ids()
                .filter(|&right_node| node_match(left_node, right_node))
                .collect()
        })
        .collect();
    let mut state = DirectEnumerationState {
        left,
        right,
        edge_match,
        embedding,
        candidates,
        used_right: bitvec![0; right.node_count()],
        matched_pairs: Vec::new(),
    };
    state.search(0, emit)
}

struct DirectEnumerationState<'g, 'm, E> {
    left: &'g Graph,
    right: &'g Graph,
    edge_match: &'m mut E,
    embedding: EmbeddingKind,
    candidates: Vec<Vec<NodeId>>,
    used_right: BitVec,
    // Doubles as the emitted pair slice at each leaf.
    matched_pairs: Vec<(NodeId, NodeId)>,
}

impl<E> DirectEnumerationState<'_, '_, E>
where
    E: FnMut(EdgeId, EdgeId) -> bool,
{
    fn search<B>(
        &mut self,
        left_index: usize,
        emit: &mut impl FnMut(&[(NodeId, NodeId)]) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        if left_index == self.left.node_count() {
            return emit(&self.matched_pairs);
        }

        let left_node = NodeId::from(left_index);
        for candidate_index in 0..self.candidates[left_index].len() {
            let right_node = self.candidates[left_index][candidate_index];
            if self.used_right[right_node.index()] || !self.compatible(left_node, right_node) {
                continue;
            }
            self.used_right.set(right_node.index(), true);
            self.matched_pairs.push((left_node, right_node));
            let result = self.search(left_index + 1, emit);
            self.matched_pairs.pop();
            self.used_right.set(right_node.index(), false);
            if let ControlFlow::Break(value) = result {
                return ControlFlow::Break(value);
            }
        }
        self.search(left_index + 1, emit)
    }

    fn compatible(&mut self, left_node: NodeId, right_node: NodeId) -> bool {
        for &(mapped_left, mapped_right) in &self.matched_pairs {
            match (
                self.left.find_edge(left_node, mapped_left),
                self.right.find_edge(right_node, mapped_right),
            ) {
                (Some(left_edge), Some(right_edge)) => {
                    if !(self.edge_match)(left_edge, right_edge) {
                        return false;
                    }
                }
                (None, None) => {}
                _ if self.embedding == EmbeddingKind::Monomorphism => {}
                _ => return false,
            }
        }
        true
    }
}

fn all_cliques<B>(
    neighbors: &[BitVec],
    pairs: &[(NodeId, NodeId)],
    clique: &mut Vec<usize>,
    candidates: BitVec,
    scratch: &mut Vec<(NodeId, NodeId)>,
    emit: &mut impl FnMut(&[(NodeId, NodeId)]) -> ControlFlow<B>,
) -> ControlFlow<B> {
    scratch.clear();
    scratch.extend(clique.iter().map(|&index| pairs[index]));
    if let ControlFlow::Break(value) = emit(scratch) {
        return ControlFlow::Break(value);
    }
    let start = clique.last().map_or(0, |&v| v + 1);
    let extend: Vec<usize> = candidates.iter_ones().filter(|&v| v >= start).collect();
    for v in extend {
        let next = candidates.clone() & neighbors[v].clone();
        clique.push(v);
        let result = all_cliques(neighbors, pairs, clique, next, scratch, emit);
        clique.pop();
        if let ControlFlow::Break(value) = result {
            return ControlFlow::Break(value);
        }
    }
    ControlFlow::Continue(())
}

fn bron_kerbosch<B>(
    neighbors: &[BitVec],
    pairs: &[(NodeId, NodeId)],
    clique: &mut Vec<usize>,
    mut candidates: BitVec,
    mut excluded: BitVec,
    scratch: &mut Vec<(NodeId, NodeId)>,
    emit: &mut impl FnMut(&[(NodeId, NodeId)]) -> ControlFlow<B>,
) -> ControlFlow<B> {
    if candidates.not_any() && excluded.not_any() {
        scratch.clear();
        scratch.extend(clique.iter().map(|&index| pairs[index]));
        return emit(scratch);
    }
    let pivot = candidates
        .iter_ones()
        .chain(excluded.iter_ones())
        .max_by_key(|&u| (candidates.clone() & neighbors[u].clone()).count_ones())
        .expect("candidates ∪ excluded is non-empty here");
    let expand: Vec<usize> = (candidates.clone() & !neighbors[pivot].clone())
        .iter_ones()
        .collect();
    for v in expand {
        let nv = neighbors[v].clone();
        clique.push(v);
        let result = bron_kerbosch(
            neighbors,
            pairs,
            clique,
            candidates.clone() & nv.clone(),
            excluded.clone() & nv,
            scratch,
            emit,
        );
        clique.pop();
        if let ControlFlow::Break(value) = result {
            return ControlFlow::Break(value);
        }
        candidates.set(v, false);
        excluded.set(v, true);
    }
    ControlFlow::Continue(())
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::correspondence::Correspondence;

    #[rstest]
    #[case::modular_product(CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking)]
    #[case::direct(CommonSubgraphEnumerationAlgorithm::DirectBacktracking)]
    fn test_graph_visit_common_subgraphs(#[case] alg: CommonSubgraphEnumerationAlgorithm) {
        let left = Graph::new(2, &[[0, 1]]);
        let right = Graph::new(2, &[[0, 1]]);
        let node_labels = [0u8, 1];

        let mut emissions: Vec<Vec<(NodeId, NodeId)>> = Vec::new();
        let flow: ControlFlow<()> = left.visit_common_subgraphs(
            &right,
            &mut |left_node, right_node| {
                node_labels[left_node.index()] == node_labels[right_node.index()]
            },
            &mut |_, _| true,
            EmbeddingKind::Induced,
            alg,
            |pairs| {
                let mut pairs = pairs.to_vec();
                pairs.sort_unstable();
                emissions.push(pairs);
                ControlFlow::Continue(())
            },
        );
        assert_eq!(flow, ControlFlow::Continue(()));
        emissions.sort();
        assert_eq!(
            emissions,
            vec![
                vec![],
                vec![(NodeId(0), NodeId(0))],
                vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))],
                vec![(NodeId(1), NodeId(1))],
            ],
        );
    }

    #[rstest]
    #[case::modular_product(CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking)]
    #[case::direct(CommonSubgraphEnumerationAlgorithm::DirectBacktracking)]
    fn test_graph_visit_common_subgraphs_termination(
        #[case] alg: CommonSubgraphEnumerationAlgorithm,
    ) {
        let left = Graph::new(2, &[[0, 1]]);
        let right = Graph::new(2, &[[0, 1]]);

        let result = left.visit_common_subgraphs(
            &right,
            &mut |_, _| true,
            &mut |_, _| true,
            EmbeddingKind::Induced,
            alg,
            |pairs| {
                let mut pairs = pairs.to_vec();
                pairs.sort_unstable();
                ControlFlow::Break(pairs)
            },
        );
        let ControlFlow::Break(first) = result else {
            panic!("expected Break on first emission");
        };
        let expected = [
            vec![],
            vec![(NodeId(0), NodeId(0))],
            vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))],
            vec![(NodeId(0), NodeId(1))],
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
            vec![(NodeId(1), NodeId(0))],
            vec![(NodeId(1), NodeId(1))],
        ];
        assert!(expected.contains(&first), "invalid emission {first:?}");
    }

    #[rstest]
    #[case::empty(
        Graph::new(0, &[]),
        Graph::new(0, &[]),
        vec![],
        vec![],
        vec![],
        vec![],
        vec![GraphCorrespondence::new(
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        )],
    )]
    #[case::isolated(
        Graph::new(1, &[]),
        Graph::new(1, &[]),
        vec![0],
        vec![0],
        vec![],
        vec![],
        vec![
            GraphCorrespondence::new(
                Correspondence::new(vec![], 1, 1).expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            ),
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(0))], 1, 1).expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            ),
        ],
    )]
    #[case::incompatible_node(
        Graph::new(1, &[]),
        Graph::new(1, &[]),
        vec![0],
        vec![1],
        vec![],
        vec![],
        vec![GraphCorrespondence::new(
            Correspondence::new(vec![], 1, 1).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        )],
    )]
    fn test_graph_enumerate_common_subgraphs(
        #[case] left: Graph,
        #[case] right: Graph,
        #[case] left_node_labels: Vec<u8>,
        #[case] right_node_labels: Vec<u8>,
        #[case] left_edge_labels: Vec<u8>,
        #[case] right_edge_labels: Vec<u8>,
        #[case] expected: Vec<GraphCorrespondence>,
    ) {
        assert_eq!(
            left.enumerate_common_subgraphs(
                &right,
                &mut |left_node, right_node| {
                    left_node_labels[left_node.index()] == right_node_labels[right_node.index()]
                },
                &mut |left_edge, right_edge| {
                    left_edge_labels[left_edge.index()] == right_edge_labels[right_edge.index()]
                },
                EmbeddingKind::Induced,
                CommonSubgraphEnumerationAlgorithm::DirectBacktracking,
            ),
            expected,
        );
    }

    #[rstest]
    #[case::injective_isolated(
        Graph::new(2, &[]),
        Graph::new(1, &[]),
        vec![0, 0],
        vec![0],
        vec![],
        vec![],
        EmbeddingKind::Induced,
    )]
    #[case::incompatible_edge(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        vec![0, 0],
        vec![0, 0],
        vec![0],
        vec![1],
        EmbeddingKind::Induced,
    )]
    #[case::disconnected(
        Graph::new(4, &[[0, 1], [2, 3]]),
        Graph::new(4, &[[0, 1], [2, 3]]),
        vec![0, 1, 0, 1],
        vec![0, 1, 0, 1],
        vec![0, 1],
        vec![0, 1],
        EmbeddingKind::Induced,
    )]
    #[case::edge_presence_induced(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[]),
        vec![0, 0],
        vec![0, 0],
        vec![0],
        vec![],
        EmbeddingKind::Induced,
    )]
    #[case::edge_presence_monomorphism(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[]),
        vec![0, 0],
        vec![0, 0],
        vec![0],
        vec![],
        EmbeddingKind::Monomorphism,
    )]
    fn test_graph_enumerate_common_subgraphs_equivalence(
        #[case] left: Graph,
        #[case] right: Graph,
        #[case] left_node_labels: Vec<u8>,
        #[case] right_node_labels: Vec<u8>,
        #[case] left_edge_labels: Vec<u8>,
        #[case] right_edge_labels: Vec<u8>,
        #[case] embedding: EmbeddingKind,
    ) {
        let direct = left.enumerate_common_subgraphs(
            &right,
            &mut |left_node, right_node| {
                left_node_labels[left_node.index()] == right_node_labels[right_node.index()]
            },
            &mut |left_edge, right_edge| {
                left_edge_labels[left_edge.index()] == right_edge_labels[right_edge.index()]
            },
            embedding,
            CommonSubgraphEnumerationAlgorithm::DirectBacktracking,
        );
        let modular_product = left.enumerate_common_subgraphs(
            &right,
            &mut |left_node, right_node| {
                left_node_labels[left_node.index()] == right_node_labels[right_node.index()]
            },
            &mut |left_edge, right_edge| {
                left_edge_labels[left_edge.index()] == right_edge_labels[right_edge.index()]
            },
            embedding,
            CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        );

        assert_eq!(direct, modular_product);
    }

    #[rstest]
    #[case::single_edge(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        vec![0, 0],
        vec![0, 0],
        vec![
            vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))],
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
        ],
    )]
    #[case::labeled_path(
        Graph::new(3, &[[0, 1], [1, 2]]),
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![0, 1, 2],
        vec![0, 1, 2],
        vec![vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1)), (NodeId(2), NodeId(2))]],
    )]
    fn test_graph_visit_maximal_common_subgraphs(
        #[case] left: Graph,
        #[case] right: Graph,
        #[case] left_node_labels: Vec<u8>,
        #[case] right_node_labels: Vec<u8>,
        #[case] expected: Vec<Vec<(NodeId, NodeId)>>,
    ) {
        let mut emissions: Vec<Vec<(NodeId, NodeId)>> = Vec::new();
        let flow: ControlFlow<()> = left.visit_maximal_common_subgraphs(
            &right,
            &mut |left_node, right_node| {
                left_node_labels[left_node.index()] == right_node_labels[right_node.index()]
            },
            &mut |_, _| true,
            EmbeddingKind::Induced,
            MaximalCommonSubgraphAlgorithm::BronKerbosch,
            |pairs| {
                let mut pairs = pairs.to_vec();
                pairs.sort_unstable();
                emissions.push(pairs);
                ControlFlow::Continue(())
            },
        );
        assert_eq!(flow, ControlFlow::Continue(()));
        emissions.sort();
        assert_eq!(emissions, expected);
    }

    #[rstest]
    #[case::bron_kerbosch(MaximalCommonSubgraphAlgorithm::BronKerbosch)]
    fn test_graph_visit_maximal_common_subgraphs_termination(
        #[case] alg: MaximalCommonSubgraphAlgorithm,
    ) {
        let left = Graph::new(2, &[[0, 1]]);
        let right = Graph::new(2, &[[0, 1]]);

        let result = left.visit_maximal_common_subgraphs(
            &right,
            &mut |_, _| true,
            &mut |_, _| true,
            EmbeddingKind::Induced,
            alg,
            |pairs| {
                let mut pairs = pairs.to_vec();
                pairs.sort_unstable();
                ControlFlow::Break(pairs)
            },
        );
        let ControlFlow::Break(first) = result else {
            panic!("expected Break on first emission");
        };
        let expected = [
            vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))],
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
        ];
        assert!(expected.contains(&first), "invalid emission {first:?}");
    }
}
