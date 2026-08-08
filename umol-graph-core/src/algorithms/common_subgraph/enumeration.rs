//! Complete and maximal common-subgraph enumeration.
//!
//! Complete enumeration supports ordered clique backtracking over the modular
//! product and direct backtracking over partial node mappings. Maximal
//! enumeration uses Bron--Kerbosch with pivoting. See
//! [Bron and Kerbosch (1973)](https://doi.org/10.1145/362342.362367).

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
        match alg {
            CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking => {
                let (pairs, neighbors) =
                    self.modular_product(other, node_match, edge_match, embedding);
                let mut cliques = Vec::new();
                all_cliques(
                    &neighbors,
                    &mut Vec::new(),
                    bitvec![1; pairs.len()],
                    &mut cliques,
                );
                subgraphs_from_cliques(cliques, &pairs, self, other)
            }
            CommonSubgraphEnumerationAlgorithm::DirectBacktracking => {
                direct_backtracking(self, other, node_match, edge_match, embedding)
            }
        }
    }

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
        let (pairs, neighbors) = self.modular_product(other, node_match, edge_match, embedding);
        let count = pairs.len();
        let mut cliques = Vec::new();
        match alg {
            MaximalCommonSubgraphAlgorithm::BronKerbosch => bron_kerbosch(
                &neighbors,
                &mut Vec::new(),
                bitvec![1; count],
                bitvec![0; count],
                &mut cliques,
            ),
        }
        subgraphs_from_cliques(cliques, &pairs, self, other)
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

fn subgraphs_from_cliques(
    cliques: Vec<Vec<usize>>,
    pairs: &[(NodeId, NodeId)],
    a: &Graph,
    b: &Graph,
) -> Vec<GraphCorrespondence> {
    let mut subgraphs: Vec<GraphCorrespondence> = cliques
        .into_iter()
        .map(|clique| {
            let mapping: Vec<(NodeId, NodeId)> = clique.iter().map(|&i| pairs[i]).collect();
            subgraph_from_mapping(mapping, a, b)
        })
        .collect();
    subgraphs.sort_by(|x, y| x.nodes().matched_pairs().cmp(y.nodes().matched_pairs()));
    subgraphs.dedup();
    subgraphs
}

fn subgraph_from_mapping(
    mapping: Vec<(NodeId, NodeId)>,
    left: &Graph,
    right: &Graph,
) -> GraphCorrespondence {
    let mut edges = Vec::new();
    for x in 0..mapping.len() {
        for y in (x + 1)..mapping.len() {
            let (left_a, right_a) = mapping[x];
            let (left_b, right_b) = mapping[y];
            if let (Some(left_edge), Some(right_edge)) = (
                left.find_edge(left_a, left_b),
                right.find_edge(right_a, right_b),
            ) {
                edges.push((left_edge, right_edge));
            }
        }
    }
    GraphCorrespondence::new(
        Correspondence::new(mapping, left.node_count(), right.node_count())
            .expect("common-subgraph node pairs form a valid correspondence"),
        Correspondence::new(edges, left.edge_count(), right.edge_count())
            .expect("common-subgraph edge pairs form a valid correspondence"),
    )
}

fn direct_backtracking<N, E>(
    left: &Graph,
    right: &Graph,
    node_match: &mut N,
    edge_match: &mut E,
    embedding: EmbeddingKind,
) -> Vec<GraphCorrespondence>
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
        subgraphs: Vec::new(),
    };
    state.search(0);
    state
        .subgraphs
        .sort_by(|a, b| a.nodes().matched_pairs().cmp(b.nodes().matched_pairs()));
    state.subgraphs
}

struct DirectEnumerationState<'g, 'm, E> {
    left: &'g Graph,
    right: &'g Graph,
    edge_match: &'m mut E,
    embedding: EmbeddingKind,
    candidates: Vec<Vec<NodeId>>,
    used_right: BitVec,
    matched_pairs: Vec<(NodeId, NodeId)>,
    subgraphs: Vec<GraphCorrespondence>,
}

impl<E> DirectEnumerationState<'_, '_, E>
where
    E: FnMut(EdgeId, EdgeId) -> bool,
{
    fn search(&mut self, left_index: usize) {
        if left_index == self.left.node_count() {
            self.subgraphs.push(subgraph_from_mapping(
                self.matched_pairs.clone(),
                self.left,
                self.right,
            ));
            return;
        }

        let left_node = NodeId::from(left_index);
        for candidate_index in 0..self.candidates[left_index].len() {
            let right_node = self.candidates[left_index][candidate_index];
            if self.used_right[right_node.index()] || !self.compatible(left_node, right_node) {
                continue;
            }
            self.used_right.set(right_node.index(), true);
            self.matched_pairs.push((left_node, right_node));
            self.search(left_index + 1);
            self.matched_pairs.pop();
            self.used_right.set(right_node.index(), false);
        }
        self.search(left_index + 1);
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

fn all_cliques(
    neighbors: &[BitVec],
    clique: &mut Vec<usize>,
    candidates: BitVec,
    out: &mut Vec<Vec<usize>>,
) {
    out.push(clique.clone());
    let start = clique.last().map_or(0, |&v| v + 1);
    let extend: Vec<usize> = candidates.iter_ones().filter(|&v| v >= start).collect();
    for v in extend {
        let next = candidates.clone() & neighbors[v].clone();
        clique.push(v);
        all_cliques(neighbors, clique, next, out);
        clique.pop();
    }
}

fn bron_kerbosch(
    neighbors: &[BitVec],
    clique: &mut Vec<usize>,
    mut candidates: BitVec,
    mut excluded: BitVec,
    out: &mut Vec<Vec<usize>>,
) {
    if candidates.not_any() && excluded.not_any() {
        out.push(clique.clone());
        return;
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
        bron_kerbosch(
            neighbors,
            clique,
            candidates.clone() & nv.clone(),
            excluded.clone() & nv,
            out,
        );
        clique.pop();
        candidates.set(v, false);
        excluded.set(v, true);
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

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
    fn test_direct_backtracking(
        #[case] left: Graph,
        #[case] right: Graph,
        #[case] left_node_labels: Vec<u8>,
        #[case] right_node_labels: Vec<u8>,
        #[case] left_edge_labels: Vec<u8>,
        #[case] right_edge_labels: Vec<u8>,
        #[case] expected: Vec<GraphCorrespondence>,
    ) {
        assert_eq!(
            direct_backtracking(
                &left,
                &right,
                &mut |left_node, right_node| {
                    left_node_labels[left_node.index()] == right_node_labels[right_node.index()]
                },
                &mut |left_edge, right_edge| {
                    left_edge_labels[left_edge.index()] == right_edge_labels[right_edge.index()]
                },
                EmbeddingKind::Induced,
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
}
