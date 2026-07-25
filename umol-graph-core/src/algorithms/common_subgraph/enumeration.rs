//! Complete and maximal common-subgraph enumeration.

use bitvec::prelude::*;

use crate::correspondence::{Correspondence, GraphCorrespondence};
use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CommonSubgraphEnumerationAlgorithm {
    Backtracking,
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
        let (pairs, neighbors) = self.modular_product(other, node_match, edge_match, embedding);
        let mut cliques = Vec::new();
        match alg {
            CommonSubgraphEnumerationAlgorithm::Backtracking => {
                all_cliques(
                    &neighbors,
                    &mut Vec::new(),
                    bitvec![1; pairs.len()],
                    &mut cliques,
                );
            }
        }
        subgraphs_from_cliques(cliques, &pairs, self, other)
    }

    pub fn maximal_common_subgraphs<N, E>(
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
            let mut edges: Vec<(EdgeId, EdgeId)> = Vec::new();
            for x in 0..clique.len() {
                for y in (x + 1)..clique.len() {
                    let (a1, b1) = pairs[clique[x]];
                    let (a2, b2) = pairs[clique[y]];
                    if let (Some(ea), Some(eb)) = (a.find_edge(a1, a2), b.find_edge(b1, b2)) {
                        edges.push((ea, eb));
                    }
                }
            }
            let mapping: Vec<(NodeId, NodeId)> = clique.iter().map(|&i| pairs[i]).collect();
            GraphCorrespondence::new(
                Correspondence::new(mapping, a.node_count(), b.node_count()),
                Correspondence::new(edges, a.edge_count(), b.edge_count()),
            )
        })
        .collect();
    subgraphs.sort_by(|x, y| x.nodes().mates().cmp(y.nodes().mates()));
    subgraphs.dedup();
    subgraphs
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
