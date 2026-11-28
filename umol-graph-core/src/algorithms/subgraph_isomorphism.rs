//! Subgraph-isomorphism search with monomorphism semantics.
//!
//! Query edges must map to target edges; no induced-subgraph reverse check is
//! applied. Selectors provide an RDKit-compatible VF2 variant (`Vf2Rdkit`),
//! Ray--Kirsch bond search (`RayKirsch`), bit-parallel RI-DS (`BitsetRi`), and
//! a Glasgow-style constraint-programming solver (`Glasgow`). Both ordinary and
//! anchored operations return the same correspondence set for every selector.
//! See [RDKit PR #2500](https://github.com/rdkit/rdkit/pull/2500),
//! [Ray and Kirsch (1957)](https://doi.org/10.1126/science.126.3278.814),
//! Bonnici and Giugno (2017), and
//! McCreesh, Prosser, and Trimble (2020).

use std::ops::ControlFlow;

use crate::correspondence::Correspondence;
use crate::graph::Graph;
use crate::id::{EdgeId, NodeId};

mod bitset;
mod bitset_ri;
mod glasgow;
mod predicate;
mod prepared;
mod ray_kirsch;
mod symmetry;
mod vf2_rdkit;

pub(crate) use predicate::ClosureSubgraphIsomorphismPredicate;
pub use predicate::SubgraphIsomorphismPredicate;
pub use prepared::{PreparedSubgraphQuery, PreparedSubgraphTarget, SubgraphIsomorphismScratch};
pub use symmetry::{QuerySymmetryConstraints, SymmetryGeneratorError, TargetSymmetry};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubgraphIsomorphismAlgorithm {
    /// RDKit's substructure engine: static BFS query order with flat CSR parent-edge
    /// lists, candidates restricted to a mapped parent's adjacency (checking the
    /// primary parent edge directly), an explicit `deg(q) <= deg(t)` bound, and the
    /// terminal-set look-ahead disabled (RDKit PR #2500).
    Vf2Rdkit,
    /// Ray-Kirsch bond-based backtracking (Ray & Kirsch 1957; Sayle `parsmart.cpp`):
    /// compiles the query into an RI-guided `Seed`/`Extend`/`Close` step program with
    /// explicit `deg(q) <= deg(t)` bounds so each extension's adjacency and edge
    /// predicate are checked directly from CSR neighbors.
    RayKirsch,
    /// Bit-parallel RI-DS (Bonnici & Giugno 2017): rarity-and-connectivity query
    /// ordering over `u64`/`u128`/wide candidate bitsets with lazy or precomputed
    /// target adjacency, direct sparse neighbor scanning for selective single-parent
    /// steps, and bit-parallel monomorphism lookahead on cyclic queries.
    BitsetRi,
    /// Glasgow-style constraint-programming solver (McCreesh, Prosser, & Trimble
    /// 2020): compact bit-parallel forward checking, dynamic minimum-remaining-values
    /// branching, unit propagation, and bit-parallel Hall/AllDifferent pruning.
    Glasgow,
}

impl Graph {
    /// Visits every subgraph isomorphism of `query` in `self` until traversal
    /// completes or the visitor returns [`ControlFlow::Break`]. The embedding
    /// slice maps each query node index to its host node and borrows search
    /// state, so it is only valid for the duration of the call. Traversal is
    /// deterministic for a fixed graph representation, but its order is not a
    /// canonical ordering contract.
    ///
    /// # Semantic properties
    ///
    /// - Every emitted embedding is an injective mapping `0..query.node_count() -> 0..self.node_count()`
    ///   satisfying `node_match` on every query node and `edge_match` on the target edge
    ///   corresponding to every query edge (monomorphism semantics).
    /// - The set of emitted embeddings is invariant across [`SubgraphIsomorphismAlgorithm`]
    ///   variants and agrees with [`Self::enumerate_subgraph_isomorphisms`],
    ///   [`Self::has_subgraph_isomorphism`], [`Self::first_subgraph_isomorphism`], and
    ///   [`Self::count_subgraph_isomorphisms`].
    pub fn visit_subgraph_isomorphisms<B, F>(
        &self,
        query: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(&[NodeId]) -> ControlFlow<B>,
    {
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(node_match, edge_match);
        match alg {
            SubgraphIsomorphismAlgorithm::Vf2Rdkit => {
                vf2_rdkit::visit_direct(self, query, None, &mut pred, &mut visitor)
            }
            SubgraphIsomorphismAlgorithm::RayKirsch => {
                ray_kirsch::visit_direct(self, query, None, &mut pred, &mut visitor)
            }
            SubgraphIsomorphismAlgorithm::BitsetRi => {
                bitset_ri::visit_direct(self, query, None, &mut pred, &mut visitor)
            }
            SubgraphIsomorphismAlgorithm::Glasgow => {
                glasgow::visit_direct(self, query, None, &mut pred, &mut visitor)
            }
        }
    }

    /// Collects every subgraph isomorphism as a [`Correspondence`] over the
    /// host's node id space by collecting [`Graph::visit_subgraph_isomorphisms`].
    pub fn enumerate_subgraph_isomorphisms(
        &self,
        query: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Correspondence<NodeId>> {
        let mut result = Vec::new();
        let _: ControlFlow<()> =
            self.visit_subgraph_isomorphisms(query, node_match, edge_match, alg, |embedding| {
                result.push(Correspondence::from_images(embedding, self.node_count()));
                ControlFlow::Continue(())
            });
        result
    }

    /// Returns whether at least one subgraph isomorphism of `query` exists in `self`.
    pub fn has_subgraph_isomorphism(
        &self,
        query: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> bool {
        self.visit_subgraph_isomorphisms(query, node_match, edge_match, alg, |_| {
            ControlFlow::Break(())
        })
        .is_break()
    }

    /// Returns the first subgraph isomorphism of `query` in `self` as a [`Correspondence`], or
    /// `None` if no match exists.
    pub fn first_subgraph_isomorphism(
        &self,
        query: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Option<Correspondence<NodeId>> {
        match self.visit_subgraph_isomorphisms(query, node_match, edge_match, alg, |embedding| {
            ControlFlow::Break(Correspondence::from_images(embedding, self.node_count()))
        }) {
            ControlFlow::Break(correspondence) => Some(correspondence),
            ControlFlow::Continue(()) => None,
        }
    }

    /// Counts every subgraph isomorphism of `query` in `self` without materializing
    /// correspondences.
    pub fn count_subgraph_isomorphisms(
        &self,
        query: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> u64 {
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(node_match, edge_match);
        match alg {
            SubgraphIsomorphismAlgorithm::Vf2Rdkit => {
                vf2_rdkit::count_direct(self, query, None, &mut pred)
            }
            SubgraphIsomorphismAlgorithm::RayKirsch => {
                ray_kirsch::count_direct(self, query, None, &mut pred)
            }
            SubgraphIsomorphismAlgorithm::BitsetRi => {
                bitset_ri::count_direct(self, query, None, &mut pred)
            }
            SubgraphIsomorphismAlgorithm::Glasgow => {
                glasgow::count_direct(self, query, None, &mut pred)
            }
        }
    }

    /// Visits every subgraph isomorphism mapping query node `anchor.0` to host
    /// node `anchor.1`; otherwise like [`Graph::visit_subgraph_isomorphisms`].
    pub fn visit_subgraph_isomorphisms_at<B, F>(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(&[NodeId]) -> ControlFlow<B>,
    {
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(node_match, edge_match);
        match alg {
            SubgraphIsomorphismAlgorithm::Vf2Rdkit => {
                vf2_rdkit::visit_direct(self, query, Some(anchor), &mut pred, &mut visitor)
            }
            SubgraphIsomorphismAlgorithm::RayKirsch => {
                ray_kirsch::visit_direct(self, query, Some(anchor), &mut pred, &mut visitor)
            }
            SubgraphIsomorphismAlgorithm::BitsetRi => {
                bitset_ri::visit_direct(self, query, Some(anchor), &mut pred, &mut visitor)
            }
            SubgraphIsomorphismAlgorithm::Glasgow => {
                glasgow::visit_direct(self, query, Some(anchor), &mut pred, &mut visitor)
            }
        }
    }

    /// Collects every anchored subgraph isomorphism as a [`Correspondence`] by
    /// collecting [`Graph::visit_subgraph_isomorphisms_at`].
    pub fn enumerate_subgraph_isomorphisms_at(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Correspondence<NodeId>> {
        let mut result = Vec::new();
        let _: ControlFlow<()> = self.visit_subgraph_isomorphisms_at(
            query,
            anchor,
            node_match,
            edge_match,
            alg,
            |embedding| {
                result.push(Correspondence::from_images(embedding, self.node_count()));
                ControlFlow::Continue(())
            },
        );
        result
    }

    /// Returns whether at least one subgraph isomorphism mapping `anchor.0` to `anchor.1` exists
    /// in `self`.
    pub fn has_subgraph_isomorphism_at(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> bool {
        self.visit_subgraph_isomorphisms_at(query, anchor, node_match, edge_match, alg, |_| {
            ControlFlow::Break(())
        })
        .is_break()
    }

    /// Returns the first subgraph isomorphism mapping `anchor.0` to `anchor.1` as a
    /// [`Correspondence`], or `None` if no anchored match exists.
    pub fn first_subgraph_isomorphism_at(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Option<Correspondence<NodeId>> {
        match self.visit_subgraph_isomorphisms_at(
            query,
            anchor,
            node_match,
            edge_match,
            alg,
            |embedding| {
                ControlFlow::Break(Correspondence::from_images(embedding, self.node_count()))
            },
        ) {
            ControlFlow::Break(correspondence) => Some(correspondence),
            ControlFlow::Continue(()) => None,
        }
    }

    /// Counts every subgraph isomorphism mapping `anchor.0` to `anchor.1` without materializing
    /// correspondences.
    pub fn count_subgraph_isomorphisms_at(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> u64 {
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(node_match, edge_match);
        match alg {
            SubgraphIsomorphismAlgorithm::Vf2Rdkit => {
                vf2_rdkit::count_direct(self, query, Some(anchor), &mut pred)
            }
            SubgraphIsomorphismAlgorithm::RayKirsch => {
                ray_kirsch::count_direct(self, query, Some(anchor), &mut pred)
            }
            SubgraphIsomorphismAlgorithm::BitsetRi => {
                bitset_ri::count_direct(self, query, Some(anchor), &mut pred)
            }
            SubgraphIsomorphismAlgorithm::Glasgow => {
                glasgow::count_direct(self, query, Some(anchor), &mut pred)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::SubgraphIsomorphismAlgorithm::{BitsetRi, Glasgow, RayKirsch, Vf2Rdkit};
    use super::*;

    const ALGORITHMS: [SubgraphIsomorphismAlgorithm; 4] = [Vf2Rdkit, RayKirsch, BitsetRi, Glasgow];

    fn any_node(_: NodeId, _: NodeId) -> bool {
        true
    }
    fn any_edge(_: EdgeId, _: EdgeId) -> bool {
        true
    }
    fn only_q0_t1(q: NodeId, t: NodeId) -> bool {
        q.0 == 0 && t.0 == 1
    }
    fn exclude_t0(_: NodeId, t: NodeId) -> bool {
        t.0 != 0
    }
    fn reject_q0_t1(q: NodeId, t: NodeId) -> bool {
        !(q.0 == 0 && t.0 == 1)
    }
    fn only_tedge1(_: EdgeId, t: EdgeId) -> bool {
        t.0 == 1
    }
    fn disc_node(q: NodeId, t: NodeId) -> bool {
        const QN: [u8; 3] = [0, 0, 1];
        const TN: [u8; 3] = [0, 1, 0];
        QN[q.index()] == TN[t.index()]
    }
    fn disc_edge(q: EdgeId, t: EdgeId) -> bool {
        const QE: [u8; 1] = [1];
        const TE: [u8; 2] = [0, 1];
        QE[q.index()] == TE[t.index()]
    }

    #[rstest]
    #[case::triangle_self(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), any_node, any_edge,
        vec![vec![0, 1, 2], vec![0, 2, 1], vec![1, 0, 2], vec![1, 2, 0], vec![2, 0, 1], vec![2, 1, 0]])]
    #[case::empty_query(Graph::new(3, &[[0, 1], [1, 2]]), Graph::default(), any_node, any_edge, vec![vec![]])]
    #[case::node_filter(Graph::new(3, &[]), Graph::new(1, &[]), only_q0_t1, any_edge, vec![vec![1]])]
    fn test_graph_visit_subgraph_isomorphisms(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        for alg in ALGORITHMS {
            let mut r: Vec<Vec<usize>> = Vec::new();
            let flow: ControlFlow<()> = target.visit_subgraph_isomorphisms(
                &query,
                &mut node_match,
                &mut edge_match,
                alg,
                |embedding| {
                    r.push(embedding.iter().map(|n| n.index()).collect());
                    ControlFlow::Continue(())
                },
            );
            assert_eq!(flow, ControlFlow::Continue(()), "algorithm {alg:?}");
            r.sort();
            assert_eq!(r, expected, "algorithm {alg:?}");
        }
    }

    #[rstest]
    #[case::vf2_rdkit(Vf2Rdkit)]
    #[case::ray_kirsch(RayKirsch)]
    #[case::bitset_ri(BitsetRi)]
    #[case::glasgow(Glasgow)]
    fn test_graph_visit_subgraph_isomorphisms_termination(
        #[case] alg: SubgraphIsomorphismAlgorithm,
    ) {
        let graph = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let expected = [
            vec![0, 1, 2],
            vec![0, 2, 1],
            vec![1, 0, 2],
            vec![1, 2, 0],
            vec![2, 0, 1],
            vec![2, 1, 0],
        ];

        let full_calls = Cell::new(0usize);
        let mut full_emissions = 0usize;
        let full: ControlFlow<Vec<usize>> = graph.visit_subgraph_isomorphisms(
            &graph,
            &mut |_, _| {
                full_calls.set(full_calls.get() + 1);
                true
            },
            &mut |_, _| {
                full_calls.set(full_calls.get() + 1);
                true
            },
            alg,
            |_| {
                full_emissions += 1;
                ControlFlow::Continue(())
            },
        );
        assert_eq!(full, ControlFlow::Continue(()));
        assert_eq!(full_emissions, 6);

        let early_calls = Cell::new(0usize);
        let early = graph.visit_subgraph_isomorphisms(
            &graph,
            &mut |_, _| {
                early_calls.set(early_calls.get() + 1);
                true
            },
            &mut |_, _| {
                early_calls.set(early_calls.get() + 1);
                true
            },
            alg,
            |embedding| ControlFlow::Break(embedding.iter().map(|n| n.index()).collect()),
        );
        let ControlFlow::Break(first) = early else {
            panic!("expected Break on first emission");
        };
        assert!(expected.contains(&first), "invalid embedding {first:?}");
        assert!(
            early_calls.get() < full_calls.get(),
            "early {} vs full {}",
            early_calls.get(),
            full_calls.get()
        );
    }

    #[rstest]
    #[case::empty_query(Graph::new(3, &[[0, 1], [1, 2]]), Graph::default(), any_node, any_edge, vec![vec![]])]
    #[case::query_larger(Graph::new(2, &[]), Graph::new(3, &[]), any_node, any_edge, vec![])]
    #[case::single_node(Graph::new(3, &[]), Graph::new(1, &[]), any_node, any_edge, vec![vec![0], vec![1], vec![2]])]
    #[case::single_edge_identity(Graph::new(2, &[[0, 1]]), Graph::new(2, &[[0, 1]]), any_node, any_edge,
        vec![vec![0, 1], vec![1, 0]])]
    #[case::edge_in_chain(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(2, &[[0, 1]]), any_node, any_edge,
        vec![vec![0, 1], vec![1, 0], vec![1, 2], vec![2, 1]])]
    #[case::triangle_self(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), any_node, any_edge,
        vec![vec![0, 1, 2], vec![0, 2, 1], vec![1, 0, 2], vec![1, 2, 0], vec![2, 0, 1], vec![2, 1, 0]])]
    #[case::no_match(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), any_node, any_edge, vec![])]
    #[case::disconnected_query(Graph::new(3, &[[0, 1]]), Graph::new(2, &[]), any_node, any_edge,
        vec![vec![0, 1], vec![0, 2], vec![1, 0], vec![1, 2], vec![2, 0], vec![2, 1]])]
    #[case::disconnected_with_edge(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(3, &[[1, 2]]),
        disc_node, disc_edge, vec![vec![0, 2, 1]])]
    #[case::path_in_cycle(Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]), Graph::new(3, &[[0, 1], [1, 2]]), any_node, any_edge,
        vec![vec![0, 1, 2], vec![0, 3, 2], vec![1, 0, 3], vec![1, 2, 3], vec![2, 1, 0], vec![2, 3, 0], vec![3, 0, 1], vec![3, 2, 1]])]
    #[case::k4_self(Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]),
        Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]), any_node, any_edge,
        vec![vec![0, 1, 2, 3], vec![0, 1, 3, 2], vec![0, 2, 1, 3], vec![0, 2, 3, 1],
        vec![0, 3, 1, 2], vec![0, 3, 2, 1], vec![1, 0, 2, 3], vec![1, 0, 3, 2],
        vec![1, 2, 0, 3], vec![1, 2, 3, 0], vec![1, 3, 0, 2], vec![1, 3, 2, 0],
        vec![2, 0, 1, 3], vec![2, 0, 3, 1], vec![2, 1, 0, 3], vec![2, 1, 3, 0],
        vec![2, 3, 0, 1], vec![2, 3, 1, 0], vec![3, 0, 1, 2], vec![3, 0, 2, 1],
        vec![3, 1, 0, 2], vec![3, 1, 2, 0], vec![3, 2, 0, 1], vec![3, 2, 1, 0]])]
    #[case::star_in_star(Graph::new(5, &[[0, 1], [0, 2], [0, 3], [0, 4]]), Graph::new(4, &[[0, 1], [0, 2], [0, 3]]), any_node, any_edge,
        vec![vec![0, 1, 2, 3], vec![0, 1, 2, 4], vec![0, 1, 3, 2], vec![0, 1, 3, 4],
        vec![0, 1, 4, 2], vec![0, 1, 4, 3], vec![0, 2, 1, 3], vec![0, 2, 1, 4],
        vec![0, 2, 3, 1], vec![0, 2, 3, 4], vec![0, 2, 4, 1], vec![0, 2, 4, 3],
        vec![0, 3, 1, 2], vec![0, 3, 1, 4], vec![0, 3, 2, 1], vec![0, 3, 2, 4],
        vec![0, 3, 4, 1], vec![0, 3, 4, 2], vec![0, 4, 1, 2], vec![0, 4, 1, 3],
        vec![0, 4, 2, 1], vec![0, 4, 2, 3], vec![0, 4, 3, 1], vec![0, 4, 3, 2]])]
    #[case::naphthalene_ring(Graph::new(10, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [3, 6], [6, 7], [7, 8], [8, 9], [9, 4]]),
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]), any_node, any_edge,
        vec![vec![0, 1, 2, 3, 4, 5], vec![0, 5, 4, 3, 2, 1], vec![1, 0, 5, 4, 3, 2],
        vec![1, 2, 3, 4, 5, 0], vec![2, 1, 0, 5, 4, 3], vec![2, 3, 4, 5, 0, 1],
        vec![3, 2, 1, 0, 5, 4], vec![3, 4, 5, 0, 1, 2], vec![3, 4, 9, 8, 7, 6],
        vec![3, 6, 7, 8, 9, 4], vec![4, 3, 2, 1, 0, 5], vec![4, 3, 6, 7, 8, 9],
        vec![4, 5, 0, 1, 2, 3], vec![4, 9, 8, 7, 6, 3], vec![5, 0, 1, 2, 3, 4],
        vec![5, 4, 3, 2, 1, 0], vec![6, 3, 4, 9, 8, 7], vec![6, 7, 8, 9, 4, 3],
        vec![7, 6, 3, 4, 9, 8], vec![7, 8, 9, 4, 3, 6], vec![8, 7, 6, 3, 4, 9],
        vec![8, 9, 4, 3, 6, 7], vec![9, 4, 3, 6, 7, 8], vec![9, 8, 7, 6, 3, 4]])]
    #[case::node_filter(Graph::new(3, &[]), Graph::new(1, &[]), only_q0_t1, any_edge, vec![vec![1]])]
    #[case::edge_filter(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(2, &[[0, 1]]), any_node, only_tedge1, vec![vec![1, 2], vec![2, 1]])]
    fn test_graph_enumerate_subgraph_isomorphisms(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        for alg in ALGORITHMS {
            let mut r: Vec<Vec<usize>> = target
                .enumerate_subgraph_isomorphisms(&query, &mut node_match, &mut edge_match, alg)
                .iter()
                .map(|c| {
                    c.matched_pairs()
                        .iter()
                        .map(|&(_, host)| host.index())
                        .collect()
                })
                .collect();
            r.sort();
            assert_eq!(r, expected, "algorithm {alg:?}");
        }
    }

    #[rstest]
    #[case::fixes_anchor(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(1)), any_node, any_edge,
        vec![vec![1, 0], vec![1, 2]])]
    #[case::empty_query(Graph::new(3, &[[0, 1], [1, 2]]), Graph::default(), (NodeId(0), NodeId(0)), any_node, any_edge, vec![])]
    #[case::out_of_bounds(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(5), NodeId(0)), any_node, any_edge, vec![])]
    fn test_graph_visit_subgraph_isomorphisms_at(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] anchor: (NodeId, NodeId),
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        for alg in ALGORITHMS {
            let mut r: Vec<Vec<usize>> = Vec::new();
            let flow: ControlFlow<()> = target.visit_subgraph_isomorphisms_at(
                &query,
                anchor,
                &mut node_match,
                &mut edge_match,
                alg,
                |embedding| {
                    r.push(embedding.iter().map(|n| n.index()).collect());
                    ControlFlow::Continue(())
                },
            );
            assert_eq!(flow, ControlFlow::Continue(()), "algorithm {alg:?}");
            r.sort();
            assert_eq!(r, expected, "algorithm {alg:?}");
        }
    }

    #[rstest]
    #[case::fixes_anchor(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(1)), any_node, any_edge,
        vec![vec![1, 0], vec![1, 2]])]
    #[case::no_adjacency(Graph::new(3, &[[1, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(0)), any_node, any_edge, vec![])]
    #[case::empty_query(Graph::new(3, &[[0, 1], [1, 2]]), Graph::default(), (NodeId(0), NodeId(0)), any_node, any_edge, vec![])]
    #[case::out_of_bounds(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(5), NodeId(0)), any_node, any_edge, vec![])]
    #[case::path_anchored(Graph::new(5, &[[0, 1], [1, 2], [2, 3], [3, 4]]), Graph::new(3, &[[0, 1], [1, 2]]), (NodeId(1), NodeId(2)), any_node, any_edge,
        vec![vec![1, 2, 3], vec![3, 2, 1]])]
    #[case::node_filter(Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(1)), exclude_t0, any_edge, vec![vec![1, 2]])]
    #[case::node_filter_rejects_anchor(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(1)), reject_q0_t1, any_edge, vec![])]
    #[case::edge_filter(Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(1)), any_node, only_tedge1, vec![vec![1, 2]])]
    fn test_graph_enumerate_subgraph_isomorphisms_at(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] anchor: (NodeId, NodeId),
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        for alg in ALGORITHMS {
            let mut r: Vec<Vec<usize>> = target
                .enumerate_subgraph_isomorphisms_at(
                    &query,
                    anchor,
                    &mut node_match,
                    &mut edge_match,
                    alg,
                )
                .iter()
                .map(|c| {
                    c.matched_pairs()
                        .iter()
                        .map(|&(_, host)| host.index())
                        .collect()
                })
                .collect();
            r.sort();
            assert_eq!(r, expected, "algorithm {alg:?}");
        }
    }

    #[rstest]
    #[case::empty_query(Graph::new(3, &[[0, 1], [1, 2]]), Graph::default(), any_node, any_edge, vec![vec![]])]
    #[case::query_larger(Graph::new(2, &[]), Graph::new(3, &[]), any_node, any_edge, vec![])]
    #[case::no_match(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), any_node, any_edge, vec![])]
    #[case::triangle_self(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        any_node,
        any_edge,
        vec![vec![0, 1, 2], vec![0, 2, 1], vec![1, 0, 2], vec![1, 2, 0], vec![2, 0, 1], vec![2, 1, 0]]
    )]
    #[case::edge_filter(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(2, &[[0, 1]]), any_node, only_tedge1, vec![vec![1, 2], vec![2, 1]])]
    fn test_graph_subgraph_isomorphism_modes(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        for alg in ALGORITHMS {
            assert_eq!(
                target.has_subgraph_isomorphism(&query, &mut node_match, &mut edge_match, alg),
                !expected.is_empty(),
                "has algorithm {alg:?}"
            );
            assert_eq!(
                target.count_subgraph_isomorphisms(&query, &mut node_match, &mut edge_match, alg),
                expected.len() as u64,
                "count algorithm {alg:?}"
            );
            match target.first_subgraph_isomorphism(&query, &mut node_match, &mut edge_match, alg) {
                None => assert!(expected.is_empty(), "first algorithm {alg:?}"),
                Some(corr) => {
                    let imgs: Vec<usize> = corr
                        .matched_pairs()
                        .iter()
                        .map(|&(_, host)| host.index())
                        .collect();
                    assert!(
                        expected.contains(&imgs),
                        "first {imgs:?} not in expected for {alg:?}"
                    );
                }
            }
        }
    }

    #[rstest]
    #[case::fixes_anchor(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        Graph::new(2, &[[0, 1]]),
        (NodeId(0), NodeId(1)),
        any_node,
        any_edge,
        vec![vec![1, 0], vec![1, 2]]
    )]
    #[case::empty_query(Graph::new(3, &[[0, 1], [1, 2]]), Graph::default(), (NodeId(0), NodeId(0)), any_node, any_edge, vec![])]
    #[case::out_of_bounds(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(5), NodeId(0)), any_node, any_edge, vec![])]
    #[case::node_filter_rejects_anchor(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(1)), reject_q0_t1, any_edge, vec![])]
    fn test_graph_subgraph_isomorphism_modes_at(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] anchor: (NodeId, NodeId),
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        for alg in ALGORITHMS {
            assert_eq!(
                target.has_subgraph_isomorphism_at(
                    &query,
                    anchor,
                    &mut node_match,
                    &mut edge_match,
                    alg
                ),
                !expected.is_empty(),
                "has_at algorithm {alg:?}"
            );
            assert_eq!(
                target.count_subgraph_isomorphisms_at(
                    &query,
                    anchor,
                    &mut node_match,
                    &mut edge_match,
                    alg
                ),
                expected.len() as u64,
                "count_at algorithm {alg:?}"
            );
            match target.first_subgraph_isomorphism_at(
                &query,
                anchor,
                &mut node_match,
                &mut edge_match,
                alg,
            ) {
                None => assert!(expected.is_empty(), "first_at algorithm {alg:?}"),
                Some(corr) => {
                    let imgs: Vec<usize> = corr
                        .matched_pairs()
                        .iter()
                        .map(|&(_, host)| host.index())
                        .collect();
                    assert!(
                        expected.contains(&imgs),
                        "first_at {imgs:?} not in expected for {alg:?}"
                    );
                }
            }
        }
    }
}
