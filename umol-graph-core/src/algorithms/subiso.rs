//! Subgraph isomorphism: find all injective node mappings from a query graph
//! into a target graph that preserve adjacency. Node and edge compatibility is
//! caller-defined via match closures. Implemented with VF2 (Cordella et al.,
//! 2004).

use crate::graph::{EdgeId, Graph, NodeId};

/// Find all subgraph isomorphisms from `query` into `target`.
///
/// Returns assignments where `result[i][q]` is the index of the target node
/// that query node `q` maps to. The mapping is injective and preserves
/// adjacency: for every edge `(u, v)` in the query, an edge exists between
/// the mapped target nodes satisfying `edge_match`.
pub fn subgraph_isomorphisms(
    query: &Graph,
    target: &Graph,
    node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
    edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
) -> Vec<Vec<usize>> {
    if query.node_count() > target.node_count() {
        return Vec::new();
    }
    if query.node_count() == 0 {
        return vec![vec![]];
    }

    let mut state = Vf2State::new(query, target);
    state.search(node_match, edge_match);
    state.results
}

/// Find all subgraph isomorphisms from `query` into `target` where `anchor.0`
/// (a query node) must map to `anchor.1` (a target node). Results have the
/// same shape as [`subgraph_isomorphisms`], with `result[i][anchor.0] ==
/// anchor.1.index()` in every returned assignment. Returns empty when query
/// is empty, when either anchor endpoint is out of bounds, or when
/// `node_match(anchor.0, anchor.1)` is false.
pub fn subgraph_isomorphisms_at(
    query: &Graph,
    target: &Graph,
    anchor: (NodeId, NodeId),
    node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
    edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
) -> Vec<Vec<usize>> {
    if query.node_count() > target.node_count() {
        return Vec::new();
    }
    if query.node_count() == 0 {
        return Vec::new();
    }
    if anchor.0.index() >= query.node_count() || anchor.1.index() >= target.node_count() {
        return Vec::new();
    }
    if !node_match(anchor.0, anchor.1) {
        return Vec::new();
    }

    let mut state = Vf2State::new(query, target);
    state.seed_anchor(anchor);
    state.search(node_match, edge_match);
    state.results
}

struct Vf2State<'g> {
    query: &'g Graph,
    target: &'g Graph,
    // query node index → target node index
    mapping: Vec<Option<u32>>,
    // target node index → already mapped?
    reverse: Vec<bool>,
    // Nonzero = depth at which node entered the terminal set.
    terminal_query: Vec<u32>,
    terminal_target: Vec<u32>,
    depth: u32,
    results: Vec<Vec<usize>>,
}

impl<'g> Vf2State<'g> {
    fn new(query: &'g Graph, target: &'g Graph) -> Self {
        Self {
            query,
            target,
            mapping: vec![None; query.node_count()],
            reverse: vec![false; target.node_count()],
            terminal_query: vec![0; query.node_count()],
            terminal_target: vec![0; target.node_count()],
            depth: 0,
            results: Vec::new(),
        }
    }

    /// Seed state with a forced pair before search. Equivalent to one level of
    /// recursion already done: pair mapped, terminal sets populated from the
    /// pair's neighbors.
    fn seed_anchor(&mut self, anchor: (NodeId, NodeId)) {
        let (q_node, t_node) = anchor;
        let q_idx = q_node.index();
        let t_idx = t_node.index();
        self.depth = 1;
        self.mapping[q_idx] = Some(t_idx as u32);
        self.reverse[t_idx] = true;
        for qn in self.query.neighbors(q_node) {
            let ni = qn.node.index();
            if self.terminal_query[ni] == 0 && self.mapping[ni].is_none() {
                self.terminal_query[ni] = self.depth;
            }
        }
        for tn in self.target.neighbors(t_node) {
            let ni = tn.node.index();
            if self.terminal_target[ni] == 0 && !self.reverse[ni] {
                self.terminal_target[ni] = self.depth;
            }
        }
    }

    fn search(
        &mut self,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) {
        if self.depth as usize == self.query.node_count() {
            let assignment = self.mapping.iter().map(|m| m.unwrap() as usize).collect();
            self.results.push(assignment);
            return;
        }

        let q_idx = self.next_query_node();
        let q_node = NodeId(q_idx as u32);
        let q_in_terminal = self.terminal_query[q_idx] > 0;

        for t_idx in 0..self.target.node_count() {
            if self.reverse[t_idx] {
                continue;
            }
            if q_in_terminal && self.terminal_target[t_idx] == 0 {
                continue;
            }

            let t_node = NodeId(t_idx as u32);
            if !node_match(q_node, t_node) {
                continue;
            }
            if !self.feasible(q_idx, t_idx, edge_match) {
                continue;
            }

            self.depth += 1;
            self.mapping[q_idx] = Some(t_idx as u32);
            self.reverse[t_idx] = true;

            let mut restore_q = Vec::new();
            let mut restore_t = Vec::new();

            for qn in self.query.neighbors(q_node) {
                let ni = qn.node.index();
                if self.terminal_query[ni] == 0 && self.mapping[ni].is_none() {
                    self.terminal_query[ni] = self.depth;
                    restore_q.push(ni);
                }
            }
            for tn in self.target.neighbors(t_node) {
                let ni = tn.node.index();
                if self.terminal_target[ni] == 0 && !self.reverse[ni] {
                    self.terminal_target[ni] = self.depth;
                    restore_t.push(ni);
                }
            }

            self.search(node_match, edge_match);

            self.depth -= 1;
            self.mapping[q_idx] = None;
            self.reverse[t_idx] = false;
            for ni in restore_q {
                self.terminal_query[ni] = 0;
            }
            for ni in restore_t {
                self.terminal_target[ni] = 0;
            }
        }
    }

    /// Smallest unmapped query node, preferring terminal nodes.
    fn next_query_node(&self) -> usize {
        for i in 0..self.query.node_count() {
            if self.mapping[i].is_none() && self.terminal_query[i] > 0 {
                return i;
            }
        }
        for i in 0..self.query.node_count() {
            if self.mapping[i].is_none() {
                return i;
            }
        }
        unreachable!()
    }

    fn feasible(
        &self,
        q_idx: usize,
        t_idx: usize,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) -> bool {
        let q_node = NodeId(q_idx as u32);
        let t_node = NodeId(t_idx as u32);

        // Consistency: each mapped query neighbor must map to a target neighbor
        // with a matching edge.
        for qn in self.query.neighbors(q_node) {
            let Some(mapped) = self.mapping[qn.node.index()] else {
                continue;
            };
            let Some(t_edge) = self.target.find_edge(t_node, NodeId(mapped)) else {
                return false;
            };
            if !edge_match(qn.edge, t_edge) {
                return false;
            }
        }

        // Look-ahead: unmapped query neighbors must not exceed unmapped target
        // neighbors in either the terminal or non-terminal category.
        let (q_term, q_new) = self.unmapped_neighbor_counts_query(q_idx);
        let (t_term, t_new) = self.unmapped_neighbor_counts_target(t_idx);

        q_term <= t_term && q_new <= t_new
    }

    fn unmapped_neighbor_counts_query(&self, idx: usize) -> (u32, u32) {
        let mut term = 0u32;
        let mut new = 0u32;
        for n in self.query.neighbors(NodeId(idx as u32)) {
            if self.mapping[n.node.index()].is_none() {
                if self.terminal_query[n.node.index()] > 0 {
                    term += 1;
                } else {
                    new += 1;
                }
            }
        }
        (term, new)
    }

    fn unmapped_neighbor_counts_target(&self, idx: usize) -> (u32, u32) {
        let mut term = 0u32;
        let mut new = 0u32;
        for n in self.target.neighbors(NodeId(idx as u32)) {
            if !self.reverse[n.node.index()] {
                if self.terminal_target[n.node.index()] > 0 {
                    term += 1;
                } else {
                    new += 1;
                }
            }
        }
        (term, new)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn accept_all_nodes(_: NodeId, _: NodeId) -> bool {
        true
    }

    fn accept_all_edges(_: EdgeId, _: EdgeId) -> bool {
        true
    }

    fn sorted(mut results: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
        results.sort();
        results
    }

    #[test]
    fn test_vf2_empty_query() {
        let q = Graph::default();
        let t = Graph::new(3, &[[0, 1], [1, 2]]);
        let r = subgraph_isomorphisms(&q, &t, &mut accept_all_nodes, &mut accept_all_edges);
        assert_eq!(r, vec![vec![]]);
    }

    #[test]
    fn test_vf2_query_larger_than_target() {
        let q = Graph::new(3, &[]);
        let t = Graph::new(2, &[]);
        let r = subgraph_isomorphisms(&q, &t, &mut accept_all_nodes, &mut accept_all_edges);
        assert!(r.is_empty());
    }

    #[test]
    fn test_vf2_single_node() {
        let q = Graph::new(1, &[]);
        let t = Graph::new(3, &[]);
        let r = sorted(subgraph_isomorphisms(
            &q,
            &t,
            &mut accept_all_nodes,
            &mut accept_all_edges,
        ));
        assert_eq!(r, vec![vec![0], vec![1], vec![2]]);
    }

    #[test]
    fn test_vf2_single_edge_identity() {
        let g = Graph::new(2, &[[0, 1]]);
        let r = sorted(subgraph_isomorphisms(
            &g,
            &g,
            &mut accept_all_nodes,
            &mut accept_all_edges,
        ));
        assert_eq!(r, vec![vec![0, 1], vec![1, 0]]);
    }

    #[test]
    fn test_vf2_edge_in_chain() {
        let q = Graph::new(2, &[[0, 1]]);
        let t = Graph::new(3, &[[0, 1], [1, 2]]);
        let r = sorted(subgraph_isomorphisms(
            &q,
            &t,
            &mut accept_all_nodes,
            &mut accept_all_edges,
        ));
        assert_eq!(r, vec![vec![0, 1], vec![1, 0], vec![1, 2], vec![2, 1]]);
    }

    #[test]
    fn test_vf2_triangle_in_triangle() {
        let g = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let r = subgraph_isomorphisms(
            &g,
            &g,
            &mut accept_all_nodes,
            &mut accept_all_edges,
        );
        // 3! = 6 automorphisms of K3
        assert_eq!(r.len(), 6);
    }

    #[test]
    fn test_vf2_no_match_missing_edge() {
        let q = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let t = Graph::new(3, &[[0, 1], [1, 2]]);
        let r = subgraph_isomorphisms(
            &q,
            &t,
            &mut accept_all_nodes,
            &mut accept_all_edges,
        );
        assert!(r.is_empty());
    }

    #[test]
    fn test_vf2_node_match_filter() {
        let q = Graph::new(1, &[]);
        let t = Graph::new(3, &[]);
        // Only match node 0 → node 1
        let mut nm = |q: NodeId, t: NodeId| q.0 == 0 && t.0 == 1;
        let r = subgraph_isomorphisms(&q, &t, &mut nm, &mut accept_all_edges);
        assert_eq!(r, vec![vec![1]]);
    }

    #[test]
    fn test_vf2_edge_match_filter() {
        let q = Graph::new(2, &[[0, 1]]);
        // Target has two edges; only edge 1 should match.
        let t = Graph::new(3, &[[0, 1], [1, 2]]);
        let mut em = |_q: EdgeId, t: EdgeId| t.0 == 1;
        let r = sorted(subgraph_isomorphisms(
            &q,
            &t,
            &mut accept_all_nodes,
            &mut em,
        ));
        assert_eq!(r, vec![vec![1, 2], vec![2, 1]]);
    }

    #[test]
    fn test_vf2_disconnected_query() {
        let q = Graph::new(2, &[]);
        let t = Graph::new(3, &[[0, 1]]);
        let r = subgraph_isomorphisms(
            &q,
            &t,
            &mut accept_all_nodes,
            &mut accept_all_edges,
        );
        // 3 * 2 = 6 injective mappings of 2 isolated nodes into 3 nodes
        assert_eq!(r.len(), 6);
    }

    #[test]
    fn test_vf2_path_in_cycle() {
        let q = Graph::new(3, &[[0, 1], [1, 2]]);
        let t = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]);
        let r = subgraph_isomorphisms(
            &q,
            &t,
            &mut accept_all_nodes,
            &mut accept_all_edges,
        );
        // 4 positions × 2 directions = 8
        assert_eq!(r.len(), 8);
    }

    #[test]
    fn test_vf2_self_isomorphism_k4() {
        let g = Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]);
        let r = subgraph_isomorphisms(
            &g,
            &g,
            &mut accept_all_nodes,
            &mut accept_all_edges,
        );
        // |Aut(K4)| = 4! = 24
        assert_eq!(r.len(), 24);
    }

    #[test]
    fn test_vf2_star_in_star() {
        // 3-star query in 4-star target
        let q = Graph::new(4, &[[0, 1], [0, 2], [0, 3]]);
        let t = Graph::new(5, &[[0, 1], [0, 2], [0, 3], [0, 4]]);
        let r = subgraph_isomorphisms(
            &q,
            &t,
            &mut accept_all_nodes,
            &mut accept_all_edges,
        );
        // Center must map to center (0→0). Leaves: P(4,3) = 24.
        assert_eq!(r.len(), 24);
    }

    #[test]
    fn test_subgraph_isomorphisms_at_fixes_anchor() {
        // Query edge 0-1 into 3-cycle target. Anchor q0→t1.
        let q = Graph::new(2, &[[0, 1]]);
        let t = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let r = sorted(subgraph_isomorphisms_at(
            &q,
            &t,
            (NodeId(0), NodeId(1)),
            &mut accept_all_nodes,
            &mut accept_all_edges,
        ));
        // q0 must be t1, so q1 can be t0 or t2.
        assert_eq!(r, vec![vec![1, 0], vec![1, 2]]);
    }

    #[test]
    fn test_subgraph_isomorphisms_at_anchor_rejects_when_adjacency_fails() {
        // Query is edge 0-1, anchor q0→t0 in a graph where t0 has no neighbors.
        let q = Graph::new(2, &[[0, 1]]);
        let t = Graph::new(3, &[[1, 2]]);
        let r = subgraph_isomorphisms_at(
            &q,
            &t,
            (NodeId(0), NodeId(0)),
            &mut accept_all_nodes,
            &mut accept_all_edges,
        );
        assert!(r.is_empty());
    }

    #[test]
    fn test_subgraph_isomorphisms_at_empty_query() {
        let q = Graph::default();
        let t = Graph::new(3, &[[0, 1], [1, 2]]);
        let r = subgraph_isomorphisms_at(
            &q,
            &t,
            (NodeId(0), NodeId(0)),
            &mut accept_all_nodes,
            &mut accept_all_edges,
        );
        assert!(r.is_empty());
    }

    #[test]
    fn test_subgraph_isomorphisms_at_out_of_bounds() {
        let q = Graph::new(2, &[[0, 1]]);
        let t = Graph::new(3, &[[0, 1], [1, 2]]);
        let r = subgraph_isomorphisms_at(
            &q,
            &t,
            (NodeId(5), NodeId(0)),
            &mut accept_all_nodes,
            &mut accept_all_edges,
        );
        assert!(r.is_empty());
    }

    #[test]
    fn test_subgraph_isomorphisms_at_node_match_gates_anchor() {
        let q = Graph::new(2, &[[0, 1]]);
        let t = Graph::new(3, &[[0, 1], [1, 2]]);
        let mut nm = |q: NodeId, t: NodeId| !(q.0 == 0 && t.0 == 1);
        let r = subgraph_isomorphisms_at(
            &q,
            &t,
            (NodeId(0), NodeId(1)),
            &mut nm,
            &mut accept_all_edges,
        );
        assert!(r.is_empty());
    }

    #[test]
    fn test_subgraph_isomorphisms_at_path_anchored_at_center() {
        // Path 0-1-2 into 5-path. Anchor q1→t2 (center to center).
        let q = Graph::new(3, &[[0, 1], [1, 2]]);
        let t = Graph::new(5, &[[0, 1], [1, 2], [2, 3], [3, 4]]);
        let r = sorted(subgraph_isomorphisms_at(
            &q,
            &t,
            (NodeId(1), NodeId(2)),
            &mut accept_all_nodes,
            &mut accept_all_edges,
        ));
        // q1=t2. Neighbors of t2 are {1,3}. Two orderings: (q0=1,q2=3), (q0=3,q2=1).
        assert_eq!(r, vec![vec![1, 2, 3], vec![3, 2, 1]]);
    }

    #[rstest]
    #[case::naphthalene_ring(
        // Query: 6-cycle
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        // Target: naphthalene (two fused 6-cycles sharing edge 3-4)
        Graph::new(10, &[
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0],
            [3, 6], [6, 7], [7, 8], [8, 9], [9, 4],
        ]),
        // Two rings, each with 2 orientations (CW/CCW) = 4
        // But wait: the 6-cycle has rotational + reflective symmetry,
        // so unique embeddings are fewer... actually no, VF2 returns
        // all distinct node mappings.
        // Ring 1: 0-1-2-3-4-5 — 6 rotations × 2 reflections = 12
        // Ring 2: 3-6-7-8-9-4 — 6 rotations × 2 reflections = 12
        // But ring 2 only has 5 unique nodes (3,6,7,8,9,4), 6 total in cycle,
        // and the cycle is 3-6-7-8-9-4 with edge 3-4 closing it.
        // Actually both are valid 6-cycles in the naphthalene graph.
        // Each 6-cycle automorphism gives a distinct assignment.
        24
    )]
    fn test_vf2_subgraph(
        #[case] query: Graph,
        #[case] target: Graph,
        #[case] expected_count: usize,
    ) {
        let r = subgraph_isomorphisms(
            &query,
            &target,
            &mut accept_all_nodes,
            &mut accept_all_edges,
        );
        assert_eq!(r.len(), expected_count);
    }
}
