//! Subgraph isomorphism via VF2.

use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubgraphIsomorphismAlgorithm {
    Vf2,
}

impl Graph {
    pub fn subgraph_isomorphisms(
        &self,
        query: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<usize>> {
        match alg {
            SubgraphIsomorphismAlgorithm::Vf2 => {
                self.subgraph_isomorphisms_vf2(query, node_match, edge_match)
            }
        }
    }

    pub fn subgraph_isomorphisms_at(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<usize>> {
        match alg {
            SubgraphIsomorphismAlgorithm::Vf2 => {
                self.subgraph_isomorphisms_at_vf2(query, anchor, node_match, edge_match)
            }
        }
    }

    // Cordella et al. 2004 "A (sub)graph isomorphism algorithm for matching large graphs".
    fn subgraph_isomorphisms_vf2(
        &self,
        query: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) -> Vec<Vec<usize>> {
        if query.node_count() > self.node_count() {
            return Vec::new();
        }
        if query.node_count() == 0 {
            return vec![vec![]];
        }

        let mut state = Vf2State::new(query, self);
        state.search(node_match, edge_match);
        state.results
    }

    // Cordella et al. 2004 "A (sub)graph isomorphism algorithm for matching large graphs".
    fn subgraph_isomorphisms_at_vf2(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) -> Vec<Vec<usize>> {
        if query.node_count() > self.node_count() {
            return Vec::new();
        }
        if query.node_count() == 0 {
            return Vec::new();
        }
        if anchor.0.index() >= query.node_count() || anchor.1.index() >= self.node_count() {
            return Vec::new();
        }
        if !node_match(anchor.0, anchor.1) {
            return Vec::new();
        }

        let mut state = Vf2State::new(query, self);
        state.seed_anchor(anchor);
        state.search(node_match, edge_match);
        state.results
    }
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

    use super::SubgraphIsomorphismAlgorithm::Vf2;
    use super::*;

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
    fn test_subgraph_isomorphisms(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        let mut r = target.subgraph_isomorphisms(&query, &mut node_match, &mut edge_match, Vf2);
        r.sort();
        assert_eq!(r, expected);
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
    fn test_subgraph_isomorphisms_at(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] anchor: (NodeId, NodeId),
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        let mut r =
            target.subgraph_isomorphisms_at(&query, anchor, &mut node_match, &mut edge_match, Vf2);
        r.sort();
        assert_eq!(r, expected);
    }
}
