//! Subgraph-isomorphism search with monomorphism semantics.
//!
//! Query edges must map to target edges; no induced-subgraph reverse check is
//! applied. Current selectors provide VF2, Ullmann, RI, ArcMatch, an
//! RDKit-compatible VF2 variant, and Ray--Kirsch bond search. Both ordinary and
//! anchored operations return the same correspondence set for every selector.
//! See [Cordella et al. (2004)](https://doi.org/10.1109/TPAMI.2004.75),
//! [Ullmann (1976)](https://doi.org/10.1145/321921.321925),
//! [Bonnici et al. (2013)](https://doi.org/10.1186/1471-2105-14-S7-S13),
//! [Bonnici et al. (2024)](https://doi.org/10.1007/s10618-024-01061-8),
//! [RDKit PR #2500](https://github.com/rdkit/rdkit/pull/2500), and
//! [Ray and Kirsch (1957)](https://doi.org/10.1126/science.126.3278.814).

use std::cmp::Reverse;
use std::collections::HashMap;

use crate::correspondence::Correspondence;
use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubgraphIsomorphismAlgorithm {
    /// VF2 — terminal-set backtracking + look-ahead (Cordella et al. 2004).
    Vf2,
    /// Ullmann — candidate-matrix refinement + backtracking (Ullmann 1976).
    Ullmann,
    /// RI — static most-constrained-first pattern ordering, light matching, no
    /// domain reduction (Bonnici et al. 2013).
    Ri,
    /// ArcMatch — arc-consistent vertex/edge-domain reduction + path-based reduction
    /// (exploits edge labels), 5-measure ordering, dynamic-parent backtracking
    /// (Bonnici et al. 2024). `path_length` is the max query-path length for the
    /// path-based reduction (Thm 1: larger = more pruning, higher cost); `< 3`
    /// disables it (arc consistency only). See `ARCMATCH_DEFAULT_PATH_LENGTH`.
    ArcMatch { path_length: usize },
    /// RDKit's substructure engine: vflib-derived VF2 with John Mayfield's
    /// chemical-graph optimizations (RDKit PR #2500) — candidates restricted to a
    /// mapped neighbor's adjacency, an explicit `deg(q) <= deg(t)` bound, and the
    /// terminal-set look-ahead disabled. A benchmark reference.
    Vf2Rdkit,
    /// Ray-Kirsch bond-based backtracking (Ray & Kirsch 1957; Sayle `parsmart.cpp`):
    /// the search advances over query *bonds*, so each extension's adjacency holds
    /// by construction and no per-neighbor re-check is needed. A benchmark reference.
    RayKirsch,
}

/// Paper's recommended ArcMatch path-length trade-off (Bonnici 2024 §3.1).
pub const ARCMATCH_DEFAULT_PATH_LENGTH: usize = 6;

impl Graph {
    // TODO: add singular `subgraph_isomorphism(...) -> Option<Vec<usize>>` that stops at
    // the first match (existence via `.is_some()`). Saves the embedding-multiplicity
    // factor on positive matches; same Vf2, `search` gains an early-exit cap.
    pub fn enumerate_subgraph_isomorphisms(
        &self,
        query: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Correspondence<NodeId>> {
        let embeddings = match alg {
            SubgraphIsomorphismAlgorithm::Vf2 => {
                self.subgraph_isomorphisms_vf2(query, node_match, edge_match)
            }
            SubgraphIsomorphismAlgorithm::Ullmann => {
                self.subgraph_isomorphisms_ullmann(query, node_match, edge_match)
            }
            SubgraphIsomorphismAlgorithm::Ri => {
                self.subgraph_isomorphisms_ri(query, node_match, edge_match)
            }
            SubgraphIsomorphismAlgorithm::ArcMatch { path_length } => {
                self.subgraph_isomorphisms_arcmatch(query, node_match, edge_match, path_length)
            }
            SubgraphIsomorphismAlgorithm::Vf2Rdkit => {
                self.subgraph_isomorphisms_vf2rdkit(query, node_match, edge_match)
            }
            SubgraphIsomorphismAlgorithm::RayKirsch => {
                self.subgraph_isomorphisms_rk(query, node_match, edge_match)
            }
        };
        self.embeddings_to_correspondences(embeddings)
    }

    /// Lift each `query→host` embedding (dense left space `0..query.node_count()`) to a
    /// [`Correspondence`] over the host's node id space.
    fn embeddings_to_correspondences(
        &self,
        embeddings: Vec<Vec<usize>>,
    ) -> Vec<Correspondence<NodeId>> {
        embeddings
            .into_iter()
            .map(|embedding| {
                let images: Vec<NodeId> = embedding.into_iter().map(NodeId::from).collect();
                Correspondence::from_images(&images, self.node_count())
            })
            .collect()
    }

    pub fn enumerate_subgraph_isomorphisms_at(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Correspondence<NodeId>> {
        let embeddings = match alg {
            SubgraphIsomorphismAlgorithm::Vf2 => {
                self.subgraph_isomorphisms_at_vf2(query, anchor, node_match, edge_match)
            }
            SubgraphIsomorphismAlgorithm::Ullmann => {
                self.subgraph_isomorphisms_at_ullmann(query, anchor, node_match, edge_match)
            }
            SubgraphIsomorphismAlgorithm::Ri => {
                self.subgraph_isomorphisms_at_ri(query, anchor, node_match, edge_match)
            }
            SubgraphIsomorphismAlgorithm::ArcMatch { path_length } => self
                .subgraph_isomorphisms_at_arcmatch(
                    query,
                    anchor,
                    node_match,
                    edge_match,
                    path_length,
                ),
            SubgraphIsomorphismAlgorithm::Vf2Rdkit => {
                self.subgraph_isomorphisms_at_vf2rdkit(query, anchor, node_match, edge_match)
            }
            SubgraphIsomorphismAlgorithm::RayKirsch => {
                self.subgraph_isomorphisms_at_rk(query, anchor, node_match, edge_match)
            }
        };
        self.embeddings_to_correspondences(embeddings)
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

    // Ullmann 1976 "An algorithm for subgraph isomorphism": candidate matrix
    // (label- and degree-compatible) refined to a fixpoint, then row-by-row
    // backtracking with re-refinement. Monomorphism (substructure) semantics.
    fn subgraph_isomorphisms_ullmann(
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
        let mut m = ullmann_matrix(query, self, node_match);
        ullmann_refine(query, self, &mut m, edge_match);
        let mut results = Vec::new();
        let mut mapping = vec![usize::MAX; query.node_count()];
        let mut used = vec![false; self.node_count()];
        ullmann_search(
            query,
            self,
            0,
            &m,
            &mut mapping,
            &mut used,
            edge_match,
            &mut results,
        );
        results
    }

    fn subgraph_isomorphisms_at_ullmann(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) -> Vec<Vec<usize>> {
        if query.node_count() > self.node_count() || query.node_count() == 0 {
            return Vec::new();
        }
        if anchor.0.index() >= query.node_count() || anchor.1.index() >= self.node_count() {
            return Vec::new();
        }
        if !node_match(anchor.0, anchor.1) {
            return Vec::new();
        }
        let n2 = self.node_count();
        let (qa, ta) = (anchor.0.index(), anchor.1.index());
        let mut m = ullmann_matrix(query, self, node_match);
        // Pin the anchor: its row keeps only the anchor target, its target column
        // keeps only the anchor row.
        for k in 0..n2 {
            if k != ta {
                m[qa * n2 + k] = false;
            }
        }
        for r in 0..query.node_count() {
            if r != qa {
                m[r * n2 + ta] = false;
            }
        }
        ullmann_refine(query, self, &mut m, edge_match);
        let mut results = Vec::new();
        let mut mapping = vec![usize::MAX; query.node_count()];
        let mut used = vec![false; n2];
        ullmann_search(
            query,
            self,
            0,
            &m,
            &mut mapping,
            &mut used,
            edge_match,
            &mut results,
        );
        results
    }

    // Bonnici et al. 2013 "A subgraph isomorphism algorithm and its application to
    // biochemical data": a static most-constrained-first pattern ordering with light
    // edge-consistency matching and no domain reduction. Monomorphism semantics.
    fn subgraph_isomorphisms_ri(
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
        let order = ri_order(query, None);
        let parents = ri_parents(query, &order);
        let mut results = Vec::new();
        let mut mapping = vec![0usize; order.len()];
        let mut used = vec![false; self.node_count()];
        ri_search(
            self,
            &order,
            &parents,
            0,
            &mut mapping,
            &mut used,
            None,
            node_match,
            edge_match,
            &mut results,
        );
        results
    }

    fn subgraph_isomorphisms_at_ri(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) -> Vec<Vec<usize>> {
        if query.node_count() > self.node_count() || query.node_count() == 0 {
            return Vec::new();
        }
        if anchor.0.index() >= query.node_count() || anchor.1.index() >= self.node_count() {
            return Vec::new();
        }
        if !node_match(anchor.0, anchor.1) {
            return Vec::new();
        }
        let order = ri_order(query, Some(anchor.0.index()));
        let parents = ri_parents(query, &order);
        let mut results = Vec::new();
        let mut mapping = vec![0usize; order.len()];
        let mut used = vec![false; self.node_count()];
        ri_search(
            self,
            &order,
            &parents,
            0,
            &mut mapping,
            &mut used,
            Some(anchor.1.index()),
            node_match,
            edge_match,
            &mut results,
        );
        results
    }

    // Bonnici et al. 2024 "ArcMatch: high-performance subgraph matching for labeled
    // graphs": vertex domains → arc consistency → edge domains → path-based reduction
    // → 5-measure ordering → dynamic-parent backtracking. Monomorphism semantics.
    fn subgraph_isomorphisms_arcmatch(
        &self,
        query: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        path_length: usize,
    ) -> Vec<Vec<usize>> {
        if query.node_count() > self.node_count() {
            return Vec::new();
        }
        if query.node_count() == 0 {
            return vec![vec![]];
        }
        let mut vertex = arcmatch_vertex_domains(query, self, node_match);
        arcmatch_arc_consistency(query, self, &mut vertex, edge_match);
        let mut domains = arcmatch_edge_domains(query, self, vertex, edge_match);
        if path_length >= 3 {
            arcmatch_path_reduction(query, &mut domains, path_length);
        }
        let order = arcmatch_variable_ordering(query, &domains);
        let mut results = Vec::new();
        let mut assigned = vec![None; query.node_count()];
        let mut used = vec![false; self.node_count()];
        arcmatch_search(
            query,
            &order,
            &domains,
            0,
            &mut assigned,
            &mut used,
            &mut results,
        );
        results
    }

    fn subgraph_isomorphisms_at_arcmatch(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        path_length: usize,
    ) -> Vec<Vec<usize>> {
        if query.node_count() > self.node_count() || query.node_count() == 0 {
            return Vec::new();
        }
        if anchor.0.index() >= query.node_count() || anchor.1.index() >= self.node_count() {
            return Vec::new();
        }
        if !node_match(anchor.0, anchor.1) {
            return Vec::new();
        }
        let n2 = self.node_count();
        let (qa, ta) = (anchor.0.index(), anchor.1.index());
        let mut vertex = arcmatch_vertex_domains(query, self, node_match);
        // Pin the anchor: its row keeps only the anchor target, its target column
        // keeps only the anchor row.
        for t in 0..n2 {
            if t != ta {
                vertex[qa * n2 + t] = false;
            }
        }
        for q in 0..query.node_count() {
            if q != qa {
                vertex[q * n2 + ta] = false;
            }
        }
        arcmatch_arc_consistency(query, self, &mut vertex, edge_match);
        let mut domains = arcmatch_edge_domains(query, self, vertex, edge_match);
        if path_length >= 3 {
            arcmatch_path_reduction(query, &mut domains, path_length);
        }
        let order = arcmatch_variable_ordering(query, &domains);
        let mut results = Vec::new();
        let mut assigned = vec![None; query.node_count()];
        let mut used = vec![false; n2];
        arcmatch_search(
            query,
            &order,
            &domains,
            0,
            &mut assigned,
            &mut used,
            &mut results,
        );
        results
    }

    // RDKit's `vf2.hpp` (vflib-derived VF2 + Mayfield PR #2500): monomorphism with
    // candidates drawn from a mapped neighbor's image adjacency, an explicit degree
    // bound, and no terminal-set look-ahead. Benchmark reference for RDKit.
    fn subgraph_isomorphisms_vf2rdkit(
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
        let mut mapping = vec![None; query.node_count()];
        let mut used = vec![false; self.node_count()];
        let mut results = Vec::new();
        vf2rdkit_search(
            query,
            self,
            &mut mapping,
            &mut used,
            node_match,
            edge_match,
            &mut results,
        );
        results
    }

    fn subgraph_isomorphisms_at_vf2rdkit(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) -> Vec<Vec<usize>> {
        if query.node_count() > self.node_count() || query.node_count() == 0 {
            return Vec::new();
        }
        if anchor.0.index() >= query.node_count() || anchor.1.index() >= self.node_count() {
            return Vec::new();
        }
        if !node_match(anchor.0, anchor.1) {
            return Vec::new();
        }
        let mut mapping = vec![None; query.node_count()];
        let mut used = vec![false; self.node_count()];
        mapping[anchor.0.index()] = Some(anchor.1.index());
        used[anchor.1.index()] = true;
        let mut results = Vec::new();
        vf2rdkit_search(
            query,
            self,
            &mut mapping,
            &mut used,
            node_match,
            edge_match,
            &mut results,
        );
        results
    }

    // Ray-Kirsch bond-based search (Sayle `parsmart.cpp`): backtrack over query
    // bonds; degree-0 query atoms are mapped as a product afterwards. Benchmark
    // reference.
    fn subgraph_isomorphisms_rk(
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
        let bond_order = rk_bond_order(query);
        let isolated: Vec<usize> = (0..query.node_count())
            .filter(|&q| query.neighbors(NodeId(q as u32)).is_empty())
            .collect();
        let mut mapping = vec![None; query.node_count()];
        let mut used = vec![false; self.node_count()];
        let mut results = Vec::new();
        rk_search(
            self,
            &bond_order,
            &isolated,
            0,
            &mut mapping,
            &mut used,
            node_match,
            edge_match,
            &mut results,
        );
        results
    }

    fn subgraph_isomorphisms_at_rk(
        &self,
        query: &Graph,
        anchor: (NodeId, NodeId),
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) -> Vec<Vec<usize>> {
        if query.node_count() > self.node_count() || query.node_count() == 0 {
            return Vec::new();
        }
        if anchor.0.index() >= query.node_count() || anchor.1.index() >= self.node_count() {
            return Vec::new();
        }
        if !node_match(anchor.0, anchor.1) {
            return Vec::new();
        }
        let bond_order = rk_bond_order(query);
        let isolated: Vec<usize> = (0..query.node_count())
            .filter(|&q| query.neighbors(NodeId(q as u32)).is_empty())
            .collect();
        let mut mapping = vec![None; query.node_count()];
        let mut used = vec![false; self.node_count()];
        mapping[anchor.0.index()] = Some(anchor.1.index());
        used[anchor.1.index()] = true;
        let mut results = Vec::new();
        rk_search(
            self,
            &bond_order,
            &isolated,
            0,
            &mut mapping,
            &mut used,
            node_match,
            edge_match,
            &mut results,
        );
        results
    }
}

mod vf2 {
    //! VF2 terminal-set backtracking with monomorphism look-ahead.
    //!
    //! See [Cordella et al.
    //! (2004)](https://doi.org/10.1109/TPAMI.2004.75).

    use super::*;

    pub(super) struct Vf2State<'g> {
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
        pub(super) results: Vec<Vec<usize>>,
    }

    impl<'g> Vf2State<'g> {
        pub(super) fn new(query: &'g Graph, target: &'g Graph) -> Self {
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
        pub(super) fn seed_anchor(&mut self, anchor: (NodeId, NodeId)) {
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

        pub(super) fn search(
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

            // VF2 monomorphism look-ahead (cf. vf2lib `VF2MonoState`): unmapped terminal
            // query neighbors must not exceed unmapped terminal target neighbors, and
            // the total unmapped-neighbor count must not either. Both are necessary for a
            // monomorphism extension. (The induced-subgraph matcher additionally requires
            // `q_new <= t_new`, which is invalid here and under-reports disconnected
            // queries.)
            let (q_term, q_new) = self.unmapped_neighbor_counts_query(q_idx);
            let (t_term, t_new) = self.unmapped_neighbor_counts_target(t_idx);

            q_term <= t_term && q_term + q_new <= t_term + t_new
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
}

use vf2::Vf2State;

mod ullmann {
    //! Ullmann candidate-matrix refinement and backtracking.
    //!
    //! See [Ullmann (1976)](https://doi.org/10.1145/321921.321925).

    use super::*;

    /// Ullmann candidate matrix `m[i * n2 + j]`: query node `i` may map to target
    /// node `j` when labels are compatible and `deg(i) <= deg(j)`.
    pub(super) fn ullmann_matrix(
        query: &Graph,
        target: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
    ) -> Vec<bool> {
        let n1 = query.node_count();
        let n2 = target.node_count();
        let target_deg: Vec<usize> = (0..n2)
            .map(|j| target.neighbors(NodeId(j as u32)).len())
            .collect();
        let mut m = vec![false; n1 * n2];
        for i in 0..n1 {
            let di = query.neighbors(NodeId(i as u32)).len();
            for j in 0..n2 {
                if di <= target_deg[j] && node_match(NodeId(i as u32), NodeId(j as u32)) {
                    m[i * n2 + j] = true;
                }
            }
        }
        m
    }

    /// Ullmann refinement to a fixpoint: clear `m[i][j]` unless every query neighbor
    /// `x` of `i` has a target neighbor `y` of `j` with `m[x][y]` and a matching edge.
    pub(super) fn ullmann_refine(
        query: &Graph,
        target: &Graph,
        m: &mut [bool],
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) {
        let n2 = target.node_count();
        let mut changed = true;
        while changed {
            changed = false;
            for i in 0..query.node_count() {
                for j in 0..n2 {
                    if !m[i * n2 + j] {
                        continue;
                    }
                    let supported = query.neighbors(NodeId(i as u32)).iter().all(|qn| {
                        let x = qn.node.index();
                        target
                            .neighbors(NodeId(j as u32))
                            .iter()
                            .any(|tn| m[x * n2 + tn.node.index()] && edge_match(qn.edge, tn.edge))
                    });
                    if !supported {
                        m[i * n2 + j] = false;
                        changed = true;
                    }
                }
            }
        }
    }

    /// Row-by-row backtracking over the refined matrix. Query nodes are assigned in
    /// index order (`depth` = next query node); each complete injective,
    /// edge-consistent mapping is recorded.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn ullmann_search(
        query: &Graph,
        target: &Graph,
        depth: usize,
        m: &[bool],
        mapping: &mut Vec<usize>,
        used: &mut [bool],
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        results: &mut Vec<Vec<usize>>,
    ) {
        let n1 = query.node_count();
        let n2 = target.node_count();
        if depth == n1 {
            results.push(mapping.clone());
            return;
        }
        let i = depth;
        for j in 0..n2 {
            if !m[i * n2 + j] || used[j] {
                continue;
            }
            // Edges from i to already-assigned query nodes must map to matching target edges.
            let consistent = query.neighbors(NodeId(i as u32)).iter().all(|qn| {
                let a = qn.node.index();
                if a >= depth {
                    return true;
                }
                match target.find_edge(NodeId(j as u32), NodeId(mapping[a] as u32)) {
                    Some(f) => edge_match(qn.edge, f),
                    None => false,
                }
            });
            if !consistent {
                continue;
            }
            let mut m2 = m.to_vec();
            for k in 0..n2 {
                if k != j {
                    m2[i * n2 + k] = false;
                }
            }
            for r in 0..n1 {
                if r != i {
                    m2[r * n2 + j] = false;
                }
            }
            ullmann_refine(query, target, &mut m2, edge_match);
            let future_ok = (depth + 1..n1).all(|r| (0..n2).any(|k| m2[r * n2 + k]));
            if future_ok {
                mapping[i] = j;
                used[j] = true;
                ullmann_search(
                    query,
                    target,
                    depth + 1,
                    &m2,
                    mapping,
                    used,
                    edge_match,
                    results,
                );
                mapping[i] = usize::MAX;
                used[j] = false;
            }
        }
    }
}

use ullmann::{ullmann_matrix, ullmann_refine, ullmann_search};

mod ri {
    //! RI static query ordering and backtracking.
    //!
    //! See [Bonnici et al.
    //! (2013)](https://doi.org/10.1186/1471-2105-14-S7-S13).

    use super::*;

    /// RI GreatestConstraintFirst ordering of the query vertices: the root is the
    /// max-degree vertex (or `first`, for an anchored search); each next vertex
    /// maximizes, lexicographically, `(V_m, V_n, V_o)` — its counts of neighbors that
    /// are already ordered, adjacent to an ordered vertex, or neither.
    pub(super) fn ri_order(query: &Graph, first: Option<usize>) -> Vec<usize> {
        let n = query.node_count();
        let mut order = Vec::with_capacity(n);
        if n == 0 {
            return order;
        }
        let mut ordered = vec![false; n];
        let mut ordered_neighbors = vec![0usize; n];
        let place = |v: usize, ordered: &mut [bool], ordered_neighbors: &mut [usize]| {
            ordered[v] = true;
            for nb in query.neighbors(NodeId(v as u32)).iter() {
                ordered_neighbors[nb.node.index()] += 1;
            }
        };
        let start = first.unwrap_or_else(|| {
            (0..n)
                .max_by_key(|&u| query.neighbors(NodeId(u as u32)).len())
                .expect("non-empty")
        });
        order.push(start);
        place(start, &mut ordered, &mut ordered_neighbors);
        while order.len() < n {
            let mut best: Option<((usize, usize, usize), usize)> = None;
            for u in 0..n {
                if ordered[u] {
                    continue;
                }
                let (mut vm, mut vn, mut vo) = (0usize, 0usize, 0usize);
                for nb in query.neighbors(NodeId(u as u32)).iter() {
                    let x = nb.node.index();
                    if ordered[x] {
                        vm += 1;
                    } else if ordered_neighbors[x] > 0 {
                        vn += 1;
                    } else {
                        vo += 1;
                    }
                }
                let key = (vm, vn, vo);
                if best.is_none_or(|(bk, _)| key > bk) {
                    best = Some((key, u));
                }
            }
            let next = best.expect("an unordered vertex remains").1;
            order.push(next);
            place(next, &mut ordered, &mut ordered_neighbors);
        }
        order
    }

    /// For each ordering position, the earlier-ordered query neighbors as
    /// `(their position, the query edge to them)` — the matcher's edge-consistency set.
    pub(super) fn ri_parents(query: &Graph, order: &[usize]) -> Vec<Vec<(usize, EdgeId)>> {
        let n = order.len();
        let mut position = vec![0usize; n];
        for (p, &u) in order.iter().enumerate() {
            position[u] = p;
        }
        let mut parents = vec![Vec::new(); n];
        for (d, &u) in order.iter().enumerate() {
            for nb in query.neighbors(NodeId(u as u32)).iter() {
                let p = position[nb.node.index()];
                if p < d {
                    parents[d].push((p, nb.edge));
                }
            }
        }
        parents
    }

    /// Backtracking in RI order. Candidates for the next vertex come from the
    /// target-neighbors of an already-mapped ordered neighbor's image (all targets at
    /// a component root); edges to every ordered neighbor must match.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn ri_search(
        target: &Graph,
        order: &[usize],
        parents: &[Vec<(usize, EdgeId)>],
        depth: usize,
        mapping: &mut [usize],
        used: &mut [bool],
        forced_root: Option<usize>,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        results: &mut Vec<Vec<usize>>,
    ) {
        let n1 = order.len();
        let n2 = target.node_count();
        if depth == n1 {
            let mut by_query = vec![0usize; n1];
            for (p, &u) in order.iter().enumerate() {
                by_query[u] = mapping[p];
            }
            results.push(by_query);
            return;
        }
        let u = order[depth];
        let candidates: Vec<usize> = match parents[depth].first() {
            _ if depth == 0 => match forced_root {
                Some(t) => vec![t],
                None => (0..n2).collect(),
            },
            Some(&(parent, _)) => target
                .neighbors(NodeId(mapping[parent] as u32))
                .iter()
                .map(|tn| tn.node.index())
                .collect(),
            None => (0..n2).collect(),
        };
        for v in candidates {
            if used[v] || !node_match(NodeId(u as u32), NodeId(v as u32)) {
                continue;
            }
            let consistent = parents[depth].iter().all(|&(p, e_q)| {
                match target.find_edge(NodeId(v as u32), NodeId(mapping[p] as u32)) {
                    Some(f) => edge_match(e_q, f),
                    None => false,
                }
            });
            if !consistent {
                continue;
            }
            mapping[depth] = v;
            used[v] = true;
            ri_search(
                target,
                order,
                parents,
                depth + 1,
                mapping,
                used,
                forced_root,
                node_match,
                edge_match,
                results,
            );
            used[v] = false;
        }
    }
}

use ri::{ri_order, ri_parents, ri_search};

mod arc_match {
    //! ArcMatch vertex- and edge-domain reduction, ordering, and search.
    //!
    //! See [Bonnici et al.
    //! (2024)](https://doi.org/10.1007/s10618-024-01061-8).

    use super::*;

    /// ArcMatch initial vertex domains `d[qi * n2 + tj]`: target vertex `tj`
    /// is in `D(qi)` when labels are compatible and `deg(qi) <= deg(tj)`.
    pub(super) fn arcmatch_vertex_domains(
        query: &Graph,
        target: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
    ) -> Vec<bool> {
        let n1 = query.node_count();
        let n2 = target.node_count();
        let target_deg: Vec<usize> = (0..n2)
            .map(|j| target.neighbors(NodeId(j as u32)).len())
            .collect();
        let mut d = vec![false; n1 * n2];
        for qi in 0..n1 {
            let dq = query.neighbors(NodeId(qi as u32)).len();
            for tj in 0..n2 {
                if dq <= target_deg[tj] && node_match(NodeId(qi as u32), NodeId(tj as u32)) {
                    d[qi * n2 + tj] = true;
                }
            }
        }
        d
    }

    /// ArcMatch arc consistency (Bonnici 2024 Algorithm 1) to a fixpoint:
    /// keep `t` in `D(a)` only if every query edge `{a, b}` has a label-matching target
    /// edge from `t` into `D(b)`. Undirected, so both endpoints of each query edge are
    /// reduced.
    pub(super) fn arcmatch_arc_consistency(
        query: &Graph,
        target: &Graph,
        d: &mut [bool],
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) {
        let n2 = target.node_count();
        let mut reduced = true;
        while reduced {
            reduced = false;
            for qi in 0..query.node_count() {
                for nb in query.neighbors(NodeId(qi as u32)).iter() {
                    let qj = nb.node.index();
                    if qj <= qi {
                        continue; // visit each undirected query edge once
                    }
                    reduced |= arcmatch_reduce_side(target, d, n2, qi, qj, nb.edge, edge_match);
                    reduced |= arcmatch_reduce_side(target, d, n2, qj, qi, nb.edge, edge_match);
                }
            }
        }
    }

    /// Remove from `D(a)` every target vertex lacking a label-matching target edge into
    /// `D(b)`; returns whether anything was removed.
    fn arcmatch_reduce_side(
        target: &Graph,
        d: &mut [bool],
        n2: usize,
        a: usize,
        b: usize,
        e_q: EdgeId,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) -> bool {
        let mut changed = false;
        for ta in 0..n2 {
            if !d[a * n2 + ta] {
                continue;
            }
            let supported = target
                .neighbors(NodeId(ta as u32))
                .iter()
                .any(|tn| d[b * n2 + tn.node.index()] && edge_match(e_q, tn.edge));
            if !supported {
                d[a * n2 + ta] = false;
                changed = true;
            }
        }
        changed
    }

    /// ArcMatch domain graph: post-arc-consistency vertex domains plus, per
    /// directed query-node pair `(a, b)`, the compatible target edges `(x, y)` with
    /// `x in D(a)`, `y in D(b)`, `{x, y}` a target edge, and matching labels. Both
    /// orientations of each query edge are stored so the matcher can pull candidates
    /// from either endpoint.
    pub(super) struct ArcMatchDomains {
        pub(super) n2: usize,
        pub(super) vertex: Vec<bool>,
        pub(super) edge: HashMap<(usize, usize), Vec<(usize, usize)>>,
    }

    /// ArcMatch, build the edge domains from the (reduced) vertex domains.
    pub(super) fn arcmatch_edge_domains(
        query: &Graph,
        target: &Graph,
        vertex: Vec<bool>,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) -> ArcMatchDomains {
        let n2 = target.node_count();
        let mut edge: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();
        for qi in 0..query.node_count() {
            for nb in query.neighbors(NodeId(qi as u32)).iter() {
                let qj = nb.node.index();
                if qj <= qi {
                    continue; // visit each undirected query edge once
                }
                let e_q = nb.edge;
                let mut fwd = Vec::new();
                for x in 0..n2 {
                    if !vertex[qi * n2 + x] {
                        continue;
                    }
                    for tn in target.neighbors(NodeId(x as u32)).iter() {
                        let y = tn.node.index();
                        if vertex[qj * n2 + y] && edge_match(e_q, tn.edge) {
                            fwd.push((x, y));
                        }
                    }
                }
                let rev: Vec<(usize, usize)> = fwd.iter().map(|&(x, y)| (y, x)).collect();
                edge.insert((qi, qj), fwd);
                edge.insert((qj, qi), rev);
            }
        }
        ArcMatchDomains { n2, vertex, edge }
    }

    /// Each undirected query edge once, as `(a, b)` with `a < b`.
    fn arcmatch_query_edges(query: &Graph) -> Vec<(usize, usize)> {
        let mut edges = Vec::new();
        for a in 0..query.node_count() {
            for nb in query.neighbors(NodeId(a as u32)).iter() {
                let b = nb.node.index();
                if a < b {
                    edges.push((a, b));
                }
            }
        }
        edges
    }

    /// ArcMatch path-based reduction (Bonnici 2024 §3.1). Runs the query-path
    /// DFS from every query vertex; provably safe (Thm 2–3: never discards a target
    /// vertex/edge that belongs to a real match). `lp` is the max path length in
    /// vertices (paper default 6).
    pub(super) fn arcmatch_path_reduction(query: &Graph, domains: &mut ArcMatchDomains, lp: usize) {
        for source in 0..query.node_count() {
            let mut omega = vec![source];
            arcmatch_path_dfs(query, &mut omega, lp, domains);
        }
    }

    /// Algorithm 2 — DFS over the query graph extracting maximal paths (`lp` vertices),
    /// rings (back to the source), and dead-ends; each (≥ 3 vertices — arc consistency
    /// already covers edges) is verified. `omega` is the path stack.
    fn arcmatch_path_dfs(
        query: &Graph,
        omega: &mut Vec<usize>,
        lp: usize,
        domains: &mut ArcMatchDomains,
    ) {
        let top = *omega.last().expect("non-empty path");
        let source = omega[0];
        let mut extended = false;
        let neighbors: Vec<usize> = query
            .neighbors(NodeId(top as u32))
            .iter()
            .map(|n| n.node.index())
            .collect();
        for v in neighbors {
            if v == source {
                if omega.len() >= 3 {
                    omega.push(v); // close the ring: source repeated at the tail
                    arcmatch_verify_path(query, omega, true, domains);
                    omega.pop();
                }
            } else if !omega.contains(&v) {
                extended = true;
                omega.push(v);
                if omega.len() == lp {
                    arcmatch_verify_path(query, omega, false, domains);
                } else {
                    arcmatch_path_dfs(query, omega, lp, domains);
                }
                omega.pop();
            }
        }
        if !extended && omega.len() >= 3 {
            arcmatch_verify_path(query, omega, false, domains);
        }
    }

    /// Algorithm 3 — discard from the source domain each target vertex that cannot
    /// reproduce the query path `omega` through the edge domains, then propagate.
    fn arcmatch_verify_path(
        query: &Graph,
        omega: &[usize],
        is_ring: bool,
        domains: &mut ArcMatchDomains,
    ) {
        let source = omega[0];
        let n2 = domains.n2;
        let candidates: Vec<usize> = (0..n2)
            .filter(|&t| domains.vertex[source * n2 + t])
            .collect();
        let mut reduced = false;
        for t in candidates {
            let mut omega_hat = vec![t];
            if !arcmatch_verify_path_dfs(omega, &mut omega_hat, is_ring, domains) {
                domains.vertex[source * n2 + t] = false;
                reduced = true;
            }
        }
        if reduced {
            arcmatch_refine_domains(query, source, domains);
        }
    }

    /// Algorithm 4 — does a target path `omega_hat` matching the query path `omega`
    /// exist, navigating the edge domains and avoiding revisits? For a ring, the final
    /// step must close back to `omega_hat[0]`.
    fn arcmatch_verify_path_dfs(
        omega: &[usize],
        omega_hat: &mut Vec<usize>,
        is_ring: bool,
        domains: &ArcMatchDomains,
    ) -> bool {
        if omega_hat.len() == omega.len() {
            return true;
        }
        let i = omega_hat.len();
        let (uq, vq) = (omega[i - 1], omega[i]);
        let current = omega_hat[i - 1];
        let Some(pairs) = domains.edge.get(&(uq, vq)) else {
            return false;
        };
        // TODO(perf): this clone (per DFS step, scales with target size) only sidesteps a
        // shared reborrow across the recursive call; removable. Defer to the stage-E benches.
        let pairs = pairs.clone();
        for (ut, vt) in pairs {
            if ut != current {
                continue;
            }
            if is_ring && omega_hat.len() == omega.len() - 1 {
                if vt == omega_hat[0] {
                    return true;
                }
            } else if !omega_hat.contains(&vt) {
                omega_hat.push(vt);
                if arcmatch_verify_path_dfs(omega, omega_hat, is_ring, domains) {
                    return true;
                }
                omega_hat.pop();
            }
        }
        false
    }

    /// Algorithm 5 — after the source domain shrank, drop edge-domain pairs whose
    /// endpoints left their vertex domains and vertex-domain entries with no supporting
    /// edge, to a fixpoint.
    fn arcmatch_refine_domains(query: &Graph, s: usize, domains: &mut ArcMatchDomains) {
        let n2 = domains.n2;
        for nb in arcmatch_query_edges(query) {
            let (a, b) = nb;
            if a == s || b == s {
                arcmatch_drop_pairs_both(domains, n2, a, b);
            }
        }
        let edges = arcmatch_query_edges(query);
        let mut reduced = true;
        while reduced {
            reduced = false;
            for &(a, b) in &edges {
                reduced |= arcmatch_drop_pairs_both(domains, n2, a, b);
            }
            for &(a, b) in &edges {
                reduced |= arcmatch_drop_unsupported(domains, n2, a, b);
                reduced |= arcmatch_drop_unsupported(domains, n2, b, a);
            }
        }
    }

    /// Keep in `D({a,b})` only pairs whose endpoints are still in both vertex domains;
    /// keeps the two directed orientations in sync. Returns whether anything changed.
    fn arcmatch_drop_pairs_both(
        domains: &mut ArcMatchDomains,
        n2: usize,
        a: usize,
        b: usize,
    ) -> bool {
        let before = domains.edge[&(a, b)].len();
        let mut kept = Vec::with_capacity(before);
        for &(x, y) in &domains.edge[&(a, b)] {
            if domains.vertex[a * n2 + x] && domains.vertex[b * n2 + y] {
                kept.push((x, y));
            }
        }
        if kept.len() == before {
            return false;
        }
        let rev: Vec<(usize, usize)> = kept.iter().map(|&(x, y)| (y, x)).collect();
        domains.edge.insert((a, b), kept);
        domains.edge.insert((b, a), rev);
        true
    }

    /// Remove from `D(a)` each target vertex with no supporting edge-domain pair (no
    /// `(a-side, _)` pair). Returns whether anything changed.
    fn arcmatch_drop_unsupported(
        domains: &mut ArcMatchDomains,
        n2: usize,
        a: usize,
        b: usize,
    ) -> bool {
        let mut changed = false;
        for t in 0..n2 {
            if !domains.vertex[a * n2 + t] {
                continue;
            }
            if !domains.edge[&(a, b)].iter().any(|&(x, _)| x == t) {
                domains.vertex[a * n2 + t] = false;
                changed = true;
            }
        }
        changed
    }

    /// ArcMatch query-vertex ordering (Bonnici 2024 §3.2). Singleton-domain
    /// vertices first, then a greedy neighbor-expansion driven by five measures, then
    /// degree-1 ("peripheral", §3.2.1) vertices last (largest domain first). Each step
    /// maximizes, lexicographically, the candidate's count of neighbors that are
    /// already ordered, that border the ordered set, and that are unrelated to it, then
    /// its degree, then minimizes its domain size, then its id.
    pub(super) fn arcmatch_variable_ordering(
        query: &Graph,
        domains: &ArcMatchDomains,
    ) -> Vec<usize> {
        let n = query.node_count();
        let n2 = domains.n2;
        let domain_size: Vec<usize> = (0..n)
            .map(|v| (0..n2).filter(|&t| domains.vertex[v * n2 + t]).count())
            .collect();
        let singleton: Vec<bool> = domain_size.iter().map(|&s| s == 1).collect();
        let peripheral: Vec<bool> = (0..n)
            .map(|v| query.neighbors(NodeId(v as u32)).len() == 1 && !singleton[v])
            .collect();

        let mut order = Vec::with_capacity(n);
        let mut placed = vec![false; n];

        arcmatch_order_phase(query, &domain_size, &mut order, &mut placed, |v| {
            singleton[v]
        });
        arcmatch_order_phase(query, &domain_size, &mut order, &mut placed, |v| {
            !singleton[v] && !peripheral[v]
        });

        let mut peris: Vec<usize> = (0..n).filter(|&v| peripheral[v] && !placed[v]).collect();
        peris.sort_by(|&a, &b| domain_size[b].cmp(&domain_size[a]).then(a.cmp(&b)));
        order.extend(peris);
        order
    }

    /// Greedy ordering of the `in_phase` vertices: repeatedly append the best candidate
    /// adjacent to the already-ordered set (or, when none, the best remaining vertex —
    /// a new component root).
    fn arcmatch_order_phase(
        query: &Graph,
        domain_size: &[usize],
        order: &mut Vec<usize>,
        placed: &mut [bool],
        in_phase: impl Fn(usize) -> bool,
    ) {
        let n = query.node_count();
        loop {
            let remaining: Vec<usize> = (0..n).filter(|&v| in_phase(v) && !placed[v]).collect();
            if remaining.is_empty() {
                break;
            }
            let frontier: Vec<usize> = remaining
                .iter()
                .copied()
                .filter(|&v| {
                    query
                        .neighbors(NodeId(v as u32))
                        .iter()
                        .any(|nb| placed[nb.node.index()])
                })
                .collect();
            let candidates = if frontier.is_empty() {
                &remaining
            } else {
                &frontier
            };
            let best = *candidates
                .iter()
                .min_by_key(|&&v| {
                    let (n1, n2, n3) = arcmatch_neighbor_split(query, placed, v);
                    (
                        Reverse(n1),
                        Reverse(n2),
                        Reverse(n3),
                        Reverse(query.neighbors(NodeId(v as u32)).len()),
                        domain_size[v],
                        v,
                    )
                })
                .expect("non-empty candidates");
            order.push(best);
            placed[best] = true;
        }
    }

    /// Split `v`'s neighbors into counts of (already-ordered, bordering-the-ordered-set,
    /// unrelated) vertices — the first three ordering measures.
    fn arcmatch_neighbor_split(query: &Graph, placed: &[bool], v: usize) -> (usize, usize, usize) {
        let (mut ordered, mut bordering, mut unrelated) = (0, 0, 0);
        for nb in query.neighbors(NodeId(v as u32)).iter() {
            let u = nb.node.index();
            if placed[u] {
                ordered += 1;
            } else if query
                .neighbors(NodeId(u as u32))
                .iter()
                .any(|w| placed[w.node.index()])
            {
                bordering += 1;
            } else {
                unrelated += 1;
            }
        }
        (ordered, bordering, unrelated)
    }

    /// ArcMatch dynamic-parent backtracking (Bonnici 2024 Algorithm 6).
    /// Places query vertices in `order`; for each, picks the already-placed neighbor
    /// (the dynamic parent) with the smallest edge domain, draws candidates from that
    /// edge domain, and keeps those consistent with every placed neighbor. A query
    /// vertex with no placed neighbor (a component root) is seeded from its vertex
    /// domain. `assigned` is indexed by query vertex.
    pub(super) fn arcmatch_search(
        query: &Graph,
        order: &[usize],
        domains: &ArcMatchDomains,
        depth: usize,
        assigned: &mut [Option<usize>],
        used: &mut [bool],
        results: &mut Vec<Vec<usize>>,
    ) {
        if depth == order.len() {
            results.push(
                assigned
                    .iter()
                    .map(|t| t.expect("complete match"))
                    .collect(),
            );
            return;
        }
        let q = order[depth];
        let n2 = domains.n2;
        let parents: Vec<usize> = query
            .neighbors(NodeId(q as u32))
            .iter()
            .map(|nb| nb.node.index())
            .filter(|&u| assigned[u].is_some())
            .collect();

        if parents.is_empty() {
            for t in 0..n2 {
                if domains.vertex[q * n2 + t] && !used[t] {
                    assigned[q] = Some(t);
                    used[t] = true;
                    arcmatch_search(query, order, domains, depth + 1, assigned, used, results);
                    used[t] = false;
                    assigned[q] = None;
                }
            }
            return;
        }

        let s = *parents
            .iter()
            .min_by_key(|&&p| domains.edge[&(p, q)].len())
            .expect("non-empty parents");
        let ms = assigned[s].expect("parent assigned");
        for &(u_t, v_t) in &domains.edge[&(s, q)] {
            if u_t != ms || used[v_t] {
                continue;
            }
            let feasible = parents.iter().all(|&p| {
                domains.edge[&(q, p)].contains(&(v_t, assigned[p].expect("parent assigned")))
            });
            if !feasible {
                continue;
            }
            assigned[q] = Some(v_t);
            used[v_t] = true;
            arcmatch_search(query, order, domains, depth + 1, assigned, used, results);
            used[v_t] = false;
            assigned[q] = None;
        }
    }
}

#[cfg(test)]
use arc_match::ArcMatchDomains;
use arc_match::{
    arcmatch_arc_consistency, arcmatch_edge_domains, arcmatch_path_reduction, arcmatch_search,
    arcmatch_variable_ordering, arcmatch_vertex_domains,
};

mod vf2_rdkit {
    //! RDKit-compatible, vflib-derived VF2 search with Mayfield's
    //! chemical-graph candidate restrictions.
    //!
    //! See [RDKit PR #2500](https://github.com/rdkit/rdkit/pull/2500).

    use super::*;

    /// Next query atom for `Vf2Rdkit`: lowest-index unmapped atom adjacent to the
    /// mapped set (terminal), else the lowest-index unmapped atom (a component root).
    fn vf2rdkit_next(query: &Graph, mapping: &[Option<usize>]) -> Option<usize> {
        let mut fallback = None;
        for q in 0..query.node_count() {
            if mapping[q].is_some() {
                continue;
            }
            if fallback.is_none() {
                fallback = Some(q);
            }
            if query
                .neighbors(NodeId(q as u32))
                .iter()
                .any(|nb| mapping[nb.node.index()].is_some())
            {
                return Some(q);
            }
        }
        fallback
    }

    /// `Vf2Rdkit` recursive search: candidates come from a mapped neighbor's image
    /// adjacency (else all atoms for a root), filtered by degree bound, node match, and
    /// edge consistency with every already-mapped query neighbor.
    pub(super) fn vf2rdkit_search(
        query: &Graph,
        target: &Graph,
        mapping: &mut [Option<usize>],
        used: &mut [bool],
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        results: &mut Vec<Vec<usize>>,
    ) {
        let Some(q) = vf2rdkit_next(query, mapping) else {
            results.push(mapping.iter().map(|m| m.expect("complete")).collect());
            return;
        };
        let q_degree = query.neighbors(NodeId(q as u32)).len();
        let parent = query
            .neighbors(NodeId(q as u32))
            .iter()
            .find_map(|nb| mapping[nb.node.index()]);
        let candidates: Vec<usize> = match parent {
            Some(image) => target
                .neighbors(NodeId(image as u32))
                .iter()
                .map(|nb| nb.node.index())
                .collect(),
            None => (0..target.node_count()).collect(),
        };
        for t in candidates {
            if used[t]
                || q_degree > target.neighbors(NodeId(t as u32)).len()
                || !node_match(NodeId(q as u32), NodeId(t as u32))
            {
                continue;
            }
            let consistent = query.neighbors(NodeId(q as u32)).iter().all(|nb| {
                let Some(image) = mapping[nb.node.index()] else {
                    return true;
                };
                match target.find_edge(NodeId(t as u32), NodeId(image as u32)) {
                    Some(te) => edge_match(nb.edge, te),
                    None => false,
                }
            });
            if !consistent {
                continue;
            }
            mapping[q] = Some(t);
            used[t] = true;
            vf2rdkit_search(
                query, target, mapping, used, node_match, edge_match, results,
            );
            used[t] = false;
            mapping[q] = None;
        }
    }
}

use vf2_rdkit::vf2rdkit_search;

mod ray_kirsch {
    //! Ray--Kirsch bond-oriented backtracking.
    //!
    //! See [Ray and Kirsch
    //! (1957)](https://doi.org/10.1126/science.126.3278.814) and Open Babel's
    //! [`parsmart.cpp`](https://github.com/openbabel/openbabel/blob/master/src/parsmart.cpp).

    use super::*;

    /// Query bond order for `RayKirsch`: a DFS edge ordering per component. Each bond is
    /// `(from, to, edge)` with `from` already on the DFS path, so a bond after a
    /// component's first shares an atom with an earlier one (tree edges then recurse,
    /// back edges close rings).
    pub(super) fn rk_bond_order(query: &Graph) -> Vec<(usize, usize, EdgeId)> {
        let mut visited = vec![false; query.node_count()];
        let mut emitted = vec![false; query.edge_count()];
        let mut order = Vec::new();
        for root in 0..query.node_count() {
            if !visited[root] {
                rk_dfs(query, root, &mut visited, &mut emitted, &mut order);
            }
        }
        order
    }

    fn rk_dfs(
        query: &Graph,
        u: usize,
        visited: &mut [bool],
        emitted: &mut [bool],
        order: &mut Vec<(usize, usize, EdgeId)>,
    ) {
        visited[u] = true;
        for nb in query.neighbors(NodeId(u as u32)).iter() {
            let v = nb.node.index();
            if !emitted[nb.edge.index()] {
                emitted[nb.edge.index()] = true;
                order.push((u, v, nb.edge));
            }
            if !visited[v] {
                rk_dfs(query, v, visited, emitted, order);
            }
        }
    }

    /// `RayKirsch` recursive search over `bond_order`. Each bond either closes a ring
    /// (both endpoints mapped), extends an unmapped endpoint along the mapped one's mol
    /// adjacency, or seeds a component root.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn rk_search(
        target: &Graph,
        bond_order: &[(usize, usize, EdgeId)],
        isolated: &[usize],
        i: usize,
        mapping: &mut [Option<usize>],
        used: &mut [bool],
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        results: &mut Vec<Vec<usize>>,
    ) {
        if i == bond_order.len() {
            rk_finish(target, isolated, 0, mapping, used, node_match, results);
            return;
        }
        let (qb, qe, qedge) = bond_order[i];
        match (mapping[qb], mapping[qe]) {
            (Some(mb), Some(me)) => {
                if let Some(te) = target.find_edge(NodeId(mb as u32), NodeId(me as u32)) {
                    if edge_match(qedge, te) {
                        rk_search(
                            target,
                            bond_order,
                            isolated,
                            i + 1,
                            mapping,
                            used,
                            node_match,
                            edge_match,
                            results,
                        );
                    }
                }
            }
            (Some(mb), None) => {
                rk_extend(
                    target, bond_order, isolated, i, qe, mb, qedge, mapping, used, node_match,
                    edge_match, results,
                );
            }
            (None, Some(me)) => {
                rk_extend(
                    target, bond_order, isolated, i, qb, me, qedge, mapping, used, node_match,
                    edge_match, results,
                );
            }
            (None, None) => {
                for atom in 0..target.node_count() {
                    if used[atom] || !node_match(NodeId(qb as u32), NodeId(atom as u32)) {
                        continue;
                    }
                    mapping[qb] = Some(atom);
                    used[atom] = true;
                    rk_search(
                        target, bond_order, isolated, i, mapping, used, node_match, edge_match,
                        results,
                    );
                    used[atom] = false;
                    mapping[qb] = None;
                }
            }
        }
    }

    /// Map the unmapped query endpoint `q_free` of bond `i` to an unused mol neighbor of
    /// the mapped endpoint's image `m_anchor`, then continue with the next bond.
    #[allow(clippy::too_many_arguments)]
    fn rk_extend(
        target: &Graph,
        bond_order: &[(usize, usize, EdgeId)],
        isolated: &[usize],
        i: usize,
        q_free: usize,
        m_anchor: usize,
        qedge: EdgeId,
        mapping: &mut [Option<usize>],
        used: &mut [bool],
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        results: &mut Vec<Vec<usize>>,
    ) {
        for k in 0..target.neighbors(NodeId(m_anchor as u32)).len() {
            let nb = target.neighbors(NodeId(m_anchor as u32))[k];
            let t = nb.node.index();
            if used[t]
                || !node_match(NodeId(q_free as u32), NodeId(t as u32))
                || !edge_match(qedge, nb.edge)
            {
                continue;
            }
            mapping[q_free] = Some(t);
            used[t] = true;
            rk_search(
                target,
                bond_order,
                isolated,
                i + 1,
                mapping,
                used,
                node_match,
                edge_match,
                results,
            );
            used[t] = false;
            mapping[q_free] = None;
        }
    }

    /// Map degree-0 query atoms (skipping any already pinned, e.g. an anchor) over the
    /// remaining unused mol atoms, recording each complete assignment.
    fn rk_finish(
        target: &Graph,
        isolated: &[usize],
        k: usize,
        mapping: &mut [Option<usize>],
        used: &mut [bool],
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        results: &mut Vec<Vec<usize>>,
    ) {
        if k == isolated.len() {
            results.push(mapping.iter().map(|m| m.expect("complete")).collect());
            return;
        }
        let q = isolated[k];
        if mapping[q].is_some() {
            rk_finish(target, isolated, k + 1, mapping, used, node_match, results);
            return;
        }
        for t in 0..target.node_count() {
            if used[t] || !node_match(NodeId(q as u32), NodeId(t as u32)) {
                continue;
            }
            mapping[q] = Some(t);
            used[t] = true;
            rk_finish(target, isolated, k + 1, mapping, used, node_match, results);
            used[t] = false;
            mapping[q] = None;
        }
    }
}

use ray_kirsch::{rk_bond_order, rk_search};

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::SubgraphIsomorphismAlgorithm::{ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit};
    use super::{ARCMATCH_DEFAULT_PATH_LENGTH, *};

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
    // Label maps for the disconnected-query regression: query labels [0,0,1] /
    // edge [1]; target labels [0,1,0] / edges [0,1].
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
    // Regression: Vf2's look-ahead under-reported disconnected queries whose isolated
    // node is mapped first (it created target terminals with no query counterpart, so
    // the now-removed `new <= new` cut spuriously pruned the only match).
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
    fn test_enumerate_subgraph_isomorphisms(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        for alg in [
            Vf2,
            Ullmann,
            Ri,
            ArcMatch {
                path_length: ARCMATCH_DEFAULT_PATH_LENGTH,
            },
            Vf2Rdkit,
            RayKirsch,
        ] {
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
    #[case::no_adjacency(Graph::new(3, &[[1, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(0)), any_node, any_edge, vec![])]
    #[case::empty_query(Graph::new(3, &[[0, 1], [1, 2]]), Graph::default(), (NodeId(0), NodeId(0)), any_node, any_edge, vec![])]
    #[case::out_of_bounds(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(5), NodeId(0)), any_node, any_edge, vec![])]
    #[case::path_anchored(Graph::new(5, &[[0, 1], [1, 2], [2, 3], [3, 4]]), Graph::new(3, &[[0, 1], [1, 2]]), (NodeId(1), NodeId(2)), any_node, any_edge,
        vec![vec![1, 2, 3], vec![3, 2, 1]])]
    #[case::node_filter(Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(1)), exclude_t0, any_edge, vec![vec![1, 2]])]
    #[case::node_filter_rejects_anchor(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(1)), reject_q0_t1, any_edge, vec![])]
    #[case::edge_filter(Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), Graph::new(2, &[[0, 1]]), (NodeId(0), NodeId(1)), any_node, only_tedge1, vec![vec![1, 2]])]
    fn test_enumerate_subgraph_isomorphisms_at(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] anchor: (NodeId, NodeId),
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        for alg in [
            Vf2,
            Ullmann,
            Ri,
            ArcMatch {
                path_length: ARCMATCH_DEFAULT_PATH_LENGTH,
            },
            Vf2Rdkit,
            RayKirsch,
        ] {
            let mut r: Vec<Vec<usize>> = target
                .enumerate_subgraph_isomorphisms_at(&query, anchor, &mut node_match, &mut edge_match, alg)
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

    // The path-length knob is a pruning/cost trade-off; the match set is invariant
    // to it (`< 3` skips path reduction, `>= 3` reduces over rings/paths).
    #[rstest]
    #[case(0)]
    #[case(2)]
    #[case(3)]
    #[case(6)]
    fn test_subgraph_isomorphisms_arcmatch_path_length(#[case] path_length: usize) {
        let graph = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let mut node_match = any_node;
        let mut edge_match = any_edge;
        let mut r: Vec<Vec<usize>> = graph
            .enumerate_subgraph_isomorphisms(
                &graph,
                &mut node_match,
                &mut edge_match,
                ArcMatch { path_length },
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
        assert_eq!(
            r,
            vec![
                vec![0, 1, 2],
                vec![0, 2, 1],
                vec![1, 0, 2],
                vec![1, 2, 0],
                vec![2, 0, 1],
                vec![2, 1, 0],
            ]
        );
    }

    #[rstest]
    fn test_arcmatch_arc_consistency() {
        // query: path 0-1-2; target: star, center 0 with leaves 1,2,3.
        let query = Graph::new(3, &[[0, 1], [1, 2]]);
        let target = Graph::new(4, &[[0, 1], [0, 2], [0, 3]]);
        let mut nm: fn(NodeId, NodeId) -> bool = any_node;
        let mut em: fn(EdgeId, EdgeId) -> bool = any_edge;
        let mut d = arcmatch_vertex_domains(&query, &target, &mut nm);
        arcmatch_arc_consistency(&query, &target, &mut d, &mut em);
        // Degree-1 path endpoints cannot map to the degree-3 center; the degree-2
        // middle must.
        assert_eq!(
            d,
            vec![
                false, true, true, true, // D(q0) = {t1, t2, t3}
                true, false, false, false, // D(q1) = {t0}
                false, true, true, true, // D(q2) = {t1, t2, t3}
            ]
        );
    }

    #[rstest]
    fn test_arcmatch_edge_domains() {
        // query: single edge 0-1; target: path 0-1-2. Degree-1 query domains admit
        // all target vertices, so the edge domain is every directed target edge.
        let query = Graph::new(2, &[[0, 1]]);
        let target = Graph::new(3, &[[0, 1], [1, 2]]);
        let mut nm: fn(NodeId, NodeId) -> bool = any_node;
        let mut em: fn(EdgeId, EdgeId) -> bool = any_edge;
        let vertex = arcmatch_vertex_domains(&query, &target, &mut nm);
        let domains = arcmatch_edge_domains(&query, &target, vertex, &mut em);
        let mut fwd = domains.edge[&(0, 1)].clone();
        fwd.sort();
        assert_eq!(fwd, vec![(0, 1), (1, 0), (1, 2), (2, 1)]);
        let mut rev = domains.edge[&(1, 0)].clone();
        rev.sort();
        assert_eq!(rev, vec![(0, 1), (1, 0), (1, 2), (2, 1)]);
    }

    // Safety (Bonnici 2024 Thm 2-3): path-based reduction must never discard a
    // target vertex/edge that belongs to a real match. Ground truth = Vf2.
    #[rstest]
    #[case::path(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(5, &[[0, 1], [1, 2], [2, 3], [3, 4]]))]
    #[case::triangle(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), Graph::new(4, &[[0, 1], [1, 2], [0, 2], [2, 3]]))]
    #[case::naphthalene(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        Graph::new(10, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [3, 6], [6, 7], [7, 8], [8, 9], [9, 4]])
    )]
    #[case::k4(
        Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]),
        Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
    )]
    fn test_arcmatch_path_reduction_safety(#[case] query: Graph, #[case] target: Graph) {
        let mut nm: fn(NodeId, NodeId) -> bool = any_node;
        let mut em: fn(EdgeId, EdgeId) -> bool = any_edge;
        let matches: Vec<Vec<usize>> = target
            .enumerate_subgraph_isomorphisms(&query, &mut nm, &mut em, Vf2)
            .iter()
            .map(|c| {
                c.matched_pairs()
                    .iter()
                    .map(|&(_, host)| host.index())
                    .collect()
            })
            .collect();
        assert!(!matches.is_empty(), "fixture should have matches");

        let mut vertex = arcmatch_vertex_domains(&query, &target, &mut nm);
        arcmatch_arc_consistency(&query, &target, &mut vertex, &mut em);
        let mut domains = arcmatch_edge_domains(&query, &target, vertex, &mut em);
        arcmatch_path_reduction(&query, &mut domains, 6);

        let n2 = target.node_count();
        for m in &matches {
            for (qi, &ti) in m.iter().enumerate() {
                assert!(
                    domains.vertex[qi * n2 + ti],
                    "match {m:?}: vertex q{qi}->t{ti} pruned"
                );
            }
            for qi in 0..query.node_count() {
                for nb in query.neighbors(NodeId(qi as u32)).iter() {
                    let qj = nb.node.index();
                    assert!(
                        domains.edge[&(qi, qj)].contains(&(m[qi], m[qj])),
                        "match {m:?}: edge q{qi}-q{qj} pruned"
                    );
                }
            }
        }
    }

    #[rstest]
    #[case::path_uniform(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        ArcMatchDomains { n2: 4, vertex: vec![true; 16], edge: HashMap::new() },
        vec![1, 2, 0, 3]
    )]
    #[case::star_center_first(
        Graph::new(4, &[[0, 1], [0, 2], [0, 3]]),
        ArcMatchDomains { n2: 4, vertex: vec![true; 16], edge: HashMap::new() },
        vec![0, 1, 2, 3]
    )]
    #[case::singleton_tail_first(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        ArcMatchDomains {
            n2: 4,
            vertex: vec![
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
                true, false, false, false,
            ],
            edge: HashMap::new(),
        },
        vec![3, 2, 1, 0]
    )]
    fn test_arcmatch_variable_ordering(
        #[case] query: Graph,
        #[case] domains: ArcMatchDomains,
        #[case] expected: Vec<usize>,
    ) {
        assert_eq!(arcmatch_variable_ordering(&query, &domains), expected);
    }
}
