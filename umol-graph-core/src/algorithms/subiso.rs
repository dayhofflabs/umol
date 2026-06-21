//! Subgraph isomorphism (monomorphism: query edges map to target edges, no
//! induced-subgraph reverse check). Multiple named algorithms behind a selector;
//! all return the same match set (a query→target index vector per occurrence).

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
}

impl Graph {
    // TODO: add singular `subgraph_isomorphism(...) -> Option<Vec<usize>>` that stops at
    // the first match (existence via `.is_some()`). Saves the embedding-multiplicity
    // factor on positive matches; same Vf2, `search` gains an early-exit cap.
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
            SubgraphIsomorphismAlgorithm::Ullmann => {
                self.subgraph_isomorphisms_ullmann(query, node_match, edge_match)
            }
            SubgraphIsomorphismAlgorithm::Ri => {
                self.subgraph_isomorphisms_ri(query, node_match, edge_match)
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
            SubgraphIsomorphismAlgorithm::Ullmann => {
                self.subgraph_isomorphisms_at_ullmann(query, anchor, node_match, edge_match)
            }
            SubgraphIsomorphismAlgorithm::Ri => {
                self.subgraph_isomorphisms_at_ri(query, anchor, node_match, edge_match)
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
        ullmann_search(query, self, 0, &m, &mut mapping, &mut used, edge_match, &mut results);
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
        ullmann_search(query, self, 0, &m, &mut mapping, &mut used, edge_match, &mut results);
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
            self, &order, &parents, 0, &mut mapping, &mut used, None, node_match, edge_match,
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

/// Ullmann candidate matrix `m[i * n2 + j]`: query node `i` may map to target
/// node `j` when labels are compatible and `deg(i) <= deg(j)`.
fn ullmann_matrix(
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
fn ullmann_refine(
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
fn ullmann_search(
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
            ullmann_search(query, target, depth + 1, &m2, mapping, used, edge_match, results);
            mapping[i] = usize::MAX;
            used[j] = false;
        }
    }
}

/// RI GreatestConstraintFirst ordering of the query vertices: the root is the
/// max-degree vertex (or `first`, for an anchored search); each next vertex
/// maximizes, lexicographically, `(V_m, V_n, V_o)` — its counts of neighbors that
/// are already ordered, adjacent to an ordered vertex, or neither.
fn ri_order(query: &Graph, first: Option<usize>) -> Vec<usize> {
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
fn ri_parents(query: &Graph, order: &[usize]) -> Vec<Vec<(usize, EdgeId)>> {
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
fn ri_search(
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
            target, order, parents, depth + 1, mapping, used, forced_root, node_match, edge_match,
            results,
        );
        used[v] = false;
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::SubgraphIsomorphismAlgorithm::{Ri, Ullmann, Vf2};
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
        for alg in [Vf2, Ullmann, Ri] {
            let mut r = target.subgraph_isomorphisms(&query, &mut node_match, &mut edge_match, alg);
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
    fn test_subgraph_isomorphisms_at(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] anchor: (NodeId, NodeId),
        #[case] mut node_match: fn(NodeId, NodeId) -> bool,
        #[case] mut edge_match: fn(EdgeId, EdgeId) -> bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        for alg in [Vf2, Ullmann, Ri] {
            let mut r = target
                .subgraph_isomorphisms_at(&query, anchor, &mut node_match, &mut edge_match, alg);
            r.sort();
            assert_eq!(r, expected, "algorithm {alg:?}");
        }
    }
}
