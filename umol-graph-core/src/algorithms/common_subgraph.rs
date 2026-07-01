//! Common subgraphs of two graphs. Three distinct tasks, each with its own entry point and
//! algorithm enum: the **maximum** (single largest) via McGregor, the **maximal**
//! (non-extendable) enumeration via Bron–Kerbosch, and the **complete** enumeration of *all*
//! common induced subgraphs by backtracking. All three treat cliques of the modular product
//! (Levi, Calcolo 9 (1973) 341–352) as the common induced subgraphs.
//!
//! Maximum (`Mcis`/`Mces`) — McGregor's vertex-mapping search (J. J. McGregor,
//! "Backtrack search algorithms and the maximal common subgraph problem", Software:
//! Practice and Experience 12 (1982) 23–34): extend a partial injective vertex
//! correspondence between two graphs, allowing vertices to be skipped, with
//! branch-and-bound on the objective. MCIS maximizes mapped vertices under an induced
//! (edge-iff-edge) constraint; MCES maximizes shared edges with no induced constraint.
//!
//! Maximal (`maximal_common_subgraphs`, `MaximalCommonSubgraphAlgorithm`) — every maximal
//! (non-extendable) common induced subgraph, as a maximal clique of the modular product.
//!
//! Complete (`enumerate_common_subgraphs`, `CommonSubgraphEnumerationAlgorithm`) — *every*
//! common induced subgraph, the empty one included, as *every* clique of the modular product.
//! Its result contains the maximal enumeration's; exponential in the worst case.

use bitvec::prelude::*;

use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum McsConnectivity {
    Connected,
    Disconnected,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum McisAlgorithm {
    McGregor,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum McesAlgorithm {
    McGregor,
}

/// Algorithm for the **complete** enumeration of all common induced subgraphs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CommonSubgraphEnumerationAlgorithm {
    /// Every clique of the modular product by recursive backtracking — each clique reported
    /// once, in ascending-vertex order, the empty clique included.
    Backtracking,
}

/// Algorithm for the **maximal** (non-extendable) common induced subgraphs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaximalCommonSubgraphAlgorithm {
    /// Bron–Kerbosch maximal-clique enumeration over the modular product (Bron &
    /// Kerbosch, CACM 16 (1973) 575–577).
    BronKerbosch,
}

/// One common subgraph: a vertex correspondence between the two graphs (sorted by
/// the first graph's node) and the number of edges it shares.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CommonSubgraph {
    mapping: Vec<(NodeId, NodeId)>,
    edge_count: usize,
}

impl CommonSubgraph {
    pub fn mapping(&self) -> &[(NodeId, NodeId)] {
        &self.mapping
    }

    pub fn node_count(&self) -> usize {
        self.mapping.len()
    }

    pub fn edge_count(&self) -> usize {
        self.edge_count
    }

    pub fn is_empty(&self) -> bool {
        self.mapping.is_empty()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum McsAlgorithmKind {
    Induced,
    Edge,
}

// Return one maximum or all maximum common subgraphs.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Enumerate {
    Single,
    All,
}

fn into_single(results: Vec<CommonSubgraph>) -> CommonSubgraph {
    results
        .into_iter()
        .next()
        .expect("a maximum common subgraph always exists (empty at worst)")
}

impl Graph {
    /// The maximum common induced subgraph (largest vertex set inducing isomorphic
    /// subgraphs). Always exists — empty at worst.
    pub fn maximum_common_induced_subgraph(
        &self,
        other: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        connectivity: McsConnectivity,
        alg: McisAlgorithm,
    ) -> CommonSubgraph {
        match alg {
            McisAlgorithm::McGregor => into_single(self.mcs_mcgregor(
                other,
                node_match,
                edge_match,
                &[],
                &[],
                connectivity,
                Enumerate::Single,
                McsAlgorithmKind::Induced,
            )),
        }
    }

    /// Every maximum common induced subgraph (all of the largest size).
    pub fn maximum_common_induced_subgraphs(
        &self,
        other: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        connectivity: McsConnectivity,
        alg: McisAlgorithm,
    ) -> Vec<CommonSubgraph> {
        match alg {
            McisAlgorithm::McGregor => self.mcs_mcgregor(
                other,
                node_match,
                edge_match,
                &[],
                &[],
                connectivity,
                Enumerate::All,
                McsAlgorithmKind::Induced,
            ),
        }
    }

    /// The maximum common edge subgraph (largest shared edge set). Always exists —
    /// empty at worst.
    pub fn maximum_common_edge_subgraph(
        &self,
        other: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        connectivity: McsConnectivity,
        alg: McesAlgorithm,
    ) -> CommonSubgraph {
        match alg {
            McesAlgorithm::McGregor => into_single(self.mcs_mcgregor(
                other,
                node_match,
                edge_match,
                &[],
                &[],
                connectivity,
                Enumerate::Single,
                McsAlgorithmKind::Edge,
            )),
        }
    }

    /// Every maximum common edge subgraph (all of the largest size).
    pub fn maximum_common_edge_subgraphs(
        &self,
        other: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        connectivity: McsConnectivity,
        alg: McesAlgorithm,
    ) -> Vec<CommonSubgraph> {
        match alg {
            McesAlgorithm::McGregor => self.mcs_mcgregor(
                other,
                node_match,
                edge_match,
                &[],
                &[],
                connectivity,
                Enumerate::All,
                McsAlgorithmKind::Edge,
            ),
        }
    }

    /// Seeded maximum common edge subgraph. `anchor` pairs are forced into the
    /// result and never skipped; `hint` pairs warm-start the incumbent bound (may be
    /// discarded). Either may be empty.
    #[allow(clippy::too_many_arguments)]
    pub fn maximum_common_edge_subgraph_seeded(
        &self,
        other: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        anchor: &[(NodeId, NodeId)],
        hint: &[(NodeId, NodeId)],
        connectivity: McsConnectivity,
        alg: McesAlgorithm,
    ) -> CommonSubgraph {
        match alg {
            McesAlgorithm::McGregor => into_single(self.mcs_mcgregor(
                other,
                node_match,
                edge_match,
                anchor,
                hint,
                connectivity,
                Enumerate::Single,
                McsAlgorithmKind::Edge,
            )),
        }
    }

    /// Every seeded maximum common edge subgraph. See
    /// [`Graph::maximum_common_edge_subgraph_seeded`] for `anchor`/`hint`.
    #[allow(clippy::too_many_arguments)]
    pub fn maximum_common_edge_subgraphs_seeded(
        &self,
        other: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        anchor: &[(NodeId, NodeId)],
        hint: &[(NodeId, NodeId)],
        connectivity: McsConnectivity,
        alg: McesAlgorithm,
    ) -> Vec<CommonSubgraph> {
        match alg {
            McesAlgorithm::McGregor => self.mcs_mcgregor(
                other,
                node_match,
                edge_match,
                anchor,
                hint,
                connectivity,
                Enumerate::All,
                McsAlgorithmKind::Edge,
            ),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn mcs_mcgregor<N, E>(
        &self,
        other: &Graph,
        node_match: &mut N,
        edge_match: &mut E,
        anchor: &[(NodeId, NodeId)],
        hint: &[(NodeId, NodeId)],
        connectivity: McsConnectivity,
        enumerate: Enumerate,
        kind: McsAlgorithmKind,
    ) -> Vec<CommonSubgraph>
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
    {
        let mut state = McsState::new(self, other, kind, connectivity, enumerate);
        if !state.place_anchor(anchor, node_match, edge_match) {
            return Vec::new();
        }
        state.seed_hint(hint, node_match, edge_match);
        state.search(node_match, edge_match);
        if state.best.is_empty() {
            let fallback = state.subgraph_from(&state.pairs, state.edges, edge_match);
            state.best.push(fallback);
        }
        let mut best = state.best;
        best.sort_by(|x, y| x.mapping.cmp(&y.mapping));
        best.dedup();
        best
    }
}

struct McsState<'g> {
    a: &'g Graph,
    b: &'g Graph,
    kind: McsAlgorithmKind,
    connectivity: McsConnectivity,
    enumerate: Enumerate,
    a_to_b: Vec<Option<NodeId>>,
    b_used: Vec<bool>,
    decided: Vec<bool>,
    anchored: Vec<bool>,
    pairs: Vec<(NodeId, NodeId)>,
    edges: usize,
    best_score: usize,
    best: Vec<CommonSubgraph>,
}

impl<'g> McsState<'g> {
    fn new(
        a: &'g Graph,
        b: &'g Graph,
        kind: McsAlgorithmKind,
        connectivity: McsConnectivity,
        enumerate: Enumerate,
    ) -> Self {
        Self {
            a,
            b,
            kind,
            connectivity,
            enumerate,
            a_to_b: vec![None; a.node_count()],
            b_used: vec![false; b.node_count()],
            decided: vec![false; a.node_count()],
            anchored: vec![false; a.node_count()],
            pairs: Vec::new(),
            edges: 0,
            best_score: 0,
            best: Vec::new(),
        }
    }

    fn score(&self) -> usize {
        match self.kind {
            McsAlgorithmKind::Induced => self.pairs.len(),
            McsAlgorithmKind::Edge => self.edges,
        }
    }

    fn bound(&self) -> usize {
        match self.kind {
            McsAlgorithmKind::Induced => {
                self.pairs.len() + self.decided.iter().filter(|&&d| !d).count()
            }
            McsAlgorithmKind::Edge => self.edges + self.remaining_possible_edges(),
        }
    }

    // Valid upper bound: only a-edges with an undecided endpoint can still be shared.
    fn remaining_possible_edges(&self) -> usize {
        self.a
            .edge_ids()
            .filter(|&e| {
                let [x, y] = self.a.edge_endpoints(e);
                !self.decided[x.index()] || !self.decided[y.index()]
            })
            .count()
    }

    // Shared edges (u,v) would add to the current mapping, or None if it violates
    // the induced edge-iff-edge constraint.
    fn pair_feasible<E: FnMut(EdgeId, EdgeId) -> bool>(
        &self,
        u: NodeId,
        v: NodeId,
        edge_match: &mut E,
    ) -> Option<usize> {
        let mut count = 0;
        for &(up, vp) in &self.pairs {
            let ea = self.a.find_edge(u, up);
            let eb = self.b.find_edge(v, vp);
            match self.kind {
                McsAlgorithmKind::Induced => match (ea, eb) {
                    (Some(ea), Some(eb)) => {
                        if !edge_match(ea, eb) {
                            return None;
                        }
                        count += 1;
                    }
                    (None, None) => {}
                    _ => return None,
                },
                McsAlgorithmKind::Edge => {
                    if let (Some(ea), Some(eb)) = (ea, eb) {
                        if edge_match(ea, eb) {
                            count += 1;
                        }
                    }
                }
            }
        }
        Some(count)
    }

    fn place_anchor<N, E>(
        &mut self,
        anchor: &[(NodeId, NodeId)],
        node_match: &mut N,
        edge_match: &mut E,
    ) -> bool
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
    {
        for &(u, v) in anchor {
            if u.index() >= self.a.node_count() || v.index() >= self.b.node_count() {
                return false;
            }
            if self.decided[u.index()] || self.b_used[v.index()] {
                return false;
            }
            if !node_match(u, v) {
                return false;
            }
            let Some(shared) = self.pair_feasible(u, v, edge_match) else {
                return false;
            };
            self.a_to_b[u.index()] = Some(v);
            self.b_used[v.index()] = true;
            self.decided[u.index()] = true;
            self.anchored[u.index()] = true;
            self.pairs.push((u, v));
            self.edges += shared;
        }
        true
    }

    // Warm start: install the largest valid sub-correspondence of `hint` as the
    // initial incumbent so branch-and-bound prunes from a good lower bound.
    fn seed_hint<N, E>(&mut self, hint: &[(NodeId, NodeId)], node_match: &mut N, edge_match: &mut E)
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
    {
        if hint.is_empty() {
            return;
        }
        let mut used_a = self.decided.clone();
        let mut used_b = self.b_used.clone();
        let mut pairs = self.pairs.clone();
        let mut edges = self.edges;
        for &(u, v) in hint {
            if u.index() >= self.a.node_count() || v.index() >= self.b.node_count() {
                continue;
            }
            if used_a[u.index()] || used_b[v.index()] || !node_match(u, v) {
                continue;
            }
            let mut shared = 0;
            let mut ok = true;
            for &(up, vp) in &pairs {
                let ea = self.a.find_edge(u, up);
                let eb = self.b.find_edge(v, vp);
                match self.kind {
                    McsAlgorithmKind::Induced => match (ea, eb) {
                        (Some(ea), Some(eb)) => {
                            if edge_match(ea, eb) {
                                shared += 1;
                            } else {
                                ok = false;
                                break;
                            }
                        }
                        (None, None) => {}
                        _ => {
                            ok = false;
                            break;
                        }
                    },
                    McsAlgorithmKind::Edge => {
                        if let (Some(ea), Some(eb)) = (ea, eb) {
                            if edge_match(ea, eb) {
                                shared += 1;
                            }
                        }
                    }
                }
            }
            if !ok {
                continue;
            }
            if self.connectivity == McsConnectivity::Connected && !pairs.is_empty() && shared == 0 {
                continue;
            }
            used_a[u.index()] = true;
            used_b[v.index()] = true;
            pairs.push((u, v));
            edges += shared;
        }
        let score = match self.kind {
            McsAlgorithmKind::Induced => pairs.len(),
            McsAlgorithmKind::Edge => edges,
        };
        if score > self.best_score {
            let cs = self.subgraph_from(&pairs, edges, edge_match);
            self.best_score = score;
            self.best = vec![cs];
        }
    }

    fn select_next(&self) -> Option<usize> {
        match self.connectivity {
            McsConnectivity::Disconnected => (0..self.a.node_count()).find(|&i| !self.decided[i]),
            McsConnectivity::Connected => {
                if self.pairs.is_empty() {
                    (0..self.a.node_count()).find(|&i| !self.decided[i])
                } else {
                    (0..self.a.node_count()).find(|&i| !self.decided[i] && self.in_a_frontier(i))
                }
            }
        }
    }

    fn in_a_frontier(&self, i: usize) -> bool {
        self.a
            .neighbors(NodeId(i as u32))
            .iter()
            .any(|n| self.a_to_b[n.node.index()].is_some())
    }

    fn a_has_undecided_neighbor(&self, u: NodeId) -> bool {
        self.a
            .neighbors(u)
            .iter()
            .any(|n| !self.decided[n.node.index()])
    }

    fn b_has_unused_neighbor(&self, v: NodeId) -> bool {
        self.b
            .neighbors(v)
            .iter()
            .any(|n| !self.b_used[n.node.index()])
    }

    fn search<N, E>(&mut self, node_match: &mut N, edge_match: &mut E)
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
    {
        let threshold = match self.enumerate {
            Enumerate::Single => self.best_score,
            Enumerate::All => self.best_score.saturating_sub(1),
        };
        if self.bound() <= threshold {
            return;
        }

        let Some(u_idx) = self.select_next() else {
            self.record(edge_match);
            return;
        };
        let u = NodeId(u_idx as u32);

        for v_idx in 0..self.b.node_count() {
            if self.b_used[v_idx] {
                continue;
            }
            let v = NodeId(v_idx as u32);
            if !node_match(u, v) {
                continue;
            }
            let Some(shared) = self.pair_feasible(u, v, edge_match) else {
                continue;
            };
            if self.connectivity == McsConnectivity::Connected
                && !self.pairs.is_empty()
                && shared == 0
            {
                continue;
            }
            // An edge-MCES vertex must be incident to a shared edge: a fresh
            // (0-shared) map only helps if it can still gain one on both sides.
            if self.kind == McsAlgorithmKind::Edge
                && shared == 0
                && (!self.a_has_undecided_neighbor(u) || !self.b_has_unused_neighbor(v))
            {
                continue;
            }

            self.a_to_b[u_idx] = Some(v);
            self.b_used[v_idx] = true;
            self.decided[u_idx] = true;
            self.pairs.push((u, v));
            self.edges += shared;
            self.search(node_match, edge_match);
            self.edges -= shared;
            self.pairs.pop();
            self.decided[u_idx] = false;
            self.b_used[v_idx] = false;
            self.a_to_b[u_idx] = None;
        }

        if !self.anchored[u_idx] {
            self.decided[u_idx] = true;
            self.search(node_match, edge_match);
            self.decided[u_idx] = false;
        }
    }

    fn record<E: FnMut(EdgeId, EdgeId) -> bool>(&mut self, edge_match: &mut E) {
        let score = self.score();
        if score == 0 {
            return;
        }
        if score > self.best_score {
            let cs = self.subgraph_from(&self.pairs, self.edges, edge_match);
            self.best_score = score;
            self.best.clear();
            self.best.push(cs);
        } else if score == self.best_score && self.enumerate == Enumerate::All {
            let cs = self.subgraph_from(&self.pairs, self.edges, edge_match);
            self.best.push(cs);
        }
    }

    fn subgraph_from<E: FnMut(EdgeId, EdgeId) -> bool>(
        &self,
        pairs: &[(NodeId, NodeId)],
        edges: usize,
        edge_match: &mut E,
    ) -> CommonSubgraph {
        let mut mapping: Vec<(NodeId, NodeId)> = match self.kind {
            McsAlgorithmKind::Induced => pairs.to_vec(),
            McsAlgorithmKind::Edge => {
                let mut incident = vec![false; pairs.len()];
                for i in 0..pairs.len() {
                    for j in (i + 1)..pairs.len() {
                        let (ui, vi) = pairs[i];
                        let (uj, vj) = pairs[j];
                        if let (Some(ea), Some(eb)) =
                            (self.a.find_edge(ui, uj), self.b.find_edge(vi, vj))
                        {
                            if edge_match(ea, eb) {
                                incident[i] = true;
                                incident[j] = true;
                            }
                        }
                    }
                }
                pairs
                    .iter()
                    .enumerate()
                    .filter(|&(i, _)| incident[i])
                    .map(|(_, &p)| p)
                    .collect()
            }
        };
        mapping.sort_unstable_by_key(|&(u, _)| u.0);
        CommonSubgraph {
            mapping,
            edge_count: edges,
        }
    }
}

impl Graph {
    /// Every maximal common induced subgraph of `self` and `other`: each a
    /// non-extendable vertex correspondence under `node_match`, with an edge present
    /// iff present (and `edge_match`-compatible) in both graphs. Unlike the maximum,
    /// this enumerates *all* maximal overlaps, disconnected ones included.
    /// Enumerate **every** common induced subgraph (the empty one included) — every clique of the
    /// modular product. Its result contains [`maximal_common_subgraphs`](Self::maximal_common_subgraphs)';
    /// exponential in the worst case.
    pub fn enumerate_common_subgraphs<N, E>(
        &self,
        other: &Graph,
        node_match: &mut N,
        edge_match: &mut E,
        alg: CommonSubgraphEnumerationAlgorithm,
    ) -> Vec<CommonSubgraph>
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
    {
        let (pairs, neighbors) = self.modular_product(other, node_match, edge_match);
        let mut cliques = Vec::new();
        match alg {
            CommonSubgraphEnumerationAlgorithm::Backtracking => {
                all_cliques(&neighbors, &mut Vec::new(), bitvec![1; pairs.len()], &mut cliques);
            }
        }
        subgraphs_from_cliques(cliques, &pairs, self, other)
    }

    /// Enumerate the **maximal** (non-extendable) common induced subgraphs — the maximal cliques of
    /// the modular product.
    pub fn maximal_common_subgraphs<N, E>(
        &self,
        other: &Graph,
        node_match: &mut N,
        edge_match: &mut E,
        alg: MaximalCommonSubgraphAlgorithm,
    ) -> Vec<CommonSubgraph>
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
    {
        let (pairs, neighbors) = self.modular_product(other, node_match, edge_match);
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

    /// The modular product of `self` and `other` under the match predicates — its cliques are the
    /// common induced subgraphs. Vertices are the `node_match`-compatible node pairs; two distinct
    /// pairs are adjacent iff their edges agree (both present and `edge_match`-compatible, or both
    /// absent). Pairs sharing a node stay non-adjacent, enforcing injectivity. Returns the pair
    /// list and its per-vertex adjacency bitsets.
    fn modular_product<N, E>(
        &self,
        other: &Graph,
        node_match: &mut N,
        edge_match: &mut E,
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
                    _ => false,
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

/// Convert modular-product cliques (index lists) into `CommonSubgraph`s: recover each pair mapping
/// (sorted by the first graph's node) and count the edges the correspondence shares. Sorted by
/// mapping and deduplicated for a stable, canonical result.
fn subgraphs_from_cliques(
    cliques: Vec<Vec<usize>>,
    pairs: &[(NodeId, NodeId)],
    a: &Graph,
    b: &Graph,
) -> Vec<CommonSubgraph> {
    let mut subgraphs: Vec<CommonSubgraph> = cliques
        .into_iter()
        .map(|clique| {
            let mut edge_count = 0;
            for x in 0..clique.len() {
                for y in (x + 1)..clique.len() {
                    let (a1, b1) = pairs[clique[x]];
                    let (a2, b2) = pairs[clique[y]];
                    if a.find_edge(a1, a2).is_some() && b.find_edge(b1, b2).is_some() {
                        edge_count += 1;
                    }
                }
            }
            let mut mapping: Vec<(NodeId, NodeId)> = clique.iter().map(|&i| pairs[i]).collect();
            mapping.sort_unstable_by_key(|&(a, _)| a.0);
            CommonSubgraph {
                mapping,
                edge_count,
            }
        })
        .collect();
    subgraphs.sort_by(|x, y| x.mapping.cmp(&y.mapping));
    subgraphs.dedup();
    subgraphs
}

/// Enumerate **every** clique of the graph given by `neighbors` (per-vertex adjacency bitsets),
/// each exactly once, appending it to `out`; the empty clique is reported at the root. The current
/// clique is extended only with candidates of higher index than its last member, so each clique is
/// generated once in ascending-vertex order; `candidates` is the common-neighbor set of the
/// current clique.
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

/// Bron–Kerbosch with pivoting (Tomita et al.): append to `out` every maximal clique
/// of the graph given by `neighbors` (per-vertex adjacency bitsets).
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
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn any_node(_: NodeId, _: NodeId) -> bool {
        true
    }
    fn any_edge(_: EdgeId, _: EdgeId) -> bool {
        true
    }
    fn cross(q: NodeId, t: NodeId) -> bool {
        q.0 != t.0
    }
    fn reject_edge(_: EdgeId, _: EdgeId) -> bool {
        false
    }

    #[rstest]
    #[case::triangle(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), 3, 3)]
    #[case::path3_triangle(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), 2, 1)]
    #[case::path4_cycle4(Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]), 3, 2)]
    #[case::edge_edge(Graph::new(2, &[[0, 1]]), Graph::new(2, &[[0, 1]]), 2, 1)]
    #[case::isolated(Graph::new(2, &[]), Graph::new(1, &[]), 1, 0)]
    fn test_graph_maximum_common_induced_subgraph(
        #[case] a: Graph,
        #[case] b: Graph,
        #[case] nodes: usize,
        #[case] edges: usize,
    ) {
        let r = a.maximum_common_induced_subgraph(
            &b,
            &mut any_node,
            &mut any_edge,
            McsConnectivity::Disconnected,
            McisAlgorithm::McGregor,
        );
        assert_eq!(r.node_count(), nodes);
        assert_eq!(r.edge_count(), edges);
    }

    #[rstest]
    #[case::path4_cycle4(Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]), McsConnectivity::Disconnected, 4, 3)]
    #[case::triangle(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), McsConnectivity::Disconnected, 3, 3)]
    #[case::two_edges_disconnected(Graph::new(4, &[[0, 1], [2, 3]]), Graph::new(4, &[[0, 1], [2, 3]]), McsConnectivity::Disconnected, 4, 2)]
    #[case::two_edges_connected(Graph::new(4, &[[0, 1], [2, 3]]), Graph::new(4, &[[0, 1], [2, 3]]), McsConnectivity::Connected, 2, 1)]
    fn test_graph_maximum_common_edge_subgraph(
        #[case] a: Graph,
        #[case] b: Graph,
        #[case] connectivity: McsConnectivity,
        #[case] nodes: usize,
        #[case] edges: usize,
    ) {
        let r = a.maximum_common_edge_subgraph(
            &b,
            &mut any_node,
            &mut any_edge,
            connectivity,
            McesAlgorithm::McGregor,
        );
        assert_eq!(r.node_count(), nodes);
        assert_eq!(r.edge_count(), edges);
    }

    #[rstest]
    fn test_graph_maximum_common_edge_subgraphs() {
        let p3 = Graph::new(3, &[[0, 1], [1, 2]]);
        let all = p3.maximum_common_edge_subgraphs(
            &p3,
            &mut any_node,
            &mut any_edge,
            McsConnectivity::Connected,
            McesAlgorithm::McGregor,
        );
        assert_eq!(
            all,
            vec![
                CommonSubgraph {
                    mapping: vec![
                        (NodeId(0), NodeId(0)),
                        (NodeId(1), NodeId(1)),
                        (NodeId(2), NodeId(2))
                    ],
                    edge_count: 2,
                },
                CommonSubgraph {
                    mapping: vec![
                        (NodeId(0), NodeId(2)),
                        (NodeId(1), NodeId(1)),
                        (NodeId(2), NodeId(0))
                    ],
                    edge_count: 2,
                },
            ]
        );
    }

    #[rstest]
    fn test_graph_maximum_common_induced_subgraph_node_filter() {
        let e = Graph::new(2, &[[0, 1]]);
        let r = e.maximum_common_induced_subgraph(
            &e,
            &mut cross,
            &mut any_edge,
            McsConnectivity::Disconnected,
            McisAlgorithm::McGregor,
        );
        assert_eq!(
            r,
            CommonSubgraph {
                mapping: vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
                edge_count: 1,
            }
        );
    }

    #[rstest]
    fn test_graph_maximum_common_induced_subgraph_edge_filter() {
        let tri = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let r = tri.maximum_common_induced_subgraph(
            &tri,
            &mut any_node,
            &mut reject_edge,
            McsConnectivity::Disconnected,
            McisAlgorithm::McGregor,
        );
        assert_eq!(r.node_count(), 1);
        assert_eq!(r.edge_count(), 0);
    }

    #[rstest]
    fn test_graph_maximum_common_edge_subgraphs_seeded() {
        let tri = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let anchor = [(NodeId(0), NodeId(1))];
        let r = tri.maximum_common_edge_subgraphs_seeded(
            &tri,
            &mut any_node,
            &mut any_edge,
            &anchor,
            &[],
            McsConnectivity::Disconnected,
            McesAlgorithm::McGregor,
        );
        assert_eq!(r.len(), 2);
        for cs in &r {
            assert!(cs.mapping().contains(&(NodeId(0), NodeId(1))));
            assert_eq!(cs.node_count(), 3);
            assert_eq!(cs.edge_count(), 3);
        }
    }

    #[rstest]
    fn test_graph_maximum_common_edge_subgraphs_seeded_hint() {
        let p4 = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        let unseeded = p4.maximum_common_edge_subgraphs(
            &p4,
            &mut any_node,
            &mut any_edge,
            McsConnectivity::Connected,
            McesAlgorithm::McGregor,
        );
        let hint = [(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))];
        let seeded = p4.maximum_common_edge_subgraphs_seeded(
            &p4,
            &mut any_node,
            &mut any_edge,
            &[],
            &hint,
            McsConnectivity::Connected,
            McesAlgorithm::McGregor,
        );
        assert_eq!(seeded, unseeded);
        assert_eq!(unseeded.len(), 2);
        assert!(unseeded
            .iter()
            .all(|cs| cs.node_count() == 4 && cs.edge_count() == 3));
    }

    #[rstest]
    fn test_graph_maximum_common_edge_subgraphs_seeded_empty() {
        // empty anchor + empty hint is exactly the unseeded edge MCS.
        let p4 = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        let plain = p4.maximum_common_edge_subgraphs(
            &p4,
            &mut any_node,
            &mut any_edge,
            McsConnectivity::Connected,
            McesAlgorithm::McGregor,
        );
        let seeded = p4.maximum_common_edge_subgraphs_seeded(
            &p4,
            &mut any_node,
            &mut any_edge,
            &[],
            &[],
            McsConnectivity::Connected,
            McesAlgorithm::McGregor,
        );
        assert_eq!(seeded, plain);
    }

    #[rstest]
    fn test_graph_maximum_common_induced_subgraph_empty() {
        let tri = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let empty = Graph::default();
        let r = tri.maximum_common_induced_subgraph(
            &empty,
            &mut any_node,
            &mut any_edge,
            McsConnectivity::Disconnected,
            McisAlgorithm::McGregor,
        );
        assert_eq!(
            r,
            CommonSubgraph {
                mapping: vec![],
                edge_count: 0,
            }
        );
        assert!(r.is_empty());
    }

    #[rstest]
    #[case::single_edge(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        vec![
            CommonSubgraph {
                mapping: vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))],
                edge_count: 1,
            },
            CommonSubgraph {
                mapping: vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
                edge_count: 1,
            },
        ]
    )]
    #[case::single_node(
        Graph::new(1, &[]),
        Graph::new(1, &[]),
        vec![CommonSubgraph {
            mapping: vec![(NodeId(0), NodeId(0))],
            edge_count: 0,
        }]
    )]
    fn test_graph_maximal_common_subgraphs(
        #[case] a: Graph,
        #[case] b: Graph,
        #[case] expected: Vec<CommonSubgraph>,
    ) {
        assert_eq!(
            a.maximal_common_subgraphs(
                &b,
                &mut any_node,
                &mut any_edge,
                MaximalCommonSubgraphAlgorithm::BronKerbosch,
            ),
            expected,
        );
    }

    #[rstest]
    #[case::path_with_d_edge(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![0u8, 1, 2],
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![0u8, 1, 2],
        vec![CommonSubgraph {
            mapping: vec![
                (NodeId(0), NodeId(0)),
                (NodeId(1), NodeId(1)),
                (NodeId(2), NodeId(2)),
            ],
            edge_count: 2,
        }]
    )]
    #[case::disjoint_edges(
        Graph::new(4, &[[0, 1], [2, 3]]),
        vec![0u8, 1, 2, 3],
        Graph::new(4, &[[0, 1], [2, 3]]),
        vec![0u8, 1, 2, 3],
        vec![CommonSubgraph {
            mapping: vec![
                (NodeId(0), NodeId(0)),
                (NodeId(1), NodeId(1)),
                (NodeId(2), NodeId(2)),
                (NodeId(3), NodeId(3)),
            ],
            edge_count: 2,
        }]
    )]
    fn test_graph_maximal_common_subgraphs_labeled(
        #[case] a: Graph,
        #[case] labels_a: Vec<u8>,
        #[case] b: Graph,
        #[case] labels_b: Vec<u8>,
        #[case] expected: Vec<CommonSubgraph>,
    ) {
        let mut node_match = |x: NodeId, y: NodeId| labels_a[x.index()] == labels_b[y.index()];
        assert_eq!(
            a.maximal_common_subgraphs(
                &b,
                &mut node_match,
                &mut any_edge,
                MaximalCommonSubgraphAlgorithm::BronKerbosch,
            ),
            expected,
        );
    }

    #[rstest]
    fn test_graph_maximal_common_subgraphs_edge_filter() {
        let edge = Graph::new(2, &[[0, 1]]);
        assert_eq!(
            edge.maximal_common_subgraphs(
                &edge,
                &mut any_node,
                &mut reject_edge,
                MaximalCommonSubgraphAlgorithm::BronKerbosch,
            ),
            vec![
                CommonSubgraph {
                    mapping: vec![(NodeId(0), NodeId(0))],
                    edge_count: 0,
                },
                CommonSubgraph {
                    mapping: vec![(NodeId(0), NodeId(1))],
                    edge_count: 0,
                },
                CommonSubgraph {
                    mapping: vec![(NodeId(1), NodeId(0))],
                    edge_count: 0,
                },
                CommonSubgraph {
                    mapping: vec![(NodeId(1), NodeId(1))],
                    edge_count: 0,
                },
            ],
        );
    }

    #[rstest]
    fn test_graph_maximal_common_subgraphs_empty() {
        let a = Graph::new(1, &[]);
        let b = Graph::new(1, &[]);
        assert_eq!(
            a.maximal_common_subgraphs(
                &b,
                &mut cross,
                &mut any_edge,
                MaximalCommonSubgraphAlgorithm::BronKerbosch,
            ),
            vec![CommonSubgraph {
                mapping: vec![],
                edge_count: 0,
            }],
        );
    }

    // Complete enumeration returns *every* clique of the modular product, empty included. For two
    // single-node graphs (`any_node`): the empty subgraph and the one singleton. For an edge
    // matched against itself: empty, four singletons, and the two 2-mappings that preserve the edge.
    #[rstest]
    #[case::single_nodes(
        Graph::new(1, &[]),
        Graph::new(1, &[]),
        vec![
            CommonSubgraph { mapping: vec![], edge_count: 0 },
            CommonSubgraph { mapping: vec![(NodeId(0), NodeId(0))], edge_count: 0 },
        ]
    )]
    #[case::edge(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        vec![
            CommonSubgraph { mapping: vec![], edge_count: 0 },
            CommonSubgraph { mapping: vec![(NodeId(0), NodeId(0))], edge_count: 0 },
            CommonSubgraph {
                mapping: vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))],
                edge_count: 1,
            },
            CommonSubgraph { mapping: vec![(NodeId(0), NodeId(1))], edge_count: 0 },
            CommonSubgraph {
                mapping: vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
                edge_count: 1,
            },
            CommonSubgraph { mapping: vec![(NodeId(1), NodeId(0))], edge_count: 0 },
            CommonSubgraph { mapping: vec![(NodeId(1), NodeId(1))], edge_count: 0 },
        ]
    )]
    fn test_graph_enumerate_common_subgraphs(
        #[case] a: Graph,
        #[case] b: Graph,
        #[case] expected: Vec<CommonSubgraph>,
    ) {
        assert_eq!(
            a.enumerate_common_subgraphs(
                &b,
                &mut any_node,
                &mut any_edge,
                CommonSubgraphEnumerationAlgorithm::Backtracking,
            ),
            expected,
        );
    }

    // The complete enumeration is a superset of the maximal one, and every maximal subgraph appears
    // in it (agreement between the two tasks on a shared modular product).
    #[rstest]
    #[case::triangles(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), Graph::new(3, &[[0, 1], [1, 2], [0, 2]]))]
    #[case::path_cycle(Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]))]
    fn test_graph_enumerate_common_subgraphs_contains_maximal(#[case] a: Graph, #[case] b: Graph) {
        let all = a.enumerate_common_subgraphs(
            &b,
            &mut any_node,
            &mut any_edge,
            CommonSubgraphEnumerationAlgorithm::Backtracking,
        );
        let maximal = a.maximal_common_subgraphs(
            &b,
            &mut any_node,
            &mut any_edge,
            MaximalCommonSubgraphAlgorithm::BronKerbosch,
        );
        assert!(all.len() > maximal.len());
        assert!(all.contains(&CommonSubgraph { mapping: vec![], edge_count: 0 }));
        for m in &maximal {
            assert!(all.contains(m), "maximal subgraph {m:?} missing from complete enumeration");
        }
    }
}
