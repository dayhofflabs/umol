//! Maximum common subgraphs.
//!
//! Current operations use McGregor backtracking for maximum common induced
//! subgraphs and maximum common edge subgraphs, including seeded edge search.
//! See [McGregor, *Backtrack Search Algorithms and the Maximal Common
//! Subgraph Problem* (1982)](https://doi.org/10.1002/spe.4380120103).

#[cfg(test)]
use super::enumeration::{
    CommonSubgraphEnumerationAlgorithm, EmbeddingKind, MaximalCommonSubgraphAlgorithm,
};
use crate::correspondence::{Correspondence, GraphCorrespondence};
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum McsAlgorithmKind {
    Induced,
    Edge,
}

impl Graph {
    /// Every maximum common induced subgraph (all of the largest size).
    pub fn maximum_common_induced_subgraphs(
        &self,
        other: &Graph,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        connectivity: McsConnectivity,
        alg: McisAlgorithm,
    ) -> Vec<GraphCorrespondence> {
        match alg {
            McisAlgorithm::McGregor => self.mcs_mcgregor(
                other,
                node_match,
                edge_match,
                &[],
                &[],
                connectivity,
                McsAlgorithmKind::Induced,
            ),
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
    ) -> Vec<GraphCorrespondence> {
        match alg {
            McesAlgorithm::McGregor => self.mcs_mcgregor(
                other,
                node_match,
                edge_match,
                &[],
                &[],
                connectivity,
                McsAlgorithmKind::Edge,
            ),
        }
    }

    /// Every seeded maximum common edge subgraph. `anchor` pairs are forced into the result and
    /// never skipped; `hint` pairs warm-start the incumbent bound (may be discarded). Either may be
    /// empty.
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
    ) -> Vec<GraphCorrespondence> {
        match alg {
            McesAlgorithm::McGregor => self.mcs_mcgregor(
                other,
                node_match,
                edge_match,
                anchor,
                hint,
                connectivity,
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
        kind: McsAlgorithmKind,
    ) -> Vec<GraphCorrespondence>
    where
        N: FnMut(NodeId, NodeId) -> bool,
        E: FnMut(EdgeId, EdgeId) -> bool,
    {
        let mut state = McsState::new(self, other, kind, connectivity);
        if !state.place_anchor(anchor, node_match, edge_match) {
            return Vec::new();
        }
        state.seed_hint(hint, node_match, edge_match);
        state.search(node_match, edge_match);
        if state.best.is_empty() {
            let fallback = state.subgraph_from(&state.pairs, edge_match);
            state.best.push(fallback);
        }
        let mut best = state.best;
        best.sort_by(|x, y| x.nodes().mates().cmp(y.nodes().mates()));
        best.dedup();
        best
    }
}

struct McsState<'g> {
    a: &'g Graph,
    b: &'g Graph,
    kind: McsAlgorithmKind,
    connectivity: McsConnectivity,
    a_to_b: Vec<Option<NodeId>>,
    b_used: Vec<bool>,
    decided: Vec<bool>,
    anchored: Vec<bool>,
    pairs: Vec<(NodeId, NodeId)>,
    edges: usize,
    best_score: usize,
    best: Vec<GraphCorrespondence>,
}

impl<'g> McsState<'g> {
    fn new(
        a: &'g Graph,
        b: &'g Graph,
        kind: McsAlgorithmKind,
        connectivity: McsConnectivity,
    ) -> Self {
        Self {
            a,
            b,
            kind,
            connectivity,
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
            let cs = self.subgraph_from(&pairs, edge_match);
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
        // Enumerate all maxima: keep ties (bound at best − 1), not just the first.
        if self.bound() <= self.best_score.saturating_sub(1) {
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
            let cs = self.subgraph_from(&self.pairs, edge_match);
            self.best_score = score;
            self.best.clear();
            self.best.push(cs);
        } else if score == self.best_score {
            let cs = self.subgraph_from(&self.pairs, edge_match);
            self.best.push(cs);
        }
    }

    fn subgraph_from<E: FnMut(EdgeId, EdgeId) -> bool>(
        &self,
        pairs: &[(NodeId, NodeId)],
        edge_match: &mut E,
    ) -> GraphCorrespondence {
        // The matched edge pairs over the correspondence (edge-predicate-filtered), and the nodes
        // they touch (an edge-subgraph mapping keeps only edge-incident nodes).
        let mut edges: Vec<(EdgeId, EdgeId)> = Vec::new();
        let mut incident = vec![false; pairs.len()];
        for i in 0..pairs.len() {
            for j in (i + 1)..pairs.len() {
                let (ui, vi) = pairs[i];
                let (uj, vj) = pairs[j];
                if let (Some(ea), Some(eb)) = (self.a.find_edge(ui, uj), self.b.find_edge(vi, vj)) {
                    if edge_match(ea, eb) {
                        incident[i] = true;
                        incident[j] = true;
                        edges.push((ea, eb));
                    }
                }
            }
        }
        let mapping: Vec<(NodeId, NodeId)> = match self.kind {
            McsAlgorithmKind::Induced => pairs.to_vec(),
            McsAlgorithmKind::Edge => pairs
                .iter()
                .enumerate()
                .filter(|&(i, _)| incident[i])
                .map(|(_, &p)| p)
                .collect(),
        };
        GraphCorrespondence::new(
            Correspondence::new(mapping, self.a.node_count(), self.b.node_count()),
            Correspondence::new(edges, self.a.edge_count(), self.b.edge_count()),
        )
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn summary(c: &GraphCorrespondence) -> (Vec<(NodeId, NodeId)>, usize) {
        (c.nodes().mates().to_vec(), c.edges().mate_count())
    }

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
        let r = a
            .maximum_common_induced_subgraphs(
                &b,
                &mut any_node,
                &mut any_edge,
                McsConnectivity::Disconnected,
                McisAlgorithm::McGregor,
            )
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(r.nodes().mate_count(), nodes);
        assert_eq!(r.edges().mate_count(), edges);
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
        let r = a
            .maximum_common_edge_subgraphs(
                &b,
                &mut any_node,
                &mut any_edge,
                connectivity,
                McesAlgorithm::McGregor,
            )
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(r.nodes().mate_count(), nodes);
        assert_eq!(r.edges().mate_count(), edges);
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
            all.iter().map(summary).collect::<Vec<_>>(),
            vec![
                (
                    vec![
                        (NodeId(0), NodeId(0)),
                        (NodeId(1), NodeId(1)),
                        (NodeId(2), NodeId(2))
                    ],
                    2
                ),
                (
                    vec![
                        (NodeId(0), NodeId(2)),
                        (NodeId(1), NodeId(1)),
                        (NodeId(2), NodeId(0))
                    ],
                    2
                ),
            ]
        );
    }

    #[rstest]
    fn test_graph_maximum_common_induced_subgraph_node_filter() {
        let e = Graph::new(2, &[[0, 1]]);
        let r = e
            .maximum_common_induced_subgraphs(
                &e,
                &mut cross,
                &mut any_edge,
                McsConnectivity::Disconnected,
                McisAlgorithm::McGregor,
            )
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(
            summary(&r),
            (vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))], 1)
        );
    }

    #[rstest]
    fn test_graph_maximum_common_induced_subgraph_edge_filter() {
        let tri = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let r = tri
            .maximum_common_induced_subgraphs(
                &tri,
                &mut any_node,
                &mut reject_edge,
                McsConnectivity::Disconnected,
                McisAlgorithm::McGregor,
            )
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(r.nodes().mate_count(), 1);
        assert_eq!(r.edges().mate_count(), 0);
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
            assert!(cs.nodes().mates().contains(&(NodeId(0), NodeId(1))));
            assert_eq!(cs.nodes().mate_count(), 3);
            assert_eq!(cs.edges().mate_count(), 3);
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
            .all(|cs| cs.nodes().mate_count() == 4 && cs.edges().mate_count() == 3));
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
        let r = tri
            .maximum_common_induced_subgraphs(
                &empty,
                &mut any_node,
                &mut any_edge,
                McsConnectivity::Disconnected,
                McisAlgorithm::McGregor,
            )
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(summary(&r), (vec![], 0));
        assert!(r.nodes().mates().is_empty());
    }

    #[rstest]
    #[case::single_edge(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        vec![
            (vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 1),
            (vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))], 1),
        ]
    )]
    #[case::single_node(
        Graph::new(1, &[]),
        Graph::new(1, &[]),
        vec![(vec![(NodeId(0), NodeId(0))], 0)]
    )]
    fn test_graph_maximal_common_subgraphs(
        #[case] a: Graph,
        #[case] b: Graph,
        #[case] expected: Vec<(Vec<(NodeId, NodeId)>, usize)>,
    ) {
        assert_eq!(
            a.maximal_common_subgraphs(
                &b,
                &mut any_node,
                &mut any_edge,
                EmbeddingKind::Induced,
                MaximalCommonSubgraphAlgorithm::BronKerbosch,
            )
            .iter()
            .map(summary)
            .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::path_with_d_edge(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![0u8, 1, 2],
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![0u8, 1, 2],
        vec![(
            vec![
                (NodeId(0), NodeId(0)),
                (NodeId(1), NodeId(1)),
                (NodeId(2), NodeId(2)),
            ],
            2
        )]
    )]
    #[case::disjoint_edges(
        Graph::new(4, &[[0, 1], [2, 3]]),
        vec![0u8, 1, 2, 3],
        Graph::new(4, &[[0, 1], [2, 3]]),
        vec![0u8, 1, 2, 3],
        vec![(
            vec![
                (NodeId(0), NodeId(0)),
                (NodeId(1), NodeId(1)),
                (NodeId(2), NodeId(2)),
                (NodeId(3), NodeId(3)),
            ],
            2
        )]
    )]
    fn test_graph_maximal_common_subgraphs_labeled(
        #[case] a: Graph,
        #[case] labels_a: Vec<u8>,
        #[case] b: Graph,
        #[case] labels_b: Vec<u8>,
        #[case] expected: Vec<(Vec<(NodeId, NodeId)>, usize)>,
    ) {
        let mut node_match = |x: NodeId, y: NodeId| labels_a[x.index()] == labels_b[y.index()];
        assert_eq!(
            a.maximal_common_subgraphs(
                &b,
                &mut node_match,
                &mut any_edge,
                EmbeddingKind::Induced,
                MaximalCommonSubgraphAlgorithm::BronKerbosch,
            )
            .iter()
            .map(summary)
            .collect::<Vec<_>>(),
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
                EmbeddingKind::Induced,
                MaximalCommonSubgraphAlgorithm::BronKerbosch,
            )
            .iter()
            .map(summary)
            .collect::<Vec<_>>(),
            vec![
                (vec![(NodeId(0), NodeId(0))], 0),
                (vec![(NodeId(0), NodeId(1))], 0),
                (vec![(NodeId(1), NodeId(0))], 0),
                (vec![(NodeId(1), NodeId(1))], 0),
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
                EmbeddingKind::Induced,
                MaximalCommonSubgraphAlgorithm::BronKerbosch,
            )
            .iter()
            .map(summary)
            .collect::<Vec<_>>(),
            vec![(vec![], 0)],
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
            (vec![], 0),
            (vec![(NodeId(0), NodeId(0))], 0),
        ]
    )]
    #[case::edge(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        vec![
            (vec![], 0),
            (vec![(NodeId(0), NodeId(0))], 0),
            (vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 1),
            (vec![(NodeId(0), NodeId(1))], 0),
            (vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))], 1),
            (vec![(NodeId(1), NodeId(0))], 0),
            (vec![(NodeId(1), NodeId(1))], 0),
        ]
    )]
    fn test_graph_enumerate_common_subgraphs(
        #[case] a: Graph,
        #[case] b: Graph,
        #[case] expected: Vec<(Vec<(NodeId, NodeId)>, usize)>,
    ) {
        assert_eq!(
            a.enumerate_common_subgraphs(
                &b,
                &mut any_node,
                &mut any_edge,
                EmbeddingKind::Induced,
                CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
            )
            .iter()
            .map(summary)
            .collect::<Vec<_>>(),
            expected,
        );
    }

    // The R1 (F–Cl) case: A's edge (1,2) is present, B's is absent. Under monomorphism the full
    // identity overlap survives (that edge is context in A); under induced it is dropped.
    #[rstest]
    fn test_graph_enumerate_common_subgraphs_monomorphism() {
        let a = Graph::new(3, &[[0, 1], [1, 2]]);
        let b = Graph::new(3, &[[0, 1]]);
        let full = vec![
            (NodeId(0), NodeId(0)),
            (NodeId(1), NodeId(1)),
            (NodeId(2), NodeId(2)),
        ];
        let monomorphism = a.enumerate_common_subgraphs(
            &b,
            &mut any_node,
            &mut any_edge,
            EmbeddingKind::Monomorphism,
            CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        );
        let induced = a.enumerate_common_subgraphs(
            &b,
            &mut any_node,
            &mut any_edge,
            EmbeddingKind::Induced,
            CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        );
        assert!(monomorphism
            .iter()
            .map(summary)
            .any(|(mapping, edges)| mapping == full && edges == 1));
        assert!(!induced
            .iter()
            .map(summary)
            .any(|(mapping, _)| mapping == full));
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
            EmbeddingKind::Induced,
            CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        );
        let maximal = a.maximal_common_subgraphs(
            &b,
            &mut any_node,
            &mut any_edge,
            EmbeddingKind::Induced,
            MaximalCommonSubgraphAlgorithm::BronKerbosch,
        );
        assert!(all.len() > maximal.len());
        assert!(all.iter().map(summary).any(|s| s == (vec![], 0)));
        for m in &maximal {
            assert!(
                all.contains(m),
                "maximal subgraph {m:?} missing from complete enumeration"
            );
        }
    }

    // A has the edge, B does not: the two "diagonal" mappings conflict on it. Induced drops them
    // (edge-iff-edge); Monomorphism keeps them adjacent (the edge is context in A).
    #[rstest]
    #[case::induced(EmbeddingKind::Induced, vec![])]
    #[case::monomorphism(EmbeddingKind::Monomorphism, vec![(0, 3), (1, 2)])]
    fn test_graph_modular_product(
        #[case] embedding: EmbeddingKind,
        #[case] expected_adjacency: Vec<(usize, usize)>,
    ) {
        let a = Graph::new(2, &[[0, 1]]);
        let b = Graph::new(2, &[]);
        let (pairs, neighbors) = a.modular_product(&b, &mut any_node, &mut any_edge, embedding);
        assert_eq!(
            pairs,
            vec![
                (NodeId(0), NodeId(0)),
                (NodeId(0), NodeId(1)),
                (NodeId(1), NodeId(0)),
                (NodeId(1), NodeId(1)),
            ],
        );
        let adjacency: Vec<(usize, usize)> = (0..pairs.len())
            .flat_map(|i| ((i + 1)..pairs.len()).map(move |j| (i, j)))
            .filter(|&(i, j)| neighbors[i][j])
            .collect();
        assert_eq!(adjacency, expected_adjacency);
    }
}
