//! Maximum matching and matching enumeration.

use std::collections::VecDeque;

use crate::graph::{EdgeId, Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaxMatchingAlgorithm {
    Edmonds,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatchingEnumerationAlgorithm {
    BranchAndBound,
}

/// A matching: a set of edges with no shared endpoints.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matching {
    edges: Vec<EdgeId>,
    mate: Vec<Option<NodeId>>,
}

impl Matching {
    fn from_mate_array(graph: &Graph, mate: &[i32]) -> Self {
        let mut edges = Vec::new();
        for eid in graph.edge_ids() {
            let [a, b] = graph.edge_endpoints(eid);
            if mate[a.index()] == b.0 as i32 {
                edges.push(eid);
            }
        }
        edges.sort_unstable();
        let mate_opt = mate
            .iter()
            .map(|&m| {
                if m >= 0 {
                    Some(NodeId(m as u32))
                } else {
                    None
                }
            })
            .collect();
        Self {
            edges,
            mate: mate_opt,
        }
    }

    pub fn edges(&self) -> &[EdgeId] {
        &self.edges
    }

    pub fn size(&self) -> usize {
        self.edges.len()
    }

    pub fn is_perfect(&self, node_count: usize) -> bool {
        self.edges.len() * 2 == node_count
    }

    pub fn mate(&self, node: NodeId) -> Option<NodeId> {
        self.mate.get(node.index()).copied().flatten()
    }

    pub fn is_matched(&self, node: NodeId) -> bool {
        self.mate(node).is_some()
    }
}

impl Graph {
    pub fn maximum_matching(&self, alg: MaxMatchingAlgorithm) -> Matching {
        match alg {
            MaxMatchingAlgorithm::Edmonds => self.maximum_matching_edmonds(),
        }
    }

    pub fn enumerate_perfect_matchings(
        &self,
        alg: MatchingEnumerationAlgorithm,
    ) -> Vec<Matching> {
        match alg {
            MatchingEnumerationAlgorithm::BranchAndBound => {
                self.enumerate_perfect_matchings_branch_and_bound()
            }
        }
    }

    pub fn enumerate_maximum_matchings(
        &self,
        alg: MatchingEnumerationAlgorithm,
    ) -> Vec<Matching> {
        match alg {
            MatchingEnumerationAlgorithm::BranchAndBound => {
                self.enumerate_maximum_matchings_branch_and_bound()
            }
        }
    }

    // Edmonds 1965, Gabow simplification 1976. Ref impl: cp-algorithms.com.
    fn maximum_matching_edmonds(&self) -> Matching {
        let n = self.node_count();
        if n == 0 {
            return Matching {
                edges: Vec::new(),
                mate: Vec::new(),
            };
        }

        let mut mate = vec![-1i32; n];

        for eid in self.edge_ids() {
            let [a, b] = self.edge_endpoints(eid);
            if mate[a.index()] < 0 && mate[b.index()] < 0 {
                mate[a.index()] = b.0 as i32;
                mate[b.index()] = a.0 as i32;
            }
        }

        for v in 0..n {
            if mate[v] < 0 {
                augment_from(self, &mut mate, v);
            }
        }

        Matching::from_mate_array(self, &mate)
    }

    // Branch-and-bound with Edmonds oracle.
    fn enumerate_perfect_matchings_branch_and_bound(&self) -> Vec<Matching> {
        let initial = self.maximum_matching_edmonds();
        if !initial.is_perfect(self.node_count()) {
            return Vec::new();
        }
        if self.node_count() == 0 {
            return vec![initial];
        }
        let mut result = Vec::new();
        let edges: Vec<EdgeId> = self.edge_ids().collect();
        let mut included = vec![false; self.edge_bound()];
        let mut excluded = vec![false; self.edge_bound()];
        enumerate_rec(
            self,
            &edges,
            &mut included,
            &mut excluded,
            self.node_count() / 2,
            &mut result,
        );
        result
    }

    // Branch-and-bound with Edmonds oracle.
    fn enumerate_maximum_matchings_branch_and_bound(&self) -> Vec<Matching> {
        let initial = self.maximum_matching_edmonds();
        let target_size = initial.size();
        if target_size == 0 {
            return vec![initial];
        }
        let mut result = Vec::new();
        let edges: Vec<EdgeId> = self.edge_ids().collect();
        let mut included = vec![false; self.edge_bound()];
        let mut excluded = vec![false; self.edge_bound()];
        enumerate_rec(
            self,
            &edges,
            &mut included,
            &mut excluded,
            target_size,
            &mut result,
        );
        result
    }
}

fn augment_from(graph: &Graph, mate: &mut [i32], root: usize) -> bool {
    let n = graph.node_bound();
    let mut used = vec![false; n];
    let mut p = vec![-1i32; n];
    let mut base: Vec<usize> = (0..n).collect();

    used[root] = true;
    let mut queue = VecDeque::new();
    queue.push_back(root);

    while let Some(v) = queue.pop_front() {
        for nbr in graph.neighbors(NodeId(v as u32)) {
            let to = nbr.node.index();

            if base[v] == base[to] || mate[v] == to as i32 {
                continue;
            }

            if to == root || (mate[to] >= 0 && p[mate[to] as usize] >= 0) {
                let curbase = blossom_lca(&base, mate, &p, root, v, to);
                let mut in_blossom = vec![false; n];
                mark_blossom_path(&mut in_blossom, &base, mate, &mut p, v, curbase, to);
                mark_blossom_path(&mut in_blossom, &base, mate, &mut p, to, curbase, v);
                for i in 0..n {
                    if in_blossom[base[i]] {
                        base[i] = curbase;
                        if !used[i] {
                            used[i] = true;
                            queue.push_back(i);
                        }
                    }
                }
            } else if p[to] < 0 {
                p[to] = v as i32;
                if mate[to] < 0 {
                    let mut u = to as i32;
                    while u >= 0 {
                        let pv = p[u as usize];
                        let ppv = mate[pv as usize];
                        mate[u as usize] = pv;
                        mate[pv as usize] = u;
                        u = ppv;
                    }
                    return true;
                }
                let matched = mate[to] as usize;
                used[matched] = true;
                queue.push_back(matched);
            }
        }
    }

    false
}

fn blossom_lca(base: &[usize], mate: &[i32], p: &[i32], root: usize, a: usize, b: usize) -> usize {
    let n = base.len();
    let mut visited = vec![false; n];
    let mut a = base[a];
    let mut b = base[b];
    loop {
        visited[a] = true;
        if a == root || mate[a] < 0 {
            break;
        }
        a = base[p[mate[a] as usize] as usize];
    }
    loop {
        if visited[b] {
            return b;
        }
        b = base[p[mate[b] as usize] as usize];
    }
}

fn mark_blossom_path(
    in_blossom: &mut [bool],
    base: &[usize],
    mate: &[i32],
    p: &mut [i32],
    mut v: usize,
    b: usize,
    mut child: usize,
) {
    while base[v] != b {
        in_blossom[base[v]] = true;
        in_blossom[base[mate[v] as usize]] = true;
        p[v] = child as i32;
        child = mate[v] as usize;
        v = p[mate[v] as usize] as usize;
    }
}

// ── Branch-and-bound enumeration ────────────────────────────────────

fn enumerate_rec(
    graph: &Graph,
    edges: &[EdgeId],
    included: &mut [bool],
    excluded: &mut [bool],
    target_size: usize,
    result: &mut Vec<Matching>,
) {
    let mut node_matched = vec![false; graph.node_bound()];
    let mut current_size = 0usize;
    for (i, &inc) in included.iter().enumerate() {
        if inc {
            let [a, b] = graph.edge_endpoints(EdgeId(i as u32));
            node_matched[a.index()] = true;
            node_matched[b.index()] = true;
            current_size += 1;
        }
    }

    if current_size == target_size {
        let mate = build_mate(graph, included);
        result.push(Matching::from_mate_array(graph, &mate));
        return;
    }

    let branch_edge = edges.iter().find(|&&eid| {
        !included[eid.index()] && !excluded[eid.index()] && {
            let [a, b] = graph.edge_endpoints(eid);
            !node_matched[a.index()] && !node_matched[b.index()]
        }
    });

    let Some(&eid) = branch_edge else {
        return;
    };

    let [a, b] = graph.edge_endpoints(eid);

    // Include branch
    included[eid.index()] = true;
    let mut newly_excluded = Vec::new();
    for nbr in graph.neighbors(a).iter().chain(graph.neighbors(b)) {
        if nbr.edge != eid && !excluded[nbr.edge.index()] && !included[nbr.edge.index()] {
            excluded[nbr.edge.index()] = true;
            newly_excluded.push(nbr.edge);
        }
    }
    if can_reach(graph, included, excluded, current_size + 1, target_size) {
        enumerate_rec(graph, edges, included, excluded, target_size, result);
    }
    included[eid.index()] = false;
    for &e in &newly_excluded {
        excluded[e.index()] = false;
    }

    // Exclude branch
    excluded[eid.index()] = true;
    if can_reach(graph, included, excluded, current_size, target_size) {
        enumerate_rec(graph, edges, included, excluded, target_size, result);
    }
    excluded[eid.index()] = false;
}

fn can_reach(
    graph: &Graph,
    included: &[bool],
    excluded: &[bool],
    current_size: usize,
    target_size: usize,
) -> bool {
    if current_size > target_size {
        return false;
    }
    let remaining = target_size - current_size;
    if remaining == 0 {
        return true;
    }

    let mut node_matched = vec![false; graph.node_bound()];
    for (i, &inc) in included.iter().enumerate() {
        if inc {
            let [a, b] = graph.edge_endpoints(EdgeId(i as u32));
            node_matched[a.index()] = true;
            node_matched[b.index()] = true;
        }
    }

    let unmatched: usize = node_matched.iter().filter(|&&m| !m).count();
    if unmatched < remaining * 2 {
        return false;
    }

    // Greedy matching on residual graph as tighter bound
    let mut residual_matched = node_matched;
    let mut greedy_count = 0usize;
    for eid in graph.edge_ids() {
        if included[eid.index()] || excluded[eid.index()] {
            continue;
        }
        let [a, b] = graph.edge_endpoints(eid);
        if !residual_matched[a.index()] && !residual_matched[b.index()] {
            residual_matched[a.index()] = true;
            residual_matched[b.index()] = true;
            greedy_count += 1;
        }
    }

    greedy_count >= remaining
}

fn build_mate(graph: &Graph, included: &[bool]) -> Vec<i32> {
    let mut mate = vec![-1i32; graph.node_bound()];
    for (i, &inc) in included.iter().enumerate() {
        if inc {
            let [a, b] = graph.edge_endpoints(EdgeId(i as u32));
            mate[a.index()] = b.0 as i32;
            mate[b.index()] = a.0 as i32;
        }
    }
    mate
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use crate::graph::Graph;

    use super::{
        MatchingEnumerationAlgorithm::BranchAndBound, MaxMatchingAlgorithm::Edmonds, Matching,
    };

    #[test]
    fn test_matching_empty() {
        let g = Graph::default();
        let m = g.maximum_matching(Edmonds);
        assert_eq!(m.size(), 0);
        assert!(m.is_perfect(0));
        assert!(m.edges().is_empty());
    }

    #[rstest]
    #[case::single_edge(2, vec![[0, 1]], 1, true)]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], 1, false)]
    #[case::square(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]], 2, true)]
    #[case::path_4(4, vec![[0, 1], [1, 2], [2, 3]], 2, true)]
    #[case::path_5(5, vec![[0, 1], [1, 2], [2, 3], [3, 4]], 2, false)]
    #[case::k4(4, vec![[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], 2, true)]
    fn test_graph_maximum_matching(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected_size: usize,
        #[case] expected_perfect: bool,
    ) {
        let g = Graph::new(node_count, &edges);
        let m = g.maximum_matching(Edmonds);
        assert_eq!(m.size(), expected_size, "matching size");
        assert_eq!(m.is_perfect(node_count), expected_perfect, "is_perfect");
        assert_matching_valid(&g, &m);
    }

    #[test]
    fn test_graph_maximum_matching_petersen() {
        let g = Graph::new(
            10,
            &[
                [0, 1], [1, 2], [2, 3], [3, 4], [4, 0],
                [5, 7], [7, 9], [9, 6], [6, 8], [8, 5],
                [0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
            ],
        );
        let m = g.maximum_matching(Edmonds);
        assert_eq!(m.size(), 5);
        assert!(m.is_perfect(10));
        assert_matching_valid(&g, &m);
    }

    #[rstest]
    #[case::single_edge(2, vec![[0, 1]], 1)]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], 0)]
    #[case::square(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]], 2)]
    #[case::c6(6, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]], 2)]
    #[case::k4(4, vec![[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], 3)]
    fn test_graph_enumerate_perfect_matchings(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected_count: usize,
    ) {
        let g = Graph::new(node_count, &edges);
        let matchings = g.enumerate_perfect_matchings(BranchAndBound);
        assert_eq!(
            matchings.len(),
            expected_count,
            "number of perfect matchings"
        );
        for m in &matchings {
            assert!(m.is_perfect(node_count));
            assert_matching_valid(&g, m);
        }
        assert_all_distinct(&matchings);
    }

    #[rstest]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], 3)]
    #[case::path_3(3, vec![[0, 1], [1, 2]], 2)]
    #[case::single_edge(2, vec![[0, 1]], 1)]
    fn test_graph_enumerate_maximum_matchings(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected_count: usize,
    ) {
        let g = Graph::new(node_count, &edges);
        let initial = g.maximum_matching(Edmonds);
        let target = initial.size();
        let matchings = g.enumerate_maximum_matchings(BranchAndBound);
        assert_eq!(
            matchings.len(),
            expected_count,
            "number of maximum matchings"
        );
        for m in &matchings {
            assert_eq!(m.size(), target);
            assert_matching_valid(&g, m);
        }
        assert_all_distinct(&matchings);
    }

    // ── helpers ─────────────────────────────────────────────────────

    fn assert_matching_valid(graph: &Graph, matching: &Matching) {
        let mut matched = vec![false; graph.node_count()];
        for &eid in matching.edges() {
            let [a, b] = graph.edge_endpoints(eid);
            assert!(!matched[a.index()], "node {} matched twice", a.0);
            assert!(!matched[b.index()], "node {} matched twice", b.0);
            matched[a.index()] = true;
            matched[b.index()] = true;
        }
    }

    fn assert_all_distinct(matchings: &[Matching]) {
        for i in 0..matchings.len() {
            for j in (i + 1)..matchings.len() {
                assert_ne!(
                    matchings[i].edges(),
                    matchings[j].edges(),
                    "matchings {} and {} are identical",
                    i,
                    j
                );
            }
        }
    }
}
