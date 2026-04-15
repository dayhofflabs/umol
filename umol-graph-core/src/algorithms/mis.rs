//! Maximum independent set enumeration.

use crate::graph::{Graph, NodeId};

impl Graph {
    pub fn maximum_independent_set(&self) -> Vec<NodeId> {
        let bound = self.node_bound();
        if bound == 0 {
            return Vec::new();
        }

        let node_ids: Vec<NodeId> = self.node_ids().collect();
        let mut available = vec![true; bound];

        let mut current = Vec::new();
        let mut best = Vec::new();
        self.mis_branch(&node_ids, &mut available, &mut current, &mut best);
        best.sort_unstable();
        best
    }

    fn mis_branch(
        &self,
        node_ids: &[NodeId],
        available: &mut [bool],
        current: &mut Vec<NodeId>,
        best: &mut Vec<NodeId>,
    ) {
        let remaining = node_ids
            .iter()
            .filter(|id| available[id.index()])
            .count();
        if current.len() + remaining <= best.len() {
            return;
        }

        let Some(&node) = node_ids.iter().find(|id| available[id.index()]) else {
            let mut candidate = current.clone();
            candidate.sort_unstable();
            if candidate.len() > best.len()
                || (candidate.len() == best.len() && candidate < *best)
            {
                *best = candidate;
            }
            return;
        };

        // Include branch
        let mut changed = Vec::new();
        if available[node.index()] {
            changed.push(node);
            available[node.index()] = false;
        }
        for neighbor in self.neighbors(node) {
            if available[neighbor.node.index()] {
                changed.push(neighbor.node);
                available[neighbor.node.index()] = false;
            }
        }
        current.push(node);
        self.mis_branch(node_ids, available, current, best);
        current.pop();
        for idx in &changed {
            available[idx.index()] = true;
        }

        // Exclude branch
        available[node.index()] = false;
        self.mis_branch(node_ids, available, current, best);
        available[node.index()] = true;
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use crate::graph::{Graph, NodeId};

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    #[test]
    fn test_graph_maximum_independent_set_empty() {
        let g = Graph::default();
        assert!(g.maximum_independent_set().is_empty());
    }

    #[rstest]
    #[case::clique(
        3,
        vec![[0, 1], [1, 2], [0, 2]],
        vec![n(0)]
    )]
    #[case::path(
        4,
        vec![[0, 1], [1, 2], [2, 3]],
        vec![n(0), n(2)]
    )]
    #[case::cycle(
        4,
        vec![[0, 1], [1, 2], [2, 3], [3, 0]],
        vec![n(0), n(2)]
    )]
    #[case::disconnected(
        4,
        vec![[0, 1], [1, 2]],
        vec![n(0), n(2), n(3)]
    )]
    fn test_graph_maximum_independent_set(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected: Vec<NodeId>,
    ) {
        let g = Graph::new(node_count, &edges);
        assert_eq!(g.maximum_independent_set(), expected);
    }
}
