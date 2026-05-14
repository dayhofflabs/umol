//! Connected component labeling.

use crate::graph::{Graph, NodeId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConnectedComponentsAlgorithm {
    Bfs,
}

impl Graph {
    pub fn connected_components(&self, alg: ConnectedComponentsAlgorithm) -> Vec<Vec<NodeId>> {
        match alg {
            ConnectedComponentsAlgorithm::Bfs => self.connected_components_bfs(),
        }
    }

    /// Connected components restricted to `nodes`. Edges to nodes outside
    /// `nodes` are ignored. Component discovery order follows `nodes`; within
    /// each component, node ids are sorted.
    pub fn connected_components_in(
        &self,
        nodes: &[NodeId],
        alg: ConnectedComponentsAlgorithm,
    ) -> Vec<Vec<NodeId>> {
        match alg {
            ConnectedComponentsAlgorithm::Bfs => self.connected_components_in_bfs(nodes),
        }
    }

    // BFS flood fill. O(V+E).
    fn connected_components_bfs(&self) -> Vec<Vec<NodeId>> {
        let mut visited = vec![false; self.node_bound()];
        let mut components = Vec::new();

        for node in self.node_ids() {
            if visited[node.index()] {
                continue;
            }
            let mut component = Vec::new();
            let mut stack = vec![node];
            visited[node.index()] = true;

            while let Some(current) = stack.pop() {
                component.push(current);
                for neighbor in self.neighbors(current) {
                    if !visited[neighbor.node.index()] {
                        visited[neighbor.node.index()] = true;
                        stack.push(neighbor.node);
                    }
                }
            }

            component.sort_unstable();
            components.push(component);
        }

        components
    }

    fn connected_components_in_bfs(&self, nodes: &[NodeId]) -> Vec<Vec<NodeId>> {
        let mut in_subset = vec![false; self.node_bound()];
        for &node in nodes {
            in_subset[node.index()] = true;
        }
        let mut visited = vec![false; self.node_bound()];
        let mut components = Vec::new();

        for &node in nodes {
            if visited[node.index()] {
                continue;
            }
            let mut component = Vec::new();
            let mut stack = vec![node];
            visited[node.index()] = true;

            while let Some(current) = stack.pop() {
                component.push(current);
                for neighbor in self.neighbors(current) {
                    if in_subset[neighbor.node.index()] && !visited[neighbor.node.index()] {
                        visited[neighbor.node.index()] = true;
                        stack.push(neighbor.node);
                    }
                }
            }

            component.sort_unstable();
            components.push(component);
        }

        components
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::ConnectedComponentsAlgorithm::Bfs;
    use crate::graph::{Graph, NodeId};

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    #[rstest]
    #[case::empty(0, vec![], vec![])]
    #[case::isolated(3, vec![], vec![vec![n(0)], vec![n(1)], vec![n(2)]])]
    #[case::single_edge(2, vec![[0, 1]], vec![vec![n(0), n(1)]])]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], vec![vec![n(0), n(1), n(2)]])]
    #[case::two_components(4, vec![[0, 1], [2, 3]], vec![vec![n(0), n(1)], vec![n(2), n(3)]])]
    fn test_graph_connected_components(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected: Vec<Vec<NodeId>>,
    ) {
        let g = Graph::new(node_count, &edges);
        assert_eq!(g.connected_components(Bfs), expected);
    }

    #[rstest]
    #[case::empty_subset(4, vec![[0, 1], [2, 3]], vec![], vec![])]
    #[case::full_subset(
        4,
        vec![[0, 1], [2, 3]],
        vec![n(0), n(1), n(2), n(3)],
        vec![vec![n(0), n(1)], vec![n(2), n(3)]]
    )]
    #[case::split_path_by_excluding_middle(
        4,
        vec![[0, 1], [1, 2], [2, 3]],
        vec![n(0), n(1), n(3)],
        vec![vec![n(0), n(1)], vec![n(3)]]
    )]
    #[case::single_component_of_two(
        4,
        vec![[0, 1], [2, 3]],
        vec![n(2), n(3)],
        vec![vec![n(2), n(3)]]
    )]
    #[case::isolated_subset_of_triangle(
        3,
        vec![[0, 1], [1, 2], [0, 2]],
        vec![n(0), n(2)],
        vec![vec![n(0), n(2)]]
    )]
    fn test_graph_connected_components_in(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] subset: Vec<NodeId>,
        #[case] expected: Vec<Vec<NodeId>>,
    ) {
        let g = Graph::new(node_count, &edges);
        assert_eq!(g.connected_components_in(&subset, Bfs), expected);
    }
}
