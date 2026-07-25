//! Bipartition and bipartiteness.
//!
//! [`Graph::bipartition`] currently computes a two-coloring by breadth-first
//! search; [`Graph::is_bipartite`] is its boolean projection. This is not a
//! general graph-coloring module.

use std::collections::VecDeque;

use crate::graph::Graph;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BipartitionAlgorithm {
    Bfs,
}

impl Graph {
    /// Returns a 2-coloring (one bool per vertex) iff the graph is bipartite,
    /// else `None`. Isolated components' colorings are independent of one
    /// another up to flipping; the function picks one canonical assignment by
    /// starting BFS at the lowest-indexed unvisited node with color `false`.
    pub fn bipartition(&self, alg: BipartitionAlgorithm) -> Option<Vec<bool>> {
        match alg {
            BipartitionAlgorithm::Bfs => self.bipartition_bfs(),
        }
    }

    pub fn is_bipartite(&self, alg: BipartitionAlgorithm) -> bool {
        self.bipartition(alg).is_some()
    }

    // BFS 2-coloring, restarting at each unvisited component. O(V+E).
    fn bipartition_bfs(&self) -> Option<Vec<bool>> {
        let n = self.node_count();
        let mut colors: Vec<Option<bool>> = vec![None; n];

        for start in self.node_ids() {
            if colors[start.index()].is_some() {
                continue;
            }
            colors[start.index()] = Some(false);
            let mut queue = VecDeque::new();
            queue.push_back(start);
            while let Some(v) = queue.pop_front() {
                let v_color = colors[v.index()].expect("queued vertex has a color");
                for nbr in self.neighbors(v) {
                    let u = nbr.node;
                    match colors[u.index()] {
                        None => {
                            colors[u.index()] = Some(!v_color);
                            queue.push_back(u);
                        }
                        Some(c) if c == v_color => return None,
                        Some(_) => {}
                    }
                }
            }
        }

        Some(colors.into_iter().map(|c| c.unwrap_or(false)).collect())
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::BipartitionAlgorithm::Bfs;
    use crate::graph::Graph;

    #[rstest]
    #[case::empty(0, vec![], true)]
    #[case::isolated_vertex(1, vec![], true)]
    #[case::single_edge(2, vec![[0, 1]], true)]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], false)]
    #[case::square(4, vec![[0, 1], [1, 2], [2, 3], [3, 0]], true)]
    #[case::pentagon(5, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]], false)]
    #[case::hexagon(6, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]], true)]
    #[case::path_4(4, vec![[0, 1], [1, 2], [2, 3]], true)]
    #[case::k4(4, vec![[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], false)]
    #[case::k_2_3(
        5,
        vec![[0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4]],
        true,
    )]
    #[case::two_components_bipartite(
        4,
        vec![[0, 1], [2, 3]],
        true,
    )]
    #[case::triangle_plus_edge(
        5,
        vec![[0, 1], [1, 2], [0, 2], [3, 4]],
        false,
    )]
    fn test_graph_is_bipartite(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected: bool,
    ) {
        let g = Graph::new(node_count, &edges);
        assert_eq!(g.is_bipartite(Bfs), expected);
    }

    #[rstest]
    fn test_graph_bipartition_bipartite_returns_valid_coloring() {
        let g = Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]);
        let coloring = g.bipartition(Bfs).expect("bipartite");
        assert_eq!(coloring.len(), 6);
        for eid in g.edge_ids() {
            let [a, b] = g.edge_endpoints(eid);
            assert_ne!(
                coloring[a.index()],
                coloring[b.index()],
                "edge {:?}-{:?} has same-color endpoints",
                a,
                b,
            );
        }
    }

    #[rstest]
    fn test_graph_bipartition_non_bipartite_returns_none() {
        let g = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        assert_eq!(g.bipartition(Bfs), None);
    }

    #[rstest]
    fn test_graph_bipartition_disconnected_components_independent() {
        // Two components: {0,1} bipartite and {2,3,4} bipartite.
        let g = Graph::new(5, &[[0, 1], [2, 3], [3, 4]]);
        let coloring = g.bipartition(Bfs).expect("bipartite");
        for eid in g.edge_ids() {
            let [a, b] = g.edge_endpoints(eid);
            assert_ne!(coloring[a.index()], coloring[b.index()]);
        }
    }
}
