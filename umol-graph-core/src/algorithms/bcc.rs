//! Biconnected components and articulation point detection.

use std::collections::HashSet;

use crate::graph::{Graph, NodeId};

impl Graph {
    /// Biconnected components with at least 3 nodes (ring-containing).
    pub fn biconnected_components(&self) -> Vec<Vec<NodeId>> {
        let bound = self.node_bound();
        let mut disc: Vec<Option<u32>> = vec![None; bound];
        let mut low: Vec<u32> = vec![0; bound];
        let mut timer: u32 = 0;
        let mut edge_stack: Vec<(NodeId, NodeId)> = Vec::new();
        let mut components: Vec<Vec<NodeId>> = Vec::new();

        for node in self.node_ids() {
            if disc[node.index()].is_some() {
                continue;
            }
            self.bcc_dfs(
                node,
                None,
                &mut disc,
                &mut low,
                &mut timer,
                &mut edge_stack,
                &mut components,
            );
        }

        components
    }

    #[allow(clippy::too_many_arguments)]
    fn bcc_dfs(
        &self,
        u: NodeId,
        parent: Option<NodeId>,
        disc: &mut [Option<u32>],
        low: &mut [u32],
        timer: &mut u32,
        edge_stack: &mut Vec<(NodeId, NodeId)>,
        components: &mut Vec<Vec<NodeId>>,
    ) {
        disc[u.index()] = Some(*timer);
        low[u.index()] = *timer;
        *timer += 1;

        let mut child_count = 0u32;
        for neighbor in self.neighbors(u) {
            let v = neighbor.node;
            if disc[v.index()].is_none() {
                child_count += 1;
                edge_stack.push((u, v));
                self.bcc_dfs(v, Some(u), disc, low, timer, edge_stack, components);

                if low[v.index()] < low[u.index()] {
                    low[u.index()] = low[v.index()];
                }

                let is_articulation = match parent {
                    None => child_count > 1,
                    Some(_) => low[v.index()] >= disc[u.index()].expect("u discovered"),
                };
                if is_articulation {
                    let mut component_nodes = HashSet::new();
                    while let Some((a, b)) = edge_stack.pop() {
                        component_nodes.insert(a);
                        component_nodes.insert(b);
                        if (a == u && b == v) || (a == v && b == u) {
                            break;
                        }
                    }
                    if component_nodes.len() >= 3 {
                        let mut component: Vec<NodeId> = component_nodes.into_iter().collect();
                        component.sort_unstable();
                        components.push(component);
                    }
                }
            } else if Some(v) != parent
                && disc[v.index()].expect("v discovered") < disc[u.index()].expect("u discovered")
            {
                edge_stack.push((u, v));
                let disc_v = disc[v.index()].expect("v discovered");
                if disc_v < low[u.index()] {
                    low[u.index()] = disc_v;
                }
            }
        }

        if parent.is_none() && !edge_stack.is_empty() {
            let mut component_nodes = HashSet::new();
            while let Some((a, b)) = edge_stack.pop() {
                component_nodes.insert(a);
                component_nodes.insert(b);
            }
            if component_nodes.len() >= 3 {
                let mut component: Vec<NodeId> = component_nodes.into_iter().collect();
                component.sort_unstable();
                components.push(component);
            }
        }
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
    fn test_graph_biconnected_components_empty() {
        let g = Graph::default();
        assert!(g.biconnected_components().is_empty());
    }

    #[test]
    fn test_graph_biconnected_components_chain() {
        let g = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        assert!(g.biconnected_components().is_empty());
    }

    #[rstest]
    #[case::single_cycle(
        4,
        vec![[0, 1], [1, 2], [2, 3], [3, 0]],
        vec![vec![n(0), n(1), n(2), n(3)]]
    )]
    #[case::hexagon(
        6,
        vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]],
        vec![vec![n(0), n(1), n(2), n(3), n(4), n(5)]]
    )]
    fn test_graph_biconnected_components(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected: Vec<Vec<NodeId>>,
    ) {
        let g = Graph::new(node_count, &edges);
        assert_eq!(g.biconnected_components(), expected);
    }

    #[test]
    fn test_graph_biconnected_components_articulation() {
        let g = Graph::new(
            5,
            &[[0, 1], [1, 2], [0, 2], [2, 3], [3, 4], [2, 4]],
        );
        let mut bcc = g.biconnected_components();
        bcc.sort();
        assert_eq!(bcc, vec![vec![n(0), n(1), n(2)], vec![n(2), n(3), n(4)]]);
    }

    #[test]
    fn test_graph_biconnected_components_disconnected() {
        let g = Graph::new(5, &[[0, 1], [1, 2], [0, 2], [3, 4]]);
        assert_eq!(
            g.biconnected_components(),
            vec![vec![n(0), n(1), n(2)]]
        );
    }
}
