//! Shortest-path state used to construct relevant-cycle candidates.

use std::collections::{HashSet, VecDeque};

use super::Cycle;
use crate::graph::{EdgeId, Graph, Neighbor, NodeId};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct ShortestPath {
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
}

impl ShortestPath {
    pub(super) fn common_prefix_len(&self, other: &Self) -> usize {
        self.nodes
            .iter()
            .zip(&other.nodes)
            .take_while(|(left, right)| left == right)
            .count()
    }

    pub(super) fn cycle_with(
        &self,
        graph: &Graph,
        other: &Self,
        closing_edge: EdgeId,
    ) -> Option<Cycle> {
        let shared_len = self.common_prefix_len(other);
        if shared_len == 0 {
            return None;
        }

        let first_tail: HashSet<_> = self.nodes[shared_len..].iter().copied().collect();
        if other.nodes[shared_len..]
            .iter()
            .any(|node| first_tail.contains(node))
        {
            return None;
        }

        let first = *self.nodes.last()?;
        let second = *other.nodes.last()?;
        let endpoints = graph.edge_endpoints(closing_edge);
        if !((endpoints == [first, second]) || (endpoints == [second, first])) {
            return None;
        }

        let mut nodes = self.nodes[shared_len - 1..].to_vec();
        nodes.extend(other.nodes[shared_len..].iter().rev());
        if nodes.len() < 3 {
            return None;
        }

        let mut edges = self.edges[shared_len - 1..].to_vec();
        edges.push(closing_edge);
        edges.extend(other.edges[shared_len - 1..].iter().rev());
        if edges.iter().copied().collect::<HashSet<_>>().len() != edges.len() {
            return None;
        }

        Some(Cycle::normalized(graph, nodes, edges))
    }
}

pub(super) struct ShortestPathDag {
    root: NodeId,
    distances: Vec<Option<usize>>,
    predecessors: Vec<Vec<Neighbor>>,
}

impl ShortestPathDag {
    pub(super) fn new(graph: &Graph, root: NodeId) -> Self {
        let mut distances = vec![None; graph.node_count()];
        let mut predecessors = vec![Vec::new(); graph.node_count()];
        let mut queue = VecDeque::from([root]);
        distances[root.index()] = Some(0);

        while let Some(current) = queue.pop_front() {
            let next_distance = distances[current.index()].expect("queued node has a distance") + 1;
            for neighbor in graph.neighbors(current) {
                match distances[neighbor.node.index()] {
                    None => {
                        distances[neighbor.node.index()] = Some(next_distance);
                        predecessors[neighbor.node.index()].push(Neighbor {
                            node: current,
                            edge: neighbor.edge,
                        });
                        queue.push_back(neighbor.node);
                    }
                    Some(distance) if distance == next_distance => {
                        predecessors[neighbor.node.index()].push(Neighbor {
                            node: current,
                            edge: neighbor.edge,
                        });
                    }
                    Some(_) => {}
                }
            }
        }

        for alternatives in &mut predecessors {
            alternatives.sort_unstable_by_key(|neighbor| (neighbor.node, neighbor.edge));
        }
        Self {
            root,
            distances,
            predecessors,
        }
    }

    pub(super) fn distance(&self, node: NodeId) -> Option<usize> {
        self.distances[node.index()]
    }

    pub(super) fn paths_to(&self, target: NodeId) -> Vec<ShortestPath> {
        if self.distance(target).is_none() {
            return Vec::new();
        }
        let mut paths = Vec::new();
        let mut nodes = vec![target];
        let mut edges = Vec::new();
        self.reconstruct(target, &mut nodes, &mut edges, &mut paths);
        paths
    }

    pub(super) fn path_to(&self, target: NodeId) -> Option<ShortestPath> {
        self.distance(target)?;
        let mut nodes = vec![target];
        let mut edges = Vec::new();
        let mut current = target;
        while current != self.root {
            let predecessor = self.predecessors[current.index()].first()?;
            nodes.push(predecessor.node);
            edges.push(predecessor.edge);
            current = predecessor.node;
        }
        nodes.reverse();
        edges.reverse();
        Some(ShortestPath { nodes, edges })
    }

    fn reconstruct(
        &self,
        current: NodeId,
        nodes: &mut Vec<NodeId>,
        edges: &mut Vec<EdgeId>,
        paths: &mut Vec<ShortestPath>,
    ) {
        if current == self.root {
            paths.push(ShortestPath {
                nodes: nodes.iter().rev().copied().collect(),
                edges: edges.iter().rev().copied().collect(),
            });
            return;
        }

        for predecessor in &self.predecessors[current.index()] {
            nodes.push(predecessor.node);
            edges.push(predecessor.edge);
            self.reconstruct(predecessor.node, nodes, edges, paths);
            edges.pop();
            nodes.pop();
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::{ShortestPath, ShortestPathDag};
    use crate::algorithms::cycles::Cycle;
    use crate::graph::{EdgeId, Graph, NodeId};

    #[rstest]
    #[case::connected(
        Graph::new(5, &[[0, 1], [0, 2], [1, 3], [2, 3]]),
        NodeId(0),
        vec![Some(0), Some(1), Some(1), Some(2), None],
    )]
    #[case::parallel(
        Graph::new(2, &[[0, 1], [0, 1]]),
        NodeId(0),
        vec![Some(0), Some(1)],
    )]
    fn test_shortest_path_dag_distance(
        #[case] graph: Graph,
        #[case] root: NodeId,
        #[case] expected: Vec<Option<usize>>,
    ) {
        let dag = ShortestPathDag::new(&graph, root);
        let actual = graph
            .node_ids()
            .map(|node| dag.distance(node))
            .collect::<Vec<_>>();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::alternatives(
        Graph::new(4, &[[0, 1], [0, 2], [1, 3], [2, 3]]),
        NodeId(0),
        NodeId(3),
        vec![
            ShortestPath {
                nodes: vec![NodeId(0), NodeId(1), NodeId(3)],
                edges: vec![EdgeId(0), EdgeId(2)],
            },
            ShortestPath {
                nodes: vec![NodeId(0), NodeId(2), NodeId(3)],
                edges: vec![EdgeId(1), EdgeId(3)],
            },
        ],
    )]
    #[case::parallel(
        Graph::new(2, &[[0, 1], [0, 1]]),
        NodeId(0),
        NodeId(1),
        vec![
            ShortestPath {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(0)],
            },
            ShortestPath {
                nodes: vec![NodeId(0), NodeId(1)],
                edges: vec![EdgeId(1)],
            },
        ],
    )]
    #[case::unreachable(
        Graph::new(2, &[]),
        NodeId(0),
        NodeId(1),
        vec![],
    )]
    fn test_shortest_path_dag_paths_to(
        #[case] graph: Graph,
        #[case] root: NodeId,
        #[case] target: NodeId,
        #[case] expected: Vec<ShortestPath>,
    ) {
        assert_eq!(
            ShortestPathDag::new(&graph, root).paths_to(target),
            expected
        );
    }

    #[rstest]
    #[case::selected(
        Graph::new(4, &[[0, 1], [0, 2], [1, 3], [2, 3]]),
        NodeId(0),
        NodeId(3),
        Some(ShortestPath {
            nodes: vec![NodeId(0), NodeId(1), NodeId(3)],
            edges: vec![EdgeId(0), EdgeId(2)],
        }),
    )]
    #[case::unreachable(Graph::new(2, &[]), NodeId(0), NodeId(1), None)]
    fn test_shortest_path_dag_path_to(
        #[case] graph: Graph,
        #[case] root: NodeId,
        #[case] target: NodeId,
        #[case] expected: Option<ShortestPath>,
    ) {
        assert_eq!(ShortestPathDag::new(&graph, root).path_to(target), expected);
    }

    #[rstest]
    #[case::triangle(
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(1)],
            edges: vec![EdgeId(0)],
        },
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(2)],
            edges: vec![EdgeId(2)],
        },
        EdgeId(1),
        Some(Cycle {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1), EdgeId(2)],
        }),
    )]
    #[case::shared_prefix(
        Graph::new(4, &[[0, 1], [1, 2], [1, 3], [2, 3]]),
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1)],
        },
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(1), NodeId(3)],
            edges: vec![EdgeId(0), EdgeId(2)],
        },
        EdgeId(3),
        Some(Cycle {
            nodes: vec![NodeId(1), NodeId(2), NodeId(3)],
            edges: vec![EdgeId(1), EdgeId(3), EdgeId(2)],
        }),
    )]
    #[case::intersecting_tails(
        Graph::new(5, &[[0, 1], [1, 2], [0, 3], [3, 1], [1, 4], [2, 4]]),
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(1), NodeId(2)],
            edges: vec![EdgeId(0), EdgeId(1)],
        },
        ShortestPath {
            nodes: vec![NodeId(0), NodeId(3), NodeId(1), NodeId(4)],
            edges: vec![EdgeId(2), EdgeId(3), EdgeId(4)],
        },
        EdgeId(5),
        None,
    )]
    fn test_shortest_path_cycle_with(
        #[case] graph: Graph,
        #[case] path: ShortestPath,
        #[case] other: ShortestPath,
        #[case] closing_edge: EdgeId,
        #[case] expected: Option<Cycle>,
    ) {
        assert_eq!(path.cycle_with(&graph, &other, closing_edge), expected);
    }
}
