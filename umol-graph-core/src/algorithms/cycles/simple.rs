//! Bounded edge-aware adaptation of Read--Tarjan simple-cycle enumeration.
//!
//! Every extension retains a bounded return path to the start while avoiding
//! the current path. Chains with only one fruitful extension are advanced
//! without recursion. The minimum node fixes cycle rotation, and the ordering
//! of the first and closing edges fixes reversal while preserving distinct
//! parallel-edge cycles.
//!
//! See [Read and Tarjan, *Bounds on Backtrack Algorithms for Listing Cycles,
//! Paths, and Spanning Trees*
//! (1975)](https://doi.org/10.1002/net.1975.5.3.237).

use std::collections::VecDeque;
use std::ops::ControlFlow;

use super::Cycle;
use crate::graph::{EdgeId, Graph, Neighbor, NodeId};

#[derive(Clone, Copy)]
enum PathExtension {
    Close(EdgeId),
    Extend(Neighbor),
}

struct Path {
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
    contains_node: Vec<bool>,
    contains_edge: Vec<bool>,
}

impl Path {
    fn new(start: NodeId, node_count: usize, edge_count: usize) -> Self {
        let mut contains_node = vec![false; node_count];
        contains_node[start.index()] = true;
        Self {
            nodes: vec![start],
            edges: Vec::new(),
            contains_node,
            contains_edge: vec![false; edge_count],
        }
    }

    fn push(&mut self, neighbor: Neighbor) {
        self.nodes.push(neighbor.node);
        self.edges.push(neighbor.edge);
        self.contains_node[neighbor.node.index()] = true;
        self.contains_edge[neighbor.edge.index()] = true;
    }

    fn truncate(&mut self, node_count: usize) {
        while self.nodes.len() > node_count {
            let node = self.nodes.pop().expect("extended path contains a node");
            let edge = self.edges.pop().expect("extended path contains an edge");
            self.contains_node[node.index()] = false;
            self.contains_edge[edge.index()] = false;
        }
    }

    fn current(&self) -> NodeId {
        *self.nodes.last().expect("a path has a current node")
    }
}

impl Graph {
    pub(super) fn visit_simple_cycles_read_tarjan<B, F>(
        &self,
        max_cycle_size: usize,
        visitor: &mut F,
    ) -> ControlFlow<B>
    where
        F: FnMut(Cycle) -> ControlFlow<B>,
    {
        if max_cycle_size == 0 {
            return ControlFlow::Continue(());
        }

        let mut loops = vec![Vec::new(); self.node_count()];
        for edge in self.edge_ids() {
            let [first, second] = self.edge_endpoints(edge);
            if first == second {
                loops[first.index()].push(edge);
            }
        }

        for start in self.node_ids() {
            for &edge in &loops[start.index()] {
                if let ControlFlow::Break(value) =
                    visitor(Cycle::normalized(self, vec![start], vec![edge]))
                {
                    return ControlFlow::Break(value);
                }
            }
            if max_cycle_size < 2 {
                continue;
            }

            let mut path = Path::new(start, self.node_count(), self.edge_count());
            if let ControlFlow::Break(value) =
                visit_from(self, start, max_cycle_size, &mut path, visitor)
            {
                return ControlFlow::Break(value);
            }
        }

        ControlFlow::Continue(())
    }
}

fn visit_from<B, F>(
    graph: &Graph,
    start: NodeId,
    max_cycle_size: usize,
    path: &mut Path,
    visitor: &mut F,
) -> ControlFlow<B>
where
    F: FnMut(Cycle) -> ControlFlow<B>,
{
    for extension in extensions(graph, start, max_cycle_size, path) {
        let node_count = path.nodes.len();
        if let ControlFlow::Break(value) =
            follow_extension(graph, start, max_cycle_size, extension, path, visitor)
        {
            return ControlFlow::Break(value);
        }
        path.truncate(node_count);
    }
    ControlFlow::Continue(())
}

fn follow_extension<B, F>(
    graph: &Graph,
    start: NodeId,
    max_cycle_size: usize,
    mut extension: PathExtension,
    path: &mut Path,
    visitor: &mut F,
) -> ControlFlow<B>
where
    F: FnMut(Cycle) -> ControlFlow<B>,
{
    loop {
        match extension {
            PathExtension::Close(edge) => {
                let mut cycle_edges = path.edges.clone();
                cycle_edges.push(edge);
                return visitor(Cycle::normalized(graph, path.nodes.clone(), cycle_edges));
            }
            PathExtension::Extend(neighbor) => {
                path.push(neighbor);

                let mut next = extensions(graph, start, max_cycle_size, path);
                if next.len() == 1 {
                    extension = next.pop().expect("one extension was found");
                } else {
                    return if next.is_empty() {
                        ControlFlow::Continue(())
                    } else {
                        visit_from(graph, start, max_cycle_size, path, visitor)
                    };
                }
            }
        }
    }
}

fn extensions(
    graph: &Graph,
    start: NodeId,
    max_cycle_size: usize,
    path: &Path,
) -> Vec<PathExtension> {
    let mut neighbors = graph.neighbors(path.current()).to_vec();
    neighbors.sort_unstable_by_key(|neighbor| (neighbor.node, neighbor.edge));

    let mut result = Vec::new();
    for neighbor in neighbors {
        if path.contains_edge[neighbor.edge.index()] {
            continue;
        }
        if neighbor.node == start {
            if path.nodes.len() >= 2
                && path
                    .edges
                    .first()
                    .is_some_and(|&first| first < neighbor.edge)
            {
                result.push(PathExtension::Close(neighbor.edge));
            }
            continue;
        }
        if neighbor.node < start
            || path.contains_node[neighbor.node.index()]
            || path.nodes.len() >= max_cycle_size
        {
            continue;
        }

        let prospective_size = path.nodes.len() + 1;
        let max_return_length = max_cycle_size
            .saturating_sub(prospective_size)
            .saturating_add(1);
        if has_return_path(
            graph,
            neighbor.node,
            start,
            max_return_length,
            path,
            neighbor.edge,
        ) {
            result.push(PathExtension::Extend(neighbor));
        }
    }
    result
}

fn has_return_path(
    graph: &Graph,
    from: NodeId,
    start: NodeId,
    max_length: usize,
    path: &Path,
    extension_edge: EdgeId,
) -> bool {
    let mut visited = vec![false; graph.node_count()];
    let mut queue = VecDeque::from([(from, 0usize)]);
    visited[from.index()] = true;

    while let Some((current, distance)) = queue.pop_front() {
        if distance == max_length {
            continue;
        }
        for neighbor in graph.neighbors(current) {
            if neighbor.edge == extension_edge || path.contains_edge[neighbor.edge.index()] {
                continue;
            }
            if neighbor.node == start {
                return true;
            }
            if path.contains_node[neighbor.node.index()] || visited[neighbor.node.index()] {
                continue;
            }
            visited[neighbor.node.index()] = true;
            queue.push_back((neighbor.node, distance + 1));
        }
    }
    false
}
