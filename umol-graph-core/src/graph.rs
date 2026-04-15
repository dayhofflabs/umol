use std::ops::{Index, IndexMut};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeId(pub u32);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct EdgeId(pub u32);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Neighbor {
    pub node: NodeId,
    pub edge: EdgeId,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct EdgeData<E> {
    endpoints: [NodeId; 2],
    data: E,
}

impl<N: PartialEq, E: PartialEq> PartialEq for Graph<N, E> {
    fn eq(&self, other: &Self) -> bool {
        self.node_count == other.node_count
            && self.edge_count == other.edge_count
            && self.nodes == other.nodes
            && self.edges == other.edges
    }
}

impl<N: Eq, E: Eq> Eq for Graph<N, E> {}

/// Undirected property graph with stable indices.
///
/// Nodes and edges are stored in slot arrays. Removed slots are recycled
/// via a free list, so identifiers remain valid across mutations.
/// Adjacency is per-node `Vec<Neighbor>`, giving cache-friendly iteration
/// for bounded-degree graphs.
#[derive(Clone, Debug)]
pub struct Graph<N, E> {
    nodes: Vec<Option<N>>,
    edges: Vec<Option<EdgeData<E>>>,
    adjacency: Vec<Vec<Neighbor>>,
    node_count: usize,
    edge_count: usize,
    free_nodes: Vec<NodeId>,
    free_edges: Vec<EdgeId>,
}

impl NodeId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

impl EdgeId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

impl<N, E> Graph<N, E> {
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
            adjacency: Vec::new(),
            node_count: 0,
            edge_count: 0,
            free_nodes: Vec::new(),
            free_edges: Vec::new(),
        }
    }

    pub fn with_capacity(node_capacity: usize, edge_capacity: usize) -> Self {
        Self {
            nodes: Vec::with_capacity(node_capacity),
            edges: Vec::with_capacity(edge_capacity),
            adjacency: Vec::with_capacity(node_capacity),
            node_count: 0,
            edge_count: 0,
            free_nodes: Vec::new(),
            free_edges: Vec::new(),
        }
    }

    pub fn add_node(&mut self, data: N) -> NodeId {
        if let Some(id) = self.free_nodes.pop() {
            self.nodes[id.index()] = Some(data);
            self.node_count += 1;
            id
        } else {
            let id = NodeId(self.nodes.len() as u32);
            self.nodes.push(Some(data));
            self.adjacency.push(Vec::new());
            self.node_count += 1;
            id
        }
    }

    pub fn add_edge(&mut self, a: NodeId, b: NodeId, data: E) -> EdgeId {
        debug_assert!(self.contains_node(a), "node {a:?} does not exist");
        debug_assert!(self.contains_node(b), "node {b:?} does not exist");

        let id = if let Some(id) = self.free_edges.pop() {
            self.edges[id.index()] = Some(EdgeData {
                endpoints: [a, b],
                data,
            });
            id
        } else {
            let id = EdgeId(self.edges.len() as u32);
            self.edges.push(Some(EdgeData {
                endpoints: [a, b],
                data,
            }));
            id
        };

        self.adjacency[a.index()].push(Neighbor { node: b, edge: id });
        self.adjacency[b.index()].push(Neighbor { node: a, edge: id });
        self.edge_count += 1;
        id
    }

    pub fn remove_node(&mut self, id: NodeId) -> Option<N> {
        let data = self.nodes.get_mut(id.index())?.take()?;
        let incident: Vec<EdgeId> = self.adjacency[id.index()]
            .drain(..)
            .map(|n| n.edge)
            .collect();
        for edge_id in incident {
            if let Some(edge_data) = self.edges[edge_id.index()].take() {
                let [a, b] = edge_data.endpoints;
                let other = if a == id { b } else { a };
                if other != id {
                    self.adjacency[other.index()].retain(|n| n.edge != edge_id);
                }
                self.free_edges.push(edge_id);
                self.edge_count -= 1;
            }
        }
        self.free_nodes.push(id);
        self.node_count -= 1;
        Some(data)
    }

    pub fn remove_edge(&mut self, id: EdgeId) -> Option<E> {
        let edge = self.edges.get_mut(id.index())?.take()?;
        let [a, b] = edge.endpoints;
        self.adjacency[a.index()].retain(|n| n.edge != id);
        if a != b {
            self.adjacency[b.index()].retain(|n| n.edge != id);
        }
        self.free_edges.push(id);
        self.edge_count -= 1;
        Some(edge.data)
    }

    pub fn contains_node(&self, id: NodeId) -> bool {
        self.nodes.get(id.index()).is_some_and(|n| n.is_some())
    }

    pub fn contains_edge(&self, id: EdgeId) -> bool {
        self.edges.get(id.index()).is_some_and(|e| e.is_some())
    }

    pub fn node(&self, id: NodeId) -> Option<&N> {
        self.nodes.get(id.index())?.as_ref()
    }

    pub fn node_mut(&mut self, id: NodeId) -> Option<&mut N> {
        self.nodes.get_mut(id.index())?.as_mut()
    }

    pub fn edge(&self, id: EdgeId) -> Option<&E> {
        self.edges
            .get(id.index())
            .and_then(|e| e.as_ref().map(|e| &e.data))
    }

    pub fn edge_mut(&mut self, id: EdgeId) -> Option<&mut E> {
        self.edges
            .get_mut(id.index())
            .and_then(|e| e.as_mut().map(|e| &mut e.data))
    }

    pub fn edge_endpoints(&self, id: EdgeId) -> Option<[NodeId; 2]> {
        self.edges
            .get(id.index())
            .and_then(|e| e.as_ref().map(|e| e.endpoints))
    }

    pub fn neighbors(&self, id: NodeId) -> &[Neighbor] {
        &self.adjacency[id.index()]
    }

    pub fn degree(&self, id: NodeId) -> usize {
        self.adjacency[id.index()].len()
    }

    pub fn find_edge(&self, a: NodeId, b: NodeId) -> Option<EdgeId> {
        self.neighbors(a)
            .iter()
            .find(|n| n.node == b)
            .map(|n| n.edge)
    }

    pub fn node_count(&self) -> usize {
        self.node_count
    }

    pub fn edge_count(&self) -> usize {
        self.edge_count
    }

    /// Upper bound on node indices. All valid NodeId values satisfy
    /// `id.index() < node_bound()`.
    pub fn node_bound(&self) -> usize {
        self.nodes.len()
    }

    /// Upper bound on edge indices. All valid EdgeId values satisfy
    /// `id.index() < edge_bound()`.
    pub fn edge_bound(&self) -> usize {
        self.edges.len()
    }

    pub fn node_ids(&self) -> impl Iterator<Item = NodeId> + '_ {
        self.nodes
            .iter()
            .enumerate()
            .filter_map(|(i, n)| n.as_ref().map(|_| NodeId(i as u32)))
    }

    pub fn edge_ids(&self) -> impl Iterator<Item = EdgeId> + '_ {
        self.edges
            .iter()
            .enumerate()
            .filter_map(|(i, e)| e.as_ref().map(|_| EdgeId(i as u32)))
    }

    /// True when the node array has no holes (all slots occupied).
    /// When true, node ids are dense `0..node_count`.
    pub fn is_dense(&self) -> bool {
        self.node_count == self.nodes.len()
    }
}

impl<N: Default, E> Graph<N, E> {
    pub fn from_edges(node_count: usize, edges: Vec<(u32, u32, E)>) -> Self {
        let mut g = Self::with_capacity(node_count, edges.len());
        for _ in 0..node_count {
            g.add_node(N::default());
        }
        for (a, b, data) in edges {
            g.add_edge(NodeId(a), NodeId(b), data);
        }
        g
    }
}

impl<N, E> Default for Graph<N, E> {
    fn default() -> Self {
        Self::new()
    }
}

impl<N, E> Index<NodeId> for Graph<N, E> {
    type Output = N;
    fn index(&self, id: NodeId) -> &N {
        self.node(id).expect("node does not exist")
    }
}

impl<N, E> IndexMut<NodeId> for Graph<N, E> {
    fn index_mut(&mut self, id: NodeId) -> &mut N {
        self.node_mut(id).expect("node does not exist")
    }
}

impl<N, E> Index<EdgeId> for Graph<N, E> {
    type Output = E;
    fn index(&self, id: EdgeId) -> &E {
        self.edge(id).expect("edge does not exist")
    }
}

impl<N, E> IndexMut<EdgeId> for Graph<N, E> {
    fn index_mut(&mut self, id: EdgeId) -> &mut E {
        self.edge_mut(id).expect("edge does not exist")
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[test]
    fn test_graph_new() {
        let g = Graph::<i32, ()>::new();
        assert_eq!(g.node_count(), 0);
        assert_eq!(g.edge_count(), 0);
        assert!(g.is_dense());
    }

    #[test]
    fn test_graph_add_node() {
        let mut g = Graph::<&str, ()>::new();
        let a = g.add_node("carbon");
        let b = g.add_node("oxygen");
        assert_eq!(g.node_count(), 2);
        assert_eq!(g[a], "carbon");
        assert_eq!(g[b], "oxygen");
        assert_eq!(a, NodeId(0));
        assert_eq!(b, NodeId(1));
    }

    #[test]
    fn test_graph_add_edge() {
        let mut g = Graph::<(), u8>::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let e = g.add_edge(a, b, 2);
        assert_eq!(g.edge_count(), 1);
        assert_eq!(g[e], 2);
        assert_eq!(g.edge_endpoints(e), Some([a, b]));
        assert_eq!(g.degree(a), 1);
        assert_eq!(g.degree(b), 1);
        assert_eq!(g.neighbors(a)[0].node, b);
        assert_eq!(g.neighbors(b)[0].node, a);
    }

    #[test]
    fn test_graph_remove_edge() {
        let mut g = Graph::<(), ()>::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let e = g.add_edge(a, b, ());
        assert_eq!(g.remove_edge(e), Some(()));
        assert_eq!(g.edge_count(), 0);
        assert_eq!(g.degree(a), 0);
        assert_eq!(g.degree(b), 0);
        assert!(!g.contains_edge(e));
    }

    #[test]
    fn test_graph_remove_node() {
        let mut g = Graph::<(), ()>::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let c = g.add_node(());
        g.add_edge(a, b, ());
        g.add_edge(b, c, ());
        g.remove_node(b);
        assert_eq!(g.node_count(), 2);
        assert_eq!(g.edge_count(), 0);
        assert!(g.contains_node(a));
        assert!(!g.contains_node(b));
        assert!(g.contains_node(c));
        assert_eq!(g.degree(a), 0);
        assert_eq!(g.degree(c), 0);
    }

    #[test]
    fn test_graph_free_list_reuse() {
        let mut g = Graph::<i32, ()>::new();
        let a = g.add_node(10);
        let b = g.add_node(20);
        g.remove_node(a);
        let c = g.add_node(30);
        // Reuses slot 0
        assert_eq!(c, NodeId(0));
        assert_eq!(g[c], 30);
        assert_eq!(g.node_count(), 2);
        assert!(g.is_dense());

        // Edge free list
        let e1 = g.add_edge(b, c, ());
        g.remove_edge(e1);
        let e2 = g.add_edge(c, b, ());
        assert_eq!(e2, EdgeId(0));
    }

    #[test]
    fn test_graph_find_edge() {
        let mut g = Graph::<(), u8>::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let c = g.add_node(());
        let e = g.add_edge(a, b, 1);
        assert_eq!(g.find_edge(a, b), Some(e));
        assert_eq!(g.find_edge(b, a), Some(e));
        assert_eq!(g.find_edge(a, c), None);
    }

    #[test]
    fn test_graph_self_loop() {
        let mut g = Graph::<(), ()>::new();
        let a = g.add_node(());
        let e = g.add_edge(a, a, ());
        assert_eq!(g.degree(a), 2);
        assert_eq!(g.edge_count(), 1);
        g.remove_edge(e);
        assert_eq!(g.degree(a), 0);
        assert_eq!(g.edge_count(), 0);
    }

    #[test]
    fn test_graph_remove_node_self_loop() {
        let mut g = Graph::<(), ()>::new();
        let a = g.add_node(());
        g.add_edge(a, a, ());
        g.remove_node(a);
        assert_eq!(g.node_count(), 0);
        assert_eq!(g.edge_count(), 0);
    }

    #[test]
    fn test_graph_node_ids() {
        let mut g = Graph::<(), ()>::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let c = g.add_node(());
        g.remove_node(b);
        let ids: Vec<NodeId> = g.node_ids().collect();
        assert_eq!(ids, vec![a, c]);
    }

    #[test]
    fn test_graph_edge_ids() {
        let mut g = Graph::<(), ()>::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let c = g.add_node(());
        let e1 = g.add_edge(a, b, ());
        let e2 = g.add_edge(b, c, ());
        g.remove_edge(e1);
        let ids: Vec<EdgeId> = g.edge_ids().collect();
        assert_eq!(ids, vec![e2]);
    }

    #[test]
    fn test_graph_node_mut() {
        let mut g = Graph::<i32, ()>::new();
        let a = g.add_node(0);
        g[a] = 99;
        assert_eq!(g[a], 99);
    }

    #[test]
    fn test_graph_edge_mut() {
        let mut g = Graph::<(), i32>::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let e = g.add_edge(a, b, 1);
        g[e] = 3;
        assert_eq!(g[e], 3);
    }

    #[rstest]
    #[case::empty(0, vec![], 0, 0)]
    #[case::single_node(1, vec![], 1, 0)]
    #[case::triangle(3, vec![(0, 1, ()), (1, 2, ()), (0, 2, ())], 3, 3)]
    fn test_graph_from_edges(
        #[case] node_count: usize,
        #[case] edges: Vec<(u32, u32, ())>,
        #[case] expected_nodes: usize,
        #[case] expected_edges: usize,
    ) {
        let g = Graph::<(), _>::from_edges(node_count, edges);
        assert_eq!(g.node_count(), expected_nodes);
        assert_eq!(g.edge_count(), expected_edges);
    }

    #[test]
    fn test_graph_parallel_edges() {
        let mut g = Graph::<(), u8>::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let e1 = g.add_edge(a, b, 1);
        let e2 = g.add_edge(a, b, 2);
        assert_eq!(g.edge_count(), 2);
        assert_eq!(g.degree(a), 2);
        assert_eq!(g[e1], 1);
        assert_eq!(g[e2], 2);
    }

    #[test]
    fn test_graph_not_dense_after_removal() {
        let mut g = Graph::<(), ()>::new();
        g.add_node(());
        let b = g.add_node(());
        g.add_node(());
        g.remove_node(b);
        assert!(!g.is_dense());
    }

    #[test]
    fn test_graph_remove_nonexistent() {
        let mut g = Graph::<(), ()>::new();
        assert_eq!(g.remove_node(NodeId(0)), None);
        assert_eq!(g.remove_edge(EdgeId(0)), None);
    }
}
