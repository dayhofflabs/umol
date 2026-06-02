//! CSR graph with Arc-based structural sharing.
//!
//! `Graph` stores only adjacency (offsets, neighbor lists, edge endpoints).
//! Node and edge data live externally in `Vec`s indexed by `NodeId`/`EdgeId`.
//! The CSR is wrapped in `Arc` for zero-cost cloning; mutations rebuild
//! it and produce a `Remapping` for reindexing external data.

use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeId(pub u32);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct EdgeId(pub u32);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Neighbor {
    pub node: NodeId,
    pub edge: EdgeId,
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

#[derive(Clone, Debug, PartialEq, Eq)]
struct Csr {
    offsets: Vec<u32>,
    neighbors: Vec<Neighbor>,
    endpoints: Vec<[NodeId; 2]>,
    node_count: usize,
    edge_count: usize,
}

/// Undirected graph stored as compressed sparse row (CSR).
///
/// Stores only adjacency structure — node and edge data live externally,
/// indexed by `NodeId` and `EdgeId` positions. The CSR is shared via
/// `Arc`; mutations trigger copy-on-write.
#[derive(Clone, Debug)]
pub struct Graph {
    csr: Arc<Csr>,
}

impl PartialEq for Graph {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.csr, &other.csr) || self.csr == other.csr
    }
}

impl Eq for Graph {}

impl Graph {
    pub fn new(node_count: usize, edges: &[[u32; 2]]) -> Self {
        Self {
            csr: Arc::new(Self::build_csr(node_count, edges)),
        }
    }

    pub fn node_count(&self) -> usize {
        self.csr.node_count
    }

    pub fn edge_count(&self) -> usize {
        self.csr.edge_count
    }

    pub fn node_bound(&self) -> usize {
        self.csr.node_count
    }

    pub fn edge_bound(&self) -> usize {
        self.csr.edge_count
    }

    /// Neighbors sorted by `NodeId`. Enables binary search in `find_edge`
    /// and set-intersection in `induced_edges`.
    pub fn neighbors(&self, id: NodeId) -> &[Neighbor] {
        let start = self.csr.offsets[id.index()] as usize;
        let end = self.csr.offsets[id.index() + 1] as usize;
        &self.csr.neighbors[start..end]
    }

    pub fn degree(&self, id: NodeId) -> usize {
        let start = self.csr.offsets[id.index()] as usize;
        let end = self.csr.offsets[id.index() + 1] as usize;
        end - start
    }

    pub fn edge_endpoints(&self, id: EdgeId) -> [NodeId; 2] {
        self.csr.endpoints[id.index()]
    }

    pub fn find_edge(&self, a: NodeId, b: NodeId) -> Option<EdgeId> {
        let nbrs = self.neighbors(a);
        nbrs.binary_search_by_key(&b, |n| n.node)
            .ok()
            .map(|i| nbrs[i].edge)
    }

    /// Edges whose both endpoints are in `nodes`. The slice must be sorted.
    /// Yields each edge exactly once; iteration order is unspecified.
    pub fn induced_edges<'a>(&'a self, nodes: &'a [NodeId]) -> impl Iterator<Item = EdgeId> + 'a {
        nodes.iter().flat_map(move |&node| {
            self.neighbors(node).iter().filter_map(move |n| {
                if nodes.binary_search(&n.node).is_ok() && self.edge_endpoints(n.edge)[0] == node {
                    Some(n.edge)
                } else {
                    None
                }
            })
        })
    }

    pub fn contains_node(&self, id: NodeId) -> bool {
        id.index() < self.csr.node_count
    }

    pub fn contains_edge(&self, id: EdgeId) -> bool {
        id.index() < self.csr.edge_count
    }

    pub fn node_ids(&self) -> impl Iterator<Item = NodeId> {
        (0..self.csr.node_count as u32).map(NodeId)
    }

    pub fn edge_ids(&self) -> impl Iterator<Item = EdgeId> {
        (0..self.csr.edge_count as u32).map(EdgeId)
    }

    pub fn is_dense(&self) -> bool {
        true
    }

    // --- Mutations (rebuild CSR, CoW via Arc) ---

    pub fn add_node(&mut self) -> NodeId {
        let old = &*self.csr;
        let new_id = NodeId(old.node_count as u32);
        let edges: Vec<[u32; 2]> = old.endpoints.iter().map(|&[a, b]| [a.0, b.0]).collect();
        self.csr = Arc::new(Self::build_csr(old.node_count + 1, &edges));
        new_id
    }

    pub fn add_edge(&mut self, a: NodeId, b: NodeId) -> EdgeId {
        let old = &*self.csr;
        let new_id = EdgeId(old.edge_count as u32);
        let mut edges: Vec<[u32; 2]> = old.endpoints.iter().map(|&[s, t]| [s.0, t.0]).collect();
        edges.push([a.0, b.0]);
        self.csr = Arc::new(Self::build_csr(old.node_count, &edges));
        new_id
    }

    pub fn remove(&mut self, nodes: &[NodeId], edges: &[EdgeId]) -> Remapping {
        let mut removed_nodes: Vec<u32> = nodes.iter().map(|n| n.0).collect();
        removed_nodes.sort_unstable();
        removed_nodes.dedup();

        let mut removed_edge_set: Vec<u32> = edges.iter().map(|e| e.0).collect();

        let old = &*self.csr;

        // Also remove edges incident to removed nodes
        for (i, &[a, b]) in old.endpoints.iter().enumerate() {
            if removed_nodes.binary_search(&a.0).is_ok()
                || removed_nodes.binary_search(&b.0).is_ok()
            {
                removed_edge_set.push(i as u32);
            }
        }
        removed_edge_set.sort_unstable();
        removed_edge_set.dedup();

        let kept_edges: Vec<[u32; 2]> = old
            .endpoints
            .iter()
            .enumerate()
            .filter(|&(i, _)| removed_edge_set.binary_search(&(i as u32)).is_err())
            .map(|(_, &[a, b])| {
                let shift_a = removed_nodes.partition_point(|&r| r < a.0) as u32;
                let shift_b = removed_nodes.partition_point(|&r| r < b.0) as u32;
                [a.0 - shift_a, b.0 - shift_b]
            })
            .collect();

        self.csr = Arc::new(Self::build_csr(
            old.node_count - removed_nodes.len(),
            &kept_edges,
        ));

        Remapping {
            removed_nodes,
            removed_edges: removed_edge_set,
        }
    }

    pub fn remove_node(&mut self, id: NodeId) -> Remapping {
        self.remove(&[id], &[])
    }

    pub fn remove_edge(&mut self, id: EdgeId) -> Remapping {
        self.remove(&[], &[id])
    }

    /// Build an induced subgraph from a subset of nodes.
    ///
    /// Returns the subgraph (with contiguous node/edge IDs starting at 0)
    /// and mappings from new IDs back to original IDs.
    /// Induced subgraph over `nodes` in caller-supplied order. Duplicates in
    /// `nodes` are deduplicated (first occurrence wins). Sub node ids
    /// 0..len(deduplicated nodes) correspond positionally to the deduplicated
    /// input. The returned [`Embedding`] borrows `self`; call
    /// [`Embedding::extract`] to materialize the sub `Graph` when needed.
    pub fn induced_subgraph(&self, nodes: &[NodeId]) -> Embedding<'_> {
        let mut host_nodes: Vec<NodeId> = Vec::with_capacity(nodes.len());
        let mut sub_nodes: HashMap<NodeId, NodeId> = HashMap::with_capacity(nodes.len());
        for &node in nodes {
            if let Entry::Vacant(entry) = sub_nodes.entry(node) {
                entry.insert(NodeId(host_nodes.len() as u32));
                host_nodes.push(node);
            }
        }

        let mut host_edges = Vec::new();
        for eid in self.edge_ids() {
            let [a, b] = self.edge_endpoints(eid);
            if sub_nodes.contains_key(&a) && sub_nodes.contains_key(&b) {
                host_edges.push(eid);
            }
        }

        Embedding {
            host_nodes,
            host_edges,
            sub_nodes,
            graph: self,
        }
    }

    fn build_csr(node_count: usize, edges: &[[u32; 2]]) -> Csr {
        let edge_count = edges.len();

        let mut degree = vec![0u32; node_count];
        for &[a, b] in edges {
            degree[a as usize] += 1;
            degree[b as usize] += 1;
        }

        let mut offsets = Vec::with_capacity(node_count + 1);
        offsets.push(0);
        for &d in &degree {
            offsets.push(offsets.last().unwrap() + d);
        }

        let total = *offsets.last().unwrap() as usize;
        let mut neighbors = vec![
            Neighbor {
                node: NodeId(0),
                edge: EdgeId(0)
            };
            total
        ];
        let mut cursor: Vec<u32> = offsets[..node_count].to_vec();

        for (i, &[a, b]) in edges.iter().enumerate() {
            let eid = EdgeId(i as u32);

            let pos = cursor[a as usize] as usize;
            neighbors[pos] = Neighbor {
                node: NodeId(b),
                edge: eid,
            };
            cursor[a as usize] += 1;

            let pos = cursor[b as usize] as usize;
            neighbors[pos] = Neighbor {
                node: NodeId(a),
                edge: eid,
            };
            cursor[b as usize] += 1;
        }

        for i in 0..node_count {
            let start = offsets[i] as usize;
            let end = offsets[i + 1] as usize;
            neighbors[start..end].sort_unstable_by_key(|n| n.node);
        }

        let endpoints: Vec<[NodeId; 2]> = edges
            .iter()
            .map(|&[a, b]| {
                let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
                [NodeId(lo), NodeId(hi)]
            })
            .collect();

        Csr {
            offsets,
            neighbors,
            endpoints,
            node_count,
            edge_count,
        }
    }
}

impl Default for Graph {
    fn default() -> Self {
        Self::new(0, &[])
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Remapping {
    pub removed_nodes: Vec<u32>,
    pub removed_edges: Vec<u32>,
}

impl Remapping {
    pub fn map_node(&self, old: NodeId) -> Option<NodeId> {
        if self.removed_nodes.binary_search(&old.0).is_ok() {
            return None;
        }
        let shift = self.removed_nodes.partition_point(|&r| r < old.0);
        Some(NodeId(old.0 - shift as u32))
    }

    pub fn map_edge(&self, old: EdgeId) -> Option<EdgeId> {
        if self.removed_edges.binary_search(&old.0).is_ok() {
            return None;
        }
        let shift = self.removed_edges.partition_point(|&r| r < old.0);
        Some(EdgeId(old.0 - shift as u32))
    }

    pub fn unmap_node(&self, post: NodeId) -> NodeId {
        NodeId(unmap_dense(&self.removed_nodes, post.0))
    }

    pub fn unmap_edge(&self, post: EdgeId) -> EdgeId {
        EdgeId(unmap_dense(&self.removed_edges, post.0))
    }

    pub fn apply_to_node_vec<T: Clone>(&self, data: &[T]) -> Vec<T> {
        data.iter()
            .enumerate()
            .filter(|(i, _)| self.removed_nodes.binary_search(&(*i as u32)).is_err())
            .map(|(_, v)| v.clone())
            .collect()
    }

    pub fn apply_to_edge_vec<T: Clone>(&self, data: &[T]) -> Vec<T> {
        data.iter()
            .enumerate()
            .filter(|(i, _)| self.removed_edges.binary_search(&(*i as u32)).is_err())
            .map(|(_, v)| v.clone())
            .collect()
    }

}

// Inverse dense shift: re-add removed ids at or below the post index (fixpoint).
fn unmap_dense(removed: &[u32], post: u32) -> u32 {
    let mut old = post;
    loop {
        let next = post + removed.partition_point(|&r| r <= old) as u32;
        if next == old {
            return old;
        }
        old = next;
    }
}

/// Subgraph induced by a node subset, borrowing the host `Graph`. Holds
/// sub→host index maps and a host→sub inverse for O(1) translation; the sub
/// `Graph` is materialized on demand via [`Embedding::extract`].
#[derive(Clone, Debug)]
pub struct Embedding<'a> {
    host_nodes: Vec<NodeId>,
    host_edges: Vec<EdgeId>,
    sub_nodes: HashMap<NodeId, NodeId>,
    graph: &'a Graph,
}

impl<'a> Embedding<'a> {
    pub fn graph(&self) -> &'a Graph {
        self.graph
    }

    pub fn host_nodes(&self) -> &[NodeId] {
        &self.host_nodes
    }

    pub fn host_edges(&self) -> &[EdgeId] {
        &self.host_edges
    }

    pub fn host_node(&self, sub: NodeId) -> NodeId {
        self.host_nodes[sub.index()]
    }

    pub fn host_edge(&self, sub: EdgeId) -> EdgeId {
        self.host_edges[sub.index()]
    }

    pub fn sub_node(&self, host: NodeId) -> Option<NodeId> {
        self.sub_nodes.get(&host).copied()
    }

    pub fn node_count(&self) -> usize {
        self.host_nodes.len()
    }

    pub fn edge_count(&self) -> usize {
        self.host_edges.len()
    }

    /// Materialize the embedded substructure as an owned `Graph` with node
    /// ids 0..node_count and edge ids 0..edge_count.
    pub fn extract(&self) -> Graph {
        let mut sub_edges = Vec::with_capacity(self.host_edges.len());
        for &eid in &self.host_edges {
            let [ha, hb] = self.graph.edge_endpoints(eid);
            let sa = self.sub_nodes[&ha];
            let sb = self.sub_nodes[&hb];
            sub_edges.push([sa.0, sb.0]);
        }
        Graph::new(self.host_nodes.len(), &sub_edges)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[test]
    fn test_graph_default() {
        let g = Graph::default();
        assert_eq!(g.node_count(), 0);
        assert_eq!(g.edge_count(), 0);
        assert!(g.is_dense());
    }

    #[rstest]
    #[case::empty(0, vec![], 0, 0)]
    #[case::isolated_nodes(3, vec![], 3, 0)]
    #[case::single_edge(2, vec![[0, 1]], 2, 1)]
    #[case::triangle(3, vec![[0, 1], [1, 2], [0, 2]], 3, 3)]
    fn test_graph_new(
        #[case] node_count: usize,
        #[case] edges: Vec<[u32; 2]>,
        #[case] expected_nodes: usize,
        #[case] expected_edges: usize,
    ) {
        let g = Graph::new(node_count, &edges);
        assert_eq!(g.node_count(), expected_nodes);
        assert_eq!(g.edge_count(), expected_edges);
        assert!(g.is_dense());
    }

    #[test]
    fn test_graph_neighbors() {
        let g = Graph::new(3, &[[0, 1], [0, 2]]);
        assert_eq!(g.degree(NodeId(0)), 2);
        assert_eq!(g.degree(NodeId(1)), 1);
        assert_eq!(g.degree(NodeId(2)), 1);
        assert_eq!(g.neighbors(NodeId(0)).len(), 2);
        assert_eq!(g.neighbors(NodeId(1))[0].node, NodeId(0));
        assert_eq!(g.neighbors(NodeId(2))[0].node, NodeId(0));
    }

    #[rstest]
    #[case::star_reverse(&[[0, 3], [0, 2], [0, 1]])]
    #[case::triangle(&[[2, 0], [1, 2], [0, 1]])]
    #[case::path_reverse(&[[3, 2], [2, 1], [1, 0]])]
    fn test_graph_neighbors_sorted(#[case] edges: &[[u32; 2]]) {
        let node_count = edges.iter().flat_map(|e| e.iter()).max().unwrap() + 1;
        let g = Graph::new(node_count as usize, edges);
        for nid in 0..g.node_count() {
            let nbrs = g.neighbors(NodeId(nid as u32));
            for w in nbrs.windows(2) {
                assert!(
                    w[0].node < w[1].node,
                    "neighbors of node {nid} not sorted: {:?}",
                    nbrs
                );
            }
        }
    }

    #[rstest]
    #[case::single_edge(&[[0, 1]], &[0, 1], vec![EdgeId(0)])]
    #[case::triangle_pair(&[[0, 1], [1, 2], [0, 2]], &[0, 1], vec![EdgeId(0)])]
    #[case::triangle_all(&[[0, 1], [1, 2], [0, 2]], &[0, 1, 2], vec![EdgeId(0), EdgeId(1), EdgeId(2)])]
    #[case::no_internal_edges(&[[0, 1], [1, 2]], &[0, 2], vec![])]
    #[case::empty_nodes(&[[0, 1]], &[], vec![])]
    fn test_graph_induced_edges(
        #[case] edges: &[[u32; 2]],
        #[case] nodes: &[u32],
        #[case] expected: Vec<EdgeId>,
    ) {
        let node_count = edges.iter().flat_map(|e| e.iter()).max().unwrap() + 1;
        let g = Graph::new(node_count as usize, edges);
        let sorted_nodes: Vec<NodeId> = nodes.iter().map(|&n| NodeId(n)).collect();
        let mut result: Vec<EdgeId> = g.induced_edges(&sorted_nodes).collect();
        result.sort_unstable();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_graph_edge_endpoints() {
        let g = Graph::new(3, &[[0, 1], [1, 2]]);
        assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
        assert_eq!(g.edge_endpoints(EdgeId(1)), [NodeId(1), NodeId(2)]);
    }

    #[test]
    fn test_graph_find_edge() {
        let g = Graph::new(3, &[[0, 1]]);
        let e = g.find_edge(NodeId(0), NodeId(1));
        assert_eq!(e, Some(EdgeId(0)));
        assert_eq!(g.find_edge(NodeId(1), NodeId(0)), Some(EdgeId(0)));
        assert_eq!(g.find_edge(NodeId(0), NodeId(2)), None);
    }

    #[test]
    fn test_graph_self_loop() {
        let g = Graph::new(1, &[[0, 0]]);
        assert_eq!(g.degree(NodeId(0)), 2);
        assert_eq!(g.edge_count(), 1);
    }

    #[test]
    fn test_graph_parallel_edges() {
        let g = Graph::new(2, &[[0, 1], [0, 1]]);
        assert_eq!(g.edge_count(), 2);
        assert_eq!(g.degree(NodeId(0)), 2);
    }

    #[test]
    fn test_graph_node_ids() {
        let g = Graph::new(3, &[]);
        let ids: Vec<NodeId> = g.node_ids().collect();
        assert_eq!(ids, vec![NodeId(0), NodeId(1), NodeId(2)]);
    }

    #[test]
    fn test_graph_edge_ids() {
        let g = Graph::new(3, &[[0, 1], [1, 2]]);
        let ids: Vec<EdgeId> = g.edge_ids().collect();
        assert_eq!(ids, vec![EdgeId(0), EdgeId(1)]);
    }

    #[test]
    fn test_graph_contains() {
        let g = Graph::new(2, &[[0, 1]]);
        assert!(g.contains_node(NodeId(0)));
        assert!(g.contains_node(NodeId(1)));
        assert!(!g.contains_node(NodeId(2)));
        assert!(g.contains_edge(EdgeId(0)));
        assert!(!g.contains_edge(EdgeId(1)));
    }

    #[test]
    fn test_graph_arc_sharing() {
        let g1 = Graph::new(3, &[[0, 1], [1, 2]]);
        let g2 = g1.clone();
        assert_eq!(g1, g2);
    }

    #[test]
    fn test_graph_add_node() {
        let mut g = Graph::new(2, &[[0, 1]]);
        let n = g.add_node();
        assert_eq!(n, NodeId(2));
        assert_eq!(g.node_count(), 3);
        assert_eq!(g.edge_count(), 1);
        assert_eq!(g.degree(NodeId(2)), 0);
        assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
    }

    #[test]
    fn test_graph_add_edge() {
        let mut g = Graph::new(3, &[[0, 1]]);
        let e = g.add_edge(NodeId(1), NodeId(2));
        assert_eq!(e, EdgeId(1));
        assert_eq!(g.edge_count(), 2);
        assert_eq!(g.degree(NodeId(1)), 2);
        assert_eq!(g.edge_endpoints(EdgeId(1)), [NodeId(1), NodeId(2)]);
    }

    #[test]
    fn test_graph_remove_node() {
        // 0--1--2, remove node 1
        let mut g = Graph::new(3, &[[0, 1], [1, 2]]);
        let remap = g.remove_node(NodeId(1));

        assert_eq!(g.node_count(), 2);
        assert_eq!(g.edge_count(), 0);
        assert_eq!(remap.removed_nodes, vec![1]);
        assert_eq!(remap.removed_edges, vec![0, 1]);

        // node 0 stays 0, node 2 becomes 1
        assert_eq!(remap.map_node(NodeId(0)), Some(NodeId(0)));
        assert_eq!(remap.map_node(NodeId(1)), None);
        assert_eq!(remap.map_node(NodeId(2)), Some(NodeId(1)));
    }

    #[test]
    fn test_graph_remove_node_partial() {
        // triangle 0-1, 1-2, 0-2; remove node 0
        let mut g = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let remap = g.remove_node(NodeId(0));

        assert_eq!(g.node_count(), 2);
        assert_eq!(g.edge_count(), 1);
        assert_eq!(remap.removed_nodes, vec![0]);
        assert_eq!(remap.removed_edges, vec![0, 2]);

        // surviving edge (old 1) maps to new 0
        assert_eq!(remap.map_edge(EdgeId(0)), None);
        assert_eq!(remap.map_edge(EdgeId(1)), Some(EdgeId(0)));
        assert_eq!(remap.map_edge(EdgeId(2)), None);

        // nodes 1,2 become 0,1
        assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
    }

    #[test]
    fn test_graph_remove_edge() {
        // triangle, remove edge 1 (1-2)
        let mut g = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let remap = g.remove_edge(EdgeId(1));

        assert_eq!(g.node_count(), 3);
        assert_eq!(g.edge_count(), 2);
        assert_eq!(remap.removed_nodes, Vec::<u32>::new());
        assert_eq!(remap.removed_edges, vec![1]);

        assert_eq!(remap.map_edge(EdgeId(0)), Some(EdgeId(0)));
        assert_eq!(remap.map_edge(EdgeId(1)), None);
        assert_eq!(remap.map_edge(EdgeId(2)), Some(EdgeId(1)));

        assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
        assert_eq!(g.edge_endpoints(EdgeId(1)), [NodeId(0), NodeId(2)]);
    }

    #[test]
    fn test_graph_cow_sharing() {
        let g1 = Graph::new(3, &[[0, 1], [1, 2]]);
        let mut g2 = g1.clone();
        g2.add_node();

        // g1 unchanged
        assert_eq!(g1.node_count(), 3);
        assert_eq!(g2.node_count(), 4);
    }

    #[test]
    fn test_graph_remove_batch() {
        // 0-1, 1-2, 2-3, 3-4; remove nodes 1 and 3
        let mut g = Graph::new(5, &[[0, 1], [1, 2], [2, 3], [3, 4]]);
        let remap = g.remove(&[NodeId(1), NodeId(3)], &[]);

        assert_eq!(g.node_count(), 3);
        // edges 0(0-1), 1(1-2), 2(2-3), 3(3-4) — all incident to 1 or 3 are removed
        // only none survive since every edge touches node 1 or 3
        assert_eq!(g.edge_count(), 0);

        assert_eq!(remap.map_node(NodeId(0)), Some(NodeId(0)));
        assert_eq!(remap.map_node(NodeId(1)), None);
        assert_eq!(remap.map_node(NodeId(2)), Some(NodeId(1)));
        assert_eq!(remap.map_node(NodeId(3)), None);
        assert_eq!(remap.map_node(NodeId(4)), Some(NodeId(2)));
    }

    #[test]
    fn test_graph_remove_nodes_and_edges() {
        // 0-1, 1-2, 2-3, 0-3; remove node 1, edge 3 (0-3)
        let mut g = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [0, 3]]);
        let remap = g.remove(&[NodeId(1)], &[EdgeId(3)]);

        assert_eq!(g.node_count(), 3);
        // edge 0(0-1) removed (incident to 1)
        // edge 1(1-2) removed (incident to 1)
        // edge 2(2-3) survives → becomes edge 0 (endpoints: 1-2 after shift)
        // edge 3(0-3) removed (explicitly)
        assert_eq!(g.edge_count(), 1);
        assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(1), NodeId(2)]);

        assert_eq!(remap.map_edge(EdgeId(2)), Some(EdgeId(0)));
    }

    #[rstest]
    #[case::identity(NodeId(0), vec![], Some(NodeId(0)))]
    #[case::before_removed(NodeId(0), vec![2], Some(NodeId(0)))]
    #[case::removed(NodeId(2), vec![2], None)]
    #[case::after_removed(NodeId(3), vec![2], Some(NodeId(2)))]
    #[case::multi_removed(NodeId(5), vec![1, 3], Some(NodeId(3)))]
    fn test_remapping_node(
        #[case] old: NodeId,
        #[case] removed: Vec<u32>,
        #[case] expected: Option<NodeId>,
    ) {
        let remap = Remapping {
            removed_nodes: removed,
            removed_edges: vec![],
        };
        assert_eq!(remap.map_node(old), expected);
    }

    #[rstest]
    #[case::identity(NodeId(0), vec![], NodeId(0))]
    #[case::before_gap(NodeId(0), vec![2], NodeId(0))]
    #[case::at_gap(NodeId(2), vec![2], NodeId(3))]
    #[case::after_gap(NodeId(3), vec![2], NodeId(4))]
    #[case::multi_removed(NodeId(3), vec![1, 3], NodeId(5))]
    fn test_remapping_unmap_node(
        #[case] post: NodeId,
        #[case] removed: Vec<u32>,
        #[case] expected: NodeId,
    ) {
        let remap = Remapping {
            removed_nodes: removed,
            removed_edges: vec![],
        };
        assert_eq!(remap.unmap_node(post), expected);
    }
}
