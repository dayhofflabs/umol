//! CSR graph with Arc-based structural sharing.
//!
//! `Graph` stores only adjacency (offsets, neighbor lists, edge endpoints).
//! Node and edge data live externally in `Vec`s indexed by `NodeId`/`EdgeId`.
//! The CSR is wrapped in `Arc` for zero-cost cloning; mutations rebuild
//! it and produce a `Compaction` for reindexing external data.

use std::collections::HashSet;
use std::sync::Arc;

use crate::correspondence::{Correspondence, GraphCorrespondence};

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

impl From<usize> for NodeId {
    fn from(index: usize) -> Self {
        Self(index as u32)
    }
}

impl EdgeId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<usize> for EdgeId {
    fn from(index: usize) -> Self {
        Self(index as u32)
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
/// `Arc`; mutations trigger copy-on-write. Each undirected edge endpoint pair
/// is stored in ascending `NodeId` order, while edge order and `EdgeId`
/// identity are preserved. Self-loops and parallel edges are retained.
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

    pub fn find_edge(&self, first: NodeId, second: NodeId) -> Option<EdgeId> {
        let nbrs = self.neighbors(first);
        nbrs.binary_search_by_key(&second, |n| n.node)
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

    /// Adjacency lists: for each node, its neighbor nodes, sorted by id and
    /// deduped (parallel edges collapse). Indexed by `NodeId::index`.
    pub fn adjacency(&self) -> Vec<Vec<NodeId>> {
        self.node_ids()
            .map(|node| {
                // `neighbors` is already sorted by `NodeId`, so dedup suffices.
                let mut nodes: Vec<NodeId> = self.neighbors(node).iter().map(|n| n.node).collect();
                nodes.dedup();
                nodes
            })
            .collect()
    }

    /// Line-graph adjacency: for each edge, the edges sharing one of its
    /// endpoints, sorted by id and deduped. Indexed by `EdgeId::index`.
    pub fn edge_adjacency(&self) -> Vec<Vec<EdgeId>> {
        let mut adjacency: Vec<Vec<EdgeId>> = vec![Vec::new(); self.edge_count()];
        for edge in self.edge_ids() {
            let [a, b] = self.edge_endpoints(edge);
            let mut neighbors: Vec<EdgeId> = Vec::new();
            for endpoint in [a, b] {
                for neighbor in self.neighbors(endpoint) {
                    if neighbor.edge != edge {
                        neighbors.push(neighbor.edge);
                    }
                }
            }
            neighbors.sort_unstable_by_key(|e| e.0);
            neighbors.dedup();
            adjacency[edge.index()] = neighbors;
        }
        adjacency
    }

    pub fn contains_node(&self, id: NodeId) -> bool {
        id.index() < self.csr.node_count
    }

    pub fn contains_edge(&self, id: EdgeId) -> bool {
        id.index() < self.csr.edge_count
    }

    pub fn node_ids(&self) -> impl ExactSizeIterator<Item = NodeId> {
        (0..self.csr.node_count as u32).map(NodeId)
    }

    pub fn edge_ids(&self) -> impl ExactSizeIterator<Item = EdgeId> {
        (0..self.csr.edge_count as u32).map(EdgeId)
    }

    pub fn is_dense(&self) -> bool {
        true
    }

    /// Returns whether this graph has neither self-loops nor parallel edges.
    pub fn is_simple(&self) -> bool {
        let mut endpoints = HashSet::with_capacity(self.edge_count());
        self.csr
            .endpoints
            .iter()
            .all(|&[first, second]| first != second && endpoints.insert([first, second]))
    }

    pub fn add_node(&mut self) -> NodeId {
        let old = &*self.csr;
        let new_id = NodeId(old.node_count as u32);
        let edges: Vec<[u32; 2]> = old.endpoints.iter().map(|&[a, b]| [a.0, b.0]).collect();
        self.csr = Arc::new(Self::build_csr(old.node_count + 1, &edges));
        new_id
    }

    pub fn add_edge(&mut self, first: NodeId, second: NodeId) -> EdgeId {
        let old = &*self.csr;
        let new_id = EdgeId(old.edge_count as u32);
        let mut edges: Vec<[u32; 2]> = old.endpoints.iter().map(|&[s, t]| [s.0, t.0]).collect();
        edges.push([first.0, second.0]);
        self.csr = Arc::new(Self::build_csr(old.node_count, &edges));
        new_id
    }

    /// SqPO-style removal: delete `nodes` and `edges`, sweeping along every edge incident to a
    /// removed node (deletion in unknown context). Always succeeds; incident edges the caller did
    /// not list are dropped too. Returns the [`Compaction`] renumbering. For the DPO discipline that
    /// rejects a stranded edge instead of sweeping it, use [`Graph::try_remove`].
    pub fn remove_cascading(&mut self, nodes: &[NodeId], edges: &[EdgeId]) -> Compaction {
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

        Compaction::new(removed_nodes, removed_edge_set)
    }

    /// DPO-style removal: delete exactly `nodes` and `edges`, or `None` when that would strand an
    /// edge — the **dangling condition**: a removed node incident to an edge the caller did not list
    /// for removal. On success the result equals [`Self::remove_cascading`] (the check guarantees
    /// there is nothing extra to sweep), so it is the pushout-complement of a matched deletion.
    pub fn try_remove(&mut self, nodes: &[NodeId], edges: &[EdgeId]) -> Option<Compaction> {
        let removed_edges: HashSet<EdgeId> = edges.iter().copied().collect();
        for &node in nodes {
            if self
                .neighbors(node)
                .iter()
                .any(|n| !removed_edges.contains(&n.edge))
            {
                return None;
            }
        }
        Some(self.remove_cascading(nodes, edges))
    }

    pub fn remove_node_cascading(&mut self, id: NodeId) -> Compaction {
        self.remove_cascading(&[id], &[])
    }

    pub fn remove_edge_cascading(&mut self, id: EdgeId) -> Compaction {
        self.remove_cascading(&[], &[id])
    }

    /// Build an induced subgraph from a subset of nodes.
    ///
    /// Induced subgraph over `nodes` (deduplicated, first occurrence wins), as an injective sub→host
    /// [`GraphCorrespondence`]: sub node `i` is `nodes[i]`, edges are *all* host edges among the node
    /// set. [`Graph::extract`] materializes it.
    pub fn induced_subgraph(&self, nodes: &[NodeId]) -> GraphCorrespondence {
        let mut host_nodes: Vec<NodeId> = Vec::with_capacity(nodes.len());
        let mut node_set: HashSet<NodeId> = HashSet::with_capacity(nodes.len());
        for &node in nodes {
            if node_set.insert(node) {
                host_nodes.push(node);
            }
        }
        let host_edges: Vec<EdgeId> = self
            .edge_ids()
            .filter(|&eid| {
                let [a, b] = self.edge_endpoints(eid);
                node_set.contains(&a) && node_set.contains(&b)
            })
            .collect();
        GraphCorrespondence::new(
            Correspondence::from_images(&host_nodes, self.node_count()),
            Correspondence::from_images(&host_edges, self.edge_count()),
        )
    }

    /// Subgraph induced by an edge subset: nodes are the endpoints of `edges`,
    /// edges are exactly `edges` (deduped, first occurrence kept). Unlike
    /// [`Graph::induced_subgraph`], chords among the endpoints are excluded, so a
    /// path and a chorded ring on the same atoms stay distinct.
    pub fn edge_induced_subgraph(&self, edges: &[EdgeId]) -> GraphCorrespondence {
        let mut host_nodes: Vec<NodeId> = Vec::new();
        let mut node_set: HashSet<NodeId> = HashSet::new();
        let mut host_edges: Vec<EdgeId> = Vec::with_capacity(edges.len());
        let mut seen_edges: HashSet<EdgeId> = HashSet::with_capacity(edges.len());
        for &eid in edges {
            if !seen_edges.insert(eid) {
                continue;
            }
            host_edges.push(eid);
            for node in self.edge_endpoints(eid) {
                if node_set.insert(node) {
                    host_nodes.push(node);
                }
            }
        }
        GraphCorrespondence::new(
            Correspondence::from_images(&host_nodes, self.node_count()),
            Correspondence::from_images(&host_edges, self.edge_count()),
        )
    }

    /// Materialize the subgraph `sub` (over `self` as host) as an owned `Graph` with node ids
    /// `0..node_count` and edge ids `0..edge_count`; edge endpoints are remapped host→sub through the
    /// node correspondence.
    pub fn extract(&self, sub: &GraphCorrespondence) -> Graph {
        let sub_edges: Vec<[u32; 2]> = sub
            .edges()
            .matched_pairs()
            .iter()
            .map(|&(_, host_edge)| {
                let [ha, hb] = self.edge_endpoints(host_edge);
                let sa = sub
                    .nodes()
                    .left_of(ha)
                    .expect("subgraph edge endpoint is a subgraph node");
                let sb = sub
                    .nodes()
                    .left_of(hb)
                    .expect("subgraph edge endpoint is a subgraph node");
                [sa.0, sb.0]
            })
            .collect();
        Graph::new(sub.nodes().matched_pair_count(), &sub_edges)
    }

    /// Subdivide every edge exactly once.
    ///
    /// Source nodes retain their ids. Each source edge receives one inserted
    /// node, and its two incident subdivision edges are consecutive.
    pub fn subdivide_edges(&self) -> SubdividedGraph {
        let source_node_count = self.node_count();
        let source_edge_count = self.edge_count();
        let edges: Vec<[u32; 2]> = self
            .edge_ids()
            .flat_map(|edge| {
                let [first, second] = self.edge_endpoints(edge);
                let inserted = NodeId((source_node_count + edge.index()) as u32);
                [[first.0, inserted.0], [inserted.0, second.0]]
            })
            .collect();

        SubdividedGraph {
            graph: Graph::new(source_node_count + source_edge_count, &edges),
            source_node_count,
            source_edge_count,
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

/// The source represented by a node in a [`SubdividedGraph`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SubdivisionNodeSource {
    Node(NodeId),
    Edge(EdgeId),
}

/// A graph in which every source edge has been replaced by a two-edge path.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubdividedGraph {
    graph: Graph,
    source_node_count: usize,
    source_edge_count: usize,
}

impl SubdividedGraph {
    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    /// Return the source node or edge represented by a subdivision node.
    ///
    /// # Panics
    ///
    /// Panics when `node` is not in the subdivided graph.
    pub fn node_source(&self, node: NodeId) -> SubdivisionNodeSource {
        assert!(
            self.graph.contains_node(node),
            "subdivision node {} out of range",
            node.0
        );
        if node.index() < self.source_node_count {
            SubdivisionNodeSource::Node(node)
        } else {
            SubdivisionNodeSource::Edge(EdgeId((node.index() - self.source_node_count) as u32))
        }
    }

    /// Return the subdivision node representing a source node or edge.
    ///
    /// # Panics
    ///
    /// Panics when the source id is not in the source graph.
    pub fn node_of(&self, source: SubdivisionNodeSource) -> NodeId {
        match source {
            SubdivisionNodeSource::Node(node) => {
                assert!(
                    node.index() < self.source_node_count,
                    "source node {} out of range",
                    node.0
                );
                node
            }
            SubdivisionNodeSource::Edge(edge) => {
                assert!(
                    edge.index() < self.source_edge_count,
                    "source edge {} out of range",
                    edge.0
                );
                NodeId((self.source_node_count + edge.index()) as u32)
            }
        }
    }

    /// Return the source edge represented by a subdivision incidence edge.
    ///
    /// # Panics
    ///
    /// Panics when `incidence` is not in the subdivided graph.
    pub fn edge_source(&self, incidence: EdgeId) -> EdgeId {
        assert!(
            self.graph.contains_edge(incidence),
            "subdivision edge {} out of range",
            incidence.0
        );
        EdgeId((incidence.index() / 2) as u32)
    }

    /// Return the two subdivision incidence edges created from a source edge.
    ///
    /// # Panics
    ///
    /// Panics when `edge` is not in the source graph.
    pub fn incidence_edges_of(&self, edge: EdgeId) -> [EdgeId; 2] {
        assert!(
            edge.index() < self.source_edge_count,
            "source edge {} out of range",
            edge.0
        );
        [
            EdgeId((2 * edge.index()) as u32),
            EdgeId((2 * edge.index() + 1) as u32),
        ]
    }
}

impl Default for Graph {
    fn default() -> Self {
        Self::new(0, &[])
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Compaction {
    removed_nodes: Vec<u32>,
    removed_edges: Vec<u32>,
}

impl Compaction {
    pub fn new(mut removed_nodes: Vec<u32>, mut removed_edges: Vec<u32>) -> Self {
        removed_nodes.sort_unstable();
        removed_nodes.dedup();
        removed_edges.sort_unstable();
        removed_edges.dedup();
        Self {
            removed_nodes,
            removed_edges,
        }
    }

    pub fn compact_node(&self, old: NodeId) -> Option<NodeId> {
        if self.removed_nodes.binary_search(&old.0).is_ok() {
            return None;
        }
        let shift = self.removed_nodes.partition_point(|&r| r < old.0);
        Some(NodeId(old.0 - shift as u32))
    }

    pub fn compact_edge(&self, old: EdgeId) -> Option<EdgeId> {
        if self.removed_edges.binary_search(&old.0).is_ok() {
            return None;
        }
        let shift = self.removed_edges.partition_point(|&r| r < old.0);
        Some(EdgeId(old.0 - shift as u32))
    }

    pub fn uncompact_node(&self, post: NodeId) -> NodeId {
        NodeId(uncompact_dense(&self.removed_nodes, post.0))
    }

    pub fn uncompact_edge(&self, post: EdgeId) -> EdgeId {
        EdgeId(uncompact_dense(&self.removed_edges, post.0))
    }
}

// Inverse dense shift: re-add removed ids at or below the post index (fixpoint).
fn uncompact_dense(removed: &[u32], post: u32) -> u32 {
    let mut old = post;
    loop {
        let next = post + removed.partition_point(|&r| r <= old) as u32;
        if next == old {
            return old;
        }
        old = next;
    }
}

/// Compact a node-indexed data column to the post-removal layout (drop removed, keep order).
pub fn compact_node_vec<T: Clone>(compaction: &Compaction, data: &[T]) -> Vec<T> {
    data.iter()
        .enumerate()
        .filter(|(i, _)| compaction.compact_node(NodeId(*i as u32)).is_some())
        .map(|(_, v)| v.clone())
        .collect()
}

/// Compact an edge-indexed data column to the post-removal layout (drop removed, keep order).
pub fn compact_edge_vec<T: Clone>(compaction: &Compaction, data: &[T]) -> Vec<T> {
    data.iter()
        .enumerate()
        .filter(|(i, _)| compaction.compact_edge(EdgeId(*i as u32)).is_some())
        .map(|(_, v)| v.clone())
        .collect()
}

/// General relabeling of node/edge ids: a **total** map old→new (no drops —
/// removal is `Compaction`). Indexed by old id, so `map_node(NodeId(i))`
/// is `nodes[i]`. The map may be an injection into a larger id space (e.g. a
/// composition's merged frame), so it is not necessarily a bijection.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Remapping {
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
}

impl Remapping {
    pub fn new(nodes: Vec<NodeId>, edges: Vec<EdgeId>) -> Self {
        Self { nodes, edges }
    }

    /// Return the image of `old`, or `None` when it lies outside the node source range.
    pub fn try_map_node(&self, old: NodeId) -> Option<NodeId> {
        self.nodes.get(old.0 as usize).copied()
    }

    /// Return the image of `old`, or `None` when it lies outside the edge source range.
    pub fn try_map_edge(&self, old: EdgeId) -> Option<EdgeId> {
        self.edges.get(old.0 as usize).copied()
    }

    /// Return the image of `old`.
    ///
    /// # Panics
    ///
    /// Panics when `old` lies outside the node source range defined at construction.
    pub fn map_node(&self, old: NodeId) -> NodeId {
        self.try_map_node(old)
            .expect("node id outside remapping source range")
    }

    /// Return the image of `old`.
    ///
    /// # Panics
    ///
    /// Panics when `old` lies outside the edge source range defined at construction.
    pub fn map_edge(&self, old: EdgeId) -> EdgeId {
        self.try_map_edge(old)
            .expect("edge id outside remapping source range")
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn assert_exact_size<T>(mut iterator: impl ExactSizeIterator<Item = T>, expected: Vec<T>)
    where
        T: Debug + PartialEq,
    {
        assert_eq!(iterator.len(), expected.len());
        assert_eq!(iterator.size_hint(), (expected.len(), Some(expected.len())));
        while let Some(expected_item) = expected.get(expected.len() - iterator.len()) {
            let previous = iterator.len();
            assert_eq!(iterator.next().as_ref(), Some(expected_item));
            let remaining = iterator.len();
            assert_eq!(remaining, previous - 1);
            assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
        }
        assert_eq!(iterator.next(), None);
        assert_eq!(iterator.len(), 0);
        assert_eq!(iterator.size_hint(), (0, Some(0)));
    }

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

    #[rstest]
    #[case::empty(0, &[], true)]
    #[case::isolated(3, &[], true)]
    #[case::simple(3, &[[0, 1], [1, 2]], true)]
    #[case::looped(1, &[[0, 0]], false)]
    #[case::parallel(2, &[[0, 1], [0, 1]], false)]
    #[case::reverse_parallel(2, &[[0, 1], [1, 0]], false)]
    #[case::mixed(4, &[[0, 1], [1, 2], [2, 2], [1, 3], [3, 1]], false)]
    fn test_graph_is_simple(
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
        #[case] expected: bool,
    ) {
        assert_eq!(Graph::new(node_count, edges).is_simple(), expected);
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

    #[rstest]
    #[case::path(&[[0, 1], [1, 2]], 3, vec![vec![NodeId(1)], vec![NodeId(0), NodeId(2)], vec![NodeId(1)]])]
    #[case::triangle(&[[0, 1], [1, 2], [0, 2]], 3, vec![vec![NodeId(1), NodeId(2)], vec![NodeId(0), NodeId(2)], vec![NodeId(0), NodeId(1)]])]
    #[case::isolated(&[[0, 1]], 3, vec![vec![NodeId(1)], vec![NodeId(0)], vec![]])]
    fn test_graph_adjacency(
        #[case] edges: &[[u32; 2]],
        #[case] node_count: usize,
        #[case] expected: Vec<Vec<NodeId>>,
    ) {
        let g = Graph::new(node_count, edges);
        assert_eq!(g.adjacency(), expected);
    }

    #[rstest]
    #[case::path(&[[0, 1], [1, 2]], vec![vec![EdgeId(1)], vec![EdgeId(0)]])]
    #[case::triangle(&[[0, 1], [1, 2], [0, 2]], vec![vec![EdgeId(1), EdgeId(2)], vec![EdgeId(0), EdgeId(2)], vec![EdgeId(0), EdgeId(1)]])]
    #[case::disjoint(&[[0, 1], [2, 3]], vec![vec![], vec![]])]
    fn test_graph_edge_adjacency(#[case] edges: &[[u32; 2]], #[case] expected: Vec<Vec<EdgeId>>) {
        let node_count = edges.iter().flat_map(|e| e.iter()).max().unwrap() + 1;
        let g = Graph::new(node_count as usize, edges);
        assert_eq!(g.edge_adjacency(), expected);
    }

    #[rstest]
    #[case::chord_excluded(&[[0, 1], [1, 2], [0, 2]], &[0, 1], vec![NodeId(0), NodeId(1), NodeId(2)], vec![EdgeId(0), EdgeId(1)])]
    #[case::full_ring(&[[0, 1], [1, 2], [0, 2]], &[0, 1, 2], vec![NodeId(0), NodeId(1), NodeId(2)], vec![EdgeId(0), EdgeId(1), EdgeId(2)])]
    #[case::deduped(&[[0, 1], [1, 2]], &[0, 0, 1], vec![NodeId(0), NodeId(1), NodeId(2)], vec![EdgeId(0), EdgeId(1)])]
    fn test_graph_edge_induced_subgraph(
        #[case] edges: &[[u32; 2]],
        #[case] subset: &[u32],
        #[case] expected_nodes: Vec<NodeId>,
        #[case] expected_edges: Vec<EdgeId>,
    ) {
        let node_count = edges.iter().flat_map(|e| e.iter()).max().unwrap() + 1;
        let g = Graph::new(node_count as usize, edges);
        let subset_edges: Vec<EdgeId> = subset.iter().map(|&e| EdgeId(e)).collect();
        let sub = g.edge_induced_subgraph(&subset_edges);
        let host_edges: Vec<EdgeId> = sub
            .edges()
            .matched_pairs()
            .iter()
            .map(|&(_, h)| h)
            .collect();
        assert_eq!(host_edges, expected_edges);
        let mut nodes: Vec<NodeId> = sub
            .nodes()
            .matched_pairs()
            .iter()
            .map(|&(_, h)| h)
            .collect();
        nodes.sort_unstable();
        assert_eq!(nodes, expected_nodes);
    }

    #[rstest]
    fn test_graph_induced_subgraph() {
        // triangle 0-1-2; induce {0, 1}: the chord to 2 is dropped, only edge 0 (0-1) survives.
        let g = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let sub = g.induced_subgraph(&[NodeId(0), NodeId(1)]);
        assert_eq!(
            sub.nodes().matched_pairs(),
            &[(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))]
        );
        assert_eq!(sub.edges().matched_pairs(), &[(EdgeId(0), EdgeId(0))]);
    }

    #[rstest]
    fn test_graph_extract() {
        // extract {0,1,2} of path 0-1-2-3 → path 0-1-2 with endpoints remapped host→sub.
        let g = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        let sub = g.induced_subgraph(&[NodeId(0), NodeId(1), NodeId(2)]);
        let extracted = g.extract(&sub);
        assert_eq!(extracted.node_count(), 3);
        assert_eq!(extracted.edge_count(), 2);
        assert_eq!(extracted.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
        assert_eq!(extracted.edge_endpoints(EdgeId(1)), [NodeId(1), NodeId(2)]);
    }

    #[rstest]
    #[case::empty(0, &[], 0, &[])]
    #[case::isolated(3, &[], 3, &[])]
    #[case::path(
        3,
        &[[0, 1], [1, 2]],
        5,
        &[[0, 3], [3, 1], [1, 4], [4, 2]],
    )]
    #[case::cycle(
        3,
        &[[0, 1], [1, 2], [2, 0]],
        6,
        &[[0, 3], [3, 1], [1, 4], [4, 2], [0, 5], [5, 2]],
    )]
    #[case::parallel(
        2,
        &[[0, 1], [0, 1]],
        4,
        &[[0, 2], [2, 1], [0, 3], [3, 1]],
    )]
    #[case::self_loop(1, &[[0, 0]], 2, &[[0, 1], [1, 0]])]
    fn test_graph_subdivide_edges(
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
        #[case] expected_node_count: usize,
        #[case] expected_edges: &[[u32; 2]],
    ) {
        let subdivision = Graph::new(node_count, edges).subdivide_edges();
        assert_eq!(
            subdivision.graph(),
            &Graph::new(expected_node_count, expected_edges)
        );
    }

    #[rstest]
    #[case::first_node(NodeId(0), SubdivisionNodeSource::Node(NodeId(0)))]
    #[case::last_node(NodeId(2), SubdivisionNodeSource::Node(NodeId(2)))]
    #[case::first_edge(NodeId(3), SubdivisionNodeSource::Edge(EdgeId(0)))]
    #[case::last_edge(NodeId(4), SubdivisionNodeSource::Edge(EdgeId(1)))]
    fn test_subdivided_graph_node_source(
        #[case] node: NodeId,
        #[case] expected: SubdivisionNodeSource,
    ) {
        let subdivision = Graph::new(3, &[[0, 1], [1, 2]]).subdivide_edges();
        assert_eq!(subdivision.node_source(node), expected);
    }

    #[rstest]
    #[case::first_node(SubdivisionNodeSource::Node(NodeId(0)), NodeId(0))]
    #[case::last_node(SubdivisionNodeSource::Node(NodeId(2)), NodeId(2))]
    #[case::first_edge(SubdivisionNodeSource::Edge(EdgeId(0)), NodeId(3))]
    #[case::last_edge(SubdivisionNodeSource::Edge(EdgeId(1)), NodeId(4))]
    fn test_subdivided_graph_node_of(
        #[case] source: SubdivisionNodeSource,
        #[case] expected: NodeId,
    ) {
        let subdivision = Graph::new(3, &[[0, 1], [1, 2]]).subdivide_edges();
        assert_eq!(subdivision.node_of(source), expected);
    }

    #[rstest]
    #[case::first_first(EdgeId(0), EdgeId(0))]
    #[case::first_second(EdgeId(1), EdgeId(0))]
    #[case::second_first(EdgeId(2), EdgeId(1))]
    #[case::second_second(EdgeId(3), EdgeId(1))]
    fn test_subdivided_graph_edge_source(#[case] incidence: EdgeId, #[case] expected: EdgeId) {
        let subdivision = Graph::new(3, &[[0, 1], [1, 2]]).subdivide_edges();
        assert_eq!(subdivision.edge_source(incidence), expected);
    }

    #[rstest]
    #[case::first(EdgeId(0), [EdgeId(0), EdgeId(1)])]
    #[case::second(EdgeId(1), [EdgeId(2), EdgeId(3)])]
    fn test_subdivided_graph_incidence_edges_of(
        #[case] edge: EdgeId,
        #[case] expected: [EdgeId; 2],
    ) {
        let subdivision = Graph::new(3, &[[0, 1], [1, 2]]).subdivide_edges();
        assert_eq!(subdivision.incidence_edges_of(edge), expected);
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

    #[rstest]
    #[case::empty(Graph::new(0, &[]), vec![])]
    #[case::populated(
        Graph::new(3, &[]),
        vec![NodeId(0), NodeId(1), NodeId(2)],
    )]
    fn test_graph_node_ids(#[case] graph: Graph, #[case] expected: Vec<NodeId>) {
        assert_exact_size(graph.node_ids(), expected);
    }

    #[rstest]
    #[case::empty(Graph::new(0, &[]), vec![])]
    #[case::populated(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![EdgeId(0), EdgeId(1)],
    )]
    fn test_graph_edge_ids(#[case] graph: Graph, #[case] expected: Vec<EdgeId>) {
        assert_exact_size(graph.edge_ids(), expected);
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
    fn test_graph_remove_node_cascading() {
        // 0--1--2, remove node 1
        let mut g = Graph::new(3, &[[0, 1], [1, 2]]);
        let compaction = g.remove_node_cascading(NodeId(1));

        assert_eq!(g.node_count(), 2);
        assert_eq!(g.edge_count(), 0);
        assert_eq!(compaction.removed_nodes, vec![1]);
        assert_eq!(compaction.removed_edges, vec![0, 1]);

        // node 0 stays 0, node 2 becomes 1
        assert_eq!(compaction.compact_node(NodeId(0)), Some(NodeId(0)));
        assert_eq!(compaction.compact_node(NodeId(1)), None);
        assert_eq!(compaction.compact_node(NodeId(2)), Some(NodeId(1)));
    }

    #[test]
    fn test_graph_remove_node_cascading_partial() {
        // triangle 0-1, 1-2, 0-2; remove node 0
        let mut g = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let compaction = g.remove_node_cascading(NodeId(0));

        assert_eq!(g.node_count(), 2);
        assert_eq!(g.edge_count(), 1);
        assert_eq!(compaction.removed_nodes, vec![0]);
        assert_eq!(compaction.removed_edges, vec![0, 2]);

        // surviving edge (old 1) maps to new 0
        assert_eq!(compaction.compact_edge(EdgeId(0)), None);
        assert_eq!(compaction.compact_edge(EdgeId(1)), Some(EdgeId(0)));
        assert_eq!(compaction.compact_edge(EdgeId(2)), None);

        // nodes 1,2 become 0,1
        assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
    }

    #[test]
    fn test_graph_remove_edge_cascading() {
        // triangle, remove edge 1 (1-2)
        let mut g = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let compaction = g.remove_edge_cascading(EdgeId(1));

        assert_eq!(g.node_count(), 3);
        assert_eq!(g.edge_count(), 2);
        assert_eq!(compaction.removed_nodes, Vec::<u32>::new());
        assert_eq!(compaction.removed_edges, vec![1]);

        assert_eq!(compaction.compact_edge(EdgeId(0)), Some(EdgeId(0)));
        assert_eq!(compaction.compact_edge(EdgeId(1)), None);
        assert_eq!(compaction.compact_edge(EdgeId(2)), Some(EdgeId(1)));

        assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
        assert_eq!(g.edge_endpoints(EdgeId(1)), [NodeId(0), NodeId(2)]);
    }

    #[fixture]
    fn remapping() -> Remapping {
        Remapping::new(
            vec![NodeId(2), NodeId(0), NodeId(5)],
            vec![EdgeId(3), EdgeId(1)],
        )
    }

    #[rstest]
    #[case::first(NodeId(0), Some(NodeId(2)))]
    #[case::last(NodeId(2), Some(NodeId(5)))]
    #[case::uncovered(NodeId(3), None)]
    fn test_remapping_try_map_node(
        remapping: Remapping,
        #[case] old: NodeId,
        #[case] expected: Option<NodeId>,
    ) {
        assert_eq!(remapping.try_map_node(old), expected);
    }

    #[rstest]
    #[case::first(EdgeId(0), Some(EdgeId(3)))]
    #[case::last(EdgeId(1), Some(EdgeId(1)))]
    #[case::uncovered(EdgeId(2), None)]
    fn test_remapping_try_map_edge(
        remapping: Remapping,
        #[case] old: EdgeId,
        #[case] expected: Option<EdgeId>,
    ) {
        assert_eq!(remapping.try_map_edge(old), expected);
    }

    #[rstest]
    #[case::first(NodeId(0), NodeId(2))]
    #[case::middle(NodeId(1), NodeId(0))]
    #[case::last(NodeId(2), NodeId(5))]
    fn test_remapping_map_node(
        remapping: Remapping,
        #[case] old: NodeId,
        #[case] expected: NodeId,
    ) {
        assert_eq!(remapping.map_node(old), expected);
    }

    #[rstest]
    #[should_panic(expected = "node id outside remapping source range")]
    fn test_remapping_map_node_error(remapping: Remapping) {
        remapping.map_node(NodeId(3));
    }

    #[rstest]
    #[case::relabel(EdgeId(0), EdgeId(3))]
    #[case::fixed(EdgeId(1), EdgeId(1))]
    fn test_remapping_map_edge(
        remapping: Remapping,
        #[case] old: EdgeId,
        #[case] expected: EdgeId,
    ) {
        assert_eq!(remapping.map_edge(old), expected);
    }

    #[rstest]
    #[should_panic(expected = "edge id outside remapping source range")]
    fn test_remapping_map_edge_error(remapping: Remapping) {
        remapping.map_edge(EdgeId(2));
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
    fn test_graph_remove_cascading_batch() {
        // 0-1, 1-2, 2-3, 3-4; remove nodes 1 and 3
        let mut g = Graph::new(5, &[[0, 1], [1, 2], [2, 3], [3, 4]]);
        let compaction = g.remove_cascading(&[NodeId(1), NodeId(3)], &[]);

        assert_eq!(g.node_count(), 3);
        // edges 0(0-1), 1(1-2), 2(2-3), 3(3-4) — all incident to 1 or 3 are removed
        // only none survive since every edge touches node 1 or 3
        assert_eq!(g.edge_count(), 0);

        assert_eq!(compaction.compact_node(NodeId(0)), Some(NodeId(0)));
        assert_eq!(compaction.compact_node(NodeId(1)), None);
        assert_eq!(compaction.compact_node(NodeId(2)), Some(NodeId(1)));
        assert_eq!(compaction.compact_node(NodeId(3)), None);
        assert_eq!(compaction.compact_node(NodeId(4)), Some(NodeId(2)));
    }

    #[test]
    fn test_graph_remove_cascading_nodes_and_edges() {
        // 0-1, 1-2, 2-3, 0-3; remove node 1, edge 3 (0-3)
        let mut g = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [0, 3]]);
        let compaction = g.remove_cascading(&[NodeId(1)], &[EdgeId(3)]);

        assert_eq!(g.node_count(), 3);
        // edge 0(0-1) removed (incident to 1)
        // edge 1(1-2) removed (incident to 1)
        // edge 2(2-3) survives → becomes edge 0 (endpoints: 1-2 after shift)
        // edge 3(0-3) removed (explicitly)
        assert_eq!(g.edge_count(), 1);
        assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(1), NodeId(2)]);

        assert_eq!(compaction.compact_edge(EdgeId(2)), Some(EdgeId(0)));
    }

    #[rstest]
    fn test_graph_try_remove_clean() {
        // path 0-1-2-3; remove node 1 together with both its incident edges (0-1, 1-2).
        let mut g = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        let compaction = g
            .try_remove(&[NodeId(1)], &[EdgeId(0), EdgeId(1)])
            .expect("no dangling: every edge on node 1 is listed");
        assert_eq!(g.node_count(), 3);
        assert_eq!(g.edge_count(), 1);
        assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(1), NodeId(2)]);
        assert_eq!(compaction.compact_node(NodeId(2)), Some(NodeId(1)));
    }

    #[rstest]
    fn test_graph_try_remove_dangling() {
        // path 0-1-2-3; removing node 1 without its incident edges strands them → rejected.
        let mut g = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        assert_eq!(g.try_remove(&[NodeId(1)], &[]), None);
        // graph is left untouched on rejection.
        assert_eq!(g.node_count(), 4);
        assert_eq!(g.edge_count(), 3);
    }

    #[rstest]
    #[case::identity(NodeId(0), vec![], Some(NodeId(0)))]
    #[case::before_removed(NodeId(0), vec![2], Some(NodeId(0)))]
    #[case::removed(NodeId(2), vec![2], None)]
    #[case::after_removed(NodeId(3), vec![2], Some(NodeId(2)))]
    #[case::multi_removed(NodeId(5), vec![1, 3], Some(NodeId(3)))]
    fn test_compaction_node(
        #[case] old: NodeId,
        #[case] removed: Vec<u32>,
        #[case] expected: Option<NodeId>,
    ) {
        let compaction = Compaction::new(removed, vec![]);
        assert_eq!(compaction.compact_node(old), expected);
    }

    #[rstest]
    #[case::identity(NodeId(0), vec![], NodeId(0))]
    #[case::before_gap(NodeId(0), vec![2], NodeId(0))]
    #[case::at_gap(NodeId(2), vec![2], NodeId(3))]
    #[case::after_gap(NodeId(3), vec![2], NodeId(4))]
    #[case::multi_removed(NodeId(3), vec![1, 3], NodeId(5))]
    fn test_uncompaction_node(
        #[case] post: NodeId,
        #[case] removed: Vec<u32>,
        #[case] expected: NodeId,
    ) {
        let compaction = Compaction::new(removed, vec![]);
        assert_eq!(compaction.uncompact_node(post), expected);
    }
}
