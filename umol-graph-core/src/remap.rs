//! Total relabeling of an id space.
//!
//! A remapping gives every source id an image and drops nothing, so it expresses alpha-renaming and
//! injection into a larger ambient space. Removal is [`crate::compact::Compaction`].

use index_vec::{Idx, IndexVec};

use crate::graph::{EdgeId, NodeId};

/// Total relabeling of one dense source id space.
///
/// The image vector defines the source domain: source id `i` maps to `images[i]`. Every source id
/// therefore has an image. Images may repeat or occupy only part of a larger target id space;
/// consumers that require an injection or a dense target establish that separately.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Remapping<Id: Idx> {
    images: IndexVec<Id, Id>,
}

impl<Id: Idx> Remapping<Id> {
    /// Construct a remapping whose source domain is `0..images.len()`.
    pub fn new(images: Vec<Id>) -> Self {
        Self {
            images: IndexVec::from_vec(images),
        }
    }

    /// Return the image of `old`, or `None` when it lies outside the source domain.
    pub fn try_map(&self, old: Id) -> Option<Id> {
        self.images.get(old).copied()
    }

    /// Return the image of `old`.
    ///
    /// # Panics
    ///
    /// Panics when `old` lies outside the source domain defined at construction.
    pub fn map(&self, old: Id) -> Id {
        self.try_map(old)
            .expect("id outside remapping source domain")
    }
}

/// General relabeling of node/edge ids: a **total** map old→new (no drops —
/// removal is `Compaction`). Indexed by old id, so `map_node(NodeId(i))`
/// is `nodes[i]`. The map may be an injection into a larger id space (e.g. a
/// composition's merged frame), so it is not necessarily a bijection.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GraphRemapping {
    nodes: Remapping<NodeId>,
    edges: Remapping<EdgeId>,
}

impl GraphRemapping {
    pub fn new(nodes: Vec<NodeId>, edges: Vec<EdgeId>) -> Self {
        Self {
            nodes: Remapping::new(nodes),
            edges: Remapping::new(edges),
        }
    }

    pub fn nodes(&self) -> &Remapping<NodeId> {
        &self.nodes
    }

    pub fn edges(&self) -> &Remapping<EdgeId> {
        &self.edges
    }

    /// Return the image of `old`, or `None` when it lies outside the node source range.
    pub fn try_map_node(&self, old: NodeId) -> Option<NodeId> {
        self.nodes.try_map(old)
    }

    /// Return the image of `old`, or `None` when it lies outside the edge source range.
    pub fn try_map_edge(&self, old: EdgeId) -> Option<EdgeId> {
        self.edges.try_map(old)
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
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::first(vec![NodeId(5), NodeId(1)], NodeId(0), Some(NodeId(5)))]
    #[case::last(vec![NodeId(5), NodeId(1)], NodeId(1), Some(NodeId(1)))]
    #[case::sparse_target(vec![NodeId(8)], NodeId(0), Some(NodeId(8)))]
    #[case::repeated_target(vec![NodeId(3), NodeId(3)], NodeId(1), Some(NodeId(3)))]
    #[case::empty(vec![], NodeId(0), None)]
    #[case::uncovered(vec![NodeId(5)], NodeId(1), None)]
    fn test_remapping_try_map(
        #[case] images: Vec<NodeId>,
        #[case] old: NodeId,
        #[case] expected: Option<NodeId>,
    ) {
        assert_eq!(Remapping::new(images).try_map(old), expected);
    }

    #[rstest]
    #[case::first(vec![NodeId(5), NodeId(1)], NodeId(0), NodeId(5))]
    #[case::last(vec![NodeId(5), NodeId(1)], NodeId(1), NodeId(1))]
    fn test_remapping_map(
        #[case] images: Vec<NodeId>,
        #[case] old: NodeId,
        #[case] expected: NodeId,
    ) {
        assert_eq!(Remapping::new(images).map(old), expected);
    }

    #[rstest]
    #[should_panic(expected = "id outside remapping source domain")]
    fn test_remapping_map_error() {
        Remapping::<NodeId>::new(vec![]).map(NodeId(0));
    }

    #[fixture]
    fn graph_remapping() -> GraphRemapping {
        GraphRemapping::new(
            vec![NodeId(2), NodeId(0), NodeId(5)],
            vec![EdgeId(3), EdgeId(1)],
        )
    }

    #[rstest]
    #[case::first(NodeId(0), Some(NodeId(2)))]
    #[case::last(NodeId(2), Some(NodeId(5)))]
    #[case::uncovered(NodeId(3), None)]
    fn test_graph_remapping_try_map_node(
        graph_remapping: GraphRemapping,
        #[case] old: NodeId,
        #[case] expected: Option<NodeId>,
    ) {
        assert_eq!(graph_remapping.try_map_node(old), expected);
    }

    #[rstest]
    #[case::first(EdgeId(0), Some(EdgeId(3)))]
    #[case::last(EdgeId(1), Some(EdgeId(1)))]
    #[case::uncovered(EdgeId(2), None)]
    fn test_graph_remapping_try_map_edge(
        graph_remapping: GraphRemapping,
        #[case] old: EdgeId,
        #[case] expected: Option<EdgeId>,
    ) {
        assert_eq!(graph_remapping.try_map_edge(old), expected);
    }

    #[rstest]
    #[case::first(NodeId(0), NodeId(2))]
    #[case::middle(NodeId(1), NodeId(0))]
    #[case::last(NodeId(2), NodeId(5))]
    fn test_graph_remapping_map_node(
        graph_remapping: GraphRemapping,
        #[case] old: NodeId,
        #[case] expected: NodeId,
    ) {
        assert_eq!(graph_remapping.map_node(old), expected);
    }

    #[rstest]
    #[should_panic(expected = "node id outside remapping source range")]
    fn test_graph_remapping_map_node_error(graph_remapping: GraphRemapping) {
        graph_remapping.map_node(NodeId(3));
    }

    #[rstest]
    #[case::relabel(EdgeId(0), EdgeId(3))]
    #[case::fixed(EdgeId(1), EdgeId(1))]
    fn test_graph_remapping_map_edge(
        graph_remapping: GraphRemapping,
        #[case] old: EdgeId,
        #[case] expected: EdgeId,
    ) {
        assert_eq!(graph_remapping.map_edge(old), expected);
    }

    #[rstest]
    #[should_panic(expected = "edge id outside remapping source range")]
    fn test_graph_remapping_map_edge_error(graph_remapping: GraphRemapping) {
        graph_remapping.map_edge(EdgeId(2));
    }
}
