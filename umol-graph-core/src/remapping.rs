//! Total relabeling of an id space.
//!
//! A remapping gives every source id an image and drops nothing, so it expresses alpha-renaming and
//! injection into a larger ambient space. Removal is [`crate::compaction::Compaction`].

use crate::graph::{EdgeId, NodeId};

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
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

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
}
