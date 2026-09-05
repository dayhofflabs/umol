//! Bijective renumbering of a dense id space.

use std::error::Error;
use std::fmt::{self, Debug, Display, Formatter};

use index_vec::{Idx, IndexVec};

use crate::graph::{EdgeId, NodeId};

/// Failure to construct a permutation of a dense id space.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RemappingError<Id> {
    ImageOutOfRange { id: Id, count: usize },
    DuplicateImage { id: Id },
}

impl<Id: Debug> Display for RemappingError<Id> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::ImageOutOfRange { id, count } => {
                write!(f, "image {id:?} is out of range for {count} entries")
            }
            Self::DuplicateImage { id } => write!(f, "image {id:?} occurs more than once"),
        }
    }
}

impl<Id: Debug> Error for RemappingError<Id> {}

/// Bijective renumbering of one dense id space.
///
/// The image vector defines the source domain: source id `i` maps to `images[i]`. Source
/// and target ids therefore range over `0..images.len()`. Each image occurs exactly once.
/// Compatibility with an independently supplied object's id space is checked by its consumer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Remapping<Id: Idx> {
    images: IndexVec<Id, Id>,
}

impl<Id: Idx> Default for Remapping<Id> {
    fn default() -> Self {
        Self {
            images: IndexVec::new(),
        }
    }
}

impl<Id: Idx> Remapping<Id> {
    /// Construct a permutation of `0..images.len()` without changing the supplied images.
    ///
    /// # Errors
    ///
    /// Returns an error for the first out-of-range or repeated image.
    pub fn new(images: Vec<Id>) -> Result<Self, RemappingError<Id>> {
        let mut seen = IndexVec::<Id, bool>::from_vec(vec![false; images.len()]);
        for &id in &images {
            let entry = seen.get_mut(id).ok_or(RemappingError::ImageOutOfRange {
                id,
                count: images.len(),
            })?;
            if *entry {
                return Err(RemappingError::DuplicateImage { id });
            }
            *entry = true;
        }
        Ok(Self {
            images: IndexVec::from_vec(images),
        })
    }

    /// Return the shared source and target size.
    pub fn len(&self) -> usize {
        self.images.len()
    }

    /// Whether the source and target spaces are empty.
    pub fn is_empty(&self) -> bool {
        self.images.is_empty()
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

/// Independent bijective renumberings of the dense node and edge id spaces.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GraphRemapping {
    nodes: Remapping<NodeId>,
    edges: Remapping<EdgeId>,
}

impl GraphRemapping {
    /// Assemble already-valid node and edge permutations.
    pub fn new(nodes: Remapping<NodeId>, edges: Remapping<EdgeId>) -> Self {
        Self { nodes, edges }
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
    #[case(RemappingError::ImageOutOfRange { id: NodeId(3), count: 2 }, "image NodeId(3) is out of range for 2 entries")]
    #[case(RemappingError::DuplicateImage { id: NodeId(1) }, "image NodeId(1) occurs more than once")]
    fn test_remapping_error_display(#[case] error: RemappingError<NodeId>, #[case] expected: &str) {
        assert_eq!(error.to_string(), expected);
    }

    #[rstest]
    #[case::empty(vec![])]
    #[case::singleton(vec![NodeId(0)])]
    #[case::identity(vec![NodeId(0), NodeId(1), NodeId(2)])]
    #[case::crossing(vec![NodeId(2), NodeId(0), NodeId(1)])]
    fn test_remapping_new(#[case] images: Vec<NodeId>) {
        let expected = Remapping {
            images: IndexVec::from_vec(images.clone()),
        };
        assert_eq!(Remapping::new(images), Ok(expected));
    }

    #[rstest]
    #[case::boundary(vec![NodeId(1)], RemappingError::ImageOutOfRange { id: NodeId(1), count: 1 })]
    #[case::maximum(vec![NodeId(u32::MAX)], RemappingError::ImageOutOfRange { id: NodeId(u32::MAX), count: 1 })]
    #[case::sparse(vec![NodeId(5), NodeId(1)], RemappingError::ImageOutOfRange { id: NodeId(5), count: 2 })]
    #[case::repeated(vec![NodeId(1), NodeId(1)], RemappingError::DuplicateImage { id: NodeId(1) })]
    #[case::first_error(vec![NodeId(0), NodeId(0), NodeId(3)], RemappingError::DuplicateImage { id: NodeId(0) })]
    fn test_remapping_new_error(
        #[case] images: Vec<NodeId>,
        #[case] expected: RemappingError<NodeId>,
    ) {
        assert_eq!(Remapping::new(images), Err(expected));
    }

    #[rstest]
    #[case(0)]
    #[case(1)]
    #[case(2)]
    #[case(3)]
    #[case(4)]
    fn test_remapping_new_exhaustive(#[case] count: u32) {
        for mut encoded in 0..(count + 1).pow(count) {
            let images: Vec<_> = (0..count)
                .map(|_| {
                    let id = NodeId(encoded % (count + 1));
                    encoded /= count + 1;
                    id
                })
                .collect();
            let mut sorted = images.clone();
            sorted.sort();
            let expected: Vec<_> = (0..count).map(NodeId).collect();
            let result = Remapping::new(images.clone());
            assert_eq!(result.is_ok(), sorted == expected, "images: {images:?}");
            if let Ok(remapping) = result {
                let actual: Vec<_> = (0..count).map(|idx| remapping.map(NodeId(idx))).collect();
                assert_eq!(actual, images);
            }
        }
    }

    #[rstest]
    #[case::empty(vec![], 0)]
    #[case::permutation(vec![NodeId(1), NodeId(0)], 2)]
    fn test_remapping_len(#[case] images: Vec<NodeId>, #[case] expected: usize) {
        assert_eq!(Remapping::new(images).unwrap().len(), expected);
    }

    #[rstest]
    #[case::empty(vec![], true)]
    #[case::permutation(vec![NodeId(1), NodeId(0)], false)]
    fn test_remapping_is_empty(#[case] images: Vec<NodeId>, #[case] expected: bool) {
        assert_eq!(Remapping::new(images).unwrap().is_empty(), expected);
    }

    #[rstest]
    #[case::first(vec![NodeId(1), NodeId(0)], NodeId(0), Some(NodeId(1)))]
    #[case::last(vec![NodeId(1), NodeId(0)], NodeId(1), Some(NodeId(0)))]
    #[case::empty(vec![], NodeId(0), None)]
    #[case::uncovered(vec![NodeId(0)], NodeId(1), None)]
    fn test_remapping_try_map(
        #[case] images: Vec<NodeId>,
        #[case] old: NodeId,
        #[case] expected: Option<NodeId>,
    ) {
        assert_eq!(Remapping::new(images).unwrap().try_map(old), expected);
    }

    #[rstest]
    #[case::first(NodeId(0), NodeId(1))]
    #[case::last(NodeId(1), NodeId(0))]
    fn test_remapping_map(#[case] old: NodeId, #[case] expected: NodeId) {
        assert_eq!(
            Remapping::new(vec![NodeId(1), NodeId(0)]).unwrap().map(old),
            expected
        );
    }

    #[rstest]
    #[should_panic(expected = "id outside remapping source domain")]
    fn test_remapping_map_error() {
        Remapping::<NodeId>::new(vec![]).unwrap().map(NodeId(0));
    }

    #[rstest]
    fn test_remapping_default() {
        assert_eq!(
            Remapping::<NodeId>::default(),
            Remapping::new(vec![]).unwrap()
        );
    }

    #[fixture]
    fn graph_remapping() -> GraphRemapping {
        GraphRemapping::new(
            Remapping::new(vec![NodeId(2), NodeId(0), NodeId(1)]).unwrap(),
            Remapping::new(vec![EdgeId(1), EdgeId(0)]).unwrap(),
        )
    }

    #[rstest]
    fn test_graph_remapping_new(graph_remapping: GraphRemapping) {
        assert_eq!(
            graph_remapping,
            GraphRemapping {
                nodes: Remapping::new(vec![NodeId(2), NodeId(0), NodeId(1)]).unwrap(),
                edges: Remapping::new(vec![EdgeId(1), EdgeId(0)]).unwrap(),
            }
        );
    }

    #[rstest]
    #[case::first(NodeId(0), Some(NodeId(2)))]
    #[case::last(NodeId(2), Some(NodeId(1)))]
    #[case::uncovered(NodeId(3), None)]
    fn test_graph_remapping_try_map_node(
        graph_remapping: GraphRemapping,
        #[case] old: NodeId,
        #[case] expected: Option<NodeId>,
    ) {
        assert_eq!(graph_remapping.try_map_node(old), expected);
    }

    #[rstest]
    #[case::first(EdgeId(0), Some(EdgeId(1)))]
    #[case::last(EdgeId(1), Some(EdgeId(0)))]
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
    #[case::last(NodeId(2), NodeId(1))]
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
    #[case::first(EdgeId(0), EdgeId(1))]
    #[case::last(EdgeId(1), EdgeId(0))]
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
