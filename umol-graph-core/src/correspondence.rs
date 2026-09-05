//! A partial bijection between two id spaces.
//!
//! Framed as a matching in the bipartite graph over the two id sets: an id is **matched** (paired
//! with a partner on the other side — the shared interface) or **unmatched** (present on only one
//! side). Generic in the id type, so the same carrier serves node correspondences (atoms) and every
//! entity kind (bonds, overlays) one layer up; `Correspondence<NodeId>` additionally exposes the
//! induced edge correspondence over the two graphs.

use std::collections::{BTreeSet, HashMap};
use std::error::Error;
use std::fmt::{self, Debug, Display, Formatter};
use std::ops::{Add, Sub};

use crate::compact::{Compaction, GraphCompaction};
use crate::graph::{EdgeId, Graph, NodeId};
use crate::remap::{GraphRemapping, Remapping};

/// A partial bijection between two `Id` spaces: the matched `(left, right)` pairs; every unmatched
/// id is reported on its side. Only the matched pairs are stored — unmatched ids are derived on
/// demand, so the carrier stays cheap to produce on the hot enumeration path.
///
/// Invariant: `matched_pairs` is a partial bijection within the declared id spaces and is sorted by
/// left id, so `right_of` is a binary search and `left_unmatched` a single merge. `new` validates
/// the partial-bijection invariant and establishes the ordering.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Correspondence<Id> {
    matched_pairs: Vec<(Id, Id)>,
    left_count: usize,
    right_count: usize,
}

impl<Id> From<&Compaction<Id>> for Correspondence<Id>
where
    Id: Copy + Ord + Into<usize> + From<usize> + Add<usize, Output = Id> + Sub<usize, Output = Id>,
{
    /// Preserve the declared counts and every surviving source-to-result pairing.
    fn from(compaction: &Compaction<Id>) -> Self {
        Self {
            matched_pairs: (0..compaction.source_count())
                .map(Id::from)
                .filter_map(|left| compaction.compact(left).map(|right| (left, right)))
                .collect(),
            left_count: compaction.source_count(),
            right_count: compaction.result_count(),
        }
    }
}

impl From<&GraphCompaction> for GraphCorrespondence {
    /// Preserve the node and edge compaction witnesses.
    fn from(compaction: &GraphCompaction) -> Self {
        Self::new(compaction.nodes().into(), compaction.edges().into())
    }
}

impl<Id: Copy + Into<usize> + From<usize>> From<&Remapping<Id>> for Correspondence<Id> {
    /// Preserve every pairing and both counts of the permutation.
    fn from(remapping: &Remapping<Id>) -> Self {
        Self {
            matched_pairs: (0..remapping.len())
                .map(|idx| {
                    let left = Id::from(idx);
                    (left, remapping.map(left))
                })
                .collect(),
            left_count: remapping.len(),
            right_count: remapping.len(),
        }
    }
}

impl<Id: Copy + Ord + From<usize>> Correspondence<Id> {
    /// The unique correspondence between two empty id spaces.
    pub const fn empty() -> Self {
        Self {
            matched_pairs: Vec::new(),
            left_count: 0,
            right_count: 0,
        }
    }

    /// Pair every id with itself in a domain of `count` entities.
    pub fn identity(count: usize) -> Self {
        Self {
            matched_pairs: (0..count)
                .map(|idx| (Id::from(idx), Id::from(idx)))
                .collect(),
            left_count: count,
            right_count: count,
        }
    }

    /// A correspondence from its matched `(left, right)` pairs over two id spaces of the given sizes.
    /// Every id of either side not appearing in `matched_pairs` is unmatched on that side. Pairs are
    /// sorted by left id to establish the lookup invariant (cheap when already sorted).
    pub fn new(
        mut matched_pairs: Vec<(Id, Id)>,
        left_count: usize,
        right_count: usize,
    ) -> Result<Self, CorrespondenceError<Id>> {
        let left_bound = Id::from(left_count);
        let right_bound = Id::from(right_count);
        let mut left_ids = BTreeSet::new();
        let mut right_ids = BTreeSet::new();
        for &(left, right) in &matched_pairs {
            if left >= left_bound {
                return Err(CorrespondenceError::LeftIdOutOfRange {
                    id: left,
                    count: left_count,
                });
            }
            if right >= right_bound {
                return Err(CorrespondenceError::RightIdOutOfRange {
                    id: right,
                    count: right_count,
                });
            }
            if !left_ids.insert(left) {
                return Err(CorrespondenceError::DuplicateLeftId { id: left });
            }
            if !right_ids.insert(right) {
                return Err(CorrespondenceError::DuplicateRightId { id: right });
            }
        }
        matched_pairs.sort_unstable_by_key(|&(left, _)| left);
        Ok(Self {
            matched_pairs,
            left_count,
            right_count,
        })
    }

    /// A correspondence whose left space is dense (`0..images.len()`), pairing each left id `i` with
    /// `images[i]` — the shape a subgraph match / embedding induces (query index → host id). Already
    /// sorted by left, so no sort. Panics if an image is outside the declared right id space or
    /// occurs more than once; algorithmic producers are expected to establish these invariants.
    pub fn from_images(images: &[Id], right_count: usize) -> Self {
        let right_bound = Id::from(right_count);
        let mut right_ids = BTreeSet::new();
        for &right in images {
            assert!(
                right < right_bound,
                "correspondence image is out of range for the declared right id space"
            );
            assert!(
                right_ids.insert(right),
                "correspondence images must be unique"
            );
        }

        Self {
            matched_pairs: images
                .iter()
                .enumerate()
                .map(|(left, &right)| (Id::from(left), right))
                .collect(),
            left_count: images.len(),
            right_count,
        }
    }

    /// The matched `(left, right)` pairs (sorted by left) — the shared interface.
    pub fn matched_pairs(&self) -> &[(Id, Id)] {
        &self.matched_pairs
    }

    /// The number of matched pairs (the interface size).
    pub fn matched_pair_count(&self) -> usize {
        self.matched_pairs.len()
    }

    /// The size of the left id space.
    pub fn left_count(&self) -> usize {
        self.left_count
    }

    /// The size of the right id space.
    pub fn right_count(&self) -> usize {
        self.right_count
    }

    /// Whether every left id is matched.
    pub fn is_total_on_left(&self) -> bool {
        self.matched_pairs.len() == self.left_count
    }

    /// Whether every right id is matched.
    pub fn is_total_on_right(&self) -> bool {
        self.matched_pairs.len() == self.right_count
    }

    /// Whether every id on both sides is matched — a total bijection with no unmatched ids. A diff
    /// through such a correspondence adds and removes nothing.
    pub fn is_total(&self) -> bool {
        self.is_total_on_left() && self.is_total_on_right()
    }

    /// The right partner of a left id, if matched. Binary search (pairs sorted by left).
    pub fn right_of(&self, left: Id) -> Option<Id> {
        self.matched_pairs
            .binary_search_by_key(&left, |&(l, _)| l)
            .ok()
            .map(|index| self.matched_pairs[index].1)
    }

    /// The left partner of a right id, if matched. Linear — the un-indexed reverse direction.
    pub fn left_of(&self, right: Id) -> Option<Id> {
        self.matched_pairs
            .iter()
            .find(|&&(_, r)| r == right)
            .map(|&(l, _)| l)
    }

    /// Left ids with no partner — the ones deleted when read as a transformation. Merge over the
    /// sorted left column.
    pub fn left_unmatched(&self) -> Vec<Id> {
        unmatched(
            self.left_count,
            self.matched_pairs.iter().map(|&(left, _)| left),
        )
    }

    /// Right ids with no partner — the ones created when read as a transformation. Sorts the right
    /// column (unindexed) once, then merges.
    pub fn right_unmatched(&self) -> Vec<Id> {
        let mut rights: Vec<Id> = self.matched_pairs.iter().map(|&(_, right)| right).collect();
        rights.sort_unstable();
        unmatched(self.right_count, rights.into_iter())
    }

    /// Append `count` unmatched ids to the right domain, reusing the pair vector.
    ///
    /// Existing pairings and the left count are unchanged.
    pub fn extend_right(mut self, count: usize) -> Self {
        self.right_count += count;
        self
    }

    /// Compact the right domain in place, discarding pairs whose right id is removed.
    ///
    /// Retains the pair-vector allocation and left count.
    ///
    /// # Errors
    ///
    /// Returns the incompatible counts when the right count differs from the compaction's
    /// source count. The receiver is consumed on failure.
    ///
    /// # Semantic properties
    ///
    /// Equivalent to composition with the compaction's correspondence.
    pub fn compact_right(
        mut self,
        compaction: &Compaction<Id>,
    ) -> Result<Self, CorrespondenceComposeError>
    where
        Id: Into<usize> + Add<usize, Output = Id> + Sub<usize, Output = Id>,
    {
        if self.right_count != compaction.source_count() {
            return Err(CorrespondenceComposeError {
                right_count: self.right_count,
                next_left_count: compaction.source_count(),
            });
        }
        self.matched_pairs.retain_mut(|(_, right)| {
            if let Some(image) = compaction.compact(*right) {
                *right = image;
                true
            } else {
                false
            }
        });
        self.right_count = compaction.result_count();
        Ok(self)
    }

    /// Expand the right domain through the inverse compaction, leaving restored ids unmatched.
    ///
    /// Retains the pair-vector allocation and left count. Does not recreate discarded pairs.
    ///
    /// # Errors
    ///
    /// Returns the incompatible counts when the right count differs from the compaction's
    /// result count. The receiver is consumed on failure.
    ///
    /// # Semantic properties
    ///
    /// Equivalent to composition with the reversed compaction correspondence.
    pub fn uncompact_right(
        mut self,
        compaction: &Compaction<Id>,
    ) -> Result<Self, CorrespondenceComposeError>
    where
        Id: Into<usize> + Add<usize, Output = Id> + Sub<usize, Output = Id>,
    {
        if self.right_count != compaction.result_count() {
            return Err(CorrespondenceComposeError {
                right_count: self.right_count,
                next_left_count: compaction.result_count(),
            });
        }
        for (_, right) in &mut self.matched_pairs {
            *right = compaction.uncompact(*right);
        }
        self.right_count = compaction.source_count();
        Ok(self)
    }

    /// Relational composition: `self` (left↔middle) followed by `other` (middle↔right), yielding a
    /// left↔right correspondence. A left id matched to a middle id that `other` leaves unmatched
    /// becomes unmatched.
    ///
    /// # Errors
    /// Returns the two intermediate counts when they disagree. Equal counts do not establish
    /// that the correspondences describe the same intermediate object.
    ///
    /// # Semantic properties
    /// Composition is associative for compatible carriers and preserves the outer counts.
    pub fn compose(
        &self,
        other: &Correspondence<Id>,
    ) -> Result<Correspondence<Id>, CorrespondenceComposeError> {
        if self.right_count != other.left_count {
            return Err(CorrespondenceComposeError {
                right_count: self.right_count,
                next_left_count: other.left_count,
            });
        }
        let matched_pairs = self
            .matched_pairs
            .iter()
            .filter_map(|&(left, middle)| other.right_of(middle).map(|right| (left, right)))
            .collect();
        Ok(Self {
            matched_pairs,
            left_count: self.left_count,
            right_count: other.right_count,
        })
    }

    /// Compose correspondences in iteration order. Returns `Ok(None)` for an empty input and
    /// `Ok(Some(value))` for a singleton.
    ///
    /// # Errors
    /// Returns the first incompatible pair of intermediate counts.
    pub fn compose_all(
        correspondences: impl IntoIterator<Item = Self>,
    ) -> Result<Option<Self>, CorrespondenceComposeError> {
        let mut correspondences = correspondences.into_iter();
        let Some(first) = correspondences.next() else {
            return Ok(None);
        };
        correspondences
            .try_fold(first, |left, right| left.compose(&right))
            .map(Some)
    }

    /// The inverse correspondence (right↔left): each matched pair swapped and the two id-space sizes
    /// exchanged. A left-unmatched id becomes right-unmatched and vice versa.
    pub fn reverse(&self) -> Correspondence<Id> {
        let matched_pairs = self
            .matched_pairs
            .iter()
            .map(|&(left, right)| (right, left))
            .collect();
        Correspondence::new(matched_pairs, self.right_count, self.left_count)
            .unwrap_or_else(|_| panic!("reversal of a valid correspondence is valid"))
    }
}

impl Correspondence<NodeId> {
    /// The induced edge correspondence: `(left_edge, right_edge)` pairs whose endpoints are matched
    /// to a unique edge on the other side. Returns `None` when the declared node spaces do not
    /// describe `left` and `right`, or when parallel edges make the induced pairing non-unique.
    pub fn edge_matched_pairs(&self, left: &Graph, right: &Graph) -> Option<Vec<(EdgeId, EdgeId)>> {
        if self.left_count != left.node_count() || self.right_count != right.node_count() {
            return None;
        }

        let mut right_edges: HashMap<[NodeId; 2], (EdgeId, bool)> = HashMap::new();
        for edge in right.edge_ids() {
            right_edges
                .entry(right.edge_endpoints(edge))
                .and_modify(|(_, unique)| *unique = false)
                .or_insert((edge, true));
        }

        let mut used_right = BTreeSet::new();
        let mut matched_pairs = Vec::new();
        for left_edge in left.edge_ids() {
            let [left_first, left_second] = left.edge_endpoints(left_edge);
            let (Some(right_first), Some(right_second)) =
                (self.right_of(left_first), self.right_of(left_second))
            else {
                continue;
            };
            let endpoints = if right_first <= right_second {
                [right_first, right_second]
            } else {
                [right_second, right_first]
            };
            let Some(&(right_edge, unique)) = right_edges.get(&endpoints) else {
                continue;
            };
            if !unique || !used_right.insert(right_edge) {
                return None;
            }
            matched_pairs.push((left_edge, right_edge));
        }
        Some(matched_pairs)
    }

    /// The number of edges shared under the node matching (the maximum-common-edge-subgraph
    /// objective), or `None` when the edge correspondence is not uniquely induced.
    pub fn shared_edge_count(&self, left: &Graph, right: &Graph) -> Option<usize> {
        self.edge_matched_pairs(left, right)
            .map(|pairs| pairs.len())
    }
}

/// A node-wise and edge-wise partial bijection between two declared graph id spaces.
///
/// Each component records matched ids and both domain sizes; either side may have unmatched ids.
/// Correspondences support matching, graph rewriting, and composition of operation witnesses.
/// Construction establishes each component's partial bijection, not compatibility with particular
/// graphs or agreement of edge pairings with node incidence. Those are contextual properties of
/// the producer or consumer; the carrier contains no graphs or entity payloads.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GraphCorrespondence {
    nodes: Correspondence<NodeId>,
    edges: Correspondence<EdgeId>,
}

impl From<&GraphRemapping> for GraphCorrespondence {
    /// Preserve the node and edge permutations as complete correspondences.
    fn from(remapping: &GraphRemapping) -> Self {
        Self::new(remapping.nodes().into(), remapping.edges().into())
    }
}

impl GraphCorrespondence {
    pub fn new(nodes: Correspondence<NodeId>, edges: Correspondence<EdgeId>) -> Self {
        Self { nodes, edges }
    }

    /// The graph correspondence a node correspondence induces over `left` / `right`: its edge component
    /// is the induced edge correspondence (an edge matched when both endpoints are matched). Exact
    /// whenever every structurally matched edge is admissible — a subiso match, a common *induced*
    /// subgraph, or a common-subgraph walk whose adjacency enforces the edge predicate; not for an
    /// edge-subgraph result whose producer selects the edge component directly. Returns `None` when the
    /// node correspondence does not describe the supplied graphs or does not induce unique edge
    /// pairs.
    pub fn induce(left: &Graph, right: &Graph, nodes: Correspondence<NodeId>) -> Option<Self> {
        let edges = Correspondence::new(
            nodes.edge_matched_pairs(left, right)?,
            left.edge_count(),
            right.edge_count(),
        )
        .expect("an induced edge correspondence is valid");
        Some(Self { nodes, edges })
    }

    pub fn nodes(&self) -> &Correspondence<NodeId> {
        &self.nodes
    }

    pub fn edges(&self) -> &Correspondence<EdgeId> {
        &self.edges
    }

    /// Append unmatched nodes and edges to the right domains without changing existing pairs.
    pub fn extend_right(self, nodes: usize, edges: usize) -> Self {
        Self::new(
            self.nodes.extend_right(nodes),
            self.edges.extend_right(edges),
        )
    }

    /// Compact both right domains, reusing their pair vectors.
    ///
    /// # Errors
    ///
    /// Returns the first source-count mismatch, nodes before edges. Consumes the receiver.
    ///
    /// # Semantic properties
    ///
    /// Equivalent to composition with the compaction's correspondence.
    pub fn compact_right(
        self,
        compaction: &GraphCompaction,
    ) -> Result<Self, GraphCorrespondenceComposeError> {
        Ok(Self::new(
            self.nodes
                .compact_right(compaction.nodes())
                .map_err(GraphCorrespondenceComposeError::Nodes)?,
            self.edges
                .compact_right(compaction.edges())
                .map_err(GraphCorrespondenceComposeError::Edges)?,
        ))
    }

    /// Expand both right domains through the inverse compaction, leaving restored ids unmatched.
    ///
    /// Reuses the pair vectors and does not recreate discarded pairs.
    ///
    /// # Errors
    ///
    /// Returns the first result-count mismatch, nodes before edges. Consumes the receiver.
    pub fn uncompact_right(
        self,
        compaction: &GraphCompaction,
    ) -> Result<Self, GraphCorrespondenceComposeError> {
        Ok(Self::new(
            self.nodes
                .uncompact_right(compaction.nodes())
                .map_err(GraphCorrespondenceComposeError::Nodes)?,
            self.edges
                .uncompact_right(compaction.edges())
                .map_err(GraphCorrespondenceComposeError::Edges)?,
        ))
    }

    /// Relational composition of the node and edge correspondences.
    ///
    /// # Errors
    /// Returns the first component with unequal intermediate counts, nodes before edges.
    /// Equal counts do not establish intermediate graph identity.
    pub fn compose(
        &self,
        other: &GraphCorrespondence,
    ) -> Result<GraphCorrespondence, GraphCorrespondenceComposeError> {
        Ok(GraphCorrespondence::new(
            self.nodes
                .compose(&other.nodes)
                .map_err(GraphCorrespondenceComposeError::Nodes)?,
            self.edges
                .compose(&other.edges)
                .map_err(GraphCorrespondenceComposeError::Edges)?,
        ))
    }

    /// Compose graph correspondences in iteration order. Returns `Ok(None)` for an empty input and
    /// `Ok(Some(value))` for a singleton.
    ///
    /// # Errors
    /// Returns the first composition error in iteration order.
    pub fn compose_all(
        correspondences: impl IntoIterator<Item = Self>,
    ) -> Result<Option<Self>, GraphCorrespondenceComposeError> {
        let mut correspondences = correspondences.into_iter();
        let Some(first) = correspondences.next() else {
            return Ok(None);
        };
        correspondences
            .try_fold(first, |left, right| left.compose(&right))
            .map(Some)
    }

    /// Whether every node and edge on the left is matched.
    pub fn is_total_on_left(&self) -> bool {
        self.nodes.is_total_on_left() && self.edges.is_total_on_left()
    }

    /// Whether every node and edge on the right is matched.
    pub fn is_total_on_right(&self) -> bool {
        self.nodes.is_total_on_right() && self.edges.is_total_on_right()
    }

    /// Whether every node and edge on both sides is matched.
    pub fn is_total(&self) -> bool {
        self.is_total_on_left() && self.is_total_on_right()
    }
}

/// The ids `0..count` absent from `sorted_matched` (which must be ascending, no duplicates) — a
/// single merge pass, no per-id search.
fn unmatched<Id: Copy + Ord + From<usize>>(
    count: usize,
    sorted_matched: impl Iterator<Item = Id>,
) -> Vec<Id> {
    let mut matched = sorted_matched.peekable();
    (0..count)
        .map(Id::from)
        .filter(|node| {
            if matched.peek() == Some(node) {
                matched.next();
                false
            } else {
                true
            }
        })
        .collect()
}

/// Failure to construct a correspondence whose pairs form a partial bijection over its declared id
/// spaces.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CorrespondenceError<Id> {
    LeftIdOutOfRange { id: Id, count: usize },
    RightIdOutOfRange { id: Id, count: usize },
    DuplicateLeftId { id: Id },
    DuplicateRightId { id: Id },
}

impl<Id: Debug> Display for CorrespondenceError<Id> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::LeftIdOutOfRange { id, count } => {
                write!(f, "left id {id:?} is out of range for {count} entries")
            }
            Self::RightIdOutOfRange { id, count } => {
                write!(f, "right id {id:?} is out of range for {count} entries")
            }
            Self::DuplicateLeftId { id } => write!(f, "left id {id:?} occurs more than once"),
            Self::DuplicateRightId { id } => write!(f, "right id {id:?} occurs more than once"),
        }
    }
}

impl<Id: Debug> Error for CorrespondenceError<Id> {}

/// The consecutive correspondences declare different intermediate sizes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CorrespondenceComposeError {
    pub right_count: usize,
    pub next_left_count: usize,
}

impl Display for CorrespondenceComposeError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "intermediate counts differ: {} and {}",
            self.right_count, self.next_left_count
        )
    }
}

impl Error for CorrespondenceComposeError {}

/// A graph correspondence component has incompatible intermediate counts.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GraphCorrespondenceComposeError {
    Nodes(CorrespondenceComposeError),
    Edges(CorrespondenceComposeError),
}

impl Display for GraphCorrespondenceComposeError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Nodes(error) => write!(f, "nodes: {error}"),
            Self::Edges(error) => write!(f, "edges: {error}"),
        }
    }
}

impl Error for GraphCorrespondenceComposeError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Nodes(error) | Self::Edges(error) => Some(error),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    #[rstest]
    #[case::empty(0, vec![])]
    #[case::one(1, vec![(NodeId(0), NodeId(0))])]
    #[case::three(3, vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1)), (NodeId(2), NodeId(2))])]
    fn test_correspondence_identity(#[case] count: usize, #[case] pairs: Vec<(NodeId, NodeId)>) {
        assert_eq!(
            Correspondence::<NodeId>::identity(count),
            Correspondence {
                matched_pairs: pairs,
                left_count: count,
                right_count: count,
            }
        );
    }

    #[rstest]
    #[case::empty(vec![], 0, 0, 2)]
    #[case::unmatched(vec![], 3, 2, 3)]
    #[case::crossing(vec![(NodeId(0), NodeId(2)), (NodeId(2), NodeId(0))], 4, 3, 2)]
    fn test_correspondence_extend_right(
        #[case] pairs: Vec<(NodeId, NodeId)>,
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] added: usize,
    ) {
        let correspondence = Correspondence {
            matched_pairs: pairs.clone(),
            left_count,
            right_count,
        };
        let ptr = correspondence.matched_pairs.as_ptr();
        let capacity = correspondence.matched_pairs.capacity();
        let result = correspondence.extend_right(added);
        assert_eq!(result.matched_pairs.as_ptr(), ptr);
        assert_eq!(result.matched_pairs.capacity(), capacity);
        assert_eq!(
            result,
            Correspondence {
                matched_pairs: pairs,
                left_count,
                right_count: right_count + added
            }
        );
    }

    #[rstest]
    #[case::empty(Correspondence::empty())]
    #[case::partial(Correspondence::new(vec![(NodeId(1), NodeId(0))], 3, 2).unwrap())]
    fn test_correspondence_extend_right_identity(#[case] correspondence: Correspondence<NodeId>) {
        assert_eq!(correspondence.clone().extend_right(0), correspondence);
    }

    #[rstest]
    #[case::unmatched(vec![NodeId(1)], vec![(NodeId(0), NodeId(1)), (NodeId(2), NodeId(0))])]
    #[case::matched(vec![NodeId(0)], vec![(NodeId(0), NodeId(1))])]
    #[case::all(vec![NodeId(0), NodeId(1), NodeId(2)], vec![])]
    fn test_correspondence_compact_right(
        #[case] removed: Vec<NodeId>,
        #[case] expected: Vec<(NodeId, NodeId)>,
    ) {
        let correspondence = Correspondence {
            matched_pairs: vec![(NodeId(0), NodeId(2)), (NodeId(2), NodeId(0))],
            left_count: 4,
            right_count: 3,
        };
        let ptr = correspondence.matched_pairs.as_ptr();
        let capacity = correspondence.matched_pairs.capacity();
        let result_count = 3 - removed.len();
        let compaction = Compaction::new(3, removed).unwrap();
        let result = correspondence.compact_right(&compaction).unwrap();
        assert_eq!(result.matched_pairs.as_ptr(), ptr);
        assert_eq!(result.matched_pairs.capacity(), capacity);
        assert_eq!(
            result,
            Correspondence {
                matched_pairs: expected,
                left_count: 4,
                right_count: result_count
            }
        );
    }

    #[rstest]
    #[case::empty(Correspondence::empty())]
    #[case::partial(Correspondence::new(vec![(NodeId(1), NodeId(0))], 3, 2).unwrap())]
    fn test_correspondence_compact_right_identity(#[case] correspondence: Correspondence<NodeId>) {
        let compaction = Compaction::identity(correspondence.right_count());
        assert_eq!(
            correspondence.clone().compact_right(&compaction),
            Ok(correspondence)
        );
    }

    #[rstest]
    #[case::smaller(2)]
    #[case::larger(4)]
    fn test_correspondence_compact_right_error(#[case] source_count: usize) {
        let correspondence = Correspondence::<NodeId>::new(vec![], 2, 3).unwrap();
        assert_eq!(
            correspondence.compact_right(&Compaction::identity(source_count)),
            Err(CorrespondenceComposeError {
                right_count: 3,
                next_left_count: source_count
            })
        );
    }

    #[rstest]
    #[case::leading(vec![NodeId(0)], vec![(NodeId(0), NodeId(3)), (NodeId(2), NodeId(1))])]
    #[case::middle(vec![NodeId(1)], vec![(NodeId(0), NodeId(3)), (NodeId(2), NodeId(0))])]
    #[case::trailing(vec![NodeId(3)], vec![(NodeId(0), NodeId(2)), (NodeId(2), NodeId(0))])]
    fn test_correspondence_uncompact_right(
        #[case] removed: Vec<NodeId>,
        #[case] expected: Vec<(NodeId, NodeId)>,
    ) {
        let correspondence = Correspondence {
            matched_pairs: vec![(NodeId(0), NodeId(2)), (NodeId(2), NodeId(0))],
            left_count: 4,
            right_count: 3,
        };
        let ptr = correspondence.matched_pairs.as_ptr();
        let capacity = correspondence.matched_pairs.capacity();
        let compaction = Compaction::new(4, removed).unwrap();
        let result = correspondence.uncompact_right(&compaction).unwrap();
        assert_eq!(result.matched_pairs.as_ptr(), ptr);
        assert_eq!(result.matched_pairs.capacity(), capacity);
        assert_eq!(
            result,
            Correspondence {
                matched_pairs: expected,
                left_count: 4,
                right_count: 4
            }
        );
    }

    #[rstest]
    #[case::empty(Correspondence::empty())]
    #[case::partial(Correspondence::new(vec![(NodeId(1), NodeId(0))], 3, 2).unwrap())]
    fn test_correspondence_uncompact_right_identity(
        #[case] correspondence: Correspondence<NodeId>,
    ) {
        let compaction = Compaction::identity(correspondence.right_count());
        assert_eq!(
            correspondence.clone().uncompact_right(&compaction),
            Ok(correspondence)
        );
    }

    #[rstest]
    #[case::smaller(3)]
    #[case::larger(5)]
    fn test_correspondence_uncompact_right_error(#[case] source_count: usize) {
        let correspondence = Correspondence::<NodeId>::new(vec![], 2, 3).unwrap();
        let compaction = Compaction::new(source_count, vec![NodeId(0)]).unwrap();
        assert_eq!(
            correspondence.uncompact_right(&compaction),
            Err(CorrespondenceComposeError {
                right_count: 3,
                next_left_count: source_count - 1
            })
        );
    }

    fn e(i: u32) -> EdgeId {
        EdgeId(i)
    }

    #[rstest]
    #[case::empty(0, vec![], vec![], 0)]
    #[case::identity(2, vec![], vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2)]
    #[case::partial(4, vec![NodeId(1), NodeId(3)], vec![(NodeId(0), NodeId(0)), (NodeId(2), NodeId(1))], 2)]
    #[case::full(2, vec![NodeId(1), NodeId(0)], vec![], 0)]
    fn test_correspondence_from_compaction(
        #[case] source_count: usize,
        #[case] removed: Vec<NodeId>,
        #[case] pairs: Vec<(NodeId, NodeId)>,
        #[case] result_count: usize,
    ) {
        let compaction = Compaction::new(source_count, removed).unwrap();
        assert_eq!(
            Correspondence::from(&compaction),
            Correspondence {
                matched_pairs: pairs,
                left_count: source_count,
                right_count: result_count,
            }
        );
    }

    #[rstest]
    fn test_graph_correspondence_from_compaction() {
        let compaction = GraphCompaction::new(
            Compaction::new(3, vec![NodeId(1)]).unwrap(),
            Compaction::new(2, vec![EdgeId(0), EdgeId(1)]).unwrap(),
        );
        assert_eq!(
            GraphCorrespondence::from(&compaction),
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(2), NodeId(1))], 3, 2)
                    .unwrap(),
                Correspondence::new(vec![], 2, 0).unwrap(),
            )
        );
    }

    #[rstest]
    #[case::empty(vec![], vec![])]
    #[case::identity(vec![NodeId(0)], vec![(NodeId(0), NodeId(0))])]
    #[case::crossing(vec![NodeId(2), NodeId(0), NodeId(1)], vec![(NodeId(0), NodeId(2)), (NodeId(1), NodeId(0)), (NodeId(2), NodeId(1))])]
    fn test_correspondence_from_remapping(
        #[case] images: Vec<NodeId>,
        #[case] pairs: Vec<(NodeId, NodeId)>,
    ) {
        let count = images.len();
        let remapping = Remapping::new(images).unwrap();
        assert_eq!(
            Correspondence::from(&remapping),
            Correspondence {
                matched_pairs: pairs,
                left_count: count,
                right_count: count,
            }
        );
    }

    #[rstest]
    #[case::empty(vec![], vec![], vec![], vec![])]
    #[case::crossing(
            vec![NodeId(1), NodeId(0)], vec![EdgeId(2), EdgeId(0), EdgeId(1)],
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
            vec![(EdgeId(0), EdgeId(2)), (EdgeId(1), EdgeId(0)), (EdgeId(2), EdgeId(1))]
        )]
    fn test_graph_correspondence_from_remapping(
        #[case] nodes: Vec<NodeId>,
        #[case] edges: Vec<EdgeId>,
        #[case] node_pairs: Vec<(NodeId, NodeId)>,
        #[case] edge_pairs: Vec<(EdgeId, EdgeId)>,
    ) {
        let node_count = nodes.len();
        let edge_count = edges.len();
        let remapping = GraphRemapping::new(
            Remapping::new(nodes).unwrap(),
            Remapping::new(edges).unwrap(),
        );
        assert_eq!(
            GraphCorrespondence::from(&remapping),
            GraphCorrespondence::new(
                Correspondence {
                    matched_pairs: node_pairs,
                    left_count: node_count,
                    right_count: node_count
                },
                Correspondence {
                    matched_pairs: edge_pairs,
                    left_count: edge_count,
                    right_count: edge_count
                },
            )
        );
    }

    #[fixture]
    fn paths() -> (Graph, Graph, Correspondence<NodeId>) {
        // left path 0-1-2 mapped onto the interior of right path 0-1-2-3; right node 0 unmatched.
        (
            Graph::new(3, &[[0, 1], [1, 2]]),
            Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
            Correspondence::new(vec![(n(0), n(1)), (n(1), n(2)), (n(2), n(3))], 3, 4)
                .expect("correspondence producer preserves partial-bijection invariants"),
        )
    }

    #[fixture]
    fn update_graph_correspondence() -> GraphCorrespondence {
        GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(2)), (NodeId(2), NodeId(0))], 3, 4)
                .unwrap(),
            Correspondence::new(vec![(EdgeId(1), EdgeId(1)), (EdgeId(2), EdgeId(0))], 3, 2)
                .unwrap(),
        )
    }

    #[rstest]
    #[case::nodes(2, 0)]
    #[case::edges(0, 1)]
    #[case::both(1, 3)]
    fn test_graph_correspondence_extend_right(
        update_graph_correspondence: GraphCorrespondence,
        #[case] nodes: usize,
        #[case] edges: usize,
    ) {
        let node_ptr = update_graph_correspondence.nodes().matched_pairs().as_ptr();
        let edge_ptr = update_graph_correspondence.edges().matched_pairs().as_ptr();
        let result = update_graph_correspondence.extend_right(nodes, edges);
        assert_eq!(result.nodes().matched_pairs().as_ptr(), node_ptr);
        assert_eq!(result.edges().matched_pairs().as_ptr(), edge_ptr);
        assert_eq!(
            result,
            GraphCorrespondence::new(
                Correspondence::new(
                    vec![(NodeId(0), NodeId(2)), (NodeId(2), NodeId(0))],
                    3,
                    4 + nodes
                )
                .unwrap(),
                Correspondence::new(
                    vec![(EdgeId(1), EdgeId(1)), (EdgeId(2), EdgeId(0))],
                    3,
                    2 + edges
                )
                .unwrap(),
            )
        );
    }

    #[rstest]
    fn test_graph_correspondence_compact_right(update_graph_correspondence: GraphCorrespondence) {
        let compaction = GraphCompaction::new(
            Compaction::new(4, vec![NodeId(1)]).unwrap(),
            Compaction::new(2, vec![EdgeId(0)]).unwrap(),
        );
        let node_ptr = update_graph_correspondence.nodes().matched_pairs().as_ptr();
        let edge_ptr = update_graph_correspondence.edges().matched_pairs().as_ptr();
        let result = update_graph_correspondence
            .compact_right(&compaction)
            .unwrap();
        assert_eq!(result.nodes().matched_pairs().as_ptr(), node_ptr);
        assert_eq!(result.edges().matched_pairs().as_ptr(), edge_ptr);
        assert_eq!(
            result,
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(1)), (NodeId(2), NodeId(0))], 3, 3)
                    .unwrap(),
                Correspondence::new(vec![(EdgeId(1), EdgeId(0))], 3, 1).unwrap(),
            )
        );
    }

    #[rstest]
    #[case::nodes(3, 2, GraphCorrespondenceComposeError::Nodes(CorrespondenceComposeError { right_count: 4, next_left_count: 3 }))]
    #[case::edges(4, 3, GraphCorrespondenceComposeError::Edges(CorrespondenceComposeError { right_count: 2, next_left_count: 3 }))]
    #[case::both(3, 3, GraphCorrespondenceComposeError::Nodes(CorrespondenceComposeError { right_count: 4, next_left_count: 3 }))]
    fn test_graph_correspondence_compact_right_error(
        update_graph_correspondence: GraphCorrespondence,
        #[case] nodes: usize,
        #[case] edges: usize,
        #[case] expected: GraphCorrespondenceComposeError,
    ) {
        let compaction =
            GraphCompaction::new(Compaction::identity(nodes), Compaction::identity(edges));
        assert_eq!(
            update_graph_correspondence.compact_right(&compaction),
            Err(expected)
        );
    }

    #[rstest]
    fn test_graph_correspondence_uncompact_right(update_graph_correspondence: GraphCorrespondence) {
        let compaction = GraphCompaction::new(
            Compaction::new(5, vec![NodeId(1)]).unwrap(),
            Compaction::new(3, vec![EdgeId(0)]).unwrap(),
        );
        let node_ptr = update_graph_correspondence.nodes().matched_pairs().as_ptr();
        let edge_ptr = update_graph_correspondence.edges().matched_pairs().as_ptr();
        let result = update_graph_correspondence
            .uncompact_right(&compaction)
            .unwrap();
        assert_eq!(result.nodes().matched_pairs().as_ptr(), node_ptr);
        assert_eq!(result.edges().matched_pairs().as_ptr(), edge_ptr);
        assert_eq!(
            result,
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(3)), (NodeId(2), NodeId(0))], 3, 5)
                    .unwrap(),
                Correspondence::new(vec![(EdgeId(1), EdgeId(2)), (EdgeId(2), EdgeId(1))], 3, 3)
                    .unwrap(),
            )
        );
    }

    #[rstest]
    #[case::nodes(3, 2, GraphCorrespondenceComposeError::Nodes(CorrespondenceComposeError { right_count: 4, next_left_count: 3 }))]
    #[case::edges(4, 3, GraphCorrespondenceComposeError::Edges(CorrespondenceComposeError { right_count: 2, next_left_count: 3 }))]
    #[case::both(3, 3, GraphCorrespondenceComposeError::Nodes(CorrespondenceComposeError { right_count: 4, next_left_count: 3 }))]
    fn test_graph_correspondence_uncompact_right_error(
        update_graph_correspondence: GraphCorrespondence,
        #[case] nodes: usize,
        #[case] edges: usize,
        #[case] expected: GraphCorrespondenceComposeError,
    ) {
        let compaction =
            GraphCompaction::new(Compaction::identity(nodes), Compaction::identity(edges));
        assert_eq!(
            update_graph_correspondence.uncompact_right(&compaction),
            Err(expected)
        );
    }

    #[rstest]
    fn test_graph_correspondence_extend_right_identity(
        update_graph_correspondence: GraphCorrespondence,
    ) {
        assert_eq!(
            update_graph_correspondence.clone().extend_right(0, 0),
            update_graph_correspondence
        );
    }

    #[rstest]
    fn test_graph_correspondence_compact_right_identity(
        update_graph_correspondence: GraphCorrespondence,
    ) {
        let compaction = GraphCompaction::new(Compaction::identity(4), Compaction::identity(2));
        assert_eq!(
            update_graph_correspondence
                .clone()
                .compact_right(&compaction),
            Ok(update_graph_correspondence)
        );
        let empty = GraphCorrespondence::new(Correspondence::empty(), Correspondence::empty());
        assert_eq!(
            empty.clone().compact_right(&GraphCompaction::empty()),
            Ok(empty)
        );
    }

    #[rstest]
    fn test_graph_correspondence_uncompact_right_identity(
        update_graph_correspondence: GraphCorrespondence,
    ) {
        let compaction = GraphCompaction::new(Compaction::identity(4), Compaction::identity(2));
        assert_eq!(
            update_graph_correspondence
                .clone()
                .uncompact_right(&compaction),
            Ok(update_graph_correspondence)
        );
        let empty = GraphCorrespondence::new(Correspondence::empty(), Correspondence::empty());
        assert_eq!(
            empty.clone().uncompact_right(&GraphCompaction::empty()),
            Ok(empty)
        );
    }

    #[fixture]
    fn graph_correspondences() -> [GraphCorrespondence; 3] {
        [
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ),
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(1), NodeId(2))], 2, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(EdgeId(0), EdgeId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ),
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(2), NodeId(0))], 3, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(EdgeId(1), EdgeId(0))], 2, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ),
        ]
    }

    #[rstest]
    fn test_correspondence_empty() {
        let correspondence = Correspondence::<NodeId>::empty();

        assert_eq!(correspondence.matched_pairs(), []);
        assert_eq!(correspondence.left_count(), 0);
        assert_eq!(correspondence.right_count(), 0);
        assert!(correspondence.is_total());
    }

    #[rstest]
    #[case::empty(Vec::new(), 0, 0, Correspondence { matched_pairs: Vec::new(), left_count: 0, right_count: 0 })]
    #[case::partial(
        vec![(NodeId(1), NodeId(3))],
        3,
        4,
        Correspondence { matched_pairs: vec![(NodeId(1), NodeId(3))], left_count: 3, right_count: 4 },
    )]
    #[case::unsorted(
        vec![(NodeId(2), NodeId(0)), (NodeId(0), NodeId(3)), (NodeId(1), NodeId(1))],
        3,
        4,
        Correspondence {
            matched_pairs: vec![(NodeId(0), NodeId(3)), (NodeId(1), NodeId(1)), (NodeId(2), NodeId(0))],
            left_count: 3,
            right_count: 4,
        },
    )]
    fn test_correspondence_new(
        #[case] matched_pairs: Vec<(NodeId, NodeId)>,
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] expected: Correspondence<NodeId>,
    ) {
        assert_eq!(
            Correspondence::new(matched_pairs, left_count, right_count),
            Ok(expected),
        );
    }

    #[rstest]
    #[case::left_out_of_range(
        vec![(NodeId(2), NodeId(0))],
        2,
        1,
        CorrespondenceError::LeftIdOutOfRange { id: NodeId(2), count: 2 },
    )]
    #[case::right_out_of_range(
        vec![(NodeId(0), NodeId(2))],
        1,
        2,
        CorrespondenceError::RightIdOutOfRange { id: NodeId(2), count: 2 },
    )]
    #[case::duplicate_left(
        vec![(NodeId(0), NodeId(0)), (NodeId(0), NodeId(1))],
        1,
        2,
        CorrespondenceError::DuplicateLeftId { id: NodeId(0) },
    )]
    #[case::duplicate_right(
        vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(1))],
        2,
        2,
        CorrespondenceError::DuplicateRightId { id: NodeId(1) },
    )]
    fn test_correspondence_new_error(
        #[case] matched_pairs: Vec<(NodeId, NodeId)>,
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] expected: CorrespondenceError<NodeId>,
    ) {
        assert_eq!(
            Correspondence::new(matched_pairs, left_count, right_count),
            Err(expected),
        );
    }

    #[rstest]
    fn test_correspondence_from_images() {
        // dense left 0,1,2 → images 3,1,0; every left matched, host id 2 right-unmatched.
        let c = Correspondence::from_images(&[n(3), n(1), n(0)], 4);
        assert_eq!(
            c.matched_pairs(),
            &[(n(0), n(3)), (n(1), n(1)), (n(2), n(0))]
        );
        assert_eq!(c.matched_pair_count(), 3);
        assert_eq!(c.left_unmatched(), Vec::<NodeId>::new());
        assert_eq!(c.right_unmatched(), vec![n(2)]);
    }

    #[rstest]
    #[case::at_bound(&[NodeId(2)], 2)]
    #[case::above_bound(&[NodeId(3)], 2)]
    #[should_panic(
        expected = "correspondence image is out of range for the declared right id space"
    )]
    fn test_correspondence_from_images_range(
        #[case] images: &[NodeId],
        #[case] right_count: usize,
    ) {
        Correspondence::from_images(images, right_count);
    }

    #[rstest]
    #[case::adjacent(&[NodeId(0), NodeId(0)])]
    #[case::separated(&[NodeId(0), NodeId(1), NodeId(0)])]
    #[should_panic(expected = "correspondence images must be unique")]
    fn test_correspondence_from_images_duplicate(#[case] images: &[NodeId]) {
        Correspondence::from_images(images, 2);
    }

    #[rstest]
    fn test_correspondence_matched_pairs() {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4)
            .expect("correspondence producer preserves partial-bijection invariants");
        assert_eq!(c.matched_pairs(), &[(n(0), n(2)), (n(1), n(3))]);
        assert_eq!(c.matched_pair_count(), 2);
    }

    #[rstest]
    fn test_graph_correspondence_new() {
        let c = GraphCorrespondence::new(
            Correspondence::from_images(&[n(1), n(0)], 3),
            Correspondence::from_images(&[e(2)], 4),
        );
        assert_eq!(c.nodes().matched_pairs(), &[(n(0), n(1)), (n(1), n(0))]);
        assert_eq!(c.edges().matched_pairs(), &[(e(0), e(2))]);
    }

    #[rstest]
    #[case::paths(
        Graph::new(3, &[[0, 1], [1, 2]]),
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        Correspondence::new(vec![(n(0), n(1)), (n(1), n(2)), (n(2), n(3))], 3, 4)
            .expect("correspondence producer preserves partial-bijection invariants"),
        GraphCorrespondence::new(
            Correspondence::new(vec![(n(0), n(1)), (n(1), n(2)), (n(2), n(3))], 3, 4)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(e(0), e(1)), (e(1), e(2))], 2, 3)
                .expect("correspondence producer preserves partial-bijection invariants"),
        ),
    )]
    #[case::empty(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        Correspondence::new(vec![], 2, 2)
            .expect("correspondence producer preserves partial-bijection invariants"),
        GraphCorrespondence::new(
            Correspondence::new(vec![], 2, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 1, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
        ),
    )]
    #[case::one_sided_edge(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[]),
        Correspondence::from_images(&[n(0), n(1)], 2),
        GraphCorrespondence::new(
            Correspondence::new(vec![(n(0), n(0)), (n(1), n(1))], 2, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 1, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
        ),
    )]
    fn test_graph_correspondence_induce(
        #[case] left: Graph,
        #[case] right: Graph,
        #[case] nodes: Correspondence<NodeId>,
        #[case] expected: GraphCorrespondence,
    ) {
        assert_eq!(
            GraphCorrespondence::induce(&left, &right, nodes),
            Some(expected)
        );
    }

    #[rstest]
    #[case::left_count(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        Correspondence::new(vec![(n(0), n(0))], 1, 2)
            .expect("correspondence producer preserves partial-bijection invariants"),
    )]
    #[case::right_count(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        Correspondence::new(vec![(n(0), n(0))], 2, 1)
            .expect("correspondence producer preserves partial-bijection invariants"),
    )]
    #[case::parallel_left(
        Graph::new(2, &[[0, 1], [0, 1]]),
        Graph::new(2, &[[0, 1]]),
        Correspondence::from_images(&[n(0), n(1)], 2),
    )]
    #[case::parallel_right(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1], [0, 1]]),
        Correspondence::from_images(&[n(0), n(1)], 2),
    )]
    fn test_graph_correspondence_induce_error(
        #[case] left: Graph,
        #[case] right: Graph,
        #[case] nodes: Correspondence<NodeId>,
    ) {
        assert_eq!(GraphCorrespondence::induce(&left, &right, nodes), None);
    }

    #[rstest]
    fn test_graph_correspondence_compose(graph_correspondences: [GraphCorrespondence; 3]) {
        let [left, right, _] = graph_correspondences;

        assert_eq!(
            left.compose(&right).unwrap(),
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(2))], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(EdgeId(0), EdgeId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ),
        );
    }

    #[rstest]
    #[case::nodes(GraphCorrespondence::new(Correspondence::new(vec![], 1, 2).unwrap(), Correspondence::empty()), GraphCorrespondenceComposeError::Nodes(CorrespondenceComposeError {right_count: 2, next_left_count: 0}))]
    #[case::edges(GraphCorrespondence::new(Correspondence::empty(), Correspondence::new(vec![], 1, 2).unwrap()), GraphCorrespondenceComposeError::Edges(CorrespondenceComposeError {right_count: 2, next_left_count: 0}))]
    fn test_graph_correspondence_compose_error(
        #[case] left: GraphCorrespondence,
        #[case] expected: GraphCorrespondenceComposeError,
    ) {
        let right = GraphCorrespondence::new(Correspondence::empty(), Correspondence::empty());
        assert_eq!(left.compose(&right), Err(expected.clone()));
        assert_eq!(
            GraphCorrespondence::compose_all([left, right]),
            Err(expected)
        );
    }

    #[rstest]
    #[case::empty(0)]
    #[case::singleton(1)]
    #[case::multiple(3)]
    fn test_graph_correspondence_compose_all(
        graph_correspondences: [GraphCorrespondence; 3],
        #[case] count: usize,
    ) {
        let expected = match count {
            0 => None,
            1 => Some(graph_correspondences[0].clone()),
            3 => Some(GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(EdgeId(0), EdgeId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            )),
            _ => unreachable!(),
        };

        assert_eq!(
            GraphCorrespondence::compose_all(graph_correspondences.into_iter().take(count))
                .unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::total(
        GraphCorrespondence::new(
            Correspondence::from_images(&[n(0)], 1),
            Correspondence::from_images(&[e(0)], 1),
        ),
        (true, true, true),
    )]
    #[case::left(
        GraphCorrespondence::new(
            Correspondence::from_images(&[n(0)], 2),
            Correspondence::from_images(&[e(0)], 2),
        ),
        (true, false, false),
    )]
    #[case::right(
        GraphCorrespondence::new(
            Correspondence::new(vec![(n(0), n(0))], 2, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(e(0), e(0))], 2, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
        ),
        (false, true, false),
    )]
    #[case::mixed_node_and_edge(
        GraphCorrespondence::new(
            Correspondence::from_images(&[n(0)], 2),
            Correspondence::new(vec![(e(0), e(0))], 2, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
        ),
        (false, false, false),
    )]
    fn test_graph_correspondence_is_total(
        #[case] correspondence: GraphCorrespondence,
        #[case] expected: (bool, bool, bool),
    ) {
        assert_eq!(
            (
                correspondence.is_total_on_left(),
                correspondence.is_total_on_right(),
                correspondence.is_total(),
            ),
            expected,
        );
    }

    #[rstest]
    #[case::matched_first(n(0), Some(n(2)))]
    #[case::matched_second(n(1), Some(n(3)))]
    #[case::unmatched(n(2), None)]
    fn test_correspondence_right_of(#[case] left: NodeId, #[case] expected: Option<NodeId>) {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4)
            .expect("correspondence producer preserves partial-bijection invariants");
        assert_eq!(c.right_of(left), expected);
    }

    #[rstest]
    #[case::matched_first(n(2), Some(n(0)))]
    #[case::matched_second(n(3), Some(n(1)))]
    #[case::unmatched(n(0), None)]
    fn test_correspondence_left_of(#[case] right: NodeId, #[case] expected: Option<NodeId>) {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4)
            .expect("correspondence producer preserves partial-bijection invariants");
        assert_eq!(c.left_of(right), expected);
    }

    #[rstest]
    fn test_correspondence_left_unmatched() {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4)
            .expect("correspondence producer preserves partial-bijection invariants");
        assert_eq!(c.left_unmatched(), vec![n(2)]);
    }

    #[rstest]
    fn test_correspondence_right_unmatched() {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4)
            .expect("correspondence producer preserves partial-bijection invariants");
        assert_eq!(c.right_unmatched(), vec![n(0), n(1)]);
    }

    #[rstest]
    fn test_correspondence_edge_matched_pairs(paths: (Graph, Graph, Correspondence<NodeId>)) {
        let (left, right, c) = paths;
        assert_eq!(
            c.edge_matched_pairs(&left, &right),
            Some(vec![(e(0), e(1)), (e(1), e(2))]),
        );
    }

    #[rstest]
    #[case::left_count(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        Correspondence::new(vec![(n(0), n(0))], 1, 2)
            .expect("correspondence producer preserves partial-bijection invariants"),
    )]
    #[case::right_count(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        Correspondence::new(vec![(n(0), n(0))], 2, 1)
            .expect("correspondence producer preserves partial-bijection invariants"),
    )]
    #[case::parallel_left(
        Graph::new(2, &[[0, 1], [0, 1]]),
        Graph::new(2, &[[0, 1]]),
        Correspondence::from_images(&[n(0), n(1)], 2),
    )]
    #[case::parallel_right(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1], [0, 1]]),
        Correspondence::from_images(&[n(0), n(1)], 2),
    )]
    fn test_correspondence_edge_matched_pairs_error(
        #[case] left: Graph,
        #[case] right: Graph,
        #[case] correspondence: Correspondence<NodeId>,
    ) {
        assert_eq!(correspondence.edge_matched_pairs(&left, &right), None);
    }

    #[rstest]
    fn test_correspondence_shared_edge_count(paths: (Graph, Graph, Correspondence<NodeId>)) {
        let (left, right, c) = paths;
        assert_eq!(c.shared_edge_count(&left, &right), Some(2));
    }

    #[rstest]
    #[case::count(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1]]),
        Correspondence::new(vec![(n(0), n(0))], 1, 2)
            .expect("correspondence producer preserves partial-bijection invariants"),
    )]
    #[case::parallel(
        Graph::new(2, &[[0, 1]]),
        Graph::new(2, &[[0, 1], [0, 1]]),
        Correspondence::from_images(&[n(0), n(1)], 2),
    )]
    fn test_correspondence_shared_edge_count_error(
        #[case] left: Graph,
        #[case] right: Graph,
        #[case] correspondence: Correspondence<NodeId>,
    ) {
        assert_eq!(correspondence.shared_edge_count(&left, &right), None);
    }

    #[rstest]
    #[case::total(vec![(n(0), n(0)), (n(1), n(1))], 2, 2, (true, true, true))]
    #[case::left_unmatched(vec![(n(0), n(0))], 2, 1, (false, true, false))]
    #[case::right_unmatched(vec![(n(0), n(0))], 1, 2, (true, false, false))]
    #[case::both_unmatched(vec![(n(0), n(0))], 2, 2, (false, false, false))]
    fn test_correspondence_is_total(
        #[case] matched_pairs: Vec<(NodeId, NodeId)>,
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] expected: (bool, bool, bool),
    ) {
        let correspondence = Correspondence::new(matched_pairs, left_count, right_count)
            .expect("correspondence producer preserves partial-bijection invariants");
        assert_eq!(
            (
                correspondence.is_total_on_left(),
                correspondence.is_total_on_right(),
                correspondence.is_total(),
            ),
            expected,
        );
    }

    #[rstest]
    #[case::same_intermediate(
        Correspondence::new(
            vec![
                (NodeId(0), NodeId(10)),
                (NodeId(1), NodeId(11)),
                (NodeId(2), NodeId(12)),
            ],
            3,
            13,
        ).expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence::new(
            vec![
                (NodeId(10), NodeId(100)),
                (NodeId(11), NodeId(101)),
            ],
            13,
            102,
        ).expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence::new(
            vec![
                (NodeId(0), NodeId(100)),
                (NodeId(1), NodeId(101)),
            ],
            3,
            102,
        ).expect("correspondence producer preserves partial-bijection invariants"),
    )]
    #[case::deletion_then_addition(
        Correspondence::new(vec![], 1, 0).expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence::new(vec![], 0, 1).expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence::new(vec![], 1, 1).expect("correspondence producer preserves partial-bijection invariants"),
    )]
    fn test_correspondence_compose(
        #[case] left: Correspondence<NodeId>,
        #[case] right: Correspondence<NodeId>,
        #[case] expected: Correspondence<NodeId>,
    ) {
        assert_eq!(left.compose(&right), Ok(expected));
    }

    #[rstest]
    #[case::shorter(3, 1)]
    #[case::longer(1, 3)]
    #[case::empty(0, 1)]
    fn test_correspondence_compose_error(
        #[case] right_count: usize,
        #[case] next_left_count: usize,
    ) {
        let left = Correspondence::<NodeId>::new(vec![], 2, right_count).unwrap();
        let right = Correspondence::new(vec![], next_left_count, 2).unwrap();
        let expected = CorrespondenceComposeError {
            right_count,
            next_left_count,
        };
        assert_eq!(left.compose(&right), Err(expected.clone()));
        assert_eq!(
            Correspondence::compose_all([left.reverse(), left, right]),
            Err(expected)
        );
    }

    #[rstest]
    #[case::empty(vec![], None)]
    #[case::singleton(
        vec![Correspondence::new(
            vec![(NodeId(0), NodeId(1))],
            1,
            2,
        ).expect("correspondence producer preserves partial-bijection invariants")],
        Some(Correspondence::new(
            vec![(NodeId(0), NodeId(1))],
            1,
            2,
        ).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    #[case::multiple(
        vec![
            Correspondence::new(
                vec![
                    (NodeId(0), NodeId(0)),
                    (NodeId(1), NodeId(2)),
                ],
                2,
                3,
            ).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(
                vec![
                    (NodeId(0), NodeId(1)),
                    (NodeId(2), NodeId(0)),
                ],
                3,
                2,
            ).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(
                vec![(NodeId(0), NodeId(2))],
                2,
                3,
            ).expect("correspondence producer preserves partial-bijection invariants"),
        ],
        Some(Correspondence::new(
            vec![(NodeId(1), NodeId(2))],
            2,
            3,
        ).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    fn test_correspondence_compose_all(
        #[case] correspondences: Vec<Correspondence<NodeId>>,
        #[case] expected: Option<Correspondence<NodeId>>,
    ) {
        assert_eq!(Correspondence::compose_all(correspondences), Ok(expected));
    }

    #[rstest]
    fn test_correspondence_reverse() {
        // pairs and counts swap; the left-unmatched id 2 becomes right-unmatched.
        let c = Correspondence::new(vec![(n(0), n(3)), (n(1), n(1))], 3, 4)
            .expect("correspondence producer preserves partial-bijection invariants");
        let reversed = c.reverse();
        assert_eq!(reversed.matched_pairs(), &[(n(1), n(1)), (n(3), n(0))]);
        assert_eq!(reversed.left_unmatched(), vec![n(0), n(2)]);
        assert_eq!(reversed.right_unmatched(), vec![n(2)]);
    }
}
