//! A partial bijection between two id spaces.
//!
//! Framed as a matching in the bipartite graph over the two id sets: an id is **matched** (paired
//! with a partner on the other side — the shared interface) or **unmatched** (present on only one
//! side). Generic in the id type, so the same carrier serves node correspondences (atoms) and every
//! entity family (bonds, overlays) one layer up; `Correspondence<NodeId>` additionally exposes the
//! induced edge correspondence over the two graphs.

use std::collections::{BTreeSet, HashMap};
use std::error::Error;
use std::fmt::{self, Debug, Display, Formatter};

use crate::graph::{EdgeId, Graph, NodeId, Remapping};

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

impl<Id: Copy + Ord + From<usize>> Correspondence<Id> {
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

    /// Relational composition: `self` (left↔middle) followed by `other` (middle↔right), yielding a
    /// left↔right correspondence. A left id matched to a middle id that `other` leaves unmatched
    /// becomes unmatched.
    ///
    /// Composition matches numerical middle ids even when the declared intermediate counts differ.
    /// Ids outside the shorter space have no pair and therefore behave as absent or unmatched.
    pub fn compose(&self, other: &Correspondence<Id>) -> Correspondence<Id> {
        let matched_pairs = self
            .matched_pairs
            .iter()
            .filter_map(|&(left, middle)| other.right_of(middle).map(|right| (left, right)))
            .collect();
        Correspondence::new(matched_pairs, self.left_count, other.right_count)
            .unwrap_or_else(|_| panic!("composition of valid correspondences is valid"))
    }

    /// Compose correspondences in iteration order. Returns `None` for an empty input and the value
    /// itself for a singleton.
    pub fn compose_all(correspondences: impl IntoIterator<Item = Self>) -> Option<Self> {
        correspondences
            .into_iter()
            .reduce(|left, right| left.compose(&right))
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

/// A subgraph↔host correspondence over a `Graph`: its node and edge families. The graph-core base
/// that the molecule-level `MoleculeCorrespondence` (atoms + bonds + overlays) extends — produced by
/// induced subgraphs, subiso matches, and common-subgraph search. The objective of each is a family
/// size: `nodes().matched_pair_count()` (induced / MCIS), `edges().matched_pair_count()` (MCES).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GraphCorrespondence {
    nodes: Correspondence<NodeId>,
    edges: Correspondence<EdgeId>,
}

impl GraphCorrespondence {
    pub fn new(nodes: Correspondence<NodeId>, edges: Correspondence<EdgeId>) -> Self {
        Self { nodes, edges }
    }

    /// The graph correspondence a node correspondence induces over `left` / `right`: its edge family
    /// is the induced edge correspondence (an edge matched when both endpoints are matched). Exact
    /// whenever every structurally matched edge is admissible — a subiso match, a common *induced*
    /// subgraph, or a common-subgraph walk whose adjacency enforces the edge predicate; not for an
    /// edge-subgraph result whose producer selects the edge family directly. Returns `None` when the
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

    /// Relational composition of the node and edge correspondences.
    pub fn compose(&self, other: &GraphCorrespondence) -> GraphCorrespondence {
        GraphCorrespondence::new(
            self.nodes.compose(&other.nodes),
            self.edges.compose(&other.edges),
        )
    }

    /// Compose graph correspondences in iteration order. Returns `None` for an empty input and the
    /// value itself for a singleton.
    pub fn compose_all(correspondences: impl IntoIterator<Item = Self>) -> Option<Self> {
        correspondences
            .into_iter()
            .reduce(|left, right| left.compose(&right))
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

    /// This correspondence as a [`Remapping`] — a dense old→new relabel of both id spaces. Requires
    /// it be **total on the left** (every left id matched), as a pushout's coprojection is: left id `i`
    /// maps to its partner.
    pub fn to_remapping(&self) -> Option<Remapping> {
        if !self.is_total_on_left() {
            return None;
        }
        Some(Remapping::new(
            self.nodes
                .matched_pairs()
                .iter()
                .map(|&(_, right)| right)
                .collect(),
            self.edges
                .matched_pairs()
                .iter()
                .map(|&(_, right)| right)
                .collect(),
        ))
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    fn e(i: u32) -> EdgeId {
        EdgeId(i)
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
            left.compose(&right),
            GraphCorrespondence::new(
                Correspondence::new(vec![(NodeId(0), NodeId(2))], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(EdgeId(0), EdgeId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ),
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
            GraphCorrespondence::compose_all(graph_correspondences.into_iter().take(count)),
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
    #[case::mixed_families(
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
    fn test_graph_correspondence_to_remapping() {
        // total-on-left node map 0→2, 1→0, 2→1 and edge map 0→1, 1→0.
        let c = GraphCorrespondence::new(
            Correspondence::from_images(&[n(2), n(0), n(1)], 3),
            Correspondence::from_images(&[e(1), e(0)], 2),
        );
        assert_eq!(
            c.to_remapping(),
            Some(Remapping::new(vec![n(2), n(0), n(1)], vec![e(1), e(0)],)),
        );
    }

    #[rstest]
    #[case::nodes(GraphCorrespondence::new(
        Correspondence::new(vec![(n(0), n(0))], 2, 1)
            .expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence::new(Vec::new(), 0, 0)
            .expect("correspondence producer preserves partial-bijection invariants"),
    ))]
    #[case::edges(GraphCorrespondence::new(
        Correspondence::new(Vec::new(), 0, 0)
            .expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence::new(vec![(e(0), e(0))], 2, 1)
            .expect("correspondence producer preserves partial-bijection invariants"),
    ))]
    fn test_graph_correspondence_to_remapping_partial(#[case] correspondence: GraphCorrespondence) {
        assert_eq!(correspondence.to_remapping(), None);
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
    #[case::mismatched_intermediate(
        Correspondence::new(
            vec![
                (NodeId(0), NodeId(0)),
                (NodeId(1), NodeId(2)),
            ],
            2,
            3,
        ).expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence::new(
            vec![(NodeId(0), NodeId(1))],
            1,
            2,
        ).expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence::new(
            vec![(NodeId(0), NodeId(1))],
            2,
            2,
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
        assert_eq!(left.compose(&right), expected);
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
                1,
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
        assert_eq!(Correspondence::compose_all(correspondences), expected);
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
