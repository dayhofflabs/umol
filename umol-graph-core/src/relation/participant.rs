//! Participant values, factor-local positions, and graph-reference transport.
//!
//! [RelationParticipant] separates a participant's complete value from its node/edge references.
//! [ParticipantRefs] routes incidence; [ParticipantPosition] addresses a stored factor position.

use std::hash::Hash;

use crate::compact::GraphCompaction;
use crate::correspondence::GraphCorrespondence;
use crate::graph::{EdgeId, NodeId};
use crate::remap::GraphRemapping;

/// Position within one relation factor, independent of node, edge, and relation ids.
///
/// Construction does not validate a factor length; the consuming operation checks the position.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ParticipantPosition(pub u32);

impl ParticipantPosition {
    /// Return the factor-local position as usize.
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Node and edge references exposed by a participant for incidence and transport coverage.
///
/// Each field is independent: a participant may report a node, an edge, both, or neither.
/// Equal references do not imply equal complete participant values. Graph membership is contextual.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ParticipantRefs {
    /// Referenced node, if any.
    pub node: Option<NodeId>,
    /// Referenced edge, if any.
    pub edge: Option<EdgeId>,
}

/// A typed relation participant with node/edge references and reference transport.
///
/// Equality, ordering, and hashing describe the complete value. Incidence uses only
/// [`Self::refs`], whose result must be stable for a given value.
///
/// # Semantic properties
///
/// Transport preserves all data other than the referenced ids. Every transported id must
/// be reported by [`Self::refs`]. A covering identity mapping preserves the value; successive
/// covered mappings agree with their composition. A covering remapping followed by its inverse
/// recovers the original value. Compaction drops a participant if any required reference has
/// no image; expanding a surviving compacted participant recovers its original value.
///
/// Relation-set transport relies on these laws. Public-API properties in
/// `tests/property/relation.rs` exercise mapping and remapping with node/edge participants;
/// `tests/property/restore.rs` checks inverse transport and combined removal/restoration sequences.
/// Structured participants exercise incidence and preservation of non-reference data.
pub trait RelationParticipant: Copy + Ord + Hash {
    /// Relabel every referenced id through a correspondence, preserving other participant data.
    ///
    /// # Panics
    /// Panics when a referenced node or edge has no image.
    fn map(self, correspondence: &GraphCorrespondence) -> Self {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Relabel every referenced id, or return `None` if any reference has no image.
    /// Unused correspondence entries need not be matched. Every id looked up must be
    /// reported by [`refs`](Self::refs); all other participant data must be preserved.
    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self>;

    /// Relabel this participant through `remapping`.
    ///
    /// Every node or edge id read from `remapping` must be reported by [`refs`](Self::refs), so
    /// checked relation-set remapping can establish coverage before calling this method.
    ///
    /// # Panics
    ///
    /// Panics if a reference is outside the corresponding source domain.
    fn remap(self, remapping: &GraphRemapping) -> Self;

    /// Compact every reference, returning `None` if any reference has no surviving image.
    ///
    /// References outside the source domain also have no image.
    fn compact(self, compaction: &GraphCompaction) -> Option<Self>;

    /// Expand references from a compaction's result space to its source space.
    ///
    /// # Panics
    ///
    /// Panics if a reference is outside the corresponding result domain.
    fn uncompact(self, compaction: &GraphCompaction) -> Self;

    /// Return every graph id used to represent this participant.
    fn refs(self) -> ParticipantRefs;
}

impl RelationParticipant for NodeId {
    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self> {
        correspondence.nodes().right_of(self)
    }

    fn remap(self, remapping: &GraphRemapping) -> Self {
        remapping.map_node(self)
    }

    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        compaction.compact_node(self)
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
        compaction.uncompact_node(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: Some(self),
            edge: None,
        }
    }
}

impl RelationParticipant for EdgeId {
    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self> {
        correspondence.edges().right_of(self)
    }

    fn remap(self, remapping: &GraphRemapping) -> Self {
        remapping.map_edge(self)
    }

    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        compaction.compact_edge(self)
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
        compaction.uncompact_edge(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: None,
            edge: Some(self),
        }
    }
}
