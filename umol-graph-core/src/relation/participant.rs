use std::hash::Hash;

use crate::compact::GraphCompaction;
use crate::correspondence::GraphCorrespondence;
use crate::graph::{EdgeId, NodeId};
use crate::remap::GraphRemapping;

/// Position of a participant within a single relation's tuple — local to the
/// relation (frame-relative), distinct from the global `NodeId`/`EdgeId`/`RelationId`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ParticipantPosition(pub u32);

impl ParticipantPosition {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// The id-space contents of a participant, surfaced for the incidence index.
/// At most one ref per space today (a node or an edge); a future port type
/// could fill both.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ParticipantRefs {
    pub node: Option<NodeId>,
    pub edge: Option<EdgeId>,
}
/// A value that can occupy a relation factor: supports compaction, remapping, and
/// correspondence-based id transport, and exposes its node/edge refs for incidence. One impl per concrete
/// id type — dispatch is static, since a factor is homogeneous.
pub trait RelationParticipant: Copy + Ord + Hash {
    fn compact(self, compaction: &GraphCompaction) -> Option<Self>;
    fn uncompact(self, compaction: &GraphCompaction) -> Self;

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
    fn remap(self, remapping: &GraphRemapping) -> Self;

    /// Return every graph id used to represent this participant.
    fn refs(self) -> ParticipantRefs;
}

impl RelationParticipant for NodeId {
    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self> {
        correspondence.nodes().right_of(self)
    }
    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        compaction.compact_node(self)
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
        compaction.uncompact_node(self)
    }

    fn remap(self, remapping: &GraphRemapping) -> Self {
        remapping.map_node(self)
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
    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        compaction.compact_edge(self)
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
        compaction.uncompact_edge(self)
    }

    fn remap(self, remapping: &GraphRemapping) -> Self {
        remapping.map_edge(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: None,
            edge: Some(self),
        }
    }
}
