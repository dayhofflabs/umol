//! Stereo ligand representation: a ligand occupying a coordination position of a stereo site.

use umol_graph_core::{
    GraphCompaction, NodeId, ParticipantAnchor, ParticipantRefs, RelationParticipant, Remapping,
};

use super::id::AtomId;

/// Stereo ligand kind: atom, or virtual ligand (implicit hydrogen or lone pair).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StereoLigandKind {
    Atom,
    ImplicitHydrogen,
    LonePair,
}

/// Stereo ligand occupying a coordination position of a stereo site.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct StereoLigand {
    pub atom_id: AtomId,
    pub kind: StereoLigandKind,
}

impl StereoLigand {
    pub fn new(atom_id: AtomId, kind: StereoLigandKind) -> Self {
        Self { atom_id, kind }
    }
}

impl RelationParticipant for StereoLigand {
    fn compact(self, compaction: &GraphCompaction) -> Option<Self> {
        Some(Self {
            atom_id: AtomId::from(compaction.compact_node(NodeId::from(self.atom_id))?),
            kind: self.kind,
        })
    }

    fn uncompact(self, compaction: &GraphCompaction) -> Self {
        Self {
            atom_id: AtomId::from(compaction.uncompact_node(NodeId::from(self.atom_id))),
            kind: self.kind,
        }
    }

    fn remap(self, remapping: &Remapping) -> Self {
        Self {
            atom_id: AtomId::from(remapping.map_node(NodeId::from(self.atom_id))),
            kind: self.kind,
        }
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: Some(NodeId::from(self.atom_id)),
            edge: None,
        }
    }

    fn anchor(self) -> Option<ParticipantAnchor> {
        Some(ParticipantAnchor::Node(NodeId::from(self.atom_id)))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    fn test_stereo_ligand_new() {
        let ligand = StereoLigand::new(AtomId(3), StereoLigandKind::Atom);
        assert_eq!(ligand.atom_id, AtomId(3));
        assert_eq!(ligand.kind, StereoLigandKind::Atom);
    }

    #[rstest]
    fn test_stereo_ligand_refs() {
        let ligand = StereoLigand::new(AtomId(2), StereoLigandKind::LonePair);
        assert_eq!(
            ligand.refs(),
            ParticipantRefs {
                node: Some(NodeId(2)),
                edge: None,
            }
        );
    }

    #[rstest]
    fn test_stereo_ligand_compact() {
        // node 1 removed ⇒ surviving node 3 densifies to 2; the kind is carried
        let compaction = GraphCompaction::new(vec![NodeId(1)], Vec::new());
        let ligand = StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen);
        assert_eq!(
            ligand.compact(&compaction),
            Some(StereoLigand::new(
                AtomId(2),
                StereoLigandKind::ImplicitHydrogen
            )),
        );
    }

    #[rstest]
    fn test_stereo_ligand_compact_removed() {
        let compaction = GraphCompaction::new(vec![NodeId(3)], Vec::new());
        let ligand = StereoLigand::new(AtomId(3), StereoLigandKind::Atom);
        assert_eq!(ligand.compact(&compaction), None);
    }

    #[rstest]
    fn test_stereo_ligand_uncompact() {
        let compaction = GraphCompaction::new(vec![NodeId(1)], Vec::new());
        let ligand = StereoLigand::new(AtomId(2), StereoLigandKind::Atom);
        assert_eq!(
            ligand.uncompact(&compaction),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        );
    }
}
