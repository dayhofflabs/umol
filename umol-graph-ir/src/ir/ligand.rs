//! Stereo ligand representation: a ligand occupying a coordination position of a stereo site.

use umol_graph_core::{
    GraphCompaction, GraphCorrespondence, GraphRemapping, NodeId, ParticipantRefs,
    RelationParticipant,
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
    fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self> {
        Some(Self {
            atom_id: AtomId::from(
                correspondence
                    .nodes()
                    .right_of(NodeId::from(self.atom_id))?,
            ),
            kind: self.kind,
        })
    }
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

    fn remap(self, remapping: &GraphRemapping) -> Self {
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
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Correspondence;

    use super::*;

    #[rstest]
    #[case::atom(StereoLigandKind::Atom)]
    #[case::hydrogen(StereoLigandKind::ImplicitHydrogen)]
    #[case::lone_pair(StereoLigandKind::LonePair)]
    fn test_stereo_ligand_map(#[case] kind: StereoLigandKind) {
        let correspondence = GraphCorrespondence::new(
            Correspondence::new(vec![(NodeId(2), NodeId(5))], 4, 6).unwrap(),
            Correspondence::empty(),
        );
        let ligand = StereoLigand::new(AtomId(2), kind);
        let expected = StereoLigand::new(AtomId(5), kind);
        assert_eq!(
            RelationParticipant::try_map(ligand, &correspondence),
            Some(expected)
        );
        assert_eq!(RelationParticipant::map(ligand, &correspondence), expected);
    }

    #[rstest]
    #[case::atom(StereoLigandKind::Atom)]
    #[case::hydrogen(StereoLigandKind::ImplicitHydrogen)]
    #[case::lone_pair(StereoLigandKind::LonePair)]
    fn test_stereo_ligand_try_map_error(#[case] kind: StereoLigandKind) {
        let correspondence = GraphCorrespondence::new(
            Correspondence::new(vec![], 4, 6).unwrap(),
            Correspondence::empty(),
        );
        assert_eq!(
            RelationParticipant::try_map(StereoLigand::new(AtomId(2), kind), &correspondence),
            None
        );
    }

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
