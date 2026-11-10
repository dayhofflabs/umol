//! Stereo ligand representation: a ligand occupying a coordination position of a stereo site.

use umol_graph_core::{
    GraphCompaction, GraphCorrespondence, GraphRemapping, NodeId, ParticipantRefs,
    RelationParticipant,
};

use super::id::AtomId;

/// Stereo ligand kind: an atom, implicit hydrogen, lone pair, or generic virtual ligand.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StereoLigandKind {
    Atom,
    ImplicitHydrogen,
    LonePair,
    /// A virtual substituent whose chemical identity has not been specified.
    Virtual,
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

    /// Map the atom reference, preserving the ligand kind, including virtual ligands.
    ///
    /// # Panics
    /// Panics if the atom reference has no image in `correspondence`.
    pub fn map(self, correspondence: &GraphCorrespondence) -> Self {
        RelationParticipant::map(self, correspondence)
    }

    /// Map the atom reference, or return `None` if it has no image.
    /// Unreferenced nodes and edges need not have images.
    pub fn try_map(self, correspondence: &GraphCorrespondence) -> Option<Self> {
        RelationParticipant::try_map(self, correspondence)
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
    use umol_graph_core::{Compaction, Correspondence, Remapping};

    use super::*;

    #[rstest]
    #[case::atom(StereoLigandKind::Atom)]
    #[case::hydrogen(StereoLigandKind::ImplicitHydrogen)]
    #[case::lone_pair(StereoLigandKind::LonePair)]
    #[case::virtual_ligand(StereoLigandKind::Virtual)]
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
        assert_eq!(ligand.try_map(&correspondence), Some(expected));
        assert_eq!(ligand.map(&correspondence), expected);
    }

    #[rstest]
    #[case::atom(StereoLigandKind::Atom)]
    #[case::hydrogen(StereoLigandKind::ImplicitHydrogen)]
    #[case::lone_pair(StereoLigandKind::LonePair)]
    #[case::virtual_ligand(StereoLigandKind::Virtual)]
    fn test_stereo_ligand_try_map_error(#[case] kind: StereoLigandKind) {
        let correspondence = GraphCorrespondence::new(
            Correspondence::new(vec![], 4, 6).unwrap(),
            Correspondence::empty(),
        );
        assert_eq!(
            RelationParticipant::try_map(StereoLigand::new(AtomId(2), kind), &correspondence),
            None
        );
        assert_eq!(
            StereoLigand::new(AtomId(2), kind).try_map(&correspondence),
            None
        );
    }

    #[rstest]
    #[case::atom(StereoLigandKind::Atom)]
    #[case::hydrogen(StereoLigandKind::ImplicitHydrogen)]
    #[case::lone_pair(StereoLigandKind::LonePair)]
    #[case::virtual_ligand(StereoLigandKind::Virtual)]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_stereo_ligand_map_error(#[case] kind: StereoLigandKind) {
        StereoLigand::new(AtomId(2), kind).map(&GraphCorrespondence::new(
            Correspondence::empty(),
            Correspondence::empty(),
        ));
    }

    #[rstest]
    fn test_stereo_ligand_new() {
        let ligand = StereoLigand::new(AtomId(3), StereoLigandKind::Atom);
        assert_eq!(ligand.atom_id, AtomId(3));
        assert_eq!(ligand.kind, StereoLigandKind::Atom);
    }

    #[rstest]
    #[case::atom(StereoLigandKind::Atom, false)]
    #[case::hydrogen(StereoLigandKind::ImplicitHydrogen, false)]
    #[case::lone_pair(StereoLigandKind::LonePair, false)]
    #[case::virtual_ligand(StereoLigandKind::Virtual, true)]
    fn test_stereo_ligand_eq(#[case] kind: StereoLigandKind, #[case] expected: bool) {
        assert_eq!(
            StereoLigand::new(AtomId(0), StereoLigandKind::Virtual)
                == StereoLigand::new(AtomId(0), kind),
            expected,
        );
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
    #[case::atom(StereoLigandKind::Atom)]
    #[case::hydrogen(StereoLigandKind::ImplicitHydrogen)]
    #[case::lone_pair(StereoLigandKind::LonePair)]
    #[case::virtual_ligand(StereoLigandKind::Virtual)]
    fn test_stereo_ligand_compact(#[case] kind: StereoLigandKind) {
        let compaction = GraphCompaction::new(
            Compaction::new(4, vec![NodeId(1)]).unwrap(),
            Compaction::empty(),
        );
        let ligand = StereoLigand::new(AtomId(3), kind);
        assert_eq!(
            ligand.compact(&compaction),
            Some(StereoLigand::new(AtomId(2), kind)),
        );
    }

    #[rstest]
    fn test_stereo_ligand_compact_removed() {
        let compaction = GraphCompaction::new(
            Compaction::new(4, vec![NodeId(3)]).unwrap(),
            Compaction::empty(),
        );
        let ligand = StereoLigand::new(AtomId(3), StereoLigandKind::Atom);
        assert_eq!(ligand.compact(&compaction), None);
    }

    #[rstest]
    #[case::atom(StereoLigandKind::Atom)]
    #[case::hydrogen(StereoLigandKind::ImplicitHydrogen)]
    #[case::lone_pair(StereoLigandKind::LonePair)]
    #[case::virtual_ligand(StereoLigandKind::Virtual)]
    fn test_stereo_ligand_uncompact(#[case] kind: StereoLigandKind) {
        let compaction = GraphCompaction::new(
            Compaction::new(4, vec![NodeId(1)]).unwrap(),
            Compaction::empty(),
        );
        let ligand = StereoLigand::new(AtomId(2), kind);
        assert_eq!(
            ligand.uncompact(&compaction),
            StereoLigand::new(AtomId(3), kind),
        );
    }

    #[rstest]
    #[case::atom(StereoLigandKind::Atom)]
    #[case::hydrogen(StereoLigandKind::ImplicitHydrogen)]
    #[case::lone_pair(StereoLigandKind::LonePair)]
    #[case::virtual_ligand(StereoLigandKind::Virtual)]
    fn test_stereo_ligand_remap(#[case] kind: StereoLigandKind) {
        let remapping = GraphRemapping::new(
            Remapping::new(vec![NodeId(2), NodeId(0), NodeId(1)]).unwrap(),
            Remapping::empty(),
        );
        assert_eq!(
            StereoLigand::new(AtomId(1), kind).remap(&remapping),
            StereoLigand::new(AtomId(0), kind),
        );
    }
}
