//! Stereo ligand AST: a ligand occupying a coordination position of a stereo site.

use umol_graph_core::{NodeId, ParticipantRefs, RelationParticipant, RemovalRemapping};

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
    fn remap(self, remapping: &RemovalRemapping) -> Option<Self> {
        Some(Self {
            atom_id: AtomId::from(remapping.map_node(NodeId::from(self.atom_id))?),
            kind: self.kind,
        })
    }

    fn unmap(self, remapping: &RemovalRemapping) -> Self {
        Self {
            atom_id: AtomId::from(remapping.unmap_node(NodeId::from(self.atom_id))),
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
    fn test_stereo_ligand_remap() {
        // node 1 removed ⇒ surviving node 3 densifies to 2; the kind is carried
        let remapping = RemovalRemapping::new(vec![1], Vec::new());
        let ligand = StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen);
        assert_eq!(
            ligand.remap(&remapping),
            Some(StereoLigand::new(
                AtomId(2),
                StereoLigandKind::ImplicitHydrogen
            )),
        );
    }

    #[rstest]
    fn test_stereo_ligand_remap_removed() {
        let remapping = RemovalRemapping::new(vec![3], Vec::new());
        let ligand = StereoLigand::new(AtomId(3), StereoLigandKind::Atom);
        assert_eq!(ligand.remap(&remapping), None);
    }

    #[rstest]
    fn test_stereo_ligand_unmap() {
        let remapping = RemovalRemapping::new(vec![1], Vec::new());
        let ligand = StereoLigand::new(AtomId(2), StereoLigandKind::Atom);
        assert_eq!(
            ligand.unmap(&remapping),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        );
    }
}
