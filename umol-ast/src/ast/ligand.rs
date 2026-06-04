//! Stereo ligand AST: a ligand occupying a coordination position of a stereo site.

use umol_graph_core::{NodeId, ParticipantRefs, RelationParticipant, Remapping};

use super::ids::AtomId;

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
    node: NodeId,
    kind: StereoLigandKind,
}

impl StereoLigand {
    pub fn new(atom: AtomId, kind: StereoLigandKind) -> Self {
        Self {
            node: atom.into(),
            kind,
        }
    }

    pub fn atom(self) -> AtomId {
        self.node.into()
    }

    pub fn kind(self) -> StereoLigandKind {
        self.kind
    }
}

impl RelationParticipant for StereoLigand {
    fn remap(self, remapping: &Remapping) -> Option<Self> {
        Some(Self {
            node: remapping.map_node(self.node)?,
            kind: self.kind,
        })
    }

    fn unmap(self, remapping: &Remapping) -> Self {
        Self {
            node: remapping.unmap_node(self.node),
            kind: self.kind,
        }
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: Some(self.node),
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
        assert_eq!(ligand.atom(), AtomId(3));
        assert_eq!(ligand.kind(), StereoLigandKind::Atom);
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
        let remapping = Remapping::new(vec![1], Vec::new());
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
        let remapping = Remapping::new(vec![3], Vec::new());
        let ligand = StereoLigand::new(AtomId(3), StereoLigandKind::Atom);
        assert_eq!(ligand.remap(&remapping), None);
    }

    #[rstest]
    fn test_stereo_ligand_unmap() {
        let remapping = Remapping::new(vec![1], Vec::new());
        let ligand = StereoLigand::new(AtomId(2), StereoLigandKind::Atom);
        assert_eq!(
            ligand.unmap(&remapping),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        );
    }
}
