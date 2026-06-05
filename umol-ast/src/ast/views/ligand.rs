//! Stereo ligand views.

use super::super::ids::AtomId;
use super::super::ligand::StereoLigandKind;
use super::super::molecule::MoleculeAst;
use super::atom::AtomView;

/// Borrowed view of a stereo ligand: its kind, bearing atom, and parent
/// molecule for resolving the atom view.
#[derive(Clone, Copy, Debug)]
pub struct StereoLigandView<'a> {
    kind: StereoLigandKind,
    atom_id: AtomId,
    molecule: &'a MoleculeAst,
}

impl<'a> StereoLigandView<'a> {
    pub(crate) fn new(kind: StereoLigandKind, atom_id: AtomId, molecule: &'a MoleculeAst) -> Self {
        Self {
            kind,
            atom_id,
            molecule,
        }
    }

    pub fn kind(&self) -> StereoLigandKind {
        self.kind
    }

    pub fn atom_id(&self) -> AtomId {
        self.atom_id
    }

    pub fn atom(&self) -> AtomView<'a> {
        self.molecule.atom(self.atom_id)
    }
}
