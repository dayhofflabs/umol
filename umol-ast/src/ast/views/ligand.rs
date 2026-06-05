//! Stereo ligand views.

use super::super::ids::AtomId;
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::molecule::MoleculeAst;
use super::atom::AtomView;

/// Borrowed view of a stereo ligand: its kind, bearing atom, and parent
/// molecule for resolving the atom view.
#[derive(Clone, Copy, Debug)]
pub struct StereoLigandView<'a> {
    ligand: StereoLigand,
    molecule: &'a MoleculeAst,
}

impl<'a> StereoLigandView<'a> {
    pub(crate) fn new(ligand: StereoLigand, molecule: &'a MoleculeAst) -> Self {
        Self { ligand, molecule }
    }

    pub fn kind(&self) -> StereoLigandKind {
        self.ligand.kind
    }

    pub fn atom_id(&self) -> AtomId {
        self.ligand.atom_id
    }

    pub fn atom(&self) -> AtomView<'a> {
        self.molecule.atom(self.ligand.atom_id)
    }
}
