//! Journal-free batch application for the molecule editor.

use umol_graph_core::Correspondence;

use super::MoleculeEditor;
use crate::ir::correspondence::MoleculeCorrespondence;
use crate::ir::edit::Edits;
use crate::ir::error::MoleculeApplyError;
use crate::ir::molecule::apply::ApplicationState;

impl MoleculeEditor {
    /// Apply an ordered [`Edits`] batch without constructing an undo journal.
    ///
    /// The editor is consumed so that a failed batch cannot expose partially applied state. On
    /// success, the returned editor remains transient. MoleculeEditor::probe and
    /// MoleculeEditor::finish check molecule integrity before publishing access or ownership.
    ///
    /// # Errors
    ///
    /// Returns MoleculeApplyError::Transaction when an edit handle, precondition, or shape is
    /// invalid for the evolving editor state. No editor is returned on failure, including any
    /// earlier direct or batch changes.
    ///
    /// # Semantic properties
    ///
    /// Applying the same batch with or without undo recording produces the same molecule on
    /// successful publication. On failure, no intermediate editor state is returned.
    pub fn apply(mut self, edits: Edits) -> Result<Self, MoleculeApplyError> {
        let mut state = ApplicationState::new(&self.molecule);
        for edit in edits {
            self.molecule.apply_edit(edit, &mut state)?;
        }
        Ok(self)
    }

    /// Consume the editor and apply a batch, returning the resulting editor and batch correspondence.
    ///
    /// The witness starts at this batch's input, not the editor's session origin. The resulting
    /// editor remains transient; publication performs the integrity check.
    ///
    /// # Errors
    ///
    /// Returns the same transaction error as [`Self::apply`], consuming the editor on failure.
    ///
    /// # Semantic properties
    ///
    /// Discarding the witness gives the same editor as the plain operation. Earlier direct changes
    /// and batches are outside this correspondence; its source is the state at this call.
    pub fn tracked_apply(
        mut self,
        edits: Edits,
    ) -> Result<(Self, MoleculeCorrespondence), MoleculeApplyError> {
        let mut correspondence = MoleculeCorrespondence::new(
            Correspondence::identity(self.atom_count()),
            Correspondence::identity(self.bond_count()),
            Correspondence::identity(self.dative_bond_count()),
            Correspondence::identity(self.aromatic_system_count()),
            Correspondence::identity(self.multicenter_bond_count()),
            Correspondence::identity(self.noncovalent_bond_count()),
            Correspondence::identity(self.stereo_atom_count()),
            Correspondence::identity(self.stereo_bond_count()),
        );
        let mut state = ApplicationState::new(&self.molecule);
        for edit in edits {
            self.molecule.apply_edit(edit, &mut state)?;
        }
        state.update_correspondence(&mut correspondence);
        Ok((self, correspondence))
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::MoleculeEditor;
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::{AtomForm, ElementForm};
    use crate::ir::bond::BondForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::edit::{AtomHandle, DativeBondHandle, Edit, Edits};
    use crate::ir::entity::EntityKind;
    use crate::ir::error::MoleculeApplyError;
    use crate::ir::id::{AtomId, BondId, DativeBondId};
    use crate::ir::ligand::{StereoLigand, StereoLigandKind};
    use crate::ir::molecule::{Molecule, TransactionError};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use crate::ir::num::NumForm;
    use crate::ir::stereo::{StereoAtomForm, StereoBondForm, StereoCoset, StereoKind};
    use crate::mol_dsl;

    #[fixture]
    fn empty() -> MoleculeEditor {
        Molecule::default().edit()
    }

    #[fixture]
    fn batched_overlays() -> MoleculeEditor {
        let mut editor = Molecule::default().edit();
        for _ in 0..6 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        for index in 0..3_u32 {
            let first = AtomId(index * 2);
            let second = AtomId(index * 2 + 1);
            let bond = editor.add_bond(first, second, BondForm::from_order(1));
            editor.add_dative_bond(&[first], second, DativeBondForm::from_order(1));
            editor.add_aromatic_system(&[first, second], AromaticSystemForm::default());
            editor.add_multicenter_bond(&[first, second], MulticenterBondForm::default());
            editor.add_noncovalent_bond(
                [first, second],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            );
            editor.add_stereo_atom(
                first,
                &[
                    StereoLigand::new(second, StereoLigandKind::Atom),
                    StereoLigand::new(first, StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(first, StereoLigandKind::LonePair),
                ],
                StereoAtomForm::default(),
            );
            editor.add_stereo_bond(
                bond,
                &[
                    StereoLigand::new(first, StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(first, StereoLigandKind::LonePair),
                    StereoLigand::new(second, StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(second, StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            );
        }
        editor
    }

    #[rstest]
    fn test_molecule_editor_apply(empty: MoleculeEditor) {
        let mut edits = Edits::new();
        edits.add_atom(AtomForm::from_element(Element::C));
        let empty = empty.apply(edits).unwrap();
        assert_eq!(empty.atom_count(), 1);
        assert_eq!(
            empty.atom(AtomId(0)).attributes().element,
            ElementForm::Lit(Element::C)
        );
    }

    #[rstest]
    fn test_molecule_editor_apply_interleaving() {
        let mut editor = mol_dsl!(r#"{:atoms ["C"]}"#).edit();
        editor.atom_mut(AtomId(0)).attributes_mut().charge = NumForm::Lit(1);
        editor.add_atom(AtomForm::from_element(Element::N));
        let mut first = Edits::new();
        let oxygen = first.add_atom(AtomForm::from_element(Element::O));
        first.add_bond(AtomHandle::Id(AtomId(0)), oxygen, BondForm::from_order(1));

        let mut editor = editor.apply(first).unwrap();
        editor.bond_mut(BondId(0)).attributes_mut().order = NumForm::Lit(2);
        let mut second = Edits::new();
        let fluorine = second.add_atom(AtomForm::from_element(Element::F));
        second.add_bond(AtomHandle::Id(AtomId(1)), fluorine, BondForm::from_order(1));

        let editor = editor.apply(second).unwrap();
        assert_eq!(
            editor.finish(),
            Ok(mol_dsl!(
                r#"{:atoms ["C#c1" "N" "O" "F"] :bonds [[0 2 "2"] [1 3 "1"]]}"#
            ))
        );
    }

    #[rstest]
    fn test_molecule_editor_apply_replace_dative_bond_acceptor(batched_overlays: MoleculeEditor) {
        let editor = batched_overlays
            .apply(Edits::from_iter([Edit::ReplaceDativeBondAcceptor {
                id: DativeBondHandle::Id(DativeBondId(0)),
                old: AtomHandle::Id(AtomId(1)),
                new: AtomHandle::Id(AtomId(3)),
            }]))
            .unwrap();
        assert_eq!(editor.dative_bond(DativeBondId(0)).acceptor_id(), AtomId(3));
    }

    #[rstest]
    fn test_molecule_editor_apply_error(empty: MoleculeEditor) {
        let mut edits = Edits::new();
        edits.remove_atom(AtomHandle::Id(AtomId(0)));
        let error = match empty.apply(edits) {
            Ok(_) => panic!("invalid edit unexpectedly applied"),
            Err(error) => error,
        };
        assert_eq!(
            error,
            MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 0,
                count: 0,
            })
        );
    }

    #[rstest]
    #[case::plain(false)]
    #[case::tracked(true)]
    fn test_molecule_editor_apply_destructive_error(
        mut empty: MoleculeEditor,
        #[case] tracked: bool,
    ) {
        empty.add_atom(AtomForm::from_element(Element::C));
        let mut first = Edits::new();
        first.add_atom(AtomForm::from_element(Element::N));
        let editor = empty.apply(first).unwrap();
        let mut second = Edits::new();
        second.add_atom(AtomForm::from_element(Element::O));
        second.remove_atom(AtomHandle::Id(AtomId(7)));

        let error = if tracked {
            editor.tracked_apply(second).err().unwrap()
        } else {
            editor.apply(second).err().unwrap()
        };

        assert_eq!(
            error,
            MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 7,
                count: 2,
            })
        );
    }
}
