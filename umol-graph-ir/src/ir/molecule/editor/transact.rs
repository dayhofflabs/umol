//! Editor batch application and detached transaction journals.
//!
//! `transact(edits)` applies each `Edit` in order, records realized `Undo`
//! entries, and either returns a rollback-capable `Transaction` or reverse-
//! replays the journal before surfacing a `TransactionError`.
//! `apply(edits)` consumes the editor and returns its modified state without
//! constructing an undo journal. A failed application drops the consumed
//! editor, so partially applied state cannot escape.

use std::mem;

use thiserror::Error;
use umol_graph_core::Correspondence;

use super::MoleculeEditor;
use crate::ir::correspondence::MoleculeCorrespondence;
use crate::ir::edit::{Edits, Undo};
use crate::ir::entity::EntityKind;
use crate::ir::molecule::apply::ApplicationState;

#[derive(Debug, Error, PartialEq, Eq, Clone)]
pub enum TransactionError {
    #[error("{kind} handle {index} is out of range for {count} entries")]
    HandleOutOfRange {
        kind: EntityKind,
        index: usize,
        count: usize,
    },

    #[error("{kind} handle {index} refers to a removed entity")]
    HandleRemoved { kind: EntityKind, index: usize },

    #[error("duplicate {kind} in removal batch")]
    DuplicateRemoval { kind: EntityKind },

    /// `Set*Field` or `Set*Constraint`: current state does not match the
    /// edit's `old` payload.
    #[error("precondition failed: old state does not match current")]
    OldStateMismatch,

    /// `Remove*Constraint` with a value that's not present.
    #[error("missing constraint entry on remove")]
    MissingEntry,

    /// Edit shape is structurally invalid.
    #[error("malformed edit: {0}")]
    MalformedEdit(&'static str),

    #[error("rollback failed after apply error: apply={apply}; rollback={rollback}")]
    RollbackFailed {
        apply: Box<TransactionError>,
        rollback: Box<TransactionError>,
    },

    /// The rollback journal cannot be structurally applied to the supplied editor state.
    ///
    /// A transaction guarantees restoration under Molecule::normalized_eq when rolled back
    /// against its post-transaction state or the end of its appended consecutive chain.
    /// Other states are rejected when a required receiver or reconstruction entry is absent;
    /// structurally compatible but unrelated states are outside that guarantee.
    #[error("rollback journal does not match editor state")]
    RollbackStateMismatch,
}

/// Detached journal of the realized undos for one successfully applied edit batch.
///
/// Detachment permits journals for consecutive transactions to be appended and rolled back as a
/// unit. Restoration under Molecule::normalized_eq is guaranteed for the post-transaction editor
/// state or the end of the consecutively appended chain. Manipulated history does not panic but
/// has no specified result. Invalid undo targets or compactions may return
/// [`TransactionError::RollbackStateMismatch`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Transaction {
    undo: Vec<Undo>,
}

impl Transaction {
    pub fn undos(&self) -> &[Undo] {
        &self.undo
    }

    /// Append the rollback journal for a transaction applied after this one.
    pub fn append(&mut self, later: Self) {
        self.undo.extend(later.undo);
    }

    /// Reverse the journal against its exact post-transaction editor state.
    ///
    /// Restores the pre-transaction state under Molecule::normalized_eq when the editor is the
    /// state produced by this transaction or its appended consecutive chain. Manipulated history
    /// does not panic but has no specified result.
    ///
    /// # Errors
    ///
    /// Returns [`TransactionError::RollbackStateMismatch`] for an incompatible undo target,
    /// field value, or compaction.
    pub fn rollback(self, editor: &mut MoleculeEditor) -> Result<(), TransactionError> {
        rollback_journal(editor, self.undo)
    }

    /// Roll back and return the correspondence from the rollback input to the restored state.
    ///
    /// The returned witness covers this rollback, not the editor's entire session. No intermediate
    /// molecule is published. Restored entities without a rollback-input partner remain unmatched.
    ///
    /// # Errors
    ///
    /// Returns the same error and leaves the same editor state as [`Self::rollback`].
    ///
    /// # Semantic properties
    ///
    /// On the transaction's exact post-state, the witness is the inverse of its forward
    /// correspondence. For appended transactions it is the inverse of their composed witness.
    pub fn tracked_rollback(
        self,
        editor: &mut MoleculeEditor,
    ) -> Result<MoleculeCorrespondence, TransactionError> {
        let identity = MoleculeCorrespondence::new(
            Correspondence::identity(editor.atom_count()),
            Correspondence::identity(editor.bond_count()),
            Correspondence::identity(editor.dative_bond_count()),
            Correspondence::identity(editor.aromatic_system_count()),
            Correspondence::identity(editor.multicenter_bond_count()),
            Correspondence::identity(editor.noncovalent_bond_count()),
            Correspondence::identity(editor.stereo_atom_count()),
            Correspondence::identity(editor.stereo_bond_count()),
        );
        let session = mem::replace(&mut editor.correspondence, identity);
        let result = self.rollback(editor);
        let correspondence =
            mem::replace(&mut editor.correspondence, MoleculeCorrespondence::empty());
        editor.correspondence = session
            .compose(&correspondence)
            .expect("rollback correspondence starts in the current editor id spaces");
        result?;
        Ok(correspondence)
    }
}

impl MoleculeEditor {
    /// Apply an ordered [`Edits`] batch atomically. On success, returns a rollback
    /// transaction. On any apply failure, reverse-replays the already-created
    /// undo journal.
    pub fn transact(&mut self, edits: Edits) -> Result<Transaction, TransactionError> {
        let correspondence = self.correspondence.clone();
        let mut journal: Vec<Undo> = Vec::with_capacity(edits.len());
        let mut state = ApplicationState::new(&self.molecule);
        for edit in edits {
            match self.molecule.apply_edit_with_undo(edit, &mut state) {
                Ok(Some(undo)) => journal.push(undo),
                Ok(None) => {}
                Err(apply) => {
                    if let Err(rollback) = rollback_journal(self, journal) {
                        return Err(TransactionError::RollbackFailed {
                            apply: Box::new(apply),
                            rollback: Box::new(rollback),
                        });
                    }
                    self.correspondence = correspondence;
                    return Err(apply);
                }
            }
        }
        Ok(Transaction { undo: journal })
    }

    /// Apply a batch atomically, returning its transaction and input-to-result correspondence.
    ///
    /// The witness starts at this batch's input, independently of the editor's session origin.
    /// The transaction retains its ordinary undo journal. No molecule is published.
    ///
    /// # Errors
    ///
    /// Returns the same error and leaves the same editor state as [`Self::transact`].
    ///
    /// # Semantic properties
    ///
    /// Discarding the witness gives the same transaction and editor state as the plain operation.
    /// Composing the previous session correspondence with this witness gives the new session.
    pub fn tracked_transact(
        &mut self,
        edits: Edits,
    ) -> Result<(Transaction, MoleculeCorrespondence), TransactionError> {
        let identity = MoleculeCorrespondence::new(
            Correspondence::identity(self.atom_count()),
            Correspondence::identity(self.bond_count()),
            Correspondence::identity(self.dative_bond_count()),
            Correspondence::identity(self.aromatic_system_count()),
            Correspondence::identity(self.multicenter_bond_count()),
            Correspondence::identity(self.noncovalent_bond_count()),
            Correspondence::identity(self.stereo_atom_count()),
            Correspondence::identity(self.stereo_bond_count()),
        );
        let session = mem::replace(&mut self.correspondence, identity);
        let result = self.transact(edits);
        let correspondence =
            mem::replace(&mut self.correspondence, MoleculeCorrespondence::empty());
        self.correspondence = session
            .compose(&correspondence)
            .expect("batch correspondence starts in the current editor id spaces");
        Ok((result?, correspondence))
    }

    /// Apply an ordered [`Edits`] batch without constructing an undo journal.
    ///
    /// The editor is consumed so that a failed batch cannot expose partially applied state. On
    /// success, the returned editor remains transient; [`MoleculeEditor::try_build`] or
    /// [`MoleculeEditor::build`] performs the molecule-integrity publication gate.
    ///
    /// # Errors
    ///
    /// Returns [`TransactionError`] when an edit handle, precondition, or shape is invalid for the
    /// evolving editor state.
    ///
    /// # Semantic properties
    ///
    /// For every edit batch accepted by [`Self::transact`] from the same initial editor,
    /// successfully building the returned editor produces the same molecule as building the
    /// post-transaction editor. On failure, no intermediate editor state is returned.
    pub fn apply(mut self, edits: Edits) -> Result<Self, TransactionError> {
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
    /// Discarding the witness gives the same editor as the plain operation. Composing the previous
    /// session correspondence with this witness gives the resulting session correspondence.
    pub fn tracked_apply(
        mut self,
        edits: Edits,
    ) -> Result<(Self, MoleculeCorrespondence), TransactionError> {
        let identity = MoleculeCorrespondence::new(
            Correspondence::identity(self.atom_count()),
            Correspondence::identity(self.bond_count()),
            Correspondence::identity(self.dative_bond_count()),
            Correspondence::identity(self.aromatic_system_count()),
            Correspondence::identity(self.multicenter_bond_count()),
            Correspondence::identity(self.noncovalent_bond_count()),
            Correspondence::identity(self.stereo_atom_count()),
            Correspondence::identity(self.stereo_bond_count()),
        );
        let session = mem::replace(&mut self.correspondence, identity);
        let mut editor = self.apply(edits)?;
        let correspondence =
            mem::replace(&mut editor.correspondence, MoleculeCorrespondence::empty());
        editor.correspondence = session
            .compose(&correspondence)
            .expect("batch correspondence starts in the current editor id spaces");
        Ok((editor, correspondence))
    }
}

fn rollback_journal(
    editor: &mut MoleculeEditor,
    journal: Vec<Undo>,
) -> Result<(), TransactionError> {
    for undo in journal.into_iter().rev() {
        editor.molecule.apply_undo(undo);
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::{Compaction, GraphCompaction};

    use super::*;
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::{AtomForm, ElementForm};
    use crate::ir::bond::BondForm;
    use crate::ir::compact::MoleculeCompaction;
    use crate::ir::constraint::{
        AromaticSystemConstraintForm, AtomConstraintForm, BondConstraintForm, Constraint,
        DativeBondConstraintForm, MoleculeConstraint, MulticenterBondConstraintForm,
        NoncovalentBondConstraintForm, RelationalConstraint, RingScope, StereoAtomConstraintForm,
        StereoBondConstraintForm, StereogenicityForm,
    };
    use crate::ir::dative::DativeBondForm;
    use crate::ir::edit::{
        AddBond, AddedAtom, AromaticSystemFieldChange, AromaticSystemHandle, AtomFieldChange,
        AtomHandle, BondFieldChange, BondHandle, CascadedConstraints, ConstraintEdit,
        DativeBondFieldChange, DativeBondHandle, Edit, EntityHandle, MulticenterBondFieldChange,
        MulticenterBondHandle, NoncovalentBondFieldChange, NoncovalentBondHandle,
        RemovedAromaticSystem, RemovedAtom, RemovedConstraint, RemovedOverlays,
        StereoAtomFieldChange, StereoAtomHandle, StereoBondFieldChange, StereoBondHandle,
    };
    use crate::ir::entity::Entity;
    use crate::ir::id::{
        AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
        StereoAtomId, StereoBondId,
    };
    use crate::ir::ligand::{StereoLigand, StereoLigandKind};
    use crate::ir::molecule::Molecule;
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{
        NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm,
    };
    use crate::ir::num::NumForm;
    use crate::ir::stereo::{
        CisTransStereoForm, StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset,
        StereoKind,
    };
    use crate::ir::traits::Normalize;
    use crate::ir::{BooleanForm, ModifiedConstraint};

    #[fixture]
    fn empty() -> MoleculeEditor {
        Molecule::default().edit()
    }

    #[fixture]
    fn one_atom() -> MoleculeEditor {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b
    }

    #[fixture]
    fn diatomic() -> MoleculeEditor {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
        b
    }

    #[rstest]
    fn test_molecule_editor_transact_add_atom(mut empty: MoleculeEditor) {
        let mut edits = Edits::new();
        edits.add_atom(AtomForm::from_element(Element::C));
        let tx = empty.transact(edits).unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedTopology { atoms, bonds }]
                if atoms.iter().map(|a| a.id).collect::<Vec<_>>() == vec![AtomId(0)]
                    && bonds.is_empty()
        ));
        let built = empty.build();
        assert_eq!(built.atoms().count(), 1);
        assert_eq!(
            built.atom(AtomId(0)).attributes().element,
            ElementForm::Lit(Element::C)
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_add_atoms(mut empty: MoleculeEditor) {
        let mut edits = Edits::new();
        edits.add_atoms([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
        ]);
        let tx = empty.transact(edits).unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedTopology { atoms, bonds }]
                if atoms.iter().map(|a| a.id).collect::<Vec<_>>() == vec![AtomId(0), AtomId(1)]
                    && bonds.is_empty()
        ));
        let built = empty.build();
        assert_eq!(built.atoms().count(), 2);
        assert_eq!(
            built.atom(AtomId(0)).attributes().element,
            ElementForm::Lit(Element::C)
        );
        assert_eq!(
            built.atom(AtomId(1)).attributes().element,
            ElementForm::Lit(Element::N)
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_add_bond_via_handle(mut empty: MoleculeEditor) {
        let mut edits = Edits::new();
        let atoms = edits.add_atoms([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ]);
        edits.add_bond(atoms[0].clone(), atoms[1].clone(), BondForm::from_order(1));
        let tx = empty.transact(edits).unwrap();
        assert!(matches!(
            tx.undos(),
            [
                Undo::RemoveAddedTopology { atoms, bonds },
                Undo::RemoveAddedTopology { atoms: bond_atoms, bonds: added_bonds },
            ] if atoms.iter().map(|a| a.id).collect::<Vec<_>>() == vec![AtomId(0), AtomId(1)]
                && bonds.is_empty()
                && bond_atoms.is_empty()
                && added_bonds.iter().map(|b| b.id).collect::<Vec<_>>() == vec![BondId(0)]
        ));
    }

    #[rstest]
    fn test_molecule_editor_transact_rollback(mut one_atom: MoleculeEditor) {
        let before = one_atom.clone().build();
        // Mid-batch failure (out-of-range id on edit 2) rolls back the
        // already-applied AddAtom on edit 1.
        let mut edits = Edits::new();
        edits.add_atom(AtomForm::from_element(Element::N));
        edits.remove_atom(AtomHandle::Id(AtomId(99)));
        let err = one_atom.transact(edits).unwrap_err();
        assert_eq!(
            err,
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 99,
                count: 1,
            }
        );
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_atom_field(mut one_atom: MoleculeEditor) {
        let tx = one_atom
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::ModifyAtomField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::default(),
                },
            }],
        );
        assert_eq!(
            one_atom.build().atom(AtomId(0)).attributes().charge,
            NumForm::Lit(1)
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_atom_field_error(mut one_atom: MoleculeEditor) {
        let err = one_atom
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(99),
                    new: NumForm::Lit(1),
                },
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    #[case::created_out_of_range(
        Edits::from_iter([Edit::AddBonds {
            bonds: vec![AddBond {
                endpoints: [AtomHandle::New(5), AtomHandle::New(6)],
                attributes: BondForm::default(),
            }],
        }]),
        TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 5,
            count: 0,
        },
    )]
    #[case::initial_out_of_range(
        Edits::from_iter([Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: Vec::new(),
        }]),
        TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 0,
            count: 0,
        },
    )]
    fn test_molecule_editor_transact_handle_error(
        mut empty: MoleculeEditor,
        #[case] edits: Edits,
        #[case] expected: TransactionError,
    ) {
        let before = empty.clone().build();
        assert_eq!(empty.transact(edits).unwrap_err(), expected);
        assert_eq!(empty.build(), before);
    }

    #[rstest]
    #[case::initial(
        1,
        Edits::from_iter([
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            },
        ]),
    )]
    #[case::created(
        0,
        Edits::from_iter([
            Edit::AddAtoms {
                atoms: vec![AtomForm::from_element(Element::C)],
            },
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::New(0)],
                bonds: Vec::new(),
            },
            Edit::ModifyAtomField {
                id: AtomHandle::New(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            },
        ]),
    )]
    fn test_molecule_editor_transact_handle_removed_error(
        #[case] initial_atom_count: usize,
        #[case] edits: Edits,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..initial_atom_count {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        let before = editor.clone().build();

        assert_eq!(
            editor.transact(edits).unwrap_err(),
            TransactionError::HandleRemoved {
                kind: EntityKind::Atom,
                index: 0,
            }
        );
        assert_eq!(editor.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_handles_initial() {
        let mut editor = Molecule::default().edit();
        editor.add_atom(AtomForm::from_element(Element::C));
        editor.add_atom(AtomForm::from_element(Element::N));
        editor.add_atom(AtomForm::from_element(Element::O));
        let edits = Edits::from_iter([
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(1)),
                change: AtomFieldChange::Element {
                    old: ElementForm::Lit(Element::N),
                    new: ElementForm::Lit(Element::F),
                },
            },
            Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(2)),
                change: AtomFieldChange::Element {
                    old: ElementForm::Lit(Element::O),
                    new: ElementForm::Lit(Element::Cl),
                },
            },
        ]);

        editor.transact(edits).unwrap();

        assert_eq!(
            (0..editor.atom_count())
                .map(|index| editor
                    .atom(AtomId(index as u32))
                    .attributes()
                    .element
                    .clone())
                .collect::<Vec<_>>(),
            vec![ElementForm::Lit(Element::F), ElementForm::Lit(Element::Cl)]
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_handles_created() {
        let mut editor = Molecule::default().edit();
        let mut edits = Edits::new();
        let atoms = edits.add_atoms([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
        ]);
        edits.remove_atom(atoms[0].clone());
        edits.push(Edit::ModifyAtomField {
            id: atoms[1].clone(),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        });

        editor.transact(edits).unwrap();

        assert_eq!(editor.atom_count(), 1);
        assert_eq!(
            (
                editor.atom(AtomId(0)).attributes().element.clone(),
                editor.atom(AtomId(0)).attributes().charge.clone(),
            ),
            (ElementForm::Lit(Element::N), NumForm::Lit(1))
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_handles_reuse() {
        let mut editor = Molecule::default().edit();
        let mut edits = Edits::new();
        let removed = edits.add_atom(AtomForm::from_element(Element::C));
        edits.remove_atom(removed);
        let surviving = edits.add_atom(AtomForm::from_element(Element::N));
        edits.push(Edit::ModifyAtomField {
            id: surviving,
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(-1),
            },
        });

        editor.transact(edits).unwrap();

        assert_eq!(editor.atom_count(), 1);
        assert_eq!(
            (
                editor.atom(AtomId(0)).attributes().element.clone(),
                editor.atom(AtomId(0)).attributes().charge.clone(),
            ),
            (ElementForm::Lit(Element::N), NumForm::Lit(-1))
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_handles_per_kind() {
        let mut editor = Molecule::default().edit();
        let before = editor.clone().build();
        let mut edits = Edits::new();
        let atoms = edits.add_atoms([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::F),
        ]);
        let bonds = edits.add_bonds([
            AddBond {
                endpoints: [atoms[0].clone(), atoms[1].clone()],
                attributes: BondForm::from_order(1),
            },
            AddBond {
                endpoints: [atoms[1].clone(), atoms[2].clone()],
                attributes: BondForm::from_order(1),
            },
            AddBond {
                endpoints: [atoms[2].clone(), atoms[3].clone()],
                attributes: BondForm::from_order(1),
            },
        ]);
        let dative = edits.add_dative_bond(
            vec![atoms[0].clone()],
            atoms[1].clone(),
            DativeBondForm::from_order(1),
        );
        let aromatic = edits.add_aromatic_system(
            vec![atoms[0].clone(), atoms[1].clone()],
            AromaticSystemForm::default(),
        );
        let multicenter = edits.add_multicenter_bond(
            vec![atoms[0].clone(), atoms[1].clone()],
            MulticenterBondForm::default(),
        );
        let noncovalent = edits.add_noncovalent_bond(
            [atoms[0].clone(), atoms[1].clone()],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        let stereo_atom = edits.add_stereo_atom(
            atoms[1].clone(),
            vec![
                (atoms[0].clone(), StereoLigandKind::Atom),
                (atoms[2].clone(), StereoLigandKind::Atom),
                (atoms[3].clone(), StereoLigandKind::Atom),
                (atoms[1].clone(), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let stereo_bond = edits.add_stereo_bond(
            bonds[1].clone(),
            vec![
                (atoms[0].clone(), StereoLigandKind::Atom),
                (atoms[1].clone(), StereoLigandKind::ImplicitHydrogen),
                (atoms[3].clone(), StereoLigandKind::Atom),
                (atoms[2].clone(), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        edits.push(Edit::ModifyAtomField {
            id: atoms[0].clone(),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        });
        edits.push(Edit::ModifyBondField {
            id: bonds[0].clone(),
            change: BondFieldChange::Order {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        });
        edits.push(Edit::ModifyDativeBondField {
            id: dative,
            change: DativeBondFieldChange::Order {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        });
        edits.push(Edit::ModifyAromaticSystemField {
            id: aromatic,
            change: AromaticSystemFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        });
        edits.push(Edit::ModifyMulticenterBondField {
            id: multicenter,
            change: MulticenterBondFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(-1),
            },
        });
        edits.push(Edit::ModifyNoncovalentBondField {
            id: noncovalent,
            change: NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
            },
        });
        edits.push(Edit::ModifyStereoAtomField {
            id: stereo_atom,
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            },
        });
        edits.push(Edit::ModifyStereoBondField {
            id: stereo_bond,
            change: StereoBondFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
            },
        });

        let transaction = editor.transact(edits).unwrap();

        assert_eq!(editor.atom(AtomId(0)).attributes().charge, NumForm::Lit(1));
        assert_eq!(editor.bond(BondId(0)).attributes().order, NumForm::Lit(2));
        assert_eq!(
            editor.dative_bond(DativeBondId(0)).attributes().order,
            NumForm::Lit(2)
        );
        assert_eq!(
            editor
                .aromatic_system(AromaticSystemId(0))
                .attributes()
                .charge,
            NumForm::Lit(1)
        );
        assert_eq!(
            editor
                .multicenter_bond(MulticenterBondId(0))
                .attributes()
                .charge,
            NumForm::Lit(-1)
        );
        assert_eq!(
            editor
                .noncovalent_bond(NoncovalentBondId(0))
                .attributes()
                .kind,
            NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic)
        );
        assert_eq!(
            editor
                .stereo_atom(StereoAtomId(0))
                .attributes()
                .configuration,
            StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0))
        );
        assert_eq!(
            editor
                .stereo_bond(StereoBondId(0))
                .attributes()
                .configuration,
            StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0))
        );

        transaction.rollback(&mut editor).unwrap();
        assert_eq!(editor.build(), before);
    }

    #[rstest]
    #[case::first(0)]
    #[case::middle(1)]
    #[case::last(2)]
    fn test_molecule_editor_transact_add_bonds_error(
        mut diatomic: MoleculeEditor,
        #[case] invalid_position: usize,
    ) {
        let before = diatomic.clone().build();
        let mut bonds = vec![
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                attributes: BondForm::from_order(1),
            },
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                attributes: BondForm::from_order(2),
            },
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                attributes: BondForm::from_order(3),
            },
        ];
        bonds[invalid_position].endpoints[1] = AtomHandle::Id(AtomId(9));

        assert_eq!(
            diatomic
                .transact(Edits::from_iter([Edit::AddBonds { bonds }]))
                .unwrap_err(),
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 9,
                count: 2,
            }
        );
        assert_eq!(diatomic.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_topology(mut diatomic: MoleculeEditor) {
        let tx = diatomic
            .transact(Edits::from_iter([Edit::RemoveTopology {
                atoms: Vec::new(),
                bonds: vec![BondHandle::Id(BondId(0))],
            }]))
            .unwrap();

        assert_eq!(diatomic.bond_count(), 0);
        let [Undo::RestoreRemovedTopology {
            atoms,
            bonds,
            overlays,
            compaction,
            ..
        }] = tx.undos()
        else {
            panic!("RemoveTopology should produce one topology-restore undo")
        };
        assert!(atoms.is_empty());
        assert_eq!(
            bonds.iter().map(|b| b.id).collect::<Vec<_>>(),
            vec![BondId(0)]
        );
        assert!(overlays.dative_bonds.is_empty());
        assert_eq!(compaction.compact_bond(BondId(0)), None);
    }

    #[rstest]
    #[case::atom(Edit::ModifyAtomConstraint { id: AtomHandle::Id(AtomId(0)), old: None, new: None })]
    #[case::bond(Edit::ModifyBondConstraint { id: BondHandle::Id(BondId(0)), old: None, new: None })]
    #[case::dative_bond(Edit::ModifyDativeBondConstraint { id: DativeBondHandle::Id(DativeBondId(0)), old: None, new: None })]
    #[case::aromatic_system(Edit::ModifyAromaticSystemConstraint { id: AromaticSystemHandle::Id(AromaticSystemId(0)), old: None, new: None })]
    #[case::multicenter_bond(Edit::ModifyMulticenterBondConstraint { id: MulticenterBondHandle::Id(MulticenterBondId(0)), old: None, new: None })]
    #[case::noncovalent_bond(Edit::ModifyNoncovalentBondConstraint { id: NoncovalentBondHandle::Id(NoncovalentBondId(0)), old: None, new: None })]
    #[case::stereo_atom(Edit::ModifyStereoAtomConstraint { id: StereoAtomHandle::Id(StereoAtomId(0)), kind: None, old: None, new: None })]
    #[case::stereo_bond(Edit::ModifyStereoBondConstraint { id: StereoBondHandle::Id(StereoBondId(0)), kind: None, old: None, new: None })]
    fn test_molecule_editor_transact_constraint_identity(
        mut batched_overlays: MoleculeEditor,
        #[case] edit: Edit,
    ) {
        let expected = batched_overlays.snapshot().unwrap();
        let transaction = batched_overlays.transact(Edits::from_iter([edit])).unwrap();

        assert_eq!(transaction.undos(), &[]);
        assert_eq!(batched_overlays.try_build(), Ok(expected));
    }

    #[rstest]
    fn test_molecule_editor_transact_add_atom_constraint(mut one_atom: MoleculeEditor) {
        one_atom
            .transact(Edits::from_iter([
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(0)),
                    old: None,
                    new: Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 1)),
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(0)),
                    old: None,
                    new: Some(AtomConstraintForm::ring_membership(RingScope::Size(6), 1)),
                },
            ]))
            .unwrap();
        let next = one_atom.build();
        let cs: Vec<_> = next
            .atom(AtomId(0))
            .attributes()
            .constraints
            .iter()
            .cloned()
            .collect();
        assert_eq!(
            cs,
            vec![
                AtomConstraintForm::ring_membership(RingScope::Size(5), 1),
                AtomConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]
        );
    }

    #[rstest]
    #[case::singleton_set(NumForm::Lit(1), NumForm::lit_set([1]))]
    fn test_molecule_editor_transact_modify_atom_field_canonical(
        mut one_atom: MoleculeEditor,
        #[case] current: NumForm,
        #[case] old: NumForm,
    ) {
        // The modify's recorded `old` is equivalent to — but structurally distinct from — the
        // stored charge, so the old-state check passes (structural `!=` would raise `OldStateMismatch`).
        one_atom.atom_mut(AtomId(0)).attributes_mut().charge = current;
        one_atom
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old,
                    new: NumForm::Lit(2),
                },
            }]))
            .unwrap();
        assert_eq!(
            one_atom.atom_mut(AtomId(0)).attributes_mut().charge,
            NumForm::Lit(2)
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_modify_atom_constraint_absent_error(
        mut one_atom: MoleculeEditor,
    ) {
        let err = one_atom
            .transact(Edits::from_iter([Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old: Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 1)),
                new: None,
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    #[case::introduce(None, Some(AtomConstraintForm::valence(4)), Some(NumForm::Lit(4)))]
    #[case::replace(
        Some(AtomConstraintForm::valence(3)),
        Some(AtomConstraintForm::valence(4)),
        Some(NumForm::Lit(4))
    )]
    #[case::remove(Some(AtomConstraintForm::valence(3)), None, None)]
    fn test_molecule_editor_transact_set_atom_constraint(
        mut one_atom: MoleculeEditor,
        #[case] old: Option<AtomConstraintForm>,
        #[case] new: Option<AtomConstraintForm>,
        #[case] expected: Option<NumForm>,
    ) {
        if let Some(c) = old.clone() {
            one_atom
                .atom_mut(AtomId(0))
                .attributes_mut()
                .constraints
                .set(c);
        }
        one_atom
            .transact(Edits::from_iter([Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old,
                new,
            }]))
            .unwrap();
        assert_eq!(
            one_atom
                .atom_mut(AtomId(0))
                .attributes_mut()
                .constraints
                .valence(),
            expected.as_ref()
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_bond_constraint(mut diatomic: MoleculeEditor) {
        diatomic
            .transact(Edits::from_iter([Edit::ModifyBondConstraint {
                id: BondHandle::Id(BondId(0)),
                old: None,
                new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
            }]))
            .unwrap();
        assert!(diatomic
            .bond_mut(BondId(0))
            .attributes_mut()
            .constraints
            .iter()
            .any(|c| *c == BondConstraintForm::Aromatic(BooleanForm::Lit(true))));
    }

    #[rstest]
    fn test_molecule_editor_transact_add_molecule_constraint(mut empty: MoleculeEditor) {
        let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        empty
            .transact(Edits::from_iter([Edit::AddMoleculeConstraint {
                constraint: c.clone().into(),
            }]))
            .unwrap();
        assert_eq!(empty.constraints_mut().as_slice(), &[c]);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_molecule_constraint(mut empty: MoleculeEditor) {
        let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        empty.constraints_mut().push(c.clone());
        empty
            .transact(Edits::from_iter([Edit::RemoveMoleculeConstraint {
                constraint: c.clone().into(),
            }]))
            .unwrap();
        assert!(empty.constraints_mut().as_slice().is_empty());
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_molecule_constraint_absent_error(
        mut empty: MoleculeEditor,
    ) {
        let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        empty.constraints_mut().push(c.clone());
        let err = empty
            .transact(Edits::from_iter([Edit::RemoveMoleculeConstraint {
                constraint: Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: None,
                    sum: NumForm::Lit(0),
                })
                .into(),
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::MissingEntry);
        assert_eq!(empty.constraints_mut().as_slice(), &[c]);
    }

    #[rstest]
    fn test_molecule_editor_transact_molecule_constraint_initial(
        mut batched_overlays: MoleculeEditor,
    ) {
        let constraint = Constraint::And(vec![
            Constraint::Atom(AtomId(1), AtomConstraintForm::valence(3_i64)),
            Constraint::Bond(BondId(1), BondConstraintForm::aromatic(true)),
            Constraint::DativeBond(DativeBondId(1), DativeBondConstraintForm::aromatic(true)),
            Constraint::AromaticSystem(
                AromaticSystemId(1),
                AromaticSystemConstraintForm::electron_count(6_i64),
            ),
            Constraint::MulticenterBond(
                MulticenterBondId(1),
                MulticenterBondConstraintForm::electron_count(2_i64),
            ),
            Constraint::NoncovalentBond(
                NoncovalentBondId(1),
                NoncovalentBondConstraintForm::intramolecular(true),
            ),
            Constraint::StereoAtom(
                StereoAtomId(1),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::StereoBond(
                StereoBondId(1),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
        ]);
        let before = batched_overlays.clone().build();
        let mut edits = Edits::new();
        edits.add_molecule_constraint(constraint.clone().into());

        let transaction = batched_overlays.transact(edits).unwrap();
        assert_eq!(batched_overlays.constraints().as_slice(), &[constraint]);

        transaction.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_molecule_constraint_created(mut empty: MoleculeEditor) {
        let before = empty.clone().build();
        let mut edits = Edits::new();
        let atoms = edits.add_atoms([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
        ]);
        let bond = edits.add_bond(atoms[0].clone(), atoms[1].clone(), BondForm::from_order(1));
        let dative = edits.add_dative_bond(
            vec![atoms[0].clone()],
            atoms[1].clone(),
            DativeBondForm::from_order(1),
        );
        let aromatic = edits.add_aromatic_system(atoms.clone(), AromaticSystemForm::default());
        let multicenter = edits.add_multicenter_bond(atoms.clone(), MulticenterBondForm::default());
        let noncovalent = edits.add_noncovalent_bond(
            [atoms[0].clone(), atoms[1].clone()],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        let stereo_atom = edits.add_stereo_atom(
            atoms[0].clone(),
            vec![
                (atoms[1].clone(), StereoLigandKind::Atom),
                (atoms[0].clone(), StereoLigandKind::ImplicitHydrogen),
                (atoms[0].clone(), StereoLigandKind::LonePair),
                (atoms[1].clone(), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let stereo_bond = edits.add_stereo_bond(
            bond.clone(),
            vec![
                (atoms[0].clone(), StereoLigandKind::Atom),
                (atoms[0].clone(), StereoLigandKind::ImplicitHydrogen),
                (atoms[1].clone(), StereoLigandKind::Atom),
                (atoms[1].clone(), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let source = Constraint::And(vec![
            Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
            Constraint::Bond(BondId(7), BondConstraintForm::aromatic(true)),
            Constraint::DativeBond(DativeBondId(7), DativeBondConstraintForm::aromatic(true)),
            Constraint::AromaticSystem(
                AromaticSystemId(7),
                AromaticSystemConstraintForm::electron_count(6_i64),
            ),
            Constraint::MulticenterBond(
                MulticenterBondId(7),
                MulticenterBondConstraintForm::electron_count(2_i64),
            ),
            Constraint::NoncovalentBond(
                NoncovalentBondId(7),
                NoncovalentBondConstraintForm::intramolecular(true),
            ),
            Constraint::StereoAtom(
                StereoAtomId(7),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::StereoBond(
                StereoBondId(7),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::Relational(RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(7),
                parallel: BondId(7),
            }),
        ]);
        let mappings = HashMap::from([
            (
                Entity::Atom(AtomId(7)),
                EntityHandle::Atom(atoms[0].clone()),
            ),
            (Entity::Bond(BondId(7)), EntityHandle::Bond(bond)),
            (
                Entity::DativeBond(DativeBondId(7)),
                EntityHandle::DativeBond(dative),
            ),
            (
                Entity::AromaticSystem(AromaticSystemId(7)),
                EntityHandle::AromaticSystem(aromatic),
            ),
            (
                Entity::MulticenterBond(MulticenterBondId(7)),
                EntityHandle::MulticenterBond(multicenter),
            ),
            (
                Entity::NoncovalentBond(NoncovalentBondId(7)),
                EntityHandle::NoncovalentBond(noncovalent),
            ),
            (
                Entity::StereoAtom(StereoAtomId(7)),
                EntityHandle::StereoAtom(stereo_atom),
            ),
            (
                Entity::StereoBond(StereoBondId(7)),
                EntityHandle::StereoBond(stereo_bond),
            ),
        ]);
        edits.add_molecule_constraint(
            ConstraintEdit::new(source, |entity| mappings.get(&entity).cloned()).unwrap(),
        );
        let expected = Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3_i64)),
            Constraint::Bond(BondId(0), BondConstraintForm::aromatic(true)),
            Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(true)),
            Constraint::AromaticSystem(
                AromaticSystemId(0),
                AromaticSystemConstraintForm::electron_count(6_i64),
            ),
            Constraint::MulticenterBond(
                MulticenterBondId(0),
                MulticenterBondConstraintForm::electron_count(2_i64),
            ),
            Constraint::NoncovalentBond(
                NoncovalentBondId(0),
                NoncovalentBondConstraintForm::intramolecular(true),
            ),
            Constraint::StereoAtom(
                StereoAtomId(0),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::StereoBond(
                StereoBondId(0),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::Relational(RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(0),
                parallel: BondId(0),
            }),
        ]);

        let transaction = empty.transact(edits).unwrap();
        assert_eq!(empty.constraints().as_slice(), &[expected]);

        transaction.rollback(&mut empty).unwrap();
        assert_eq!(empty.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_molecule_constraint_compaction(
        mut batched_overlays: MoleculeEditor,
    ) {
        let removed = Constraint::AromaticSystem(
            AromaticSystemId(1),
            AromaticSystemConstraintForm::electron_count(6_i64),
        );
        let added = Constraint::AromaticSystem(
            AromaticSystemId(1),
            AromaticSystemConstraintForm::electron_count(4_i64),
        );
        batched_overlays.constraints_mut().push(removed.clone());
        let before = batched_overlays.clone().build();
        let mut edits = Edits::from_iter([Edit::RemoveAromaticSystems {
            removes: vec![(
                AromaticSystemHandle::Id(AromaticSystemId(0)),
                vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                AromaticSystemForm::default(),
            )],
        }]);
        edits.add_molecule_constraint(added.into());
        edits.remove_molecule_constraint(removed.into());

        let transaction = batched_overlays.transact(edits).unwrap();
        assert_eq!(
            batched_overlays.constraints().as_slice(),
            &[Constraint::AromaticSystem(
                AromaticSystemId(0),
                AromaticSystemConstraintForm::electron_count(4_i64),
            )],
        );

        transaction.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    #[case::forward(
        0,
        Edits::from_iter([Edit::AddMoleculeConstraint {
            constraint: ConstraintEdit::new(
                Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
                |_| Some(EntityHandle::Atom(AtomHandle::New(0))),
            ).unwrap(),
        }]),
        TransactionError::HandleOutOfRange { kind: EntityKind::Atom, index: 0, count: 0 },
    )]
    #[case::removed(
        1,
        Edits::from_iter([
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            Edit::AddMoleculeConstraint {
                constraint: ConstraintEdit::new(
                    Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
                    |_| Some(EntityHandle::Atom(AtomHandle::Id(AtomId(0)))),
                ).unwrap(),
            },
        ]),
        TransactionError::HandleRemoved { kind: EntityKind::Atom, index: 0 },
    )]
    fn test_molecule_editor_transact_molecule_constraint_error(
        #[case] initial_atom_count: usize,
        #[case] edits: Edits,
        #[case] expected: TransactionError,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..initial_atom_count {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        let before = editor.clone().build();

        assert_eq!(editor.transact(edits), Err(expected));
        assert_eq!(editor.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_topology_atom_error(mut one_atom: MoleculeEditor) {
        let before = one_atom.clone().build();
        let mut edits = Edits::new();
        edits.remove_atom(AtomHandle::Id(AtomId(9)));
        let err = one_atom.transact(edits).unwrap_err();
        assert_eq!(
            err,
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 9,
                count: 1,
            }
        );
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_topology_bond_error(mut diatomic: MoleculeEditor) {
        let before = diatomic.clone().build();
        let mut edits = Edits::new();
        edits.remove_bond(BondHandle::Id(BondId(9)));
        let err = diatomic.transact(edits).unwrap_err();
        assert_eq!(
            err,
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Bond,
                index: 9,
                count: 1,
            }
        );
        assert_eq!(diatomic.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_add_dative_bond_acceptor_error(mut one_atom: MoleculeEditor) {
        let err = one_atom
            .transact(Edits::from_iter([Edit::AddDativeBond {
                donors: vec![],
                acceptor: AtomHandle::Id(AtomId(9)),
                attributes: DativeBondForm::from_order(1),
            }]))
            .unwrap_err();
        assert_eq!(
            err,
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 9,
                count: 1,
            }
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_bond_field(mut diatomic: MoleculeEditor) {
        diatomic
            .transact(Edits::from_iter([Edit::ModifyBondField {
                id: BondHandle::Id(BondId(0)),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            }]))
            .unwrap();
        assert_eq!(diatomic.bond(BondId(0)).attributes().order, NumForm::Lit(2));
    }

    #[rstest]
    fn test_molecule_editor_transact_set_bond_field_error(mut diatomic: MoleculeEditor) {
        let err = diatomic
            .transact(Edits::from_iter([Edit::ModifyBondField {
                id: BondHandle::Id(BondId(0)),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(99),
                    new: NumForm::Lit(2),
                },
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    // Stereo elements — transactional add/remove + undo, and topology-cascade
    // restore. (D3i behavior, tested against D3j's view-based read path.)

    #[fixture]
    fn stereo_atom_skeleton() -> MoleculeEditor {
        let mut b = Molecule::default().edit();
        for el in [Element::C, Element::F, Element::Cl, Element::Br, Element::I] {
            b.add_atom(AtomForm::from_element(el));
        }
        for t in 1u32..=4 {
            b.add_bond(AtomId(0), AtomId(t), BondForm::from_order(1));
        }
        b
    }

    fn tetrahedral_ligands() -> Vec<StereoLigand> {
        (1u32..=4)
            .map(|t| StereoLigand::new(AtomId(t), StereoLigandKind::Atom))
            .collect()
    }

    #[rstest]
    fn test_molecule_editor_transact_add_stereo_atom(mut stereo_atom_skeleton: MoleculeEditor) {
        let before = stereo_atom_skeleton.clone().build();
        let tx = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::AddStereoAtom {
                site: AtomHandle::Id(AtomId(0)),
                ligands: (1u32..=4)
                    .map(|t| (AtomHandle::Id(AtomId(t)), StereoLigandKind::Atom))
                    .collect(),
                attributes: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            }]))
            .unwrap();
        assert_eq!(stereo_atom_skeleton.stereo_atom_count(), 1);
        tx.rollback(&mut stereo_atom_skeleton).unwrap();
        assert_eq!(stereo_atom_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_stereo_atom(mut stereo_atom_skeleton: MoleculeEditor) {
        stereo_atom_skeleton.add_stereo_atom(
            AtomId(0),
            &tetrahedral_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let before = stereo_atom_skeleton.clone().build();
        let tx = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    (1u32..=4)
                        .map(|t| (AtomHandle::Id(AtomId(t)), StereoLigandKind::Atom))
                        .collect(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                )],
            }]))
            .unwrap();
        assert_eq!(stereo_atom_skeleton.stereo_atom_count(), 0);
        tx.rollback(&mut stereo_atom_skeleton).unwrap();
        assert_eq!(stereo_atom_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_stereo_atom_error(
        mut stereo_atom_skeleton: MoleculeEditor,
    ) {
        stereo_atom_skeleton.add_stereo_atom(
            AtomId(0),
            &tetrahedral_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let err = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    (1u32..=4)
                        .map(|t| (AtomHandle::Id(AtomId(t)), StereoLigandKind::Atom))
                        .collect(),
                    // Wrong recorded coset (Th0 vs the stored Th1).
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                )],
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_editor_transact_topology_removal_restores_stereo_atom(
        mut stereo_atom_skeleton: MoleculeEditor,
    ) {
        stereo_atom_skeleton.add_stereo_atom(
            AtomId(0),
            &tetrahedral_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let before = stereo_atom_skeleton.clone().build();
        // Removing a ligand atom cascades the stereo element away.
        let tx = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(1))],
                bonds: Vec::new(),
            }]))
            .unwrap();
        assert_eq!(stereo_atom_skeleton.stereo_atom_count(), 0);
        tx.rollback(&mut stereo_atom_skeleton).unwrap();
        assert_eq!(stereo_atom_skeleton.build(), before);
    }

    #[fixture]
    fn stereo_bond_skeleton() -> MoleculeEditor {
        let mut b = Molecule::default().edit();
        for _ in 0..4 {
            b.add_atom(AtomForm::from_element(Element::C));
        }
        b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
        b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(2));
        b.add_bond(AtomId(2), AtomId(3), BondForm::from_order(1));
        b
    }

    #[rstest]
    fn test_molecule_editor_transact_add_stereo_bond(mut stereo_bond_skeleton: MoleculeEditor) {
        let before = stereo_bond_skeleton.clone().build();
        let tx = stereo_bond_skeleton
            .transact(Edits::from_iter([Edit::AddStereoBond {
                site: BondHandle::Id(BondId(1)),
                ligands: vec![
                    (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(1)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                    (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(2)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                ],
                attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            }]))
            .unwrap();
        assert_eq!(stereo_bond_skeleton.stereo_bond_count(), 1);
        tx.rollback(&mut stereo_bond_skeleton).unwrap();
        assert_eq!(stereo_bond_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_stereo_bond(mut stereo_bond_skeleton: MoleculeEditor) {
        stereo_bond_skeleton.add_stereo_bond(
            BondId(1),
            &[
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let before = stereo_bond_skeleton.clone().build();
        let tx = stereo_bond_skeleton
            .transact(Edits::from_iter([Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(1)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(2)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                )],
            }]))
            .unwrap();
        assert_eq!(stereo_bond_skeleton.stereo_bond_count(), 0);
        tx.rollback(&mut stereo_bond_skeleton).unwrap();
        assert_eq!(stereo_bond_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_stereo_bond_error(
        mut stereo_bond_skeleton: MoleculeEditor,
    ) {
        stereo_bond_skeleton.add_stereo_bond(
            BondId(1),
            &[
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let error = stereo_bond_skeleton
            .transact(Edits::from_iter([Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(1)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(2)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                    ],
                    // Wrong recorded coset (Ct0 vs the stored Ct1).
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                )],
            }]))
            .unwrap_err();

        assert_eq!(error, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_stereo_atom_field(
        mut stereo_atom_skeleton: MoleculeEditor,
    ) {
        stereo_atom_skeleton.add_stereo_atom(
            AtomId(0),
            &tetrahedral_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let before = stereo_atom_skeleton.clone().build();
        let tx = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::ModifyStereoAtomField {
                id: StereoAtomHandle::Id(StereoAtomId(0)),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1),
                    ),
                    new: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0),
                    ),
                },
            }]))
            .unwrap();
        assert_eq!(
            stereo_atom_skeleton
                .stereo_atom(StereoAtomId(0))
                .attributes()
                .configuration,
            StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0),),
        );
        tx.rollback(&mut stereo_atom_skeleton).unwrap();
        assert_eq!(stereo_atom_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_stereo_atom_field_error(
        mut stereo_atom_skeleton: MoleculeEditor,
    ) {
        stereo_atom_skeleton.add_stereo_atom(
            AtomId(0),
            &tetrahedral_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let err = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::ModifyStereoAtomField {
                id: StereoAtomHandle::Id(StereoAtomId(0)),
                change: StereoAtomFieldChange::Configuration {
                    // Wrong recorded coset (Th0 vs the stored Th1).
                    old: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0),
                    ),
                    new: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1),
                    ),
                },
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_stereo_bond_field(
        mut stereo_bond_skeleton: MoleculeEditor,
    ) {
        stereo_bond_skeleton.add_stereo_bond(
            BondId(1),
            &[
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let before = stereo_bond_skeleton.clone().build();
        let tx = stereo_bond_skeleton
            .transact(Edits::from_iter([Edit::ModifyStereoBondField {
                id: StereoBondHandle::Id(StereoBondId(0)),
                change: StereoBondFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                    new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
                },
            }]))
            .unwrap();
        assert_eq!(
            stereo_bond_skeleton
                .stereo_bond(StereoBondId(0))
                .attributes()
                .configuration,
            StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
        );
        tx.rollback(&mut stereo_bond_skeleton).unwrap();
        assert_eq!(stereo_bond_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_stereo_bond_field_error(
        mut stereo_bond_skeleton: MoleculeEditor,
    ) {
        stereo_bond_skeleton.add_stereo_bond(
            BondId(1),
            &[
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let err = stereo_bond_skeleton
            .transact(Edits::from_iter([Edit::ModifyStereoBondField {
                id: StereoBondHandle::Id(StereoBondId(0)),
                change: StereoBondFieldChange::Configuration {
                    // Wrong recorded coset (vs the stored 1).
                    old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                },
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[fixture]
    fn diatomic_with_overlays() -> MoleculeEditor {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::N));
        b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
        b.add_dative_bond(&[AtomId(0)], AtomId(1), DativeBondForm::from_order(1));
        b.add_aromatic_system(&[AtomId(0), AtomId(1)], AromaticSystemForm::default());
        b.add_multicenter_bond(&[AtomId(0), AtomId(1)], MulticenterBondForm::default());
        b.add_noncovalent_bond(
            [AtomId(0), AtomId(1)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        b
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
    fn test_molecule_editor_transact_replace_dative_bond_donors(
        mut batched_overlays: MoleculeEditor,
    ) {
        let before = batched_overlays.clone().build();
        let donors = vec![AtomId(2), AtomId(4)];
        let tx = batched_overlays
            .transact(Edits::from_iter([Edit::ReplaceDativeBondDonors {
                id: DativeBondHandle::Id(DativeBondId(0)),
                old: vec![AtomHandle::Id(AtomId(0))],
                new: donors.iter().copied().map(AtomHandle::Id).collect(),
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::RestoreDativeBondDonors {
                id: DativeBondId(0),
                donors: vec![AtomId(0)],
            }]
        );
        assert_eq!(
            batched_overlays
                .dative_bond(DativeBondId(0))
                .donor_ids()
                .collect::<Vec<_>>(),
            donors
        );
        tx.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_replace_dative_bond_acceptor(
        mut batched_overlays: MoleculeEditor,
    ) {
        let before = batched_overlays.clone().build();
        let tx = batched_overlays
            .transact(Edits::from_iter([Edit::ReplaceDativeBondAcceptor {
                id: DativeBondHandle::Id(DativeBondId(0)),
                old: AtomHandle::Id(AtomId(1)),
                new: AtomHandle::Id(AtomId(3)),
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::RestoreDativeBondAcceptor {
                id: DativeBondId(0),
                acceptor: AtomId(1),
            }]
        );
        assert_eq!(
            batched_overlays.dative_bond(DativeBondId(0)).acceptor_id(),
            AtomId(3)
        );
        tx.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_replace_aromatic_system_atoms(
        mut batched_overlays: MoleculeEditor,
    ) {
        let before = batched_overlays.clone().build();
        let tx = batched_overlays
            .transact(Edits::from_iter([Edit::ReplaceAromaticSystemAtoms {
                id: AromaticSystemHandle::Id(AromaticSystemId(0)),
                old: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                new: vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::RestoreAromaticSystemAtoms {
                id: AromaticSystemId(0),
                atoms: vec![AtomId(0), AtomId(1)],
            }]
        );
        assert_eq!(
            batched_overlays
                .aromatic_system(AromaticSystemId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(0)]
        );
        tx.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_replace_multicenter_bond_atoms(
        mut batched_overlays: MoleculeEditor,
    ) {
        let before = batched_overlays.clone().build();
        let tx = batched_overlays
            .transact(Edits::from_iter([Edit::ReplaceMulticenterBondAtoms {
                id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                old: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                new: vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::RestoreMulticenterBondAtoms {
                id: MulticenterBondId(0),
                atoms: vec![AtomId(0), AtomId(1)],
            }]
        );
        assert_eq!(
            batched_overlays
                .multicenter_bond(MulticenterBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(0)]
        );
        tx.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_replace_noncovalent_bond_atoms(
        mut batched_overlays: MoleculeEditor,
    ) {
        let before = batched_overlays.clone().build();
        let tx = batched_overlays
            .transact(Edits::from_iter([Edit::ReplaceNoncovalentBondAtoms {
                id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                old: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                new: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::RestoreNoncovalentBondAtoms {
                id: NoncovalentBondId(0),
                atoms: [AtomId(0), AtomId(1)],
            }]
        );
        assert_eq!(
            batched_overlays
                .noncovalent_bond(NoncovalentBondId(0))
                .atom_ids(),
            [AtomId(1), AtomId(0)]
        );
        tx.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_replace_stereo_atom_site(
        mut batched_overlays: MoleculeEditor,
    ) {
        let before = batched_overlays.clone().build();
        let tx = batched_overlays
            .transact(Edits::from_iter([Edit::ReplaceStereoAtomSite {
                id: StereoAtomHandle::Id(StereoAtomId(0)),
                old: AtomHandle::Id(AtomId(0)),
                new: AtomHandle::Id(AtomId(2)),
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::RestoreStereoAtomSite {
                id: StereoAtomId(0),
                site: AtomId(0),
            }]
        );
        assert_eq!(
            batched_overlays.stereo_atom(StereoAtomId(0)).site_id(),
            AtomId(2)
        );
        tx.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_replace_stereo_atom_ligands(
        mut batched_overlays: MoleculeEditor,
    ) {
        let before = batched_overlays.clone().build();
        let old = vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
        ];
        let new = vec![old[2], old[0], old[1]];
        let tx = batched_overlays
            .transact(Edits::from_iter([Edit::ReplaceStereoAtomLigands {
                id: StereoAtomHandle::Id(StereoAtomId(0)),
                old: old
                    .iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                new: new
                    .iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::RestoreStereoAtomLigands {
                id: StereoAtomId(0),
                ligands: old,
            }]
        );
        assert_eq!(
            batched_overlays.stereo_atom(StereoAtomId(0)).ligand_ids(),
            new
        );
        tx.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_replace_stereo_bond_site(
        mut batched_overlays: MoleculeEditor,
    ) {
        let before = batched_overlays.clone().build();
        let tx = batched_overlays
            .transact(Edits::from_iter([Edit::ReplaceStereoBondSite {
                id: StereoBondHandle::Id(StereoBondId(0)),
                old: BondHandle::Id(BondId(0)),
                new: BondHandle::Id(BondId(1)),
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::RestoreStereoBondSite {
                id: StereoBondId(0),
                site: BondId(0),
            }]
        );
        assert_eq!(
            batched_overlays.stereo_bond(StereoBondId(0)).site_id(),
            BondId(1)
        );
        tx.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_replace_stereo_bond_ligands(
        mut batched_overlays: MoleculeEditor,
    ) {
        let before = batched_overlays.clone().build();
        let old = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
        ];
        let new = vec![old[1], old[0], old[3], old[2]];
        let tx = batched_overlays
            .transact(Edits::from_iter([Edit::ReplaceStereoBondLigands {
                id: StereoBondHandle::Id(StereoBondId(0)),
                old: old
                    .iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                new: new
                    .iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::RestoreStereoBondLigands {
                id: StereoBondId(0),
                ligands: old,
            }]
        );
        assert_eq!(
            batched_overlays.stereo_bond(StereoBondId(0)).ligand_ids(),
            new
        );
        tx.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    #[case::dative_donors(Edit::ReplaceDativeBondDonors {
        id: DativeBondHandle::Id(DativeBondId(0)),
        old: vec![],
        new: vec![AtomHandle::Id(AtomId(2))],
    })]
    #[case::dative_acceptor(Edit::ReplaceDativeBondAcceptor {
        id: DativeBondHandle::Id(DativeBondId(0)),
        old: AtomHandle::Id(AtomId(2)),
        new: AtomHandle::Id(AtomId(3)),
    })]
    #[case::aromatic_atoms(Edit::ReplaceAromaticSystemAtoms {
        id: AromaticSystemHandle::Id(AromaticSystemId(0)),
        old: vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
        new: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(2))],
    })]
    #[case::multicenter_atoms(Edit::ReplaceMulticenterBondAtoms {
        id: MulticenterBondHandle::Id(MulticenterBondId(0)),
        old: vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
        new: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(2))],
    })]
    #[case::noncovalent_atoms(Edit::ReplaceNoncovalentBondAtoms {
        id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
        old: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
        new: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(2))],
    })]
    #[case::stereo_atom_site(Edit::ReplaceStereoAtomSite {
        id: StereoAtomHandle::Id(StereoAtomId(0)),
        old: AtomHandle::Id(AtomId(2)),
        new: AtomHandle::Id(AtomId(3)),
    })]
    #[case::stereo_atom_ligands(Edit::ReplaceStereoAtomLigands {
        id: StereoAtomHandle::Id(StereoAtomId(0)),
        old: vec![],
        new: vec![(AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom)],
    })]
    #[case::stereo_bond_site(Edit::ReplaceStereoBondSite {
        id: StereoBondHandle::Id(StereoBondId(0)),
        old: BondHandle::Id(BondId(1)),
        new: BondHandle::Id(BondId(2)),
    })]
    #[case::stereo_bond_ligands(Edit::ReplaceStereoBondLigands {
        id: StereoBondHandle::Id(StereoBondId(0)),
        old: vec![],
        new: vec![(AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom)],
    })]
    fn test_molecule_editor_transact_replace_old_state_error(
        mut batched_overlays: MoleculeEditor,
        #[case] edit: Edit,
    ) {
        let before = batched_overlays.clone().build();
        assert_eq!(
            batched_overlays.transact(Edits::from_iter([edit])),
            Err(TransactionError::OldStateMismatch),
        );
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    #[case::dative_donors(Undo::RestoreDativeBondDonors {
        id: DativeBondId(9), donors: vec![AtomId(0)],
    })]
    #[case::dative_acceptor(Undo::RestoreDativeBondAcceptor {
        id: DativeBondId(9), acceptor: AtomId(0),
    })]
    #[case::aromatic_atoms(Undo::RestoreAromaticSystemAtoms {
        id: AromaticSystemId(9), atoms: vec![AtomId(0)],
    })]
    #[case::multicenter_atoms(Undo::RestoreMulticenterBondAtoms {
        id: MulticenterBondId(9), atoms: vec![AtomId(0)],
    })]
    #[case::noncovalent_atoms(Undo::RestoreNoncovalentBondAtoms {
        id: NoncovalentBondId(9), atoms: [AtomId(0), AtomId(1)],
    })]
    #[case::stereo_atom_site(Undo::RestoreStereoAtomSite {
        id: StereoAtomId(9), site: AtomId(0),
    })]
    #[case::stereo_atom_ligands(Undo::RestoreStereoAtomLigands {
        id: StereoAtomId(9), ligands: vec![],
    })]
    #[case::stereo_bond_site(Undo::RestoreStereoBondSite {
        id: StereoBondId(9), site: BondId(0),
    })]
    #[case::stereo_bond_ligands(Undo::RestoreStereoBondLigands {
        id: StereoBondId(9), ligands: vec![],
    })]
    fn test_transaction_rollback_restore_target_error(
        mut batched_overlays: MoleculeEditor,
        #[case] undo: Undo,
    ) {
        let before = batched_overlays.clone().build();
        assert_eq!(
            (Transaction { undo: vec![undo] }).rollback(&mut batched_overlays),
            Err(TransactionError::RollbackStateMismatch),
        );
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    #[case::dative_first(EntityKind::DativeBond, 0)]
    #[case::dative_middle(EntityKind::DativeBond, 1)]
    #[case::dative_last(EntityKind::DativeBond, 2)]
    #[case::aromatic_first(EntityKind::AromaticSystem, 0)]
    #[case::aromatic_middle(EntityKind::AromaticSystem, 1)]
    #[case::aromatic_last(EntityKind::AromaticSystem, 2)]
    #[case::multicenter_first(EntityKind::MulticenterBond, 0)]
    #[case::multicenter_middle(EntityKind::MulticenterBond, 1)]
    #[case::multicenter_last(EntityKind::MulticenterBond, 2)]
    #[case::noncovalent_first(EntityKind::NoncovalentBond, 0)]
    #[case::noncovalent_middle(EntityKind::NoncovalentBond, 1)]
    #[case::noncovalent_last(EntityKind::NoncovalentBond, 2)]
    #[case::stereo_atom_first(EntityKind::StereoAtom, 0)]
    #[case::stereo_atom_middle(EntityKind::StereoAtom, 1)]
    #[case::stereo_atom_last(EntityKind::StereoAtom, 2)]
    #[case::stereo_bond_first(EntityKind::StereoBond, 0)]
    #[case::stereo_bond_middle(EntityKind::StereoBond, 1)]
    #[case::stereo_bond_last(EntityKind::StereoBond, 2)]
    fn test_molecule_editor_transact_remove_overlays_error(
        mut batched_overlays: MoleculeEditor,
        #[case] kind: EntityKind,
        #[case] invalid_position: usize,
    ) {
        let before = batched_overlays.clone().build();
        let edit = match kind {
            EntityKind::DativeBond => Edit::RemoveDativeBonds {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            DativeBondHandle::Id(DativeBondId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            vec![AtomHandle::Id(AtomId(index * 2))],
                            AtomHandle::Id(AtomId(index * 2 + 1)),
                            DativeBondForm::from_order(1),
                        )
                    })
                    .collect(),
            },
            EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            AromaticSystemHandle::Id(AromaticSystemId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            vec![
                                AtomHandle::Id(AtomId(index * 2)),
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                            ],
                            AromaticSystemForm::default(),
                        )
                    })
                    .collect(),
            },
            EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            MulticenterBondHandle::Id(MulticenterBondId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            vec![
                                AtomHandle::Id(AtomId(index * 2)),
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                            ],
                            MulticenterBondForm::default(),
                        )
                    })
                    .collect(),
            },
            EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            NoncovalentBondHandle::Id(NoncovalentBondId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            [
                                AtomHandle::Id(AtomId(index * 2)),
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                            ],
                            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                        )
                    })
                    .collect(),
            },
            EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            StereoAtomHandle::Id(StereoAtomId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            AtomHandle::Id(AtomId(index * 2)),
                            vec![
                                (
                                    AtomHandle::Id(AtomId(index * 2 + 1)),
                                    StereoLigandKind::Atom,
                                ),
                                (
                                    AtomHandle::Id(AtomId(index * 2)),
                                    StereoLigandKind::ImplicitHydrogen,
                                ),
                                (
                                    AtomHandle::Id(AtomId(index * 2)),
                                    StereoLigandKind::LonePair,
                                ),
                            ],
                            StereoAtomForm::default(),
                        )
                    })
                    .collect(),
            },
            EntityKind::StereoBond => Edit::RemoveStereoBonds {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            StereoBondHandle::Id(StereoBondId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            BondHandle::Id(BondId(index)),
                            vec![
                                (
                                    AtomHandle::Id(AtomId(index * 2)),
                                    StereoLigandKind::ImplicitHydrogen,
                                ),
                                (
                                    AtomHandle::Id(AtomId(index * 2)),
                                    StereoLigandKind::LonePair,
                                ),
                                (
                                    AtomHandle::Id(AtomId(index * 2 + 1)),
                                    StereoLigandKind::ImplicitHydrogen,
                                ),
                                (
                                    AtomHandle::Id(AtomId(index * 2 + 1)),
                                    StereoLigandKind::LonePair,
                                ),
                            ],
                            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                        )
                    })
                    .collect(),
            },
            EntityKind::Atom | EntityKind::Bond => unreachable!(),
        };

        assert_eq!(
            batched_overlays
                .transact(Edits::from_iter([edit]))
                .unwrap_err(),
            TransactionError::HandleOutOfRange {
                kind,
                index: 9,
                count: 3,
            }
        );
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    #[case::atom(EntityKind::Atom)]
    #[case::bond(EntityKind::Bond)]
    #[case::dative_bond(EntityKind::DativeBond)]
    #[case::aromatic_system(EntityKind::AromaticSystem)]
    #[case::multicenter_bond(EntityKind::MulticenterBond)]
    #[case::noncovalent_bond(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_molecule_editor_transact_duplicate_removal_error(
        mut batched_overlays: MoleculeEditor,
        #[case] kind: EntityKind,
    ) {
        let before = batched_overlays.clone().build();
        let edit = match kind {
            EntityKind::Atom => Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            EntityKind::Bond => Edit::RemoveTopology {
                atoms: Vec::new(),
                bonds: vec![BondHandle::Id(BondId(0)), BondHandle::Id(BondId(0))],
            },
            EntityKind::DativeBond => Edit::RemoveDativeBonds {
                removes: vec![
                    (
                        DativeBondHandle::Id(DativeBondId(0)),
                        vec![AtomHandle::Id(AtomId(0))],
                        AtomHandle::Id(AtomId(1)),
                        DativeBondForm::from_order(1),
                    ),
                    (
                        DativeBondHandle::Id(DativeBondId(0)),
                        vec![AtomHandle::Id(AtomId(0))],
                        AtomHandle::Id(AtomId(1)),
                        DativeBondForm::from_order(1),
                    ),
                ],
            },
            EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
                removes: vec![
                    (
                        AromaticSystemHandle::Id(AromaticSystemId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        AromaticSystemForm::default(),
                    ),
                    (
                        AromaticSystemHandle::Id(AromaticSystemId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        AromaticSystemForm::default(),
                    ),
                ],
            },
            EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
                removes: vec![
                    (
                        MulticenterBondHandle::Id(MulticenterBondId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        MulticenterBondForm::default(),
                    ),
                    (
                        MulticenterBondHandle::Id(MulticenterBondId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        MulticenterBondForm::default(),
                    ),
                ],
            },
            EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
                removes: vec![
                    (
                        NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                        [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    ),
                    (
                        NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                        [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    ),
                ],
            },
            EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
                removes: vec![
                    (
                        StereoAtomHandle::Id(StereoAtomId(0)),
                        AtomHandle::Id(AtomId(0)),
                        vec![
                            (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                            (
                                AtomHandle::Id(AtomId(0)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                        ],
                        StereoAtomForm::default(),
                    ),
                    (
                        StereoAtomHandle::Id(StereoAtomId(0)),
                        AtomHandle::Id(AtomId(0)),
                        vec![
                            (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                            (
                                AtomHandle::Id(AtomId(0)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                        ],
                        StereoAtomForm::default(),
                    ),
                ],
            },
            EntityKind::StereoBond => Edit::RemoveStereoBonds {
                removes: vec![
                    (
                        StereoBondHandle::Id(StereoBondId(0)),
                        BondHandle::Id(BondId(0)),
                        vec![
                            (
                                AtomHandle::Id(AtomId(0)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                            (
                                AtomHandle::Id(AtomId(1)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (AtomHandle::Id(AtomId(1)), StereoLigandKind::LonePair),
                        ],
                        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                    ),
                    (
                        StereoBondHandle::Id(StereoBondId(0)),
                        BondHandle::Id(BondId(0)),
                        vec![
                            (
                                AtomHandle::Id(AtomId(0)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                            (
                                AtomHandle::Id(AtomId(1)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (AtomHandle::Id(AtomId(1)), StereoLigandKind::LonePair),
                        ],
                        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                    ),
                ],
            },
        };

        assert_eq!(
            batched_overlays
                .transact(Edits::from_iter([edit]))
                .unwrap_err(),
            TransactionError::DuplicateRemoval { kind }
        );
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    #[case::dative_bond(EntityKind::DativeBond)]
    #[case::aromatic_system(EntityKind::AromaticSystem)]
    #[case::multicenter_bond(EntityKind::MulticenterBond)]
    #[case::noncovalent_bond(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_molecule_editor_transact_handle_removed_error_cascade(
        mut batched_overlays: MoleculeEditor,
        #[case] kind: EntityKind,
    ) {
        let before = batched_overlays.clone().build();
        let mut edits = Edits::from_iter([Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: Vec::new(),
        }]);
        edits.push(match kind {
            EntityKind::DativeBond => Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(0))],
                    AtomHandle::Id(AtomId(1)),
                    DativeBondForm::from_order(1),
                )],
            },
            EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    AromaticSystemForm::default(),
                )],
            },
            EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
                removes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    MulticenterBondForm::default(),
                )],
            },
            EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
            },
            EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    vec![(AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom)],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                )],
            },
            EntityKind::StereoBond => Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                )],
            },
            EntityKind::Atom | EntityKind::Bond => unreachable!(),
        });

        assert_eq!(
            batched_overlays.transact(edits).unwrap_err(),
            TransactionError::HandleRemoved { kind, index: 0 }
        );
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_dative_bond_field(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyDativeBondField {
                id: DativeBondHandle::Id(DativeBondId(0)),
                change: DativeBondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .dative_bond(DativeBondId(0))
                .attributes()
                .order,
            NumForm::Lit(2),
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_aromatic_system_field(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyAromaticSystemField {
                id: AromaticSystemHandle::Id(AromaticSystemId(0)),
                change: AromaticSystemFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .aromatic_system(AromaticSystemId(0))
                .attributes()
                .charge,
            NumForm::Lit(1),
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_multicenter_bond_field(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                change: MulticenterBondFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(-1),
                },
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .multicenter_bond(MulticenterBondId(0))
                .attributes()
                .charge,
            NumForm::Lit(-1),
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_noncovalent_bond_field(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyNoncovalentBondField {
                id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                    new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                },
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .noncovalent_bond(NoncovalentBondId(0))
                .attributes()
                .kind,
            NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_add_dative_bond(mut diatomic: MoleculeEditor) {
        let tx = diatomic
            .transact(Edits::from_iter([Edit::AddDativeBond {
                donors: vec![AtomHandle::Id(AtomId(0))],
                acceptor: AtomHandle::Id(AtomId(1)),
                attributes: DativeBondForm::from_order(1),
            }]))
            .unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedDativeBond(added)] if added.id == DativeBondId(0)
        ));
        assert_eq!(diatomic.dative_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_editor_transact_add_aromatic_system(mut diatomic: MoleculeEditor) {
        let tx = diatomic
            .transact(Edits::from_iter([Edit::AddAromaticSystem {
                atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                attributes: AromaticSystemForm::default(),
            }]))
            .unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedAromaticSystem(added)] if added.id == AromaticSystemId(0)
        ));
        assert_eq!(diatomic.aromatic_system_count(), 1);
    }

    #[rstest]
    fn test_molecule_editor_transact_add_multicenter_bond(mut diatomic: MoleculeEditor) {
        let tx = diatomic
            .transact(Edits::from_iter([Edit::AddMulticenterBond {
                atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                attributes: MulticenterBondForm::default(),
            }]))
            .unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedMulticenterBond(added)] if added.id == MulticenterBondId(0)
        ));
        assert_eq!(diatomic.multicenter_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_editor_transact_add_noncovalent_bond(mut diatomic: MoleculeEditor) {
        let tx = diatomic
            .transact(Edits::from_iter([Edit::AddNoncovalentBond {
                atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            }]))
            .unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedNoncovalentBond(added)] if added.id == NoncovalentBondId(0)
        ));
        assert_eq!(diatomic.noncovalent_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_dative_bond(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let before = diatomic_with_overlays.clone().build();
        let transaction = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(0))],
                    AtomHandle::Id(AtomId(1)),
                    DativeBondForm {
                        order: NumForm::Lit(1),
                        constraints: Default::default(),
                    },
                )],
            }]))
            .unwrap();
        assert_eq!(diatomic_with_overlays.dative_bond_count(), 0);
        transaction.rollback(&mut diatomic_with_overlays).unwrap();
        assert_eq!(diatomic_with_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_dative_bond_roles_error(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let err = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(1))],
                    AtomHandle::Id(AtomId(0)),
                    DativeBondForm {
                        order: NumForm::Lit(1),
                        constraints: Default::default(),
                    },
                )],
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
        assert_eq!(diatomic_with_overlays.dative_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_aromatic_system(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let before = diatomic_with_overlays.clone().build();
        let transaction = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    AromaticSystemForm::default(),
                )],
            }]))
            .unwrap();
        assert_eq!(diatomic_with_overlays.aromatic_system_count(), 0);
        transaction.rollback(&mut diatomic_with_overlays).unwrap();
        assert_eq!(diatomic_with_overlays.build(), before);
    }

    // Batch removal of non-contiguous same-kind ids (0 and 2) in one edit: ids resolve against the
    // pre-removal state and compact once, so the survivor (former id 1) remaps to id 0. A single-id
    // sequence would stale id 2 after removing id 0.
    #[rstest]
    fn test_molecule_editor_transact_remove_aromatic_systems() {
        let mut b = Molecule::default().edit();
        for _ in 0..6 {
            b.add_atom(AtomForm::from_element(Element::C));
        }
        b.add_aromatic_system(&[AtomId(0), AtomId(1)], AromaticSystemForm::default());
        b.add_aromatic_system(&[AtomId(2), AtomId(3)], AromaticSystemForm::default());
        b.add_aromatic_system(&[AtomId(4), AtomId(5)], AromaticSystemForm::default());
        b.transact(Edits::from_iter([Edit::RemoveAromaticSystems {
            removes: vec![
                (
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    AromaticSystemForm::default(),
                ),
                (
                    AromaticSystemHandle::Id(AromaticSystemId(2)),
                    vec![AtomHandle::Id(AtomId(4)), AtomHandle::Id(AtomId(5))],
                    AromaticSystemForm::default(),
                ),
            ],
        }]))
        .unwrap();
        assert_eq!(b.aromatic_system_count(), 1);
        assert_eq!(
            b.aromatic_system(AromaticSystemId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2), AtomId(3)],
        );
    }

    // Rolling back an aromatic-system removal restores a molecule constraint the removal dropped
    // (`dropped`) or remapped (`remapped`) — the overlay-remove undo captures the constraint cascade.
    #[rstest]
    #[case::dropped(AromaticSystemId(0), 0)]
    #[case::remapped(AromaticSystemId(1), 1)]
    fn test_molecule_editor_transact_remove_aromatic_system_rollback(
        #[case] constrained: AromaticSystemId,
        #[case] forward_constraint_count: usize,
    ) {
        let mut b = Molecule::default().edit();
        for _ in 0..6 {
            b.add_atom(AtomForm::from_element(Element::C));
        }
        b.add_aromatic_system(
            &[AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::default(),
        );
        b.add_aromatic_system(
            &[AtomId(3), AtomId(4), AtomId(5)],
            AromaticSystemForm::default(),
        );
        b.transact(Edits::from_iter([Edit::AddMoleculeConstraint {
            constraint: Constraint::AromaticSystem(
                constrained,
                AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6)),
            )
            .into(),
        }]))
        .unwrap();
        let before = b.clone().build();

        let tx = b
            .transact(Edits::from_iter([Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![
                        AtomHandle::Id(AtomId(0)),
                        AtomHandle::Id(AtomId(1)),
                        AtomHandle::Id(AtomId(2)),
                    ],
                    AromaticSystemForm::default(),
                )],
            }]))
            .unwrap();
        assert_eq!(b.aromatic_system_count(), 1);
        assert_eq!(b.constraints().iter().count(), forward_constraint_count);

        tx.rollback(&mut b).unwrap();
        assert_eq!(b.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_multicenter_bond(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let before = diatomic_with_overlays.clone().build();
        let transaction = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveMulticenterBonds {
                removes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    MulticenterBondForm::default(),
                )],
            }]))
            .unwrap();
        assert_eq!(diatomic_with_overlays.multicenter_bond_count(), 0);
        transaction.rollback(&mut diatomic_with_overlays).unwrap();
        assert_eq!(diatomic_with_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_noncovalent_bond(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let before = diatomic_with_overlays.clone().build();
        let transaction = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
            }]))
            .unwrap();
        assert_eq!(diatomic_with_overlays.noncovalent_bond_count(), 0);
        transaction.rollback(&mut diatomic_with_overlays).unwrap();
        assert_eq!(diatomic_with_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_noncovalent_bond_form_mismatch_error(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let err = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic), // wrong
                )],
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_bond_constraint_value_bearing(
        mut diatomic: MoleculeEditor,
    ) {
        diatomic
            .transact(Edits::from_iter([Edit::ModifyBondConstraint {
                id: BondHandle::Id(BondId(0)),
                old: None,
                new: Some(BondConstraintForm::cis_trans_stereo(
                    CisTransStereoForm::NotStereo,
                )),
            }]))
            .unwrap();
        assert_eq!(
            diatomic
                .bond(BondId(0))
                .attributes()
                .constraints
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            vec![BondConstraintForm::cis_trans_stereo(
                CisTransStereoForm::NotStereo
            )],
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_add_bond_constraint(mut diatomic: MoleculeEditor) {
        diatomic
            .transact(Edits::from_iter([Edit::ModifyBondConstraint {
                id: BondHandle::Id(BondId(0)),
                old: None,
                new: Some(BondConstraintForm::ring_membership(RingScope::Size(5), 1)),
            }]))
            .unwrap();
        assert!(diatomic
            .bond(BondId(0))
            .attributes()
            .constraints
            .iter()
            .any(|c| *c == BondConstraintForm::ring_membership(RingScope::Size(5), 1)));
    }

    #[rstest]
    fn test_molecule_editor_transact_modify_bond_constraint_absent_error(
        mut diatomic: MoleculeEditor,
    ) {
        let err = diatomic
            .transact(Edits::from_iter([Edit::ModifyBondConstraint {
                id: BondHandle::Id(BondId(0)),
                old: Some(BondConstraintForm::ring_membership(RingScope::Size(5), 1)),
                new: None,
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_dative_bond_constraint(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyDativeBondConstraint {
                id: DativeBondHandle::Id(DativeBondId(0)),
                old: None,
                new: Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
            }]))
            .unwrap();
        assert!(diatomic_with_overlays
            .dative_bond(DativeBondId(0))
            .attributes()
            .constraints
            .iter()
            .any(|c| *c == DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))));
    }

    #[rstest]
    fn test_molecule_editor_transact_set_aromatic_system_constraint(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyAromaticSystemConstraint {
                id: AromaticSystemHandle::Id(AromaticSystemId(0)),
                old: None,
                new: Some(AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))),
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .aromatic_system(AromaticSystemId(0))
                .attributes()
                .constraints
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            vec![AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))],
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_multicenter_bond_constraint(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyMulticenterBondConstraint {
                id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                old: None,
                new: Some(MulticenterBondConstraintForm::ElectronCount(NumForm::Lit(
                    2,
                ))),
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .multicenter_bond(MulticenterBondId(0))
                .attributes()
                .constraints
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            vec![MulticenterBondConstraintForm::ElectronCount(NumForm::Lit(
                2
            ))],
        );
    }

    #[fixture]
    fn triatomic_with_overlays() -> MoleculeEditor {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::N));
        b.add_atom(AtomForm::from_element(Element::O));
        b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
        b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(1));
        b.add_dative_bond(&[AtomId(0)], AtomId(1), DativeBondForm::from_order(1));
        b.add_aromatic_system(
            &[AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::default(),
        );
        b.add_multicenter_bond(
            &[AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        );
        b.add_noncovalent_bond(
            [AtomId(0), AtomId(2)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        b
    }

    #[rstest]
    fn test_transaction_append(mut diatomic: MoleculeEditor) {
        let before = diatomic.clone().build();
        let first = diatomic
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }]))
            .unwrap();
        let second = diatomic
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            }]))
            .unwrap();
        let expected_undos = [first.undos(), second.undos()].concat();

        let mut combined = Transaction::default();
        combined.append(first);
        combined.append(second);

        assert_eq!(combined.undos(), expected_undos);
        combined.rollback(&mut diatomic).unwrap();
        assert_eq!(diatomic.build(), before);
    }

    #[rstest]
    fn test_transaction_append_error(mut diatomic: MoleculeEditor) {
        let before = diatomic.clone().build();
        let first = diatomic
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }]))
            .unwrap();
        let second = diatomic
            .transact(Edits::from_iter([Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old: None,
                new: Some(AtomConstraintForm::degree(1)),
            }]))
            .unwrap();
        let third = diatomic
            .transact(Edits::from_iter([Edit::AddDativeBond {
                donors: vec![AtomHandle::Id(AtomId(0))],
                acceptor: AtomHandle::Id(AtomId(1)),
                attributes: DativeBondForm::from_order(1),
            }]))
            .unwrap();
        let expected_undos = [first.undos(), second.undos(), third.undos()].concat();

        let mut combined = Transaction::default();
        combined.append(first);
        combined.append(second);
        combined.append(third);

        let materialized = diatomic.clone().build();
        let mut rejected = Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        }]);
        rejected.remove_atom(AtomHandle::Id(AtomId(99)));
        let error = diatomic.transact(rejected).unwrap_err();
        assert_eq!(
            error,
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 99,
                count: 2,
            }
        );
        assert_eq!(diatomic.clone().build(), materialized);
        assert_eq!(combined.undos(), expected_undos);

        combined.rollback(&mut diatomic).unwrap();
        assert_eq!(diatomic.build(), before);
    }

    enum RollbackCase {
        RemoveTopology,
        RemoveBond,
        AddTopology,
        Field,
        AddOverlay,
        RemoveOverlay,
        Constraint,
        CascadedConstraints,
    }

    #[rstest]
    #[case::remove_topology(RollbackCase::RemoveTopology)]
    #[case::remove_bond(RollbackCase::RemoveBond)]
    #[case::add_topology(RollbackCase::AddTopology)]
    #[case::field(RollbackCase::Field)]
    #[case::add_overlay(RollbackCase::AddOverlay)]
    #[case::remove_overlay(RollbackCase::RemoveOverlay)]
    #[case::constraint(RollbackCase::Constraint)]
    #[case::cascade(RollbackCase::CascadedConstraints)]
    fn test_transaction_rollback(#[case] case: RollbackCase) {
        let mut editor = match case {
            RollbackCase::AddTopology => Molecule::default().edit(),
            RollbackCase::Field | RollbackCase::Constraint => {
                let mut b = Molecule::default().edit();
                b.add_atom(AtomForm::from_element(Element::C));
                b
            }
            RollbackCase::AddOverlay => {
                let mut b = Molecule::default().edit();
                b.add_atom(AtomForm::from_element(Element::C));
                b.add_atom(AtomForm::from_element(Element::C));
                b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
                b
            }
            RollbackCase::RemoveBond => {
                let mut b = Molecule::default().edit();
                b.add_atom(AtomForm::from_element(Element::C));
                b.add_atom(AtomForm::from_element(Element::N));
                b.add_atom(AtomForm::from_element(Element::O));
                b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
                b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(2));
                b
            }
            RollbackCase::CascadedConstraints => {
                let mut b = Molecule::default().edit();
                b.add_atom(AtomForm::from_element(Element::C));
                b.add_atom(AtomForm::from_element(Element::N));
                b.constraints_mut()
                    .push(Constraint::Atom(AtomId(1), AtomConstraintForm::degree(3)));
                b
            }
            RollbackCase::RemoveTopology | RollbackCase::RemoveOverlay => triatomic_with_overlays(),
        };
        let before = editor.clone().build();
        let edits = match case {
            RollbackCase::RemoveTopology => Edits::from_iter([Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(1))],
                bonds: vec![],
            }]),
            RollbackCase::RemoveBond => Edits::from_iter([Edit::RemoveTopology {
                atoms: Vec::new(),
                bonds: vec![BondHandle::Id(BondId(0))],
            }]),
            RollbackCase::AddTopology => Edits::from_iter([
                Edit::AddAtoms {
                    atoms: vec![
                        AtomForm::from_element(Element::C),
                        AtomForm::from_element(Element::O),
                    ],
                },
                Edit::AddBonds {
                    bonds: vec![AddBond {
                        endpoints: [AtomHandle::New(0), AtomHandle::New(1)],
                        attributes: BondForm::from_order(2),
                    }],
                },
            ]),
            RollbackCase::Field => Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }]),
            RollbackCase::AddOverlay => Edits::from_iter([Edit::AddDativeBond {
                donors: vec![AtomHandle::Id(AtomId(0))],
                acceptor: AtomHandle::Id(AtomId(1)),
                attributes: DativeBondForm::from_order(1),
            }]),
            RollbackCase::RemoveOverlay => Edits::from_iter([Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(0))],
                    AtomHandle::Id(AtomId(1)),
                    DativeBondForm {
                        order: NumForm::Lit(1),
                        constraints: Default::default(),
                    },
                )],
            }]),
            RollbackCase::Constraint => Edits::from_iter([Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old: None,
                new: Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 1)),
            }]),
            RollbackCase::CascadedConstraints => Edits::from_iter([Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(1))],
                bonds: Vec::new(),
            }]),
        };
        let tx = editor.transact(edits).unwrap();
        tx.rollback(&mut editor).unwrap();
        assert_eq!(editor.build(), before);
    }

    #[rstest]
    #[case::first(vec![AtomHandle::Id(AtomId(0))], 2)]
    #[case::middle(vec![AtomHandle::Id(AtomId(2))], 2)]
    #[case::separated(vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(4))], 1)]
    fn test_transaction_rollback_topology(
        mut batched_overlays: MoleculeEditor,
        #[case] atoms: Vec<AtomHandle>,
        #[case] remaining: usize,
    ) {
        batched_overlays.constraints_mut().extend(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1)),
            Constraint::Atom(AtomId(5), AtomConstraintForm::degree(1)),
            Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1)),
        ]);
        let before = batched_overlays.clone().build();
        let transaction = batched_overlays
            .transact(Edits::from_iter([Edit::RemoveTopology {
                atoms,
                bonds: vec![],
            }]))
            .unwrap();
        assert_eq!(
            [
                batched_overlays.dative_bond_count(),
                batched_overlays.aromatic_system_count(),
                batched_overlays.multicenter_bond_count(),
                batched_overlays.noncovalent_bond_count(),
                batched_overlays.stereo_atom_count(),
                batched_overlays.stereo_bond_count(),
            ],
            [remaining; 6],
        );
        transaction.rollback(&mut batched_overlays).unwrap();
        let restored = batched_overlays.build();
        assert!(restored.normalized_eq(&before));
        assert_eq!(restored, before);
    }

    #[rstest]
    #[case::missing(7)]
    #[case::overlapping(0)]
    fn test_transaction_rollback_constraint_history(
        mut one_atom: MoleculeEditor,
        #[case] position: usize,
    ) {
        let constraint = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1));
        one_atom.constraints_mut().push(constraint.clone());
        let transaction = Transaction {
            undo: vec![Undo::RestoreMoleculeConstraints(CascadedConstraints {
                removed: vec![
                    RemovedConstraint {
                        position,
                        constraint: constraint.clone(),
                    },
                    RemovedConstraint {
                        position,
                        constraint: constraint.clone(),
                    },
                ],
                modified: vec![ModifiedConstraint {
                    position,
                    old: constraint,
                    new: Constraint::Or(vec![]),
                }],
            })],
        };
        let _ = transaction.rollback(&mut one_atom);
    }

    #[rstest]
    fn test_transaction_rollback_empty(mut one_atom: MoleculeEditor) {
        let before = one_atom.clone().build();
        Transaction::default().rollback(&mut one_atom).unwrap();
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    #[case::atom(EntityKind::Atom)]
    #[case::bond(EntityKind::Bond)]
    #[case::dative_bond(EntityKind::DativeBond)]
    #[case::aromatic_system(EntityKind::AromaticSystem)]
    #[case::multicenter_bond(EntityKind::MulticenterBond)]
    #[case::noncovalent_bond(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_transaction_rollback_field_receiver(
        mut empty: MoleculeEditor,
        #[case] kind: EntityKind,
    ) {
        let undo = match kind {
            EntityKind::Atom => Undo::ModifyAtomField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::default(),
                },
            },
            EntityKind::Bond => Undo::ModifyBondField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(2),
                    new: NumForm::Lit(1),
                },
            },
            EntityKind::DativeBond => Undo::ModifyDativeBondField {
                id: DativeBondId(0),
                change: DativeBondFieldChange::Order {
                    old: NumForm::Lit(2),
                    new: NumForm::Lit(1),
                },
            },
            EntityKind::AromaticSystem => Undo::ModifyAromaticSystemField {
                id: AromaticSystemId(0),
                change: AromaticSystemFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::default(),
                },
            },
            EntityKind::MulticenterBond => Undo::ModifyMulticenterBondField {
                id: MulticenterBondId(0),
                change: MulticenterBondFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::default(),
                },
            },
            EntityKind::NoncovalentBond => Undo::ModifyNoncovalentBondField {
                id: NoncovalentBondId(0),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                    new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                },
            },
            EntityKind::StereoAtom => Undo::ModifyStereoAtomField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0),
                    ),
                    new: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1),
                    ),
                },
            },
            EntityKind::StereoBond => Undo::ModifyStereoBondField {
                id: StereoBondId(0),
                change: StereoBondFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                },
            },
        };

        assert_eq!(
            (Transaction { undo: vec![undo] }).rollback(&mut empty),
            Err(TransactionError::RollbackStateMismatch),
        );
    }

    #[rstest]
    fn test_transaction_rollback_added_topology_duplicate(mut one_atom: MoleculeEditor) {
        let before = one_atom.clone().build();
        let added = AddedAtom {
            id: AtomId(0),
            attributes: AtomForm::from_element(Element::C),
        };
        let transaction = Transaction {
            undo: vec![Undo::RemoveAddedTopology {
                atoms: vec![added.clone(), added],
                bonds: Vec::new(),
            }],
        };

        assert_eq!(
            transaction.rollback(&mut one_atom),
            Err(TransactionError::RollbackStateMismatch),
        );
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    fn test_transaction_rollback_reconstruction_entry(mut empty: MoleculeEditor) {
        let transaction = Transaction {
            undo: vec![Undo::RestoreRemovedAromaticSystems {
                removed: vec![RemovedAromaticSystem {
                    id: AromaticSystemId(1),
                    atoms: Vec::new(),
                    attributes: AromaticSystemForm::default(),
                }],
                undo_compaction: MoleculeCompaction::empty().undo_compaction(),
                cascade: CascadedConstraints::default(),
            }],
        };

        assert_eq!(
            transaction.rollback(&mut empty),
            Err(TransactionError::RollbackStateMismatch),
        );
    }

    #[rstest]
    #[case::atoms_short([0, 0, 0, 0, 0, 0, 0, 0])]
    #[case::atoms_long([2, 0, 0, 0, 0, 0, 0, 0])]
    #[case::bonds([1, 1, 0, 0, 0, 0, 0, 0])]
    #[case::dative([1, 0, 1, 0, 0, 0, 0, 0])]
    #[case::aromatic([1, 0, 0, 1, 0, 0, 0, 0])]
    #[case::multicenter([1, 0, 0, 0, 1, 0, 0, 0])]
    #[case::noncovalent([1, 0, 0, 0, 0, 1, 0, 0])]
    #[case::stereo_atoms([1, 0, 0, 0, 0, 0, 1, 0])]
    #[case::stereo_bonds([1, 0, 0, 0, 0, 0, 0, 1])]
    fn test_transaction_rollback_compaction_counts(
        mut one_atom: MoleculeEditor,
        #[case] counts: [usize; 8],
    ) {
        let before = one_atom.clone().build();
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(counts[0]),
                Compaction::identity(counts[1]),
            ),
            Compaction::identity(counts[2]),
            Compaction::identity(counts[3]),
            Compaction::identity(counts[4]),
            Compaction::identity(counts[5]),
            Compaction::identity(counts[6]),
            Compaction::identity(counts[7]),
        );
        let transaction = Transaction {
            undo: vec![Undo::RestoreRemovedTopology {
                atoms: Vec::new(),
                bonds: Vec::new(),
                overlays: RemovedOverlays::default(),
                undo_compaction: compaction.undo_compaction(),
                compaction,
                cascade: CascadedConstraints::default(),
            }],
        };
        assert_eq!(
            transaction.rollback(&mut one_atom),
            Err(TransactionError::RollbackStateMismatch)
        );
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    fn test_transaction_rollback_compaction_dimension(mut one_atom: MoleculeEditor) {
        let compaction = MoleculeCompaction::empty();
        let transaction = Transaction {
            undo: vec![Undo::RestoreRemovedTopology {
                atoms: vec![RemovedAtom {
                    id: AtomId(0),
                    attributes: AtomForm::from_element(Element::N),
                }],
                bonds: Vec::new(),
                overlays: RemovedOverlays::default(),
                undo_compaction: compaction.undo_compaction(),
                compaction,
                cascade: CascadedConstraints::default(),
            }],
        };

        assert_eq!(
            transaction.rollback(&mut one_atom),
            Err(TransactionError::RollbackStateMismatch),
        );
    }

    #[rstest]
    fn test_transaction_rollback_molecule_constraint_order(mut one_atom: MoleculeEditor) {
        let repeated = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1));
        let middle = Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4));
        one_atom.constraints_mut().push(repeated.clone());
        one_atom.constraints_mut().push(middle.clone());
        one_atom.constraints_mut().push(repeated.clone());
        let before = one_atom.clone().build();

        let transaction = one_atom
            .transact(Edits::from_iter([Edit::RemoveMoleculeConstraint {
                constraint: repeated.clone().into(),
            }]))
            .unwrap();
        assert_eq!(one_atom.constraints().as_slice(), &[repeated, middle]);

        transaction.rollback(&mut one_atom).unwrap();
        assert_eq!(one_atom.build(), before);
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
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 0,
                count: 0,
            }
        );
    }
}
