//! Scoped molecule transactions with private undo-based restoration.

use std::mem;

use thiserror::Error;
use umol_graph_core::Correspondence;

use super::apply::ApplicationState;
use super::{Molecule, MoleculeIntegrityError};
use crate::ir::correspondence::MoleculeCorrespondence;
use crate::ir::edit::{Edits, Undo};
use crate::ir::entity::EntityKind;
use crate::ir::error::MoleculeApplyError;

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

    /// The current state does not match the edit's old value.
    #[error("precondition failed: old state does not match current")]
    OldStateMismatch,

    /// `Remove*Constraint` with a value that's not present.
    #[error("missing constraint entry on remove")]
    MissingEntry,

    /// Edit shape is structurally invalid.
    #[error("malformed edit: {0}")]
    MalformedEdit(&'static str),

    /// An earlier application or commit failure aborted the transaction.
    #[error("transaction has been aborted")]
    Aborted,
}

/// Scoped access to a molecule and its private undo journal.
///
/// Transaction::run owns the restoration guard and lends this handle to its callback.
/// Completed changes are accepted only when commit succeeds and the callback returns Ok.
/// Otherwise the guard restores the molecule, including when the handle is forgotten or
/// the callback unwinds between completed edits.
///
/// # Semantic properties
///
/// Rollback restores transaction-entry state under Molecule::normalized_eq. Participant
/// sequences retain their original order. Recovery from allocation failure or an internal
/// panic during an Edit is outside this contract.
pub struct Transaction<'scope> {
    molecule: &'scope mut Molecule,
    journal: &'scope mut Vec<Undo>,
    correspondence: &'scope mut MoleculeCorrespondence,
    status: &'scope mut TransactionStatus,
}

impl Transaction<'_> {
    /// Run a callback with exclusive transaction access to the molecule.
    ///
    /// A successful commit requests acceptance. The callback must also return Ok for its
    /// changes to be retained. Returning without commit restores the molecule and returns
    /// the callback's value. Callback errors and unwinding restore unaccepted completed
    /// changes; forgetting the handle does not disable restoration.
    ///
    /// # Errors
    ///
    /// Returns the callback's error unchanged. If the callback returns Ok after the
    /// transaction was aborted, returns TransactionError::Aborted through MoleculeApplyError
    /// and E.
    pub fn run<T, E>(
        molecule: &mut Molecule,
        f: impl for<'scope> FnOnce(Transaction<'scope>) -> Result<T, E>,
    ) -> Result<T, E>
    where
        E: From<MoleculeApplyError>,
    {
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::identity(molecule.atoms().count()),
            Correspondence::identity(molecule.bonds().count()),
            Correspondence::identity(molecule.dative_bonds().count()),
            Correspondence::identity(molecule.aromatic_systems().count()),
            Correspondence::identity(molecule.multicenter_bonds().count()),
            Correspondence::identity(molecule.noncovalent_bonds().count()),
            Correspondence::identity(molecule.stereo_atoms().count()),
            Correspondence::identity(molecule.stereo_bonds().count()),
        );
        let mut guard = TransactionGuard {
            molecule,
            correspondence,
            journal: Vec::new(),
            status: TransactionStatus::Active,
        };
        let value = f(Transaction {
            molecule: guard.molecule,
            journal: &mut guard.journal,
            correspondence: &mut guard.correspondence,
            status: &mut guard.status,
        })?;
        match guard.status {
            TransactionStatus::CommitRequested => guard.status = TransactionStatus::Accepted,
            TransactionStatus::Aborted => {
                return Err(MoleculeApplyError::Transaction(TransactionError::Aborted).into());
            }
            TransactionStatus::Active
            | TransactionStatus::Accepted
            | TransactionStatus::RolledBack => {}
        }
        Ok(value)
    }

    /// Apply a batch immediately, recording undo entries for completed edits.
    ///
    /// Each call starts a separate handle namespace. Id handles name entities at this
    /// batch's entry; New handles name entities created within this batch. Intermediate
    /// state may violate integrity until a later batch repairs it.
    ///
    /// # Errors
    ///
    /// An invalid handle, old value, or edit shape restores the whole transaction and
    /// aborts it. Later apply or commit calls return TransactionError::Aborted.
    pub fn apply(&mut self, edits: Edits) -> Result<(), MoleculeApplyError> {
        if *self.status == TransactionStatus::Aborted {
            return Err(TransactionError::Aborted.into());
        }
        let mut state = ApplicationState::new(self.molecule);
        for edit in edits {
            match self.molecule.apply_edit_with_undo(edit, &mut state) {
                Ok(Some(undo)) => self.journal.push(undo),
                Ok(None) => {}
                Err(error) => {
                    while let Some(undo) = self.journal.pop() {
                        self.molecule.apply_undo(undo);
                    }
                    *self.status = TransactionStatus::Aborted;
                    return Err(error.into());
                }
            }
        }
        state.update_correspondence(self.correspondence);
        Ok(())
    }

    /// Borrow the molecule after checking its representation integrity.
    ///
    /// A failed check leaves the transaction active. After an application failure,
    /// the borrow observes the restored transaction-entry state.
    ///
    /// # Errors
    ///
    /// Returns the same integrity errors as molecule construction, without changing state.
    pub fn probe(&self) -> Result<&Molecule, MoleculeIntegrityError> {
        self.molecule.check_integrity()?;
        Ok(self.molecule)
    }

    /// Check integrity and request acceptance of the transaction's changes.
    ///
    /// Changes are retained only if the enclosing Transaction::run callback also returns Ok.
    ///
    /// # Errors
    ///
    /// Returns Aborted after an earlier application failure. Failed integrity restores
    /// transaction-entry state, aborts the transaction, and returns the integrity error.
    pub fn commit(self) -> Result<(), MoleculeApplyError> {
        if *self.status == TransactionStatus::Aborted {
            return Err(TransactionError::Aborted.into());
        }
        if let Err(error) = self.molecule.check_integrity() {
            while let Some(undo) = self.journal.pop() {
                self.molecule.apply_undo(undo);
            }
            *self.status = TransactionStatus::Aborted;
            return Err(error.into());
        }
        *self.status = TransactionStatus::CommitRequested;
        Ok(())
    }

    /// Request acceptance and return the transaction-entry to result correspondence.
    ///
    /// Entity ids are paired across all batches. Added and removed entities remain
    /// unmatched; attribute and participant changes preserve the entity's pairing.
    /// The correspondence is accumulated after each successful batch, independently of rollback.
    ///
    /// # Errors
    ///
    /// Returns the same errors and restores the same state as commit.
    ///
    /// # Semantic properties
    ///
    /// Discarding the correspondence gives commit's result and final molecule state.
    /// Regrouping the same ordered edits into batches preserves the final molecule and
    /// correspondence, provided handles name the same entities in each batch's namespace.
    pub fn tracked_commit(self) -> Result<MoleculeCorrespondence, MoleculeApplyError> {
        if *self.status == TransactionStatus::Aborted {
            return Err(TransactionError::Aborted.into());
        }
        if let Err(error) = self.molecule.check_integrity() {
            while let Some(undo) = self.journal.pop() {
                self.molecule.apply_undo(undo);
            }
            *self.status = TransactionStatus::Aborted;
            return Err(error.into());
        }
        *self.status = TransactionStatus::CommitRequested;
        Ok(mem::replace(
            self.correspondence,
            MoleculeCorrespondence::empty(),
        ))
    }

    /// Restore transaction-entry state and consume the handle.
    ///
    /// An earlier application failure remains latched as Aborted, so the enclosing
    /// run cannot accept a callback's Ok result after that failure.
    pub fn rollback(self) {
        while let Some(undo) = self.journal.pop() {
            self.molecule.apply_undo(undo);
        }
        if *self.status != TransactionStatus::Aborted {
            *self.status = TransactionStatus::RolledBack;
        }
    }
}

/// Scoped molecule transaction with an accumulated entity correspondence.
///
/// TrackedTransaction::run owns the restoration guard and lends this handle to its callback.
/// The correspondence pairs transaction-entry entity ids with their surviving ids at commit
/// across all batches. Commit requests acceptance; changes are retained only when the callback
/// also returns Ok. Forgetting the handle does not disable restoration.
///
/// # Semantic properties
///
/// For the same batches and callback outcome, discarding the commit correspondence gives the
/// same final molecule and errors as Transaction. Rollback restores transaction-entry state
/// under Molecule::normalized_eq and preserves participant order. Recovery from allocation
/// failure or an internal panic during an Edit is outside this contract.
pub struct TrackedTransaction<'scope> {
    transaction: Transaction<'scope>,
}

impl TrackedTransaction<'_> {
    /// Run a callback with exclusive tracked transaction access to the molecule.
    ///
    /// A successful commit requests acceptance and returns the correspondence to the callback.
    /// The callback must return Ok for its changes to be retained. Returning without commit,
    /// returning an error, or unwinding restores unaccepted completed changes. Returns the
    /// callback's value after acceptance or cancellation.
    ///
    /// # Errors
    ///
    /// Returns the callback's error unchanged. If the callback returns Ok after the transaction
    /// was aborted, returns TransactionError::Aborted through MoleculeApplyError and E.
    pub fn run<T, E>(
        molecule: &mut Molecule,
        f: impl for<'scope> FnOnce(TrackedTransaction<'scope>) -> Result<T, E>,
    ) -> Result<T, E>
    where
        E: From<MoleculeApplyError>,
    {
        Transaction::run(molecule, |transaction| {
            f(TrackedTransaction { transaction })
        })
    }

    /// Apply one batch immediately and update the accumulated correspondence.
    ///
    /// Each call starts a separate handle namespace. Id handles name entities at batch entry;
    /// New handles name creations within that batch. Intermediate state may violate integrity
    /// until a later batch repairs it.
    ///
    /// # Errors
    ///
    /// An invalid handle, old value or edit shape restores the whole transaction and aborts it.
    /// Later apply or commit calls return TransactionError::Aborted.
    pub fn apply(&mut self, edits: Edits) -> Result<(), MoleculeApplyError> {
        self.transaction.apply(edits)
    }

    /// Borrow the molecule after checking its representation integrity.
    ///
    /// Failed checks leave the transaction active. After an application failure, the borrow
    /// observes the restored transaction-entry state.
    ///
    /// # Errors
    ///
    /// Returns MoleculeIntegrityError without changing state when integrity does not hold.
    pub fn probe(&self) -> Result<&Molecule, MoleculeIntegrityError> {
        self.transaction.probe()
    }

    /// Check integrity, request acceptance and return the transaction-entry to commit correspondence.
    ///
    /// The correspondence covers all eight entity kinds. Additions and removals remain unmatched;
    /// attribute and participant changes preserve pairing. Participant positions are not mapped.
    /// Changes are retained only if the enclosing run callback also returns Ok.
    ///
    /// # Errors
    ///
    /// Returns TransactionError::Aborted after an earlier application failure. Failed integrity
    /// restores transaction-entry state, aborts the transaction and returns the integrity error
    /// through MoleculeApplyError.
    ///
    /// # Semantic properties
    ///
    /// The correspondence equals the composition of consecutive successful batches' mappings.
    /// Regrouping the same ordered edits into batches preserves the final molecule and
    /// correspondence, provided handles name the same entities in each batch's namespace.
    pub fn commit(self) -> Result<MoleculeCorrespondence, MoleculeApplyError> {
        self.transaction.tracked_commit()
    }

    /// Restore transaction-entry state and consume the handle.
    ///
    /// An earlier application failure remains latched as Aborted, so run rejects a subsequent
    /// callback Ok result even after explicit rollback.
    pub fn rollback(self) {
        self.transaction.rollback();
    }
}

impl Molecule {
    /// Apply prepared batches in one transaction and commit their combined result.
    ///
    /// Each batch has its own handle namespace. Integrity is checked once, after all
    /// batches have been applied. No molecule copy is made for recovery.
    ///
    /// # Errors
    ///
    /// Application and integrity failures restore the receiver to its entry state under
    /// Molecule::normalized_eq and return the original failure. Participant order is preserved.
    pub fn transact(
        &mut self,
        batches: impl IntoIterator<Item = Edits>,
    ) -> Result<(), MoleculeApplyError> {
        Transaction::run(self, |mut transaction| {
            for edits in batches {
                transaction.apply(edits)?;
            }
            transaction.commit()
        })
    }

    /// Apply prepared batches and return the receiver-entry to result correspondence.
    ///
    /// # Errors
    ///
    /// Returns the same errors and restores the same state as transact.
    ///
    /// # Semantic properties
    ///
    /// Discarding the correspondence gives transact's result and final molecule state.
    pub fn tracked_transact(
        &mut self,
        batches: impl IntoIterator<Item = Edits>,
    ) -> Result<MoleculeCorrespondence, MoleculeApplyError> {
        Transaction::run(self, |mut transaction| {
            for edits in batches {
                transaction.apply(edits)?;
            }
            transaction.tracked_commit()
        })
    }
}

struct TransactionGuard<'a> {
    molecule: &'a mut Molecule,
    journal: Vec<Undo>,
    correspondence: MoleculeCorrespondence,
    status: TransactionStatus,
}

impl Drop for TransactionGuard<'_> {
    fn drop(&mut self) {
        if matches!(
            self.status,
            TransactionStatus::Active | TransactionStatus::CommitRequested
        ) {
            while let Some(undo) = self.journal.pop() {
                self.molecule.apply_undo(undo);
            }
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum TransactionStatus {
    Active,
    CommitRequested,
    Accepted,
    RolledBack,
    Aborted,
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod lifecycle_tests {
    use std::mem::ManuallyDrop;
    use std::panic::{catch_unwind, AssertUnwindSafe};

    use rstest::{fixture, rstest};
    use umol_chem::element::Element;
    use umol_graph_core::Correspondence;

    use super::{TrackedTransaction, Transaction, TransactionError};
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::correspondence::MoleculeCorrespondence;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::edit::{
        AromaticSystemHandle, AtomFieldChange, AtomHandle, BondHandle, DativeBondHandle, Edit,
        Edits, MulticenterBondHandle, NoncovalentBondHandle, StereoAtomHandle, StereoBondHandle,
    };
    use crate::ir::entity::{Entity, EntityKind};
    use crate::ir::error::MoleculeApplyError;
    use crate::ir::id::{
        AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
        StereoAtomId, StereoBondId,
    };
    use crate::ir::ligand::{StereoLigand, StereoLigandKind};
    use crate::ir::molecule::{Molecule, MoleculeEntries, MoleculeIntegrityError};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::NoncovalentBondForm;
    use crate::ir::num::NumForm;
    use crate::ir::stereo::{StereoAtomForm, StereoBondForm, StereoCoset, StereoKind};

    #[fixture]
    fn molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..MoleculeEntries::default()
        })
    }

    #[fixture]
    fn all_entities() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::default(),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::default(),
            )],
            noncovalent: vec![([AtomId(0), AtomId(3)], NoncovalentBondForm::default())],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..MoleculeEntries::default()
        })
    }

    #[rstest]
    #[case::no_commit(false)]
    #[case::commit(true)]
    fn test_transaction_run(mut molecule: Molecule, #[case] commit: bool) {
        let original = molecule.clone();
        let result = Transaction::run(&mut molecule, |mut transaction| {
            for element in [Element::N, Element::O] {
                let mut edits = Edits::new();
                edits.add_atom(AtomForm::from_element(element));
                transaction.apply(edits)?;
            }
            let count = transaction.probe()?.atoms().count();
            if commit {
                transaction.commit()?;
            }
            Ok::<_, MoleculeApplyError>(count)
        });
        let expected = if commit {
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::N),
                    AtomForm::from_element(Element::O),
                ],
                ..MoleculeEntries::default()
            })
        } else {
            original
        };
        assert_eq!(result, Ok(3));
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_transaction_run_forgotten(mut molecule: Molecule) {
        let original = molecule.clone();
        let result = Transaction::run(&mut molecule, |mut transaction| {
            for element in [Element::N, Element::O] {
                let mut edits = Edits::new();
                edits.add_atom(AtomForm::from_element(element));
                transaction.apply(edits)?;
            }
            let _ = ManuallyDrop::new(transaction);
            Ok::<_, MoleculeApplyError>(())
        });
        assert_eq!(result, Ok(()));
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::active(None)]
    #[case::commit_requested(Some(false))]
    #[case::tracked_commit_requested(Some(true))]
    fn test_transaction_run_error(mut molecule: Molecule, #[case] commit: Option<bool>) {
        let original = molecule.clone();
        let result = Transaction::run(&mut molecule, |mut transaction| {
            for element in [Element::N, Element::O] {
                let mut edits = Edits::new();
                edits.add_atom(AtomForm::from_element(element));
                transaction.apply(edits)?;
            }
            match commit {
                Some(false) => transaction.commit()?,
                Some(true) => {
                    transaction.tracked_commit()?;
                }
                None => {}
            }
            Err::<(), _>(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch,
            ))
        });
        assert_eq!(
            result,
            Err(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch
            ))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::active(None)]
    #[case::commit_requested(Some(false))]
    #[case::tracked_commit_requested(Some(true))]
    fn test_transaction_run_unwind(mut molecule: Molecule, #[case] commit: Option<bool>) {
        let original = molecule.clone();
        let result = catch_unwind(AssertUnwindSafe(|| {
            Transaction::run(
                &mut molecule,
                |mut transaction| -> Result<(), MoleculeApplyError> {
                    for element in [Element::N, Element::O] {
                        let mut edits = Edits::new();
                        edits.add_atom(AtomForm::from_element(element));
                        transaction.apply(edits)?;
                    }
                    match commit {
                        Some(false) => transaction.commit()?,
                        Some(true) => {
                            transaction.tracked_commit()?;
                        }
                        None => {}
                    }
                    panic!("callback failure");
                },
            )
        }));
        assert_eq!(
            *result.unwrap_err().downcast::<&str>().unwrap(),
            "callback failure"
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_transaction_apply(mut molecule: Molecule) {
        let mut first = Edits::new();
        let nitrogen = first.add_atom(AtomForm::from_element(Element::N));
        assert_eq!(nitrogen, AtomHandle::New(0));
        first.add_bond(AtomHandle::Id(AtomId(0)), nitrogen, BondForm::from_order(1));
        let mut second = Edits::new();
        let oxygen = second.add_atom(AtomForm::from_element(Element::O));
        assert_eq!(oxygen, AtomHandle::New(0));
        second.add_bond(AtomHandle::Id(AtomId(1)), oxygen, BondForm::from_order(2));
        Transaction::run(&mut molecule, |mut transaction| {
            transaction.apply(first)?;
            transaction.apply(second)?;
            transaction.commit()
        })
        .unwrap();
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
            ],
            ..MoleculeEntries::default()
        });
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::drop(false)]
    #[case::rollback(true)]
    fn test_transaction_apply_error(mut molecule: Molecule, #[case] rollback: bool) {
        let original = molecule.clone();
        let result = Transaction::run(&mut molecule, |mut transaction| {
            let mut first = Edits::new();
            first.add_atom(AtomForm::from_element(Element::N));
            transaction.apply(first)?;
            let mut second = Edits::new();
            second.add_atom(AtomForm::from_element(Element::O));
            second.remove_atom(AtomHandle::Id(AtomId(2)));
            assert_eq!(
                transaction.apply(second),
                Err(MoleculeApplyError::Transaction(
                    TransactionError::HandleOutOfRange {
                        kind: EntityKind::Atom,
                        index: 2,
                        count: 2
                    }
                ))
            );
            assert_eq!(transaction.probe()?, &original);
            assert_eq!(
                transaction.apply(Edits::new()),
                Err(TransactionError::Aborted.into())
            );
            if rollback {
                transaction.rollback();
            }
            Ok::<_, MoleculeApplyError>(())
        });
        assert_eq!(result, Err(TransactionError::Aborted.into()));
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_transaction_probe(mut molecule: Molecule) {
        let original = molecule.clone();
        Transaction::run(&mut molecule, |mut transaction| {
            let mut first = Edits::new();
            first.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(0)),
                BondForm::from_order(1),
            );
            transaction.apply(first)?;
            assert_eq!(
                transaction.probe(),
                Err(MoleculeIntegrityError::DuplicateAtom {
                    entity: Entity::Bond(BondId(0)),
                    atom: AtomId(0),
                })
            );
            let mut second = Edits::new();
            second.remove_bond(BondHandle::Id(BondId(0)));
            transaction.apply(second)?;
            assert_eq!(transaction.probe()?, &original);
            transaction.commit()
        })
        .unwrap();
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::commit(false)]
    #[case::tracked_commit(true)]
    fn test_transaction_commit_error(mut molecule: Molecule, #[case] tracked: bool) {
        let original = molecule.clone();
        let constructor_error = Molecule::try_from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            bonds: vec![(AtomId(0), AtomId(0), BondForm::from_order(1))],
            ..MoleculeEntries::default()
        })
        .unwrap_err();
        let result = Transaction::run(&mut molecule, |mut transaction| {
            let mut edits = Edits::new();
            edits.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(0)),
                BondForm::from_order(1),
            );
            transaction.apply(edits)?;
            let commit = if tracked {
                transaction.tracked_commit().map(|_| ())
            } else {
                transaction.commit()
            };
            assert_eq!(
                commit,
                Err(MoleculeApplyError::Integrity(constructor_error))
            );
            Ok::<_, MoleculeApplyError>(())
        });
        assert_eq!(result, Err(TransactionError::Aborted.into()));
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::commit(false)]
    #[case::tracked_commit(true)]
    fn test_transaction_commit_aborted(mut molecule: Molecule, #[case] tracked: bool) {
        let original = molecule.clone();
        let result = Transaction::run(&mut molecule, |mut transaction| {
            let mut edits = Edits::new();
            edits.remove_atom(AtomHandle::Id(AtomId(1)));
            assert_eq!(
                transaction.apply(edits),
                Err(MoleculeApplyError::Transaction(
                    TransactionError::HandleOutOfRange {
                        kind: EntityKind::Atom,
                        index: 1,
                        count: 1
                    }
                ))
            );
            if tracked {
                transaction.tracked_commit().map(|_| ())
            } else {
                transaction.commit()
            }
        });
        assert_eq!(result, Err(TransactionError::Aborted.into()));
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_transaction_rollback(mut all_entities: Molecule) {
        let original = all_entities.clone();
        let result = Transaction::run(&mut all_entities, |mut transaction| {
            let mut first = Edits::new();
            first.remove_atom(AtomHandle::Id(AtomId(0)));
            transaction.apply(first)?;
            let mut second = Edits::new();
            second.add_atom(AtomForm::from_element(Element::N));
            transaction.apply(second)?;
            transaction.rollback();
            Ok::<_, MoleculeApplyError>(())
        });
        assert_eq!(result, Ok(()));
        assert_eq!(all_entities, original);
    }

    #[rstest]
    #[case::no_commit(false)]
    #[case::commit(true)]
    fn test_tracked_transaction_run(mut molecule: Molecule, #[case] commit: bool) {
        let original = molecule.clone();
        let result = TrackedTransaction::run(&mut molecule, |mut transaction| {
            for element in [Element::N, Element::O] {
                let mut edits = Edits::new();
                edits.add_atom(AtomForm::from_element(element));
                transaction.apply(edits)?;
            }
            let count = transaction.probe()?.atoms().count();
            if commit {
                transaction.commit()?;
            }
            Ok::<_, MoleculeApplyError>(count)
        });
        let expected = if commit {
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::N),
                    AtomForm::from_element(Element::O),
                ],
                ..MoleculeEntries::default()
            })
        } else {
            original
        };
        assert_eq!(result, Ok(3));
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_tracked_transaction_run_forgotten(mut molecule: Molecule) {
        let original = molecule.clone();
        let result = TrackedTransaction::run(&mut molecule, |mut transaction| {
            for element in [Element::N, Element::O] {
                let mut edits = Edits::new();
                edits.add_atom(AtomForm::from_element(element));
                transaction.apply(edits)?;
            }
            let _ = ManuallyDrop::new(transaction);
            Ok::<_, MoleculeApplyError>(())
        });
        assert_eq!(result, Ok(()));
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::active(false)]
    #[case::commit_requested(true)]
    fn test_tracked_transaction_run_error(mut molecule: Molecule, #[case] commit: bool) {
        let original = molecule.clone();
        let result = TrackedTransaction::run(&mut molecule, |mut transaction| {
            for element in [Element::N, Element::O] {
                let mut edits = Edits::new();
                edits.add_atom(AtomForm::from_element(element));
                transaction.apply(edits)?;
            }
            if commit {
                transaction.commit()?;
            }
            Err::<(), _>(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch,
            ))
        });
        assert_eq!(
            result,
            Err(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch
            ))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::active(false)]
    #[case::commit_requested(true)]
    fn test_tracked_transaction_run_unwind(mut molecule: Molecule, #[case] commit: bool) {
        let original = molecule.clone();
        let result = catch_unwind(AssertUnwindSafe(|| {
            TrackedTransaction::run(
                &mut molecule,
                |mut transaction| -> Result<(), MoleculeApplyError> {
                    for element in [Element::N, Element::O] {
                        let mut edits = Edits::new();
                        edits.add_atom(AtomForm::from_element(element));
                        transaction.apply(edits)?;
                    }
                    if commit {
                        transaction.commit()?;
                    }
                    panic!("callback failure");
                },
            )
        }));
        assert_eq!(
            *result.unwrap_err().downcast::<&str>().unwrap(),
            "callback failure"
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_tracked_transaction_apply(mut molecule: Molecule) {
        let mut first = Edits::new();
        let nitrogen = first.add_atom(AtomForm::from_element(Element::N));
        assert_eq!(nitrogen, AtomHandle::New(0));
        first.add_bond(AtomHandle::Id(AtomId(0)), nitrogen, BondForm::from_order(1));
        let mut second = Edits::new();
        let oxygen = second.add_atom(AtomForm::from_element(Element::O));
        assert_eq!(oxygen, AtomHandle::New(0));
        second.add_bond(AtomHandle::Id(AtomId(1)), oxygen, BondForm::from_order(2));
        TrackedTransaction::run(&mut molecule, |mut transaction| {
            transaction.apply(first)?;
            transaction.apply(second)?;
            transaction.commit()
        })
        .unwrap();
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
            ],
            ..MoleculeEntries::default()
        });
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::drop(false)]
    #[case::rollback(true)]
    fn test_tracked_transaction_apply_error(
        mut molecule: Molecule,
        #[case] rollback: bool,
        #[values(0, 1)] successful_batches: usize,
    ) {
        let original = molecule.clone();
        let result = TrackedTransaction::run(&mut molecule, |mut transaction| {
            for _ in 0..successful_batches {
                let mut first = Edits::new();
                first.add_atom(AtomForm::from_element(Element::N));
                transaction.apply(first)?;
            }
            let count = 1 + successful_batches;
            let mut second = Edits::new();
            second.add_atom(AtomForm::from_element(Element::O));
            second.remove_atom(AtomHandle::Id(AtomId(count as u32)));
            assert_eq!(
                transaction.apply(second),
                Err(MoleculeApplyError::Transaction(
                    TransactionError::HandleOutOfRange {
                        kind: EntityKind::Atom,
                        index: count,
                        count,
                    }
                ))
            );
            assert_eq!(transaction.probe()?, &original);
            assert_eq!(
                transaction.apply(Edits::new()),
                Err(TransactionError::Aborted.into())
            );
            if rollback {
                transaction.rollback();
            }
            Ok::<_, MoleculeApplyError>(())
        });
        assert_eq!(result, Err(TransactionError::Aborted.into()));
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_tracked_transaction_probe(mut molecule: Molecule) {
        let original = molecule.clone();
        TrackedTransaction::run(&mut molecule, |mut transaction| {
            let mut first = Edits::new();
            first.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(0)),
                BondForm::from_order(1),
            );
            transaction.apply(first)?;
            assert_eq!(
                transaction.probe(),
                Err(MoleculeIntegrityError::DuplicateAtom {
                    entity: Entity::Bond(BondId(0)),
                    atom: AtomId(0),
                })
            );
            let mut second = Edits::new();
            second.remove_bond(BondHandle::Id(BondId(0)));
            transaction.apply(second)?;
            assert_eq!(transaction.probe()?, &original);
            transaction.commit()
        })
        .unwrap();
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_tracked_transaction_commit() {
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: [Element::C, Element::N, Element::O, Element::F, Element::S]
                .into_iter()
                .map(AtomForm::from_element)
                .collect(),
            bonds: (0..4)
                .map(|id| (AtomId(id), AtomId(id + 1), BondForm::from_order(1)))
                .collect(),
            ..MoleculeEntries::default()
        });
        let mut first = Edits::new();
        first.remove_atom(AtomHandle::Id(AtomId(0)));
        let hydrogen = first.add_atom(AtomForm::from_element(Element::H));
        first.add_bond(AtomHandle::Id(AtomId(4)), hydrogen, BondForm::from_order(1));
        let mut second = Edits::new();
        second.remove_atom(AtomHandle::Id(AtomId(1)));
        let chlorine = second.add_atom(AtomForm::from_element(Element::Cl));
        second.add_bond(AtomHandle::Id(AtomId(2)), chlorine, BondForm::from_order(1));
        second.remove_atom(AtomHandle::Id(AtomId(4)));
        let mut ordinary = molecule.clone();
        ordinary.transact([first.clone(), second.clone()]).unwrap();
        let correspondence = TrackedTransaction::run(&mut molecule, |mut transaction| {
            transaction.apply(first)?;
            transaction.apply(second)?;
            transaction.commit()
        })
        .unwrap();
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: [Element::N, Element::F, Element::S, Element::Cl]
                .into_iter()
                .map(AtomForm::from_element)
                .collect(),
            bonds: vec![
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(1), AtomId(3), BondForm::from_order(1)),
            ],
            ..MoleculeEntries::default()
        });
        assert_eq!(molecule, expected);
        assert_eq!(molecule, ordinary);
        assert_eq!(
            correspondence,
            MoleculeCorrespondence::new(
                Correspondence::new(
                    vec![
                        (AtomId(1), AtomId(0)),
                        (AtomId(3), AtomId(1)),
                        (AtomId(4), AtomId(2))
                    ],
                    5,
                    4
                )
                .unwrap(),
                Correspondence::new(vec![(BondId(3), BondId(0))], 4, 2).unwrap(),
                Correspondence::empty(),
                Correspondence::empty(),
                Correspondence::empty(),
                Correspondence::empty(),
                Correspondence::empty(),
                Correspondence::empty(),
            )
        );
    }

    #[rstest]
    fn test_tracked_transaction_commit_additions(mut molecule: Molecule) {
        let mut edits = Edits::new();
        let nitrogen = edits.add_atom(AtomForm::from_element(Element::N));
        let oxygen = edits.add_atom(AtomForm::from_element(Element::O));
        edits.add_bond(
            AtomHandle::Id(AtomId(0)),
            nitrogen.clone(),
            BondForm::from_order(1),
        );
        edits.add_bond(nitrogen.clone(), oxygen.clone(), BondForm::from_order(1));
        edits.remove_atom(nitrogen);
        let fluorine = edits.add_atom(AtomForm::from_element(Element::F));
        edits.add_bond(oxygen, fluorine, BondForm::from_order(1));
        let correspondence = TrackedTransaction::run(&mut molecule, |mut transaction| {
            transaction.apply(edits)?;
            transaction.commit()
        })
        .unwrap();
        assert_eq!(
            molecule,
            Molecule::from_entries(MoleculeEntries {
                atoms: [Element::C, Element::O, Element::F]
                    .into_iter()
                    .map(AtomForm::from_element)
                    .collect(),
                bonds: vec![(AtomId(1), AtomId(2), BondForm::from_order(1))],
                ..MoleculeEntries::default()
            })
        );
        assert_eq!(
            correspondence,
            MoleculeCorrespondence::new(
                Correspondence::new(vec![(AtomId(0), AtomId(0))], 1, 3).unwrap(),
                Correspondence::new(Vec::new(), 0, 1).unwrap(),
                Correspondence::empty(),
                Correspondence::empty(),
                Correspondence::empty(),
                Correspondence::empty(),
                Correspondence::empty(),
                Correspondence::empty(),
            )
        );
    }

    #[rstest]
    #[case::dative_bond(EntityKind::DativeBond)]
    #[case::aromatic_system(EntityKind::AromaticSystem)]
    #[case::multicenter_bond(EntityKind::MulticenterBond)]
    #[case::noncovalent_bond(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_tracked_transaction_commit_overlays(
        mut all_entities: Molecule,
        #[case] kind: EntityKind,
    ) {
        let original = all_entities.clone();
        let mut edits = Edits::new();
        match kind {
            EntityKind::DativeBond => {
                let donors = vec![AtomHandle::Id(AtomId(0))];
                let acceptor = AtomHandle::Id(AtomId(1));
                let attributes = DativeBondForm::from_order(1);
                edits.remove_dative_bonds(vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    donors.clone(),
                    acceptor.clone(),
                    attributes.clone(),
                )]);
                edits.add_dative_bond(donors, acceptor, attributes);
            }
            EntityKind::AromaticSystem => {
                let atoms = vec![
                    AtomHandle::Id(AtomId(0)),
                    AtomHandle::Id(AtomId(1)),
                    AtomHandle::Id(AtomId(2)),
                ];
                let attributes = AromaticSystemForm::default();
                edits.remove_aromatic_systems(vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    atoms.clone(),
                    attributes.clone(),
                )]);
                edits.add_aromatic_system(atoms, attributes);
            }
            EntityKind::MulticenterBond => {
                let atoms = vec![
                    AtomHandle::Id(AtomId(0)),
                    AtomHandle::Id(AtomId(1)),
                    AtomHandle::Id(AtomId(2)),
                ];
                let attributes = MulticenterBondForm::default();
                edits.remove_multicenter_bonds(vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    atoms.clone(),
                    attributes.clone(),
                )]);
                edits.add_multicenter_bond(atoms, attributes);
            }
            EntityKind::NoncovalentBond => {
                let atoms = [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(3))];
                let attributes = NoncovalentBondForm::default();
                edits.remove_noncovalent_bonds(vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    atoms.clone(),
                    attributes.clone(),
                )]);
                edits.add_noncovalent_bond(atoms, attributes);
            }
            EntityKind::StereoAtom => {
                let site = AtomHandle::Id(AtomId(0));
                let ligands: Vec<_> = original
                    .stereo_atom(StereoAtomId(0))
                    .ligand_ids()
                    .iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect();
                let attributes = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1));
                edits.remove_stereo_atoms(vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    site.clone(),
                    ligands.clone(),
                    attributes.clone(),
                )]);
                edits.add_stereo_atom(site, ligands, attributes);
            }
            EntityKind::StereoBond => {
                let site = BondHandle::Id(BondId(0));
                let ligands: Vec<_> = original
                    .stereo_bond(StereoBondId(0))
                    .ligand_ids()
                    .iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect();
                let attributes = StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1));
                edits.remove_stereo_bonds(vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    site.clone(),
                    ligands.clone(),
                    attributes.clone(),
                )]);
                edits.add_stereo_bond(site, ligands, attributes);
            }
            EntityKind::Atom | EntityKind::Bond => unreachable!(),
        }
        let correspondence = TrackedTransaction::run(&mut all_entities, |mut transaction| {
            transaction.apply(edits)?;
            transaction.commit()
        })
        .unwrap();
        assert_eq!(all_entities, original);
        assert_eq!(
            correspondence,
            MoleculeCorrespondence::new(
                Correspondence::identity(4),
                Correspondence::identity(4),
                if kind == EntityKind::DativeBond {
                    Correspondence::new(Vec::new(), 1, 1).unwrap()
                } else {
                    Correspondence::identity(1)
                },
                if kind == EntityKind::AromaticSystem {
                    Correspondence::new(Vec::new(), 1, 1).unwrap()
                } else {
                    Correspondence::identity(1)
                },
                if kind == EntityKind::MulticenterBond {
                    Correspondence::new(Vec::new(), 1, 1).unwrap()
                } else {
                    Correspondence::identity(1)
                },
                if kind == EntityKind::NoncovalentBond {
                    Correspondence::new(Vec::new(), 1, 1).unwrap()
                } else {
                    Correspondence::identity(1)
                },
                if kind == EntityKind::StereoAtom {
                    Correspondence::new(Vec::new(), 1, 1).unwrap()
                } else {
                    Correspondence::identity(1)
                },
                if kind == EntityKind::StereoBond {
                    Correspondence::new(Vec::new(), 1, 1).unwrap()
                } else {
                    Correspondence::identity(1)
                },
            )
        );
    }

    #[rstest]
    fn test_tracked_transaction_commit_modification(mut all_entities: Molecule) {
        let edits = Edits::from_iter([
            Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            },
            Edit::ReplaceNoncovalentBondAtoms {
                id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                old: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(3))],
                new: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(3))],
            },
        ]);
        let expected = all_entities.clone().apply(edits.clone()).unwrap();
        let correspondence = TrackedTransaction::run(&mut all_entities, |mut transaction| {
            transaction.apply(edits)?;
            transaction.commit()
        })
        .unwrap();
        assert_eq!(all_entities, expected);
        assert_eq!(
            correspondence,
            MoleculeCorrespondence::new(
                Correspondence::identity(4),
                Correspondence::identity(4),
                Correspondence::identity(1),
                Correspondence::identity(1),
                Correspondence::identity(1),
                Correspondence::identity(1),
                Correspondence::identity(1),
                Correspondence::identity(1),
            )
        );
    }

    #[rstest]
    fn test_tracked_transaction_commit_cascade(mut all_entities: Molecule) {
        let mut edits = Edits::new();
        edits.remove_atom(AtomHandle::Id(AtomId(0)));
        let correspondence = TrackedTransaction::run(&mut all_entities, |mut transaction| {
            transaction.apply(edits)?;
            transaction.commit()
        })
        .unwrap();
        assert_eq!(
            all_entities,
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 3],
                bonds: vec![(AtomId(1), AtomId(2), BondForm::from_order(1))],
                ..MoleculeEntries::default()
            })
        );
        assert_eq!(
            correspondence,
            MoleculeCorrespondence::new(
                Correspondence::new(
                    vec![
                        (AtomId(1), AtomId(0)),
                        (AtomId(2), AtomId(1)),
                        (AtomId(3), AtomId(2))
                    ],
                    4,
                    3
                )
                .unwrap(),
                Correspondence::new(vec![(BondId(1), BondId(0))], 4, 1).unwrap(),
                Correspondence::new(Vec::new(), 1, 0).unwrap(),
                Correspondence::new(Vec::new(), 1, 0).unwrap(),
                Correspondence::new(Vec::new(), 1, 0).unwrap(),
                Correspondence::new(Vec::new(), 1, 0).unwrap(),
                Correspondence::new(Vec::new(), 1, 0).unwrap(),
                Correspondence::new(Vec::new(), 1, 0).unwrap(),
            )
        );
    }

    #[rstest]
    #[case::no_batches(0)]
    #[case::empty_batches(3)]
    fn test_tracked_transaction_commit_identity(
        mut all_entities: Molecule,
        #[case] batches: usize,
    ) {
        let original = all_entities.clone();
        let correspondence = TrackedTransaction::run(&mut all_entities, |mut transaction| {
            for _ in 0..batches {
                transaction.apply(Edits::new())?;
            }
            transaction.commit()
        })
        .unwrap();
        assert_eq!(all_entities, original);
        assert_eq!(
            correspondence,
            MoleculeCorrespondence::new(
                Correspondence::identity(4),
                Correspondence::identity(4),
                Correspondence::identity(1),
                Correspondence::identity(1),
                Correspondence::identity(1),
                Correspondence::identity(1),
                Correspondence::identity(1),
                Correspondence::identity(1),
            )
        );
    }

    #[rstest]
    fn test_tracked_transaction_commit_error(mut molecule: Molecule) {
        let original = molecule.clone();
        let result = TrackedTransaction::run(&mut molecule, |mut transaction| {
            let mut edits = Edits::new();
            edits.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(0)),
                BondForm::from_order(1),
            );
            transaction.apply(edits)?;
            assert_eq!(
                transaction.commit(),
                Err(MoleculeApplyError::Integrity(
                    MoleculeIntegrityError::DuplicateAtom {
                        entity: Entity::Bond(BondId(0)),
                        atom: AtomId(0),
                    }
                ))
            );
            Ok::<_, MoleculeApplyError>(())
        });
        assert_eq!(result, Err(TransactionError::Aborted.into()));
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_tracked_transaction_commit_aborted(mut molecule: Molecule) {
        let original = molecule.clone();
        let result = TrackedTransaction::run(&mut molecule, |mut transaction| {
            let mut edits = Edits::new();
            edits.remove_atom(AtomHandle::Id(AtomId(1)));
            assert_eq!(
                transaction.apply(edits),
                Err(MoleculeApplyError::Transaction(
                    TransactionError::HandleOutOfRange {
                        kind: EntityKind::Atom,
                        index: 1,
                        count: 1
                    }
                ))
            );
            transaction.commit()
        });
        assert_eq!(result, Err(TransactionError::Aborted.into()));
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_tracked_transaction_rollback(mut all_entities: Molecule) {
        let original = all_entities.clone();
        let result = TrackedTransaction::run(&mut all_entities, |mut transaction| {
            let mut first = Edits::new();
            first.remove_atom(AtomHandle::Id(AtomId(0)));
            transaction.apply(first)?;
            let mut second = Edits::new();
            second.add_atom(AtomForm::from_element(Element::N));
            transaction.apply(second)?;
            transaction.rollback();
            Ok::<_, MoleculeApplyError>(())
        });
        assert_eq!(result, Ok(()));
        assert_eq!(all_entities, original);
    }

    #[rstest]
    #[case::application(
        Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::Lit(0) },
        }]),
        MoleculeApplyError::Transaction(TransactionError::OldStateMismatch),
    )]
    #[case::integrity(
        Edits::from_iter([Edit::ReplaceNoncovalentBondAtoms {
            id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            old: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(3))],
            new: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(0))],
        }]),
        MoleculeApplyError::Integrity(MoleculeIntegrityError::DuplicateAtom {
            entity: Entity::NoncovalentBond(NoncovalentBondId(0)), atom: AtomId(0),
        }),
    )]
    fn test_molecule_transact_error(
        mut all_entities: Molecule,
        #[case] failing: Edits,
        #[case] error: MoleculeApplyError,
        #[values(false, true)] tracked: bool,
    ) {
        let original = all_entities.clone();
        let first = Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(1)),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        }]);
        let result = if tracked {
            all_entities.tracked_transact([first, failing]).map(|_| ())
        } else {
            all_entities.transact([first, failing])
        };
        assert_eq!(result, Err(error));
        assert_eq!(all_entities, original);
    }
}
