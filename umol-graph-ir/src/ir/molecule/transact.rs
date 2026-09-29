//! Scoped molecule transactions with private undo-based restoration.

use thiserror::Error;

use super::Molecule;
use crate::ir::edit::Undo;
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

    /// An earlier application or commit failure aborted the transaction.
    #[error("transaction has been aborted")]
    Aborted,

    #[error("rollback failed after apply error: apply={apply}; rollback={rollback}")]
    RollbackFailed {
        apply: Box<TransactionError>,
        rollback: Box<TransactionError>,
    },

    /// The rollback journal cannot be structurally applied to the supplied editor state.
    #[error("rollback journal does not match editor state")]
    RollbackStateMismatch,
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
        let mut guard = TransactionGuard {
            molecule,
            journal: Vec::new(),
            status: TransactionStatus::Active,
        };
        let value = f(Transaction {
            molecule: guard.molecule,
            journal: &mut guard.journal,
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
}

struct TransactionGuard<'a> {
    molecule: &'a mut Molecule,
    journal: Vec<Undo>,
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
mod tests {
    use std::mem::ManuallyDrop;
    use std::panic::{catch_unwind, AssertUnwindSafe};

    use rstest::{fixture, rstest};
    use umol_chem::element::Element;

    use super::{Transaction, TransactionError, TransactionStatus};
    use crate::ir::atom::AtomForm;
    use crate::ir::edit::Edit;
    use crate::ir::error::MoleculeApplyError;
    use crate::ir::molecule::apply::ApplicationState;
    use crate::ir::molecule::{Molecule, MoleculeEntries};

    #[fixture]
    fn molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..MoleculeEntries::default()
        })
    }

    #[rstest]
    #[case::no_commit(false)]
    #[case::commit_requested(true)]
    fn test_transaction_run(
        mut molecule: Molecule,
        #[case] commit_requested: bool,
        #[values(false, true)] forgotten: bool,
    ) {
        let expected = if commit_requested {
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::N),
                    AtomForm::from_element(Element::O),
                ],
                ..MoleculeEntries::default()
            })
        } else {
            molecule.clone()
        };
        let result = Transaction::run(&mut molecule, |transaction| {
            let mut state = ApplicationState::new(transaction.molecule);
            for element in [Element::N, Element::O] {
                let undo = transaction
                    .molecule
                    .apply_edit_with_undo(
                        Edit::AddAtoms {
                            atoms: vec![AtomForm::from_element(element)],
                        },
                        &mut state,
                    )
                    .unwrap()
                    .unwrap();
                transaction.journal.push(undo);
            }
            if commit_requested {
                *transaction.status = TransactionStatus::CommitRequested;
            }
            let count = transaction.molecule.atoms().count();
            if forgotten {
                let _ = ManuallyDrop::new(transaction);
            }
            Ok::<_, MoleculeApplyError>(count)
        });

        assert_eq!(result, Ok(3));
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::active(false)]
    #[case::commit_requested(true)]
    fn test_transaction_run_error(mut molecule: Molecule, #[case] commit_requested: bool) {
        let original = molecule.clone();
        let result = Transaction::run(&mut molecule, |transaction| {
            let mut state = ApplicationState::new(transaction.molecule);
            for element in [Element::N, Element::O] {
                let undo = transaction
                    .molecule
                    .apply_edit_with_undo(
                        Edit::AddAtoms {
                            atoms: vec![AtomForm::from_element(element)],
                        },
                        &mut state,
                    )
                    .unwrap()
                    .unwrap();
                transaction.journal.push(undo);
            }
            if commit_requested {
                *transaction.status = TransactionStatus::CommitRequested;
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
    fn test_transaction_run_unwind(mut molecule: Molecule, #[case] commit_requested: bool) {
        let original = molecule.clone();
        let result = catch_unwind(AssertUnwindSafe(|| {
            Transaction::run(
                &mut molecule,
                |transaction| -> Result<(), MoleculeApplyError> {
                    let mut state = ApplicationState::new(transaction.molecule);
                    for element in [Element::N, Element::O] {
                        let undo = transaction
                            .molecule
                            .apply_edit_with_undo(
                                Edit::AddAtoms {
                                    atoms: vec![AtomForm::from_element(element)],
                                },
                                &mut state,
                            )
                            .unwrap()
                            .unwrap();
                        transaction.journal.push(undo);
                    }
                    if commit_requested {
                        *transaction.status = TransactionStatus::CommitRequested;
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
    #[case::rolled_back(TransactionStatus::RolledBack, Ok(()))]
    #[case::aborted(
        TransactionStatus::Aborted,
        Err(MoleculeApplyError::Transaction(TransactionError::Aborted))
    )]
    fn test_transaction_run_status(
        mut molecule: Molecule,
        #[case] status: TransactionStatus,
        #[case] expected: Result<(), MoleculeApplyError>,
    ) {
        let original = molecule.clone();
        let result = Transaction::run(&mut molecule, |transaction| {
            *transaction.status = status;
            Ok::<_, MoleculeApplyError>(())
        });

        assert_eq!(result, expected);
        assert_eq!(molecule, original);
    }
}
