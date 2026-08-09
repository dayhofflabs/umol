//! Python ownership wrappers for transactional molecule editing.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use umol_graph_ir::ir::{
    MoleculeEditor as GraphIrMoleculeEditor, Transaction as GraphIrTransaction,
};

use crate::edit::Edits;
use crate::error::transaction_error;
use crate::molecule::MoleculeAst;

/// A mutable molecule editor that can be inspected before it is finalized.
#[pyclass]
pub struct MoleculeEditor {
    inner: Option<GraphIrMoleculeEditor>,
}

impl MoleculeEditor {
    pub(crate) fn from_rust(editor: GraphIrMoleculeEditor) -> Self {
        Self {
            inner: Some(editor),
        }
    }
}

#[pymethods]
impl MoleculeEditor {
    /// Materialize the editor's current state without consuming it.
    fn snapshot(&self) -> PyResult<MoleculeAst> {
        self.inner
            .as_ref()
            .map(|editor| MoleculeAst::from_rust(editor.snapshot()))
            .ok_or_else(consumed_editor_error)
    }

    /// Finalize the editor and consume its mutable state.
    fn build(&mut self) -> PyResult<MoleculeAst> {
        self.inner
            .take()
            .map(|editor| MoleculeAst::from_rust(editor.build()))
            .ok_or_else(consumed_editor_error)
    }

    /// Apply a checked edit batch atomically and return its rollback journal.
    fn transact(&mut self, py: Python<'_>, edits: Py<Edits>) -> PyResult<Transaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        editor
            .transact(edits.bind(py).borrow().to_rust())
            .map(|transaction| Transaction {
                inner: Some(transaction),
            })
            .map_err(transaction_error)
    }
}

/// A detached, one-shot rollback journal for a successful transaction.
#[pyclass]
#[derive(Debug)]
pub struct Transaction {
    inner: Option<GraphIrTransaction>,
}

#[pymethods]
impl Transaction {
    /// Roll back against the editor state produced by this transaction.
    fn rollback(&mut self, py: Python<'_>, editor: Py<MoleculeEditor>) -> PyResult<()> {
        if self.inner.is_none() {
            return Err(consumed_transaction_error());
        }

        let mut editor = editor.bind(py).try_borrow_mut()?;
        let editor = editor.inner.as_mut().ok_or_else(consumed_editor_error)?;
        let transaction = self.inner.take().ok_or_else(consumed_transaction_error)?;
        transaction.rollback(editor).map_err(transaction_error)
    }
}

fn consumed_editor_error() -> PyErr {
    PyRuntimeError::new_err("molecule editor has been consumed")
}

fn consumed_transaction_error() -> PyErr {
    PyRuntimeError::new_err("transaction has been consumed")
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::PyRuntimeError;
    use pyo3::prelude::*;
    use rstest::{fixture, rstest};
    use umol_chem::element::Element as ChemElement;
    use umol_graph_ir::ir::{
        AtomFieldChange as GraphIrAtomFieldChange, AtomForm as GraphIrAtomForm,
        AtomHandle as GraphIrAtomHandle, AtomId as GraphIrAtomId, Edit as GraphIrEdit,
        Edits as GraphIrEdits, NumForm as GraphIrNumForm,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::error::TransactionError;

    #[fixture]
    fn carbon_editor() -> MoleculeEditor {
        MoleculeEditor {
            inner: Some(mol_dsl!(r#"{:atoms ["C"]}"#).edit()),
        }
    }

    #[fixture]
    fn add_nitrogen() -> Edits {
        let mut edits = GraphIrEdits::new();
        edits.add_atom(GraphIrAtomForm::from_element(ChemElement::N));
        Edits::from_rust(edits)
    }

    #[rstest]
    fn test_molecule_editor_snapshot() {
        let initial = mol_dsl!(r#"{:atoms ["C"]}"#);
        let mut editor = MoleculeEditor {
            inner: Some(initial.edit()),
        };

        let first = editor.snapshot().unwrap();
        editor
            .inner
            .as_mut()
            .unwrap()
            .add_atom(GraphIrAtomForm::from_element(ChemElement::N));
        let second = editor.snapshot().unwrap();

        assert_eq!(first.inner(), &initial);
        assert_eq!(second.inner(), &mol_dsl!(r#"{:atoms ["C" "N"]}"#));
        assert_eq!(editor.snapshot().unwrap(), second);
    }

    #[rstest]
    fn test_molecule_editor_build() {
        let initial = mol_dsl!(r#"{:atoms ["C"]}"#);
        let mut editor = MoleculeEditor {
            inner: Some(initial.edit()),
        };
        let snapshot = editor.snapshot().unwrap();

        let mut built = editor.build().unwrap();
        *built.inner_mut().atom_mut(GraphIrAtomId(0)).ast =
            GraphIrAtomForm::from_element(ChemElement::N);
        let snapshot_error = editor.snapshot().unwrap_err();
        let build_error = editor.build().unwrap_err();

        assert_eq!(snapshot.inner(), &initial);
        assert_eq!(built.inner(), &mol_dsl!(r#"{:atoms ["N"]}"#));
        Python::attach(|py| {
            assert!(snapshot_error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                snapshot_error
                    .value(py)
                    .str()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "molecule editor has been consumed"
            );
            assert!(build_error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                build_error
                    .value(py)
                    .str()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "molecule editor has been consumed"
            );
        });
    }

    #[rstest]
    fn test_molecule_editor_transact(mut carbon_editor: MoleculeEditor, add_nitrogen: Edits) {
        Python::attach(|py| {
            let edits = Py::new(py, add_nitrogen).unwrap();

            let _transaction = carbon_editor.transact(py, edits).unwrap();

            assert_eq!(
                carbon_editor.snapshot().unwrap().inner(),
                &mol_dsl!(r#"{:atoms ["C" "N"]}"#)
            );
        });
    }

    #[rstest]
    fn test_molecule_editor_transact_error(mut carbon_editor: MoleculeEditor) {
        let initial = carbon_editor.snapshot().unwrap();
        let mut edits = GraphIrEdits::new();
        edits.add_atom(GraphIrAtomForm::from_element(ChemElement::N));
        edits.push(GraphIrEdit::ModifyAtomField {
            id: GraphIrAtomHandle::Id(GraphIrAtomId(7)),
            change: GraphIrAtomFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(1),
            },
        });

        Python::attach(|py| {
            let edits = Py::new(py, Edits::from_rust(edits)).unwrap();

            let error = carbon_editor.transact(py, edits).unwrap_err();

            assert!(error.is_instance_of::<TransactionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "atom handle 7 is out of range for 1 entries"
            );
            assert_eq!(carbon_editor.snapshot().unwrap(), initial);
        });
    }

    #[rstest]
    fn test_transaction_rollback(carbon_editor: MoleculeEditor, add_nitrogen: Edits) {
        Python::attach(|py| {
            let editor = Py::new(py, carbon_editor).unwrap();
            let edits = Py::new(py, add_nitrogen).unwrap();
            let mut transaction = editor.bind(py).borrow_mut().transact(py, edits).unwrap();
            assert_eq!(
                editor.bind(py).borrow().snapshot().unwrap().inner(),
                &mol_dsl!(r#"{:atoms ["C" "N"]}"#)
            );

            transaction.rollback(py, editor.clone_ref(py)).unwrap();
            let second_error = transaction.rollback(py, editor.clone_ref(py)).unwrap_err();

            assert_eq!(
                editor.bind(py).borrow().snapshot().unwrap().inner(),
                &mol_dsl!(r#"{:atoms ["C"]}"#)
            );
            assert!(second_error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                second_error
                    .value(py)
                    .str()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "transaction has been consumed"
            );
        });
    }

    #[rstest]
    fn test_transaction_rollback_error(carbon_editor: MoleculeEditor, add_nitrogen: Edits) {
        Python::attach(|py| {
            let source = Py::new(py, carbon_editor).unwrap();
            let edits = Py::new(py, add_nitrogen).unwrap();
            let mut transaction = source.bind(py).borrow_mut().transact(py, edits).unwrap();
            let incompatible = Py::new(
                py,
                MoleculeEditor {
                    inner: Some(mol_dsl!(r#"{:atoms []}"#).edit()),
                },
            )
            .unwrap();

            let error = transaction
                .rollback(py, incompatible.clone_ref(py))
                .unwrap_err();
            let second_error = transaction.rollback(py, incompatible).unwrap_err();

            assert!(error.is_instance_of::<TransactionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "rollback journal does not match editor state"
            );
            assert!(second_error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                second_error
                    .value(py)
                    .str()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "transaction has been consumed"
            );
        });
    }

    #[rstest]
    fn test_molecule_editor_transact_consumed(
        mut carbon_editor: MoleculeEditor,
        add_nitrogen: Edits,
    ) {
        carbon_editor.build().unwrap();

        Python::attach(|py| {
            let edits = Py::new(py, add_nitrogen).unwrap();

            let error = carbon_editor.transact(py, edits).unwrap_err();

            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "molecule editor has been consumed"
            );
        });
    }

    #[rstest]
    fn test_transaction_rollback_consumed_editor(
        carbon_editor: MoleculeEditor,
        add_nitrogen: Edits,
    ) {
        Python::attach(|py| {
            let editor = Py::new(py, carbon_editor).unwrap();
            let edits = Py::new(py, add_nitrogen).unwrap();
            let mut transaction = editor.bind(py).borrow_mut().transact(py, edits).unwrap();
            editor.bind(py).borrow_mut().build().unwrap();

            let error = transaction.rollback(py, editor).unwrap_err();

            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "molecule editor has been consumed"
            );
        });
    }
}
