//! Python ownership wrappers for transactional molecule editing.

use pyo3::exceptions::{PyIndexError, PyRuntimeError};
use pyo3::prelude::*;
use umol_graph_ir::ir::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MoleculeEditor as GraphIrMoleculeEditor,
    MulticenterBondId, NoncovalentBondId, StereoAtomId, StereoBondId,
    Transaction as GraphIrTransaction,
};

use crate::compact::MoleculeCompaction;
use crate::correspondence::MoleculeCorrespondence;
use crate::edit::Edits;
use crate::error::{molecule_integrity_error, transaction_error};
use crate::molecule::Molecule;

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
    fn snapshot(&self) -> PyResult<Molecule> {
        self.inner
            .as_ref()
            .ok_or_else(consumed_editor_error)?
            .snapshot()
            .map(Molecule::from_rust)
            .map_err(molecule_integrity_error)
    }

    /// Materialize the current state and its initial-to-current correspondence.
    fn tracked_snapshot(&self) -> PyResult<(Molecule, MoleculeCorrespondence)> {
        self.inner
            .as_ref()
            .ok_or_else(consumed_editor_error)?
            .tracked_snapshot()
            .map(|(molecule, correspondence)| {
                (
                    Molecule::from_rust(molecule),
                    MoleculeCorrespondence::from_rust(correspondence),
                )
            })
            .map_err(molecule_integrity_error)
    }

    /// Finalize the editor and consume its mutable state.
    fn build(&mut self) -> PyResult<Molecule> {
        self.inner
            .take()
            .ok_or_else(consumed_editor_error)?
            .try_build()
            .map(Molecule::from_rust)
            .map_err(molecule_integrity_error)
    }

    /// Finalize the editor and return its initial-to-result correspondence.
    fn tracked_build(&mut self) -> PyResult<(Molecule, MoleculeCorrespondence)> {
        self.inner
            .take()
            .ok_or_else(consumed_editor_error)?
            .try_tracked_build()
            .map(|(molecule, correspondence)| {
                (
                    Molecule::from_rust(molecule),
                    MoleculeCorrespondence::from_rust(correspondence),
                )
            })
            .map_err(molecule_integrity_error)
    }

    /// Consume this editor and apply a checked edit batch without constructing a rollback journal.
    fn apply(&mut self, py: Python<'_>, edits: Py<Edits>) -> PyResult<Self> {
        self.inner
            .take()
            .ok_or_else(consumed_editor_error)?
            .apply(edits.bind(py).borrow().to_rust().clone())
            .map(Self::from_rust)
            .map_err(transaction_error)
    }

    /// Apply the same consuming batch and return its input-to-result correspondence.
    fn tracked_apply(
        &mut self,
        py: Python<'_>,
        edits: Py<Edits>,
    ) -> PyResult<(Self, MoleculeCorrespondence)> {
        self.inner
            .take()
            .ok_or_else(consumed_editor_error)?
            .tracked_apply(edits.bind(py).borrow().to_rust().clone())
            .map(|(editor, correspondence)| {
                (
                    Self::from_rust(editor),
                    MoleculeCorrespondence::from_rust(correspondence),
                )
            })
            .map_err(transaction_error)
    }

    /// Apply a checked edit batch atomically and return its rollback journal.
    fn transact(&mut self, py: Python<'_>, edits: Py<Edits>) -> PyResult<Transaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        editor
            .transact(edits.bind(py).borrow().to_rust().clone())
            .map(|transaction| Transaction {
                inner: Some(transaction),
            })
            .map_err(transaction_error)
    }

    /// Apply the same atomic batch and return its input-to-result correspondence.
    fn tracked_transact(
        &mut self,
        py: Python<'_>,
        edits: Py<Edits>,
    ) -> PyResult<(Transaction, MoleculeCorrespondence)> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        editor
            .tracked_transact(edits.bind(py).borrow().to_rust().clone())
            .map(|(transaction, correspondence)| {
                (
                    Transaction {
                        inner: Some(transaction),
                    },
                    MoleculeCorrespondence::from_rust(correspondence),
                )
            })
            .map_err(transaction_error)
    }

    /// Remove atoms and bonds, cascading dependent entities.
    fn remove(&mut self, atoms: Vec<u32>, bonds: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&atoms, editor.atom_count(), "atom")?;
        ensure_in_range(&bonds, editor.bond_count(), "bond")?;
        let atoms = atoms.into_iter().map(AtomId).collect::<Vec<_>>();
        let bonds = bonds.into_iter().map(BondId).collect::<Vec<_>>();
        editor.remove(&atoms, &bonds);
        Ok(())
    }

    /// Remove atoms and bonds, returning the source-to-result compaction.
    fn tracked_remove(&mut self, atoms: Vec<u32>, bonds: Vec<u32>) -> PyResult<MoleculeCompaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&atoms, editor.atom_count(), "atom")?;
        ensure_in_range(&bonds, editor.bond_count(), "bond")?;
        let atoms = atoms.into_iter().map(AtomId).collect::<Vec<_>>();
        let bonds = bonds.into_iter().map(BondId).collect::<Vec<_>>();
        Ok(MoleculeCompaction::from_rust(
            editor.tracked_remove(&atoms, &bonds),
        ))
    }

    /// Remove dative bonds and compact that entity space.
    fn remove_dative_bonds(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.dative_bond_count(), "dative bond")?;
        editor.remove_dative_bonds(&ids.into_iter().map(DativeBondId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove dative bonds and return the source-to-result compaction.
    fn tracked_remove_dative_bonds(&mut self, ids: Vec<u32>) -> PyResult<MoleculeCompaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.dative_bond_count(), "dative bond")?;
        Ok(MoleculeCompaction::from_rust(
            editor.tracked_remove_dative_bonds(
                &ids.into_iter().map(DativeBondId).collect::<Vec<_>>(),
            ),
        ))
    }

    /// Remove aromatic systems and compact that entity space.
    fn remove_aromatic_systems(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.aromatic_system_count(), "aromatic system")?;
        editor.remove_aromatic_systems(&ids.into_iter().map(AromaticSystemId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove aromatic systems and return the source-to-result compaction.
    fn tracked_remove_aromatic_systems(&mut self, ids: Vec<u32>) -> PyResult<MoleculeCompaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.aromatic_system_count(), "aromatic system")?;
        Ok(MoleculeCompaction::from_rust(
            editor.tracked_remove_aromatic_systems(
                &ids.into_iter().map(AromaticSystemId).collect::<Vec<_>>(),
            ),
        ))
    }

    /// Remove multicenter bonds and compact that entity space.
    fn remove_multicenter_bonds(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.multicenter_bond_count(), "multicenter bond")?;
        editor
            .remove_multicenter_bonds(&ids.into_iter().map(MulticenterBondId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove multicenter bonds and return the source-to-result compaction.
    fn tracked_remove_multicenter_bonds(&mut self, ids: Vec<u32>) -> PyResult<MoleculeCompaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.multicenter_bond_count(), "multicenter bond")?;
        Ok(MoleculeCompaction::from_rust(
            editor.tracked_remove_multicenter_bonds(
                &ids.into_iter().map(MulticenterBondId).collect::<Vec<_>>(),
            ),
        ))
    }

    /// Remove noncovalent bonds and compact that entity space.
    fn remove_noncovalent_bonds(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.noncovalent_bond_count(), "noncovalent bond")?;
        editor
            .remove_noncovalent_bonds(&ids.into_iter().map(NoncovalentBondId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove noncovalent bonds and return the source-to-result compaction.
    fn tracked_remove_noncovalent_bonds(&mut self, ids: Vec<u32>) -> PyResult<MoleculeCompaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.noncovalent_bond_count(), "noncovalent bond")?;
        Ok(MoleculeCompaction::from_rust(
            editor.tracked_remove_noncovalent_bonds(
                &ids.into_iter().map(NoncovalentBondId).collect::<Vec<_>>(),
            ),
        ))
    }

    /// Remove stereo atoms and compact that entity space.
    fn remove_stereo_atoms(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.stereo_atom_count(), "stereo atom")?;
        editor.remove_stereo_atoms(&ids.into_iter().map(StereoAtomId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove stereo atoms and return the source-to-result compaction.
    fn tracked_remove_stereo_atoms(&mut self, ids: Vec<u32>) -> PyResult<MoleculeCompaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.stereo_atom_count(), "stereo atom")?;
        Ok(MoleculeCompaction::from_rust(
            editor.tracked_remove_stereo_atoms(
                &ids.into_iter().map(StereoAtomId).collect::<Vec<_>>(),
            ),
        ))
    }

    /// Remove stereo bonds and compact that entity space.
    fn remove_stereo_bonds(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.stereo_bond_count(), "stereo bond")?;
        editor.remove_stereo_bonds(&ids.into_iter().map(StereoBondId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove stereo bonds and return the source-to-result compaction.
    fn tracked_remove_stereo_bonds(&mut self, ids: Vec<u32>) -> PyResult<MoleculeCompaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.stereo_bond_count(), "stereo bond")?;
        Ok(MoleculeCompaction::from_rust(
            editor.tracked_remove_stereo_bonds(
                &ids.into_iter().map(StereoBondId).collect::<Vec<_>>(),
            ),
        ))
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

    /// Roll back and return the rollback-input-to-restored-state correspondence.
    fn tracked_rollback(
        &mut self,
        py: Python<'_>,
        editor: Py<MoleculeEditor>,
    ) -> PyResult<MoleculeCorrespondence> {
        if self.inner.is_none() {
            return Err(consumed_transaction_error());
        }

        let mut editor = editor.bind(py).try_borrow_mut()?;
        let editor = editor.inner.as_mut().ok_or_else(consumed_editor_error)?;
        let transaction = self.inner.take().ok_or_else(consumed_transaction_error)?;
        transaction
            .tracked_rollback(editor)
            .map(MoleculeCorrespondence::from_rust)
            .map_err(transaction_error)
    }
}

fn ensure_in_range(ids: &[u32], count: usize, entity: &str) -> PyResult<()> {
    if ids.iter().any(|&id| id as usize >= count) {
        return Err(PyIndexError::new_err(format!("{entity} id out of range")));
    }
    Ok(())
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
    use rstest::{fixture, rstest};
    use umol_chem::element::Element as ChemElement;
    use umol_graph_ir::ir::{
        AtomFieldChange as GraphIrAtomFieldChange, AtomForm as GraphIrAtomForm,
        AtomHandle as GraphIrAtomHandle, AtomId as GraphIrAtomId, BondForm as GraphIrBondForm,
        Edit as GraphIrEdit, Edits as GraphIrEdits, NumForm as GraphIrNumForm,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::error::{InvalidStructureError, TransactionError};

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

        assert_eq!(first.to_rust(), &initial);
        assert_eq!(second.to_rust(), &mol_dsl!(r#"{:atoms ["C" "N"]}"#));
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
        *built.to_rust_mut().atom_mut(GraphIrAtomId(0)).attributes =
            GraphIrAtomForm::from_element(ChemElement::N);
        let snapshot_error = editor.snapshot().unwrap_err();
        let build_error = editor.build().unwrap_err();

        assert_eq!(snapshot.to_rust(), &initial);
        assert_eq!(built.to_rust(), &mol_dsl!(r#"{:atoms ["N"]}"#));
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
    fn test_molecule_editor_publication_error() {
        let molecule = mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
        let mut snapshot_editor = MoleculeEditor {
            inner: Some(molecule.clone().edit()),
        };
        snapshot_editor.inner.as_mut().unwrap().add_bond(
            GraphIrAtomId(0),
            GraphIrAtomId(1),
            GraphIrBondForm::from_order(1),
        );
        let mut build_editor = MoleculeEditor {
            inner: Some(molecule.edit()),
        };
        build_editor.inner.as_mut().unwrap().add_bond(
            GraphIrAtomId(0),
            GraphIrAtomId(1),
            GraphIrBondForm::from_order(1),
        );

        let snapshot_error = snapshot_editor.snapshot().unwrap_err();
        let build_error = build_editor.build().unwrap_err();

        Python::attach(|py| {
            assert!(snapshot_error.is_instance_of::<InvalidStructureError>(py));
            assert_eq!(
                snapshot_error
                    .value(py)
                    .str()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "bond: parallel bonds on atoms [AtomId(0), AtomId(1)]"
            );
            assert!(build_error.is_instance_of::<InvalidStructureError>(py));
            assert_eq!(
                build_error
                    .value(py)
                    .str()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "bond: parallel bonds on atoms [AtomId(0), AtomId(1)]"
            );
        });
    }

    #[rstest]
    fn test_molecule_editor_transact(mut carbon_editor: MoleculeEditor, add_nitrogen: Edits) {
        Python::attach(|py| {
            let edits = Py::new(py, add_nitrogen).unwrap();

            let _transaction = carbon_editor.transact(py, edits).unwrap();

            assert_eq!(
                carbon_editor.snapshot().unwrap().to_rust(),
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
                editor.bind(py).borrow().snapshot().unwrap().to_rust(),
                &mol_dsl!(r#"{:atoms ["C" "N"]}"#)
            );

            transaction.rollback(py, editor.clone_ref(py)).unwrap();
            let second_error = transaction.rollback(py, editor.clone_ref(py)).unwrap_err();

            assert_eq!(
                editor.bind(py).borrow().snapshot().unwrap().to_rust(),
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
