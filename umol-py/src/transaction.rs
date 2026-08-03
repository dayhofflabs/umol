//! Python ownership wrappers for transactional molecule editing.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use umol_ast::ast::MoleculeEditor as AstMoleculeEditor;

use crate::molecule::MoleculeAst;

/// A mutable molecule editor that can be inspected before it is finalized.
#[pyclass]
pub struct MoleculeEditor {
    inner: Option<AstMoleculeEditor>,
}

#[pymethods]
impl MoleculeEditor {
    /// Materialize the editor's current state without consuming it.
    fn snapshot(&self) -> PyResult<MoleculeAst> {
        self.inner
            .as_ref()
            .map(|editor| MoleculeAst::from_inner(editor.snapshot()))
            .ok_or_else(consumed_editor_error)
    }

    /// Finalize the editor and consume its mutable state.
    fn build(&mut self) -> PyResult<MoleculeAst> {
        self.inner
            .take()
            .map(|editor| MoleculeAst::from_inner(editor.build()))
            .ok_or_else(consumed_editor_error)
    }
}

fn consumed_editor_error() -> PyErr {
    PyRuntimeError::new_err("molecule editor has been consumed")
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::PyRuntimeError;
    use rstest::rstest;
    use umol_ast::ast::{AtomAst as AstAtomAst, AtomId as AstAtomId};
    use umol_ast::mol_dsl;
    use umol_chem::element::Element as ChemElement;

    use super::*;

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
            .add_atom(AstAtomAst::from_element(ChemElement::N));
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
        *built.inner_mut().atom_mut(AstAtomId(0)).ast = AstAtomAst::from_element(ChemElement::N);
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
}
