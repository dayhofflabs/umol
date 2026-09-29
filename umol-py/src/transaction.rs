//! Python ownership wrapper for the molecule editor.

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use umol_graph_ir::ir::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MoleculeEditor as GraphIrMoleculeEditor,
    MulticenterBondId, NoncovalentBondId, StereoAtomId, StereoBondId,
};

use crate::correspondence::MoleculeCorrespondence;
use crate::edit::Edits;
use crate::error::{molecule_integrity_error, transaction_error, ConsumedError};
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
            .apply(edits.try_borrow_mut(py)?.take()?)
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
            .tracked_apply(edits.try_borrow_mut(py)?.take()?)
            .map(|(editor, correspondence)| {
                (
                    Self::from_rust(editor),
                    MoleculeCorrespondence::from_rust(correspondence),
                )
            })
            .map_err(transaction_error)
    }

    /// Remove atoms and bonds, cascading dependent entities.
    fn remove_topology(&mut self, atoms: Vec<u32>, bonds: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&atoms, editor.atom_count(), "atom")?;
        ensure_in_range(&bonds, editor.bond_count(), "bond")?;
        let atoms = atoms.into_iter().map(AtomId).collect::<Vec<_>>();
        let bonds = bonds.into_iter().map(BondId).collect::<Vec<_>>();
        editor.remove_topology(&atoms, &bonds);
        Ok(())
    }

    /// Remove dative bonds and compact that entity space.
    fn remove_dative_bonds(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.dative_bond_count(), "dative bond")?;
        editor.remove_dative_bonds(&ids.into_iter().map(DativeBondId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove aromatic systems and compact that entity space.
    fn remove_aromatic_systems(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.aromatic_system_count(), "aromatic system")?;
        editor.remove_aromatic_systems(&ids.into_iter().map(AromaticSystemId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove multicenter bonds and compact that entity space.
    fn remove_multicenter_bonds(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.multicenter_bond_count(), "multicenter bond")?;
        editor
            .remove_multicenter_bonds(&ids.into_iter().map(MulticenterBondId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove noncovalent bonds and compact that entity space.
    fn remove_noncovalent_bonds(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.noncovalent_bond_count(), "noncovalent bond")?;
        editor
            .remove_noncovalent_bonds(&ids.into_iter().map(NoncovalentBondId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove stereo atoms and compact that entity space.
    fn remove_stereo_atoms(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.stereo_atom_count(), "stereo atom")?;
        editor.remove_stereo_atoms(&ids.into_iter().map(StereoAtomId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove stereo bonds and compact that entity space.
    fn remove_stereo_bonds(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.stereo_bond_count(), "stereo bond")?;
        editor.remove_stereo_bonds(&ids.into_iter().map(StereoBondId).collect::<Vec<_>>());
        Ok(())
    }
}

fn ensure_in_range(ids: &[u32], count: usize, entity: &str) -> PyResult<()> {
    if ids.iter().any(|&id| id as usize >= count) {
        return Err(PyIndexError::new_err(format!("{entity} id out of range")));
    }
    Ok(())
}

fn consumed_editor_error() -> PyErr {
    ConsumedError::new_err("MoleculeEditor has been consumed")
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_chem::element::Element as ChemElement;
    use umol_graph_ir::ir::{
        AtomConstraintForm as GraphIrAtomConstraintForm, AtomForm as GraphIrAtomForm,
        AtomId as GraphIrAtomId, BondForm as GraphIrBondForm, Constraint as GraphIrConstraint,
        Entity as GraphIrEntity, MoleculeIntegrityError as GraphIrMoleculeIntegrityError,
        StereoAtomConstraintForm as GraphIrStereoAtomConstraintForm,
        StereoAtomId as GraphIrStereoAtomId, StereoKind as GraphIrStereoKind,
        StereogenicityForm as GraphIrStereogenicityForm,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::error::InvalidStructureError;

    #[fixture]
    fn carbon_editor() -> MoleculeEditor {
        MoleculeEditor {
            inner: Some(mol_dsl!(r#"{:atoms ["C"]}"#).edit()),
        }
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
        *built
            .to_rust_mut()
            .atom_mut(GraphIrAtomId(0))
            .attributes_mut() = GraphIrAtomForm::from_element(ChemElement::N);
        let snapshot_error = editor.snapshot().unwrap_err();
        let build_error = editor.build().unwrap_err();

        assert_eq!(snapshot.to_rust(), &initial);
        assert_eq!(built.to_rust(), &mol_dsl!(r#"{:atoms ["N"]}"#));
        Python::attach(|py| {
            assert!(snapshot_error.is_instance_of::<ConsumedError>(py));
            assert_eq!(
                snapshot_error
                    .value(py)
                    .str()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "MoleculeEditor has been consumed"
            );
            assert!(build_error.is_instance_of::<ConsumedError>(py));
            assert_eq!(
                build_error
                    .value(py)
                    .str()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "MoleculeEditor has been consumed"
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
    #[case::reference(
        GraphIrConstraint::Atom(GraphIrAtomId(5), GraphIrAtomConstraintForm::degree(1)),
        GraphIrMoleculeIntegrityError::InvalidReference { entity: GraphIrEntity::Atom(GraphIrAtomId(5)) }
    )]
    #[case::stereo_frame(
        GraphIrConstraint::StereoAtom(
            GraphIrStereoAtomId(0), GraphIrStereoKind::Octahedral,
            GraphIrStereoAtomConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Undetermined)
        ),
        GraphIrMoleculeIntegrityError::StereoLigandArity {
            entity: GraphIrEntity::StereoAtom(GraphIrStereoAtomId(0)),
            kind: GraphIrStereoKind::Octahedral, expected: 6, actual: 4
        }
    )]
    fn test_molecule_editor_build_constraint_error(
        #[case] constraint: GraphIrConstraint,
        #[case] expected: GraphIrMoleculeIntegrityError,
    ) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
        );
        let mut editor = MoleculeEditor {
            inner: Some(molecule.edit()),
        };
        editor
            .inner
            .as_mut()
            .unwrap()
            .constraints_mut()
            .push(constraint);
        let error = editor.build().unwrap_err();
        Python::attach(|py| {
            assert!(error.is_instance_of::<InvalidStructureError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected.to_string()
            );
        });
    }
}
