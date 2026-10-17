//! Molecule editor ownership and transaction recovery selection.

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use umol_graph_ir::ir::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MoleculeEditor as GraphIrMoleculeEditor,
    MulticenterBondId, NoncovalentBondId, RollbackMethod as GraphIrRollbackMethod, StereoAtomId,
    StereoBondId,
};

use crate::compact::MoleculeCompaction;
use crate::correspondence::MoleculeCorrespondence;
use crate::edit::Edits;
use crate::error::{
    molecule_apply_error, molecule_integrity_error, transaction_error, ConsumedError,
};
use crate::molecule::Molecule;

/// Recovery method for molecule transactions and in-place resolution.
///
/// Clone edits a working copy; Journal records undo while editing the receiver.
/// Rejection or failure restores the receiver under normalized equality, preserving
/// participant order.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RollbackMethod {
    Clone,
    Journal,
}

impl RollbackMethod {
    pub(crate) fn to_rust(self) -> GraphIrRollbackMethod {
        match self {
            Self::Clone => GraphIrRollbackMethod::Clone,
            Self::Journal => GraphIrRollbackMethod::Journal,
        }
    }
}

/// An owning molecule editor, consumed by publication.
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
    /// Consume the editor and publish its molecule after checking integrity.
    ///
    /// Raises InvalidStructureError if publication fails. The editor remains consumed
    /// on either outcome; later access through any alias raises ConsumedError.
    fn finish(&mut self) -> PyResult<Molecule> {
        self.inner
            .take()
            .ok_or_else(consumed_editor_error)?
            .finish()
            .map(Molecule::from_rust)
            .map_err(molecule_integrity_error)
    }

    /// Apply an edit batch to this editor, checking each edit's preconditions.
    ///
    /// Returns None without checking molecule integrity or recording undo. Failure retains
    /// completed edits and leaves the editor available; a transferred batch stays consumed.
    fn apply(&mut self, py: Python<'_>, edits: Py<Edits>) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(consumed_editor_error)?
            .apply(edits.try_borrow_mut(py)?.take()?)
            .map_err(molecule_apply_error)
    }

    /// Apply a batch as in apply, returning its input-to-result correspondence.
    fn tracked_apply(
        &mut self,
        py: Python<'_>,
        edits: Py<Edits>,
    ) -> PyResult<MoleculeCorrespondence> {
        self.inner
            .as_mut()
            .ok_or_else(consumed_editor_error)?
            .tracked_apply(edits.try_borrow_mut(py)?.take()?)
            .map(MoleculeCorrespondence::from_rust)
            .map_err(molecule_apply_error)
    }

    /// Remove exactly the supplied atoms and bonds.
    ///
    /// Raises TransactionError before mutation for unavailable stored references or references
    /// that would be broken by removal. The editor remains available after rejection.
    fn remove_topology(&mut self, atoms: Vec<u32>, bonds: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&atoms, editor.atom_count(), "atom")?;
        ensure_in_range(&bonds, editor.bond_count(), "bond")?;
        let atoms = atoms.into_iter().map(AtomId).collect::<Vec<_>>();
        let bonds = bonds.into_iter().map(BondId).collect::<Vec<_>>();
        editor
            .remove_topology(&atoms, &bonds)
            .map_err(transaction_error)
    }

    /// Exact removal as in remove_topology, returning all eight source-to-result compactions.
    ///
    /// Unchanged entity kinds retain their source counts. Reference errors leave the editor
    /// unchanged and available; out-of-range removal IDs raise IndexError.
    fn tracked_remove_topology(
        &mut self,
        atoms: Vec<u32>,
        bonds: Vec<u32>,
    ) -> PyResult<MoleculeCompaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&atoms, editor.atom_count(), "atom")?;
        ensure_in_range(&bonds, editor.bond_count(), "bond")?;
        let atoms = atoms.into_iter().map(AtomId).collect::<Vec<_>>();
        let bonds = bonds.into_iter().map(BondId).collect::<Vec<_>>();
        editor
            .tracked_remove_topology(&atoms, &bonds)
            .map(MoleculeCompaction::from_rust)
            .map_err(transaction_error)
    }

    /// Remove topology, incident bonds, and overlays with unavailable stored participants or sites.
    ///
    /// Compacts surviving references and molecule constraints. The editor remains available;
    /// out-of-range removal IDs raise IndexError before mutation.
    fn remove_topology_cascading(&mut self, atoms: Vec<u32>, bonds: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&atoms, editor.atom_count(), "atom")?;
        ensure_in_range(&bonds, editor.bond_count(), "bond")?;
        let atoms = atoms.into_iter().map(AtomId).collect::<Vec<_>>();
        let bonds = bonds.into_iter().map(BondId).collect::<Vec<_>>();
        editor.remove_topology_cascading(&atoms, &bonds);
        Ok(())
    }

    /// Cascading removal as in remove_topology_cascading, returning all eight compactions.
    ///
    /// Unchanged entity kinds retain their source counts. The editor remains available;
    /// out-of-range removal IDs raise IndexError before mutation.
    fn tracked_remove_topology_cascading(
        &mut self,
        atoms: Vec<u32>,
        bonds: Vec<u32>,
    ) -> PyResult<MoleculeCompaction> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&atoms, editor.atom_count(), "atom")?;
        ensure_in_range(&bonds, editor.bond_count(), "bond")?;
        let atoms = atoms.into_iter().map(AtomId).collect::<Vec<_>>();
        let bonds = bonds.into_iter().map(BondId).collect::<Vec<_>>();
        Ok(MoleculeCompaction::from_rust(
            editor.tracked_remove_topology_cascading(&atoms, &bonds),
        ))
    }

    /// Remove dative bonds and compact that entity space.
    fn remove_dative_bonds(&mut self, ids: Vec<u32>) -> PyResult<()> {
        let editor = self.inner.as_mut().ok_or_else(consumed_editor_error)?;
        ensure_in_range(&ids, editor.dative_bond_count(), "dative bond")?;
        editor.remove_dative_bonds(&ids.into_iter().map(DativeBondId).collect::<Vec<_>>());
        Ok(())
    }

    /// Remove dative bonds and return all eight compactions, retaining unchanged source counts.
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

    /// Remove aromatic systems and return all eight compactions, retaining unchanged source counts.
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

    /// Remove multicenter bonds and return all eight compactions, retaining unchanged source counts.
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

    /// Remove noncovalent bonds and return all eight compactions, retaining unchanged source counts.
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

    /// Remove stereo atoms and return all eight compactions, retaining unchanged source counts.
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

    /// Remove stereo bonds and return all eight compactions, retaining unchanged source counts.
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
    use pyo3::exceptions::PyRuntimeError;
    use rstest::{fixture, rstest};
    use umol_chem::element::Element as ChemElement;
    use umol_graph_ir::ir::{
        AtomConstraintForm as GraphIrAtomConstraintForm, AtomForm as GraphIrAtomForm,
        AtomId as GraphIrAtomId, BondForm as GraphIrBondForm, Constraint as GraphIrConstraint,
        DativeBondForm as GraphIrDativeBondForm, Edits as GraphIrEdits, Entity as GraphIrEntity,
        Molecule as GraphIrMolecule, MoleculeEntries as GraphIrMoleculeEntries,
        MoleculeIntegrityError as GraphIrMoleculeIntegrityError,
        StereoAtomConstraintForm as GraphIrStereoAtomConstraintForm,
        StereoAtomId as GraphIrStereoAtomId, StereoKind as GraphIrStereoKind,
        StereogenicityForm as GraphIrStereogenicityForm,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::error::{InvalidStructureError, TransactionError};

    #[rstest]
    #[case(RollbackMethod::Clone, GraphIrRollbackMethod::Clone)]
    #[case(RollbackMethod::Journal, GraphIrRollbackMethod::Journal)]
    fn test_rollback_method_to_rust(
        #[case] method: RollbackMethod,
        #[case] expected: GraphIrRollbackMethod,
    ) {
        assert_eq!(method.to_rust(), expected);
    }

    #[fixture]
    fn carbon_editor() -> MoleculeEditor {
        MoleculeEditor {
            inner: Some(mol_dsl!(r#"{:atoms ["C"]}"#).edit()),
        }
    }

    #[rstest]
    fn test_molecule_editor_finish(mut carbon_editor: MoleculeEditor) {
        let mut molecule = carbon_editor.finish().unwrap();
        assert_eq!(molecule.to_rust(), &mol_dsl!(r#"{:atoms ["C"]}"#));
        *molecule
            .to_rust_mut()
            .atom_mut(GraphIrAtomId(0))
            .attributes_mut() = GraphIrAtomForm::from_element(ChemElement::N);
        assert_eq!(molecule.to_rust(), &mol_dsl!(r#"{:atoms ["N"]}"#));
        Python::attach(|py| {
            let error = carbon_editor.finish().unwrap_err();
            assert!(error.is_instance_of::<ConsumedError>(py));
            assert_eq!(
                error.to_string(),
                "ConsumedError: MoleculeEditor has been consumed"
            );
        });
    }

    #[rstest]
    fn test_molecule_editor_finish_error() {
        let molecule = mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
        let mut editor = MoleculeEditor::from_rust(molecule.edit());
        editor.inner.as_mut().unwrap().add_bond(
            GraphIrAtomId(0),
            GraphIrAtomId(1),
            GraphIrBondForm::from_order(1),
        );
        Python::attach(|py| {
            let error = editor.finish().unwrap_err();
            assert!(error.is_instance_of::<InvalidStructureError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "bond: parallel bonds on atoms [AtomId(0), AtomId(1)]"
            );
            assert!(editor
                .finish()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
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
    fn test_molecule_editor_finish_constraint_error(
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
        let error = editor.finish().unwrap_err();
        Python::attach(|py| {
            assert!(error.is_instance_of::<InvalidStructureError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected.to_string()
            );
        });
    }

    #[rstest]
    #[case::plain(false)]
    #[case::tracked(true)]
    fn test_molecule_editor_apply_borrowed_edits(
        mut carbon_editor: MoleculeEditor,
        #[case] tracked: bool,
    ) {
        Python::attach(|py| {
            let mut changes = GraphIrEdits::new();
            changes.add_atom(GraphIrAtomForm::from_element(ChemElement::N));
            let edits = Py::new(py, Edits::from_rust(changes.clone())).unwrap();
            let borrowed = edits.borrow(py);
            let result = if tracked {
                carbon_editor
                    .tracked_apply(py, edits.clone_ref(py))
                    .map(|_| ())
            } else {
                carbon_editor.apply(py, edits.clone_ref(py))
            };
            let error = result.unwrap_err();
            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "Already borrowed"
            );
            assert_eq!(borrowed.to_rust().unwrap(), &changes);
            assert_eq!(
                carbon_editor.inner.as_ref().unwrap().probe(),
                Ok(&mol_dsl!(r#"{:atoms ["C"]}"#))
            );
            drop(borrowed);

            carbon_editor.apply(py, edits.clone_ref(py)).unwrap();
            assert_eq!(
                carbon_editor.finish().unwrap().to_rust(),
                &mol_dsl!(r#"{:atoms ["C" "N"]}"#)
            );
            assert!(edits
                .borrow(py)
                .to_rust()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
        });
    }

    #[rstest]
    fn test_molecule_editor_remove_topology_reference(
        mut carbon_editor: MoleculeEditor,
        #[values(false, true)] tracked: bool,
    ) {
        carbon_editor.inner.as_mut().unwrap().add_dative_bond(
            &[GraphIrAtomId(2)],
            GraphIrAtomId(0),
            GraphIrDativeBondForm::from_order(1),
        );
        let result = if tracked {
            carbon_editor
                .tracked_remove_topology(vec![0], vec![])
                .map(|_| ())
        } else {
            carbon_editor.remove_topology(vec![0], vec![])
        };
        Python::attach(|py| {
            let error = result.unwrap_err();
            assert!(error.is_instance_of::<TransactionError>(py));
            assert_eq!(
                error.to_string(),
                "TransactionError: invalid reference to atom 2"
            );
        });
        carbon_editor.inner.as_mut().unwrap().add_atoms(vec![
            GraphIrAtomForm::from_element(ChemElement::N),
            GraphIrAtomForm::from_element(ChemElement::O),
        ]);
        let expected = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::N),
                GraphIrAtomForm::from_element(ChemElement::O),
            ],
            dative: vec![(
                vec![GraphIrAtomId(2)],
                GraphIrAtomId(0),
                GraphIrDativeBondForm::from_order(1),
            )],
            ..Default::default()
        });
        assert_eq!(carbon_editor.finish().unwrap().to_rust(), &expected);
    }
}
