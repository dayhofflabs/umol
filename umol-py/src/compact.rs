//! Count-bearing compaction witnesses at the Python boundary.

use std::fmt::Debug;
use std::ops::{Add, Sub};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_graph_core::{
    Compaction as GraphCoreCompaction, Correspondence as GraphCoreCorrespondence, GraphCompaction,
    NodeId,
};
use umol_graph_ir::ir::{
    AromaticSystemId, AtomId, BondId, DativeBondId,
    MoleculeCompaction as GraphIrMoleculeCompaction,
    MoleculeCorrespondence as GraphIrMoleculeCorrespondence, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

use crate::correspondence::{Correspondence, MoleculeCorrespondence};

trait CompactionId:
    Copy
    + Debug
    + Ord
    + Into<usize>
    + From<usize>
    + Add<usize, Output = Self>
    + Sub<usize, Output = Self>
{
}

impl<T> CompactionId for T where
    T: Copy
        + Debug
        + Ord
        + Into<usize>
        + From<usize>
        + Add<usize, Output = Self>
        + Sub<usize, Output = Self>
{
}

/// An order-preserving removal from one finite integer id space.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Compaction(GraphCoreCompaction<NodeId>);

#[pymethods]
impl Compaction {
    /// Construct from a source count and the removed ids.
    ///
    /// Input order and repetitions do not matter. Raises `ValueError` when a removed id is outside
    /// the declared source space.
    #[new]
    fn new(source_count: usize, removed: Vec<u32>) -> PyResult<Self> {
        GraphCoreCompaction::new(source_count, removed.into_iter().map(NodeId).collect())
            .map(Self)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// No removals from an empty source space.
    #[staticmethod]
    fn empty() -> Self {
        Self(GraphCoreCompaction::empty())
    }

    /// No removals from a source space of `source_count` ids.
    #[staticmethod]
    fn identity(source_count: usize) -> Self {
        Self(GraphCoreCompaction::identity(source_count))
    }

    /// Number of ids before removal.
    #[getter]
    fn source_count(&self) -> usize {
        self.0.source_count()
    }

    /// Number of surviving ids.
    #[getter]
    fn result_count(&self) -> usize {
        self.0.result_count()
    }

    /// Removed ids in ascending order, without repetitions.
    #[getter]
    fn removed(&self) -> Vec<u32> {
        self.0.removed().iter().map(|id| id.0).collect()
    }

    /// Return the survivor's result id, or `None` when it was removed or is out of range.
    fn compact(&self, old: u32) -> Option<u32> {
        self.0.compact(NodeId(old)).map(|id| id.0)
    }

    /// Return the survivor's source id, or `None` outside the result space.
    fn try_uncompact(&self, post: u32) -> Option<u32> {
        self.0.try_uncompact(NodeId(post)).map(|id| id.0)
    }

    /// Return the survivor's source id; raises `ValueError` outside the result space.
    fn uncompact(&self, post: u32) -> PyResult<u32> {
        self.try_uncompact(post)
            .ok_or_else(|| PyValueError::new_err("id outside compaction result domain"))
    }

    /// Preserve every surviving pairing and both counts as a correspondence.
    fn to_correspondence(&self) -> Correspondence {
        Correspondence::from_rust(&GraphCoreCorrespondence::from(&self.0))
    }

    fn __repr__(&self) -> String {
        format!(
            "Compaction(source_count={}, removed={:?})",
            self.source_count(),
            self.removed()
        )
    }
}

impl Compaction {
    fn from_rust<Id: CompactionId>(compaction: &GraphCoreCompaction<Id>) -> Self {
        Self(
            GraphCoreCompaction::new(
                compaction.source_count(),
                compaction
                    .removed()
                    .iter()
                    .map(|&id| NodeId::from(Into::<usize>::into(id)))
                    .collect(),
            )
            .expect("typed compaction preserves its source domain"),
        )
    }

    fn to_rust<Id: CompactionId>(&self) -> GraphCoreCompaction<Id> {
        GraphCoreCompaction::new(
            self.0.source_count(),
            self.0
                .removed()
                .iter()
                .map(|&id| Id::from(id.index()))
                .collect(),
        )
        .expect("Python compaction preserves its source domain")
    }
}

/// Read-only compactions for all eight molecule entity kinds.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeCompaction(GraphIrMoleculeCompaction);

#[pymethods]
impl MoleculeCompaction {
    /// Assemble eight validated compactions without binding them to a molecule.
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        atoms: &Compaction,
        bonds: &Compaction,
        dative_bonds: &Compaction,
        aromatic_systems: &Compaction,
        multicenter_bonds: &Compaction,
        noncovalent_bonds: &Compaction,
        stereo_atoms: &Compaction,
        stereo_bonds: &Compaction,
    ) -> Self {
        Self(GraphIrMoleculeCompaction::new(
            GraphCompaction::new(atoms.to_rust(), bonds.to_rust()),
            dative_bonds.to_rust(),
            aromatic_systems.to_rust(),
            multicenter_bonds.to_rust(),
            noncovalent_bonds.to_rust(),
            stereo_atoms.to_rust(),
            stereo_bonds.to_rust(),
        ))
    }

    /// Empty compactions for all eight entity kinds.
    #[staticmethod]
    fn empty() -> Self {
        Self(GraphIrMoleculeCompaction::empty())
    }

    #[getter]
    fn atoms(&self) -> Compaction {
        Compaction::from_rust(self.0.graph().nodes())
    }

    #[getter]
    fn bonds(&self) -> Compaction {
        Compaction::from_rust(self.0.graph().edges())
    }

    #[getter]
    fn dative_bonds(&self) -> Compaction {
        Compaction::from_rust(self.0.dative_bonds())
    }

    #[getter]
    fn aromatic_systems(&self) -> Compaction {
        Compaction::from_rust(self.0.aromatic_systems())
    }

    #[getter]
    fn multicenter_bonds(&self) -> Compaction {
        Compaction::from_rust(self.0.multicenter_bonds())
    }

    #[getter]
    fn noncovalent_bonds(&self) -> Compaction {
        Compaction::from_rust(self.0.noncovalent_bonds())
    }

    #[getter]
    fn stereo_atoms(&self) -> Compaction {
        Compaction::from_rust(self.0.stereo_atoms())
    }

    #[getter]
    fn stereo_bonds(&self) -> Compaction {
        Compaction::from_rust(self.0.stereo_bonds())
    }

    fn compact_atom(&self, id: u32) -> Option<u32> {
        self.0.compact_atom(AtomId(id)).map(|id| id.0)
    }

    fn compact_bond(&self, id: u32) -> Option<u32> {
        self.0.compact_bond(BondId(id)).map(|id| id.0)
    }

    fn compact_dative_bond(&self, id: u32) -> Option<u32> {
        self.0.compact_dative_bond(DativeBondId(id)).map(|id| id.0)
    }

    fn compact_aromatic_system(&self, id: u32) -> Option<u32> {
        self.0
            .compact_aromatic_system(AromaticSystemId(id))
            .map(|id| id.0)
    }

    fn compact_multicenter_bond(&self, id: u32) -> Option<u32> {
        self.0
            .compact_multicenter_bond(MulticenterBondId(id))
            .map(|id| id.0)
    }

    fn compact_noncovalent_bond(&self, id: u32) -> Option<u32> {
        self.0
            .compact_noncovalent_bond(NoncovalentBondId(id))
            .map(|id| id.0)
    }

    fn compact_stereo_atom(&self, id: u32) -> Option<u32> {
        self.0.compact_stereo_atom(StereoAtomId(id)).map(|id| id.0)
    }

    fn compact_stereo_bond(&self, id: u32) -> Option<u32> {
        self.0.compact_stereo_bond(StereoBondId(id)).map(|id| id.0)
    }

    /// Preserve all eight survivor mappings and counts as a molecule correspondence.
    fn to_correspondence(&self) -> MoleculeCorrespondence {
        MoleculeCorrespondence::from_rust(GraphIrMoleculeCorrespondence::from(&self.0))
    }

    fn __repr__(&self) -> String {
        format!(
            concat!(
                "MoleculeCompaction(",
                "atoms={}, bonds={}, dative_bonds={}, aromatic_systems={}, ",
                "multicenter_bonds={}, noncovalent_bonds={}, stereo_atoms={}, ",
                "stereo_bonds={})"
            ),
            self.atoms().__repr__(),
            self.bonds().__repr__(),
            self.dative_bonds().__repr__(),
            self.aromatic_systems().__repr__(),
            self.multicenter_bonds().__repr__(),
            self.noncovalent_bonds().__repr__(),
            self.stereo_atoms().__repr__(),
            self.stereo_bonds().__repr__(),
        )
    }
}

impl MoleculeCompaction {
    pub(crate) fn from_rust(compaction: GraphIrMoleculeCompaction) -> Self {
        Self(compaction)
    }
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};

    use super::*;

    #[fixture]
    fn compaction() -> Compaction {
        Compaction::new(5, vec![3, 1, 3]).unwrap()
    }

    #[rstest]
    fn test_compaction_new(compaction: Compaction) {
        assert_eq!(compaction.source_count(), 5);
        assert_eq!(compaction.result_count(), 3);
        assert_eq!(compaction.removed(), vec![1, 3]);
        assert_eq!(
            compaction.__repr__(),
            "Compaction(source_count=5, removed=[1, 3])"
        );
    }

    #[rstest]
    fn test_compaction_new_error() {
        let error = Compaction::new(2, vec![2]).unwrap_err();
        Python::attach(|py| {
            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().to_string(),
                "removed id NodeId(2) is out of range for 2 entries"
            );
        });
    }

    #[rstest]
    fn test_compaction_maps_survivors(compaction: Compaction) {
        assert_eq!(compaction.compact(0), Some(0));
        assert_eq!(compaction.compact(1), None);
        assert_eq!(compaction.compact(4), Some(2));
        assert_eq!(compaction.compact(5), None);
        assert_eq!(compaction.try_uncompact(2), Some(4));
        assert_eq!(compaction.try_uncompact(3), None);
        assert_eq!(compaction.uncompact(2).unwrap(), 4);
    }

    #[rstest]
    fn test_compaction_uncompact_error(compaction: Compaction) {
        let error = compaction.uncompact(3).unwrap_err();
        Python::attach(|py| {
            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().to_string(),
                "id outside compaction result domain"
            );
        });
    }

    #[rstest]
    fn test_compaction_to_correspondence(compaction: Compaction) {
        let expected = GraphCoreCorrespondence::new(
            vec![
                (NodeId(0), NodeId(0)),
                (NodeId(2), NodeId(1)),
                (NodeId(4), NodeId(2)),
            ],
            5,
            3,
        )
        .unwrap();

        assert_eq!(
            compaction.to_correspondence(),
            Correspondence::from_rust(&expected)
        );
    }

    #[rstest]
    fn test_molecule_compaction_components() {
        let atoms = Compaction::new(3, vec![1]).unwrap();
        let identity = Compaction::identity(2);
        let empty = Compaction::empty();
        let compaction = MoleculeCompaction::new(
            &atoms, &identity, &empty, &empty, &empty, &empty, &empty, &empty,
        );

        assert_eq!(compaction.atoms(), atoms);
        assert_eq!(compaction.bonds(), identity);
        assert_eq!(compaction.compact_atom(2), Some(1));
        assert_eq!(compaction.compact_atom(1), None);
        let correspondence = compaction.to_correspondence();
        assert_eq!(
            correspondence.to_rust().atoms().matched_pairs(),
            &[(AtomId(0), AtomId(0)), (AtomId(2), AtomId(1))]
        );
    }
}
