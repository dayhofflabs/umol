//! Dense permutation witnesses at the Python boundary.

use std::fmt::Debug;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_graph_core::{
    Correspondence as GraphCoreCorrespondence, GraphRemapping, NodeId,
    Remapping as GraphCoreRemapping,
};
use umol_graph_ir::ir::MoleculeRemapping as GraphIrMoleculeRemapping;

use crate::correspondence::{Correspondence, MoleculeCorrespondence};

/// A read-only permutation of a dense integer id space.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Remapping(GraphCoreRemapping<NodeId>);

#[pymethods]
impl Remapping {
    /// Construct from images in source-id order.
    /// Raises ValueError for an out-of-range or repeated image.
    #[new]
    fn new(images: Vec<u32>) -> PyResult<Self> {
        GraphCoreRemapping::new(images.into_iter().map(NodeId).collect())
            .map(Self)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Images in source-id order.
    #[getter]
    fn images(&self) -> Vec<u32> {
        (0..self.0.len())
            .map(|idx| self.0.map(NodeId::from(idx)).0)
            .collect()
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Return the image, or None outside the source domain.
    fn try_map(&self, old: u32) -> Option<u32> {
        self.0.try_map(NodeId(old)).map(|id| id.0)
    }

    /// Return the image; raises ValueError outside the source domain.
    fn map(&self, old: u32) -> PyResult<u32> {
        self.try_map(old)
            .ok_or_else(|| PyValueError::new_err("id outside remapping source domain"))
    }

    /// Preserve all pairings and both counts as a correspondence.
    fn to_correspondence(&self) -> Correspondence {
        Correspondence::from_rust(&GraphCoreCorrespondence::from(&self.0))
    }

    fn __repr__(&self) -> String {
        format!("Remapping(images={:?})", self.images())
    }
}

impl Remapping {
    pub(crate) fn from_rust<Id: Copy + Into<usize> + From<usize>>(
        remapping: &GraphCoreRemapping<Id>,
    ) -> Self {
        Self(
            GraphCoreRemapping::new(
                (0..remapping.len())
                    .map(|idx| NodeId::from(Into::<usize>::into(remapping.map(Id::from(idx)))))
                    .collect(),
            )
            .expect("typed permutation preserves images"),
        )
    }

    pub(crate) fn to_rust<Id: Copy + Debug + Into<usize> + From<usize>>(
        &self,
    ) -> GraphCoreRemapping<Id> {
        GraphCoreRemapping::new(
            self.images()
                .into_iter()
                .map(|idx| Id::from(idx as usize))
                .collect(),
        )
        .expect("Python permutation preserves images")
    }
}

/// Read-only permutations for all eight molecule entity kinds.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeRemapping(GraphIrMoleculeRemapping);

#[pymethods]
impl MoleculeRemapping {
    /// Assemble eight validated permutations without binding them to a molecule.
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        atoms: &Remapping,
        bonds: &Remapping,
        dative_bonds: &Remapping,
        aromatic_systems: &Remapping,
        multicenter_bonds: &Remapping,
        noncovalent_bonds: &Remapping,
        stereo_atoms: &Remapping,
        stereo_bonds: &Remapping,
    ) -> Self {
        Self(GraphIrMoleculeRemapping::new(
            GraphRemapping::new(atoms.to_rust(), bonds.to_rust()),
            dative_bonds.to_rust(),
            aromatic_systems.to_rust(),
            multicenter_bonds.to_rust(),
            noncovalent_bonds.to_rust(),
            stereo_atoms.to_rust(),
            stereo_bonds.to_rust(),
        ))
    }

    #[getter]
    fn atoms(&self) -> Remapping {
        Remapping::from_rust(self.0.graph().nodes())
    }

    #[getter]
    fn bonds(&self) -> Remapping {
        Remapping::from_rust(self.0.graph().edges())
    }

    #[getter]
    fn dative_bonds(&self) -> Remapping {
        Remapping::from_rust(self.0.dative_bonds())
    }

    #[getter]
    fn aromatic_systems(&self) -> Remapping {
        Remapping::from_rust(self.0.aromatic_systems())
    }

    #[getter]
    fn multicenter_bonds(&self) -> Remapping {
        Remapping::from_rust(self.0.multicenter_bonds())
    }

    #[getter]
    fn noncovalent_bonds(&self) -> Remapping {
        Remapping::from_rust(self.0.noncovalent_bonds())
    }

    #[getter]
    fn stereo_atoms(&self) -> Remapping {
        Remapping::from_rust(self.0.stereo_atoms())
    }

    #[getter]
    fn stereo_bonds(&self) -> Remapping {
        Remapping::from_rust(self.0.stereo_bonds())
    }

    /// Preserve all eight pairings and counts as a molecule correspondence.
    fn to_correspondence(&self) -> MoleculeCorrespondence {
        MoleculeCorrespondence::from_rust(self.to_rust().into())
    }

    fn __repr__(&self) -> String {
        format!("MoleculeRemapping(atoms={}, bonds={}, dative_bonds={}, aromatic_systems={}, multicenter_bonds={}, noncovalent_bonds={}, stereo_atoms={}, stereo_bonds={})",
            self.atoms().__repr__(), self.bonds().__repr__(), self.dative_bonds().__repr__(), self.aromatic_systems().__repr__(), self.multicenter_bonds().__repr__(), self.noncovalent_bonds().__repr__(), self.stereo_atoms().__repr__(), self.stereo_bonds().__repr__()
        )
    }
}

impl MoleculeRemapping {
    pub(crate) fn from_rust(remapping: GraphIrMoleculeRemapping) -> Self {
        Self(remapping)
    }
    pub(crate) fn to_rust(&self) -> &GraphIrMoleculeRemapping {
        &self.0
    }
}
