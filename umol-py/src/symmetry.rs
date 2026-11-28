//! Graph-symmetry auxiliary data and molecule-backed readings.

use pyo3::prelude::*;
use umol_graph_ir::ir::{
    AtomId as GraphIrAtomId, GraphSymmetry as GraphIrGraphSymmetry,
    GraphSymmetryConfig as GraphIrGraphSymmetryConfig, Molecule as GraphIrMolecule,
};

use crate::algorithm::AutomorphismAlgorithm;
use crate::error::ConsumedError;
use crate::metadata::Entity;
use crate::molecule::Molecule;
use crate::stereo::StereogenicSiteKind;

/// Owned graph-symmetry data, without a molecule or stored query results.
/// Obtain it by consuming GraphSymmetryView with into_graph_symmetry().
#[pyclass(frozen)]
pub struct GraphSymmetry(GraphIrGraphSymmetry);

impl GraphSymmetry {
    pub(crate) fn from_rust(symmetry: GraphIrGraphSymmetry) -> Self {
        Self(symmetry)
    }

    pub(crate) fn to_rust(&self) -> &GraphIrGraphSymmetry {
        &self.0
    }
}

#[pymethods]
impl GraphSymmetry {
    /// Explicitly copy the auxiliary data for the same molecular representation.
    fn copy(&self) -> Self {
        Self::from_rust(self.0.clone())
    }
}

/// Graph-symmetry readings retaining a live molecule and owned or borrowed auxiliary data.
/// Whole-molecule mutation invalidates this view. Extraction consumes it.
#[pyclass]
pub struct GraphSymmetryView {
    molecule: Py<Molecule>,
    counter: u64,
    storage: Option<GraphSymmetryStorage>,
}

enum GraphSymmetryStorage {
    Owned(GraphIrGraphSymmetry),
    Borrowed(Py<GraphSymmetry>),
}

impl GraphSymmetryView {
    pub(crate) fn from_rust(
        py: Python<'_>,
        symmetry: GraphIrGraphSymmetry,
        molecule: Py<Molecule>,
    ) -> PyResult<Self> {
        let counter = molecule.try_borrow(py)?.view_counter();
        Ok(Self {
            molecule,
            counter,
            storage: Some(GraphSymmetryStorage::Owned(symmetry)),
        })
    }

    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&GraphIrGraphSymmetry, &GraphIrMolecule) -> PyResult<R>,
    ) -> PyResult<R> {
        let storage = self
            .storage
            .as_ref()
            .ok_or_else(|| ConsumedError::new_err("GraphSymmetryView has been consumed"))?;
        let molecule = self.molecule.try_borrow(py)?;
        molecule.check_access(self.counter, "GraphSymmetryView")?;
        match storage {
            GraphSymmetryStorage::Owned(symmetry) => f(symmetry, molecule.to_rust()),
            GraphSymmetryStorage::Borrowed(owner) => {
                f(owner.try_borrow(py)?.to_rust(), molecule.to_rust())
            }
        }
    }
}

#[pymethods]
impl GraphSymmetryView {
    /// Pair existing data with a molecule without copying or recomputing either.
    /// The data must describe the molecule; that relationship is not validated.
    #[staticmethod]
    fn from_graph_symmetry(
        py: Python<'_>,
        symmetry: Py<GraphSymmetry>,
        molecule: Py<Molecule>,
    ) -> PyResult<Self> {
        let counter = molecule.try_borrow(py)?.view_counter();
        Ok(Self {
            molecule,
            counter,
            storage: Some(GraphSymmetryStorage::Borrowed(symmetry)),
        })
    }

    /// Consume this view, moving owned data or cloning borrowed data.
    /// The molecule and any source GraphSymmetry remain usable.
    #[allow(
        clippy::wrong_self_convention,
        reason = "Consumes storage while retaining the Python owner in its consumed state"
    )]
    fn into_graph_symmetry(&mut self, py: Python<'_>) -> PyResult<GraphSymmetry> {
        self.read(py, |_, _| Ok(()))?;
        match self
            .storage
            .take()
            .ok_or_else(|| ConsumedError::new_err("GraphSymmetryView has been consumed"))?
        {
            GraphSymmetryStorage::Owned(symmetry) => Ok(GraphSymmetry::from_rust(symmetry)),
            GraphSymmetryStorage::Borrowed(owner) => Ok(GraphSymmetry::from_rust(
                owner.try_borrow(py)?.to_rust().clone(),
            )),
        }
    }

    /// Whether the molecule has no remapping to its mirror image.
    fn is_chiral(&self, py: Python<'_>) -> PyResult<bool> {
        self.read(py, |symmetry, molecule| Ok(molecule.is_chiral(symmetry)))
    }

    /// Whether a preserving remapping sends a to b; false for unavailable ids.
    fn same_proper_orbit(&self, py: Python<'_>, a: u32, b: u32) -> PyResult<bool> {
        self.read(py, |symmetry, molecule| {
            Ok(molecule.same_proper_orbit(symmetry, GraphIrAtomId(a), GraphIrAtomId(b)))
        })
    }

    /// Whether a preserving or mirror remapping sends a to b; false for unavailable ids.
    fn same_star_orbit(&self, py: Python<'_>, a: u32, b: u32) -> PyResult<bool> {
        self.read(py, |symmetry, molecule| {
            Ok(molecule.same_star_orbit(symmetry, GraphIrAtomId(a), GraphIrAtomId(b)))
        })
    }

    /// The proper orbit in atom-id order, empty for an unavailable id.
    fn proper_orbit_of(&self, py: Python<'_>, a: u32) -> PyResult<Vec<u32>> {
        self.read(py, |symmetry, molecule| {
            Ok(molecule
                .proper_orbit_of(symmetry, GraphIrAtomId(a))
                .into_iter()
                .map(|id| id.0)
                .collect())
        })
    }

    /// The proper-or-mirror orbit in atom-id order, empty for an unavailable id.
    fn star_orbit_of(&self, py: Python<'_>, a: u32) -> PyResult<Vec<u32>> {
        self.read(py, |symmetry, molecule| {
            Ok(molecule
                .star_orbit_of(symmetry, GraphIrAtomId(a))
                .into_iter()
                .map(|id| id.0)
                .collect())
        })
    }

    /// Classify all kinded stereo records in stereo-record id order.
    /// Kindless records are omitted. The returned list belongs to the caller.
    fn stereogenic_site_kinds(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<(Entity, StereogenicSiteKind)>> {
        self.read(py, |symmetry, molecule| {
            Ok(molecule
                .stereogenic_site_kinds(symmetry)
                .into_iter()
                .map(|(entity, kind)| {
                    (
                        Entity::from_rust(entity),
                        StereogenicSiteKind::from_rust(kind),
                    )
                })
                .collect())
        })
    }
}

/// Algorithm, classification-trial limit and constraint participation for graph symmetry.
/// The three fields are independently writable, matching the Rust config.
#[pyclass(eq, from_py_object)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GraphSymmetryConfig(GraphIrGraphSymmetryConfig);

#[pymethods]
impl GraphSymmetryConfig {
    #[new]
    #[pyo3(signature = (*, automorphism_algorithm=AutomorphismAlgorithm::TrailIr(), max_classification_trials=4096, include_constraints=false))]
    fn new(
        automorphism_algorithm: AutomorphismAlgorithm,
        max_classification_trials: u32,
        include_constraints: bool,
    ) -> Self {
        Self::from_rust(GraphIrGraphSymmetryConfig {
            automorphism_algorithm: automorphism_algorithm.to_rust(),
            max_classification_trials,
            include_constraints,
        })
    }

    #[staticmethod]
    fn default() -> Self {
        Default::default()
    }

    #[getter]
    fn automorphism_algorithm(&self) -> AutomorphismAlgorithm {
        AutomorphismAlgorithm::from_rust(self.0.automorphism_algorithm)
    }

    #[setter]
    fn set_automorphism_algorithm(&mut self, algorithm: AutomorphismAlgorithm) {
        self.0.automorphism_algorithm = algorithm.to_rust();
    }

    #[getter]
    fn max_classification_trials(&self) -> u32 {
        self.0.max_classification_trials
    }

    #[setter]
    fn set_max_classification_trials(&mut self, limit: u32) {
        self.0.max_classification_trials = limit;
    }

    #[getter]
    fn include_constraints(&self) -> bool {
        self.0.include_constraints
    }

    #[setter]
    fn set_include_constraints(&mut self, include: bool) {
        self.0.include_constraints = include;
    }

    fn copy(&self) -> Self {
        *self
    }

    fn __repr__(&self) -> String {
        format!(
            "GraphSymmetryConfig(automorphism_algorithm={}, max_classification_trials={}, include_constraints={})",
            self.automorphism_algorithm().repr(), self.0.max_classification_trials,
            if self.0.include_constraints { "True" } else { "False" },
        )
    }
}

impl GraphSymmetryConfig {
    pub(crate) fn from_rust(config: GraphIrGraphSymmetryConfig) -> Self {
        Self(config)
    }

    pub(crate) fn to_rust(self) -> GraphIrGraphSymmetryConfig {
        self.0
    }
}

#[cfg(test)]
mod tests;
