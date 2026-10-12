//! Ring data, enumeration configuration, and owner-backed ring views.

use pyo3::exceptions::{PyIndexError, PyTypeError};
use pyo3::prelude::*;
use umol_graph_ir::ir::{
    AtomId as GraphIrAtomId, BondId as GraphIrBondId, Molecule as GraphIrMolecule,
    RingAtomView as GraphIrRingAtomView, RingBondView as GraphIrRingBondView,
    RingConfig as GraphIrRingConfig, RingConnection as GraphIrRingConnection,
    RingGraph as GraphIrRingGraph, RingId as GraphIrRingId, RingModel as GraphIrRingModel,
    RingRelation as GraphIrRingRelation, RingSet as GraphIrRingSet,
    RingSetKind as GraphIrRingSetKind, RingView as GraphIrRingView, RingViews as GraphIrRingViews,
};

use crate::algorithm::{RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm};
use crate::constraint::ring::RingScope;
use crate::convert::hash_rust;
use crate::error::{ConsumedError, InvalidatedViewError};
use crate::molecule::Molecule;
use crate::num::NumForm;

/// Which cycles constitute the ring set.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RingSetKind {
    Simple,
    Relevant,
}

#[pymethods]
impl RingSetKind {
    fn __hash__(&self) -> u64 {
        hash_rust(self)
    }
}

impl RingSetKind {
    fn from_rust(kind: GraphIrRingSetKind) -> Self {
        match kind {
            GraphIrRingSetKind::Simple => Self::Simple,
            GraphIrRingSetKind::Relevant => Self::Relevant,
        }
    }

    fn to_rust(self) -> GraphIrRingSetKind {
        match self {
            Self::Simple => GraphIrRingSetKind::Simple,
            Self::Relevant => GraphIrRingSetKind::Relevant,
        }
    }
}

/// Ring-set semantics and maximum ring size. Fields are independently writable.
#[pyclass(eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RingModel {
    #[pyo3(get, set)]
    kind: RingSetKind,
    #[pyo3(get, set)]
    max_ring_size: usize,
}

#[pymethods]
impl RingModel {
    #[new]
    #[pyo3(signature = (*, kind=RingSetKind::Relevant, max_ring_size=22))]
    fn new(kind: RingSetKind, max_ring_size: usize) -> Self {
        Self {
            kind,
            max_ring_size,
        }
    }

    fn copy(&self) -> Self {
        *self
    }

    fn __repr__(&self) -> String {
        format!(
            "RingModel(kind=RingSetKind.{:?}, max_ring_size={})",
            self.kind, self.max_ring_size
        )
    }
}

impl RingModel {
    pub(crate) fn to_rust(self) -> GraphIrRingModel {
        GraphIrRingModel {
            kind: self.kind.to_rust(),
            max_ring_size: self.max_ring_size,
        }
    }
}

/// Algorithms used to compute each supported ring-set kind.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RingConfig(GraphIrRingConfig);

#[pymethods]
impl RingConfig {
    #[new]
    #[pyo3(signature = (*, simple_cycle_algorithm=None, relevant_cycle_algorithm=None))]
    fn new(
        simple_cycle_algorithm: Option<SimpleCycleEnumerationAlgorithm>,
        relevant_cycle_algorithm: Option<RelevantCycleEnumerationAlgorithm>,
    ) -> Self {
        let defaults = Self::default().0;
        Self(GraphIrRingConfig {
            simple_cycle_algorithm: simple_cycle_algorithm
                .map_or(defaults.simple_cycle_algorithm, |algorithm| {
                    algorithm.to_rust()
                }),
            relevant_cycle_algorithm: relevant_cycle_algorithm
                .map_or(defaults.relevant_cycle_algorithm, |algorithm| {
                    algorithm.to_rust()
                }),
        })
    }

    #[getter]
    fn simple_cycle_algorithm(&self) -> SimpleCycleEnumerationAlgorithm {
        SimpleCycleEnumerationAlgorithm::from_rust(self.0.simple_cycle_algorithm)
    }

    #[getter]
    fn relevant_cycle_algorithm(&self) -> RelevantCycleEnumerationAlgorithm {
        RelevantCycleEnumerationAlgorithm::from_rust(self.0.relevant_cycle_algorithm)
    }

    pub(crate) fn __repr__(&self) -> String {
        format!(
            "RingConfig(simple_cycle_algorithm={}, relevant_cycle_algorithm={})",
            self.simple_cycle_algorithm().repr(),
            self.relevant_cycle_algorithm().repr(),
        )
    }
}

impl RingConfig {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for ring-consuming configurations"
    )]
    pub(crate) fn from_rust(config: GraphIrRingConfig) -> Self {
        Self(config)
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for ring-consuming configurations"
    )]
    pub(crate) fn to_rust(self) -> GraphIrRingConfig {
        self.0
    }
}

/// Topological relation between two rings.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RingRelation {
    Identical,
    Disjoint,
    Spiro,
    Fused,
    Bridged,
    MultiSpiro,
    Noncontiguous,
}

#[pymethods]
impl RingRelation {
    fn __hash__(&self) -> u64 {
        hash_rust(self)
    }
}

impl RingRelation {
    fn from_rust(relation: GraphIrRingRelation) -> Self {
        match relation {
            GraphIrRingRelation::Identical => Self::Identical,
            GraphIrRingRelation::Disjoint => Self::Disjoint,
            GraphIrRingRelation::Spiro => Self::Spiro,
            GraphIrRingRelation::Fused => Self::Fused,
            GraphIrRingRelation::Bridged => Self::Bridged,
            GraphIrRingRelation::MultiSpiro => Self::MultiSpiro,
            GraphIrRingRelation::Noncontiguous => Self::Noncontiguous,
        }
    }

    fn to_rust(self) -> GraphIrRingRelation {
        match self {
            Self::Identical => GraphIrRingRelation::Identical,
            Self::Disjoint => GraphIrRingRelation::Disjoint,
            Self::Spiro => GraphIrRingRelation::Spiro,
            Self::Fused => GraphIrRingRelation::Fused,
            Self::Bridged => GraphIrRingRelation::Bridged,
            Self::MultiSpiro => GraphIrRingRelation::MultiSpiro,
            Self::Noncontiguous => GraphIrRingRelation::Noncontiguous,
        }
    }
}

/// An owned ring connection, or a read-only edge borrowed from a RingGraph.
/// Use copy() to obtain an independent writable connection.
#[pyclass]
pub struct RingConnection {
    storage: RingConnectionStorage,
}

enum RingConnectionStorage {
    Owned(GraphIrRingConnection),
    Graph {
        owner: Py<RingGraph>,
        position: usize,
    },
}

impl RingConnection {
    fn read<R>(&self, py: Python<'_>, f: impl FnOnce(&GraphIrRingConnection) -> R) -> PyResult<R> {
        match &self.storage {
            RingConnectionStorage::Owned(value) => Ok(f(value)),
            RingConnectionStorage::Graph { owner, position } => {
                let graph = owner.try_borrow(py)?;
                Ok(f(&graph.0.edges()[*position]))
            }
        }
    }

    fn to_rust_mut(&mut self) -> PyResult<&mut GraphIrRingConnection> {
        match &mut self.storage {
            RingConnectionStorage::Owned(value) => Ok(value),
            RingConnectionStorage::Graph { .. } => {
                Err(PyTypeError::new_err("RingGraph edges are read-only"))
            }
        }
    }
}

#[pymethods]
impl RingConnection {
    #[new]
    fn new(first: u32, second: u32, relation: RingRelation) -> Self {
        Self {
            storage: RingConnectionStorage::Owned(GraphIrRingConnection {
                first: GraphIrRingId(first),
                second: GraphIrRingId(second),
                relation: relation.to_rust(),
            }),
        }
    }

    #[getter]
    fn first(&self, py: Python<'_>) -> PyResult<u32> {
        self.read(py, |value| value.first.0)
    }

    #[setter]
    fn set_first(&mut self, first: u32) -> PyResult<()> {
        self.to_rust_mut()?.first = GraphIrRingId(first);
        Ok(())
    }

    #[getter]
    fn second(&self, py: Python<'_>) -> PyResult<u32> {
        self.read(py, |value| value.second.0)
    }

    #[setter]
    fn set_second(&mut self, second: u32) -> PyResult<()> {
        self.to_rust_mut()?.second = GraphIrRingId(second);
        Ok(())
    }

    #[getter]
    fn relation(&self, py: Python<'_>) -> PyResult<RingRelation> {
        self.read(py, |value| RingRelation::from_rust(value.relation))
    }

    #[setter]
    fn set_relation(&mut self, relation: RingRelation) -> PyResult<()> {
        self.to_rust_mut()?.relation = relation.to_rust();
        Ok(())
    }

    #[getter]
    fn readonly(&self) -> bool {
        matches!(self.storage, RingConnectionStorage::Graph { .. })
    }

    fn copy(&self, py: Python<'_>) -> PyResult<Self> {
        self.read(py, |value| Self {
            storage: RingConnectionStorage::Owned(*value),
        })
    }

    fn __eq__(&self, py: Python<'_>, other: &Self) -> PyResult<bool> {
        let value = self.read(py, |value| *value)?;
        other.read(py, |other| value == *other)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |value| {
            format!(
                "RingConnection(first={}, second={}, relation=RingRelation.{:?})",
                value.first.0, value.second.0, value.relation
            )
        })
    }
}

/// An independently owned graph of ring relations, produced by RingSet.graph().
#[pyclass(frozen)]
pub struct RingGraph(GraphIrRingGraph);

#[pymethods]
impl RingGraph {
    /// Borrow the graph's connections as a read-only sequence.
    #[getter]
    fn edges(slf: Py<Self>) -> RingConnections {
        RingConnections { owner: slf }
    }

    fn neighbors(&self, ring: u32) -> Vec<(u32, RingRelation)> {
        self.0
            .neighbors(GraphIrRingId(ring))
            .into_iter()
            .map(|(id, relation)| (id.0, RingRelation::from_rust(relation)))
            .collect()
    }

    fn relation(&self, first: u32, second: u32) -> RingRelation {
        RingRelation::from_rust(self.0.relation(GraphIrRingId(first), GraphIrRingId(second)))
    }

    fn copy(&self) -> Self {
        Self(self.0.clone())
    }

    fn __repr__(&self) -> String {
        format!("RingGraph(edges={})", self.0.edges().len())
    }
}

#[pyclass(sequence, frozen)]
struct RingConnections {
    owner: Py<RingGraph>,
}

#[pymethods]
impl RingConnections {
    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        Ok(self.owner.try_borrow(py)?.0.edges().len())
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<RingConnection> {
        let position = resolve_index(self.__len__(py)?, index, "ring connection")?;
        Ok(RingConnection {
            storage: RingConnectionStorage::Graph {
                owner: self.owner.clone_ref(py),
                position,
            },
        })
    }

    fn __iter__(slf: Py<Self>) -> RingConnectionIter {
        RingConnectionIter {
            edges: slf,
            position: 0,
        }
    }
}

#[pyclass]
struct RingConnectionIter {
    edges: Py<RingConnections>,
    position: usize,
}

#[pymethods]
impl RingConnectionIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<RingConnection>> {
        let edges = self.edges.try_borrow(py)?;
        if self.position == edges.__len__(py)? {
            return Ok(None);
        }
        let value = edges.__getitem__(py, self.position as isize)?;
        self.position += 1;
        Ok(Some(value))
    }
}

/// Owned ring topology, independent of the molecule used to perceive it.
/// Obtain a set by consuming RingViews with into_ring_set().
#[pyclass(frozen)]
pub struct RingSet(GraphIrRingSet);

#[pymethods]
impl RingSet {
    #[getter]
    fn kind(&self) -> RingSetKind {
        RingSetKind::from_rust(self.0.kind())
    }

    #[getter]
    fn max_ring_size(&self) -> usize {
        self.0.max_ring_size()
    }

    fn __len__(&self) -> usize {
        self.0.count()
    }

    fn __getitem__(slf: Py<Self>, py: Python<'_>, index: isize) -> PyResult<RingView> {
        let position = resolve_index(slf.try_borrow(py)?.0.count(), index, "ring")?;
        Ok(RingView::new(
            RingOwner::Set(slf),
            GraphIrRingId(position as u32),
        ))
    }

    fn get(slf: Py<Self>, py: Python<'_>, id: u32) -> PyResult<Option<RingView>> {
        let present = slf.try_borrow(py)?.0.get(GraphIrRingId(id)).is_some();
        Ok(present.then(|| RingView::new(RingOwner::Set(slf), GraphIrRingId(id))))
    }

    fn __contains__(&self, id: u32) -> bool {
        self.0.get(GraphIrRingId(id)).is_some()
    }

    fn ids(slf: Py<Self>) -> RingIdIter {
        RingIdIter {
            owner: RingOwner::Set(slf),
            position: 0,
        }
    }

    fn __iter__(slf: Py<Self>) -> RingViewIter {
        RingViewIter {
            source: RingIteration::Rings(RingOwner::Set(slf)),
            position: 0,
        }
    }

    fn shared_atoms(&self, first: u32, second: u32) -> Vec<u32> {
        self.0
            .shared_atoms(GraphIrRingId(first), GraphIrRingId(second))
            .into_iter()
            .map(|id| id.0)
            .collect()
    }

    fn shared_bonds(&self, first: u32, second: u32) -> Vec<u32> {
        self.0
            .shared_bonds(GraphIrRingId(first), GraphIrRingId(second))
            .into_iter()
            .map(|id| id.0)
            .collect()
    }

    fn relation(&self, first: u32, second: u32) -> RingRelation {
        RingRelation::from_rust(self.0.relation(GraphIrRingId(first), GraphIrRingId(second)))
    }

    fn are_spiro(&self, first: u32, second: u32) -> bool {
        self.0
            .are_spiro(GraphIrRingId(first), GraphIrRingId(second))
    }

    fn are_fused(&self, first: u32, second: u32) -> bool {
        self.0
            .are_fused(GraphIrRingId(first), GraphIrRingId(second))
    }

    fn are_bridged(&self, first: u32, second: u32) -> bool {
        self.0
            .are_bridged(GraphIrRingId(first), GraphIrRingId(second))
    }

    fn spiro_neighbors(&self, ring: u32) -> Vec<u32> {
        self.0
            .spiro_neighbors(GraphIrRingId(ring))
            .into_iter()
            .map(|id| id.0)
            .collect()
    }

    fn fused_neighbors(&self, ring: u32) -> Vec<u32> {
        self.0
            .fused_neighbors(GraphIrRingId(ring))
            .into_iter()
            .map(|id| id.0)
            .collect()
    }

    fn bridged_neighbors(&self, ring: u32) -> Vec<u32> {
        self.0
            .bridged_neighbors(GraphIrRingId(ring))
            .into_iter()
            .map(|id| id.0)
            .collect()
    }

    fn fused_component(&self, ring: u32) -> Vec<u32> {
        self.0
            .fused_component(GraphIrRingId(ring))
            .into_iter()
            .map(|id| id.0)
            .collect()
    }

    fn fused_components(&self) -> Vec<Vec<u32>> {
        self.0
            .fused_components()
            .into_iter()
            .map(|component| component.into_iter().map(|id| id.0).collect())
            .collect()
    }

    fn contains_atom(&self, atom: u32) -> bool {
        self.0.contains_atom(GraphIrAtomId(atom))
    }

    fn atom_smallest_ring_size(&self, atom: u32) -> Option<usize> {
        self.0.atom_smallest_ring_size(GraphIrAtomId(atom))
    }

    fn contains_bond(&self, bond: u32) -> bool {
        self.0.contains_bond(GraphIrBondId(bond))
    }

    fn bond_smallest_ring_size(&self, bond: u32) -> Option<usize> {
        self.0.bond_smallest_ring_size(GraphIrBondId(bond))
    }

    /// Return an independent copy of the ring-relation graph.
    fn graph(&self) -> RingGraph {
        RingGraph(self.0.graph())
    }

    fn copy(&self) -> Self {
        Self(self.0.clone())
    }

    fn __repr__(&self) -> String {
        format!(
            "RingSet(len={}, kind={:?}, max_ring_size={})",
            self.0.count(),
            self.0.kind(),
            self.0.max_ring_size()
        )
    }
}

/// Ring topology with a live molecule context.
/// Consuming into_ring_set invalidates this view and all dependent accessors.
#[pyclass]
pub struct RingViews {
    molecule: Py<Molecule>,
    counter: u64,
    storage: Option<RingStorage>,
}

enum RingStorage {
    Owned(GraphIrRingSet),
    Borrowed(Py<RingSet>),
}

impl RingViews {
    pub(crate) fn from_rust(
        py: Python<'_>,
        rings: GraphIrRingSet,
        molecule: Py<Molecule>,
    ) -> PyResult<Self> {
        let counter = molecule.try_borrow(py)?.view_counter()?;
        Ok(Self {
            molecule,
            counter,
            storage: Some(RingStorage::Owned(rings)),
        })
    }

    fn check_access(&self, accessor: &str) -> PyResult<()> {
        if self.storage.is_none() {
            return Err(InvalidatedViewError::new_err(format!(
                "{accessor} was invalidated by RingViews consumption"
            )));
        }
        Ok(())
    }

    fn read<R>(
        &self,
        py: Python<'_>,
        accessor: &'static str,
        f: impl FnOnce(&GraphIrRingSet, &GraphIrMolecule) -> PyResult<R>,
    ) -> PyResult<R> {
        let storage = self
            .storage
            .as_ref()
            .ok_or_else(|| ConsumedError::new_err("RingViews has been consumed"))?;
        let molecule = self.molecule.try_borrow(py)?;
        molecule.check_access(self.counter, accessor)?;
        let molecule = molecule.to_rust()?;
        match storage {
            RingStorage::Owned(rings) => f(rings, molecule),
            RingStorage::Borrowed(owner) => f(&owner.try_borrow(py)?.0, molecule),
        }
    }
}

#[pymethods]
impl RingViews {
    /// Borrow a ring set and molecule without copying or repeating perception.
    /// Topology and ids must agree for meaningful readings; this is not validated.
    #[staticmethod]
    fn from_ring_set(py: Python<'_>, rings: Py<RingSet>, molecule: Py<Molecule>) -> PyResult<Self> {
        let counter = molecule.try_borrow(py)?.view_counter()?;
        Ok(Self {
            molecule,
            counter,
            storage: Some(RingStorage::Borrowed(rings)),
        })
    }

    #[getter]
    fn kind(&self, py: Python<'_>) -> PyResult<RingSetKind> {
        self.read(py, "RingViews", |rings, _| {
            Ok(RingSetKind::from_rust(rings.kind()))
        })
    }

    #[getter]
    fn max_ring_size(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, "RingViews", |rings, _| Ok(rings.max_ring_size()))
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, "RingViews", |rings, _| Ok(rings.count()))
    }

    fn __getitem__(slf: Py<Self>, py: Python<'_>, index: isize) -> PyResult<RingView> {
        let position = resolve_index(slf.try_borrow(py)?.__len__(py)?, index, "ring")?;
        Ok(RingView::new(
            RingOwner::Views(slf),
            GraphIrRingId(position as u32),
        ))
    }

    fn get(slf: Py<Self>, py: Python<'_>, id: u32) -> PyResult<Option<RingView>> {
        let present = slf.try_borrow(py)?.read(py, "RingViews", |rings, _| {
            Ok(rings.get(GraphIrRingId(id)).is_some())
        })?;
        Ok(present.then(|| RingView::new(RingOwner::Views(slf), GraphIrRingId(id))))
    }

    fn __contains__(&self, py: Python<'_>, id: u32) -> PyResult<bool> {
        self.read(py, "RingViews", |rings, _| {
            Ok(rings.get(GraphIrRingId(id)).is_some())
        })
    }

    fn ids(slf: Py<Self>, py: Python<'_>) -> PyResult<RingIdIter> {
        slf.try_borrow(py)?.read(py, "RingViews", |_, _| Ok(()))?;
        Ok(RingIdIter {
            owner: RingOwner::Views(slf),
            position: 0,
        })
    }

    fn __iter__(slf: Py<Self>, py: Python<'_>) -> PyResult<RingViewIter> {
        slf.try_borrow(py)?.read(py, "RingViews", |_, _| Ok(()))?;
        Ok(RingViewIter {
            source: RingIteration::Rings(RingOwner::Views(slf)),
            position: 0,
        })
    }

    fn atom(slf: Py<Self>, py: Python<'_>, id: u32) -> PyResult<RingAtomView> {
        slf.try_borrow(py)?.read(py, "RingViews", |_, _| Ok(()))?;
        Ok(RingAtomView {
            owner: slf,
            id: GraphIrAtomId(id),
        })
    }

    fn bond(slf: Py<Self>, py: Python<'_>, id: u32) -> PyResult<RingBondView> {
        slf.try_borrow(py)?.read(py, "RingViews", |_, _| Ok(()))?;
        Ok(RingBondView {
            owner: slf,
            id: GraphIrBondId(id),
        })
    }

    /// Consume this view, moving owned ring data or cloning borrowed ring data.
    /// The returned set is independent of both the molecule and the source set.
    #[allow(
        clippy::wrong_self_convention,
        reason = "Consumes storage while retaining the Python owner in its consumed state"
    )]
    fn into_ring_set(&mut self, py: Python<'_>) -> PyResult<RingSet> {
        self.read(py, "RingViews", |_, _| Ok(()))?;
        match self
            .storage
            .take()
            .ok_or_else(|| ConsumedError::new_err("RingViews has been consumed"))?
        {
            RingStorage::Owned(rings) => Ok(RingSet(rings)),
            RingStorage::Borrowed(owner) => Ok(RingSet(owner.try_borrow(py)?.0.clone())),
        }
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, "RingViews", |rings, _| {
            Ok(format!("RingViews(len={})", rings.count()))
        })
    }
}

enum RingOwner {
    Set(Py<RingSet>),
    Views(Py<RingViews>),
}

impl RingOwner {
    fn clone_ref(&self, py: Python<'_>) -> Self {
        match self {
            Self::Set(owner) => Self::Set(owner.clone_ref(py)),
            Self::Views(owner) => Self::Views(owner.clone_ref(py)),
        }
    }

    fn read<R>(
        &self,
        py: Python<'_>,
        accessor: &'static str,
        f: impl FnOnce(&GraphIrRingSet) -> PyResult<R>,
    ) -> PyResult<R> {
        match self {
            Self::Set(owner) => f(&owner.try_borrow(py)?.0),
            Self::Views(owner) => {
                let owner = owner.try_borrow(py)?;
                owner.check_access(accessor)?;
                owner.read(py, accessor, |rings, _| f(rings))
            }
        }
    }
}

/// A single ring with read-only member sequences.
/// The writable id labels this view; changing it does not change its members.
#[pyclass]
pub struct RingView {
    owner: RingOwner,
    position: GraphIrRingId,
    id: GraphIrRingId,
}

impl RingView {
    fn new(owner: RingOwner, id: GraphIrRingId) -> Self {
        Self {
            owner,
            position: id,
            id,
        }
    }

    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(GraphIrRingView<'_>) -> PyResult<R>,
    ) -> PyResult<R> {
        self.owner.read(py, "RingView", |rings| {
            f(rings
                .get(self.position)
                .ok_or_else(|| PyIndexError::new_err("ring id out of range"))?)
        })
    }
}

#[pymethods]
impl RingView {
    #[getter]
    fn id(&self, py: Python<'_>) -> PyResult<u32> {
        self.read(py, |_| Ok(self.id.0))
    }

    #[setter]
    fn set_id(&mut self, py: Python<'_>, id: u32) -> PyResult<()> {
        self.read(py, |_| Ok(()))?;
        self.id = GraphIrRingId(id);
        Ok(())
    }

    #[getter]
    fn atoms(slf: Py<Self>, py: Python<'_>) -> PyResult<RingMembers> {
        slf.try_borrow(py)?.read(py, |_| Ok(()))?;
        Ok(RingMembers {
            owner: slf,
            kind: RingMemberKind::Atoms,
        })
    }

    #[getter]
    fn bonds(slf: Py<Self>, py: Python<'_>) -> PyResult<RingMembers> {
        slf.try_borrow(py)?.read(py, |_| Ok(()))?;
        Ok(RingMembers {
            owner: slf,
            kind: RingMemberKind::Bonds,
        })
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |ring| Ok(ring.len()))
    }

    fn is_empty(&self, py: Python<'_>) -> PyResult<bool> {
        self.read(py, |ring| Ok(ring.is_empty()))
    }

    fn copy(&self, py: Python<'_>) -> PyResult<Self> {
        self.read(py, |_| {
            Ok(Self {
                owner: self.owner.clone_ref(py),
                position: self.position,
                id: self.id,
            })
        })
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |ring| {
            Ok(format!("RingView(id={}, len={})", self.id.0, ring.len()))
        })
    }
}

enum RingMemberKind {
    Atoms,
    Bonds,
}

#[pyclass(sequence, frozen)]
struct RingMembers {
    owner: Py<RingView>,
    kind: RingMemberKind,
}

#[pymethods]
impl RingMembers {
    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.owner.try_borrow(py)?.__len__(py)
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<u32> {
        self.owner.try_borrow(py)?.read(py, |ring| {
            let position = resolve_index(ring.len(), index, "ring member")?;
            Ok(match self.kind {
                RingMemberKind::Atoms => ring.atoms()[position].0,
                RingMemberKind::Bonds => ring.bonds()[position].0,
            })
        })
    }

    fn __iter__(slf: Py<Self>, py: Python<'_>) -> PyResult<RingMemberIter> {
        slf.try_borrow(py)?.__len__(py)?;
        Ok(RingMemberIter {
            members: slf,
            position: 0,
        })
    }
}

#[pyclass]
struct RingMemberIter {
    members: Py<RingMembers>,
    position: usize,
}

#[pymethods]
impl RingMemberIter {
    fn __iter__<'py>(slf: PyRef<'py, Self>, py: Python<'py>) -> PyResult<PyRef<'py, Self>> {
        slf.members.try_borrow(py)?.__len__(py)?;
        Ok(slf)
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<u32>> {
        let members = self.members.try_borrow(py)?;
        if self.position == members.__len__(py)? {
            return Ok(None);
        }
        let value = members.__getitem__(py, self.position as isize)?;
        self.position += 1;
        Ok(Some(value))
    }
}

#[pyclass]
struct RingIdIter {
    owner: RingOwner,
    position: usize,
}

#[pymethods]
impl RingIdIter {
    fn __iter__<'py>(slf: PyRef<'py, Self>, py: Python<'py>) -> PyResult<PyRef<'py, Self>> {
        slf.owner.read(py, "RingIdIter", |_| Ok(()))?;
        Ok(slf)
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<u32>> {
        let value = self
            .owner
            .read(py, "RingIdIter", |rings| Ok(rings.ids().nth(self.position)))?;
        if value.is_some() {
            self.position += 1;
        }
        Ok(value.map(|id| id.0))
    }
}

enum RingIteration {
    Rings(RingOwner),
    Atom {
        owner: Py<RingViews>,
        id: GraphIrAtomId,
    },
    Bond {
        owner: Py<RingViews>,
        id: GraphIrBondId,
    },
}

#[pyclass]
struct RingViewIter {
    source: RingIteration,
    position: usize,
}

impl RingViewIter {
    fn read(&self, py: Python<'_>, position: Option<usize>) -> PyResult<Option<GraphIrRingId>> {
        match &self.source {
            RingIteration::Rings(owner) => owner.read(py, "RingViewIter", |rings| {
                Ok(position.and_then(|position| rings.ids().nth(position)))
            }),
            RingIteration::Atom { owner, id } => {
                let owner = owner.try_borrow(py)?;
                owner.check_access("RingViewIter")?;
                owner.read(py, "RingViewIter", |rings, molecule| {
                    Ok(position.and_then(|position| {
                        GraphIrRingViews::from_ring_set(rings, molecule)
                            .atom(*id)
                            .rings()
                            .nth(position)
                            .map(|ring| ring.id)
                    }))
                })
            }
            RingIteration::Bond { owner, id } => {
                let owner = owner.try_borrow(py)?;
                owner.check_access("RingViewIter")?;
                owner.read(py, "RingViewIter", |rings, molecule| {
                    Ok(position.and_then(|position| {
                        GraphIrRingViews::from_ring_set(rings, molecule)
                            .bond(*id)
                            .rings()
                            .nth(position)
                            .map(|ring| ring.id)
                    }))
                })
            }
        }
    }
}

#[pymethods]
impl RingViewIter {
    fn __iter__<'py>(slf: PyRef<'py, Self>, py: Python<'py>) -> PyResult<PyRef<'py, Self>> {
        slf.read(py, None)?;
        Ok(slf)
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<RingView>> {
        let Some(id) = self.read(py, Some(self.position))? else {
            return Ok(None);
        };
        self.position += 1;
        let owner = match &self.source {
            RingIteration::Rings(owner) => owner.clone_ref(py),
            RingIteration::Atom { owner, .. } | RingIteration::Bond { owner, .. } => {
                RingOwner::Views(owner.clone_ref(py))
            }
        };
        Ok(Some(RingView::new(owner, id)))
    }
}

/// Ring readings for one atom, retaining the ring and molecule owners.
#[pyclass]
pub struct RingAtomView {
    owner: Py<RingViews>,
    id: GraphIrAtomId,
}

impl RingAtomView {
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(GraphIrRingAtomView<'_>) -> PyResult<R>,
    ) -> PyResult<R> {
        let owner = self.owner.try_borrow(py)?;
        owner.check_access("RingAtomView")?;
        owner.read(py, "RingAtomView", |rings, molecule| {
            f(GraphIrRingViews::from_ring_set(rings, molecule).atom(self.id))
        })
    }
}

#[pymethods]
impl RingAtomView {
    fn is_in_ring(&self, py: Python<'_>) -> PyResult<bool> {
        self.read(py, |view| Ok(view.is_in_ring()))
    }

    fn rings(&self, py: Python<'_>) -> PyResult<RingViewIter> {
        self.read(py, |_| Ok(()))?;
        Ok(RingViewIter {
            source: RingIteration::Atom {
                owner: self.owner.clone_ref(py),
                id: self.id,
            },
            position: 0,
        })
    }

    fn ring_membership(&self, py: Python<'_>, scope: &RingScope) -> PyResult<NumForm> {
        self.read(py, |view| {
            NumForm::from_rust(py, &view.ring_membership(scope.to_rust()))
        })
    }

    #[getter]
    fn ring_count(&self, py: Python<'_>) -> PyResult<NumForm> {
        self.read(py, |view| NumForm::from_rust(py, &view.ring_count()))
    }

    fn ring_size_count(&self, py: Python<'_>, size: u8) -> PyResult<NumForm> {
        self.read(py, |view| {
            NumForm::from_rust(py, &view.ring_size_count(size))
        })
    }

    #[getter]
    fn smallest_ring_size(&self, py: Python<'_>) -> PyResult<Option<usize>> {
        self.read(py, |view| Ok(view.smallest_ring_size()))
    }

    #[getter]
    fn ring_degree(&self, py: Python<'_>) -> PyResult<NumForm> {
        self.read(py, |view| NumForm::from_rust(py, &view.ring_degree()))
    }

    /// Sum of incident ring-bond orders: None denotes bottom, Undetermined denotes top.
    #[getter]
    fn ring_valence(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |view| {
            view.ring_valence()
                .map(|value| NumForm::from_rust(py, &value))
                .transpose()
        })
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |_| Ok(format!("RingAtomView(atom={})", self.id.0)))
    }
}

/// Ring readings for one bond, retaining the ring and molecule owners.
#[pyclass]
pub struct RingBondView {
    owner: Py<RingViews>,
    id: GraphIrBondId,
}

impl RingBondView {
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(GraphIrRingBondView<'_>) -> PyResult<R>,
    ) -> PyResult<R> {
        let owner = self.owner.try_borrow(py)?;
        owner.check_access("RingBondView")?;
        owner.read(py, "RingBondView", |rings, molecule| {
            f(GraphIrRingViews::from_ring_set(rings, molecule).bond(self.id))
        })
    }
}

#[pymethods]
impl RingBondView {
    fn is_in_ring(&self, py: Python<'_>) -> PyResult<bool> {
        self.read(py, |view| Ok(view.is_in_ring()))
    }

    fn rings(&self, py: Python<'_>) -> PyResult<RingViewIter> {
        self.read(py, |_| Ok(()))?;
        Ok(RingViewIter {
            source: RingIteration::Bond {
                owner: self.owner.clone_ref(py),
                id: self.id,
            },
            position: 0,
        })
    }

    fn ring_membership(&self, py: Python<'_>, scope: &RingScope) -> PyResult<NumForm> {
        self.read(py, |view| {
            NumForm::from_rust(py, &view.ring_membership(scope.to_rust()))
        })
    }

    #[getter]
    fn ring_count(&self, py: Python<'_>) -> PyResult<NumForm> {
        self.read(py, |view| NumForm::from_rust(py, &view.ring_count()))
    }

    fn ring_size_count(&self, py: Python<'_>, size: u8) -> PyResult<NumForm> {
        self.read(py, |view| {
            NumForm::from_rust(py, &view.ring_size_count(size))
        })
    }

    #[getter]
    fn smallest_ring_size(&self, py: Python<'_>) -> PyResult<Option<usize>> {
        self.read(py, |view| Ok(view.smallest_ring_size()))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |_| Ok(format!("RingBondView(bond={})", self.id.0)))
    }
}

fn resolve_index(len: usize, index: isize, name: &str) -> PyResult<usize> {
    let position = if index < 0 {
        len as isize + index
    } else {
        index
    };
    if position < 0 || position as usize >= len {
        return Err(PyIndexError::new_err(format!("{name} index out of range")));
    }
    Ok(position as usize)
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_chem::element::Element as ChemElement;
    use umol_graph_core::{
        RelevantCycleEnumerationAlgorithm as GraphCoreRelevantCycleEnumerationAlgorithm,
        SimpleCycleEnumerationAlgorithm as GraphCoreSimpleCycleEnumerationAlgorithm,
    };
    use umol_graph_ir::ir::{
        AtomForm as GraphIrAtomForm, BondForm as GraphIrBondForm,
        MoleculeEntries as GraphIrMoleculeEntries,
    };

    use super::*;

    #[fixture]
    fn triangle() -> GraphIrMolecule {
        GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C); 3],
            bonds: vec![
                (
                    GraphIrAtomId(0),
                    GraphIrAtomId(1),
                    GraphIrBondForm::from_order(1),
                ),
                (
                    GraphIrAtomId(1),
                    GraphIrAtomId(2),
                    GraphIrBondForm::from_order(1),
                ),
                (
                    GraphIrAtomId(2),
                    GraphIrAtomId(0),
                    GraphIrBondForm::from_order(1),
                ),
            ],
            ..Default::default()
        })
    }

    #[rstest]
    #[case::default(
        None,
        None,
        GraphIrRingConfig::default(),
        "RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara())",
    )]
    #[case::explicit(
        Some(SimpleCycleEnumerationAlgorithm::ReadTarjan()),
        Some(RelevantCycleEnumerationAlgorithm::Vismara()),
        GraphIrRingConfig {
            simple_cycle_algorithm: GraphCoreSimpleCycleEnumerationAlgorithm::ReadTarjan,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        },
        "RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara())",
    )]
    fn test_ring_config_new(
        #[case] simple_cycle_algorithm: Option<SimpleCycleEnumerationAlgorithm>,
        #[case] relevant_cycle_algorithm: Option<RelevantCycleEnumerationAlgorithm>,
        #[case] expected: GraphIrRingConfig,
        #[case] expected_repr: &str,
    ) {
        let config = RingConfig::new(simple_cycle_algorithm, relevant_cycle_algorithm);

        assert_eq!(config.to_rust(), expected);
        assert_eq!(config.__repr__(), expected_repr);
    }

    #[rstest]
    #[case::default(GraphIrRingConfig::default())]
    fn test_ring_config_from_rust(#[case] config: GraphIrRingConfig) {
        assert_eq!(RingConfig::from_rust(config).to_rust(), config);
    }

    #[rstest]
    #[case::owned(false)]
    #[case::borrowed(true)]
    fn test_ring_views_into_ring_set(triangle: GraphIrMolecule, #[case] borrowed: bool) {
        let rings = triangle
            .rings(GraphIrRingModel::default(), GraphIrRingConfig::default())
            .into_ring_set();
        let atoms = rings.get(GraphIrRingId(0)).unwrap().atoms().as_ptr();
        let bonds = rings.get(GraphIrRingId(0)).unwrap().bonds().as_ptr();
        Python::attach(|py| {
            let molecule = Py::new(py, Molecule::from_rust(triangle)).unwrap();
            let mut views = if borrowed {
                let source = Py::new(py, RingSet(rings)).unwrap();
                RingViews::from_ring_set(py, source, molecule).unwrap()
            } else {
                RingViews::from_rust(py, rings, molecule).unwrap()
            };
            views
                .read(py, "RingViews", |rings, _| {
                    let ring = rings.get(GraphIrRingId(0)).unwrap();
                    assert_eq!(ring.atoms().as_ptr(), atoms);
                    assert_eq!(ring.bonds().as_ptr(), bonds);
                    Ok(())
                })
                .unwrap();
            let extracted = views.into_ring_set(py).unwrap();
            let ring = extracted.0.get(GraphIrRingId(0)).unwrap();
            assert_eq!(ring.atoms().as_ptr() == atoms, !borrowed);
            assert_eq!(ring.bonds().as_ptr() == bonds, !borrowed);
            assert_eq!(
                ring.atoms(),
                &[GraphIrAtomId(0), GraphIrAtomId(1), GraphIrAtomId(2)]
            );
            assert_eq!(
                ring.bonds(),
                &[GraphIrBondId(0), GraphIrBondId(1), GraphIrBondId(2)]
            );
            assert!(views
                .__len__(py)
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
        });
    }
}
