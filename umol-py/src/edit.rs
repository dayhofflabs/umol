//! Python edit values and symbolic creation handles.

use std::collections::HashMap;
use std::str::FromStr;

use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use umol_graph_ir::dsl::EditsDsl as GraphIrEditsDsl;
use umol_graph_ir::ir::{
    AddAromaticSystem as GraphIrAddAromaticSystem, AddBond as GraphIrAddBond,
    AddMulticenterBond as GraphIrAddMulticenterBond, AddStereoAtom as GraphIrAddStereoAtom,
    AddStereoBond as GraphIrAddStereoBond, AromaticSystemHandle as GraphIrAromaticSystemHandle,
    AromaticSystemId as GraphIrAromaticSystemId, AtomHandle as GraphIrAtomHandle,
    AtomId as GraphIrAtomId, BondHandle as GraphIrBondHandle, BondId as GraphIrBondId,
    ConstraintEdit as GraphIrConstraintEdit, DativeBondHandle as GraphIrDativeBondHandle,
    DativeBondId as GraphIrDativeBondId, Edit as GraphIrEdit, Edits as GraphIrEdits,
    Entity as GraphIrEntity, EntityHandle as GraphIrEntityHandle, FromIr, IntoIr,
    MulticenterBondHandle as GraphIrMulticenterBondHandle,
    MulticenterBondId as GraphIrMulticenterBondId,
    NoncovalentBondHandle as GraphIrNoncovalentBondHandle,
    NoncovalentBondId as GraphIrNoncovalentBondId, StereoAtomHandle as GraphIrStereoAtomHandle,
    StereoAtomId as GraphIrStereoAtomId, StereoBondHandle as GraphIrStereoBondHandle,
    StereoBondId as GraphIrStereoBondId,
};

use crate::aromatic::{AromaticSystemForm, AromaticSystemUpdate};
use crate::atom::{AtomForm, AtomUpdate, ElementForm, IsotopeMassForm};
use crate::bond::{BondForm, BondUpdate};
use crate::constraint::aromatic::{AromaticSystemConstraintForm, AromaticSystemConstraintKey};
use crate::constraint::atom::{AtomConstraintForm, AtomConstraintKey};
use crate::constraint::bond::{BondConstraintForm, BondConstraintKey};
use crate::constraint::dative::{DativeBondConstraintForm, DativeBondConstraintKey};
use crate::constraint::molecule::Constraint;
use crate::constraint::multicenter::{MulticenterBondConstraintForm, MulticenterBondConstraintKey};
use crate::constraint::noncovalent::{NoncovalentBondConstraintForm, NoncovalentBondConstraintKey};
use crate::constraint::stereo::{
    StereoAtomConstraintForm, StereoAtomConstraintKey, StereoBondConstraintForm,
    StereoBondConstraintKey,
};
use crate::convert::into_py_variant;
use crate::dative::{DativeBondForm, DativeBondUpdate};
use crate::defaults::MoleculeDefaults;
use crate::delta::{
    AromaticSystemFieldChange, AtomFieldChange, BondFieldChange, DativeBondFieldChange,
    MulticenterBondFieldChange, NoncovalentBondFieldChange, StereoAtomFieldChange,
    StereoBondFieldChange,
};
use crate::electrons::ElectronCountsForm;
use crate::entity::Readonly;
use crate::error::{parse_error, ConsumedError, InvalidatedViewError};
use crate::metadata::Entity;
use crate::multicenter::{MulticenterBondForm, MulticenterBondUpdate};
use crate::noncovalent::{NoncovalentBondForm, NoncovalentBondKindForm, NoncovalentBondUpdate};
use crate::num::NumForm;
use crate::spin::UnpairedElectronsForm;
use crate::stereo::{
    StereoAtomForm, StereoAtomUpdate, StereoBondForm, StereoBondUpdate, StereoConfigurationForm,
    StereoKind, StereoLigandKind,
};

/// A same-kind creation ordinal in an edit sequence.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct New {
    index: usize,
}

#[pymethods]
impl New {
    #[new]
    fn from_index(index: usize) -> Self {
        Self { index }
    }

    /// Zero-based creation ordinal within the surrounding entity kind.
    #[getter]
    fn index(&self) -> usize {
        self.index
    }

    fn __repr__(&self) -> String {
        format!("New({})", self.index)
    }
}

impl New {
    fn from_rust(handle: GraphIrEntityHandle) -> Self {
        let index = match handle {
            GraphIrEntityHandle::Atom(GraphIrAtomHandle::New(index))
            | GraphIrEntityHandle::Bond(GraphIrBondHandle::New(index))
            | GraphIrEntityHandle::DativeBond(GraphIrDativeBondHandle::New(index))
            | GraphIrEntityHandle::AromaticSystem(GraphIrAromaticSystemHandle::New(index))
            | GraphIrEntityHandle::MulticenterBond(GraphIrMulticenterBondHandle::New(index))
            | GraphIrEntityHandle::NoncovalentBond(GraphIrNoncovalentBondHandle::New(index))
            | GraphIrEntityHandle::StereoAtom(GraphIrStereoAtomHandle::New(index))
            | GraphIrEntityHandle::StereoBond(GraphIrStereoBondHandle::New(index)) => index,
            _ => unreachable!("Edits addition methods always return creation handles"),
        };
        Self { index }
    }
}

/// A Python host id or same-kind creation ordinal. The argument position supplies the entity kind.
#[derive(Clone, FromPyObject)]
#[allow(
    dead_code,
    reason = "Python-to-Rust handle input for the Edit and Edits bindings"
)]
pub(crate) enum HandleLike {
    New(New),
    Id(u32),
}

#[allow(
    dead_code,
    reason = "Python-to-Rust handle conversions for the Edit and Edits bindings"
)]
impl HandleLike {
    fn from_atom_handle(handle: &GraphIrAtomHandle) -> Self {
        match handle {
            GraphIrAtomHandle::Id(id) => Self::Id(id.0),
            GraphIrAtomHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_bond_handle(handle: &GraphIrBondHandle) -> Self {
        match handle {
            GraphIrBondHandle::Id(id) => Self::Id(id.0),
            GraphIrBondHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_dative_bond_handle(handle: &GraphIrDativeBondHandle) -> Self {
        match handle {
            GraphIrDativeBondHandle::Id(id) => Self::Id(id.0),
            GraphIrDativeBondHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_aromatic_system_handle(handle: &GraphIrAromaticSystemHandle) -> Self {
        match handle {
            GraphIrAromaticSystemHandle::Id(id) => Self::Id(id.0),
            GraphIrAromaticSystemHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_multicenter_bond_handle(handle: &GraphIrMulticenterBondHandle) -> Self {
        match handle {
            GraphIrMulticenterBondHandle::Id(id) => Self::Id(id.0),
            GraphIrMulticenterBondHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_noncovalent_bond_handle(handle: &GraphIrNoncovalentBondHandle) -> Self {
        match handle {
            GraphIrNoncovalentBondHandle::Id(id) => Self::Id(id.0),
            GraphIrNoncovalentBondHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_stereo_atom_handle(handle: &GraphIrStereoAtomHandle) -> Self {
        match handle {
            GraphIrStereoAtomHandle::Id(id) => Self::Id(id.0),
            GraphIrStereoAtomHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_stereo_bond_handle(handle: &GraphIrStereoBondHandle) -> Self {
        match handle {
            GraphIrStereoBondHandle::Id(id) => Self::Id(id.0),
            GraphIrStereoBondHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    pub(crate) fn to_atom_handle(&self) -> GraphIrAtomHandle {
        match self {
            Self::Id(index) => GraphIrAtomHandle::Id(GraphIrAtomId(*index)),
            Self::New(new) => GraphIrAtomHandle::New(new.index),
        }
    }

    pub(crate) fn to_bond_handle(&self) -> GraphIrBondHandle {
        match self {
            Self::Id(index) => GraphIrBondHandle::Id(GraphIrBondId(*index)),
            Self::New(new) => GraphIrBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_dative_bond_handle(&self) -> GraphIrDativeBondHandle {
        match self {
            Self::Id(index) => GraphIrDativeBondHandle::Id(GraphIrDativeBondId(*index)),
            Self::New(new) => GraphIrDativeBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_aromatic_system_handle(&self) -> GraphIrAromaticSystemHandle {
        match self {
            Self::Id(index) => GraphIrAromaticSystemHandle::Id(GraphIrAromaticSystemId(*index)),
            Self::New(new) => GraphIrAromaticSystemHandle::New(new.index),
        }
    }

    pub(crate) fn to_multicenter_bond_handle(&self) -> GraphIrMulticenterBondHandle {
        match self {
            Self::Id(index) => GraphIrMulticenterBondHandle::Id(GraphIrMulticenterBondId(*index)),
            Self::New(new) => GraphIrMulticenterBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_noncovalent_bond_handle(&self) -> GraphIrNoncovalentBondHandle {
        match self {
            Self::Id(index) => GraphIrNoncovalentBondHandle::Id(GraphIrNoncovalentBondId(*index)),
            Self::New(new) => GraphIrNoncovalentBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_stereo_atom_handle(&self) -> GraphIrStereoAtomHandle {
        match self {
            Self::Id(index) => GraphIrStereoAtomHandle::Id(GraphIrStereoAtomId(*index)),
            Self::New(new) => GraphIrStereoAtomHandle::New(new.index),
        }
    }

    pub(crate) fn to_stereo_bond_handle(&self) -> GraphIrStereoBondHandle {
        match self {
            Self::Id(index) => GraphIrStereoBondHandle::Id(GraphIrStereoBondId(*index)),
            Self::New(new) => GraphIrStereoBondHandle::New(new.index),
        }
    }

    fn to_entity_handle(&self, entity: GraphIrEntity) -> GraphIrEntityHandle {
        match entity {
            GraphIrEntity::Atom(_) => GraphIrEntityHandle::Atom(self.to_atom_handle()),
            GraphIrEntity::Bond(_) => GraphIrEntityHandle::Bond(self.to_bond_handle()),
            GraphIrEntity::DativeBond(_) => {
                GraphIrEntityHandle::DativeBond(self.to_dative_bond_handle())
            }
            GraphIrEntity::AromaticSystem(_) => {
                GraphIrEntityHandle::AromaticSystem(self.to_aromatic_system_handle())
            }
            GraphIrEntity::MulticenterBond(_) => {
                GraphIrEntityHandle::MulticenterBond(self.to_multicenter_bond_handle())
            }
            GraphIrEntity::NoncovalentBond(_) => {
                GraphIrEntityHandle::NoncovalentBond(self.to_noncovalent_bond_handle())
            }
            GraphIrEntity::StereoAtom(_) => {
                GraphIrEntityHandle::StereoAtom(self.to_stereo_atom_handle())
            }
            GraphIrEntity::StereoBond(_) => {
                GraphIrEntityHandle::StereoBond(self.to_stereo_bond_handle())
            }
        }
    }
}

impl<'py> IntoPyObject<'py> for &HandleLike {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        match self {
            HandleLike::Id(index) => Ok(index.into_pyobject(py)?.into_any()),
            HandleLike::New(new) => Ok(Bound::new(py, *new)?.into_any()),
        }
    }
}

impl<'py> IntoPyObject<'py> for HandleLike {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        (&self).into_pyobject(py)
    }
}

/// A handle-aware molecule constraint used by standalone edit entries.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConstraintEdit(GraphIrConstraintEdit);

#[pymethods]
impl ConstraintEdit {
    #[new]
    #[pyo3(signature = (constraint, *, handles=None))]
    fn new(
        py: Python<'_>,
        constraint: Py<Constraint>,
        handles: Option<HashMap<Entity, HandleLike>>,
    ) -> PyResult<Self> {
        let constraint = constraint.bind(py).borrow().to_rust(py)?;
        let Some(handles) = handles else {
            return Ok(Self(GraphIrConstraintEdit::from(constraint)));
        };
        GraphIrConstraintEdit::new(constraint, |entity| {
            handles
                .get(&Entity::from_rust(entity))
                .map(|handle| handle.to_entity_handle(entity))
        })
        .map(Self)
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn __repr__(&self) -> &'static str {
        "ConstraintEdit(...)"
    }
}

impl ConstraintEdit {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API used by Edit snapshots"
    )]
    pub(crate) fn from_rust(constraint: &GraphIrConstraintEdit) -> Self {
        Self(constraint.clone())
    }

    pub(crate) fn to_rust(&self) -> &GraphIrConstraintEdit {
        &self.0
    }
}

type ReadonlyBondAddition = ((HandleLike, HandleLike), Readonly<BondForm>);
type ReadonlyDativeBondRemoval = (
    HandleLike,
    Vec<HandleLike>,
    HandleLike,
    Readonly<DativeBondForm>,
);
type ReadonlyAromaticSystemRemoval = (HandleLike, Vec<HandleLike>, Readonly<AromaticSystemForm>);
type ReadonlyMulticenterBondRemoval = (HandleLike, Vec<HandleLike>, Readonly<MulticenterBondForm>);
type ReadonlyNoncovalentBondRemoval = (
    HandleLike,
    (HandleLike, HandleLike),
    Readonly<NoncovalentBondForm>,
);
type StereoLigandInput = (HandleLike, StereoLigandKind);
type ReadonlyStereoAtomRemoval = (
    HandleLike,
    HandleLike,
    Vec<StereoLigandInput>,
    Readonly<StereoAtomForm>,
);
type ReadonlyStereoBondRemoval = (
    HandleLike,
    HandleLike,
    Vec<StereoLigandInput>,
    Readonly<StereoBondForm>,
);

type BondAddition = ((HandleLike, HandleLike), Py<BondForm>);
type DativeBondAddition = (Vec<HandleLike>, HandleLike, Py<DativeBondForm>);
type AromaticSystemAddition = (Vec<HandleLike>, Py<AromaticSystemForm>);
type MulticenterBondAddition = (Vec<HandleLike>, Py<MulticenterBondForm>);
type NoncovalentBondAddition = ((HandleLike, HandleLike), Py<NoncovalentBondForm>);
type StereoAtomAddition = (HandleLike, Vec<StereoLigandInput>, Py<StereoAtomForm>);
type StereoBondAddition = (HandleLike, Vec<StereoLigandInput>, Py<StereoBondForm>);
type DativeBondRemoval = (HandleLike, Vec<HandleLike>, HandleLike, Py<DativeBondForm>);
type AromaticSystemRemoval = (HandleLike, Vec<HandleLike>, Py<AromaticSystemForm>);
type MulticenterBondRemoval = (HandleLike, Vec<HandleLike>, Py<MulticenterBondForm>);
type NoncovalentBondRemoval = (
    HandleLike,
    (HandleLike, HandleLike),
    Py<NoncovalentBondForm>,
);
type StereoAtomRemoval = (
    HandleLike,
    HandleLike,
    Vec<StereoLigandInput>,
    Py<StereoAtomForm>,
);
type StereoBondRemoval = (
    HandleLike,
    HandleLike,
    Vec<StereoLigandInput>,
    Py<StereoBondForm>,
);

struct ReadonlyStereoAtomRemovals(Vec<ReadonlyStereoAtomRemoval>);

impl FromPyObject<'_, '_> for ReadonlyStereoAtomRemovals {
    type Error = PyErr;

    fn extract(object: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        object.extract().map(Self)
    }
}

impl<'py> IntoPyObject<'py> for &ReadonlyStereoAtomRemovals {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        let entries = self
            .0
            .iter()
            .map(|(id, site, ligands, attributes)| {
                Ok((
                    id.clone(),
                    site.clone(),
                    ligands.clone(),
                    attributes.into_pyobject(py)?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        entries.into_pyobject(py).map(Bound::into_any)
    }
}

struct ReadonlyStereoBondRemovals(Vec<ReadonlyStereoBondRemoval>);

impl FromPyObject<'_, '_> for ReadonlyStereoBondRemovals {
    type Error = PyErr;

    fn extract(object: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        object.extract().map(Self)
    }
}

impl<'py> IntoPyObject<'py> for &ReadonlyStereoBondRemovals {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        let entries = self
            .0
            .iter()
            .map(|(id, site, ligands, attributes)| {
                Ok((
                    id.clone(),
                    site.clone(),
                    ligands.clone(),
                    attributes.into_pyobject(py)?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        entries.into_pyobject(py).map(Bound::into_any)
    }
}

type AtomConstraintChange = (
    HandleLike,
    Option<Py<AtomConstraintForm>>,
    Option<Py<AtomConstraintForm>>,
);

type BondConstraintChange = (
    HandleLike,
    Option<Py<BondConstraintForm>>,
    Option<Py<BondConstraintForm>>,
);

type DativeBondConstraintChange = (
    HandleLike,
    Option<Py<DativeBondConstraintForm>>,
    Option<Py<DativeBondConstraintForm>>,
);

type AromaticSystemConstraintChange = (
    HandleLike,
    Option<Py<AromaticSystemConstraintForm>>,
    Option<Py<AromaticSystemConstraintForm>>,
);

type MulticenterBondConstraintChange = (
    HandleLike,
    Option<Py<MulticenterBondConstraintForm>>,
    Option<Py<MulticenterBondConstraintForm>>,
);

type NoncovalentBondConstraintChange = (
    HandleLike,
    Option<Py<NoncovalentBondConstraintForm>>,
    Option<Py<NoncovalentBondConstraintForm>>,
);

type StereoAtomConstraintChange = (
    HandleLike,
    Option<Py<StereoAtomConstraintForm>>,
    Option<Py<StereoAtomConstraintForm>>,
);

type StereoBondConstraintChange = (
    HandleLike,
    Option<Py<StereoBondConstraintForm>>,
    Option<Py<StereoBondConstraintForm>>,
);

struct TupleItems<T>(Vec<T>);

impl<'a, 'py, T> FromPyObject<'a, 'py> for TupleItems<T>
where
    Vec<T>: FromPyObject<'a, 'py, Error = PyErr>,
{
    type Error = PyErr;

    fn extract(object: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        object.extract().map(Self)
    }
}

impl<'a, 'py, T> IntoPyObject<'py> for &'a TupleItems<T>
where
    &'a T: IntoPyObject<'py>,
{
    type Target = PyTuple;
    type Output = Bound<'py, PyTuple>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        PyTuple::new(py, &self.0)
    }
}

/// One immutable host-specific molecule edit.
#[pyclass(frozen)]
#[allow(
    private_interfaces,
    reason = "PyO3 exposes Python field values through private coercion adapters"
)]
pub enum Edit {
    AddAtoms {
        atoms: Vec<Readonly<AtomForm>>,
    },
    AddBonds {
        bonds: Vec<ReadonlyBondAddition>,
    },
    RemoveTopology {
        atoms: Vec<HandleLike>,
        bonds: Vec<HandleLike>,
    },
    ModifyAtomElements {
        changes: TupleItems<(HandleLike, Py<ElementForm>, Py<ElementForm>)>,
    },
    ModifyAtomIsotopeMasses {
        changes: TupleItems<(HandleLike, Py<IsotopeMassForm>, Py<IsotopeMassForm>)>,
    },
    ModifyAtomCharges {
        changes: TupleItems<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    },
    ModifyAtomImplicitHydrogens {
        changes: TupleItems<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    },
    ModifyAtomLonePairs {
        changes: TupleItems<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    },
    ModifyAtomUnpairedElectrons {
        changes: TupleItems<(
            HandleLike,
            Py<UnpairedElectronsForm>,
            Py<UnpairedElectronsForm>,
        )>,
    },
    ModifyAtomConstraints {
        key: Py<AtomConstraintKey>,
        changes: TupleItems<AtomConstraintChange>,
    },
    ModifyAtomField {
        id: HandleLike,
        change: Py<AtomFieldChange>,
    },
    ModifyBondOrders {
        changes: TupleItems<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    },
    ModifyBondCharges {
        changes: TupleItems<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    },
    ModifyBondUnpairedElectrons {
        changes: TupleItems<(
            HandleLike,
            Py<UnpairedElectronsForm>,
            Py<UnpairedElectronsForm>,
        )>,
    },
    ModifyBondField {
        id: HandleLike,
        change: Py<BondFieldChange>,
    },
    AddDativeBond {
        donors: Vec<HandleLike>,
        acceptor: HandleLike,
        attributes: Readonly<DativeBondForm>,
    },
    RemoveDativeBonds {
        removes: Vec<ReadonlyDativeBondRemoval>,
    },
    ReplaceDativeBondDonors {
        id: HandleLike,
        old: Vec<HandleLike>,
        new: Vec<HandleLike>,
    },
    ReplaceDativeBondAcceptor {
        id: HandleLike,
        old: HandleLike,
        new: HandleLike,
    },
    ModifyDativeBondOrders {
        changes: TupleItems<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    },
    ModifyDativeBondField {
        id: HandleLike,
        change: Py<DativeBondFieldChange>,
    },
    AddAromaticSystem {
        atoms: Vec<HandleLike>,
        attributes: Readonly<AromaticSystemForm>,
    },
    RemoveAromaticSystems {
        removes: Vec<ReadonlyAromaticSystemRemoval>,
    },
    ReplaceAromaticSystemAtoms {
        id: HandleLike,
        old: Vec<HandleLike>,
        new: Vec<HandleLike>,
    },
    ModifyAromaticSystemElectrons {
        changes: TupleItems<(HandleLike, Py<ElectronCountsForm>, Py<ElectronCountsForm>)>,
    },
    ModifyAromaticSystemCharges {
        changes: TupleItems<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    },
    ModifyAromaticSystemUnpairedElectrons {
        changes: TupleItems<(
            HandleLike,
            Py<UnpairedElectronsForm>,
            Py<UnpairedElectronsForm>,
        )>,
    },
    ModifyAromaticSystemField {
        id: HandleLike,
        change: Py<AromaticSystemFieldChange>,
    },
    AddMulticenterBond {
        atoms: Vec<HandleLike>,
        attributes: Readonly<MulticenterBondForm>,
    },
    RemoveMulticenterBonds {
        removes: Vec<ReadonlyMulticenterBondRemoval>,
    },
    ReplaceMulticenterBondAtoms {
        id: HandleLike,
        old: Vec<HandleLike>,
        new: Vec<HandleLike>,
    },
    ModifyMulticenterBondElectrons {
        changes: TupleItems<(HandleLike, Py<ElectronCountsForm>, Py<ElectronCountsForm>)>,
    },
    ModifyMulticenterBondCharges {
        changes: TupleItems<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    },
    ModifyMulticenterBondUnpairedElectrons {
        changes: TupleItems<(
            HandleLike,
            Py<UnpairedElectronsForm>,
            Py<UnpairedElectronsForm>,
        )>,
    },
    ModifyMulticenterBondField {
        id: HandleLike,
        change: Py<MulticenterBondFieldChange>,
    },
    AddNoncovalentBond {
        atoms: (HandleLike, HandleLike),
        attributes: Readonly<NoncovalentBondForm>,
    },
    RemoveNoncovalentBonds {
        removes: Vec<ReadonlyNoncovalentBondRemoval>,
    },
    ReplaceNoncovalentBondAtoms {
        id: HandleLike,
        old: (HandleLike, HandleLike),
        new: (HandleLike, HandleLike),
    },
    ModifyNoncovalentBondKinds {
        changes: TupleItems<(
            HandleLike,
            Py<NoncovalentBondKindForm>,
            Py<NoncovalentBondKindForm>,
        )>,
    },
    ModifyNoncovalentBondField {
        id: HandleLike,
        change: Py<NoncovalentBondFieldChange>,
    },
    AddStereoAtom {
        site: HandleLike,
        ligands: Vec<StereoLigandInput>,
        attributes: Readonly<StereoAtomForm>,
    },
    RemoveStereoAtoms {
        removes: ReadonlyStereoAtomRemovals,
    },
    ReplaceStereoAtomSite {
        id: HandleLike,
        old: HandleLike,
        new: HandleLike,
    },
    ReplaceStereoAtomLigands {
        id: HandleLike,
        old: Vec<StereoLigandInput>,
        new: Vec<StereoLigandInput>,
    },
    ModifyStereoAtomConfigurations {
        changes: TupleItems<(
            HandleLike,
            Py<StereoConfigurationForm>,
            Py<StereoConfigurationForm>,
        )>,
    },
    ModifyStereoAtomField {
        id: HandleLike,
        change: Py<StereoAtomFieldChange>,
    },
    AddStereoBond {
        site: HandleLike,
        ligands: Vec<StereoLigandInput>,
        attributes: Readonly<StereoBondForm>,
    },
    RemoveStereoBonds {
        removes: ReadonlyStereoBondRemovals,
    },
    ReplaceStereoBondSite {
        id: HandleLike,
        old: HandleLike,
        new: HandleLike,
    },
    ReplaceStereoBondLigands {
        id: HandleLike,
        old: Vec<StereoLigandInput>,
        new: Vec<StereoLigandInput>,
    },
    ModifyStereoBondConfigurations {
        changes: TupleItems<(
            HandleLike,
            Py<StereoConfigurationForm>,
            Py<StereoConfigurationForm>,
        )>,
    },
    ModifyStereoBondField {
        id: HandleLike,
        change: Py<StereoBondFieldChange>,
    },
    ModifyAtomConstraint {
        id: HandleLike,
        old: Option<Py<AtomConstraintForm>>,
        new: Option<Py<AtomConstraintForm>>,
    },
    ModifyBondConstraints {
        key: Py<BondConstraintKey>,
        changes: TupleItems<BondConstraintChange>,
    },
    ModifyBondConstraint {
        id: HandleLike,
        old: Option<Py<BondConstraintForm>>,
        new: Option<Py<BondConstraintForm>>,
    },
    ModifyDativeBondConstraints {
        key: Py<DativeBondConstraintKey>,
        changes: TupleItems<DativeBondConstraintChange>,
    },
    ModifyDativeBondConstraint {
        id: HandleLike,
        old: Option<Py<DativeBondConstraintForm>>,
        new: Option<Py<DativeBondConstraintForm>>,
    },
    ModifyAromaticSystemConstraints {
        key: Py<AromaticSystemConstraintKey>,
        changes: TupleItems<AromaticSystemConstraintChange>,
    },
    ModifyAromaticSystemConstraint {
        id: HandleLike,
        old: Option<Py<AromaticSystemConstraintForm>>,
        new: Option<Py<AromaticSystemConstraintForm>>,
    },
    ModifyMulticenterBondConstraints {
        key: Py<MulticenterBondConstraintKey>,
        changes: TupleItems<MulticenterBondConstraintChange>,
    },
    ModifyMulticenterBondConstraint {
        id: HandleLike,
        old: Option<Py<MulticenterBondConstraintForm>>,
        new: Option<Py<MulticenterBondConstraintForm>>,
    },
    ModifyNoncovalentBondConstraints {
        key: Py<NoncovalentBondConstraintKey>,
        changes: TupleItems<NoncovalentBondConstraintChange>,
    },
    ModifyNoncovalentBondConstraint {
        id: HandleLike,
        old: Option<Py<NoncovalentBondConstraintForm>>,
        new: Option<Py<NoncovalentBondConstraintForm>>,
    },
    #[pyo3(constructor = (kind, key, changes))]
    ModifyStereoAtomConstraints {
        kind: Option<StereoKind>,
        key: Py<StereoAtomConstraintKey>,
        changes: TupleItems<StereoAtomConstraintChange>,
    },
    ModifyStereoAtomConstraint {
        id: HandleLike,
        kind: Option<StereoKind>,
        old: Option<Py<StereoAtomConstraintForm>>,
        new: Option<Py<StereoAtomConstraintForm>>,
    },
    #[pyo3(constructor = (kind, key, changes))]
    ModifyStereoBondConstraints {
        kind: Option<StereoKind>,
        key: Py<StereoBondConstraintKey>,
        changes: TupleItems<StereoBondConstraintChange>,
    },
    ModifyStereoBondConstraint {
        id: HandleLike,
        kind: Option<StereoKind>,
        old: Option<Py<StereoBondConstraintForm>>,
        new: Option<Py<StereoBondConstraintForm>>,
    },
    AddMoleculeConstraint {
        constraint: Py<ConstraintEdit>,
    },
    RemoveMoleculeConstraint {
        constraint: Py<ConstraintEdit>,
    },
}

fn edit_variant_repr(
    object: &Bound<'_, PyAny>,
    variant: &str,
    fields: &[&str],
) -> PyResult<String> {
    let mut parts = Vec::with_capacity(fields.len());
    for field in fields {
        let value = object.getattr(*field)?.repr()?.extract::<String>()?;
        parts.push(format!("{field}={value}"));
    }
    Ok(format!("Edit.{variant}({})", parts.join(", ")))
}

#[pymethods]
impl Edit {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::AddAtoms { .. } => ("AddAtoms", &["atoms"]),
            Self::AddBonds { .. } => ("AddBonds", &["bonds"]),
            Self::RemoveTopology { .. } => ("RemoveTopology", &["atoms", "bonds"]),
            Self::ModifyAtomElements { .. } => ("ModifyAtomElements", &["changes"]),
            Self::ModifyAtomIsotopeMasses { .. } => ("ModifyAtomIsotopeMasses", &["changes"]),
            Self::ModifyAtomCharges { .. } => ("ModifyAtomCharges", &["changes"]),
            Self::ModifyAtomImplicitHydrogens { .. } => {
                ("ModifyAtomImplicitHydrogens", &["changes"])
            }
            Self::ModifyAtomLonePairs { .. } => ("ModifyAtomLonePairs", &["changes"]),
            Self::ModifyAtomUnpairedElectrons { .. } => {
                ("ModifyAtomUnpairedElectrons", &["changes"])
            }
            Self::ModifyAtomConstraints { .. } => ("ModifyAtomConstraints", &["key", "changes"]),
            Self::ModifyAtomField { .. } => ("ModifyAtomField", &["id", "change"]),
            Self::ModifyBondOrders { .. } => ("ModifyBondOrders", &["changes"]),
            Self::ModifyBondCharges { .. } => ("ModifyBondCharges", &["changes"]),
            Self::ModifyBondUnpairedElectrons { .. } => {
                ("ModifyBondUnpairedElectrons", &["changes"])
            }
            Self::ModifyBondField { .. } => ("ModifyBondField", &["id", "change"]),
            Self::AddDativeBond { .. } => ("AddDativeBond", &["donors", "acceptor", "attributes"]),
            Self::RemoveDativeBonds { .. } => ("RemoveDativeBonds", &["removes"]),
            Self::ReplaceDativeBondDonors { .. } => {
                ("ReplaceDativeBondDonors", &["id", "old", "new"])
            }
            Self::ReplaceDativeBondAcceptor { .. } => {
                ("ReplaceDativeBondAcceptor", &["id", "old", "new"])
            }
            Self::ModifyDativeBondOrders { .. } => ("ModifyDativeBondOrders", &["changes"]),
            Self::ModifyDativeBondField { .. } => ("ModifyDativeBondField", &["id", "change"]),
            Self::AddAromaticSystem { .. } => ("AddAromaticSystem", &["atoms", "attributes"]),
            Self::RemoveAromaticSystems { .. } => ("RemoveAromaticSystems", &["removes"]),
            Self::ReplaceAromaticSystemAtoms { .. } => {
                ("ReplaceAromaticSystemAtoms", &["id", "old", "new"])
            }
            Self::ModifyAromaticSystemElectrons { .. } => {
                ("ModifyAromaticSystemElectrons", &["changes"])
            }
            Self::ModifyAromaticSystemCharges { .. } => {
                ("ModifyAromaticSystemCharges", &["changes"])
            }
            Self::ModifyAromaticSystemUnpairedElectrons { .. } => {
                ("ModifyAromaticSystemUnpairedElectrons", &["changes"])
            }
            Self::ModifyAromaticSystemField { .. } => {
                ("ModifyAromaticSystemField", &["id", "change"])
            }
            Self::AddMulticenterBond { .. } => ("AddMulticenterBond", &["atoms", "attributes"]),
            Self::RemoveMulticenterBonds { .. } => ("RemoveMulticenterBonds", &["removes"]),
            Self::ReplaceMulticenterBondAtoms { .. } => {
                ("ReplaceMulticenterBondAtoms", &["id", "old", "new"])
            }
            Self::ModifyMulticenterBondElectrons { .. } => {
                ("ModifyMulticenterBondElectrons", &["changes"])
            }
            Self::ModifyMulticenterBondCharges { .. } => {
                ("ModifyMulticenterBondCharges", &["changes"])
            }
            Self::ModifyMulticenterBondUnpairedElectrons { .. } => {
                ("ModifyMulticenterBondUnpairedElectrons", &["changes"])
            }
            Self::ModifyMulticenterBondField { .. } => {
                ("ModifyMulticenterBondField", &["id", "change"])
            }
            Self::AddNoncovalentBond { .. } => ("AddNoncovalentBond", &["atoms", "attributes"]),
            Self::RemoveNoncovalentBonds { .. } => ("RemoveNoncovalentBonds", &["removes"]),
            Self::ReplaceNoncovalentBondAtoms { .. } => {
                ("ReplaceNoncovalentBondAtoms", &["id", "old", "new"])
            }
            Self::ModifyNoncovalentBondKinds { .. } => ("ModifyNoncovalentBondKinds", &["changes"]),
            Self::ModifyNoncovalentBondField { .. } => {
                ("ModifyNoncovalentBondField", &["id", "change"])
            }
            Self::AddStereoAtom { .. } => ("AddStereoAtom", &["site", "ligands", "attributes"]),
            Self::RemoveStereoAtoms { .. } => ("RemoveStereoAtoms", &["removes"]),
            Self::ReplaceStereoAtomSite { .. } => ("ReplaceStereoAtomSite", &["id", "old", "new"]),
            Self::ReplaceStereoAtomLigands { .. } => {
                ("ReplaceStereoAtomLigands", &["id", "old", "new"])
            }
            Self::ModifyStereoAtomConfigurations { .. } => {
                ("ModifyStereoAtomConfigurations", &["changes"])
            }
            Self::ModifyStereoAtomField { .. } => ("ModifyStereoAtomField", &["id", "change"]),
            Self::AddStereoBond { .. } => ("AddStereoBond", &["site", "ligands", "attributes"]),
            Self::RemoveStereoBonds { .. } => ("RemoveStereoBonds", &["removes"]),
            Self::ReplaceStereoBondSite { .. } => ("ReplaceStereoBondSite", &["id", "old", "new"]),
            Self::ReplaceStereoBondLigands { .. } => {
                ("ReplaceStereoBondLigands", &["id", "old", "new"])
            }
            Self::ModifyStereoBondConfigurations { .. } => {
                ("ModifyStereoBondConfigurations", &["changes"])
            }
            Self::ModifyStereoBondField { .. } => ("ModifyStereoBondField", &["id", "change"]),
            Self::ModifyAtomConstraint { .. } => ("ModifyAtomConstraint", &["id", "old", "new"]),
            Self::ModifyBondConstraints { .. } => ("ModifyBondConstraints", &["key", "changes"]),
            Self::ModifyBondConstraint { .. } => ("ModifyBondConstraint", &["id", "old", "new"]),
            Self::ModifyDativeBondConstraints { .. } => {
                ("ModifyDativeBondConstraints", &["key", "changes"])
            }
            Self::ModifyDativeBondConstraint { .. } => {
                ("ModifyDativeBondConstraint", &["id", "old", "new"])
            }
            Self::ModifyAromaticSystemConstraints { .. } => {
                ("ModifyAromaticSystemConstraints", &["key", "changes"])
            }
            Self::ModifyAromaticSystemConstraint { .. } => {
                ("ModifyAromaticSystemConstraint", &["id", "old", "new"])
            }
            Self::ModifyMulticenterBondConstraints { .. } => {
                ("ModifyMulticenterBondConstraints", &["key", "changes"])
            }
            Self::ModifyMulticenterBondConstraint { .. } => {
                ("ModifyMulticenterBondConstraint", &["id", "old", "new"])
            }
            Self::ModifyNoncovalentBondConstraints { .. } => {
                ("ModifyNoncovalentBondConstraints", &["key", "changes"])
            }
            Self::ModifyNoncovalentBondConstraint { .. } => {
                ("ModifyNoncovalentBondConstraint", &["id", "old", "new"])
            }
            Self::ModifyStereoAtomConstraints { .. } => {
                ("ModifyStereoAtomConstraints", &["kind", "key", "changes"])
            }
            Self::ModifyStereoAtomConstraint { .. } => {
                ("ModifyStereoAtomConstraint", &["id", "kind", "old", "new"])
            }
            Self::ModifyStereoBondConstraints { .. } => {
                ("ModifyStereoBondConstraints", &["kind", "key", "changes"])
            }
            Self::ModifyStereoBondConstraint { .. } => {
                ("ModifyStereoBondConstraint", &["id", "kind", "old", "new"])
            }
            Self::AddMoleculeConstraint { .. } => ("AddMoleculeConstraint", &["constraint"]),
            Self::RemoveMoleculeConstraint { .. } => ("RemoveMoleculeConstraint", &["constraint"]),
        };
        edit_variant_repr(slf.bind(py).as_any(), variant, fields)
    }
}

impl Edit {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API used by Edits snapshots"
    )]
    pub(crate) fn from_rust(py: Python<'_>, edit: &GraphIrEdit) -> PyResult<Self> {
        Ok(match edit {
            GraphIrEdit::AddAtoms { atoms } => Self::AddAtoms {
                atoms: atoms
                    .iter()
                    .map(|atom| Readonly::<AtomForm>::from_rust(py, atom))
                    .collect::<PyResult<_>>()?,
            },
            GraphIrEdit::AddBonds { bonds } => Self::AddBonds {
                bonds: bonds
                    .iter()
                    .map(|bond| {
                        Ok((
                            (
                                HandleLike::from_atom_handle(&bond.endpoints[0]),
                                HandleLike::from_atom_handle(&bond.endpoints[1]),
                            ),
                            Readonly::<BondForm>::from_rust(py, &bond.attributes)?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            GraphIrEdit::RemoveTopology { atoms, bonds } => Self::RemoveTopology {
                atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                bonds: bonds.iter().map(HandleLike::from_bond_handle).collect(),
            },
            GraphIrEdit::ModifyAtomElements { changes } => Self::ModifyAtomElements {
                changes: TupleItems(
                    changes
                        .iter()
                        .map(|(id, old, new)| {
                            Ok((
                                HandleLike::from_atom_handle(id),
                                into_py_variant(py, ElementForm::from_rust(old))?,
                                into_py_variant(py, ElementForm::from_rust(new))?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ModifyAtomIsotopeMasses { changes } => Self::ModifyAtomIsotopeMasses {
                changes: TupleItems(
                    changes
                        .iter()
                        .map(|(id, old, new)| {
                            Ok((
                                HandleLike::from_atom_handle(id),
                                into_py_variant(py, IsotopeMassForm::from_rust(old))?,
                                into_py_variant(py, IsotopeMassForm::from_rust(new))?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ModifyAtomCharges { changes } => Self::ModifyAtomCharges {
                changes: TupleItems(
                    changes
                        .iter()
                        .map(|(id, old, new)| {
                            Ok((
                                HandleLike::from_atom_handle(id),
                                into_py_variant(py, NumForm::from_rust(py, old)?)?,
                                into_py_variant(py, NumForm::from_rust(py, new)?)?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ModifyAtomImplicitHydrogens { changes } => {
                Self::ModifyAtomImplicitHydrogens {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_atom_handle(id),
                                    into_py_variant(py, NumForm::from_rust(py, old)?)?,
                                    into_py_variant(py, NumForm::from_rust(py, new)?)?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyAtomLonePairs { changes } => Self::ModifyAtomLonePairs {
                changes: TupleItems(
                    changes
                        .iter()
                        .map(|(id, old, new)| {
                            Ok((
                                HandleLike::from_atom_handle(id),
                                into_py_variant(py, NumForm::from_rust(py, old)?)?,
                                into_py_variant(py, NumForm::from_rust(py, new)?)?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ModifyAtomUnpairedElectrons { changes } => {
                Self::ModifyAtomUnpairedElectrons {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_atom_handle(id),
                                    into_py_variant(
                                        py,
                                        UnpairedElectronsForm::from_rust(py, old)?,
                                    )?,
                                    into_py_variant(
                                        py,
                                        UnpairedElectronsForm::from_rust(py, new)?,
                                    )?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyAtomConstraints { key, changes } => Self::ModifyAtomConstraints {
                key: into_py_variant(py, AtomConstraintKey::from_rust(py, key)?)?,
                changes: TupleItems(
                    changes
                        .iter()
                        .map(|(id, old, new)| {
                            Ok((
                                HandleLike::from_atom_handle(id),
                                old.as_ref()
                                    .map(|value| {
                                        into_py_variant(
                                            py,
                                            AtomConstraintForm::from_rust(py, value)?,
                                        )
                                    })
                                    .transpose()?,
                                new.as_ref()
                                    .map(|value| {
                                        into_py_variant(
                                            py,
                                            AtomConstraintForm::from_rust(py, value)?,
                                        )
                                    })
                                    .transpose()?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ModifyAtomField { id, change } => Self::ModifyAtomField {
                id: HandleLike::from_atom_handle(id),
                change: into_py_variant(py, AtomFieldChange::from_rust(py, change)?)?,
            },
            GraphIrEdit::ModifyBondOrders { changes } => Self::ModifyBondOrders {
                changes: TupleItems(
                    changes
                        .iter()
                        .map(|(id, old, new)| {
                            Ok((
                                HandleLike::from_bond_handle(id),
                                into_py_variant(py, NumForm::from_rust(py, old)?)?,
                                into_py_variant(py, NumForm::from_rust(py, new)?)?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ModifyBondCharges { changes } => Self::ModifyBondCharges {
                changes: TupleItems(
                    changes
                        .iter()
                        .map(|(id, old, new)| {
                            Ok((
                                HandleLike::from_bond_handle(id),
                                into_py_variant(py, NumForm::from_rust(py, old)?)?,
                                into_py_variant(py, NumForm::from_rust(py, new)?)?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ModifyBondUnpairedElectrons { changes } => {
                Self::ModifyBondUnpairedElectrons {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_bond_handle(id),
                                    into_py_variant(
                                        py,
                                        UnpairedElectronsForm::from_rust(py, old)?,
                                    )?,
                                    into_py_variant(
                                        py,
                                        UnpairedElectronsForm::from_rust(py, new)?,
                                    )?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyBondField { id, change } => Self::ModifyBondField {
                id: HandleLike::from_bond_handle(id),
                change: into_py_variant(py, BondFieldChange::from_rust(py, change)?)?,
            },
            GraphIrEdit::AddDativeBond {
                donors,
                acceptor,
                attributes,
            } => Self::AddDativeBond {
                donors: donors.iter().map(HandleLike::from_atom_handle).collect(),
                acceptor: HandleLike::from_atom_handle(acceptor),
                attributes: Readonly::<DativeBondForm>::from_rust(py, attributes)?,
            },
            GraphIrEdit::RemoveDativeBonds { removes } => Self::RemoveDativeBonds {
                removes: removes
                    .iter()
                    .map(|(id, donors, acceptor, attributes)| {
                        Ok((
                            HandleLike::from_dative_bond_handle(id),
                            donors.iter().map(HandleLike::from_atom_handle).collect(),
                            HandleLike::from_atom_handle(acceptor),
                            Readonly::<DativeBondForm>::from_rust(py, attributes)?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            GraphIrEdit::ReplaceDativeBondDonors { id, old, new } => {
                Self::ReplaceDativeBondDonors {
                    id: HandleLike::from_dative_bond_handle(id),
                    old: old.iter().map(HandleLike::from_atom_handle).collect(),
                    new: new.iter().map(HandleLike::from_atom_handle).collect(),
                }
            }
            GraphIrEdit::ReplaceDativeBondAcceptor { id, old, new } => {
                Self::ReplaceDativeBondAcceptor {
                    id: HandleLike::from_dative_bond_handle(id),
                    old: HandleLike::from_atom_handle(old),
                    new: HandleLike::from_atom_handle(new),
                }
            }
            GraphIrEdit::ModifyDativeBondOrders { changes } => Self::ModifyDativeBondOrders {
                changes: TupleItems(
                    changes
                        .iter()
                        .map(|(id, old, new)| {
                            Ok((
                                HandleLike::from_dative_bond_handle(id),
                                into_py_variant(py, NumForm::from_rust(py, old)?)?,
                                into_py_variant(py, NumForm::from_rust(py, new)?)?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ModifyDativeBondField { id, change } => Self::ModifyDativeBondField {
                id: HandleLike::from_dative_bond_handle(id),
                change: into_py_variant(py, DativeBondFieldChange::from_rust(py, change)?)?,
            },
            GraphIrEdit::AddAromaticSystem(addition) => {
                let GraphIrAddAromaticSystem { atoms, attributes } = addition.as_ref();
                Self::AddAromaticSystem {
                    atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                    attributes: Readonly::<AromaticSystemForm>::from_rust(py, attributes)?,
                }
            }
            GraphIrEdit::RemoveAromaticSystems { removes } => Self::RemoveAromaticSystems {
                removes: removes
                    .iter()
                    .map(|(id, atoms, attributes)| {
                        Ok((
                            HandleLike::from_aromatic_system_handle(id),
                            atoms.iter().map(HandleLike::from_atom_handle).collect(),
                            Readonly::<AromaticSystemForm>::from_rust(py, attributes)?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            GraphIrEdit::ReplaceAromaticSystemAtoms { id, old, new } => {
                Self::ReplaceAromaticSystemAtoms {
                    id: HandleLike::from_aromatic_system_handle(id),
                    old: old.iter().map(HandleLike::from_atom_handle).collect(),
                    new: new.iter().map(HandleLike::from_atom_handle).collect(),
                }
            }
            GraphIrEdit::ModifyAromaticSystemElectrons { changes } => {
                Self::ModifyAromaticSystemElectrons {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_aromatic_system_handle(id),
                                    into_py_variant(py, ElectronCountsForm::from_rust(old))?,
                                    into_py_variant(py, ElectronCountsForm::from_rust(new))?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyAromaticSystemCharges { changes } => {
                Self::ModifyAromaticSystemCharges {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_aromatic_system_handle(id),
                                    into_py_variant(py, NumForm::from_rust(py, old)?)?,
                                    into_py_variant(py, NumForm::from_rust(py, new)?)?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyAromaticSystemUnpairedElectrons { changes } => {
                Self::ModifyAromaticSystemUnpairedElectrons {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_aromatic_system_handle(id),
                                    into_py_variant(
                                        py,
                                        UnpairedElectronsForm::from_rust(py, old)?,
                                    )?,
                                    into_py_variant(
                                        py,
                                        UnpairedElectronsForm::from_rust(py, new)?,
                                    )?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyAromaticSystemField { id, change } => {
                Self::ModifyAromaticSystemField {
                    id: HandleLike::from_aromatic_system_handle(id),
                    change: into_py_variant(py, AromaticSystemFieldChange::from_rust(py, change)?)?,
                }
            }
            GraphIrEdit::AddMulticenterBond(addition) => {
                let GraphIrAddMulticenterBond { atoms, attributes } = addition.as_ref();
                Self::AddMulticenterBond {
                    atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                    attributes: Readonly::<MulticenterBondForm>::from_rust(py, attributes)?,
                }
            }
            GraphIrEdit::RemoveMulticenterBonds { removes } => Self::RemoveMulticenterBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, attributes)| {
                        Ok((
                            HandleLike::from_multicenter_bond_handle(id),
                            atoms.iter().map(HandleLike::from_atom_handle).collect(),
                            Readonly::<MulticenterBondForm>::from_rust(py, attributes)?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            GraphIrEdit::ReplaceMulticenterBondAtoms { id, old, new } => {
                Self::ReplaceMulticenterBondAtoms {
                    id: HandleLike::from_multicenter_bond_handle(id),
                    old: old.iter().map(HandleLike::from_atom_handle).collect(),
                    new: new.iter().map(HandleLike::from_atom_handle).collect(),
                }
            }
            GraphIrEdit::ModifyMulticenterBondElectrons { changes } => {
                Self::ModifyMulticenterBondElectrons {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_multicenter_bond_handle(id),
                                    into_py_variant(py, ElectronCountsForm::from_rust(old))?,
                                    into_py_variant(py, ElectronCountsForm::from_rust(new))?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyMulticenterBondCharges { changes } => {
                Self::ModifyMulticenterBondCharges {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_multicenter_bond_handle(id),
                                    into_py_variant(py, NumForm::from_rust(py, old)?)?,
                                    into_py_variant(py, NumForm::from_rust(py, new)?)?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyMulticenterBondUnpairedElectrons { changes } => {
                Self::ModifyMulticenterBondUnpairedElectrons {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_multicenter_bond_handle(id),
                                    into_py_variant(
                                        py,
                                        UnpairedElectronsForm::from_rust(py, old)?,
                                    )?,
                                    into_py_variant(
                                        py,
                                        UnpairedElectronsForm::from_rust(py, new)?,
                                    )?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyMulticenterBondField { id, change } => {
                Self::ModifyMulticenterBondField {
                    id: HandleLike::from_multicenter_bond_handle(id),
                    change: into_py_variant(
                        py,
                        MulticenterBondFieldChange::from_rust(py, change)?,
                    )?,
                }
            }
            GraphIrEdit::AddNoncovalentBond { atoms, attributes } => Self::AddNoncovalentBond {
                atoms: (
                    HandleLike::from_atom_handle(&atoms[0]),
                    HandleLike::from_atom_handle(&atoms[1]),
                ),
                attributes: Readonly::<NoncovalentBondForm>::from_rust(py, attributes)?,
            },
            GraphIrEdit::RemoveNoncovalentBonds { removes } => Self::RemoveNoncovalentBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, attributes)| {
                        Ok((
                            HandleLike::from_noncovalent_bond_handle(id),
                            (
                                HandleLike::from_atom_handle(&atoms[0]),
                                HandleLike::from_atom_handle(&atoms[1]),
                            ),
                            Readonly::<NoncovalentBondForm>::from_rust(py, attributes)?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            GraphIrEdit::ReplaceNoncovalentBondAtoms { id, old, new } => {
                Self::ReplaceNoncovalentBondAtoms {
                    id: HandleLike::from_noncovalent_bond_handle(id),
                    old: (
                        HandleLike::from_atom_handle(&old[0]),
                        HandleLike::from_atom_handle(&old[1]),
                    ),
                    new: (
                        HandleLike::from_atom_handle(&new[0]),
                        HandleLike::from_atom_handle(&new[1]),
                    ),
                }
            }
            GraphIrEdit::ModifyNoncovalentBondKinds { changes } => {
                Self::ModifyNoncovalentBondKinds {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_noncovalent_bond_handle(id),
                                    into_py_variant(py, NoncovalentBondKindForm::from_rust(old))?,
                                    into_py_variant(py, NoncovalentBondKindForm::from_rust(new))?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyNoncovalentBondField { id, change } => {
                Self::ModifyNoncovalentBondField {
                    id: HandleLike::from_noncovalent_bond_handle(id),
                    change: into_py_variant(
                        py,
                        NoncovalentBondFieldChange::from_rust(py, change)?,
                    )?,
                }
            }
            GraphIrEdit::AddStereoAtom(addition) => {
                let GraphIrAddStereoAtom {
                    site,
                    ligands,
                    attributes,
                } = addition.as_ref();
                Self::AddStereoAtom {
                    site: HandleLike::from_atom_handle(site),
                    ligands: ligands
                        .iter()
                        .map(|(atom, kind)| {
                            (
                                HandleLike::from_atom_handle(atom),
                                StereoLigandKind::from_rust(*kind),
                            )
                        })
                        .collect(),
                    attributes: Readonly::<StereoAtomForm>::from_rust(py, attributes)?,
                }
            }
            GraphIrEdit::RemoveStereoAtoms { removes } => Self::RemoveStereoAtoms {
                removes: ReadonlyStereoAtomRemovals(
                    removes
                        .iter()
                        .map(|(id, site, ligands, attributes)| {
                            Ok((
                                HandleLike::from_stereo_atom_handle(id),
                                HandleLike::from_atom_handle(site),
                                ligands
                                    .iter()
                                    .map(|(atom, kind)| {
                                        (
                                            HandleLike::from_atom_handle(atom),
                                            StereoLigandKind::from_rust(*kind),
                                        )
                                    })
                                    .collect(),
                                Readonly::<StereoAtomForm>::from_rust(py, attributes)?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ReplaceStereoAtomSite { id, old, new } => Self::ReplaceStereoAtomSite {
                id: HandleLike::from_stereo_atom_handle(id),
                old: HandleLike::from_atom_handle(old),
                new: HandleLike::from_atom_handle(new),
            },
            GraphIrEdit::ReplaceStereoAtomLigands { id, old, new } => {
                Self::ReplaceStereoAtomLigands {
                    id: HandleLike::from_stereo_atom_handle(id),
                    old: old
                        .iter()
                        .map(|(atom, kind)| {
                            (
                                HandleLike::from_atom_handle(atom),
                                StereoLigandKind::from_rust(*kind),
                            )
                        })
                        .collect(),
                    new: new
                        .iter()
                        .map(|(atom, kind)| {
                            (
                                HandleLike::from_atom_handle(atom),
                                StereoLigandKind::from_rust(*kind),
                            )
                        })
                        .collect(),
                }
            }
            GraphIrEdit::ModifyStereoAtomConfigurations { changes } => {
                Self::ModifyStereoAtomConfigurations {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_stereo_atom_handle(id),
                                    into_py_variant(
                                        py,
                                        StereoConfigurationForm::from_rust(py, old)?,
                                    )?,
                                    into_py_variant(
                                        py,
                                        StereoConfigurationForm::from_rust(py, new)?,
                                    )?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyStereoAtomField { id, change } => Self::ModifyStereoAtomField {
                id: HandleLike::from_stereo_atom_handle(id),
                change: into_py_variant(py, StereoAtomFieldChange::from_rust(py, change)?)?,
            },
            GraphIrEdit::AddStereoBond(addition) => {
                let GraphIrAddStereoBond {
                    site,
                    ligands,
                    attributes,
                } = addition.as_ref();
                Self::AddStereoBond {
                    site: HandleLike::from_bond_handle(site),
                    ligands: ligands
                        .iter()
                        .map(|(atom, kind)| {
                            (
                                HandleLike::from_atom_handle(atom),
                                StereoLigandKind::from_rust(*kind),
                            )
                        })
                        .collect(),
                    attributes: Readonly::<StereoBondForm>::from_rust(py, attributes)?,
                }
            }
            GraphIrEdit::RemoveStereoBonds { removes } => Self::RemoveStereoBonds {
                removes: ReadonlyStereoBondRemovals(
                    removes
                        .iter()
                        .map(|(id, site, ligands, attributes)| {
                            Ok((
                                HandleLike::from_stereo_bond_handle(id),
                                HandleLike::from_bond_handle(site),
                                ligands
                                    .iter()
                                    .map(|(atom, kind)| {
                                        (
                                            HandleLike::from_atom_handle(atom),
                                            StereoLigandKind::from_rust(*kind),
                                        )
                                    })
                                    .collect(),
                                Readonly::<StereoBondForm>::from_rust(py, attributes)?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ReplaceStereoBondSite { id, old, new } => Self::ReplaceStereoBondSite {
                id: HandleLike::from_stereo_bond_handle(id),
                old: HandleLike::from_bond_handle(old),
                new: HandleLike::from_bond_handle(new),
            },
            GraphIrEdit::ReplaceStereoBondLigands { id, old, new } => {
                Self::ReplaceStereoBondLigands {
                    id: HandleLike::from_stereo_bond_handle(id),
                    old: old
                        .iter()
                        .map(|(atom, kind)| {
                            (
                                HandleLike::from_atom_handle(atom),
                                StereoLigandKind::from_rust(*kind),
                            )
                        })
                        .collect(),
                    new: new
                        .iter()
                        .map(|(atom, kind)| {
                            (
                                HandleLike::from_atom_handle(atom),
                                StereoLigandKind::from_rust(*kind),
                            )
                        })
                        .collect(),
                }
            }
            GraphIrEdit::ModifyStereoBondConfigurations { changes } => {
                Self::ModifyStereoBondConfigurations {
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_stereo_bond_handle(id),
                                    into_py_variant(
                                        py,
                                        StereoConfigurationForm::from_rust(py, old)?,
                                    )?,
                                    into_py_variant(
                                        py,
                                        StereoConfigurationForm::from_rust(py, new)?,
                                    )?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyStereoBondField { id, change } => Self::ModifyStereoBondField {
                id: HandleLike::from_stereo_bond_handle(id),
                change: into_py_variant(py, StereoBondFieldChange::from_rust(py, change)?)?,
            },
            GraphIrEdit::ModifyAtomConstraint { id, old, new } => Self::ModifyAtomConstraint {
                id: HandleLike::from_atom_handle(id),
                old: old
                    .as_ref()
                    .map(|value| into_py_variant(py, AtomConstraintForm::from_rust(py, value)?))
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|value| into_py_variant(py, AtomConstraintForm::from_rust(py, value)?))
                    .transpose()?,
            },
            GraphIrEdit::ModifyBondConstraints { key, changes } => Self::ModifyBondConstraints {
                key: into_py_variant(py, BondConstraintKey::from_rust(py, key)?)?,
                changes: TupleItems(
                    changes
                        .iter()
                        .map(|(id, old, new)| {
                            Ok((
                                HandleLike::from_bond_handle(id),
                                old.as_ref()
                                    .map(|value| {
                                        into_py_variant(
                                            py,
                                            BondConstraintForm::from_rust(py, value)?,
                                        )
                                    })
                                    .transpose()?,
                                new.as_ref()
                                    .map(|value| {
                                        into_py_variant(
                                            py,
                                            BondConstraintForm::from_rust(py, value)?,
                                        )
                                    })
                                    .transpose()?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ModifyBondConstraint { id, old, new } => Self::ModifyBondConstraint {
                id: HandleLike::from_bond_handle(id),
                old: old
                    .as_ref()
                    .map(|value| into_py_variant(py, BondConstraintForm::from_rust(py, value)?))
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|value| into_py_variant(py, BondConstraintForm::from_rust(py, value)?))
                    .transpose()?,
            },
            GraphIrEdit::ModifyDativeBondConstraints { key, changes } => {
                Self::ModifyDativeBondConstraints {
                    key: into_py_variant(py, DativeBondConstraintKey::from_rust(py, key)?)?,
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_dative_bond_handle(id),
                                    old.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                DativeBondConstraintForm::from_rust(py, value)?,
                                            )
                                        })
                                        .transpose()?,
                                    new.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                DativeBondConstraintForm::from_rust(py, value)?,
                                            )
                                        })
                                        .transpose()?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyDativeBondConstraint { id, old, new } => {
                Self::ModifyDativeBondConstraint {
                    id: HandleLike::from_dative_bond_handle(id),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, DativeBondConstraintForm::from_rust(py, value)?)
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, DativeBondConstraintForm::from_rust(py, value)?)
                        })
                        .transpose()?,
                }
            }
            GraphIrEdit::ModifyAromaticSystemConstraints { key, changes } => {
                Self::ModifyAromaticSystemConstraints {
                    key: into_py_variant(py, AromaticSystemConstraintKey::from_rust(key))?,
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_aromatic_system_handle(id),
                                    old.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                AromaticSystemConstraintForm::from_rust(py, value)?,
                                            )
                                        })
                                        .transpose()?,
                                    new.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                AromaticSystemConstraintForm::from_rust(py, value)?,
                                            )
                                        })
                                        .transpose()?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyAromaticSystemConstraint { id, old, new } => {
                Self::ModifyAromaticSystemConstraint {
                    id: HandleLike::from_aromatic_system_handle(id),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, AromaticSystemConstraintForm::from_rust(py, value)?)
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, AromaticSystemConstraintForm::from_rust(py, value)?)
                        })
                        .transpose()?,
                }
            }
            GraphIrEdit::ModifyMulticenterBondConstraints { key, changes } => {
                Self::ModifyMulticenterBondConstraints {
                    key: into_py_variant(py, MulticenterBondConstraintKey::from_rust(key))?,
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_multicenter_bond_handle(id),
                                    old.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                MulticenterBondConstraintForm::from_rust(
                                                    py, value,
                                                )?,
                                            )
                                        })
                                        .transpose()?,
                                    new.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                MulticenterBondConstraintForm::from_rust(
                                                    py, value,
                                                )?,
                                            )
                                        })
                                        .transpose()?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyMulticenterBondConstraint { id, old, new } => {
                Self::ModifyMulticenterBondConstraint {
                    id: HandleLike::from_multicenter_bond_handle(id),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(
                                py,
                                MulticenterBondConstraintForm::from_rust(py, value)?,
                            )
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(
                                py,
                                MulticenterBondConstraintForm::from_rust(py, value)?,
                            )
                        })
                        .transpose()?,
                }
            }
            GraphIrEdit::ModifyNoncovalentBondConstraints { key, changes } => {
                Self::ModifyNoncovalentBondConstraints {
                    key: into_py_variant(py, NoncovalentBondConstraintKey::from_rust(key))?,
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_noncovalent_bond_handle(id),
                                    old.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                NoncovalentBondConstraintForm::from_rust(
                                                    py, value,
                                                )?,
                                            )
                                        })
                                        .transpose()?,
                                    new.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                NoncovalentBondConstraintForm::from_rust(
                                                    py, value,
                                                )?,
                                            )
                                        })
                                        .transpose()?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyNoncovalentBondConstraint { id, old, new } => {
                Self::ModifyNoncovalentBondConstraint {
                    id: HandleLike::from_noncovalent_bond_handle(id),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(
                                py,
                                NoncovalentBondConstraintForm::from_rust(py, value)?,
                            )
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(
                                py,
                                NoncovalentBondConstraintForm::from_rust(py, value)?,
                            )
                        })
                        .transpose()?,
                }
            }
            GraphIrEdit::ModifyStereoAtomConstraints { kind, key, changes } => {
                Self::ModifyStereoAtomConstraints {
                    kind: kind.map(StereoKind::from_rust),
                    key: into_py_variant(py, StereoAtomConstraintKey::from_rust(py, key)?)?,
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_stereo_atom_handle(id),
                                    old.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                StereoAtomConstraintForm::from_rust(py, value)?,
                                            )
                                        })
                                        .transpose()?,
                                    new.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                StereoAtomConstraintForm::from_rust(py, value)?,
                                            )
                                        })
                                        .transpose()?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyStereoAtomConstraint { id, kind, old, new } => {
                Self::ModifyStereoAtomConstraint {
                    id: HandleLike::from_stereo_atom_handle(id),
                    kind: kind.map(StereoKind::from_rust),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, StereoAtomConstraintForm::from_rust(py, value)?)
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, StereoAtomConstraintForm::from_rust(py, value)?)
                        })
                        .transpose()?,
                }
            }
            GraphIrEdit::ModifyStereoBondConstraints { kind, key, changes } => {
                Self::ModifyStereoBondConstraints {
                    kind: kind.map(StereoKind::from_rust),
                    key: into_py_variant(py, StereoBondConstraintKey::from_rust(py, key)?)?,
                    changes: TupleItems(
                        changes
                            .iter()
                            .map(|(id, old, new)| {
                                Ok((
                                    HandleLike::from_stereo_bond_handle(id),
                                    old.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                StereoBondConstraintForm::from_rust(py, value)?,
                                            )
                                        })
                                        .transpose()?,
                                    new.as_ref()
                                        .map(|value| {
                                            into_py_variant(
                                                py,
                                                StereoBondConstraintForm::from_rust(py, value)?,
                                            )
                                        })
                                        .transpose()?,
                                ))
                            })
                            .collect::<PyResult<_>>()?,
                    ),
                }
            }
            GraphIrEdit::ModifyStereoBondConstraint { id, kind, old, new } => {
                Self::ModifyStereoBondConstraint {
                    id: HandleLike::from_stereo_bond_handle(id),
                    kind: kind.map(StereoKind::from_rust),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, StereoBondConstraintForm::from_rust(py, value)?)
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, StereoBondConstraintForm::from_rust(py, value)?)
                        })
                        .transpose()?,
                }
            }
            GraphIrEdit::AddMoleculeConstraint { constraint } => Self::AddMoleculeConstraint {
                constraint: Py::new(py, ConstraintEdit::from_rust(constraint))?,
            },
            GraphIrEdit::RemoveMoleculeConstraint { constraint } => {
                Self::RemoveMoleculeConstraint {
                    constraint: Py::new(py, ConstraintEdit::from_rust(constraint))?,
                }
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrEdit {
        match self {
            Self::AddAtoms { atoms } => GraphIrEdit::AddAtoms {
                atoms: atoms.iter().map(|atom| atom.to_rust(py)).collect(),
            },
            Self::AddBonds { bonds } => GraphIrEdit::AddBonds {
                bonds: bonds
                    .iter()
                    .map(|((first, second), attributes)| GraphIrAddBond {
                        endpoints: [first.to_atom_handle(), second.to_atom_handle()],
                        attributes: attributes.to_rust(py),
                    })
                    .collect(),
            },
            Self::RemoveTopology { atoms, bonds } => GraphIrEdit::RemoveTopology {
                atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                bonds: bonds.iter().map(HandleLike::to_bond_handle).collect(),
            },
            Self::ModifyAtomElements { changes } => GraphIrEdit::ModifyAtomElements {
                changes: changes
                    .0
                    .iter()
                    .map(|(id, old, new)| {
                        (
                            id.to_atom_handle(),
                            old.bind(py).borrow().to_rust(),
                            new.bind(py).borrow().to_rust(),
                        )
                    })
                    .collect(),
            },
            Self::ModifyAtomIsotopeMasses { changes } => GraphIrEdit::ModifyAtomIsotopeMasses {
                changes: changes
                    .0
                    .iter()
                    .map(|(id, old, new)| {
                        (
                            id.to_atom_handle(),
                            old.bind(py).borrow().to_rust(),
                            new.bind(py).borrow().to_rust(),
                        )
                    })
                    .collect(),
            },
            Self::ModifyAtomCharges { changes } => GraphIrEdit::ModifyAtomCharges {
                changes: changes
                    .0
                    .iter()
                    .map(|(id, old, new)| {
                        (
                            id.to_atom_handle(),
                            old.bind(py).borrow().to_rust(py),
                            new.bind(py).borrow().to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ModifyAtomImplicitHydrogens { changes } => {
                GraphIrEdit::ModifyAtomImplicitHydrogens {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_atom_handle(),
                                old.bind(py).borrow().to_rust(py),
                                new.bind(py).borrow().to_rust(py),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyAtomLonePairs { changes } => GraphIrEdit::ModifyAtomLonePairs {
                changes: changes
                    .0
                    .iter()
                    .map(|(id, old, new)| {
                        (
                            id.to_atom_handle(),
                            old.bind(py).borrow().to_rust(py),
                            new.bind(py).borrow().to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ModifyAtomUnpairedElectrons { changes } => {
                GraphIrEdit::ModifyAtomUnpairedElectrons {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_atom_handle(),
                                old.bind(py).borrow().to_rust(py),
                                new.bind(py).borrow().to_rust(py),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyAtomConstraints { key, changes } => GraphIrEdit::ModifyAtomConstraints {
                key: key.bind(py).borrow().to_rust(py),
                changes: changes
                    .0
                    .iter()
                    .map(|(id, old, new)| {
                        (
                            id.to_atom_handle(),
                            old.as_ref()
                                .map(|value| value.bind(py).borrow().to_rust(py)),
                            new.as_ref()
                                .map(|value| value.bind(py).borrow().to_rust(py)),
                        )
                    })
                    .collect(),
            },
            Self::ModifyAtomField { id, change } => GraphIrEdit::ModifyAtomField {
                id: id.to_atom_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyBondOrders { changes } => GraphIrEdit::ModifyBondOrders {
                changes: changes
                    .0
                    .iter()
                    .map(|(id, old, new)| {
                        (
                            id.to_bond_handle(),
                            old.bind(py).borrow().to_rust(py),
                            new.bind(py).borrow().to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ModifyBondCharges { changes } => GraphIrEdit::ModifyBondCharges {
                changes: changes
                    .0
                    .iter()
                    .map(|(id, old, new)| {
                        (
                            id.to_bond_handle(),
                            old.bind(py).borrow().to_rust(py),
                            new.bind(py).borrow().to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ModifyBondUnpairedElectrons { changes } => {
                GraphIrEdit::ModifyBondUnpairedElectrons {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_bond_handle(),
                                old.bind(py).borrow().to_rust(py),
                                new.bind(py).borrow().to_rust(py),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyBondField { id, change } => GraphIrEdit::ModifyBondField {
                id: id.to_bond_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::AddDativeBond {
                donors,
                acceptor,
                attributes,
            } => GraphIrEdit::AddDativeBond {
                donors: donors.iter().map(HandleLike::to_atom_handle).collect(),
                acceptor: acceptor.to_atom_handle(),
                attributes: attributes.to_rust(py),
            },
            Self::RemoveDativeBonds { removes } => GraphIrEdit::RemoveDativeBonds {
                removes: removes
                    .iter()
                    .map(|(id, donors, acceptor, attributes)| {
                        (
                            id.to_dative_bond_handle(),
                            donors.iter().map(HandleLike::to_atom_handle).collect(),
                            acceptor.to_atom_handle(),
                            attributes.to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ReplaceDativeBondDonors { id, old, new } => {
                GraphIrEdit::ReplaceDativeBondDonors {
                    id: id.to_dative_bond_handle(),
                    old: old.iter().map(HandleLike::to_atom_handle).collect(),
                    new: new.iter().map(HandleLike::to_atom_handle).collect(),
                }
            }
            Self::ReplaceDativeBondAcceptor { id, old, new } => {
                GraphIrEdit::ReplaceDativeBondAcceptor {
                    id: id.to_dative_bond_handle(),
                    old: old.to_atom_handle(),
                    new: new.to_atom_handle(),
                }
            }
            Self::ModifyDativeBondOrders { changes } => GraphIrEdit::ModifyDativeBondOrders {
                changes: changes
                    .0
                    .iter()
                    .map(|(id, old, new)| {
                        (
                            id.to_dative_bond_handle(),
                            old.bind(py).borrow().to_rust(py),
                            new.bind(py).borrow().to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ModifyDativeBondField { id, change } => GraphIrEdit::ModifyDativeBondField {
                id: id.to_dative_bond_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::AddAromaticSystem { atoms, attributes } => {
                GraphIrEdit::AddAromaticSystem(Box::new(GraphIrAddAromaticSystem {
                    atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                    attributes: attributes.to_rust(py),
                }))
            }
            Self::RemoveAromaticSystems { removes } => GraphIrEdit::RemoveAromaticSystems {
                removes: removes
                    .iter()
                    .map(|(id, atoms, attributes)| {
                        (
                            id.to_aromatic_system_handle(),
                            atoms.iter().map(HandleLike::to_atom_handle).collect(),
                            attributes.to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ReplaceAromaticSystemAtoms { id, old, new } => {
                GraphIrEdit::ReplaceAromaticSystemAtoms {
                    id: id.to_aromatic_system_handle(),
                    old: old.iter().map(HandleLike::to_atom_handle).collect(),
                    new: new.iter().map(HandleLike::to_atom_handle).collect(),
                }
            }
            Self::ModifyAromaticSystemElectrons { changes } => {
                GraphIrEdit::ModifyAromaticSystemElectrons {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_aromatic_system_handle(),
                                old.bind(py).borrow().to_rust(),
                                new.bind(py).borrow().to_rust(),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyAromaticSystemCharges { changes } => {
                GraphIrEdit::ModifyAromaticSystemCharges {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_aromatic_system_handle(),
                                old.bind(py).borrow().to_rust(py),
                                new.bind(py).borrow().to_rust(py),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyAromaticSystemUnpairedElectrons { changes } => {
                GraphIrEdit::ModifyAromaticSystemUnpairedElectrons {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_aromatic_system_handle(),
                                old.bind(py).borrow().to_rust(py),
                                new.bind(py).borrow().to_rust(py),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyAromaticSystemField { id, change } => {
                GraphIrEdit::ModifyAromaticSystemField {
                    id: id.to_aromatic_system_handle(),
                    change: change.bind(py).borrow().to_rust(py),
                }
            }
            Self::AddMulticenterBond { atoms, attributes } => {
                GraphIrEdit::AddMulticenterBond(Box::new(GraphIrAddMulticenterBond {
                    atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                    attributes: attributes.to_rust(py),
                }))
            }
            Self::RemoveMulticenterBonds { removes } => GraphIrEdit::RemoveMulticenterBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, attributes)| {
                        (
                            id.to_multicenter_bond_handle(),
                            atoms.iter().map(HandleLike::to_atom_handle).collect(),
                            attributes.to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ReplaceMulticenterBondAtoms { id, old, new } => {
                GraphIrEdit::ReplaceMulticenterBondAtoms {
                    id: id.to_multicenter_bond_handle(),
                    old: old.iter().map(HandleLike::to_atom_handle).collect(),
                    new: new.iter().map(HandleLike::to_atom_handle).collect(),
                }
            }
            Self::ModifyMulticenterBondElectrons { changes } => {
                GraphIrEdit::ModifyMulticenterBondElectrons {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_multicenter_bond_handle(),
                                old.bind(py).borrow().to_rust(),
                                new.bind(py).borrow().to_rust(),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyMulticenterBondCharges { changes } => {
                GraphIrEdit::ModifyMulticenterBondCharges {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_multicenter_bond_handle(),
                                old.bind(py).borrow().to_rust(py),
                                new.bind(py).borrow().to_rust(py),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyMulticenterBondUnpairedElectrons { changes } => {
                GraphIrEdit::ModifyMulticenterBondUnpairedElectrons {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_multicenter_bond_handle(),
                                old.bind(py).borrow().to_rust(py),
                                new.bind(py).borrow().to_rust(py),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyMulticenterBondField { id, change } => {
                GraphIrEdit::ModifyMulticenterBondField {
                    id: id.to_multicenter_bond_handle(),
                    change: change.bind(py).borrow().to_rust(py),
                }
            }
            Self::AddNoncovalentBond { atoms, attributes } => GraphIrEdit::AddNoncovalentBond {
                atoms: [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                attributes: attributes.to_rust(py),
            },
            Self::RemoveNoncovalentBonds { removes } => GraphIrEdit::RemoveNoncovalentBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, attributes)| {
                        (
                            id.to_noncovalent_bond_handle(),
                            [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                            attributes.to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ReplaceNoncovalentBondAtoms { id, old, new } => {
                GraphIrEdit::ReplaceNoncovalentBondAtoms {
                    id: id.to_noncovalent_bond_handle(),
                    old: [old.0.to_atom_handle(), old.1.to_atom_handle()],
                    new: [new.0.to_atom_handle(), new.1.to_atom_handle()],
                }
            }
            Self::ModifyNoncovalentBondKinds { changes } => {
                GraphIrEdit::ModifyNoncovalentBondKinds {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_noncovalent_bond_handle(),
                                old.bind(py).borrow().to_rust(),
                                new.bind(py).borrow().to_rust(),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyNoncovalentBondField { id, change } => {
                GraphIrEdit::ModifyNoncovalentBondField {
                    id: id.to_noncovalent_bond_handle(),
                    change: change.bind(py).borrow().to_rust(py),
                }
            }
            Self::AddStereoAtom {
                site,
                ligands,
                attributes,
            } => GraphIrEdit::AddStereoAtom(Box::new(GraphIrAddStereoAtom {
                site: site.to_atom_handle(),
                ligands: ligands
                    .iter()
                    .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                    .collect(),
                attributes: attributes.to_rust(py),
            })),
            Self::RemoveStereoAtoms { removes } => GraphIrEdit::RemoveStereoAtoms {
                removes: removes
                    .0
                    .iter()
                    .map(|(id, site, ligands, attributes)| {
                        (
                            id.to_stereo_atom_handle(),
                            site.to_atom_handle(),
                            ligands
                                .iter()
                                .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                                .collect(),
                            attributes.to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ReplaceStereoAtomSite { id, old, new } => GraphIrEdit::ReplaceStereoAtomSite {
                id: id.to_stereo_atom_handle(),
                old: old.to_atom_handle(),
                new: new.to_atom_handle(),
            },
            Self::ReplaceStereoAtomLigands { id, old, new } => {
                GraphIrEdit::ReplaceStereoAtomLigands {
                    id: id.to_stereo_atom_handle(),
                    old: old
                        .iter()
                        .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                        .collect(),
                    new: new
                        .iter()
                        .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                        .collect(),
                }
            }
            Self::ModifyStereoAtomConfigurations { changes } => {
                GraphIrEdit::ModifyStereoAtomConfigurations {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_stereo_atom_handle(),
                                old.bind(py).borrow().to_rust(py),
                                new.bind(py).borrow().to_rust(py),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyStereoAtomField { id, change } => GraphIrEdit::ModifyStereoAtomField {
                id: id.to_stereo_atom_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::AddStereoBond {
                site,
                ligands,
                attributes,
            } => GraphIrEdit::AddStereoBond(Box::new(GraphIrAddStereoBond {
                site: site.to_bond_handle(),
                ligands: ligands
                    .iter()
                    .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                    .collect(),
                attributes: attributes.to_rust(py),
            })),
            Self::RemoveStereoBonds { removes } => GraphIrEdit::RemoveStereoBonds {
                removes: removes
                    .0
                    .iter()
                    .map(|(id, site, ligands, attributes)| {
                        (
                            id.to_stereo_bond_handle(),
                            site.to_bond_handle(),
                            ligands
                                .iter()
                                .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                                .collect(),
                            attributes.to_rust(py),
                        )
                    })
                    .collect(),
            },
            Self::ReplaceStereoBondSite { id, old, new } => GraphIrEdit::ReplaceStereoBondSite {
                id: id.to_stereo_bond_handle(),
                old: old.to_bond_handle(),
                new: new.to_bond_handle(),
            },
            Self::ReplaceStereoBondLigands { id, old, new } => {
                GraphIrEdit::ReplaceStereoBondLigands {
                    id: id.to_stereo_bond_handle(),
                    old: old
                        .iter()
                        .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                        .collect(),
                    new: new
                        .iter()
                        .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                        .collect(),
                }
            }
            Self::ModifyStereoBondConfigurations { changes } => {
                GraphIrEdit::ModifyStereoBondConfigurations {
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_stereo_bond_handle(),
                                old.bind(py).borrow().to_rust(py),
                                new.bind(py).borrow().to_rust(py),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyStereoBondField { id, change } => GraphIrEdit::ModifyStereoBondField {
                id: id.to_stereo_bond_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyAtomConstraint { id, old, new } => GraphIrEdit::ModifyAtomConstraint {
                id: id.to_atom_handle(),
                old: old
                    .as_ref()
                    .map(|value| value.bind(py).borrow().to_rust(py)),
                new: new
                    .as_ref()
                    .map(|value| value.bind(py).borrow().to_rust(py)),
            },
            Self::ModifyBondConstraints { key, changes } => GraphIrEdit::ModifyBondConstraints {
                key: key.bind(py).borrow().to_rust(py),
                changes: changes
                    .0
                    .iter()
                    .map(|(id, old, new)| {
                        (
                            id.to_bond_handle(),
                            old.as_ref()
                                .map(|value| value.bind(py).borrow().to_rust(py)),
                            new.as_ref()
                                .map(|value| value.bind(py).borrow().to_rust(py)),
                        )
                    })
                    .collect(),
            },
            Self::ModifyBondConstraint { id, old, new } => GraphIrEdit::ModifyBondConstraint {
                id: id.to_bond_handle(),
                old: old
                    .as_ref()
                    .map(|value| value.bind(py).borrow().to_rust(py)),
                new: new
                    .as_ref()
                    .map(|value| value.bind(py).borrow().to_rust(py)),
            },
            Self::ModifyDativeBondConstraints { key, changes } => {
                GraphIrEdit::ModifyDativeBondConstraints {
                    key: key.bind(py).borrow().to_rust(py),
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_dative_bond_handle(),
                                old.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                                new.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyDativeBondConstraint { id, old, new } => {
                GraphIrEdit::ModifyDativeBondConstraint {
                    id: id.to_dative_bond_handle(),
                    old: old
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                    new: new
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                }
            }
            Self::ModifyAromaticSystemConstraints { key, changes } => {
                GraphIrEdit::ModifyAromaticSystemConstraints {
                    key: key.bind(py).borrow().to_rust(),
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_aromatic_system_handle(),
                                old.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                                new.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyAromaticSystemConstraint { id, old, new } => {
                GraphIrEdit::ModifyAromaticSystemConstraint {
                    id: id.to_aromatic_system_handle(),
                    old: old
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                    new: new
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                }
            }
            Self::ModifyMulticenterBondConstraints { key, changes } => {
                GraphIrEdit::ModifyMulticenterBondConstraints {
                    key: key.bind(py).borrow().to_rust(),
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_multicenter_bond_handle(),
                                old.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                                new.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyMulticenterBondConstraint { id, old, new } => {
                GraphIrEdit::ModifyMulticenterBondConstraint {
                    id: id.to_multicenter_bond_handle(),
                    old: old
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                    new: new
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                }
            }
            Self::ModifyNoncovalentBondConstraints { key, changes } => {
                GraphIrEdit::ModifyNoncovalentBondConstraints {
                    key: key.bind(py).borrow().to_rust(),
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_noncovalent_bond_handle(),
                                old.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                                new.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyNoncovalentBondConstraint { id, old, new } => {
                GraphIrEdit::ModifyNoncovalentBondConstraint {
                    id: id.to_noncovalent_bond_handle(),
                    old: old
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                    new: new
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                }
            }
            Self::ModifyStereoAtomConstraints { kind, key, changes } => {
                GraphIrEdit::ModifyStereoAtomConstraints {
                    kind: kind.map(StereoKind::to_rust),
                    key: key.bind(py).borrow().to_rust(py),
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_stereo_atom_handle(),
                                old.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                                new.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyStereoAtomConstraint { id, kind, old, new } => {
                GraphIrEdit::ModifyStereoAtomConstraint {
                    id: id.to_stereo_atom_handle(),
                    kind: kind.map(StereoKind::to_rust),
                    old: old
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                    new: new
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                }
            }
            Self::ModifyStereoBondConstraints { kind, key, changes } => {
                GraphIrEdit::ModifyStereoBondConstraints {
                    kind: kind.map(StereoKind::to_rust),
                    key: key.bind(py).borrow().to_rust(py),
                    changes: changes
                        .0
                        .iter()
                        .map(|(id, old, new)| {
                            (
                                id.to_stereo_bond_handle(),
                                old.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                                new.as_ref()
                                    .map(|value| value.bind(py).borrow().to_rust(py)),
                            )
                        })
                        .collect(),
                }
            }
            Self::ModifyStereoBondConstraint { id, kind, old, new } => {
                GraphIrEdit::ModifyStereoBondConstraint {
                    id: id.to_stereo_bond_handle(),
                    kind: kind.map(StereoKind::to_rust),
                    old: old
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                    new: new
                        .as_ref()
                        .map(|value| value.bind(py).borrow().to_rust(py)),
                }
            }
            Self::AddMoleculeConstraint { constraint } => GraphIrEdit::AddMoleculeConstraint {
                constraint: Box::new(constraint.bind(py).borrow().to_rust().clone()),
            },
            Self::RemoveMoleculeConstraint { constraint } => {
                GraphIrEdit::RemoveMoleculeConstraint {
                    constraint: Box::new(constraint.bind(py).borrow().to_rust().clone()),
                }
            }
        }
    }
}

fn resolve_edit_index(len: usize, index: isize) -> PyResult<usize> {
    let resolved = if index < 0 {
        index + len as isize
    } else {
        index
    };
    if resolved < 0 || resolved as usize >= len {
        Err(PyIndexError::new_err("edit index out of range"))
    } else {
        Ok(resolved as usize)
    }
}

/// A lazy iterator yielding independent copies of the batch's initial entries.
///
/// Retains the batch until dropped, but excludes entries appended after creation.
#[pyclass]
pub(crate) struct EditIter {
    owner: Py<Edits>,
    position: usize,
    end: usize,
}

#[pymethods]
impl EditIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<Edit>>> {
        let owner = self.owner.try_borrow(py)?;
        let edits = owner.value.as_ref().ok_or_else(|| {
            InvalidatedViewError::new_err(
                "Edit iterator is invalid because Edits has been consumed",
            )
        })?;
        if self.position == self.end {
            return Ok(None);
        }
        let entry = into_py_variant(py, Edit::from_rust(py, &edits.as_slice()[self.position])?)?;
        self.position += 1;
        Ok(Some(entry))
    }
}

/// An ordered, append-only batch of host-specific molecule edits.
///
/// Passing the batch to an operation that consumes Edits moves its contents.
/// Subsequent access raises ConsumedError; its iterators raise InvalidatedViewError.
/// Entries already returned by indexing or iteration are independent values.
#[pyclass]
#[derive(Debug)]
pub struct Edits {
    value: Option<GraphIrEdits>,
}

#[pymethods]
impl Edits {
    #[new]
    #[pyo3(signature = (entries=Vec::new()))]
    fn new(py: Python<'_>, entries: Vec<Py<Edit>>) -> Self {
        Self::from_rust(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py))
                .collect(),
        )
    }

    #[staticmethod]
    #[pyo3(signature = (text, *, defaults=None))]
    fn parse(text: &str, defaults: Option<MoleculeDefaults>) -> PyResult<Self> {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new);
        let edits = GraphIrEditsDsl::from_str(text)
            .map_err(parse_error)?
            .into_ir(defaults.to_rust());
        Ok(Self::from_rust(edits))
    }

    /// Render edits, omitting whole stereo constraint operations with no kind context.
    #[pyo3(signature = (*, defaults=None))]
    fn render(&self, defaults: Option<MoleculeDefaults>) -> PyResult<String> {
        let edits = self.to_rust()?;
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new);
        Ok(GraphIrEditsDsl::from_ir(edits, defaults.to_rust()).to_string())
    }

    fn __eq__(&self, other: &Self) -> PyResult<bool> {
        Ok(self.to_rust()? == other.to_rust()?)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let edits = self.to_rust()?;
        let mut parts = Vec::with_capacity(edits.len());
        for entry in edits.iter() {
            let value = into_py_variant(py, Edit::from_rust(py, entry)?)?;
            parts.push(value.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("Edits([{}])", parts.join(", ")))
    }

    /// Append one detached raw edit and account for every entity it creates.
    fn append(&mut self, py: Python<'_>, edit: Py<Edit>) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.push(edit.bind(py).borrow().to_rust(py));
        Ok(())
    }

    fn __len__(&self) -> PyResult<usize> {
        Ok(self.to_rust()?.len())
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<Edit> {
        let edits = self.to_rust()?;
        let index = resolve_edit_index(edits.len(), index)?;
        Edit::from_rust(py, &edits.as_slice()[index])
    }

    /// Lazily copy entries present now; later appends are excluded.
    fn __iter__(slf: Py<Self>, py: Python<'_>) -> PyResult<EditIter> {
        let end = slf.try_borrow(py)?.to_rust()?.len();
        Ok(EditIter {
            owner: slf,
            position: 0,
            end,
        })
    }

    fn add_atom(&mut self, py: Python<'_>, attributes: Py<AtomForm>) -> PyResult<New> {
        let edits = self.to_rust_mut()?;
        Ok(New::from_rust(GraphIrEntityHandle::Atom(
            edits.add_atom(attributes.bind(py).borrow().to_rust().clone()),
        )))
    }

    fn add_atoms(&mut self, py: Python<'_>, atoms: Vec<Py<AtomForm>>) -> PyResult<Vec<New>> {
        let edits = self.to_rust_mut()?;
        Ok(edits
            .add_atoms(
                atoms
                    .into_iter()
                    .map(|attributes| attributes.bind(py).borrow().to_rust().clone()),
            )
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::Atom(handle)))
            .collect())
    }

    fn add_bond(
        &mut self,
        py: Python<'_>,
        first: HandleLike,
        second: HandleLike,
        attributes: Py<BondForm>,
    ) -> PyResult<New> {
        let edits = self.to_rust_mut()?;
        Ok(New::from_rust(GraphIrEntityHandle::Bond(edits.add_bond(
            first.to_atom_handle(),
            second.to_atom_handle(),
            attributes.bind(py).borrow().to_rust().clone(),
        ))))
    }

    fn add_bonds(&mut self, py: Python<'_>, bonds: Vec<BondAddition>) -> PyResult<Vec<New>> {
        let edits = self.to_rust_mut()?;
        Ok(edits
            .add_bonds(
                bonds
                    .into_iter()
                    .map(|((first, second), attributes)| GraphIrAddBond {
                        endpoints: [first.to_atom_handle(), second.to_atom_handle()],
                        attributes: attributes.bind(py).borrow().to_rust().clone(),
                    }),
            )
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::Bond(handle)))
            .collect())
    }

    fn add_dative_bond(
        &mut self,
        py: Python<'_>,
        donors: Vec<HandleLike>,
        acceptor: HandleLike,
        attributes: Py<DativeBondForm>,
    ) -> PyResult<New> {
        let edits = self.to_rust_mut()?;
        Ok(New::from_rust(GraphIrEntityHandle::DativeBond(
            edits.add_dative_bond(
                donors.iter().map(HandleLike::to_atom_handle).collect(),
                acceptor.to_atom_handle(),
                attributes.bind(py).borrow().to_rust().clone(),
            ),
        )))
    }

    fn add_dative_bonds(
        &mut self,
        py: Python<'_>,
        bonds: Vec<DativeBondAddition>,
    ) -> PyResult<Vec<New>> {
        let edits = self.to_rust_mut()?;
        Ok(edits
            .add_dative_bonds(bonds.into_iter().map(|(donors, acceptor, attributes)| {
                (
                    donors.iter().map(HandleLike::to_atom_handle).collect(),
                    acceptor.to_atom_handle(),
                    attributes.bind(py).borrow().to_rust().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::DativeBond(handle)))
            .collect())
    }

    fn add_aromatic_system(
        &mut self,
        py: Python<'_>,
        atoms: Vec<HandleLike>,
        attributes: Py<AromaticSystemForm>,
    ) -> PyResult<New> {
        let edits = self.to_rust_mut()?;
        Ok(New::from_rust(GraphIrEntityHandle::AromaticSystem(
            edits.add_aromatic_system(
                atoms.iter().map(HandleLike::to_atom_handle).collect(),
                attributes.bind(py).borrow().to_rust().clone(),
            ),
        )))
    }

    fn add_aromatic_systems(
        &mut self,
        py: Python<'_>,
        systems: Vec<AromaticSystemAddition>,
    ) -> PyResult<Vec<New>> {
        let edits = self.to_rust_mut()?;
        Ok(edits
            .add_aromatic_systems(systems.into_iter().map(|(atoms, attributes)| {
                (
                    atoms.iter().map(HandleLike::to_atom_handle).collect(),
                    attributes.bind(py).borrow().to_rust().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::AromaticSystem(handle)))
            .collect())
    }

    fn add_multicenter_bond(
        &mut self,
        py: Python<'_>,
        atoms: Vec<HandleLike>,
        attributes: Py<MulticenterBondForm>,
    ) -> PyResult<New> {
        let edits = self.to_rust_mut()?;
        Ok(New::from_rust(GraphIrEntityHandle::MulticenterBond(
            edits.add_multicenter_bond(
                atoms.iter().map(HandleLike::to_atom_handle).collect(),
                attributes.bind(py).borrow().to_rust().clone(),
            ),
        )))
    }

    fn add_multicenter_bonds(
        &mut self,
        py: Python<'_>,
        bonds: Vec<MulticenterBondAddition>,
    ) -> PyResult<Vec<New>> {
        let edits = self.to_rust_mut()?;
        Ok(edits
            .add_multicenter_bonds(bonds.into_iter().map(|(atoms, attributes)| {
                (
                    atoms.iter().map(HandleLike::to_atom_handle).collect(),
                    attributes.bind(py).borrow().to_rust().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::MulticenterBond(handle)))
            .collect())
    }

    fn add_noncovalent_bond(
        &mut self,
        py: Python<'_>,
        atoms: (HandleLike, HandleLike),
        attributes: Py<NoncovalentBondForm>,
    ) -> PyResult<New> {
        let edits = self.to_rust_mut()?;
        Ok(New::from_rust(GraphIrEntityHandle::NoncovalentBond(
            edits.add_noncovalent_bond(
                [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                attributes.bind(py).borrow().to_rust().clone(),
            ),
        )))
    }

    fn add_noncovalent_bonds(
        &mut self,
        py: Python<'_>,
        bonds: Vec<NoncovalentBondAddition>,
    ) -> PyResult<Vec<New>> {
        let edits = self.to_rust_mut()?;
        Ok(edits
            .add_noncovalent_bonds(bonds.into_iter().map(|(atoms, attributes)| {
                (
                    [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                    attributes.bind(py).borrow().to_rust().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::NoncovalentBond(handle)))
            .collect())
    }

    fn add_stereo_atom(
        &mut self,
        py: Python<'_>,
        site: HandleLike,
        ligands: Vec<StereoLigandInput>,
        attributes: Py<StereoAtomForm>,
    ) -> PyResult<New> {
        let edits = self.to_rust_mut()?;
        Ok(New::from_rust(GraphIrEntityHandle::StereoAtom(
            edits.add_stereo_atom(
                site.to_atom_handle(),
                ligands
                    .iter()
                    .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                    .collect(),
                attributes.bind(py).borrow().to_rust().clone(),
            ),
        )))
    }

    fn add_stereo_atoms(
        &mut self,
        py: Python<'_>,
        atoms: Vec<StereoAtomAddition>,
    ) -> PyResult<Vec<New>> {
        let edits = self.to_rust_mut()?;
        Ok(edits
            .add_stereo_atoms(atoms.into_iter().map(|(site, ligands, attributes)| {
                (
                    site.to_atom_handle(),
                    ligands
                        .iter()
                        .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                        .collect(),
                    attributes.bind(py).borrow().to_rust().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::StereoAtom(handle)))
            .collect())
    }

    fn add_stereo_bond(
        &mut self,
        py: Python<'_>,
        site: HandleLike,
        ligands: Vec<StereoLigandInput>,
        attributes: Py<StereoBondForm>,
    ) -> PyResult<New> {
        let edits = self.to_rust_mut()?;
        Ok(New::from_rust(GraphIrEntityHandle::StereoBond(
            edits.add_stereo_bond(
                site.to_bond_handle(),
                ligands
                    .iter()
                    .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                    .collect(),
                attributes.bind(py).borrow().to_rust().clone(),
            ),
        )))
    }

    fn add_stereo_bonds(
        &mut self,
        py: Python<'_>,
        bonds: Vec<StereoBondAddition>,
    ) -> PyResult<Vec<New>> {
        let edits = self.to_rust_mut()?;
        Ok(edits
            .add_stereo_bonds(bonds.into_iter().map(|(site, ligands, attributes)| {
                (
                    site.to_bond_handle(),
                    ligands
                        .iter()
                        .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                        .collect(),
                    attributes.bind(py).borrow().to_rust().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::StereoBond(handle)))
            .collect())
    }

    fn remove_topology(&mut self, atoms: Vec<HandleLike>, bonds: Vec<HandleLike>) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.remove_topology(
            atoms.iter().map(HandleLike::to_atom_handle).collect(),
            bonds.iter().map(HandleLike::to_bond_handle).collect(),
        );
        Ok(())
    }

    fn remove_dative_bonds(
        &mut self,
        py: Python<'_>,
        removes: Vec<DativeBondRemoval>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.remove_dative_bonds(
            removes
                .into_iter()
                .map(|(id, donors, acceptor, attributes)| {
                    (
                        id.to_dative_bond_handle(),
                        donors.iter().map(HandleLike::to_atom_handle).collect(),
                        acceptor.to_atom_handle(),
                        attributes.bind(py).borrow().to_rust().clone(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn remove_aromatic_systems(
        &mut self,
        py: Python<'_>,
        removes: Vec<AromaticSystemRemoval>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.remove_aromatic_systems(
            removes
                .into_iter()
                .map(|(id, atoms, attributes)| {
                    (
                        id.to_aromatic_system_handle(),
                        atoms.iter().map(HandleLike::to_atom_handle).collect(),
                        attributes.bind(py).borrow().to_rust().clone(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn remove_multicenter_bonds(
        &mut self,
        py: Python<'_>,
        removes: Vec<MulticenterBondRemoval>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.remove_multicenter_bonds(
            removes
                .into_iter()
                .map(|(id, atoms, attributes)| {
                    (
                        id.to_multicenter_bond_handle(),
                        atoms.iter().map(HandleLike::to_atom_handle).collect(),
                        attributes.bind(py).borrow().to_rust().clone(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn remove_noncovalent_bonds(
        &mut self,
        py: Python<'_>,
        removes: Vec<NoncovalentBondRemoval>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.remove_noncovalent_bonds(
            removes
                .into_iter()
                .map(|(id, atoms, attributes)| {
                    (
                        id.to_noncovalent_bond_handle(),
                        [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                        attributes.bind(py).borrow().to_rust().clone(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn remove_stereo_atoms(
        &mut self,
        py: Python<'_>,
        removes: Vec<StereoAtomRemoval>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.remove_stereo_atoms(
            removes
                .into_iter()
                .map(|(id, site, ligands, attributes)| {
                    (
                        id.to_stereo_atom_handle(),
                        site.to_atom_handle(),
                        ligands
                            .iter()
                            .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                            .collect(),
                        attributes.bind(py).borrow().to_rust().clone(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn remove_stereo_bonds(
        &mut self,
        py: Python<'_>,
        removes: Vec<StereoBondRemoval>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.remove_stereo_bonds(
            removes
                .into_iter()
                .map(|(id, site, ligands, attributes)| {
                    (
                        id.to_stereo_bond_handle(),
                        site.to_bond_handle(),
                        ligands
                            .iter()
                            .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                            .collect(),
                        attributes.bind(py).borrow().to_rust().clone(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn add_molecule_constraint(
        &mut self,
        py: Python<'_>,
        constraint: Py<ConstraintEdit>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.add_molecule_constraint(constraint.bind(py).borrow().to_rust().clone());
        Ok(())
    }

    fn remove_molecule_constraint(
        &mut self,
        py: Python<'_>,
        constraint: Py<ConstraintEdit>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.remove_molecule_constraint(constraint.bind(py).borrow().to_rust().clone());
        Ok(())
    }

    fn modify_atom_elements(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<ElementForm>, Py<ElementForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_atom_elements(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_atom_handle(),
                        old.bind(py).borrow().to_rust(),
                        new.bind(py).borrow().to_rust(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_atom_isotope_masses(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<IsotopeMassForm>, Py<IsotopeMassForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_atom_isotope_masses(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_atom_handle(),
                        old.bind(py).borrow().to_rust(),
                        new.bind(py).borrow().to_rust(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_atom_charges(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_atom_charges(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_atom_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_atom_implicit_hydrogens(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_atom_implicit_hydrogens(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_atom_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_atom_lone_pairs(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_atom_lone_pairs(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_atom_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_atom_unpaired_electrons(
        &mut self,
        py: Python<'_>,
        changes: Vec<(
            HandleLike,
            Py<UnpairedElectronsForm>,
            Py<UnpairedElectronsForm>,
        )>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_atom_unpaired_electrons(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_atom_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_atom_constraints(
        &mut self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
        changes: Vec<AtomConstraintChange>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_atom_constraints(
            key.bind(py).borrow().to_rust(py),
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_atom_handle(),
                        old.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                        new.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_bond_orders(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_bond_orders(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_bond_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_bond_charges(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_bond_charges(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_bond_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_bond_unpaired_electrons(
        &mut self,
        py: Python<'_>,
        changes: Vec<(
            HandleLike,
            Py<UnpairedElectronsForm>,
            Py<UnpairedElectronsForm>,
        )>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_bond_unpaired_electrons(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_bond_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_bond_constraints(
        &mut self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
        changes: Vec<BondConstraintChange>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_bond_constraints(
            key.bind(py).borrow().to_rust(py),
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_bond_handle(),
                        old.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                        new.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_dative_bond_orders(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_dative_bond_orders(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_dative_bond_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_dative_bond_constraints(
        &mut self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
        changes: Vec<DativeBondConstraintChange>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_dative_bond_constraints(
            key.bind(py).borrow().to_rust(py),
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_dative_bond_handle(),
                        old.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                        new.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_aromatic_system_electrons(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<ElectronCountsForm>, Py<ElectronCountsForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_aromatic_system_electrons(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_aromatic_system_handle(),
                        old.bind(py).borrow().to_rust(),
                        new.bind(py).borrow().to_rust(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_aromatic_system_charges(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_aromatic_system_charges(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_aromatic_system_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_aromatic_system_unpaired_electrons(
        &mut self,
        py: Python<'_>,
        changes: Vec<(
            HandleLike,
            Py<UnpairedElectronsForm>,
            Py<UnpairedElectronsForm>,
        )>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_aromatic_system_unpaired_electrons(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_aromatic_system_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_aromatic_system_constraints(
        &mut self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
        changes: Vec<AromaticSystemConstraintChange>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_aromatic_system_constraints(
            key.bind(py).borrow().to_rust(),
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_aromatic_system_handle(),
                        old.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                        new.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_multicenter_bond_electrons(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<ElectronCountsForm>, Py<ElectronCountsForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_multicenter_bond_electrons(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_multicenter_bond_handle(),
                        old.bind(py).borrow().to_rust(),
                        new.bind(py).borrow().to_rust(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_multicenter_bond_charges(
        &mut self,
        py: Python<'_>,
        changes: Vec<(HandleLike, Py<NumForm>, Py<NumForm>)>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_multicenter_bond_charges(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_multicenter_bond_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_multicenter_bond_unpaired_electrons(
        &mut self,
        py: Python<'_>,
        changes: Vec<(
            HandleLike,
            Py<UnpairedElectronsForm>,
            Py<UnpairedElectronsForm>,
        )>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_multicenter_bond_unpaired_electrons(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_multicenter_bond_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_multicenter_bond_constraints(
        &mut self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
        changes: Vec<MulticenterBondConstraintChange>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_multicenter_bond_constraints(
            key.bind(py).borrow().to_rust(),
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_multicenter_bond_handle(),
                        old.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                        new.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_noncovalent_bond_kinds(
        &mut self,
        py: Python<'_>,
        changes: Vec<(
            HandleLike,
            Py<NoncovalentBondKindForm>,
            Py<NoncovalentBondKindForm>,
        )>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_noncovalent_bond_kinds(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_noncovalent_bond_handle(),
                        old.bind(py).borrow().to_rust(),
                        new.bind(py).borrow().to_rust(),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_noncovalent_bond_constraints(
        &mut self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
        changes: Vec<NoncovalentBondConstraintChange>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_noncovalent_bond_constraints(
            key.bind(py).borrow().to_rust(),
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_noncovalent_bond_handle(),
                        old.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                        new.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_stereo_atom_configurations(
        &mut self,
        py: Python<'_>,
        changes: Vec<(
            HandleLike,
            Py<StereoConfigurationForm>,
            Py<StereoConfigurationForm>,
        )>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_stereo_atom_configurations(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_stereo_atom_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    #[pyo3(signature = (kind, key, changes))]
    fn modify_stereo_atom_constraints(
        &mut self,
        py: Python<'_>,
        kind: Option<StereoKind>,
        key: Py<StereoAtomConstraintKey>,
        changes: Vec<StereoAtomConstraintChange>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_stereo_atom_constraints(
            kind.map(StereoKind::to_rust),
            key.bind(py).borrow().to_rust(py),
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_stereo_atom_handle(),
                        old.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                        new.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn modify_stereo_bond_configurations(
        &mut self,
        py: Python<'_>,
        changes: Vec<(
            HandleLike,
            Py<StereoConfigurationForm>,
            Py<StereoConfigurationForm>,
        )>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_stereo_bond_configurations(
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_stereo_bond_handle(),
                        old.bind(py).borrow().to_rust(py),
                        new.bind(py).borrow().to_rust(py),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    #[pyo3(signature = (kind, key, changes))]
    fn modify_stereo_bond_constraints(
        &mut self,
        py: Python<'_>,
        kind: Option<StereoKind>,
        key: Py<StereoBondConstraintKey>,
        changes: Vec<StereoBondConstraintChange>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        edits.modify_stereo_bond_constraints(
            kind.map(StereoKind::to_rust),
            key.bind(py).borrow().to_rust(py),
            changes
                .into_iter()
                .map(|(id, old, new)| {
                    (
                        id.to_stereo_bond_handle(),
                        old.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                        new.as_ref()
                            .map(|value| value.bind(py).borrow().to_rust(py)),
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn update_atom(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<AtomForm>,
        update: Py<AtomUpdate>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow();
        edits.update_atom(id.to_atom_handle(), current.to_rust(), update.to_rust());
        Ok(())
    }

    fn update_bond(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<BondForm>,
        update: Py<BondUpdate>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow();
        edits.update_bond(id.to_bond_handle(), current.to_rust(), update.to_rust());
        Ok(())
    }

    fn update_dative_bond(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<DativeBondForm>,
        update: Py<DativeBondUpdate>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow();
        edits.update_dative_bond(
            id.to_dative_bond_handle(),
            current.to_rust(),
            update.to_rust(),
        );
        Ok(())
    }

    fn update_aromatic_system(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<AromaticSystemForm>,
        update: Py<AromaticSystemUpdate>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow();
        edits.update_aromatic_system(
            id.to_aromatic_system_handle(),
            current.to_rust(),
            update.to_rust(),
        );
        Ok(())
    }

    fn update_multicenter_bond(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<MulticenterBondForm>,
        update: Py<MulticenterBondUpdate>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow();
        edits.update_multicenter_bond(
            id.to_multicenter_bond_handle(),
            current.to_rust(),
            update.to_rust(),
        );
        Ok(())
    }

    fn update_noncovalent_bond(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<NoncovalentBondForm>,
        update: Py<NoncovalentBondUpdate>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow();
        edits.update_noncovalent_bond(
            id.to_noncovalent_bond_handle(),
            current.to_rust(),
            update.to_rust(),
        );
        Ok(())
    }

    fn update_stereo_atom(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<StereoAtomForm>,
        update: Py<StereoAtomUpdate>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow();
        edits.update_stereo_atom(
            id.to_stereo_atom_handle(),
            current.to_rust(),
            update.to_rust(),
        );
        Ok(())
    }

    fn update_stereo_bond(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<StereoBondForm>,
        update: Py<StereoBondUpdate>,
    ) -> PyResult<()> {
        let edits = self.to_rust_mut()?;
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow();
        edits.update_stereo_bond(
            id.to_stereo_bond_handle(),
            current.to_rust(),
            update.to_rust(),
        );
        Ok(())
    }
}

impl Edits {
    pub(crate) fn from_rust(edits: GraphIrEdits) -> Self {
        Self { value: Some(edits) }
    }

    pub(crate) fn to_rust(&self) -> PyResult<&GraphIrEdits> {
        self.value
            .as_ref()
            .ok_or_else(|| ConsumedError::new_err("Edits has been consumed"))
    }

    pub(crate) fn to_rust_mut(&mut self) -> PyResult<&mut GraphIrEdits> {
        self.value
            .as_mut()
            .ok_or_else(|| ConsumedError::new_err("Edits has been consumed"))
    }

    pub(crate) fn take(&mut self) -> PyResult<GraphIrEdits> {
        self.value
            .take()
            .ok_or_else(|| ConsumedError::new_err("Edits has been consumed"))
    }
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::PyRuntimeError;
    use rstest::rstest;
    use umol_graph_ir::ir::{
        AromaticSystemFieldChange as GraphIrAromaticSystemFieldChange,
        AromaticSystemForm as GraphIrAromaticSystemForm, AtomFieldChange as GraphIrAtomFieldChange,
        AtomForm as GraphIrAtomForm, BondFieldChange as GraphIrBondFieldChange,
        BondForm as GraphIrBondForm, Constraint as GraphIrConstraint,
        DativeBondFieldChange as GraphIrDativeBondFieldChange,
        DativeBondForm as GraphIrDativeBondForm, MoleculeConstraint as GraphIrMoleculeConstraint,
        MulticenterBondFieldChange as GraphIrMulticenterBondFieldChange,
        MulticenterBondForm as GraphIrMulticenterBondForm,
        NoncovalentBondFieldChange as GraphIrNoncovalentBondFieldChange,
        NoncovalentBondForm as GraphIrNoncovalentBondForm,
        NoncovalentBondKind as GraphIrNoncovalentBondKind,
        NoncovalentBondKindForm as GraphIrNoncovalentBondKindForm, NumForm as GraphIrNumForm,
        StereoAtomFieldChange as GraphIrStereoAtomFieldChange,
        StereoAtomForm as GraphIrStereoAtomForm,
        StereoBondFieldChange as GraphIrStereoBondFieldChange,
        StereoBondForm as GraphIrStereoBondForm,
        StereoConfigurationForm as GraphIrStereoConfigurationForm, StereoKind as GraphIrStereoKind,
        StereoLigandKind as GraphIrStereoLigandKind,
    };

    use super::*;

    #[rstest]
    #[case::id(HandleLike::Id(7), GraphIrAtomHandle::Id(GraphIrAtomId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), GraphIrAtomHandle::New(7))]
    fn test_handle_like_to_atom_handle(
        #[case] input: HandleLike,
        #[case] expected: GraphIrAtomHandle,
    ) {
        assert_eq!(input.to_atom_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), GraphIrBondHandle::Id(GraphIrBondId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), GraphIrBondHandle::New(7))]
    fn test_handle_like_to_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: GraphIrBondHandle,
    ) {
        assert_eq!(input.to_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), GraphIrDativeBondHandle::Id(GraphIrDativeBondId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), GraphIrDativeBondHandle::New(7))]
    fn test_handle_like_to_dative_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: GraphIrDativeBondHandle,
    ) {
        assert_eq!(input.to_dative_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(
        HandleLike::Id(7),
        GraphIrAromaticSystemHandle::Id(GraphIrAromaticSystemId(7))
    )]
    #[case::new(HandleLike::New(New { index: 7 }), GraphIrAromaticSystemHandle::New(7))]
    fn test_handle_like_to_aromatic_system_handle(
        #[case] input: HandleLike,
        #[case] expected: GraphIrAromaticSystemHandle,
    ) {
        assert_eq!(input.to_aromatic_system_handle(), expected);
    }

    #[rstest]
    #[case::id(
        HandleLike::Id(7),
        GraphIrMulticenterBondHandle::Id(GraphIrMulticenterBondId(7))
    )]
    #[case::new(HandleLike::New(New { index: 7 }), GraphIrMulticenterBondHandle::New(7))]
    fn test_handle_like_to_multicenter_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: GraphIrMulticenterBondHandle,
    ) {
        assert_eq!(input.to_multicenter_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(
        HandleLike::Id(7),
        GraphIrNoncovalentBondHandle::Id(GraphIrNoncovalentBondId(7))
    )]
    #[case::new(HandleLike::New(New { index: 7 }), GraphIrNoncovalentBondHandle::New(7))]
    fn test_handle_like_to_noncovalent_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: GraphIrNoncovalentBondHandle,
    ) {
        assert_eq!(input.to_noncovalent_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), GraphIrStereoAtomHandle::Id(GraphIrStereoAtomId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), GraphIrStereoAtomHandle::New(7))]
    fn test_handle_like_to_stereo_atom_handle(
        #[case] input: HandleLike,
        #[case] expected: GraphIrStereoAtomHandle,
    ) {
        assert_eq!(input.to_stereo_atom_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), GraphIrStereoBondHandle::Id(GraphIrStereoBondId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), GraphIrStereoBondHandle::New(7))]
    fn test_handle_like_to_stereo_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: GraphIrStereoBondHandle,
    ) {
        assert_eq!(input.to_stereo_bond_handle(), expected);
    }

    #[rstest]
    #[case::identity(GraphIrConstraintEdit::from(GraphIrConstraint::Molecule(
        GraphIrMoleculeConstraint::Connected { atoms: None },
    )))]
    fn test_constraint_edit_roundtrip(#[case] expected: GraphIrConstraintEdit) {
        let constraint = ConstraintEdit::from_rust(&expected);

        assert_eq!(constraint.to_rust(), &expected);
    }

    #[rstest]
    #[case::inventory(vec![
        GraphIrEdit::AddAtoms { atoms: vec![GraphIrAtomForm::default()] },
        GraphIrEdit::AddBonds { bonds: vec![GraphIrAddBond {
            endpoints: [GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrAtomHandle::New(0)],
            attributes: GraphIrBondForm::default(),
        }] },
        GraphIrEdit::RemoveTopology {
            atoms: vec![GraphIrAtomHandle::New(0)],
            bonds: vec![GraphIrBondHandle::Id(GraphIrBondId(1))],
        },
        GraphIrEdit::ModifyAtomField {
            id: GraphIrAtomHandle::Id(GraphIrAtomId(0)),
            change: GraphIrAtomFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(1),
            },
        },
        GraphIrEdit::ModifyBondField {
            id: GraphIrBondHandle::New(0),
            change: GraphIrBondFieldChange::Order {
                old: GraphIrNumForm::Lit(1),
                new: GraphIrNumForm::Lit(2),
            },
        },
        GraphIrEdit::AddDativeBond {
            donors: vec![GraphIrAtomHandle::Id(GraphIrAtomId(0))],
            acceptor: GraphIrAtomHandle::New(0),
            attributes: GraphIrDativeBondForm::default(),
        },
        GraphIrEdit::RemoveDativeBonds { removes: vec![(
            GraphIrDativeBondHandle::New(0),
            vec![],
            GraphIrAtomHandle::Id(GraphIrAtomId(0)),
            GraphIrDativeBondForm::default(),
        )] },
        GraphIrEdit::ReplaceDativeBondDonors {
            id: GraphIrDativeBondHandle::New(0),
            old: vec![GraphIrAtomHandle::Id(GraphIrAtomId(0))],
            new: vec![GraphIrAtomHandle::New(1), GraphIrAtomHandle::Id(GraphIrAtomId(2))],
        },
        GraphIrEdit::ReplaceDativeBondAcceptor {
            id: GraphIrDativeBondHandle::New(0),
            old: GraphIrAtomHandle::Id(GraphIrAtomId(0)),
            new: GraphIrAtomHandle::New(1),
        },
        GraphIrEdit::ModifyDativeBondField {
            id: GraphIrDativeBondHandle::Id(GraphIrDativeBondId(0)),
            change: GraphIrDativeBondFieldChange::Order {
                old: GraphIrNumForm::Lit(1),
                new: GraphIrNumForm::Lit(2),
            },
        },
        GraphIrEdit::AddAromaticSystem(Box::new(GraphIrAddAromaticSystem {
            atoms: vec![GraphIrAtomHandle::New(0)],
            attributes: GraphIrAromaticSystemForm::default(),
        })),
        GraphIrEdit::RemoveAromaticSystems { removes: vec![(
            GraphIrAromaticSystemHandle::Id(GraphIrAromaticSystemId(0)),
            vec![GraphIrAtomHandle::New(0)],
            GraphIrAromaticSystemForm::default(),
        )] },
        GraphIrEdit::ReplaceAromaticSystemAtoms {
            id: GraphIrAromaticSystemHandle::New(0),
            old: vec![GraphIrAtomHandle::Id(GraphIrAtomId(0))],
            new: vec![GraphIrAtomHandle::New(1), GraphIrAtomHandle::Id(GraphIrAtomId(2))],
        },
        GraphIrEdit::ModifyAromaticSystemField {
            id: GraphIrAromaticSystemHandle::New(0),
            change: GraphIrAromaticSystemFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(-1),
            },
        },
        GraphIrEdit::AddMulticenterBond(Box::new(GraphIrAddMulticenterBond {
            atoms: vec![GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrAtomHandle::New(0)],
            attributes: GraphIrMulticenterBondForm::default(),
        })),
        GraphIrEdit::RemoveMulticenterBonds { removes: vec![(
            GraphIrMulticenterBondHandle::New(0),
            vec![GraphIrAtomHandle::Id(GraphIrAtomId(0))],
            GraphIrMulticenterBondForm::default(),
        )] },
        GraphIrEdit::ReplaceMulticenterBondAtoms {
            id: GraphIrMulticenterBondHandle::New(0),
            old: vec![GraphIrAtomHandle::Id(GraphIrAtomId(0))],
            new: vec![GraphIrAtomHandle::New(1), GraphIrAtomHandle::Id(GraphIrAtomId(2))],
        },
        GraphIrEdit::ModifyMulticenterBondField {
            id: GraphIrMulticenterBondHandle::Id(GraphIrMulticenterBondId(0)),
            change: GraphIrMulticenterBondFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(1),
            },
        },
        GraphIrEdit::AddNoncovalentBond {
            atoms: [GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrAtomHandle::New(0)],
            attributes: GraphIrNoncovalentBondForm::default(),
        },
        GraphIrEdit::RemoveNoncovalentBonds { removes: vec![(
            GraphIrNoncovalentBondHandle::New(0),
            [GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrAtomHandle::New(0)],
            GraphIrNoncovalentBondForm::default(),
        )] },
        GraphIrEdit::ReplaceNoncovalentBondAtoms {
            id: GraphIrNoncovalentBondHandle::New(0),
            old: [GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrAtomHandle::New(0)],
            new: [GraphIrAtomHandle::New(1), GraphIrAtomHandle::Id(GraphIrAtomId(2))],
        },
        GraphIrEdit::ModifyNoncovalentBondField {
            id: GraphIrNoncovalentBondHandle::Id(GraphIrNoncovalentBondId(0)),
            change: GraphIrNoncovalentBondFieldChange::Kind {
                old: GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond),
                new: GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::Ionic),
            },
        },
        GraphIrEdit::AddStereoAtom(Box::new(GraphIrAddStereoAtom {
            site: GraphIrAtomHandle::New(0),
            ligands: vec![(GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrStereoLigandKind::Atom)],
            attributes: GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, 0_u32),
        })),
        GraphIrEdit::RemoveStereoAtoms { removes: vec![(
            GraphIrStereoAtomHandle::Id(GraphIrStereoAtomId(0)),
            GraphIrAtomHandle::New(0),
            vec![(GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrStereoLigandKind::ImplicitHydrogen)],
            GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, 1_u32),
        )] },
        GraphIrEdit::ReplaceStereoAtomSite {
            id: GraphIrStereoAtomHandle::New(0),
            old: GraphIrAtomHandle::Id(GraphIrAtomId(0)),
            new: GraphIrAtomHandle::New(1),
        },
        GraphIrEdit::ReplaceStereoAtomLigands {
            id: GraphIrStereoAtomHandle::New(0),
            old: vec![(GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrStereoLigandKind::Atom)],
            new: vec![(GraphIrAtomHandle::New(1), GraphIrStereoLigandKind::ImplicitHydrogen), (GraphIrAtomHandle::New(0), GraphIrStereoLigandKind::LonePair)],
        },
        GraphIrEdit::ModifyStereoAtomField {
            id: GraphIrStereoAtomHandle::New(0),
            change: GraphIrStereoAtomFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::kinded(GraphIrStereoKind::Tetrahedral, 0_u32),
                new: GraphIrStereoConfigurationForm::kinded(GraphIrStereoKind::Tetrahedral, 1_u32),
            },
        },
        GraphIrEdit::AddStereoBond(Box::new(GraphIrAddStereoBond {
            site: GraphIrBondHandle::Id(GraphIrBondId(0)),
            ligands: vec![(GraphIrAtomHandle::New(0), GraphIrStereoLigandKind::LonePair)],
            attributes: GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, 0_u32),
        })),
        GraphIrEdit::RemoveStereoBonds { removes: vec![(
            GraphIrStereoBondHandle::New(0),
            GraphIrBondHandle::Id(GraphIrBondId(0)),
            vec![(GraphIrAtomHandle::New(0), GraphIrStereoLigandKind::Atom)],
            GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, 1_u32),
        )] },
        GraphIrEdit::ReplaceStereoBondSite {
            id: GraphIrStereoBondHandle::New(0),
            old: GraphIrBondHandle::Id(GraphIrBondId(0)),
            new: GraphIrBondHandle::New(1),
        },
        GraphIrEdit::ReplaceStereoBondLigands {
            id: GraphIrStereoBondHandle::New(0),
            old: vec![(GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrStereoLigandKind::Atom)],
            new: vec![(GraphIrAtomHandle::New(1), GraphIrStereoLigandKind::ImplicitHydrogen), (GraphIrAtomHandle::New(0), GraphIrStereoLigandKind::LonePair)],
        },
        GraphIrEdit::ModifyStereoBondField {
            id: GraphIrStereoBondHandle::Id(GraphIrStereoBondId(0)),
            change: GraphIrStereoBondFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::kinded(GraphIrStereoKind::CisTrans, 0_u32),
                new: GraphIrStereoConfigurationForm::kinded(GraphIrStereoKind::CisTrans, 1_u32),
            },
        },
        GraphIrEdit::ModifyAtomConstraint {
            id: GraphIrAtomHandle::Id(GraphIrAtomId(0)), old: None, new: None,
        },
        GraphIrEdit::ModifyBondConstraint {
            id: GraphIrBondHandle::New(0), old: None, new: None,
        },
        GraphIrEdit::ModifyDativeBondConstraint {
            id: GraphIrDativeBondHandle::Id(GraphIrDativeBondId(0)), old: None, new: None,
        },
        GraphIrEdit::ModifyAromaticSystemConstraint {
            id: GraphIrAromaticSystemHandle::New(0), old: None, new: None,
        },
        GraphIrEdit::ModifyMulticenterBondConstraint {
            id: GraphIrMulticenterBondHandle::Id(GraphIrMulticenterBondId(0)), old: None, new: None,
        },
        GraphIrEdit::ModifyNoncovalentBondConstraint {
            id: GraphIrNoncovalentBondHandle::New(0), old: None, new: None,
        },
        GraphIrEdit::ModifyStereoAtomConstraint {
            id: GraphIrStereoAtomHandle::Id(GraphIrStereoAtomId(0)),
            kind: Some(GraphIrStereoKind::Tetrahedral), old: None, new: None,
        },
        GraphIrEdit::ModifyStereoBondConstraint {
            id: GraphIrStereoBondHandle::New(0),
            kind: Some(GraphIrStereoKind::CisTrans), old: None, new: None,
        },
        GraphIrEdit::AddMoleculeConstraint {
            constraint: Box::new(GraphIrConstraintEdit::from(GraphIrConstraint::Molecule(
                GraphIrMoleculeConstraint::Connected { atoms: None },
            ))),
        },
        GraphIrEdit::RemoveMoleculeConstraint {
            constraint: Box::new(GraphIrConstraintEdit::from(GraphIrConstraint::Molecule(
                GraphIrMoleculeConstraint::Connected {
                    atoms: Some(vec![GraphIrAtomId(0)]),
                },
            ))),
        },
    ])]
    fn test_edit_roundtrip(#[case] edits: Vec<GraphIrEdit>) {
        Python::attach(|py| {
            for expected in edits {
                let edit = Edit::from_rust(py, &expected).unwrap();
                let rust = edit.to_rust(py);
                let rebuilt = Edit::from_rust(py, &rust).unwrap();

                assert_eq!(rust, expected);
                assert!(edit.__eq__(&rebuilt, py));
            }
        });
    }

    #[rstest]
    #[case::empty(GraphIrEdits::new())]
    #[case::atom(GraphIrEdits::from_iter([GraphIrEdit::AddAtoms {
        atoms: vec![GraphIrAtomForm::default()],
    }]))]
    fn test_edits_take(#[case] expected: GraphIrEdits) {
        Python::attach(|py| {
            let mut edits = Edits::from_rust(expected.clone());
            assert_eq!(edits.to_rust().unwrap(), &expected);
            assert_eq!(edits.take().unwrap(), expected);
            for error in [
                edits.to_rust().unwrap_err(),
                edits.to_rust_mut().unwrap_err(),
                edits.take().unwrap_err(),
            ] {
                assert!(error.is_instance_of::<ConsumedError>(py));
                assert_eq!(
                    error.value(py).str().unwrap().extract::<String>().unwrap(),
                    "Edits has been consumed"
                );
            }
        });
    }

    #[rstest]
    #[case::live(1, 0)]
    #[case::exhausted(1, 1)]
    #[case::empty(0, 0)]
    fn test_edit_iter_next_consumed(#[case] count: usize, #[case] consumed: usize) {
        Python::attach(|py| {
            let expected = GraphIrEdit::AddAtoms {
                atoms: vec![GraphIrAtomForm::default()],
            };
            let owner = Py::new(
                py,
                Edits::from_rust(GraphIrEdits::from_iter(vec![expected.clone(); count])),
            )
            .unwrap();
            let mut iter = Edits::__iter__(owner.clone_ref(py), py).unwrap();
            for _ in 0..consumed {
                assert_eq!(
                    iter.__next__(py).unwrap().unwrap().borrow(py).to_rust(py),
                    expected
                );
            }
            if consumed == count {
                assert!(iter.__next__(py).unwrap().is_none());
            }
            assert_eq!(owner.borrow_mut(py).take().unwrap().len(), count);
            let error = iter.__next__(py).unwrap_err();
            assert!(error.is_instance_of::<InvalidatedViewError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "Edit iterator is invalid because Edits has been consumed"
            );
        });
    }

    #[rstest]
    fn test_edit_iter_next_error() {
        Python::attach(|py| {
            let expected = GraphIrEdit::AddAtoms {
                atoms: vec![GraphIrAtomForm::default()],
            };
            let owner = Py::new(
                py,
                Edits::from_rust([expected.clone()].into_iter().collect()),
            )
            .unwrap();
            let mut iter = Edits::__iter__(owner.clone_ref(py), py).unwrap();
            let borrow = owner.borrow_mut(py);

            let error = iter.__next__(py).unwrap_err();
            assert!(error.is_instance_of::<PyRuntimeError>(py));
            drop(borrow);

            let entry = iter.__next__(py).unwrap().unwrap();
            assert_eq!(entry.borrow(py).to_rust(py), expected);
            assert!(owner.try_borrow_mut(py).is_ok());
            assert!(iter.__next__(py).unwrap().is_none());
        });
    }
}
