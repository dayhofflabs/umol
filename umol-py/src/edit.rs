//! Python edit values and symbolic creation handles.

use std::collections::HashMap;
use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use umol_graph_ir::dsl::EditsDsl as GraphIrEditsDsl;
use umol_graph_ir::ir::{
    AddBond as GraphIrAddBond, AromaticSystemHandle as GraphIrAromaticSystemHandle,
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
use crate::atom::{AtomForm, AtomUpdate};
use crate::bond::{BondForm, BondUpdate};
use crate::constraint::aromatic::AromaticSystemConstraintForm;
use crate::constraint::atom::AtomConstraintForm;
use crate::constraint::bond::BondConstraintForm;
use crate::constraint::dative::DativeBondConstraintForm;
use crate::constraint::molecule::Constraint;
use crate::constraint::multicenter::MulticenterBondConstraintForm;
use crate::constraint::noncovalent::NoncovalentBondConstraintForm;
use crate::constraint::stereo::{StereoAtomConstraintForm, StereoBondConstraintForm};
use crate::convert::into_py_variant;
use crate::dative::{DativeBondForm, DativeBondUpdate};
use crate::defaults::MoleculeDefaults;
use crate::delta::{
    AromaticSystemFieldChange, AtomFieldChange, BondFieldChange, DativeBondFieldChange,
    MulticenterBondFieldChange, NoncovalentBondFieldChange, StereoAtomFieldChange,
    StereoBondFieldChange,
};
use crate::error::parse_error;
use crate::metadata::Entity;
use crate::multicenter::{MulticenterBondForm, MulticenterBondUpdate};
use crate::noncovalent::{NoncovalentBondForm, NoncovalentBondUpdate};
use crate::stereo::{
    StereoAtomForm, StereoAtomUpdate, StereoBondForm, StereoBondUpdate, StereoKind,
    StereoLigandKind,
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
        let constraint = constraint.bind(py).borrow().to_rust(py);
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

    pub(crate) fn to_rust(&self) -> GraphIrConstraintEdit {
        self.0.clone()
    }
}

type BondAddition = ((HandleLike, HandleLike), Py<BondForm>);
type DativeBondAddition = (Vec<HandleLike>, Py<DativeBondForm>);
type AromaticSystemAddition = (Vec<HandleLike>, Py<AromaticSystemForm>);
type MulticenterBondAddition = (Vec<HandleLike>, Py<MulticenterBondForm>);
type NoncovalentBondAddition = ((HandleLike, HandleLike), Py<NoncovalentBondForm>);
type StereoAtomAddition = (HandleLike, Vec<StereoLigandInput>, Py<StereoAtomForm>);
type StereoBondAddition = (HandleLike, Vec<StereoLigandInput>, Py<StereoBondForm>);
type DativeBondRemoval = (HandleLike, Vec<HandleLike>, Py<DativeBondForm>);
type AromaticSystemRemoval = (HandleLike, Vec<HandleLike>, Py<AromaticSystemForm>);
type MulticenterBondRemoval = (HandleLike, Vec<HandleLike>, Py<MulticenterBondForm>);
type NoncovalentBondRemoval = (
    HandleLike,
    (HandleLike, HandleLike),
    Py<NoncovalentBondForm>,
);
type StereoLigandInput = (HandleLike, StereoLigandKind);
type StereoAtomRemovalEntry = (
    HandleLike,
    HandleLike,
    Vec<StereoLigandInput>,
    Py<StereoAtomForm>,
);
type StereoBondRemovalEntry = (
    HandleLike,
    HandleLike,
    Vec<StereoLigandInput>,
    Py<StereoBondForm>,
);

struct StereoAtomRemovals(Vec<StereoAtomRemovalEntry>);

impl FromPyObject<'_, '_> for StereoAtomRemovals {
    type Error = PyErr;

    fn extract(object: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        object.extract().map(Self)
    }
}

impl<'py> IntoPyObject<'py> for &StereoAtomRemovals {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        let entries = self
            .0
            .iter()
            .map(|(id, site, ligands, ast)| {
                (id.clone(), site.clone(), ligands.clone(), ast.clone_ref(py))
            })
            .collect::<Vec<_>>();
        entries.into_pyobject(py).map(Bound::into_any)
    }
}

struct StereoBondRemovals(Vec<StereoBondRemovalEntry>);

impl FromPyObject<'_, '_> for StereoBondRemovals {
    type Error = PyErr;

    fn extract(object: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        object.extract().map(Self)
    }
}

impl<'py> IntoPyObject<'py> for &StereoBondRemovals {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        let entries = self
            .0
            .iter()
            .map(|(id, site, ligands, ast)| {
                (id.clone(), site.clone(), ligands.clone(), ast.clone_ref(py))
            })
            .collect::<Vec<_>>();
        entries.into_pyobject(py).map(Bound::into_any)
    }
}

/// One raw host-specific molecule edit.
#[pyclass]
#[allow(
    private_interfaces,
    reason = "PyO3 exposes Python field values through private coercion adapters"
)]
pub enum Edit {
    AddAtoms {
        atoms: Vec<Py<AtomForm>>,
    },
    AddBonds {
        bonds: Vec<BondAddition>,
    },
    RemoveTopology {
        atoms: Vec<HandleLike>,
        bonds: Vec<HandleLike>,
    },
    ModifyAtomField {
        id: HandleLike,
        change: Py<AtomFieldChange>,
    },
    ModifyBondField {
        id: HandleLike,
        change: Py<BondFieldChange>,
    },
    AddDativeBond {
        atoms: Vec<HandleLike>,
        ast: Py<DativeBondForm>,
    },
    RemoveDativeBonds {
        removes: Vec<DativeBondRemoval>,
    },
    ModifyDativeBondField {
        id: HandleLike,
        change: Py<DativeBondFieldChange>,
    },
    AddAromaticSystem {
        atoms: Vec<HandleLike>,
        ast: Py<AromaticSystemForm>,
    },
    RemoveAromaticSystems {
        removes: Vec<AromaticSystemRemoval>,
    },
    ModifyAromaticSystemField {
        id: HandleLike,
        change: Py<AromaticSystemFieldChange>,
    },
    AddMulticenterBond {
        atoms: Vec<HandleLike>,
        ast: Py<MulticenterBondForm>,
    },
    RemoveMulticenterBonds {
        removes: Vec<MulticenterBondRemoval>,
    },
    ModifyMulticenterBondField {
        id: HandleLike,
        change: Py<MulticenterBondFieldChange>,
    },
    AddNoncovalentBond {
        atoms: (HandleLike, HandleLike),
        ast: Py<NoncovalentBondForm>,
    },
    RemoveNoncovalentBonds {
        removes: Vec<NoncovalentBondRemoval>,
    },
    ModifyNoncovalentBondField {
        id: HandleLike,
        change: Py<NoncovalentBondFieldChange>,
    },
    AddStereoAtom {
        site: HandleLike,
        ligands: Vec<StereoLigandInput>,
        ast: Py<StereoAtomForm>,
    },
    RemoveStereoAtoms {
        removes: StereoAtomRemovals,
    },
    ModifyStereoAtomField {
        id: HandleLike,
        change: Py<StereoAtomFieldChange>,
    },
    AddStereoBond {
        site: HandleLike,
        ligands: Vec<StereoLigandInput>,
        ast: Py<StereoBondForm>,
    },
    RemoveStereoBonds {
        removes: StereoBondRemovals,
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
    ModifyBondConstraint {
        id: HandleLike,
        old: Option<Py<BondConstraintForm>>,
        new: Option<Py<BondConstraintForm>>,
    },
    ModifyDativeBondConstraint {
        id: HandleLike,
        old: Option<Py<DativeBondConstraintForm>>,
        new: Option<Py<DativeBondConstraintForm>>,
    },
    ModifyAromaticSystemConstraint {
        id: HandleLike,
        old: Option<Py<AromaticSystemConstraintForm>>,
        new: Option<Py<AromaticSystemConstraintForm>>,
    },
    ModifyMulticenterBondConstraint {
        id: HandleLike,
        old: Option<Py<MulticenterBondConstraintForm>>,
        new: Option<Py<MulticenterBondConstraintForm>>,
    },
    ModifyNoncovalentBondConstraint {
        id: HandleLike,
        old: Option<Py<NoncovalentBondConstraintForm>>,
        new: Option<Py<NoncovalentBondConstraintForm>>,
    },
    ModifyStereoAtomConstraint {
        id: HandleLike,
        kind: Option<StereoKind>,
        old: Option<Py<StereoAtomConstraintForm>>,
        new: Option<Py<StereoAtomConstraintForm>>,
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
            Self::ModifyAtomField { .. } => ("ModifyAtomField", &["id", "change"]),
            Self::ModifyBondField { .. } => ("ModifyBondField", &["id", "change"]),
            Self::AddDativeBond { .. } => ("AddDativeBond", &["atoms", "ast"]),
            Self::RemoveDativeBonds { .. } => ("RemoveDativeBonds", &["removes"]),
            Self::ModifyDativeBondField { .. } => ("ModifyDativeBondField", &["id", "change"]),
            Self::AddAromaticSystem { .. } => ("AddAromaticSystem", &["atoms", "ast"]),
            Self::RemoveAromaticSystems { .. } => ("RemoveAromaticSystems", &["removes"]),
            Self::ModifyAromaticSystemField { .. } => {
                ("ModifyAromaticSystemField", &["id", "change"])
            }
            Self::AddMulticenterBond { .. } => ("AddMulticenterBond", &["atoms", "ast"]),
            Self::RemoveMulticenterBonds { .. } => ("RemoveMulticenterBonds", &["removes"]),
            Self::ModifyMulticenterBondField { .. } => {
                ("ModifyMulticenterBondField", &["id", "change"])
            }
            Self::AddNoncovalentBond { .. } => ("AddNoncovalentBond", &["atoms", "ast"]),
            Self::RemoveNoncovalentBonds { .. } => ("RemoveNoncovalentBonds", &["removes"]),
            Self::ModifyNoncovalentBondField { .. } => {
                ("ModifyNoncovalentBondField", &["id", "change"])
            }
            Self::AddStereoAtom { .. } => ("AddStereoAtom", &["site", "ligands", "ast"]),
            Self::RemoveStereoAtoms { .. } => ("RemoveStereoAtoms", &["removes"]),
            Self::ModifyStereoAtomField { .. } => ("ModifyStereoAtomField", &["id", "change"]),
            Self::AddStereoBond { .. } => ("AddStereoBond", &["site", "ligands", "ast"]),
            Self::RemoveStereoBonds { .. } => ("RemoveStereoBonds", &["removes"]),
            Self::ModifyStereoBondField { .. } => ("ModifyStereoBondField", &["id", "change"]),
            Self::ModifyAtomConstraint { .. } => ("ModifyAtomConstraint", &["id", "old", "new"]),
            Self::ModifyBondConstraint { .. } => ("ModifyBondConstraint", &["id", "old", "new"]),
            Self::ModifyDativeBondConstraint { .. } => {
                ("ModifyDativeBondConstraint", &["id", "old", "new"])
            }
            Self::ModifyAromaticSystemConstraint { .. } => {
                ("ModifyAromaticSystemConstraint", &["id", "old", "new"])
            }
            Self::ModifyMulticenterBondConstraint { .. } => {
                ("ModifyMulticenterBondConstraint", &["id", "old", "new"])
            }
            Self::ModifyNoncovalentBondConstraint { .. } => {
                ("ModifyNoncovalentBondConstraint", &["id", "old", "new"])
            }
            Self::ModifyStereoAtomConstraint { .. } => {
                ("ModifyStereoAtomConstraint", &["id", "kind", "old", "new"])
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
                    .cloned()
                    .map(|atom| Py::new(py, AtomForm::from_inner(atom)))
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
                            Py::new(py, BondForm::from_inner(bond.attributes.clone()))?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            GraphIrEdit::RemoveTopology { atoms, bonds } => Self::RemoveTopology {
                atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                bonds: bonds.iter().map(HandleLike::from_bond_handle).collect(),
            },
            GraphIrEdit::ModifyAtomField { id, change } => Self::ModifyAtomField {
                id: HandleLike::from_atom_handle(id),
                change: into_py_variant(py, AtomFieldChange::from_rust(py, change)?)?,
            },
            GraphIrEdit::ModifyBondField { id, change } => Self::ModifyBondField {
                id: HandleLike::from_bond_handle(id),
                change: into_py_variant(py, BondFieldChange::from_rust(py, change)?)?,
            },
            GraphIrEdit::AddDativeBond { atoms, attributes } => Self::AddDativeBond {
                atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                ast: Py::new(py, DativeBondForm::from_inner(attributes.clone()))?,
            },
            GraphIrEdit::RemoveDativeBonds { removes } => Self::RemoveDativeBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        Ok((
                            HandleLike::from_dative_bond_handle(id),
                            atoms.iter().map(HandleLike::from_atom_handle).collect(),
                            Py::new(py, DativeBondForm::from_inner(ast.clone()))?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            GraphIrEdit::ModifyDativeBondField { id, change } => Self::ModifyDativeBondField {
                id: HandleLike::from_dative_bond_handle(id),
                change: into_py_variant(py, DativeBondFieldChange::from_rust(py, change)?)?,
            },
            GraphIrEdit::AddAromaticSystem { atoms, attributes } => Self::AddAromaticSystem {
                atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                ast: Py::new(py, AromaticSystemForm::from_inner(attributes.clone()))?,
            },
            GraphIrEdit::RemoveAromaticSystems { removes } => Self::RemoveAromaticSystems {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        Ok((
                            HandleLike::from_aromatic_system_handle(id),
                            atoms.iter().map(HandleLike::from_atom_handle).collect(),
                            Py::new(py, AromaticSystemForm::from_inner(ast.clone()))?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            GraphIrEdit::ModifyAromaticSystemField { id, change } => {
                Self::ModifyAromaticSystemField {
                    id: HandleLike::from_aromatic_system_handle(id),
                    change: into_py_variant(py, AromaticSystemFieldChange::from_rust(py, change)?)?,
                }
            }
            GraphIrEdit::AddMulticenterBond { atoms, attributes } => Self::AddMulticenterBond {
                atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                ast: Py::new(py, MulticenterBondForm::from_inner(attributes.clone()))?,
            },
            GraphIrEdit::RemoveMulticenterBonds { removes } => Self::RemoveMulticenterBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        Ok((
                            HandleLike::from_multicenter_bond_handle(id),
                            atoms.iter().map(HandleLike::from_atom_handle).collect(),
                            Py::new(py, MulticenterBondForm::from_inner(ast.clone()))?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
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
                ast: Py::new(py, NoncovalentBondForm::from_inner(attributes.clone()))?,
            },
            GraphIrEdit::RemoveNoncovalentBonds { removes } => Self::RemoveNoncovalentBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        Ok((
                            HandleLike::from_noncovalent_bond_handle(id),
                            (
                                HandleLike::from_atom_handle(&atoms[0]),
                                HandleLike::from_atom_handle(&atoms[1]),
                            ),
                            Py::new(py, NoncovalentBondForm::from_inner(ast.clone()))?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            GraphIrEdit::ModifyNoncovalentBondField { id, change } => {
                Self::ModifyNoncovalentBondField {
                    id: HandleLike::from_noncovalent_bond_handle(id),
                    change: into_py_variant(
                        py,
                        NoncovalentBondFieldChange::from_rust(py, change)?,
                    )?,
                }
            }
            GraphIrEdit::AddStereoAtom {
                site,
                ligands,
                attributes,
            } => Self::AddStereoAtom {
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
                ast: Py::new(py, StereoAtomForm::from_inner(attributes.clone()))?,
            },
            GraphIrEdit::RemoveStereoAtoms { removes } => Self::RemoveStereoAtoms {
                removes: StereoAtomRemovals(
                    removes
                        .iter()
                        .map(|(id, site, ligands, ast)| {
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
                                Py::new(py, StereoAtomForm::from_inner(ast.clone()))?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            GraphIrEdit::ModifyStereoAtomField { id, change } => Self::ModifyStereoAtomField {
                id: HandleLike::from_stereo_atom_handle(id),
                change: into_py_variant(py, StereoAtomFieldChange::from_rust(py, change)?)?,
            },
            GraphIrEdit::AddStereoBond {
                site,
                ligands,
                attributes,
            } => Self::AddStereoBond {
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
                ast: Py::new(py, StereoBondForm::from_inner(attributes.clone()))?,
            },
            GraphIrEdit::RemoveStereoBonds { removes } => Self::RemoveStereoBonds {
                removes: StereoBondRemovals(
                    removes
                        .iter()
                        .map(|(id, site, ligands, ast)| {
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
                                Py::new(py, StereoBondForm::from_inner(ast.clone()))?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
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
                atoms: atoms
                    .iter()
                    .map(|atom| atom.bind(py).borrow().inner().clone())
                    .collect(),
            },
            Self::AddBonds { bonds } => GraphIrEdit::AddBonds {
                bonds: bonds
                    .iter()
                    .map(|((first, second), ast)| GraphIrAddBond {
                        endpoints: [first.to_atom_handle(), second.to_atom_handle()],
                        attributes: ast.bind(py).borrow().inner().clone(),
                    })
                    .collect(),
            },
            Self::RemoveTopology { atoms, bonds } => GraphIrEdit::RemoveTopology {
                atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                bonds: bonds.iter().map(HandleLike::to_bond_handle).collect(),
            },
            Self::ModifyAtomField { id, change } => GraphIrEdit::ModifyAtomField {
                id: id.to_atom_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyBondField { id, change } => GraphIrEdit::ModifyBondField {
                id: id.to_bond_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::AddDativeBond { atoms, ast } => GraphIrEdit::AddDativeBond {
                atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                attributes: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveDativeBonds { removes } => GraphIrEdit::RemoveDativeBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        (
                            id.to_dative_bond_handle(),
                            atoms.iter().map(HandleLike::to_atom_handle).collect(),
                            ast.bind(py).borrow().inner().clone(),
                        )
                    })
                    .collect(),
            },
            Self::ModifyDativeBondField { id, change } => GraphIrEdit::ModifyDativeBondField {
                id: id.to_dative_bond_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::AddAromaticSystem { atoms, ast } => GraphIrEdit::AddAromaticSystem {
                atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                attributes: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveAromaticSystems { removes } => GraphIrEdit::RemoveAromaticSystems {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        (
                            id.to_aromatic_system_handle(),
                            atoms.iter().map(HandleLike::to_atom_handle).collect(),
                            ast.bind(py).borrow().inner().clone(),
                        )
                    })
                    .collect(),
            },
            Self::ModifyAromaticSystemField { id, change } => {
                GraphIrEdit::ModifyAromaticSystemField {
                    id: id.to_aromatic_system_handle(),
                    change: change.bind(py).borrow().to_rust(py),
                }
            }
            Self::AddMulticenterBond { atoms, ast } => GraphIrEdit::AddMulticenterBond {
                atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                attributes: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveMulticenterBonds { removes } => GraphIrEdit::RemoveMulticenterBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        (
                            id.to_multicenter_bond_handle(),
                            atoms.iter().map(HandleLike::to_atom_handle).collect(),
                            ast.bind(py).borrow().inner().clone(),
                        )
                    })
                    .collect(),
            },
            Self::ModifyMulticenterBondField { id, change } => {
                GraphIrEdit::ModifyMulticenterBondField {
                    id: id.to_multicenter_bond_handle(),
                    change: change.bind(py).borrow().to_rust(py),
                }
            }
            Self::AddNoncovalentBond { atoms, ast } => GraphIrEdit::AddNoncovalentBond {
                atoms: [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                attributes: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveNoncovalentBonds { removes } => GraphIrEdit::RemoveNoncovalentBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        (
                            id.to_noncovalent_bond_handle(),
                            [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                            ast.bind(py).borrow().inner().clone(),
                        )
                    })
                    .collect(),
            },
            Self::ModifyNoncovalentBondField { id, change } => {
                GraphIrEdit::ModifyNoncovalentBondField {
                    id: id.to_noncovalent_bond_handle(),
                    change: change.bind(py).borrow().to_rust(py),
                }
            }
            Self::AddStereoAtom { site, ligands, ast } => GraphIrEdit::AddStereoAtom {
                site: site.to_atom_handle(),
                ligands: ligands
                    .iter()
                    .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                    .collect(),
                attributes: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveStereoAtoms { removes } => GraphIrEdit::RemoveStereoAtoms {
                removes: removes
                    .0
                    .iter()
                    .map(|(id, site, ligands, ast)| {
                        (
                            id.to_stereo_atom_handle(),
                            site.to_atom_handle(),
                            ligands
                                .iter()
                                .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                                .collect(),
                            ast.bind(py).borrow().inner().clone(),
                        )
                    })
                    .collect(),
            },
            Self::ModifyStereoAtomField { id, change } => GraphIrEdit::ModifyStereoAtomField {
                id: id.to_stereo_atom_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::AddStereoBond { site, ligands, ast } => GraphIrEdit::AddStereoBond {
                site: site.to_bond_handle(),
                ligands: ligands
                    .iter()
                    .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                    .collect(),
                attributes: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveStereoBonds { removes } => GraphIrEdit::RemoveStereoBonds {
                removes: removes
                    .0
                    .iter()
                    .map(|(id, site, ligands, ast)| {
                        (
                            id.to_stereo_bond_handle(),
                            site.to_bond_handle(),
                            ligands
                                .iter()
                                .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                                .collect(),
                            ast.bind(py).borrow().inner().clone(),
                        )
                    })
                    .collect(),
            },
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
            Self::ModifyBondConstraint { id, old, new } => GraphIrEdit::ModifyBondConstraint {
                id: id.to_bond_handle(),
                old: old
                    .as_ref()
                    .map(|value| value.bind(py).borrow().to_rust(py)),
                new: new
                    .as_ref()
                    .map(|value| value.bind(py).borrow().to_rust(py)),
            },
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
                constraint: constraint.bind(py).borrow().to_rust(),
            },
            Self::RemoveMoleculeConstraint { constraint } => {
                GraphIrEdit::RemoveMoleculeConstraint {
                    constraint: constraint.bind(py).borrow().to_rust(),
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

fn edit_iter(py: Python<'_>, edits: &GraphIrEdits) -> PyResult<EditIter> {
    let entries = edits
        .iter()
        .map(|edit| into_py_variant(py, Edit::from_rust(py, edit)?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(EditIter {
        entries: entries.into_iter(),
    })
}

/// A snapshot iterator over host-specific edits.
#[pyclass]
pub(crate) struct EditIter {
    entries: IntoIter<Py<Edit>>,
}

#[pymethods]
impl EditIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<Edit>> {
        self.entries.next()
    }
}

/// An ordered, append-only batch of host-specific molecule edits.
#[pyclass(eq)]
#[derive(Debug, PartialEq)]
pub struct Edits(GraphIrEdits);

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
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new).to_rust();
        let edits = GraphIrEditsDsl::from_str(text)
            .map_err(parse_error)?
            .into_ir(&defaults);
        Ok(Self::from_rust(edits))
    }

    #[pyo3(signature = (*, defaults=None))]
    fn render(&self, defaults: Option<MoleculeDefaults>) -> String {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new).to_rust();
        GraphIrEditsDsl::from_ir(&self.0, &defaults).to_string()
    }

    /// Append one detached raw edit and account for every entity it creates.
    fn append(&mut self, py: Python<'_>, edit: Py<Edit>) {
        self.0.push(edit.bind(py).borrow().to_rust(py));
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<Edit> {
        let index = resolve_edit_index(self.0.len(), index)?;
        Edit::from_rust(py, &self.0.as_slice()[index])
    }

    fn __iter__(&self, py: Python<'_>) -> PyResult<EditIter> {
        edit_iter(py, &self.0)
    }

    fn add_atom(&mut self, py: Python<'_>, ast: Py<AtomForm>) -> New {
        New::from_rust(GraphIrEntityHandle::Atom(
            self.0.add_atom(ast.bind(py).borrow().inner().clone()),
        ))
    }

    fn add_atoms(&mut self, py: Python<'_>, atoms: Vec<Py<AtomForm>>) -> Vec<New> {
        self.0
            .add_atoms(
                atoms
                    .into_iter()
                    .map(|ast| ast.bind(py).borrow().inner().clone()),
            )
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::Atom(handle)))
            .collect()
    }

    fn add_bond(
        &mut self,
        py: Python<'_>,
        first: HandleLike,
        second: HandleLike,
        ast: Py<BondForm>,
    ) -> New {
        New::from_rust(GraphIrEntityHandle::Bond(self.0.add_bond(
            first.to_atom_handle(),
            second.to_atom_handle(),
            ast.bind(py).borrow().inner().clone(),
        )))
    }

    fn add_bonds(&mut self, py: Python<'_>, bonds: Vec<BondAddition>) -> Vec<New> {
        self.0
            .add_bonds(
                bonds
                    .into_iter()
                    .map(|((first, second), ast)| GraphIrAddBond {
                        endpoints: [first.to_atom_handle(), second.to_atom_handle()],
                        attributes: ast.bind(py).borrow().inner().clone(),
                    }),
            )
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::Bond(handle)))
            .collect()
    }

    fn add_dative_bond(
        &mut self,
        py: Python<'_>,
        atoms: Vec<HandleLike>,
        ast: Py<DativeBondForm>,
    ) -> New {
        New::from_rust(GraphIrEntityHandle::DativeBond(self.0.add_dative_bond(
            atoms.iter().map(HandleLike::to_atom_handle).collect(),
            ast.bind(py).borrow().inner().clone(),
        )))
    }

    fn add_dative_bonds(&mut self, py: Python<'_>, bonds: Vec<DativeBondAddition>) -> Vec<New> {
        self.0
            .add_dative_bonds(bonds.into_iter().map(|(atoms, ast)| {
                (
                    atoms.iter().map(HandleLike::to_atom_handle).collect(),
                    ast.bind(py).borrow().inner().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::DativeBond(handle)))
            .collect()
    }

    fn add_aromatic_system(
        &mut self,
        py: Python<'_>,
        atoms: Vec<HandleLike>,
        ast: Py<AromaticSystemForm>,
    ) -> New {
        New::from_rust(GraphIrEntityHandle::AromaticSystem(
            self.0.add_aromatic_system(
                atoms.iter().map(HandleLike::to_atom_handle).collect(),
                ast.bind(py).borrow().inner().clone(),
            ),
        ))
    }

    fn add_aromatic_systems(
        &mut self,
        py: Python<'_>,
        systems: Vec<AromaticSystemAddition>,
    ) -> Vec<New> {
        self.0
            .add_aromatic_systems(systems.into_iter().map(|(atoms, ast)| {
                (
                    atoms.iter().map(HandleLike::to_atom_handle).collect(),
                    ast.bind(py).borrow().inner().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::AromaticSystem(handle)))
            .collect()
    }

    fn add_multicenter_bond(
        &mut self,
        py: Python<'_>,
        atoms: Vec<HandleLike>,
        ast: Py<MulticenterBondForm>,
    ) -> New {
        New::from_rust(GraphIrEntityHandle::MulticenterBond(
            self.0.add_multicenter_bond(
                atoms.iter().map(HandleLike::to_atom_handle).collect(),
                ast.bind(py).borrow().inner().clone(),
            ),
        ))
    }

    fn add_multicenter_bonds(
        &mut self,
        py: Python<'_>,
        bonds: Vec<MulticenterBondAddition>,
    ) -> Vec<New> {
        self.0
            .add_multicenter_bonds(bonds.into_iter().map(|(atoms, ast)| {
                (
                    atoms.iter().map(HandleLike::to_atom_handle).collect(),
                    ast.bind(py).borrow().inner().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::MulticenterBond(handle)))
            .collect()
    }

    fn add_noncovalent_bond(
        &mut self,
        py: Python<'_>,
        atoms: (HandleLike, HandleLike),
        ast: Py<NoncovalentBondForm>,
    ) -> New {
        New::from_rust(GraphIrEntityHandle::NoncovalentBond(
            self.0.add_noncovalent_bond(
                [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                ast.bind(py).borrow().inner().clone(),
            ),
        ))
    }

    fn add_noncovalent_bonds(
        &mut self,
        py: Python<'_>,
        bonds: Vec<NoncovalentBondAddition>,
    ) -> Vec<New> {
        self.0
            .add_noncovalent_bonds(bonds.into_iter().map(|(atoms, ast)| {
                (
                    [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                    ast.bind(py).borrow().inner().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::NoncovalentBond(handle)))
            .collect()
    }

    fn add_stereo_atom(
        &mut self,
        py: Python<'_>,
        site: HandleLike,
        ligands: Vec<StereoLigandInput>,
        ast: Py<StereoAtomForm>,
    ) -> New {
        New::from_rust(GraphIrEntityHandle::StereoAtom(
            self.0.add_stereo_atom(
                site.to_atom_handle(),
                ligands
                    .iter()
                    .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                    .collect(),
                ast.bind(py).borrow().inner().clone(),
            ),
        ))
    }

    fn add_stereo_atoms(&mut self, py: Python<'_>, atoms: Vec<StereoAtomAddition>) -> Vec<New> {
        self.0
            .add_stereo_atoms(atoms.into_iter().map(|(site, ligands, ast)| {
                (
                    site.to_atom_handle(),
                    ligands
                        .iter()
                        .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                        .collect(),
                    ast.bind(py).borrow().inner().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::StereoAtom(handle)))
            .collect()
    }

    fn add_stereo_bond(
        &mut self,
        py: Python<'_>,
        site: HandleLike,
        ligands: Vec<StereoLigandInput>,
        ast: Py<StereoBondForm>,
    ) -> New {
        New::from_rust(GraphIrEntityHandle::StereoBond(
            self.0.add_stereo_bond(
                site.to_bond_handle(),
                ligands
                    .iter()
                    .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                    .collect(),
                ast.bind(py).borrow().inner().clone(),
            ),
        ))
    }

    fn add_stereo_bonds(&mut self, py: Python<'_>, bonds: Vec<StereoBondAddition>) -> Vec<New> {
        self.0
            .add_stereo_bonds(bonds.into_iter().map(|(site, ligands, ast)| {
                (
                    site.to_bond_handle(),
                    ligands
                        .iter()
                        .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                        .collect(),
                    ast.bind(py).borrow().inner().clone(),
                )
            }))
            .into_iter()
            .map(|handle| New::from_rust(GraphIrEntityHandle::StereoBond(handle)))
            .collect()
    }

    fn remove_topology(&mut self, atoms: Vec<HandleLike>, bonds: Vec<HandleLike>) {
        self.0.remove_topology(
            atoms.iter().map(HandleLike::to_atom_handle).collect(),
            bonds.iter().map(HandleLike::to_bond_handle).collect(),
        );
    }

    fn remove_dative_bonds(&mut self, py: Python<'_>, removes: Vec<DativeBondRemoval>) {
        self.0.remove_dative_bonds(
            removes
                .into_iter()
                .map(|(id, atoms, ast)| {
                    (
                        id.to_dative_bond_handle(),
                        atoms.iter().map(HandleLike::to_atom_handle).collect(),
                        ast.bind(py).borrow().inner().clone(),
                    )
                })
                .collect(),
        );
    }

    fn remove_aromatic_systems(&mut self, py: Python<'_>, removes: Vec<AromaticSystemRemoval>) {
        self.0.remove_aromatic_systems(
            removes
                .into_iter()
                .map(|(id, atoms, ast)| {
                    (
                        id.to_aromatic_system_handle(),
                        atoms.iter().map(HandleLike::to_atom_handle).collect(),
                        ast.bind(py).borrow().inner().clone(),
                    )
                })
                .collect(),
        );
    }

    fn remove_multicenter_bonds(&mut self, py: Python<'_>, removes: Vec<MulticenterBondRemoval>) {
        self.0.remove_multicenter_bonds(
            removes
                .into_iter()
                .map(|(id, atoms, ast)| {
                    (
                        id.to_multicenter_bond_handle(),
                        atoms.iter().map(HandleLike::to_atom_handle).collect(),
                        ast.bind(py).borrow().inner().clone(),
                    )
                })
                .collect(),
        );
    }

    fn remove_noncovalent_bonds(&mut self, py: Python<'_>, removes: Vec<NoncovalentBondRemoval>) {
        self.0.remove_noncovalent_bonds(
            removes
                .into_iter()
                .map(|(id, atoms, ast)| {
                    (
                        id.to_noncovalent_bond_handle(),
                        [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                        ast.bind(py).borrow().inner().clone(),
                    )
                })
                .collect(),
        );
    }

    fn remove_stereo_atoms(&mut self, py: Python<'_>, removes: StereoAtomRemovals) {
        self.0.remove_stereo_atoms(
            removes
                .0
                .into_iter()
                .map(|(id, site, ligands, ast)| {
                    (
                        id.to_stereo_atom_handle(),
                        site.to_atom_handle(),
                        ligands
                            .iter()
                            .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                            .collect(),
                        ast.bind(py).borrow().inner().clone(),
                    )
                })
                .collect(),
        );
    }

    fn remove_stereo_bonds(&mut self, py: Python<'_>, removes: StereoBondRemovals) {
        self.0.remove_stereo_bonds(
            removes
                .0
                .into_iter()
                .map(|(id, site, ligands, ast)| {
                    (
                        id.to_stereo_bond_handle(),
                        site.to_bond_handle(),
                        ligands
                            .iter()
                            .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                            .collect(),
                        ast.bind(py).borrow().inner().clone(),
                    )
                })
                .collect(),
        );
    }

    fn add_molecule_constraint(&mut self, py: Python<'_>, constraint: Py<ConstraintEdit>) {
        self.0
            .add_molecule_constraint(constraint.bind(py).borrow().to_rust());
    }

    fn remove_molecule_constraint(&mut self, py: Python<'_>, constraint: Py<ConstraintEdit>) {
        self.0
            .remove_molecule_constraint(constraint.bind(py).borrow().to_rust());
    }

    fn update_atom(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<AtomForm>,
        update: Py<AtomUpdate>,
    ) {
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow().to_rust();
        self.0
            .update_atom(id.to_atom_handle(), current.inner(), &update);
    }

    fn update_bond(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<BondForm>,
        update: Py<BondUpdate>,
    ) {
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow().to_rust();
        self.0
            .update_bond(id.to_bond_handle(), current.inner(), &update);
    }

    fn update_dative_bond(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<DativeBondForm>,
        update: Py<DativeBondUpdate>,
    ) {
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow().to_rust();
        self.0
            .update_dative_bond(id.to_dative_bond_handle(), current.inner(), &update);
    }

    fn update_aromatic_system(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<AromaticSystemForm>,
        update: Py<AromaticSystemUpdate>,
    ) {
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow().to_rust();
        self.0
            .update_aromatic_system(id.to_aromatic_system_handle(), current.inner(), &update);
    }

    fn update_multicenter_bond(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<MulticenterBondForm>,
        update: Py<MulticenterBondUpdate>,
    ) {
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow().to_rust();
        self.0
            .update_multicenter_bond(id.to_multicenter_bond_handle(), current.inner(), &update);
    }

    fn update_noncovalent_bond(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<NoncovalentBondForm>,
        update: Py<NoncovalentBondUpdate>,
    ) {
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow().to_rust();
        self.0
            .update_noncovalent_bond(id.to_noncovalent_bond_handle(), current.inner(), &update);
    }

    fn update_stereo_atom(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<StereoAtomForm>,
        update: Py<StereoAtomUpdate>,
    ) {
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow().to_rust();
        self.0
            .update_stereo_atom(id.to_stereo_atom_handle(), current.inner(), &update);
    }

    fn update_stereo_bond(
        &mut self,
        py: Python<'_>,
        id: HandleLike,
        current: Py<StereoBondForm>,
        update: Py<StereoBondUpdate>,
    ) {
        let current = current.bind(py).borrow();
        let update = update.bind(py).borrow().to_rust();
        self.0
            .update_stereo_bond(id.to_stereo_bond_handle(), current.inner(), &update);
    }
}

impl Edits {
    pub(crate) fn from_rust(edits: GraphIrEdits) -> Self {
        Self(edits)
    }

    pub(crate) fn to_rust(&self) -> GraphIrEdits {
        self.0.clone()
    }
}

#[cfg(test)]
mod tests {
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

        assert_eq!(constraint.to_rust(), expected);
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
            atoms: vec![GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrAtomHandle::New(0)],
            attributes: GraphIrDativeBondForm::default(),
        },
        GraphIrEdit::RemoveDativeBonds { removes: vec![(
            GraphIrDativeBondHandle::New(0),
            vec![GraphIrAtomHandle::Id(GraphIrAtomId(0))],
            GraphIrDativeBondForm::default(),
        )] },
        GraphIrEdit::ModifyDativeBondField {
            id: GraphIrDativeBondHandle::Id(GraphIrDativeBondId(0)),
            change: GraphIrDativeBondFieldChange::Order {
                old: GraphIrNumForm::Lit(1),
                new: GraphIrNumForm::Lit(2),
            },
        },
        GraphIrEdit::AddAromaticSystem {
            atoms: vec![GraphIrAtomHandle::New(0)],
            attributes: GraphIrAromaticSystemForm::default(),
        },
        GraphIrEdit::RemoveAromaticSystems { removes: vec![(
            GraphIrAromaticSystemHandle::Id(GraphIrAromaticSystemId(0)),
            vec![GraphIrAtomHandle::New(0)],
            GraphIrAromaticSystemForm::default(),
        )] },
        GraphIrEdit::ModifyAromaticSystemField {
            id: GraphIrAromaticSystemHandle::New(0),
            change: GraphIrAromaticSystemFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(-1),
            },
        },
        GraphIrEdit::AddMulticenterBond {
            atoms: vec![GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrAtomHandle::New(0)],
            attributes: GraphIrMulticenterBondForm::default(),
        },
        GraphIrEdit::RemoveMulticenterBonds { removes: vec![(
            GraphIrMulticenterBondHandle::New(0),
            vec![GraphIrAtomHandle::Id(GraphIrAtomId(0))],
            GraphIrMulticenterBondForm::default(),
        )] },
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
        GraphIrEdit::ModifyNoncovalentBondField {
            id: GraphIrNoncovalentBondHandle::Id(GraphIrNoncovalentBondId(0)),
            change: GraphIrNoncovalentBondFieldChange::Kind {
                old: GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond),
                new: GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::Ionic),
            },
        },
        GraphIrEdit::AddStereoAtom {
            site: GraphIrAtomHandle::New(0),
            ligands: vec![(GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrStereoLigandKind::Atom)],
            attributes: GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, 0_u32),
        },
        GraphIrEdit::RemoveStereoAtoms { removes: vec![(
            GraphIrStereoAtomHandle::Id(GraphIrStereoAtomId(0)),
            GraphIrAtomHandle::New(0),
            vec![(GraphIrAtomHandle::Id(GraphIrAtomId(0)), GraphIrStereoLigandKind::ImplicitHydrogen)],
            GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, 1_u32),
        )] },
        GraphIrEdit::ModifyStereoAtomField {
            id: GraphIrStereoAtomHandle::New(0),
            change: GraphIrStereoAtomFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::kinded(GraphIrStereoKind::Tetrahedral, 0_u32),
                new: GraphIrStereoConfigurationForm::kinded(GraphIrStereoKind::Tetrahedral, 1_u32),
            },
        },
        GraphIrEdit::AddStereoBond {
            site: GraphIrBondHandle::Id(GraphIrBondId(0)),
            ligands: vec![(GraphIrAtomHandle::New(0), GraphIrStereoLigandKind::LonePair)],
            attributes: GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, 0_u32),
        },
        GraphIrEdit::RemoveStereoBonds { removes: vec![(
            GraphIrStereoBondHandle::New(0),
            GraphIrBondHandle::Id(GraphIrBondId(0)),
            vec![(GraphIrAtomHandle::New(0), GraphIrStereoLigandKind::Atom)],
            GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, 1_u32),
        )] },
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
            constraint: GraphIrConstraintEdit::from(GraphIrConstraint::Molecule(
                GraphIrMoleculeConstraint::Connected { atoms: None },
            )),
        },
        GraphIrEdit::RemoveMoleculeConstraint {
            constraint: GraphIrConstraintEdit::from(GraphIrConstraint::Molecule(
                GraphIrMoleculeConstraint::Connected { atoms: Some(vec![GraphIrAtomId(0)]) },
            )),
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
}
