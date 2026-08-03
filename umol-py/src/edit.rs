//! Python edit values and symbolic creation handles.

use std::collections::HashMap;
use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use umol_ast::ast::{
    AddBond as AstAddBond, AromaticSystemHandle as AstAromaticSystemHandle,
    AromaticSystemId as AstAromaticSystemId, AtomHandle as AstAtomHandle, AtomId as AstAtomId,
    BondHandle as AstBondHandle, BondId as AstBondId, ConstraintEdit as AstConstraintEdit,
    DativeBondHandle as AstDativeBondHandle, DativeBondId as AstDativeBondId, Edit as AstEdit,
    Edits as AstEdits, Entity as AstEntity, EntityHandle as AstEntityHandle,
    MulticenterBondHandle as AstMulticenterBondHandle, MulticenterBondId as AstMulticenterBondId,
    NoncovalentBondHandle as AstNoncovalentBondHandle, NoncovalentBondId as AstNoncovalentBondId,
    StereoAtomHandle as AstStereoAtomHandle, StereoAtomId as AstStereoAtomId,
    StereoBondHandle as AstStereoBondHandle, StereoBondId as AstStereoBondId,
};

use crate::aromatic::AromaticSystemAst;
use crate::atom::AtomAst;
use crate::bond::BondAst;
use crate::constraint::aromatic::AromaticSystemConstraintAst;
use crate::constraint::atom::AtomConstraintAst;
use crate::constraint::bond::BondConstraintAst;
use crate::constraint::dative::DativeBondConstraintAst;
use crate::constraint::molecule::Constraint;
use crate::constraint::multicenter::MulticenterBondConstraintAst;
use crate::constraint::noncovalent::NoncovalentBondConstraintAst;
use crate::constraint::stereo::{StereoAtomConstraintAst, StereoBondConstraintAst};
use crate::convert::into_py_variant;
use crate::dative::DativeBondAst;
use crate::delta::{
    AromaticSystemFieldChange, AtomFieldChange, BondFieldChange, DativeBondFieldChange,
    MulticenterBondFieldChange, NoncovalentBondFieldChange, StereoAtomFieldChange,
    StereoBondFieldChange,
};
use crate::metadata::Entity;
use crate::multicenter::MulticenterBondAst;
use crate::noncovalent::NoncovalentBondAst;
use crate::stereo::{StereoAtomAst, StereoBondAst, StereoKind, StereoLigandKind};

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
    fn from_rust(handle: AstEntityHandle) -> Self {
        let index = match handle {
            AstEntityHandle::Atom(AstAtomHandle::New(index))
            | AstEntityHandle::Bond(AstBondHandle::New(index))
            | AstEntityHandle::DativeBond(AstDativeBondHandle::New(index))
            | AstEntityHandle::AromaticSystem(AstAromaticSystemHandle::New(index))
            | AstEntityHandle::MulticenterBond(AstMulticenterBondHandle::New(index))
            | AstEntityHandle::NoncovalentBond(AstNoncovalentBondHandle::New(index))
            | AstEntityHandle::StereoAtom(AstStereoAtomHandle::New(index))
            | AstEntityHandle::StereoBond(AstStereoBondHandle::New(index)) => index,
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
    fn from_atom_handle(handle: &AstAtomHandle) -> Self {
        match handle {
            AstAtomHandle::Id(id) => Self::Id(id.0),
            AstAtomHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_bond_handle(handle: &AstBondHandle) -> Self {
        match handle {
            AstBondHandle::Id(id) => Self::Id(id.0),
            AstBondHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_dative_bond_handle(handle: &AstDativeBondHandle) -> Self {
        match handle {
            AstDativeBondHandle::Id(id) => Self::Id(id.0),
            AstDativeBondHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_aromatic_system_handle(handle: &AstAromaticSystemHandle) -> Self {
        match handle {
            AstAromaticSystemHandle::Id(id) => Self::Id(id.0),
            AstAromaticSystemHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_multicenter_bond_handle(handle: &AstMulticenterBondHandle) -> Self {
        match handle {
            AstMulticenterBondHandle::Id(id) => Self::Id(id.0),
            AstMulticenterBondHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_noncovalent_bond_handle(handle: &AstNoncovalentBondHandle) -> Self {
        match handle {
            AstNoncovalentBondHandle::Id(id) => Self::Id(id.0),
            AstNoncovalentBondHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_stereo_atom_handle(handle: &AstStereoAtomHandle) -> Self {
        match handle {
            AstStereoAtomHandle::Id(id) => Self::Id(id.0),
            AstStereoAtomHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    fn from_stereo_bond_handle(handle: &AstStereoBondHandle) -> Self {
        match handle {
            AstStereoBondHandle::Id(id) => Self::Id(id.0),
            AstStereoBondHandle::New(index) => Self::New(New { index: *index }),
        }
    }

    pub(crate) fn to_atom_handle(&self) -> AstAtomHandle {
        match self {
            Self::Id(index) => AstAtomHandle::Id(AstAtomId(*index)),
            Self::New(new) => AstAtomHandle::New(new.index),
        }
    }

    pub(crate) fn to_bond_handle(&self) -> AstBondHandle {
        match self {
            Self::Id(index) => AstBondHandle::Id(AstBondId(*index)),
            Self::New(new) => AstBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_dative_bond_handle(&self) -> AstDativeBondHandle {
        match self {
            Self::Id(index) => AstDativeBondHandle::Id(AstDativeBondId(*index)),
            Self::New(new) => AstDativeBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_aromatic_system_handle(&self) -> AstAromaticSystemHandle {
        match self {
            Self::Id(index) => AstAromaticSystemHandle::Id(AstAromaticSystemId(*index)),
            Self::New(new) => AstAromaticSystemHandle::New(new.index),
        }
    }

    pub(crate) fn to_multicenter_bond_handle(&self) -> AstMulticenterBondHandle {
        match self {
            Self::Id(index) => AstMulticenterBondHandle::Id(AstMulticenterBondId(*index)),
            Self::New(new) => AstMulticenterBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_noncovalent_bond_handle(&self) -> AstNoncovalentBondHandle {
        match self {
            Self::Id(index) => AstNoncovalentBondHandle::Id(AstNoncovalentBondId(*index)),
            Self::New(new) => AstNoncovalentBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_stereo_atom_handle(&self) -> AstStereoAtomHandle {
        match self {
            Self::Id(index) => AstStereoAtomHandle::Id(AstStereoAtomId(*index)),
            Self::New(new) => AstStereoAtomHandle::New(new.index),
        }
    }

    pub(crate) fn to_stereo_bond_handle(&self) -> AstStereoBondHandle {
        match self {
            Self::Id(index) => AstStereoBondHandle::Id(AstStereoBondId(*index)),
            Self::New(new) => AstStereoBondHandle::New(new.index),
        }
    }

    fn to_entity_handle(&self, entity: AstEntity) -> AstEntityHandle {
        match entity {
            AstEntity::Atom(_) => AstEntityHandle::Atom(self.to_atom_handle()),
            AstEntity::Bond(_) => AstEntityHandle::Bond(self.to_bond_handle()),
            AstEntity::DativeBond(_) => AstEntityHandle::DativeBond(self.to_dative_bond_handle()),
            AstEntity::AromaticSystem(_) => {
                AstEntityHandle::AromaticSystem(self.to_aromatic_system_handle())
            }
            AstEntity::MulticenterBond(_) => {
                AstEntityHandle::MulticenterBond(self.to_multicenter_bond_handle())
            }
            AstEntity::NoncovalentBond(_) => {
                AstEntityHandle::NoncovalentBond(self.to_noncovalent_bond_handle())
            }
            AstEntity::StereoAtom(_) => AstEntityHandle::StereoAtom(self.to_stereo_atom_handle()),
            AstEntity::StereoBond(_) => AstEntityHandle::StereoBond(self.to_stereo_bond_handle()),
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
pub struct ConstraintEdit(AstConstraintEdit);

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
            return Ok(Self(AstConstraintEdit::from(constraint)));
        };
        AstConstraintEdit::new(constraint, |entity| {
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
    pub(crate) fn from_rust(constraint: &AstConstraintEdit) -> Self {
        Self(constraint.clone())
    }

    pub(crate) fn to_rust(&self) -> AstConstraintEdit {
        self.0.clone()
    }
}

type BondAddition = ((HandleLike, HandleLike), Py<BondAst>);
type DativeBondAddition = (Vec<HandleLike>, Py<DativeBondAst>);
type AromaticSystemAddition = (Vec<HandleLike>, Py<AromaticSystemAst>);
type MulticenterBondAddition = (Vec<HandleLike>, Py<MulticenterBondAst>);
type NoncovalentBondAddition = ((HandleLike, HandleLike), Py<NoncovalentBondAst>);
type StereoAtomAddition = (HandleLike, Vec<StereoLigandInput>, Py<StereoAtomAst>);
type StereoBondAddition = (HandleLike, Vec<StereoLigandInput>, Py<StereoBondAst>);
type DativeBondRemoval = (HandleLike, Vec<HandleLike>, Py<DativeBondAst>);
type AromaticSystemRemoval = (HandleLike, Vec<HandleLike>, Py<AromaticSystemAst>);
type MulticenterBondRemoval = (HandleLike, Vec<HandleLike>, Py<MulticenterBondAst>);
type NoncovalentBondRemoval = (HandleLike, (HandleLike, HandleLike), Py<NoncovalentBondAst>);
type StereoLigandInput = (HandleLike, StereoLigandKind);
type StereoAtomRemovalEntry = (
    HandleLike,
    HandleLike,
    Vec<StereoLigandInput>,
    Py<StereoAtomAst>,
);
type StereoBondRemovalEntry = (
    HandleLike,
    HandleLike,
    Vec<StereoLigandInput>,
    Py<StereoBondAst>,
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
        atoms: Vec<Py<AtomAst>>,
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
        ast: Py<DativeBondAst>,
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
        ast: Py<AromaticSystemAst>,
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
        ast: Py<MulticenterBondAst>,
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
        ast: Py<NoncovalentBondAst>,
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
        ast: Py<StereoAtomAst>,
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
        ast: Py<StereoBondAst>,
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
        old: Option<Py<AtomConstraintAst>>,
        new: Option<Py<AtomConstraintAst>>,
    },
    ModifyBondConstraint {
        id: HandleLike,
        old: Option<Py<BondConstraintAst>>,
        new: Option<Py<BondConstraintAst>>,
    },
    ModifyDativeBondConstraint {
        id: HandleLike,
        old: Option<Py<DativeBondConstraintAst>>,
        new: Option<Py<DativeBondConstraintAst>>,
    },
    ModifyAromaticSystemConstraint {
        id: HandleLike,
        old: Option<Py<AromaticSystemConstraintAst>>,
        new: Option<Py<AromaticSystemConstraintAst>>,
    },
    ModifyMulticenterBondConstraint {
        id: HandleLike,
        old: Option<Py<MulticenterBondConstraintAst>>,
        new: Option<Py<MulticenterBondConstraintAst>>,
    },
    ModifyNoncovalentBondConstraint {
        id: HandleLike,
        old: Option<Py<NoncovalentBondConstraintAst>>,
        new: Option<Py<NoncovalentBondConstraintAst>>,
    },
    ModifyStereoAtomConstraint {
        id: HandleLike,
        kind: Option<StereoKind>,
        old: Option<Py<StereoAtomConstraintAst>>,
        new: Option<Py<StereoAtomConstraintAst>>,
    },
    ModifyStereoBondConstraint {
        id: HandleLike,
        kind: Option<StereoKind>,
        old: Option<Py<StereoBondConstraintAst>>,
        new: Option<Py<StereoBondConstraintAst>>,
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
    pub(crate) fn from_rust(py: Python<'_>, edit: &AstEdit) -> PyResult<Self> {
        Ok(match edit {
            AstEdit::AddAtoms { atoms } => Self::AddAtoms {
                atoms: atoms
                    .iter()
                    .cloned()
                    .map(|atom| Py::new(py, AtomAst::from_inner(atom)))
                    .collect::<PyResult<_>>()?,
            },
            AstEdit::AddBonds { bonds } => Self::AddBonds {
                bonds: bonds
                    .iter()
                    .map(|bond| {
                        Ok((
                            (
                                HandleLike::from_atom_handle(&bond.endpoints[0]),
                                HandleLike::from_atom_handle(&bond.endpoints[1]),
                            ),
                            Py::new(py, BondAst::from_inner(bond.ast.clone()))?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            AstEdit::RemoveTopology { atoms, bonds } => Self::RemoveTopology {
                atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                bonds: bonds.iter().map(HandleLike::from_bond_handle).collect(),
            },
            AstEdit::ModifyAtomField { id, change } => Self::ModifyAtomField {
                id: HandleLike::from_atom_handle(id),
                change: into_py_variant(py, AtomFieldChange::from_rust(py, change)?)?,
            },
            AstEdit::ModifyBondField { id, change } => Self::ModifyBondField {
                id: HandleLike::from_bond_handle(id),
                change: into_py_variant(py, BondFieldChange::from_rust(py, change)?)?,
            },
            AstEdit::AddDativeBond { atoms, ast } => Self::AddDativeBond {
                atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                ast: Py::new(py, DativeBondAst::from_inner(ast.clone()))?,
            },
            AstEdit::RemoveDativeBonds { removes } => Self::RemoveDativeBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        Ok((
                            HandleLike::from_dative_bond_handle(id),
                            atoms.iter().map(HandleLike::from_atom_handle).collect(),
                            Py::new(py, DativeBondAst::from_inner(ast.clone()))?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            AstEdit::ModifyDativeBondField { id, change } => Self::ModifyDativeBondField {
                id: HandleLike::from_dative_bond_handle(id),
                change: into_py_variant(py, DativeBondFieldChange::from_rust(py, change)?)?,
            },
            AstEdit::AddAromaticSystem { atoms, ast } => Self::AddAromaticSystem {
                atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                ast: Py::new(py, AromaticSystemAst::from_inner(ast.clone()))?,
            },
            AstEdit::RemoveAromaticSystems { removes } => Self::RemoveAromaticSystems {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        Ok((
                            HandleLike::from_aromatic_system_handle(id),
                            atoms.iter().map(HandleLike::from_atom_handle).collect(),
                            Py::new(py, AromaticSystemAst::from_inner(ast.clone()))?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            AstEdit::ModifyAromaticSystemField { id, change } => Self::ModifyAromaticSystemField {
                id: HandleLike::from_aromatic_system_handle(id),
                change: into_py_variant(py, AromaticSystemFieldChange::from_rust(py, change)?)?,
            },
            AstEdit::AddMulticenterBond { atoms, ast } => Self::AddMulticenterBond {
                atoms: atoms.iter().map(HandleLike::from_atom_handle).collect(),
                ast: Py::new(py, MulticenterBondAst::from_inner(ast.clone()))?,
            },
            AstEdit::RemoveMulticenterBonds { removes } => Self::RemoveMulticenterBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        Ok((
                            HandleLike::from_multicenter_bond_handle(id),
                            atoms.iter().map(HandleLike::from_atom_handle).collect(),
                            Py::new(py, MulticenterBondAst::from_inner(ast.clone()))?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            AstEdit::ModifyMulticenterBondField { id, change } => {
                Self::ModifyMulticenterBondField {
                    id: HandleLike::from_multicenter_bond_handle(id),
                    change: into_py_variant(
                        py,
                        MulticenterBondFieldChange::from_rust(py, change)?,
                    )?,
                }
            }
            AstEdit::AddNoncovalentBond { atoms, ast } => Self::AddNoncovalentBond {
                atoms: (
                    HandleLike::from_atom_handle(&atoms[0]),
                    HandleLike::from_atom_handle(&atoms[1]),
                ),
                ast: Py::new(py, NoncovalentBondAst::from_inner(ast.clone()))?,
            },
            AstEdit::RemoveNoncovalentBonds { removes } => Self::RemoveNoncovalentBonds {
                removes: removes
                    .iter()
                    .map(|(id, atoms, ast)| {
                        Ok((
                            HandleLike::from_noncovalent_bond_handle(id),
                            (
                                HandleLike::from_atom_handle(&atoms[0]),
                                HandleLike::from_atom_handle(&atoms[1]),
                            ),
                            Py::new(py, NoncovalentBondAst::from_inner(ast.clone()))?,
                        ))
                    })
                    .collect::<PyResult<_>>()?,
            },
            AstEdit::ModifyNoncovalentBondField { id, change } => {
                Self::ModifyNoncovalentBondField {
                    id: HandleLike::from_noncovalent_bond_handle(id),
                    change: into_py_variant(
                        py,
                        NoncovalentBondFieldChange::from_rust(py, change)?,
                    )?,
                }
            }
            AstEdit::AddStereoAtom { site, ligands, ast } => Self::AddStereoAtom {
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
                ast: Py::new(py, StereoAtomAst::from_inner(ast.clone()))?,
            },
            AstEdit::RemoveStereoAtoms { removes } => Self::RemoveStereoAtoms {
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
                                Py::new(py, StereoAtomAst::from_inner(ast.clone()))?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            AstEdit::ModifyStereoAtomField { id, change } => Self::ModifyStereoAtomField {
                id: HandleLike::from_stereo_atom_handle(id),
                change: into_py_variant(py, StereoAtomFieldChange::from_rust(py, change)?)?,
            },
            AstEdit::AddStereoBond { site, ligands, ast } => Self::AddStereoBond {
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
                ast: Py::new(py, StereoBondAst::from_inner(ast.clone()))?,
            },
            AstEdit::RemoveStereoBonds { removes } => Self::RemoveStereoBonds {
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
                                Py::new(py, StereoBondAst::from_inner(ast.clone()))?,
                            ))
                        })
                        .collect::<PyResult<_>>()?,
                ),
            },
            AstEdit::ModifyStereoBondField { id, change } => Self::ModifyStereoBondField {
                id: HandleLike::from_stereo_bond_handle(id),
                change: into_py_variant(py, StereoBondFieldChange::from_rust(py, change)?)?,
            },
            AstEdit::ModifyAtomConstraint { id, old, new } => Self::ModifyAtomConstraint {
                id: HandleLike::from_atom_handle(id),
                old: old
                    .as_ref()
                    .map(|value| into_py_variant(py, AtomConstraintAst::from_rust(py, value)?))
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|value| into_py_variant(py, AtomConstraintAst::from_rust(py, value)?))
                    .transpose()?,
            },
            AstEdit::ModifyBondConstraint { id, old, new } => Self::ModifyBondConstraint {
                id: HandleLike::from_bond_handle(id),
                old: old
                    .as_ref()
                    .map(|value| into_py_variant(py, BondConstraintAst::from_rust(py, value)?))
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|value| into_py_variant(py, BondConstraintAst::from_rust(py, value)?))
                    .transpose()?,
            },
            AstEdit::ModifyDativeBondConstraint { id, old, new } => {
                Self::ModifyDativeBondConstraint {
                    id: HandleLike::from_dative_bond_handle(id),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, DativeBondConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, DativeBondConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                }
            }
            AstEdit::ModifyAromaticSystemConstraint { id, old, new } => {
                Self::ModifyAromaticSystemConstraint {
                    id: HandleLike::from_aromatic_system_handle(id),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, AromaticSystemConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, AromaticSystemConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                }
            }
            AstEdit::ModifyMulticenterBondConstraint { id, old, new } => {
                Self::ModifyMulticenterBondConstraint {
                    id: HandleLike::from_multicenter_bond_handle(id),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, MulticenterBondConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, MulticenterBondConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                }
            }
            AstEdit::ModifyNoncovalentBondConstraint { id, old, new } => {
                Self::ModifyNoncovalentBondConstraint {
                    id: HandleLike::from_noncovalent_bond_handle(id),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, NoncovalentBondConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, NoncovalentBondConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                }
            }
            AstEdit::ModifyStereoAtomConstraint { id, kind, old, new } => {
                Self::ModifyStereoAtomConstraint {
                    id: HandleLike::from_stereo_atom_handle(id),
                    kind: kind.map(StereoKind::from_rust),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, StereoAtomConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, StereoAtomConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                }
            }
            AstEdit::ModifyStereoBondConstraint { id, kind, old, new } => {
                Self::ModifyStereoBondConstraint {
                    id: HandleLike::from_stereo_bond_handle(id),
                    kind: kind.map(StereoKind::from_rust),
                    old: old
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, StereoBondConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|value| {
                            into_py_variant(py, StereoBondConstraintAst::from_rust(py, value)?)
                        })
                        .transpose()?,
                }
            }
            AstEdit::AddMoleculeConstraint { constraint } => Self::AddMoleculeConstraint {
                constraint: Py::new(py, ConstraintEdit::from_rust(constraint))?,
            },
            AstEdit::RemoveMoleculeConstraint { constraint } => Self::RemoveMoleculeConstraint {
                constraint: Py::new(py, ConstraintEdit::from_rust(constraint))?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstEdit {
        match self {
            Self::AddAtoms { atoms } => AstEdit::AddAtoms {
                atoms: atoms
                    .iter()
                    .map(|atom| atom.bind(py).borrow().inner().clone())
                    .collect(),
            },
            Self::AddBonds { bonds } => AstEdit::AddBonds {
                bonds: bonds
                    .iter()
                    .map(|((first, second), ast)| AstAddBond {
                        endpoints: [first.to_atom_handle(), second.to_atom_handle()],
                        ast: ast.bind(py).borrow().inner().clone(),
                    })
                    .collect(),
            },
            Self::RemoveTopology { atoms, bonds } => AstEdit::RemoveTopology {
                atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                bonds: bonds.iter().map(HandleLike::to_bond_handle).collect(),
            },
            Self::ModifyAtomField { id, change } => AstEdit::ModifyAtomField {
                id: id.to_atom_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyBondField { id, change } => AstEdit::ModifyBondField {
                id: id.to_bond_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::AddDativeBond { atoms, ast } => AstEdit::AddDativeBond {
                atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                ast: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveDativeBonds { removes } => AstEdit::RemoveDativeBonds {
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
            Self::ModifyDativeBondField { id, change } => AstEdit::ModifyDativeBondField {
                id: id.to_dative_bond_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::AddAromaticSystem { atoms, ast } => AstEdit::AddAromaticSystem {
                atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                ast: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveAromaticSystems { removes } => AstEdit::RemoveAromaticSystems {
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
            Self::ModifyAromaticSystemField { id, change } => AstEdit::ModifyAromaticSystemField {
                id: id.to_aromatic_system_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::AddMulticenterBond { atoms, ast } => AstEdit::AddMulticenterBond {
                atoms: atoms.iter().map(HandleLike::to_atom_handle).collect(),
                ast: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveMulticenterBonds { removes } => AstEdit::RemoveMulticenterBonds {
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
                AstEdit::ModifyMulticenterBondField {
                    id: id.to_multicenter_bond_handle(),
                    change: change.bind(py).borrow().to_rust(py),
                }
            }
            Self::AddNoncovalentBond { atoms, ast } => AstEdit::AddNoncovalentBond {
                atoms: [atoms.0.to_atom_handle(), atoms.1.to_atom_handle()],
                ast: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveNoncovalentBonds { removes } => AstEdit::RemoveNoncovalentBonds {
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
                AstEdit::ModifyNoncovalentBondField {
                    id: id.to_noncovalent_bond_handle(),
                    change: change.bind(py).borrow().to_rust(py),
                }
            }
            Self::AddStereoAtom { site, ligands, ast } => AstEdit::AddStereoAtom {
                site: site.to_atom_handle(),
                ligands: ligands
                    .iter()
                    .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                    .collect(),
                ast: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveStereoAtoms { removes } => AstEdit::RemoveStereoAtoms {
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
            Self::ModifyStereoAtomField { id, change } => AstEdit::ModifyStereoAtomField {
                id: id.to_stereo_atom_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::AddStereoBond { site, ligands, ast } => AstEdit::AddStereoBond {
                site: site.to_bond_handle(),
                ligands: ligands
                    .iter()
                    .map(|(atom, kind)| (atom.to_atom_handle(), kind.to_rust()))
                    .collect(),
                ast: ast.bind(py).borrow().inner().clone(),
            },
            Self::RemoveStereoBonds { removes } => AstEdit::RemoveStereoBonds {
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
            Self::ModifyStereoBondField { id, change } => AstEdit::ModifyStereoBondField {
                id: id.to_stereo_bond_handle(),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyAtomConstraint { id, old, new } => AstEdit::ModifyAtomConstraint {
                id: id.to_atom_handle(),
                old: old
                    .as_ref()
                    .map(|value| value.bind(py).borrow().to_rust(py)),
                new: new
                    .as_ref()
                    .map(|value| value.bind(py).borrow().to_rust(py)),
            },
            Self::ModifyBondConstraint { id, old, new } => AstEdit::ModifyBondConstraint {
                id: id.to_bond_handle(),
                old: old
                    .as_ref()
                    .map(|value| value.bind(py).borrow().to_rust(py)),
                new: new
                    .as_ref()
                    .map(|value| value.bind(py).borrow().to_rust(py)),
            },
            Self::ModifyDativeBondConstraint { id, old, new } => {
                AstEdit::ModifyDativeBondConstraint {
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
                AstEdit::ModifyAromaticSystemConstraint {
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
                AstEdit::ModifyMulticenterBondConstraint {
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
                AstEdit::ModifyNoncovalentBondConstraint {
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
                AstEdit::ModifyStereoAtomConstraint {
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
                AstEdit::ModifyStereoBondConstraint {
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
            Self::AddMoleculeConstraint { constraint } => AstEdit::AddMoleculeConstraint {
                constraint: constraint.bind(py).borrow().to_rust(),
            },
            Self::RemoveMoleculeConstraint { constraint } => AstEdit::RemoveMoleculeConstraint {
                constraint: constraint.bind(py).borrow().to_rust(),
            },
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

fn edit_iter(py: Python<'_>, edits: &AstEdits) -> PyResult<EditIter> {
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
pub struct Edits(AstEdits);

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

    fn add_atom(&mut self, py: Python<'_>, ast: Py<AtomAst>) -> New {
        New::from_rust(AstEntityHandle::Atom(
            self.0.add_atom(ast.bind(py).borrow().inner().clone()),
        ))
    }

    fn add_atoms(&mut self, py: Python<'_>, atoms: Vec<Py<AtomAst>>) -> Vec<New> {
        self.0
            .add_atoms(
                atoms
                    .into_iter()
                    .map(|ast| ast.bind(py).borrow().inner().clone()),
            )
            .into_iter()
            .map(|handle| New::from_rust(AstEntityHandle::Atom(handle)))
            .collect()
    }

    fn add_bond(
        &mut self,
        py: Python<'_>,
        first: HandleLike,
        second: HandleLike,
        ast: Py<BondAst>,
    ) -> New {
        New::from_rust(AstEntityHandle::Bond(self.0.add_bond(
            first.to_atom_handle(),
            second.to_atom_handle(),
            ast.bind(py).borrow().inner().clone(),
        )))
    }

    fn add_bonds(&mut self, py: Python<'_>, bonds: Vec<BondAddition>) -> Vec<New> {
        self.0
            .add_bonds(bonds.into_iter().map(|((first, second), ast)| AstAddBond {
                endpoints: [first.to_atom_handle(), second.to_atom_handle()],
                ast: ast.bind(py).borrow().inner().clone(),
            }))
            .into_iter()
            .map(|handle| New::from_rust(AstEntityHandle::Bond(handle)))
            .collect()
    }

    fn add_dative_bond(
        &mut self,
        py: Python<'_>,
        atoms: Vec<HandleLike>,
        ast: Py<DativeBondAst>,
    ) -> New {
        New::from_rust(AstEntityHandle::DativeBond(self.0.add_dative_bond(
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
            .map(|handle| New::from_rust(AstEntityHandle::DativeBond(handle)))
            .collect()
    }

    fn add_aromatic_system(
        &mut self,
        py: Python<'_>,
        atoms: Vec<HandleLike>,
        ast: Py<AromaticSystemAst>,
    ) -> New {
        New::from_rust(AstEntityHandle::AromaticSystem(self.0.add_aromatic_system(
            atoms.iter().map(HandleLike::to_atom_handle).collect(),
            ast.bind(py).borrow().inner().clone(),
        )))
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
            .map(|handle| New::from_rust(AstEntityHandle::AromaticSystem(handle)))
            .collect()
    }

    fn add_multicenter_bond(
        &mut self,
        py: Python<'_>,
        atoms: Vec<HandleLike>,
        ast: Py<MulticenterBondAst>,
    ) -> New {
        New::from_rust(AstEntityHandle::MulticenterBond(
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
            .map(|handle| New::from_rust(AstEntityHandle::MulticenterBond(handle)))
            .collect()
    }

    fn add_noncovalent_bond(
        &mut self,
        py: Python<'_>,
        atoms: (HandleLike, HandleLike),
        ast: Py<NoncovalentBondAst>,
    ) -> New {
        New::from_rust(AstEntityHandle::NoncovalentBond(
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
            .map(|handle| New::from_rust(AstEntityHandle::NoncovalentBond(handle)))
            .collect()
    }

    fn add_stereo_atom(
        &mut self,
        py: Python<'_>,
        site: HandleLike,
        ligands: Vec<StereoLigandInput>,
        ast: Py<StereoAtomAst>,
    ) -> New {
        New::from_rust(AstEntityHandle::StereoAtom(
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
            .map(|handle| New::from_rust(AstEntityHandle::StereoAtom(handle)))
            .collect()
    }

    fn add_stereo_bond(
        &mut self,
        py: Python<'_>,
        site: HandleLike,
        ligands: Vec<StereoLigandInput>,
        ast: Py<StereoBondAst>,
    ) -> New {
        New::from_rust(AstEntityHandle::StereoBond(
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
            .map(|handle| New::from_rust(AstEntityHandle::StereoBond(handle)))
            .collect()
    }
}

impl Edits {
    pub(crate) fn from_rust(edits: AstEdits) -> Self {
        Self(edits)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AromaticSystemAst as AstAromaticSystemAst,
        AromaticSystemFieldChange as AstAromaticSystemFieldChange, AtomAst as AstAtomAst,
        AtomFieldChange as AstAtomFieldChange, BondAst as AstBondAst,
        BondFieldChange as AstBondFieldChange, Constraint as AstConstraint,
        DativeBondAst as AstDativeBondAst, DativeBondFieldChange as AstDativeBondFieldChange,
        MoleculeConstraint as AstMoleculeConstraint, MulticenterBondAst as AstMulticenterBondAst,
        MulticenterBondFieldChange as AstMulticenterBondFieldChange,
        NoncovalentBondAst as AstNoncovalentBondAst,
        NoncovalentBondFieldChange as AstNoncovalentBondFieldChange,
        NoncovalentBondKind as AstNoncovalentBondKind,
        NoncovalentBondKindAst as AstNoncovalentBondKindAst, StereoAtomAst as AstStereoAtomAst,
        StereoAtomFieldChange as AstStereoAtomFieldChange, StereoBondAst as AstStereoBondAst,
        StereoBondFieldChange as AstStereoBondFieldChange,
        StereoConfigurationAst as AstStereoConfigurationAst, StereoKind as AstStereoKind,
        StereoLigandKind as AstStereoLigandKind, ValueAst as AstValueAst,
    };

    use super::*;

    #[rstest]
    #[case::id(HandleLike::Id(7), AstAtomHandle::Id(AstAtomId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstAtomHandle::New(7))]
    fn test_handle_like_to_atom_handle(#[case] input: HandleLike, #[case] expected: AstAtomHandle) {
        assert_eq!(input.to_atom_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), AstBondHandle::Id(AstBondId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstBondHandle::New(7))]
    fn test_handle_like_to_bond_handle(#[case] input: HandleLike, #[case] expected: AstBondHandle) {
        assert_eq!(input.to_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), AstDativeBondHandle::Id(AstDativeBondId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstDativeBondHandle::New(7))]
    fn test_handle_like_to_dative_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: AstDativeBondHandle,
    ) {
        assert_eq!(input.to_dative_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), AstAromaticSystemHandle::Id(AstAromaticSystemId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstAromaticSystemHandle::New(7))]
    fn test_handle_like_to_aromatic_system_handle(
        #[case] input: HandleLike,
        #[case] expected: AstAromaticSystemHandle,
    ) {
        assert_eq!(input.to_aromatic_system_handle(), expected);
    }

    #[rstest]
    #[case::id(
        HandleLike::Id(7),
        AstMulticenterBondHandle::Id(AstMulticenterBondId(7))
    )]
    #[case::new(HandleLike::New(New { index: 7 }), AstMulticenterBondHandle::New(7))]
    fn test_handle_like_to_multicenter_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: AstMulticenterBondHandle,
    ) {
        assert_eq!(input.to_multicenter_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(
        HandleLike::Id(7),
        AstNoncovalentBondHandle::Id(AstNoncovalentBondId(7))
    )]
    #[case::new(HandleLike::New(New { index: 7 }), AstNoncovalentBondHandle::New(7))]
    fn test_handle_like_to_noncovalent_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: AstNoncovalentBondHandle,
    ) {
        assert_eq!(input.to_noncovalent_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), AstStereoAtomHandle::Id(AstStereoAtomId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstStereoAtomHandle::New(7))]
    fn test_handle_like_to_stereo_atom_handle(
        #[case] input: HandleLike,
        #[case] expected: AstStereoAtomHandle,
    ) {
        assert_eq!(input.to_stereo_atom_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), AstStereoBondHandle::Id(AstStereoBondId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstStereoBondHandle::New(7))]
    fn test_handle_like_to_stereo_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: AstStereoBondHandle,
    ) {
        assert_eq!(input.to_stereo_bond_handle(), expected);
    }

    #[rstest]
    #[case::identity(AstConstraintEdit::from(AstConstraint::Molecule(
        AstMoleculeConstraint::Connected { atoms: None },
    )))]
    fn test_constraint_edit_roundtrip(#[case] expected: AstConstraintEdit) {
        let constraint = ConstraintEdit::from_rust(&expected);

        assert_eq!(constraint.to_rust(), expected);
    }

    #[rstest]
    #[case::inventory(vec![
        AstEdit::AddAtoms { atoms: vec![AstAtomAst::default()] },
        AstEdit::AddBonds { bonds: vec![AstAddBond {
            endpoints: [AstAtomHandle::Id(AstAtomId(0)), AstAtomHandle::New(0)],
            ast: AstBondAst::default(),
        }] },
        AstEdit::RemoveTopology {
            atoms: vec![AstAtomHandle::New(0)],
            bonds: vec![AstBondHandle::Id(AstBondId(1))],
        },
        AstEdit::ModifyAtomField {
            id: AstAtomHandle::Id(AstAtomId(0)),
            change: AstAtomFieldChange::Charge {
                old: AstValueAst::Lit(0),
                new: AstValueAst::Lit(1),
            },
        },
        AstEdit::ModifyBondField {
            id: AstBondHandle::New(0),
            change: AstBondFieldChange::Order {
                old: AstValueAst::Lit(1),
                new: AstValueAst::Lit(2),
            },
        },
        AstEdit::AddDativeBond {
            atoms: vec![AstAtomHandle::Id(AstAtomId(0)), AstAtomHandle::New(0)],
            ast: AstDativeBondAst::default(),
        },
        AstEdit::RemoveDativeBonds { removes: vec![(
            AstDativeBondHandle::New(0),
            vec![AstAtomHandle::Id(AstAtomId(0))],
            AstDativeBondAst::default(),
        )] },
        AstEdit::ModifyDativeBondField {
            id: AstDativeBondHandle::Id(AstDativeBondId(0)),
            change: AstDativeBondFieldChange::Order {
                old: AstValueAst::Lit(1),
                new: AstValueAst::Lit(2),
            },
        },
        AstEdit::AddAromaticSystem {
            atoms: vec![AstAtomHandle::New(0)],
            ast: AstAromaticSystemAst::default(),
        },
        AstEdit::RemoveAromaticSystems { removes: vec![(
            AstAromaticSystemHandle::Id(AstAromaticSystemId(0)),
            vec![AstAtomHandle::New(0)],
            AstAromaticSystemAst::default(),
        )] },
        AstEdit::ModifyAromaticSystemField {
            id: AstAromaticSystemHandle::New(0),
            change: AstAromaticSystemFieldChange::Charge {
                old: AstValueAst::Lit(0),
                new: AstValueAst::Lit(-1),
            },
        },
        AstEdit::AddMulticenterBond {
            atoms: vec![AstAtomHandle::Id(AstAtomId(0)), AstAtomHandle::New(0)],
            ast: AstMulticenterBondAst::default(),
        },
        AstEdit::RemoveMulticenterBonds { removes: vec![(
            AstMulticenterBondHandle::New(0),
            vec![AstAtomHandle::Id(AstAtomId(0))],
            AstMulticenterBondAst::default(),
        )] },
        AstEdit::ModifyMulticenterBondField {
            id: AstMulticenterBondHandle::Id(AstMulticenterBondId(0)),
            change: AstMulticenterBondFieldChange::Charge {
                old: AstValueAst::Lit(0),
                new: AstValueAst::Lit(1),
            },
        },
        AstEdit::AddNoncovalentBond {
            atoms: [AstAtomHandle::Id(AstAtomId(0)), AstAtomHandle::New(0)],
            ast: AstNoncovalentBondAst::default(),
        },
        AstEdit::RemoveNoncovalentBonds { removes: vec![(
            AstNoncovalentBondHandle::New(0),
            [AstAtomHandle::Id(AstAtomId(0)), AstAtomHandle::New(0)],
            AstNoncovalentBondAst::default(),
        )] },
        AstEdit::ModifyNoncovalentBondField {
            id: AstNoncovalentBondHandle::Id(AstNoncovalentBondId(0)),
            change: AstNoncovalentBondFieldChange::Kind {
                old: AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond),
                new: AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::Ionic),
            },
        },
        AstEdit::AddStereoAtom {
            site: AstAtomHandle::New(0),
            ligands: vec![(AstAtomHandle::Id(AstAtomId(0)), AstStereoLigandKind::Atom)],
            ast: AstStereoAtomAst::new(AstStereoKind::Tetrahedral, 0_u32),
        },
        AstEdit::RemoveStereoAtoms { removes: vec![(
            AstStereoAtomHandle::Id(AstStereoAtomId(0)),
            AstAtomHandle::New(0),
            vec![(AstAtomHandle::Id(AstAtomId(0)), AstStereoLigandKind::ImplicitHydrogen)],
            AstStereoAtomAst::new(AstStereoKind::Tetrahedral, 1_u32),
        )] },
        AstEdit::ModifyStereoAtomField {
            id: AstStereoAtomHandle::New(0),
            change: AstStereoAtomFieldChange::Configuration {
                old: AstStereoConfigurationAst::kinded(AstStereoKind::Tetrahedral, 0_u32),
                new: AstStereoConfigurationAst::kinded(AstStereoKind::Tetrahedral, 1_u32),
            },
        },
        AstEdit::AddStereoBond {
            site: AstBondHandle::Id(AstBondId(0)),
            ligands: vec![(AstAtomHandle::New(0), AstStereoLigandKind::LonePair)],
            ast: AstStereoBondAst::new(AstStereoKind::CisTrans, 0_u32),
        },
        AstEdit::RemoveStereoBonds { removes: vec![(
            AstStereoBondHandle::New(0),
            AstBondHandle::Id(AstBondId(0)),
            vec![(AstAtomHandle::New(0), AstStereoLigandKind::Atom)],
            AstStereoBondAst::new(AstStereoKind::CisTrans, 1_u32),
        )] },
        AstEdit::ModifyStereoBondField {
            id: AstStereoBondHandle::Id(AstStereoBondId(0)),
            change: AstStereoBondFieldChange::Configuration {
                old: AstStereoConfigurationAst::kinded(AstStereoKind::CisTrans, 0_u32),
                new: AstStereoConfigurationAst::kinded(AstStereoKind::CisTrans, 1_u32),
            },
        },
        AstEdit::ModifyAtomConstraint {
            id: AstAtomHandle::Id(AstAtomId(0)), old: None, new: None,
        },
        AstEdit::ModifyBondConstraint {
            id: AstBondHandle::New(0), old: None, new: None,
        },
        AstEdit::ModifyDativeBondConstraint {
            id: AstDativeBondHandle::Id(AstDativeBondId(0)), old: None, new: None,
        },
        AstEdit::ModifyAromaticSystemConstraint {
            id: AstAromaticSystemHandle::New(0), old: None, new: None,
        },
        AstEdit::ModifyMulticenterBondConstraint {
            id: AstMulticenterBondHandle::Id(AstMulticenterBondId(0)), old: None, new: None,
        },
        AstEdit::ModifyNoncovalentBondConstraint {
            id: AstNoncovalentBondHandle::New(0), old: None, new: None,
        },
        AstEdit::ModifyStereoAtomConstraint {
            id: AstStereoAtomHandle::Id(AstStereoAtomId(0)),
            kind: Some(AstStereoKind::Tetrahedral), old: None, new: None,
        },
        AstEdit::ModifyStereoBondConstraint {
            id: AstStereoBondHandle::New(0),
            kind: Some(AstStereoKind::CisTrans), old: None, new: None,
        },
        AstEdit::AddMoleculeConstraint {
            constraint: AstConstraintEdit::from(AstConstraint::Molecule(
                AstMoleculeConstraint::Connected { atoms: None },
            )),
        },
        AstEdit::RemoveMoleculeConstraint {
            constraint: AstConstraintEdit::from(AstConstraint::Molecule(
                AstMoleculeConstraint::Connected { atoms: Some(vec![AstAtomId(0)]) },
            )),
        },
    ])]
    fn test_edit_roundtrip(#[case] edits: Vec<AstEdit>) {
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
