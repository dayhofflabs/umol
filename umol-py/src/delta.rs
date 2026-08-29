//! Python bindings for resolved molecule deltas and their field-change payloads.

use std::vec::IntoIter;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use umol_graph_ir::ir::{
    AromaticSystemDelta as GraphIrAromaticSystemDelta,
    AromaticSystemFieldChange as GraphIrAromaticSystemFieldChange,
    AromaticSystemId as GraphIrAromaticSystemId, AtomDelta as GraphIrAtomDelta,
    AtomFieldChange as GraphIrAtomFieldChange, AtomId as GraphIrAtomId,
    BondDelta as GraphIrBondDelta, BondFieldChange as GraphIrBondFieldChange,
    BondId as GraphIrBondId, ConstraintDelta as GraphIrConstraintDelta,
    DativeBondDelta as GraphIrDativeBondDelta,
    DativeBondFieldChange as GraphIrDativeBondFieldChange, DativeBondId as GraphIrDativeBondId,
    Delta as GraphIrDelta, Deltas as GraphIrDeltas,
    MulticenterBondDelta as GraphIrMulticenterBondDelta,
    MulticenterBondFieldChange as GraphIrMulticenterBondFieldChange,
    MulticenterBondId as GraphIrMulticenterBondId,
    NoncovalentBondDelta as GraphIrNoncovalentBondDelta,
    NoncovalentBondFieldChange as GraphIrNoncovalentBondFieldChange,
    NoncovalentBondId as GraphIrNoncovalentBondId, StereoAtomDelta as GraphIrStereoAtomDelta,
    StereoAtomFieldChange as GraphIrStereoAtomFieldChange, StereoAtomId as GraphIrStereoAtomId,
    StereoBondDelta as GraphIrStereoBondDelta,
    StereoBondFieldChange as GraphIrStereoBondFieldChange, StereoBondId as GraphIrStereoBondId,
};
#[cfg(test)]
use umol_graph_ir::ir::{
    AromaticSystemForm as GraphIrAromaticSystemForm, AtomForm as GraphIrAtomForm,
    BondForm as GraphIrBondForm, Constraint as GraphIrConstraint,
    DativeBondForm as GraphIrDativeBondForm, MulticenterBondForm as GraphIrMulticenterBondForm,
    NoncovalentBondForm as GraphIrNoncovalentBondForm, StereoAtomForm as GraphIrStereoAtomForm,
    StereoBondForm as GraphIrStereoBondForm,
};

use crate::aromatic::AromaticSystemForm;
use crate::atom::{AtomForm, ElementForm, IsotopeMassForm};
use crate::bond::BondForm;
use crate::constraint::aromatic::AromaticSystemConstraintForm;
use crate::constraint::atom::AtomConstraintForm;
use crate::constraint::bond::BondConstraintForm;
use crate::constraint::dative::DativeBondConstraintForm;
use crate::constraint::molecule::Constraint;
use crate::constraint::multicenter::MulticenterBondConstraintForm;
use crate::constraint::noncovalent::NoncovalentBondConstraintForm;
use crate::constraint::stereo::{StereoAtomConstraintForm, StereoBondConstraintForm};
use crate::convert::{into_py_variant, variant_repr};
use crate::dative::DativeBondForm;
use crate::electrons::ElectronCountsForm;
use crate::entity::Readonly;
use crate::lattice::impl_py_normalize;
use crate::multicenter::MulticenterBondForm;
use crate::noncovalent::{NoncovalentBondForm, NoncovalentBondKindForm};
use crate::num::NumForm;
use crate::spin::UnpairedElectronsForm;
use crate::stereo::{
    StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoKind, StereoLigand,
};

/// Render a named old/new complex-enum variant using the child objects' reprs.
fn field_change_repr(obj: &Bound<'_, PyAny>, type_name: &str, variant: &str) -> PyResult<String> {
    let old = obj.getattr("old")?.repr()?.extract::<String>()?;
    let new = obj.getattr("new")?.repr()?.extract::<String>()?;
    Ok(format!("{type_name}.{variant}(old={old}, new={new})"))
}

/// Render a named complex-enum variant from its declared field names.
fn entity_delta_repr(
    obj: &Bound<'_, PyAny>,
    type_name: &str,
    variant: &str,
    fields: &[&str],
) -> PyResult<String> {
    let mut parts = Vec::with_capacity(fields.len());
    for field in fields {
        let value = obj.getattr(*field)?.repr()?.extract::<String>()?;
        parts.push(format!("{field}={value}"));
    }
    Ok(format!("{type_name}.{variant}({})", parts.join(", ")))
}

macro_rules! field_change {
    (
        $(#[$meta:meta])*
        $name:ident {
            $(
                $variant:ident($value:ty)
            ),+ $(,)?
        }
    ) => {
        $(#[$meta])*
        #[pyclass(frozen)]
        pub enum $name {
            $(
                $variant {
                    old: Py<$value>,
                    new: Py<$value>,
                },
            )+
        }

        #[pymethods]
        impl $name {
            fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
                self.to_rust(py) == other.to_rust(py)
            }

            fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
                let variant = match &*slf.bind(py).borrow() {
                    $(Self::$variant { .. } => stringify!($variant),)+
                };
                field_change_repr(
                    slf.bind(py).as_any(),
                    stringify!($name),
                    variant,
                )
            }

            /// Return the same field change with its old and new values exchanged.
            fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
                into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
            }
        }

    };
}

field_change! {
    /// An atom attribute change carrying the field's old and new forms.
    AtomFieldChange {
        Element(ElementForm),
        IsotopeMass(IsotopeMassForm),
        Charge(NumForm),
        ImplicitHydrogens(NumForm),
        LonePairs(NumForm),
        UnpairedElectrons(UnpairedElectronsForm),
    }
}

field_change! {
    /// A covalent-bond attribute change carrying the field's old and new forms.
    BondFieldChange {
        Order(NumForm),
        Charge(NumForm),
        UnpairedElectrons(UnpairedElectronsForm),
    }
}

field_change! {
    /// A dative-bond attribute change carrying the field's old and new forms.
    DativeBondFieldChange {
        Order(NumForm),
    }
}

field_change! {
    /// An aromatic-system attribute change carrying the field's old and new forms.
    AromaticSystemFieldChange {
        Electrons(ElectronCountsForm),
        Charge(NumForm),
        UnpairedElectrons(UnpairedElectronsForm),
    }
}

field_change! {
    /// A multicenter-bond attribute change carrying the field's old and new forms.
    MulticenterBondFieldChange {
        Electrons(ElectronCountsForm),
        Charge(NumForm),
        UnpairedElectrons(UnpairedElectronsForm),
    }
}

field_change! {
    /// A noncovalent-bond kind change carrying the field's old and new forms.
    NoncovalentBondFieldChange {
        Kind(NoncovalentBondKindForm),
    }
}

field_change! {
    /// A stereo-atom configuration change carrying the field's old and new forms.
    StereoAtomFieldChange {
        Configuration(StereoConfigurationForm),
    }
}

field_change! {
    /// A stereo-bond configuration change carrying the field's old and new forms.
    StereoBondFieldChange {
        Configuration(StereoConfigurationForm),
    }
}

impl AtomFieldChange {
    pub(crate) fn from_rust(py: Python<'_>, change: &GraphIrAtomFieldChange) -> PyResult<Self> {
        Ok(match change {
            GraphIrAtomFieldChange::Element { old, new } => Self::Element {
                old: into_py_variant(py, ElementForm::from_rust(old))?,
                new: into_py_variant(py, ElementForm::from_rust(new))?,
            },
            GraphIrAtomFieldChange::IsotopeMass { old, new } => Self::IsotopeMass {
                old: into_py_variant(py, IsotopeMassForm::from_rust(old))?,
                new: into_py_variant(py, IsotopeMassForm::from_rust(new))?,
            },
            GraphIrAtomFieldChange::Charge { old, new } => Self::Charge {
                old: into_py_variant(py, NumForm::from_rust(py, old)?)?,
                new: into_py_variant(py, NumForm::from_rust(py, new)?)?,
            },
            GraphIrAtomFieldChange::ImplicitHydrogens { old, new } => Self::ImplicitHydrogens {
                old: into_py_variant(py, NumForm::from_rust(py, old)?)?,
                new: into_py_variant(py, NumForm::from_rust(py, new)?)?,
            },
            GraphIrAtomFieldChange::LonePairs { old, new } => Self::LonePairs {
                old: into_py_variant(py, NumForm::from_rust(py, old)?)?,
                new: into_py_variant(py, NumForm::from_rust(py, new)?)?,
            },
            GraphIrAtomFieldChange::UnpairedElectrons { old, new } => Self::UnpairedElectrons {
                old: Py::new(py, UnpairedElectronsForm::from_rust(py, old)?)?,
                new: Py::new(py, UnpairedElectronsForm::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrAtomFieldChange {
        match self {
            Self::Element { old, new } => GraphIrAtomFieldChange::Element {
                old: old.bind(py).borrow().to_rust(),
                new: new.bind(py).borrow().to_rust(),
            },
            Self::IsotopeMass { old, new } => GraphIrAtomFieldChange::IsotopeMass {
                old: old.bind(py).borrow().to_rust(),
                new: new.bind(py).borrow().to_rust(),
            },
            Self::Charge { old, new } => GraphIrAtomFieldChange::Charge {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::ImplicitHydrogens { old, new } => GraphIrAtomFieldChange::ImplicitHydrogens {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::LonePairs { old, new } => GraphIrAtomFieldChange::LonePairs {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::UnpairedElectrons { old, new } => GraphIrAtomFieldChange::UnpairedElectrons {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

impl BondFieldChange {
    pub(crate) fn from_rust(py: Python<'_>, change: &GraphIrBondFieldChange) -> PyResult<Self> {
        Ok(match change {
            GraphIrBondFieldChange::Order { old, new } => Self::Order {
                old: into_py_variant(py, NumForm::from_rust(py, old)?)?,
                new: into_py_variant(py, NumForm::from_rust(py, new)?)?,
            },
            GraphIrBondFieldChange::Charge { old, new } => Self::Charge {
                old: into_py_variant(py, NumForm::from_rust(py, old)?)?,
                new: into_py_variant(py, NumForm::from_rust(py, new)?)?,
            },
            GraphIrBondFieldChange::UnpairedElectrons { old, new } => Self::UnpairedElectrons {
                old: Py::new(py, UnpairedElectronsForm::from_rust(py, old)?)?,
                new: Py::new(py, UnpairedElectronsForm::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrBondFieldChange {
        match self {
            Self::Order { old, new } => GraphIrBondFieldChange::Order {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::Charge { old, new } => GraphIrBondFieldChange::Charge {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::UnpairedElectrons { old, new } => GraphIrBondFieldChange::UnpairedElectrons {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

impl DativeBondFieldChange {
    pub(crate) fn from_rust(
        py: Python<'_>,
        change: &GraphIrDativeBondFieldChange,
    ) -> PyResult<Self> {
        Ok(match change {
            GraphIrDativeBondFieldChange::Order { old, new } => Self::Order {
                old: into_py_variant(py, NumForm::from_rust(py, old)?)?,
                new: into_py_variant(py, NumForm::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrDativeBondFieldChange {
        match self {
            Self::Order { old, new } => GraphIrDativeBondFieldChange::Order {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

impl AromaticSystemFieldChange {
    pub(crate) fn from_rust(
        py: Python<'_>,
        change: &GraphIrAromaticSystemFieldChange,
    ) -> PyResult<Self> {
        Ok(match change {
            GraphIrAromaticSystemFieldChange::Electrons { old, new } => Self::Electrons {
                old: into_py_variant(py, ElectronCountsForm::from_rust(old))?,
                new: into_py_variant(py, ElectronCountsForm::from_rust(new))?,
            },
            GraphIrAromaticSystemFieldChange::Charge { old, new } => Self::Charge {
                old: into_py_variant(py, NumForm::from_rust(py, old)?)?,
                new: into_py_variant(py, NumForm::from_rust(py, new)?)?,
            },
            GraphIrAromaticSystemFieldChange::UnpairedElectrons { old, new } => {
                Self::UnpairedElectrons {
                    old: Py::new(py, UnpairedElectronsForm::from_rust(py, old)?)?,
                    new: Py::new(py, UnpairedElectronsForm::from_rust(py, new)?)?,
                }
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrAromaticSystemFieldChange {
        match self {
            Self::Electrons { old, new } => GraphIrAromaticSystemFieldChange::Electrons {
                old: old.bind(py).borrow().to_rust(),
                new: new.bind(py).borrow().to_rust(),
            },
            Self::Charge { old, new } => GraphIrAromaticSystemFieldChange::Charge {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::UnpairedElectrons { old, new } => {
                GraphIrAromaticSystemFieldChange::UnpairedElectrons {
                    old: old.bind(py).borrow().to_rust(py),
                    new: new.bind(py).borrow().to_rust(py),
                }
            }
        }
    }
}

impl MulticenterBondFieldChange {
    pub(crate) fn from_rust(
        py: Python<'_>,
        change: &GraphIrMulticenterBondFieldChange,
    ) -> PyResult<Self> {
        Ok(match change {
            GraphIrMulticenterBondFieldChange::Electrons { old, new } => Self::Electrons {
                old: into_py_variant(py, ElectronCountsForm::from_rust(old))?,
                new: into_py_variant(py, ElectronCountsForm::from_rust(new))?,
            },
            GraphIrMulticenterBondFieldChange::Charge { old, new } => Self::Charge {
                old: into_py_variant(py, NumForm::from_rust(py, old)?)?,
                new: into_py_variant(py, NumForm::from_rust(py, new)?)?,
            },
            GraphIrMulticenterBondFieldChange::UnpairedElectrons { old, new } => {
                Self::UnpairedElectrons {
                    old: Py::new(py, UnpairedElectronsForm::from_rust(py, old)?)?,
                    new: Py::new(py, UnpairedElectronsForm::from_rust(py, new)?)?,
                }
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrMulticenterBondFieldChange {
        match self {
            Self::Electrons { old, new } => GraphIrMulticenterBondFieldChange::Electrons {
                old: old.bind(py).borrow().to_rust(),
                new: new.bind(py).borrow().to_rust(),
            },
            Self::Charge { old, new } => GraphIrMulticenterBondFieldChange::Charge {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::UnpairedElectrons { old, new } => {
                GraphIrMulticenterBondFieldChange::UnpairedElectrons {
                    old: old.bind(py).borrow().to_rust(py),
                    new: new.bind(py).borrow().to_rust(py),
                }
            }
        }
    }
}

impl NoncovalentBondFieldChange {
    pub(crate) fn from_rust(
        py: Python<'_>,
        change: &GraphIrNoncovalentBondFieldChange,
    ) -> PyResult<Self> {
        Ok(match change {
            GraphIrNoncovalentBondFieldChange::Kind { old, new } => Self::Kind {
                old: into_py_variant(py, NoncovalentBondKindForm::from_rust(old))?,
                new: into_py_variant(py, NoncovalentBondKindForm::from_rust(new))?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrNoncovalentBondFieldChange {
        match self {
            Self::Kind { old, new } => GraphIrNoncovalentBondFieldChange::Kind {
                old: old.bind(py).borrow().to_rust(),
                new: new.bind(py).borrow().to_rust(),
            },
        }
    }
}

impl StereoAtomFieldChange {
    pub(crate) fn from_rust(
        py: Python<'_>,
        change: &GraphIrStereoAtomFieldChange,
    ) -> PyResult<Self> {
        Ok(match change {
            GraphIrStereoAtomFieldChange::Configuration { old, new } => Self::Configuration {
                old: into_py_variant(py, StereoConfigurationForm::from_rust(py, old)?)?,
                new: into_py_variant(py, StereoConfigurationForm::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrStereoAtomFieldChange {
        match self {
            Self::Configuration { old, new } => GraphIrStereoAtomFieldChange::Configuration {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

impl StereoBondFieldChange {
    pub(crate) fn from_rust(
        py: Python<'_>,
        change: &GraphIrStereoBondFieldChange,
    ) -> PyResult<Self> {
        Ok(match change {
            GraphIrStereoBondFieldChange::Configuration { old, new } => Self::Configuration {
                old: into_py_variant(py, StereoConfigurationForm::from_rust(py, old)?)?,
                new: into_py_variant(py, StereoConfigurationForm::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrStereoBondFieldChange {
        match self {
            Self::Configuration { old, new } => GraphIrStereoBondFieldChange::Configuration {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

/// A resolved edit to one atom.
#[pyclass(frozen)]
pub enum AtomDelta {
    Add {
        id: u32,
        attributes: Readonly<AtomForm>,
    },
    Remove {
        id: u32,
        attributes: Readonly<AtomForm>,
    },
    ModifyField {
        id: u32,
        change: Py<AtomFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<AtomConstraintForm>>,
        new: Option<Py<AtomConstraintForm>>,
    },
}

#[pymethods]
impl AtomDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "attributes"]),
            Self::Remove { .. } => ("Remove", &["id", "attributes"]),
            Self::ModifyField { .. } => ("ModifyField", &["id", "change"]),
            Self::ModifyConstraint { .. } => ("ModifyConstraint", &["id", "old", "new"]),
        };
        entity_delta_repr(slf.bind(py).as_any(), "AtomDelta", variant, fields)
    }

    /// Return the inverse resolved edit.
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl AtomDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &GraphIrAtomDelta) -> PyResult<Self> {
        Ok(match delta {
            GraphIrAtomDelta::Add { id, attributes } => Self::Add {
                id: id.0,
                attributes: Readonly::<AtomForm>::from_rust(py, attributes)?,
            },
            GraphIrAtomDelta::Remove { id, attributes } => Self::Remove {
                id: id.0,
                attributes: Readonly::<AtomForm>::from_rust(py, attributes)?,
            },
            GraphIrAtomDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, AtomFieldChange::from_rust(py, change)?)?,
            },
            GraphIrAtomDelta::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id: id.0,
                old: old
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, AtomConstraintForm::from_rust(py, constraint)?)
                    })
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, AtomConstraintForm::from_rust(py, constraint)?)
                    })
                    .transpose()?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrAtomDelta {
        match self {
            Self::Add { id, attributes } => GraphIrAtomDelta::Add {
                id: GraphIrAtomId(*id),
                attributes: attributes.to_rust(py),
            },
            Self::Remove { id, attributes } => GraphIrAtomDelta::Remove {
                id: GraphIrAtomId(*id),
                attributes: attributes.to_rust(py),
            },
            Self::ModifyField { id, change } => GraphIrAtomDelta::ModifyField {
                id: GraphIrAtomId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => GraphIrAtomDelta::ModifyConstraint {
                id: GraphIrAtomId(*id),
                old: old
                    .as_ref()
                    .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                new: new
                    .as_ref()
                    .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
            },
        }
    }
}

/// A resolved edit to one covalent bond.
#[pyclass(frozen)]
pub enum BondDelta {
    Add {
        id: u32,
        atoms: (u32, u32),
        attributes: Readonly<BondForm>,
    },
    Remove {
        id: u32,
        atoms: (u32, u32),
        attributes: Readonly<BondForm>,
    },
    ModifyField {
        id: u32,
        change: Py<BondFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<BondConstraintForm>>,
        new: Option<Py<BondConstraintForm>>,
    },
}

#[pymethods]
impl BondDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "atoms", "attributes"]),
            Self::Remove { .. } => ("Remove", &["id", "atoms", "attributes"]),
            Self::ModifyField { .. } => ("ModifyField", &["id", "change"]),
            Self::ModifyConstraint { .. } => ("ModifyConstraint", &["id", "old", "new"]),
        };
        entity_delta_repr(slf.bind(py).as_any(), "BondDelta", variant, fields)
    }

    /// Return the inverse resolved edit.
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl BondDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &GraphIrBondDelta) -> PyResult<Self> {
        Ok(match delta {
            GraphIrBondDelta::Add {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id: id.0,
                atoms: (atoms[0].0, atoms[1].0),
                attributes: Readonly::<BondForm>::from_rust(py, attributes)?,
            },
            GraphIrBondDelta::Remove {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id: id.0,
                atoms: (atoms[0].0, atoms[1].0),
                attributes: Readonly::<BondForm>::from_rust(py, attributes)?,
            },
            GraphIrBondDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, BondFieldChange::from_rust(py, change)?)?,
            },
            GraphIrBondDelta::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id: id.0,
                old: old
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, BondConstraintForm::from_rust(py, constraint)?)
                    })
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, BondConstraintForm::from_rust(py, constraint)?)
                    })
                    .transpose()?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrBondDelta {
        match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => GraphIrBondDelta::Add {
                id: GraphIrBondId(*id),
                atoms: [GraphIrAtomId(atoms.0), GraphIrAtomId(atoms.1)],
                attributes: attributes.to_rust(py),
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => GraphIrBondDelta::Remove {
                id: GraphIrBondId(*id),
                atoms: [GraphIrAtomId(atoms.0), GraphIrAtomId(atoms.1)],
                attributes: attributes.to_rust(py),
            },
            Self::ModifyField { id, change } => GraphIrBondDelta::ModifyField {
                id: GraphIrBondId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => GraphIrBondDelta::ModifyConstraint {
                id: GraphIrBondId(*id),
                old: old
                    .as_ref()
                    .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                new: new
                    .as_ref()
                    .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
            },
        }
    }
}

/// A resolved edit to one dative bond.
#[pyclass(frozen)]
pub enum DativeBondDelta {
    Add {
        id: u32,
        donors: Vec<u32>,
        acceptor: u32,
        attributes: Readonly<DativeBondForm>,
    },
    Remove {
        id: u32,
        donors: Vec<u32>,
        acceptor: u32,
        attributes: Readonly<DativeBondForm>,
    },
    ModifyField {
        id: u32,
        change: Py<DativeBondFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<DativeBondConstraintForm>>,
        new: Option<Py<DativeBondConstraintForm>>,
    },
}

#[pymethods]
impl DativeBondDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "donors", "acceptor", "attributes"]),
            Self::Remove { .. } => ("Remove", &["id", "donors", "acceptor", "attributes"]),
            Self::ModifyField { .. } => ("ModifyField", &["id", "change"]),
            Self::ModifyConstraint { .. } => ("ModifyConstraint", &["id", "old", "new"]),
        };
        entity_delta_repr(slf.bind(py).as_any(), "DativeBondDelta", variant, fields)
    }

    /// Return the inverse resolved edit.
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl DativeBondDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &GraphIrDativeBondDelta) -> PyResult<Self> {
        Ok(match delta {
            GraphIrDativeBondDelta::Add {
                id,
                donors,
                acceptor,
                attributes,
            } => Self::Add {
                id: id.0,
                donors: donors.iter().map(|atom| atom.0).collect(),
                acceptor: acceptor.0,
                attributes: Readonly::<DativeBondForm>::from_rust(py, attributes)?,
            },
            GraphIrDativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                attributes,
            } => Self::Remove {
                id: id.0,
                donors: donors.iter().map(|atom| atom.0).collect(),
                acceptor: acceptor.0,
                attributes: Readonly::<DativeBondForm>::from_rust(py, attributes)?,
            },
            GraphIrDativeBondDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, DativeBondFieldChange::from_rust(py, change)?)?,
            },
            GraphIrDativeBondDelta::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id: id.0,
                old: old
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, DativeBondConstraintForm::from_rust(py, constraint)?)
                    })
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, DativeBondConstraintForm::from_rust(py, constraint)?)
                    })
                    .transpose()?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrDativeBondDelta {
        match self {
            Self::Add {
                id,
                donors,
                acceptor,
                attributes,
            } => GraphIrDativeBondDelta::Add {
                id: GraphIrDativeBondId(*id),
                donors: donors.iter().copied().map(GraphIrAtomId).collect(),
                acceptor: GraphIrAtomId(*acceptor),
                attributes: attributes.to_rust(py),
            },
            Self::Remove {
                id,
                donors,
                acceptor,
                attributes,
            } => GraphIrDativeBondDelta::Remove {
                id: GraphIrDativeBondId(*id),
                donors: donors.iter().copied().map(GraphIrAtomId).collect(),
                acceptor: GraphIrAtomId(*acceptor),
                attributes: attributes.to_rust(py),
            },
            Self::ModifyField { id, change } => GraphIrDativeBondDelta::ModifyField {
                id: GraphIrDativeBondId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => GraphIrDativeBondDelta::ModifyConstraint {
                id: GraphIrDativeBondId(*id),
                old: old
                    .as_ref()
                    .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                new: new
                    .as_ref()
                    .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
            },
        }
    }
}

/// A resolved edit to one aromatic system.
#[pyclass(frozen)]
pub enum AromaticSystemDelta {
    Add {
        id: u32,
        atoms: Vec<u32>,
        attributes: Readonly<AromaticSystemForm>,
    },
    Remove {
        id: u32,
        atoms: Vec<u32>,
        attributes: Readonly<AromaticSystemForm>,
    },
    ModifyField {
        id: u32,
        change: Py<AromaticSystemFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<AromaticSystemConstraintForm>>,
        new: Option<Py<AromaticSystemConstraintForm>>,
    },
}

#[pymethods]
impl AromaticSystemDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "atoms", "attributes"]),
            Self::Remove { .. } => ("Remove", &["id", "atoms", "attributes"]),
            Self::ModifyField { .. } => ("ModifyField", &["id", "change"]),
            Self::ModifyConstraint { .. } => ("ModifyConstraint", &["id", "old", "new"]),
        };
        entity_delta_repr(
            slf.bind(py).as_any(),
            "AromaticSystemDelta",
            variant,
            fields,
        )
    }

    /// Return the inverse resolved edit.
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl AromaticSystemDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &GraphIrAromaticSystemDelta) -> PyResult<Self> {
        Ok(match delta {
            GraphIrAromaticSystemDelta::Add {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id: id.0,
                atoms: atoms.iter().map(|atom| atom.0).collect(),
                attributes: Readonly::<AromaticSystemForm>::from_rust(py, attributes)?,
            },
            GraphIrAromaticSystemDelta::Remove {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id: id.0,
                atoms: atoms.iter().map(|atom| atom.0).collect(),
                attributes: Readonly::<AromaticSystemForm>::from_rust(py, attributes)?,
            },
            GraphIrAromaticSystemDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, AromaticSystemFieldChange::from_rust(py, change)?)?,
            },
            GraphIrAromaticSystemDelta::ModifyConstraint { id, old, new } => {
                Self::ModifyConstraint {
                    id: id.0,
                    old: old
                        .as_ref()
                        .map(|constraint| {
                            into_py_variant(
                                py,
                                AromaticSystemConstraintForm::from_rust(py, constraint)?,
                            )
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|constraint| {
                            into_py_variant(
                                py,
                                AromaticSystemConstraintForm::from_rust(py, constraint)?,
                            )
                        })
                        .transpose()?,
                }
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrAromaticSystemDelta {
        match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => GraphIrAromaticSystemDelta::Add {
                id: GraphIrAromaticSystemId(*id),
                atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                attributes: attributes.to_rust(py),
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => GraphIrAromaticSystemDelta::Remove {
                id: GraphIrAromaticSystemId(*id),
                atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                attributes: attributes.to_rust(py),
            },
            Self::ModifyField { id, change } => GraphIrAromaticSystemDelta::ModifyField {
                id: GraphIrAromaticSystemId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => {
                GraphIrAromaticSystemDelta::ModifyConstraint {
                    id: GraphIrAromaticSystemId(*id),
                    old: old
                        .as_ref()
                        .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                    new: new
                        .as_ref()
                        .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                }
            }
        }
    }
}

/// A resolved edit to one multicenter bond.
#[pyclass(frozen)]
pub enum MulticenterBondDelta {
    Add {
        id: u32,
        atoms: Vec<u32>,
        attributes: Readonly<MulticenterBondForm>,
    },
    Remove {
        id: u32,
        atoms: Vec<u32>,
        attributes: Readonly<MulticenterBondForm>,
    },
    ModifyField {
        id: u32,
        change: Py<MulticenterBondFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<MulticenterBondConstraintForm>>,
        new: Option<Py<MulticenterBondConstraintForm>>,
    },
}

#[pymethods]
impl MulticenterBondDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "atoms", "attributes"]),
            Self::Remove { .. } => ("Remove", &["id", "atoms", "attributes"]),
            Self::ModifyField { .. } => ("ModifyField", &["id", "change"]),
            Self::ModifyConstraint { .. } => ("ModifyConstraint", &["id", "old", "new"]),
        };
        entity_delta_repr(
            slf.bind(py).as_any(),
            "MulticenterBondDelta",
            variant,
            fields,
        )
    }

    /// Return the inverse resolved edit.
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl MulticenterBondDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &GraphIrMulticenterBondDelta) -> PyResult<Self> {
        Ok(match delta {
            GraphIrMulticenterBondDelta::Add {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id: id.0,
                atoms: atoms.iter().map(|atom| atom.0).collect(),
                attributes: Readonly::<MulticenterBondForm>::from_rust(py, attributes)?,
            },
            GraphIrMulticenterBondDelta::Remove {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id: id.0,
                atoms: atoms.iter().map(|atom| atom.0).collect(),
                attributes: Readonly::<MulticenterBondForm>::from_rust(py, attributes)?,
            },
            GraphIrMulticenterBondDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, MulticenterBondFieldChange::from_rust(py, change)?)?,
            },
            GraphIrMulticenterBondDelta::ModifyConstraint { id, old, new } => {
                Self::ModifyConstraint {
                    id: id.0,
                    old: old
                        .as_ref()
                        .map(|constraint| {
                            into_py_variant(
                                py,
                                MulticenterBondConstraintForm::from_rust(py, constraint)?,
                            )
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|constraint| {
                            into_py_variant(
                                py,
                                MulticenterBondConstraintForm::from_rust(py, constraint)?,
                            )
                        })
                        .transpose()?,
                }
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrMulticenterBondDelta {
        match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => GraphIrMulticenterBondDelta::Add {
                id: GraphIrMulticenterBondId(*id),
                atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                attributes: attributes.to_rust(py),
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => GraphIrMulticenterBondDelta::Remove {
                id: GraphIrMulticenterBondId(*id),
                atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                attributes: attributes.to_rust(py),
            },
            Self::ModifyField { id, change } => GraphIrMulticenterBondDelta::ModifyField {
                id: GraphIrMulticenterBondId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => {
                GraphIrMulticenterBondDelta::ModifyConstraint {
                    id: GraphIrMulticenterBondId(*id),
                    old: old
                        .as_ref()
                        .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                    new: new
                        .as_ref()
                        .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                }
            }
        }
    }
}

/// A resolved edit to one noncovalent bond.
#[pyclass(frozen)]
pub enum NoncovalentBondDelta {
    Add {
        id: u32,
        atoms: (u32, u32),
        attributes: Readonly<NoncovalentBondForm>,
    },
    Remove {
        id: u32,
        atoms: (u32, u32),
        attributes: Readonly<NoncovalentBondForm>,
    },
    ModifyField {
        id: u32,
        change: Py<NoncovalentBondFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<NoncovalentBondConstraintForm>>,
        new: Option<Py<NoncovalentBondConstraintForm>>,
    },
}

#[pymethods]
impl NoncovalentBondDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "atoms", "attributes"]),
            Self::Remove { .. } => ("Remove", &["id", "atoms", "attributes"]),
            Self::ModifyField { .. } => ("ModifyField", &["id", "change"]),
            Self::ModifyConstraint { .. } => ("ModifyConstraint", &["id", "old", "new"]),
        };
        entity_delta_repr(
            slf.bind(py).as_any(),
            "NoncovalentBondDelta",
            variant,
            fields,
        )
    }

    /// Return the inverse resolved edit.
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl NoncovalentBondDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &GraphIrNoncovalentBondDelta) -> PyResult<Self> {
        Ok(match delta {
            GraphIrNoncovalentBondDelta::Add {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id: id.0,
                atoms: (atoms[0].0, atoms[1].0),
                attributes: Readonly::<NoncovalentBondForm>::from_rust(py, attributes)?,
            },
            GraphIrNoncovalentBondDelta::Remove {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id: id.0,
                atoms: (atoms[0].0, atoms[1].0),
                attributes: Readonly::<NoncovalentBondForm>::from_rust(py, attributes)?,
            },
            GraphIrNoncovalentBondDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, NoncovalentBondFieldChange::from_rust(py, change)?)?,
            },
            GraphIrNoncovalentBondDelta::ModifyConstraint { id, old, new } => {
                Self::ModifyConstraint {
                    id: id.0,
                    old: old
                        .as_ref()
                        .map(|constraint| {
                            into_py_variant(
                                py,
                                NoncovalentBondConstraintForm::from_rust(py, constraint)?,
                            )
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|constraint| {
                            into_py_variant(
                                py,
                                NoncovalentBondConstraintForm::from_rust(py, constraint)?,
                            )
                        })
                        .transpose()?,
                }
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrNoncovalentBondDelta {
        match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => GraphIrNoncovalentBondDelta::Add {
                id: GraphIrNoncovalentBondId(*id),
                atoms: [GraphIrAtomId(atoms.0), GraphIrAtomId(atoms.1)],
                attributes: attributes.to_rust(py),
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => GraphIrNoncovalentBondDelta::Remove {
                id: GraphIrNoncovalentBondId(*id),
                atoms: [GraphIrAtomId(atoms.0), GraphIrAtomId(atoms.1)],
                attributes: attributes.to_rust(py),
            },
            Self::ModifyField { id, change } => GraphIrNoncovalentBondDelta::ModifyField {
                id: GraphIrNoncovalentBondId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => {
                GraphIrNoncovalentBondDelta::ModifyConstraint {
                    id: GraphIrNoncovalentBondId(*id),
                    old: old
                        .as_ref()
                        .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                    new: new
                        .as_ref()
                        .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                }
            }
        }
    }
}

/// A resolved edit to one atom-centered stereo element.
#[pyclass(frozen)]
pub enum StereoAtomDelta {
    Add {
        id: u32,
        site: u32,
        ligands: Vec<StereoLigand>,
        attributes: Readonly<StereoAtomForm>,
    },
    Remove {
        id: u32,
        site: u32,
        ligands: Vec<StereoLigand>,
        attributes: Readonly<StereoAtomForm>,
    },
    ModifyField {
        id: u32,
        change: Py<StereoAtomFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        kind: Option<StereoKind>,
        old: Option<Py<StereoAtomConstraintForm>>,
        new: Option<Py<StereoAtomConstraintForm>>,
    },
}

#[pymethods]
impl StereoAtomDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "site", "ligands", "attributes"]),
            Self::Remove { .. } => ("Remove", &["id", "site", "ligands", "attributes"]),
            Self::ModifyField { .. } => ("ModifyField", &["id", "change"]),
            Self::ModifyConstraint { .. } => ("ModifyConstraint", &["id", "kind", "old", "new"]),
        };
        entity_delta_repr(slf.bind(py).as_any(), "StereoAtomDelta", variant, fields)
    }

    /// Return the inverse resolved edit.
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl StereoAtomDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &GraphIrStereoAtomDelta) -> PyResult<Self> {
        Ok(match delta {
            GraphIrStereoAtomDelta::Add {
                id,
                site,
                ligands,
                attributes,
            } => Self::Add {
                id: id.0,
                site: site.0,
                ligands: ligands
                    .iter()
                    .copied()
                    .map(StereoLigand::from_rust)
                    .collect(),
                attributes: Readonly::<StereoAtomForm>::from_rust(py, attributes)?,
            },
            GraphIrStereoAtomDelta::Remove {
                id,
                site,
                ligands,
                attributes,
            } => Self::Remove {
                id: id.0,
                site: site.0,
                ligands: ligands
                    .iter()
                    .copied()
                    .map(StereoLigand::from_rust)
                    .collect(),
                attributes: Readonly::<StereoAtomForm>::from_rust(py, attributes)?,
            },
            GraphIrStereoAtomDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, StereoAtomFieldChange::from_rust(py, change)?)?,
            },
            GraphIrStereoAtomDelta::ModifyConstraint { id, kind, old, new } => {
                Self::ModifyConstraint {
                    id: id.0,
                    kind: kind.map(StereoKind::from_rust),
                    old: old
                        .as_ref()
                        .map(|constraint| {
                            into_py_variant(
                                py,
                                StereoAtomConstraintForm::from_rust(py, constraint)?,
                            )
                        })
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|constraint| {
                            into_py_variant(
                                py,
                                StereoAtomConstraintForm::from_rust(py, constraint)?,
                            )
                        })
                        .transpose()?,
                }
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrStereoAtomDelta {
        match self {
            Self::Add {
                id,
                site,
                ligands,
                attributes,
            } => GraphIrStereoAtomDelta::Add {
                id: GraphIrStereoAtomId(*id),
                site: GraphIrAtomId(*site),
                ligands: ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                attributes: attributes.to_rust(py),
            },
            Self::Remove {
                id,
                site,
                ligands,
                attributes,
            } => GraphIrStereoAtomDelta::Remove {
                id: GraphIrStereoAtomId(*id),
                site: GraphIrAtomId(*site),
                ligands: ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                attributes: attributes.to_rust(py),
            },
            Self::ModifyField { id, change } => GraphIrStereoAtomDelta::ModifyField {
                id: GraphIrStereoAtomId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, kind, old, new } => {
                GraphIrStereoAtomDelta::ModifyConstraint {
                    id: GraphIrStereoAtomId(*id),
                    kind: kind.map(StereoKind::to_rust),
                    old: old
                        .as_ref()
                        .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                    new: new
                        .as_ref()
                        .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                }
            }
        }
    }
}

/// A resolved edit to one bond-centered stereo element.
#[pyclass(frozen)]
pub enum StereoBondDelta {
    Add {
        id: u32,
        site: u32,
        ligands: Vec<StereoLigand>,
        attributes: Readonly<StereoBondForm>,
    },
    Remove {
        id: u32,
        site: u32,
        ligands: Vec<StereoLigand>,
        attributes: Readonly<StereoBondForm>,
    },
    ModifyField {
        id: u32,
        change: Py<StereoBondFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        kind: Option<StereoKind>,
        old: Option<Py<StereoBondConstraintForm>>,
        new: Option<Py<StereoBondConstraintForm>>,
    },
}

#[pymethods]
impl StereoBondDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }
    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "site", "ligands", "attributes"]),
            Self::Remove { .. } => ("Remove", &["id", "site", "ligands", "attributes"]),
            Self::ModifyField { .. } => ("ModifyField", &["id", "change"]),
            Self::ModifyConstraint { .. } => ("ModifyConstraint", &["id", "kind", "old", "new"]),
        };
        entity_delta_repr(slf.bind(py).as_any(), "StereoBondDelta", variant, fields)
    }
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl StereoBondDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &GraphIrStereoBondDelta) -> PyResult<Self> {
        Ok(match delta {
            GraphIrStereoBondDelta::Add {
                id,
                site,
                ligands,
                attributes,
            } => Self::Add {
                id: id.0,
                site: site.0,
                ligands: ligands
                    .iter()
                    .copied()
                    .map(StereoLigand::from_rust)
                    .collect(),
                attributes: Readonly::<StereoBondForm>::from_rust(py, attributes)?,
            },
            GraphIrStereoBondDelta::Remove {
                id,
                site,
                ligands,
                attributes,
            } => Self::Remove {
                id: id.0,
                site: site.0,
                ligands: ligands
                    .iter()
                    .copied()
                    .map(StereoLigand::from_rust)
                    .collect(),
                attributes: Readonly::<StereoBondForm>::from_rust(py, attributes)?,
            },
            GraphIrStereoBondDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, StereoBondFieldChange::from_rust(py, change)?)?,
            },
            GraphIrStereoBondDelta::ModifyConstraint { id, kind, old, new } => {
                Self::ModifyConstraint {
                    id: id.0,
                    kind: kind.map(StereoKind::from_rust),
                    old: old
                        .as_ref()
                        .map(|c| into_py_variant(py, StereoBondConstraintForm::from_rust(py, c)?))
                        .transpose()?,
                    new: new
                        .as_ref()
                        .map(|c| into_py_variant(py, StereoBondConstraintForm::from_rust(py, c)?))
                        .transpose()?,
                }
            }
        })
    }
    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrStereoBondDelta {
        match self {
            Self::Add {
                id,
                site,
                ligands,
                attributes,
            } => GraphIrStereoBondDelta::Add {
                id: GraphIrStereoBondId(*id),
                site: GraphIrBondId(*site),
                ligands: ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                attributes: attributes.to_rust(py),
            },
            Self::Remove {
                id,
                site,
                ligands,
                attributes,
            } => GraphIrStereoBondDelta::Remove {
                id: GraphIrStereoBondId(*id),
                site: GraphIrBondId(*site),
                ligands: ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                attributes: attributes.to_rust(py),
            },
            Self::ModifyField { id, change } => GraphIrStereoBondDelta::ModifyField {
                id: GraphIrStereoBondId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, kind, old, new } => {
                GraphIrStereoBondDelta::ModifyConstraint {
                    id: GraphIrStereoBondId(*id),
                    kind: kind.map(StereoKind::to_rust),
                    old: old.as_ref().map(|c| c.bind(py).borrow().to_rust(py)),
                    new: new.as_ref().map(|c| c.bind(py).borrow().to_rust(py)),
                }
            }
        }
    }
}

/// A resolved edit adding or removing a molecule constraint.
#[pyclass(frozen)]
pub enum ConstraintDelta {
    Add { constraint: Py<Constraint> },
    Remove { constraint: Py<Constraint> },
}

#[pymethods]
impl ConstraintDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            Self::Add { .. } => "Add",
            Self::Remove { .. } => "Remove",
        };
        entity_delta_repr(
            slf.bind(py).as_any(),
            "ConstraintDelta",
            variant,
            &["constraint"],
        )
    }

    /// Return the inverse resolved edit.
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl ConstraintDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &GraphIrConstraintDelta) -> PyResult<Self> {
        Ok(match delta {
            GraphIrConstraintDelta::Add(constraint) => Self::Add {
                constraint: into_py_variant(py, Constraint::from_rust(py, constraint)?)?,
            },
            GraphIrConstraintDelta::Remove(constraint) => Self::Remove {
                constraint: into_py_variant(py, Constraint::from_rust(py, constraint)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrConstraintDelta {
        match self {
            Self::Add { constraint } => {
                GraphIrConstraintDelta::Add(constraint.bind(py).borrow().to_rust(py))
            }
            Self::Remove { constraint } => {
                GraphIrConstraintDelta::Remove(constraint.bind(py).borrow().to_rust(py))
            }
        }
    }
}

/// One resolved edit from any localized-topology entity kind.
#[pyclass(frozen)]
pub enum Delta {
    Atom(Py<AtomDelta>),
    Bond(Py<BondDelta>),
    DativeBond(Py<DativeBondDelta>),
    AromaticSystem(Py<AromaticSystemDelta>),
    MulticenterBond(Py<MulticenterBondDelta>),
    NoncovalentBond(Py<NoncovalentBondDelta>),
    StereoAtom(Py<StereoAtomDelta>),
    StereoBond(Py<StereoBondDelta>),
    Constraint(Py<ConstraintDelta>),
}

#[pymethods]
impl Delta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            Self::Atom(_) => "Atom",
            Self::Bond(_) => "Bond",
            Self::DativeBond(_) => "DativeBond",
            Self::AromaticSystem(_) => "AromaticSystem",
            Self::MulticenterBond(_) => "MulticenterBond",
            Self::NoncovalentBond(_) => "NoncovalentBond",
            Self::StereoAtom(_) => "StereoAtom",
            Self::StereoBond(_) => "StereoBond",
            Self::Constraint(_) => "Constraint",
        };
        variant_repr(slf.bind(py).as_any(), "Delta", variant, 1)
    }

    /// Return the inverse resolved edit.
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl Delta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &GraphIrDelta) -> PyResult<Self> {
        Ok(match delta {
            GraphIrDelta::Atom(delta) => {
                Self::Atom(into_py_variant(py, AtomDelta::from_rust(py, delta)?)?)
            }
            GraphIrDelta::Bond(delta) => {
                Self::Bond(into_py_variant(py, BondDelta::from_rust(py, delta)?)?)
            }
            GraphIrDelta::DativeBond(delta) => {
                Self::DativeBond(into_py_variant(py, DativeBondDelta::from_rust(py, delta)?)?)
            }
            GraphIrDelta::AromaticSystem(delta) => Self::AromaticSystem(into_py_variant(
                py,
                AromaticSystemDelta::from_rust(py, delta)?,
            )?),
            GraphIrDelta::MulticenterBond(delta) => Self::MulticenterBond(into_py_variant(
                py,
                MulticenterBondDelta::from_rust(py, delta)?,
            )?),
            GraphIrDelta::NoncovalentBond(delta) => Self::NoncovalentBond(into_py_variant(
                py,
                NoncovalentBondDelta::from_rust(py, delta)?,
            )?),
            GraphIrDelta::StereoAtom(delta) => {
                Self::StereoAtom(into_py_variant(py, StereoAtomDelta::from_rust(py, delta)?)?)
            }
            GraphIrDelta::StereoBond(delta) => {
                Self::StereoBond(into_py_variant(py, StereoBondDelta::from_rust(py, delta)?)?)
            }
            GraphIrDelta::Constraint(delta) => {
                Self::Constraint(into_py_variant(py, ConstraintDelta::from_rust(py, delta)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrDelta {
        match self {
            Self::Atom(delta) => GraphIrDelta::Atom(delta.bind(py).borrow().to_rust(py)),
            Self::Bond(delta) => GraphIrDelta::Bond(delta.bind(py).borrow().to_rust(py)),
            Self::DativeBond(delta) => {
                GraphIrDelta::DativeBond(delta.bind(py).borrow().to_rust(py))
            }
            Self::AromaticSystem(delta) => {
                GraphIrDelta::AromaticSystem(delta.bind(py).borrow().to_rust(py))
            }
            Self::MulticenterBond(delta) => {
                GraphIrDelta::MulticenterBond(delta.bind(py).borrow().to_rust(py))
            }
            Self::NoncovalentBond(delta) => {
                GraphIrDelta::NoncovalentBond(delta.bind(py).borrow().to_rust(py))
            }
            Self::StereoAtom(delta) => {
                GraphIrDelta::StereoAtom(delta.bind(py).borrow().to_rust(py))
            }
            Self::StereoBond(delta) => {
                GraphIrDelta::StereoBond(delta.bind(py).borrow().to_rust(py))
            }
            Self::Constraint(delta) => {
                GraphIrDelta::Constraint(delta.bind(py).borrow().to_rust(py))
            }
        }
    }
}

/// Resolve a possibly-negative Python index into an existing delta position.
fn resolve_delta_index(len: usize, index: isize) -> PyResult<usize> {
    let resolved = if index < 0 {
        index + len as isize
    } else {
        index
    };
    if resolved < 0 || resolved as usize >= len {
        Err(PyIndexError::new_err("delta index out of range"))
    } else {
        Ok(resolved as usize)
    }
}

/// Build a detached iterator of concrete Python delta variants.
fn delta_iter(py: Python<'_>, deltas: &GraphIrDeltas) -> PyResult<DeltaIter> {
    let entries = deltas
        .iter()
        .map(|delta| into_py_variant(py, Delta::from_rust(py, delta)?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(DeltaIter {
        entries: entries.into_iter(),
    })
}

/// A snapshot iterator over resolved deltas.
#[pyclass]
pub(crate) struct DeltaIter {
    entries: IntoIter<Py<Delta>>,
}

#[pymethods]
impl DeltaIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<Delta>> {
        self.entries.next()
    }
}

/// The argument to `Deltas.extend`: another container or delta entries.
#[derive(FromPyObject)]
pub(crate) enum DeltasExtend {
    Container(Py<Deltas>),
    Entries(Vec<Py<Delta>>),
}

impl DeltasExtend {
    /// Snapshot every Python input before the target takes a write borrow.
    fn resolve(&self, py: Python<'_>) -> ResolvedDeltasExtend {
        let entries = match self {
            Self::Container(container) => container.bind(py).borrow().to_rust().as_slice().to_vec(),
            Self::Entries(entries) => entries
                .iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py))
                .collect(),
        };
        ResolvedDeltasExtend(entries)
    }
}

/// An extend input containing no Python references that need to be read.
pub(crate) struct ResolvedDeltasExtend(Vec<GraphIrDelta>);

impl ResolvedDeltasExtend {
    /// Append the resolved deltas in order, preserving duplicates.
    fn apply(self, target: &mut GraphIrDeltas) {
        for delta in self.0 {
            target.push(delta);
        }
    }
}

/// Resolved deltas in insertion order. Mutable, value-equal, and unhashable.
#[pyclass(eq)]
#[derive(Debug, PartialEq)]
pub struct Deltas(GraphIrDeltas);

#[pymethods]
impl Deltas {
    /// Build an owned container from delta entries, preserving order and duplicates.
    #[new]
    #[pyo3(signature = (entries=Vec::new()))]
    fn new(py: Python<'_>, entries: Vec<Py<Delta>>) -> Self {
        Self::from_rust(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py))
                .collect(),
        )
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let value = into_py_variant(py, Delta::from_rust(py, entry)?)?;
            parts.push(value.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("Deltas([{}])", parts.join(", ")))
    }

    /// Append one detached delta snapshot.
    fn append(&mut self, py: Python<'_>, delta: Py<Delta>) {
        self.0.push(delta.bind(py).borrow().to_rust(py));
    }

    /// Append another container or iterable after snapshotting the complete RHS.
    fn extend(slf: Py<Self>, py: Python<'_>, other: DeltasExtend) {
        let resolved = other.resolve(py);
        resolved.apply(slf.borrow_mut(py).to_rust_mut());
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<Delta> {
        let index = resolve_delta_index(self.0.len(), index)?;
        Delta::from_rust(py, &self.0.as_slice()[index])
    }

    fn __iter__(&self, py: Python<'_>) -> PyResult<DeltaIter> {
        delta_iter(py, &self.0)
    }
}

impl_py_normalize!(
    Deltas,
    GraphIrDeltas,
    |value: &Deltas, _py: Python<'_>| -> PyResult<GraphIrDeltas> { Ok(value.to_rust().clone()) },
    |_py: Python<'_>, value: GraphIrDeltas| -> PyResult<Deltas> { Ok(Deltas::from_rust(value)) }
);

impl Deltas {
    pub(crate) fn from_rust(deltas: GraphIrDeltas) -> Self {
        Self(deltas)
    }

    pub(crate) fn to_rust(&self) -> &GraphIrDeltas {
        &self.0
    }

    fn to_rust_mut(&mut self) -> &mut GraphIrDeltas {
        &mut self.0
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element as ChemElement;
    use umol_graph_ir::ir::{
        AromaticSystemConstraintForm as GraphIrAromaticSystemConstraintForm,
        AtomConstraintForm as GraphIrAtomConstraintForm,
        BondConstraintForm as GraphIrBondConstraintForm, BooleanForm as GraphIrBooleanForm,
        DativeBondConstraintForm as GraphIrDativeBondConstraintForm,
        ElectronCountsForm as GraphIrElectronCountsForm, ElementForm as GraphIrElementForm,
        IsotopeMassForm as GraphIrIsotopeMassForm,
        MulticenterBondConstraintForm as GraphIrMulticenterBondConstraintForm,
        NoncovalentBondConstraintForm as GraphIrNoncovalentBondConstraintForm,
        NoncovalentBondKind as GraphIrNoncovalentBondKind,
        NoncovalentBondKindForm as GraphIrNoncovalentBondKindForm, NumForm as GraphIrNumForm,
        StereoAtomConstraintForm as GraphIrStereoAtomConstraintForm,
        StereoBondConstraintForm as GraphIrStereoBondConstraintForm,
        StereoConfigurationForm as GraphIrStereoConfigurationForm,
        StereoCoset as GraphIrStereoCoset, StereoKind as GraphIrStereoKind,
        StereoLigand as GraphIrStereoLigand, StereoLigandKind as GraphIrStereoLigandKind,
        Stereogenicity as GraphIrStereogenicity, StereogenicityForm as GraphIrStereogenicityForm,
        UnpairedElectronsForm as GraphIrUnpairedElectronsForm,
    };

    use super::*;
    use crate::error::ContradictionError;

    #[rstest]
    #[case::element(GraphIrAtomFieldChange::Element {
        old: GraphIrElementForm::Lit(ChemElement::C),
        new: GraphIrElementForm::Lit(ChemElement::N),
    })]
    #[case::isotope_mass(GraphIrAtomFieldChange::IsotopeMass {
        old: GraphIrIsotopeMassForm::Lit(12),
        new: GraphIrIsotopeMassForm::Lit(13),
    })]
    #[case::charge(GraphIrAtomFieldChange::Charge {
        old: GraphIrNumForm::Lit(0),
        new: GraphIrNumForm::Lit(-1),
    })]
    #[case::implicit_hydrogens(GraphIrAtomFieldChange::ImplicitHydrogens {
        old: GraphIrNumForm::Lit(3),
        new: GraphIrNumForm::Lit(2),
    })]
    #[case::lone_pairs(GraphIrAtomFieldChange::LonePairs {
        old: GraphIrNumForm::Lit(1),
        new: GraphIrNumForm::Lit(2),
    })]
    #[case::unpaired_electrons(GraphIrAtomFieldChange::UnpairedElectrons {
        old: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(0),
            multiplicity: GraphIrNumForm::Lit(1),
        },
        new: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(1),
            multiplicity: GraphIrNumForm::Lit(2),
        },
    })]
    fn test_atom_field_change_roundtrip(#[case] change: GraphIrAtomFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                AtomFieldChange::from_rust(py, &change).unwrap().to_rust(py),
                change
            );
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrAtomFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(-1),
        },
        GraphIrAtomFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(-1),
        },
        true,
    )]
    #[case::different(
        GraphIrAtomFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(-1),
        },
        GraphIrAtomFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(1),
        },
        false,
    )]
    fn test_atom_field_change_eq(
        #[case] lhs: GraphIrAtomFieldChange,
        #[case] rhs: GraphIrAtomFieldChange,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = AtomFieldChange::from_rust(py, &lhs).unwrap();
            let rhs = AtomFieldChange::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::element(
        GraphIrAtomFieldChange::Element {
            old: GraphIrElementForm::Lit(ChemElement::C),
            new: GraphIrElementForm::Lit(ChemElement::N),
        },
        "ElementForm.Lit(Element('C'))",
        "ElementForm.Lit(Element('N'))",
        "AtomFieldChange.Element(old=ElementForm.Lit(Element('C')), new=ElementForm.Lit(Element('N')))"
    )]
    #[case::isotope_mass(
        GraphIrAtomFieldChange::IsotopeMass {
            old: GraphIrIsotopeMassForm::Lit(12),
            new: GraphIrIsotopeMassForm::Lit(13),
        },
        "IsotopeMassForm.Lit(12)",
        "IsotopeMassForm.Lit(13)",
        "AtomFieldChange.IsotopeMass(old=IsotopeMassForm.Lit(12), new=IsotopeMassForm.Lit(13))"
    )]
    #[case::charge(
        GraphIrAtomFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(-1),
        },
        "NumForm.Lit(0)",
        "NumForm.Lit(-1)",
        "AtomFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(-1))"
    )]
    #[case::implicit_hydrogens(
        GraphIrAtomFieldChange::ImplicitHydrogens {
            old: GraphIrNumForm::Lit(3),
            new: GraphIrNumForm::Lit(2),
        },
        "NumForm.Lit(3)",
        "NumForm.Lit(2)",
        "AtomFieldChange.ImplicitHydrogens(old=NumForm.Lit(3), new=NumForm.Lit(2))"
    )]
    #[case::lone_pairs(
        GraphIrAtomFieldChange::LonePairs {
            old: GraphIrNumForm::Lit(1),
            new: GraphIrNumForm::Lit(2),
        },
        "NumForm.Lit(1)",
        "NumForm.Lit(2)",
        "AtomFieldChange.LonePairs(old=NumForm.Lit(1), new=NumForm.Lit(2))"
    )]
    #[case::unpaired_electrons(
        GraphIrAtomFieldChange::UnpairedElectrons {
            old: GraphIrUnpairedElectronsForm {
                count: GraphIrNumForm::Lit(0),
                multiplicity: GraphIrNumForm::Lit(1),
            },
            new: GraphIrUnpairedElectronsForm {
                count: GraphIrNumForm::Lit(1),
                multiplicity: GraphIrNumForm::Lit(2),
            },
        },
        "UnpairedElectronsForm(NumForm.Lit(0), NumForm.Lit(1))",
        "UnpairedElectronsForm(NumForm.Lit(1), NumForm.Lit(2))",
        "AtomFieldChange.UnpairedElectrons(old=UnpairedElectronsForm(NumForm.Lit(0), NumForm.Lit(1)), new=UnpairedElectronsForm(NumForm.Lit(1), NumForm.Lit(2)))"
    )]
    fn test_atom_field_change_repr(
        #[case] change: GraphIrAtomFieldChange,
        #[case] old: &str,
        #[case] new: &str,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let change =
                into_py_variant(py, AtomFieldChange::from_rust(py, &change).unwrap()).unwrap();
            let bound = change.bind(py).as_any();
            assert_eq!(
                bound
                    .getattr("old")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                old
            );
            assert_eq!(
                bound
                    .getattr("new")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                new
            );
            assert_eq!(bound.repr().unwrap().extract::<String>().unwrap(), expected);
        });
    }

    #[rstest]
    #[case::element(GraphIrAtomFieldChange::Element {
        old: GraphIrElementForm::Lit(ChemElement::C),
        new: GraphIrElementForm::Lit(ChemElement::N),
    })]
    #[case::isotope_mass(GraphIrAtomFieldChange::IsotopeMass {
        old: GraphIrIsotopeMassForm::Lit(12),
        new: GraphIrIsotopeMassForm::Lit(13),
    })]
    #[case::charge(GraphIrAtomFieldChange::Charge {
        old: GraphIrNumForm::Lit(0),
        new: GraphIrNumForm::Lit(-1),
    })]
    #[case::implicit_hydrogens(GraphIrAtomFieldChange::ImplicitHydrogens {
        old: GraphIrNumForm::Lit(3),
        new: GraphIrNumForm::Lit(2),
    })]
    #[case::lone_pairs(GraphIrAtomFieldChange::LonePairs {
        old: GraphIrNumForm::Lit(1),
        new: GraphIrNumForm::Lit(2),
    })]
    #[case::unpaired_electrons(GraphIrAtomFieldChange::UnpairedElectrons {
        old: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(0),
            multiplicity: GraphIrNumForm::Lit(1),
        },
        new: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(1),
            multiplicity: GraphIrNumForm::Lit(2),
        },
    })]
    fn test_atom_field_change_inverse(#[case] change: GraphIrAtomFieldChange) {
        Python::attach(|py| {
            let binding = AtomFieldChange::from_rust(py, &change).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), change);
        });
    }

    #[rstest]
    #[case::order(GraphIrBondFieldChange::Order {
        old: GraphIrNumForm::Lit(1),
        new: GraphIrNumForm::Lit(2),
    })]
    #[case::charge(GraphIrBondFieldChange::Charge {
        old: GraphIrNumForm::Lit(0),
        new: GraphIrNumForm::Lit(1),
    })]
    #[case::unpaired_electrons(GraphIrBondFieldChange::UnpairedElectrons {
        old: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(0),
            multiplicity: GraphIrNumForm::Lit(1),
        },
        new: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(1),
            multiplicity: GraphIrNumForm::Lit(2),
        },
    })]
    fn test_bond_field_change_roundtrip(#[case] change: GraphIrBondFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                BondFieldChange::from_rust(py, &change).unwrap().to_rust(py),
                change
            );
        });
    }

    #[rstest]
    #[case::order(GraphIrBondFieldChange::Order {
        old: GraphIrNumForm::Lit(1),
        new: GraphIrNumForm::Lit(2),
    })]
    #[case::charge(GraphIrBondFieldChange::Charge {
        old: GraphIrNumForm::Lit(0),
        new: GraphIrNumForm::Lit(1),
    })]
    #[case::unpaired_electrons(GraphIrBondFieldChange::UnpairedElectrons {
        old: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(0),
            multiplicity: GraphIrNumForm::Lit(1),
        },
        new: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(1),
            multiplicity: GraphIrNumForm::Lit(2),
        },
    })]
    fn test_bond_field_change_inverse(#[case] change: GraphIrBondFieldChange) {
        Python::attach(|py| {
            let binding = BondFieldChange::from_rust(py, &change).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), change);
        });
    }

    #[rstest]
    #[case::order(GraphIrDativeBondFieldChange::Order {
        old: GraphIrNumForm::Lit(1),
        new: GraphIrNumForm::Lit(2),
    })]
    fn test_dative_bond_field_change_roundtrip(#[case] change: GraphIrDativeBondFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                DativeBondFieldChange::from_rust(py, &change)
                    .unwrap()
                    .to_rust(py),
                change
            );
        });
    }

    #[rstest]
    #[case::order(GraphIrDativeBondFieldChange::Order {
        old: GraphIrNumForm::Lit(1),
        new: GraphIrNumForm::Lit(2),
    })]
    fn test_dative_bond_field_change_inverse(#[case] change: GraphIrDativeBondFieldChange) {
        Python::attach(|py| {
            let binding = DativeBondFieldChange::from_rust(py, &change).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), change);
        });
    }

    #[rstest]
    #[case::electrons(GraphIrAromaticSystemFieldChange::Electrons {
        old: GraphIrElectronCountsForm::Undetermined,
        new: GraphIrElectronCountsForm::Lit(vec![1, 1, 1]),
    })]
    #[case::charge(GraphIrAromaticSystemFieldChange::Charge {
        old: GraphIrNumForm::Lit(0),
        new: GraphIrNumForm::Lit(-1),
    })]
    #[case::unpaired_electrons(GraphIrAromaticSystemFieldChange::UnpairedElectrons {
        old: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(0),
            multiplicity: GraphIrNumForm::Lit(1),
        },
        new: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(1),
            multiplicity: GraphIrNumForm::Lit(2),
        },
    })]
    fn test_aromatic_system_field_change_roundtrip(
        #[case] change: GraphIrAromaticSystemFieldChange,
    ) {
        Python::attach(|py| {
            assert_eq!(
                AromaticSystemFieldChange::from_rust(py, &change)
                    .unwrap()
                    .to_rust(py),
                change
            );
        });
    }

    #[rstest]
    #[case::electrons(GraphIrAromaticSystemFieldChange::Electrons {
        old: GraphIrElectronCountsForm::Undetermined,
        new: GraphIrElectronCountsForm::Lit(vec![1, 1, 1]),
    })]
    #[case::charge(GraphIrAromaticSystemFieldChange::Charge {
        old: GraphIrNumForm::Lit(0),
        new: GraphIrNumForm::Lit(-1),
    })]
    #[case::unpaired_electrons(GraphIrAromaticSystemFieldChange::UnpairedElectrons {
        old: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(0),
            multiplicity: GraphIrNumForm::Lit(1),
        },
        new: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(1),
            multiplicity: GraphIrNumForm::Lit(2),
        },
    })]
    fn test_aromatic_system_field_change_inverse(#[case] change: GraphIrAromaticSystemFieldChange) {
        Python::attach(|py| {
            let binding = AromaticSystemFieldChange::from_rust(py, &change).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), change);
        });
    }

    #[rstest]
    #[case::electrons(GraphIrMulticenterBondFieldChange::Electrons {
        old: GraphIrElectronCountsForm::Lit(vec![1, 0, 1]),
        new: GraphIrElectronCountsForm::Lit(vec![2, 0, 1]),
    })]
    #[case::charge(GraphIrMulticenterBondFieldChange::Charge {
        old: GraphIrNumForm::Lit(0),
        new: GraphIrNumForm::Lit(1),
    })]
    #[case::unpaired_electrons(GraphIrMulticenterBondFieldChange::UnpairedElectrons {
        old: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(0),
            multiplicity: GraphIrNumForm::Lit(1),
        },
        new: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(2),
            multiplicity: GraphIrNumForm::Lit(3),
        },
    })]
    fn test_multicenter_bond_field_change_roundtrip(
        #[case] change: GraphIrMulticenterBondFieldChange,
    ) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterBondFieldChange::from_rust(py, &change)
                    .unwrap()
                    .to_rust(py),
                change
            );
        });
    }

    #[rstest]
    #[case::electrons(GraphIrMulticenterBondFieldChange::Electrons {
        old: GraphIrElectronCountsForm::Lit(vec![1, 0, 1]),
        new: GraphIrElectronCountsForm::Lit(vec![2, 0, 1]),
    })]
    #[case::charge(GraphIrMulticenterBondFieldChange::Charge {
        old: GraphIrNumForm::Lit(0),
        new: GraphIrNumForm::Lit(1),
    })]
    #[case::unpaired_electrons(GraphIrMulticenterBondFieldChange::UnpairedElectrons {
        old: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(0),
            multiplicity: GraphIrNumForm::Lit(1),
        },
        new: GraphIrUnpairedElectronsForm {
            count: GraphIrNumForm::Lit(2),
            multiplicity: GraphIrNumForm::Lit(3),
        },
    })]
    fn test_multicenter_bond_field_change_inverse(
        #[case] change: GraphIrMulticenterBondFieldChange,
    ) {
        Python::attach(|py| {
            let binding = MulticenterBondFieldChange::from_rust(py, &change).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), change);
        });
    }

    #[rstest]
    #[case::kind(GraphIrNoncovalentBondFieldChange::Kind {
        old: GraphIrNoncovalentBondKindForm::Undetermined,
        new: GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond),
    })]
    fn test_noncovalent_bond_field_change_roundtrip(
        #[case] change: GraphIrNoncovalentBondFieldChange,
    ) {
        Python::attach(|py| {
            assert_eq!(
                NoncovalentBondFieldChange::from_rust(py, &change)
                    .unwrap()
                    .to_rust(py),
                change
            );
        });
    }

    #[rstest]
    #[case::kind(GraphIrNoncovalentBondFieldChange::Kind {
        old: GraphIrNoncovalentBondKindForm::Undetermined,
        new: GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond),
    })]
    fn test_noncovalent_bond_field_change_inverse(
        #[case] change: GraphIrNoncovalentBondFieldChange,
    ) {
        Python::attach(|py| {
            let binding = NoncovalentBondFieldChange::from_rust(py, &change).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), change);
        });
    }

    #[rstest]
    #[case::geometry_unknown(GraphIrStereoAtomFieldChange::Configuration {
        old: GraphIrStereoConfigurationForm::Undetermined,
        new: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::Tetrahedral,
            GraphIrStereoCoset::Undetermined,
        ),
    })]
    #[case::coset_resolved(GraphIrStereoAtomFieldChange::Configuration {
        old: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::Tetrahedral,
            GraphIrStereoCoset::Undetermined,
        ),
        new: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::Tetrahedral,
            GraphIrStereoCoset::Lit(1),
        ),
    })]
    fn test_stereo_atom_field_change_roundtrip(#[case] change: GraphIrStereoAtomFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                StereoAtomFieldChange::from_rust(py, &change)
                    .unwrap()
                    .to_rust(py),
                change
            );
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrStereoAtomFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Undetermined,
            ),
        },
        GraphIrStereoAtomFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Undetermined,
            ),
        },
        true,
    )]
    #[case::different(
        GraphIrStereoAtomFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Undetermined,
            ),
        },
        GraphIrStereoAtomFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(1),
            ),
        },
        false,
    )]
    fn test_stereo_atom_field_change_eq(
        #[case] lhs: GraphIrStereoAtomFieldChange,
        #[case] rhs: GraphIrStereoAtomFieldChange,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = StereoAtomFieldChange::from_rust(py, &lhs).unwrap();
            let rhs = StereoAtomFieldChange::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::geometry_unknown(
        GraphIrStereoAtomFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Undetermined,
            ),
        },
        "StereoConfigurationForm.Undetermined()",
        "StereoConfigurationForm.Kinded(StereoKind.Tetrahedral, StereoCoset.Undetermined())",
        "StereoAtomFieldChange.Configuration(old=StereoConfigurationForm.Undetermined(), new=StereoConfigurationForm.Kinded(StereoKind.Tetrahedral, StereoCoset.Undetermined()))",
    )]
    #[case::coset_resolved(
        GraphIrStereoAtomFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Undetermined,
            ),
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(1),
            ),
        },
        "StereoConfigurationForm.Kinded(StereoKind.Tetrahedral, StereoCoset.Undetermined())",
        "StereoConfigurationForm.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(1))",
        "StereoAtomFieldChange.Configuration(old=StereoConfigurationForm.Kinded(StereoKind.Tetrahedral, StereoCoset.Undetermined()), new=StereoConfigurationForm.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(1)))",
    )]
    fn test_stereo_atom_field_change_repr(
        #[case] change: GraphIrStereoAtomFieldChange,
        #[case] old: &str,
        #[case] new: &str,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let change =
                into_py_variant(py, StereoAtomFieldChange::from_rust(py, &change).unwrap())
                    .unwrap();
            let bound = change.bind(py).as_any();
            assert_eq!(
                bound
                    .getattr("old")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                old
            );
            assert_eq!(
                bound
                    .getattr("new")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                new
            );
            assert_eq!(bound.repr().unwrap().extract::<String>().unwrap(), expected);
        });
    }

    #[rstest]
    #[case::geometry_unknown(GraphIrStereoAtomFieldChange::Configuration {
        old: GraphIrStereoConfigurationForm::Undetermined,
        new: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::Tetrahedral,
            GraphIrStereoCoset::Undetermined,
        ),
    })]
    #[case::coset_resolved(GraphIrStereoAtomFieldChange::Configuration {
        old: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::Tetrahedral,
            GraphIrStereoCoset::Undetermined,
        ),
        new: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::Tetrahedral,
            GraphIrStereoCoset::Lit(1),
        ),
    })]
    fn test_stereo_atom_field_change_inverse(#[case] change: GraphIrStereoAtomFieldChange) {
        Python::attach(|py| {
            let binding = StereoAtomFieldChange::from_rust(py, &change).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), change);
        });
    }

    #[rstest]
    #[case::geometry_unknown(GraphIrStereoBondFieldChange::Configuration {
        old: GraphIrStereoConfigurationForm::Undetermined,
        new: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::CisTrans,
            GraphIrStereoCoset::Undetermined,
        ),
    })]
    #[case::coset_resolved(GraphIrStereoBondFieldChange::Configuration {
        old: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::CisTrans,
            GraphIrStereoCoset::Undetermined,
        ),
        new: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::CisTrans,
            GraphIrStereoCoset::Lit(1),
        ),
    })]
    fn test_stereo_bond_field_change_roundtrip(#[case] change: GraphIrStereoBondFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                StereoBondFieldChange::from_rust(py, &change)
                    .unwrap()
                    .to_rust(py),
                change
            );
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrStereoBondFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Undetermined,
            ),
        },
        GraphIrStereoBondFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Undetermined,
            ),
        },
        true,
    )]
    #[case::different(
        GraphIrStereoBondFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Undetermined,
            ),
        },
        GraphIrStereoBondFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Lit(1),
            ),
        },
        false,
    )]
    fn test_stereo_bond_field_change_eq(
        #[case] lhs: GraphIrStereoBondFieldChange,
        #[case] rhs: GraphIrStereoBondFieldChange,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = StereoBondFieldChange::from_rust(py, &lhs).unwrap();
            let rhs = StereoBondFieldChange::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::geometry_unknown(
        GraphIrStereoBondFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Undetermined,
            ),
        },
        "StereoConfigurationForm.Undetermined()",
        "StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Undetermined())",
        "StereoBondFieldChange.Configuration(old=StereoConfigurationForm.Undetermined(), new=StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Undetermined()))",
    )]
    #[case::coset_resolved(
        GraphIrStereoBondFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Undetermined,
            ),
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Lit(1),
            ),
        },
        "StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Undetermined())",
        "StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Lit(1))",
        "StereoBondFieldChange.Configuration(old=StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Undetermined()), new=StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Lit(1)))",
    )]
    fn test_stereo_bond_field_change_repr(
        #[case] change: GraphIrStereoBondFieldChange,
        #[case] old: &str,
        #[case] new: &str,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let change =
                into_py_variant(py, StereoBondFieldChange::from_rust(py, &change).unwrap())
                    .unwrap();
            let bound = change.bind(py).as_any();
            assert_eq!(
                bound
                    .getattr("old")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                old
            );
            assert_eq!(
                bound
                    .getattr("new")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                new
            );
            assert_eq!(bound.repr().unwrap().extract::<String>().unwrap(), expected);
        });
    }

    #[rstest]
    #[case::geometry_unknown(GraphIrStereoBondFieldChange::Configuration {
        old: GraphIrStereoConfigurationForm::Undetermined,
        new: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::CisTrans,
            GraphIrStereoCoset::Undetermined,
        ),
    })]
    #[case::coset_resolved(GraphIrStereoBondFieldChange::Configuration {
        old: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::CisTrans,
            GraphIrStereoCoset::Undetermined,
        ),
        new: GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::CisTrans,
            GraphIrStereoCoset::Lit(1),
        ),
    })]
    fn test_stereo_bond_field_change_inverse(#[case] change: GraphIrStereoBondFieldChange) {
        Python::attach(|py| {
            let binding = StereoBondFieldChange::from_rust(py, &change).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), change);
        });
    }

    #[rstest]
    #[case::add(GraphIrAtomDelta::Add {
        id: GraphIrAtomId(3),
        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
    })]
    #[case::remove(GraphIrAtomDelta::Remove {
        id: GraphIrAtomId(3),
        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
    })]
    #[case::modify_field(GraphIrAtomDelta::ModifyField {
        id: GraphIrAtomId(3),
        change: GraphIrAtomFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(-1),
        },
    })]
    #[case::constraint_added(GraphIrAtomDelta::ModifyConstraint {
        id: GraphIrAtomId(3),
        old: None,
        new: Some(GraphIrAtomConstraintForm::Valence(GraphIrNumForm::Lit(4))),
    })]
    #[case::constraint_removed(GraphIrAtomDelta::ModifyConstraint {
        id: GraphIrAtomId(3),
        old: Some(GraphIrAtomConstraintForm::Valence(GraphIrNumForm::Lit(4))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrAtomDelta::ModifyConstraint {
        id: GraphIrAtomId(3),
        old: Some(GraphIrAtomConstraintForm::Valence(GraphIrNumForm::Lit(3))),
        new: Some(GraphIrAtomConstraintForm::Valence(GraphIrNumForm::Lit(4))),
    })]
    fn test_atom_delta_roundtrip(#[case] delta: GraphIrAtomDelta) {
        Python::attach(|py| {
            assert_eq!(AtomDelta::from_rust(py, &delta).unwrap().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        },
        GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        },
        true,
    )]
    #[case::different(
        GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        },
        GraphIrAtomDelta::Add {
            id: GraphIrAtomId(4),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        },
        false,
    )]
    fn test_atom_delta_eq(
        #[case] lhs: GraphIrAtomDelta,
        #[case] rhs: GraphIrAtomDelta,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = AtomDelta::from_rust(py, &lhs).unwrap();
            let rhs = AtomDelta::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::add(
        GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        },
        "AtomDelta.Add(id=3, attributes=AtomForm.parse('C'))",
    )]
    #[case::remove(
        GraphIrAtomDelta::Remove {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        },
        "AtomDelta.Remove(id=3, attributes=AtomForm.parse('C'))",
    )]
    #[case::modify_field(
        GraphIrAtomDelta::ModifyField {
            id: GraphIrAtomId(3),
            change: GraphIrAtomFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(-1),
            },
        },
        "AtomDelta.ModifyField(id=3, change=AtomFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(-1)))",
    )]
    #[case::modify_constraint(
        GraphIrAtomDelta::ModifyConstraint {
            id: GraphIrAtomId(3),
            old: None,
            new: Some(GraphIrAtomConstraintForm::Valence(GraphIrNumForm::Lit(4))),
        },
        "AtomDelta.ModifyConstraint(id=3, old=None, new=AtomConstraintForm.Valence(NumForm.Lit(4)))",
    )]
    fn test_atom_delta_repr(#[case] delta: GraphIrAtomDelta, #[case] expected: &str) {
        Python::attach(|py| {
            let delta = into_py_variant(py, AtomDelta::from_rust(py, &delta).unwrap()).unwrap();
            assert_eq!(
                delta
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::add(GraphIrAtomDelta::Add {
        id: GraphIrAtomId(3),
        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
    })]
    #[case::remove(GraphIrAtomDelta::Remove {
        id: GraphIrAtomId(3),
        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
    })]
    #[case::modify_field(GraphIrAtomDelta::ModifyField {
        id: GraphIrAtomId(3),
        change: GraphIrAtomFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(-1),
        },
    })]
    #[case::constraint_added(GraphIrAtomDelta::ModifyConstraint {
        id: GraphIrAtomId(3),
        old: None,
        new: Some(GraphIrAtomConstraintForm::Valence(GraphIrNumForm::Lit(4))),
    })]
    #[case::constraint_removed(GraphIrAtomDelta::ModifyConstraint {
        id: GraphIrAtomId(3),
        old: Some(GraphIrAtomConstraintForm::Valence(GraphIrNumForm::Lit(4))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrAtomDelta::ModifyConstraint {
        id: GraphIrAtomId(3),
        old: Some(GraphIrAtomConstraintForm::Valence(GraphIrNumForm::Lit(3))),
        new: Some(GraphIrAtomConstraintForm::Valence(GraphIrNumForm::Lit(4))),
    })]
    fn test_atom_delta_inverse(#[case] delta: GraphIrAtomDelta) {
        Python::attach(|py| {
            let binding = AtomDelta::from_rust(py, &delta).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                delta.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::add(GraphIrBondDelta::Add {
        id: GraphIrBondId(2),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
        attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
    })]
    #[case::remove(GraphIrBondDelta::Remove {
        id: GraphIrBondId(2),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
        attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
    })]
    #[case::modify_field(GraphIrBondDelta::ModifyField {
        id: GraphIrBondId(2),
        change: GraphIrBondFieldChange::Order {
            old: GraphIrNumForm::Lit(1),
            new: GraphIrNumForm::Lit(2),
        },
    })]
    #[case::constraint_added(GraphIrBondDelta::ModifyConstraint {
        id: GraphIrBondId(2),
        old: None,
        new: Some(GraphIrBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
    })]
    #[case::constraint_removed(GraphIrBondDelta::ModifyConstraint {
        id: GraphIrBondId(2),
        old: Some(GraphIrBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrBondDelta::ModifyConstraint {
        id: GraphIrBondId(2),
        old: Some(GraphIrBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(false))),
        new: Some(GraphIrBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
    })]
    fn test_bond_delta_roundtrip(#[case] delta: GraphIrBondDelta) {
        Python::attach(|py| {
            assert_eq!(BondDelta::from_rust(py, &delta).unwrap().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrBondDelta::Add {
            id: GraphIrBondId(2),
            atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
            attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
        },
        GraphIrBondDelta::Add {
            id: GraphIrBondId(2),
            atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
            attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
        },
        true,
    )]
    #[case::different_order(
        GraphIrBondDelta::Add {
            id: GraphIrBondId(2),
            atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
            attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
        },
        GraphIrBondDelta::Add {
            id: GraphIrBondId(2),
            atoms: [GraphIrAtomId(1), GraphIrAtomId(5)],
            attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
        },
        false,
    )]
    fn test_bond_delta_eq(
        #[case] lhs: GraphIrBondDelta,
        #[case] rhs: GraphIrBondDelta,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = BondDelta::from_rust(py, &lhs).unwrap();
            let rhs = BondDelta::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::add(
        GraphIrBondDelta::Add {
            id: GraphIrBondId(2),
            atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
            attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
        },
        "BondDelta.Add(id=2, atoms=(5, 1), attributes=BondForm.parse('1'))",
    )]
    #[case::remove(
        GraphIrBondDelta::Remove {
            id: GraphIrBondId(2),
            atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
            attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
        },
        "BondDelta.Remove(id=2, atoms=(5, 1), attributes=BondForm.parse('1'))",
    )]
    #[case::modify_field(
        GraphIrBondDelta::ModifyField {
            id: GraphIrBondId(2),
            change: GraphIrBondFieldChange::Order {
                old: GraphIrNumForm::Lit(1),
                new: GraphIrNumForm::Lit(2),
            },
        },
        "BondDelta.ModifyField(id=2, change=BondFieldChange.Order(old=NumForm.Lit(1), new=NumForm.Lit(2)))",
    )]
    #[case::modify_constraint(
        GraphIrBondDelta::ModifyConstraint {
            id: GraphIrBondId(2),
            old: None,
            new: Some(GraphIrBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
        },
        "BondDelta.ModifyConstraint(id=2, old=None, new=BondConstraintForm.Aromatic(BooleanForm.Lit(True)))",
    )]
    fn test_bond_delta_repr(#[case] delta: GraphIrBondDelta, #[case] expected: &str) {
        Python::attach(|py| {
            let delta = into_py_variant(py, BondDelta::from_rust(py, &delta).unwrap()).unwrap();
            assert_eq!(
                delta
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::add(GraphIrBondDelta::Add {
        id: GraphIrBondId(2),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
        attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
    })]
    #[case::remove(GraphIrBondDelta::Remove {
        id: GraphIrBondId(2),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
        attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
    })]
    #[case::modify_field(GraphIrBondDelta::ModifyField {
        id: GraphIrBondId(2),
        change: GraphIrBondFieldChange::Order {
            old: GraphIrNumForm::Lit(1),
            new: GraphIrNumForm::Lit(2),
        },
    })]
    #[case::constraint_added(GraphIrBondDelta::ModifyConstraint {
        id: GraphIrBondId(2),
        old: None,
        new: Some(GraphIrBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
    })]
    #[case::constraint_removed(GraphIrBondDelta::ModifyConstraint {
        id: GraphIrBondId(2),
        old: Some(GraphIrBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrBondDelta::ModifyConstraint {
        id: GraphIrBondId(2),
        old: Some(GraphIrBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(false))),
        new: Some(GraphIrBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
    })]
    fn test_bond_delta_inverse(#[case] delta: GraphIrBondDelta) {
        Python::attach(|py| {
            let binding = BondDelta::from_rust(py, &delta).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                delta.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::add(GraphIrDativeBondDelta::Add {
        id: GraphIrDativeBondId(1),
        donors: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        acceptor: GraphIrAtomId(3),
        attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
    })]
    #[case::remove(GraphIrDativeBondDelta::Remove {
        id: GraphIrDativeBondId(1),
        donors: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        acceptor: GraphIrAtomId(3),
        attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
    })]
    #[case::modify_field(GraphIrDativeBondDelta::ModifyField {
        id: GraphIrDativeBondId(1),
        change: GraphIrDativeBondFieldChange::Order {
            old: GraphIrNumForm::Lit(1),
            new: GraphIrNumForm::Lit(2),
        },
    })]
    #[case::constraint_added(GraphIrDativeBondDelta::ModifyConstraint {
        id: GraphIrDativeBondId(1),
        old: None,
        new: Some(GraphIrDativeBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
    })]
    #[case::constraint_removed(GraphIrDativeBondDelta::ModifyConstraint {
        id: GraphIrDativeBondId(1),
        old: Some(GraphIrDativeBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrDativeBondDelta::ModifyConstraint {
        id: GraphIrDativeBondId(1),
        old: Some(GraphIrDativeBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(false))),
        new: Some(GraphIrDativeBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
    })]
    fn test_dative_bond_delta_roundtrip(#[case] delta: GraphIrDativeBondDelta) {
        Python::attach(|py| {
            assert_eq!(
                DativeBondDelta::from_rust(py, &delta).unwrap().to_rust(py),
                delta
            );
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrDativeBondDelta::Add {
            id: GraphIrDativeBondId(1),
            donors: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            acceptor: GraphIrAtomId(3),
            attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
        },
        GraphIrDativeBondDelta::Add {
            id: GraphIrDativeBondId(1),
            donors: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            acceptor: GraphIrAtomId(3),
            attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
        },
        true,
    )]
    #[case::different_donor_order(
        GraphIrDativeBondDelta::Add {
            id: GraphIrDativeBondId(1),
            donors: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            acceptor: GraphIrAtomId(3),
            attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
        },
        GraphIrDativeBondDelta::Add {
            id: GraphIrDativeBondId(1),
            donors: vec![GraphIrAtomId(2), GraphIrAtomId(4), GraphIrAtomId(4)],
            acceptor: GraphIrAtomId(3),
            attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
        },
        false,
    )]
    fn test_dative_bond_delta_eq(
        #[case] lhs: GraphIrDativeBondDelta,
        #[case] rhs: GraphIrDativeBondDelta,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = DativeBondDelta::from_rust(py, &lhs).unwrap();
            let rhs = DativeBondDelta::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::add(
        GraphIrDativeBondDelta::Add {
            id: GraphIrDativeBondId(1),
            donors: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            acceptor: GraphIrAtomId(3),
            attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
        },
        "DativeBondDelta.Add(id=1, donors=[4, 2, 4], acceptor=3, attributes=DativeBondForm.parse('1'))",
    )]
    #[case::remove(
        GraphIrDativeBondDelta::Remove {
            id: GraphIrDativeBondId(1),
            donors: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            acceptor: GraphIrAtomId(3),
            attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
        },
        "DativeBondDelta.Remove(id=1, donors=[4, 2, 4], acceptor=3, attributes=DativeBondForm.parse('1'))",
    )]
    #[case::modify_field(
        GraphIrDativeBondDelta::ModifyField {
            id: GraphIrDativeBondId(1),
            change: GraphIrDativeBondFieldChange::Order {
                old: GraphIrNumForm::Lit(1),
                new: GraphIrNumForm::Lit(2),
            },
        },
        "DativeBondDelta.ModifyField(id=1, change=DativeBondFieldChange.Order(old=NumForm.Lit(1), new=NumForm.Lit(2)))",
    )]
    #[case::modify_constraint(
        GraphIrDativeBondDelta::ModifyConstraint {
            id: GraphIrDativeBondId(1),
            old: None,
            new: Some(GraphIrDativeBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
        },
        "DativeBondDelta.ModifyConstraint(id=1, old=None, new=DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)))",
    )]
    fn test_dative_bond_delta_repr(#[case] delta: GraphIrDativeBondDelta, #[case] expected: &str) {
        Python::attach(|py| {
            let delta =
                into_py_variant(py, DativeBondDelta::from_rust(py, &delta).unwrap()).unwrap();
            assert_eq!(
                delta
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::add(GraphIrDativeBondDelta::Add {
        id: GraphIrDativeBondId(1),
        donors: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        acceptor: GraphIrAtomId(3),
        attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
    })]
    #[case::remove(GraphIrDativeBondDelta::Remove {
        id: GraphIrDativeBondId(1),
        donors: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        acceptor: GraphIrAtomId(3),
        attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
    })]
    #[case::modify_field(GraphIrDativeBondDelta::ModifyField {
        id: GraphIrDativeBondId(1),
        change: GraphIrDativeBondFieldChange::Order {
            old: GraphIrNumForm::Lit(1),
            new: GraphIrNumForm::Lit(2),
        },
    })]
    #[case::constraint_added(GraphIrDativeBondDelta::ModifyConstraint {
        id: GraphIrDativeBondId(1),
        old: None,
        new: Some(GraphIrDativeBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
    })]
    #[case::constraint_removed(GraphIrDativeBondDelta::ModifyConstraint {
        id: GraphIrDativeBondId(1),
        old: Some(GraphIrDativeBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrDativeBondDelta::ModifyConstraint {
        id: GraphIrDativeBondId(1),
        old: Some(GraphIrDativeBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(false))),
        new: Some(GraphIrDativeBondConstraintForm::Aromatic(GraphIrBooleanForm::Lit(true))),
    })]
    fn test_dative_bond_delta_inverse(#[case] delta: GraphIrDativeBondDelta) {
        Python::attach(|py| {
            let binding = DativeBondDelta::from_rust(py, &delta).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                delta.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::add(GraphIrAromaticSystemDelta::Add {
        id: GraphIrAromaticSystemId(2),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1, 1]),
    })]
    #[case::remove(GraphIrAromaticSystemDelta::Remove {
        id: GraphIrAromaticSystemId(2),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1, 1]),
    })]
    #[case::modify_field(GraphIrAromaticSystemDelta::ModifyField {
        id: GraphIrAromaticSystemId(2),
        change: GraphIrAromaticSystemFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(-1),
        },
    })]
    #[case::constraint_added(GraphIrAromaticSystemDelta::ModifyConstraint {
        id: GraphIrAromaticSystemId(2),
        old: None,
        new: Some(GraphIrAromaticSystemConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
    })]
    #[case::constraint_removed(GraphIrAromaticSystemDelta::ModifyConstraint {
        id: GraphIrAromaticSystemId(2),
        old: Some(GraphIrAromaticSystemConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrAromaticSystemDelta::ModifyConstraint {
        id: GraphIrAromaticSystemId(2),
        old: Some(GraphIrAromaticSystemConstraintForm::ElectronCount(GraphIrNumForm::Lit(5))),
        new: Some(GraphIrAromaticSystemConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
    })]
    fn test_aromatic_system_delta_roundtrip(#[case] delta: GraphIrAromaticSystemDelta) {
        Python::attach(|py| {
            assert_eq!(
                AromaticSystemDelta::from_rust(py, &delta)
                    .unwrap()
                    .to_rust(py),
                delta
            );
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrAromaticSystemDelta::Add {
            id: GraphIrAromaticSystemId(2),
            atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1, 1]),
        },
        GraphIrAromaticSystemDelta::Add {
            id: GraphIrAromaticSystemId(2),
            atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1, 1]),
        },
        true,
    )]
    #[case::different_atom_order(
        GraphIrAromaticSystemDelta::Add {
            id: GraphIrAromaticSystemId(2),
            atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1, 1]),
        },
        GraphIrAromaticSystemDelta::Add {
            id: GraphIrAromaticSystemId(2),
            atoms: vec![GraphIrAtomId(2), GraphIrAtomId(4), GraphIrAtomId(4)],
            attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1, 1]),
        },
        false,
    )]
    fn test_aromatic_system_delta_eq(
        #[case] lhs: GraphIrAromaticSystemDelta,
        #[case] rhs: GraphIrAromaticSystemDelta,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = AromaticSystemDelta::from_rust(py, &lhs).unwrap();
            let rhs = AromaticSystemDelta::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::add(
        GraphIrAromaticSystemDelta::Add {
            id: GraphIrAromaticSystemId(2),
            atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1, 1]),
        },
        "AromaticSystemDelta.Add(id=2, atoms=[4, 2, 4], attributes=AromaticSystemForm.parse('[1,1,1]'))",
    )]
    #[case::remove(
        GraphIrAromaticSystemDelta::Remove {
            id: GraphIrAromaticSystemId(2),
            atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1, 1]),
        },
        "AromaticSystemDelta.Remove(id=2, atoms=[4, 2, 4], attributes=AromaticSystemForm.parse('[1,1,1]'))",
    )]
    #[case::modify_field(
        GraphIrAromaticSystemDelta::ModifyField {
            id: GraphIrAromaticSystemId(2),
            change: GraphIrAromaticSystemFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(-1),
            },
        },
        "AromaticSystemDelta.ModifyField(id=2, change=AromaticSystemFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(-1)))",
    )]
    #[case::modify_constraint(
        GraphIrAromaticSystemDelta::ModifyConstraint {
            id: GraphIrAromaticSystemId(2),
            old: None,
            new: Some(GraphIrAromaticSystemConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
        },
        "AromaticSystemDelta.ModifyConstraint(id=2, old=None, new=AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(6)))",
    )]
    fn test_aromatic_system_delta_repr(
        #[case] delta: GraphIrAromaticSystemDelta,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let delta =
                into_py_variant(py, AromaticSystemDelta::from_rust(py, &delta).unwrap()).unwrap();
            assert_eq!(
                delta
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::add(GraphIrAromaticSystemDelta::Add {
        id: GraphIrAromaticSystemId(2),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1, 1]),
    })]
    #[case::remove(GraphIrAromaticSystemDelta::Remove {
        id: GraphIrAromaticSystemId(2),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1, 1]),
    })]
    #[case::modify_field(GraphIrAromaticSystemDelta::ModifyField {
        id: GraphIrAromaticSystemId(2),
        change: GraphIrAromaticSystemFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(-1),
        },
    })]
    #[case::constraint_added(GraphIrAromaticSystemDelta::ModifyConstraint {
        id: GraphIrAromaticSystemId(2),
        old: None,
        new: Some(GraphIrAromaticSystemConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
    })]
    #[case::constraint_removed(GraphIrAromaticSystemDelta::ModifyConstraint {
        id: GraphIrAromaticSystemId(2),
        old: Some(GraphIrAromaticSystemConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrAromaticSystemDelta::ModifyConstraint {
        id: GraphIrAromaticSystemId(2),
        old: Some(GraphIrAromaticSystemConstraintForm::ElectronCount(GraphIrNumForm::Lit(5))),
        new: Some(GraphIrAromaticSystemConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
    })]
    fn test_aromatic_system_delta_inverse(#[case] delta: GraphIrAromaticSystemDelta) {
        Python::attach(|py| {
            let binding = AromaticSystemDelta::from_rust(py, &delta).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                delta.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::add(GraphIrMulticenterBondDelta::Add {
        id: GraphIrMulticenterBondId(3),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
    })]
    #[case::remove(GraphIrMulticenterBondDelta::Remove {
        id: GraphIrMulticenterBondId(3),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
    })]
    #[case::modify_field(GraphIrMulticenterBondDelta::ModifyField {
        id: GraphIrMulticenterBondId(3),
        change: GraphIrMulticenterBondFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(-1),
        },
    })]
    #[case::constraint_added(GraphIrMulticenterBondDelta::ModifyConstraint {
        id: GraphIrMulticenterBondId(3),
        old: None,
        new: Some(GraphIrMulticenterBondConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
    })]
    #[case::constraint_removed(GraphIrMulticenterBondDelta::ModifyConstraint {
        id: GraphIrMulticenterBondId(3),
        old: Some(GraphIrMulticenterBondConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrMulticenterBondDelta::ModifyConstraint {
        id: GraphIrMulticenterBondId(3),
        old: Some(GraphIrMulticenterBondConstraintForm::ElectronCount(GraphIrNumForm::Lit(5))),
        new: Some(GraphIrMulticenterBondConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
    })]
    fn test_multicenter_bond_delta_roundtrip(#[case] delta: GraphIrMulticenterBondDelta) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterBondDelta::from_rust(py, &delta)
                    .unwrap()
                    .to_rust(py),
                delta
            );
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrMulticenterBondDelta::Add {
            id: GraphIrMulticenterBondId(3),
            atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
        },
        GraphIrMulticenterBondDelta::Add {
            id: GraphIrMulticenterBondId(3),
            atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
        },
        true,
    )]
    #[case::different_atom_order(
        GraphIrMulticenterBondDelta::Add {
            id: GraphIrMulticenterBondId(3),
            atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
        },
        GraphIrMulticenterBondDelta::Add {
            id: GraphIrMulticenterBondId(3),
            atoms: vec![GraphIrAtomId(2), GraphIrAtomId(4), GraphIrAtomId(4)],
            attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
        },
        false,
    )]
    fn test_multicenter_bond_delta_eq(
        #[case] lhs: GraphIrMulticenterBondDelta,
        #[case] rhs: GraphIrMulticenterBondDelta,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = MulticenterBondDelta::from_rust(py, &lhs).unwrap();
            let rhs = MulticenterBondDelta::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::add(
        GraphIrMulticenterBondDelta::Add {
            id: GraphIrMulticenterBondId(3),
            atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
        },
        "MulticenterBondDelta.Add(id=3, atoms=[4, 2, 4], attributes=MulticenterBondForm.parse('[1,1,1]'))",
    )]
    #[case::remove(
        GraphIrMulticenterBondDelta::Remove {
            id: GraphIrMulticenterBondId(3),
            atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
            attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
        },
        "MulticenterBondDelta.Remove(id=3, atoms=[4, 2, 4], attributes=MulticenterBondForm.parse('[1,1,1]'))",
    )]
    #[case::modify_field(
        GraphIrMulticenterBondDelta::ModifyField {
            id: GraphIrMulticenterBondId(3),
            change: GraphIrMulticenterBondFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(-1),
            },
        },
        "MulticenterBondDelta.ModifyField(id=3, change=MulticenterBondFieldChange.Charge(old=NumForm.Lit(0), new=NumForm.Lit(-1)))",
    )]
    #[case::modify_constraint(
        GraphIrMulticenterBondDelta::ModifyConstraint {
            id: GraphIrMulticenterBondId(3),
            old: None,
            new: Some(GraphIrMulticenterBondConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
        },
        "MulticenterBondDelta.ModifyConstraint(id=3, old=None, new=MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(6)))",
    )]
    fn test_multicenter_bond_delta_repr(
        #[case] delta: GraphIrMulticenterBondDelta,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let delta =
                into_py_variant(py, MulticenterBondDelta::from_rust(py, &delta).unwrap()).unwrap();
            assert_eq!(
                delta
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::add(GraphIrMulticenterBondDelta::Add {
        id: GraphIrMulticenterBondId(3),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
    })]
    #[case::remove(GraphIrMulticenterBondDelta::Remove {
        id: GraphIrMulticenterBondId(3),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2), GraphIrAtomId(4)],
        attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
    })]
    #[case::modify_field(GraphIrMulticenterBondDelta::ModifyField {
        id: GraphIrMulticenterBondId(3),
        change: GraphIrMulticenterBondFieldChange::Charge {
            old: GraphIrNumForm::Lit(0),
            new: GraphIrNumForm::Lit(-1),
        },
    })]
    #[case::constraint_added(GraphIrMulticenterBondDelta::ModifyConstraint {
        id: GraphIrMulticenterBondId(3),
        old: None,
        new: Some(GraphIrMulticenterBondConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
    })]
    #[case::constraint_removed(GraphIrMulticenterBondDelta::ModifyConstraint {
        id: GraphIrMulticenterBondId(3),
        old: Some(GraphIrMulticenterBondConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrMulticenterBondDelta::ModifyConstraint {
        id: GraphIrMulticenterBondId(3),
        old: Some(GraphIrMulticenterBondConstraintForm::ElectronCount(GraphIrNumForm::Lit(5))),
        new: Some(GraphIrMulticenterBondConstraintForm::ElectronCount(GraphIrNumForm::Lit(6))),
    })]
    fn test_multicenter_bond_delta_inverse(#[case] delta: GraphIrMulticenterBondDelta) {
        Python::attach(|py| {
            let binding = MulticenterBondDelta::from_rust(py, &delta).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                delta.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::add(GraphIrNoncovalentBondDelta::Add {
        id: GraphIrNoncovalentBondId(4),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
        attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
    })]
    #[case::remove(GraphIrNoncovalentBondDelta::Remove {
        id: GraphIrNoncovalentBondId(4),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
        attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
    })]
    #[case::modify_field(GraphIrNoncovalentBondDelta::ModifyField {
        id: GraphIrNoncovalentBondId(4),
        change: GraphIrNoncovalentBondFieldChange::Kind {
            old: GraphIrNoncovalentBondKindForm::Undetermined,
            new: GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond),
        },
    })]
    #[case::constraint_added(GraphIrNoncovalentBondDelta::ModifyConstraint {
        id: GraphIrNoncovalentBondId(4),
        old: None,
        new: Some(GraphIrNoncovalentBondConstraintForm::Intramolecular(GraphIrBooleanForm::Lit(true))),
    })]
    #[case::constraint_removed(GraphIrNoncovalentBondDelta::ModifyConstraint {
        id: GraphIrNoncovalentBondId(4),
        old: Some(GraphIrNoncovalentBondConstraintForm::Intramolecular(GraphIrBooleanForm::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrNoncovalentBondDelta::ModifyConstraint {
        id: GraphIrNoncovalentBondId(4),
        old: Some(GraphIrNoncovalentBondConstraintForm::Intramolecular(GraphIrBooleanForm::Lit(false))),
        new: Some(GraphIrNoncovalentBondConstraintForm::Intramolecular(GraphIrBooleanForm::Lit(true))),
    })]
    fn test_noncovalent_bond_delta_roundtrip(#[case] delta: GraphIrNoncovalentBondDelta) {
        Python::attach(|py| {
            assert_eq!(
                NoncovalentBondDelta::from_rust(py, &delta)
                    .unwrap()
                    .to_rust(py),
                delta
            );
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrNoncovalentBondDelta::Add {
            id: GraphIrNoncovalentBondId(4),
            atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
            attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
        },
        GraphIrNoncovalentBondDelta::Add {
            id: GraphIrNoncovalentBondId(4),
            atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
            attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
        },
        true,
    )]
    #[case::different_atom_order(
        GraphIrNoncovalentBondDelta::Add {
            id: GraphIrNoncovalentBondId(4),
            atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
            attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
        },
        GraphIrNoncovalentBondDelta::Add {
            id: GraphIrNoncovalentBondId(4),
            atoms: [GraphIrAtomId(2), GraphIrAtomId(5)],
            attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
        },
        false,
    )]
    fn test_noncovalent_bond_delta_eq(
        #[case] lhs: GraphIrNoncovalentBondDelta,
        #[case] rhs: GraphIrNoncovalentBondDelta,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = NoncovalentBondDelta::from_rust(py, &lhs).unwrap();
            let rhs = NoncovalentBondDelta::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::add(
        GraphIrNoncovalentBondDelta::Add {
            id: GraphIrNoncovalentBondId(4),
            atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
            attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
        },
        "NoncovalentBondDelta.Add(id=4, atoms=(5, 2), attributes=NoncovalentBondForm.parse('Hbd'))",
    )]
    #[case::remove(
        GraphIrNoncovalentBondDelta::Remove {
            id: GraphIrNoncovalentBondId(4),
            atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
            attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
        },
        "NoncovalentBondDelta.Remove(id=4, atoms=(5, 2), attributes=NoncovalentBondForm.parse('Hbd'))",
    )]
    #[case::modify_field(
        GraphIrNoncovalentBondDelta::ModifyField {
            id: GraphIrNoncovalentBondId(4),
            change: GraphIrNoncovalentBondFieldChange::Kind {
                old: GraphIrNoncovalentBondKindForm::Undetermined,
                new: GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond),
            },
        },
        "NoncovalentBondDelta.ModifyField(id=4, change=NoncovalentBondFieldChange.Kind(old=NoncovalentBondKindForm.Undetermined(), new=NoncovalentBondKindForm.Lit(NoncovalentBondKind.HydrogenBond)))",
    )]
    #[case::modify_constraint(
        GraphIrNoncovalentBondDelta::ModifyConstraint {
            id: GraphIrNoncovalentBondId(4),
            old: None,
            new: Some(GraphIrNoncovalentBondConstraintForm::Intramolecular(GraphIrBooleanForm::Lit(true))),
        },
        "NoncovalentBondDelta.ModifyConstraint(id=4, old=None, new=NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Lit(True)))",
    )]
    fn test_noncovalent_bond_delta_repr(
        #[case] delta: GraphIrNoncovalentBondDelta,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let delta =
                into_py_variant(py, NoncovalentBondDelta::from_rust(py, &delta).unwrap()).unwrap();
            assert_eq!(
                delta
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::add(GraphIrNoncovalentBondDelta::Add {
        id: GraphIrNoncovalentBondId(4),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
        attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
    })]
    #[case::remove(GraphIrNoncovalentBondDelta::Remove {
        id: GraphIrNoncovalentBondId(4),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
        attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
    })]
    #[case::modify_field(GraphIrNoncovalentBondDelta::ModifyField {
        id: GraphIrNoncovalentBondId(4),
        change: GraphIrNoncovalentBondFieldChange::Kind {
            old: GraphIrNoncovalentBondKindForm::Undetermined,
            new: GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond),
        },
    })]
    #[case::constraint_added(GraphIrNoncovalentBondDelta::ModifyConstraint {
        id: GraphIrNoncovalentBondId(4),
        old: None,
        new: Some(GraphIrNoncovalentBondConstraintForm::Intramolecular(GraphIrBooleanForm::Lit(true))),
    })]
    #[case::constraint_removed(GraphIrNoncovalentBondDelta::ModifyConstraint {
        id: GraphIrNoncovalentBondId(4),
        old: Some(GraphIrNoncovalentBondConstraintForm::Intramolecular(GraphIrBooleanForm::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(GraphIrNoncovalentBondDelta::ModifyConstraint {
        id: GraphIrNoncovalentBondId(4),
        old: Some(GraphIrNoncovalentBondConstraintForm::Intramolecular(GraphIrBooleanForm::Lit(false))),
        new: Some(GraphIrNoncovalentBondConstraintForm::Intramolecular(GraphIrBooleanForm::Lit(true))),
    })]
    fn test_noncovalent_bond_delta_inverse(#[case] delta: GraphIrNoncovalentBondDelta) {
        Python::attach(|py| {
            let binding = NoncovalentBondDelta::from_rust(py, &delta).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                delta.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::add(GraphIrStereoAtomDelta::Add {
        id: GraphIrStereoAtomId(5),
        site: GraphIrAtomId(3),
        ligands: vec![
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
        ],
        attributes: GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, GraphIrStereoCoset::Lit(0)),
    })]
    #[case::remove(GraphIrStereoAtomDelta::Remove {
        id: GraphIrStereoAtomId(5),
        site: GraphIrAtomId(3),
        ligands: vec![
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
        ],
        attributes: GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, GraphIrStereoCoset::Lit(0)),
    })]
    #[case::modify_field(GraphIrStereoAtomDelta::ModifyField {
        id: GraphIrStereoAtomId(5),
        change: GraphIrStereoAtomFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(0),
            ),
        },
    })]
    #[case::constraint_added_with_kind(GraphIrStereoAtomDelta::ModifyConstraint {
        id: GraphIrStereoAtomId(5),
        kind: Some(GraphIrStereoKind::Tetrahedral),
        old: None,
        new: Some(GraphIrStereoAtomConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
        )),
    })]
    #[case::constraint_removed_without_kind(GraphIrStereoAtomDelta::ModifyConstraint {
        id: GraphIrStereoAtomId(5),
        kind: None,
        old: Some(GraphIrStereoAtomConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Undetermined,
        )),
        new: None,
    })]
    #[case::constraint_modified(GraphIrStereoAtomDelta::ModifyConstraint {
        id: GraphIrStereoAtomId(5),
        kind: Some(GraphIrStereoKind::Tetrahedral),
        old: Some(GraphIrStereoAtomConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Undetermined,
        )),
        new: Some(GraphIrStereoAtomConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
        )),
    })]
    fn test_stereo_atom_delta_roundtrip(#[case] delta: GraphIrStereoAtomDelta) {
        Python::attach(|py| {
            assert_eq!(
                StereoAtomDelta::from_rust(py, &delta).unwrap().to_rust(py),
                delta
            );
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrStereoAtomDelta::ModifyField {
            id: GraphIrStereoAtomId(5),
            change: GraphIrStereoAtomFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                ),
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(1),
                ),
            },
        },
        GraphIrStereoAtomDelta::ModifyField {
            id: GraphIrStereoAtomId(5),
            change: GraphIrStereoAtomFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                ),
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(1),
                ),
            },
        },
        true,
    )]
    #[case::different_ligand_order(
        GraphIrStereoAtomDelta::Add {
            id: GraphIrStereoAtomId(5),
            site: GraphIrAtomId(3),
            ligands: vec![
                GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
            ],
            attributes: GraphIrStereoAtomForm::new(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(0),
            ),
        },
        GraphIrStereoAtomDelta::Add {
            id: GraphIrStereoAtomId(5),
            site: GraphIrAtomId(3),
            ligands: vec![
                GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            ],
            attributes: GraphIrStereoAtomForm::new(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(0),
            ),
        },
        false,
    )]
    #[case::different_configuration(
        GraphIrStereoAtomDelta::ModifyField {
            id: GraphIrStereoAtomId(5),
            change: GraphIrStereoAtomFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                ),
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(1),
                ),
            },
        },
        GraphIrStereoAtomDelta::ModifyField {
            id: GraphIrStereoAtomId(5),
            change: GraphIrStereoAtomFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                ),
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                ),
            },
        },
        false,
    )]
    fn test_stereo_atom_delta_eq(
        #[case] lhs: GraphIrStereoAtomDelta,
        #[case] rhs: GraphIrStereoAtomDelta,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = StereoAtomDelta::from_rust(py, &lhs).unwrap();
            let rhs = StereoAtomDelta::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::add(
        GraphIrStereoAtomDelta::Add {
            id: GraphIrStereoAtomId(5),
            site: GraphIrAtomId(3),
            ligands: vec![
                GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::LonePair),
            ],
            attributes: GraphIrStereoAtomForm::new(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(0),
            ),
        },
        "StereoAtomDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair)], attributes=StereoAtomForm.parse('Th0'))",
    )]
    #[case::remove(
        GraphIrStereoAtomDelta::Remove {
            id: GraphIrStereoAtomId(5),
            site: GraphIrAtomId(3),
            ligands: vec![
                GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::LonePair),
            ],
            attributes: GraphIrStereoAtomForm::new(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(0),
            ),
        },
        "StereoAtomDelta.Remove(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair)], attributes=StereoAtomForm.parse('Th0'))",
    )]
    #[case::modify_field(
        GraphIrStereoAtomDelta::ModifyField {
            id: GraphIrStereoAtomId(5),
            change: GraphIrStereoAtomFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Undetermined,
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                ),
            },
        },
        "StereoAtomDelta.ModifyField(id=5, change=StereoAtomFieldChange.Configuration(old=StereoConfigurationForm.Undetermined(), new=StereoConfigurationForm.Kinded(StereoKind.Tetrahedral, StereoCoset.Lit(0))))",
    )]
    #[case::modify_constraint(
        GraphIrStereoAtomDelta::ModifyConstraint {
            id: GraphIrStereoAtomId(5),
            kind: Some(GraphIrStereoKind::Tetrahedral),
            old: None,
            new: Some(GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Undetermined,
            )),
        },
        "StereoAtomDelta.ModifyConstraint(id=5, kind=StereoKind.Tetrahedral, old=None, new=StereoAtomConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()))",
    )]
    fn test_stereo_atom_delta_repr(#[case] delta: GraphIrStereoAtomDelta, #[case] expected: &str) {
        Python::attach(|py| {
            let delta =
                into_py_variant(py, StereoAtomDelta::from_rust(py, &delta).unwrap()).unwrap();
            assert_eq!(
                delta
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::add(GraphIrStereoAtomDelta::Add {
        id: GraphIrStereoAtomId(5),
        site: GraphIrAtomId(3),
        ligands: vec![
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
        ],
        attributes: GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, GraphIrStereoCoset::Lit(0)),
    })]
    #[case::remove(GraphIrStereoAtomDelta::Remove {
        id: GraphIrStereoAtomId(5),
        site: GraphIrAtomId(3),
        ligands: vec![
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
        ],
        attributes: GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, GraphIrStereoCoset::Lit(0)),
    })]
    #[case::modify_field(GraphIrStereoAtomDelta::ModifyField {
        id: GraphIrStereoAtomId(5),
        change: GraphIrStereoAtomFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(0),
            ),
        },
    })]
    #[case::constraint_added_with_kind(GraphIrStereoAtomDelta::ModifyConstraint {
        id: GraphIrStereoAtomId(5),
        kind: Some(GraphIrStereoKind::Tetrahedral),
        old: None,
        new: Some(GraphIrStereoAtomConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
        )),
    })]
    #[case::constraint_removed_without_kind(GraphIrStereoAtomDelta::ModifyConstraint {
        id: GraphIrStereoAtomId(5),
        kind: None,
        old: Some(GraphIrStereoAtomConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Undetermined,
        )),
        new: None,
    })]
    #[case::constraint_modified(GraphIrStereoAtomDelta::ModifyConstraint {
        id: GraphIrStereoAtomId(5),
        kind: Some(GraphIrStereoKind::Tetrahedral),
        old: Some(GraphIrStereoAtomConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Undetermined,
        )),
        new: Some(GraphIrStereoAtomConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
        )),
    })]
    fn test_stereo_atom_delta_inverse(#[case] delta: GraphIrStereoAtomDelta) {
        Python::attach(|py| {
            let binding = StereoAtomDelta::from_rust(py, &delta).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                delta.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), delta);
        });
    }
    #[rstest]
    #[case::add(GraphIrStereoBondDelta::Add {
        id: GraphIrStereoBondId(5),
        site: GraphIrBondId(3),
        ligands: vec![
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
        ],
        attributes: GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(0)),
    })]
    #[case::remove(GraphIrStereoBondDelta::Remove {
        id: GraphIrStereoBondId(5),
        site: GraphIrBondId(3),
        ligands: vec![
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
        ],
        attributes: GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(0)),
    })]
    #[case::modify_field(GraphIrStereoBondDelta::ModifyField {
        id: GraphIrStereoBondId(5),
        change: GraphIrStereoBondFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Lit(0),
            ),
        },
    })]
    #[case::constraint_added_with_kind(GraphIrStereoBondDelta::ModifyConstraint {
        id: GraphIrStereoBondId(5),
        kind: Some(GraphIrStereoKind::CisTrans),
        old: None,
        new: Some(GraphIrStereoBondConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
        )),
    })]
    #[case::constraint_removed_without_kind(GraphIrStereoBondDelta::ModifyConstraint {
        id: GraphIrStereoBondId(5),
        kind: None,
        old: Some(GraphIrStereoBondConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Undetermined,
        )),
        new: None,
    })]
    #[case::constraint_modified(GraphIrStereoBondDelta::ModifyConstraint {
        id: GraphIrStereoBondId(5),
        kind: Some(GraphIrStereoKind::CisTrans),
        old: Some(GraphIrStereoBondConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Undetermined,
        )),
        new: Some(GraphIrStereoBondConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
        )),
    })]
    fn test_stereo_bond_delta_roundtrip(#[case] delta: GraphIrStereoBondDelta) {
        Python::attach(|py| {
            assert_eq!(
                StereoBondDelta::from_rust(py, &delta).unwrap().to_rust(py),
                delta
            );
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrStereoBondDelta::ModifyField {
            id: GraphIrStereoBondId(5),
            change: GraphIrStereoBondFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(0),
                ),
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(1),
                ),
            },
        },
        GraphIrStereoBondDelta::ModifyField {
            id: GraphIrStereoBondId(5),
            change: GraphIrStereoBondFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(0),
                ),
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(1),
                ),
            },
        },
        true,
    )]
    #[case::different_ligand_order(
        GraphIrStereoBondDelta::Add {
            id: GraphIrStereoBondId(5),
            site: GraphIrBondId(3),
            ligands: vec![
                GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
            ],
            attributes: GraphIrStereoBondForm::new(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Lit(0),
            ),
        },
        GraphIrStereoBondDelta::Add {
            id: GraphIrStereoBondId(5),
            site: GraphIrBondId(3),
            ligands: vec![
                GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            ],
            attributes: GraphIrStereoBondForm::new(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Lit(0),
            ),
        },
        false,
    )]
    #[case::different_configuration(
        GraphIrStereoBondDelta::ModifyField {
            id: GraphIrStereoBondId(5),
            change: GraphIrStereoBondFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(0),
                ),
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(1),
                ),
            },
        },
        GraphIrStereoBondDelta::ModifyField {
            id: GraphIrStereoBondId(5),
            change: GraphIrStereoBondFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(0),
                ),
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(0),
                ),
            },
        },
        false,
    )]
    fn test_stereo_bond_delta_eq(
        #[case] lhs: GraphIrStereoBondDelta,
        #[case] rhs: GraphIrStereoBondDelta,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = StereoBondDelta::from_rust(py, &lhs).unwrap();
            let rhs = StereoBondDelta::from_rust(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::add(
        GraphIrStereoBondDelta::Add {
            id: GraphIrStereoBondId(5),
            site: GraphIrBondId(3),
            ligands: vec![
                GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::LonePair),
            ],
            attributes: GraphIrStereoBondForm::new(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Lit(0),
            ),
        },
        "StereoBondDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair)], attributes=StereoBondForm.parse('Ct0'))",
    )]
    #[case::remove(
        GraphIrStereoBondDelta::Remove {
            id: GraphIrStereoBondId(5),
            site: GraphIrBondId(3),
            ligands: vec![
                GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::LonePair),
            ],
            attributes: GraphIrStereoBondForm::new(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Lit(0),
            ),
        },
        "StereoBondDelta.Remove(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair)], attributes=StereoBondForm.parse('Ct0'))",
    )]
    #[case::modify_field(
        GraphIrStereoBondDelta::ModifyField {
            id: GraphIrStereoBondId(5),
            change: GraphIrStereoBondFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Undetermined,
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(0),
                ),
            },
        },
        "StereoBondDelta.ModifyField(id=5, change=StereoBondFieldChange.Configuration(old=StereoConfigurationForm.Undetermined(), new=StereoConfigurationForm.Kinded(StereoKind.CisTrans, StereoCoset.Lit(0))))",
    )]
    #[case::modify_constraint(
        GraphIrStereoBondDelta::ModifyConstraint {
            id: GraphIrStereoBondId(5),
            kind: Some(GraphIrStereoKind::CisTrans),
            old: None,
            new: Some(GraphIrStereoBondConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Undetermined,
            )),
        },
        "StereoBondDelta.ModifyConstraint(id=5, kind=StereoKind.CisTrans, old=None, new=StereoBondConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()))",
    )]
    fn test_stereo_bond_delta_repr(#[case] delta: GraphIrStereoBondDelta, #[case] expected: &str) {
        Python::attach(|py| {
            let delta =
                into_py_variant(py, StereoBondDelta::from_rust(py, &delta).unwrap()).unwrap();
            assert_eq!(
                delta
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::add(GraphIrStereoBondDelta::Add {
        id: GraphIrStereoBondId(5),
        site: GraphIrBondId(3),
        ligands: vec![
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
        ],
        attributes: GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(0)),
    })]
    #[case::remove(GraphIrStereoBondDelta::Remove {
        id: GraphIrStereoBondId(5),
        site: GraphIrBondId(3),
        ligands: vec![
            GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
        ],
        attributes: GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(0)),
    })]
    #[case::modify_field(GraphIrStereoBondDelta::ModifyField {
        id: GraphIrStereoBondId(5),
        change: GraphIrStereoBondFieldChange::Configuration {
            old: GraphIrStereoConfigurationForm::Undetermined,
            new: GraphIrStereoConfigurationForm::Kinded(
                GraphIrStereoKind::CisTrans,
                GraphIrStereoCoset::Lit(0),
            ),
        },
    })]
    #[case::constraint_added_with_kind(GraphIrStereoBondDelta::ModifyConstraint {
        id: GraphIrStereoBondId(5),
        kind: Some(GraphIrStereoKind::CisTrans),
        old: None,
        new: Some(GraphIrStereoBondConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
        )),
    })]
    #[case::constraint_removed_without_kind(GraphIrStereoBondDelta::ModifyConstraint {
        id: GraphIrStereoBondId(5),
        kind: None,
        old: Some(GraphIrStereoBondConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Undetermined,
        )),
        new: None,
    })]
    #[case::constraint_modified(GraphIrStereoBondDelta::ModifyConstraint {
        id: GraphIrStereoBondId(5),
        kind: Some(GraphIrStereoKind::CisTrans),
        old: Some(GraphIrStereoBondConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Undetermined,
        )),
        new: Some(GraphIrStereoBondConstraintForm::Stereogenicity(
            GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
        )),
    })]
    fn test_stereo_bond_delta_inverse(#[case] delta: GraphIrStereoBondDelta) {
        Python::attach(|py| {
            let binding = StereoBondDelta::from_rust(py, &delta).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                delta.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::add_leaf(GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
        GraphIrAtomId(3),
        GraphIrAtomConstraintForm::degree(2),
    )))]
    #[case::remove_recursive(GraphIrConstraintDelta::Remove(GraphIrConstraint::And(vec![
        GraphIrConstraint::Atom(GraphIrAtomId(7), GraphIrAtomConstraintForm::valence(4)),
        GraphIrConstraint::Not(Box::new(GraphIrConstraint::Or(Vec::new()))),
    ])))]
    fn test_constraint_delta_roundtrip(#[case] delta: GraphIrConstraintDelta) {
        Python::attach(|py| {
            let binding = ConstraintDelta::from_rust(py, &delta).unwrap();
            assert_eq!(binding.to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
            GraphIrAtomId(3),
            GraphIrAtomConstraintForm::degree(2),
        )),
        GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
            GraphIrAtomId(3),
            GraphIrAtomConstraintForm::degree(2),
        )),
        true
    )]
    #[case::variant(
        GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
            GraphIrAtomId(3),
            GraphIrAtomConstraintForm::degree(2),
        )),
        GraphIrConstraintDelta::Remove(GraphIrConstraint::Atom(
            GraphIrAtomId(3),
            GraphIrAtomConstraintForm::degree(2),
        )),
        false
    )]
    #[case::constraint(
        GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
            GraphIrAtomId(3),
            GraphIrAtomConstraintForm::degree(2),
        )),
        GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
            GraphIrAtomId(3),
            GraphIrAtomConstraintForm::valence(2),
        )),
        false
    )]
    fn test_constraint_delta_eq(
        #[case] left: GraphIrConstraintDelta,
        #[case] right: GraphIrConstraintDelta,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let left = ConstraintDelta::from_rust(py, &left).unwrap();
            let right = ConstraintDelta::from_rust(py, &right).unwrap();
            assert_eq!(left.__eq__(&right, py), expected);
        });
    }

    #[rstest]
    #[case::add_leaf(
        GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
            GraphIrAtomId(3),
            GraphIrAtomConstraintForm::degree(2),
        )),
        "ConstraintDelta.Add(constraint=Constraint.Atom(3, AtomConstraintForm.Degree(NumForm.Lit(2))))",
    )]
    #[case::remove_recursive(
        GraphIrConstraintDelta::Remove(GraphIrConstraint::And(vec![
            GraphIrConstraint::Atom(GraphIrAtomId(7), GraphIrAtomConstraintForm::valence(4)),
            GraphIrConstraint::Not(Box::new(GraphIrConstraint::Or(Vec::new()))),
        ])),
        "ConstraintDelta.Remove(constraint=Constraint.And([Constraint.Atom(7, AtomConstraintForm.Valence(NumForm.Lit(4))), Constraint.Not(Constraint.Or([]))]))",
    )]
    fn test_constraint_delta_repr(#[case] delta: GraphIrConstraintDelta, #[case] expected: &str) {
        Python::attach(|py| {
            let binding =
                into_py_variant(py, ConstraintDelta::from_rust(py, &delta).unwrap()).unwrap();
            assert_eq!(
                binding
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::add_leaf(GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
        GraphIrAtomId(3),
        GraphIrAtomConstraintForm::degree(2),
    )))]
    #[case::remove_recursive(GraphIrConstraintDelta::Remove(GraphIrConstraint::And(vec![
        GraphIrConstraint::Atom(GraphIrAtomId(7), GraphIrAtomConstraintForm::valence(4)),
        GraphIrConstraint::Not(Box::new(GraphIrConstraint::Or(Vec::new()))),
    ])))]
    fn test_constraint_delta_inverse(#[case] delta: GraphIrConstraintDelta) {
        Python::attach(|py| {
            let binding = ConstraintDelta::from_rust(py, &delta).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                delta.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::atom(GraphIrDelta::Atom(GraphIrAtomDelta::Add {
        id: GraphIrAtomId(3),
        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
    }))]
    #[case::bond(GraphIrDelta::Bond(GraphIrBondDelta::Add {
        id: GraphIrBondId(2),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
        attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
    }))]
    #[case::dative_bond(GraphIrDelta::DativeBond(GraphIrDativeBondDelta::Add {
        id: GraphIrDativeBondId(1),
        donors: vec![GraphIrAtomId(4), GraphIrAtomId(2)],
        acceptor: GraphIrAtomId(3),
        attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
    }))]
    #[case::aromatic_system(GraphIrDelta::AromaticSystem(GraphIrAromaticSystemDelta::Add {
        id: GraphIrAromaticSystemId(2),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2)],
        attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1]),
    }))]
    #[case::multicenter_bond(GraphIrDelta::MulticenterBond(GraphIrMulticenterBondDelta::Add {
        id: GraphIrMulticenterBondId(3),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2)],
        attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1]),
    }))]
    #[case::noncovalent_bond(GraphIrDelta::NoncovalentBond(GraphIrNoncovalentBondDelta::Add {
        id: GraphIrNoncovalentBondId(4),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
        attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
    }))]
    #[case::stereo_atom(GraphIrDelta::StereoAtom(GraphIrStereoAtomDelta::Add {
        id: GraphIrStereoAtomId(5),
        site: GraphIrAtomId(3),
        ligands: vec![GraphIrStereoLigand::new(
            GraphIrAtomId(4),
            GraphIrStereoLigandKind::Atom,
        )],
        attributes: GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, GraphIrStereoCoset::Lit(0)),
    }))]
    #[case::stereo_bond(GraphIrDelta::StereoBond(GraphIrStereoBondDelta::Add {
        id: GraphIrStereoBondId(5),
        site: GraphIrBondId(3),
        ligands: vec![GraphIrStereoLigand::new(
            GraphIrAtomId(4),
            GraphIrStereoLigandKind::Atom,
        )],
        attributes: GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(0)),
    }))]
    #[case::constraint(GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(
        GraphIrConstraint::Atom(GraphIrAtomId(3), GraphIrAtomConstraintForm::degree(2)),
    )))]
    fn test_delta_roundtrip(#[case] delta: GraphIrDelta) {
        Python::attach(|py| {
            let binding = Delta::from_rust(py, &delta).unwrap();
            assert_eq!(binding.to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        }),
        GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        }),
        true,
    )]
    #[case::outer_variant(
        GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        }),
        GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
            GraphIrAtomId(3),
            GraphIrAtomConstraintForm::degree(2),
        ))),
        false,
    )]
    #[case::child(
        GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        }),
        GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(4),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        }),
        false,
    )]
    fn test_delta_eq(
        #[case] left: GraphIrDelta,
        #[case] right: GraphIrDelta,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let left = Delta::from_rust(py, &left).unwrap();
            let right = Delta::from_rust(py, &right).unwrap();
            assert_eq!(left.__eq__(&right, py), expected);
        });
    }

    #[rstest]
    #[case::atom(
        GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        }),
        "Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm.parse('C')))"
    )]
    #[case::stereo_atom(
        GraphIrDelta::StereoAtom(GraphIrStereoAtomDelta::Add {
            id: GraphIrStereoAtomId(5),
            site: GraphIrAtomId(3),
            ligands: vec![GraphIrStereoLigand::new(
                GraphIrAtomId(4),
                GraphIrStereoLigandKind::Atom,
            )],
            attributes: GraphIrStereoAtomForm::new(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(0),
            ),
        }),
        "Delta.StereoAtom(StereoAtomDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], attributes=StereoAtomForm.parse('Th0')))"
    )]
    #[case::constraint(
        GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
            GraphIrAtomId(3),
            GraphIrAtomConstraintForm::degree(2),
        ))),
        "Delta.Constraint(ConstraintDelta.Add(constraint=Constraint.Atom(3, AtomConstraintForm.Degree(NumForm.Lit(2)))))"
    )]
    fn test_delta_repr(#[case] delta: GraphIrDelta, #[case] expected: &str) {
        Python::attach(|py| {
            let binding = into_py_variant(py, Delta::from_rust(py, &delta).unwrap()).unwrap();
            assert_eq!(
                binding
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::atom(GraphIrDelta::Atom(GraphIrAtomDelta::Add {
        id: GraphIrAtomId(3),
        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
    }))]
    #[case::bond(GraphIrDelta::Bond(GraphIrBondDelta::Add {
        id: GraphIrBondId(2),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(1)],
        attributes: GraphIrBondForm::new(GraphIrNumForm::Lit(1)),
    }))]
    #[case::dative_bond(GraphIrDelta::DativeBond(GraphIrDativeBondDelta::Add {
        id: GraphIrDativeBondId(1),
        donors: vec![GraphIrAtomId(4), GraphIrAtomId(2)],
        acceptor: GraphIrAtomId(3),
        attributes: GraphIrDativeBondForm::new(GraphIrNumForm::Lit(1)),
    }))]
    #[case::aromatic_system(GraphIrDelta::AromaticSystem(GraphIrAromaticSystemDelta::Add {
        id: GraphIrAromaticSystemId(2),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2)],
        attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1]),
    }))]
    #[case::multicenter_bond(GraphIrDelta::MulticenterBond(GraphIrMulticenterBondDelta::Add {
        id: GraphIrMulticenterBondId(3),
        atoms: vec![GraphIrAtomId(4), GraphIrAtomId(2)],
        attributes: GraphIrMulticenterBondForm::from_electrons(vec![1, 1]),
    }))]
    #[case::noncovalent_bond(GraphIrDelta::NoncovalentBond(GraphIrNoncovalentBondDelta::Add {
        id: GraphIrNoncovalentBondId(4),
        atoms: [GraphIrAtomId(5), GraphIrAtomId(2)],
        attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
    }))]
    #[case::stereo_atom(GraphIrDelta::StereoAtom(GraphIrStereoAtomDelta::Add {
        id: GraphIrStereoAtomId(5),
        site: GraphIrAtomId(3),
        ligands: vec![GraphIrStereoLigand::new(
            GraphIrAtomId(4),
            GraphIrStereoLigandKind::Atom,
        )],
        attributes: GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, GraphIrStereoCoset::Lit(0)),
    }))]
    #[case::stereo_bond(GraphIrDelta::StereoBond(GraphIrStereoBondDelta::Add {
        id: GraphIrStereoBondId(5),
        site: GraphIrBondId(3),
        ligands: vec![GraphIrStereoLigand::new(
            GraphIrAtomId(4),
            GraphIrStereoLigandKind::Atom,
        )],
        attributes: GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(0)),
    }))]
    #[case::constraint(GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(
        GraphIrConstraint::Atom(GraphIrAtomId(3), GraphIrAtomConstraintForm::degree(2)),
    )))]
    fn test_delta_inverse(#[case] delta: GraphIrDelta) {
        Python::attach(|py| {
            let binding = Delta::from_rust(py, &delta).unwrap();
            let inverse = binding.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_rust(py),
                delta.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::positive_first(3, 0, 0)]
    #[case::positive_last(3, 2, 2)]
    #[case::negative_last(3, -1, 2)]
    #[case::negative_first(3, -3, 0)]
    fn test_resolve_delta_index(#[case] len: usize, #[case] index: isize, #[case] expected: usize) {
        assert_eq!(resolve_delta_index(len, index).unwrap(), expected);
    }

    #[rstest]
    #[case::empty(0, 0)]
    #[case::positive(2, 2)]
    #[case::negative(2, -3)]
    fn test_resolve_delta_index_error(#[case] len: usize, #[case] index: isize) {
        Python::attach(|py| {
            let error = resolve_delta_index(len, index).unwrap_err();
            assert!(error.is_instance_of::<PyIndexError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "delta index out of range"
            );
        });
    }

    #[rstest]
    #[case::empty(Vec::new())]
    #[case::populated(vec![
        GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        }),
        GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
            GraphIrAtomId(3),
            GraphIrAtomConstraintForm::degree(2),
        ))),
    ])]
    fn test_deltas_new(#[case] entries: Vec<GraphIrDelta>) {
        Python::attach(|py| {
            let python_entries = entries
                .iter()
                .map(|entry| into_py_variant(py, Delta::from_rust(py, entry).unwrap()).unwrap())
                .collect();
            let expected: GraphIrDeltas = entries.into_iter().collect();
            assert_eq!(Deltas::new(py, python_entries).to_rust(), &expected);
        });
    }

    #[rstest]
    #[case::equal(
        vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        })],
        vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        })],
        true,
    )]
    #[case::different(
        vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        })],
        vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(4),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        })],
        false,
    )]
    fn test_deltas_eq(
        #[case] left: Vec<GraphIrDelta>,
        #[case] right: Vec<GraphIrDelta>,
        #[case] expected: bool,
    ) {
        assert_eq!(
            Deltas::from_rust(left.into_iter().collect())
                == Deltas::from_rust(right.into_iter().collect()),
            expected
        );
    }

    #[rstest]
    #[case::empty(Vec::new(), "Deltas([])")]
    #[case::populated(
        vec![
            GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                id: GraphIrAtomId(3),
                attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
            }),
            GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
                GraphIrAtomId(3),
                GraphIrAtomConstraintForm::degree(2),
            ))),
        ],
        "Deltas([Delta.Atom(AtomDelta.Add(id=3, attributes=AtomForm.parse('C'))), Delta.Constraint(ConstraintDelta.Add(constraint=Constraint.Atom(3, AtomConstraintForm.Degree(NumForm.Lit(2)))))])",
    )]
    fn test_deltas_repr(#[case] entries: Vec<GraphIrDelta>, #[case] expected: &str) {
        Python::attach(|py| {
            let deltas = Deltas::from_rust(entries.into_iter().collect());
            assert_eq!(deltas.__repr__(py).unwrap(), expected);
        });
    }

    #[rstest]
    fn test_deltas_append() {
        let appended = GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(
            GraphIrConstraint::Atom(GraphIrAtomId(3), GraphIrAtomConstraintForm::degree(2)),
        ));
        Python::attach(|py| {
            let mut deltas = Deltas::from_rust(
                vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                    id: GraphIrAtomId(3),
                    attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
                })]
                .into_iter()
                .collect(),
            );
            let value = into_py_variant(py, Delta::from_rust(py, &appended).unwrap()).unwrap();

            deltas.append(py, value);

            assert_eq!(
                deltas.to_rust().as_slice(),
                &[
                    GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(3),
                        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
                    }),
                    appended,
                ]
            );
        });
    }

    #[rstest]
    fn test_deltas_extend_container() {
        Python::attach(|py| {
            let target = Py::new(
                py,
                Deltas::from_rust(
                    vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(3),
                        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let source = Py::new(
                py,
                Deltas::from_rust(
                    vec![GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(
                        GraphIrConstraint::Atom(
                            GraphIrAtomId(3),
                            GraphIrAtomConstraintForm::degree(2),
                        ),
                    ))]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();

            Deltas::extend(target.clone_ref(py), py, DeltasExtend::Container(source));

            assert_eq!(
                target.bind(py).borrow().to_rust().as_slice(),
                &[
                    GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(3),
                        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
                    }),
                    GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
                        GraphIrAtomId(3),
                        GraphIrAtomConstraintForm::degree(2),
                    ))),
                ]
            );
        });
    }

    #[rstest]
    fn test_deltas_extend_entries() {
        Python::attach(|py| {
            let target = Py::new(py, Deltas::from_rust(GraphIrDeltas::new())).unwrap();
            let atom = GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                id: GraphIrAtomId(3),
                attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
            });
            let constraint = GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(
                GraphIrConstraint::Atom(GraphIrAtomId(3), GraphIrAtomConstraintForm::degree(2)),
            ));
            let entries = vec![
                into_py_variant(py, Delta::from_rust(py, &atom).unwrap()).unwrap(),
                into_py_variant(py, Delta::from_rust(py, &constraint).unwrap()).unwrap(),
            ];

            Deltas::extend(target.clone_ref(py), py, DeltasExtend::Entries(entries));

            assert_eq!(
                target.bind(py).borrow().to_rust().as_slice(),
                &[atom, constraint]
            );
        });
    }

    #[rstest]
    fn test_deltas_extend_self() {
        Python::attach(|py| {
            let atom = GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                id: GraphIrAtomId(3),
                attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
            });
            let constraint = GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(
                GraphIrConstraint::Atom(GraphIrAtomId(3), GraphIrAtomConstraintForm::degree(2)),
            ));
            let target = Py::new(
                py,
                Deltas::from_rust(vec![atom.clone(), constraint.clone()].into_iter().collect()),
            )
            .unwrap();

            Deltas::extend(
                target.clone_ref(py),
                py,
                DeltasExtend::Container(target.clone_ref(py)),
            );

            assert_eq!(
                target.bind(py).borrow().to_rust().as_slice(),
                &[atom.clone(), constraint.clone(), atom, constraint]
            );
        });
    }

    #[rstest]
    #[case::field_fusion(
        vec![
            GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                id: GraphIrAtomId(0),
                change: GraphIrAtomFieldChange::Charge {
                    old: GraphIrNumForm::Lit(0),
                    new: GraphIrNumForm::Lit(1),
                },
            }),
            GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                id: GraphIrAtomId(0),
                change: GraphIrAtomFieldChange::Charge {
                    old: GraphIrNumForm::Lit(1),
                    new: GraphIrNumForm::Lit(2),
                },
            }),
        ],
        vec![GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
            id: GraphIrAtomId(0),
            change: GraphIrAtomFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(2),
            },
        })],
    )]
    #[case::add_remove_cancellation(
        vec![
            GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                id: GraphIrAtomId(0),
                attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
            }),
            GraphIrDelta::Atom(GraphIrAtomDelta::Remove {
                id: GraphIrAtomId(0),
                attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
            }),
        ],
        Vec::new(),
    )]
    #[case::entity_kind_order(
        vec![
            GraphIrDelta::Bond(GraphIrBondDelta::ModifyField {
                id: GraphIrBondId(0),
                change: GraphIrBondFieldChange::Order {
                    old: GraphIrNumForm::Lit(1),
                    new: GraphIrNumForm::Lit(2),
                },
            }),
            GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                id: GraphIrAtomId(0),
                change: GraphIrAtomFieldChange::Charge {
                    old: GraphIrNumForm::Lit(0),
                    new: GraphIrNumForm::Lit(1),
                },
            }),
        ],
        vec![
            GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                id: GraphIrAtomId(0),
                change: GraphIrAtomFieldChange::Charge {
                    old: GraphIrNumForm::Lit(0),
                    new: GraphIrNumForm::Lit(1),
                },
            }),
            GraphIrDelta::Bond(GraphIrBondDelta::ModifyField {
                id: GraphIrBondId(0),
                change: GraphIrBondFieldChange::Order {
                    old: GraphIrNumForm::Lit(1),
                    new: GraphIrNumForm::Lit(2),
                },
            }),
        ],
    )]
    fn test_deltas_normalize(
        #[case] input: Vec<GraphIrDelta>,
        #[case] expected: Vec<GraphIrDelta>,
    ) {
        let source = Deltas::from_rust(input.into_iter().collect());
        let before = source.to_rust().clone();

        Python::attach(|py| {
            let normalized = source.normalize(py).unwrap();

            let expected: GraphIrDeltas = expected.into_iter().collect();
            assert_eq!(normalized.to_rust(), &expected);
            assert_eq!(source.to_rust(), &before);
            assert_eq!(normalized.normalize(py).unwrap(), normalized);
        });
    }

    #[rstest]
    fn test_deltas_normalize_error() {
        let source = Deltas::from_rust(
            vec![
                GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                    id: GraphIrAtomId(0),
                    change: GraphIrAtomFieldChange::Charge {
                        old: GraphIrNumForm::Lit(0),
                        new: GraphIrNumForm::Lit(1),
                    },
                }),
                GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                    id: GraphIrAtomId(0),
                    change: GraphIrAtomFieldChange::Charge {
                        old: GraphIrNumForm::Lit(2),
                        new: GraphIrNumForm::Lit(3),
                    },
                }),
            ]
            .into_iter()
            .collect(),
        );
        let before = source.to_rust();

        Python::attach(|py| {
            let error = source.normalize(py).err().unwrap();
            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reached a contradiction"
            );
        });
        assert_eq!(source.to_rust(), before);
    }

    #[rstest]
    #[case::empty(Vec::new(), 0)]
    #[case::populated(
        vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(3),
            attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
        })],
        1,
    )]
    fn test_deltas_len(#[case] entries: Vec<GraphIrDelta>, #[case] expected: usize) {
        assert_eq!(
            Deltas::from_rust(entries.into_iter().collect()).__len__(),
            expected
        );
    }

    #[rstest]
    #[case::positive(0, GraphIrDelta::Atom(GraphIrAtomDelta::Add {
        id: GraphIrAtomId(3),
        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
    }))]
    #[case::negative(-1, GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(
        GraphIrConstraint::Atom(GraphIrAtomId(3), GraphIrAtomConstraintForm::degree(2)),
    )))]
    fn test_deltas_getitem(#[case] index: isize, #[case] expected: GraphIrDelta) {
        Python::attach(|py| {
            let deltas = Deltas::from_rust(
                vec![
                    GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(3),
                        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
                    }),
                    GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
                        GraphIrAtomId(3),
                        GraphIrAtomConstraintForm::degree(2),
                    ))),
                ]
                .into_iter()
                .collect(),
            );
            assert_eq!(deltas.__getitem__(py, index).unwrap().to_rust(py), expected);
        });
    }

    #[rstest]
    #[case::positive(2)]
    #[case::negative(-3)]
    fn test_deltas_getitem_error(#[case] index: isize) {
        Python::attach(|py| {
            let deltas = Deltas::from_rust(
                vec![
                    GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(3),
                        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
                    }),
                    GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
                        GraphIrAtomId(3),
                        GraphIrAtomConstraintForm::degree(2),
                    ))),
                ]
                .into_iter()
                .collect(),
            );
            let error = deltas.__getitem__(py, index).err().unwrap();
            assert!(error.is_instance_of::<PyIndexError>(py));
        });
    }

    #[rstest]
    fn test_deltas_iter() {
        let expected = vec![
            GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                id: GraphIrAtomId(3),
                attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
            }),
            GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(GraphIrConstraint::Atom(
                GraphIrAtomId(3),
                GraphIrAtomConstraintForm::degree(2),
            ))),
        ];
        Python::attach(|py| {
            let deltas = Deltas::from_rust(expected.clone().into_iter().collect());
            let mut iter = deltas.__iter__(py).unwrap();
            assert_eq!(
                iter.__next__().unwrap().bind(py).borrow().to_rust(py),
                expected[0]
            );
            assert_eq!(
                iter.__next__().unwrap().bind(py).borrow().to_rust(py),
                expected[1]
            );
            assert!(iter.__next__().is_none());
        });
    }

    #[rstest]
    #[case::empty(Vec::new())]
    #[case::populated(vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
        id: GraphIrAtomId(3),
        attributes: GraphIrAtomForm::new(GraphIrElementForm::Lit(ChemElement::C)),
    })])]
    fn test_deltas_roundtrip(#[case] entries: Vec<GraphIrDelta>) {
        let rust: GraphIrDeltas = entries.into_iter().collect();
        assert_eq!(Deltas::from_rust(rust.clone()).to_rust(), &rust);
    }
}
