//! Python bindings for resolved molecule deltas and their field-change payloads.

use std::vec::IntoIter;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticSystemAst as AstAromaticSystemAst, AromaticSystemDelta as AstAromaticSystemDelta,
    AromaticSystemFieldChange as AstAromaticSystemFieldChange,
    AromaticSystemId as AstAromaticSystemId, AtomAst as AstAtomAst, AtomDelta as AstAtomDelta,
    AtomFieldChange as AstAtomFieldChange, AtomId as AstAtomId, BondAst as AstBondAst,
    BondDelta as AstBondDelta, BondFieldChange as AstBondFieldChange, BondId as AstBondId,
    Constraint as AstConstraint, ConstraintDelta as AstConstraintDelta,
    DativeBondAst as AstDativeBondAst, DativeBondDelta as AstDativeBondDelta,
    DativeBondFieldChange as AstDativeBondFieldChange, DativeBondId as AstDativeBondId,
    Delta as AstDelta, Deltas as AstDeltas, MulticenterBondAst as AstMulticenterBondAst,
    MulticenterBondDelta as AstMulticenterBondDelta,
    MulticenterBondFieldChange as AstMulticenterBondFieldChange,
    MulticenterBondId as AstMulticenterBondId, NoncovalentBondAst as AstNoncovalentBondAst,
    NoncovalentBondDelta as AstNoncovalentBondDelta,
    NoncovalentBondFieldChange as AstNoncovalentBondFieldChange,
    NoncovalentBondId as AstNoncovalentBondId, StereoAtomAst as AstStereoAtomAst,
    StereoAtomDelta as AstStereoAtomDelta, StereoAtomFieldChange as AstStereoAtomFieldChange,
    StereoAtomId as AstStereoAtomId, StereoBondAst as AstStereoBondAst,
    StereoBondDelta as AstStereoBondDelta, StereoBondFieldChange as AstStereoBondFieldChange,
    StereoBondId as AstStereoBondId,
};

use crate::aromatic::AromaticSystemAst;
use crate::atom::{AtomAst, ElementAst, IsotopeMassAst};
use crate::bond::BondAst;
use crate::constraint::aromatic::AromaticSystemConstraintAst;
use crate::constraint::atom::AtomConstraintAst;
use crate::constraint::bond::BondConstraintAst;
use crate::constraint::dative::DativeBondConstraintAst;
use crate::constraint::molecule::Constraint;
use crate::constraint::multicenter::MulticenterBondConstraintAst;
use crate::constraint::noncovalent::NoncovalentBondConstraintAst;
use crate::constraint::stereo::{StereoAtomConstraintAst, StereoBondConstraintAst};
use crate::convert::{into_py_variant, variant_repr};
use crate::dative::DativeBondAst;
use crate::electrons::ElectronCountsAst;
use crate::multicenter::MulticenterBondAst;
use crate::noncovalent::{NoncovalentBondAst, NoncovalentBondKindAst};
use crate::spin::SpinStateAst;
use crate::stereo::{
    Permutation, StereoAtomAst, StereoBondAst, StereoConfigurationAst, StereoKind, StereoLigand,
};
use crate::value::ValueAst;

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
        #[pyclass]
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
    /// An atom attribute change carrying the field's old and new AST values.
    AtomFieldChange {
        Element(ElementAst),
        IsotopeMass(IsotopeMassAst),
        Charge(ValueAst),
        ImplicitHydrogens(ValueAst),
        LonePairs(ValueAst),
        Spin(SpinStateAst),
    }
}

field_change! {
    /// A covalent-bond attribute change carrying the field's old and new AST values.
    BondFieldChange {
        Order(ValueAst),
        Charge(ValueAst),
        Spin(SpinStateAst),
    }
}

field_change! {
    /// A dative-bond attribute change carrying the field's old and new AST values.
    DativeBondFieldChange {
        Order(ValueAst),
    }
}

field_change! {
    /// An aromatic-system attribute change carrying the field's old and new AST values.
    AromaticSystemFieldChange {
        Electrons(ElectronCountsAst),
        Charge(ValueAst),
        Spin(SpinStateAst),
    }
}

field_change! {
    /// A multicenter-bond attribute change carrying the field's old and new AST values.
    MulticenterBondFieldChange {
        Electrons(ElectronCountsAst),
        Charge(ValueAst),
        Spin(SpinStateAst),
    }
}

field_change! {
    /// A noncovalent-bond kind change carrying the field's old and new AST values.
    NoncovalentBondFieldChange {
        Kind(NoncovalentBondKindAst),
    }
}

field_change! {
    /// A stereo-atom configuration change carrying the field's old and new AST values.
    StereoAtomFieldChange {
        Configuration(StereoConfigurationAst),
    }
}

field_change! {
    /// A stereo-bond configuration change carrying the field's old and new AST values.
    StereoBondFieldChange {
        Configuration(StereoConfigurationAst),
    }
}

impl AtomFieldChange {
    pub(crate) fn from_rust(py: Python<'_>, change: &AstAtomFieldChange) -> PyResult<Self> {
        Ok(match change {
            AstAtomFieldChange::Element { old, new } => Self::Element {
                old: into_py_variant(py, ElementAst::from_rust(old))?,
                new: into_py_variant(py, ElementAst::from_rust(new))?,
            },
            AstAtomFieldChange::IsotopeMass { old, new } => Self::IsotopeMass {
                old: into_py_variant(py, IsotopeMassAst::from_rust(old))?,
                new: into_py_variant(py, IsotopeMassAst::from_rust(new))?,
            },
            AstAtomFieldChange::Charge { old, new } => Self::Charge {
                old: into_py_variant(py, ValueAst::from_rust(py, old)?)?,
                new: into_py_variant(py, ValueAst::from_rust(py, new)?)?,
            },
            AstAtomFieldChange::ImplicitHydrogens { old, new } => Self::ImplicitHydrogens {
                old: into_py_variant(py, ValueAst::from_rust(py, old)?)?,
                new: into_py_variant(py, ValueAst::from_rust(py, new)?)?,
            },
            AstAtomFieldChange::LonePairs { old, new } => Self::LonePairs {
                old: into_py_variant(py, ValueAst::from_rust(py, old)?)?,
                new: into_py_variant(py, ValueAst::from_rust(py, new)?)?,
            },
            AstAtomFieldChange::Spin { old, new } => Self::Spin {
                old: Py::new(py, SpinStateAst::from_rust(py, old)?)?,
                new: Py::new(py, SpinStateAst::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstAtomFieldChange {
        match self {
            Self::Element { old, new } => AstAtomFieldChange::Element {
                old: old.bind(py).borrow().to_rust(),
                new: new.bind(py).borrow().to_rust(),
            },
            Self::IsotopeMass { old, new } => AstAtomFieldChange::IsotopeMass {
                old: old.bind(py).borrow().to_rust(),
                new: new.bind(py).borrow().to_rust(),
            },
            Self::Charge { old, new } => AstAtomFieldChange::Charge {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::ImplicitHydrogens { old, new } => AstAtomFieldChange::ImplicitHydrogens {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::LonePairs { old, new } => AstAtomFieldChange::LonePairs {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::Spin { old, new } => AstAtomFieldChange::Spin {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

impl BondFieldChange {
    pub(crate) fn from_rust(py: Python<'_>, change: &AstBondFieldChange) -> PyResult<Self> {
        Ok(match change {
            AstBondFieldChange::Order { old, new } => Self::Order {
                old: into_py_variant(py, ValueAst::from_rust(py, old)?)?,
                new: into_py_variant(py, ValueAst::from_rust(py, new)?)?,
            },
            AstBondFieldChange::Charge { old, new } => Self::Charge {
                old: into_py_variant(py, ValueAst::from_rust(py, old)?)?,
                new: into_py_variant(py, ValueAst::from_rust(py, new)?)?,
            },
            AstBondFieldChange::Spin { old, new } => Self::Spin {
                old: Py::new(py, SpinStateAst::from_rust(py, old)?)?,
                new: Py::new(py, SpinStateAst::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstBondFieldChange {
        match self {
            Self::Order { old, new } => AstBondFieldChange::Order {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::Charge { old, new } => AstBondFieldChange::Charge {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::Spin { old, new } => AstBondFieldChange::Spin {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

impl DativeBondFieldChange {
    pub(crate) fn from_rust(py: Python<'_>, change: &AstDativeBondFieldChange) -> PyResult<Self> {
        Ok(match change {
            AstDativeBondFieldChange::Order { old, new } => Self::Order {
                old: into_py_variant(py, ValueAst::from_rust(py, old)?)?,
                new: into_py_variant(py, ValueAst::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstDativeBondFieldChange {
        match self {
            Self::Order { old, new } => AstDativeBondFieldChange::Order {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

impl AromaticSystemFieldChange {
    pub(crate) fn from_rust(
        py: Python<'_>,
        change: &AstAromaticSystemFieldChange,
    ) -> PyResult<Self> {
        Ok(match change {
            AstAromaticSystemFieldChange::Electrons { old, new } => Self::Electrons {
                old: into_py_variant(py, ElectronCountsAst::from_rust(old))?,
                new: into_py_variant(py, ElectronCountsAst::from_rust(new))?,
            },
            AstAromaticSystemFieldChange::Charge { old, new } => Self::Charge {
                old: into_py_variant(py, ValueAst::from_rust(py, old)?)?,
                new: into_py_variant(py, ValueAst::from_rust(py, new)?)?,
            },
            AstAromaticSystemFieldChange::Spin { old, new } => Self::Spin {
                old: Py::new(py, SpinStateAst::from_rust(py, old)?)?,
                new: Py::new(py, SpinStateAst::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstAromaticSystemFieldChange {
        match self {
            Self::Electrons { old, new } => AstAromaticSystemFieldChange::Electrons {
                old: old.bind(py).borrow().to_rust(),
                new: new.bind(py).borrow().to_rust(),
            },
            Self::Charge { old, new } => AstAromaticSystemFieldChange::Charge {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::Spin { old, new } => AstAromaticSystemFieldChange::Spin {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

impl MulticenterBondFieldChange {
    pub(crate) fn from_rust(
        py: Python<'_>,
        change: &AstMulticenterBondFieldChange,
    ) -> PyResult<Self> {
        Ok(match change {
            AstMulticenterBondFieldChange::Electrons { old, new } => Self::Electrons {
                old: into_py_variant(py, ElectronCountsAst::from_rust(old))?,
                new: into_py_variant(py, ElectronCountsAst::from_rust(new))?,
            },
            AstMulticenterBondFieldChange::Charge { old, new } => Self::Charge {
                old: into_py_variant(py, ValueAst::from_rust(py, old)?)?,
                new: into_py_variant(py, ValueAst::from_rust(py, new)?)?,
            },
            AstMulticenterBondFieldChange::Spin { old, new } => Self::Spin {
                old: Py::new(py, SpinStateAst::from_rust(py, old)?)?,
                new: Py::new(py, SpinStateAst::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstMulticenterBondFieldChange {
        match self {
            Self::Electrons { old, new } => AstMulticenterBondFieldChange::Electrons {
                old: old.bind(py).borrow().to_rust(),
                new: new.bind(py).borrow().to_rust(),
            },
            Self::Charge { old, new } => AstMulticenterBondFieldChange::Charge {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
            Self::Spin { old, new } => AstMulticenterBondFieldChange::Spin {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

impl NoncovalentBondFieldChange {
    pub(crate) fn from_rust(
        py: Python<'_>,
        change: &AstNoncovalentBondFieldChange,
    ) -> PyResult<Self> {
        Ok(match change {
            AstNoncovalentBondFieldChange::Kind { old, new } => Self::Kind {
                old: into_py_variant(py, NoncovalentBondKindAst::from_rust(old))?,
                new: into_py_variant(py, NoncovalentBondKindAst::from_rust(new))?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstNoncovalentBondFieldChange {
        match self {
            Self::Kind { old, new } => AstNoncovalentBondFieldChange::Kind {
                old: old.bind(py).borrow().to_rust(),
                new: new.bind(py).borrow().to_rust(),
            },
        }
    }
}

impl StereoAtomFieldChange {
    pub(crate) fn from_rust(py: Python<'_>, change: &AstStereoAtomFieldChange) -> PyResult<Self> {
        Ok(match change {
            AstStereoAtomFieldChange::Configuration { old, new } => Self::Configuration {
                old: into_py_variant(py, StereoConfigurationAst::from_rust(py, old)?)?,
                new: into_py_variant(py, StereoConfigurationAst::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstStereoAtomFieldChange {
        match self {
            Self::Configuration { old, new } => AstStereoAtomFieldChange::Configuration {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

impl StereoBondFieldChange {
    pub(crate) fn from_rust(py: Python<'_>, change: &AstStereoBondFieldChange) -> PyResult<Self> {
        Ok(match change {
            AstStereoBondFieldChange::Configuration { old, new } => Self::Configuration {
                old: into_py_variant(py, StereoConfigurationAst::from_rust(py, old)?)?,
                new: into_py_variant(py, StereoConfigurationAst::from_rust(py, new)?)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstStereoBondFieldChange {
        match self {
            Self::Configuration { old, new } => AstStereoBondFieldChange::Configuration {
                old: old.bind(py).borrow().to_rust(py),
                new: new.bind(py).borrow().to_rust(py),
            },
        }
    }
}

pub struct AtomDeltaAstValue(Py<AtomAst>);

impl FromPyObject<'_, '_> for AtomDeltaAstValue {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, AtomAst>>()?;
        let ast = source.inner().clone();
        drop(source);
        Ok(Self(Py::new(obj.py(), AtomAst::from_inner(ast))?))
    }
}

impl<'py> IntoPyObject<'py> for &AtomDeltaAstValue {
    type Target = AtomAst;
    type Output = Bound<'py, AtomAst>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Ok(self.0.clone_ref(py).into_bound(py))
    }
}

impl AtomDeltaAstValue {
    fn from_rust(py: Python<'_>, ast: &AstAtomAst) -> PyResult<Self> {
        Ok(Self(Py::new(py, AtomAst::from_inner(ast.clone()))?))
    }

    fn to_rust(&self, py: Python<'_>) -> AstAtomAst {
        self.0.bind(py).borrow().inner().clone()
    }
}

/// A resolved edit to one atom.
#[pyclass]
pub enum AtomDelta {
    Add {
        id: u32,
        ast: AtomDeltaAstValue,
    },
    Remove {
        id: u32,
        ast: AtomDeltaAstValue,
    },
    ModifyField {
        id: u32,
        change: Py<AtomFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<AtomConstraintAst>>,
        new: Option<Py<AtomConstraintAst>>,
    },
}

#[pymethods]
impl AtomDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "ast"]),
            Self::Remove { .. } => ("Remove", &["id", "ast"]),
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
    pub(crate) fn from_rust(py: Python<'_>, delta: &AstAtomDelta) -> PyResult<Self> {
        Ok(match delta {
            AstAtomDelta::Add { id, ast } => Self::Add {
                id: id.0,
                ast: AtomDeltaAstValue::from_rust(py, ast)?,
            },
            AstAtomDelta::Remove { id, ast } => Self::Remove {
                id: id.0,
                ast: AtomDeltaAstValue::from_rust(py, ast)?,
            },
            AstAtomDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, AtomFieldChange::from_rust(py, change)?)?,
            },
            AstAtomDelta::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id: id.0,
                old: old
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, AtomConstraintAst::from_rust(py, constraint)?)
                    })
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, AtomConstraintAst::from_rust(py, constraint)?)
                    })
                    .transpose()?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstAtomDelta {
        match self {
            Self::Add { id, ast } => AstAtomDelta::Add {
                id: AstAtomId(*id),
                ast: ast.to_rust(py),
            },
            Self::Remove { id, ast } => AstAtomDelta::Remove {
                id: AstAtomId(*id),
                ast: ast.to_rust(py),
            },
            Self::ModifyField { id, change } => AstAtomDelta::ModifyField {
                id: AstAtomId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => AstAtomDelta::ModifyConstraint {
                id: AstAtomId(*id),
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

pub struct BondDeltaAstValue(Py<BondAst>);

impl FromPyObject<'_, '_> for BondDeltaAstValue {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, BondAst>>()?;
        let ast = source.inner().clone();
        drop(source);
        Ok(Self(Py::new(obj.py(), BondAst::from_inner(ast))?))
    }
}

impl<'py> IntoPyObject<'py> for &BondDeltaAstValue {
    type Target = BondAst;
    type Output = Bound<'py, BondAst>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Ok(self.0.clone_ref(py).into_bound(py))
    }
}

impl BondDeltaAstValue {
    fn from_rust(py: Python<'_>, ast: &AstBondAst) -> PyResult<Self> {
        Ok(Self(Py::new(py, BondAst::from_inner(ast.clone()))?))
    }

    fn to_rust(&self, py: Python<'_>) -> AstBondAst {
        self.0.bind(py).borrow().inner().clone()
    }
}

/// A resolved edit to one covalent bond.
#[pyclass]
pub enum BondDelta {
    Add {
        id: u32,
        atoms: (u32, u32),
        ast: BondDeltaAstValue,
    },
    Remove {
        id: u32,
        atoms: (u32, u32),
        ast: BondDeltaAstValue,
    },
    ModifyField {
        id: u32,
        change: Py<BondFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<BondConstraintAst>>,
        new: Option<Py<BondConstraintAst>>,
    },
}

#[pymethods]
impl BondDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "atoms", "ast"]),
            Self::Remove { .. } => ("Remove", &["id", "atoms", "ast"]),
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
    pub(crate) fn from_rust(py: Python<'_>, delta: &AstBondDelta) -> PyResult<Self> {
        Ok(match delta {
            AstBondDelta::Add { id, atoms, ast } => Self::Add {
                id: id.0,
                atoms: (atoms[0].0, atoms[1].0),
                ast: BondDeltaAstValue::from_rust(py, ast)?,
            },
            AstBondDelta::Remove { id, atoms, ast } => Self::Remove {
                id: id.0,
                atoms: (atoms[0].0, atoms[1].0),
                ast: BondDeltaAstValue::from_rust(py, ast)?,
            },
            AstBondDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, BondFieldChange::from_rust(py, change)?)?,
            },
            AstBondDelta::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id: id.0,
                old: old
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, BondConstraintAst::from_rust(py, constraint)?)
                    })
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, BondConstraintAst::from_rust(py, constraint)?)
                    })
                    .transpose()?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstBondDelta {
        match self {
            Self::Add { id, atoms, ast } => AstBondDelta::Add {
                id: AstBondId(*id),
                atoms: [AstAtomId(atoms.0), AstAtomId(atoms.1)],
                ast: ast.to_rust(py),
            },
            Self::Remove { id, atoms, ast } => AstBondDelta::Remove {
                id: AstBondId(*id),
                atoms: [AstAtomId(atoms.0), AstAtomId(atoms.1)],
                ast: ast.to_rust(py),
            },
            Self::ModifyField { id, change } => AstBondDelta::ModifyField {
                id: AstBondId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => AstBondDelta::ModifyConstraint {
                id: AstBondId(*id),
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

pub struct DativeBondDeltaAstValue(Py<DativeBondAst>);

impl FromPyObject<'_, '_> for DativeBondDeltaAstValue {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, DativeBondAst>>()?;
        let ast = source.inner().clone();
        drop(source);
        Ok(Self(Py::new(obj.py(), DativeBondAst::from_inner(ast))?))
    }
}

impl<'py> IntoPyObject<'py> for &DativeBondDeltaAstValue {
    type Target = DativeBondAst;
    type Output = Bound<'py, DativeBondAst>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Ok(self.0.clone_ref(py).into_bound(py))
    }
}

impl DativeBondDeltaAstValue {
    fn from_rust(py: Python<'_>, ast: &AstDativeBondAst) -> PyResult<Self> {
        Ok(Self(Py::new(py, DativeBondAst::from_inner(ast.clone()))?))
    }

    fn to_rust(&self, py: Python<'_>) -> AstDativeBondAst {
        self.0.bind(py).borrow().inner().clone()
    }
}

/// A resolved edit to one dative bond.
#[pyclass]
pub enum DativeBondDelta {
    Add {
        id: u32,
        donors: Vec<u32>,
        acceptor: u32,
        ast: DativeBondDeltaAstValue,
    },
    Remove {
        id: u32,
        donors: Vec<u32>,
        acceptor: u32,
        ast: DativeBondDeltaAstValue,
    },
    ModifyField {
        id: u32,
        change: Py<DativeBondFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<DativeBondConstraintAst>>,
        new: Option<Py<DativeBondConstraintAst>>,
    },
}

#[pymethods]
impl DativeBondDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "donors", "acceptor", "ast"]),
            Self::Remove { .. } => ("Remove", &["id", "donors", "acceptor", "ast"]),
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
    pub(crate) fn from_rust(py: Python<'_>, delta: &AstDativeBondDelta) -> PyResult<Self> {
        Ok(match delta {
            AstDativeBondDelta::Add {
                id,
                donors,
                acceptor,
                ast,
            } => Self::Add {
                id: id.0,
                donors: donors.iter().map(|atom| atom.0).collect(),
                acceptor: acceptor.0,
                ast: DativeBondDeltaAstValue::from_rust(py, ast)?,
            },
            AstDativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                ast,
            } => Self::Remove {
                id: id.0,
                donors: donors.iter().map(|atom| atom.0).collect(),
                acceptor: acceptor.0,
                ast: DativeBondDeltaAstValue::from_rust(py, ast)?,
            },
            AstDativeBondDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, DativeBondFieldChange::from_rust(py, change)?)?,
            },
            AstDativeBondDelta::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id: id.0,
                old: old
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, DativeBondConstraintAst::from_rust(py, constraint)?)
                    })
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, DativeBondConstraintAst::from_rust(py, constraint)?)
                    })
                    .transpose()?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstDativeBondDelta {
        match self {
            Self::Add {
                id,
                donors,
                acceptor,
                ast,
            } => AstDativeBondDelta::Add {
                id: AstDativeBondId(*id),
                donors: donors.iter().copied().map(AstAtomId).collect(),
                acceptor: AstAtomId(*acceptor),
                ast: ast.to_rust(py),
            },
            Self::Remove {
                id,
                donors,
                acceptor,
                ast,
            } => AstDativeBondDelta::Remove {
                id: AstDativeBondId(*id),
                donors: donors.iter().copied().map(AstAtomId).collect(),
                acceptor: AstAtomId(*acceptor),
                ast: ast.to_rust(py),
            },
            Self::ModifyField { id, change } => AstDativeBondDelta::ModifyField {
                id: AstDativeBondId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => AstDativeBondDelta::ModifyConstraint {
                id: AstDativeBondId(*id),
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

pub struct AromaticSystemDeltaAstValue(Py<AromaticSystemAst>);

impl FromPyObject<'_, '_> for AromaticSystemDeltaAstValue {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, AromaticSystemAst>>()?;
        let ast = source.inner().clone();
        drop(source);
        Ok(Self(Py::new(obj.py(), AromaticSystemAst::from_inner(ast))?))
    }
}

impl<'py> IntoPyObject<'py> for &AromaticSystemDeltaAstValue {
    type Target = AromaticSystemAst;
    type Output = Bound<'py, AromaticSystemAst>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Ok(self.0.clone_ref(py).into_bound(py))
    }
}

impl AromaticSystemDeltaAstValue {
    fn from_rust(py: Python<'_>, ast: &AstAromaticSystemAst) -> PyResult<Self> {
        Ok(Self(Py::new(
            py,
            AromaticSystemAst::from_inner(ast.clone()),
        )?))
    }

    fn to_rust(&self, py: Python<'_>) -> AstAromaticSystemAst {
        self.0.bind(py).borrow().inner().clone()
    }
}

/// A resolved edit to one aromatic system.
#[pyclass]
pub enum AromaticSystemDelta {
    Add {
        id: u32,
        atoms: Vec<u32>,
        ast: AromaticSystemDeltaAstValue,
    },
    Remove {
        id: u32,
        atoms: Vec<u32>,
        ast: AromaticSystemDeltaAstValue,
    },
    ModifyField {
        id: u32,
        change: Py<AromaticSystemFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<AromaticSystemConstraintAst>>,
        new: Option<Py<AromaticSystemConstraintAst>>,
    },
}

#[pymethods]
impl AromaticSystemDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "atoms", "ast"]),
            Self::Remove { .. } => ("Remove", &["id", "atoms", "ast"]),
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
    pub(crate) fn from_rust(py: Python<'_>, delta: &AstAromaticSystemDelta) -> PyResult<Self> {
        Ok(match delta {
            AstAromaticSystemDelta::Add { id, atoms, ast } => Self::Add {
                id: id.0,
                atoms: atoms.iter().map(|atom| atom.0).collect(),
                ast: AromaticSystemDeltaAstValue::from_rust(py, ast)?,
            },
            AstAromaticSystemDelta::Remove { id, atoms, ast } => Self::Remove {
                id: id.0,
                atoms: atoms.iter().map(|atom| atom.0).collect(),
                ast: AromaticSystemDeltaAstValue::from_rust(py, ast)?,
            },
            AstAromaticSystemDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, AromaticSystemFieldChange::from_rust(py, change)?)?,
            },
            AstAromaticSystemDelta::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id: id.0,
                old: old
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, AromaticSystemConstraintAst::from_rust(py, constraint)?)
                    })
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, AromaticSystemConstraintAst::from_rust(py, constraint)?)
                    })
                    .transpose()?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstAromaticSystemDelta {
        match self {
            Self::Add { id, atoms, ast } => AstAromaticSystemDelta::Add {
                id: AstAromaticSystemId(*id),
                atoms: atoms.iter().copied().map(AstAtomId).collect(),
                ast: ast.to_rust(py),
            },
            Self::Remove { id, atoms, ast } => AstAromaticSystemDelta::Remove {
                id: AstAromaticSystemId(*id),
                atoms: atoms.iter().copied().map(AstAtomId).collect(),
                ast: ast.to_rust(py),
            },
            Self::ModifyField { id, change } => AstAromaticSystemDelta::ModifyField {
                id: AstAromaticSystemId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => AstAromaticSystemDelta::ModifyConstraint {
                id: AstAromaticSystemId(*id),
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

pub struct MulticenterBondDeltaAstValue(Py<MulticenterBondAst>);

impl FromPyObject<'_, '_> for MulticenterBondDeltaAstValue {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, MulticenterBondAst>>()?;
        let ast = source.inner().clone();
        drop(source);
        Ok(Self(Py::new(
            obj.py(),
            MulticenterBondAst::from_inner(ast),
        )?))
    }
}

impl<'py> IntoPyObject<'py> for &MulticenterBondDeltaAstValue {
    type Target = MulticenterBondAst;
    type Output = Bound<'py, MulticenterBondAst>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Ok(self.0.clone_ref(py).into_bound(py))
    }
}

impl MulticenterBondDeltaAstValue {
    fn from_rust(py: Python<'_>, ast: &AstMulticenterBondAst) -> PyResult<Self> {
        Ok(Self(Py::new(
            py,
            MulticenterBondAst::from_inner(ast.clone()),
        )?))
    }

    fn to_rust(&self, py: Python<'_>) -> AstMulticenterBondAst {
        self.0.bind(py).borrow().inner().clone()
    }
}

/// A resolved edit to one multicenter bond.
#[pyclass]
pub enum MulticenterBondDelta {
    Add {
        id: u32,
        atoms: Vec<u32>,
        ast: MulticenterBondDeltaAstValue,
    },
    Remove {
        id: u32,
        atoms: Vec<u32>,
        ast: MulticenterBondDeltaAstValue,
    },
    ModifyField {
        id: u32,
        change: Py<MulticenterBondFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<MulticenterBondConstraintAst>>,
        new: Option<Py<MulticenterBondConstraintAst>>,
    },
}

#[pymethods]
impl MulticenterBondDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "atoms", "ast"]),
            Self::Remove { .. } => ("Remove", &["id", "atoms", "ast"]),
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
    pub(crate) fn from_rust(py: Python<'_>, delta: &AstMulticenterBondDelta) -> PyResult<Self> {
        Ok(match delta {
            AstMulticenterBondDelta::Add { id, atoms, ast } => Self::Add {
                id: id.0,
                atoms: atoms.iter().map(|atom| atom.0).collect(),
                ast: MulticenterBondDeltaAstValue::from_rust(py, ast)?,
            },
            AstMulticenterBondDelta::Remove { id, atoms, ast } => Self::Remove {
                id: id.0,
                atoms: atoms.iter().map(|atom| atom.0).collect(),
                ast: MulticenterBondDeltaAstValue::from_rust(py, ast)?,
            },
            AstMulticenterBondDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, MulticenterBondFieldChange::from_rust(py, change)?)?,
            },
            AstMulticenterBondDelta::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id: id.0,
                old: old
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(
                            py,
                            MulticenterBondConstraintAst::from_rust(py, constraint)?,
                        )
                    })
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(
                            py,
                            MulticenterBondConstraintAst::from_rust(py, constraint)?,
                        )
                    })
                    .transpose()?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstMulticenterBondDelta {
        match self {
            Self::Add { id, atoms, ast } => AstMulticenterBondDelta::Add {
                id: AstMulticenterBondId(*id),
                atoms: atoms.iter().copied().map(AstAtomId).collect(),
                ast: ast.to_rust(py),
            },
            Self::Remove { id, atoms, ast } => AstMulticenterBondDelta::Remove {
                id: AstMulticenterBondId(*id),
                atoms: atoms.iter().copied().map(AstAtomId).collect(),
                ast: ast.to_rust(py),
            },
            Self::ModifyField { id, change } => AstMulticenterBondDelta::ModifyField {
                id: AstMulticenterBondId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => AstMulticenterBondDelta::ModifyConstraint {
                id: AstMulticenterBondId(*id),
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

pub struct NoncovalentBondDeltaAstValue(Py<NoncovalentBondAst>);

impl FromPyObject<'_, '_> for NoncovalentBondDeltaAstValue {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, NoncovalentBondAst>>()?;
        let ast = source.inner().clone();
        drop(source);
        Ok(Self(Py::new(
            obj.py(),
            NoncovalentBondAst::from_inner(ast),
        )?))
    }
}

impl<'py> IntoPyObject<'py> for &NoncovalentBondDeltaAstValue {
    type Target = NoncovalentBondAst;
    type Output = Bound<'py, NoncovalentBondAst>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Ok(self.0.clone_ref(py).into_bound(py))
    }
}

impl NoncovalentBondDeltaAstValue {
    fn from_rust(py: Python<'_>, ast: &AstNoncovalentBondAst) -> PyResult<Self> {
        Ok(Self(Py::new(
            py,
            NoncovalentBondAst::from_inner(ast.clone()),
        )?))
    }

    fn to_rust(&self, py: Python<'_>) -> AstNoncovalentBondAst {
        self.0.bind(py).borrow().inner().clone()
    }
}

/// A resolved edit to one noncovalent bond.
#[pyclass]
pub enum NoncovalentBondDelta {
    Add {
        id: u32,
        atoms: (u32, u32),
        ast: NoncovalentBondDeltaAstValue,
    },
    Remove {
        id: u32,
        atoms: (u32, u32),
        ast: NoncovalentBondDeltaAstValue,
    },
    ModifyField {
        id: u32,
        change: Py<NoncovalentBondFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        old: Option<Py<NoncovalentBondConstraintAst>>,
        new: Option<Py<NoncovalentBondConstraintAst>>,
    },
}

#[pymethods]
impl NoncovalentBondDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "atoms", "ast"]),
            Self::Remove { .. } => ("Remove", &["id", "atoms", "ast"]),
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
    pub(crate) fn from_rust(py: Python<'_>, delta: &AstNoncovalentBondDelta) -> PyResult<Self> {
        Ok(match delta {
            AstNoncovalentBondDelta::Add { id, atoms, ast } => Self::Add {
                id: id.0,
                atoms: (atoms[0].0, atoms[1].0),
                ast: NoncovalentBondDeltaAstValue::from_rust(py, ast)?,
            },
            AstNoncovalentBondDelta::Remove { id, atoms, ast } => Self::Remove {
                id: id.0,
                atoms: (atoms[0].0, atoms[1].0),
                ast: NoncovalentBondDeltaAstValue::from_rust(py, ast)?,
            },
            AstNoncovalentBondDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, NoncovalentBondFieldChange::from_rust(py, change)?)?,
            },
            AstNoncovalentBondDelta::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id: id.0,
                old: old
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(
                            py,
                            NoncovalentBondConstraintAst::from_rust(py, constraint)?,
                        )
                    })
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(
                            py,
                            NoncovalentBondConstraintAst::from_rust(py, constraint)?,
                        )
                    })
                    .transpose()?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstNoncovalentBondDelta {
        match self {
            Self::Add { id, atoms, ast } => AstNoncovalentBondDelta::Add {
                id: AstNoncovalentBondId(*id),
                atoms: [AstAtomId(atoms.0), AstAtomId(atoms.1)],
                ast: ast.to_rust(py),
            },
            Self::Remove { id, atoms, ast } => AstNoncovalentBondDelta::Remove {
                id: AstNoncovalentBondId(*id),
                atoms: [AstAtomId(atoms.0), AstAtomId(atoms.1)],
                ast: ast.to_rust(py),
            },
            Self::ModifyField { id, change } => AstNoncovalentBondDelta::ModifyField {
                id: AstNoncovalentBondId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, old, new } => AstNoncovalentBondDelta::ModifyConstraint {
                id: AstNoncovalentBondId(*id),
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

pub struct StereoAtomDeltaAstValue(Py<StereoAtomAst>);

impl FromPyObject<'_, '_> for StereoAtomDeltaAstValue {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, StereoAtomAst>>()?;
        let ast = source.inner().clone();
        drop(source);
        Ok(Self(Py::new(obj.py(), StereoAtomAst::from_inner(ast))?))
    }
}

impl<'py> IntoPyObject<'py> for &StereoAtomDeltaAstValue {
    type Target = StereoAtomAst;
    type Output = Bound<'py, StereoAtomAst>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Ok(self.0.clone_ref(py).into_bound(py))
    }
}

impl StereoAtomDeltaAstValue {
    fn from_rust(py: Python<'_>, ast: &AstStereoAtomAst) -> PyResult<Self> {
        Ok(Self(Py::new(py, StereoAtomAst::from_inner(ast.clone()))?))
    }

    fn to_rust(&self, py: Python<'_>) -> AstStereoAtomAst {
        self.0.bind(py).borrow().inner().clone()
    }
}

/// A resolved edit to one atom-centered stereo element.
#[pyclass]
pub enum StereoAtomDelta {
    Add {
        id: u32,
        site: u32,
        ligands: Vec<StereoLigand>,
        ast: StereoAtomDeltaAstValue,
    },
    Remove {
        id: u32,
        site: u32,
        ligands: Vec<StereoLigand>,
        ast: StereoAtomDeltaAstValue,
    },
    ModifyField {
        id: u32,
        change: Py<StereoAtomFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        kind: Option<StereoKind>,
        old: Option<Py<StereoAtomConstraintAst>>,
        new: Option<Py<StereoAtomConstraintAst>>,
    },
    Apply {
        id: u32,
        kind: StereoKind,
        permutation: Permutation,
    },
    Swap {
        id: u32,
        kind: StereoKind,
    },
    Mirror {
        id: u32,
        kind: StereoKind,
    },
}

#[pymethods]
impl StereoAtomDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "site", "ligands", "ast"]),
            Self::Remove { .. } => ("Remove", &["id", "site", "ligands", "ast"]),
            Self::ModifyField { .. } => ("ModifyField", &["id", "change"]),
            Self::ModifyConstraint { .. } => ("ModifyConstraint", &["id", "kind", "old", "new"]),
            Self::Apply { .. } => ("Apply", &["id", "kind", "permutation"]),
            Self::Swap { .. } => ("Swap", &["id", "kind"]),
            Self::Mirror { .. } => ("Mirror", &["id", "kind"]),
        };
        entity_delta_repr(slf.bind(py).as_any(), "StereoAtomDelta", variant, fields)
    }

    /// Return the inverse resolved edit.
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl StereoAtomDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &AstStereoAtomDelta) -> PyResult<Self> {
        Ok(match delta {
            AstStereoAtomDelta::Add {
                id,
                site,
                ligands,
                ast,
            } => Self::Add {
                id: id.0,
                site: site.0,
                ligands: ligands
                    .iter()
                    .copied()
                    .map(StereoLigand::from_rust)
                    .collect(),
                ast: StereoAtomDeltaAstValue::from_rust(py, ast)?,
            },
            AstStereoAtomDelta::Remove {
                id,
                site,
                ligands,
                ast,
            } => Self::Remove {
                id: id.0,
                site: site.0,
                ligands: ligands
                    .iter()
                    .copied()
                    .map(StereoLigand::from_rust)
                    .collect(),
                ast: StereoAtomDeltaAstValue::from_rust(py, ast)?,
            },
            AstStereoAtomDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, StereoAtomFieldChange::from_rust(py, change)?)?,
            },
            AstStereoAtomDelta::ModifyConstraint { id, kind, old, new } => Self::ModifyConstraint {
                id: id.0,
                kind: kind.map(StereoKind::from_rust),
                old: old
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, StereoAtomConstraintAst::from_rust(py, constraint)?)
                    })
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|constraint| {
                        into_py_variant(py, StereoAtomConstraintAst::from_rust(py, constraint)?)
                    })
                    .transpose()?,
            },
            AstStereoAtomDelta::Apply {
                id,
                kind,
                permutation,
            } => Self::Apply {
                id: id.0,
                kind: StereoKind::from_rust(*kind),
                permutation: Permutation::from_inner(*permutation),
            },
            AstStereoAtomDelta::Swap { id, kind } => Self::Swap {
                id: id.0,
                kind: StereoKind::from_rust(*kind),
            },
            AstStereoAtomDelta::Mirror { id, kind } => Self::Mirror {
                id: id.0,
                kind: StereoKind::from_rust(*kind),
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstStereoAtomDelta {
        match self {
            Self::Add {
                id,
                site,
                ligands,
                ast,
            } => AstStereoAtomDelta::Add {
                id: AstStereoAtomId(*id),
                site: AstAtomId(*site),
                ligands: ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                ast: ast.to_rust(py),
            },
            Self::Remove {
                id,
                site,
                ligands,
                ast,
            } => AstStereoAtomDelta::Remove {
                id: AstStereoAtomId(*id),
                site: AstAtomId(*site),
                ligands: ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                ast: ast.to_rust(py),
            },
            Self::ModifyField { id, change } => AstStereoAtomDelta::ModifyField {
                id: AstStereoAtomId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, kind, old, new } => AstStereoAtomDelta::ModifyConstraint {
                id: AstStereoAtomId(*id),
                kind: kind.map(StereoKind::to_rust),
                old: old
                    .as_ref()
                    .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
                new: new
                    .as_ref()
                    .map(|constraint| constraint.bind(py).borrow().to_rust(py)),
            },
            Self::Apply {
                id,
                kind,
                permutation,
            } => AstStereoAtomDelta::Apply {
                id: AstStereoAtomId(*id),
                kind: kind.to_rust(),
                permutation: permutation.inner(),
            },
            Self::Swap { id, kind } => AstStereoAtomDelta::Swap {
                id: AstStereoAtomId(*id),
                kind: kind.to_rust(),
            },
            Self::Mirror { id, kind } => AstStereoAtomDelta::Mirror {
                id: AstStereoAtomId(*id),
                kind: kind.to_rust(),
            },
        }
    }
}

pub struct StereoBondDeltaAstValue(Py<StereoBondAst>);

impl FromPyObject<'_, '_> for StereoBondDeltaAstValue {
    type Error = PyErr;
    fn extract(obj: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, StereoBondAst>>()?;
        let ast = source.inner().clone();
        drop(source);
        Ok(Self(Py::new(obj.py(), StereoBondAst::from_inner(ast))?))
    }
}

impl<'py> IntoPyObject<'py> for &StereoBondDeltaAstValue {
    type Target = StereoBondAst;
    type Output = Bound<'py, StereoBondAst>;
    type Error = PyErr;
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Ok(self.0.clone_ref(py).into_bound(py))
    }
}

impl StereoBondDeltaAstValue {
    fn from_rust(py: Python<'_>, ast: &AstStereoBondAst) -> PyResult<Self> {
        Ok(Self(Py::new(py, StereoBondAst::from_inner(ast.clone()))?))
    }
    fn to_rust(&self, py: Python<'_>) -> AstStereoBondAst {
        self.0.bind(py).borrow().inner().clone()
    }
}

/// A resolved edit to one bond-centered stereo element.
#[pyclass]
pub enum StereoBondDelta {
    Add {
        id: u32,
        site: u32,
        ligands: Vec<StereoLigand>,
        ast: StereoBondDeltaAstValue,
    },
    Remove {
        id: u32,
        site: u32,
        ligands: Vec<StereoLigand>,
        ast: StereoBondDeltaAstValue,
    },
    ModifyField {
        id: u32,
        change: Py<StereoBondFieldChange>,
    },
    ModifyConstraint {
        id: u32,
        kind: Option<StereoKind>,
        old: Option<Py<StereoBondConstraintAst>>,
        new: Option<Py<StereoBondConstraintAst>>,
    },
    Apply {
        id: u32,
        kind: StereoKind,
        permutation: Permutation,
    },
    Swap {
        id: u32,
        kind: StereoKind,
    },
    Mirror {
        id: u32,
        kind: StereoKind,
    },
}

#[pymethods]
impl StereoBondDelta {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }
    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, fields): (&str, &[&str]) = match &*slf.bind(py).borrow() {
            Self::Add { .. } => ("Add", &["id", "site", "ligands", "ast"]),
            Self::Remove { .. } => ("Remove", &["id", "site", "ligands", "ast"]),
            Self::ModifyField { .. } => ("ModifyField", &["id", "change"]),
            Self::ModifyConstraint { .. } => ("ModifyConstraint", &["id", "kind", "old", "new"]),
            Self::Apply { .. } => ("Apply", &["id", "kind", "permutation"]),
            Self::Swap { .. } => ("Swap", &["id", "kind"]),
            Self::Mirror { .. } => ("Mirror", &["id", "kind"]),
        };
        entity_delta_repr(slf.bind(py).as_any(), "StereoBondDelta", variant, fields)
    }
    fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        into_py_variant(py, Self::from_rust(py, &self.to_rust(py).inverse())?)
    }
}

impl StereoBondDelta {
    pub(crate) fn from_rust(py: Python<'_>, delta: &AstStereoBondDelta) -> PyResult<Self> {
        Ok(match delta {
            AstStereoBondDelta::Add {
                id,
                site,
                ligands,
                ast,
            } => Self::Add {
                id: id.0,
                site: site.0,
                ligands: ligands
                    .iter()
                    .copied()
                    .map(StereoLigand::from_rust)
                    .collect(),
                ast: StereoBondDeltaAstValue::from_rust(py, ast)?,
            },
            AstStereoBondDelta::Remove {
                id,
                site,
                ligands,
                ast,
            } => Self::Remove {
                id: id.0,
                site: site.0,
                ligands: ligands
                    .iter()
                    .copied()
                    .map(StereoLigand::from_rust)
                    .collect(),
                ast: StereoBondDeltaAstValue::from_rust(py, ast)?,
            },
            AstStereoBondDelta::ModifyField { id, change } => Self::ModifyField {
                id: id.0,
                change: into_py_variant(py, StereoBondFieldChange::from_rust(py, change)?)?,
            },
            AstStereoBondDelta::ModifyConstraint { id, kind, old, new } => Self::ModifyConstraint {
                id: id.0,
                kind: kind.map(StereoKind::from_rust),
                old: old
                    .as_ref()
                    .map(|c| into_py_variant(py, StereoBondConstraintAst::from_rust(py, c)?))
                    .transpose()?,
                new: new
                    .as_ref()
                    .map(|c| into_py_variant(py, StereoBondConstraintAst::from_rust(py, c)?))
                    .transpose()?,
            },
            AstStereoBondDelta::Apply {
                id,
                kind,
                permutation,
            } => Self::Apply {
                id: id.0,
                kind: StereoKind::from_rust(*kind),
                permutation: Permutation::from_inner(*permutation),
            },
            AstStereoBondDelta::Swap { id, kind } => Self::Swap {
                id: id.0,
                kind: StereoKind::from_rust(*kind),
            },
            AstStereoBondDelta::Mirror { id, kind } => Self::Mirror {
                id: id.0,
                kind: StereoKind::from_rust(*kind),
            },
        })
    }
    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstStereoBondDelta {
        match self {
            Self::Add {
                id,
                site,
                ligands,
                ast,
            } => AstStereoBondDelta::Add {
                id: AstStereoBondId(*id),
                site: AstBondId(*site),
                ligands: ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                ast: ast.to_rust(py),
            },
            Self::Remove {
                id,
                site,
                ligands,
                ast,
            } => AstStereoBondDelta::Remove {
                id: AstStereoBondId(*id),
                site: AstBondId(*site),
                ligands: ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                ast: ast.to_rust(py),
            },
            Self::ModifyField { id, change } => AstStereoBondDelta::ModifyField {
                id: AstStereoBondId(*id),
                change: change.bind(py).borrow().to_rust(py),
            },
            Self::ModifyConstraint { id, kind, old, new } => AstStereoBondDelta::ModifyConstraint {
                id: AstStereoBondId(*id),
                kind: kind.map(StereoKind::to_rust),
                old: old.as_ref().map(|c| c.bind(py).borrow().to_rust(py)),
                new: new.as_ref().map(|c| c.bind(py).borrow().to_rust(py)),
            },
            Self::Apply {
                id,
                kind,
                permutation,
            } => AstStereoBondDelta::Apply {
                id: AstStereoBondId(*id),
                kind: kind.to_rust(),
                permutation: permutation.inner(),
            },
            Self::Swap { id, kind } => AstStereoBondDelta::Swap {
                id: AstStereoBondId(*id),
                kind: kind.to_rust(),
            },
            Self::Mirror { id, kind } => AstStereoBondDelta::Mirror {
                id: AstStereoBondId(*id),
                kind: kind.to_rust(),
            },
        }
    }
}

pub struct ConstraintDeltaValue(Py<Constraint>);

impl FromPyObject<'_, '_> for ConstraintDeltaValue {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, '_, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, Constraint>>()?;
        let constraint = source.to_rust(obj.py());
        drop(source);
        Ok(Self(into_py_variant(
            obj.py(),
            Constraint::from_rust(obj.py(), &constraint)?,
        )?))
    }
}

impl<'py> IntoPyObject<'py> for &ConstraintDeltaValue {
    type Target = Constraint;
    type Output = Bound<'py, Constraint>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Ok(self.0.clone_ref(py).into_bound(py))
    }
}

impl ConstraintDeltaValue {
    fn from_rust(py: Python<'_>, constraint: &AstConstraint) -> PyResult<Self> {
        Ok(Self(into_py_variant(
            py,
            Constraint::from_rust(py, constraint)?,
        )?))
    }

    fn to_rust(&self, py: Python<'_>) -> AstConstraint {
        self.0.bind(py).borrow().to_rust(py)
    }
}

/// A resolved edit adding or removing a molecule constraint.
#[pyclass]
pub enum ConstraintDelta {
    Add { constraint: ConstraintDeltaValue },
    Remove { constraint: ConstraintDeltaValue },
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
    pub(crate) fn from_rust(py: Python<'_>, delta: &AstConstraintDelta) -> PyResult<Self> {
        Ok(match delta {
            AstConstraintDelta::Add(constraint) => Self::Add {
                constraint: ConstraintDeltaValue::from_rust(py, constraint)?,
            },
            AstConstraintDelta::Remove(constraint) => Self::Remove {
                constraint: ConstraintDeltaValue::from_rust(py, constraint)?,
            },
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstConstraintDelta {
        match self {
            Self::Add { constraint } => AstConstraintDelta::Add(constraint.to_rust(py)),
            Self::Remove { constraint } => AstConstraintDelta::Remove(constraint.to_rust(py)),
        }
    }
}

/// One resolved edit from any localized-topology family.
#[pyclass]
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
    pub(crate) fn from_rust(py: Python<'_>, delta: &AstDelta) -> PyResult<Self> {
        Ok(match delta {
            AstDelta::Atom(delta) => {
                Self::Atom(into_py_variant(py, AtomDelta::from_rust(py, delta)?)?)
            }
            AstDelta::Bond(delta) => {
                Self::Bond(into_py_variant(py, BondDelta::from_rust(py, delta)?)?)
            }
            AstDelta::DativeBond(delta) => {
                Self::DativeBond(into_py_variant(py, DativeBondDelta::from_rust(py, delta)?)?)
            }
            AstDelta::AromaticSystem(delta) => Self::AromaticSystem(into_py_variant(
                py,
                AromaticSystemDelta::from_rust(py, delta)?,
            )?),
            AstDelta::MulticenterBond(delta) => Self::MulticenterBond(into_py_variant(
                py,
                MulticenterBondDelta::from_rust(py, delta)?,
            )?),
            AstDelta::NoncovalentBond(delta) => Self::NoncovalentBond(into_py_variant(
                py,
                NoncovalentBondDelta::from_rust(py, delta)?,
            )?),
            AstDelta::StereoAtom(delta) => {
                Self::StereoAtom(into_py_variant(py, StereoAtomDelta::from_rust(py, delta)?)?)
            }
            AstDelta::StereoBond(delta) => {
                Self::StereoBond(into_py_variant(py, StereoBondDelta::from_rust(py, delta)?)?)
            }
            AstDelta::Constraint(delta) => {
                Self::Constraint(into_py_variant(py, ConstraintDelta::from_rust(py, delta)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstDelta {
        match self {
            Self::Atom(delta) => AstDelta::Atom(delta.bind(py).borrow().to_rust(py)),
            Self::Bond(delta) => AstDelta::Bond(delta.bind(py).borrow().to_rust(py)),
            Self::DativeBond(delta) => AstDelta::DativeBond(delta.bind(py).borrow().to_rust(py)),
            Self::AromaticSystem(delta) => {
                AstDelta::AromaticSystem(delta.bind(py).borrow().to_rust(py))
            }
            Self::MulticenterBond(delta) => {
                AstDelta::MulticenterBond(delta.bind(py).borrow().to_rust(py))
            }
            Self::NoncovalentBond(delta) => {
                AstDelta::NoncovalentBond(delta.bind(py).borrow().to_rust(py))
            }
            Self::StereoAtom(delta) => AstDelta::StereoAtom(delta.bind(py).borrow().to_rust(py)),
            Self::StereoBond(delta) => AstDelta::StereoBond(delta.bind(py).borrow().to_rust(py)),
            Self::Constraint(delta) => AstDelta::Constraint(delta.bind(py).borrow().to_rust(py)),
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
fn delta_iter(py: Python<'_>, deltas: &AstDeltas) -> PyResult<DeltaIter> {
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
            Self::Container(container) => container
                .bind(py)
                .borrow()
                .to_rust()
                .as_slice()
                .to_vec(),
            Self::Entries(entries) => entries
                .iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py))
                .collect(),
        };
        ResolvedDeltasExtend(entries)
    }
}

/// An extend input containing no Python references that need to be read.
pub(crate) struct ResolvedDeltasExtend(Vec<AstDelta>);

impl ResolvedDeltasExtend {
    /// Append the resolved deltas in order, preserving duplicates.
    fn apply(self, target: &mut AstDeltas) {
        for delta in self.0 {
            target.push(delta);
        }
    }
}

/// Resolved deltas in insertion order. Mutable, value-equal, and unhashable.
#[pyclass(eq)]
#[derive(Debug, PartialEq)]
pub struct Deltas(AstDeltas);

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
        resolved.apply(slf.borrow_mut(py).inner_mut());
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

impl Deltas {
    pub(crate) fn from_rust(deltas: AstDeltas) -> Self {
        Self(deltas)
    }

    #[allow(
        dead_code,
        reason = "component snapshot conversion for owning wrappers"
    )]
    pub(crate) fn to_rust(&self) -> AstDeltas {
        self.0.clone()
    }

    fn inner_mut(&mut self) -> &mut AstDeltas {
        &mut self.0
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AromaticSystemConstraintAst as AstAromaticSystemConstraintAst,
        AtomConstraintAst as AstAtomConstraintAst, BondConstraintAst as AstBondConstraintAst,
        BooleanAst as AstBooleanAst, DativeBondConstraintAst as AstDativeBondConstraintAst,
        ElectronCountsAst as AstElectronCountsAst, ElementAst as AstElementAst,
        IsotopeMassAst as AstIsotopeMassAst,
        MulticenterBondConstraintAst as AstMulticenterBondConstraintAst,
        NoncovalentBondConstraintAst as AstNoncovalentBondConstraintAst,
        NoncovalentBondKind as AstNoncovalentBondKind,
        NoncovalentBondKindAst as AstNoncovalentBondKindAst, SpinStateAst as AstSpinStateAst,
        StereoAtomConstraintAst as AstStereoAtomConstraintAst,
        StereoBondConstraintAst as AstStereoBondConstraintAst,
        StereoConfigurationAst as AstStereoConfigurationAst, StereoCosetAst as AstStereoCosetAst,
        StereoKind as AstStereoKind, StereoLigand as AstStereoLigand,
        StereoLigandKind as AstStereoLigandKind, Stereogenicity as AstStereogenicity,
        StereogenicityAst as AstStereogenicityAst, ValueAst as AstValueAst,
    };
    use umol_chem::element::Element as ChemElement;
    use umol_perm::Permutation as PermPermutation;

    use super::*;

    #[rstest]
    #[case::element(AstAtomFieldChange::Element {
        old: AstElementAst::Lit(ChemElement::C),
        new: AstElementAst::Lit(ChemElement::N),
    })]
    #[case::isotope_mass(AstAtomFieldChange::IsotopeMass {
        old: AstIsotopeMassAst::Lit(12),
        new: AstIsotopeMassAst::Lit(13),
    })]
    #[case::charge(AstAtomFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(-1),
    })]
    #[case::implicit_hydrogens(AstAtomFieldChange::ImplicitHydrogens {
        old: AstValueAst::Lit(3),
        new: AstValueAst::Lit(2),
    })]
    #[case::lone_pairs(AstAtomFieldChange::LonePairs {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    #[case::spin(AstAtomFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_atom_field_change_roundtrip(#[case] change: AstAtomFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                AtomFieldChange::from_rust(py, &change).unwrap().to_rust(py),
                change
            );
        });
    }

    #[rstest]
    #[case::equal(
        AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
        AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
        true,
    )]
    #[case::different(
        AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
        AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(1),
        },
        false,
    )]
    fn test_atom_field_change_eq(
        #[case] lhs: AstAtomFieldChange,
        #[case] rhs: AstAtomFieldChange,
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
        AstAtomFieldChange::Element {
            old: AstElementAst::Lit(ChemElement::C),
            new: AstElementAst::Lit(ChemElement::N),
        },
        "ElementAst.Lit(Element('C'))",
        "ElementAst.Lit(Element('N'))",
        "AtomFieldChange.Element(old=ElementAst.Lit(Element('C')), new=ElementAst.Lit(Element('N')))"
    )]
    #[case::isotope_mass(
        AstAtomFieldChange::IsotopeMass {
            old: AstIsotopeMassAst::Lit(12),
            new: AstIsotopeMassAst::Lit(13),
        },
        "IsotopeMassAst.Lit(12)",
        "IsotopeMassAst.Lit(13)",
        "AtomFieldChange.IsotopeMass(old=IsotopeMassAst.Lit(12), new=IsotopeMassAst.Lit(13))"
    )]
    #[case::charge(
        AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
        "ValueAst.Lit(0)",
        "ValueAst.Lit(-1)",
        "AtomFieldChange.Charge(old=ValueAst.Lit(0), new=ValueAst.Lit(-1))"
    )]
    #[case::implicit_hydrogens(
        AstAtomFieldChange::ImplicitHydrogens {
            old: AstValueAst::Lit(3),
            new: AstValueAst::Lit(2),
        },
        "ValueAst.Lit(3)",
        "ValueAst.Lit(2)",
        "AtomFieldChange.ImplicitHydrogens(old=ValueAst.Lit(3), new=ValueAst.Lit(2))"
    )]
    #[case::lone_pairs(
        AstAtomFieldChange::LonePairs {
            old: AstValueAst::Lit(1),
            new: AstValueAst::Lit(2),
        },
        "ValueAst.Lit(1)",
        "ValueAst.Lit(2)",
        "AtomFieldChange.LonePairs(old=ValueAst.Lit(1), new=ValueAst.Lit(2))"
    )]
    #[case::spin(
        AstAtomFieldChange::Spin {
            old: AstSpinStateAst {
                unpaired: AstValueAst::Lit(0),
                multiplicity: AstValueAst::Lit(1),
            },
            new: AstSpinStateAst {
                unpaired: AstValueAst::Lit(1),
                multiplicity: AstValueAst::Lit(2),
            },
        },
        "SpinStateAst(ValueAst.Lit(0), ValueAst.Lit(1))",
        "SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2))",
        "AtomFieldChange.Spin(old=SpinStateAst(ValueAst.Lit(0), ValueAst.Lit(1)), new=SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2)))"
    )]
    fn test_atom_field_change_repr(
        #[case] change: AstAtomFieldChange,
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
    #[case::element(AstAtomFieldChange::Element {
        old: AstElementAst::Lit(ChemElement::C),
        new: AstElementAst::Lit(ChemElement::N),
    })]
    #[case::isotope_mass(AstAtomFieldChange::IsotopeMass {
        old: AstIsotopeMassAst::Lit(12),
        new: AstIsotopeMassAst::Lit(13),
    })]
    #[case::charge(AstAtomFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(-1),
    })]
    #[case::implicit_hydrogens(AstAtomFieldChange::ImplicitHydrogens {
        old: AstValueAst::Lit(3),
        new: AstValueAst::Lit(2),
    })]
    #[case::lone_pairs(AstAtomFieldChange::LonePairs {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    #[case::spin(AstAtomFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_atom_field_change_inverse(#[case] change: AstAtomFieldChange) {
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
    #[case::order(AstBondFieldChange::Order {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    #[case::charge(AstBondFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(1),
    })]
    #[case::spin(AstBondFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_bond_field_change_roundtrip(#[case] change: AstBondFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                BondFieldChange::from_rust(py, &change).unwrap().to_rust(py),
                change
            );
        });
    }

    #[rstest]
    #[case::order(AstBondFieldChange::Order {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    #[case::charge(AstBondFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(1),
    })]
    #[case::spin(AstBondFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_bond_field_change_inverse(#[case] change: AstBondFieldChange) {
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
    #[case::order(AstDativeBondFieldChange::Order {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    fn test_dative_bond_field_change_roundtrip(#[case] change: AstDativeBondFieldChange) {
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
    #[case::order(AstDativeBondFieldChange::Order {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    fn test_dative_bond_field_change_inverse(#[case] change: AstDativeBondFieldChange) {
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
    #[case::electrons(AstAromaticSystemFieldChange::Electrons {
        old: AstElectronCountsAst::Undetermined,
        new: AstElectronCountsAst::Lit(vec![1, 1, 1]),
    })]
    #[case::charge(AstAromaticSystemFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(-1),
    })]
    #[case::spin(AstAromaticSystemFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_aromatic_system_field_change_roundtrip(#[case] change: AstAromaticSystemFieldChange) {
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
    #[case::electrons(AstAromaticSystemFieldChange::Electrons {
        old: AstElectronCountsAst::Undetermined,
        new: AstElectronCountsAst::Lit(vec![1, 1, 1]),
    })]
    #[case::charge(AstAromaticSystemFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(-1),
    })]
    #[case::spin(AstAromaticSystemFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_aromatic_system_field_change_inverse(#[case] change: AstAromaticSystemFieldChange) {
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
    #[case::electrons(AstMulticenterBondFieldChange::Electrons {
        old: AstElectronCountsAst::Lit(vec![1, 0, 1]),
        new: AstElectronCountsAst::Lit(vec![2, 0, 1]),
    })]
    #[case::charge(AstMulticenterBondFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(1),
    })]
    #[case::spin(AstMulticenterBondFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(2),
            multiplicity: AstValueAst::Lit(3),
        },
    })]
    fn test_multicenter_bond_field_change_roundtrip(#[case] change: AstMulticenterBondFieldChange) {
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
    #[case::electrons(AstMulticenterBondFieldChange::Electrons {
        old: AstElectronCountsAst::Lit(vec![1, 0, 1]),
        new: AstElectronCountsAst::Lit(vec![2, 0, 1]),
    })]
    #[case::charge(AstMulticenterBondFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(1),
    })]
    #[case::spin(AstMulticenterBondFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(2),
            multiplicity: AstValueAst::Lit(3),
        },
    })]
    fn test_multicenter_bond_field_change_inverse(#[case] change: AstMulticenterBondFieldChange) {
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
    #[case::kind(AstNoncovalentBondFieldChange::Kind {
        old: AstNoncovalentBondKindAst::Undetermined,
        new: AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond),
    })]
    fn test_noncovalent_bond_field_change_roundtrip(#[case] change: AstNoncovalentBondFieldChange) {
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
    #[case::kind(AstNoncovalentBondFieldChange::Kind {
        old: AstNoncovalentBondKindAst::Undetermined,
        new: AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond),
    })]
    fn test_noncovalent_bond_field_change_inverse(#[case] change: AstNoncovalentBondFieldChange) {
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
    #[case::geometry_unknown(AstStereoAtomFieldChange::Configuration {
        old: AstStereoConfigurationAst::Undetermined,
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Undetermined,
        ),
    })]
    #[case::coset_resolved(AstStereoAtomFieldChange::Configuration {
        old: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Undetermined,
        ),
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Lit(1),
        ),
    })]
    fn test_stereo_atom_field_change_roundtrip(#[case] change: AstStereoAtomFieldChange) {
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
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Undetermined,
            ),
        },
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Undetermined,
            ),
        },
        true,
    )]
    #[case::different(
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Undetermined,
            ),
        },
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(1),
            ),
        },
        false,
    )]
    fn test_stereo_atom_field_change_eq(
        #[case] lhs: AstStereoAtomFieldChange,
        #[case] rhs: AstStereoAtomFieldChange,
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
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Undetermined,
            ),
        },
        "StereoConfigurationAst.Undetermined()",
        "StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Undetermined())",
        "StereoAtomFieldChange.Configuration(old=StereoConfigurationAst.Undetermined(), new=StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Undetermined()))",
    )]
    #[case::coset_resolved(
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Undetermined,
            ),
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(1),
            ),
        },
        "StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Undetermined())",
        "StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Lit(1))",
        "StereoAtomFieldChange.Configuration(old=StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Undetermined()), new=StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Lit(1)))",
    )]
    fn test_stereo_atom_field_change_repr(
        #[case] change: AstStereoAtomFieldChange,
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
    #[case::geometry_unknown(AstStereoAtomFieldChange::Configuration {
        old: AstStereoConfigurationAst::Undetermined,
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Undetermined,
        ),
    })]
    #[case::coset_resolved(AstStereoAtomFieldChange::Configuration {
        old: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Undetermined,
        ),
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Lit(1),
        ),
    })]
    fn test_stereo_atom_field_change_inverse(#[case] change: AstStereoAtomFieldChange) {
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
    #[case::geometry_unknown(AstStereoBondFieldChange::Configuration {
        old: AstStereoConfigurationAst::Undetermined,
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Undetermined,
        ),
    })]
    #[case::coset_resolved(AstStereoBondFieldChange::Configuration {
        old: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Undetermined,
        ),
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Lit(1),
        ),
    })]
    fn test_stereo_bond_field_change_roundtrip(#[case] change: AstStereoBondFieldChange) {
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
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Undetermined,
            ),
        },
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Undetermined,
            ),
        },
        true,
    )]
    #[case::different(
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Undetermined,
            ),
        },
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Lit(1),
            ),
        },
        false,
    )]
    fn test_stereo_bond_field_change_eq(
        #[case] lhs: AstStereoBondFieldChange,
        #[case] rhs: AstStereoBondFieldChange,
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
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Undetermined,
            ),
        },
        "StereoConfigurationAst.Undetermined()",
        "StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Undetermined())",
        "StereoBondFieldChange.Configuration(old=StereoConfigurationAst.Undetermined(), new=StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Undetermined()))",
    )]
    #[case::coset_resolved(
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Undetermined,
            ),
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Lit(1),
            ),
        },
        "StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Undetermined())",
        "StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Lit(1))",
        "StereoBondFieldChange.Configuration(old=StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Undetermined()), new=StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Lit(1)))",
    )]
    fn test_stereo_bond_field_change_repr(
        #[case] change: AstStereoBondFieldChange,
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
    #[case::geometry_unknown(AstStereoBondFieldChange::Configuration {
        old: AstStereoConfigurationAst::Undetermined,
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Undetermined,
        ),
    })]
    #[case::coset_resolved(AstStereoBondFieldChange::Configuration {
        old: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Undetermined,
        ),
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Lit(1),
        ),
    })]
    fn test_stereo_bond_field_change_inverse(#[case] change: AstStereoBondFieldChange) {
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
    #[case::add(AstAtomDelta::Add {
        id: AstAtomId(3),
        ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
    })]
    #[case::remove(AstAtomDelta::Remove {
        id: AstAtomId(3),
        ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
    })]
    #[case::modify_field(AstAtomDelta::ModifyField {
        id: AstAtomId(3),
        change: AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
    })]
    #[case::constraint_added(AstAtomDelta::ModifyConstraint {
        id: AstAtomId(3),
        old: None,
        new: Some(AstAtomConstraintAst::Valence(AstValueAst::Lit(4))),
    })]
    #[case::constraint_removed(AstAtomDelta::ModifyConstraint {
        id: AstAtomId(3),
        old: Some(AstAtomConstraintAst::Valence(AstValueAst::Lit(4))),
        new: None,
    })]
    #[case::constraint_modified(AstAtomDelta::ModifyConstraint {
        id: AstAtomId(3),
        old: Some(AstAtomConstraintAst::Valence(AstValueAst::Lit(3))),
        new: Some(AstAtomConstraintAst::Valence(AstValueAst::Lit(4))),
    })]
    fn test_atom_delta_roundtrip(#[case] delta: AstAtomDelta) {
        Python::attach(|py| {
            assert_eq!(AtomDelta::from_rust(py, &delta).unwrap().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::equal(
        AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        },
        AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        },
        true,
    )]
    #[case::different(
        AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        },
        AstAtomDelta::Add {
            id: AstAtomId(4),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        },
        false,
    )]
    fn test_atom_delta_eq(
        #[case] lhs: AstAtomDelta,
        #[case] rhs: AstAtomDelta,
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
        AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        },
        "AtomDelta.Add(id=3, ast=AtomAst.parse('C'))",
    )]
    #[case::remove(
        AstAtomDelta::Remove {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        },
        "AtomDelta.Remove(id=3, ast=AtomAst.parse('C'))",
    )]
    #[case::modify_field(
        AstAtomDelta::ModifyField {
            id: AstAtomId(3),
            change: AstAtomFieldChange::Charge {
                old: AstValueAst::Lit(0),
                new: AstValueAst::Lit(-1),
            },
        },
        "AtomDelta.ModifyField(id=3, change=AtomFieldChange.Charge(old=ValueAst.Lit(0), new=ValueAst.Lit(-1)))",
    )]
    #[case::modify_constraint(
        AstAtomDelta::ModifyConstraint {
            id: AstAtomId(3),
            old: None,
            new: Some(AstAtomConstraintAst::Valence(AstValueAst::Lit(4))),
        },
        "AtomDelta.ModifyConstraint(id=3, old=None, new=AtomConstraintAst.Valence(ValueAst.Lit(4)))",
    )]
    fn test_atom_delta_repr(#[case] delta: AstAtomDelta, #[case] expected: &str) {
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
    #[case::add(AstAtomDelta::Add {
        id: AstAtomId(3),
        ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
    })]
    #[case::remove(AstAtomDelta::Remove {
        id: AstAtomId(3),
        ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
    })]
    #[case::modify_field(AstAtomDelta::ModifyField {
        id: AstAtomId(3),
        change: AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
    })]
    #[case::constraint_added(AstAtomDelta::ModifyConstraint {
        id: AstAtomId(3),
        old: None,
        new: Some(AstAtomConstraintAst::Valence(AstValueAst::Lit(4))),
    })]
    #[case::constraint_removed(AstAtomDelta::ModifyConstraint {
        id: AstAtomId(3),
        old: Some(AstAtomConstraintAst::Valence(AstValueAst::Lit(4))),
        new: None,
    })]
    #[case::constraint_modified(AstAtomDelta::ModifyConstraint {
        id: AstAtomId(3),
        old: Some(AstAtomConstraintAst::Valence(AstValueAst::Lit(3))),
        new: Some(AstAtomConstraintAst::Valence(AstValueAst::Lit(4))),
    })]
    fn test_atom_delta_inverse(#[case] delta: AstAtomDelta) {
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
    #[case::add(AstBondDelta::Add {
        id: AstBondId(2),
        atoms: [AstAtomId(5), AstAtomId(1)],
        ast: AstBondAst::new(AstValueAst::Lit(1)),
    })]
    #[case::remove(AstBondDelta::Remove {
        id: AstBondId(2),
        atoms: [AstAtomId(5), AstAtomId(1)],
        ast: AstBondAst::new(AstValueAst::Lit(1)),
    })]
    #[case::modify_field(AstBondDelta::ModifyField {
        id: AstBondId(2),
        change: AstBondFieldChange::Order {
            old: AstValueAst::Lit(1),
            new: AstValueAst::Lit(2),
        },
    })]
    #[case::constraint_added(AstBondDelta::ModifyConstraint {
        id: AstBondId(2),
        old: None,
        new: Some(AstBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
    })]
    #[case::constraint_removed(AstBondDelta::ModifyConstraint {
        id: AstBondId(2),
        old: Some(AstBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(AstBondDelta::ModifyConstraint {
        id: AstBondId(2),
        old: Some(AstBondConstraintAst::Aromatic(AstBooleanAst::Lit(false))),
        new: Some(AstBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
    })]
    fn test_bond_delta_roundtrip(#[case] delta: AstBondDelta) {
        Python::attach(|py| {
            assert_eq!(BondDelta::from_rust(py, &delta).unwrap().to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::equal(
        AstBondDelta::Add {
            id: AstBondId(2),
            atoms: [AstAtomId(5), AstAtomId(1)],
            ast: AstBondAst::new(AstValueAst::Lit(1)),
        },
        AstBondDelta::Add {
            id: AstBondId(2),
            atoms: [AstAtomId(5), AstAtomId(1)],
            ast: AstBondAst::new(AstValueAst::Lit(1)),
        },
        true,
    )]
    #[case::different_order(
        AstBondDelta::Add {
            id: AstBondId(2),
            atoms: [AstAtomId(5), AstAtomId(1)],
            ast: AstBondAst::new(AstValueAst::Lit(1)),
        },
        AstBondDelta::Add {
            id: AstBondId(2),
            atoms: [AstAtomId(1), AstAtomId(5)],
            ast: AstBondAst::new(AstValueAst::Lit(1)),
        },
        false,
    )]
    fn test_bond_delta_eq(
        #[case] lhs: AstBondDelta,
        #[case] rhs: AstBondDelta,
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
        AstBondDelta::Add {
            id: AstBondId(2),
            atoms: [AstAtomId(5), AstAtomId(1)],
            ast: AstBondAst::new(AstValueAst::Lit(1)),
        },
        "BondDelta.Add(id=2, atoms=(5, 1), ast=BondAst.parse('1'))",
    )]
    #[case::remove(
        AstBondDelta::Remove {
            id: AstBondId(2),
            atoms: [AstAtomId(5), AstAtomId(1)],
            ast: AstBondAst::new(AstValueAst::Lit(1)),
        },
        "BondDelta.Remove(id=2, atoms=(5, 1), ast=BondAst.parse('1'))",
    )]
    #[case::modify_field(
        AstBondDelta::ModifyField {
            id: AstBondId(2),
            change: AstBondFieldChange::Order {
                old: AstValueAst::Lit(1),
                new: AstValueAst::Lit(2),
            },
        },
        "BondDelta.ModifyField(id=2, change=BondFieldChange.Order(old=ValueAst.Lit(1), new=ValueAst.Lit(2)))",
    )]
    #[case::modify_constraint(
        AstBondDelta::ModifyConstraint {
            id: AstBondId(2),
            old: None,
            new: Some(AstBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
        },
        "BondDelta.ModifyConstraint(id=2, old=None, new=BondConstraintAst.Aromatic(BooleanAst.Lit(True)))",
    )]
    fn test_bond_delta_repr(#[case] delta: AstBondDelta, #[case] expected: &str) {
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
    #[case::add(AstBondDelta::Add {
        id: AstBondId(2),
        atoms: [AstAtomId(5), AstAtomId(1)],
        ast: AstBondAst::new(AstValueAst::Lit(1)),
    })]
    #[case::remove(AstBondDelta::Remove {
        id: AstBondId(2),
        atoms: [AstAtomId(5), AstAtomId(1)],
        ast: AstBondAst::new(AstValueAst::Lit(1)),
    })]
    #[case::modify_field(AstBondDelta::ModifyField {
        id: AstBondId(2),
        change: AstBondFieldChange::Order {
            old: AstValueAst::Lit(1),
            new: AstValueAst::Lit(2),
        },
    })]
    #[case::constraint_added(AstBondDelta::ModifyConstraint {
        id: AstBondId(2),
        old: None,
        new: Some(AstBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
    })]
    #[case::constraint_removed(AstBondDelta::ModifyConstraint {
        id: AstBondId(2),
        old: Some(AstBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(AstBondDelta::ModifyConstraint {
        id: AstBondId(2),
        old: Some(AstBondConstraintAst::Aromatic(AstBooleanAst::Lit(false))),
        new: Some(AstBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
    })]
    fn test_bond_delta_inverse(#[case] delta: AstBondDelta) {
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
    #[case::add(AstDativeBondDelta::Add {
        id: AstDativeBondId(1),
        donors: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        acceptor: AstAtomId(3),
        ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
    })]
    #[case::remove(AstDativeBondDelta::Remove {
        id: AstDativeBondId(1),
        donors: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        acceptor: AstAtomId(3),
        ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
    })]
    #[case::modify_field(AstDativeBondDelta::ModifyField {
        id: AstDativeBondId(1),
        change: AstDativeBondFieldChange::Order {
            old: AstValueAst::Lit(1),
            new: AstValueAst::Lit(2),
        },
    })]
    #[case::constraint_added(AstDativeBondDelta::ModifyConstraint {
        id: AstDativeBondId(1),
        old: None,
        new: Some(AstDativeBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
    })]
    #[case::constraint_removed(AstDativeBondDelta::ModifyConstraint {
        id: AstDativeBondId(1),
        old: Some(AstDativeBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(AstDativeBondDelta::ModifyConstraint {
        id: AstDativeBondId(1),
        old: Some(AstDativeBondConstraintAst::Aromatic(AstBooleanAst::Lit(false))),
        new: Some(AstDativeBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
    })]
    fn test_dative_bond_delta_roundtrip(#[case] delta: AstDativeBondDelta) {
        Python::attach(|py| {
            assert_eq!(
                DativeBondDelta::from_rust(py, &delta).unwrap().to_rust(py),
                delta
            );
        });
    }

    #[rstest]
    #[case::equal(
        AstDativeBondDelta::Add {
            id: AstDativeBondId(1),
            donors: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            acceptor: AstAtomId(3),
            ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
        },
        AstDativeBondDelta::Add {
            id: AstDativeBondId(1),
            donors: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            acceptor: AstAtomId(3),
            ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
        },
        true,
    )]
    #[case::different_donor_order(
        AstDativeBondDelta::Add {
            id: AstDativeBondId(1),
            donors: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            acceptor: AstAtomId(3),
            ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
        },
        AstDativeBondDelta::Add {
            id: AstDativeBondId(1),
            donors: vec![AstAtomId(2), AstAtomId(4), AstAtomId(4)],
            acceptor: AstAtomId(3),
            ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
        },
        false,
    )]
    fn test_dative_bond_delta_eq(
        #[case] lhs: AstDativeBondDelta,
        #[case] rhs: AstDativeBondDelta,
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
        AstDativeBondDelta::Add {
            id: AstDativeBondId(1),
            donors: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            acceptor: AstAtomId(3),
            ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
        },
        "DativeBondDelta.Add(id=1, donors=[4, 2, 4], acceptor=3, ast=DativeBondAst.parse('1'))",
    )]
    #[case::remove(
        AstDativeBondDelta::Remove {
            id: AstDativeBondId(1),
            donors: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            acceptor: AstAtomId(3),
            ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
        },
        "DativeBondDelta.Remove(id=1, donors=[4, 2, 4], acceptor=3, ast=DativeBondAst.parse('1'))",
    )]
    #[case::modify_field(
        AstDativeBondDelta::ModifyField {
            id: AstDativeBondId(1),
            change: AstDativeBondFieldChange::Order {
                old: AstValueAst::Lit(1),
                new: AstValueAst::Lit(2),
            },
        },
        "DativeBondDelta.ModifyField(id=1, change=DativeBondFieldChange.Order(old=ValueAst.Lit(1), new=ValueAst.Lit(2)))",
    )]
    #[case::modify_constraint(
        AstDativeBondDelta::ModifyConstraint {
            id: AstDativeBondId(1),
            old: None,
            new: Some(AstDativeBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
        },
        "DativeBondDelta.ModifyConstraint(id=1, old=None, new=DativeBondConstraintAst.Aromatic(BooleanAst.Lit(True)))",
    )]
    fn test_dative_bond_delta_repr(#[case] delta: AstDativeBondDelta, #[case] expected: &str) {
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
    #[case::add(AstDativeBondDelta::Add {
        id: AstDativeBondId(1),
        donors: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        acceptor: AstAtomId(3),
        ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
    })]
    #[case::remove(AstDativeBondDelta::Remove {
        id: AstDativeBondId(1),
        donors: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        acceptor: AstAtomId(3),
        ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
    })]
    #[case::modify_field(AstDativeBondDelta::ModifyField {
        id: AstDativeBondId(1),
        change: AstDativeBondFieldChange::Order {
            old: AstValueAst::Lit(1),
            new: AstValueAst::Lit(2),
        },
    })]
    #[case::constraint_added(AstDativeBondDelta::ModifyConstraint {
        id: AstDativeBondId(1),
        old: None,
        new: Some(AstDativeBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
    })]
    #[case::constraint_removed(AstDativeBondDelta::ModifyConstraint {
        id: AstDativeBondId(1),
        old: Some(AstDativeBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(AstDativeBondDelta::ModifyConstraint {
        id: AstDativeBondId(1),
        old: Some(AstDativeBondConstraintAst::Aromatic(AstBooleanAst::Lit(false))),
        new: Some(AstDativeBondConstraintAst::Aromatic(AstBooleanAst::Lit(true))),
    })]
    fn test_dative_bond_delta_inverse(#[case] delta: AstDativeBondDelta) {
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
    #[case::add(AstAromaticSystemDelta::Add {
        id: AstAromaticSystemId(2),
        atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        ast: AstAromaticSystemAst::from_electrons(vec![1, 1, 1]),
    })]
    #[case::remove(AstAromaticSystemDelta::Remove {
        id: AstAromaticSystemId(2),
        atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        ast: AstAromaticSystemAst::from_electrons(vec![1, 1, 1]),
    })]
    #[case::modify_field(AstAromaticSystemDelta::ModifyField {
        id: AstAromaticSystemId(2),
        change: AstAromaticSystemFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
    })]
    #[case::constraint_added(AstAromaticSystemDelta::ModifyConstraint {
        id: AstAromaticSystemId(2),
        old: None,
        new: Some(AstAromaticSystemConstraintAst::ElectronCount(AstValueAst::Lit(6))),
    })]
    #[case::constraint_removed(AstAromaticSystemDelta::ModifyConstraint {
        id: AstAromaticSystemId(2),
        old: Some(AstAromaticSystemConstraintAst::ElectronCount(AstValueAst::Lit(6))),
        new: None,
    })]
    #[case::constraint_modified(AstAromaticSystemDelta::ModifyConstraint {
        id: AstAromaticSystemId(2),
        old: Some(AstAromaticSystemConstraintAst::ElectronCount(AstValueAst::Lit(5))),
        new: Some(AstAromaticSystemConstraintAst::ElectronCount(AstValueAst::Lit(6))),
    })]
    fn test_aromatic_system_delta_roundtrip(#[case] delta: AstAromaticSystemDelta) {
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
        AstAromaticSystemDelta::Add {
            id: AstAromaticSystemId(2),
            atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            ast: AstAromaticSystemAst::from_electrons(vec![1, 1, 1]),
        },
        AstAromaticSystemDelta::Add {
            id: AstAromaticSystemId(2),
            atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            ast: AstAromaticSystemAst::from_electrons(vec![1, 1, 1]),
        },
        true,
    )]
    #[case::different_atom_order(
        AstAromaticSystemDelta::Add {
            id: AstAromaticSystemId(2),
            atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            ast: AstAromaticSystemAst::from_electrons(vec![1, 1, 1]),
        },
        AstAromaticSystemDelta::Add {
            id: AstAromaticSystemId(2),
            atoms: vec![AstAtomId(2), AstAtomId(4), AstAtomId(4)],
            ast: AstAromaticSystemAst::from_electrons(vec![1, 1, 1]),
        },
        false,
    )]
    fn test_aromatic_system_delta_eq(
        #[case] lhs: AstAromaticSystemDelta,
        #[case] rhs: AstAromaticSystemDelta,
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
        AstAromaticSystemDelta::Add {
            id: AstAromaticSystemId(2),
            atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            ast: AstAromaticSystemAst::from_electrons(vec![1, 1, 1]),
        },
        "AromaticSystemDelta.Add(id=2, atoms=[4, 2, 4], ast=AromaticSystemAst.parse('[1,1,1]'))",
    )]
    #[case::remove(
        AstAromaticSystemDelta::Remove {
            id: AstAromaticSystemId(2),
            atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            ast: AstAromaticSystemAst::from_electrons(vec![1, 1, 1]),
        },
        "AromaticSystemDelta.Remove(id=2, atoms=[4, 2, 4], ast=AromaticSystemAst.parse('[1,1,1]'))",
    )]
    #[case::modify_field(
        AstAromaticSystemDelta::ModifyField {
            id: AstAromaticSystemId(2),
            change: AstAromaticSystemFieldChange::Charge {
                old: AstValueAst::Lit(0),
                new: AstValueAst::Lit(-1),
            },
        },
        "AromaticSystemDelta.ModifyField(id=2, change=AromaticSystemFieldChange.Charge(old=ValueAst.Lit(0), new=ValueAst.Lit(-1)))",
    )]
    #[case::modify_constraint(
        AstAromaticSystemDelta::ModifyConstraint {
            id: AstAromaticSystemId(2),
            old: None,
            new: Some(AstAromaticSystemConstraintAst::ElectronCount(AstValueAst::Lit(6))),
        },
        "AromaticSystemDelta.ModifyConstraint(id=2, old=None, new=AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6)))",
    )]
    fn test_aromatic_system_delta_repr(
        #[case] delta: AstAromaticSystemDelta,
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
    #[case::add(AstAromaticSystemDelta::Add {
        id: AstAromaticSystemId(2),
        atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        ast: AstAromaticSystemAst::from_electrons(vec![1, 1, 1]),
    })]
    #[case::remove(AstAromaticSystemDelta::Remove {
        id: AstAromaticSystemId(2),
        atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        ast: AstAromaticSystemAst::from_electrons(vec![1, 1, 1]),
    })]
    #[case::modify_field(AstAromaticSystemDelta::ModifyField {
        id: AstAromaticSystemId(2),
        change: AstAromaticSystemFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
    })]
    #[case::constraint_added(AstAromaticSystemDelta::ModifyConstraint {
        id: AstAromaticSystemId(2),
        old: None,
        new: Some(AstAromaticSystemConstraintAst::ElectronCount(AstValueAst::Lit(6))),
    })]
    #[case::constraint_removed(AstAromaticSystemDelta::ModifyConstraint {
        id: AstAromaticSystemId(2),
        old: Some(AstAromaticSystemConstraintAst::ElectronCount(AstValueAst::Lit(6))),
        new: None,
    })]
    #[case::constraint_modified(AstAromaticSystemDelta::ModifyConstraint {
        id: AstAromaticSystemId(2),
        old: Some(AstAromaticSystemConstraintAst::ElectronCount(AstValueAst::Lit(5))),
        new: Some(AstAromaticSystemConstraintAst::ElectronCount(AstValueAst::Lit(6))),
    })]
    fn test_aromatic_system_delta_inverse(#[case] delta: AstAromaticSystemDelta) {
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
    #[case::add(AstMulticenterBondDelta::Add {
        id: AstMulticenterBondId(3),
        atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        ast: AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
    })]
    #[case::remove(AstMulticenterBondDelta::Remove {
        id: AstMulticenterBondId(3),
        atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        ast: AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
    })]
    #[case::modify_field(AstMulticenterBondDelta::ModifyField {
        id: AstMulticenterBondId(3),
        change: AstMulticenterBondFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
    })]
    #[case::constraint_added(AstMulticenterBondDelta::ModifyConstraint {
        id: AstMulticenterBondId(3),
        old: None,
        new: Some(AstMulticenterBondConstraintAst::ElectronCount(AstValueAst::Lit(6))),
    })]
    #[case::constraint_removed(AstMulticenterBondDelta::ModifyConstraint {
        id: AstMulticenterBondId(3),
        old: Some(AstMulticenterBondConstraintAst::ElectronCount(AstValueAst::Lit(6))),
        new: None,
    })]
    #[case::constraint_modified(AstMulticenterBondDelta::ModifyConstraint {
        id: AstMulticenterBondId(3),
        old: Some(AstMulticenterBondConstraintAst::ElectronCount(AstValueAst::Lit(5))),
        new: Some(AstMulticenterBondConstraintAst::ElectronCount(AstValueAst::Lit(6))),
    })]
    fn test_multicenter_bond_delta_roundtrip(#[case] delta: AstMulticenterBondDelta) {
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
        AstMulticenterBondDelta::Add {
            id: AstMulticenterBondId(3),
            atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            ast: AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
        },
        AstMulticenterBondDelta::Add {
            id: AstMulticenterBondId(3),
            atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            ast: AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
        },
        true,
    )]
    #[case::different_atom_order(
        AstMulticenterBondDelta::Add {
            id: AstMulticenterBondId(3),
            atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            ast: AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
        },
        AstMulticenterBondDelta::Add {
            id: AstMulticenterBondId(3),
            atoms: vec![AstAtomId(2), AstAtomId(4), AstAtomId(4)],
            ast: AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
        },
        false,
    )]
    fn test_multicenter_bond_delta_eq(
        #[case] lhs: AstMulticenterBondDelta,
        #[case] rhs: AstMulticenterBondDelta,
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
        AstMulticenterBondDelta::Add {
            id: AstMulticenterBondId(3),
            atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            ast: AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
        },
        "MulticenterBondDelta.Add(id=3, atoms=[4, 2, 4], ast=MulticenterBondAst.parse('[1,1,1]'))",
    )]
    #[case::remove(
        AstMulticenterBondDelta::Remove {
            id: AstMulticenterBondId(3),
            atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
            ast: AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
        },
        "MulticenterBondDelta.Remove(id=3, atoms=[4, 2, 4], ast=MulticenterBondAst.parse('[1,1,1]'))",
    )]
    #[case::modify_field(
        AstMulticenterBondDelta::ModifyField {
            id: AstMulticenterBondId(3),
            change: AstMulticenterBondFieldChange::Charge {
                old: AstValueAst::Lit(0),
                new: AstValueAst::Lit(-1),
            },
        },
        "MulticenterBondDelta.ModifyField(id=3, change=MulticenterBondFieldChange.Charge(old=ValueAst.Lit(0), new=ValueAst.Lit(-1)))",
    )]
    #[case::modify_constraint(
        AstMulticenterBondDelta::ModifyConstraint {
            id: AstMulticenterBondId(3),
            old: None,
            new: Some(AstMulticenterBondConstraintAst::ElectronCount(AstValueAst::Lit(6))),
        },
        "MulticenterBondDelta.ModifyConstraint(id=3, old=None, new=MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6)))",
    )]
    fn test_multicenter_bond_delta_repr(
        #[case] delta: AstMulticenterBondDelta,
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
    #[case::add(AstMulticenterBondDelta::Add {
        id: AstMulticenterBondId(3),
        atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        ast: AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
    })]
    #[case::remove(AstMulticenterBondDelta::Remove {
        id: AstMulticenterBondId(3),
        atoms: vec![AstAtomId(4), AstAtomId(2), AstAtomId(4)],
        ast: AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
    })]
    #[case::modify_field(AstMulticenterBondDelta::ModifyField {
        id: AstMulticenterBondId(3),
        change: AstMulticenterBondFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
    })]
    #[case::constraint_added(AstMulticenterBondDelta::ModifyConstraint {
        id: AstMulticenterBondId(3),
        old: None,
        new: Some(AstMulticenterBondConstraintAst::ElectronCount(AstValueAst::Lit(6))),
    })]
    #[case::constraint_removed(AstMulticenterBondDelta::ModifyConstraint {
        id: AstMulticenterBondId(3),
        old: Some(AstMulticenterBondConstraintAst::ElectronCount(AstValueAst::Lit(6))),
        new: None,
    })]
    #[case::constraint_modified(AstMulticenterBondDelta::ModifyConstraint {
        id: AstMulticenterBondId(3),
        old: Some(AstMulticenterBondConstraintAst::ElectronCount(AstValueAst::Lit(5))),
        new: Some(AstMulticenterBondConstraintAst::ElectronCount(AstValueAst::Lit(6))),
    })]
    fn test_multicenter_bond_delta_inverse(#[case] delta: AstMulticenterBondDelta) {
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
    #[case::add(AstNoncovalentBondDelta::Add {
        id: AstNoncovalentBondId(4),
        atoms: [AstAtomId(5), AstAtomId(2)],
        ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
    })]
    #[case::remove(AstNoncovalentBondDelta::Remove {
        id: AstNoncovalentBondId(4),
        atoms: [AstAtomId(5), AstAtomId(2)],
        ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
    })]
    #[case::modify_field(AstNoncovalentBondDelta::ModifyField {
        id: AstNoncovalentBondId(4),
        change: AstNoncovalentBondFieldChange::Kind {
            old: AstNoncovalentBondKindAst::Undetermined,
            new: AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond),
        },
    })]
    #[case::constraint_added(AstNoncovalentBondDelta::ModifyConstraint {
        id: AstNoncovalentBondId(4),
        old: None,
        new: Some(AstNoncovalentBondConstraintAst::Intramolecular(AstBooleanAst::Lit(true))),
    })]
    #[case::constraint_removed(AstNoncovalentBondDelta::ModifyConstraint {
        id: AstNoncovalentBondId(4),
        old: Some(AstNoncovalentBondConstraintAst::Intramolecular(AstBooleanAst::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(AstNoncovalentBondDelta::ModifyConstraint {
        id: AstNoncovalentBondId(4),
        old: Some(AstNoncovalentBondConstraintAst::Intramolecular(AstBooleanAst::Lit(false))),
        new: Some(AstNoncovalentBondConstraintAst::Intramolecular(AstBooleanAst::Lit(true))),
    })]
    fn test_noncovalent_bond_delta_roundtrip(#[case] delta: AstNoncovalentBondDelta) {
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
        AstNoncovalentBondDelta::Add {
            id: AstNoncovalentBondId(4),
            atoms: [AstAtomId(5), AstAtomId(2)],
            ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
        },
        AstNoncovalentBondDelta::Add {
            id: AstNoncovalentBondId(4),
            atoms: [AstAtomId(5), AstAtomId(2)],
            ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
        },
        true,
    )]
    #[case::different_atom_order(
        AstNoncovalentBondDelta::Add {
            id: AstNoncovalentBondId(4),
            atoms: [AstAtomId(5), AstAtomId(2)],
            ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
        },
        AstNoncovalentBondDelta::Add {
            id: AstNoncovalentBondId(4),
            atoms: [AstAtomId(2), AstAtomId(5)],
            ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
        },
        false,
    )]
    fn test_noncovalent_bond_delta_eq(
        #[case] lhs: AstNoncovalentBondDelta,
        #[case] rhs: AstNoncovalentBondDelta,
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
        AstNoncovalentBondDelta::Add {
            id: AstNoncovalentBondId(4),
            atoms: [AstAtomId(5), AstAtomId(2)],
            ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
        },
        "NoncovalentBondDelta.Add(id=4, atoms=(5, 2), ast=NoncovalentBondAst.parse('Hbd'))",
    )]
    #[case::remove(
        AstNoncovalentBondDelta::Remove {
            id: AstNoncovalentBondId(4),
            atoms: [AstAtomId(5), AstAtomId(2)],
            ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
        },
        "NoncovalentBondDelta.Remove(id=4, atoms=(5, 2), ast=NoncovalentBondAst.parse('Hbd'))",
    )]
    #[case::modify_field(
        AstNoncovalentBondDelta::ModifyField {
            id: AstNoncovalentBondId(4),
            change: AstNoncovalentBondFieldChange::Kind {
                old: AstNoncovalentBondKindAst::Undetermined,
                new: AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond),
            },
        },
        "NoncovalentBondDelta.ModifyField(id=4, change=NoncovalentBondFieldChange.Kind(old=NoncovalentBondKindAst.Undetermined(), new=NoncovalentBondKindAst.Lit(NoncovalentBondKind.HydrogenBond)))",
    )]
    #[case::modify_constraint(
        AstNoncovalentBondDelta::ModifyConstraint {
            id: AstNoncovalentBondId(4),
            old: None,
            new: Some(AstNoncovalentBondConstraintAst::Intramolecular(AstBooleanAst::Lit(true))),
        },
        "NoncovalentBondDelta.ModifyConstraint(id=4, old=None, new=NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True)))",
    )]
    fn test_noncovalent_bond_delta_repr(
        #[case] delta: AstNoncovalentBondDelta,
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
    #[case::add(AstNoncovalentBondDelta::Add {
        id: AstNoncovalentBondId(4),
        atoms: [AstAtomId(5), AstAtomId(2)],
        ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
    })]
    #[case::remove(AstNoncovalentBondDelta::Remove {
        id: AstNoncovalentBondId(4),
        atoms: [AstAtomId(5), AstAtomId(2)],
        ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
    })]
    #[case::modify_field(AstNoncovalentBondDelta::ModifyField {
        id: AstNoncovalentBondId(4),
        change: AstNoncovalentBondFieldChange::Kind {
            old: AstNoncovalentBondKindAst::Undetermined,
            new: AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond),
        },
    })]
    #[case::constraint_added(AstNoncovalentBondDelta::ModifyConstraint {
        id: AstNoncovalentBondId(4),
        old: None,
        new: Some(AstNoncovalentBondConstraintAst::Intramolecular(AstBooleanAst::Lit(true))),
    })]
    #[case::constraint_removed(AstNoncovalentBondDelta::ModifyConstraint {
        id: AstNoncovalentBondId(4),
        old: Some(AstNoncovalentBondConstraintAst::Intramolecular(AstBooleanAst::Lit(true))),
        new: None,
    })]
    #[case::constraint_modified(AstNoncovalentBondDelta::ModifyConstraint {
        id: AstNoncovalentBondId(4),
        old: Some(AstNoncovalentBondConstraintAst::Intramolecular(AstBooleanAst::Lit(false))),
        new: Some(AstNoncovalentBondConstraintAst::Intramolecular(AstBooleanAst::Lit(true))),
    })]
    fn test_noncovalent_bond_delta_inverse(#[case] delta: AstNoncovalentBondDelta) {
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
    #[case::add(AstStereoAtomDelta::Add {
        id: AstStereoAtomId(5),
        site: AstAtomId(3),
        ligands: vec![
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
        ],
        ast: AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(0)),
    })]
    #[case::remove(AstStereoAtomDelta::Remove {
        id: AstStereoAtomId(5),
        site: AstAtomId(3),
        ligands: vec![
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
        ],
        ast: AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(0)),
    })]
    #[case::modify_field(AstStereoAtomDelta::ModifyField {
        id: AstStereoAtomId(5),
        change: AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(0),
            ),
        },
    })]
    #[case::constraint_added_with_kind(AstStereoAtomDelta::ModifyConstraint {
        id: AstStereoAtomId(5),
        kind: Some(AstStereoKind::Tetrahedral),
        old: None,
        new: Some(AstStereoAtomConstraintAst::Stereogenicity(
            AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
        )),
    })]
    #[case::constraint_removed_without_kind(AstStereoAtomDelta::ModifyConstraint {
        id: AstStereoAtomId(5),
        kind: None,
        old: Some(AstStereoAtomConstraintAst::Stereogenicity(
            AstStereogenicityAst::Undetermined,
        )),
        new: None,
    })]
    #[case::constraint_modified(AstStereoAtomDelta::ModifyConstraint {
        id: AstStereoAtomId(5),
        kind: Some(AstStereoKind::Tetrahedral),
        old: Some(AstStereoAtomConstraintAst::Stereogenicity(
            AstStereogenicityAst::Undetermined,
        )),
        new: Some(AstStereoAtomConstraintAst::Stereogenicity(
            AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
        )),
    })]
    #[case::apply(AstStereoAtomDelta::Apply {
        id: AstStereoAtomId(5),
        kind: AstStereoKind::Tetrahedral,
        permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
    })]
    #[case::swap(AstStereoAtomDelta::Swap {
        id: AstStereoAtomId(5),
        kind: AstStereoKind::Tetrahedral,
    })]
    #[case::mirror(AstStereoAtomDelta::Mirror {
        id: AstStereoAtomId(5),
        kind: AstStereoKind::Tetrahedral,
    })]
    fn test_stereo_atom_delta_roundtrip(#[case] delta: AstStereoAtomDelta) {
        Python::attach(|py| {
            assert_eq!(
                StereoAtomDelta::from_rust(py, &delta).unwrap().to_rust(py),
                delta
            );
        });
    }

    #[rstest]
    #[case::equal(
        AstStereoAtomDelta::Apply {
            id: AstStereoAtomId(5),
            kind: AstStereoKind::Tetrahedral,
            permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
        },
        AstStereoAtomDelta::Apply {
            id: AstStereoAtomId(5),
            kind: AstStereoKind::Tetrahedral,
            permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
        },
        true,
    )]
    #[case::different_ligand_order(
        AstStereoAtomDelta::Add {
            id: AstStereoAtomId(5),
            site: AstAtomId(3),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
            ],
            ast: AstStereoAtomAst::new(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(0),
            ),
        },
        AstStereoAtomDelta::Add {
            id: AstStereoAtomId(5),
            site: AstAtomId(3),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            ],
            ast: AstStereoAtomAst::new(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(0),
            ),
        },
        false,
    )]
    #[case::different_permutation(
        AstStereoAtomDelta::Apply {
            id: AstStereoAtomId(5),
            kind: AstStereoKind::Tetrahedral,
            permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
        },
        AstStereoAtomDelta::Apply {
            id: AstStereoAtomId(5),
            kind: AstStereoKind::Tetrahedral,
            permutation: PermPermutation::from_image(4, &[2, 0, 1, 3]),
        },
        false,
    )]
    fn test_stereo_atom_delta_eq(
        #[case] lhs: AstStereoAtomDelta,
        #[case] rhs: AstStereoAtomDelta,
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
        AstStereoAtomDelta::Add {
            id: AstStereoAtomId(5),
            site: AstAtomId(3),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::LonePair),
            ],
            ast: AstStereoAtomAst::new(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(0),
            ),
        },
        "StereoAtomDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair)], ast=StereoAtomAst.parse('Th0'))",
    )]
    #[case::remove(
        AstStereoAtomDelta::Remove {
            id: AstStereoAtomId(5),
            site: AstAtomId(3),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::LonePair),
            ],
            ast: AstStereoAtomAst::new(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(0),
            ),
        },
        "StereoAtomDelta.Remove(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair)], ast=StereoAtomAst.parse('Th0'))",
    )]
    #[case::modify_field(
        AstStereoAtomDelta::ModifyField {
            id: AstStereoAtomId(5),
            change: AstStereoAtomFieldChange::Configuration {
                old: AstStereoConfigurationAst::Undetermined,
                new: AstStereoConfigurationAst::Kinded(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(0),
                ),
            },
        },
        "StereoAtomDelta.ModifyField(id=5, change=StereoAtomFieldChange.Configuration(old=StereoConfigurationAst.Undetermined(), new=StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Lit(0))))",
    )]
    #[case::modify_constraint(
        AstStereoAtomDelta::ModifyConstraint {
            id: AstStereoAtomId(5),
            kind: Some(AstStereoKind::Tetrahedral),
            old: None,
            new: Some(AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Undetermined,
            )),
        },
        "StereoAtomDelta.ModifyConstraint(id=5, kind=StereoKind.Tetrahedral, old=None, new=StereoAtomConstraintAst.Stereogenicity(StereogenicityAst.Undetermined()))",
    )]
    #[case::apply(
        AstStereoAtomDelta::Apply {
            id: AstStereoAtomId(5),
            kind: AstStereoKind::Tetrahedral,
            permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
        },
        "StereoAtomDelta.Apply(id=5, kind=StereoKind.Tetrahedral, permutation=Permutation([1, 2, 0, 3]))",
    )]
    #[case::swap(
        AstStereoAtomDelta::Swap {
            id: AstStereoAtomId(5),
            kind: AstStereoKind::Tetrahedral,
        },
        "StereoAtomDelta.Swap(id=5, kind=StereoKind.Tetrahedral)",
    )]
    #[case::mirror(
        AstStereoAtomDelta::Mirror {
            id: AstStereoAtomId(5),
            kind: AstStereoKind::Tetrahedral,
        },
        "StereoAtomDelta.Mirror(id=5, kind=StereoKind.Tetrahedral)",
    )]
    fn test_stereo_atom_delta_repr(#[case] delta: AstStereoAtomDelta, #[case] expected: &str) {
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
    #[case::add(AstStereoAtomDelta::Add {
        id: AstStereoAtomId(5),
        site: AstAtomId(3),
        ligands: vec![
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
        ],
        ast: AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(0)),
    })]
    #[case::remove(AstStereoAtomDelta::Remove {
        id: AstStereoAtomId(5),
        site: AstAtomId(3),
        ligands: vec![
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
        ],
        ast: AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(0)),
    })]
    #[case::modify_field(AstStereoAtomDelta::ModifyField {
        id: AstStereoAtomId(5),
        change: AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(0),
            ),
        },
    })]
    #[case::constraint_added_with_kind(AstStereoAtomDelta::ModifyConstraint {
        id: AstStereoAtomId(5),
        kind: Some(AstStereoKind::Tetrahedral),
        old: None,
        new: Some(AstStereoAtomConstraintAst::Stereogenicity(
            AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
        )),
    })]
    #[case::constraint_removed_without_kind(AstStereoAtomDelta::ModifyConstraint {
        id: AstStereoAtomId(5),
        kind: None,
        old: Some(AstStereoAtomConstraintAst::Stereogenicity(
            AstStereogenicityAst::Undetermined,
        )),
        new: None,
    })]
    #[case::constraint_modified(AstStereoAtomDelta::ModifyConstraint {
        id: AstStereoAtomId(5),
        kind: Some(AstStereoKind::Tetrahedral),
        old: Some(AstStereoAtomConstraintAst::Stereogenicity(
            AstStereogenicityAst::Undetermined,
        )),
        new: Some(AstStereoAtomConstraintAst::Stereogenicity(
            AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
        )),
    })]
    #[case::apply(AstStereoAtomDelta::Apply {
        id: AstStereoAtomId(5),
        kind: AstStereoKind::Tetrahedral,
        permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
    })]
    #[case::swap(AstStereoAtomDelta::Swap {
        id: AstStereoAtomId(5),
        kind: AstStereoKind::Tetrahedral,
    })]
    #[case::mirror(AstStereoAtomDelta::Mirror {
        id: AstStereoAtomId(5),
        kind: AstStereoKind::Tetrahedral,
    })]
    fn test_stereo_atom_delta_inverse(#[case] delta: AstStereoAtomDelta) {
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
    #[case::add(AstStereoBondDelta::Add {
        id: AstStereoBondId(5),
        site: AstBondId(3),
        ligands: vec![
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
        ],
        ast: AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(0)),
    })]
    #[case::remove(AstStereoBondDelta::Remove {
        id: AstStereoBondId(5),
        site: AstBondId(3),
        ligands: vec![
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
        ],
        ast: AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(0)),
    })]
    #[case::modify_field(AstStereoBondDelta::ModifyField {
        id: AstStereoBondId(5),
        change: AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Lit(0),
            ),
        },
    })]
    #[case::constraint_added_with_kind(AstStereoBondDelta::ModifyConstraint {
        id: AstStereoBondId(5),
        kind: Some(AstStereoKind::CisTrans),
        old: None,
        new: Some(AstStereoBondConstraintAst::Stereogenicity(
            AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
        )),
    })]
    #[case::constraint_removed_without_kind(AstStereoBondDelta::ModifyConstraint {
        id: AstStereoBondId(5),
        kind: None,
        old: Some(AstStereoBondConstraintAst::Stereogenicity(
            AstStereogenicityAst::Undetermined,
        )),
        new: None,
    })]
    #[case::constraint_modified(AstStereoBondDelta::ModifyConstraint {
        id: AstStereoBondId(5),
        kind: Some(AstStereoKind::CisTrans),
        old: Some(AstStereoBondConstraintAst::Stereogenicity(
            AstStereogenicityAst::Undetermined,
        )),
        new: Some(AstStereoBondConstraintAst::Stereogenicity(
            AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
        )),
    })]
    #[case::apply(AstStereoBondDelta::Apply {
        id: AstStereoBondId(5),
        kind: AstStereoKind::CisTrans,
        permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
    })]
    #[case::swap(AstStereoBondDelta::Swap {
        id: AstStereoBondId(5),
        kind: AstStereoKind::CisTrans,
    })]
    #[case::mirror(AstStereoBondDelta::Mirror {
        id: AstStereoBondId(5),
        kind: AstStereoKind::CisTrans,
    })]
    fn test_stereo_bond_delta_roundtrip(#[case] delta: AstStereoBondDelta) {
        Python::attach(|py| {
            assert_eq!(
                StereoBondDelta::from_rust(py, &delta).unwrap().to_rust(py),
                delta
            );
        });
    }

    #[rstest]
    #[case::equal(
        AstStereoBondDelta::Apply {
            id: AstStereoBondId(5),
            kind: AstStereoKind::CisTrans,
            permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
        },
        AstStereoBondDelta::Apply {
            id: AstStereoBondId(5),
            kind: AstStereoKind::CisTrans,
            permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
        },
        true,
    )]
    #[case::different_ligand_order(
        AstStereoBondDelta::Add {
            id: AstStereoBondId(5),
            site: AstBondId(3),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
            ],
            ast: AstStereoBondAst::new(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Lit(0),
            ),
        },
        AstStereoBondDelta::Add {
            id: AstStereoBondId(5),
            site: AstBondId(3),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            ],
            ast: AstStereoBondAst::new(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Lit(0),
            ),
        },
        false,
    )]
    #[case::different_permutation(
        AstStereoBondDelta::Apply {
            id: AstStereoBondId(5),
            kind: AstStereoKind::CisTrans,
            permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
        },
        AstStereoBondDelta::Apply {
            id: AstStereoBondId(5),
            kind: AstStereoKind::CisTrans,
            permutation: PermPermutation::from_image(4, &[2, 0, 1, 3]),
        },
        false,
    )]
    fn test_stereo_bond_delta_eq(
        #[case] lhs: AstStereoBondDelta,
        #[case] rhs: AstStereoBondDelta,
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
        AstStereoBondDelta::Add {
            id: AstStereoBondId(5),
            site: AstBondId(3),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::LonePair),
            ],
            ast: AstStereoBondAst::new(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Lit(0),
            ),
        },
        "StereoBondDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair)], ast=StereoBondAst.parse('Ct0'))",
    )]
    #[case::remove(
        AstStereoBondDelta::Remove {
            id: AstStereoBondId(5),
            site: AstBondId(3),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::LonePair),
            ],
            ast: AstStereoBondAst::new(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Lit(0),
            ),
        },
        "StereoBondDelta.Remove(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom), StereoLigand(atom_id=2, kind=StereoLigandKind.LonePair)], ast=StereoBondAst.parse('Ct0'))",
    )]
    #[case::modify_field(
        AstStereoBondDelta::ModifyField {
            id: AstStereoBondId(5),
            change: AstStereoBondFieldChange::Configuration {
                old: AstStereoConfigurationAst::Undetermined,
                new: AstStereoConfigurationAst::Kinded(
                    AstStereoKind::CisTrans,
                    AstStereoCosetAst::Lit(0),
                ),
            },
        },
        "StereoBondDelta.ModifyField(id=5, change=StereoBondFieldChange.Configuration(old=StereoConfigurationAst.Undetermined(), new=StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Lit(0))))",
    )]
    #[case::modify_constraint(
        AstStereoBondDelta::ModifyConstraint {
            id: AstStereoBondId(5),
            kind: Some(AstStereoKind::CisTrans),
            old: None,
            new: Some(AstStereoBondConstraintAst::Stereogenicity(
                AstStereogenicityAst::Undetermined,
            )),
        },
        "StereoBondDelta.ModifyConstraint(id=5, kind=StereoKind.CisTrans, old=None, new=StereoBondConstraintAst.Stereogenicity(StereogenicityAst.Undetermined()))",
    )]
    #[case::apply(
        AstStereoBondDelta::Apply {
            id: AstStereoBondId(5),
            kind: AstStereoKind::CisTrans,
            permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
        },
        "StereoBondDelta.Apply(id=5, kind=StereoKind.CisTrans, permutation=Permutation([1, 2, 0, 3]))",
    )]
    #[case::swap(
        AstStereoBondDelta::Swap {
            id: AstStereoBondId(5),
            kind: AstStereoKind::CisTrans,
        },
        "StereoBondDelta.Swap(id=5, kind=StereoKind.CisTrans)",
    )]
    #[case::mirror(
        AstStereoBondDelta::Mirror {
            id: AstStereoBondId(5),
            kind: AstStereoKind::CisTrans,
        },
        "StereoBondDelta.Mirror(id=5, kind=StereoKind.CisTrans)",
    )]
    fn test_stereo_bond_delta_repr(#[case] delta: AstStereoBondDelta, #[case] expected: &str) {
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
    #[case::add(AstStereoBondDelta::Add {
        id: AstStereoBondId(5),
        site: AstBondId(3),
        ligands: vec![
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
        ],
        ast: AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(0)),
    })]
    #[case::remove(AstStereoBondDelta::Remove {
        id: AstStereoBondId(5),
        site: AstBondId(3),
        ligands: vec![
            AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
        ],
        ast: AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(0)),
    })]
    #[case::modify_field(AstStereoBondDelta::ModifyField {
        id: AstStereoBondId(5),
        change: AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Lit(0),
            ),
        },
    })]
    #[case::constraint_added_with_kind(AstStereoBondDelta::ModifyConstraint {
        id: AstStereoBondId(5),
        kind: Some(AstStereoKind::CisTrans),
        old: None,
        new: Some(AstStereoBondConstraintAst::Stereogenicity(
            AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
        )),
    })]
    #[case::constraint_removed_without_kind(AstStereoBondDelta::ModifyConstraint {
        id: AstStereoBondId(5),
        kind: None,
        old: Some(AstStereoBondConstraintAst::Stereogenicity(
            AstStereogenicityAst::Undetermined,
        )),
        new: None,
    })]
    #[case::constraint_modified(AstStereoBondDelta::ModifyConstraint {
        id: AstStereoBondId(5),
        kind: Some(AstStereoKind::CisTrans),
        old: Some(AstStereoBondConstraintAst::Stereogenicity(
            AstStereogenicityAst::Undetermined,
        )),
        new: Some(AstStereoBondConstraintAst::Stereogenicity(
            AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
        )),
    })]
    #[case::apply(AstStereoBondDelta::Apply {
        id: AstStereoBondId(5),
        kind: AstStereoKind::CisTrans,
        permutation: PermPermutation::from_image(4, &[1, 2, 0, 3]),
    })]
    #[case::swap(AstStereoBondDelta::Swap {
        id: AstStereoBondId(5),
        kind: AstStereoKind::CisTrans,
    })]
    #[case::mirror(AstStereoBondDelta::Mirror {
        id: AstStereoBondId(5),
        kind: AstStereoKind::CisTrans,
    })]
    fn test_stereo_bond_delta_inverse(#[case] delta: AstStereoBondDelta) {
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
    #[case::add_leaf(AstConstraintDelta::Add(AstConstraint::Atom(
        AstAtomId(3),
        AstAtomConstraintAst::degree(2),
    )))]
    #[case::remove_recursive(AstConstraintDelta::Remove(AstConstraint::And(vec![
        AstConstraint::Atom(AstAtomId(7), AstAtomConstraintAst::valence(4)),
        AstConstraint::Not(Box::new(AstConstraint::Or(Vec::new()))),
    ])))]
    fn test_constraint_delta_roundtrip(#[case] delta: AstConstraintDelta) {
        Python::attach(|py| {
            let binding = ConstraintDelta::from_rust(py, &delta).unwrap();
            assert_eq!(binding.to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::equal(
        AstConstraintDelta::Add(AstConstraint::Atom(
            AstAtomId(3),
            AstAtomConstraintAst::degree(2),
        )),
        AstConstraintDelta::Add(AstConstraint::Atom(
            AstAtomId(3),
            AstAtomConstraintAst::degree(2),
        )),
        true
    )]
    #[case::variant(
        AstConstraintDelta::Add(AstConstraint::Atom(
            AstAtomId(3),
            AstAtomConstraintAst::degree(2),
        )),
        AstConstraintDelta::Remove(AstConstraint::Atom(
            AstAtomId(3),
            AstAtomConstraintAst::degree(2),
        )),
        false
    )]
    #[case::constraint(
        AstConstraintDelta::Add(AstConstraint::Atom(
            AstAtomId(3),
            AstAtomConstraintAst::degree(2),
        )),
        AstConstraintDelta::Add(AstConstraint::Atom(
            AstAtomId(3),
            AstAtomConstraintAst::valence(2),
        )),
        false
    )]
    fn test_constraint_delta_eq(
        #[case] left: AstConstraintDelta,
        #[case] right: AstConstraintDelta,
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
        AstConstraintDelta::Add(AstConstraint::Atom(
            AstAtomId(3),
            AstAtomConstraintAst::degree(2),
        )),
        "ConstraintDelta.Add(constraint=Constraint.Atom(3, AtomConstraintAst.Degree(ValueAst.Lit(2))))",
    )]
    #[case::remove_recursive(
        AstConstraintDelta::Remove(AstConstraint::And(vec![
            AstConstraint::Atom(AstAtomId(7), AstAtomConstraintAst::valence(4)),
            AstConstraint::Not(Box::new(AstConstraint::Or(Vec::new()))),
        ])),
        "ConstraintDelta.Remove(constraint=Constraint.And([Constraint.Atom(7, AtomConstraintAst.Valence(ValueAst.Lit(4))), Constraint.Not(Constraint.Or([]))]))",
    )]
    fn test_constraint_delta_repr(#[case] delta: AstConstraintDelta, #[case] expected: &str) {
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
    #[case::add_leaf(AstConstraintDelta::Add(AstConstraint::Atom(
        AstAtomId(3),
        AstAtomConstraintAst::degree(2),
    )))]
    #[case::remove_recursive(AstConstraintDelta::Remove(AstConstraint::And(vec![
        AstConstraint::Atom(AstAtomId(7), AstAtomConstraintAst::valence(4)),
        AstConstraint::Not(Box::new(AstConstraint::Or(Vec::new()))),
    ])))]
    fn test_constraint_delta_inverse(#[case] delta: AstConstraintDelta) {
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
    #[case::atom(AstDelta::Atom(AstAtomDelta::Add {
        id: AstAtomId(3),
        ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
    }))]
    #[case::bond(AstDelta::Bond(AstBondDelta::Add {
        id: AstBondId(2),
        atoms: [AstAtomId(5), AstAtomId(1)],
        ast: AstBondAst::new(AstValueAst::Lit(1)),
    }))]
    #[case::dative_bond(AstDelta::DativeBond(AstDativeBondDelta::Add {
        id: AstDativeBondId(1),
        donors: vec![AstAtomId(4), AstAtomId(2)],
        acceptor: AstAtomId(3),
        ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
    }))]
    #[case::aromatic_system(AstDelta::AromaticSystem(AstAromaticSystemDelta::Add {
        id: AstAromaticSystemId(2),
        atoms: vec![AstAtomId(4), AstAtomId(2)],
        ast: AstAromaticSystemAst::from_electrons(vec![1, 1]),
    }))]
    #[case::multicenter_bond(AstDelta::MulticenterBond(AstMulticenterBondDelta::Add {
        id: AstMulticenterBondId(3),
        atoms: vec![AstAtomId(4), AstAtomId(2)],
        ast: AstMulticenterBondAst::from_electrons(vec![1, 1]),
    }))]
    #[case::noncovalent_bond(AstDelta::NoncovalentBond(AstNoncovalentBondDelta::Add {
        id: AstNoncovalentBondId(4),
        atoms: [AstAtomId(5), AstAtomId(2)],
        ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
    }))]
    #[case::stereo_atom(AstDelta::StereoAtom(AstStereoAtomDelta::Add {
        id: AstStereoAtomId(5),
        site: AstAtomId(3),
        ligands: vec![AstStereoLigand::new(
            AstAtomId(4),
            AstStereoLigandKind::Atom,
        )],
        ast: AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(0)),
    }))]
    #[case::stereo_bond(AstDelta::StereoBond(AstStereoBondDelta::Add {
        id: AstStereoBondId(5),
        site: AstBondId(3),
        ligands: vec![AstStereoLigand::new(
            AstAtomId(4),
            AstStereoLigandKind::Atom,
        )],
        ast: AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(0)),
    }))]
    #[case::constraint(AstDelta::Constraint(AstConstraintDelta::Add(AstConstraint::Atom(
        AstAtomId(3),
        AstAtomConstraintAst::degree(2)
    ),)))]
    fn test_delta_roundtrip(#[case] delta: AstDelta) {
        Python::attach(|py| {
            let binding = Delta::from_rust(py, &delta).unwrap();
            assert_eq!(binding.to_rust(py), delta);
        });
    }

    #[rstest]
    #[case::equal(
        AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        }),
        AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        }),
        true,
    )]
    #[case::outer_variant(
        AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        }),
        AstDelta::Constraint(AstConstraintDelta::Add(AstConstraint::Atom(
            AstAtomId(3),
            AstAtomConstraintAst::degree(2),
        ))),
        false,
    )]
    #[case::child(
        AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        }),
        AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(4),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        }),
        false,
    )]
    fn test_delta_eq(#[case] left: AstDelta, #[case] right: AstDelta, #[case] expected: bool) {
        Python::attach(|py| {
            let left = Delta::from_rust(py, &left).unwrap();
            let right = Delta::from_rust(py, &right).unwrap();
            assert_eq!(left.__eq__(&right, py), expected);
        });
    }

    #[rstest]
    #[case::atom(
        AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        }),
        "Delta.Atom(AtomDelta.Add(id=3, ast=AtomAst.parse('C')))"
    )]
    #[case::stereo_atom(
        AstDelta::StereoAtom(AstStereoAtomDelta::Add {
            id: AstStereoAtomId(5),
            site: AstAtomId(3),
            ligands: vec![AstStereoLigand::new(
                AstAtomId(4),
                AstStereoLigandKind::Atom,
            )],
            ast: AstStereoAtomAst::new(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(0),
            ),
        }),
        "Delta.StereoAtom(StereoAtomDelta.Add(id=5, site=3, ligands=[StereoLigand(atom_id=4, kind=StereoLigandKind.Atom)], ast=StereoAtomAst.parse('Th0')))"
    )]
    #[case::constraint(
        AstDelta::Constraint(AstConstraintDelta::Add(AstConstraint::Atom(
            AstAtomId(3),
            AstAtomConstraintAst::degree(2),
        ))),
        "Delta.Constraint(ConstraintDelta.Add(constraint=Constraint.Atom(3, AtomConstraintAst.Degree(ValueAst.Lit(2)))))"
    )]
    fn test_delta_repr(#[case] delta: AstDelta, #[case] expected: &str) {
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
    #[case::atom(AstDelta::Atom(AstAtomDelta::Add {
        id: AstAtomId(3),
        ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
    }))]
    #[case::bond(AstDelta::Bond(AstBondDelta::Add {
        id: AstBondId(2),
        atoms: [AstAtomId(5), AstAtomId(1)],
        ast: AstBondAst::new(AstValueAst::Lit(1)),
    }))]
    #[case::dative_bond(AstDelta::DativeBond(AstDativeBondDelta::Add {
        id: AstDativeBondId(1),
        donors: vec![AstAtomId(4), AstAtomId(2)],
        acceptor: AstAtomId(3),
        ast: AstDativeBondAst::new(AstValueAst::Lit(1)),
    }))]
    #[case::aromatic_system(AstDelta::AromaticSystem(AstAromaticSystemDelta::Add {
        id: AstAromaticSystemId(2),
        atoms: vec![AstAtomId(4), AstAtomId(2)],
        ast: AstAromaticSystemAst::from_electrons(vec![1, 1]),
    }))]
    #[case::multicenter_bond(AstDelta::MulticenterBond(AstMulticenterBondDelta::Add {
        id: AstMulticenterBondId(3),
        atoms: vec![AstAtomId(4), AstAtomId(2)],
        ast: AstMulticenterBondAst::from_electrons(vec![1, 1]),
    }))]
    #[case::noncovalent_bond(AstDelta::NoncovalentBond(AstNoncovalentBondDelta::Add {
        id: AstNoncovalentBondId(4),
        atoms: [AstAtomId(5), AstAtomId(2)],
        ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
    }))]
    #[case::stereo_atom(AstDelta::StereoAtom(AstStereoAtomDelta::Add {
        id: AstStereoAtomId(5),
        site: AstAtomId(3),
        ligands: vec![AstStereoLigand::new(
            AstAtomId(4),
            AstStereoLigandKind::Atom,
        )],
        ast: AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(0)),
    }))]
    #[case::stereo_bond(AstDelta::StereoBond(AstStereoBondDelta::Add {
        id: AstStereoBondId(5),
        site: AstBondId(3),
        ligands: vec![AstStereoLigand::new(
            AstAtomId(4),
            AstStereoLigandKind::Atom,
        )],
        ast: AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(0)),
    }))]
    #[case::constraint(AstDelta::Constraint(AstConstraintDelta::Add(AstConstraint::Atom(
        AstAtomId(3),
        AstAtomConstraintAst::degree(2)
    ),)))]
    fn test_delta_inverse(#[case] delta: AstDelta) {
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
        AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        }),
        AstDelta::Constraint(AstConstraintDelta::Add(AstConstraint::Atom(
            AstAtomId(3),
            AstAtomConstraintAst::degree(2),
        ))),
    ])]
    fn test_deltas_new(#[case] entries: Vec<AstDelta>) {
        Python::attach(|py| {
            let python_entries = entries
                .iter()
                .map(|entry| into_py_variant(py, Delta::from_rust(py, entry).unwrap()).unwrap())
                .collect();
            assert_eq!(
                Deltas::new(py, python_entries).to_rust(),
                entries.into_iter().collect()
            );
        });
    }

    #[rstest]
    #[case::equal(
        vec![AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        })],
        vec![AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        })],
        true,
    )]
    #[case::different(
        vec![AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        })],
        vec![AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(4),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        })],
        false,
    )]
    fn test_deltas_eq(
        #[case] left: Vec<AstDelta>,
        #[case] right: Vec<AstDelta>,
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
            AstDelta::Atom(AstAtomDelta::Add {
                id: AstAtomId(3),
                ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
            }),
            AstDelta::Constraint(AstConstraintDelta::Add(AstConstraint::Atom(
                AstAtomId(3),
                AstAtomConstraintAst::degree(2),
            ))),
        ],
        "Deltas([Delta.Atom(AtomDelta.Add(id=3, ast=AtomAst.parse('C'))), Delta.Constraint(ConstraintDelta.Add(constraint=Constraint.Atom(3, AtomConstraintAst.Degree(ValueAst.Lit(2)))))])",
    )]
    fn test_deltas_repr(#[case] entries: Vec<AstDelta>, #[case] expected: &str) {
        Python::attach(|py| {
            let deltas = Deltas::from_rust(entries.into_iter().collect());
            assert_eq!(deltas.__repr__(py).unwrap(), expected);
        });
    }

    #[rstest]
    #[case::empty(Vec::new(), 0)]
    #[case::populated(
        vec![AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(3),
            ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
        })],
        1,
    )]
    fn test_deltas_len(#[case] entries: Vec<AstDelta>, #[case] expected: usize) {
        assert_eq!(
            Deltas::from_rust(entries.into_iter().collect()).__len__(),
            expected
        );
    }

    #[rstest]
    #[case::positive(0, AstDelta::Atom(AstAtomDelta::Add {
        id: AstAtomId(3),
        ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
    }))]
    #[case::negative(-1, AstDelta::Constraint(AstConstraintDelta::Add(
        AstConstraint::Atom(AstAtomId(3), AstAtomConstraintAst::degree(2)),
    )))]
    fn test_deltas_getitem(#[case] index: isize, #[case] expected: AstDelta) {
        Python::attach(|py| {
            let deltas = Deltas::from_rust(
                vec![
                    AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(3),
                        ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
                    }),
                    AstDelta::Constraint(AstConstraintDelta::Add(AstConstraint::Atom(
                        AstAtomId(3),
                        AstAtomConstraintAst::degree(2),
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
                    AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(3),
                        ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
                    }),
                    AstDelta::Constraint(AstConstraintDelta::Add(AstConstraint::Atom(
                        AstAtomId(3),
                        AstAtomConstraintAst::degree(2),
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
            AstDelta::Atom(AstAtomDelta::Add {
                id: AstAtomId(3),
                ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
            }),
            AstDelta::Constraint(AstConstraintDelta::Add(AstConstraint::Atom(
                AstAtomId(3),
                AstAtomConstraintAst::degree(2),
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
    #[case::populated(vec![AstDelta::Atom(AstAtomDelta::Add {
        id: AstAtomId(3),
        ast: AstAtomAst::new(AstElementAst::Lit(ChemElement::C)),
    })])]
    fn test_deltas_roundtrip(#[case] entries: Vec<AstDelta>) {
        let rust: AstDeltas = entries.into_iter().collect();
        assert_eq!(Deltas::from_rust(rust.clone()).to_rust(), rust);
    }
}
