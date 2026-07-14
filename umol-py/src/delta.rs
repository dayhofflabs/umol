//! Python bindings for resolved molecule deltas and their field-change payloads.

use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticSystemFieldChange as AstAromaticSystemFieldChange, AtomAst as AstAtomAst,
    AtomDelta as AstAtomDelta, AtomFieldChange as AstAtomFieldChange, AtomId as AstAtomId,
    BondAst as AstBondAst, BondDelta as AstBondDelta, BondFieldChange as AstBondFieldChange,
    BondId as AstBondId, DativeBondAst as AstDativeBondAst, DativeBondDelta as AstDativeBondDelta,
    DativeBondFieldChange as AstDativeBondFieldChange, DativeBondId as AstDativeBondId,
    MulticenterBondFieldChange as AstMulticenterBondFieldChange,
    NoncovalentBondFieldChange as AstNoncovalentBondFieldChange,
    StereoAtomFieldChange as AstStereoAtomFieldChange,
    StereoBondFieldChange as AstStereoBondFieldChange,
};

use crate::atom::{AtomAst, ElementAst, IsotopeMassAst};
use crate::bond::BondAst;
use crate::constraint::atom::AtomConstraintAst;
use crate::constraint::bond::BondConstraintAst;
use crate::constraint::dative::DativeBondConstraintAst;
use crate::convert::into_py_variant;
use crate::dative::DativeBondAst;
use crate::electrons::ElectronCountsAst;
use crate::noncovalent::NoncovalentBondKindAst;
use crate::spin::SpinStateAst;
use crate::stereo::StereoConfigurationAst;
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

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AtomConstraintAst as AstAtomConstraintAst, BondConstraintAst as AstBondConstraintAst,
        BooleanAst as AstBooleanAst, DativeBondConstraintAst as AstDativeBondConstraintAst,
        ElectronCountsAst as AstElectronCountsAst, ElementAst as AstElementAst,
        IsotopeMassAst as AstIsotopeMassAst, NoncovalentBondKind as AstNoncovalentBondKind,
        NoncovalentBondKindAst as AstNoncovalentBondKindAst, SpinStateAst as AstSpinStateAst,
        StereoConfigurationAst as AstStereoConfigurationAst, StereoCosetAst as AstStereoCosetAst,
        StereoKind as AstStereoKind, ValueAst as AstValueAst,
    };
    use umol_chem::element::Element as ChemElement;

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
}
