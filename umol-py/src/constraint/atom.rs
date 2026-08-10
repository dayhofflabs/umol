//! Constraint payloads specific to atoms.

use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use umol_graph_ir::ir::{
    AromaticValence as GraphIrAromaticValence, AromaticValenceForm as GraphIrAromaticValenceForm,
    AsLit, AtomConstraintForm as GraphIrAtomConstraintForm,
    AtomConstraintKey as GraphIrAtomConstraintKey,
    AtomConstraintsForm as GraphIrAtomConstraintsForm, AtomId as GraphIrAtomId,
    MulticenterValence as GraphIrMulticenterValence,
    MulticenterValenceForm as GraphIrMulticenterValenceForm, RingScope as GraphIrRingScope,
    TetrahedralStereoForm as GraphIrTetrahedralStereoForm,
};

use super::ring::{RingMembershipForm, RingScope};
use crate::atom::AtomForm;
use crate::convert::{hash_rust, into_py_variant, variant_repr};
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;
use crate::stereo::{TetrahedralConfiguration, TetrahedralStereoForm};
use crate::value::{NumForm, NumLike};

/// Exact ground aromatic-valence state.
#[pyclass(from_py_object)]
#[derive(Clone, Copy)]
pub enum AromaticValence {
    NotAromatic(),
    Aromatic(i64),
}

#[pymethods]
impl AromaticValence {
    /// Aromatic-valence count, identifying explicit non-aromaticity with zero.
    fn valence_count(&self) -> i64 {
        self.to_rust().valence_count()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::NotAromatic() => ("NotAromatic", 0),
            Self::Aromatic(_) => ("Aromatic", 1),
        };
        variant_repr(slf.bind(py).as_any(), "AromaticValence", variant, arity)
    }
}

impl AromaticValence {
    pub(crate) fn from_rust(valence: GraphIrAromaticValence) -> Self {
        match valence {
            GraphIrAromaticValence::NotAromatic => Self::NotAromatic(),
            GraphIrAromaticValence::Aromatic(valence) => Self::Aromatic(valence),
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrAromaticValence {
        match self {
            Self::NotAromatic() => GraphIrAromaticValence::NotAromatic,
            Self::Aromatic(valence) => GraphIrAromaticValence::Aromatic(valence),
        }
    }
}

/// Aromatic-valence state: undetermined, explicitly not aromatic, or aromatic with
/// an aromatic-valence count. `Aromatic` coerces `int | NumForm` on construction.
#[pyclass]
pub enum AromaticValenceForm {
    Undetermined(),
    NotAromatic(),
    Aromatic(NumLike),
}

#[pymethods]
impl AromaticValenceForm {
    /// The exact aromatic-valence state, or `None` when this expression is not ground.
    pub(crate) fn as_lit(&self, py: Python<'_>) -> Option<AromaticValence> {
        self.to_rust(py).as_lit().map(AromaticValence::from_rust)
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::Undetermined() => ("Undetermined", 0),
            Self::NotAromatic() => ("NotAromatic", 0),
            Self::Aromatic(_) => ("Aromatic", 1),
        };
        variant_repr(slf.bind(py).as_any(), "AromaticValenceForm", variant, arity)
    }
}

impl_py_lattice!(
    AromaticValenceForm,
    GraphIrAromaticValenceForm,
    |value: &AromaticValenceForm, py: Python<'_>| -> PyResult<GraphIrAromaticValenceForm> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrAromaticValenceForm| -> PyResult<AromaticValenceForm> {
        AromaticValenceForm::from_rust(py, &value)
    }
);

impl AromaticValenceForm {
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrAromaticValenceForm) -> PyResult<Self> {
        Ok(match ast {
            GraphIrAromaticValenceForm::Undetermined => Self::Undetermined(),
            GraphIrAromaticValenceForm::NotAromatic => Self::NotAromatic(),
            GraphIrAromaticValenceForm::Aromatic(v) => Self::Aromatic(NumLike::Ast(
                into_py_variant(py, NumForm::from_rust(py, v)?)?,
            )),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrAromaticValenceForm {
        match self {
            Self::Undetermined() => GraphIrAromaticValenceForm::Undetermined,
            Self::NotAromatic() => GraphIrAromaticValenceForm::NotAromatic,
            Self::Aromatic(v) => GraphIrAromaticValenceForm::Aromatic(v.to_rust(py)),
        }
    }
}

#[derive(FromPyObject)]
pub(crate) enum AromaticValenceLike {
    Flag(bool),
    Value(NumLike),
    Ast(Py<AromaticValenceForm>),
}

impl AromaticValenceLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<GraphIrAromaticValenceForm> {
        Ok(match self {
            Self::Flag(false) => GraphIrAromaticValenceForm::NotAromatic,
            Self::Flag(true) => {
                return Err(PyValueError::new_err(
                    "aromatic_valence = True is not meaningful; use an int count or False",
                ))
            }
            Self::Value(v) => GraphIrAromaticValenceForm::Aromatic(v.to_rust(py)),
            Self::Ast(a) => a.bind(py).borrow().to_rust(py),
        })
    }
}

/// Exact ground multicenter-valence state.
#[pyclass(from_py_object)]
#[derive(Clone, Copy)]
pub enum MulticenterValence {
    NotMulticenter(),
    Multicenter(i64),
}

#[pymethods]
impl MulticenterValence {
    /// Multicenter-valence count, identifying explicit non-multicenter participation with zero.
    fn valence_count(&self) -> i64 {
        self.to_rust().valence_count()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::NotMulticenter() => ("NotMulticenter", 0),
            Self::Multicenter(_) => ("Multicenter", 1),
        };
        variant_repr(slf.bind(py).as_any(), "MulticenterValence", variant, arity)
    }
}

impl MulticenterValence {
    pub(crate) fn from_rust(valence: GraphIrMulticenterValence) -> Self {
        match valence {
            GraphIrMulticenterValence::NotMulticenter => Self::NotMulticenter(),
            GraphIrMulticenterValence::Multicenter(valence) => Self::Multicenter(valence),
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrMulticenterValence {
        match self {
            Self::NotMulticenter() => GraphIrMulticenterValence::NotMulticenter,
            Self::Multicenter(valence) => GraphIrMulticenterValence::Multicenter(valence),
        }
    }
}

/// Multicenter-valence state: undetermined, explicitly not multicenter, or
/// multicenter with a multicenter-valence count.
#[pyclass]
pub enum MulticenterValenceForm {
    Undetermined(),
    NotMulticenter(),
    Multicenter(NumLike),
}

#[pymethods]
impl MulticenterValenceForm {
    /// The exact multicenter-valence state, or `None` when this expression is not ground.
    pub(crate) fn as_lit(&self, py: Python<'_>) -> Option<MulticenterValence> {
        self.to_rust(py).as_lit().map(MulticenterValence::from_rust)
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::Undetermined() => ("Undetermined", 0),
            Self::NotMulticenter() => ("NotMulticenter", 0),
            Self::Multicenter(_) => ("Multicenter", 1),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "MulticenterValenceForm",
            variant,
            arity,
        )
    }
}

impl MulticenterValenceForm {
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrMulticenterValenceForm) -> PyResult<Self> {
        Ok(match ast {
            GraphIrMulticenterValenceForm::Undetermined => Self::Undetermined(),
            GraphIrMulticenterValenceForm::NotMulticenter => Self::NotMulticenter(),
            GraphIrMulticenterValenceForm::Multicenter(v) => Self::Multicenter(NumLike::Ast(
                into_py_variant(py, NumForm::from_rust(py, v)?)?,
            )),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrMulticenterValenceForm {
        match self {
            Self::Undetermined() => GraphIrMulticenterValenceForm::Undetermined,
            Self::NotMulticenter() => GraphIrMulticenterValenceForm::NotMulticenter,
            Self::Multicenter(v) => GraphIrMulticenterValenceForm::Multicenter(v.to_rust(py)),
        }
    }
}

impl_py_lattice!(
    MulticenterValenceForm,
    GraphIrMulticenterValenceForm,
    |value: &MulticenterValenceForm, py: Python<'_>| -> PyResult<GraphIrMulticenterValenceForm> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrMulticenterValenceForm| -> PyResult<MulticenterValenceForm> {
        MulticenterValenceForm::from_rust(py, &value)
    }
);

#[derive(FromPyObject)]
pub(crate) enum MulticenterValenceLike {
    Flag(bool),
    Value(NumLike),
    Ast(Py<MulticenterValenceForm>),
}

impl MulticenterValenceLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<GraphIrMulticenterValenceForm> {
        Ok(match self {
            Self::Flag(false) => GraphIrMulticenterValenceForm::NotMulticenter,
            Self::Flag(true) => {
                return Err(PyValueError::new_err(
                    "multicenter_valence = True is not meaningful; use an int count or False",
                ))
            }
            Self::Value(v) => GraphIrMulticenterValenceForm::Multicenter(v.to_rust(py)),
            Self::Ast(a) => a.bind(py).borrow().to_rust(py),
        })
    }
}

#[derive(FromPyObject)]
pub(crate) enum TetrahedralStereoLike {
    Flag(bool),
    Config(TetrahedralConfiguration),
    Ast(Py<TetrahedralStereoForm>),
}

impl TetrahedralStereoLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<GraphIrTetrahedralStereoForm> {
        Ok(match self {
            Self::Flag(false) => GraphIrTetrahedralStereoForm::NotStereo,
            Self::Flag(true) => {
                return Err(PyValueError::new_err(
                    "tetrahedral_stereo = True is not meaningful; use TetrahedralConfiguration.Ccw/Cw or False",
                ))
            }
            Self::Config(ts) => ts.to_rust().into(),
            Self::Ast(a) => a.bind(py).borrow().to_rust(py),
        })
    }
}

/// The key (identity) of an atom constraint, for keyed lookup. The ring-membership
/// key carries its ring scope; all other keys are the bare discriminant.
#[pyclass]
pub enum AtomConstraintKey {
    Valence(),
    DonatedPairs(),
    AcceptedPairs(),
    AromaticValence(),
    MulticenterValence(),
    TetrahedralStereo(),
    Degree(),
    TotalDegree(),
    TotalValence(),
    RingDegree(),
    RingValence(),
    TotalHydrogens(),
    RingMembership(Py<RingScope>),
}

#[pymethods]
impl AtomConstraintKey {
    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            AtomConstraintKey::Valence() => ("Valence", 0),
            AtomConstraintKey::DonatedPairs() => ("DonatedPairs", 0),
            AtomConstraintKey::AcceptedPairs() => ("AcceptedPairs", 0),
            AtomConstraintKey::AromaticValence() => ("AromaticValence", 0),
            AtomConstraintKey::MulticenterValence() => ("MulticenterValence", 0),
            AtomConstraintKey::TetrahedralStereo() => ("TetrahedralStereo", 0),
            AtomConstraintKey::Degree() => ("Degree", 0),
            AtomConstraintKey::TotalDegree() => ("TotalDegree", 0),
            AtomConstraintKey::TotalValence() => ("TotalValence", 0),
            AtomConstraintKey::RingDegree() => ("RingDegree", 0),
            AtomConstraintKey::RingValence() => ("RingValence", 0),
            AtomConstraintKey::TotalHydrogens() => ("TotalHydrogens", 0),
            AtomConstraintKey::RingMembership(_) => ("RingMembership", 1),
        };
        variant_repr(slf.bind(py).as_any(), "AtomConstraintKey", variant, arity)
    }
}

impl AtomConstraintKey {
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrAtomConstraintKey) -> PyResult<Self> {
        Ok(match ast {
            GraphIrAtomConstraintKey::Valence => Self::Valence(),
            GraphIrAtomConstraintKey::DonatedPairs => Self::DonatedPairs(),
            GraphIrAtomConstraintKey::AcceptedPairs => Self::AcceptedPairs(),
            GraphIrAtomConstraintKey::AromaticValence => Self::AromaticValence(),
            GraphIrAtomConstraintKey::MulticenterValence => Self::MulticenterValence(),
            GraphIrAtomConstraintKey::TetrahedralStereo => Self::TetrahedralStereo(),
            GraphIrAtomConstraintKey::Degree => Self::Degree(),
            GraphIrAtomConstraintKey::TotalDegree => Self::TotalDegree(),
            GraphIrAtomConstraintKey::TotalValence => Self::TotalValence(),
            GraphIrAtomConstraintKey::RingDegree => Self::RingDegree(),
            GraphIrAtomConstraintKey::RingValence => Self::RingValence(),
            GraphIrAtomConstraintKey::TotalHydrogens => Self::TotalHydrogens(),
            GraphIrAtomConstraintKey::RingMembership(scope) => {
                Self::RingMembership(into_py_variant(py, RingScope::from_rust(scope))?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrAtomConstraintKey {
        match self {
            Self::Valence() => GraphIrAtomConstraintKey::Valence,
            Self::DonatedPairs() => GraphIrAtomConstraintKey::DonatedPairs,
            Self::AcceptedPairs() => GraphIrAtomConstraintKey::AcceptedPairs,
            Self::AromaticValence() => GraphIrAtomConstraintKey::AromaticValence,
            Self::MulticenterValence() => GraphIrAtomConstraintKey::MulticenterValence,
            Self::TetrahedralStereo() => GraphIrAtomConstraintKey::TetrahedralStereo,
            Self::Degree() => GraphIrAtomConstraintKey::Degree,
            Self::TotalDegree() => GraphIrAtomConstraintKey::TotalDegree,
            Self::TotalValence() => GraphIrAtomConstraintKey::TotalValence,
            Self::RingDegree() => GraphIrAtomConstraintKey::RingDegree,
            Self::RingValence() => GraphIrAtomConstraintKey::RingValence,
            Self::TotalHydrogens() => GraphIrAtomConstraintKey::TotalHydrogens,
            Self::RingMembership(scope) => {
                GraphIrAtomConstraintKey::RingMembership(scope.bind(py).borrow().to_rust())
            }
        }
    }
}

/// An atom-scope constraint: a predicate on a valence, degree, ring, or stereo
/// property of a single atom.
#[pyclass(frozen)]
pub enum AtomConstraintForm {
    Valence(Py<NumForm>),
    TotalValence(Py<NumForm>),
    AromaticValence(Py<AromaticValenceForm>),
    MulticenterValence(Py<MulticenterValenceForm>),
    DonatedPairs(Py<NumForm>),
    AcceptedPairs(Py<NumForm>),
    Degree(Py<NumForm>),
    TotalDegree(Py<NumForm>),
    RingDegree(Py<NumForm>),
    RingValence(Py<NumForm>),
    TotalHydrogens(Py<NumForm>),
    RingMembership(Py<RingMembershipForm>),
    TetrahedralStereo(Py<TetrahedralStereoForm>),
}

#[pymethods]
impl AtomConstraintForm {
    /// The constraint's key (identity).
    #[getter]
    pub(crate) fn key(&self, py: Python<'_>) -> PyResult<AtomConstraintKey> {
        AtomConstraintKey::from_rust(py, &self.to_rust(py).key())
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            AtomConstraintForm::Valence(_) => "Valence",
            AtomConstraintForm::TotalValence(_) => "TotalValence",
            AtomConstraintForm::AromaticValence(_) => "AromaticValence",
            AtomConstraintForm::MulticenterValence(_) => "MulticenterValence",
            AtomConstraintForm::DonatedPairs(_) => "DonatedPairs",
            AtomConstraintForm::AcceptedPairs(_) => "AcceptedPairs",
            AtomConstraintForm::Degree(_) => "Degree",
            AtomConstraintForm::TotalDegree(_) => "TotalDegree",
            AtomConstraintForm::RingDegree(_) => "RingDegree",
            AtomConstraintForm::RingValence(_) => "RingValence",
            AtomConstraintForm::TotalHydrogens(_) => "TotalHydrogens",
            AtomConstraintForm::RingMembership(_) => "RingMembership",
            AtomConstraintForm::TetrahedralStereo(_) => "TetrahedralStereo",
        };
        variant_repr(slf.bind(py).as_any(), "AtomConstraintForm", variant, 1)
    }
}

impl_py_lattice!(
    AtomConstraintForm,
    GraphIrAtomConstraintForm,
    |value: &AtomConstraintForm, py: Python<'_>| -> PyResult<GraphIrAtomConstraintForm> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrAtomConstraintForm| -> PyResult<AtomConstraintForm> {
        AtomConstraintForm::from_rust(py, &value)
    }
);

impl AtomConstraintForm {
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrAtomConstraintForm) -> PyResult<Self> {
        Ok(match ast {
            GraphIrAtomConstraintForm::Valence(v) => {
                Self::Valence(into_py_variant(py, NumForm::from_rust(py, v)?)?)
            }
            GraphIrAtomConstraintForm::TotalValence(v) => {
                Self::TotalValence(into_py_variant(py, NumForm::from_rust(py, v)?)?)
            }
            GraphIrAtomConstraintForm::AromaticValence(c) => {
                Self::AromaticValence(into_py_variant(py, AromaticValenceForm::from_rust(py, c)?)?)
            }
            GraphIrAtomConstraintForm::MulticenterValence(c) => Self::MulticenterValence(
                into_py_variant(py, MulticenterValenceForm::from_rust(py, c)?)?,
            ),
            GraphIrAtomConstraintForm::DonatedPairs(v) => {
                Self::DonatedPairs(into_py_variant(py, NumForm::from_rust(py, v)?)?)
            }
            GraphIrAtomConstraintForm::AcceptedPairs(v) => {
                Self::AcceptedPairs(into_py_variant(py, NumForm::from_rust(py, v)?)?)
            }
            GraphIrAtomConstraintForm::Degree(v) => {
                Self::Degree(into_py_variant(py, NumForm::from_rust(py, v)?)?)
            }
            GraphIrAtomConstraintForm::TotalDegree(v) => {
                Self::TotalDegree(into_py_variant(py, NumForm::from_rust(py, v)?)?)
            }
            GraphIrAtomConstraintForm::RingDegree(v) => {
                Self::RingDegree(into_py_variant(py, NumForm::from_rust(py, v)?)?)
            }
            GraphIrAtomConstraintForm::RingValence(v) => {
                Self::RingValence(into_py_variant(py, NumForm::from_rust(py, v)?)?)
            }
            GraphIrAtomConstraintForm::TotalHydrogens(v) => {
                Self::TotalHydrogens(into_py_variant(py, NumForm::from_rust(py, v)?)?)
            }
            GraphIrAtomConstraintForm::RingMembership(m) => {
                Self::RingMembership(into_py_variant(py, RingMembershipForm::from_rust(py, m)?)?)
            }
            GraphIrAtomConstraintForm::TetrahedralStereo(c) => Self::TetrahedralStereo(
                into_py_variant(py, TetrahedralStereoForm::from_rust(py, c)?)?,
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrAtomConstraintForm {
        match self {
            Self::Valence(v) => GraphIrAtomConstraintForm::Valence(v.bind(py).borrow().to_rust(py)),
            Self::TotalValence(v) => {
                GraphIrAtomConstraintForm::TotalValence(v.bind(py).borrow().to_rust(py))
            }
            Self::AromaticValence(c) => {
                GraphIrAtomConstraintForm::AromaticValence(c.bind(py).borrow().to_rust(py))
            }
            Self::MulticenterValence(c) => {
                GraphIrAtomConstraintForm::MulticenterValence(c.bind(py).borrow().to_rust(py))
            }
            Self::DonatedPairs(v) => {
                GraphIrAtomConstraintForm::DonatedPairs(v.bind(py).borrow().to_rust(py))
            }
            Self::AcceptedPairs(v) => {
                GraphIrAtomConstraintForm::AcceptedPairs(v.bind(py).borrow().to_rust(py))
            }
            Self::Degree(v) => GraphIrAtomConstraintForm::Degree(v.bind(py).borrow().to_rust(py)),
            Self::TotalDegree(v) => {
                GraphIrAtomConstraintForm::TotalDegree(v.bind(py).borrow().to_rust(py))
            }
            Self::RingDegree(v) => {
                GraphIrAtomConstraintForm::RingDegree(v.bind(py).borrow().to_rust(py))
            }
            Self::RingValence(v) => {
                GraphIrAtomConstraintForm::RingValence(v.bind(py).borrow().to_rust(py))
            }
            Self::TotalHydrogens(v) => {
                GraphIrAtomConstraintForm::TotalHydrogens(v.bind(py).borrow().to_rust(py))
            }
            Self::RingMembership(m) => {
                GraphIrAtomConstraintForm::RingMembership(m.bind(py).borrow().to_rust(py))
            }
            Self::TetrahedralStereo(c) => {
                GraphIrAtomConstraintForm::TetrahedralStereo(c.bind(py).borrow().to_rust(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container (value or live view) or
/// an iterable of `AtomConstraintForm` (each `set`, last-wins).
#[derive(FromPyObject)]
pub(crate) enum AtomConstraintsUpdate {
    Container(Py<AtomConstraintsForm>),
    View(Py<AtomConstraintsView>),
    Entries(Vec<Py<AtomConstraintForm>>),
}

impl AtomConstraintsUpdate {
    /// Read every Python object into owned data — no write target is touched. Callers
    /// resolve *before* taking the write borrow so a view (or container) that aliases the
    /// same atom is read while nothing is borrowed (otherwise
    /// `atom.constraints.update(atom.constraints)` self-aliases into a double-borrow panic).
    pub(crate) fn resolve(&self, py: Python<'_>) -> PyResult<ResolvedAtomConstraintsUpdate> {
        Ok(match self {
            AtomConstraintsUpdate::Container(c) => {
                ResolvedAtomConstraintsUpdate::Overlay(c.bind(py).borrow().inner().clone())
            }
            AtomConstraintsUpdate::View(v) => ResolvedAtomConstraintsUpdate::Overlay(
                v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?,
            ),
            AtomConstraintsUpdate::Entries(entries) => ResolvedAtomConstraintsUpdate::Entries(
                entries
                    .iter()
                    .map(|entry| entry.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
        })
    }
}

/// A `AtomConstraintsUpdate` with all Python-object reads already done, so it can be applied
/// under a write borrow without re-entering Python.
pub(crate) enum ResolvedAtomConstraintsUpdate {
    /// A whole container (from another container or a live view): overlaid via `update`
    /// (last-wins per key; undetermined entries remove).
    Overlay(GraphIrAtomConstraintsForm),
    /// Loose entries: `set` each (last-wins; undetermined entries stored, not removed).
    Entries(Vec<GraphIrAtomConstraintForm>),
}

impl ResolvedAtomConstraintsUpdate {
    /// Overlay onto `target` in place. No Python reads.
    pub(crate) fn apply(self, target: &mut GraphIrAtomConstraintsForm) {
        match self {
            ResolvedAtomConstraintsUpdate::Overlay(overlay) => target.update(&overlay),
            ResolvedAtomConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry);
                }
            }
        }
    }
}

/// A whole-container argument that snapshots either a value container or a live
/// view — for the atom `constraints` setter, which accepts either.
#[derive(FromPyObject)]
pub(crate) enum AtomConstraintsLike {
    Container(Py<AtomConstraintsForm>),
    View(Py<AtomConstraintsView>),
}

impl AtomConstraintsLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<GraphIrAtomConstraintsForm> {
        match self {
            AtomConstraintsLike::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
            AtomConstraintsLike::View(v) => v.bind(py).borrow().read(py, |cs| Ok(cs.clone())),
        }
    }
}

/// The atom-scope constraints on an atom, in kind-sorted order. Mutable, hence
/// value-equal but unhashable (matching `AtomForm`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct AtomConstraintsForm(GraphIrAtomConstraintsForm);

#[pymethods]
impl AtomConstraintsForm {
    /// Build from a sequence of constraints (kind-sorted; a unique kind replaces
    /// an earlier one, ring memberships accumulate per scope).
    #[new]
    pub(crate) fn new(py: Python<'_>, entries: Vec<Py<AtomConstraintForm>>) -> Self {
        let mut constraints = GraphIrAtomConstraintsForm::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py)),
        );
        AtomConstraintsForm(constraints)
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let value = into_py_variant(py, AtomConstraintForm::from_rust(py, entry)?)?;
            parts.push(value.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("AtomConstraintsForm([{}])", parts.join(", ")))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    pub(crate) fn set(&mut self, py: Python<'_>, c: Py<AtomConstraintForm>) {
        self.0.set(c.bind(py).borrow().to_rust(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    pub(crate) fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<Option<AtomConstraintForm>> {
        self.0
            .remove(key.bind(py).borrow().to_rust(py))
            .map(|c| AtomConstraintForm::from_rust(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container, a live view, or an
    /// iterable of `AtomConstraintForm` (last-wins per key; undetermined entries remove).
    /// Takes `slf` by handle so `other` is fully read *before* the write borrow —
    /// `cs.update(cs)` on the same container is then a no-op, not a double-borrow panic.
    pub(crate) fn update(
        slf: Py<Self>,
        py: Python<'_>,
        other: AtomConstraintsUpdate,
    ) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        resolved.apply(&mut slf.borrow_mut(py).0);
        Ok(())
    }

    pub(crate) fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<AtomConstraintKeyIter> {
        atom_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<AtomConstraintKeyIter> {
        atom_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<AtomConstraintIter> {
        atom_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<AtomConstraintItemsIter> {
        atom_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_rust(py)) {
            Some(constraint) => {
                Ok(into_py_variant(py, AtomConstraintForm::from_rust(py, constraint)?)?.into_any())
            }
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    pub(crate) fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<AtomConstraintForm> {
        match self.0.get(key.bind(py).borrow().to_rust(py)) {
            Some(constraint) => AtomConstraintForm::from_rust(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_rust(py)).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    pub(crate) fn __contains__(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_rust(py))
    }

    /// The valence value, or `None`.
    #[getter]
    pub(crate) fn valence(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .valence()
            .map(|v| NumForm::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_valence(&mut self, py: Python<'_>, value: NumLike) {
        self.0
            .set(GraphIrAtomConstraintForm::valence(value.to_rust(py)));
    }

    /// The donated-pairs value, or `None`.
    #[getter]
    pub(crate) fn donated_pairs(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .donated_pairs()
            .map(|v| NumForm::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_donated_pairs(&mut self, py: Python<'_>, value: NumLike) {
        self.0
            .set(GraphIrAtomConstraintForm::donated_pairs(value.to_rust(py)));
    }

    /// The accepted-pairs value, or `None`.
    #[getter]
    pub(crate) fn accepted_pairs(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .accepted_pairs()
            .map(|v| NumForm::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_accepted_pairs(&mut self, py: Python<'_>, value: NumLike) {
        self.0
            .set(GraphIrAtomConstraintForm::accepted_pairs(value.to_rust(py)));
    }

    /// The aromatic-valence state, or `None`.
    #[getter]
    pub(crate) fn aromatic_valence(&self, py: Python<'_>) -> PyResult<Option<AromaticValenceForm>> {
        self.0
            .aromatic_valence()
            .map(|c| AromaticValenceForm::from_rust(py, c))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_aromatic_valence(
        &mut self,
        py: Python<'_>,
        value: AromaticValenceLike,
    ) -> PyResult<()> {
        self.0.set(GraphIrAtomConstraintForm::aromatic_valence(
            value.to_rust(py)?,
        ));
        Ok(())
    }

    /// The multicenter-valence state, or `None`.
    #[getter]
    pub(crate) fn multicenter_valence(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<MulticenterValenceForm>> {
        self.0
            .multicenter_valence()
            .map(|c| MulticenterValenceForm::from_rust(py, c))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_multicenter_valence(
        &mut self,
        py: Python<'_>,
        value: MulticenterValenceLike,
    ) -> PyResult<()> {
        self.0.set(GraphIrAtomConstraintForm::multicenter_valence(
            value.to_rust(py)?,
        ));
        Ok(())
    }

    /// The tetrahedral-stereo state, or `None`.
    #[getter]
    pub(crate) fn tetrahedral_stereo(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<TetrahedralStereoForm>> {
        self.0
            .tetrahedral_stereo()
            .map(|c| TetrahedralStereoForm::from_rust(py, c))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_tetrahedral_stereo(
        &mut self,
        py: Python<'_>,
        value: TetrahedralStereoLike,
    ) -> PyResult<()> {
        self.0.set(GraphIrAtomConstraintForm::tetrahedral_stereo(
            value.to_rust(py)?,
        ));
        Ok(())
    }

    /// The degree value, or `None`.
    #[getter]
    pub(crate) fn degree(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .degree()
            .map(|v| NumForm::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_degree(&mut self, py: Python<'_>, value: NumLike) {
        self.0
            .set(GraphIrAtomConstraintForm::degree(value.to_rust(py)));
    }

    /// The total-degree value, or `None`.
    #[getter]
    pub(crate) fn total_degree(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .total_degree()
            .map(|v| NumForm::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_total_degree(&mut self, py: Python<'_>, value: NumLike) {
        self.0
            .set(GraphIrAtomConstraintForm::total_degree(value.to_rust(py)));
    }

    /// The total-valence value, or `None`.
    #[getter]
    pub(crate) fn total_valence(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .total_valence()
            .map(|v| NumForm::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_total_valence(&mut self, py: Python<'_>, value: NumLike) {
        self.0
            .set(GraphIrAtomConstraintForm::total_valence(value.to_rust(py)));
    }

    /// The ring-degree value, or `None`.
    #[getter]
    pub(crate) fn ring_degree(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .ring_degree()
            .map(|v| NumForm::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_ring_degree(&mut self, py: Python<'_>, value: NumLike) {
        self.0
            .set(GraphIrAtomConstraintForm::ring_degree(value.to_rust(py)));
    }

    /// The ring-valence value, or `None`.
    #[getter]
    pub(crate) fn ring_valence(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .ring_valence()
            .map(|v| NumForm::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_ring_valence(&mut self, py: Python<'_>, value: NumLike) {
        self.0
            .set(GraphIrAtomConstraintForm::ring_valence(value.to_rust(py)));
    }

    /// The total-hydrogens value, or `None`.
    #[getter]
    pub(crate) fn total_hydrogens(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .total_hydrogens()
            .map(|v| NumForm::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_total_hydrogens(&mut self, py: Python<'_>, value: NumLike) {
        self.0.set(GraphIrAtomConstraintForm::total_hydrogens(
            value.to_rust(py),
        ));
    }

    /// The all-rings membership count, or `None`.
    #[getter]
    pub(crate) fn ring_count(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .ring_count()
            .map(|v| NumForm::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_ring_count(&mut self, py: Python<'_>, value: NumLike) {
        self.0.set(GraphIrAtomConstraintForm::ring_membership(
            GraphIrRingScope::All,
            value.to_rust(py),
        ));
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    pub(crate) fn ring_size_count(slf: Py<Self>) -> AtomRingSizeCounts {
        AtomRingSizeCounts {
            backing: AtomRingSizeBacking::Value(slf),
        }
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// Python values. Ring memberships key by scope: `ring_count` for the
    /// all-rings scope, `ring_size_count_<n>` for a specific ring size.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        atom_constraints_asdict(py, &self.0)
    }
}

impl AtomConstraintsForm {
    /// The wrapped AST constraints — read access for atom construction.
    pub(crate) fn inner(&self) -> &GraphIrAtomConstraintsForm {
        &self.0
    }

    /// Mutable access to the wrapped AST constraints — for the value-backed proxy.
    pub(crate) fn inner_mut(&mut self) -> &mut GraphIrAtomConstraintsForm {
        &mut self.0
    }

    /// Wrap owned AST constraints.
    pub(crate) fn from_inner(constraints: GraphIrAtomConstraintsForm) -> Self {
        AtomConstraintsForm(constraints)
    }
}

impl_py_lattice!(
    AtomConstraintsForm,
    GraphIrAtomConstraintsForm,
    |value: &AtomConstraintsForm, _py: Python<'_>| -> PyResult<GraphIrAtomConstraintsForm> {
        Ok(value.inner().clone())
    },
    |_py: Python<'_>, value: GraphIrAtomConstraintsForm| -> PyResult<AtomConstraintsForm> {
        Ok(AtomConstraintsForm(value))
    }
);

/// Build the per-constraint iterator handle from a borrowed container.
pub(crate) fn atom_constraints_iter(
    py: Python<'_>,
    constraints: &GraphIrAtomConstraintsForm,
) -> PyResult<AtomConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| into_py_variant(py, AtomConstraintForm::from_rust(py, constraint)?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AtomConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
pub(crate) fn atom_constraint_keys(
    py: Python<'_>,
    constraints: &GraphIrAtomConstraintsForm,
) -> PyResult<AtomConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| into_py_variant(py, AtomConstraintKey::from_rust(py, &constraint.key())?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AtomConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
pub(crate) fn atom_constraint_items(
    py: Python<'_>,
    constraints: &GraphIrAtomConstraintsForm,
) -> PyResult<AtomConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(py, AtomConstraintKey::from_rust(py, &constraint.key())?)?,
                into_py_variant(py, AtomConstraintForm::from_rust(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AtomConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// Python values. Ring memberships key by scope: `ring_count` for the
/// all-rings scope, `ring_size_count_<n>` for a specific ring size.
pub(crate) fn atom_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &GraphIrAtomConstraintsForm,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            GraphIrAtomConstraintForm::Valence(v) => {
                dict.set_item("valence", NumForm::from_rust(py, v)?)?
            }
            GraphIrAtomConstraintForm::DonatedPairs(v) => {
                dict.set_item("donated_pairs", NumForm::from_rust(py, v)?)?
            }
            GraphIrAtomConstraintForm::AcceptedPairs(v) => {
                dict.set_item("accepted_pairs", NumForm::from_rust(py, v)?)?
            }
            GraphIrAtomConstraintForm::AromaticValence(c) => {
                dict.set_item("aromatic_valence", AromaticValenceForm::from_rust(py, c)?)?
            }
            GraphIrAtomConstraintForm::MulticenterValence(c) => dict.set_item(
                "multicenter_valence",
                MulticenterValenceForm::from_rust(py, c)?,
            )?,
            GraphIrAtomConstraintForm::TetrahedralStereo(c) => dict.set_item(
                "tetrahedral_stereo",
                TetrahedralStereoForm::from_rust(py, c)?,
            )?,
            GraphIrAtomConstraintForm::Degree(v) => {
                dict.set_item("degree", NumForm::from_rust(py, v)?)?
            }
            GraphIrAtomConstraintForm::TotalDegree(v) => {
                dict.set_item("total_degree", NumForm::from_rust(py, v)?)?
            }
            GraphIrAtomConstraintForm::TotalValence(v) => {
                dict.set_item("total_valence", NumForm::from_rust(py, v)?)?
            }
            GraphIrAtomConstraintForm::RingDegree(v) => {
                dict.set_item("ring_degree", NumForm::from_rust(py, v)?)?
            }
            GraphIrAtomConstraintForm::RingValence(v) => {
                dict.set_item("ring_valence", NumForm::from_rust(py, v)?)?
            }
            GraphIrAtomConstraintForm::TotalHydrogens(v) => {
                dict.set_item("total_hydrogens", NumForm::from_rust(py, v)?)?
            }
            GraphIrAtomConstraintForm::RingMembership(m) => {
                let key = match m.scope {
                    GraphIrRingScope::All => "ring_count".to_string(),
                    GraphIrRingScope::Size(size) => format!("ring_size_count_{size}"),
                };
                dict.set_item(key, NumForm::from_rust(py, &m.count)?)?
            }
        }
    }
    Ok(dict)
}

/// What an `AtomConstraintsView` writes through to: an atom within a molecule
/// (by index) or a standalone `AtomForm`.
pub(crate) enum AtomConstraintsBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: GraphIrAtomId,
    },
    Atom(Py<AtomForm>),
}

/// A live handle onto one atom's constraints, backed by either a molecule-atom or
/// a standalone `AtomForm`. Reads borrow the atom's constraints and read only the
/// item they need (no whole-container clone); mutators write through to the atom in
/// place, without a clone-and-writeback.
#[pyclass]
pub struct AtomConstraintsView {
    pub(crate) backing: AtomConstraintsBacking,
}

impl AtomConstraintsView {
    /// Borrow the backing atom's constraints and read one item through `f` — no clone.
    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&GraphIrAtomConstraintsForm) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            AtomConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .atoms()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("atom id out of range"))?;
                f(&view.attributes.constraints)
            }
            AtomConstraintsBacking::Atom(atom) => {
                let atom = atom.bind(py).borrow();
                f(&atom.inner().constraints)
            }
        }
    }

    /// Mutate the backing atom's constraints in place through `f`.
    pub(crate) fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut GraphIrAtomConstraintsForm) -> R,
    ) -> PyResult<R> {
        match &self.backing {
            AtomConstraintsBacking::Molecule { owner, id } => Ok(f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .atom_mut(*id)
                .attributes
                .constraints)),
            AtomConstraintsBacking::Atom(atom) => {
                Ok(f(&mut atom.borrow_mut(py).try_inner_mut()?.constraints))
            }
        }
    }

    /// Set one constraint on the backing atom in place (last-wins per key).
    pub(crate) fn set_ast(
        &self,
        py: Python<'_>,
        constraint: GraphIrAtomConstraintForm,
    ) -> PyResult<()> {
        self.with_mut(py, |cs| cs.set(constraint))
    }

    /// Remove one key from the backing atom in place, returning the removed entry.
    pub(crate) fn remove_ast(
        &self,
        py: Python<'_>,
        key: GraphIrAtomConstraintKey,
    ) -> PyResult<Option<GraphIrAtomConstraintForm>> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl AtomConstraintsView {
    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("AtomConstraintsView({count} entries)"))
    }

    /// Insert `c` on the atom in place, replacing any existing entry of the same
    /// key (last-wins).
    pub(crate) fn set(&self, py: Python<'_>, c: Py<AtomConstraintForm>) -> PyResult<()> {
        self.set_ast(py, c.bind(py).borrow().to_rust(py))
    }

    /// Remove the entry with the given key from the atom in place, returning it if
    /// present (dict `pop`).
    pub(crate) fn pop(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<Option<AtomConstraintForm>> {
        self.remove_ast(py, key.bind(py).borrow().to_rust(py))?
            .map(|c| AtomConstraintForm::from_rust(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> PyResult<()> {
        if self
            .remove_ast(py, key.bind(py).borrow().to_rust(py))?
            .is_some()
        {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    /// Overlay `other` onto the atom's constraints in place — another container, a live
    /// view, or an iterable of `AtomConstraintForm` (last-wins per key; undetermined
    /// entries remove). Resolves `other` to owned data *before* the write borrow, so a
    /// view aliasing the same atom is not a double-borrow panic.
    pub(crate) fn update(&self, py: Python<'_>, other: AtomConstraintsUpdate) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        self.with_mut(py, |cs| resolved.apply(cs))
    }

    pub(crate) fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<AtomConstraintKeyIter> {
        self.read(py, |cs| atom_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<AtomConstraintKeyIter> {
        self.read(py, |cs| atom_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<AtomConstraintIter> {
        self.read(py, |cs| atom_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<AtomConstraintItemsIter> {
        self.read(py, |cs| atom_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_rust(py);
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| AtomConstraintForm::from_rust(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(into_py_variant(py, constraint)?.into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    pub(crate) fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<AtomConstraintForm> {
        let ast_key = key.bind(py).borrow().to_rust(py);
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| AtomConstraintForm::from_rust(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    pub(crate) fn __contains__(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_rust(py);
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// The valence value, or `None`.
    #[getter]
    pub(crate) fn valence(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.valence().map(|v| NumForm::from_rust(py, v)).transpose()
        })
    }

    #[setter]
    pub(crate) fn set_valence(&self, py: Python<'_>, value: NumLike) -> PyResult<()> {
        self.set_ast(py, GraphIrAtomConstraintForm::valence(value.to_rust(py)))
    }

    /// The donated-pairs value, or `None`.
    #[getter]
    pub(crate) fn donated_pairs(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.donated_pairs()
                .map(|v| NumForm::from_rust(py, v))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_donated_pairs(&self, py: Python<'_>, value: NumLike) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::donated_pairs(value.to_rust(py)),
        )
    }

    /// The accepted-pairs value, or `None`.
    #[getter]
    pub(crate) fn accepted_pairs(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.accepted_pairs()
                .map(|v| NumForm::from_rust(py, v))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_accepted_pairs(&self, py: Python<'_>, value: NumLike) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::accepted_pairs(value.to_rust(py)),
        )
    }

    /// The aromatic-valence state, or `None`.
    #[getter]
    pub(crate) fn aromatic_valence(&self, py: Python<'_>) -> PyResult<Option<AromaticValenceForm>> {
        self.read(py, |cs| {
            cs.aromatic_valence()
                .map(|c| AromaticValenceForm::from_rust(py, c))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_aromatic_valence(
        &self,
        py: Python<'_>,
        value: AromaticValenceLike,
    ) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::aromatic_valence(value.to_rust(py)?),
        )
    }

    /// The multicenter-valence state, or `None`.
    #[getter]
    pub(crate) fn multicenter_valence(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<MulticenterValenceForm>> {
        self.read(py, |cs| {
            cs.multicenter_valence()
                .map(|c| MulticenterValenceForm::from_rust(py, c))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_multicenter_valence(
        &self,
        py: Python<'_>,
        value: MulticenterValenceLike,
    ) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::multicenter_valence(value.to_rust(py)?),
        )
    }

    /// The tetrahedral-stereo state, or `None`.
    #[getter]
    pub(crate) fn tetrahedral_stereo(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<TetrahedralStereoForm>> {
        self.read(py, |cs| {
            cs.tetrahedral_stereo()
                .map(|c| TetrahedralStereoForm::from_rust(py, c))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_tetrahedral_stereo(
        &self,
        py: Python<'_>,
        value: TetrahedralStereoLike,
    ) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::tetrahedral_stereo(value.to_rust(py)?),
        )
    }

    /// The degree value, or `None`.
    #[getter]
    pub(crate) fn degree(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.degree().map(|v| NumForm::from_rust(py, v)).transpose()
        })
    }

    #[setter]
    pub(crate) fn set_degree(&self, py: Python<'_>, value: NumLike) -> PyResult<()> {
        self.set_ast(py, GraphIrAtomConstraintForm::degree(value.to_rust(py)))
    }

    /// The total-degree value, or `None`.
    #[getter]
    pub(crate) fn total_degree(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.total_degree()
                .map(|v| NumForm::from_rust(py, v))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_total_degree(&self, py: Python<'_>, value: NumLike) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::total_degree(value.to_rust(py)),
        )
    }

    /// The total-valence value, or `None`.
    #[getter]
    pub(crate) fn total_valence(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.total_valence()
                .map(|v| NumForm::from_rust(py, v))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_total_valence(&self, py: Python<'_>, value: NumLike) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::total_valence(value.to_rust(py)),
        )
    }

    /// The ring-degree value, or `None`.
    #[getter]
    pub(crate) fn ring_degree(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.ring_degree()
                .map(|v| NumForm::from_rust(py, v))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_ring_degree(&self, py: Python<'_>, value: NumLike) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::ring_degree(value.to_rust(py)),
        )
    }

    /// The ring-valence value, or `None`.
    #[getter]
    pub(crate) fn ring_valence(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.ring_valence()
                .map(|v| NumForm::from_rust(py, v))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_ring_valence(&self, py: Python<'_>, value: NumLike) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::ring_valence(value.to_rust(py)),
        )
    }

    /// The total-hydrogens value, or `None`.
    #[getter]
    pub(crate) fn total_hydrogens(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.total_hydrogens()
                .map(|v| NumForm::from_rust(py, v))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_total_hydrogens(&self, py: Python<'_>, value: NumLike) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::total_hydrogens(value.to_rust(py)),
        )
    }

    /// The all-rings membership count, or `None`.
    #[getter]
    pub(crate) fn ring_count(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.ring_count()
                .map(|v| NumForm::from_rust(py, v))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_ring_count(&self, py: Python<'_>, value: NumLike) -> PyResult<()> {
        self.set_ast(
            py,
            GraphIrAtomConstraintForm::ring_membership(GraphIrRingScope::All, value.to_rust(py)),
        )
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    pub(crate) fn ring_size_count(&self, py: Python<'_>) -> AtomRingSizeCounts {
        let backing = match &self.backing {
            AtomConstraintsBacking::Molecule { owner, id } => AtomRingSizeBacking::Molecule {
                owner: owner.clone_ref(py),
                id: *id,
            },
            AtomConstraintsBacking::Atom(atom) => AtomRingSizeBacking::Atom(atom.clone_ref(py)),
        };
        AtomRingSizeCounts { backing }
    }

    /// The present constraints as a dict keyed by snake_case name.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| atom_constraints_asdict(py, cs))
    }
}

/// What a `AtomRingSizeCounts` proxy reads/writes through to: an atom within a molecule,
/// a standalone `AtomForm`, or a standalone `AtomConstraintsForm` value.
pub(crate) enum AtomRingSizeBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: GraphIrAtomId,
    },
    Atom(Py<AtomForm>),
    Value(Py<AtomConstraintsForm>),
}

/// A subscriptable proxy over the sized-ring membership counts of an atom, keyed by
/// ring size: `proxy[size]` reads, `proxy[size] = count` sets, `del proxy[size]`
/// removes. Backs onto whichever container produced it (dual-backing, like
/// `AtomConstraintsView`).
#[pyclass]
pub struct AtomRingSizeCounts {
    pub(crate) backing: AtomRingSizeBacking,
}

impl AtomRingSizeCounts {
    /// Borrow the backing constraints and read through `f` — no clone.
    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&GraphIrAtomConstraintsForm) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            AtomRingSizeBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .atoms()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("atom id out of range"))?;
                f(&view.attributes.constraints)
            }
            AtomRingSizeBacking::Atom(atom) => f(&atom.bind(py).borrow().inner().constraints),
            AtomRingSizeBacking::Value(value) => f(value.bind(py).borrow().inner()),
        }
    }

    /// Mutate the backing constraints in place through `f`.
    pub(crate) fn write(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut GraphIrAtomConstraintsForm),
    ) -> PyResult<()> {
        match &self.backing {
            AtomRingSizeBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .atom_mut(*id)
                .attributes
                .constraints),
            AtomRingSizeBacking::Atom(atom) => {
                f(&mut atom.borrow_mut(py).try_inner_mut()?.constraints)
            }
            AtomRingSizeBacking::Value(value) => f(value.borrow_mut(py).inner_mut()),
        }
        Ok(())
    }
}

#[pymethods]
impl AtomRingSizeCounts {
    /// The membership count for rings of `size`, or `None`.
    pub(crate) fn __getitem__(&self, py: Python<'_>, size: u8) -> PyResult<Option<NumForm>> {
        self.read(py, |cs| {
            cs.ring_size_count(size)
                .map(|v| NumForm::from_rust(py, v))
                .transpose()
        })
    }

    /// The number of distinct ring sizes with a membership constraint.
    pub(crate) fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(ring_sizes(cs).count()))
    }

    pub(crate) fn __contains__(&self, py: Python<'_>, size: u8) -> PyResult<bool> {
        self.read(py, |cs| Ok(cs.ring_size_count(size).is_some()))
    }

    /// Iterate the present ring sizes (as ints).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<AtomRingSizeIter> {
        let sizes = self.read(py, |cs| Ok(ring_sizes(cs).collect::<Vec<u8>>()))?;
        Ok(AtomRingSizeIter {
            sizes: sizes.into_iter(),
        })
    }

    /// Set the membership count for rings of `size` in place.
    pub(crate) fn __setitem__(&self, py: Python<'_>, size: u8, count: NumLike) -> PyResult<()> {
        let constraint = GraphIrAtomConstraintForm::ring_membership(
            GraphIrRingScope::Size(size),
            count.to_rust(py),
        );
        self.write(py, |cs| cs.set(constraint))
    }

    /// Remove the sized-ring membership for `size` in place.
    pub(crate) fn __delitem__(&self, py: Python<'_>, size: u8) -> PyResult<()> {
        self.write(py, |cs| {
            cs.remove(GraphIrAtomConstraintKey::RingMembership(
                GraphIrRingScope::Size(size),
            ));
        })
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |cs| {
            let mut parts = Vec::new();
            for entry in cs.iter() {
                if let GraphIrAtomConstraintForm::RingMembership(m) = entry {
                    if let GraphIrRingScope::Size(size) = m.scope {
                        let count = into_py_variant(py, NumForm::from_rust(py, &m.count)?)?;
                        parts.push(format!(
                            "{size}: {}",
                            count.bind(py).as_any().repr()?.extract::<String>()?
                        ));
                    }
                }
            }
            Ok(format!("AtomRingSizeCounts({{{}}})", parts.join(", ")))
        })
    }
}

/// The ring sizes with a membership constraint, in kind-sorted order.
pub(crate) fn ring_sizes(
    constraints: &GraphIrAtomConstraintsForm,
) -> impl Iterator<Item = u8> + '_ {
    constraints.iter().filter_map(|entry| match entry {
        GraphIrAtomConstraintForm::RingMembership(m) => match m.scope {
            GraphIrRingScope::Size(size) => Some(size),
            GraphIrRingScope::All => None,
        },
        _ => None,
    })
}

#[pyclass]
pub(crate) struct AtomRingSizeIter {
    sizes: IntoIter<u8>,
}

#[pymethods]
impl AtomRingSizeIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<u8> {
        self.sizes.next()
    }
}

#[pyclass]
pub(crate) struct AtomConstraintIter {
    entries: IntoIter<Py<AtomConstraintForm>>,
}

#[pymethods]
impl AtomConstraintIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<AtomConstraintForm>> {
        self.entries.next()
    }
}

#[pyclass]
pub(crate) struct AtomConstraintKeyIter {
    keys: IntoIter<Py<AtomConstraintKey>>,
}

#[pymethods]
impl AtomConstraintKeyIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<AtomConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
pub(crate) struct AtomConstraintItemsIter {
    items: IntoIter<(Py<AtomConstraintKey>, Py<AtomConstraintForm>)>,
}

#[pymethods]
impl AtomConstraintItemsIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<(Py<AtomConstraintKey>, Py<AtomConstraintForm>)> {
        self.items.next()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(GraphIrAromaticValenceForm::Undetermined)]
    #[case(GraphIrAromaticValenceForm::NotAromatic)]
    #[case(GraphIrAromaticValenceForm::aromatic(1))]
    pub(crate) fn test_aromatic_valence_form_roundtrip(#[case] ast: GraphIrAromaticValenceForm) {
        Python::attach(|py| {
            assert_eq!(
                AromaticValenceForm::from_rust(py, &ast)
                    .unwrap()
                    .to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(
        GraphIrAromaticValenceForm::NotAromatic,
        Some(GraphIrAromaticValence::NotAromatic)
    )]
    #[case(
        GraphIrAromaticValenceForm::aromatic(2),
        Some(GraphIrAromaticValence::Aromatic(2))
    )]
    #[case(GraphIrAromaticValenceForm::Undetermined, None)]
    pub(crate) fn test_aromatic_valence_form_as_lit(
        #[case] ast: GraphIrAromaticValenceForm,
        #[case] expected: Option<GraphIrAromaticValence>,
    ) {
        Python::attach(|py| {
            assert_eq!(
                AromaticValenceForm::from_rust(py, &ast)
                    .unwrap()
                    .as_lit(py)
                    .map(AromaticValence::to_rust),
                expected
            );
        });
    }

    #[rstest]
    #[case(GraphIrMulticenterValenceForm::Undetermined)]
    #[case(GraphIrMulticenterValenceForm::NotMulticenter)]
    #[case(GraphIrMulticenterValenceForm::multicenter(2))]
    pub(crate) fn test_multicenter_valence_form_roundtrip(
        #[case] ast: GraphIrMulticenterValenceForm,
    ) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterValenceForm::from_rust(py, &ast)
                    .unwrap()
                    .to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(
        GraphIrMulticenterValenceForm::NotMulticenter,
        Some(GraphIrMulticenterValence::NotMulticenter)
    )]
    #[case(
        GraphIrMulticenterValenceForm::multicenter(3),
        Some(GraphIrMulticenterValence::Multicenter(3))
    )]
    #[case(GraphIrMulticenterValenceForm::Undetermined, None)]
    pub(crate) fn test_multicenter_valence_form_as_lit(
        #[case] ast: GraphIrMulticenterValenceForm,
        #[case] expected: Option<GraphIrMulticenterValence>,
    ) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterValenceForm::from_rust(py, &ast)
                    .unwrap()
                    .as_lit(py)
                    .map(MulticenterValence::to_rust),
                expected
            );
        });
    }
}
