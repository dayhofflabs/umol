//! Atom-constraint sub-ASTs mirroring `umol_ast::ast::constraint` (S5a): the
//! aromatic/multicenter valence states, ring scope, and ring membership. The
//! `AtomConstraintAst` enum and `AtomConstraintsAst` container follow at S5b.

use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use umol_ast::ast::{
    AromaticValenceAst as AstAromaticValenceAst, AtomConstraintAst as AstAtomConstraintAst,
    AtomConstraintKey as AstAtomConstraintKey, AtomConstraintsAst as AstAtomConstraintsAst,
    AtomId as AstAtomId, MulticenterValenceAst as AstMulticenterValenceAst,
    RingMembershipAst as AstRingMembershipAst, RingScope as AstRingScope,
    TetrahedralStereoAst as AstTetrahedralStereoAst,
};

use crate::atom::AtomAst;
use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::molecule::MoleculeAst;
use crate::stereo::{TetrahedralStereo, TetrahedralStereoAst};
use crate::value::{ValueArg, ValueAst};

/// Aromatic-valence state: undetermined, explicitly not aromatic, or aromatic with
/// an aromatic-valence count. `Aromatic` coerces `int | ValueAst` on construction.
#[pyclass]
pub enum AromaticValenceAst {
    Undetermined(),
    NotAromatic(),
    Aromatic(ValueArg),
}

#[pymethods]
impl AromaticValenceAst {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            AromaticValenceAst::Undetermined() => ("Undetermined", 0),
            AromaticValenceAst::NotAromatic() => ("NotAromatic", 0),
            AromaticValenceAst::Aromatic(_) => ("Aromatic", 1),
        };
        variant_repr(slf.bind(py).as_any(), "AromaticValenceAst", variant, arity)
    }
}

impl AromaticValenceAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAromaticValenceAst) -> PyResult<Self> {
        Ok(match ast {
            AstAromaticValenceAst::Undetermined => Self::Undetermined(),
            AstAromaticValenceAst::NotAromatic => Self::NotAromatic(),
            AstAromaticValenceAst::Aromatic(v) => {
                Self::Aromatic(ValueArg::Ast(into_py_variant(py, ValueAst::from_ast(py, v)?)?))
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAromaticValenceAst {
        match self {
            Self::Undetermined() => AstAromaticValenceAst::Undetermined,
            Self::NotAromatic() => AstAromaticValenceAst::NotAromatic,
            Self::Aromatic(v) => AstAromaticValenceAst::Aromatic(v.to_ast(py)),
        }
    }
}

/// Setter coercion for `aromatic_valence`: `False` → not aromatic, an `int`/`ValueAst`
/// → aromatic with that valence, or an `AromaticValenceAst` passthrough.
#[derive(FromPyObject)]
enum AromaticValenceArg {
    Flag(bool),
    Value(ValueArg),
    Ast(Py<AromaticValenceAst>),
}

impl AromaticValenceArg {
    fn to_ast(&self, py: Python<'_>) -> PyResult<AstAromaticValenceAst> {
        Ok(match self {
            AromaticValenceArg::Flag(false) => AstAromaticValenceAst::NotAromatic,
            AromaticValenceArg::Flag(true) => {
                return Err(PyValueError::new_err(
                    "aromatic_valence = True is not meaningful; use an int count or False",
                ))
            }
            AromaticValenceArg::Value(v) => AstAromaticValenceAst::Aromatic(v.to_ast(py)),
            AromaticValenceArg::Ast(a) => a.bind(py).borrow().to_ast(py),
        })
    }
}

/// Multicenter-valence state: undetermined, explicitly not multicenter, or
/// multicenter with a multicenter-valence count. `Multicenter` coerces
/// `int | ValueAst` on construction.
#[pyclass]
pub enum MulticenterValenceAst {
    Undetermined(),
    NotMulticenter(),
    Multicenter(ValueArg),
}

#[pymethods]
impl MulticenterValenceAst {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            MulticenterValenceAst::Undetermined() => ("Undetermined", 0),
            MulticenterValenceAst::NotMulticenter() => ("NotMulticenter", 0),
            MulticenterValenceAst::Multicenter(_) => ("Multicenter", 1),
        };
        variant_repr(slf.bind(py).as_any(), "MulticenterValenceAst", variant, arity)
    }
}

impl MulticenterValenceAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstMulticenterValenceAst) -> PyResult<Self> {
        Ok(match ast {
            AstMulticenterValenceAst::Undetermined => Self::Undetermined(),
            AstMulticenterValenceAst::NotMulticenter => Self::NotMulticenter(),
            AstMulticenterValenceAst::Multicenter(v) => {
                Self::Multicenter(ValueArg::Ast(into_py_variant(py, ValueAst::from_ast(py, v)?)?))
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstMulticenterValenceAst {
        match self {
            Self::Undetermined() => AstMulticenterValenceAst::Undetermined,
            Self::NotMulticenter() => AstMulticenterValenceAst::NotMulticenter,
            Self::Multicenter(v) => AstMulticenterValenceAst::Multicenter(v.to_ast(py)),
        }
    }
}

/// Setter coercion for `multicenter_valence`: `False` → not multicenter, an
/// `int`/`ValueAst` → multicenter with that valence, or a `MulticenterValenceAst`
/// passthrough.
#[derive(FromPyObject)]
enum MulticenterValenceArg {
    Flag(bool),
    Value(ValueArg),
    Ast(Py<MulticenterValenceAst>),
}

impl MulticenterValenceArg {
    fn to_ast(&self, py: Python<'_>) -> PyResult<AstMulticenterValenceAst> {
        Ok(match self {
            MulticenterValenceArg::Flag(false) => AstMulticenterValenceAst::NotMulticenter,
            MulticenterValenceArg::Flag(true) => {
                return Err(PyValueError::new_err(
                    "multicenter_valence = True is not meaningful; use an int count or False",
                ))
            }
            MulticenterValenceArg::Value(v) => AstMulticenterValenceAst::Multicenter(v.to_ast(py)),
            MulticenterValenceArg::Ast(a) => a.bind(py).borrow().to_ast(py),
        })
    }
}

/// Setter coercion for `tetrahedral_stereo`: `False` → not stereogenic, a
/// `TetrahedralStereo` (`Ccw`/`Cw`) → that coset, or a `TetrahedralStereoAst`
/// passthrough.
#[derive(FromPyObject)]
enum TetrahedralStereoArg {
    Flag(bool),
    Config(TetrahedralStereo),
    Ast(Py<TetrahedralStereoAst>),
}

impl TetrahedralStereoArg {
    fn to_ast(&self, py: Python<'_>) -> PyResult<AstTetrahedralStereoAst> {
        Ok(match self {
            TetrahedralStereoArg::Flag(false) => AstTetrahedralStereoAst::NotStereo,
            TetrahedralStereoArg::Flag(true) => {
                return Err(PyValueError::new_err(
                    "tetrahedral_stereo = True is not meaningful; use TetrahedralStereo.Ccw/Cw or False",
                ))
            }
            TetrahedralStereoArg::Config(ts) => ts.to_ast(),
            TetrahedralStereoArg::Ast(a) => a.bind(py).borrow().to_ast(py),
        })
    }
}

/// Ring scope: all rings, or rings of a given size.
#[pyclass]
pub enum RingScope {
    All(),
    Size(u8),
}

#[pymethods]
impl RingScope {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            RingScope::All() => ("All", 0),
            RingScope::Size(_) => ("Size", 1),
        };
        variant_repr(slf.bind(py).as_any(), "RingScope", variant, arity)
    }
}

impl RingScope {
    pub(crate) fn from_ast(ast: &AstRingScope) -> Self {
        match ast {
            AstRingScope::All => Self::All(),
            AstRingScope::Size(size) => Self::Size(*size),
        }
    }

    pub(crate) fn to_ast(&self) -> AstRingScope {
        match self {
            Self::All() => AstRingScope::All,
            Self::Size(size) => AstRingScope::Size(*size),
        }
    }
}

/// Ring-membership fact: a ring scope and a membership count.
#[pyclass]
pub struct RingMembershipAst {
    #[pyo3(get)]
    scope: Py<RingScope>,
    #[pyo3(get)]
    count: Py<ValueAst>,
}

#[pymethods]
impl RingMembershipAst {
    #[new]
    fn new(py: Python<'_>, scope: Py<RingScope>, count: ValueArg) -> PyResult<Self> {
        Ok(RingMembershipAst {
            scope,
            count: count.to_py(py)?,
        })
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "RingMembershipAst({}, {})",
            self.scope.bind(py).as_any().repr()?.extract::<String>()?,
            self.count.bind(py).as_any().repr()?.extract::<String>()?,
        ))
    }
}

impl RingMembershipAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstRingMembershipAst) -> PyResult<Self> {
        Ok(RingMembershipAst {
            scope: into_py_variant(py, RingScope::from_ast(&ast.scope))?,
            count: into_py_variant(py, ValueAst::from_ast(py, &ast.count)?)?,
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstRingMembershipAst {
        AstRingMembershipAst::new(
            self.scope.bind(py).borrow().to_ast(),
            self.count.bind(py).borrow().to_ast(py),
        )
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
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
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
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAtomConstraintKey) -> PyResult<Self> {
        Ok(match ast {
            AstAtomConstraintKey::Valence => Self::Valence(),
            AstAtomConstraintKey::DonatedPairs => Self::DonatedPairs(),
            AstAtomConstraintKey::AcceptedPairs => Self::AcceptedPairs(),
            AstAtomConstraintKey::AromaticValence => Self::AromaticValence(),
            AstAtomConstraintKey::MulticenterValence => Self::MulticenterValence(),
            AstAtomConstraintKey::TetrahedralStereo => Self::TetrahedralStereo(),
            AstAtomConstraintKey::Degree => Self::Degree(),
            AstAtomConstraintKey::TotalDegree => Self::TotalDegree(),
            AstAtomConstraintKey::TotalValence => Self::TotalValence(),
            AstAtomConstraintKey::RingDegree => Self::RingDegree(),
            AstAtomConstraintKey::RingValence => Self::RingValence(),
            AstAtomConstraintKey::TotalHydrogens => Self::TotalHydrogens(),
            AstAtomConstraintKey::RingMembership(scope) => {
                Self::RingMembership(into_py_variant(py, RingScope::from_ast(scope))?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAtomConstraintKey {
        match self {
            Self::Valence() => AstAtomConstraintKey::Valence,
            Self::DonatedPairs() => AstAtomConstraintKey::DonatedPairs,
            Self::AcceptedPairs() => AstAtomConstraintKey::AcceptedPairs,
            Self::AromaticValence() => AstAtomConstraintKey::AromaticValence,
            Self::MulticenterValence() => AstAtomConstraintKey::MulticenterValence,
            Self::TetrahedralStereo() => AstAtomConstraintKey::TetrahedralStereo,
            Self::Degree() => AstAtomConstraintKey::Degree,
            Self::TotalDegree() => AstAtomConstraintKey::TotalDegree,
            Self::TotalValence() => AstAtomConstraintKey::TotalValence,
            Self::RingDegree() => AstAtomConstraintKey::RingDegree,
            Self::RingValence() => AstAtomConstraintKey::RingValence,
            Self::TotalHydrogens() => AstAtomConstraintKey::TotalHydrogens,
            Self::RingMembership(scope) => {
                AstAtomConstraintKey::RingMembership(scope.bind(py).borrow().to_ast())
            }
        }
    }
}

/// An atom-scope constraint: a predicate on a valence, degree, ring, or stereo
/// property of a single atom.
#[pyclass]
pub enum AtomConstraintAst {
    Valence(Py<ValueAst>),
    TotalValence(Py<ValueAst>),
    AromaticValence(Py<AromaticValenceAst>),
    MulticenterValence(Py<MulticenterValenceAst>),
    DonatedPairs(Py<ValueAst>),
    AcceptedPairs(Py<ValueAst>),
    Degree(Py<ValueAst>),
    TotalDegree(Py<ValueAst>),
    RingDegree(Py<ValueAst>),
    RingValence(Py<ValueAst>),
    TotalHydrogens(Py<ValueAst>),
    RingMembership(Py<RingMembershipAst>),
    TetrahedralStereo(Py<TetrahedralStereoAst>),
}

#[pymethods]
impl AtomConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    fn key(&self, py: Python<'_>) -> PyResult<AtomConstraintKey> {
        AtomConstraintKey::from_ast(py, &self.to_ast(py).key())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            AtomConstraintAst::Valence(_) => "Valence",
            AtomConstraintAst::TotalValence(_) => "TotalValence",
            AtomConstraintAst::AromaticValence(_) => "AromaticValence",
            AtomConstraintAst::MulticenterValence(_) => "MulticenterValence",
            AtomConstraintAst::DonatedPairs(_) => "DonatedPairs",
            AtomConstraintAst::AcceptedPairs(_) => "AcceptedPairs",
            AtomConstraintAst::Degree(_) => "Degree",
            AtomConstraintAst::TotalDegree(_) => "TotalDegree",
            AtomConstraintAst::RingDegree(_) => "RingDegree",
            AtomConstraintAst::RingValence(_) => "RingValence",
            AtomConstraintAst::TotalHydrogens(_) => "TotalHydrogens",
            AtomConstraintAst::RingMembership(_) => "RingMembership",
            AtomConstraintAst::TetrahedralStereo(_) => "TetrahedralStereo",
        };
        variant_repr(slf.bind(py).as_any(), "AtomConstraintAst", variant, 1)
    }
}

impl AtomConstraintAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAtomConstraintAst) -> PyResult<Self> {
        Ok(match ast {
            AstAtomConstraintAst::Valence(v) => {
                Self::Valence(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::TotalValence(v) => {
                Self::TotalValence(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::AromaticValence(c) => {
                Self::AromaticValence(into_py_variant(py, AromaticValenceAst::from_ast(py, c)?)?)
            }
            AstAtomConstraintAst::MulticenterValence(c) => Self::MulticenterValence(
                into_py_variant(py, MulticenterValenceAst::from_ast(py, c)?)?,
            ),
            AstAtomConstraintAst::DonatedPairs(v) => {
                Self::DonatedPairs(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::AcceptedPairs(v) => {
                Self::AcceptedPairs(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::Degree(v) => {
                Self::Degree(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::TotalDegree(v) => {
                Self::TotalDegree(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::RingDegree(v) => {
                Self::RingDegree(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::RingValence(v) => {
                Self::RingValence(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::TotalHydrogens(v) => {
                Self::TotalHydrogens(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::RingMembership(m) => {
                Self::RingMembership(into_py_variant(py, RingMembershipAst::from_ast(py, m)?)?)
            }
            AstAtomConstraintAst::TetrahedralStereo(c) => Self::TetrahedralStereo(into_py_variant(
                py,
                TetrahedralStereoAst::from_ast(py, c)?,
            )?),
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAtomConstraintAst {
        match self {
            Self::Valence(v) => AstAtomConstraintAst::Valence(v.bind(py).borrow().to_ast(py)),
            Self::TotalValence(v) => {
                AstAtomConstraintAst::TotalValence(v.bind(py).borrow().to_ast(py))
            }
            Self::AromaticValence(c) => {
                AstAtomConstraintAst::AromaticValence(c.bind(py).borrow().to_ast(py))
            }
            Self::MulticenterValence(c) => {
                AstAtomConstraintAst::MulticenterValence(c.bind(py).borrow().to_ast(py))
            }
            Self::DonatedPairs(v) => {
                AstAtomConstraintAst::DonatedPairs(v.bind(py).borrow().to_ast(py))
            }
            Self::AcceptedPairs(v) => {
                AstAtomConstraintAst::AcceptedPairs(v.bind(py).borrow().to_ast(py))
            }
            Self::Degree(v) => AstAtomConstraintAst::Degree(v.bind(py).borrow().to_ast(py)),
            Self::TotalDegree(v) => {
                AstAtomConstraintAst::TotalDegree(v.bind(py).borrow().to_ast(py))
            }
            Self::RingDegree(v) => AstAtomConstraintAst::RingDegree(v.bind(py).borrow().to_ast(py)),
            Self::RingValence(v) => {
                AstAtomConstraintAst::RingValence(v.bind(py).borrow().to_ast(py))
            }
            Self::TotalHydrogens(v) => {
                AstAtomConstraintAst::TotalHydrogens(v.bind(py).borrow().to_ast(py))
            }
            Self::RingMembership(m) => {
                AstAtomConstraintAst::RingMembership(m.bind(py).borrow().to_ast(py))
            }
            Self::TetrahedralStereo(c) => {
                AstAtomConstraintAst::TetrahedralStereo(c.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container (value or live view) or
/// an iterable of `AtomConstraintAst` (each `set`, last-wins).
#[derive(FromPyObject)]
enum ConstraintsUpdate {
    Container(Py<AtomConstraintsAst>),
    View(Py<AtomConstraintsView>),
    Entries(Vec<Py<AtomConstraintAst>>),
}

impl ConstraintsUpdate {
    /// Overlay this update onto `target` in place.
    fn apply(&self, py: Python<'_>, target: &mut AstAtomConstraintsAst) -> PyResult<()> {
        match self {
            ConstraintsUpdate::Container(c) => target.update(c.bind(py).borrow().inner()),
            ConstraintsUpdate::View(v) => {
                let snapshot = v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?;
                target.update(&snapshot);
            }
            ConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry.bind(py).borrow().to_ast(py));
                }
            }
        }
        Ok(())
    }
}

/// The atom-scope constraints on an atom, in kind-sorted order.
#[pyclass]
pub struct AtomConstraintsAst(AstAtomConstraintsAst);

#[pymethods]
impl AtomConstraintsAst {
    /// Build from a sequence of constraints (kind-sorted; a unique kind replaces
    /// an earlier one, ring memberships accumulate per scope).
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<AtomConstraintAst>>) -> Self {
        let mut constraints = AstAtomConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_ast(py)),
        );
        AtomConstraintsAst(constraints)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.0)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let mirror = into_py_variant(py, AtomConstraintAst::from_ast(py, entry)?)?;
            parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("AtomConstraintsAst([{}])", parts.join(", ")))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    fn set(&mut self, py: Python<'_>, c: Py<AtomConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key, returning it if present.
    fn remove(
        &mut self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<Option<AtomConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_ast(py))
            .map(|c| AtomConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container or an iterable of
    /// `AtomConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&mut self, py: Python<'_>, other: ConstraintsUpdate) -> PyResult<()> {
        other.apply(py, &mut self.0)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __iter__(&self, py: Python<'_>) -> PyResult<AtomConstraintIter> {
        atom_constraints_iter(py, &self.0)
    }

    /// The constraint with the given key, or `None`.
    fn get(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<Option<AtomConstraintAst>> {
        self.0
            .get(key.bind(py).borrow().to_ast(py))
            .map(|constraint| AtomConstraintAst::from_ast(py, constraint))
            .transpose()
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> PyResult<AtomConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_ast(py)) {
            Some(constraint) => AtomConstraintAst::from_ast(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&mut self, py: Python<'_>, key: Py<AtomConstraintKey>) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_ast(py)).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast(py))
    }

    /// The valence value, or `None`.
    #[getter]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .valence()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_valence(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstAtomConstraintAst::valence(value.to_ast(py)));
    }

    /// The donated-pairs value, or `None`.
    #[getter]
    fn donated_pairs(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .donated_pairs()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_donated_pairs(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::donated_pairs(value.to_ast(py)));
    }

    /// The accepted-pairs value, or `None`.
    #[getter]
    fn accepted_pairs(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .accepted_pairs()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_accepted_pairs(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::accepted_pairs(value.to_ast(py)));
    }

    /// The aromatic-valence state, or `None`.
    #[getter]
    fn aromatic_valence(&self, py: Python<'_>) -> PyResult<Option<AromaticValenceAst>> {
        self.0
            .aromatic_valence()
            .map(|c| AromaticValenceAst::from_ast(py, c))
            .transpose()
    }

    #[setter]
    fn set_aromatic_valence(&mut self, py: Python<'_>, value: AromaticValenceArg) -> PyResult<()> {
        self.0
            .set(AstAtomConstraintAst::aromatic_valence(value.to_ast(py)?));
        Ok(())
    }

    /// The multicenter-valence state, or `None`.
    #[getter]
    fn multicenter_valence(&self, py: Python<'_>) -> PyResult<Option<MulticenterValenceAst>> {
        self.0
            .multicenter_valence()
            .map(|c| MulticenterValenceAst::from_ast(py, c))
            .transpose()
    }

    #[setter]
    fn set_multicenter_valence(
        &mut self,
        py: Python<'_>,
        value: MulticenterValenceArg,
    ) -> PyResult<()> {
        self.0
            .set(AstAtomConstraintAst::multicenter_valence(value.to_ast(py)?));
        Ok(())
    }

    /// The tetrahedral-stereo state, or `None`.
    #[getter]
    fn tetrahedral_stereo(&self, py: Python<'_>) -> PyResult<Option<TetrahedralStereoAst>> {
        self.0
            .tetrahedral_stereo()
            .map(|c| TetrahedralStereoAst::from_ast(py, c))
            .transpose()
    }

    #[setter]
    fn set_tetrahedral_stereo(
        &mut self,
        py: Python<'_>,
        value: TetrahedralStereoArg,
    ) -> PyResult<()> {
        self.0
            .set(AstAtomConstraintAst::tetrahedral_stereo(value.to_ast(py)?));
        Ok(())
    }

    /// The degree value, or `None`.
    #[getter]
    fn degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .degree()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_degree(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstAtomConstraintAst::degree(value.to_ast(py)));
    }

    /// The total-degree value, or `None`.
    #[getter]
    fn total_degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .total_degree()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_total_degree(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::total_degree(value.to_ast(py)));
    }

    /// The total-valence value, or `None`.
    #[getter]
    fn total_valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .total_valence()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_total_valence(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::total_valence(value.to_ast(py)));
    }

    /// The ring-degree value, or `None`.
    #[getter]
    fn ring_degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_degree()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_ring_degree(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::ring_degree(value.to_ast(py)));
    }

    /// The ring-valence value, or `None`.
    #[getter]
    fn ring_valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_valence()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_ring_valence(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::ring_valence(value.to_ast(py)));
    }

    /// The total-hydrogens value, or `None`.
    #[getter]
    fn total_hydrogens(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .total_hydrogens()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_total_hydrogens(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::total_hydrogens(value.to_ast(py)));
    }

    /// The all-rings membership count, or `None`.
    #[getter]
    fn ring_count(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_count()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_ring_count(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstAtomConstraintAst::ring_membership(
            AstRingScope::All,
            value.to_ast(py),
        ));
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    fn ring_size_count(slf: Py<Self>) -> RingSizeCounts {
        RingSizeCounts {
            backing: RingSizeBacking::Value(slf),
        }
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors. Ring memberships key by scope: `ring_count` for the
    /// all-rings scope, `ring_size_count_<n>` for a specific ring size.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        atom_constraints_asdict(py, &self.0)
    }
}

impl AtomConstraintsAst {
    /// The wrapped AST constraints — read access for atom construction.
    pub(crate) fn inner(&self) -> &AstAtomConstraintsAst {
        &self.0
    }

    /// Mutable access to the wrapped AST constraints — for the value-backed proxy.
    pub(crate) fn inner_mut(&mut self) -> &mut AstAtomConstraintsAst {
        &mut self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge). Test-only —
    /// in-crate construction wraps `AtomConstraintsAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(constraints: AstAtomConstraintsAst) -> Self {
        AtomConstraintsAst(constraints)
    }
}

/// Build the per-constraint iterator handle from a borrowed container.
fn atom_constraints_iter(
    py: Python<'_>,
    constraints: &AstAtomConstraintsAst,
) -> PyResult<AtomConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| into_py_variant(py, AtomConstraintAst::from_ast(py, constraint)?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AtomConstraintIter {
        entries: entries.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// inner-value mirrors. Ring memberships key by scope: `ring_count` for the
/// all-rings scope, `ring_size_count_<n>` for a specific ring size.
pub(crate) fn atom_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &AstAtomConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            AstAtomConstraintAst::Valence(v) => {
                dict.set_item("valence", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::DonatedPairs(v) => {
                dict.set_item("donated_pairs", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::AcceptedPairs(v) => {
                dict.set_item("accepted_pairs", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::AromaticValence(c) => {
                dict.set_item("aromatic_valence", AromaticValenceAst::from_ast(py, c)?)?
            }
            AstAtomConstraintAst::MulticenterValence(c) => {
                dict.set_item("multicenter_valence", MulticenterValenceAst::from_ast(py, c)?)?
            }
            AstAtomConstraintAst::TetrahedralStereo(c) => {
                dict.set_item("tetrahedral_stereo", TetrahedralStereoAst::from_ast(py, c)?)?
            }
            AstAtomConstraintAst::Degree(v) => dict.set_item("degree", ValueAst::from_ast(py, v)?)?,
            AstAtomConstraintAst::TotalDegree(v) => {
                dict.set_item("total_degree", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::TotalValence(v) => {
                dict.set_item("total_valence", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::RingDegree(v) => {
                dict.set_item("ring_degree", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::RingValence(v) => {
                dict.set_item("ring_valence", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::TotalHydrogens(v) => {
                dict.set_item("total_hydrogens", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::RingMembership(m) => {
                let key = match m.scope {
                    AstRingScope::All => "ring_count".to_string(),
                    AstRingScope::Size(size) => format!("ring_size_count_{size}"),
                };
                dict.set_item(key, ValueAst::from_ast(py, &m.count)?)?
            }
        }
    }
    Ok(dict)
}

/// What an `AtomConstraintsView` writes through to: an atom within a molecule
/// (by index) or a standalone `AtomAst`.
pub(crate) enum ConstraintsBacking {
    Molecule { owner: Py<MoleculeAst>, id: AstAtomId },
    Atom(Py<AtomAst>),
}

/// A live handle onto one atom's constraints, backed by either a molecule-atom or
/// a standalone `AtomAst`. Reads borrow the atom's constraints and read only the
/// item they need (no whole-container clone); mutators write through to the atom in
/// place, without a clone-and-writeback.
#[pyclass]
pub struct AtomConstraintsView {
    pub(crate) backing: ConstraintsBacking,
}

impl AtomConstraintsView {
    /// Borrow the backing atom's constraints and read one item through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstAtomConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            ConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .atoms()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("atom id out of range"))?;
                f(&view.ast.constraints)
            }
            ConstraintsBacking::Atom(atom) => {
                let atom = atom.bind(py).borrow();
                f(&atom.inner().constraints)
            }
        }
    }

    /// Mutate the backing atom's constraints in place through `f`.
    fn with_mut<R>(&self, py: Python<'_>, f: impl FnOnce(&mut AstAtomConstraintsAst) -> R) -> R {
        match &self.backing {
            ConstraintsBacking::Molecule { owner, id } => {
                f(&mut owner.borrow_mut(py).inner_mut().atom_mut(*id).ast.constraints)
            }
            ConstraintsBacking::Atom(atom) => f(&mut atom.borrow_mut(py).inner_mut().constraints),
        }
    }

    /// Set one constraint on the backing atom in place (last-wins per key).
    fn set_ast(&self, py: Python<'_>, constraint: AstAtomConstraintAst) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing atom in place, returning the removed entry.
    fn remove_ast(
        &self,
        py: Python<'_>,
        key: AstAtomConstraintKey,
    ) -> Option<AstAtomConstraintAst> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl AtomConstraintsView {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("AtomConstraintsView({count} entries)"))
    }

    /// Insert `c` on the atom in place, replacing any existing entry of the same
    /// key (last-wins).
    fn set(&self, py: Python<'_>, c: Py<AtomConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key from the atom in place, returning it if
    /// present.
    fn remove(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<Option<AtomConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_ast(py))
            .map(|c| AtomConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> PyResult<()> {
        if self
            .remove_ast(py, key.bind(py).borrow().to_ast(py))
            .is_some()
        {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    /// Overlay `other` onto the atom's constraints in place — another container or an
    /// iterable of `AtomConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&self, py: Python<'_>, other: ConstraintsUpdate) -> PyResult<()> {
        self.with_mut(py, |cs| other.apply(py, cs))
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    fn __iter__(&self, py: Python<'_>) -> PyResult<AtomConstraintIter> {
        self.read(py, |cs| atom_constraints_iter(py, cs))
    }

    /// The constraint with the given key, or `None`.
    fn get(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<Option<AtomConstraintAst>> {
        let key = key.bind(py).borrow().to_ast(py);
        self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| AtomConstraintAst::from_ast(py, constraint))
                .transpose()
        })
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> PyResult<AtomConstraintAst> {
        let ast_key = key.bind(py).borrow().to_ast(py);
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| AtomConstraintAst::from_ast(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_ast(py);
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// The valence value, or `None`.
    #[getter]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.valence().map(|v| ValueAst::from_ast(py, v)).transpose()
        })
    }

    #[setter]
    fn set_valence(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::valence(value.to_ast(py)));
    }

    /// The donated-pairs value, or `None`.
    #[getter]
    fn donated_pairs(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.donated_pairs()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_donated_pairs(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::donated_pairs(value.to_ast(py)));
    }

    /// The accepted-pairs value, or `None`.
    #[getter]
    fn accepted_pairs(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.accepted_pairs()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_accepted_pairs(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::accepted_pairs(value.to_ast(py)));
    }

    /// The aromatic-valence state, or `None`.
    #[getter]
    fn aromatic_valence(&self, py: Python<'_>) -> PyResult<Option<AromaticValenceAst>> {
        self.read(py, |cs| {
            cs.aromatic_valence()
                .map(|c| AromaticValenceAst::from_ast(py, c))
                .transpose()
        })
    }

    #[setter]
    fn set_aromatic_valence(&self, py: Python<'_>, value: AromaticValenceArg) -> PyResult<()> {
        self.set_ast(py, AstAtomConstraintAst::aromatic_valence(value.to_ast(py)?));
        Ok(())
    }

    /// The multicenter-valence state, or `None`.
    #[getter]
    fn multicenter_valence(&self, py: Python<'_>) -> PyResult<Option<MulticenterValenceAst>> {
        self.read(py, |cs| {
            cs.multicenter_valence()
                .map(|c| MulticenterValenceAst::from_ast(py, c))
                .transpose()
        })
    }

    #[setter]
    fn set_multicenter_valence(
        &self,
        py: Python<'_>,
        value: MulticenterValenceArg,
    ) -> PyResult<()> {
        self.set_ast(py, AstAtomConstraintAst::multicenter_valence(value.to_ast(py)?));
        Ok(())
    }

    /// The tetrahedral-stereo state, or `None`.
    #[getter]
    fn tetrahedral_stereo(&self, py: Python<'_>) -> PyResult<Option<TetrahedralStereoAst>> {
        self.read(py, |cs| {
            cs.tetrahedral_stereo()
                .map(|c| TetrahedralStereoAst::from_ast(py, c))
                .transpose()
        })
    }

    #[setter]
    fn set_tetrahedral_stereo(&self, py: Python<'_>, value: TetrahedralStereoArg) -> PyResult<()> {
        self.set_ast(py, AstAtomConstraintAst::tetrahedral_stereo(value.to_ast(py)?));
        Ok(())
    }

    /// The degree value, or `None`.
    #[getter]
    fn degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.degree().map(|v| ValueAst::from_ast(py, v)).transpose()
        })
    }

    #[setter]
    fn set_degree(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::degree(value.to_ast(py)));
    }

    /// The total-degree value, or `None`.
    #[getter]
    fn total_degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.total_degree()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_total_degree(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::total_degree(value.to_ast(py)));
    }

    /// The total-valence value, or `None`.
    #[getter]
    fn total_valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.total_valence()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_total_valence(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::total_valence(value.to_ast(py)));
    }

    /// The ring-degree value, or `None`.
    #[getter]
    fn ring_degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_degree()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_ring_degree(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::ring_degree(value.to_ast(py)));
    }

    /// The ring-valence value, or `None`.
    #[getter]
    fn ring_valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_valence()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_ring_valence(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::ring_valence(value.to_ast(py)));
    }

    /// The total-hydrogens value, or `None`.
    #[getter]
    fn total_hydrogens(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.total_hydrogens()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_total_hydrogens(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::total_hydrogens(value.to_ast(py)));
    }

    /// The all-rings membership count, or `None`.
    #[getter]
    fn ring_count(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_count()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_ring_count(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(
            py,
            AstAtomConstraintAst::ring_membership(AstRingScope::All, value.to_ast(py)),
        );
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    fn ring_size_count(&self, py: Python<'_>) -> RingSizeCounts {
        let backing = match &self.backing {
            ConstraintsBacking::Molecule { owner, id } => RingSizeBacking::Molecule {
                owner: owner.clone_ref(py),
                id: *id,
            },
            ConstraintsBacking::Atom(atom) => RingSizeBacking::Atom(atom.clone_ref(py)),
        };
        RingSizeCounts { backing }
    }

    /// The present constraints as a dict keyed by snake_case name.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| atom_constraints_asdict(py, cs))
    }
}

/// What a `RingSizeCounts` proxy reads/writes through to: an atom within a molecule,
/// a standalone `AtomAst`, or a standalone `AtomConstraintsAst` value.
pub(crate) enum RingSizeBacking {
    Molecule { owner: Py<MoleculeAst>, id: AstAtomId },
    Atom(Py<AtomAst>),
    Value(Py<AtomConstraintsAst>),
}

/// A subscriptable proxy over the sized-ring membership counts of an atom, keyed by
/// ring size: `proxy[size]` reads, `proxy[size] = count` sets, `del proxy[size]`
/// removes. Backs onto whichever container produced it (dual-backing, like
/// `AtomConstraintsView`).
#[pyclass]
pub struct RingSizeCounts {
    backing: RingSizeBacking,
}

impl RingSizeCounts {
    /// Borrow the backing constraints and read through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstAtomConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            RingSizeBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .atoms()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("atom id out of range"))?;
                f(&view.ast.constraints)
            }
            RingSizeBacking::Atom(atom) => f(&atom.bind(py).borrow().inner().constraints),
            RingSizeBacking::Value(value) => f(value.bind(py).borrow().inner()),
        }
    }

    /// Mutate the backing constraints in place through `f`.
    fn write(&self, py: Python<'_>, f: impl FnOnce(&mut AstAtomConstraintsAst)) {
        match &self.backing {
            RingSizeBacking::Molecule { owner, id } => {
                f(&mut owner.borrow_mut(py).inner_mut().atom_mut(*id).ast.constraints)
            }
            RingSizeBacking::Atom(atom) => f(&mut atom.borrow_mut(py).inner_mut().constraints),
            RingSizeBacking::Value(value) => f(value.borrow_mut(py).inner_mut()),
        }
    }
}

#[pymethods]
impl RingSizeCounts {
    /// The membership count for rings of `size`, or `None`.
    fn __getitem__(&self, py: Python<'_>, size: u8) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_size_count(size)
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    /// The number of distinct ring sizes with a membership constraint.
    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(ring_sizes(cs).count()))
    }

    fn __contains__(&self, py: Python<'_>, size: u8) -> PyResult<bool> {
        self.read(py, |cs| Ok(cs.ring_size_count(size).is_some()))
    }

    /// Iterate the present ring sizes (as ints).
    fn __iter__(&self, py: Python<'_>) -> PyResult<RingSizeIter> {
        let sizes = self.read(py, |cs| Ok(ring_sizes(cs).collect::<Vec<u8>>()))?;
        Ok(RingSizeIter {
            sizes: sizes.into_iter(),
        })
    }

    /// Set the membership count for rings of `size` in place.
    fn __setitem__(&self, py: Python<'_>, size: u8, count: ValueArg) {
        let constraint =
            AstAtomConstraintAst::ring_membership(AstRingScope::Size(size), count.to_ast(py));
        self.write(py, |cs| cs.set(constraint));
    }

    /// Remove the sized-ring membership for `size` in place.
    fn __delitem__(&self, py: Python<'_>, size: u8) {
        self.write(py, |cs| {
            cs.remove(AstAtomConstraintKey::RingMembership(AstRingScope::Size(size)));
        });
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |cs| {
            let mut parts = Vec::new();
            for entry in cs.iter() {
                if let AstAtomConstraintAst::RingMembership(m) = entry {
                    if let AstRingScope::Size(size) = m.scope {
                        let count = into_py_variant(py, ValueAst::from_ast(py, &m.count)?)?;
                        parts.push(format!(
                            "{size}: {}",
                            count.bind(py).as_any().repr()?.extract::<String>()?
                        ));
                    }
                }
            }
            Ok(format!("RingSizeCounts({{{}}})", parts.join(", ")))
        })
    }
}

/// The ring sizes with a membership constraint, in kind-sorted order.
fn ring_sizes(constraints: &AstAtomConstraintsAst) -> impl Iterator<Item = u8> + '_ {
    constraints.iter().filter_map(|entry| match entry {
        AstAtomConstraintAst::RingMembership(m) => match m.scope {
            AstRingScope::Size(size) => Some(size),
            AstRingScope::All => None,
        },
        _ => None,
    })
}

#[pyclass]
struct RingSizeIter {
    sizes: IntoIter<u8>,
}

#[pymethods]
impl RingSizeIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<u8> {
        self.sizes.next()
    }
}

#[pyclass]
struct AtomConstraintIter {
    entries: IntoIter<Py<AtomConstraintAst>>,
}

#[pymethods]
impl AtomConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<AtomConstraintAst>> {
        self.entries.next()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AtomAst as AstAtomAst, MoleculeAst as AstMoleculeAst, StereoCosetAst as AstStereoCosetAst,
        TetrahedralStereoAst as AstTetrahedralStereoAst, ValueAst as AstValueAst,
    };
    use umol_chem::element::Element as ChemElement;

    use super::*;

    #[rstest]
    #[case(AstAromaticValenceAst::Undetermined)]
    #[case(AstAromaticValenceAst::NotAromatic)]
    #[case(AstAromaticValenceAst::aromatic(1))]
    fn test_aromatic_valence_ast_roundtrip(#[case] ast: AstAromaticValenceAst) {
        Python::attach(|py| {
            assert_eq!(
                AromaticValenceAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(AstMulticenterValenceAst::Undetermined)]
    #[case(AstMulticenterValenceAst::NotMulticenter)]
    #[case(AstMulticenterValenceAst::multicenter(2))]
    fn test_multicenter_valence_ast_roundtrip(#[case] ast: AstMulticenterValenceAst) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterValenceAst::from_ast(py, &ast)
                    .unwrap()
                    .to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(AstRingScope::All)]
    #[case(AstRingScope::Size(6))]
    fn test_ring_scope_roundtrip(#[case] ast: AstRingScope) {
        assert_eq!(RingScope::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstRingMembershipAst::new(AstRingScope::All, 2))]
    #[case(AstRingMembershipAst::new(AstRingScope::Size(6), 1))]
    fn test_ring_membership_ast_roundtrip(#[case] ast: AstRingMembershipAst) {
        Python::attach(|py| {
            assert_eq!(
                RingMembershipAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(AstAtomConstraintAst::valence(4))]
    #[case(AstAtomConstraintAst::aromatic_valence(AstAromaticValenceAst::aromatic(1)))]
    #[case(AstAtomConstraintAst::ring_membership(AstRingScope::All, 2))]
    #[case(AstAtomConstraintAst::tetrahedral_stereo(AstTetrahedralStereoAst::not_stereo()))]
    fn test_atom_constraint_roundtrip(#[case] ast: AstAtomConstraintAst) {
        Python::attach(|py| {
            assert_eq!(
                AtomConstraintAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_len_contains() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::degree(3)).unwrap(),
            )
            .unwrap();
            let constraints = AtomConstraintsAst::new(py, vec![valence, degree]);
            assert_eq!(constraints.__len__(), 2);
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, AtomConstraintKey::Valence()).unwrap()
            ));
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, AtomConstraintKey::Degree()).unwrap()
            ));
            assert!(!constraints.__contains__(
                py,
                into_py_variant(py, AtomConstraintKey::TotalHydrogens()).unwrap()
            ));
        });
    }

    #[rstest]
    fn test_atom_constraints_valence() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::degree(3)).unwrap(),
            )
            .unwrap();
            let constraints = AtomConstraintsAst::new(py, vec![valence, degree]);
            assert_eq!(
                constraints.valence(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(4)
            );
            assert_eq!(
                constraints.degree(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(3)
            );
            assert!(constraints.total_valence(py).unwrap().is_none());
            assert!(constraints.aromatic_valence(py).unwrap().is_none());
        });
    }

    #[rstest]
    fn test_atom_constraints_ring_size_count() {
        Python::attach(|py| {
            let membership = into_py_variant(
                py,
                AtomConstraintAst::from_ast(
                    py,
                    &AstAtomConstraintAst::ring_membership(AstRingScope::Size(6), 1),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![membership])).unwrap();
            let proxy = AtomConstraintsAst::ring_size_count(constraints.clone_ref(py));
            assert_eq!(
                proxy.__getitem__(py, 6).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(1)
            );
            assert!(proxy.__getitem__(py, 5).unwrap().is_none());
            assert!(constraints
                .bind(py)
                .borrow()
                .ring_count(py)
                .unwrap()
                .is_none());
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            constraints.set(py, valence);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.valence(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(4)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_remove() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            let mut constraints = AtomConstraintsAst::new(py, vec![valence]);
            let removed = constraints
                .remove(py, into_py_variant(py, AtomConstraintKey::Valence()).unwrap())
                .unwrap();
            match removed {
                Some(AtomConstraintAst::Valence(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(4))
                }
                _ => panic!("expected removed Valence(Lit(4))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_update() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            let mut other = AstAtomConstraintsAst::new();
            other.set(AstAtomConstraintAst::valence(4));
            other.set(AstAtomConstraintAst::degree(3));
            constraints
                .update(
                    py,
                    ConstraintsUpdate::Container(
                        Py::new(py, AtomConstraintsAst::from_inner(other)).unwrap(),
                    ),
                )
                .unwrap();
            assert_eq!(constraints.__len__(), 2);
            assert_eq!(
                constraints.valence(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(4)
            );
            assert_eq!(
                constraints.degree(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(3)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_view_set() {
        Python::attach(|py| {
            let owner = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_atoms_and_bonds(
                    vec![AstAtomAst::from_element(ChemElement::C)],
                    vec![],
                )),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: ConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAtomId(0),
                },
            };
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            view.set(py, valence);
            let fresh = AtomConstraintsView {
                backing: ConstraintsBacking::Molecule {
                    owner,
                    id: AstAtomId(0),
                },
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            match fresh
                .get(py, into_py_variant(py, AtomConstraintKey::Valence()).unwrap())
                .unwrap()
            {
                Some(AtomConstraintAst::Valence(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(4))
                }
                _ => panic!("expected Valence(Lit(4))"),
            }
        });
    }

    #[rstest]
    fn test_atom_constraints_view_remove() {
        Python::attach(|py| {
            let atom = AstAtomAst::from_element(ChemElement::C)
                .with_constraint(AstAtomConstraintAst::valence(4));
            let owner = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_atoms_and_bonds(vec![atom], vec![])),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: ConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAtomId(0),
                },
            };
            let removed = view
                .remove(py, into_py_variant(py, AtomConstraintKey::Valence()).unwrap())
                .unwrap();
            match removed {
                Some(AtomConstraintAst::Valence(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(4))
                }
                _ => panic!("expected removed Valence(Lit(4))"),
            }
            let fresh = AtomConstraintsView {
                backing: ConstraintsBacking::Molecule {
                    owner,
                    id: AstAtomId(0),
                },
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_atom_constraints_view_update() {
        Python::attach(|py| {
            let owner = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_atoms_and_bonds(
                    vec![AstAtomAst::from_element(ChemElement::C)],
                    vec![],
                )),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: ConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAtomId(0),
                },
            };
            let mut other = AstAtomConstraintsAst::new();
            other.set(AstAtomConstraintAst::valence(4));
            other.set(AstAtomConstraintAst::degree(3));
            view
                .update(
                    py,
                    ConstraintsUpdate::Container(
                        Py::new(py, AtomConstraintsAst::from_inner(other)).unwrap(),
                    ),
                )
                .unwrap();
            let fresh = AtomConstraintsView {
                backing: ConstraintsBacking::Molecule {
                    owner,
                    id: AstAtomId(0),
                },
            };
            assert_eq!(fresh.__len__(py).unwrap(), 2);
        });
    }

    #[rstest]
    fn test_atom_constraints_view_set_atom_backed() {
        Python::attach(|py| {
            let atom =
                Py::new(py, AtomAst::from_inner(AstAtomAst::from_element(ChemElement::C))).unwrap();
            let view = AtomConstraintsView {
                backing: ConstraintsBacking::Atom(atom.clone_ref(py)),
            };
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            view.set(py, valence);
            // a fresh view proves the write hit the standalone atom, not a copy
            let fresh = AtomConstraintsView {
                backing: ConstraintsBacking::Atom(atom),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            match fresh
                .get(py, into_py_variant(py, AtomConstraintKey::Valence()).unwrap())
                .unwrap()
            {
                Some(AtomConstraintAst::Valence(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(4))
                }
                _ => panic!("expected Valence(Lit(4))"),
            }
        });
    }

    #[rstest]
    fn test_atom_constraints_view_remove_atom_backed() {
        Python::attach(|py| {
            let atom = Py::new(
                py,
                AtomAst::from_inner(
                    AstAtomAst::from_element(ChemElement::C)
                        .with_constraint(AstAtomConstraintAst::valence(4)),
                ),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: ConstraintsBacking::Atom(atom.clone_ref(py)),
            };
            let removed = view
                .remove(py, into_py_variant(py, AtomConstraintKey::Valence()).unwrap())
                .unwrap();
            match removed {
                Some(AtomConstraintAst::Valence(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(4))
                }
                _ => panic!("expected removed Valence(Lit(4))"),
            }
            let fresh = AtomConstraintsView {
                backing: ConstraintsBacking::Atom(atom),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_atom_constraints_view_update_atom_backed() {
        Python::attach(|py| {
            let atom =
                Py::new(py, AtomAst::from_inner(AstAtomAst::from_element(ChemElement::C))).unwrap();
            let view = AtomConstraintsView {
                backing: ConstraintsBacking::Atom(atom.clone_ref(py)),
            };
            let mut other = AstAtomConstraintsAst::new();
            other.set(AstAtomConstraintAst::valence(4));
            other.set(AstAtomConstraintAst::degree(3));
            view
                .update(
                    py,
                    ConstraintsUpdate::Container(
                        Py::new(py, AtomConstraintsAst::from_inner(other)).unwrap(),
                    ),
                )
                .unwrap();
            let fresh = AtomConstraintsView {
                backing: ConstraintsBacking::Atom(atom),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 2);
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_valence() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints.set_valence(py, ValueArg::Lit(4));
            assert_eq!(
                constraints.valence(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(4)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_ring_count() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints.set_ring_count(py, ValueArg::Lit(2));
            assert_eq!(
                constraints.ring_count(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(2)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_aromatic_valence() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints
                .set_aromatic_valence(py, AromaticValenceArg::Value(ValueArg::Lit(1)))
                .unwrap();
            match constraints.aromatic_valence(py).unwrap().unwrap() {
                AromaticValenceAst::Aromatic(v) => assert_eq!(v.to_ast(py), AstValueAst::Lit(1)),
                _ => panic!("expected Aromatic"),
            }
            constraints
                .set_aromatic_valence(py, AromaticValenceArg::Flag(false))
                .unwrap();
            match constraints.aromatic_valence(py).unwrap().unwrap() {
                AromaticValenceAst::NotAromatic() => {}
                _ => panic!("expected NotAromatic"),
            }
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_aromatic_valence_error() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            assert!(constraints
                .set_aromatic_valence(py, AromaticValenceArg::Flag(true))
                .is_err());
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_tetrahedral_stereo() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints
                .set_tetrahedral_stereo(py, TetrahedralStereoArg::Config(TetrahedralStereo::Cw))
                .unwrap();
            match constraints.tetrahedral_stereo(py).unwrap().unwrap() {
                TetrahedralStereoAst::Stereo(coset) => {
                    assert_eq!(coset.bind(py).borrow().to_ast(py), AstStereoCosetAst::Lit(1))
                }
                _ => panic!("expected Stereo"),
            }
        });
    }

    #[rstest]
    fn test_atom_constraints_view_set_aromatic_valence() {
        Python::attach(|py| {
            let owner = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_atoms_and_bonds(
                    vec![AstAtomAst::from_element(ChemElement::C)],
                    vec![],
                )),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: ConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAtomId(0),
                },
            };
            view.set_aromatic_valence(py, AromaticValenceArg::Value(ValueArg::Lit(1)))
                .unwrap();
            let fresh = AtomConstraintsView {
                backing: ConstraintsBacking::Molecule {
                    owner,
                    id: AstAtomId(0),
                },
            };
            match fresh.aromatic_valence(py).unwrap().unwrap() {
                AromaticValenceAst::Aromatic(v) => assert_eq!(v.to_ast(py), AstValueAst::Lit(1)),
                _ => panic!("expected Aromatic"),
            }
        });
    }

    #[rstest]
    fn test_ring_size_counts_value_backed() {
        Python::attach(|py| {
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![])).unwrap();
            let proxy = AtomConstraintsAst::ring_size_count(constraints.clone_ref(py));
            proxy.__setitem__(py, 6, ValueArg::Lit(3));
            assert_eq!(
                proxy.__getitem__(py, 6).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(3)
            );
            proxy.__delitem__(py, 6);
            assert!(proxy.__getitem__(py, 6).unwrap().is_none());
        });
    }

    #[rstest]
    fn test_ring_size_counts_molecule_backed() {
        Python::attach(|py| {
            let owner = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_atoms_and_bonds(
                    vec![AstAtomAst::from_element(ChemElement::C)],
                    vec![],
                )),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: ConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAtomId(0),
                },
            };
            view.ring_size_count(py).__setitem__(py, 5, ValueArg::Lit(1));
            let fresh = AtomConstraintsView {
                backing: ConstraintsBacking::Molecule {
                    owner,
                    id: AstAtomId(0),
                },
            };
            assert_eq!(
                fresh
                    .ring_size_count(py)
                    .__getitem__(py, 5)
                    .unwrap()
                    .unwrap()
                    .to_ast(py),
                AstValueAst::Lit(1)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_update_entries() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::degree(3)).unwrap(),
            )
            .unwrap();
            constraints
                .update(py, ConstraintsUpdate::Entries(vec![valence, degree]))
                .unwrap();
            assert_eq!(constraints.__len__(), 2);
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = AtomConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, AtomConstraintKey::Valence()).unwrap();
            assert!(constraints.__getitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, AtomConstraintKey::Valence()).unwrap();
            assert!(constraints.__delitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_ring_size_counts_len_iter_contains() {
        Python::attach(|py| {
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![])).unwrap();
            let proxy = AtomConstraintsAst::ring_size_count(constraints.clone_ref(py));
            proxy.__setitem__(py, 6, ValueArg::Lit(3));
            proxy.__setitem__(py, 5, ValueArg::Lit(1));
            assert_eq!(proxy.__len__(py).unwrap(), 2);
            assert!(proxy.__contains__(py, 6).unwrap());
            assert!(!proxy.__contains__(py, 4).unwrap());
            let mut iter = proxy.__iter__(py).unwrap();
            let mut sizes = Vec::new();
            while let Some(size) = iter.__next__() {
                sizes.push(size);
            }
            sizes.sort_unstable();
            assert_eq!(sizes, vec![5, 6]);
        });
    }
}
