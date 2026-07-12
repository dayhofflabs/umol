//! Dative bond value type and dative-bond-constraint surface mirroring
//! `umol_ast::ast`: `DativeBondAst`, the `DativeBondConstraintAst` enum, the
//! `DativeBondConstraintsAst` container, and the `DativeBondConstraintsView` live
//! handle. A dative bond carries only an order and bond-scope constraints; the
//! acceptor and donor atoms are the participants of the owning molecule's dative
//! relation, so they are topology (the view half) rather than part of the value.

use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use umol_ast::ast::{
    DativeBondAst as AstDativeBondAst, DativeBondConstraintAst as AstDativeBondConstraintAst,
    DativeBondConstraintKey as AstDativeBondConstraintKey,
    DativeBondConstraintsAst as AstDativeBondConstraintsAst, RingScope as AstRingScope,
};

use crate::boolean::{BooleanArg, BooleanAst};
use crate::constraint::{RingMembershipAst, RingScope};
use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::error::parse_error;
use crate::value::{ValueArg, ValueAst};

/// A dative bond: order and bond-scope constraints.
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct DativeBondAst(AstDativeBondAst);

#[pymethods]
impl DativeBondAst {
    /// Construct from an order — an `int` or a `ValueAst` expression — optionally
    /// setting constraints.
    #[new]
    #[pyo3(signature = (order, *, constraints=None))]
    fn new(
        py: Python<'_>,
        order: ValueArg,
        constraints: Option<Py<DativeBondConstraintsAst>>,
    ) -> Self {
        let mut bond = AstDativeBondAst::new(order.to_ast(py));
        if let Some(constraints) = constraints {
            bond.constraints = constraints.bind(py).borrow().inner().clone();
        }
        DativeBondAst(bond)
    }

    /// Parse a dative-bond-DSL string (e.g. `"1#R(6)"`) into a `DativeBondAst`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        AstDativeBondAst::from_str(s).map(Self).map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("DativeBondAst.parse('{}')", self.0)
    }

    #[getter]
    fn order(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.order)
    }

    #[setter]
    fn set_order(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.order = value.to_ast(py);
    }

    /// The dative bond's constraints as a live handle onto this bond: reads borrow
    /// the current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(slf: Py<Self>) -> DativeBondConstraintsView {
        DativeBondConstraintsView {
            backing: DativeBondConstraintsBacking::DativeBond(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or
    /// a live view.
    #[setter]
    fn set_constraints(&mut self, py: Python<'_>, value: DativeBondConstraintsArg) -> PyResult<()> {
        self.0.constraints = value.to_ast(py)?;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are the field mirrors.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("order", self.order(py)?)?;
        dict.set_item(
            "constraints",
            dative_bond_constraints_asdict(py, &self.0.constraints)?,
        )?;
        Ok(dict)
    }
}

impl DativeBondAst {
    /// The wrapped AST bond — read access for the bond-backed constraints view.
    pub(crate) fn inner(&self) -> &AstDativeBondAst {
        &self.0
    }

    /// Mutable access to the wrapped AST bond — write access for the bond-backed
    /// constraints view.
    pub(crate) fn inner_mut(&mut self) -> &mut AstDativeBondAst {
        &mut self.0
    }

    /// Wrap an AST bond (the hold-the-value `from_inner` bridge, paired with
    /// `inner`). Test-only — in-crate construction wraps `DativeBondAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(bond: AstDativeBondAst) -> Self {
        DativeBondAst(bond)
    }
}

/// The key (identity) of a dative bond constraint, for keyed lookup. The
/// ring-membership key carries its ring scope; all other keys are the bare
/// discriminant.
#[pyclass]
pub enum DativeBondConstraintKey {
    Aromatic(),
    RingMembership(Py<RingScope>),
}

#[pymethods]
impl DativeBondConstraintKey {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            DativeBondConstraintKey::Aromatic() => ("Aromatic", 0),
            DativeBondConstraintKey::RingMembership(_) => ("RingMembership", 1),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "DativeBondConstraintKey",
            variant,
            arity,
        )
    }
}

impl DativeBondConstraintKey {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstDativeBondConstraintKey) -> PyResult<Self> {
        Ok(match ast {
            AstDativeBondConstraintKey::Aromatic => Self::Aromatic(),
            AstDativeBondConstraintKey::RingMembership(scope) => {
                Self::RingMembership(into_py_variant(py, RingScope::from_ast(scope))?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstDativeBondConstraintKey {
        match self {
            Self::Aromatic() => AstDativeBondConstraintKey::Aromatic,
            Self::RingMembership(scope) => {
                AstDativeBondConstraintKey::RingMembership(scope.bind(py).borrow().to_ast())
            }
        }
    }
}

/// A dative-bond-scope constraint: the aromatic flag or a ring membership of a
/// single dative bond.
#[pyclass]
pub enum DativeBondConstraintAst {
    Aromatic(Py<BooleanAst>),
    RingMembership(Py<RingMembershipAst>),
}

#[pymethods]
impl DativeBondConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    fn key(&self, py: Python<'_>) -> PyResult<DativeBondConstraintKey> {
        DativeBondConstraintKey::from_ast(py, &self.to_ast(py).key())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            DativeBondConstraintAst::Aromatic(_) => "Aromatic",
            DativeBondConstraintAst::RingMembership(_) => "RingMembership",
        };
        variant_repr(slf.bind(py).as_any(), "DativeBondConstraintAst", variant, 1)
    }
}

impl DativeBondConstraintAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstDativeBondConstraintAst) -> PyResult<Self> {
        Ok(match ast {
            AstDativeBondConstraintAst::Aromatic(b) => {
                Self::Aromatic(into_py_variant(py, BooleanAst::from_ast(b))?)
            }
            AstDativeBondConstraintAst::RingMembership(m) => {
                Self::RingMembership(into_py_variant(py, RingMembershipAst::from_ast(py, m)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstDativeBondConstraintAst {
        match self {
            Self::Aromatic(b) => AstDativeBondConstraintAst::Aromatic(b.bind(py).borrow().to_ast()),
            Self::RingMembership(m) => {
                AstDativeBondConstraintAst::RingMembership(m.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container (value or live view) or
/// an iterable of `DativeBondConstraintAst` (each `set`, last-wins).
#[derive(FromPyObject)]
enum DativeBondConstraintsUpdate {
    Container(Py<DativeBondConstraintsAst>),
    View(Py<DativeBondConstraintsView>),
    Entries(Vec<Py<DativeBondConstraintAst>>),
}

impl DativeBondConstraintsUpdate {
    /// Overlay this update onto `target` in place.
    fn apply(&self, py: Python<'_>, target: &mut AstDativeBondConstraintsAst) -> PyResult<()> {
        match self {
            DativeBondConstraintsUpdate::Container(c) => target.update(c.bind(py).borrow().inner()),
            DativeBondConstraintsUpdate::View(v) => {
                let snapshot = v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?;
                target.update(&snapshot);
            }
            DativeBondConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry.bind(py).borrow().to_ast(py));
                }
            }
        }
        Ok(())
    }
}

/// A whole-container argument that snapshots either a value container or a live
/// view — for the dative bond `constraints` setter, which accepts either.
#[derive(FromPyObject)]
enum DativeBondConstraintsArg {
    Container(Py<DativeBondConstraintsAst>),
    View(Py<DativeBondConstraintsView>),
}

impl DativeBondConstraintsArg {
    fn to_ast(&self, py: Python<'_>) -> PyResult<AstDativeBondConstraintsAst> {
        match self {
            DativeBondConstraintsArg::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
            DativeBondConstraintsArg::View(v) => v.bind(py).borrow().read(py, |cs| Ok(cs.clone())),
        }
    }
}

/// The dative-bond-scope constraints on a dative bond, in kind-sorted order.
/// Mutable, hence value-equal but unhashable (matching `DativeBondAst`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct DativeBondConstraintsAst(AstDativeBondConstraintsAst);

#[pymethods]
impl DativeBondConstraintsAst {
    /// Build from a sequence of constraints (kind-sorted; a unique kind replaces
    /// an earlier one, ring memberships accumulate per scope).
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<DativeBondConstraintAst>>) -> Self {
        let mut constraints = AstDativeBondConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_ast(py)),
        );
        DativeBondConstraintsAst(constraints)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let mirror = into_py_variant(py, DativeBondConstraintAst::from_ast(py, entry)?)?;
            parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("DativeBondConstraintsAst([{}])", parts.join(", ")))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    fn set(&mut self, py: Python<'_>, c: Py<DativeBondConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<Option<DativeBondConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_ast(py))
            .map(|c| DativeBondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container or an iterable of
    /// `DativeBondConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&mut self, py: Python<'_>, other: DativeBondConstraintsUpdate) -> PyResult<()> {
        other.apply(py, &mut self.0)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<DativeBondConstraintKeyIter> {
        dative_bond_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<DativeBondConstraintKeyIter> {
        dative_bond_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<DativeBondConstraintIter> {
        dative_bond_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<DativeBondConstraintItemsIter> {
        dative_bond_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_ast(py)) {
            Some(constraint) => Ok(into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(py, constraint)?,
            )?
            .into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<DativeBondConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_ast(py)) {
            Some(constraint) => DativeBondConstraintAst::from_ast(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&mut self, py: Python<'_>, key: Py<DativeBondConstraintKey>) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_ast(py)).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<DativeBondConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast(py))
    }

    /// The aromatic value; `Undetermined` when no `Aromatic` constraint is present
    /// (mirroring the non-optional Rust accessor).
    #[getter]
    fn aromatic(&self) -> BooleanAst {
        BooleanAst::from_ast(&self.0.aromatic())
    }

    #[setter]
    fn set_aromatic(&mut self, py: Python<'_>, value: BooleanArg) {
        self.0
            .set(AstDativeBondConstraintAst::aromatic(value.to_ast(py)));
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
        self.0.set(AstDativeBondConstraintAst::ring_membership(
            AstRingScope::All,
            value.to_ast(py),
        ));
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    fn ring_size_count(slf: Py<Self>) -> DativeBondRingSizeCounts {
        DativeBondRingSizeCounts {
            backing: DativeBondRingSizeBacking::Value(slf),
        }
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors. Ring memberships key by scope: `ring_count` for the
    /// all-rings scope, `ring_size_count_<n>` for a specific ring size.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        dative_bond_constraints_asdict(py, &self.0)
    }
}

impl DativeBondConstraintsAst {
    /// The wrapped AST constraints — read access for dative bond construction.
    pub(crate) fn inner(&self) -> &AstDativeBondConstraintsAst {
        &self.0
    }

    /// Mutable access to the wrapped AST constraints — for the value-backed proxy.
    pub(crate) fn inner_mut(&mut self) -> &mut AstDativeBondConstraintsAst {
        &mut self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge). Test-only —
    /// in-crate construction wraps `DativeBondConstraintsAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(constraints: AstDativeBondConstraintsAst) -> Self {
        DativeBondConstraintsAst(constraints)
    }
}

/// Build the per-constraint iterator handle from a borrowed container.
fn dative_bond_constraints_iter(
    py: Python<'_>,
    constraints: &AstDativeBondConstraintsAst,
) -> PyResult<DativeBondConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| into_py_variant(py, DativeBondConstraintAst::from_ast(py, constraint)?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(DativeBondConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
fn dative_bond_constraint_keys(
    py: Python<'_>,
    constraints: &AstDativeBondConstraintsAst,
) -> PyResult<DativeBondConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(
                py,
                DativeBondConstraintKey::from_ast(py, &constraint.key())?,
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(DativeBondConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
fn dative_bond_constraint_items(
    py: Python<'_>,
    constraints: &AstDativeBondConstraintsAst,
) -> PyResult<DativeBondConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(
                    py,
                    DativeBondConstraintKey::from_ast(py, &constraint.key())?,
                )?,
                into_py_variant(py, DativeBondConstraintAst::from_ast(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(DativeBondConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// inner-value mirrors. Ring memberships key by scope: `ring_count` for the
/// all-rings scope, `ring_size_count_<n>` for a specific ring size.
fn dative_bond_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &AstDativeBondConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            AstDativeBondConstraintAst::Aromatic(b) => {
                dict.set_item("aromatic", BooleanAst::from_ast(b))?
            }
            AstDativeBondConstraintAst::RingMembership(m) => {
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

/// What a `DativeBondConstraintsView` writes through to. Only the standalone
/// `DativeBondAst` backing exists at this stage; the molecule-bond backing lands
/// with `DativeBondView`, which constructs it.
enum DativeBondConstraintsBacking {
    DativeBond(Py<DativeBondAst>),
}

/// A live handle onto one dative bond's constraints, backed by a standalone
/// `DativeBondAst`. Reads borrow the bond's constraints and read only the item
/// they need (no whole-container clone); mutators write through to the bond in
/// place, without a clone-and-writeback.
#[pyclass]
pub struct DativeBondConstraintsView {
    backing: DativeBondConstraintsBacking,
}

impl DativeBondConstraintsView {
    /// Borrow the backing bond's constraints and read one item through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstDativeBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            DativeBondConstraintsBacking::DativeBond(bond) => {
                let bond = bond.bind(py).borrow();
                f(&bond.inner().constraints)
            }
        }
    }

    /// Mutate the backing bond's constraints in place through `f`.
    fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut AstDativeBondConstraintsAst) -> R,
    ) -> R {
        match &self.backing {
            DativeBondConstraintsBacking::DativeBond(bond) => {
                f(&mut bond.borrow_mut(py).inner_mut().constraints)
            }
        }
    }

    /// Set one constraint on the backing bond in place (last-wins per key).
    fn set_ast(&self, py: Python<'_>, constraint: AstDativeBondConstraintAst) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing bond in place, returning the removed entry.
    fn remove_ast(
        &self,
        py: Python<'_>,
        key: AstDativeBondConstraintKey,
    ) -> Option<AstDativeBondConstraintAst> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl DativeBondConstraintsView {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("DativeBondConstraintsView({count} entries)"))
    }

    /// Insert `c` on the bond in place, replacing any existing entry of the same
    /// key (last-wins).
    fn set(&self, py: Python<'_>, c: Py<DativeBondConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key from the bond in place, returning it if
    /// present (dict `pop`).
    fn pop(
        &self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<Option<DativeBondConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_ast(py))
            .map(|c| DativeBondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&self, py: Python<'_>, key: Py<DativeBondConstraintKey>) -> PyResult<()> {
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

    /// Overlay `other` onto the bond's constraints in place — another container or an
    /// iterable of `DativeBondConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&self, py: Python<'_>, other: DativeBondConstraintsUpdate) -> PyResult<()> {
        self.with_mut(py, |cs| other.apply(py, cs))
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<DativeBondConstraintKeyIter> {
        self.read(py, |cs| dative_bond_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<DativeBondConstraintKeyIter> {
        self.read(py, |cs| dative_bond_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<DativeBondConstraintIter> {
        self.read(py, |cs| dative_bond_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<DativeBondConstraintItemsIter> {
        self.read(py, |cs| dative_bond_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_ast(py);
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| DativeBondConstraintAst::from_ast(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(into_py_variant(py, constraint)?.into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<DativeBondConstraintAst> {
        let ast_key = key.bind(py).borrow().to_ast(py);
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| DativeBondConstraintAst::from_ast(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<DativeBondConstraintKey>) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_ast(py);
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// The aromatic value; `Undetermined` when no `Aromatic` constraint is present
    /// (mirroring the non-optional Rust accessor).
    #[getter]
    fn aromatic(&self, py: Python<'_>) -> PyResult<BooleanAst> {
        self.read(py, |cs| Ok(BooleanAst::from_ast(&cs.aromatic())))
    }

    #[setter]
    fn set_aromatic(&self, py: Python<'_>, value: BooleanArg) {
        self.set_ast(py, AstDativeBondConstraintAst::aromatic(value.to_ast(py)));
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
            AstDativeBondConstraintAst::ring_membership(AstRingScope::All, value.to_ast(py)),
        );
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    fn ring_size_count(&self, py: Python<'_>) -> DativeBondRingSizeCounts {
        let backing = match &self.backing {
            DativeBondConstraintsBacking::DativeBond(bond) => {
                DativeBondRingSizeBacking::DativeBond(bond.clone_ref(py))
            }
        };
        DativeBondRingSizeCounts { backing }
    }

    /// The present constraints as a dict keyed by snake_case name.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| dative_bond_constraints_asdict(py, cs))
    }
}

/// What a `DativeBondRingSizeCounts` proxy reads/writes through to: a standalone
/// `DativeBondAst` or a standalone `DativeBondConstraintsAst` value. The
/// molecule-bond backing lands with `DativeBondView`.
enum DativeBondRingSizeBacking {
    DativeBond(Py<DativeBondAst>),
    Value(Py<DativeBondConstraintsAst>),
}

/// A subscriptable proxy over the sized-ring membership counts of a dative bond,
/// keyed by ring size: `proxy[size]` reads, `proxy[size] = count` sets, `del
/// proxy[size]` removes. Backs onto whichever container produced it (dual-backing,
/// like `DativeBondConstraintsView`).
#[pyclass]
pub struct DativeBondRingSizeCounts {
    backing: DativeBondRingSizeBacking,
}

impl DativeBondRingSizeCounts {
    /// Borrow the backing constraints and read through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstDativeBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            DativeBondRingSizeBacking::DativeBond(bond) => {
                f(&bond.bind(py).borrow().inner().constraints)
            }
            DativeBondRingSizeBacking::Value(value) => f(value.bind(py).borrow().inner()),
        }
    }

    /// Mutate the backing constraints in place through `f`.
    fn write(&self, py: Python<'_>, f: impl FnOnce(&mut AstDativeBondConstraintsAst)) {
        match &self.backing {
            DativeBondRingSizeBacking::DativeBond(bond) => {
                f(&mut bond.borrow_mut(py).inner_mut().constraints)
            }
            DativeBondRingSizeBacking::Value(value) => f(value.borrow_mut(py).inner_mut()),
        }
    }
}

#[pymethods]
impl DativeBondRingSizeCounts {
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
    fn __iter__(&self, py: Python<'_>) -> PyResult<DativeBondRingSizeIter> {
        let sizes = self.read(py, |cs| Ok(ring_sizes(cs).collect::<Vec<u8>>()))?;
        Ok(DativeBondRingSizeIter {
            sizes: sizes.into_iter(),
        })
    }

    /// Set the membership count for rings of `size` in place.
    fn __setitem__(&self, py: Python<'_>, size: u8, count: ValueArg) {
        let constraint =
            AstDativeBondConstraintAst::ring_membership(AstRingScope::Size(size), count.to_ast(py));
        self.write(py, |cs| cs.set(constraint));
    }

    /// Remove the sized-ring membership for `size` in place.
    fn __delitem__(&self, py: Python<'_>, size: u8) {
        self.write(py, |cs| {
            cs.remove(AstDativeBondConstraintKey::RingMembership(
                AstRingScope::Size(size),
            ));
        });
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |cs| {
            let mut parts = Vec::new();
            for entry in cs.iter() {
                if let AstDativeBondConstraintAst::RingMembership(m) = entry {
                    if let AstRingScope::Size(size) = m.scope {
                        let count = into_py_variant(py, ValueAst::from_ast(py, &m.count)?)?;
                        parts.push(format!(
                            "{size}: {}",
                            count.bind(py).as_any().repr()?.extract::<String>()?
                        ));
                    }
                }
            }
            Ok(format!(
                "DativeBondRingSizeCounts({{{}}})",
                parts.join(", ")
            ))
        })
    }
}

/// The ring sizes with a membership constraint, in kind-sorted order.
fn ring_sizes(constraints: &AstDativeBondConstraintsAst) -> impl Iterator<Item = u8> + '_ {
    constraints.iter().filter_map(|entry| match entry {
        AstDativeBondConstraintAst::RingMembership(m) => match m.scope {
            AstRingScope::Size(size) => Some(size),
            AstRingScope::All => None,
        },
        _ => None,
    })
}

#[pyclass]
struct DativeBondRingSizeIter {
    sizes: IntoIter<u8>,
}

#[pymethods]
impl DativeBondRingSizeIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<u8> {
        self.sizes.next()
    }
}

#[pyclass]
struct DativeBondConstraintIter {
    entries: IntoIter<Py<DativeBondConstraintAst>>,
}

#[pymethods]
impl DativeBondConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<DativeBondConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
struct DativeBondConstraintKeyIter {
    keys: IntoIter<Py<DativeBondConstraintKey>>,
}

#[pymethods]
impl DativeBondConstraintKeyIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<DativeBondConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
struct DativeBondConstraintItemsIter {
    items: IntoIter<(Py<DativeBondConstraintKey>, Py<DativeBondConstraintAst>)>,
}

#[pymethods]
impl DativeBondConstraintItemsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<(Py<DativeBondConstraintKey>, Py<DativeBondConstraintAst>)> {
        self.items.next()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{BooleanAst as AstBooleanAst, ValueAst as AstValueAst};

    use super::*;

    #[rstest]
    #[case::single("1")]
    #[case::aromatic("1#a")]
    #[case::ring_size("1#R(6)")]
    fn test_dative_bond_ast_parse(#[case] dsl: &str) {
        let bond = DativeBondAst::parse(dsl).unwrap();
        assert_eq!(bond.__str__(), dsl);
        assert_eq!(bond.__repr__(), format!("DativeBondAst.parse('{dsl}')"));
    }

    #[rstest]
    fn test_dative_bond_ast_parse_error() {
        assert!(DativeBondAst::parse("x#").is_err());
    }

    #[rstest]
    fn test_dative_bond_ast_constraints() {
        let bond = DativeBondAst(AstDativeBondAst::from_order(1).with_constraint(
            AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
        ));
        assert_eq!(bond.inner().constraints.len(), 1);
    }

    #[rstest]
    fn test_dative_bond_ast_set_constraints_from_view() {
        Python::attach(|py| {
            let src = Py::new(
                py,
                DativeBondAst::from_inner(AstDativeBondAst::from_order(1).with_constraint(
                    AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )),
            )
            .unwrap();
            let view = Py::new(
                py,
                DativeBondConstraintsView {
                    backing: DativeBondConstraintsBacking::DativeBond(src),
                },
            )
            .unwrap();
            let mut dst = DativeBondAst::from_inner(AstDativeBondAst::from_order(2));
            dst.set_constraints(py, DativeBondConstraintsArg::View(view))
                .unwrap();
            assert_eq!(dst.inner().constraints.aromatic(), AstBooleanAst::Lit(true));
        });
    }

    #[rstest]
    #[case(AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)))]
    #[case(AstDativeBondConstraintAst::ring_membership(AstRingScope::All, 2))]
    #[case(AstDativeBondConstraintAst::ring_membership(AstRingScope::Size(6), 1))]
    fn test_dative_bond_constraint_ast_roundtrip(#[case] ast: AstDativeBondConstraintAst) {
        Python::attach(|py| {
            assert_eq!(
                DativeBondConstraintAst::from_ast(py, &ast)
                    .unwrap()
                    .to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_len_contains() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::ring_membership(AstRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = DativeBondConstraintsAst::new(py, vec![aromatic, ring]);
            assert_eq!(constraints.__len__(), 2);
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap()
            ));
            assert!(!constraints.__contains__(
                py,
                into_py_variant(
                    py,
                    DativeBondConstraintKey::RingMembership(
                        into_py_variant(py, RingScope::Size(5)).unwrap()
                    ),
                )
                .unwrap()
            ));
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_keys_values_items() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::ring_membership(AstRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = DativeBondConstraintsAst::new(py, vec![aromatic, ring]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstDativeBondConstraintKey::Aromatic
            );
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstDativeBondConstraintKey::RingMembership(AstRingScope::All)
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true))
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_ast(py),
                AstDativeBondConstraintKey::Aromatic
            );
            assert_eq!(
                value.bind(py).borrow().to_ast(py),
                AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true))
            );
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_get() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = DativeBondConstraintsAst::new(py, vec![aromatic]);
            let present = constraints
                .get(
                    py,
                    into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap(),
                    None,
                )
                .unwrap();
            let expected = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap()
            .into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());
            let sentinel_key = into_py_variant(
                py,
                DativeBondConstraintKey::RingMembership(
                    into_py_variant(py, RingScope::All()).unwrap(),
                ),
            )
            .unwrap();
            let absent = constraints
                .get(py, sentinel_key.clone_ref(py), None)
                .unwrap();
            assert!(absent.bind(py).is_none());
            let sentinel = sentinel_key.into_any();
            let defaulted = constraints
                .get(
                    py,
                    into_py_variant(
                        py,
                        DativeBondConstraintKey::RingMembership(
                            into_py_variant(py, RingScope::All()).unwrap(),
                        ),
                    )
                    .unwrap(),
                    Some(sentinel.clone_ref(py)),
                )
                .unwrap();
            assert_eq!(defaulted.as_ptr(), sentinel.as_ptr());
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_aromatic() {
        Python::attach(|py| {
            let empty = DativeBondConstraintsAst::new(py, vec![]);
            assert_eq!(empty.aromatic().to_ast(), AstBooleanAst::Undetermined);
            assert!(empty.ring_count(py).unwrap().is_none());
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = DativeBondConstraintsAst::new(py, vec![aromatic]);
            assert_eq!(constraints.aromatic().to_ast(), AstBooleanAst::Lit(true));
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_ring_size_count() {
        Python::attach(|py| {
            let membership = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::ring_membership(AstRingScope::Size(6), 1),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints =
                Py::new(py, DativeBondConstraintsAst::new(py, vec![membership])).unwrap();
            let proxy = DativeBondConstraintsAst::ring_size_count(constraints.clone_ref(py));
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
    fn test_dative_bond_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = DativeBondConstraintsAst::new(py, vec![]);
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            constraints.set(py, aromatic);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(constraints.aromatic().to_ast(), AstBooleanAst::Lit(true));
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_pop() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let mut constraints = DativeBondConstraintsAst::new(py, vec![aromatic]);
            let removed = constraints
                .pop(
                    py,
                    into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(DativeBondConstraintAst::Aromatic(b)) => {
                    assert_eq!(b.bind(py).borrow().to_ast(), AstBooleanAst::Lit(true))
                }
                _ => panic!("expected removed Aromatic(Lit(true))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_update() {
        Python::attach(|py| {
            let mut constraints = DativeBondConstraintsAst::new(py, vec![]);
            let mut other = AstDativeBondConstraintsAst::new();
            other.set(AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(
                true,
            )));
            other.set(AstDativeBondConstraintAst::ring_membership(
                AstRingScope::All,
                2,
            ));
            constraints
                .update(
                    py,
                    DativeBondConstraintsUpdate::Container(
                        Py::new(py, DativeBondConstraintsAst::from_inner(other)).unwrap(),
                    ),
                )
                .unwrap();
            assert_eq!(constraints.__len__(), 2);
            assert_eq!(constraints.aromatic().to_ast(), AstBooleanAst::Lit(true));
            assert_eq!(
                constraints.ring_count(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(2)
            );
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_update_entries() {
        Python::attach(|py| {
            let mut constraints = DativeBondConstraintsAst::new(py, vec![]);
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::ring_membership(AstRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            constraints
                .update(
                    py,
                    DativeBondConstraintsUpdate::Entries(vec![aromatic, ring]),
                )
                .unwrap();
            assert_eq!(constraints.__len__(), 2);
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_set_aromatic() {
        Python::attach(|py| {
            let mut constraints = DativeBondConstraintsAst::new(py, vec![]);
            constraints.set_aromatic(py, BooleanArg::Lit(true));
            assert_eq!(constraints.aromatic().to_ast(), AstBooleanAst::Lit(true));
            constraints.set_aromatic(py, BooleanArg::Lit(false));
            assert_eq!(constraints.aromatic().to_ast(), AstBooleanAst::Lit(false));
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_set_ring_count() {
        Python::attach(|py| {
            let mut constraints = DativeBondConstraintsAst::new(py, vec![]);
            constraints.set_ring_count(py, ValueArg::Lit(2));
            assert_eq!(
                constraints.ring_count(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(2)
            );
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = DativeBondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap();
            assert!(constraints.__getitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = DativeBondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap();
            assert!(constraints.__delitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_view_set() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(AstDativeBondAst::from_order(1)),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_ast(
                    py,
                    &AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, aromatic);
            // a fresh view proves the write hit the standalone bond, not a copy
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            match fresh
                .__getitem__(
                    py,
                    into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap()
            {
                DativeBondConstraintAst::Aromatic(b) => {
                    assert_eq!(b.bind(py).borrow().to_ast(), AstBooleanAst::Lit(true))
                }
                _ => panic!("expected Aromatic(Lit(true))"),
            }
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_view_pop() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(AstDativeBondAst::from_order(1).with_constraint(
                    AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            let removed = view
                .pop(
                    py,
                    into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(DativeBondConstraintAst::Aromatic(b)) => {
                    assert_eq!(b.bind(py).borrow().to_ast(), AstBooleanAst::Lit(true))
                }
                _ => panic!("expected removed Aromatic(Lit(true))"),
            }
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_view_update() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(AstDativeBondAst::from_order(1)),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            let mut other = AstDativeBondConstraintsAst::new();
            other.set(AstDativeBondConstraintAst::aromatic(AstBooleanAst::Lit(
                true,
            )));
            other.set(AstDativeBondConstraintAst::ring_membership(
                AstRingScope::All,
                2,
            ));
            view.update(
                py,
                DativeBondConstraintsUpdate::Container(
                    Py::new(py, DativeBondConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 2);
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_view_set_aromatic() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(AstDativeBondAst::from_order(1)),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            assert_eq!(
                view.aromatic(py).unwrap().to_ast(),
                AstBooleanAst::Undetermined
            );
            view.set_aromatic(py, BooleanArg::Lit(true));
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond),
            };
            assert_eq!(
                fresh.aromatic(py).unwrap().to_ast(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_dative_bond_ring_size_counts_value_backed() {
        Python::attach(|py| {
            let constraints = Py::new(py, DativeBondConstraintsAst::new(py, vec![])).unwrap();
            let proxy = DativeBondConstraintsAst::ring_size_count(constraints.clone_ref(py));
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
    fn test_dative_bond_ring_size_counts_bond_backed() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(AstDativeBondAst::from_order(1)),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            view.ring_size_count(py)
                .__setitem__(py, 5, ValueArg::Lit(1));
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond),
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
    fn test_dative_bond_ring_size_counts_len_iter_contains() {
        Python::attach(|py| {
            let constraints = Py::new(py, DativeBondConstraintsAst::new(py, vec![])).unwrap();
            let proxy = DativeBondConstraintsAst::ring_size_count(constraints.clone_ref(py));
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
