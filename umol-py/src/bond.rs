//! Bond value type and bond-constraint surface mirroring `umol_ast::ast`:
//! `BondAst`, the `BondConstraintAst` enum, the `BondConstraintsAst` container,
//! and the `BondConstraintsView` live handle.

use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use umol_ast::ast::{
    BondAst as AstBondAst, BondConstraintAst as AstBondConstraintAst,
    BondConstraintKey as AstBondConstraintKey, BondConstraintsAst as AstBondConstraintsAst,
    RingScope as AstRingScope,
};

use crate::atom::SpinStateAst;
use crate::boolean::{BooleanArg, BooleanAst};
use crate::constraint::{RingMembershipAst, RingScope};
use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::error::parse_error;
use crate::stereo::{CisTransStereoArg, CisTransStereoAst};
use crate::value::{ValueArg, ValueAst};

/// A bond: order, charge, spin, and bond-scope constraints.
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct BondAst(AstBondAst);

#[pymethods]
impl BondAst {
    /// Construct from an order — an `int` or a `ValueAst` expression — optionally
    /// setting fields.
    #[new]
    #[pyo3(signature = (order, *, charge=None, spin=None, constraints=None))]
    fn new(
        py: Python<'_>,
        order: ValueArg,
        charge: Option<ValueArg>,
        spin: Option<PyRef<'_, SpinStateAst>>,
        constraints: Option<Py<BondConstraintsAst>>,
    ) -> Self {
        let mut bond = AstBondAst::new(order.to_ast(py));
        if let Some(charge) = charge {
            bond = bond.with_charge(charge.to_ast(py));
        }
        if let Some(spin) = spin {
            bond = bond.with_spin(spin.to_ast(py));
        }
        if let Some(constraints) = constraints {
            bond.constraints = constraints.bind(py).borrow().inner().clone();
        }
        BondAst(bond)
    }

    /// Parse a bond-DSL string (e.g. `"2#c-1"`) into a `BondAst`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        AstBondAst::from_str(s).map(Self).map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("BondAst.parse('{}')", self.0)
    }

    #[getter]
    fn order(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.order)
    }

    #[setter]
    fn set_order(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.order = value.to_ast(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.charge)
    }

    #[setter]
    fn set_charge(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.charge = value.to_ast(py);
    }

    #[getter]
    fn spin(&self, py: Python<'_>) -> PyResult<SpinStateAst> {
        SpinStateAst::from_ast(py, &self.0.spin)
    }

    #[setter]
    fn set_spin(&mut self, py: Python<'_>, value: PyRef<'_, SpinStateAst>) {
        self.0.spin = value.to_ast(py);
    }

    /// The bond's constraints as a live handle onto this bond: reads borrow the
    /// current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(slf: Py<Self>) -> BondConstraintsView {
        BondConstraintsView {
            backing: BondConstraintsBacking::Bond(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or
    /// a live view.
    #[setter]
    fn set_constraints(&mut self, py: Python<'_>, value: BondConstraintsArg) -> PyResult<()> {
        self.0.constraints = value.to_ast(py)?;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are the field mirrors.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("order", self.order(py)?)?;
        dict.set_item("charge", self.charge(py)?)?;
        dict.set_item("spin", self.spin(py)?)?;
        dict.set_item(
            "constraints",
            bond_constraints_asdict(py, &self.0.constraints)?,
        )?;
        Ok(dict)
    }
}

impl BondAst {
    /// The wrapped AST bond — read access for the bond-backed constraints view.
    pub(crate) fn inner(&self) -> &AstBondAst {
        &self.0
    }

    /// Mutable access to the wrapped AST bond — write access for the bond-backed
    /// constraints view.
    pub(crate) fn inner_mut(&mut self) -> &mut AstBondAst {
        &mut self.0
    }

    /// Wrap an AST bond (the hold-the-value `from_inner` bridge, paired with
    /// `inner`). Test-only — in-crate construction wraps `BondAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(bond: AstBondAst) -> Self {
        BondAst(bond)
    }
}

/// The key (identity) of a bond constraint, for keyed lookup. The ring-membership
/// key carries its ring scope; all other keys are the bare discriminant.
#[pyclass]
pub enum BondConstraintKey {
    Aromatic(),
    CisTransStereo(),
    RingMembership(Py<RingScope>),
}

#[pymethods]
impl BondConstraintKey {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            BondConstraintKey::Aromatic() => ("Aromatic", 0),
            BondConstraintKey::CisTransStereo() => ("CisTransStereo", 0),
            BondConstraintKey::RingMembership(_) => ("RingMembership", 1),
        };
        variant_repr(slf.bind(py).as_any(), "BondConstraintKey", variant, arity)
    }
}

impl BondConstraintKey {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstBondConstraintKey) -> PyResult<Self> {
        Ok(match ast {
            AstBondConstraintKey::Aromatic => Self::Aromatic(),
            AstBondConstraintKey::CisTransStereo => Self::CisTransStereo(),
            AstBondConstraintKey::RingMembership(scope) => {
                Self::RingMembership(into_py_variant(py, RingScope::from_ast(scope))?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstBondConstraintKey {
        match self {
            Self::Aromatic() => AstBondConstraintKey::Aromatic,
            Self::CisTransStereo() => AstBondConstraintKey::CisTransStereo,
            Self::RingMembership(scope) => {
                AstBondConstraintKey::RingMembership(scope.bind(py).borrow().to_ast())
            }
        }
    }
}

/// A bond-scope constraint: the aromatic flag, cis/trans stereo, or a ring
/// membership of a single bond.
#[pyclass]
pub enum BondConstraintAst {
    Aromatic(Py<BooleanAst>),
    CisTransStereo(Py<CisTransStereoAst>),
    RingMembership(Py<RingMembershipAst>),
}

#[pymethods]
impl BondConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    fn key(&self, py: Python<'_>) -> PyResult<BondConstraintKey> {
        BondConstraintKey::from_ast(py, &self.to_ast(py).key())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            BondConstraintAst::Aromatic(_) => "Aromatic",
            BondConstraintAst::CisTransStereo(_) => "CisTransStereo",
            BondConstraintAst::RingMembership(_) => "RingMembership",
        };
        variant_repr(slf.bind(py).as_any(), "BondConstraintAst", variant, 1)
    }
}

impl BondConstraintAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstBondConstraintAst) -> PyResult<Self> {
        Ok(match ast {
            AstBondConstraintAst::Aromatic(b) => {
                Self::Aromatic(into_py_variant(py, BooleanAst::from_ast(b))?)
            }
            AstBondConstraintAst::CisTransStereo(c) => Self::CisTransStereo(into_py_variant(
                py,
                CisTransStereoAst::from_ast(py, c)?,
            )?),
            AstBondConstraintAst::RingMembership(m) => {
                Self::RingMembership(into_py_variant(py, RingMembershipAst::from_ast(py, m)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstBondConstraintAst {
        match self {
            Self::Aromatic(b) => AstBondConstraintAst::Aromatic(b.bind(py).borrow().to_ast()),
            Self::CisTransStereo(c) => {
                AstBondConstraintAst::CisTransStereo(c.bind(py).borrow().to_ast(py))
            }
            Self::RingMembership(m) => {
                AstBondConstraintAst::RingMembership(m.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container (value or live view) or
/// an iterable of `BondConstraintAst` (each `set`, last-wins).
#[derive(FromPyObject)]
enum BondConstraintsUpdate {
    Container(Py<BondConstraintsAst>),
    View(Py<BondConstraintsView>),
    Entries(Vec<Py<BondConstraintAst>>),
}

impl BondConstraintsUpdate {
    /// Overlay this update onto `target` in place.
    fn apply(&self, py: Python<'_>, target: &mut AstBondConstraintsAst) -> PyResult<()> {
        match self {
            BondConstraintsUpdate::Container(c) => target.update(c.bind(py).borrow().inner()),
            BondConstraintsUpdate::View(v) => {
                let snapshot = v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?;
                target.update(&snapshot);
            }
            BondConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry.bind(py).borrow().to_ast(py));
                }
            }
        }
        Ok(())
    }
}

/// A whole-container argument that snapshots either a value container or a live
/// view — for the bond `constraints` setter, which accepts either.
#[derive(FromPyObject)]
enum BondConstraintsArg {
    Container(Py<BondConstraintsAst>),
    View(Py<BondConstraintsView>),
}

impl BondConstraintsArg {
    fn to_ast(&self, py: Python<'_>) -> PyResult<AstBondConstraintsAst> {
        match self {
            BondConstraintsArg::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
            BondConstraintsArg::View(v) => v.bind(py).borrow().read(py, |cs| Ok(cs.clone())),
        }
    }
}

/// The bond-scope constraints on a bond, in kind-sorted order. Mutable, hence
/// value-equal but unhashable (matching `BondAst`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct BondConstraintsAst(AstBondConstraintsAst);

#[pymethods]
impl BondConstraintsAst {
    /// Build from a sequence of constraints (kind-sorted; a unique kind replaces
    /// an earlier one, ring memberships accumulate per scope).
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<BondConstraintAst>>) -> Self {
        let mut constraints = AstBondConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_ast(py)),
        );
        BondConstraintsAst(constraints)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let mirror = into_py_variant(py, BondConstraintAst::from_ast(py, entry)?)?;
            parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("BondConstraintsAst([{}])", parts.join(", ")))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    fn set(&mut self, py: Python<'_>, c: Py<BondConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
    ) -> PyResult<Option<BondConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_ast(py))
            .map(|c| BondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container or an iterable of
    /// `BondConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&mut self, py: Python<'_>, other: BondConstraintsUpdate) -> PyResult<()> {
        other.apply(py, &mut self.0)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<BondConstraintKeyIter> {
        bond_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<BondConstraintKeyIter> {
        bond_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<BondConstraintIter> {
        bond_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<BondConstraintItemsIter> {
        bond_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_ast(py)) {
            Some(constraint) => {
                Ok(into_py_variant(py, BondConstraintAst::from_ast(py, constraint)?)?.into_any())
            }
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(&self, py: Python<'_>, key: Py<BondConstraintKey>) -> PyResult<BondConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_ast(py)) {
            Some(constraint) => BondConstraintAst::from_ast(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&mut self, py: Python<'_>, key: Py<BondConstraintKey>) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_ast(py)).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<BondConstraintKey>) -> bool {
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
        self.0.set(AstBondConstraintAst::aromatic(value.to_ast(py)));
    }

    /// The cis/trans-stereo state, or `None`.
    #[getter]
    fn cis_trans_stereo(&self, py: Python<'_>) -> PyResult<Option<CisTransStereoAst>> {
        self.0
            .cis_trans_stereo()
            .map(|c| CisTransStereoAst::from_ast(py, c))
            .transpose()
    }

    #[setter]
    fn set_cis_trans_stereo(&mut self, py: Python<'_>, value: CisTransStereoArg) -> PyResult<()> {
        self.0
            .set(AstBondConstraintAst::cis_trans_stereo(value.to_ast(py)?));
        Ok(())
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
        self.0.set(AstBondConstraintAst::ring_membership(
            AstRingScope::All,
            value.to_ast(py),
        ));
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    fn ring_size_count(slf: Py<Self>) -> BondRingSizeCounts {
        BondRingSizeCounts {
            backing: BondRingSizeBacking::Value(slf),
        }
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors. Ring memberships key by scope: `ring_count` for the
    /// all-rings scope, `ring_size_count_<n>` for a specific ring size.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        bond_constraints_asdict(py, &self.0)
    }
}

impl BondConstraintsAst {
    /// The wrapped AST constraints — read access for bond construction.
    pub(crate) fn inner(&self) -> &AstBondConstraintsAst {
        &self.0
    }

    /// Mutable access to the wrapped AST constraints — for the value-backed proxy.
    pub(crate) fn inner_mut(&mut self) -> &mut AstBondConstraintsAst {
        &mut self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge). Test-only —
    /// in-crate construction wraps `BondConstraintsAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(constraints: AstBondConstraintsAst) -> Self {
        BondConstraintsAst(constraints)
    }
}

/// Build the per-constraint iterator handle from a borrowed container.
fn bond_constraints_iter(
    py: Python<'_>,
    constraints: &AstBondConstraintsAst,
) -> PyResult<BondConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| into_py_variant(py, BondConstraintAst::from_ast(py, constraint)?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(BondConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
fn bond_constraint_keys(
    py: Python<'_>,
    constraints: &AstBondConstraintsAst,
) -> PyResult<BondConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| into_py_variant(py, BondConstraintKey::from_ast(py, &constraint.key())?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(BondConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
fn bond_constraint_items(
    py: Python<'_>,
    constraints: &AstBondConstraintsAst,
) -> PyResult<BondConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(py, BondConstraintKey::from_ast(py, &constraint.key())?)?,
                into_py_variant(py, BondConstraintAst::from_ast(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(BondConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// inner-value mirrors. Ring memberships key by scope: `ring_count` for the
/// all-rings scope, `ring_size_count_<n>` for a specific ring size.
fn bond_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &AstBondConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            AstBondConstraintAst::Aromatic(b) => {
                dict.set_item("aromatic", BooleanAst::from_ast(b))?
            }
            AstBondConstraintAst::CisTransStereo(c) => {
                dict.set_item("cis_trans_stereo", CisTransStereoAst::from_ast(py, c)?)?
            }
            AstBondConstraintAst::RingMembership(m) => {
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

/// What a `BondConstraintsView` writes through to: a standalone `BondAst`. (The
/// molecule-bond arm arrives with the `BondView` half.)
enum BondConstraintsBacking {
    Bond(Py<BondAst>),
}

/// A live handle onto one bond's constraints, backed by a standalone `BondAst`.
/// Reads borrow the bond's constraints and read only the item they need (no
/// whole-container clone); mutators write through to the bond in place, without a
/// clone-and-writeback.
#[pyclass]
pub struct BondConstraintsView {
    backing: BondConstraintsBacking,
}

impl BondConstraintsView {
    /// Borrow the backing bond's constraints and read one item through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            BondConstraintsBacking::Bond(bond) => {
                let bond = bond.bind(py).borrow();
                f(&bond.inner().constraints)
            }
        }
    }

    /// Mutate the backing bond's constraints in place through `f`.
    fn with_mut<R>(&self, py: Python<'_>, f: impl FnOnce(&mut AstBondConstraintsAst) -> R) -> R {
        match &self.backing {
            BondConstraintsBacking::Bond(bond) => {
                f(&mut bond.borrow_mut(py).inner_mut().constraints)
            }
        }
    }

    /// Set one constraint on the backing bond in place (last-wins per key).
    fn set_ast(&self, py: Python<'_>, constraint: AstBondConstraintAst) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing bond in place, returning the removed entry.
    fn remove_ast(
        &self,
        py: Python<'_>,
        key: AstBondConstraintKey,
    ) -> Option<AstBondConstraintAst> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl BondConstraintsView {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("BondConstraintsView({count} entries)"))
    }

    /// Insert `c` on the bond in place, replacing any existing entry of the same
    /// key (last-wins).
    fn set(&self, py: Python<'_>, c: Py<BondConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key from the bond in place, returning it if
    /// present (dict `pop`).
    fn pop(
        &self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
    ) -> PyResult<Option<BondConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_ast(py))
            .map(|c| BondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&self, py: Python<'_>, key: Py<BondConstraintKey>) -> PyResult<()> {
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
    /// iterable of `BondConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&self, py: Python<'_>, other: BondConstraintsUpdate) -> PyResult<()> {
        self.with_mut(py, |cs| other.apply(py, cs))
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<BondConstraintKeyIter> {
        self.read(py, |cs| bond_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<BondConstraintKeyIter> {
        self.read(py, |cs| bond_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<BondConstraintIter> {
        self.read(py, |cs| bond_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<BondConstraintItemsIter> {
        self.read(py, |cs| bond_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_ast(py);
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| BondConstraintAst::from_ast(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(into_py_variant(py, constraint)?.into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(&self, py: Python<'_>, key: Py<BondConstraintKey>) -> PyResult<BondConstraintAst> {
        let ast_key = key.bind(py).borrow().to_ast(py);
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| BondConstraintAst::from_ast(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<BondConstraintKey>) -> PyResult<bool> {
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
        self.set_ast(py, AstBondConstraintAst::aromatic(value.to_ast(py)));
    }

    /// The cis/trans-stereo state, or `None`.
    #[getter]
    fn cis_trans_stereo(&self, py: Python<'_>) -> PyResult<Option<CisTransStereoAst>> {
        self.read(py, |cs| {
            cs.cis_trans_stereo()
                .map(|c| CisTransStereoAst::from_ast(py, c))
                .transpose()
        })
    }

    #[setter]
    fn set_cis_trans_stereo(&self, py: Python<'_>, value: CisTransStereoArg) -> PyResult<()> {
        self.set_ast(py, AstBondConstraintAst::cis_trans_stereo(value.to_ast(py)?));
        Ok(())
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
            AstBondConstraintAst::ring_membership(AstRingScope::All, value.to_ast(py)),
        );
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    fn ring_size_count(&self, py: Python<'_>) -> BondRingSizeCounts {
        let backing = match &self.backing {
            BondConstraintsBacking::Bond(bond) => BondRingSizeBacking::Bond(bond.clone_ref(py)),
        };
        BondRingSizeCounts { backing }
    }

    /// The present constraints as a dict keyed by snake_case name.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| bond_constraints_asdict(py, cs))
    }
}

/// What a `BondRingSizeCounts` proxy reads/writes through to: a standalone
/// `BondAst` or a standalone `BondConstraintsAst` value.
enum BondRingSizeBacking {
    Bond(Py<BondAst>),
    Value(Py<BondConstraintsAst>),
}

/// A subscriptable proxy over the sized-ring membership counts of a bond, keyed by
/// ring size: `proxy[size]` reads, `proxy[size] = count` sets, `del proxy[size]`
/// removes. Backs onto whichever container produced it (dual-backing, like
/// `BondConstraintsView`).
#[pyclass]
pub struct BondRingSizeCounts {
    backing: BondRingSizeBacking,
}

impl BondRingSizeCounts {
    /// Borrow the backing constraints and read through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            BondRingSizeBacking::Bond(bond) => f(&bond.bind(py).borrow().inner().constraints),
            BondRingSizeBacking::Value(value) => f(value.bind(py).borrow().inner()),
        }
    }

    /// Mutate the backing constraints in place through `f`.
    fn write(&self, py: Python<'_>, f: impl FnOnce(&mut AstBondConstraintsAst)) {
        match &self.backing {
            BondRingSizeBacking::Bond(bond) => f(&mut bond.borrow_mut(py).inner_mut().constraints),
            BondRingSizeBacking::Value(value) => f(value.borrow_mut(py).inner_mut()),
        }
    }
}

#[pymethods]
impl BondRingSizeCounts {
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
    fn __iter__(&self, py: Python<'_>) -> PyResult<BondRingSizeIter> {
        let sizes = self.read(py, |cs| Ok(ring_sizes(cs).collect::<Vec<u8>>()))?;
        Ok(BondRingSizeIter {
            sizes: sizes.into_iter(),
        })
    }

    /// Set the membership count for rings of `size` in place.
    fn __setitem__(&self, py: Python<'_>, size: u8, count: ValueArg) {
        let constraint =
            AstBondConstraintAst::ring_membership(AstRingScope::Size(size), count.to_ast(py));
        self.write(py, |cs| cs.set(constraint));
    }

    /// Remove the sized-ring membership for `size` in place.
    fn __delitem__(&self, py: Python<'_>, size: u8) {
        self.write(py, |cs| {
            cs.remove(AstBondConstraintKey::RingMembership(AstRingScope::Size(size)));
        });
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |cs| {
            let mut parts = Vec::new();
            for entry in cs.iter() {
                if let AstBondConstraintAst::RingMembership(m) = entry {
                    if let AstRingScope::Size(size) = m.scope {
                        let count = into_py_variant(py, ValueAst::from_ast(py, &m.count)?)?;
                        parts.push(format!(
                            "{size}: {}",
                            count.bind(py).as_any().repr()?.extract::<String>()?
                        ));
                    }
                }
            }
            Ok(format!("BondRingSizeCounts({{{}}})", parts.join(", ")))
        })
    }
}

/// The ring sizes with a membership constraint, in kind-sorted order.
fn ring_sizes(constraints: &AstBondConstraintsAst) -> impl Iterator<Item = u8> + '_ {
    constraints.iter().filter_map(|entry| match entry {
        AstBondConstraintAst::RingMembership(m) => match m.scope {
            AstRingScope::Size(size) => Some(size),
            AstRingScope::All => None,
        },
        _ => None,
    })
}

#[pyclass]
struct BondRingSizeIter {
    sizes: IntoIter<u8>,
}

#[pymethods]
impl BondRingSizeIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<u8> {
        self.sizes.next()
    }
}

#[pyclass]
struct BondConstraintIter {
    entries: IntoIter<Py<BondConstraintAst>>,
}

#[pymethods]
impl BondConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<BondConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
struct BondConstraintKeyIter {
    keys: IntoIter<Py<BondConstraintKey>>,
}

#[pymethods]
impl BondConstraintKeyIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<BondConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
struct BondConstraintItemsIter {
    items: IntoIter<(Py<BondConstraintKey>, Py<BondConstraintAst>)>,
}

#[pymethods]
impl BondConstraintItemsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<(Py<BondConstraintKey>, Py<BondConstraintAst>)> {
        self.items.next()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        BooleanAst as AstBooleanAst, CisTransStereoAst as AstCisTransStereoAst,
        StereoCosetAst as AstStereoCosetAst, ValueAst as AstValueAst,
    };

    use super::*;
    use crate::stereo::CisTransStereo;

    #[rstest]
    #[case::single("1")]
    #[case::charge("2#c-")]
    #[case::aromatic("1#a")]
    #[case::ring_size("1#R(6)")]
    fn test_bond_ast_parse(#[case] dsl: &str) {
        let bond = BondAst::parse(dsl).unwrap();
        assert_eq!(bond.__str__(), dsl);
        assert_eq!(bond.__repr__(), format!("BondAst.parse('{dsl}')"));
    }

    #[rstest]
    fn test_bond_ast_parse_error() {
        assert!(BondAst::parse("x#").is_err());
    }

    #[rstest]
    fn test_bond_ast_constraints() {
        let bond = BondAst(
            AstBondAst::from_order(1)
                .with_constraint(AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true))),
        );
        assert_eq!(bond.inner().constraints.len(), 1);
    }

    #[rstest]
    fn test_bond_ast_set_constraints_from_view() {
        Python::attach(|py| {
            let src = Py::new(
                py,
                BondAst::from_inner(
                    AstBondAst::from_order(1)
                        .with_constraint(AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true))),
                ),
            )
            .unwrap();
            let view = Py::new(
                py,
                BondConstraintsView {
                    backing: BondConstraintsBacking::Bond(src),
                },
            )
            .unwrap();
            let mut dst = BondAst::from_inner(AstBondAst::from_order(2));
            dst.set_constraints(py, BondConstraintsArg::View(view))
                .unwrap();
            assert_eq!(dst.inner().constraints.aromatic(), AstBooleanAst::Lit(true));
        });
    }

    #[rstest]
    #[case(AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)))]
    #[case(AstBondConstraintAst::cis_trans_stereo(AstCisTransStereoAst::NotStereo))]
    #[case(AstBondConstraintAst::cis_trans_stereo(AstCisTransStereoAst::Stereo(AstStereoCosetAst::Lit(1))))]
    #[case(AstBondConstraintAst::ring_membership(AstRingScope::All, 2))]
    #[case(AstBondConstraintAst::ring_membership(AstRingScope::Size(6), 1))]
    fn test_bond_constraint_ast_roundtrip(#[case] ast: AstBondConstraintAst) {
        Python::attach(|py| {
            assert_eq!(
                BondConstraintAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_len_contains() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::ring_membership(AstRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = BondConstraintsAst::new(py, vec![aromatic, ring]);
            assert_eq!(constraints.__len__(), 2);
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, BondConstraintKey::Aromatic()).unwrap()
            ));
            assert!(!constraints.__contains__(
                py,
                into_py_variant(py, BondConstraintKey::CisTransStereo()).unwrap()
            ));
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_keys_values_items() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::ring_membership(AstRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = BondConstraintsAst::new(py, vec![aromatic, ring]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstBondConstraintKey::Aromatic
            );
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstBondConstraintKey::RingMembership(AstRingScope::All)
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true))
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_ast(py),
                AstBondConstraintKey::Aromatic
            );
            assert_eq!(
                value.bind(py).borrow().to_ast(py),
                AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true))
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_get() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = BondConstraintsAst::new(py, vec![aromatic]);
            let present = constraints
                .get(
                    py,
                    into_py_variant(py, BondConstraintKey::Aromatic()).unwrap(),
                    None,
                )
                .unwrap();
            let expected = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap()
            .into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());
            let absent = constraints
                .get(
                    py,
                    into_py_variant(py, BondConstraintKey::CisTransStereo()).unwrap(),
                    None,
                )
                .unwrap();
            assert!(absent.bind(py).is_none());
            let sentinel = into_py_variant(py, BondConstraintKey::CisTransStereo())
                .unwrap()
                .into_any();
            let defaulted = constraints
                .get(
                    py,
                    into_py_variant(py, BondConstraintKey::CisTransStereo()).unwrap(),
                    Some(sentinel.clone_ref(py)),
                )
                .unwrap();
            assert_eq!(defaulted.as_ptr(), sentinel.as_ptr());
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_aromatic() {
        Python::attach(|py| {
            let empty = BondConstraintsAst::new(py, vec![]);
            assert_eq!(empty.aromatic().to_ast(), AstBooleanAst::Undetermined);
            assert!(empty.cis_trans_stereo(py).unwrap().is_none());
            assert!(empty.ring_count(py).unwrap().is_none());
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = BondConstraintsAst::new(py, vec![aromatic]);
            assert_eq!(constraints.aromatic().to_ast(), AstBooleanAst::Lit(true));
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_ring_size_count() {
        Python::attach(|py| {
            let membership = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::ring_membership(AstRingScope::Size(6), 1),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = Py::new(py, BondConstraintsAst::new(py, vec![membership])).unwrap();
            let proxy = BondConstraintsAst::ring_size_count(constraints.clone_ref(py));
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
    fn test_bond_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
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
    fn test_bond_constraints_ast_pop() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let mut constraints = BondConstraintsAst::new(py, vec![aromatic]);
            let removed = constraints
                .pop(
                    py,
                    into_py_variant(py, BondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(BondConstraintAst::Aromatic(b)) => {
                    assert_eq!(b.bind(py).borrow().to_ast(), AstBooleanAst::Lit(true))
                }
                _ => panic!("expected removed Aromatic(Lit(true))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_update() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            let mut other = AstBondConstraintsAst::new();
            other.set(AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)));
            other.set(AstBondConstraintAst::ring_membership(AstRingScope::All, 2));
            constraints
                .update(
                    py,
                    BondConstraintsUpdate::Container(
                        Py::new(py, BondConstraintsAst::from_inner(other)).unwrap(),
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
    fn test_bond_constraints_ast_update_entries() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::ring_membership(AstRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            constraints
                .update(py, BondConstraintsUpdate::Entries(vec![aromatic, ring]))
                .unwrap();
            assert_eq!(constraints.__len__(), 2);
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_set_aromatic() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            constraints.set_aromatic(py, BooleanArg::Lit(true));
            assert_eq!(constraints.aromatic().to_ast(), AstBooleanAst::Lit(true));
            constraints.set_aromatic(py, BooleanArg::Lit(false));
            assert_eq!(constraints.aromatic().to_ast(), AstBooleanAst::Lit(false));
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_set_cis_trans_stereo() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            constraints
                .set_cis_trans_stereo(py, CisTransStereoArg::Config(CisTransStereo::E))
                .unwrap();
            match constraints.cis_trans_stereo(py).unwrap().unwrap() {
                CisTransStereoAst::Stereo(coset) => {
                    assert_eq!(coset.bind(py).borrow().to_ast(py), AstStereoCosetAst::Lit(1))
                }
                _ => panic!("expected Stereo"),
            }
            constraints
                .set_cis_trans_stereo(py, CisTransStereoArg::Flag(false))
                .unwrap();
            match constraints.cis_trans_stereo(py).unwrap().unwrap() {
                CisTransStereoAst::NotStereo() => {}
                _ => panic!("expected NotStereo"),
            }
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_set_cis_trans_stereo_error() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            assert!(constraints
                .set_cis_trans_stereo(py, CisTransStereoArg::Flag(true))
                .is_err());
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_set_ring_count() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            constraints.set_ring_count(py, ValueArg::Lit(2));
            assert_eq!(
                constraints.ring_count(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(2)
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = BondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, BondConstraintKey::Aromatic()).unwrap();
            assert!(constraints.__getitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, BondConstraintKey::Aromatic()).unwrap();
            assert!(constraints.__delitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_bond_constraints_view_set() {
        Python::attach(|py| {
            let bond = Py::new(py, BondAst::from_inner(AstBondAst::from_order(1))).unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_ast(
                    py,
                    &AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, aromatic);
            // a fresh view proves the write hit the standalone bond, not a copy
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            match fresh
                .__getitem__(
                    py,
                    into_py_variant(py, BondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap()
            {
                BondConstraintAst::Aromatic(b) => {
                    assert_eq!(b.bind(py).borrow().to_ast(), AstBooleanAst::Lit(true))
                }
                _ => panic!("expected Aromatic(Lit(true))"),
            }
        });
    }

    #[rstest]
    fn test_bond_constraints_view_pop() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                BondAst::from_inner(
                    AstBondAst::from_order(1)
                        .with_constraint(AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true))),
                ),
            )
            .unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            let removed = view
                .pop(
                    py,
                    into_py_variant(py, BondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(BondConstraintAst::Aromatic(b)) => {
                    assert_eq!(b.bind(py).borrow().to_ast(), AstBooleanAst::Lit(true))
                }
                _ => panic!("expected removed Aromatic(Lit(true))"),
            }
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_bond_constraints_view_update() {
        Python::attach(|py| {
            let bond = Py::new(py, BondAst::from_inner(AstBondAst::from_order(1))).unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            let mut other = AstBondConstraintsAst::new();
            other.set(AstBondConstraintAst::aromatic(AstBooleanAst::Lit(true)));
            other.set(AstBondConstraintAst::ring_membership(AstRingScope::All, 2));
            view.update(
                py,
                BondConstraintsUpdate::Container(
                    Py::new(py, BondConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 2);
        });
    }

    #[rstest]
    fn test_bond_constraints_view_set_aromatic() {
        Python::attach(|py| {
            let bond = Py::new(py, BondAst::from_inner(AstBondAst::from_order(1))).unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            assert_eq!(
                view.aromatic(py).unwrap().to_ast(),
                AstBooleanAst::Undetermined
            );
            view.set_aromatic(py, BooleanArg::Lit(true));
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond),
            };
            assert_eq!(fresh.aromatic(py).unwrap().to_ast(), AstBooleanAst::Lit(true));
        });
    }

    #[rstest]
    fn test_bond_ring_size_counts_value_backed() {
        Python::attach(|py| {
            let constraints = Py::new(py, BondConstraintsAst::new(py, vec![])).unwrap();
            let proxy = BondConstraintsAst::ring_size_count(constraints.clone_ref(py));
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
    fn test_bond_ring_size_counts_bond_backed() {
        Python::attach(|py| {
            let bond = Py::new(py, BondAst::from_inner(AstBondAst::from_order(1))).unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            view.ring_size_count(py).__setitem__(py, 5, ValueArg::Lit(1));
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond),
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
    fn test_bond_ring_size_counts_len_iter_contains() {
        Python::attach(|py| {
            let constraints = Py::new(py, BondConstraintsAst::new(py, vec![])).unwrap();
            let proxy = BondConstraintsAst::ring_size_count(constraints.clone_ref(py));
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
