//! Dative-bond constraint values, containers, and live views.

use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use umol_graph_ir::ir::{
    DativeBondConstraintAst as GraphIrDativeBondConstraintAst,
    DativeBondConstraintKey as GraphIrDativeBondConstraintKey,
    DativeBondConstraintsAst as GraphIrDativeBondConstraintsAst,
    DativeBondId as GraphIrDativeBondId, RingScope as GraphIrRingScope,
};

use super::ring::{RingMembershipAst, RingScope};
use crate::boolean::{BooleanAst, BooleanLike};
use crate::convert::{hash_rust, into_py_variant, variant_repr};
use crate::dative::DativeBondAst;
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;
use crate::value::{NumLike, ValueAst};

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
    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
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
    pub(crate) fn from_rust(
        py: Python<'_>,
        ast: &GraphIrDativeBondConstraintKey,
    ) -> PyResult<Self> {
        Ok(match ast {
            GraphIrDativeBondConstraintKey::Aromatic => Self::Aromatic(),
            GraphIrDativeBondConstraintKey::RingMembership(scope) => {
                Self::RingMembership(into_py_variant(py, RingScope::from_rust(scope))?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrDativeBondConstraintKey {
        match self {
            Self::Aromatic() => GraphIrDativeBondConstraintKey::Aromatic,
            Self::RingMembership(scope) => {
                GraphIrDativeBondConstraintKey::RingMembership(scope.bind(py).borrow().to_rust())
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
    pub(crate) fn key(&self, py: Python<'_>) -> PyResult<DativeBondConstraintKey> {
        DativeBondConstraintKey::from_rust(py, &self.to_rust(py).key())
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            DativeBondConstraintAst::Aromatic(_) => "Aromatic",
            DativeBondConstraintAst::RingMembership(_) => "RingMembership",
        };
        variant_repr(slf.bind(py).as_any(), "DativeBondConstraintAst", variant, 1)
    }
}

impl_py_lattice!(
    DativeBondConstraintAst,
    GraphIrDativeBondConstraintAst,
    |value: &DativeBondConstraintAst, py: Python<'_>| -> PyResult<GraphIrDativeBondConstraintAst> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrDativeBondConstraintAst| -> PyResult<DativeBondConstraintAst> {
        DativeBondConstraintAst::from_rust(py, &value)
    }
);

impl DativeBondConstraintAst {
    pub(crate) fn from_rust(
        py: Python<'_>,
        ast: &GraphIrDativeBondConstraintAst,
    ) -> PyResult<Self> {
        Ok(match ast {
            GraphIrDativeBondConstraintAst::Aromatic(b) => {
                Self::Aromatic(into_py_variant(py, BooleanAst::from_rust(b))?)
            }
            GraphIrDativeBondConstraintAst::RingMembership(m) => {
                Self::RingMembership(into_py_variant(py, RingMembershipAst::from_rust(py, m)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrDativeBondConstraintAst {
        match self {
            Self::Aromatic(b) => {
                GraphIrDativeBondConstraintAst::Aromatic(b.bind(py).borrow().to_rust())
            }
            Self::RingMembership(m) => {
                GraphIrDativeBondConstraintAst::RingMembership(m.bind(py).borrow().to_rust(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container (value or live view) or
/// an iterable of `DativeBondConstraintAst` (each `set`, last-wins).
#[derive(FromPyObject)]
pub(crate) enum DativeBondConstraintsUpdate {
    Container(Py<DativeBondConstraintsAst>),
    View(Py<DativeBondConstraintsView>),
    Entries(Vec<Py<DativeBondConstraintAst>>),
}

impl DativeBondConstraintsUpdate {
    /// Read every Python object into owned data — no write target is touched. Callers
    /// resolve *before* taking the write borrow so a view (or container) that aliases
    /// the same bond is read while nothing is borrowed (otherwise
    /// `bond.constraints.update(bond.constraints)` self-aliases into a double-borrow panic).
    pub(crate) fn resolve(&self, py: Python<'_>) -> PyResult<ResolvedDativeBondConstraintsUpdate> {
        Ok(match self {
            DativeBondConstraintsUpdate::Container(c) => {
                ResolvedDativeBondConstraintsUpdate::Overlay(c.bind(py).borrow().inner().clone())
            }
            DativeBondConstraintsUpdate::View(v) => ResolvedDativeBondConstraintsUpdate::Overlay(
                v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?,
            ),
            DativeBondConstraintsUpdate::Entries(entries) => {
                ResolvedDativeBondConstraintsUpdate::Entries(
                    entries
                        .iter()
                        .map(|entry| entry.bind(py).borrow().to_rust(py))
                        .collect(),
                )
            }
        })
    }
}

/// A `DativeBondConstraintsUpdate` with all Python-object reads already done, so it can
/// be applied under a write borrow without re-entering Python.
pub(crate) enum ResolvedDativeBondConstraintsUpdate {
    /// A whole container (from another container or a live view): overlaid via `update`
    /// (last-wins per key; undetermined entries remove).
    Overlay(GraphIrDativeBondConstraintsAst),
    /// Loose entries: `set` each (last-wins; undetermined entries stored, not removed).
    Entries(Vec<GraphIrDativeBondConstraintAst>),
}

impl ResolvedDativeBondConstraintsUpdate {
    /// Overlay onto `target` in place. No Python reads.
    pub(crate) fn apply(self, target: &mut GraphIrDativeBondConstraintsAst) {
        match self {
            ResolvedDativeBondConstraintsUpdate::Overlay(overlay) => target.update(&overlay),
            ResolvedDativeBondConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry);
                }
            }
        }
    }
}

/// A whole-container argument that snapshots either a value container or a live
/// view — for the dative bond `constraints` setter, which accepts either.
#[derive(FromPyObject)]
pub(crate) enum DativeBondConstraintsLike {
    Container(Py<DativeBondConstraintsAst>),
    View(Py<DativeBondConstraintsView>),
}

impl DativeBondConstraintsLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<GraphIrDativeBondConstraintsAst> {
        match self {
            DativeBondConstraintsLike::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
            DativeBondConstraintsLike::View(v) => v.bind(py).borrow().read(py, |cs| Ok(cs.clone())),
        }
    }
}

/// The dative-bond-scope constraints on a dative bond, in kind-sorted order.
/// Mutable, hence value-equal but unhashable (matching `DativeBondAst`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct DativeBondConstraintsAst(GraphIrDativeBondConstraintsAst);

#[pymethods]
impl DativeBondConstraintsAst {
    /// Build from a sequence of constraints (kind-sorted; a unique kind replaces
    /// an earlier one, ring memberships accumulate per scope).
    #[new]
    pub(crate) fn new(py: Python<'_>, entries: Vec<Py<DativeBondConstraintAst>>) -> Self {
        let mut constraints = GraphIrDativeBondConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py)),
        );
        DativeBondConstraintsAst(constraints)
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let value = into_py_variant(py, DativeBondConstraintAst::from_rust(py, entry)?)?;
            parts.push(value.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("DativeBondConstraintsAst([{}])", parts.join(", ")))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    pub(crate) fn set(&mut self, py: Python<'_>, c: Py<DativeBondConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_rust(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    pub(crate) fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<Option<DativeBondConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_rust(py))
            .map(|c| DativeBondConstraintAst::from_rust(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container, a live view, or an
    /// iterable of `DativeBondConstraintAst` (last-wins per key; undetermined entries
    /// remove). Takes `slf` by handle so `other` is fully read *before* the write borrow
    /// — `cs.update(cs)` on the same container is then a no-op, not a double-borrow panic.
    pub(crate) fn update(
        slf: Py<Self>,
        py: Python<'_>,
        other: DativeBondConstraintsUpdate,
    ) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        resolved.apply(&mut slf.borrow_mut(py).0);
        Ok(())
    }

    pub(crate) fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<DativeBondConstraintKeyIter> {
        dative_bond_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<DativeBondConstraintKeyIter> {
        dative_bond_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<DativeBondConstraintIter> {
        dative_bond_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<DativeBondConstraintItemsIter> {
        dative_bond_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
        &self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_rust(py)) {
            Some(constraint) => Ok(into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(py, constraint)?,
            )?
            .into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    pub(crate) fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<DativeBondConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_rust(py)) {
            Some(constraint) => DativeBondConstraintAst::from_rust(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_rust(py)).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    pub(crate) fn __contains__(&self, py: Python<'_>, key: Py<DativeBondConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_rust(py))
    }

    /// The aromatic value; `Undetermined` when no `Aromatic` constraint is present
    /// (matching the non-optional Rust accessor).
    #[getter]
    pub(crate) fn aromatic(&self) -> BooleanAst {
        BooleanAst::from_rust(&self.0.aromatic())
    }

    #[setter]
    pub(crate) fn set_aromatic(&mut self, py: Python<'_>, value: BooleanLike) {
        self.0
            .set(GraphIrDativeBondConstraintAst::aromatic(value.to_rust(py)));
    }

    /// The all-rings membership count, or `None`.
    #[getter]
    pub(crate) fn ring_count(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_count()
            .map(|v| ValueAst::from_rust(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_ring_count(&mut self, py: Python<'_>, value: NumLike) {
        self.0.set(GraphIrDativeBondConstraintAst::ring_membership(
            GraphIrRingScope::All,
            value.to_rust(py),
        ));
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    pub(crate) fn ring_size_count(slf: Py<Self>) -> DativeBondRingSizeCounts {
        DativeBondRingSizeCounts {
            backing: DativeBondRingSizeBacking::Value(slf),
        }
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// Python values. Ring memberships key by scope: `ring_count` for the
    /// all-rings scope, `ring_size_count_<n>` for a specific ring size.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        dative_bond_constraints_asdict(py, &self.0)
    }
}

impl DativeBondConstraintsAst {
    /// The wrapped AST constraints — read access for dative bond construction.
    pub(crate) fn inner(&self) -> &GraphIrDativeBondConstraintsAst {
        &self.0
    }

    /// Mutable access to the wrapped AST constraints — for the value-backed proxy.
    pub(crate) fn inner_mut(&mut self) -> &mut GraphIrDativeBondConstraintsAst {
        &mut self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge). Test-only —
    /// in-crate construction wraps `DativeBondConstraintsAst(..)` directly.
    pub(crate) fn from_inner(constraints: GraphIrDativeBondConstraintsAst) -> Self {
        DativeBondConstraintsAst(constraints)
    }
}

impl_py_lattice!(
    DativeBondConstraintsAst,
    GraphIrDativeBondConstraintsAst,
    |value: &DativeBondConstraintsAst,
     _py: Python<'_>|
     -> PyResult<GraphIrDativeBondConstraintsAst> { Ok(value.inner().clone()) },
    |_py: Python<'_>,
     value: GraphIrDativeBondConstraintsAst|
     -> PyResult<DativeBondConstraintsAst> { Ok(DativeBondConstraintsAst(value)) }
);

/// Build the per-constraint iterator handle from a borrowed container.
pub(crate) fn dative_bond_constraints_iter(
    py: Python<'_>,
    constraints: &GraphIrDativeBondConstraintsAst,
) -> PyResult<DativeBondConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| into_py_variant(py, DativeBondConstraintAst::from_rust(py, constraint)?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(DativeBondConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
pub(crate) fn dative_bond_constraint_keys(
    py: Python<'_>,
    constraints: &GraphIrDativeBondConstraintsAst,
) -> PyResult<DativeBondConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(
                py,
                DativeBondConstraintKey::from_rust(py, &constraint.key())?,
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(DativeBondConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
pub(crate) fn dative_bond_constraint_items(
    py: Python<'_>,
    constraints: &GraphIrDativeBondConstraintsAst,
) -> PyResult<DativeBondConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(
                    py,
                    DativeBondConstraintKey::from_rust(py, &constraint.key())?,
                )?,
                into_py_variant(py, DativeBondConstraintAst::from_rust(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(DativeBondConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// Python values. Ring memberships key by scope: `ring_count` for the
/// all-rings scope, `ring_size_count_<n>` for a specific ring size.
pub(crate) fn dative_bond_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &GraphIrDativeBondConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            GraphIrDativeBondConstraintAst::Aromatic(b) => {
                dict.set_item("aromatic", BooleanAst::from_rust(b))?
            }
            GraphIrDativeBondConstraintAst::RingMembership(m) => {
                let key = match m.scope {
                    GraphIrRingScope::All => "ring_count".to_string(),
                    GraphIrRingScope::Size(size) => format!("ring_size_count_{size}"),
                };
                dict.set_item(key, ValueAst::from_rust(py, &m.count)?)?
            }
        }
    }
    Ok(dict)
}

/// What a `DativeBondConstraintsView` writes through to. Only the standalone
/// `DativeBondAst` backing or a dative bond within a molecule (by index).
pub(crate) enum DativeBondConstraintsBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: GraphIrDativeBondId,
    },
    DativeBond(Py<DativeBondAst>),
}

/// A live handle onto one dative bond's constraints, backed by either a
/// molecule-bond or a standalone `DativeBondAst`. Reads borrow the bond's
/// constraints and read only the item they need (no whole-container clone);
/// mutators write through to the bond in place, without a clone-and-writeback.
#[pyclass]
pub struct DativeBondConstraintsView {
    pub(crate) backing: DativeBondConstraintsBacking,
}

impl DativeBondConstraintsView {
    /// Borrow the backing bond's constraints and read one item through `f` — no clone.
    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&GraphIrDativeBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            DativeBondConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .dative_bonds()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("dative bond id out of range"))?;
                f(&view.ast.constraints)
            }
            DativeBondConstraintsBacking::DativeBond(bond) => {
                let bond = bond.bind(py).borrow();
                f(&bond.inner().constraints)
            }
        }
    }

    /// Mutate the backing bond's constraints in place through `f`.
    pub(crate) fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut GraphIrDativeBondConstraintsAst) -> R,
    ) -> R {
        match &self.backing {
            DativeBondConstraintsBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .dative_bond_mut(*id)
                .ast
                .constraints),
            DativeBondConstraintsBacking::DativeBond(bond) => {
                f(&mut bond.borrow_mut(py).inner_mut().constraints)
            }
        }
    }

    /// Set one constraint on the backing bond in place (last-wins per key).
    pub(crate) fn set_ast(&self, py: Python<'_>, constraint: GraphIrDativeBondConstraintAst) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing bond in place, returning the removed entry.
    pub(crate) fn remove_ast(
        &self,
        py: Python<'_>,
        key: GraphIrDativeBondConstraintKey,
    ) -> Option<GraphIrDativeBondConstraintAst> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl DativeBondConstraintsView {
    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("DativeBondConstraintsView({count} entries)"))
    }

    /// Insert `c` on the bond in place, replacing any existing entry of the same
    /// key (last-wins).
    pub(crate) fn set(&self, py: Python<'_>, c: Py<DativeBondConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_rust(py));
    }

    /// Remove the entry with the given key from the bond in place, returning it if
    /// present (dict `pop`).
    pub(crate) fn pop(
        &self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<Option<DativeBondConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_rust(py))
            .map(|c| DativeBondConstraintAst::from_rust(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(
        &self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<()> {
        if self
            .remove_ast(py, key.bind(py).borrow().to_rust(py))
            .is_some()
        {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    /// Overlay `other` onto the bond's constraints in place — another container, a live
    /// view, or an iterable of `DativeBondConstraintAst` (last-wins per key; undetermined
    /// entries remove). Resolves `other` to owned data *before* the write borrow, so a
    /// view aliasing the same bond is not a double-borrow panic.
    pub(crate) fn update(
        &self,
        py: Python<'_>,
        other: DativeBondConstraintsUpdate,
    ) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        self.with_mut(py, |cs| resolved.apply(cs));
        Ok(())
    }

    pub(crate) fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<DativeBondConstraintKeyIter> {
        self.read(py, |cs| dative_bond_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<DativeBondConstraintKeyIter> {
        self.read(py, |cs| dative_bond_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<DativeBondConstraintIter> {
        self.read(py, |cs| dative_bond_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<DativeBondConstraintItemsIter> {
        self.read(py, |cs| dative_bond_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
        &self,
        py: Python<'_>,
        key: Py<DativeBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_rust(py);
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| DativeBondConstraintAst::from_rust(py, constraint))
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
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<DativeBondConstraintAst> {
        let ast_key = key.bind(py).borrow().to_rust(py);
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| DativeBondConstraintAst::from_rust(py, constraint))
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
        key: Py<DativeBondConstraintKey>,
    ) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_rust(py);
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// The aromatic value; `Undetermined` when no `Aromatic` constraint is present
    /// (matching the non-optional Rust accessor).
    #[getter]
    pub(crate) fn aromatic(&self, py: Python<'_>) -> PyResult<BooleanAst> {
        self.read(py, |cs| Ok(BooleanAst::from_rust(&cs.aromatic())))
    }

    #[setter]
    pub(crate) fn set_aromatic(&self, py: Python<'_>, value: BooleanLike) {
        self.set_ast(
            py,
            GraphIrDativeBondConstraintAst::aromatic(value.to_rust(py)),
        );
    }

    /// The all-rings membership count, or `None`.
    #[getter]
    pub(crate) fn ring_count(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_count()
                .map(|v| ValueAst::from_rust(py, v))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_ring_count(&self, py: Python<'_>, value: NumLike) {
        self.set_ast(
            py,
            GraphIrDativeBondConstraintAst::ring_membership(
                GraphIrRingScope::All,
                value.to_rust(py),
            ),
        );
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    pub(crate) fn ring_size_count(&self, py: Python<'_>) -> DativeBondRingSizeCounts {
        let backing = match &self.backing {
            DativeBondConstraintsBacking::Molecule { owner, id } => {
                DativeBondRingSizeBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: *id,
                }
            }
            DativeBondConstraintsBacking::DativeBond(bond) => {
                DativeBondRingSizeBacking::DativeBond(bond.clone_ref(py))
            }
        };
        DativeBondRingSizeCounts { backing }
    }

    /// The present constraints as a dict keyed by snake_case name.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| dative_bond_constraints_asdict(py, cs))
    }
}

/// What a `DativeBondRingSizeCounts` proxy reads/writes through to: a dative bond
/// within a molecule, a standalone `DativeBondAst`, or a standalone
/// `DativeBondConstraintsAst` value.
pub(crate) enum DativeBondRingSizeBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: GraphIrDativeBondId,
    },
    DativeBond(Py<DativeBondAst>),
    Value(Py<DativeBondConstraintsAst>),
}

/// A subscriptable proxy over the sized-ring membership counts of a dative bond,
/// keyed by ring size: `proxy[size]` reads, `proxy[size] = count` sets, `del
/// proxy[size]` removes. Backs onto whichever container produced it (dual-backing,
/// like `DativeBondConstraintsView`).
#[pyclass]
pub struct DativeBondRingSizeCounts {
    pub(crate) backing: DativeBondRingSizeBacking,
}

impl DativeBondRingSizeCounts {
    /// Borrow the backing constraints and read through `f` — no clone.
    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&GraphIrDativeBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            DativeBondRingSizeBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .dative_bonds()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("dative bond id out of range"))?;
                f(&view.ast.constraints)
            }
            DativeBondRingSizeBacking::DativeBond(bond) => {
                f(&bond.bind(py).borrow().inner().constraints)
            }
            DativeBondRingSizeBacking::Value(value) => f(value.bind(py).borrow().inner()),
        }
    }

    /// Mutate the backing constraints in place through `f`.
    pub(crate) fn write(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut GraphIrDativeBondConstraintsAst),
    ) {
        match &self.backing {
            DativeBondRingSizeBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .dative_bond_mut(*id)
                .ast
                .constraints),
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
    pub(crate) fn __getitem__(&self, py: Python<'_>, size: u8) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_size_count(size)
                .map(|v| ValueAst::from_rust(py, v))
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
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<DativeBondRingSizeIter> {
        let sizes = self.read(py, |cs| Ok(ring_sizes(cs).collect::<Vec<u8>>()))?;
        Ok(DativeBondRingSizeIter {
            sizes: sizes.into_iter(),
        })
    }

    /// Set the membership count for rings of `size` in place.
    pub(crate) fn __setitem__(&self, py: Python<'_>, size: u8, count: NumLike) {
        let constraint = GraphIrDativeBondConstraintAst::ring_membership(
            GraphIrRingScope::Size(size),
            count.to_rust(py),
        );
        self.write(py, |cs| cs.set(constraint));
    }

    /// Remove the sized-ring membership for `size` in place.
    pub(crate) fn __delitem__(&self, py: Python<'_>, size: u8) {
        self.write(py, |cs| {
            cs.remove(GraphIrDativeBondConstraintKey::RingMembership(
                GraphIrRingScope::Size(size),
            ));
        });
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |cs| {
            let mut parts = Vec::new();
            for entry in cs.iter() {
                if let GraphIrDativeBondConstraintAst::RingMembership(m) = entry {
                    if let GraphIrRingScope::Size(size) = m.scope {
                        let count = into_py_variant(py, ValueAst::from_rust(py, &m.count)?)?;
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
pub(crate) fn ring_sizes(
    constraints: &GraphIrDativeBondConstraintsAst,
) -> impl Iterator<Item = u8> + '_ {
    constraints.iter().filter_map(|entry| match entry {
        GraphIrDativeBondConstraintAst::RingMembership(m) => match m.scope {
            GraphIrRingScope::Size(size) => Some(size),
            GraphIrRingScope::All => None,
        },
        _ => None,
    })
}

#[pyclass]
pub(crate) struct DativeBondRingSizeIter {
    sizes: IntoIter<u8>,
}

#[pymethods]
impl DativeBondRingSizeIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<u8> {
        self.sizes.next()
    }
}

#[pyclass]
pub(crate) struct DativeBondConstraintIter {
    entries: IntoIter<Py<DativeBondConstraintAst>>,
}

#[pymethods]
impl DativeBondConstraintIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<DativeBondConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
pub(crate) struct DativeBondConstraintKeyIter {
    keys: IntoIter<Py<DativeBondConstraintKey>>,
}

#[pymethods]
impl DativeBondConstraintKeyIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<DativeBondConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
pub(crate) struct DativeBondConstraintItemsIter {
    items: IntoIter<(Py<DativeBondConstraintKey>, Py<DativeBondConstraintAst>)>,
}

#[pymethods]
impl DativeBondConstraintItemsIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(
        &mut self,
    ) -> Option<(Py<DativeBondConstraintKey>, Py<DativeBondConstraintAst>)> {
        self.items.next()
    }
}
