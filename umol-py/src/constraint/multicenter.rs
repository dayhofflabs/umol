//! Multicenter-bond constraint values, containers, and live views.

use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use umol_ast::ast::{
    MulticenterBondConstraintAst as AstMulticenterBondConstraintAst,
    MulticenterBondConstraintKey as AstMulticenterBondConstraintKey,
    MulticenterBondConstraintsAst as AstMulticenterBondConstraintsAst,
    MulticenterBondId as AstMulticenterBondId,
};

use crate::convert::{hash_rust, into_py_variant, variant_repr};
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;
use crate::multicenter::MulticenterBondAst;
use crate::value::{ValueArg, ValueAst};

/// The key (identity) of a multicenter-bond constraint, for keyed lookup. The
/// single key `ElectronCount` is the bare discriminant (no sub-key).
#[pyclass]
pub enum MulticenterBondConstraintKey {
    ElectronCount(),
}

#[pymethods]
impl MulticenterBondConstraintKey {
    pub(crate) fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    pub(crate) fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            MulticenterBondConstraintKey::ElectronCount() => ("ElectronCount", 0),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "MulticenterBondConstraintKey",
            variant,
            arity,
        )
    }
}

impl MulticenterBondConstraintKey {
    pub(crate) fn from_rust(ast: &AstMulticenterBondConstraintKey) -> Self {
        match ast {
            AstMulticenterBondConstraintKey::ElectronCount => Self::ElectronCount(),
        }
    }

    pub(crate) fn to_rust(&self) -> AstMulticenterBondConstraintKey {
        match self {
            Self::ElectronCount() => AstMulticenterBondConstraintKey::ElectronCount,
        }
    }
}

/// A multicenter-bond-scope constraint: the asserted total electron count of the
/// bond (cross-checked against `sum(MulticenterBondAst::electrons)`).
#[pyclass]
pub enum MulticenterBondConstraintAst {
    ElectronCount(Py<ValueAst>),
}

#[pymethods]
impl MulticenterBondConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    pub(crate) fn key(&self, py: Python<'_>) -> MulticenterBondConstraintKey {
        MulticenterBondConstraintKey::from_rust(&self.to_rust(py).key())
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            MulticenterBondConstraintAst::ElectronCount(_) => "ElectronCount",
        };
        variant_repr(
            slf.bind(py).as_any(),
            "MulticenterBondConstraintAst",
            variant,
            1,
        )
    }
}

impl_py_lattice!(
    MulticenterBondConstraintAst,
    AstMulticenterBondConstraintAst,
    |value: &MulticenterBondConstraintAst,
     py: Python<'_>|
     -> PyResult<AstMulticenterBondConstraintAst> { Ok(value.to_rust(py)) },
    |py: Python<'_>,
     value: AstMulticenterBondConstraintAst|
     -> PyResult<MulticenterBondConstraintAst> {
        MulticenterBondConstraintAst::from_rust(py, &value)
    }
);

impl MulticenterBondConstraintAst {
    pub(crate) fn from_rust(
        py: Python<'_>,
        ast: &AstMulticenterBondConstraintAst,
    ) -> PyResult<Self> {
        Ok(match ast {
            AstMulticenterBondConstraintAst::ElectronCount(v) => {
                Self::ElectronCount(into_py_variant(py, ValueAst::from_rust(py, v)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstMulticenterBondConstraintAst {
        match self {
            Self::ElectronCount(v) => {
                AstMulticenterBondConstraintAst::ElectronCount(v.bind(py).borrow().to_rust(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container (value or live view) or an
/// iterable of `MulticenterBondConstraintAst` (each `set`, last-wins).
#[derive(FromPyObject)]
pub(crate) enum MulticenterBondConstraintsUpdate {
    Container(Py<MulticenterBondConstraintsAst>),
    View(Py<MulticenterBondConstraintsView>),
    Entries(Vec<Py<MulticenterBondConstraintAst>>),
}

impl MulticenterBondConstraintsUpdate {
    /// Read every Python object into owned data — no write target is touched. Callers
    /// resolve *before* taking the write borrow so a view (or container) that aliases the
    /// same bond is read while nothing is borrowed (otherwise
    /// `bond.constraints.update(bond.constraints)` self-aliases into a double-borrow panic).
    pub(crate) fn resolve(
        &self,
        py: Python<'_>,
    ) -> PyResult<ResolvedMulticenterBondConstraintsUpdate> {
        Ok(match self {
            MulticenterBondConstraintsUpdate::Container(c) => {
                ResolvedMulticenterBondConstraintsUpdate::Overlay(
                    c.bind(py).borrow().inner().clone(),
                )
            }
            MulticenterBondConstraintsUpdate::View(v) => {
                ResolvedMulticenterBondConstraintsUpdate::Overlay(
                    v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?,
                )
            }
            MulticenterBondConstraintsUpdate::Entries(entries) => {
                ResolvedMulticenterBondConstraintsUpdate::Entries(
                    entries
                        .iter()
                        .map(|entry| entry.bind(py).borrow().to_rust(py))
                        .collect(),
                )
            }
        })
    }
}

/// A `MulticenterBondConstraintsUpdate` with all Python-object reads already done, so it
/// can be applied under a write borrow without re-entering Python.
pub(crate) enum ResolvedMulticenterBondConstraintsUpdate {
    /// A whole container (from another container or a live view): overlaid via `update`
    /// (last-wins per key; undetermined entries remove).
    Overlay(AstMulticenterBondConstraintsAst),
    /// Loose entries: `set` each (last-wins; undetermined entries stored, not removed).
    Entries(Vec<AstMulticenterBondConstraintAst>),
}

impl ResolvedMulticenterBondConstraintsUpdate {
    /// Overlay onto `target` in place. No Python reads.
    pub(crate) fn apply(self, target: &mut AstMulticenterBondConstraintsAst) {
        match self {
            ResolvedMulticenterBondConstraintsUpdate::Overlay(overlay) => target.update(&overlay),
            ResolvedMulticenterBondConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry);
                }
            }
        }
    }
}

/// A whole-container argument that snapshots either a value container or a live view
/// — for the multicenter bond `constraints` setter, which accepts either.
#[derive(FromPyObject)]
pub(crate) enum MulticenterBondConstraintsArg {
    Container(Py<MulticenterBondConstraintsAst>),
    View(Py<MulticenterBondConstraintsView>),
}

impl MulticenterBondConstraintsArg {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<AstMulticenterBondConstraintsAst> {
        match self {
            MulticenterBondConstraintsArg::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
            MulticenterBondConstraintsArg::View(v) => {
                v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))
            }
        }
    }
}

/// The multicenter-bond-scope constraints on a multicenter bond, in kind-sorted order.
/// Mutable, hence value-equal but unhashable (matching `MulticenterBondAst`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct MulticenterBondConstraintsAst(AstMulticenterBondConstraintsAst);

#[pymethods]
impl MulticenterBondConstraintsAst {
    /// Build from a sequence of constraints (a later entry of the same key replaces
    /// an earlier one, last-wins).
    #[new]
    pub(crate) fn new(py: Python<'_>, entries: Vec<Py<MulticenterBondConstraintAst>>) -> Self {
        let mut constraints = AstMulticenterBondConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py)),
        );
        MulticenterBondConstraintsAst(constraints)
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let value = into_py_variant(py, MulticenterBondConstraintAst::from_rust(py, entry)?)?;
            parts.push(value.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!(
            "MulticenterBondConstraintsAst([{}])",
            parts.join(", ")
        ))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    pub(crate) fn set(&mut self, py: Python<'_>, c: Py<MulticenterBondConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_rust(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    pub(crate) fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<Option<MulticenterBondConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_rust())
            .map(|c| MulticenterBondConstraintAst::from_rust(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container, a live view, or an
    /// iterable of `MulticenterBondConstraintAst` (last-wins per key; undetermined entries
    /// remove). Takes `slf` by handle so `other` is fully read *before* the write borrow —
    /// `cs.update(cs)` on the same container is then a no-op, not a double-borrow panic.
    pub(crate) fn update(
        slf: Py<Self>,
        py: Python<'_>,
        other: MulticenterBondConstraintsUpdate,
    ) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        resolved.apply(&mut slf.borrow_mut(py).0);
        Ok(())
    }

    pub(crate) fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintKeyIter> {
        multicenter_bond_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintKeyIter> {
        multicenter_bond_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintIter> {
        multicenter_bond_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintItemsIter> {
        multicenter_bond_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_rust()) {
            Some(constraint) => Ok(into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(py, constraint)?,
            )?
            .into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    pub(crate) fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<MulticenterBondConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_rust()) {
            Some(constraint) => MulticenterBondConstraintAst::from_rust(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_rust()).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    pub(crate) fn __contains__(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> bool {
        self.0.contains(key.bind(py).borrow().to_rust())
    }

    /// The asserted total electron count; `Undetermined` when no `ElectronCount`
    /// constraint is present (matching the non-optional Rust accessor).
    #[getter]
    pub(crate) fn electron_count(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_rust(py, &self.0.electron_count())
    }

    #[setter]
    pub(crate) fn set_electron_count(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstMulticenterBondConstraintAst::electron_count(
            value.to_rust(py),
        ));
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// Python values.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        multicenter_bond_constraints_asdict(py, &self.0)
    }
}

impl MulticenterBondConstraintsAst {
    /// The wrapped AST constraints — read access for multicenter bond construction.
    pub(crate) fn inner(&self) -> &AstMulticenterBondConstraintsAst {
        &self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge). Test-only —
    /// in-crate construction wraps `MulticenterBondConstraintsAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(constraints: AstMulticenterBondConstraintsAst) -> Self {
        MulticenterBondConstraintsAst(constraints)
    }
}

impl_py_lattice!(
    MulticenterBondConstraintsAst,
    AstMulticenterBondConstraintsAst,
    |value: &MulticenterBondConstraintsAst,
     _py: Python<'_>|
     -> PyResult<AstMulticenterBondConstraintsAst> { Ok(value.inner().clone()) },
    |_py: Python<'_>,
     value: AstMulticenterBondConstraintsAst|
     -> PyResult<MulticenterBondConstraintsAst> { Ok(MulticenterBondConstraintsAst(value)) }
);

/// Build the per-constraint iterator handle from a borrowed container.
pub(crate) fn multicenter_bond_constraints_iter(
    py: Python<'_>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<MulticenterBondConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(py, MulticenterBondConstraintAst::from_rust(py, constraint)?)
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(MulticenterBondConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
pub(crate) fn multicenter_bond_constraint_keys(
    py: Python<'_>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<MulticenterBondConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(
                py,
                MulticenterBondConstraintKey::from_rust(&constraint.key()),
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(MulticenterBondConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
pub(crate) fn multicenter_bond_constraint_items(
    py: Python<'_>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<MulticenterBondConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(
                    py,
                    MulticenterBondConstraintKey::from_rust(&constraint.key()),
                )?,
                into_py_variant(py, MulticenterBondConstraintAst::from_rust(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(MulticenterBondConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// Python values.
pub(crate) fn multicenter_bond_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            AstMulticenterBondConstraintAst::ElectronCount(v) => {
                dict.set_item("electron_count", ValueAst::from_rust(py, v)?)?
            }
        }
    }
    Ok(dict)
}

/// What a `MulticenterBondConstraintsView` writes through to: a multicenter bond
/// within a molecule (by index) or a standalone `MulticenterBondAst`.
pub(crate) enum MulticenterBondConstraintsBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: AstMulticenterBondId,
    },
    MulticenterBond(Py<MulticenterBondAst>),
}

/// A live handle onto one multicenter bond's constraints, backed by either a
/// molecule-bond or a standalone `MulticenterBondAst`. Reads borrow the constraints
/// and read only the item they need (no whole-container clone); mutators write through
/// to the bond in place, without a clone-and-writeback.
#[pyclass]
pub struct MulticenterBondConstraintsView {
    pub(crate) backing: MulticenterBondConstraintsBacking,
}

impl MulticenterBondConstraintsView {
    /// Borrow the backing bond's constraints and read one item through `f` — no clone.
    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstMulticenterBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            MulticenterBondConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .multicenter_bonds()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("multicenter bond id out of range"))?;
                f(&view.ast.constraints)
            }
            MulticenterBondConstraintsBacking::MulticenterBond(bond) => {
                let bond = bond.bind(py).borrow();
                f(&bond.inner().constraints)
            }
        }
    }

    /// Mutate the backing bond's constraints in place through `f`.
    pub(crate) fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut AstMulticenterBondConstraintsAst) -> R,
    ) -> R {
        match &self.backing {
            MulticenterBondConstraintsBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .multicenter_bond_mut(*id)
                .ast
                .constraints),
            MulticenterBondConstraintsBacking::MulticenterBond(bond) => {
                f(&mut bond.borrow_mut(py).inner_mut().constraints)
            }
        }
    }

    /// Set one constraint on the backing bond in place (last-wins per key).
    pub(crate) fn set_ast(&self, py: Python<'_>, constraint: AstMulticenterBondConstraintAst) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing bond in place, returning the removed entry.
    pub(crate) fn remove_ast(
        &self,
        py: Python<'_>,
        key: AstMulticenterBondConstraintKey,
    ) -> Option<AstMulticenterBondConstraintAst> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl MulticenterBondConstraintsView {
    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("MulticenterBondConstraintsView({count} entries)"))
    }

    /// Insert `c` on the bond in place, replacing any existing entry of the same key
    /// (last-wins).
    pub(crate) fn set(&self, py: Python<'_>, c: Py<MulticenterBondConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_rust(py));
    }

    /// Remove the entry with the given key from the bond in place, returning it if
    /// present (dict `pop`).
    pub(crate) fn pop(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<Option<MulticenterBondConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_rust())
            .map(|c| MulticenterBondConstraintAst::from_rust(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<()> {
        if self
            .remove_ast(py, key.bind(py).borrow().to_rust())
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
    /// view, or an iterable of `MulticenterBondConstraintAst` (last-wins per key;
    /// undetermined entries remove). Resolves `other` to owned data *before* the write
    /// borrow, so a view aliasing the same bond is not a double-borrow panic.
    pub(crate) fn update(
        &self,
        py: Python<'_>,
        other: MulticenterBondConstraintsUpdate,
    ) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        self.with_mut(py, |cs| resolved.apply(cs));
        Ok(())
    }

    pub(crate) fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintKeyIter> {
        self.read(py, |cs| multicenter_bond_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintKeyIter> {
        self.read(py, |cs| multicenter_bond_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintIter> {
        self.read(py, |cs| multicenter_bond_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintItemsIter> {
        self.read(py, |cs| multicenter_bond_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_rust();
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| MulticenterBondConstraintAst::from_rust(py, constraint))
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
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<MulticenterBondConstraintAst> {
        let ast_key = key.bind(py).borrow().to_rust();
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| MulticenterBondConstraintAst::from_rust(py, constraint))
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
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_rust();
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// The asserted total electron count; `Undetermined` when no `ElectronCount`
    /// constraint is present (matching the non-optional Rust accessor).
    #[getter]
    pub(crate) fn electron_count(&self, py: Python<'_>) -> PyResult<ValueAst> {
        self.read(py, |cs| ValueAst::from_rust(py, &cs.electron_count()))
    }

    #[setter]
    pub(crate) fn set_electron_count(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(
            py,
            AstMulticenterBondConstraintAst::electron_count(value.to_rust(py)),
        );
    }

    /// The present constraints as a dict keyed by snake_case name.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| multicenter_bond_constraints_asdict(py, cs))
    }
}

#[pyclass]
pub(crate) struct MulticenterBondConstraintIter {
    entries: IntoIter<Py<MulticenterBondConstraintAst>>,
}

#[pymethods]
impl MulticenterBondConstraintIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<MulticenterBondConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
pub(crate) struct MulticenterBondConstraintKeyIter {
    keys: IntoIter<Py<MulticenterBondConstraintKey>>,
}

#[pymethods]
impl MulticenterBondConstraintKeyIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<MulticenterBondConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
pub(crate) struct MulticenterBondConstraintItemsIter {
    items: IntoIter<(
        Py<MulticenterBondConstraintKey>,
        Py<MulticenterBondConstraintAst>,
    )>,
}

#[pymethods]
impl MulticenterBondConstraintItemsIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(
        &mut self,
    ) -> Option<(
        Py<MulticenterBondConstraintKey>,
        Py<MulticenterBondConstraintAst>,
    )> {
        self.items.next()
    }
}
