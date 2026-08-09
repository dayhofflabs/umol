//! Aromatic-system constraint values, containers, and live views.

use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use umol_graph_ir::ir::{
    AromaticSystemConstraintForm as GraphIrAromaticSystemConstraintForm,
    AromaticSystemConstraintKey as GraphIrAromaticSystemConstraintKey,
    AromaticSystemConstraintsForm as GraphIrAromaticSystemConstraintsForm,
    AromaticSystemId as GraphIrAromaticSystemId,
};

use crate::aromatic::AromaticSystemForm;
use crate::convert::{hash_rust, into_py_variant, variant_repr};
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;
use crate::value::{NumForm, NumLike};

/// The key (identity) of an aromatic-system constraint, for keyed lookup. The
/// single key `ElectronCount` is the bare discriminant (no sub-key).
#[pyclass]
pub enum AromaticSystemConstraintKey {
    ElectronCount(),
}

#[pymethods]
impl AromaticSystemConstraintKey {
    pub(crate) fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    pub(crate) fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            AromaticSystemConstraintKey::ElectronCount() => ("ElectronCount", 0),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "AromaticSystemConstraintKey",
            variant,
            arity,
        )
    }
}

impl AromaticSystemConstraintKey {
    pub(crate) fn from_rust(ast: &GraphIrAromaticSystemConstraintKey) -> Self {
        match ast {
            GraphIrAromaticSystemConstraintKey::ElectronCount => Self::ElectronCount(),
        }
    }

    pub(crate) fn to_rust(&self) -> GraphIrAromaticSystemConstraintKey {
        match self {
            Self::ElectronCount() => GraphIrAromaticSystemConstraintKey::ElectronCount,
        }
    }
}

/// An aromatic-system-scope constraint: the asserted total π-electron count of the
/// system (cross-checked against `sum(AromaticSystemForm::electrons)`).
#[pyclass]
pub enum AromaticSystemConstraintAst {
    ElectronCount(Py<NumForm>),
}

#[pymethods]
impl AromaticSystemConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    pub(crate) fn key(&self, py: Python<'_>) -> AromaticSystemConstraintKey {
        AromaticSystemConstraintKey::from_rust(&self.to_rust(py).key())
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            AromaticSystemConstraintAst::ElectronCount(_) => "ElectronCount",
        };
        variant_repr(
            slf.bind(py).as_any(),
            "AromaticSystemConstraintAst",
            variant,
            1,
        )
    }
}

impl_py_lattice!(
    AromaticSystemConstraintAst,
    GraphIrAromaticSystemConstraintForm,
    |value: &AromaticSystemConstraintAst,
     py: Python<'_>|
     -> PyResult<GraphIrAromaticSystemConstraintForm> { Ok(value.to_rust(py)) },
    |py: Python<'_>,
     value: GraphIrAromaticSystemConstraintForm|
     -> PyResult<AromaticSystemConstraintAst> {
        AromaticSystemConstraintAst::from_rust(py, &value)
    }
);

impl AromaticSystemConstraintAst {
    pub(crate) fn from_rust(
        py: Python<'_>,
        ast: &GraphIrAromaticSystemConstraintForm,
    ) -> PyResult<Self> {
        Ok(match ast {
            GraphIrAromaticSystemConstraintForm::ElectronCount(v) => {
                Self::ElectronCount(into_py_variant(py, NumForm::from_rust(py, v)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrAromaticSystemConstraintForm {
        match self {
            Self::ElectronCount(v) => {
                GraphIrAromaticSystemConstraintForm::ElectronCount(v.bind(py).borrow().to_rust(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container (value or live view) or an
/// iterable of `AromaticSystemConstraintAst` (each `set`, last-wins).
#[derive(FromPyObject)]
pub(crate) enum AromaticSystemConstraintsUpdate {
    Container(Py<AromaticSystemConstraintsAst>),
    View(Py<AromaticSystemConstraintsView>),
    Entries(Vec<Py<AromaticSystemConstraintAst>>),
}

impl AromaticSystemConstraintsUpdate {
    /// Read every Python object into owned data — no write target is touched. Callers
    /// resolve *before* taking the write borrow so a view (or container) that aliases the
    /// same system is read while nothing is borrowed (otherwise
    /// `sys.constraints.update(sys.constraints)` self-aliases into a double-borrow panic).
    pub(crate) fn resolve(
        &self,
        py: Python<'_>,
    ) -> PyResult<ResolvedAromaticSystemConstraintsUpdate> {
        Ok(match self {
            AromaticSystemConstraintsUpdate::Container(c) => {
                ResolvedAromaticSystemConstraintsUpdate::Overlay(
                    c.bind(py).borrow().inner().clone(),
                )
            }
            AromaticSystemConstraintsUpdate::View(v) => {
                ResolvedAromaticSystemConstraintsUpdate::Overlay(
                    v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?,
                )
            }
            AromaticSystemConstraintsUpdate::Entries(entries) => {
                ResolvedAromaticSystemConstraintsUpdate::Entries(
                    entries
                        .iter()
                        .map(|entry| entry.bind(py).borrow().to_rust(py))
                        .collect(),
                )
            }
        })
    }
}

/// An `AromaticSystemConstraintsUpdate` with all Python-object reads already done, so it
/// can be applied under a write borrow without re-entering Python.
pub(crate) enum ResolvedAromaticSystemConstraintsUpdate {
    /// A whole container (from another container or a live view): overlaid via `update`
    /// (last-wins per key; undetermined entries remove).
    Overlay(GraphIrAromaticSystemConstraintsForm),
    /// Loose entries: `set` each (last-wins; undetermined entries stored, not removed).
    Entries(Vec<GraphIrAromaticSystemConstraintForm>),
}

impl ResolvedAromaticSystemConstraintsUpdate {
    /// Overlay onto `target` in place. No Python reads.
    pub(crate) fn apply(self, target: &mut GraphIrAromaticSystemConstraintsForm) {
        match self {
            ResolvedAromaticSystemConstraintsUpdate::Overlay(overlay) => target.update(&overlay),
            ResolvedAromaticSystemConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry);
                }
            }
        }
    }
}

/// A whole-container argument that snapshots either a value container or a live view
/// — for the aromatic system `constraints` setter, which accepts either.
#[derive(FromPyObject)]
pub(crate) enum AromaticSystemConstraintsLike {
    Container(Py<AromaticSystemConstraintsAst>),
    View(Py<AromaticSystemConstraintsView>),
}

impl AromaticSystemConstraintsLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<GraphIrAromaticSystemConstraintsForm> {
        match self {
            AromaticSystemConstraintsLike::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
            AromaticSystemConstraintsLike::View(v) => {
                v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))
            }
        }
    }
}

/// The aromatic-system-scope constraints on an aromatic system, in kind-sorted order.
/// Mutable, hence value-equal but unhashable (matching `AromaticSystemForm`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct AromaticSystemConstraintsAst(GraphIrAromaticSystemConstraintsForm);

#[pymethods]
impl AromaticSystemConstraintsAst {
    /// Build from a sequence of constraints (a later entry of the same key replaces
    /// an earlier one, last-wins).
    #[new]
    pub(crate) fn new(py: Python<'_>, entries: Vec<Py<AromaticSystemConstraintAst>>) -> Self {
        let mut constraints = GraphIrAromaticSystemConstraintsForm::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py)),
        );
        AromaticSystemConstraintsAst(constraints)
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let value = into_py_variant(py, AromaticSystemConstraintAst::from_rust(py, entry)?)?;
            parts.push(value.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!(
            "AromaticSystemConstraintsAst([{}])",
            parts.join(", ")
        ))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    pub(crate) fn set(&mut self, py: Python<'_>, c: Py<AromaticSystemConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_rust(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    pub(crate) fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<Option<AromaticSystemConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_rust())
            .map(|c| AromaticSystemConstraintAst::from_rust(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container, a live view, or an
    /// iterable of `AromaticSystemConstraintAst` (last-wins per key; undetermined entries
    /// remove). Takes `slf` by handle so `other` is fully read *before* the write borrow —
    /// `cs.update(cs)` on the same container is then a no-op, not a double-borrow panic.
    pub(crate) fn update(
        slf: Py<Self>,
        py: Python<'_>,
        other: AromaticSystemConstraintsUpdate,
    ) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        resolved.apply(&mut slf.borrow_mut(py).0);
        Ok(())
    }

    pub(crate) fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintKeyIter> {
        aromatic_system_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintKeyIter> {
        aromatic_system_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintIter> {
        aromatic_system_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintItemsIter> {
        aromatic_system_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_rust()) {
            Some(constraint) => Ok(into_py_variant(
                py,
                AromaticSystemConstraintAst::from_rust(py, constraint)?,
            )?
            .into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    pub(crate) fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<AromaticSystemConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_rust()) {
            Some(constraint) => AromaticSystemConstraintAst::from_rust(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
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
        key: Py<AromaticSystemConstraintKey>,
    ) -> bool {
        self.0.contains(key.bind(py).borrow().to_rust())
    }

    /// The asserted total π-electron count; `Undetermined` when no `ElectronCount`
    /// constraint is present (matching the non-optional Rust accessor).
    #[getter]
    pub(crate) fn electron_count(&self, py: Python<'_>) -> PyResult<NumForm> {
        NumForm::from_rust(py, &self.0.electron_count())
    }

    #[setter]
    pub(crate) fn set_electron_count(&mut self, py: Python<'_>, value: NumLike) {
        self.0
            .set(GraphIrAromaticSystemConstraintForm::electron_count(
                value.to_rust(py),
            ));
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// Python values.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        aromatic_system_constraints_asdict(py, &self.0)
    }
}

impl AromaticSystemConstraintsAst {
    /// The wrapped AST constraints — read access for aromatic system construction.
    pub(crate) fn inner(&self) -> &GraphIrAromaticSystemConstraintsForm {
        &self.0
    }

    /// Wrap owned AST constraints.
    pub(crate) fn from_inner(constraints: GraphIrAromaticSystemConstraintsForm) -> Self {
        AromaticSystemConstraintsAst(constraints)
    }
}

impl_py_lattice!(
    AromaticSystemConstraintsAst,
    GraphIrAromaticSystemConstraintsForm,
    |value: &AromaticSystemConstraintsAst,
     _py: Python<'_>|
     -> PyResult<GraphIrAromaticSystemConstraintsForm> { Ok(value.inner().clone()) },
    |_py: Python<'_>,
     value: GraphIrAromaticSystemConstraintsForm|
     -> PyResult<AromaticSystemConstraintsAst> { Ok(AromaticSystemConstraintsAst(value)) }
);

/// Build the per-constraint iterator handle from a borrowed container.
pub(crate) fn aromatic_system_constraints_iter(
    py: Python<'_>,
    constraints: &GraphIrAromaticSystemConstraintsForm,
) -> PyResult<AromaticSystemConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(py, AromaticSystemConstraintAst::from_rust(py, constraint)?)
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AromaticSystemConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
pub(crate) fn aromatic_system_constraint_keys(
    py: Python<'_>,
    constraints: &GraphIrAromaticSystemConstraintsForm,
) -> PyResult<AromaticSystemConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(
                py,
                AromaticSystemConstraintKey::from_rust(&constraint.key()),
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AromaticSystemConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
pub(crate) fn aromatic_system_constraint_items(
    py: Python<'_>,
    constraints: &GraphIrAromaticSystemConstraintsForm,
) -> PyResult<AromaticSystemConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(
                    py,
                    AromaticSystemConstraintKey::from_rust(&constraint.key()),
                )?,
                into_py_variant(py, AromaticSystemConstraintAst::from_rust(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AromaticSystemConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// Python values.
pub(crate) fn aromatic_system_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &GraphIrAromaticSystemConstraintsForm,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            GraphIrAromaticSystemConstraintForm::ElectronCount(v) => {
                dict.set_item("electron_count", NumForm::from_rust(py, v)?)?
            }
        }
    }
    Ok(dict)
}

/// What an `AromaticSystemConstraintsView` writes through to: an aromatic system
/// within a molecule (by index) or a standalone `AromaticSystemForm`.
pub(crate) enum AromaticSystemConstraintsBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: GraphIrAromaticSystemId,
    },
    AromaticSystem(Py<AromaticSystemForm>),
}

/// A live handle onto one aromatic system's constraints, backed by either a
/// molecule-system or a standalone `AromaticSystemForm`. Reads borrow the constraints
/// and read only the item they need (no whole-container clone); mutators write through
/// to the system in place, without a clone-and-writeback.
#[pyclass]
pub struct AromaticSystemConstraintsView {
    pub(crate) backing: AromaticSystemConstraintsBacking,
}

impl AromaticSystemConstraintsView {
    /// Borrow the backing system's constraints and read one item through `f` — no clone.
    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&GraphIrAromaticSystemConstraintsForm) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            AromaticSystemConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .aromatic_systems()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("aromatic system id out of range"))?;
                f(&view.ast.constraints)
            }
            AromaticSystemConstraintsBacking::AromaticSystem(system) => {
                let system = system.bind(py).borrow();
                f(&system.inner().constraints)
            }
        }
    }

    /// Mutate the backing system's constraints in place through `f`.
    pub(crate) fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut GraphIrAromaticSystemConstraintsForm) -> R,
    ) -> R {
        match &self.backing {
            AromaticSystemConstraintsBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .aromatic_system_mut(*id)
                .ast
                .constraints),
            AromaticSystemConstraintsBacking::AromaticSystem(system) => {
                f(&mut system.borrow_mut(py).inner_mut().constraints)
            }
        }
    }

    /// Set one constraint on the backing system in place (last-wins per key).
    pub(crate) fn set_ast(&self, py: Python<'_>, constraint: GraphIrAromaticSystemConstraintForm) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing system in place, returning the removed entry.
    pub(crate) fn remove_ast(
        &self,
        py: Python<'_>,
        key: GraphIrAromaticSystemConstraintKey,
    ) -> Option<GraphIrAromaticSystemConstraintForm> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl AromaticSystemConstraintsView {
    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("AromaticSystemConstraintsView({count} entries)"))
    }

    /// Insert `c` on the system in place, replacing any existing entry of the same key
    /// (last-wins).
    pub(crate) fn set(&self, py: Python<'_>, c: Py<AromaticSystemConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_rust(py));
    }

    /// Remove the entry with the given key from the system in place, returning it if
    /// present (dict `pop`).
    pub(crate) fn pop(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<Option<AromaticSystemConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_rust())
            .map(|c| AromaticSystemConstraintAst::from_rust(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
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

    /// Overlay `other` onto the system's constraints in place — another container, a live
    /// view, or an iterable of `AromaticSystemConstraintAst` (last-wins per key;
    /// undetermined entries remove). Resolves `other` to owned data *before* the write
    /// borrow, so a view aliasing the same system is not a double-borrow panic.
    pub(crate) fn update(
        &self,
        py: Python<'_>,
        other: AromaticSystemConstraintsUpdate,
    ) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        self.with_mut(py, |cs| resolved.apply(cs));
        Ok(())
    }

    pub(crate) fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintKeyIter> {
        self.read(py, |cs| aromatic_system_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintKeyIter> {
        self.read(py, |cs| aromatic_system_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintIter> {
        self.read(py, |cs| aromatic_system_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintItemsIter> {
        self.read(py, |cs| aromatic_system_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_rust();
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| AromaticSystemConstraintAst::from_rust(py, constraint))
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
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<AromaticSystemConstraintAst> {
        let ast_key = key.bind(py).borrow().to_rust();
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| AromaticSystemConstraintAst::from_rust(py, constraint))
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
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_rust();
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// The asserted total π-electron count; `Undetermined` when no `ElectronCount`
    /// constraint is present (matching the non-optional Rust accessor).
    #[getter]
    pub(crate) fn electron_count(&self, py: Python<'_>) -> PyResult<NumForm> {
        self.read(py, |cs| NumForm::from_rust(py, &cs.electron_count()))
    }

    #[setter]
    pub(crate) fn set_electron_count(&self, py: Python<'_>, value: NumLike) {
        self.set_ast(
            py,
            GraphIrAromaticSystemConstraintForm::electron_count(value.to_rust(py)),
        );
    }

    /// The present constraints as a dict keyed by snake_case name.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| aromatic_system_constraints_asdict(py, cs))
    }
}

#[pyclass]
pub(crate) struct AromaticSystemConstraintIter {
    entries: IntoIter<Py<AromaticSystemConstraintAst>>,
}

#[pymethods]
impl AromaticSystemConstraintIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<AromaticSystemConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
pub(crate) struct AromaticSystemConstraintKeyIter {
    keys: IntoIter<Py<AromaticSystemConstraintKey>>,
}

#[pymethods]
impl AromaticSystemConstraintKeyIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<AromaticSystemConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
pub(crate) struct AromaticSystemConstraintItemsIter {
    items: IntoIter<(
        Py<AromaticSystemConstraintKey>,
        Py<AromaticSystemConstraintAst>,
    )>,
}

#[pymethods]
impl AromaticSystemConstraintItemsIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(
        &mut self,
    ) -> Option<(
        Py<AromaticSystemConstraintKey>,
        Py<AromaticSystemConstraintAst>,
    )> {
        self.items.next()
    }
}
