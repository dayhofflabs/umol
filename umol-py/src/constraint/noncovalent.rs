//! Noncovalent-bond constraint values, containers, and live views.

use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use umol_graph_ir::ir::{
    NoncovalentBondConstraintForm as GraphIrNoncovalentBondConstraintForm,
    NoncovalentBondConstraintKey as GraphIrNoncovalentBondConstraintKey,
    NoncovalentBondConstraintsForm as GraphIrNoncovalentBondConstraintsForm,
    NoncovalentBondId as GraphIrNoncovalentBondId,
};

use crate::boolean::{BooleanForm, BooleanLike};
use crate::convert::{hash_rust, into_py_variant, variant_repr};
use crate::lattice::impl_py_lattice;
use crate::molecule::Molecule;
use crate::noncovalent::NoncovalentBondForm;

/// The key (identity) of a noncovalent-bond constraint, for keyed lookup. The single
/// key `Intramolecular` is the bare discriminant (no sub-key).
#[pyclass]
pub enum NoncovalentBondConstraintKey {
    Intramolecular(),
}

#[pymethods]
impl NoncovalentBondConstraintKey {
    pub(crate) fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    pub(crate) fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            NoncovalentBondConstraintKey::Intramolecular() => ("Intramolecular", 0),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "NoncovalentBondConstraintKey",
            variant,
            arity,
        )
    }
}

impl NoncovalentBondConstraintKey {
    pub(crate) fn from_rust(key: &GraphIrNoncovalentBondConstraintKey) -> Self {
        match key {
            GraphIrNoncovalentBondConstraintKey::Intramolecular => Self::Intramolecular(),
        }
    }

    pub(crate) fn to_rust(&self) -> GraphIrNoncovalentBondConstraintKey {
        match self {
            Self::Intramolecular() => GraphIrNoncovalentBondConstraintKey::Intramolecular,
        }
    }
}

/// A noncovalent-bond-scope constraint: whether the bond is intramolecular (a boolean
/// value; `Undetermined` when unspecified).
#[pyclass(frozen)]
pub enum NoncovalentBondConstraintForm {
    Intramolecular(Py<BooleanForm>),
}

#[pymethods]
impl NoncovalentBondConstraintForm {
    /// The constraint's key (identity).
    #[getter]
    pub(crate) fn key(&self, py: Python<'_>) -> NoncovalentBondConstraintKey {
        NoncovalentBondConstraintKey::from_rust(&self.to_rust(py).key())
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            NoncovalentBondConstraintForm::Intramolecular(_) => "Intramolecular",
        };
        variant_repr(
            slf.bind(py).as_any(),
            "NoncovalentBondConstraintForm",
            variant,
            1,
        )
    }
}

impl_py_lattice!(
    NoncovalentBondConstraintForm,
    GraphIrNoncovalentBondConstraintForm,
    |value: &NoncovalentBondConstraintForm,
     py: Python<'_>|
     -> PyResult<GraphIrNoncovalentBondConstraintForm> { Ok(value.to_rust(py)) },
    |py: Python<'_>,
     value: GraphIrNoncovalentBondConstraintForm|
     -> PyResult<NoncovalentBondConstraintForm> {
        NoncovalentBondConstraintForm::from_rust(py, &value)
    }
);

impl NoncovalentBondConstraintForm {
    pub(crate) fn from_rust(
        py: Python<'_>,
        form: &GraphIrNoncovalentBondConstraintForm,
    ) -> PyResult<Self> {
        Ok(match form {
            GraphIrNoncovalentBondConstraintForm::Intramolecular(b) => {
                Self::Intramolecular(into_py_variant(py, BooleanForm::from_rust(b))?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrNoncovalentBondConstraintForm {
        match self {
            Self::Intramolecular(b) => {
                GraphIrNoncovalentBondConstraintForm::Intramolecular(b.bind(py).borrow().to_rust())
            }
        }
    }
}

/// The argument to `update`: another constraint container, a live view, or an
/// iterable of `NoncovalentBondConstraintForm` (each `set`, last-wins).
#[derive(FromPyObject)]
pub(crate) enum NoncovalentBondConstraintsUpdate {
    Container(Py<NoncovalentBondConstraintsForm>),
    View(Py<NoncovalentBondConstraintsView>),
    Entries(Vec<Py<NoncovalentBondConstraintForm>>),
}

impl NoncovalentBondConstraintsUpdate {
    /// Read every Python object into owned data — no write target is touched. Callers
    /// resolve *before* taking the write borrow so the re-entrant read of a view (or
    /// container) that aliases the same bond happens while nothing is borrowed
    /// (otherwise `bond.constraints.update(bond.constraints)` self-aliases into a
    /// RefCell double-borrow panic).
    pub(crate) fn resolve(
        &self,
        py: Python<'_>,
    ) -> PyResult<ResolvedNoncovalentBondConstraintsUpdate> {
        Ok(match self {
            NoncovalentBondConstraintsUpdate::Container(c) => {
                ResolvedNoncovalentBondConstraintsUpdate::Overlay(
                    c.bind(py).borrow().to_rust().clone(),
                )
            }
            NoncovalentBondConstraintsUpdate::View(v) => {
                ResolvedNoncovalentBondConstraintsUpdate::Overlay(
                    v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?,
                )
            }
            NoncovalentBondConstraintsUpdate::Entries(entries) => {
                ResolvedNoncovalentBondConstraintsUpdate::Entries(
                    entries
                        .iter()
                        .map(|entry| entry.bind(py).borrow().to_rust(py))
                        .collect(),
                )
            }
        })
    }
}

/// A `NoncovalentBondConstraintsUpdate` with all Python-object reads already done, so
/// it can be applied under a write borrow without re-entering Python.
pub(crate) enum ResolvedNoncovalentBondConstraintsUpdate {
    /// A whole container (from another container or a live view): overlaid via
    /// `update` (last-wins per key; undetermined entries remove).
    Overlay(GraphIrNoncovalentBondConstraintsForm),
    /// Loose entries: `set` each (last-wins; undetermined entries stored, not removed).
    Entries(Vec<GraphIrNoncovalentBondConstraintForm>),
}

impl ResolvedNoncovalentBondConstraintsUpdate {
    /// Overlay onto `target` in place. No Python reads.
    pub(crate) fn apply(self, target: &mut GraphIrNoncovalentBondConstraintsForm) {
        match self {
            ResolvedNoncovalentBondConstraintsUpdate::Overlay(overlay) => target.update(&overlay),
            ResolvedNoncovalentBondConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry);
                }
            }
        }
    }
}

/// The noncovalent-bond-scope constraints on a noncovalent bond, in kind-sorted order.
/// Mutable, hence value-equal but unhashable (matching `NoncovalentBondForm`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct NoncovalentBondConstraintsForm(GraphIrNoncovalentBondConstraintsForm);

#[pymethods]
impl NoncovalentBondConstraintsForm {
    /// Build from a sequence of constraints (a later entry of the same key replaces
    /// an earlier one, last-wins).
    #[new]
    pub(crate) fn new(py: Python<'_>, entries: Vec<Py<NoncovalentBondConstraintForm>>) -> Self {
        let mut constraints = GraphIrNoncovalentBondConstraintsForm::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py)),
        );
        NoncovalentBondConstraintsForm(constraints)
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let value = into_py_variant(py, NoncovalentBondConstraintForm::from_rust(py, entry)?)?;
            parts.push(value.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!(
            "NoncovalentBondConstraintsForm([{}])",
            parts.join(", ")
        ))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    pub(crate) fn set(&mut self, py: Python<'_>, c: Py<NoncovalentBondConstraintForm>) {
        self.0.set(c.bind(py).borrow().to_rust(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    pub(crate) fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<Option<NoncovalentBondConstraintForm>> {
        self.0
            .remove(key.bind(py).borrow().to_rust())
            .map(|c| NoncovalentBondConstraintForm::from_rust(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container, a live view, or an
    /// iterable of `NoncovalentBondConstraintForm` (last-wins per key; undetermined
    /// entries remove). Takes `slf` by handle so `other` is fully read *before* the
    /// write borrow — `cs.update(cs)` on the same container is then a no-op, not a
    /// double-borrow panic.
    pub(crate) fn update(
        slf: Py<Self>,
        py: Python<'_>,
        other: NoncovalentBondConstraintsUpdate,
    ) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        resolved.apply(&mut slf.borrow_mut(py).0);
        Ok(())
    }

    pub(crate) fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintKeyIter> {
        noncovalent_bond_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintKeyIter> {
        noncovalent_bond_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintIter> {
        noncovalent_bond_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintItemsIter> {
        noncovalent_bond_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
        &self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_rust()) {
            Some(constraint) => Ok(into_py_variant(
                py,
                NoncovalentBondConstraintForm::from_rust(py, constraint)?,
            )?
            .into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    pub(crate) fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<NoncovalentBondConstraintForm> {
        match self.0.get(key.bind(py).borrow().to_rust()) {
            Some(constraint) => NoncovalentBondConstraintForm::from_rust(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
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
        key: Py<NoncovalentBondConstraintKey>,
    ) -> bool {
        self.0.contains(key.bind(py).borrow().to_rust())
    }

    /// Whether the bond is intramolecular; `Undetermined` when no `Intramolecular`
    /// constraint is present (matching the non-optional Rust accessor).
    #[getter]
    pub(crate) fn intramolecular(&self) -> BooleanForm {
        BooleanForm::from_rust(&self.0.intramolecular())
    }

    #[setter]
    pub(crate) fn set_intramolecular(&mut self, py: Python<'_>, value: BooleanLike) {
        self.0
            .set(GraphIrNoncovalentBondConstraintForm::intramolecular(
                value.to_rust(py),
            ));
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// Python values.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        noncovalent_bond_constraints_asdict(py, &self.0)
    }
}

impl NoncovalentBondConstraintsForm {
    /// The wrapped IR constraints — read access for noncovalent bond construction.
    pub(crate) fn to_rust(&self) -> &GraphIrNoncovalentBondConstraintsForm {
        &self.0
    }

    /// Wrap IR constraints (the hold-the-value `from_rust` bridge). Test-only —
    /// in-crate construction wraps `NoncovalentBondConstraintsForm(..)` directly.
    pub(crate) fn from_rust(constraints: GraphIrNoncovalentBondConstraintsForm) -> Self {
        NoncovalentBondConstraintsForm(constraints)
    }
}

impl_py_lattice!(
    NoncovalentBondConstraintsForm,
    GraphIrNoncovalentBondConstraintsForm,
    |value: &NoncovalentBondConstraintsForm,
     _py: Python<'_>|
     -> PyResult<GraphIrNoncovalentBondConstraintsForm> { Ok(value.to_rust().clone()) },
    |_py: Python<'_>,
     value: GraphIrNoncovalentBondConstraintsForm|
     -> PyResult<NoncovalentBondConstraintsForm> { Ok(NoncovalentBondConstraintsForm(value)) }
);

/// Build the per-constraint iterator handle from a borrowed container.
pub(crate) fn noncovalent_bond_constraints_iter(
    py: Python<'_>,
    constraints: &GraphIrNoncovalentBondConstraintsForm,
) -> PyResult<NoncovalentBondConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(
                py,
                NoncovalentBondConstraintForm::from_rust(py, constraint)?,
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(NoncovalentBondConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
pub(crate) fn noncovalent_bond_constraint_keys(
    py: Python<'_>,
    constraints: &GraphIrNoncovalentBondConstraintsForm,
) -> PyResult<NoncovalentBondConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(
                py,
                NoncovalentBondConstraintKey::from_rust(&constraint.key()),
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(NoncovalentBondConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
pub(crate) fn noncovalent_bond_constraint_items(
    py: Python<'_>,
    constraints: &GraphIrNoncovalentBondConstraintsForm,
) -> PyResult<NoncovalentBondConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(
                    py,
                    NoncovalentBondConstraintKey::from_rust(&constraint.key()),
                )?,
                into_py_variant(
                    py,
                    NoncovalentBondConstraintForm::from_rust(py, constraint)?,
                )?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(NoncovalentBondConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// Python values.
pub(crate) fn noncovalent_bond_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &GraphIrNoncovalentBondConstraintsForm,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            GraphIrNoncovalentBondConstraintForm::Intramolecular(b) => {
                dict.set_item("intramolecular", BooleanForm::from_rust(b))?
            }
        }
    }
    Ok(dict)
}

#[pyclass]
pub(crate) struct NoncovalentBondConstraintIter {
    entries: IntoIter<Py<NoncovalentBondConstraintForm>>,
}

#[pymethods]
impl NoncovalentBondConstraintIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<NoncovalentBondConstraintForm>> {
        self.entries.next()
    }
}

#[pyclass]
pub(crate) struct NoncovalentBondConstraintKeyIter {
    keys: IntoIter<Py<NoncovalentBondConstraintKey>>,
}

#[pymethods]
impl NoncovalentBondConstraintKeyIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<NoncovalentBondConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
pub(crate) struct NoncovalentBondConstraintItemsIter {
    items: IntoIter<(
        Py<NoncovalentBondConstraintKey>,
        Py<NoncovalentBondConstraintForm>,
    )>,
}

#[pymethods]
impl NoncovalentBondConstraintItemsIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(
        &mut self,
    ) -> Option<(
        Py<NoncovalentBondConstraintKey>,
        Py<NoncovalentBondConstraintForm>,
    )> {
        self.items.next()
    }
}

/// A whole-container argument that snapshots either a value container or a live view
/// — for the noncovalent bond `constraints` setter, which accepts either.
#[derive(FromPyObject)]
pub(crate) enum NoncovalentBondConstraintsLike {
    Container(Py<NoncovalentBondConstraintsForm>),
    View(Py<NoncovalentBondConstraintsView>),
}

impl NoncovalentBondConstraintsLike {
    pub(crate) fn to_rust(
        &self,
        py: Python<'_>,
    ) -> PyResult<GraphIrNoncovalentBondConstraintsForm> {
        match self {
            NoncovalentBondConstraintsLike::Container(c) => {
                Ok(c.bind(py).borrow().to_rust().clone())
            }
            NoncovalentBondConstraintsLike::View(v) => {
                v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))
            }
        }
    }
}

/// What a `NoncovalentBondConstraintsView` writes through to: a noncovalent bond within
/// a molecule (by index) or a standalone `NoncovalentBondForm`.
pub(crate) enum NoncovalentBondConstraintsBacking {
    Molecule {
        owner: Py<Molecule>,
        id: GraphIrNoncovalentBondId,
    },
    Noncovalent(Py<NoncovalentBondForm>),
}

/// A live handle onto one noncovalent bond's constraints, backed by either a
/// molecule-bond or a standalone `NoncovalentBondForm`. Reads borrow the constraints and
/// read only the item they need (no whole-container clone); mutators write through to the
/// bond in place, without a clone-and-writeback.
#[pyclass]
pub struct NoncovalentBondConstraintsView {
    pub(crate) backing: NoncovalentBondConstraintsBacking,
}

impl NoncovalentBondConstraintsView {
    /// Borrow the backing bond's constraints and read one item through `f` — no clone.
    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&GraphIrNoncovalentBondConstraintsForm) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            NoncovalentBondConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .to_rust()
                    .noncovalent_bonds()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("noncovalent bond id out of range"))?;
                f(&view.attributes.constraints)
            }
            NoncovalentBondConstraintsBacking::Noncovalent(bond) => {
                let bond = bond.bind(py).borrow();
                f(&bond.to_rust().constraints)
            }
        }
    }

    /// Mutate the backing bond's constraints in place through `f`.
    pub(crate) fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut GraphIrNoncovalentBondConstraintsForm) -> R,
    ) -> PyResult<R> {
        match &self.backing {
            NoncovalentBondConstraintsBacking::Molecule { owner, id } => Ok(f(&mut owner
                .borrow_mut(py)
                .to_rust_mut()
                .noncovalent_bond_mut(*id)
                .attributes
                .constraints)),
            NoncovalentBondConstraintsBacking::Noncovalent(bond) => {
                Ok(f(&mut bond.borrow_mut(py).to_rust_mut()?.constraints))
            }
        }
    }

    /// Set one constraint on the backing bond in place (last-wins per key).
    pub(crate) fn set_form(
        &self,
        py: Python<'_>,
        constraint: GraphIrNoncovalentBondConstraintForm,
    ) -> PyResult<()> {
        self.with_mut(py, |cs| cs.set(constraint))
    }

    /// Remove one key from the backing bond in place, returning the removed entry.
    pub(crate) fn remove_form(
        &self,
        py: Python<'_>,
        key: GraphIrNoncovalentBondConstraintKey,
    ) -> PyResult<Option<GraphIrNoncovalentBondConstraintForm>> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl NoncovalentBondConstraintsView {
    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("NoncovalentBondConstraintsView({count} entries)"))
    }

    /// Insert `c` on the bond in place, replacing any existing entry of the same key
    /// (last-wins).
    pub(crate) fn set(&self, py: Python<'_>, c: Py<NoncovalentBondConstraintForm>) -> PyResult<()> {
        self.set_form(py, c.bind(py).borrow().to_rust(py))
    }

    /// Remove the entry with the given key from the bond in place, returning it if
    /// present (dict `pop`).
    pub(crate) fn pop(
        &self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<Option<NoncovalentBondConstraintForm>> {
        self.remove_form(py, key.bind(py).borrow().to_rust())?
            .map(|c| NoncovalentBondConstraintForm::from_rust(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(
        &self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<()> {
        if self
            .remove_form(py, key.bind(py).borrow().to_rust())?
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
    /// view, or an iterable of `NoncovalentBondConstraintForm` (last-wins per key;
    /// undetermined entries remove).
    pub(crate) fn update(
        &self,
        py: Python<'_>,
        other: NoncovalentBondConstraintsUpdate,
    ) -> PyResult<()> {
        // Resolve `other` to owned data *before* the write borrow, so a view aliasing
        // the same bond (`bond.constraints.update(bond.constraints)`) reads while the
        // bond is unborrowed instead of self-aliasing into a double-borrow panic.
        let resolved = other.resolve(py)?;
        self.with_mut(py, |cs| resolved.apply(cs))
    }

    pub(crate) fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintKeyIter> {
        self.read(py, |cs| noncovalent_bond_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintKeyIter> {
        self.read(py, |cs| noncovalent_bond_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintIter> {
        self.read(py, |cs| noncovalent_bond_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintItemsIter> {
        self.read(py, |cs| noncovalent_bond_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
        &self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_rust();
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| NoncovalentBondConstraintForm::from_rust(py, constraint))
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
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<NoncovalentBondConstraintForm> {
        let rust_key = key.bind(py).borrow().to_rust();
        let found = self.read(py, |cs| {
            cs.get(rust_key)
                .map(|constraint| NoncovalentBondConstraintForm::from_rust(py, constraint))
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
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_rust();
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// Whether the bond is intramolecular; `Undetermined` when no `Intramolecular`
    /// constraint is present (matching the non-optional Rust accessor).
    #[getter]
    pub(crate) fn intramolecular(&self, py: Python<'_>) -> PyResult<BooleanForm> {
        self.read(py, |cs| Ok(BooleanForm::from_rust(&cs.intramolecular())))
    }

    #[setter]
    pub(crate) fn set_intramolecular(&self, py: Python<'_>, value: BooleanLike) -> PyResult<()> {
        self.set_form(
            py,
            GraphIrNoncovalentBondConstraintForm::intramolecular(value.to_rust(py)),
        )
    }

    /// The present constraints as a dict keyed by snake_case name.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| noncovalent_bond_constraints_asdict(py, cs))
    }
}
