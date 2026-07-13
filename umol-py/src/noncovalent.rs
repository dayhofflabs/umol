//! Noncovalent-bond bindings. This file grows over the B5 slice; it opens with the
//! kind leaf: `NoncovalentBondKind` (the interaction kind) and `NoncovalentBondKindAst`
//! (`Undetermined | Lit(kind)`), mirroring `umol_ast::ast::{NoncovalentBondKind,
//! NoncovalentBondKindAst}` — the noncovalent analog of `atom.element: ElementAst`
//! over the `Element` value enum.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use umol_ast::ast::{
    AsLit, AtomId as AstAtomId, MoleculeAst as AstMoleculeAst,
    NoncovalentBondAst as AstNoncovalentBondAst,
    NoncovalentBondConstraintAst as AstNoncovalentBondConstraintAst,
    NoncovalentBondConstraintKey as AstNoncovalentBondConstraintKey,
    NoncovalentBondConstraintsAst as AstNoncovalentBondConstraintsAst,
    NoncovalentBondId as AstNoncovalentBondId, NoncovalentBondKind as AstNoncovalentBondKind,
    NoncovalentBondKindAst as AstNoncovalentBondKindAst,
    NoncovalentBondView as AstNoncovalentBondView,
};

use crate::boolean::{BooleanArg, BooleanAst};
use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::error::parse_error;
use crate::molecule::MoleculeAst;

/// A noncovalent interaction kind. A fieldless, hashable value enum whose members
/// mirror the Rust `NoncovalentBondKind` exactly.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondKind {
    HydrogenBond,
    HalogenBond,
    ChalcogenBond,
    Ionic,
    VanDerWaals,
}

impl NoncovalentBondKind {
    pub(crate) fn from_ast(ast: AstNoncovalentBondKind) -> Self {
        match ast {
            AstNoncovalentBondKind::HydrogenBond => Self::HydrogenBond,
            AstNoncovalentBondKind::HalogenBond => Self::HalogenBond,
            AstNoncovalentBondKind::ChalcogenBond => Self::ChalcogenBond,
            AstNoncovalentBondKind::Ionic => Self::Ionic,
            AstNoncovalentBondKind::VanDerWaals => Self::VanDerWaals,
        }
    }

    pub(crate) fn to_ast(self) -> AstNoncovalentBondKind {
        match self {
            Self::HydrogenBond => AstNoncovalentBondKind::HydrogenBond,
            Self::HalogenBond => AstNoncovalentBondKind::HalogenBond,
            Self::ChalcogenBond => AstNoncovalentBondKind::ChalcogenBond,
            Self::Ionic => AstNoncovalentBondKind::Ionic,
            Self::VanDerWaals => AstNoncovalentBondKind::VanDerWaals,
        }
    }
}

/// A noncovalent bond's interaction kind: undetermined, or a concrete
/// `NoncovalentBondKind`. Mirrors `NoncovalentBondKindAst`.
#[pyclass]
pub enum NoncovalentBondKindAst {
    Undetermined(),
    Lit(NoncovalentBondKind),
}

#[pymethods]
impl NoncovalentBondKindAst {
    /// The concrete interaction kind, or `None` when undetermined.
    fn as_lit(&self) -> Option<NoncovalentBondKind> {
        self.to_ast().as_lit().map(NoncovalentBondKind::from_ast)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            NoncovalentBondKindAst::Undetermined() => ("Undetermined", 0),
            NoncovalentBondKindAst::Lit(_) => ("Lit", 1),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "NoncovalentBondKindAst",
            variant,
            arity,
        )
    }
}

impl NoncovalentBondKindAst {
    pub(crate) fn from_ast(ast: &AstNoncovalentBondKindAst) -> Self {
        match ast {
            AstNoncovalentBondKindAst::Undetermined => Self::Undetermined(),
            AstNoncovalentBondKindAst::Lit(k) => Self::Lit(NoncovalentBondKind::from_ast(*k)),
        }
    }

    pub(crate) fn to_ast(&self) -> AstNoncovalentBondKindAst {
        match self {
            Self::Undetermined() => AstNoncovalentBondKindAst::Undetermined,
            Self::Lit(k) => AstNoncovalentBondKindAst::Lit(k.to_ast()),
        }
    }
}

/// Setter coercion for a noncovalent `kind` field: a bare `NoncovalentBondKind` →
/// `Lit`, or a `NoncovalentBondKindAst` passthrough (mirroring the `Undetermined |
/// Lit` structure).
#[derive(FromPyObject)]
pub(crate) enum NoncovalentBondKindArg {
    Kind(NoncovalentBondKind),
    Ast(Py<NoncovalentBondKindAst>),
}

impl NoncovalentBondKindArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstNoncovalentBondKindAst {
        match self {
            NoncovalentBondKindArg::Kind(k) => AstNoncovalentBondKindAst::Lit(k.to_ast()),
            NoncovalentBondKindArg::Ast(a) => a.bind(py).borrow().to_ast(),
        }
    }
}

/// The key (identity) of a noncovalent-bond constraint, for keyed lookup. The single
/// key `Intramolecular` is the bare discriminant (no sub-key).
#[pyclass]
pub enum NoncovalentBondConstraintKey {
    Intramolecular(),
}

#[pymethods]
impl NoncovalentBondConstraintKey {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
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
    pub(crate) fn from_ast(ast: &AstNoncovalentBondConstraintKey) -> Self {
        match ast {
            AstNoncovalentBondConstraintKey::Intramolecular => Self::Intramolecular(),
        }
    }

    pub(crate) fn to_ast(&self) -> AstNoncovalentBondConstraintKey {
        match self {
            Self::Intramolecular() => AstNoncovalentBondConstraintKey::Intramolecular,
        }
    }
}

/// A noncovalent-bond-scope constraint: whether the bond is intramolecular (a boolean
/// value; `Undetermined` when unspecified).
#[pyclass]
pub enum NoncovalentBondConstraintAst {
    Intramolecular(Py<BooleanAst>),
}

#[pymethods]
impl NoncovalentBondConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    fn key(&self, py: Python<'_>) -> NoncovalentBondConstraintKey {
        NoncovalentBondConstraintKey::from_ast(&self.to_ast(py).key())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            NoncovalentBondConstraintAst::Intramolecular(_) => "Intramolecular",
        };
        variant_repr(
            slf.bind(py).as_any(),
            "NoncovalentBondConstraintAst",
            variant,
            1,
        )
    }
}

impl NoncovalentBondConstraintAst {
    pub(crate) fn from_ast(
        py: Python<'_>,
        ast: &AstNoncovalentBondConstraintAst,
    ) -> PyResult<Self> {
        Ok(match ast {
            AstNoncovalentBondConstraintAst::Intramolecular(b) => {
                Self::Intramolecular(into_py_variant(py, BooleanAst::from_ast(b))?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstNoncovalentBondConstraintAst {
        match self {
            Self::Intramolecular(b) => {
                AstNoncovalentBondConstraintAst::Intramolecular(b.bind(py).borrow().to_ast())
            }
        }
    }
}

/// The argument to `update`: another constraint container, a live view, or an
/// iterable of `NoncovalentBondConstraintAst` (each `set`, last-wins).
#[derive(FromPyObject)]
enum NoncovalentBondConstraintsUpdate {
    Container(Py<NoncovalentBondConstraintsAst>),
    View(Py<NoncovalentBondConstraintsView>),
    Entries(Vec<Py<NoncovalentBondConstraintAst>>),
}

impl NoncovalentBondConstraintsUpdate {
    /// Read every Python object into owned data — no write target is touched. Callers
    /// resolve *before* taking the write borrow so the re-entrant read of a view (or
    /// container) that aliases the same bond happens while nothing is borrowed
    /// (otherwise `bond.constraints.update(bond.constraints)` self-aliases into a
    /// RefCell double-borrow panic).
    fn resolve(&self, py: Python<'_>) -> PyResult<ResolvedNoncovalentBondConstraintsUpdate> {
        Ok(match self {
            NoncovalentBondConstraintsUpdate::Container(c) => {
                ResolvedNoncovalentBondConstraintsUpdate::Overlay(
                    c.bind(py).borrow().inner().clone(),
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
                        .map(|entry| entry.bind(py).borrow().to_ast(py))
                        .collect(),
                )
            }
        })
    }
}

/// A `NoncovalentBondConstraintsUpdate` with all Python-object reads already done, so
/// it can be applied under a write borrow without re-entering Python.
enum ResolvedNoncovalentBondConstraintsUpdate {
    /// A whole container (from another container or a live view): overlaid via
    /// `update` (last-wins per key; undetermined entries remove).
    Overlay(AstNoncovalentBondConstraintsAst),
    /// Loose entries: `set` each (last-wins; undetermined entries stored, not removed).
    Entries(Vec<AstNoncovalentBondConstraintAst>),
}

impl ResolvedNoncovalentBondConstraintsUpdate {
    /// Overlay onto `target` in place. No Python reads.
    fn apply(self, target: &mut AstNoncovalentBondConstraintsAst) {
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
/// Mutable, hence value-equal but unhashable (matching `NoncovalentBondAst`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct NoncovalentBondConstraintsAst(AstNoncovalentBondConstraintsAst);

#[pymethods]
impl NoncovalentBondConstraintsAst {
    /// Build from a sequence of constraints (a later entry of the same key replaces
    /// an earlier one, last-wins).
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<NoncovalentBondConstraintAst>>) -> Self {
        let mut constraints = AstNoncovalentBondConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_ast(py)),
        );
        NoncovalentBondConstraintsAst(constraints)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let mirror = into_py_variant(py, NoncovalentBondConstraintAst::from_ast(py, entry)?)?;
            parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!(
            "NoncovalentBondConstraintsAst([{}])",
            parts.join(", ")
        ))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    fn set(&mut self, py: Python<'_>, c: Py<NoncovalentBondConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<Option<NoncovalentBondConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_ast())
            .map(|c| NoncovalentBondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container, a live view, or an
    /// iterable of `NoncovalentBondConstraintAst` (last-wins per key; undetermined
    /// entries remove). Takes `slf` by handle so `other` is fully read *before* the
    /// write borrow — `cs.update(cs)` on the same container is then a no-op, not a
    /// double-borrow panic.
    fn update(
        slf: Py<Self>,
        py: Python<'_>,
        other: NoncovalentBondConstraintsUpdate,
    ) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        resolved.apply(&mut slf.borrow_mut(py).0);
        Ok(())
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintKeyIter> {
        noncovalent_bond_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintKeyIter> {
        noncovalent_bond_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintIter> {
        noncovalent_bond_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintItemsIter> {
        noncovalent_bond_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_ast()) {
            Some(constraint) => Ok(into_py_variant(
                py,
                NoncovalentBondConstraintAst::from_ast(py, constraint)?,
            )?
            .into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<NoncovalentBondConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_ast()) {
            Some(constraint) => NoncovalentBondConstraintAst::from_ast(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_ast()).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<NoncovalentBondConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast())
    }

    /// Whether the bond is intramolecular; `Undetermined` when no `Intramolecular`
    /// constraint is present (mirroring the non-optional Rust accessor).
    #[getter]
    fn intramolecular(&self) -> BooleanAst {
        BooleanAst::from_ast(&self.0.intramolecular())
    }

    #[setter]
    fn set_intramolecular(&mut self, py: Python<'_>, value: BooleanArg) {
        self.0.set(AstNoncovalentBondConstraintAst::intramolecular(
            value.to_ast(py),
        ));
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        noncovalent_bond_constraints_asdict(py, &self.0)
    }
}

impl NoncovalentBondConstraintsAst {
    /// The wrapped AST constraints — read access for noncovalent bond construction.
    pub(crate) fn inner(&self) -> &AstNoncovalentBondConstraintsAst {
        &self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge). Test-only —
    /// in-crate construction wraps `NoncovalentBondConstraintsAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(constraints: AstNoncovalentBondConstraintsAst) -> Self {
        NoncovalentBondConstraintsAst(constraints)
    }
}

/// Build the per-constraint iterator handle from a borrowed container.
fn noncovalent_bond_constraints_iter(
    py: Python<'_>,
    constraints: &AstNoncovalentBondConstraintsAst,
) -> PyResult<NoncovalentBondConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(py, NoncovalentBondConstraintAst::from_ast(py, constraint)?)
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(NoncovalentBondConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
fn noncovalent_bond_constraint_keys(
    py: Python<'_>,
    constraints: &AstNoncovalentBondConstraintsAst,
) -> PyResult<NoncovalentBondConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(
                py,
                NoncovalentBondConstraintKey::from_ast(&constraint.key()),
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(NoncovalentBondConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
fn noncovalent_bond_constraint_items(
    py: Python<'_>,
    constraints: &AstNoncovalentBondConstraintsAst,
) -> PyResult<NoncovalentBondConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(
                    py,
                    NoncovalentBondConstraintKey::from_ast(&constraint.key()),
                )?,
                into_py_variant(py, NoncovalentBondConstraintAst::from_ast(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(NoncovalentBondConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// inner-value mirrors.
fn noncovalent_bond_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &AstNoncovalentBondConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            AstNoncovalentBondConstraintAst::Intramolecular(b) => {
                dict.set_item("intramolecular", BooleanAst::from_ast(b))?
            }
        }
    }
    Ok(dict)
}

#[pyclass]
struct NoncovalentBondConstraintIter {
    entries: IntoIter<Py<NoncovalentBondConstraintAst>>,
}

#[pymethods]
impl NoncovalentBondConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<NoncovalentBondConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
struct NoncovalentBondConstraintKeyIter {
    keys: IntoIter<Py<NoncovalentBondConstraintKey>>,
}

#[pymethods]
impl NoncovalentBondConstraintKeyIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<NoncovalentBondConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
struct NoncovalentBondConstraintItemsIter {
    items: IntoIter<(
        Py<NoncovalentBondConstraintKey>,
        Py<NoncovalentBondConstraintAst>,
    )>,
}

#[pymethods]
impl NoncovalentBondConstraintItemsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(
        &mut self,
    ) -> Option<(
        Py<NoncovalentBondConstraintKey>,
        Py<NoncovalentBondConstraintAst>,
    )> {
        self.items.next()
    }
}

/// A noncovalent bond: an interaction `kind` plus noncovalent-bond-scope constraints.
/// No bond order, charge, or spin — these do not apply to noncovalent interactions.
/// The endpoint atom pair is the owning molecule's relation topology (the view half),
/// not part of the value.
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct NoncovalentBondAst(AstNoncovalentBondAst);

#[pymethods]
impl NoncovalentBondAst {
    /// Construct from an interaction kind — a `NoncovalentBondKind` or a
    /// `NoncovalentBondKindAst` — optionally setting constraints.
    #[new]
    #[pyo3(signature = (kind, *, constraints=None))]
    fn new(
        py: Python<'_>,
        kind: NoncovalentBondKindArg,
        constraints: Option<Py<NoncovalentBondConstraintsAst>>,
    ) -> Self {
        let mut bond = AstNoncovalentBondAst::new(kind.to_ast(py));
        if let Some(constraints) = constraints {
            bond.constraints = constraints.bind(py).borrow().inner().clone();
        }
        NoncovalentBondAst(bond)
    }

    /// Parse a noncovalent-bond-DSL string (e.g. `"Hbd#I"`) into a `NoncovalentBondAst`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        AstNoncovalentBondAst::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("NoncovalentBondAst.parse('{}')", self.0)
    }

    /// The interaction kind.
    #[getter]
    fn kind(&self) -> NoncovalentBondKindAst {
        NoncovalentBondKindAst::from_ast(&self.0.kind)
    }

    #[setter]
    fn set_kind(&mut self, py: Python<'_>, value: NoncovalentBondKindArg) {
        self.0.kind = value.to_ast(py);
    }

    /// The bond's constraints as a live handle onto this bond: reads borrow the
    /// current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(slf: Py<Self>) -> NoncovalentBondConstraintsView {
        NoncovalentBondConstraintsView {
            backing: NoncovalentBondConstraintsBacking::Noncovalent(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or a
    /// live view. Takes `slf` by handle and snapshots `value` *before* the write borrow,
    /// so `bond.constraints = bond.constraints` (a view over the same bond) reads while
    /// the bond is unborrowed instead of self-aliasing into a double-borrow panic.
    #[setter]
    fn set_constraints(
        slf: Py<Self>,
        py: Python<'_>,
        value: NoncovalentBondConstraintsArg,
    ) -> PyResult<()> {
        let snapshot = value.to_ast(py)?;
        slf.borrow_mut(py).0.constraints = snapshot;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are the field mirrors.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("kind", self.kind())?;
        dict.set_item(
            "constraints",
            noncovalent_bond_constraints_asdict(py, &self.0.constraints)?,
        )?;
        Ok(dict)
    }
}

impl NoncovalentBondAst {
    /// The wrapped AST bond — read access for the bond-backed constraints view.
    pub(crate) fn inner(&self) -> &AstNoncovalentBondAst {
        &self.0
    }

    /// Mutable access to the wrapped AST bond — write access for the bond-backed
    /// constraints view.
    pub(crate) fn inner_mut(&mut self) -> &mut AstNoncovalentBondAst {
        &mut self.0
    }

    /// Wrap an AST bond (the hold-the-value `from_inner` bridge, paired with `inner`).
    /// Test-only — in-crate construction wraps `NoncovalentBondAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(bond: AstNoncovalentBondAst) -> Self {
        NoncovalentBondAst(bond)
    }
}

/// A whole-container argument that snapshots either a value container or a live view
/// — for the noncovalent bond `constraints` setter, which accepts either.
#[derive(FromPyObject)]
enum NoncovalentBondConstraintsArg {
    Container(Py<NoncovalentBondConstraintsAst>),
    View(Py<NoncovalentBondConstraintsView>),
}

impl NoncovalentBondConstraintsArg {
    fn to_ast(&self, py: Python<'_>) -> PyResult<AstNoncovalentBondConstraintsAst> {
        match self {
            NoncovalentBondConstraintsArg::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
            NoncovalentBondConstraintsArg::View(v) => {
                v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))
            }
        }
    }
}

/// What a `NoncovalentBondConstraintsView` writes through to: a noncovalent bond within
/// a molecule (by index) or a standalone `NoncovalentBondAst`.
enum NoncovalentBondConstraintsBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: AstNoncovalentBondId,
    },
    Noncovalent(Py<NoncovalentBondAst>),
}

/// A live handle onto one noncovalent bond's constraints, backed by either a
/// molecule-bond or a standalone `NoncovalentBondAst`. Reads borrow the constraints and
/// read only the item they need (no whole-container clone); mutators write through to the
/// bond in place, without a clone-and-writeback.
#[pyclass]
pub struct NoncovalentBondConstraintsView {
    backing: NoncovalentBondConstraintsBacking,
}

impl NoncovalentBondConstraintsView {
    /// Borrow the backing bond's constraints and read one item through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstNoncovalentBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            NoncovalentBondConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .noncovalent_bonds()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("noncovalent bond id out of range"))?;
                f(&view.ast.constraints)
            }
            NoncovalentBondConstraintsBacking::Noncovalent(bond) => {
                let bond = bond.bind(py).borrow();
                f(&bond.inner().constraints)
            }
        }
    }

    /// Mutate the backing bond's constraints in place through `f`.
    fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut AstNoncovalentBondConstraintsAst) -> R,
    ) -> R {
        match &self.backing {
            NoncovalentBondConstraintsBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .noncovalent_bond_mut(*id)
                .ast
                .constraints),
            NoncovalentBondConstraintsBacking::Noncovalent(bond) => {
                f(&mut bond.borrow_mut(py).inner_mut().constraints)
            }
        }
    }

    /// Set one constraint on the backing bond in place (last-wins per key).
    fn set_ast(&self, py: Python<'_>, constraint: AstNoncovalentBondConstraintAst) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing bond in place, returning the removed entry.
    fn remove_ast(
        &self,
        py: Python<'_>,
        key: AstNoncovalentBondConstraintKey,
    ) -> Option<AstNoncovalentBondConstraintAst> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl NoncovalentBondConstraintsView {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("NoncovalentBondConstraintsView({count} entries)"))
    }

    /// Insert `c` on the bond in place, replacing any existing entry of the same key
    /// (last-wins).
    fn set(&self, py: Python<'_>, c: Py<NoncovalentBondConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key from the bond in place, returning it if
    /// present (dict `pop`).
    fn pop(
        &self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<Option<NoncovalentBondConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_ast())
            .map(|c| NoncovalentBondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&self, py: Python<'_>, key: Py<NoncovalentBondConstraintKey>) -> PyResult<()> {
        if self
            .remove_ast(py, key.bind(py).borrow().to_ast())
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
    /// view, or an iterable of `NoncovalentBondConstraintAst` (last-wins per key;
    /// undetermined entries remove).
    fn update(&self, py: Python<'_>, other: NoncovalentBondConstraintsUpdate) -> PyResult<()> {
        // Resolve `other` to owned data *before* the write borrow, so a view aliasing
        // the same bond (`bond.constraints.update(bond.constraints)`) reads while the
        // bond is unborrowed instead of self-aliasing into a double-borrow panic.
        let resolved = other.resolve(py)?;
        self.with_mut(py, |cs| resolved.apply(cs));
        Ok(())
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintKeyIter> {
        self.read(py, |cs| noncovalent_bond_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintKeyIter> {
        self.read(py, |cs| noncovalent_bond_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintIter> {
        self.read(py, |cs| noncovalent_bond_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<NoncovalentBondConstraintItemsIter> {
        self.read(py, |cs| noncovalent_bond_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_ast();
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| NoncovalentBondConstraintAst::from_ast(py, constraint))
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
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<NoncovalentBondConstraintAst> {
        let ast_key = key.bind(py).borrow().to_ast();
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| NoncovalentBondConstraintAst::from_ast(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    fn __contains__(
        &self,
        py: Python<'_>,
        key: Py<NoncovalentBondConstraintKey>,
    ) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_ast();
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// Whether the bond is intramolecular; `Undetermined` when no `Intramolecular`
    /// constraint is present (mirroring the non-optional Rust accessor).
    #[getter]
    fn intramolecular(&self, py: Python<'_>) -> PyResult<BooleanAst> {
        self.read(py, |cs| Ok(BooleanAst::from_ast(&cs.intramolecular())))
    }

    #[setter]
    fn set_intramolecular(&self, py: Python<'_>, value: BooleanArg) {
        self.set_ast(
            py,
            AstNoncovalentBondConstraintAst::intramolecular(value.to_ast(py)),
        );
    }

    /// The present constraints as a dict keyed by snake_case name.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| noncovalent_bond_constraints_asdict(py, cs))
    }
}

/// A view of one noncovalent bond within a molecule: a handle to the molecule plus the
/// bond's index. Field reads rebuild the transient Rust view; the molecule is never
/// copied. The two endpoint atom indices are read-only topology; the kind and
/// constraints are the mutable bond value.
#[pyclass]
pub struct NoncovalentBondView {
    owner: Py<MoleculeAst>,
    id: AstNoncovalentBondId,
}

impl NoncovalentBondView {
    fn noncovalent_bond<'a>(
        &self,
        molecule: &'a AstMoleculeAst,
    ) -> PyResult<AstNoncovalentBondView<'a>> {
        molecule
            .noncovalent_bonds()
            .get(self.id)
            .ok_or_else(|| PyIndexError::new_err("noncovalent bond id out of range"))
    }
}

#[pymethods]
impl NoncovalentBondView {
    #[getter]
    fn id(&self) -> u32 {
        self.id.0
    }

    /// The two endpoint atom indices (read-only — participants are topology, not part of
    /// the bond value; the pair is unordered).
    #[getter]
    fn atom_ids(&self, py: Python<'_>) -> PyResult<(u32, u32)> {
        let molecule = self.owner.bind(py).borrow();
        let [first, second] = self.noncovalent_bond(molecule.inner())?.atom_ids();
        Ok((first.0, second.0))
    }

    fn __repr__(&self) -> String {
        format!("NoncovalentBondView(id={})", self.id.0)
    }

    /// The interaction kind.
    #[getter]
    fn kind(&self, py: Python<'_>) -> PyResult<NoncovalentBondKindAst> {
        let molecule = self.owner.bind(py).borrow();
        Ok(NoncovalentBondKindAst::from_ast(
            &self.noncovalent_bond(molecule.inner())?.ast.kind,
        ))
    }

    #[setter]
    fn set_kind(&self, py: Python<'_>, value: NoncovalentBondKindArg) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .noncovalent_bond_mut(self.id)
            .ast
            .kind = value.to_ast(py);
    }

    /// The bond's constraints as a live handle onto the molecule: reads borrow the
    /// current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(&self, py: Python<'_>) -> NoncovalentBondConstraintsView {
        NoncovalentBondConstraintsView {
            backing: NoncovalentBondConstraintsBacking::Molecule {
                owner: self.owner.clone_ref(py),
                id: self.id,
            },
        }
    }

    /// Replace the whole constraint set of the backing bond in place (wipe-and-set)
    /// from a value container or a live view.
    #[setter]
    fn set_constraints(
        &self,
        py: Python<'_>,
        value: NoncovalentBondConstraintsArg,
    ) -> PyResult<()> {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .noncovalent_bond_mut(self.id)
            .ast
            .constraints = value.to_ast(py)?;
        Ok(())
    }

    /// The value fields as a dict keyed by field name; values are the field mirrors —
    /// symmetric with `NoncovalentBondAst.asdict`, read through the view.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let molecule = self.owner.bind(py).borrow();
        let bond = self.noncovalent_bond(molecule.inner())?.ast;
        let dict = PyDict::new(py);
        dict.set_item("kind", NoncovalentBondKindAst::from_ast(&bond.kind))?;
        dict.set_item(
            "constraints",
            noncovalent_bond_constraints_asdict(py, &bond.constraints)?,
        )?;
        Ok(dict)
    }
}

/// Resolve a possibly-negative Python index (negative counts from the end) into an
/// existing noncovalent bond id, or `IndexError`. `NoncovalentBondId` is `RelationId`-
/// backed but contiguous for fresh molecules, so integer positions address it directly.
fn resolve_noncovalent_bond_index(
    molecule: &AstMoleculeAst,
    index: isize,
) -> PyResult<AstNoncovalentBondId> {
    let count = molecule.noncovalent_bonds().count();
    let resolved = if index < 0 {
        index + count as isize
    } else {
        index
    };
    if resolved < 0 {
        return Err(PyIndexError::new_err("noncovalent bond id out of range"));
    }
    let id = AstNoncovalentBondId(resolved as u32);
    if molecule.noncovalent_bonds().contains(id) {
        Ok(id)
    } else {
        Err(PyIndexError::new_err("noncovalent bond id out of range"))
    }
}

/// The noncovalent bonds of a molecule, indexed by integer position.
#[pyclass]
pub struct NoncovalentBondViews {
    owner: Py<MoleculeAst>,
}

#[pymethods]
impl NoncovalentBondViews {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner
            .bind(py)
            .borrow()
            .inner()
            .noncovalent_bonds()
            .count()
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "NoncovalentBondViews(len={})",
            self.owner
                .bind(py)
                .borrow()
                .inner()
                .noncovalent_bonds()
                .count()
        )
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<NoncovalentBondView> {
        let molecule = self.owner.bind(py).borrow();
        let id = resolve_noncovalent_bond_index(molecule.inner(), index)?;
        Ok(NoncovalentBondView {
            owner: self.owner.clone_ref(py),
            id,
        })
    }

    /// Replace the whole noncovalent bond value at `index` in place (endpoints unchanged).
    fn __setitem__(
        &self,
        py: Python<'_>,
        index: isize,
        bond: PyRef<'_, NoncovalentBondAst>,
    ) -> PyResult<()> {
        let mut molecule = self.owner.borrow_mut(py);
        let id = resolve_noncovalent_bond_index(molecule.inner(), index)?;
        *molecule.inner_mut().noncovalent_bond_mut(id).ast = bond.inner().clone();
        Ok(())
    }

    /// The noncovalent bond between atoms `first` and `second`, or `None`.
    fn connecting(&self, py: Python<'_>, first: u32, second: u32) -> Option<NoncovalentBondView> {
        let molecule = self.owner.bind(py).borrow();
        molecule
            .inner()
            .noncovalent_bonds()
            .connecting_id(AstAtomId(first), AstAtomId(second))
            .map(|id| NoncovalentBondView {
                owner: self.owner.clone_ref(py),
                id,
            })
    }

    /// The noncovalent bonds `atom` is an endpoint of.
    fn incident(&self, py: Python<'_>, atom: u32) -> Vec<NoncovalentBondView> {
        let molecule = self.owner.bind(py).borrow();
        molecule
            .inner()
            .noncovalent_bonds()
            .incident_ids(AstAtomId(atom))
            .map(|id| NoncovalentBondView {
                owner: self.owner.clone_ref(py),
                id,
            })
            .collect()
    }

    fn __iter__(&self, py: Python<'_>) -> NoncovalentBondViewIter {
        let ids = self
            .owner
            .bind(py)
            .borrow()
            .inner()
            .noncovalent_bonds()
            .ids()
            .collect::<Vec<_>>();
        NoncovalentBondViewIter {
            owner: self.owner.clone_ref(py),
            ids: ids.into_iter(),
        }
    }
}

#[pyclass]
struct NoncovalentBondViewIter {
    owner: Py<MoleculeAst>,
    ids: IntoIter<AstNoncovalentBondId>,
}

#[pymethods]
impl NoncovalentBondViewIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> Option<NoncovalentBondView> {
        self.ids.next().map(|id| NoncovalentBondView {
            owner: self.owner.clone_ref(py),
            id,
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{AtomAst as AstAtomAst, BooleanAst as AstBooleanAst, MoleculeParts};
    use umol_chem::element::Element as ChemElement;

    use super::*;

    #[rstest]
    #[case(AstNoncovalentBondKindAst::Undetermined)]
    #[case(AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond))]
    #[case(AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::VanDerWaals))]
    fn test_noncovalent_bond_kind_ast_roundtrip(#[case] ast: AstNoncovalentBondKindAst) {
        assert_eq!(NoncovalentBondKindAst::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstNoncovalentBondKindAst::Undetermined, None)]
    #[case(
        AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::Ionic),
        Some(NoncovalentBondKind::Ionic)
    )]
    fn test_noncovalent_bond_kind_ast_as_lit(
        #[case] ast: AstNoncovalentBondKindAst,
        #[case] expected: Option<NoncovalentBondKind>,
    ) {
        assert_eq!(NoncovalentBondKindAst::from_ast(&ast).as_lit(), expected);
    }

    #[rstest]
    #[case(AstNoncovalentBondKind::HydrogenBond)]
    #[case(AstNoncovalentBondKind::ChalcogenBond)]
    fn test_noncovalent_bond_kind_roundtrip(#[case] ast: AstNoncovalentBondKind) {
        assert_eq!(NoncovalentBondKind::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    fn test_noncovalent_bond_kind_arg_to_ast() {
        Python::attach(|py| {
            // a bare kind coerces to Lit
            assert_eq!(
                NoncovalentBondKindArg::Kind(NoncovalentBondKind::HydrogenBond).to_ast(py),
                AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond)
            );
            // a NoncovalentBondKindAst passes through
            let ast = Py::new(py, NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic)).unwrap();
            assert_eq!(
                NoncovalentBondKindArg::Ast(ast).to_ast(py),
                AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::Ionic)
            );
        });
    }

    /// A `Py<NoncovalentBondConstraintAst>` for `Intramolecular(b)`.
    fn intramolecular(py: Python<'_>, b: bool) -> Py<NoncovalentBondConstraintAst> {
        into_py_variant(
            py,
            NoncovalentBondConstraintAst::from_ast(
                py,
                &AstNoncovalentBondConstraintAst::intramolecular(b),
            )
            .unwrap(),
        )
        .unwrap()
    }

    fn intramolecular_key(py: Python<'_>) -> Py<NoncovalentBondConstraintKey> {
        into_py_variant(py, NoncovalentBondConstraintKey::Intramolecular()).unwrap()
    }

    #[rstest]
    fn test_noncovalent_bond_constraint_key_roundtrip() {
        let key = NoncovalentBondConstraintKey::from_ast(
            &AstNoncovalentBondConstraintKey::Intramolecular,
        );
        assert_eq!(
            key.to_ast(),
            AstNoncovalentBondConstraintKey::Intramolecular
        );
    }

    #[rstest]
    fn test_noncovalent_bond_constraint_ast_key() {
        Python::attach(|py| {
            let constraint = AstNoncovalentBondConstraintAst::intramolecular(true);
            let key = NoncovalentBondConstraintAst::from_ast(py, &constraint)
                .unwrap()
                .key(py);
            assert_eq!(
                key.to_ast(),
                AstNoncovalentBondConstraintKey::Intramolecular
            );
        });
    }

    #[rstest]
    #[case(AstNoncovalentBondConstraintAst::intramolecular(true))]
    #[case(AstNoncovalentBondConstraintAst::intramolecular(false))]
    #[case(AstNoncovalentBondConstraintAst::Intramolecular(AstBooleanAst::Undetermined))]
    fn test_noncovalent_bond_constraint_ast_roundtrip(
        #[case] ast: AstNoncovalentBondConstraintAst,
    ) {
        Python::attach(|py| {
            assert_eq!(
                NoncovalentBondConstraintAst::from_ast(py, &ast)
                    .unwrap()
                    .to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_new() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.intramolecular().to_ast(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_repr() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]);
            assert_eq!(
                constraints.__repr__(py).unwrap(),
                "NoncovalentBondConstraintsAst([NoncovalentBondConstraintAst.Intramolecular(BooleanAst.Lit(True))])"
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = NoncovalentBondConstraintsAst::new(py, vec![]);
            constraints.set(py, intramolecular(py, false));
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.intramolecular().to_ast(),
                AstBooleanAst::Lit(false)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_pop() {
        Python::attach(|py| {
            let mut constraints =
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]);
            let removed = constraints.pop(py, intramolecular_key(py)).unwrap();
            match removed {
                Some(NoncovalentBondConstraintAst::Intramolecular(b)) => {
                    assert_eq!(b.bind(py).borrow().to_ast(), AstBooleanAst::Lit(true))
                }
                _ => panic!("expected removed Intramolecular(Lit(true))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_update_container() {
        Python::attach(|py| {
            let constraints = Py::new(py, NoncovalentBondConstraintsAst::new(py, vec![])).unwrap();
            let mut other = AstNoncovalentBondConstraintsAst::new();
            other.set(AstNoncovalentBondConstraintAst::intramolecular(true));
            NoncovalentBondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                NoncovalentBondConstraintsUpdate::Container(
                    Py::new(py, NoncovalentBondConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            assert_eq!(
                constraints.bind(py).borrow().intramolecular().to_ast(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_update_entries() {
        Python::attach(|py| {
            let constraints = Py::new(py, NoncovalentBondConstraintsAst::new(py, vec![])).unwrap();
            NoncovalentBondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                NoncovalentBondConstraintsUpdate::Entries(vec![intramolecular(py, false)]),
            )
            .unwrap();
            assert_eq!(
                constraints.bind(py).borrow().intramolecular().to_ast(),
                AstBooleanAst::Lit(false)
            );
        });
    }

    /// Regression: a container updating itself resolves `other` before the write borrow,
    /// so it is an idempotent no-op, not a RefCell double-borrow panic.
    #[rstest]
    fn test_noncovalent_bond_constraints_ast_update_self() {
        Python::attach(|py| {
            let constraints = Py::new(
                py,
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]),
            )
            .unwrap();
            NoncovalentBondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                NoncovalentBondConstraintsUpdate::Container(constraints.clone_ref(py)),
            )
            .unwrap();
            assert_eq!(
                constraints.bind(py).borrow().intramolecular().to_ast(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_len_contains() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]);
            assert_eq!(constraints.__len__(), 1);
            assert!(constraints.__contains__(py, intramolecular_key(py)));
            let empty = NoncovalentBondConstraintsAst::new(py, vec![]);
            assert!(!empty.__contains__(py, intramolecular_key(py)));
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_keys_values_items() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]);
            let mut keys = constraints.keys(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(),
                AstNoncovalentBondConstraintKey::Intramolecular
            );
            assert!(keys.__next__().is_none());
            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstNoncovalentBondConstraintAst::intramolecular(true)
            );
            let mut items = constraints.items(py).unwrap();
            let (k, v) = items.__next__().unwrap();
            assert_eq!(
                k.bind(py).borrow().to_ast(),
                AstNoncovalentBondConstraintKey::Intramolecular
            );
            assert_eq!(
                v.bind(py).borrow().to_ast(py),
                AstNoncovalentBondConstraintAst::intramolecular(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_get() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]);
            let present = constraints.get(py, intramolecular_key(py), None).unwrap();
            let expected = intramolecular(py, true).into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());
            // absent → None
            let empty = NoncovalentBondConstraintsAst::new(py, vec![]);
            let absent = empty.get(py, intramolecular_key(py), None).unwrap();
            assert!(absent.bind(py).is_none());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_intramolecular() {
        Python::attach(|py| {
            let present = NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]);
            assert_eq!(present.intramolecular().to_ast(), AstBooleanAst::Lit(true));
            // absent → Undetermined (non-optional accessor)
            let empty = NoncovalentBondConstraintsAst::new(py, vec![]);
            assert_eq!(empty.intramolecular().to_ast(), AstBooleanAst::Undetermined);
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_set_intramolecular() {
        Python::attach(|py| {
            let mut constraints = NoncovalentBondConstraintsAst::new(py, vec![]);
            constraints.set_intramolecular(py, BooleanArg::Lit(false));
            assert_eq!(
                constraints.intramolecular().to_ast(),
                AstBooleanAst::Lit(false)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = NoncovalentBondConstraintsAst::new(py, vec![]);
            assert!(constraints.__getitem__(py, intramolecular_key(py)).is_err());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = NoncovalentBondConstraintsAst::new(py, vec![]);
            assert!(constraints.__delitem__(py, intramolecular_key(py)).is_err());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_asdict() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]);
            let dict = constraints.asdict(py).unwrap();
            let value = dict.get_item("intramolecular").unwrap().unwrap();
            let expected = into_py_variant(py, BooleanAst::Lit(true)).unwrap();
            assert!(value.eq(expected.bind(py)).unwrap());
        });
    }

    /// A standalone `NoncovalentBondAst` value pyclass (hydrogen bond, no constraints).
    fn hbond(py: Python<'_>) -> Py<NoncovalentBondAst> {
        Py::new(
            py,
            NoncovalentBondAst::from_inner(AstNoncovalentBondAst::from_kind(
                AstNoncovalentBondKind::HydrogenBond,
            )),
        )
        .unwrap()
    }

    #[rstest]
    fn test_noncovalent_bond_ast_new() {
        Python::attach(|py| {
            let bond = NoncovalentBondAst::new(
                py,
                NoncovalentBondKindArg::Kind(NoncovalentBondKind::HydrogenBond),
                None,
            );
            assert_eq!(
                bond.inner().kind,
                AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond)
            );
            assert_eq!(bond.inner().constraints.len(), 0);
        });
    }

    #[rstest]
    fn test_noncovalent_bond_ast_new_constraints() {
        Python::attach(|py| {
            let constraints = Py::new(
                py,
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]),
            )
            .unwrap();
            let bond = NoncovalentBondAst::new(
                py,
                NoncovalentBondKindArg::Kind(NoncovalentBondKind::HalogenBond),
                Some(constraints),
            );
            assert_eq!(
                bond.inner().constraints.intramolecular(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    #[case("Hbd")]
    #[case("Hbd#I")]
    #[case("Hbd#I!")]
    #[case("*")]
    fn test_noncovalent_bond_ast_parse(#[case] dsl: &str) {
        let bond = NoncovalentBondAst::parse(dsl).unwrap();
        assert_eq!(bond.__str__(), dsl);
        assert_eq!(
            bond.__repr__(),
            format!("NoncovalentBondAst.parse('{dsl}')")
        );
    }

    #[rstest]
    fn test_noncovalent_bond_ast_parse_error() {
        assert!(NoncovalentBondAst::parse("z").is_err());
    }

    #[rstest]
    fn test_noncovalent_bond_ast_kind() {
        Python::attach(|py| {
            let mut bond = NoncovalentBondAst::from_inner(AstNoncovalentBondAst::from_kind(
                AstNoncovalentBondKind::HydrogenBond,
            ));
            assert_eq!(
                bond.kind().to_ast(),
                AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond)
            );
            bond.set_kind(py, NoncovalentBondKindArg::Kind(NoncovalentBondKind::Ionic));
            assert_eq!(
                bond.kind().to_ast(),
                AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::Ionic)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_ast_set_constraints() {
        Python::attach(|py| {
            let bond = hbond(py);
            let constraints = Py::new(
                py,
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, true)]),
            )
            .unwrap();
            NoncovalentBondAst::set_constraints(
                bond.clone_ref(py),
                py,
                NoncovalentBondConstraintsArg::Container(constraints),
            )
            .unwrap();
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.intramolecular(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_ast_set_constraints_from_view() {
        Python::attach(|py| {
            // source bond carrying a constraint, exposed as a live view
            let source = Py::new(
                py,
                NoncovalentBondAst::from_inner(
                    AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond)
                        .with_constraint(AstNoncovalentBondConstraintAst::intramolecular(true)),
                ),
            )
            .unwrap();
            let view = NoncovalentBondAst::constraints(source);
            let dest = hbond(py);
            NoncovalentBondAst::set_constraints(
                dest.clone_ref(py),
                py,
                NoncovalentBondConstraintsArg::View(Py::new(py, view).unwrap()),
            )
            .unwrap();
            assert_eq!(
                dest.bind(py).borrow().inner().constraints.intramolecular(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    /// Regression: assigning a bond's own constraints view back to it snapshots before
    /// the write borrow, so it is a no-op, not a RefCell double-borrow panic
    /// (`bond.constraints = bond.constraints`).
    #[rstest]
    fn test_noncovalent_bond_ast_set_constraints_self() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                NoncovalentBondAst::from_inner(
                    AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond)
                        .with_constraint(AstNoncovalentBondConstraintAst::intramolecular(true)),
                ),
            )
            .unwrap();
            let own_view = NoncovalentBondAst::constraints(bond.clone_ref(py));
            NoncovalentBondAst::set_constraints(
                bond.clone_ref(py),
                py,
                NoncovalentBondConstraintsArg::View(Py::new(py, own_view).unwrap()),
            )
            .unwrap();
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.intramolecular(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_ast_constraints_write_through() {
        Python::attach(|py| {
            let bond = hbond(py);
            let view = NoncovalentBondAst::constraints(bond.clone_ref(py));
            view.set(py, intramolecular(py, true));
            // the write hit the standalone bond, not a copy
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.intramolecular(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_ast_asdict() {
        Python::attach(|py| {
            let bond = NoncovalentBondAst::from_inner(
                AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond)
                    .with_constraint(AstNoncovalentBondConstraintAst::intramolecular(true)),
            );
            let dict = bond.asdict(py).unwrap();
            let kind = dict.get_item("kind").unwrap().unwrap();
            let expected_kind = into_py_variant(
                py,
                NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
            )
            .unwrap();
            assert!(kind.eq(expected_kind.bind(py)).unwrap());
            let constraints = dict.get_item("constraints").unwrap().unwrap();
            assert_eq!(constraints.len().unwrap(), 1);
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_view_set() {
        Python::attach(|py| {
            let bond = hbond(py);
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            view.set(py, intramolecular(py, true));
            // a fresh view proves the write hit the standalone bond, not a copy
            let fresh = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            assert_eq!(
                fresh.intramolecular(py).unwrap().to_ast(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_view_pop() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                NoncovalentBondAst::from_inner(
                    AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond)
                        .with_constraint(AstNoncovalentBondConstraintAst::intramolecular(true)),
                ),
            )
            .unwrap();
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            let removed = view.pop(py, intramolecular_key(py)).unwrap();
            match removed {
                Some(NoncovalentBondConstraintAst::Intramolecular(b)) => {
                    assert_eq!(b.bind(py).borrow().to_ast(), AstBooleanAst::Lit(true))
                }
                _ => panic!("expected removed Intramolecular(Lit(true))"),
            }
            let fresh = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_view_update() {
        Python::attach(|py| {
            let bond = hbond(py);
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            view.update(
                py,
                NoncovalentBondConstraintsUpdate::Entries(vec![intramolecular(py, false)]),
            )
            .unwrap();
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.intramolecular(),
                AstBooleanAst::Lit(false)
            );
        });
    }

    /// Regression: a view updating from a view over the same bond resolves `other`
    /// before the write borrow, so it is an idempotent no-op, not a double-borrow panic
    /// (`bond.constraints.update(bond.constraints)`).
    #[rstest]
    fn test_noncovalent_bond_constraints_view_update_self() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                NoncovalentBondAst::from_inner(
                    AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond)
                        .with_constraint(AstNoncovalentBondConstraintAst::intramolecular(true)),
                ),
            )
            .unwrap();
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            let other = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            view.update(
                py,
                NoncovalentBondConstraintsUpdate::View(Py::new(py, other).unwrap()),
            )
            .unwrap();
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.intramolecular(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_view_set_intramolecular() {
        Python::attach(|py| {
            let bond = hbond(py);
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            view.set_intramolecular(py, BooleanArg::Lit(true));
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.intramolecular(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    /// A molecule of two oxygens with one hydrogen bond over atoms (0, 1), noncovalent id 0.
    fn molecule_with_hbond(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = AstMoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::O),
                AstAtomAst::from_element(ChemElement::O),
            ],
            noncovalent: vec![(
                AstAtomId(0),
                AstAtomId(1),
                AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_inner(molecule)).unwrap()
    }

    #[rstest]
    fn test_noncovalent_bond_view_id_atom_ids() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let view = NoncovalentBondView {
                owner,
                id: AstNoncovalentBondId(0),
            };
            assert_eq!(view.id(), 0);
            assert_eq!(view.atom_ids(py).unwrap(), (0, 1));
            assert_eq!(view.__repr__(), "NoncovalentBondView(id=0)");
        });
    }

    #[rstest]
    fn test_noncovalent_bond_view_kind() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let view = NoncovalentBondView {
                owner: owner.clone_ref(py),
                id: AstNoncovalentBondId(0),
            };
            assert_eq!(
                view.kind(py).unwrap().to_ast(),
                AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond)
            );
            view.set_kind(py, NoncovalentBondKindArg::Kind(NoncovalentBondKind::Ionic));
            // a fresh read proves the write hit the molecule
            let fresh = NoncovalentBondView {
                owner,
                id: AstNoncovalentBondId(0),
            };
            assert_eq!(
                fresh.kind(py).unwrap().to_ast(),
                AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::Ionic)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_view_constraints_write_through() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let view = NoncovalentBondView {
                owner: owner.clone_ref(py),
                id: AstNoncovalentBondId(0),
            };
            // the constraints handle is molecule-backed; a write goes through to the bond
            view.constraints(py)
                .set_intramolecular(py, BooleanArg::Lit(true));
            assert_eq!(
                owner
                    .bind(py)
                    .borrow()
                    .inner()
                    .noncovalent_bond(AstNoncovalentBondId(0))
                    .ast
                    .constraints
                    .intramolecular(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_view_set_constraints() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let view = NoncovalentBondView {
                owner: owner.clone_ref(py),
                id: AstNoncovalentBondId(0),
            };
            let constraints = Py::new(
                py,
                NoncovalentBondConstraintsAst::new(py, vec![intramolecular(py, false)]),
            )
            .unwrap();
            view.set_constraints(py, NoncovalentBondConstraintsArg::Container(constraints))
                .unwrap();
            assert_eq!(
                owner
                    .bind(py)
                    .borrow()
                    .inner()
                    .noncovalent_bond(AstNoncovalentBondId(0))
                    .ast
                    .constraints
                    .intramolecular(),
                AstBooleanAst::Lit(false)
            );
        });
    }

    /// Regression: `mol.noncovalent_bonds[i].constraints.update(same-view)` resolves before
    /// the molecule write borrow, so the molecule-backed self-alias is a no-op, not a panic.
    #[rstest]
    fn test_noncovalent_bond_view_constraints_update_self() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            owner
                .bind(py)
                .borrow_mut()
                .inner_mut()
                .noncovalent_bond_mut(AstNoncovalentBondId(0))
                .ast
                .constraints
                .set(AstNoncovalentBondConstraintAst::intramolecular(true));
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstNoncovalentBondId(0),
                },
            };
            let other = Py::new(
                py,
                NoncovalentBondConstraintsView {
                    backing: NoncovalentBondConstraintsBacking::Molecule {
                        owner: owner.clone_ref(py),
                        id: AstNoncovalentBondId(0),
                    },
                },
            )
            .unwrap();
            view.update(py, NoncovalentBondConstraintsUpdate::View(other))
                .unwrap();
            assert_eq!(
                owner
                    .bind(py)
                    .borrow()
                    .inner()
                    .noncovalent_bond(AstNoncovalentBondId(0))
                    .ast
                    .constraints
                    .intramolecular(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_view_asdict() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let view = NoncovalentBondView {
                owner,
                id: AstNoncovalentBondId(0),
            };
            let dict = view.asdict(py).unwrap();
            assert_eq!(
                dict.keys()
                    .iter()
                    .map(|k| k.extract::<String>().unwrap())
                    .collect::<Vec<_>>(),
                vec!["kind".to_string(), "constraints".to_string()]
            );
            let kind = dict.get_item("kind").unwrap().unwrap();
            let expected = into_py_variant(
                py,
                NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
            )
            .unwrap();
            assert!(kind.eq(expected.bind(py)).unwrap());
        });
    }

    /// Three atoms, one hydrogen bond over (0, 1), atom 2 isolated. For the collection
    /// negative cases (`connecting` / `incident` with no bond).
    fn molecule_with_hbond_and_isolated(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = AstMoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::O),
                AstAtomAst::from_element(ChemElement::O),
                AstAtomAst::from_element(ChemElement::O),
            ],
            noncovalent: vec![(
                AstAtomId(0),
                AstAtomId(1),
                AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_inner(molecule)).unwrap()
    }

    #[rstest]
    fn test_noncovalent_bond_views_len_getitem() {
        Python::attach(|py| {
            let views = NoncovalentBondViews {
                owner: molecule_with_hbond(py),
            };
            assert_eq!(views.__len__(py), 1);
            assert_eq!(views.__repr__(py), "NoncovalentBondViews(len=1)");
            assert_eq!(views.__getitem__(py, 0).unwrap().id(), 0);
            // negative index counts from the end
            assert_eq!(views.__getitem__(py, -1).unwrap().id(), 0);
            assert!(views.__getitem__(py, 1).is_err());
            assert!(views.__getitem__(py, -2).is_err());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_views_setitem() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let views = NoncovalentBondViews {
                owner: owner.clone_ref(py),
            };
            let replacement = Py::new(
                py,
                NoncovalentBondAst::from_inner(AstNoncovalentBondAst::from_kind(
                    AstNoncovalentBondKind::Ionic,
                )),
            )
            .unwrap();
            views
                .__setitem__(py, 0, replacement.bind(py).borrow())
                .unwrap();
            let view = views.__getitem__(py, 0).unwrap();
            // value replaced, endpoints preserved
            assert_eq!(
                view.kind(py).unwrap().to_ast(),
                AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::Ionic)
            );
            assert_eq!(view.atom_ids(py).unwrap(), (0, 1));
        });
    }

    #[rstest]
    fn test_noncovalent_bond_views_setitem_out_of_range() {
        Python::attach(|py| {
            let views = NoncovalentBondViews {
                owner: molecule_with_hbond(py),
            };
            let bond = Py::new(
                py,
                NoncovalentBondAst::from_inner(AstNoncovalentBondAst::from_kind(
                    AstNoncovalentBondKind::Ionic,
                )),
            )
            .unwrap();
            assert!(views.__setitem__(py, 5, bond.bind(py).borrow()).is_err());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_views_iter() {
        Python::attach(|py| {
            let views = NoncovalentBondViews {
                owner: molecule_with_hbond(py),
            };
            let mut iter = views.__iter__(py);
            assert_eq!(iter.__next__(py).unwrap().id(), 0);
            assert!(iter.__next__(py).is_none());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_views_connecting() {
        Python::attach(|py| {
            let views = NoncovalentBondViews {
                owner: molecule_with_hbond_and_isolated(py),
            };
            // unordered pair — both orders find the same bond
            assert_eq!(views.connecting(py, 0, 1).unwrap().id(), 0);
            assert_eq!(views.connecting(py, 1, 0).unwrap().id(), 0);
            // no bond between 0 and the isolated atom 2
            assert!(views.connecting(py, 0, 2).is_none());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_views_incident() {
        Python::attach(|py| {
            let views = NoncovalentBondViews {
                owner: molecule_with_hbond_and_isolated(py),
            };
            assert_eq!(
                views
                    .incident(py, 0)
                    .iter()
                    .map(|v| v.id())
                    .collect::<Vec<_>>(),
                vec![0]
            );
            assert!(views.incident(py, 2).is_empty());
        });
    }
}
