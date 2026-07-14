//! Bond constraint mirrors, containers, and live views.

use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use umol_ast::ast::{
    BondConstraintAst as AstBondConstraintAst, BondConstraintKey as AstBondConstraintKey,
    BondConstraintsAst as AstBondConstraintsAst, BondId as AstBondId, RingScope as AstRingScope,
};

use super::ring::{RingMembershipAst, RingScope};
use crate::bond::BondAst;
use crate::boolean::{BooleanArg, BooleanAst};
use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::molecule::MoleculeAst;
use crate::stereo::{CisTransStereoArg, CisTransStereoAst};
use crate::value::{ValueArg, ValueAst};

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
    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
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
    pub(crate) fn key(&self, py: Python<'_>) -> PyResult<BondConstraintKey> {
        BondConstraintKey::from_ast(py, &self.to_ast(py).key())
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
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
            AstBondConstraintAst::CisTransStereo(c) => {
                Self::CisTransStereo(into_py_variant(py, CisTransStereoAst::from_ast(py, c)?)?)
            }
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
pub(crate) enum BondConstraintsUpdate {
    Container(Py<BondConstraintsAst>),
    View(Py<BondConstraintsView>),
    Entries(Vec<Py<BondConstraintAst>>),
}

impl BondConstraintsUpdate {
    /// Read every Python object into owned data — no write target is touched. Callers
    /// resolve *before* taking the write borrow so a view (or container) that aliases
    /// the same bond is read while nothing is borrowed (otherwise
    /// `bond.constraints.update(bond.constraints)` self-aliases into a double-borrow panic).
    pub(crate) fn resolve(&self, py: Python<'_>) -> PyResult<ResolvedBondConstraintsUpdate> {
        Ok(match self {
            BondConstraintsUpdate::Container(c) => {
                ResolvedBondConstraintsUpdate::Overlay(c.bind(py).borrow().inner().clone())
            }
            BondConstraintsUpdate::View(v) => ResolvedBondConstraintsUpdate::Overlay(
                v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?,
            ),
            BondConstraintsUpdate::Entries(entries) => ResolvedBondConstraintsUpdate::Entries(
                entries
                    .iter()
                    .map(|entry| entry.bind(py).borrow().to_ast(py))
                    .collect(),
            ),
        })
    }
}

/// A `BondConstraintsUpdate` with all Python-object reads already done, so it can be
/// applied under a write borrow without re-entering Python.
pub(crate) enum ResolvedBondConstraintsUpdate {
    /// A whole container (from another container or a live view): overlaid via `update`
    /// (last-wins per key; undetermined entries remove).
    Overlay(AstBondConstraintsAst),
    /// Loose entries: `set` each (last-wins; undetermined entries stored, not removed).
    Entries(Vec<AstBondConstraintAst>),
}

impl ResolvedBondConstraintsUpdate {
    /// Overlay onto `target` in place. No Python reads.
    pub(crate) fn apply(self, target: &mut AstBondConstraintsAst) {
        match self {
            ResolvedBondConstraintsUpdate::Overlay(overlay) => target.update(&overlay),
            ResolvedBondConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry);
                }
            }
        }
    }
}

/// A whole-container argument that snapshots either a value container or a live
/// view — for the bond `constraints` setter, which accepts either.
#[derive(FromPyObject)]
pub(crate) enum BondConstraintsArg {
    Container(Py<BondConstraintsAst>),
    View(Py<BondConstraintsView>),
}

impl BondConstraintsArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> PyResult<AstBondConstraintsAst> {
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
    pub(crate) fn new(py: Python<'_>, entries: Vec<Py<BondConstraintAst>>) -> Self {
        let mut constraints = AstBondConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_ast(py)),
        );
        BondConstraintsAst(constraints)
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let mirror = into_py_variant(py, BondConstraintAst::from_ast(py, entry)?)?;
            parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("BondConstraintsAst([{}])", parts.join(", ")))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    pub(crate) fn set(&mut self, py: Python<'_>, c: Py<BondConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    pub(crate) fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
    ) -> PyResult<Option<BondConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_ast(py))
            .map(|c| BondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container, a live view, or an
    /// iterable of `BondConstraintAst` (last-wins per key; undetermined entries remove).
    /// Takes `slf` by handle so `other` is fully read *before* the write borrow —
    /// `cs.update(cs)` on the same container is then a no-op, not a double-borrow panic.
    pub(crate) fn update(
        slf: Py<Self>,
        py: Python<'_>,
        other: BondConstraintsUpdate,
    ) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        resolved.apply(&mut slf.borrow_mut(py).0);
        Ok(())
    }

    pub(crate) fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<BondConstraintKeyIter> {
        bond_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<BondConstraintKeyIter> {
        bond_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<BondConstraintIter> {
        bond_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<BondConstraintItemsIter> {
        bond_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
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
    pub(crate) fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
    ) -> PyResult<BondConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_ast(py)) {
            Some(constraint) => BondConstraintAst::from_ast(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
    ) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_ast(py)).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    pub(crate) fn __contains__(&self, py: Python<'_>, key: Py<BondConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast(py))
    }

    /// The aromatic value; `Undetermined` when no `Aromatic` constraint is present
    /// (mirroring the non-optional Rust accessor).
    #[getter]
    pub(crate) fn aromatic(&self) -> BooleanAst {
        BooleanAst::from_ast(&self.0.aromatic())
    }

    #[setter]
    pub(crate) fn set_aromatic(&mut self, py: Python<'_>, value: BooleanArg) {
        self.0.set(AstBondConstraintAst::aromatic(value.to_ast(py)));
    }

    /// The cis/trans-stereo state, or `None`.
    #[getter]
    pub(crate) fn cis_trans_stereo(&self, py: Python<'_>) -> PyResult<Option<CisTransStereoAst>> {
        self.0
            .cis_trans_stereo()
            .map(|c| CisTransStereoAst::from_ast(py, c))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_cis_trans_stereo(
        &mut self,
        py: Python<'_>,
        value: CisTransStereoArg,
    ) -> PyResult<()> {
        self.0
            .set(AstBondConstraintAst::cis_trans_stereo(value.to_ast(py)?));
        Ok(())
    }

    /// The all-rings membership count, or `None`.
    #[getter]
    pub(crate) fn ring_count(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_count()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    pub(crate) fn set_ring_count(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstBondConstraintAst::ring_membership(
            AstRingScope::All,
            value.to_ast(py),
        ));
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    pub(crate) fn ring_size_count(slf: Py<Self>) -> BondRingSizeCounts {
        BondRingSizeCounts {
            backing: BondRingSizeBacking::Value(slf),
        }
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors. Ring memberships key by scope: `ring_count` for the
    /// all-rings scope, `ring_size_count_<n>` for a specific ring size.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
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
pub(crate) fn bond_constraints_iter(
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
pub(crate) fn bond_constraint_keys(
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
pub(crate) fn bond_constraint_items(
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
pub(crate) fn bond_constraints_asdict<'py>(
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

/// What a `BondConstraintsView` writes through to: a bond within a molecule (by
/// index) or a standalone `BondAst`.
pub(crate) enum BondConstraintsBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: AstBondId,
    },
    Bond(Py<BondAst>),
}

/// A live handle onto one bond's constraints, backed by either a molecule-bond or a
/// standalone `BondAst`. Reads borrow the bond's constraints and read only the item
/// they need (no whole-container clone); mutators write through to the bond in place,
/// without a clone-and-writeback.
#[pyclass]
pub struct BondConstraintsView {
    pub(crate) backing: BondConstraintsBacking,
}

impl BondConstraintsView {
    /// Borrow the backing bond's constraints and read one item through `f` — no clone.
    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            BondConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .bonds()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("bond id out of range"))?;
                f(&view.ast.constraints)
            }
            BondConstraintsBacking::Bond(bond) => {
                let bond = bond.bind(py).borrow();
                f(&bond.inner().constraints)
            }
        }
    }

    /// Mutate the backing bond's constraints in place through `f`.
    pub(crate) fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut AstBondConstraintsAst) -> R,
    ) -> R {
        match &self.backing {
            BondConstraintsBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .bond_mut(*id)
                .ast
                .constraints),
            BondConstraintsBacking::Bond(bond) => {
                f(&mut bond.borrow_mut(py).inner_mut().constraints)
            }
        }
    }

    /// Set one constraint on the backing bond in place (last-wins per key).
    pub(crate) fn set_ast(&self, py: Python<'_>, constraint: AstBondConstraintAst) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing bond in place, returning the removed entry.
    pub(crate) fn remove_ast(
        &self,
        py: Python<'_>,
        key: AstBondConstraintKey,
    ) -> Option<AstBondConstraintAst> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl BondConstraintsView {
    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("BondConstraintsView({count} entries)"))
    }

    /// Insert `c` on the bond in place, replacing any existing entry of the same
    /// key (last-wins).
    pub(crate) fn set(&self, py: Python<'_>, c: Py<BondConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key from the bond in place, returning it if
    /// present (dict `pop`).
    pub(crate) fn pop(
        &self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
    ) -> PyResult<Option<BondConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_ast(py))
            .map(|c| BondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    pub(crate) fn __delitem__(&self, py: Python<'_>, key: Py<BondConstraintKey>) -> PyResult<()> {
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

    /// Overlay `other` onto the bond's constraints in place — another container, a live
    /// view, or an iterable of `BondConstraintAst` (last-wins per key; undetermined
    /// entries remove). Resolves `other` to owned data *before* the write borrow, so a
    /// view aliasing the same bond is not a double-borrow panic.
    pub(crate) fn update(&self, py: Python<'_>, other: BondConstraintsUpdate) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        self.with_mut(py, |cs| resolved.apply(cs));
        Ok(())
    }

    pub(crate) fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<BondConstraintKeyIter> {
        self.read(py, |cs| bond_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<BondConstraintKeyIter> {
        self.read(py, |cs| bond_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    pub(crate) fn values(&self, py: Python<'_>) -> PyResult<BondConstraintIter> {
        self.read(py, |cs| bond_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    pub(crate) fn items(&self, py: Python<'_>) -> PyResult<BondConstraintItemsIter> {
        self.read(py, |cs| bond_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    pub(crate) fn get(
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
    pub(crate) fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
    ) -> PyResult<BondConstraintAst> {
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

    pub(crate) fn __contains__(
        &self,
        py: Python<'_>,
        key: Py<BondConstraintKey>,
    ) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_ast(py);
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// The aromatic value; `Undetermined` when no `Aromatic` constraint is present
    /// (mirroring the non-optional Rust accessor).
    #[getter]
    pub(crate) fn aromatic(&self, py: Python<'_>) -> PyResult<BooleanAst> {
        self.read(py, |cs| Ok(BooleanAst::from_ast(&cs.aromatic())))
    }

    #[setter]
    pub(crate) fn set_aromatic(&self, py: Python<'_>, value: BooleanArg) {
        self.set_ast(py, AstBondConstraintAst::aromatic(value.to_ast(py)));
    }

    /// The cis/trans-stereo state, or `None`.
    #[getter]
    pub(crate) fn cis_trans_stereo(&self, py: Python<'_>) -> PyResult<Option<CisTransStereoAst>> {
        self.read(py, |cs| {
            cs.cis_trans_stereo()
                .map(|c| CisTransStereoAst::from_ast(py, c))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_cis_trans_stereo(
        &self,
        py: Python<'_>,
        value: CisTransStereoArg,
    ) -> PyResult<()> {
        self.set_ast(
            py,
            AstBondConstraintAst::cis_trans_stereo(value.to_ast(py)?),
        );
        Ok(())
    }

    /// The all-rings membership count, or `None`.
    #[getter]
    pub(crate) fn ring_count(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_count()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    pub(crate) fn set_ring_count(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(
            py,
            AstBondConstraintAst::ring_membership(AstRingScope::All, value.to_ast(py)),
        );
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    pub(crate) fn ring_size_count(&self, py: Python<'_>) -> BondRingSizeCounts {
        let backing = match &self.backing {
            BondConstraintsBacking::Molecule { owner, id } => BondRingSizeBacking::Molecule {
                owner: owner.clone_ref(py),
                id: *id,
            },
            BondConstraintsBacking::Bond(bond) => BondRingSizeBacking::Bond(bond.clone_ref(py)),
        };
        BondRingSizeCounts { backing }
    }

    /// The present constraints as a dict keyed by snake_case name.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| bond_constraints_asdict(py, cs))
    }
}

/// What a `BondRingSizeCounts` proxy reads/writes through to: a bond within a
/// molecule, a standalone `BondAst`, or a standalone `BondConstraintsAst` value.
pub(crate) enum BondRingSizeBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: AstBondId,
    },
    Bond(Py<BondAst>),
    Value(Py<BondConstraintsAst>),
}

/// A subscriptable proxy over the sized-ring membership counts of a bond, keyed by
/// ring size: `proxy[size]` reads, `proxy[size] = count` sets, `del proxy[size]`
/// removes. Backs onto whichever container produced it (dual-backing, like
/// `BondConstraintsView`).
#[pyclass]
pub struct BondRingSizeCounts {
    pub(crate) backing: BondRingSizeBacking,
}

impl BondRingSizeCounts {
    /// Borrow the backing constraints and read through `f` — no clone.
    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            BondRingSizeBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .bonds()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("bond id out of range"))?;
                f(&view.ast.constraints)
            }
            BondRingSizeBacking::Bond(bond) => f(&bond.bind(py).borrow().inner().constraints),
            BondRingSizeBacking::Value(value) => f(value.bind(py).borrow().inner()),
        }
    }

    /// Mutate the backing constraints in place through `f`.
    pub(crate) fn write(&self, py: Python<'_>, f: impl FnOnce(&mut AstBondConstraintsAst)) {
        match &self.backing {
            BondRingSizeBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .bond_mut(*id)
                .ast
                .constraints),
            BondRingSizeBacking::Bond(bond) => f(&mut bond.borrow_mut(py).inner_mut().constraints),
            BondRingSizeBacking::Value(value) => f(value.borrow_mut(py).inner_mut()),
        }
    }
}

#[pymethods]
impl BondRingSizeCounts {
    /// The membership count for rings of `size`, or `None`.
    pub(crate) fn __getitem__(&self, py: Python<'_>, size: u8) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_size_count(size)
                .map(|v| ValueAst::from_ast(py, v))
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
    pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<BondRingSizeIter> {
        let sizes = self.read(py, |cs| Ok(ring_sizes(cs).collect::<Vec<u8>>()))?;
        Ok(BondRingSizeIter {
            sizes: sizes.into_iter(),
        })
    }

    /// Set the membership count for rings of `size` in place.
    pub(crate) fn __setitem__(&self, py: Python<'_>, size: u8, count: ValueArg) {
        let constraint =
            AstBondConstraintAst::ring_membership(AstRingScope::Size(size), count.to_ast(py));
        self.write(py, |cs| cs.set(constraint));
    }

    /// Remove the sized-ring membership for `size` in place.
    pub(crate) fn __delitem__(&self, py: Python<'_>, size: u8) {
        self.write(py, |cs| {
            cs.remove(AstBondConstraintKey::RingMembership(AstRingScope::Size(
                size,
            )));
        });
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
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
pub(crate) fn ring_sizes(constraints: &AstBondConstraintsAst) -> impl Iterator<Item = u8> + '_ {
    constraints.iter().filter_map(|entry| match entry {
        AstBondConstraintAst::RingMembership(m) => match m.scope {
            AstRingScope::Size(size) => Some(size),
            AstRingScope::All => None,
        },
        _ => None,
    })
}

#[pyclass]
pub(crate) struct BondRingSizeIter {
    sizes: IntoIter<u8>,
}

#[pymethods]
impl BondRingSizeIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<u8> {
        self.sizes.next()
    }
}

#[pyclass]
pub(crate) struct BondConstraintIter {
    entries: IntoIter<Py<BondConstraintAst>>,
}

#[pymethods]
impl BondConstraintIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<BondConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
pub(crate) struct BondConstraintKeyIter {
    keys: IntoIter<Py<BondConstraintKey>>,
}

#[pymethods]
impl BondConstraintKeyIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<Py<BondConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
pub(crate) struct BondConstraintItemsIter {
    items: IntoIter<(Py<BondConstraintKey>, Py<BondConstraintAst>)>,
}

#[pymethods]
impl BondConstraintItemsIter {
    pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub(crate) fn __next__(&mut self) -> Option<(Py<BondConstraintKey>, Py<BondConstraintAst>)> {
        self.items.next()
    }
}
