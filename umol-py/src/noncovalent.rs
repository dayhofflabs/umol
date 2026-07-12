//! Noncovalent-bond bindings. This file grows over the B5 slice; it opens with the
//! kind leaf: `NoncovalentBondKind` (the interaction kind) and `NoncovalentBondKindAst`
//! (`Undetermined | Lit(kind)`), mirroring `umol_ast::ast::{NoncovalentBondKind,
//! NoncovalentBondKindAst}` — the noncovalent analog of `atom.element: ElementAst`
//! over the `Element` value enum.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::vec::IntoIter;

use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use umol_ast::ast::{
    AsLit, NoncovalentBondConstraintAst as AstNoncovalentBondConstraintAst,
    NoncovalentBondConstraintKey as AstNoncovalentBondConstraintKey,
    NoncovalentBondConstraintsAst as AstNoncovalentBondConstraintsAst,
    NoncovalentBondKind as AstNoncovalentBondKind,
    NoncovalentBondKindAst as AstNoncovalentBondKindAst,
};

use crate::boolean::{BooleanArg, BooleanAst};
use crate::convert::{hash_ast, into_py_variant, variant_repr};

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
    #[cfg(test)]
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
#[cfg(test)]
#[derive(FromPyObject)]
pub(crate) enum NoncovalentBondKindArg {
    Kind(NoncovalentBondKind),
    Ast(Py<NoncovalentBondKindAst>),
}

#[cfg(test)]
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

/// The argument to `update`: another constraint container or an iterable of
/// `NoncovalentBondConstraintAst` (each `set`, last-wins). The live-view variant is
/// added with the view.
#[derive(FromPyObject)]
enum NoncovalentBondConstraintsUpdate {
    Container(Py<NoncovalentBondConstraintsAst>),
    Entries(Vec<Py<NoncovalentBondConstraintAst>>),
}

impl NoncovalentBondConstraintsUpdate {
    /// Overlay this update onto `target` in place.
    fn apply(&self, py: Python<'_>, target: &mut AstNoncovalentBondConstraintsAst) -> PyResult<()> {
        match self {
            NoncovalentBondConstraintsUpdate::Container(c) => {
                target.update(c.bind(py).borrow().inner())
            }
            NoncovalentBondConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry.bind(py).borrow().to_ast(py));
                }
            }
        }
        Ok(())
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

    /// Overlay `other` onto self in place — another container or an iterable of
    /// `NoncovalentBondConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&mut self, py: Python<'_>, other: NoncovalentBondConstraintsUpdate) -> PyResult<()> {
        other.apply(py, &mut self.0)
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

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::BooleanAst as AstBooleanAst;

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
            let mut constraints = NoncovalentBondConstraintsAst::new(py, vec![]);
            let mut other = AstNoncovalentBondConstraintsAst::new();
            other.set(AstNoncovalentBondConstraintAst::intramolecular(true));
            constraints
                .update(
                    py,
                    NoncovalentBondConstraintsUpdate::Container(
                        Py::new(py, NoncovalentBondConstraintsAst::from_inner(other)).unwrap(),
                    ),
                )
                .unwrap();
            assert_eq!(
                constraints.intramolecular().to_ast(),
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_update_entries() {
        Python::attach(|py| {
            let mut constraints = NoncovalentBondConstraintsAst::new(py, vec![]);
            constraints
                .update(
                    py,
                    NoncovalentBondConstraintsUpdate::Entries(vec![intramolecular(py, false)]),
                )
                .unwrap();
            assert_eq!(
                constraints.intramolecular().to_ast(),
                AstBooleanAst::Lit(false)
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
}
