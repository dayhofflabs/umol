//! Multicenter bond value type and multicenter-bond-constraint surface mirroring
//! `umol_ast::ast`: `MulticenterBondAst`, the `MulticenterBondConstraintAst` enum, the
//! `MulticenterBondConstraintsAst` container, and the `MulticenterBondConstraintsView`
//! live handle. A multicenter bond is a single unordered set of member atoms; the
//! value carries a positional per-atom `electrons` vector plus charge, spin, and
//! constraints. The member atoms are the participants of the owning molecule's
//! multicenter relation, so they are topology (the view half) rather than value.

use std::vec::IntoIter;

use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use umol_ast::ast::{
    MulticenterBondConstraintAst as AstMulticenterBondConstraintAst,
    MulticenterBondConstraintKey as AstMulticenterBondConstraintKey,
    MulticenterBondConstraintsAst as AstMulticenterBondConstraintsAst,
};

use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::value::{ValueArg, ValueAst};

/// The key (identity) of a multicenter-bond constraint, for keyed lookup. The
/// single key `ElectronCount` is the bare discriminant (no sub-key).
#[pyclass]
pub enum MulticenterBondConstraintKey {
    ElectronCount(),
}

#[pymethods]
impl MulticenterBondConstraintKey {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
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
    pub(crate) fn from_ast(ast: &AstMulticenterBondConstraintKey) -> Self {
        match ast {
            AstMulticenterBondConstraintKey::ElectronCount => Self::ElectronCount(),
        }
    }

    pub(crate) fn to_ast(&self) -> AstMulticenterBondConstraintKey {
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
    fn key(&self, py: Python<'_>) -> MulticenterBondConstraintKey {
        MulticenterBondConstraintKey::from_ast(&self.to_ast(py).key())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
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

impl MulticenterBondConstraintAst {
    pub(crate) fn from_ast(
        py: Python<'_>,
        ast: &AstMulticenterBondConstraintAst,
    ) -> PyResult<Self> {
        Ok(match ast {
            AstMulticenterBondConstraintAst::ElectronCount(v) => {
                Self::ElectronCount(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstMulticenterBondConstraintAst {
        match self {
            Self::ElectronCount(v) => {
                AstMulticenterBondConstraintAst::ElectronCount(v.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container or an iterable of
/// `MulticenterBondConstraintAst` (each `set`, last-wins). The live-view variant
/// lands with `MulticenterBondConstraintsView` (S1b).
#[derive(FromPyObject)]
enum MulticenterBondConstraintsUpdate {
    Container(Py<MulticenterBondConstraintsAst>),
    Entries(Vec<Py<MulticenterBondConstraintAst>>),
}

impl MulticenterBondConstraintsUpdate {
    /// Overlay this update onto `target` in place.
    fn apply(&self, py: Python<'_>, target: &mut AstMulticenterBondConstraintsAst) -> PyResult<()> {
        match self {
            MulticenterBondConstraintsUpdate::Container(c) => {
                target.update(c.bind(py).borrow().inner())
            }
            MulticenterBondConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry.bind(py).borrow().to_ast(py));
                }
            }
        }
        Ok(())
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
    fn new(py: Python<'_>, entries: Vec<Py<MulticenterBondConstraintAst>>) -> Self {
        let mut constraints = AstMulticenterBondConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_ast(py)),
        );
        MulticenterBondConstraintsAst(constraints)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let mirror = into_py_variant(py, MulticenterBondConstraintAst::from_ast(py, entry)?)?;
            parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!(
            "MulticenterBondConstraintsAst([{}])",
            parts.join(", ")
        ))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    fn set(&mut self, py: Python<'_>, c: Py<MulticenterBondConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<Option<MulticenterBondConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_ast())
            .map(|c| MulticenterBondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container or an iterable of
    /// `MulticenterBondConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&mut self, py: Python<'_>, other: MulticenterBondConstraintsUpdate) -> PyResult<()> {
        other.apply(py, &mut self.0)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintKeyIter> {
        multicenter_bond_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintKeyIter> {
        multicenter_bond_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintIter> {
        multicenter_bond_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintItemsIter> {
        multicenter_bond_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_ast()) {
            Some(constraint) => Ok(into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(py, constraint)?,
            )?
            .into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<MulticenterBondConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_ast()) {
            Some(constraint) => MulticenterBondConstraintAst::from_ast(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_ast()).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<MulticenterBondConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast())
    }

    /// The asserted total electron count; `Undetermined` when no `ElectronCount`
    /// constraint is present (mirroring the non-optional Rust accessor).
    #[getter]
    fn electron_count(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.electron_count())
    }

    #[setter]
    fn set_electron_count(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstMulticenterBondConstraintAst::electron_count(
            value.to_ast(py),
        ));
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
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

/// Build the per-constraint iterator handle from a borrowed container.
fn multicenter_bond_constraints_iter(
    py: Python<'_>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<MulticenterBondConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(py, MulticenterBondConstraintAst::from_ast(py, constraint)?)
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(MulticenterBondConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
fn multicenter_bond_constraint_keys(
    py: Python<'_>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<MulticenterBondConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(
                py,
                MulticenterBondConstraintKey::from_ast(&constraint.key()),
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(MulticenterBondConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
fn multicenter_bond_constraint_items(
    py: Python<'_>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<MulticenterBondConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(
                    py,
                    MulticenterBondConstraintKey::from_ast(&constraint.key()),
                )?,
                into_py_variant(py, MulticenterBondConstraintAst::from_ast(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(MulticenterBondConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// inner-value mirrors.
fn multicenter_bond_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            AstMulticenterBondConstraintAst::ElectronCount(v) => {
                dict.set_item("electron_count", ValueAst::from_ast(py, v)?)?
            }
        }
    }
    Ok(dict)
}

#[pyclass]
struct MulticenterBondConstraintIter {
    entries: IntoIter<Py<MulticenterBondConstraintAst>>,
}

#[pymethods]
impl MulticenterBondConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<MulticenterBondConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
struct MulticenterBondConstraintKeyIter {
    keys: IntoIter<Py<MulticenterBondConstraintKey>>,
}

#[pymethods]
impl MulticenterBondConstraintKeyIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<MulticenterBondConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
struct MulticenterBondConstraintItemsIter {
    items: IntoIter<(
        Py<MulticenterBondConstraintKey>,
        Py<MulticenterBondConstraintAst>,
    )>,
}

#[pymethods]
impl MulticenterBondConstraintItemsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(
        &mut self,
    ) -> Option<(
        Py<MulticenterBondConstraintKey>,
        Py<MulticenterBondConstraintAst>,
    )> {
        self.items.next()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::ValueAst as AstValueAst;

    use super::*;

    #[rstest]
    fn test_multicenter_bond_constraint_key_roundtrip() {
        let key =
            MulticenterBondConstraintKey::from_ast(&AstMulticenterBondConstraintKey::ElectronCount);
        assert_eq!(key.to_ast(), AstMulticenterBondConstraintKey::ElectronCount);
    }

    #[rstest]
    fn test_multicenter_bond_constraint_ast_key() {
        Python::attach(|py| {
            let constraint = AstMulticenterBondConstraintAst::electron_count(6);
            let key = MulticenterBondConstraintAst::from_ast(py, &constraint)
                .unwrap()
                .key(py);
            assert_eq!(key.to_ast(), AstMulticenterBondConstraintKey::ElectronCount);
        });
    }

    #[rstest]
    #[case(AstMulticenterBondConstraintAst::electron_count(6))]
    #[case(AstMulticenterBondConstraintAst::electron_count(AstValueAst::Undetermined))]
    fn test_multicenter_bond_constraint_ast_roundtrip(
        #[case] ast: AstMulticenterBondConstraintAst,
    ) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterBondConstraintAst::from_ast(py, &ast)
                    .unwrap()
                    .to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_new() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_repr() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            assert_eq!(
                constraints.__repr__(py).unwrap(),
                "MulticenterBondConstraintsAst([MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))])"
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            constraints.set(py, ec);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_pop() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            let key = into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap();
            let removed = constraints.pop(py, key).unwrap();
            match removed {
                Some(MulticenterBondConstraintAst::ElectronCount(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(6))
                }
                _ => panic!("expected removed ElectronCount(Lit(6))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_update() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            let mut other = AstMulticenterBondConstraintsAst::new();
            other.set(AstMulticenterBondConstraintAst::electron_count(6));
            constraints
                .update(
                    py,
                    MulticenterBondConstraintsUpdate::Container(
                        Py::new(py, MulticenterBondConstraintsAst::from_inner(other)).unwrap(),
                    ),
                )
                .unwrap();
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_update_entries() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            constraints
                .update(py, MulticenterBondConstraintsUpdate::Entries(vec![ec]))
                .unwrap();
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_len_contains() {
        Python::attach(|py| {
            let empty = MulticenterBondConstraintsAst::new(py, vec![]);
            assert_eq!(empty.__len__(), 0);
            assert!(!empty.__contains__(
                py,
                into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap()
            ));
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap()
            ));
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_keys_values_items() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(),
                AstMulticenterBondConstraintKey::ElectronCount
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstMulticenterBondConstraintAst::electron_count(6)
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_ast(),
                AstMulticenterBondConstraintKey::ElectronCount
            );
            assert_eq!(
                value.bind(py).borrow().to_ast(py),
                AstMulticenterBondConstraintAst::electron_count(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_get() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            let present = constraints
                .get(
                    py,
                    into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap(),
                    None,
                )
                .unwrap();
            let expected = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap()
            .into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());

            let empty = MulticenterBondConstraintsAst::new(py, vec![]);
            let absent = empty
                .get(
                    py,
                    into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap(),
                    None,
                )
                .unwrap();
            assert!(absent.bind(py).is_none());

            // a caller-supplied default is returned verbatim when the key is absent
            let sentinel = into_py_variant(py, MulticenterBondConstraintKey::ElectronCount())
                .unwrap()
                .into_any();
            let defaulted = empty
                .get(
                    py,
                    into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap(),
                    Some(sentinel.clone_ref(py)),
                )
                .unwrap();
            assert_eq!(defaulted.as_ptr(), sentinel.as_ptr());
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_electron_count() {
        Python::attach(|py| {
            let empty = MulticenterBondConstraintsAst::new(py, vec![]);
            assert_eq!(
                empty.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Undetermined
            );
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_set_electron_count() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            constraints.set_electron_count(py, ValueArg::Lit(6));
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap();
            assert!(constraints.__getitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap();
            assert!(constraints.__delitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_asdict() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            let dict = constraints.asdict(py).unwrap();
            assert_eq!(dict.len(), 1);
            let value = dict.get_item("electron_count").unwrap().unwrap();
            let expected =
                into_py_variant(py, ValueAst::from_ast(py, &AstValueAst::Lit(6)).unwrap()).unwrap();
            assert!(value.eq(expected.bind(py)).unwrap());
        });
    }
}
