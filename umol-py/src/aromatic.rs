//! Aromatic system value type and aromatic-constraint surface mirroring
//! `umol_ast::ast`: `AromaticSystemAst`, the `AromaticSystemConstraintAst` enum, the
//! `AromaticSystemConstraintsAst` container, and the `AromaticSystemConstraintsView`
//! live handle. An aromatic system is a single unordered set of member atoms; the
//! value carries a positional per-atom `electrons` vector plus charge, spin, and
//! constraints. The member atoms are the participants of the owning molecule's
//! aromatic relation, so they are topology (the view half) rather than value.

use std::vec::IntoIter;

use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use umol_ast::ast::{
    AromaticSystemConstraintAst as AstAromaticSystemConstraintAst,
    AromaticSystemConstraintKey as AstAromaticSystemConstraintKey,
    AromaticSystemConstraintsAst as AstAromaticSystemConstraintsAst,
};

use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::value::{ValueArg, ValueAst};

/// The key (identity) of an aromatic-system constraint, for keyed lookup. The
/// single key `ElectronCount` is the bare discriminant (no sub-key).
#[pyclass]
pub enum AromaticSystemConstraintKey {
    ElectronCount(),
}

#[pymethods]
impl AromaticSystemConstraintKey {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
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
    pub(crate) fn from_ast(ast: &AstAromaticSystemConstraintKey) -> Self {
        match ast {
            AstAromaticSystemConstraintKey::ElectronCount => Self::ElectronCount(),
        }
    }

    pub(crate) fn to_ast(&self) -> AstAromaticSystemConstraintKey {
        match self {
            Self::ElectronCount() => AstAromaticSystemConstraintKey::ElectronCount,
        }
    }
}

/// An aromatic-system-scope constraint: the asserted total π-electron count of the
/// system (cross-checked against `sum(AromaticSystemAst::electrons)`).
#[pyclass]
pub enum AromaticSystemConstraintAst {
    ElectronCount(Py<ValueAst>),
}

#[pymethods]
impl AromaticSystemConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    fn key(&self, py: Python<'_>) -> AromaticSystemConstraintKey {
        AromaticSystemConstraintKey::from_ast(&self.to_ast(py).key())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
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

impl AromaticSystemConstraintAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAromaticSystemConstraintAst) -> PyResult<Self> {
        Ok(match ast {
            AstAromaticSystemConstraintAst::ElectronCount(v) => {
                Self::ElectronCount(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAromaticSystemConstraintAst {
        match self {
            Self::ElectronCount(v) => {
                AstAromaticSystemConstraintAst::ElectronCount(v.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container or an iterable of
/// `AromaticSystemConstraintAst` (each `set`, last-wins). The live-view variant lands
/// with `AromaticSystemConstraintsView` (S1c).
#[derive(FromPyObject)]
enum AromaticSystemConstraintsUpdate {
    Container(Py<AromaticSystemConstraintsAst>),
    Entries(Vec<Py<AromaticSystemConstraintAst>>),
}

impl AromaticSystemConstraintsUpdate {
    /// Overlay this update onto `target` in place.
    fn apply(&self, py: Python<'_>, target: &mut AstAromaticSystemConstraintsAst) -> PyResult<()> {
        match self {
            AromaticSystemConstraintsUpdate::Container(c) => {
                target.update(c.bind(py).borrow().inner())
            }
            AromaticSystemConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry.bind(py).borrow().to_ast(py));
                }
            }
        }
        Ok(())
    }
}

/// The aromatic-system-scope constraints on an aromatic system, in kind-sorted order.
/// Mutable, hence value-equal but unhashable (matching `AromaticSystemAst`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct AromaticSystemConstraintsAst(AstAromaticSystemConstraintsAst);

#[pymethods]
impl AromaticSystemConstraintsAst {
    /// Build from a sequence of constraints (a later entry of the same key replaces
    /// an earlier one, last-wins).
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<AromaticSystemConstraintAst>>) -> Self {
        let mut constraints = AstAromaticSystemConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_ast(py)),
        );
        AromaticSystemConstraintsAst(constraints)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let mirror = into_py_variant(py, AromaticSystemConstraintAst::from_ast(py, entry)?)?;
            parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!(
            "AromaticSystemConstraintsAst([{}])",
            parts.join(", ")
        ))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    fn set(&mut self, py: Python<'_>, c: Py<AromaticSystemConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<Option<AromaticSystemConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_ast())
            .map(|c| AromaticSystemConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container or an iterable of
    /// `AromaticSystemConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&mut self, py: Python<'_>, other: AromaticSystemConstraintsUpdate) -> PyResult<()> {
        other.apply(py, &mut self.0)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintKeyIter> {
        aromatic_system_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintKeyIter> {
        aromatic_system_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintIter> {
        aromatic_system_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintItemsIter> {
        aromatic_system_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_ast()) {
            Some(constraint) => Ok(into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(py, constraint)?,
            )?
            .into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<AromaticSystemConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_ast()) {
            Some(constraint) => AromaticSystemConstraintAst::from_ast(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_ast()).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<AromaticSystemConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast())
    }

    /// The asserted total π-electron count; `Undetermined` when no `ElectronCount`
    /// constraint is present (mirroring the non-optional Rust accessor).
    #[getter]
    fn electron_count(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.electron_count())
    }

    #[setter]
    fn set_electron_count(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstAromaticSystemConstraintAst::electron_count(
            value.to_ast(py),
        ));
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        aromatic_system_constraints_asdict(py, &self.0)
    }
}

impl AromaticSystemConstraintsAst {
    /// The wrapped AST constraints — read access for aromatic system construction.
    pub(crate) fn inner(&self) -> &AstAromaticSystemConstraintsAst {
        &self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge). Test-only —
    /// in-crate construction wraps `AromaticSystemConstraintsAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(constraints: AstAromaticSystemConstraintsAst) -> Self {
        AromaticSystemConstraintsAst(constraints)
    }
}

/// Build the per-constraint iterator handle from a borrowed container.
fn aromatic_system_constraints_iter(
    py: Python<'_>,
    constraints: &AstAromaticSystemConstraintsAst,
) -> PyResult<AromaticSystemConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(py, AromaticSystemConstraintAst::from_ast(py, constraint)?)
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AromaticSystemConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
fn aromatic_system_constraint_keys(
    py: Python<'_>,
    constraints: &AstAromaticSystemConstraintsAst,
) -> PyResult<AromaticSystemConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(py, AromaticSystemConstraintKey::from_ast(&constraint.key()))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AromaticSystemConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
fn aromatic_system_constraint_items(
    py: Python<'_>,
    constraints: &AstAromaticSystemConstraintsAst,
) -> PyResult<AromaticSystemConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(py, AromaticSystemConstraintKey::from_ast(&constraint.key()))?,
                into_py_variant(py, AromaticSystemConstraintAst::from_ast(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AromaticSystemConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// inner-value mirrors.
fn aromatic_system_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &AstAromaticSystemConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            AstAromaticSystemConstraintAst::ElectronCount(v) => {
                dict.set_item("electron_count", ValueAst::from_ast(py, v)?)?
            }
        }
    }
    Ok(dict)
}

#[pyclass]
struct AromaticSystemConstraintIter {
    entries: IntoIter<Py<AromaticSystemConstraintAst>>,
}

#[pymethods]
impl AromaticSystemConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<AromaticSystemConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
struct AromaticSystemConstraintKeyIter {
    keys: IntoIter<Py<AromaticSystemConstraintKey>>,
}

#[pymethods]
impl AromaticSystemConstraintKeyIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<AromaticSystemConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
struct AromaticSystemConstraintItemsIter {
    items: IntoIter<(
        Py<AromaticSystemConstraintKey>,
        Py<AromaticSystemConstraintAst>,
    )>,
}

#[pymethods]
impl AromaticSystemConstraintItemsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(
        &mut self,
    ) -> Option<(
        Py<AromaticSystemConstraintKey>,
        Py<AromaticSystemConstraintAst>,
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
    fn test_aromatic_system_constraint_key_roundtrip() {
        let key =
            AromaticSystemConstraintKey::from_ast(&AstAromaticSystemConstraintKey::ElectronCount);
        assert_eq!(key.to_ast(), AstAromaticSystemConstraintKey::ElectronCount);
    }

    #[rstest]
    fn test_aromatic_system_constraint_ast_key() {
        Python::attach(|py| {
            let constraint = AstAromaticSystemConstraintAst::electron_count(6);
            let key = AromaticSystemConstraintAst::from_ast(py, &constraint)
                .unwrap()
                .key(py);
            assert_eq!(key.to_ast(), AstAromaticSystemConstraintKey::ElectronCount);
        });
    }

    #[rstest]
    #[case(AstAromaticSystemConstraintAst::electron_count(6))]
    #[case(AstAromaticSystemConstraintAst::electron_count(AstValueAst::Undetermined))]
    fn test_aromatic_system_constraint_ast_roundtrip(#[case] ast: AstAromaticSystemConstraintAst) {
        Python::attach(|py| {
            assert_eq!(
                AromaticSystemConstraintAst::from_ast(py, &ast)
                    .unwrap()
                    .to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_new() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            assert_eq!(constraints.__len__(), 1);
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_repr() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            assert_eq!(
                constraints.__repr__(py).unwrap(),
                "AromaticSystemConstraintsAst([AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))])"
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
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
    fn test_aromatic_system_constraints_ast_pop() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            let key = into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap();
            let removed = constraints.pop(py, key).unwrap();
            match removed {
                Some(AromaticSystemConstraintAst::ElectronCount(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(6))
                }
                _ => panic!("expected removed ElectronCount(Lit(6))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_update() {
        Python::attach(|py| {
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            let mut other = AstAromaticSystemConstraintsAst::new();
            other.set(AstAromaticSystemConstraintAst::electron_count(6));
            constraints
                .update(
                    py,
                    AromaticSystemConstraintsUpdate::Container(
                        Py::new(py, AromaticSystemConstraintsAst::from_inner(other)).unwrap(),
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
    fn test_aromatic_system_constraints_ast_update_entries() {
        Python::attach(|py| {
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            constraints
                .update(py, AromaticSystemConstraintsUpdate::Entries(vec![ec]))
                .unwrap();
            assert_eq!(constraints.__len__(), 1);
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_len_contains() {
        Python::attach(|py| {
            let empty = AromaticSystemConstraintsAst::new(py, vec![]);
            assert_eq!(empty.__len__(), 0);
            assert!(!empty.__contains__(
                py,
                into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap()
            ));
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap()
            ));
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_keys_values_items() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(),
                AstAromaticSystemConstraintKey::ElectronCount
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstAromaticSystemConstraintAst::electron_count(6)
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_ast(),
                AstAromaticSystemConstraintKey::ElectronCount
            );
            assert_eq!(
                value.bind(py).borrow().to_ast(py),
                AstAromaticSystemConstraintAst::electron_count(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_get() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            let present = constraints
                .get(
                    py,
                    into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap(),
                    None,
                )
                .unwrap();
            let expected = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap()
            .into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());

            let empty = AromaticSystemConstraintsAst::new(py, vec![]);
            let absent = empty
                .get(
                    py,
                    into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap(),
                    None,
                )
                .unwrap();
            assert!(absent.bind(py).is_none());
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_electron_count() {
        Python::attach(|py| {
            let empty = AromaticSystemConstraintsAst::new(py, vec![]);
            assert_eq!(
                empty.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Undetermined
            );
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_set_electron_count() {
        Python::attach(|py| {
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            constraints.set_electron_count(py, ValueArg::Lit(6));
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap();
            assert!(constraints.__getitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap();
            assert!(constraints.__delitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_asdict() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            let dict = constraints.asdict(py).unwrap();
            assert_eq!(dict.len(), 1);
            let value = dict.get_item("electron_count").unwrap().unwrap();
            let expected =
                into_py_variant(py, ValueAst::from_ast(py, &AstValueAst::Lit(6)).unwrap()).unwrap();
            assert!(value.eq(expected.bind(py)).unwrap());
        });
    }
}
