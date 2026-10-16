//! Conversion helpers shared by the Python binding types.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::ops::Deref;

use pyo3::prelude::*;
use pyo3::types::{PyFrozenSet, PyTuple};

/// A collection extracted normally and exposed as a Python tuple.
#[derive(Clone, Debug, PartialEq, Eq, Hash, FromPyObject)]
#[pyo3(transparent)]
pub struct TupleItems<T>(pub T);

impl<T> Deref for TupleItems<T> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.0
    }
}

impl<T, C: FromIterator<T>> FromIterator<T> for TupleItems<C> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<'a, 'py, T> IntoPyObject<'py> for &'a TupleItems<T>
where
    &'a T: IntoIterator,
    <&'a T as IntoIterator>::IntoIter: ExactSizeIterator,
    <&'a T as IntoIterator>::Item: IntoPyObject<'py>,
{
    type Target = PyTuple;
    type Output = Bound<'py, PyTuple>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        PyTuple::new(py, &self.0)
    }
}

impl<'py, T> IntoPyObject<'py> for TupleItems<T>
where
    T: IntoIterator,
    T::IntoIter: ExactSizeIterator,
    T::Item: IntoPyObject<'py>,
{
    type Target = PyTuple;
    type Output = Bound<'py, PyTuple>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        PyTuple::new(py, self.0)
    }
}

/// A collection extracted normally and exposed as a Python frozenset.
#[derive(Clone, Debug, PartialEq, Eq, Hash, FromPyObject)]
#[pyo3(transparent)]
pub struct FrozenSetItems<T>(pub T);

impl<T> Deref for FrozenSetItems<T> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.0
    }
}

impl<T, C: FromIterator<T>> FromIterator<T> for FrozenSetItems<C> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<'a, 'py, T, E: 'a> IntoPyObject<'py> for &'a FrozenSetItems<T>
where
    &'a T: IntoIterator<Item = &'a E>,
    E: Clone + IntoPyObject<'py>,
{
    type Target = PyFrozenSet;
    type Output = Bound<'py, PyFrozenSet>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        PyFrozenSet::new(py, self.0.into_iter().cloned())
    }
}

impl<'py, T> IntoPyObject<'py> for FrozenSetItems<T>
where
    T: IntoIterator,
    T::Item: IntoPyObject<'py>,
{
    type Target = PyFrozenSet;
    type Output = Bound<'py, PyFrozenSet>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        PyFrozenSet::new(py, self.0)
    }
}

/// Hash a Rust value with the default hasher for a binding type's `__hash__`.
pub(crate) fn hash_rust<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

/// The eval-able `__repr__` of a complex-enum variant instance: `Type.Variant(a, b)`,
/// recursing through Python `repr()` on each tuple field (`_0`, `_1`, …).
pub(crate) fn variant_repr(
    obj: &Bound<'_, PyAny>,
    type_name: &str,
    variant: &str,
    arity: usize,
) -> PyResult<String> {
    let mut parts = Vec::with_capacity(arity);
    for i in 0..arity {
        parts.push(
            obj.getattr(format!("_{i}").as_str())?
                .repr()?
                .extract::<String>()?,
        );
    }
    Ok(format!("{type_name}.{variant}({})", parts.join(", ")))
}

/// The eval-able representation of a variant with Rust-named fields.
pub(crate) fn named_variant_repr(
    obj: &Bound<'_, PyAny>,
    type_name: &str,
    variant: &str,
    fields: &[&str],
) -> PyResult<String> {
    let mut parts = Vec::with_capacity(fields.len());
    for field in fields {
        let value = obj.getattr(*field)?.repr()?.extract::<String>()?;
        parts.push(format!("{field}={value}"));
    }
    Ok(format!("{type_name}.{variant}({})", parts.join(", ")))
}

/// Wrap a complex-enum value into a `Py<T>` pointing at the **variant** instance.
///
/// `Py::new(py, value)` creates a *base*-type instance whose variant fields (`_0`,
/// …) and `match` support are absent from Python; `IntoPyObject` creates the proper
/// variant subtype. Use this for every nested `Py<…>` child in a `from_rust`.
pub(crate) fn into_py_variant<'py, T>(py: Python<'py>, value: T) -> PyResult<Py<T>>
where
    T: IntoPyObject<'py, Output = Bound<'py, T>>,
    T::Error: Into<PyErr>,
{
    value
        .into_pyobject(py)
        .map_err(Into::into)
        .map(Bound::unbind)
}
