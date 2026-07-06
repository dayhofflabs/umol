//! Conversion helpers shared by the mirror types.

use pyo3::prelude::*;

/// Wrap a complex-enum value into a `Py<T>` pointing at the **variant** instance.
///
/// `Py::new(py, value)` creates a *base*-type instance whose variant fields (`_0`,
/// …) and `match` support are absent from Python; `IntoPyObject` creates the proper
/// variant subtype. Use this for every nested `Py<…>` child in a `from_ast`.
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
