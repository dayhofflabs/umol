//! Conversion helpers shared by the Python binding types.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use pyo3::prelude::*;

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
