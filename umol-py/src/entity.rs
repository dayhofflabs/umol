//! Cross-entity support for owned Python entity forms.

use pyo3::prelude::*;
use pyo3::PyClass;

/// Operations needed by entity forms retained in immutable mutation values.
pub trait EntityForm: PyClass {
    type RustForm: Clone;

    fn to_rust(&self) -> &Self::RustForm;
    fn new_readonly(py: Python<'_>, value: Self::RustForm) -> PyResult<Py<Self>>;
}

/// A retained read-only Python form used as a field of an immutable delta or edit.
pub struct Readonly<T: EntityForm>(Py<T>);

impl<'a, 'py, T> FromPyObject<'a, 'py> for Readonly<T>
where
    T: EntityForm,
{
    type Error = PyErr;

    fn extract(obj: Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, T>>()?;
        let value = source.to_rust().clone();
        drop(source);
        Ok(Self(T::new_readonly(obj.py(), value)?))
    }
}

impl<'py, T> IntoPyObject<'py> for &Readonly<T>
where
    T: EntityForm,
{
    type Target = T;
    type Output = Bound<'py, T>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Ok(self.0.clone_ref(py).into_bound(py))
    }
}

impl<T> Readonly<T>
where
    T: EntityForm,
{
    pub(crate) fn from_rust(py: Python<'_>, value: &T::RustForm) -> PyResult<Self> {
        Ok(Self(T::new_readonly(py, value.clone())?))
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> T::RustForm {
        self.0.bind(py).borrow().to_rust().clone()
    }
}
