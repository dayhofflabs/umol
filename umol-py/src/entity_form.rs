//! Shared storage for owned Python entity-form values.

use std::fmt;
use std::ops::Deref;

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::PyClass;

/// An owned entity form with Python-side write access metadata.
///
/// Read-only state is deliberately not part of value equality: it describes how
/// this particular Python object may be used, not the represented form.
pub(crate) struct EntityFormValue<T> {
    value: T,
    readonly: bool,
}

impl<T> EntityFormValue<T> {
    pub(crate) fn writable(value: T) -> Self {
        Self {
            value,
            readonly: false,
        }
    }

    pub(crate) fn readonly(value: T) -> Self {
        Self {
            value,
            readonly: true,
        }
    }

    pub(crate) fn is_readonly(&self) -> bool {
        self.readonly
    }

    pub(crate) fn value(&self) -> &T {
        &self.value
    }

    pub(crate) fn value_mut(&mut self) -> PyResult<&mut T> {
        if self.readonly {
            Err(PyTypeError::new_err("read-only entity form"))
        } else {
            Ok(&mut self.value)
        }
    }
}

impl<T> Deref for EntityFormValue<T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.value
    }
}

impl<T: PartialEq> PartialEq for EntityFormValue<T> {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<T: fmt::Display> fmt::Display for EntityFormValue<T> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.value.fmt(formatter)
    }
}

/// Operations needed by the generic read-only form holder used in deltas.
pub trait EntityForm: PyClass {
    type RustForm: Clone;

    fn clone_rust(&self) -> Self::RustForm;
    fn new_readonly(py: Python<'_>, value: Self::RustForm) -> PyResult<Py<Self>>;
}

/// A retained read-only Python form used as a field of an immutable delta.
pub struct ReadonlyForm<T: EntityForm>(Py<T>);

impl<'a, 'py, T> FromPyObject<'a, 'py> for ReadonlyForm<T>
where
    T: EntityForm,
{
    type Error = PyErr;

    fn extract(obj: Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        let source = obj.extract::<PyRef<'_, T>>()?;
        let value = source.clone_rust();
        drop(source);
        Ok(Self(T::new_readonly(obj.py(), value)?))
    }
}

impl<'py, T> IntoPyObject<'py> for &ReadonlyForm<T>
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

impl<T> ReadonlyForm<T>
where
    T: EntityForm,
{
    pub(crate) fn from_rust(py: Python<'_>, value: &T::RustForm) -> PyResult<Self> {
        Ok(Self(T::new_readonly(py, value.clone())?))
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> T::RustForm {
        self.0.bind(py).borrow().clone_rust()
    }
}
