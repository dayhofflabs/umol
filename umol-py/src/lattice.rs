//! Shared generation of Python methods for form lattice and normalization operations.

pub(crate) use umol_graph_ir::ir::{Lattice, Normalize};

pub(crate) use crate::error::contradiction_error;

macro_rules! impl_py_normalize {
    ($py_type:ty, $rust_type:ty, $to_rust:expr, $from_rust:expr) => {
        #[pymethods]
        impl $py_type {
            /// Return the normal form without mutating the receiver.
            fn normalize(&self, py: Python<'_>) -> PyResult<Self> {
                let to_rust = $to_rust;
                let from_rust = $from_rust;
                let value: $rust_type = to_rust(self, py)?;
                let normalized = $crate::lattice::Normalize::normalize(value)
                    .map_err($crate::lattice::contradiction_error)?;
                from_rust(py, normalized)
            }

            /// Compare normal forms while leaving structural equality unchanged.
            fn normalized_eq(&self, py: Python<'_>, other: &$py_type) -> PyResult<bool> {
                let to_rust = $to_rust;
                let lhs: $rust_type = to_rust(self, py)?;
                let rhs: $rust_type = to_rust(other, py)?;
                Ok($crate::lattice::Normalize::normalized_eq(&lhs, &rhs))
            }
        }
    };
}

macro_rules! impl_py_lattice {
    ($py_type:ty, $rust_type:ty, $to_rust:expr, $from_rust:expr) => {
        $crate::lattice::impl_py_normalize!($py_type, $rust_type, $to_rust, $from_rust);

        #[pymethods]
        impl $py_type {
            /// Whether this value is the undetermined top of its lattice.
            fn is_undetermined(&self, py: Python<'_>) -> PyResult<bool> {
                let to_rust = $to_rust;
                let value: $rust_type = to_rust(self, py)?;
                Ok($crate::lattice::Lattice::is_undetermined(&value))
            }

            /// Whether this value resolves to one concrete value.
            fn is_ground(&self, py: Python<'_>) -> PyResult<bool> {
                let to_rust = $to_rust;
                let value: $rust_type = to_rust(self, py)?;
                Ok($crate::lattice::Lattice::is_ground(&value))
            }

            /// Return the greatest lower bound, or `None` when the values are incompatible.
            fn meet(&self, py: Python<'_>, other: &$py_type) -> PyResult<Option<Self>> {
                let to_rust = $to_rust;
                let from_rust = $from_rust;
                let lhs: $rust_type = to_rust(self, py)?;
                let rhs: $rust_type = to_rust(other, py)?;
                match $crate::lattice::Lattice::meet(&lhs, &rhs) {
                    Some(value) => from_rust(py, value).map(Some),
                    None => Ok(None),
                }
            }

            /// Return the least upper bound, or `None` when no join exists.
            fn join(&self, py: Python<'_>, other: &$py_type) -> PyResult<Option<Self>> {
                let to_rust = $to_rust;
                let from_rust = $from_rust;
                let lhs: $rust_type = to_rust(self, py)?;
                let rhs: $rust_type = to_rust(other, py)?;
                match $crate::lattice::Lattice::join(&lhs, &rhs) {
                    Ok(value) => from_rust(py, value).map(Some),
                    Err(_) => Ok(None),
                }
            }

            /// Whether `target` refines this pattern.
            pub(crate) fn matches(&self, py: Python<'_>, target: &$py_type) -> PyResult<bool> {
                let to_rust = $to_rust;
                let pattern: $rust_type = to_rust(self, py)?;
                let target: $rust_type = to_rust(target, py)?;
                Ok($crate::lattice::Lattice::matches(&pattern, &target))
            }

            /// Whether the two values admit a common ground refinement.
            fn is_compatible(&self, py: Python<'_>, other: &$py_type) -> PyResult<bool> {
                let to_rust = $to_rust;
                let lhs: $rust_type = to_rust(self, py)?;
                let rhs: $rust_type = to_rust(other, py)?;
                Ok($crate::lattice::Lattice::is_compatible(&lhs, &rhs))
            }
        }
    };
}

pub(crate) use impl_py_lattice;
pub(crate) use impl_py_normalize;
