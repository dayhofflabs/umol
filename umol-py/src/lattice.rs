//! Shared generation of Python methods for form lattice and normalization operations.

pub(crate) use umol_graph_ir::ir::{Lattice, Normalize};

pub(crate) use crate::error::{contradiction_error, NoJoinError};

macro_rules! impl_py_normalize {
    (borrowed, $py_type:ty, $rust_type:ty, $from_rust:expr) => {
        #[pymethods]
        impl $py_type {
            /// Return the normal form without mutating the receiver.
            fn normalize(&self, py: Python<'_>) -> PyResult<Self> {
                let from_rust = $from_rust;
                let value: &$rust_type = self.to_rust();
                let normalized = $crate::lattice::Normalize::normalize(value.clone())
                    .map_err($crate::lattice::contradiction_error)?;
                from_rust(py, normalized)
            }

            /// Compare normal forms while leaving structural equality unchanged.
            fn normalized_eq(&self, _py: Python<'_>, other: &$py_type) -> PyResult<bool> {
                let lhs: &$rust_type = self.to_rust();
                let rhs: &$rust_type = other.to_rust();
                Ok($crate::lattice::Normalize::normalized_eq(lhs, rhs))
            }
        }
    };
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
    (borrowed, $py_type:ty, $rust_type:ty, $from_rust:expr) => {
        $crate::lattice::impl_py_normalize!(borrowed, $py_type, $rust_type, $from_rust);

        #[pymethods]
        impl $py_type {
            /// Whether this value is the undetermined top of its lattice.
            fn is_undetermined(&self, _py: Python<'_>) -> PyResult<bool> {
                let value: &$rust_type = self.to_rust();
                Ok($crate::lattice::Lattice::is_undetermined(value))
            }

            /// Whether this value resolves to one concrete value.
            fn is_ground(&self, _py: Python<'_>) -> PyResult<bool> {
                let value: &$rust_type = self.to_rust();
                Ok($crate::lattice::Lattice::is_ground(value))
            }

            /// Return the greatest lower bound, or `None` when the values are incompatible.
            fn meet(&self, py: Python<'_>, other: &$py_type) -> PyResult<Option<Self>> {
                let from_rust = $from_rust;
                let lhs: &$rust_type = self.to_rust();
                let rhs: &$rust_type = other.to_rust();
                match $crate::lattice::Lattice::meet(lhs, rhs) {
                    Some(value) => from_rust(py, value).map(Some),
                    None => Ok(None),
                }
            }

            /// Return the least upper bound, raising `NoJoinError` when no join exists.
            fn join(&self, py: Python<'_>, other: &$py_type) -> PyResult<Self> {
                let from_rust = $from_rust;
                let lhs: &$rust_type = self.to_rust();
                let rhs: &$rust_type = other.to_rust();
                let value = $crate::lattice::Lattice::join(lhs, rhs)
                    .map_err(|error| $crate::lattice::NoJoinError::new_err(error.to_string()))?;
                from_rust(py, value)
            }

            /// Whether `target` refines this pattern.
            pub(crate) fn matches(&self, _py: Python<'_>, target: &$py_type) -> PyResult<bool> {
                let pattern: &$rust_type = self.to_rust();
                let target: &$rust_type = target.to_rust();
                Ok($crate::lattice::Lattice::matches(pattern, target))
            }

            /// Whether the two values admit a common ground refinement.
            fn is_compatible(&self, _py: Python<'_>, other: &$py_type) -> PyResult<bool> {
                let lhs: &$rust_type = self.to_rust();
                let rhs: &$rust_type = other.to_rust();
                Ok($crate::lattice::Lattice::is_compatible(lhs, rhs))
            }
        }
    };
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

            /// Return the least upper bound, raising `NoJoinError` when no join exists.
            fn join(&self, py: Python<'_>, other: &$py_type) -> PyResult<Self> {
                let to_rust = $to_rust;
                let from_rust = $from_rust;
                let lhs: $rust_type = to_rust(self, py)?;
                let rhs: $rust_type = to_rust(other, py)?;
                let value = $crate::lattice::Lattice::join(&lhs, &rhs)
                    .map_err(|error| $crate::lattice::NoJoinError::new_err(error.to_string()))?;
                from_rust(py, value)
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
