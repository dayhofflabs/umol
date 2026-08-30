//! Ring scope and membership constraint payloads used by atom and bond kinds.

use pyo3::prelude::*;
use umol_graph_ir::ir::{
    RingMembershipForm as GraphIrRingMembershipForm, RingScope as GraphIrRingScope,
};

use crate::convert::{hash_rust, into_py_variant, variant_repr};
use crate::lattice::impl_py_lattice;
use crate::num::{NumForm, NumLike};

#[pyclass]
pub enum RingScope {
    All(),
    Size(u8),
}

#[pymethods]
impl RingScope {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::All() => ("All", 0),
            Self::Size(_) => ("Size", 1),
        };
        variant_repr(slf.bind(py).as_any(), "RingScope", variant, arity)
    }
}

impl RingScope {
    pub(crate) fn from_rust(scope: &GraphIrRingScope) -> Self {
        match scope {
            GraphIrRingScope::All => Self::All(),
            GraphIrRingScope::Size(size) => Self::Size(*size),
        }
    }

    pub(crate) fn to_rust(&self) -> GraphIrRingScope {
        match self {
            Self::All() => GraphIrRingScope::All,
            Self::Size(size) => GraphIrRingScope::Size(*size),
        }
    }
}

#[pyclass]
pub struct RingMembershipForm {
    #[pyo3(get)]
    scope: Py<RingScope>,
    #[pyo3(get)]
    count: Py<NumForm>,
}

#[pymethods]
impl RingMembershipForm {
    #[new]
    fn new(py: Python<'_>, scope: Py<RingScope>, count: NumLike) -> PyResult<Self> {
        Ok(Self {
            scope,
            count: count.to_py(py)?,
        })
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "RingMembershipForm({}, {})",
            self.scope.bind(py).as_any().repr()?.extract::<String>()?,
            self.count.bind(py).as_any().repr()?.extract::<String>()?,
        ))
    }
}

impl_py_lattice!(
    RingMembershipForm,
    GraphIrRingMembershipForm,
    |value: &RingMembershipForm, py: Python<'_>| -> PyResult<GraphIrRingMembershipForm> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrRingMembershipForm| -> PyResult<RingMembershipForm> {
        RingMembershipForm::from_rust(py, &value)
    }
);

impl RingMembershipForm {
    pub(crate) fn from_rust(py: Python<'_>, form: &GraphIrRingMembershipForm) -> PyResult<Self> {
        Ok(Self {
            scope: into_py_variant(py, RingScope::from_rust(&form.scope))?,
            count: into_py_variant(py, NumForm::from_rust(py, &form.count)?)?,
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrRingMembershipForm {
        GraphIrRingMembershipForm::new(
            self.scope.bind(py).borrow().to_rust(),
            self.count.bind(py).borrow().to_rust(py),
        )
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(GraphIrRingScope::All)]
    #[case(GraphIrRingScope::Size(6))]
    fn test_ring_scope_roundtrip(#[case] form: GraphIrRingScope) {
        assert_eq!(RingScope::from_rust(&form).to_rust(), form);
    }

    #[rstest]
    #[case(GraphIrRingMembershipForm::new(GraphIrRingScope::All, 2))]
    #[case(GraphIrRingMembershipForm::new(GraphIrRingScope::Size(6), 1))]
    fn test_ring_membership_form_roundtrip(#[case] form: GraphIrRingMembershipForm) {
        Python::attach(|py| {
            assert_eq!(
                RingMembershipForm::from_rust(py, &form)
                    .unwrap()
                    .to_rust(py),
                form
            );
        });
    }
}
