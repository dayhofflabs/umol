//! Boolean AST value: `Undetermined` or a boolean literal. A leaf of the bond and
//! dative `Aromatic` constraint.

use pyo3::prelude::*;
use umol_graph_ir::ir::BooleanForm as GraphIrBooleanForm;

use crate::convert::{hash_rust, variant_repr};
use crate::lattice::impl_py_lattice;

/// A boolean expression: undetermined, or a literal `True` / `False`.
#[pyclass]
pub enum BooleanForm {
    Undetermined(),
    Lit(bool),
}

#[pymethods]
impl BooleanForm {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            BooleanForm::Undetermined() => ("Undetermined", 0),
            BooleanForm::Lit(_) => ("Lit", 1),
        };
        variant_repr(slf.bind(py).as_any(), "BooleanForm", variant, arity)
    }
}

impl BooleanForm {
    pub(crate) fn from_rust(ast: &GraphIrBooleanForm) -> Self {
        match ast {
            GraphIrBooleanForm::Undetermined => Self::Undetermined(),
            GraphIrBooleanForm::Lit(b) => Self::Lit(*b),
        }
    }

    pub(crate) fn to_rust(&self) -> GraphIrBooleanForm {
        match self {
            Self::Undetermined() => GraphIrBooleanForm::Undetermined,
            Self::Lit(b) => GraphIrBooleanForm::Lit(*b),
        }
    }
}

impl_py_lattice!(
    BooleanForm,
    GraphIrBooleanForm,
    |value: &BooleanForm, _py: Python<'_>| -> PyResult<GraphIrBooleanForm> { Ok(value.to_rust()) },
    |_py: Python<'_>, value: GraphIrBooleanForm| -> PyResult<BooleanForm> {
        Ok(BooleanForm::from_rust(&value))
    }
);

/// Setter coercion for a boolean field: a Python `bool` → `Lit`, or a `BooleanForm`
/// passthrough (matching `impl Into<BooleanForm>`).
#[derive(FromPyObject)]
pub(crate) enum BooleanLike {
    Lit(bool),
    Ast(Py<BooleanForm>),
}

impl BooleanLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrBooleanForm {
        match self {
            BooleanLike::Lit(b) => GraphIrBooleanForm::Lit(*b),
            BooleanLike::Ast(a) => a.bind(py).borrow().to_rust(),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(GraphIrBooleanForm::Undetermined)]
    #[case(GraphIrBooleanForm::Lit(true))]
    #[case(GraphIrBooleanForm::Lit(false))]
    fn test_boolean_ast_roundtrip(#[case] ast: GraphIrBooleanForm) {
        assert_eq!(BooleanForm::from_rust(&ast).to_rust(), ast);
    }
}
