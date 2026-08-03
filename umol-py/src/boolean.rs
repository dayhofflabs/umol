//! Boolean AST value: `Undetermined` or a boolean literal. A leaf of the bond and
//! dative `Aromatic` constraint.

use pyo3::prelude::*;
use umol_ast::ast::BooleanAst as AstBooleanAst;

use crate::convert::{hash_rust, variant_repr};
use crate::lattice::impl_py_lattice;

/// A boolean expression: undetermined, or a literal `True` / `False`.
#[pyclass]
pub enum BooleanAst {
    Undetermined(),
    Lit(bool),
}

#[pymethods]
impl BooleanAst {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            BooleanAst::Undetermined() => ("Undetermined", 0),
            BooleanAst::Lit(_) => ("Lit", 1),
        };
        variant_repr(slf.bind(py).as_any(), "BooleanAst", variant, arity)
    }
}

impl BooleanAst {
    pub(crate) fn from_rust(ast: &AstBooleanAst) -> Self {
        match ast {
            AstBooleanAst::Undetermined => Self::Undetermined(),
            AstBooleanAst::Lit(b) => Self::Lit(*b),
        }
    }

    pub(crate) fn to_rust(&self) -> AstBooleanAst {
        match self {
            Self::Undetermined() => AstBooleanAst::Undetermined,
            Self::Lit(b) => AstBooleanAst::Lit(*b),
        }
    }
}

impl_py_lattice!(
    BooleanAst,
    AstBooleanAst,
    |value: &BooleanAst, _py: Python<'_>| -> PyResult<AstBooleanAst> { Ok(value.to_rust()) },
    |_py: Python<'_>, value: AstBooleanAst| -> PyResult<BooleanAst> {
        Ok(BooleanAst::from_rust(&value))
    }
);

/// Setter coercion for a boolean field: a Python `bool` → `Lit`, or a `BooleanAst`
/// passthrough (matching `impl Into<BooleanAst>`).
#[derive(FromPyObject)]
pub(crate) enum BooleanLike {
    Lit(bool),
    Ast(Py<BooleanAst>),
}

impl BooleanLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstBooleanAst {
        match self {
            BooleanLike::Lit(b) => AstBooleanAst::Lit(*b),
            BooleanLike::Ast(a) => a.bind(py).borrow().to_rust(),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(AstBooleanAst::Undetermined)]
    #[case(AstBooleanAst::Lit(true))]
    #[case(AstBooleanAst::Lit(false))]
    fn test_boolean_ast_roundtrip(#[case] ast: AstBooleanAst) {
        assert_eq!(BooleanAst::from_rust(&ast).to_rust(), ast);
    }
}
