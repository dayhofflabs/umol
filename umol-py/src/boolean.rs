//! Boolean AST mirror: `Undetermined` or a boolean literal. A leaf of the bond and
//! dative `Aromatic` constraint.

use pyo3::prelude::*;
use umol_ast::ast::BooleanAst as AstBooleanAst;

use crate::convert::{hash_ast, variant_repr};

/// A boolean expression: undetermined, or a literal `True` / `False`.
#[pyclass]
pub enum BooleanAst {
    Undetermined(),
    Lit(bool),
}

#[pymethods]
impl BooleanAst {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
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
    pub(crate) fn from_ast(ast: &AstBooleanAst) -> Self {
        match ast {
            AstBooleanAst::Undetermined => Self::Undetermined(),
            AstBooleanAst::Lit(b) => Self::Lit(*b),
        }
    }

    pub(crate) fn to_ast(&self) -> AstBooleanAst {
        match self {
            Self::Undetermined() => AstBooleanAst::Undetermined,
            Self::Lit(b) => AstBooleanAst::Lit(*b),
        }
    }
}

/// Setter coercion for a boolean field: a Python `bool` → `Lit`, or a `BooleanAst`
/// passthrough (mirroring `impl Into<BooleanAst>`).
#[derive(FromPyObject)]
pub(crate) enum BooleanArg {
    Lit(bool),
    Ast(Py<BooleanAst>),
}

impl BooleanArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstBooleanAst {
        match self {
            BooleanArg::Lit(b) => AstBooleanAst::Lit(*b),
            BooleanArg::Ast(a) => a.bind(py).borrow().to_ast(),
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
        assert_eq!(BooleanAst::from_ast(&ast).to_ast(), ast);
    }
}
