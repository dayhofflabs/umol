//! `ValueTerm` — arithmetic term over integers, mirroring `umol_ast::ast::ValueTerm`
//! as a native PyO3 complex enum. AST recursion (`Box<Self>` / `Vec<Self>`) becomes
//! `Py<Self>` / `Vec<Py<Self>>`; per-variant construction and `match` work natively
//! on the Python side.

use pyo3::prelude::*;
use umol_ast::ast::{MemOp as AstMemOp, RelOp as AstRelOp, ValueTerm as AstValueTerm};

/// Relational operator in a value predicate (`<=`, `>=`, `==`, `<`, `>`, `!=`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub enum RelOp {
    Le,
    Ge,
    Eq,
    Lt,
    Gt,
    Ne,
}

// The AST bridge; consumed by `ValuePredicate` at S1c (unused in the lib until then).
#[allow(dead_code)]
impl RelOp {
    pub(crate) fn from_ast(ast: AstRelOp) -> RelOp {
        match ast {
            AstRelOp::Le => RelOp::Le,
            AstRelOp::Ge => RelOp::Ge,
            AstRelOp::Eq => RelOp::Eq,
            AstRelOp::Lt => RelOp::Lt,
            AstRelOp::Gt => RelOp::Gt,
            AstRelOp::Ne => RelOp::Ne,
        }
    }

    pub(crate) fn to_ast(&self) -> AstRelOp {
        match self {
            RelOp::Le => AstRelOp::Le,
            RelOp::Ge => AstRelOp::Ge,
            RelOp::Eq => AstRelOp::Eq,
            RelOp::Lt => AstRelOp::Lt,
            RelOp::Gt => AstRelOp::Gt,
            RelOp::Ne => AstRelOp::Ne,
        }
    }
}

/// Membership operator in a value predicate (`in`, `not in`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub enum MemOp {
    In,
    NotIn,
}

// The AST bridge; consumed by `ValuePredicate` at S1c (unused in the lib until then).
#[allow(dead_code)]
impl MemOp {
    pub(crate) fn from_ast(ast: AstMemOp) -> MemOp {
        match ast {
            AstMemOp::In => MemOp::In,
            AstMemOp::NotIn => MemOp::NotIn,
        }
    }

    pub(crate) fn to_ast(&self) -> AstMemOp {
        match self {
            MemOp::In => AstMemOp::In,
            MemOp::NotIn => AstMemOp::NotIn,
        }
    }
}

/// Arithmetic term over integers.
#[pyclass]
pub enum ValueTerm {
    Lit(i64),
    Var(String),
    Neg(Py<ValueTerm>),
    Sum(Vec<Py<ValueTerm>>),
    Product(Vec<Py<ValueTerm>>),
    Div(Py<ValueTerm>, Py<ValueTerm>),
    Rem(Py<ValueTerm>, Py<ValueTerm>),
}

// The AST bridge; consumed by `ValueAst` at S1d (unused in the lib until then).
#[allow(dead_code)]
impl ValueTerm {
    /// Build the Python mirror from the AST term (one Python object per node).
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstValueTerm) -> PyResult<ValueTerm> {
        Ok(match ast {
            AstValueTerm::Lit(n) => ValueTerm::Lit(*n),
            AstValueTerm::Var(name) => ValueTerm::Var(name.clone()),
            AstValueTerm::Neg(t) => ValueTerm::Neg(Py::new(py, Self::from_ast(py, t)?)?),
            AstValueTerm::Sum(terms) => ValueTerm::Sum(
                terms
                    .iter()
                    .map(|t| Py::new(py, Self::from_ast(py, t)?))
                    .collect::<PyResult<_>>()?,
            ),
            AstValueTerm::Product(terms) => ValueTerm::Product(
                terms
                    .iter()
                    .map(|t| Py::new(py, Self::from_ast(py, t)?))
                    .collect::<PyResult<_>>()?,
            ),
            AstValueTerm::Div(a, b) => ValueTerm::Div(
                Py::new(py, Self::from_ast(py, a)?)?,
                Py::new(py, Self::from_ast(py, b)?)?,
            ),
            AstValueTerm::Rem(a, b) => ValueTerm::Rem(
                Py::new(py, Self::from_ast(py, a)?)?,
                Py::new(py, Self::from_ast(py, b)?)?,
            ),
        })
    }

    /// Lower the Python mirror back to the AST term.
    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstValueTerm {
        match self {
            ValueTerm::Lit(n) => AstValueTerm::Lit(*n),
            ValueTerm::Var(name) => AstValueTerm::Var(name.clone()),
            ValueTerm::Neg(t) => AstValueTerm::Neg(Box::new(t.bind(py).borrow().to_ast(py))),
            ValueTerm::Sum(terms) => {
                AstValueTerm::Sum(terms.iter().map(|t| t.bind(py).borrow().to_ast(py)).collect())
            }
            ValueTerm::Product(terms) => {
                AstValueTerm::Product(terms.iter().map(|t| t.bind(py).borrow().to_ast(py)).collect())
            }
            ValueTerm::Div(a, b) => AstValueTerm::Div(
                Box::new(a.bind(py).borrow().to_ast(py)),
                Box::new(b.bind(py).borrow().to_ast(py)),
            ),
            ValueTerm::Rem(a, b) => AstValueTerm::Rem(
                Box::new(a.bind(py).borrow().to_ast(py)),
                Box::new(b.bind(py).borrow().to_ast(py)),
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(AstRelOp::Le)]
    #[case(AstRelOp::Ge)]
    #[case(AstRelOp::Eq)]
    #[case(AstRelOp::Lt)]
    #[case(AstRelOp::Gt)]
    #[case(AstRelOp::Ne)]
    fn test_rel_op_roundtrip(#[case] ast: AstRelOp) {
        assert_eq!(RelOp::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstMemOp::In)]
    #[case(AstMemOp::NotIn)]
    fn test_mem_op_roundtrip(#[case] ast: AstMemOp) {
        assert_eq!(MemOp::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstValueTerm::Lit(5))]
    #[case(AstValueTerm::Var("x".to_string()))]
    #[case(AstValueTerm::Neg(Box::new(AstValueTerm::Lit(3))))]
    #[case(AstValueTerm::Sum(vec![AstValueTerm::Lit(1), AstValueTerm::Lit(2)]))]
    #[case(AstValueTerm::Div(Box::new(AstValueTerm::Lit(6)), Box::new(AstValueTerm::Lit(2))))]
    fn test_value_term_roundtrip(#[case] ast: AstValueTerm) {
        Python::attach(|py| {
            let mirror = ValueTerm::from_ast(py, &ast).unwrap();
            assert_eq!(mirror.to_ast(py), ast);
        });
    }
}
