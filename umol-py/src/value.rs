//! `ValueTerm` — arithmetic term over integers, mirroring `umol_ast::ast::ValueTerm`
//! as a native PyO3 complex enum. AST recursion (`Box<Self>` / `Vec<Self>`) becomes
//! `Py<Self>` / `Vec<Py<Self>>`; per-variant construction and `match` work natively
//! on the Python side.

use std::collections::BTreeSet;

use pyo3::prelude::*;
use umol_ast::ast::{
    MemOp as AstMemOp, RelOp as AstRelOp, ValueAst as AstValueAst,
    ValuePredicate as AstValuePredicate, ValueTerm as AstValueTerm,
};

/// Relational operator in a value predicate (`<=`, `>=`, `==`, `<`, `>`, `!=`).
#[pyclass(eq, from_py_object)]
#[derive(Clone, PartialEq)]
pub enum RelOp {
    Le,
    Ge,
    Eq,
    Lt,
    Gt,
    Ne,
}

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
#[pyclass(eq, from_py_object)]
#[derive(Clone, PartialEq)]
pub enum MemOp {
    In,
    NotIn,
}

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

/// Boolean predicate over value terms.
#[pyclass]
pub enum ValuePredicate {
    Rel(Py<ValueTerm>, RelOp, Py<ValueTerm>),
    Mem(Py<ValueTerm>, MemOp, BTreeSet<i64>),
    Not(Py<ValuePredicate>),
    And(Vec<Py<ValuePredicate>>),
    Or(Vec<Py<ValuePredicate>>),
}

impl ValuePredicate {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstValuePredicate) -> PyResult<ValuePredicate> {
        Ok(match ast {
            AstValuePredicate::Rel(a, op, b) => ValuePredicate::Rel(
                Py::new(py, ValueTerm::from_ast(py, a)?)?,
                RelOp::from_ast(*op),
                Py::new(py, ValueTerm::from_ast(py, b)?)?,
            ),
            AstValuePredicate::Mem(t, op, members) => ValuePredicate::Mem(
                Py::new(py, ValueTerm::from_ast(py, t)?)?,
                MemOp::from_ast(*op),
                members.clone(),
            ),
            AstValuePredicate::Not(p) => ValuePredicate::Not(Py::new(py, Self::from_ast(py, p)?)?),
            AstValuePredicate::And(ps) => ValuePredicate::And(
                ps.iter()
                    .map(|p| Py::new(py, Self::from_ast(py, p)?))
                    .collect::<PyResult<_>>()?,
            ),
            AstValuePredicate::Or(ps) => ValuePredicate::Or(
                ps.iter()
                    .map(|p| Py::new(py, Self::from_ast(py, p)?))
                    .collect::<PyResult<_>>()?,
            ),
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstValuePredicate {
        match self {
            ValuePredicate::Rel(a, op, b) => AstValuePredicate::Rel(
                a.bind(py).borrow().to_ast(py),
                op.to_ast(),
                b.bind(py).borrow().to_ast(py),
            ),
            ValuePredicate::Mem(t, op, members) => AstValuePredicate::Mem(
                t.bind(py).borrow().to_ast(py),
                op.to_ast(),
                members.clone(),
            ),
            ValuePredicate::Not(p) => {
                AstValuePredicate::Not(Box::new(p.bind(py).borrow().to_ast(py)))
            }
            ValuePredicate::And(ps) => {
                AstValuePredicate::And(ps.iter().map(|p| p.bind(py).borrow().to_ast(py)).collect())
            }
            ValuePredicate::Or(ps) => {
                AstValuePredicate::Or(ps.iter().map(|p| p.bind(py).borrow().to_ast(py)).collect())
            }
        }
    }
}

/// Integer-valued atom/bond field: the undetermined wildcard, a literal, a literal
/// set, a range, an arithmetic term, or a boolean predicate.
#[pyclass]
pub enum ValueAst {
    Undetermined(),
    Lit(i64),
    LitSet(BTreeSet<i64>),
    RangeFrom(i64),
    RangeTo(i64),
    Term(Py<ValueTerm>),
    Predicate(Py<ValuePredicate>),
}

// The AST bridge; consumed by `AtomAst` at S3 (unused in the lib until then).
#[allow(dead_code)]
impl ValueAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstValueAst) -> PyResult<ValueAst> {
        Ok(match ast {
            AstValueAst::Undetermined => ValueAst::Undetermined(),
            AstValueAst::Lit(n) => ValueAst::Lit(*n),
            AstValueAst::LitSet(members) => ValueAst::LitSet((**members).clone()),
            AstValueAst::RangeFrom(n) => ValueAst::RangeFrom(*n),
            AstValueAst::RangeTo(n) => ValueAst::RangeTo(*n),
            AstValueAst::Term(t) => ValueAst::Term(Py::new(py, ValueTerm::from_ast(py, t)?)?),
            AstValueAst::Predicate(p) => {
                ValueAst::Predicate(Py::new(py, ValuePredicate::from_ast(py, p)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstValueAst {
        match self {
            ValueAst::Undetermined() => AstValueAst::Undetermined,
            ValueAst::Lit(n) => AstValueAst::Lit(*n),
            ValueAst::LitSet(members) => AstValueAst::LitSet(Box::new(members.clone())),
            ValueAst::RangeFrom(n) => AstValueAst::RangeFrom(*n),
            ValueAst::RangeTo(n) => AstValueAst::RangeTo(*n),
            ValueAst::Term(t) => AstValueAst::Term(Box::new(t.bind(py).borrow().to_ast(py))),
            ValueAst::Predicate(p) => {
                AstValueAst::Predicate(Box::new(p.bind(py).borrow().to_ast(py)))
            }
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

    #[rstest]
    #[case(AstValuePredicate::Rel(AstValueTerm::Var("h".into()), AstRelOp::Le, AstValueTerm::Lit(3)))]
    #[case(AstValuePredicate::Mem(AstValueTerm::Lit(0), AstMemOp::In, BTreeSet::from([1, 2, 3])))]
    #[case(AstValuePredicate::Not(Box::new(AstValuePredicate::Rel(
        AstValueTerm::Lit(1),
        AstRelOp::Eq,
        AstValueTerm::Lit(1),
    ))))]
    #[case(AstValuePredicate::And(vec![
        AstValuePredicate::Rel(AstValueTerm::Lit(1), AstRelOp::Lt, AstValueTerm::Lit(2)),
        AstValuePredicate::Rel(AstValueTerm::Lit(3), AstRelOp::Gt, AstValueTerm::Lit(2)),
    ]))]
    fn test_value_predicate_roundtrip(#[case] ast: AstValuePredicate) {
        Python::attach(|py| {
            let mirror = ValuePredicate::from_ast(py, &ast).unwrap();
            assert_eq!(mirror.to_ast(py), ast);
        });
    }

    #[rstest]
    #[case(AstValueAst::Undetermined)]
    #[case(AstValueAst::Lit(7))]
    #[case(AstValueAst::LitSet(Box::new(BTreeSet::from([1, 2, 3]))))]
    #[case(AstValueAst::RangeFrom(1))]
    #[case(AstValueAst::RangeTo(9))]
    #[case(AstValueAst::Term(Box::new(AstValueTerm::Var("x".into()))))]
    #[case(AstValueAst::Predicate(Box::new(AstValuePredicate::Rel(
        AstValueTerm::Var("h".into()),
        AstRelOp::Le,
        AstValueTerm::Lit(3),
    ))))]
    fn test_value_ast_roundtrip(#[case] ast: AstValueAst) {
        Python::attach(|py| {
            let mirror = ValueAst::from_ast(py, &ast).unwrap();
            assert_eq!(mirror.to_ast(py), ast);
        });
    }
}
