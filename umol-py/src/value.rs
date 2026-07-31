//! `ValueTerm` — arithmetic term over integers, matching `umol_ast::ast::ValueTerm`
//! as a native PyO3 complex enum. AST recursion (`Box<Self>` / `Vec<Self>`) becomes
//! `Py<Self>` / `Vec<Py<Self>>`; per-variant construction and `match` work natively
//! on the Python side.
// Blanket-allow the `absolute_paths` false positives from pyo3's `hash` derive
// (hygienic `::std::…` paths). Hand-written code here imports at top.
#![allow(clippy::absolute_paths)]

use std::collections::BTreeSet;

use pyo3::prelude::*;
use umol_ast::ast::{
    AsLit, MemOp as AstMemOp, RelOp as AstRelOp, ValueAst as AstValueAst,
    ValuePredicate as AstValuePredicate, ValueTerm as AstValueTerm,
};

use crate::convert::{hash_rust, into_py_variant, variant_repr};

/// Relational operator in a value predicate (`<=`, `>=`, `==`, `<`, `>`, `!=`).
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, PartialEq, Eq, Hash)]
pub enum RelOp {
    Le,
    Ge,
    Eq,
    Lt,
    Gt,
    Ne,
}

impl RelOp {
    pub(crate) fn from_rust(ast: AstRelOp) -> RelOp {
        match ast {
            AstRelOp::Le => RelOp::Le,
            AstRelOp::Ge => RelOp::Ge,
            AstRelOp::Eq => RelOp::Eq,
            AstRelOp::Lt => RelOp::Lt,
            AstRelOp::Gt => RelOp::Gt,
            AstRelOp::Ne => RelOp::Ne,
        }
    }

    pub(crate) fn to_rust(&self) -> AstRelOp {
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
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, PartialEq, Eq, Hash)]
pub enum MemOp {
    In,
    NotIn,
}

impl MemOp {
    pub(crate) fn from_rust(ast: AstMemOp) -> MemOp {
        match ast {
            AstMemOp::In => MemOp::In,
            AstMemOp::NotIn => MemOp::NotIn,
        }
    }

    pub(crate) fn to_rust(&self) -> AstMemOp {
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

#[pymethods]
impl ValueTerm {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            ValueTerm::Lit(_) => ("Lit", 1),
            ValueTerm::Var(_) => ("Var", 1),
            ValueTerm::Neg(_) => ("Neg", 1),
            ValueTerm::Sum(_) => ("Sum", 1),
            ValueTerm::Product(_) => ("Product", 1),
            ValueTerm::Div(_, _) => ("Div", 2),
            ValueTerm::Rem(_, _) => ("Rem", 2),
        };
        variant_repr(slf.bind(py).as_any(), "ValueTerm", variant, arity)
    }
}

impl ValueTerm {
    /// Build the Python value from the AST term (one Python object per node).
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstValueTerm) -> PyResult<ValueTerm> {
        Ok(match ast {
            AstValueTerm::Lit(n) => ValueTerm::Lit(*n),
            AstValueTerm::Var(name) => ValueTerm::Var(name.clone()),
            AstValueTerm::Neg(t) => ValueTerm::Neg(into_py_variant(py, Self::from_rust(py, t)?)?),
            AstValueTerm::Sum(terms) => ValueTerm::Sum(
                terms
                    .iter()
                    .map(|t| into_py_variant(py, Self::from_rust(py, t)?))
                    .collect::<PyResult<_>>()?,
            ),
            AstValueTerm::Product(terms) => ValueTerm::Product(
                terms
                    .iter()
                    .map(|t| into_py_variant(py, Self::from_rust(py, t)?))
                    .collect::<PyResult<_>>()?,
            ),
            AstValueTerm::Div(a, b) => ValueTerm::Div(
                into_py_variant(py, Self::from_rust(py, a)?)?,
                into_py_variant(py, Self::from_rust(py, b)?)?,
            ),
            AstValueTerm::Rem(a, b) => ValueTerm::Rem(
                into_py_variant(py, Self::from_rust(py, a)?)?,
                into_py_variant(py, Self::from_rust(py, b)?)?,
            ),
        })
    }

    /// Lower the Python value back to the AST term.
    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstValueTerm {
        match self {
            ValueTerm::Lit(n) => AstValueTerm::Lit(*n),
            ValueTerm::Var(name) => AstValueTerm::Var(name.clone()),
            ValueTerm::Neg(t) => AstValueTerm::Neg(Box::new(t.bind(py).borrow().to_rust(py))),
            ValueTerm::Sum(terms) => AstValueTerm::Sum(
                terms
                    .iter()
                    .map(|t| t.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
            ValueTerm::Product(terms) => AstValueTerm::Product(
                terms
                    .iter()
                    .map(|t| t.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
            ValueTerm::Div(a, b) => AstValueTerm::Div(
                Box::new(a.bind(py).borrow().to_rust(py)),
                Box::new(b.bind(py).borrow().to_rust(py)),
            ),
            ValueTerm::Rem(a, b) => AstValueTerm::Rem(
                Box::new(a.bind(py).borrow().to_rust(py)),
                Box::new(b.bind(py).borrow().to_rust(py)),
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

#[pymethods]
impl ValuePredicate {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            ValuePredicate::Rel(_, _, _) => ("Rel", 3),
            ValuePredicate::Mem(_, _, _) => ("Mem", 3),
            ValuePredicate::Not(_) => ("Not", 1),
            ValuePredicate::And(_) => ("And", 1),
            ValuePredicate::Or(_) => ("Or", 1),
        };
        variant_repr(slf.bind(py).as_any(), "ValuePredicate", variant, arity)
    }
}

impl ValuePredicate {
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstValuePredicate) -> PyResult<ValuePredicate> {
        Ok(match ast {
            AstValuePredicate::Rel(a, op, b) => ValuePredicate::Rel(
                into_py_variant(py, ValueTerm::from_rust(py, a)?)?,
                RelOp::from_rust(*op),
                into_py_variant(py, ValueTerm::from_rust(py, b)?)?,
            ),
            AstValuePredicate::Mem(t, op, members) => ValuePredicate::Mem(
                into_py_variant(py, ValueTerm::from_rust(py, t)?)?,
                MemOp::from_rust(*op),
                members.clone(),
            ),
            AstValuePredicate::Not(p) => {
                ValuePredicate::Not(into_py_variant(py, Self::from_rust(py, p)?)?)
            }
            AstValuePredicate::And(ps) => ValuePredicate::And(
                ps.iter()
                    .map(|p| into_py_variant(py, Self::from_rust(py, p)?))
                    .collect::<PyResult<_>>()?,
            ),
            AstValuePredicate::Or(ps) => ValuePredicate::Or(
                ps.iter()
                    .map(|p| into_py_variant(py, Self::from_rust(py, p)?))
                    .collect::<PyResult<_>>()?,
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstValuePredicate {
        match self {
            ValuePredicate::Rel(a, op, b) => AstValuePredicate::Rel(
                a.bind(py).borrow().to_rust(py),
                op.to_rust(),
                b.bind(py).borrow().to_rust(py),
            ),
            ValuePredicate::Mem(t, op, members) => AstValuePredicate::Mem(
                t.bind(py).borrow().to_rust(py),
                op.to_rust(),
                members.clone(),
            ),
            ValuePredicate::Not(p) => {
                AstValuePredicate::Not(Box::new(p.bind(py).borrow().to_rust(py)))
            }
            ValuePredicate::And(ps) => {
                AstValuePredicate::And(ps.iter().map(|p| p.bind(py).borrow().to_rust(py)).collect())
            }
            ValuePredicate::Or(ps) => {
                AstValuePredicate::Or(ps.iter().map(|p| p.bind(py).borrow().to_rust(py)).collect())
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

#[pymethods]
impl ValueAst {
    /// The concrete integer this resolves to, or `None` when it is not a bare
    /// literal (undetermined, a set, a range, or an expression).
    fn as_lit(&self, py: Python<'_>) -> Option<i64> {
        self.to_rust(py).as_lit()
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            ValueAst::Undetermined() => ("Undetermined", 0),
            ValueAst::Lit(_) => ("Lit", 1),
            ValueAst::LitSet(_) => ("LitSet", 1),
            ValueAst::RangeFrom(_) => ("RangeFrom", 1),
            ValueAst::RangeTo(_) => ("RangeTo", 1),
            ValueAst::Term(_) => ("Term", 1),
            ValueAst::Predicate(_) => ("Predicate", 1),
        };
        variant_repr(slf.bind(py).as_any(), "ValueAst", variant, arity)
    }
}

impl ValueAst {
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstValueAst) -> PyResult<ValueAst> {
        Ok(match ast {
            AstValueAst::Undetermined => ValueAst::Undetermined(),
            AstValueAst::Lit(n) => ValueAst::Lit(*n),
            AstValueAst::LitSet(members) => ValueAst::LitSet((**members).clone()),
            AstValueAst::RangeFrom(n) => ValueAst::RangeFrom(*n),
            AstValueAst::RangeTo(n) => ValueAst::RangeTo(*n),
            AstValueAst::Term(t) => {
                ValueAst::Term(into_py_variant(py, ValueTerm::from_rust(py, t)?)?)
            }
            AstValueAst::Predicate(p) => {
                ValueAst::Predicate(into_py_variant(py, ValuePredicate::from_rust(py, p)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstValueAst {
        match self {
            ValueAst::Undetermined() => AstValueAst::Undetermined,
            ValueAst::Lit(n) => AstValueAst::Lit(*n),
            ValueAst::LitSet(members) => AstValueAst::LitSet(Box::new(members.clone())),
            ValueAst::RangeFrom(n) => AstValueAst::RangeFrom(*n),
            ValueAst::RangeTo(n) => AstValueAst::RangeTo(*n),
            ValueAst::Term(t) => AstValueAst::Term(Box::new(t.bind(py).borrow().to_rust(py))),
            ValueAst::Predicate(p) => {
                AstValueAst::Predicate(Box::new(p.bind(py).borrow().to_rust(py)))
            }
        }
    }
}

/// A `ValueAst` or a Python `int` (→ `ValueAst::Lit`), matching `impl Into<ValueAst>`
/// on the Rust builders. The `*Arg` convention for binding coercion inputs (`*Input`
/// is the DSL side); shared by the atom fields, unpaired-electron components, and
/// ring-membership count.
#[derive(FromPyObject)]
pub enum ValueArg {
    Ast(Py<ValueAst>),
    Lit(i64),
}

impl ValueArg {
    /// Coerce to the value AST (for `impl Into<ValueAst>` Rust builders).
    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstValueAst {
        match self {
            ValueArg::Ast(value) => value.bind(py).borrow().to_rust(py),
            ValueArg::Lit(number) => AstValueAst::Lit(*number),
        }
    }

    /// Coerce to a `Py<ValueAst>` (for value structs that store the value field).
    pub(crate) fn to_py(&self, py: Python<'_>) -> PyResult<Py<ValueAst>> {
        match self {
            ValueArg::Ast(value) => Ok(value.clone_ref(py)),
            ValueArg::Lit(number) => into_py_variant(py, ValueAst::Lit(*number)),
        }
    }
}

/// `IntoPyObject` for `&ValueArg` so it can be a complex-enum field: constructors
/// (`AromaticValenceAst.Aromatic(1)`) coerce `int | ValueAst` in, and the field
/// reads back as a `ValueAst`.
impl<'py> IntoPyObject<'py> for &ValueArg {
    type Target = ValueAst;
    type Output = Bound<'py, ValueAst>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Bound<'py, ValueAst>> {
        Ok(self.to_py(py)?.into_bound(py))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(AstRelOp::Le)]
    #[case(AstRelOp::Ge)]
    #[case(AstRelOp::Eq)]
    #[case(AstRelOp::Lt)]
    #[case(AstRelOp::Gt)]
    #[case(AstRelOp::Ne)]
    fn test_rel_op_roundtrip(#[case] ast: AstRelOp) {
        assert_eq!(RelOp::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case(AstMemOp::In)]
    #[case(AstMemOp::NotIn)]
    fn test_mem_op_roundtrip(#[case] ast: AstMemOp) {
        assert_eq!(MemOp::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case(AstValueTerm::Lit(5))]
    #[case(AstValueTerm::Var("x".to_string()))]
    #[case(AstValueTerm::Neg(Box::new(AstValueTerm::Lit(3))))]
    #[case(AstValueTerm::Sum(vec![AstValueTerm::Lit(1), AstValueTerm::Lit(2)]))]
    #[case(AstValueTerm::Div(Box::new(AstValueTerm::Lit(6)), Box::new(AstValueTerm::Lit(2))))]
    fn test_value_term_roundtrip(#[case] ast: AstValueTerm) {
        Python::attach(|py| {
            let value = ValueTerm::from_rust(py, &ast).unwrap();
            assert_eq!(value.to_rust(py), ast);
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
            let value = ValuePredicate::from_rust(py, &ast).unwrap();
            assert_eq!(value.to_rust(py), ast);
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
            let value = ValueAst::from_rust(py, &ast).unwrap();
            assert_eq!(value.to_rust(py), ast);
        });
    }

    #[rstest]
    #[case(AstValueAst::Lit(4), Some(4))]
    #[case(AstValueAst::Lit(-1), Some(-1))]
    #[case(AstValueAst::Undetermined, None)]
    #[case(AstValueAst::RangeFrom(1), None)]
    #[case(AstValueAst::LitSet(Box::new(BTreeSet::from([1, 2]))), None)]
    fn test_value_ast_as_lit(#[case] ast: AstValueAst, #[case] expected: Option<i64>) {
        Python::attach(|py| {
            assert_eq!(ValueAst::from_rust(py, &ast).unwrap().as_lit(py), expected);
        });
    }
}
