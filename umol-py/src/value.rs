//! `ValueTerm` — arithmetic term over integers, matching `umol_graph_ir::ir::ValueTerm`
//! as a native PyO3 complex enum. AST recursion (`Box<Self>` / `Vec<Self>`) becomes
//! `Py<Self>` / `Vec<Py<Self>>`; per-variant construction and `match` work natively
//! on the Python side.
// Blanket-allow the `absolute_paths` false positives from pyo3's `hash` derive
// (hygienic `::std::…` paths). Hand-written code here imports at top.
#![allow(clippy::absolute_paths)]

use std::collections::BTreeSet;

use pyo3::prelude::*;
use umol_graph_ir::ir::{
    AsLit, MemOp as GraphIrMemOp, RelOp as GraphIrRelOp, ValueAst as GraphIrValueAst,
    ValuePredicate as GraphIrValuePredicate, ValueTerm as GraphIrValueTerm,
};

use crate::convert::{hash_rust, into_py_variant, variant_repr};
use crate::lattice::{impl_py_canonicalize, impl_py_lattice};

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
    pub(crate) fn from_rust(ast: GraphIrRelOp) -> RelOp {
        match ast {
            GraphIrRelOp::Le => RelOp::Le,
            GraphIrRelOp::Ge => RelOp::Ge,
            GraphIrRelOp::Eq => RelOp::Eq,
            GraphIrRelOp::Lt => RelOp::Lt,
            GraphIrRelOp::Gt => RelOp::Gt,
            GraphIrRelOp::Ne => RelOp::Ne,
        }
    }

    pub(crate) fn to_rust(&self) -> GraphIrRelOp {
        match self {
            RelOp::Le => GraphIrRelOp::Le,
            RelOp::Ge => GraphIrRelOp::Ge,
            RelOp::Eq => GraphIrRelOp::Eq,
            RelOp::Lt => GraphIrRelOp::Lt,
            RelOp::Gt => GraphIrRelOp::Gt,
            RelOp::Ne => GraphIrRelOp::Ne,
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
    pub(crate) fn from_rust(ast: GraphIrMemOp) -> MemOp {
        match ast {
            GraphIrMemOp::In => MemOp::In,
            GraphIrMemOp::NotIn => MemOp::NotIn,
        }
    }

    pub(crate) fn to_rust(&self) -> GraphIrMemOp {
        match self {
            MemOp::In => GraphIrMemOp::In,
            MemOp::NotIn => GraphIrMemOp::NotIn,
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
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrValueTerm) -> PyResult<ValueTerm> {
        Ok(match ast {
            GraphIrValueTerm::Lit(n) => ValueTerm::Lit(*n),
            GraphIrValueTerm::Var(name) => ValueTerm::Var(name.clone()),
            GraphIrValueTerm::Neg(t) => {
                ValueTerm::Neg(into_py_variant(py, Self::from_rust(py, t)?)?)
            }
            GraphIrValueTerm::Sum(terms) => ValueTerm::Sum(
                terms
                    .iter()
                    .map(|t| into_py_variant(py, Self::from_rust(py, t)?))
                    .collect::<PyResult<_>>()?,
            ),
            GraphIrValueTerm::Product(terms) => ValueTerm::Product(
                terms
                    .iter()
                    .map(|t| into_py_variant(py, Self::from_rust(py, t)?))
                    .collect::<PyResult<_>>()?,
            ),
            GraphIrValueTerm::Div(a, b) => ValueTerm::Div(
                into_py_variant(py, Self::from_rust(py, a)?)?,
                into_py_variant(py, Self::from_rust(py, b)?)?,
            ),
            GraphIrValueTerm::Rem(a, b) => ValueTerm::Rem(
                into_py_variant(py, Self::from_rust(py, a)?)?,
                into_py_variant(py, Self::from_rust(py, b)?)?,
            ),
        })
    }

    /// Lower the Python value back to the AST term.
    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrValueTerm {
        match self {
            ValueTerm::Lit(n) => GraphIrValueTerm::Lit(*n),
            ValueTerm::Var(name) => GraphIrValueTerm::Var(name.clone()),
            ValueTerm::Neg(t) => GraphIrValueTerm::Neg(Box::new(t.bind(py).borrow().to_rust(py))),
            ValueTerm::Sum(terms) => GraphIrValueTerm::Sum(
                terms
                    .iter()
                    .map(|t| t.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
            ValueTerm::Product(terms) => GraphIrValueTerm::Product(
                terms
                    .iter()
                    .map(|t| t.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
            ValueTerm::Div(a, b) => GraphIrValueTerm::Div(
                Box::new(a.bind(py).borrow().to_rust(py)),
                Box::new(b.bind(py).borrow().to_rust(py)),
            ),
            ValueTerm::Rem(a, b) => GraphIrValueTerm::Rem(
                Box::new(a.bind(py).borrow().to_rust(py)),
                Box::new(b.bind(py).borrow().to_rust(py)),
            ),
        }
    }
}

impl_py_canonicalize!(
    ValueTerm,
    GraphIrValueTerm,
    |value: &ValueTerm, py: Python<'_>| -> PyResult<GraphIrValueTerm> { Ok(value.to_rust(py)) },
    |py: Python<'_>, value: GraphIrValueTerm| -> PyResult<ValueTerm> {
        ValueTerm::from_rust(py, &value)
    }
);

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
    pub(crate) fn from_rust(
        py: Python<'_>,
        ast: &GraphIrValuePredicate,
    ) -> PyResult<ValuePredicate> {
        Ok(match ast {
            GraphIrValuePredicate::Rel(a, op, b) => ValuePredicate::Rel(
                into_py_variant(py, ValueTerm::from_rust(py, a)?)?,
                RelOp::from_rust(*op),
                into_py_variant(py, ValueTerm::from_rust(py, b)?)?,
            ),
            GraphIrValuePredicate::Mem(t, op, members) => ValuePredicate::Mem(
                into_py_variant(py, ValueTerm::from_rust(py, t)?)?,
                MemOp::from_rust(*op),
                members.clone(),
            ),
            GraphIrValuePredicate::Not(p) => {
                ValuePredicate::Not(into_py_variant(py, Self::from_rust(py, p)?)?)
            }
            GraphIrValuePredicate::And(ps) => ValuePredicate::And(
                ps.iter()
                    .map(|p| into_py_variant(py, Self::from_rust(py, p)?))
                    .collect::<PyResult<_>>()?,
            ),
            GraphIrValuePredicate::Or(ps) => ValuePredicate::Or(
                ps.iter()
                    .map(|p| into_py_variant(py, Self::from_rust(py, p)?))
                    .collect::<PyResult<_>>()?,
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrValuePredicate {
        match self {
            ValuePredicate::Rel(a, op, b) => GraphIrValuePredicate::Rel(
                a.bind(py).borrow().to_rust(py),
                op.to_rust(),
                b.bind(py).borrow().to_rust(py),
            ),
            ValuePredicate::Mem(t, op, members) => GraphIrValuePredicate::Mem(
                t.bind(py).borrow().to_rust(py),
                op.to_rust(),
                members.clone(),
            ),
            ValuePredicate::Not(p) => {
                GraphIrValuePredicate::Not(Box::new(p.bind(py).borrow().to_rust(py)))
            }
            ValuePredicate::And(ps) => GraphIrValuePredicate::And(
                ps.iter().map(|p| p.bind(py).borrow().to_rust(py)).collect(),
            ),
            ValuePredicate::Or(ps) => GraphIrValuePredicate::Or(
                ps.iter().map(|p| p.bind(py).borrow().to_rust(py)).collect(),
            ),
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
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrValueAst) -> PyResult<ValueAst> {
        Ok(match ast {
            GraphIrValueAst::Undetermined => ValueAst::Undetermined(),
            GraphIrValueAst::Lit(n) => ValueAst::Lit(*n),
            GraphIrValueAst::LitSet(members) => ValueAst::LitSet((**members).clone()),
            GraphIrValueAst::RangeFrom(n) => ValueAst::RangeFrom(*n),
            GraphIrValueAst::RangeTo(n) => ValueAst::RangeTo(*n),
            GraphIrValueAst::Term(t) => {
                ValueAst::Term(into_py_variant(py, ValueTerm::from_rust(py, t)?)?)
            }
            GraphIrValueAst::Predicate(p) => {
                ValueAst::Predicate(into_py_variant(py, ValuePredicate::from_rust(py, p)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrValueAst {
        match self {
            ValueAst::Undetermined() => GraphIrValueAst::Undetermined,
            ValueAst::Lit(n) => GraphIrValueAst::Lit(*n),
            ValueAst::LitSet(members) => GraphIrValueAst::LitSet(Box::new(members.clone())),
            ValueAst::RangeFrom(n) => GraphIrValueAst::RangeFrom(*n),
            ValueAst::RangeTo(n) => GraphIrValueAst::RangeTo(*n),
            ValueAst::Term(t) => GraphIrValueAst::Term(Box::new(t.bind(py).borrow().to_rust(py))),
            ValueAst::Predicate(p) => {
                GraphIrValueAst::Predicate(Box::new(p.bind(py).borrow().to_rust(py)))
            }
        }
    }
}

impl_py_lattice!(
    ValueAst,
    GraphIrValueAst,
    |value: &ValueAst, py: Python<'_>| -> PyResult<GraphIrValueAst> { Ok(value.to_rust(py)) },
    |py: Python<'_>, value: GraphIrValueAst| -> PyResult<ValueAst> {
        ValueAst::from_rust(py, &value)
    }
);

/// A `ValueAst` or a Python `int` (→ `ValueAst::Lit`), matching `impl Into<ValueAst>`
/// on the Rust builders. The `*Like` convention for binding coercion inputs (`*Input`
/// is the DSL side); shared by the atom fields, unpaired-electron components, and
/// ring-membership count.
#[derive(FromPyObject)]
pub enum ValueLike {
    Ast(Py<ValueAst>),
    Lit(i64),
}

impl ValueLike {
    /// Coerce to the value AST (for `impl Into<ValueAst>` Rust builders).
    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrValueAst {
        match self {
            ValueLike::Ast(value) => value.bind(py).borrow().to_rust(py),
            ValueLike::Lit(number) => GraphIrValueAst::Lit(*number),
        }
    }

    /// Coerce to a `Py<ValueAst>` (for value structs that store the value field).
    pub(crate) fn to_py(&self, py: Python<'_>) -> PyResult<Py<ValueAst>> {
        match self {
            ValueLike::Ast(value) => Ok(value.clone_ref(py)),
            ValueLike::Lit(number) => into_py_variant(py, ValueAst::Lit(*number)),
        }
    }
}

/// `IntoPyObject` for `&ValueLike` so it can be a complex-enum field: constructors
/// (`AromaticValenceAst.Aromatic(1)`) coerce `int | ValueAst` in, and the field
/// reads back as a `ValueAst`.
impl<'py> IntoPyObject<'py> for &ValueLike {
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
    #[case(GraphIrRelOp::Le)]
    #[case(GraphIrRelOp::Ge)]
    #[case(GraphIrRelOp::Eq)]
    #[case(GraphIrRelOp::Lt)]
    #[case(GraphIrRelOp::Gt)]
    #[case(GraphIrRelOp::Ne)]
    fn test_rel_op_roundtrip(#[case] ast: GraphIrRelOp) {
        assert_eq!(RelOp::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case(GraphIrMemOp::In)]
    #[case(GraphIrMemOp::NotIn)]
    fn test_mem_op_roundtrip(#[case] ast: GraphIrMemOp) {
        assert_eq!(MemOp::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case(GraphIrValueTerm::Lit(5))]
    #[case(GraphIrValueTerm::Var("x".to_string()))]
    #[case(GraphIrValueTerm::Neg(Box::new(GraphIrValueTerm::Lit(3))))]
    #[case(GraphIrValueTerm::Sum(vec![GraphIrValueTerm::Lit(1), GraphIrValueTerm::Lit(2)]))]
    #[case(GraphIrValueTerm::Div(
        Box::new(GraphIrValueTerm::Lit(6)),
        Box::new(GraphIrValueTerm::Lit(2))
    ))]
    fn test_value_term_roundtrip(#[case] ast: GraphIrValueTerm) {
        Python::attach(|py| {
            let value = ValueTerm::from_rust(py, &ast).unwrap();
            assert_eq!(value.to_rust(py), ast);
        });
    }

    #[rstest]
    #[case(GraphIrValuePredicate::Rel(GraphIrValueTerm::Var("h".into()), GraphIrRelOp::Le, GraphIrValueTerm::Lit(3)))]
    #[case(GraphIrValuePredicate::Mem(GraphIrValueTerm::Lit(0), GraphIrMemOp::In, BTreeSet::from([1, 2, 3])))]
    #[case(GraphIrValuePredicate::Not(Box::new(GraphIrValuePredicate::Rel(
        GraphIrValueTerm::Lit(1),
        GraphIrRelOp::Eq,
        GraphIrValueTerm::Lit(1),
    ))))]
    #[case(GraphIrValuePredicate::And(vec![
        GraphIrValuePredicate::Rel(GraphIrValueTerm::Lit(1), GraphIrRelOp::Lt, GraphIrValueTerm::Lit(2)),
        GraphIrValuePredicate::Rel(GraphIrValueTerm::Lit(3), GraphIrRelOp::Gt, GraphIrValueTerm::Lit(2)),
    ]))]
    fn test_value_predicate_roundtrip(#[case] ast: GraphIrValuePredicate) {
        Python::attach(|py| {
            let value = ValuePredicate::from_rust(py, &ast).unwrap();
            assert_eq!(value.to_rust(py), ast);
        });
    }

    #[rstest]
    #[case(GraphIrValueAst::Undetermined)]
    #[case(GraphIrValueAst::Lit(7))]
    #[case(GraphIrValueAst::LitSet(Box::new(BTreeSet::from([1, 2, 3]))))]
    #[case(GraphIrValueAst::RangeFrom(1))]
    #[case(GraphIrValueAst::RangeTo(9))]
    #[case(GraphIrValueAst::Term(Box::new(GraphIrValueTerm::Var("x".into()))))]
    #[case(GraphIrValueAst::Predicate(Box::new(GraphIrValuePredicate::Rel(
        GraphIrValueTerm::Var("h".into()),
        GraphIrRelOp::Le,
        GraphIrValueTerm::Lit(3),
    ))))]
    fn test_value_ast_roundtrip(#[case] ast: GraphIrValueAst) {
        Python::attach(|py| {
            let value = ValueAst::from_rust(py, &ast).unwrap();
            assert_eq!(value.to_rust(py), ast);
        });
    }

    #[rstest]
    #[case(GraphIrValueAst::Lit(4), Some(4))]
    #[case(GraphIrValueAst::Lit(-1), Some(-1))]
    #[case(GraphIrValueAst::Undetermined, None)]
    #[case(GraphIrValueAst::RangeFrom(1), None)]
    #[case(GraphIrValueAst::LitSet(Box::new(BTreeSet::from([1, 2]))), None)]
    fn test_value_ast_as_lit(#[case] ast: GraphIrValueAst, #[case] expected: Option<i64>) {
        Python::attach(|py| {
            assert_eq!(ValueAst::from_rust(py, &ast).unwrap().as_lit(py), expected);
        });
    }
}
