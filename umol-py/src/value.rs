//! `ArithExpr` — arithmetic expression over integers, matching `umol_graph_ir::ir::ArithExpr`
//! as a native PyO3 complex enum. IR recursion (`Box<Self>` / `Vec<Self>`) becomes
//! `Py<Self>` / `Vec<Py<Self>>`; per-variant construction and `match` work natively
//! on the Python side.
// Blanket-allow the `absolute_paths` false positives from pyo3's `hash` derive
// (hygienic `::std::…` paths). Hand-written code here imports at top.
#![allow(clippy::absolute_paths)]

use std::collections::BTreeSet;

use pyo3::prelude::*;
use umol_graph_ir::ir::{
    ArithExpr as GraphIrArithExpr, AsLit, MemOp as GraphIrMemOp, NumForm as GraphIrNumForm,
    PredExpr as GraphIrPredExpr, RelOp as GraphIrRelOp,
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

/// Arithmetic expression over integers.
#[pyclass]
pub enum ArithExpr {
    Lit(i64),
    Var(String),
    Neg(Py<ArithExpr>),
    Sum(Vec<Py<ArithExpr>>),
    Product(Vec<Py<ArithExpr>>),
    Div(Py<ArithExpr>, Py<ArithExpr>),
    Rem(Py<ArithExpr>, Py<ArithExpr>),
}

#[pymethods]
impl ArithExpr {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            ArithExpr::Lit(_) => ("Lit", 1),
            ArithExpr::Var(_) => ("Var", 1),
            ArithExpr::Neg(_) => ("Neg", 1),
            ArithExpr::Sum(_) => ("Sum", 1),
            ArithExpr::Product(_) => ("Product", 1),
            ArithExpr::Div(_, _) => ("Div", 2),
            ArithExpr::Rem(_, _) => ("Rem", 2),
        };
        variant_repr(slf.bind(py).as_any(), "ArithExpr", variant, arity)
    }
}

impl ArithExpr {
    /// Build the Python value from the IR expression (one Python object per node).
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrArithExpr) -> PyResult<ArithExpr> {
        Ok(match ast {
            GraphIrArithExpr::Lit(n) => ArithExpr::Lit(*n),
            GraphIrArithExpr::Var(name) => ArithExpr::Var(name.clone()),
            GraphIrArithExpr::Neg(t) => {
                ArithExpr::Neg(into_py_variant(py, Self::from_rust(py, t)?)?)
            }
            GraphIrArithExpr::Sum(terms) => ArithExpr::Sum(
                terms
                    .iter()
                    .map(|t| into_py_variant(py, Self::from_rust(py, t)?))
                    .collect::<PyResult<_>>()?,
            ),
            GraphIrArithExpr::Product(terms) => ArithExpr::Product(
                terms
                    .iter()
                    .map(|t| into_py_variant(py, Self::from_rust(py, t)?))
                    .collect::<PyResult<_>>()?,
            ),
            GraphIrArithExpr::Div(a, b) => ArithExpr::Div(
                into_py_variant(py, Self::from_rust(py, a)?)?,
                into_py_variant(py, Self::from_rust(py, b)?)?,
            ),
            GraphIrArithExpr::Rem(a, b) => ArithExpr::Rem(
                into_py_variant(py, Self::from_rust(py, a)?)?,
                into_py_variant(py, Self::from_rust(py, b)?)?,
            ),
        })
    }

    /// Convert the Python value back to the IR expression.
    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrArithExpr {
        match self {
            ArithExpr::Lit(n) => GraphIrArithExpr::Lit(*n),
            ArithExpr::Var(name) => GraphIrArithExpr::Var(name.clone()),
            ArithExpr::Neg(t) => GraphIrArithExpr::Neg(Box::new(t.bind(py).borrow().to_rust(py))),
            ArithExpr::Sum(terms) => GraphIrArithExpr::Sum(
                terms
                    .iter()
                    .map(|t| t.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
            ArithExpr::Product(terms) => GraphIrArithExpr::Product(
                terms
                    .iter()
                    .map(|t| t.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
            ArithExpr::Div(a, b) => GraphIrArithExpr::Div(
                Box::new(a.bind(py).borrow().to_rust(py)),
                Box::new(b.bind(py).borrow().to_rust(py)),
            ),
            ArithExpr::Rem(a, b) => GraphIrArithExpr::Rem(
                Box::new(a.bind(py).borrow().to_rust(py)),
                Box::new(b.bind(py).borrow().to_rust(py)),
            ),
        }
    }
}

impl_py_canonicalize!(
    ArithExpr,
    GraphIrArithExpr,
    |value: &ArithExpr, py: Python<'_>| -> PyResult<GraphIrArithExpr> { Ok(value.to_rust(py)) },
    |py: Python<'_>, value: GraphIrArithExpr| -> PyResult<ArithExpr> {
        ArithExpr::from_rust(py, &value)
    }
);

/// Predicate expression over arithmetic expressions.
#[pyclass]
pub enum PredExpr {
    Rel(Py<ArithExpr>, RelOp, Py<ArithExpr>),
    Mem(Py<ArithExpr>, MemOp, BTreeSet<i64>),
    Not(Py<PredExpr>),
    And(Vec<Py<PredExpr>>),
    Or(Vec<Py<PredExpr>>),
}

#[pymethods]
impl PredExpr {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            PredExpr::Rel(_, _, _) => ("Rel", 3),
            PredExpr::Mem(_, _, _) => ("Mem", 3),
            PredExpr::Not(_) => ("Not", 1),
            PredExpr::And(_) => ("And", 1),
            PredExpr::Or(_) => ("Or", 1),
        };
        variant_repr(slf.bind(py).as_any(), "PredExpr", variant, arity)
    }
}

impl PredExpr {
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrPredExpr) -> PyResult<PredExpr> {
        Ok(match ast {
            GraphIrPredExpr::Rel(a, op, b) => PredExpr::Rel(
                into_py_variant(py, ArithExpr::from_rust(py, a)?)?,
                RelOp::from_rust(*op),
                into_py_variant(py, ArithExpr::from_rust(py, b)?)?,
            ),
            GraphIrPredExpr::Mem(t, op, members) => PredExpr::Mem(
                into_py_variant(py, ArithExpr::from_rust(py, t)?)?,
                MemOp::from_rust(*op),
                members.clone(),
            ),
            GraphIrPredExpr::Not(p) => PredExpr::Not(into_py_variant(py, Self::from_rust(py, p)?)?),
            GraphIrPredExpr::And(ps) => PredExpr::And(
                ps.iter()
                    .map(|p| into_py_variant(py, Self::from_rust(py, p)?))
                    .collect::<PyResult<_>>()?,
            ),
            GraphIrPredExpr::Or(ps) => PredExpr::Or(
                ps.iter()
                    .map(|p| into_py_variant(py, Self::from_rust(py, p)?))
                    .collect::<PyResult<_>>()?,
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrPredExpr {
        match self {
            PredExpr::Rel(a, op, b) => GraphIrPredExpr::Rel(
                a.bind(py).borrow().to_rust(py),
                op.to_rust(),
                b.bind(py).borrow().to_rust(py),
            ),
            PredExpr::Mem(t, op, members) => GraphIrPredExpr::Mem(
                t.bind(py).borrow().to_rust(py),
                op.to_rust(),
                members.clone(),
            ),
            PredExpr::Not(p) => GraphIrPredExpr::Not(Box::new(p.bind(py).borrow().to_rust(py))),
            PredExpr::And(ps) => {
                GraphIrPredExpr::And(ps.iter().map(|p| p.bind(py).borrow().to_rust(py)).collect())
            }
            PredExpr::Or(ps) => {
                GraphIrPredExpr::Or(ps.iter().map(|p| p.bind(py).borrow().to_rust(py)).collect())
            }
        }
    }
}

/// Integer-valued atom/bond field: the undetermined wildcard, a literal, a literal
/// set, a range, an arithmetic expression, or a predicate expression.
#[pyclass]
pub enum NumForm {
    Undetermined(),
    Lit(i64),
    LitSet(BTreeSet<i64>),
    RangeFrom(i64),
    RangeTo(i64),
    ArithExpr(Py<ArithExpr>),
    PredExpr(Py<PredExpr>),
}

#[pymethods]
impl NumForm {
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
            NumForm::Undetermined() => ("Undetermined", 0),
            NumForm::Lit(_) => ("Lit", 1),
            NumForm::LitSet(_) => ("LitSet", 1),
            NumForm::RangeFrom(_) => ("RangeFrom", 1),
            NumForm::RangeTo(_) => ("RangeTo", 1),
            NumForm::ArithExpr(_) => ("ArithExpr", 1),
            NumForm::PredExpr(_) => ("PredExpr", 1),
        };
        variant_repr(slf.bind(py).as_any(), "NumForm", variant, arity)
    }
}

impl NumForm {
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrNumForm) -> PyResult<NumForm> {
        Ok(match ast {
            GraphIrNumForm::Undetermined => NumForm::Undetermined(),
            GraphIrNumForm::Lit(n) => NumForm::Lit(*n),
            GraphIrNumForm::LitSet(members) => NumForm::LitSet((**members).clone()),
            GraphIrNumForm::RangeFrom(n) => NumForm::RangeFrom(*n),
            GraphIrNumForm::RangeTo(n) => NumForm::RangeTo(*n),
            GraphIrNumForm::ArithExpr(t) => {
                NumForm::ArithExpr(into_py_variant(py, ArithExpr::from_rust(py, t)?)?)
            }
            GraphIrNumForm::PredExpr(p) => {
                NumForm::PredExpr(into_py_variant(py, PredExpr::from_rust(py, p)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrNumForm {
        match self {
            NumForm::Undetermined() => GraphIrNumForm::Undetermined,
            NumForm::Lit(n) => GraphIrNumForm::Lit(*n),
            NumForm::LitSet(members) => GraphIrNumForm::LitSet(Box::new(members.clone())),
            NumForm::RangeFrom(n) => GraphIrNumForm::RangeFrom(*n),
            NumForm::RangeTo(n) => GraphIrNumForm::RangeTo(*n),
            NumForm::ArithExpr(t) => {
                GraphIrNumForm::ArithExpr(Box::new(t.bind(py).borrow().to_rust(py)))
            }
            NumForm::PredExpr(p) => {
                GraphIrNumForm::PredExpr(Box::new(p.bind(py).borrow().to_rust(py)))
            }
        }
    }
}

impl_py_lattice!(
    NumForm,
    GraphIrNumForm,
    |value: &NumForm, py: Python<'_>| -> PyResult<GraphIrNumForm> { Ok(value.to_rust(py)) },
    |py: Python<'_>, value: GraphIrNumForm| -> PyResult<NumForm> { NumForm::from_rust(py, &value) }
);

/// A Python `NumForm` or `int` (→ `NumForm::Lit`), matching `impl Into<NumForm>`
/// on the Rust builders. The `*Like` convention for binding coercion inputs (`*Input`
/// is the DSL side); shared by the atom fields, unpaired-electron components, and
/// ring-membership count.
#[derive(FromPyObject)]
pub enum NumLike {
    Ast(Py<NumForm>),
    Lit(i64),
}

impl NumLike {
    /// Coerce to the numeric form used by `impl Into<NumForm>` Rust builders.
    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrNumForm {
        match self {
            NumLike::Ast(value) => value.bind(py).borrow().to_rust(py),
            NumLike::Lit(number) => GraphIrNumForm::Lit(*number),
        }
    }

    /// Coerce to a `Py<NumForm>` (for value structs that store the value field).
    pub(crate) fn to_py(&self, py: Python<'_>) -> PyResult<Py<NumForm>> {
        match self {
            NumLike::Ast(value) => Ok(value.clone_ref(py)),
            NumLike::Lit(number) => into_py_variant(py, NumForm::Lit(*number)),
        }
    }
}

/// `IntoPyObject` for `&NumLike` so it can be a complex-enum field: constructors
/// (`AromaticValenceAst.Aromatic(1)`) coerce `int | NumForm` in, and the field
/// reads back as a `NumForm`.
impl<'py> IntoPyObject<'py> for &NumLike {
    type Target = NumForm;
    type Output = Bound<'py, NumForm>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Bound<'py, NumForm>> {
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
    #[case(GraphIrArithExpr::Lit(5))]
    #[case(GraphIrArithExpr::Var("x".to_string()))]
    #[case(GraphIrArithExpr::Neg(Box::new(GraphIrArithExpr::Lit(3))))]
    #[case(GraphIrArithExpr::Sum(vec![GraphIrArithExpr::Lit(1), GraphIrArithExpr::Lit(2)]))]
    #[case(GraphIrArithExpr::Div(
        Box::new(GraphIrArithExpr::Lit(6)),
        Box::new(GraphIrArithExpr::Lit(2))
    ))]
    fn test_arith_expr_roundtrip(#[case] ast: GraphIrArithExpr) {
        Python::attach(|py| {
            let value = ArithExpr::from_rust(py, &ast).unwrap();
            assert_eq!(value.to_rust(py), ast);
        });
    }

    #[rstest]
    #[case(GraphIrPredExpr::Rel(GraphIrArithExpr::Var("h".into()), GraphIrRelOp::Le, GraphIrArithExpr::Lit(3)))]
    #[case(GraphIrPredExpr::Mem(GraphIrArithExpr::Lit(0), GraphIrMemOp::In, BTreeSet::from([1, 2, 3])))]
    #[case(GraphIrPredExpr::Not(Box::new(GraphIrPredExpr::Rel(
        GraphIrArithExpr::Lit(1),
        GraphIrRelOp::Eq,
        GraphIrArithExpr::Lit(1),
    ))))]
    #[case(GraphIrPredExpr::And(vec![
        GraphIrPredExpr::Rel(GraphIrArithExpr::Lit(1), GraphIrRelOp::Lt, GraphIrArithExpr::Lit(2)),
        GraphIrPredExpr::Rel(GraphIrArithExpr::Lit(3), GraphIrRelOp::Gt, GraphIrArithExpr::Lit(2)),
    ]))]
    fn test_pred_expr_roundtrip(#[case] ast: GraphIrPredExpr) {
        Python::attach(|py| {
            let value = PredExpr::from_rust(py, &ast).unwrap();
            assert_eq!(value.to_rust(py), ast);
        });
    }

    #[rstest]
    #[case(GraphIrNumForm::Undetermined)]
    #[case(GraphIrNumForm::Lit(7))]
    #[case(GraphIrNumForm::LitSet(Box::new(BTreeSet::from([1, 2, 3]))))]
    #[case(GraphIrNumForm::RangeFrom(1))]
    #[case(GraphIrNumForm::RangeTo(9))]
    #[case(GraphIrNumForm::ArithExpr(Box::new(GraphIrArithExpr::Var("x".into()))))]
    #[case(GraphIrNumForm::PredExpr(Box::new(GraphIrPredExpr::Rel(
        GraphIrArithExpr::Var("h".into()),
        GraphIrRelOp::Le,
        GraphIrArithExpr::Lit(3),
    ))))]
    fn test_value_ast_roundtrip(#[case] ast: GraphIrNumForm) {
        Python::attach(|py| {
            let value = NumForm::from_rust(py, &ast).unwrap();
            assert_eq!(value.to_rust(py), ast);
        });
    }

    #[rstest]
    #[case(GraphIrNumForm::Lit(4), Some(4))]
    #[case(GraphIrNumForm::Lit(-1), Some(-1))]
    #[case(GraphIrNumForm::Undetermined, None)]
    #[case(GraphIrNumForm::RangeFrom(1), None)]
    #[case(GraphIrNumForm::LitSet(Box::new(BTreeSet::from([1, 2]))), None)]
    fn test_value_ast_as_lit(#[case] ast: GraphIrNumForm, #[case] expected: Option<i64>) {
        Python::attach(|py| {
            assert_eq!(NumForm::from_rust(py, &ast).unwrap().as_lit(py), expected);
        });
    }
}
