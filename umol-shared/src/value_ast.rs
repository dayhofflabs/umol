//! Value AST: integer literals, sets, and arithmetic/boolean expressions.

use std::collections::HashMap;

use thiserror::Error;

/// Variable bindings used by [`Expr::evaluate`] and [`Expr::evaluate_bool`].
pub type Bindings = HashMap<String, i64>;

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum ValueAst {
    #[default]
    Undetermined,
    LitSet(Vec<i64>),
    Lit(i64),
    Expr(Expr),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Expr {
    Lit(i64),
    Var(String),
    Neg(Box<Expr>),
    BinOp(Box<Expr>, ArithOp, Box<Expr>),
    Mem(Box<Expr>, Vec<i64>),
    Rel(Box<Expr>, RelOp, Box<Expr>),
    Not(Box<Expr>),
    And(Vec<Expr>),
    Or(Vec<Expr>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ArithOp {
    Add,
    Sub,
    Mul,
    Div,
    Rem,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RelOp {
    Le,
    Ge,
    Eq,
    Lt,
    Gt,
}

#[derive(Clone, Debug, PartialEq, Error)]
pub enum EvaluationError {
    #[error("Unbound variable: {0}")]
    UnboundVariable(String),
    #[error("Division by zero")]
    DivisionByZero,
    #[error("Type mismatch")]
    TypeMismatch,
}

impl ValueAst {
    pub fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    /// Match a concrete integer value against this pattern
    pub fn matches(&self, value: i64) -> bool {
        self.capture(value).is_some()
    }

    /// Match a concrete integer value against this pattern, returning variable bindings
    ///
    /// Variables in the pattern are bound to `value`. For boolean expressions the
    /// predicate is evaluated with those bindings; for arithmetic expressions the
    /// result is compared to `value`
    pub fn capture(&self, value: i64) -> Option<Bindings> {
        match self {
            ValueAst::Undetermined => Some(Bindings::new()),
            ValueAst::Lit(n) => {
                if *n == value {
                    Some(Bindings::new())
                } else {
                    None
                }
            }
            ValueAst::LitSet(s) => {
                if s.contains(&value) {
                    Some(Bindings::new())
                } else {
                    None
                }
            }
            ValueAst::Expr(e) => {
                let mut bindings = Bindings::new();
                collect_bindings(e, value, &mut bindings);
                if e.is_arithmetic() {
                    match e.evaluate(&bindings) {
                        Ok(v) if v == value => Some(bindings),
                        _ => None,
                    }
                } else {
                    match e.evaluate_bool(&bindings) {
                        Ok(true) => Some(bindings),
                        _ => None,
                    }
                }
            }
        }
    }
}

/// Recursively bind every variable in `expr` to `value`
fn collect_bindings(expr: &Expr, value: i64, bindings: &mut Bindings) {
    match expr {
        Expr::Var(name) => {
            bindings.insert(name.clone(), value);
        }
        Expr::Neg(e) => collect_bindings(e, value, bindings),
        Expr::BinOp(l, _, r) => {
            collect_bindings(l, value, bindings);
            collect_bindings(r, value, bindings);
        }
        Expr::Mem(e, _) => collect_bindings(e, value, bindings),
        Expr::Rel(l, _, r) => {
            collect_bindings(l, value, bindings);
            collect_bindings(r, value, bindings);
        }
        Expr::Not(e) => collect_bindings(e, value, bindings),
        Expr::And(exprs) | Expr::Or(exprs) => {
            for e in exprs {
                collect_bindings(e, value, bindings);
            }
        }
        Expr::Lit(_) => {}
    }
}

impl Expr {
    pub fn is_arithmetic(&self) -> bool {
        matches!(
            self,
            Expr::Lit(..) | Expr::Var(..) | Expr::Neg(..) | Expr::BinOp(..)
        )
    }

    /// Evaluate an arithmetic expression to an `i64`
    ///
    /// Returns [`EvaluationError::TypeMismatch`] if called on a boolean-domain
    /// expression (`Rel`, `Mem`, `Not`, `And`, `Or`)
    pub fn evaluate(&self, vars: &Bindings) -> Result<i64, EvaluationError> {
        match self {
            Expr::Lit(n) => Ok(*n),
            Expr::Var(name) => vars
                .get(name)
                .copied()
                .ok_or_else(|| EvaluationError::UnboundVariable(name.clone())),
            Expr::Neg(e) => Ok(-e.evaluate(vars)?),
            Expr::BinOp(l, op, r) => {
                let l = l.evaluate(vars)?;
                let r = r.evaluate(vars)?;
                match op {
                    ArithOp::Add => Ok(l + r),
                    ArithOp::Sub => Ok(l - r),
                    ArithOp::Mul => Ok(l * r),
                    ArithOp::Div => {
                        if r == 0 {
                            Err(EvaluationError::DivisionByZero)
                        } else {
                            Ok(l / r)
                        }
                    }
                    ArithOp::Rem => {
                        if r == 0 {
                            Err(EvaluationError::DivisionByZero)
                        } else {
                            Ok(l % r)
                        }
                    }
                }
            }
            Expr::Rel(..) | Expr::Mem(..) | Expr::Not(..) | Expr::And(..) | Expr::Or(..) => {
                Err(EvaluationError::TypeMismatch)
            }
        }
    }

    /// Evaluate a boolean expression to a `bool`
    ///
    /// Returns [`EvaluationError::TypeMismatch`] if called on an arithmetic-domain
    /// expression (`Lit`, `Var`, `Neg`, `BinOp`)
    pub fn evaluate_bool(&self, vars: &Bindings) -> Result<bool, EvaluationError> {
        match self {
            Expr::Rel(l, op, r) => {
                let l = l.evaluate(vars)?;
                let r = r.evaluate(vars)?;
                Ok(match op {
                    RelOp::Le => l <= r,
                    RelOp::Ge => l >= r,
                    RelOp::Eq => l == r,
                    RelOp::Lt => l < r,
                    RelOp::Gt => l > r,
                })
            }
            Expr::Mem(e, set) => Ok(set.contains(&e.evaluate(vars)?)),
            Expr::Not(e) => Ok(!e.evaluate_bool(vars)?),
            Expr::And(exprs) => {
                for e in exprs {
                    if !e.evaluate_bool(vars)? {
                        return Ok(false);
                    }
                }
                Ok(true)
            }
            Expr::Or(exprs) => {
                for e in exprs {
                    if e.evaluate_bool(vars)? {
                        return Ok(true);
                    }
                }
                Ok(false)
            }
            Expr::Lit(..) | Expr::Var(..) | Expr::Neg(..) | Expr::BinOp(..) => {
                Err(EvaluationError::TypeMismatch)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(Expr::Lit(5), Bindings::new(), 5)]
    #[case::var_bound(Expr::Var("x".to_string()), Bindings::from([("x".to_string(), 3)]), 3)]
    #[case::neg(Expr::Neg(Box::new(Expr::Lit(3))), Bindings::new(), -3)]
    #[case::add(Expr::BinOp(Box::new(Expr::Lit(2)), ArithOp::Add, Box::new(Expr::Lit(3))), Bindings::new(), 5)]
    #[case::sub(Expr::BinOp(Box::new(Expr::Lit(5)), ArithOp::Sub, Box::new(Expr::Lit(3))), Bindings::new(), 2)]
    #[case::mul(Expr::BinOp(Box::new(Expr::Lit(3)), ArithOp::Mul, Box::new(Expr::Lit(4))), Bindings::new(), 12)]
    #[case::div(Expr::BinOp(Box::new(Expr::Lit(10)), ArithOp::Div, Box::new(Expr::Lit(3))), Bindings::new(), 3)]
    #[case::rem(Expr::BinOp(Box::new(Expr::Lit(10)), ArithOp::Rem, Box::new(Expr::Lit(3))), Bindings::new(), 1)]
    fn test_evaluate(#[case] expr: Expr, #[case] vars: Bindings, #[case] expected: i64) {
        let result = expr.evaluate(&vars);
        assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", expr, result.clone().unwrap_err());
        assert_eq!(result.unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::var_unbound(Expr::Var("x".to_string()), Bindings::new(), EvaluationError::UnboundVariable("x".to_string()))]
    #[case::div_zero(Expr::BinOp(Box::new(Expr::Lit(10)), ArithOp::Div, Box::new(Expr::Lit(0))), Bindings::new(), EvaluationError::DivisionByZero)]
    #[case::rem_zero(Expr::BinOp(Box::new(Expr::Lit(10)), ArithOp::Rem, Box::new(Expr::Lit(0))), Bindings::new(), EvaluationError::DivisionByZero)]
    #[case::type_mismatch(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(1))), Bindings::new(), EvaluationError::TypeMismatch)]
    fn test_evaluate_invalid(#[case] expr: Expr, #[case] vars: Bindings, #[case] expected: EvaluationError) {
        let result = expr.evaluate(&vars);
        assert!(result.is_err(), "{:?} should have failed, error: {:?}", expr, result.clone().unwrap_err());
        assert_eq!(result.unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::rel_eq_true(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(1))), Bindings::new(), true)]
    #[case::rel_eq_false(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(2))), Bindings::new(), false)]
    #[case::rel_lt_true(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Lt, Box::new(Expr::Lit(2))), Bindings::new(), true)]
    #[case::rel_lt_false(Expr::Rel(Box::new(Expr::Lit(2)), RelOp::Lt, Box::new(Expr::Lit(1))), Bindings::new(), false)]
    #[case::rel_le_true(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Le, Box::new(Expr::Lit(1))), Bindings::new(), true)]
    #[case::rel_le_false(Expr::Rel(Box::new(Expr::Lit(2)), RelOp::Le, Box::new(Expr::Lit(1))), Bindings::new(), false)]
    #[case::rel_gt_true(Expr::Rel(Box::new(Expr::Lit(2)), RelOp::Gt, Box::new(Expr::Lit(1))), Bindings::new(), true)]
    #[case::rel_gt_false(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Gt, Box::new(Expr::Lit(2))), Bindings::new(), false)]
    #[case::rel_ge_true(Expr::Rel(Box::new(Expr::Lit(2)), RelOp::Ge, Box::new(Expr::Lit(2))), Bindings::new(), true)]
    #[case::rel_ge_false(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Ge, Box::new(Expr::Lit(2))), Bindings::new(), false)]
    #[case::mem_true(Expr::Mem(Box::new(Expr::Lit(2)), vec![1, 2, 3]), Bindings::new(), true)]
    #[case::mem_false(Expr::Mem(Box::new(Expr::Lit(4)), vec![1, 2, 3]), Bindings::new(), false)]
    #[case::not_true(Expr::Not(Box::new(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(2))))), Bindings::new(), true)]
    #[case::not_false(Expr::Not(Box::new(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(1))))), Bindings::new(), false)]
    #[case::and_true(Expr::And(vec![Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Lt, Box::new(Expr::Lit(2))), Expr::Rel(Box::new(Expr::Lit(3)), RelOp::Gt, Box::new(Expr::Lit(2)))]), Bindings::new(), true)]
    #[case::and_false(Expr::And(vec![Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Lt, Box::new(Expr::Lit(2))), Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Gt, Box::new(Expr::Lit(2)))]), Bindings::new(), false)]
    #[case::or_true(Expr::Or(vec![Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(2))), Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Lt, Box::new(Expr::Lit(2)))]), Bindings::new(), true)]
    #[case::or_false(Expr::Or(vec![Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(2))), Expr::Rel(Box::new(Expr::Lit(3)), RelOp::Lt, Box::new(Expr::Lit(2)))]), Bindings::new(), false)]
    #[case::var_in_rel(Expr::Rel(Box::new(Expr::Var("x".to_string())), RelOp::Gt, Box::new(Expr::Lit(0))), Bindings::from([("x".to_string(), 5)]), true)]
    fn test_evaluate_bool(#[case] expr: Expr, #[case] vars: Bindings, #[case] expected: bool) {
        let result = expr.evaluate_bool(&vars);
        assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", expr, result.clone().unwrap_err());
        assert_eq!(result.unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unbound_in_rel(Expr::Rel(Box::new(Expr::Var("x".to_string())), RelOp::Gt, Box::new(Expr::Lit(0))), Bindings::new(), EvaluationError::UnboundVariable("x".to_string()))]
    #[case::type_mismatch(Expr::Lit(1), Bindings::new(), EvaluationError::TypeMismatch)]
    fn test_evaluate_bool_invalid(#[case] expr: Expr, #[case] vars: Bindings, #[case] expected: EvaluationError) {
        let result = expr.evaluate_bool(&vars);
        assert!(result.is_err(), "{:?} should have failed, error: {:?}", expr, result.clone().unwrap_err());
        assert_eq!(result.unwrap_err(), expected);
    }

    #[rstest]
    #[case::lit(ValueAst::Lit(3), true)]
    #[case::wildcard(ValueAst::Undetermined, false)]
    #[case::lit_set(ValueAst::LitSet(vec![1, 2]), false)]
    #[case::expr(ValueAst::Expr(Expr::Var("x".to_string())), false)]
    fn test_value_ast_is_ground(#[case] ast: ValueAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(ValueAst::Undetermined, 3, true)]
    #[case::lit_match(ValueAst::Lit(3), 3,  true)]
    #[case::lit_set_match(ValueAst::LitSet(vec![1, 2, 3]), 2, true)]
    #[case::expr_var(ValueAst::Expr(Expr::Var("h".to_string())), 5, true)]
    #[case::expr_lit_match(ValueAst::Expr(Expr::Lit(3)), 3, true)]
    #[case::expr_rel_match(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), 3, true)]
    #[case::expr_mem_match(ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![0, 1])), 1, true)]
    #[case::lit_no_match(ValueAst::Lit(3), 4, false)]
    #[case::expr_lit_no_match(ValueAst::Expr(Expr::Lit(3)), 4, false)]
    #[case::expr_rel_no_match(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), 0, false)]
    fn test_matches(#[case] pattern: ValueAst, #[case] value: i64, #[case] expected: bool) {
        assert_eq!(pattern.matches(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(ValueAst::Undetermined, 3, Bindings::new())]
    #[case::lit_match(ValueAst::Lit(3), 3, Bindings::new())]
    #[case::lit_set_match(ValueAst::LitSet(vec![1, 2, 3]), 2, Bindings::new())]
    #[case::expr_var(ValueAst::Expr(Expr::Var("h".to_string())), 5, Bindings::from([("h".to_string(), 5)]))]
    #[case::expr_lit_match(ValueAst::Expr(Expr::Lit(3)), 3, Bindings::new())]
    #[case::expr_rel_match(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), 3, Bindings::from([("h".to_string(), 3)]))]
    #[case::expr_mem_match(ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![0, 1])), 1, Bindings::from([("h".to_string(), 1)]))]
    fn test_capture(#[case] pattern: ValueAst, #[case] value: i64, #[case] expected: Bindings) {
        assert_eq!(pattern.capture(value), Some(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_no_match(ValueAst::Lit(3), 4)]
    #[case::expr_lit_no_match(ValueAst::Expr(Expr::Lit(3)), 4)]
    #[case::expr_rel_no_match(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), 0)]
    fn test_capture_no_match(#[case] pattern: ValueAst, #[case] value: i64) {
        assert_eq!(pattern.capture(value), None);
    }
}
