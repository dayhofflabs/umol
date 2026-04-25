//! Value AST.

use std::collections::HashMap;

use super::error::EvaluationError;

/// Variable bindings used by [`Expr::evaluate`] and [`Expr::evaluate_bool`].
pub type Bindings = HashMap<String, i64>;

/// Integer-valued atom/bond field: undetermined (pattern wildcard), a
/// literal, a finite literal set, or an arithmetic/boolean expression
/// pattern. Used for charge, hydrogen count, isotope mass, valence, etc.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum ValueAst {
    #[default]
    Undetermined,
    LitSet(Vec<i64>),
    Lit(i64),
    Expr(Expr),
}

impl ValueAst {
    pub fn new(value: i64) -> Self {
        Self::Lit(value)
    }

    /// The pattern denotes a single concrete integer. Semantic, not
    /// syntactic: `Expr` that folds to a constant is ground, and a
    /// `LitSet` of a single value (regardless of duplicates) is ground.
    ///
    /// Fast path — `Lit` and `Undetermined` dispatch with two tag compares
    /// so the common case (a fully-lowered ground molecule) doesn't pay
    /// for the `LitSet`/`Expr` logic.
    #[inline]
    pub fn is_ground(&self) -> bool {
        match self {
            Self::Lit(_) => true,
            Self::Undetermined => false,
            _ => self.is_ground_slow(),
        }
    }

    #[inline(never)]
    #[cold]
    fn is_ground_slow(&self) -> bool {
        match self {
            Self::LitSet(s) => litset_is_ground(s),
            Self::Expr(e) => e.is_ground(),
            Self::Lit(_) | Self::Undetermined => unreachable!(),
        }
    }

    pub fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    /// Pattern matches target iff every integer target admits is also
    /// admitted by pattern (superset semantics). `Expr` targets cannot
    /// be certified generically and are rejected here; pattern `Expr`
    /// is evaluated pointwise on ground / set targets.
    pub fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (_, Self::Expr(_)) => false,
            (pattern, Self::Lit(n)) => pattern.matches_value(*n),
            (pattern, Self::LitSet(ns)) => ns.iter().all(|n| pattern.matches_value(*n)),
        }
    }

    /// Match a concrete integer value against this pattern.
    pub fn matches_value(&self, value: i64) -> bool {
        self.capture(value).is_some()
    }

    /// Match a concrete integer value against this pattern, returning variable bindings
    ///
    /// Variables in the pattern are bound to `value`. For boolean expressions the
    /// predicate is evaluated with those bindings; for arithmetic expressions the
    /// result is compared to `value`
    /// Reduce to canonical form by lifting trivial `Expr` wrappers and
    /// recursively simplifying any inner expression. Specifically:
    ///
    /// - `Expr(Expr::Lit(n))` → `Lit(n)`
    /// - `Expr(Expr::Neg(Expr::Lit(n)))` → `Lit(-n)` when `-n` does not
    ///   overflow `i64` (otherwise the wrapped form is preserved)
    /// - `Expr(e)` for any other shape → `Expr(e.simplify())`
    /// - `Lit` / `LitSet` / `Undetermined` → unchanged
    pub fn simplify(self) -> Self {
        match self {
            Self::Expr(e) => match e.simplify() {
                Expr::Lit(n) => Self::Lit(n),
                Expr::Neg(inner) => match *inner {
                    Expr::Lit(n) => match n.checked_neg() {
                        Some(neg) => Self::Lit(neg),
                        None => Self::Expr(Expr::Neg(Box::new(Expr::Lit(n)))),
                    },
                    other => Self::Expr(Expr::Neg(Box::new(other))),
                },
                other => Self::Expr(other),
            },
            other => other,
        }
    }

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

/// Arithmetic / boolean expression tree over `ValueAst`. Captures
/// atom/bond field constraints that can't be expressed as a literal or
/// literal set, including bound variables (`Var`), membership tests
/// (`Mem`), relational comparisons (`Rel`), and boolean combinators.
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

impl Expr {
    /// Reduce to canonical form. Recursively simplifies children, then
    /// applies one round of structural folding:
    ///
    /// - `Lit(n)` for `n < 0` → `Neg(Lit(-n))`. Inside an `Expr` context
    ///   the parser always produces `Neg(Lit(n))` for `-n`; only the
    ///   ValueAst-level `Lit` slot reads signed integers directly.
    ///   `Lit(i64::MIN)` is left as-is (cannot be negated without overflow).
    /// - `Neg(Neg(e))` → `e`
    /// - `Or(... Or(inner) ...)` flattens the inner `Or` one level
    ///   (recursively, so the result has no `Or` direct child)
    /// - `And(... And(inner) ...)` flattens identically
    /// - `And([single])` / `Or([single])` → `single`. The renderer emits a
    ///   single-child And/Or as just the child, and the parser reads it
    ///   back without the wrapper.
    ///
    /// Idempotent. Mirrors the parser's normalization.
    pub fn simplify(self) -> Self {
        match self {
            Expr::Lit(n) if n < 0 => match n.checked_neg() {
                Some(pos) => Expr::Neg(Box::new(Expr::Lit(pos))),
                None => Expr::Lit(n),
            },
            Expr::Lit(_) | Expr::Var(_) => self,
            Expr::Neg(inner) => match inner.simplify() {
                Expr::Neg(grand) => *grand,
                other => Expr::Neg(Box::new(other)),
            },
            Expr::BinOp(l, op, r) => {
                Expr::BinOp(Box::new(l.simplify()), op, Box::new(r.simplify()))
            }
            Expr::Mem(e, set) => Expr::Mem(Box::new(e.simplify()), set),
            Expr::Rel(l, op, r) => Expr::Rel(Box::new(l.simplify()), op, Box::new(r.simplify())),
            Expr::Not(e) => Expr::Not(Box::new(e.simplify())),
            Expr::And(exprs) => {
                let mut flat = flatten_simplified(exprs, |e| matches!(e, Expr::And(_)));
                if flat.len() == 1 {
                    flat.pop().unwrap()
                } else {
                    Expr::And(flat)
                }
            }
            Expr::Or(exprs) => {
                let mut flat = flatten_simplified(exprs, |e| matches!(e, Expr::Or(_)));
                if flat.len() == 1 {
                    flat.pop().unwrap()
                } else {
                    Expr::Or(flat)
                }
            }
        }
    }

    pub fn is_arithmetic(&self) -> bool {
        matches!(
            self,
            Expr::Lit(..) | Expr::Var(..) | Expr::Neg(..) | Expr::BinOp(..)
        )
    }

    /// The expression denotes a single concrete integer: it is arithmetic
    /// (not a boolean-domain predicate) and evaluates without free variables
    /// or error (including i64 overflow via the checked evaluator). A
    /// `ValueAst::Expr` containing a ground expression is itself ground.
    pub fn is_ground(&self) -> bool {
        self.is_arithmetic() && self.evaluate_checked(&Bindings::new()).is_some()
    }

    /// Overflow-safe arithmetic evaluation. Returns `None` for free
    /// variables, division/remainder by zero, type mismatch (boolean-domain
    /// Expr), and any `i64` overflow in `Neg`/`BinOp`. Intended as the
    /// foundation of [`Expr::is_ground`]; for error-reporting callers use
    /// [`Expr::evaluate`].
    pub fn evaluate_checked(&self, vars: &Bindings) -> Option<i64> {
        match self {
            Expr::Lit(n) => Some(*n),
            Expr::Var(name) => vars.get(name).copied(),
            Expr::Neg(e) => e.evaluate_checked(vars)?.checked_neg(),
            Expr::BinOp(l, op, r) => {
                let l = l.evaluate_checked(vars)?;
                let r = r.evaluate_checked(vars)?;
                match op {
                    ArithOp::Add => l.checked_add(r),
                    ArithOp::Sub => l.checked_sub(r),
                    ArithOp::Mul => l.checked_mul(r),
                    ArithOp::Div => l.checked_div(r),
                    ArithOp::Rem => l.checked_rem(r),
                }
            }
            Expr::Rel(..) | Expr::Mem(..) | Expr::Not(..) | Expr::And(..) | Expr::Or(..) => None,
        }
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

/// Simplify each child, then flatten any whose simplified form satisfies
/// `is_same_op` (i.e. an `Or` child of an `Or` parent, or `And` of `And`).
/// One level of flattening is enough: each child has already been simplified,
/// so its own children are themselves not same-op containers.
fn flatten_simplified(exprs: Vec<Expr>, is_same_op: impl Fn(&Expr) -> bool) -> Vec<Expr> {
    let mut out = Vec::with_capacity(exprs.len());
    for child in exprs {
        let simplified = child.simplify();
        if is_same_op(&simplified) {
            match simplified {
                Expr::And(inner) | Expr::Or(inner) => out.extend(inner),
                _ => unreachable!("is_same_op rejects non-And/Or"),
            }
        } else {
            out.push(simplified);
        }
    }
    out
}

/// A `LitSet` is ground iff non-empty and all elements are equal (semantic
/// singleton). Shared by `ValueAst::is_ground` and the atom-field types that
/// embed a `LitSet` directly (`IsotopeAst`, `ImplicitHydrogensAst`), so they
/// avoid cloning the Vec just to delegate.
#[inline(never)]
pub(crate) fn litset_is_ground(s: &[i64]) -> bool {
    match s {
        [] => false,
        [first, rest @ ..] => rest.iter().all(|x| x == first),
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

/// Arithmetic operators for `Expr::BinOp`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ArithOp {
    Add,
    Sub,
    Mul,
    Div,
    Rem,
}

/// Relational operators for `Expr::Rel`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RelOp {
    Le,
    Ge,
    Eq,
    Lt,
    Gt,
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
    #[case::undetermined(ValueAst::Undetermined, false)]
    #[case::lit_set(ValueAst::LitSet(vec![1, 2]), false)]
    #[case::expr(ValueAst::Expr(Expr::Var("x".to_string())), false)]
    fn test_value_ast_is_ground(#[case] ast: ValueAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(ValueAst::Undetermined, 3, true)]
    #[case::lit_match(ValueAst::Lit(3), 3,  true)]
    #[case::lit_set_match(ValueAst::LitSet(vec![1, 2, 3]), 2, true)]
    #[case::expr_var(ValueAst::Expr(Expr::Var("h".to_string())), 5, true)]
    #[case::expr_lit_match(ValueAst::Expr(Expr::Lit(3)), 3, true)]
    #[case::expr_rel_match(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), 3, true)]
    #[case::expr_mem_match(ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![0, 1])), 1, true)]
    #[case::lit_no_match(ValueAst::Lit(3), 4, false)]
    #[case::expr_lit_no_match(ValueAst::Expr(Expr::Lit(3)), 4, false)]
    #[case::expr_rel_no_match(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), 0, false)]
    fn test_matches_value(#[case] pattern: ValueAst, #[case] value: i64, #[case] expected: bool) {
        assert_eq!(pattern.matches_value(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(ValueAst::Undetermined, 3, Bindings::new())]
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

    // region: simplify

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(Expr::Lit(5), Expr::Lit(5))]
    #[case::var(Expr::Var("x".into()), Expr::Var("x".into()))]
    #[case::neg_lit(Expr::Neg(Box::new(Expr::Lit(3))), Expr::Neg(Box::new(Expr::Lit(3))))]
    #[case::neg_neg_collapses(
        Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Lit(3))))),
        Expr::Lit(3),
    )]
    #[case::neg_neg_neg_collapses_to_one(
        Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Var("x".into()))))))),
        Expr::Neg(Box::new(Expr::Var("x".into()))),
    )]
    #[case::or_flattens_or_child(
        Expr::Or(vec![
            Expr::Var("a".into()),
            Expr::Or(vec![Expr::Var("b".into()), Expr::Var("c".into())]),
        ]),
        Expr::Or(vec![
            Expr::Var("a".into()),
            Expr::Var("b".into()),
            Expr::Var("c".into()),
        ]),
    )]
    #[case::and_flattens_and_child(
        Expr::And(vec![
            Expr::And(vec![Expr::Var("a".into()), Expr::Var("b".into())]),
            Expr::Var("c".into()),
        ]),
        Expr::And(vec![
            Expr::Var("a".into()),
            Expr::Var("b".into()),
            Expr::Var("c".into()),
        ]),
    )]
    #[case::or_does_not_flatten_and(
        Expr::Or(vec![
            Expr::Var("a".into()),
            Expr::And(vec![Expr::Var("b".into()), Expr::Var("c".into())]),
        ]),
        Expr::Or(vec![
            Expr::Var("a".into()),
            Expr::And(vec![Expr::Var("b".into()), Expr::Var("c".into())]),
        ]),
    )]
    #[case::recursive_into_binop(
        Expr::BinOp(
            Box::new(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Lit(2)))))),
            ArithOp::Add,
            Box::new(Expr::Lit(3)),
        ),
        Expr::BinOp(Box::new(Expr::Lit(2)), ArithOp::Add, Box::new(Expr::Lit(3))),
    )]
    #[case::recursive_into_rel(
        Expr::Rel(
            Box::new(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Var("h".into())))))),
            RelOp::Ge,
            Box::new(Expr::Lit(1)),
        ),
        Expr::Rel(Box::new(Expr::Var("h".into())), RelOp::Ge, Box::new(Expr::Lit(1))),
    )]
    fn test_expr_simplify(#[case] input: Expr, #[case] expected: Expr) {
        assert_eq!(input.simplify(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ValueAst::Lit(5), ValueAst::Lit(5))]
    #[case::undetermined(ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::lit_set(ValueAst::LitSet(vec![1, 2]), ValueAst::LitSet(vec![1, 2]))]
    #[case::expr_lit_lifts(ValueAst::Expr(Expr::Lit(5)), ValueAst::Lit(5))]
    #[case::expr_neg_lit_lifts(
        ValueAst::Expr(Expr::Neg(Box::new(Expr::Lit(7)))),
        ValueAst::Lit(-7),
    )]
    #[case::expr_neg_neg_lit_lifts(
        ValueAst::Expr(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Lit(4)))))),
        ValueAst::Lit(4),
    )]
    #[case::expr_var_stays(
        ValueAst::Expr(Expr::Var("x".into())),
        ValueAst::Expr(Expr::Var("x".into())),
    )]
    #[case::expr_neg_var_stays(
        ValueAst::Expr(Expr::Neg(Box::new(Expr::Var("x".into())))),
        ValueAst::Expr(Expr::Neg(Box::new(Expr::Var("x".into())))),
    )]
    #[case::expr_neg_lit_min_overflow_keeps_form(
        ValueAst::Expr(Expr::Neg(Box::new(Expr::Lit(i64::MIN)))),
        ValueAst::Expr(Expr::Neg(Box::new(Expr::Lit(i64::MIN)))),
    )]
    fn test_value_ast_simplify(#[case] input: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::neg_neg(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Lit(3))))))]
    #[case::nested_or(
        Expr::Or(vec![
            Expr::Or(vec![Expr::Var("a".into()), Expr::Var("b".into())]),
            Expr::Or(vec![Expr::Var("c".into()), Expr::Var("d".into())]),
        ])
    )]
    #[case::deep_neg(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Neg(Box::new(
        Expr::Neg(Box::new(Expr::Lit(1)))
    )))))))]
    fn test_expr_simplify_idempotent(#[case] input: Expr) {
        let once = input.clone().simplify();
        let twice = once.clone().simplify();
        assert_eq!(once, twice);
    }

    // endregion: simplify
}
