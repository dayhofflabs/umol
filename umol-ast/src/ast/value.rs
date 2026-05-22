//! Value AST.

use std::collections::HashMap;
use std::ops::{Add, Div, Mul, Sub};

use umol_shared::spin::SpinMultiplicity;

use super::error::EvaluationError;
use super::traits::{AsLit, Lattice};

/// Variable bindings used by [`Expr::evaluate`] and [`Expr::evaluate_bool`].
pub type Bindings = HashMap<String, i64>;

/// Integer-valued atom/bond field: undetermined (pattern wildcard), a
/// literal, a finite literal set, an arithmetic/boolean expression
/// pattern, or a named bind / reference for cross-field joint constraints.
/// Used for charge, hydrogen count, isotope mass, valence, etc.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ValueAst {
    #[default]
    Undetermined,
    Lit(i64),
    LitSet(Box<Vec<i64>>),
    Expr(Box<Expr>),
    Bind {
        id: String,
        set: Box<Vec<i64>>,
    },
    Ref(String),
}

impl ValueAst {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn lit(n: i64) -> Self {
        Self::Lit(n)
    }

    pub fn lit_set(values: Vec<i64>) -> Self {
        Self::LitSet(Box::new(values))
    }

    pub fn expr(e: Expr) -> Self {
        Self::Expr(Box::new(e))
    }

    pub fn bind(id: impl Into<String>, set: Vec<i64>) -> Self {
        Self::Bind {
            id: id.into(),
            set: Box::new(set),
        }
    }

    pub fn reference(id: impl Into<String>) -> Self {
        Self::Ref(id.into())
    }

    #[inline(never)]
    #[cold]
    fn is_ground_slow(&self) -> bool {
        match self {
            Self::LitSet(s) => litset_is_ground(s),
            Self::Expr(e) => e.is_ground(),
            Self::Bind { .. } | Self::Ref(_) => false,
            Self::Lit(_) | Self::Undetermined => unreachable!(),
        }
    }

    #[inline(never)]
    #[cold]
    fn as_lit_slow(&self) -> Option<i64> {
        match self {
            Self::LitSet(s) => litset_is_ground(s).then(|| s[0]),
            Self::Expr(e) => e.evaluate_checked(&Bindings::new()),
            Self::Bind { .. } | Self::Ref(_) => None,
            Self::Lit(_) | Self::Undetermined => unreachable!(),
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
    /// - `Expr(Expr::Var(id))` → `Ref(id)` (top-level bind reference)
    /// - `Expr(Expr::Mem(Var(id), set))` → `Bind { id, set }` (named domain
    ///   constraint at top level)
    /// - `Expr(e)` for any other shape → `Expr(e.simplify())`
    /// - `Lit` / `LitSet` / `Undetermined` / `Bind` / `Ref` → unchanged
    pub fn simplify(self) -> Self {
        match self {
            Self::Expr(e) => match e.simplify() {
                Expr::Lit(n) => Self::Lit(n),
                Expr::Neg(inner) => match *inner {
                    Expr::Lit(n) => Self::Lit(-n),
                    other => Self::Expr(Box::new(Expr::Neg(Box::new(other)))),
                },
                Expr::Var(id) => Self::Ref(id),
                Expr::Mem(inner, set) => match *inner {
                    Expr::Var(id) => Self::Bind {
                        id,
                        set: Box::new(set),
                    },
                    other => Self::Expr(Box::new(Expr::Mem(Box::new(other), set))),
                },
                other => Self::Expr(Box::new(other)),
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
            ValueAst::Bind { id, set } => {
                if set.contains(&value) {
                    let mut bindings = Bindings::new();
                    bindings.insert(id.clone(), value);
                    Some(bindings)
                } else {
                    None
                }
            }
            ValueAst::Ref(_) => None,
        }
    }
}

/// A `LitSet` is ground iff non-empty and all elements are equal (semantic
/// singleton). Shared by `ValueAst::is_ground` and the atom-field types that
/// embed a `LitSet` directly (`IsotopeAst`, ``), so they
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

impl AsLit for ValueAst {
    type Lit = i64;

    /// The single integer this value denotes when ground; `None` otherwise.
    /// Aligned with [`Lattice::is_ground`]: `is_ground() == as_lit().is_some()`.
    /// Non-destructive — does not mutate or simplify in place.
    #[inline]
    fn as_lit(&self) -> Option<i64> {
        match self {
            Self::Lit(n) => Some(*n),
            Self::Undetermined => None,
            _ => self.as_lit_slow(),
        }
    }
}

impl Lattice for ValueAst {
    #[inline]
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    /// The pattern denotes a single concrete integer. Semantic, not
    /// syntactic: `Expr` that folds to a constant is ground, and a
    /// `LitSet` of a single value (regardless of duplicates) is ground.
    ///
    /// Fast path — `Lit` and `Undetermined` dispatch with two tag compares
    /// so the common case (a fully-lowered ground molecule) doesn't pay
    /// for the `LitSet`/`Expr` logic.
    #[inline]
    fn is_ground(&self) -> bool {
        match self {
            Self::Lit(_) => true,
            Self::Undetermined => false,
            _ => self.is_ground_slow(),
        }
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::Lit(a), Self::Lit(b)) => (a == b).then_some(Self::Lit(*a)),
            (Self::Lit(a), Self::LitSet(s)) | (Self::LitSet(s), Self::Lit(a)) => {
                s.contains(a).then_some(Self::Lit(*a))
            }
            (Self::LitSet(s), Self::LitSet(t)) => {
                let intersection: Vec<i64> = s.iter().filter(|x| t.contains(x)).copied().collect();
                match intersection.len() {
                    0 => None,
                    1 => Some(Self::Lit(intersection[0])),
                    _ => Some(Self::LitSet(Box::new(intersection))),
                }
            }
            (Self::Expr(e), Self::Expr(f)) => (e == f).then(|| Self::Expr(e.clone())),
            (Self::Ref(a), Self::Ref(b)) if a == b => Some(Self::Ref(a.clone())),
            (Self::Bind { id: a, set: s }, Self::Bind { id: b, set: t }) if a == b && s == t => {
                Some(self.clone())
            }
            (Self::Expr(_), _) | (_, Self::Expr(_)) => None,
            (Self::Bind { .. } | Self::Ref(_), _) | (_, Self::Bind { .. } | Self::Ref(_)) => None,
        }
    }

    fn join(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Lit(a), Self::Lit(b)) => {
                if a == b {
                    Self::Lit(*a)
                } else {
                    Self::LitSet(Box::new(vec![*a, *b]))
                }
            }
            (Self::Lit(a), Self::LitSet(s)) => {
                let mut v: Vec<i64> = Vec::with_capacity(s.len() + 1);
                v.push(*a);
                for &x in s.iter() {
                    if x != *a {
                        v.push(x);
                    }
                }
                if v.len() == 1 {
                    Self::Lit(v[0])
                } else {
                    Self::LitSet(Box::new(v))
                }
            }
            (Self::LitSet(s), Self::Lit(a)) => {
                let mut v: Vec<i64> = s.to_vec();
                if !v.contains(a) {
                    v.push(*a);
                }
                if v.len() == 1 {
                    Self::Lit(v[0])
                } else {
                    Self::LitSet(Box::new(v))
                }
            }
            (Self::LitSet(s), Self::LitSet(t)) => {
                let mut v: Vec<i64> = s.to_vec();
                for &x in t.iter() {
                    if !v.contains(&x) {
                        v.push(x);
                    }
                }
                if v.len() == 1 {
                    Self::Lit(v[0])
                } else {
                    Self::LitSet(Box::new(v))
                }
            }
            (Self::Expr(e), Self::Expr(f)) if e == f => Self::Expr(e.clone()),
            (Self::Ref(a), Self::Ref(b)) if a == b => Self::Ref(a.clone()),
            (Self::Bind { id: a, set: s }, Self::Bind { id: b, set: t }) if a == b && s == t => {
                self.clone()
            }
            _ => Self::Undetermined,
        }
    }

    /// Pattern matches target iff every integer target admits is also
    /// admitted by pattern (superset semantics). `Expr` and `Ref` targets
    /// cannot be certified generically and are rejected. `Bind` is treated
    /// as a `LitSet` of its admissible values on either side.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Ref(_), _) | (_, Self::Ref(_)) => false,
            (_, Self::Expr(_)) => false,
            (pattern, Self::Lit(n)) => pattern.matches_value(*n),
            (pattern, Self::LitSet(ns) | Self::Bind { set: ns, .. }) => {
                ns.iter().all(|n| pattern.matches_value(*n))
            }
        }
    }
}

// Arithmetic on `ValueAst` propagates `Undetermined`. Every binop has impls
// for all four `(owned|ref) × (owned|ref)` combinations; the owned forms
// delegate to the ref-ref form so the match is written once. Each binop
// additionally accepts a bare `i64` on either side.
macro_rules! impl_value_binop {
    ($Op:ident, $op:ident, $lit_op:tt) => {
        impl $Op<&ValueAst> for &ValueAst {
            type Output = ValueAst;
            fn $op(self, rhs: &ValueAst) -> ValueAst {
                match (self, rhs) {
                    (ValueAst::Lit(a), ValueAst::Lit(b)) => ValueAst::Lit(a $lit_op b),
                    _ => ValueAst::Undetermined,
                }
            }
        }
        impl $Op<ValueAst> for &ValueAst {
            type Output = ValueAst;
            fn $op(self, rhs: ValueAst) -> ValueAst { self.$op(&rhs) }
        }
        impl $Op<&ValueAst> for ValueAst {
            type Output = ValueAst;
            fn $op(self, rhs: &ValueAst) -> ValueAst { (&self).$op(rhs) }
        }
        impl $Op<ValueAst> for ValueAst {
            type Output = ValueAst;
            fn $op(self, rhs: ValueAst) -> ValueAst { (&self).$op(&rhs) }
        }
        impl $Op<i64> for &ValueAst {
            type Output = ValueAst;
            fn $op(self, rhs: i64) -> ValueAst { self.$op(&ValueAst::Lit(rhs)) }
        }
        impl $Op<i64> for ValueAst {
            type Output = ValueAst;
            fn $op(self, rhs: i64) -> ValueAst { (&self).$op(&ValueAst::Lit(rhs)) }
        }
        impl $Op<&ValueAst> for i64 {
            type Output = ValueAst;
            fn $op(self, rhs: &ValueAst) -> ValueAst { (&ValueAst::Lit(self)).$op(rhs) }
        }
        impl $Op<ValueAst> for i64 {
            type Output = ValueAst;
            fn $op(self, rhs: ValueAst) -> ValueAst { (&ValueAst::Lit(self)).$op(&rhs) }
        }
    };
}

impl_value_binop!(Add, add, +);
impl_value_binop!(Sub, sub, -);
impl_value_binop!(Mul, mul, *);
impl_value_binop!(Div, div, /);

impl From<i64> for ValueAst {
    fn from(value: i64) -> Self {
        Self::Lit(value)
    }
}

impl From<SpinMultiplicity> for ValueAst {
    fn from(m: SpinMultiplicity) -> Self {
        Self::Lit(u8::from(m) as i64)
    }
}

impl From<Vec<i64>> for ValueAst {
    fn from(values: Vec<i64>) -> Self {
        Self::LitSet(Box::new(values))
    }
}

impl From<Expr> for ValueAst {
    fn from(e: Expr) -> Self {
        Self::Expr(Box::new(e))
    }
}

/// Arithmetic / boolean expression tree over `ValueAst`. Captures
/// atom/bond field constraints that can't be expressed as a literal or
/// literal set, including bound variables (`Var`), membership tests
/// (`Mem`), relational comparisons (`Rel`), and boolean combinators.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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
            Expr::Lit(n) if n < 0 => Expr::Neg(Box::new(Expr::Lit(-n))),
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

    /// Arithmetic evaluation. Returns `None` for free variables and type
    /// mismatch (boolean-domain Expr). Uses standard Rust arithmetic — `i64`
    /// overflow follows debug-panic / release-wrap semantics, and division
    /// or remainder by zero panics. Intended as the foundation of
    /// [`Expr::is_ground`]; for error-reporting callers use [`Expr::evaluate`].
    pub fn evaluate_checked(&self, vars: &Bindings) -> Option<i64> {
        match self {
            Expr::Lit(n) => Some(*n),
            Expr::Var(name) => vars.get(name).copied(),
            Expr::Neg(e) => Some(-e.evaluate_checked(vars)?),
            Expr::BinOp(l, op, r) => {
                let l = l.evaluate_checked(vars)?;
                let r = r.evaluate_checked(vars)?;
                Some(match op {
                    ArithOp::Add => l + r,
                    ArithOp::Sub => l - r,
                    ArithOp::Mul => l * r,
                    ArithOp::Div => l / r,
                    ArithOp::Rem => l % r,
                })
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

/// Arithmetic operators for `Expr::BinOp`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ArithOp {
    Add,
    Sub,
    Mul,
    Div,
    Rem,
}

/// Relational operators for `Expr::Rel`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ValueAst::Lit(3), Some(3))]
    #[case::lit_zero(ValueAst::Lit(0), Some(0))]
    #[case::lit_neg(ValueAst::Lit(-5), Some(-5))]
    #[case::undetermined(ValueAst::Undetermined, None)]
    #[case::lit_set_empty(ValueAst::LitSet(Box::default()), None)]
    #[case::lit_set_singleton(ValueAst::LitSet(Box::new(vec![7])), Some(7))]
    #[case::lit_set_all_equal(ValueAst::LitSet(Box::new(vec![3, 3, 3])), Some(3))]
    #[case::lit_set_multi(ValueAst::LitSet(Box::new(vec![1, 2])), None)]
    #[case::expr_lit(ValueAst::Expr(Box::new(Expr::Lit(5))), Some(5))]
    #[case::expr_neg_lit(ValueAst::Expr(Box::new(Expr::Neg(Box::new(Expr::Lit(3))))), Some(-3))]
    #[case::expr_const_add(
        ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(2)), ArithOp::Add, Box::new(Expr::Lit(3))))),
        Some(5),
    )]
    #[case::expr_var(ValueAst::Expr(Box::new(Expr::Var("x".to_string()))), None)]
    #[case::expr_boolean(
        ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(1))))),
        None,
    )]
    #[case::bind(ValueAst::bind("n", vec![1, 2]), None)]
    #[case::bind_singleton(ValueAst::bind("n", vec![3]), None)]
    #[case::reference(ValueAst::reference("n"), None)]
    fn test_value_ast_literal_and_is_ground(
        #[case] ast: ValueAst,
        #[case] expected_literal: Option<i64>,
    ) {
        assert_eq!(ast.as_lit(), expected_literal);
        assert_eq!(ast.is_ground(), expected_literal.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(ValueAst::Undetermined, 3, true)]
    #[case::lit_match(ValueAst::Lit(3), 3,  true)]
    #[case::lit_set_match(ValueAst::LitSet(Box::new(vec![1, 2, 3])), 2, true)]
    #[case::expr_var(ValueAst::Expr(Box::new(Expr::Var("h".to_string()))), 5, true)]
    #[case::expr_lit_match(ValueAst::Expr(Box::new(Expr::Lit(3))), 3, true)]
    #[case::expr_rel_match(ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1))))), 3, true)]
    #[case::expr_mem_match(ValueAst::Expr(Box::new(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![0, 1]))), 1, true)]
    #[case::bind_in_set(ValueAst::bind("n", vec![1, 2, 3]), 2, true)]
    #[case::lit_no_match(ValueAst::Lit(3), 4, false)]
    #[case::expr_lit_no_match(ValueAst::Expr(Box::new(Expr::Lit(3))), 4, false)]
    #[case::expr_rel_no_match(ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1))))), 0, false)]
    #[case::bind_not_in_set(ValueAst::bind("n", vec![1, 2]), 5, false)]
    #[case::reference_no_capture(ValueAst::reference("n"), 3, false)]
    fn test_matches_value(#[case] pattern: ValueAst, #[case] value: i64, #[case] expected: bool) {
        assert_eq!(pattern.matches_value(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(ValueAst::Undetermined, 3, Bindings::new())]
    #[case::lit_match(ValueAst::Lit(3), 3, Bindings::new())]
    #[case::lit_set_match(ValueAst::LitSet(Box::new(vec![1, 2, 3])), 2, Bindings::new())]
    #[case::expr_var(ValueAst::Expr(Box::new(Expr::Var("h".to_string()))), 5, Bindings::from([("h".to_string(), 5)]))]
    #[case::expr_lit_match(ValueAst::Expr(Box::new(Expr::Lit(3))), 3, Bindings::new())]
    #[case::expr_rel_match(ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1))))), 3, Bindings::from([("h".to_string(), 3)]))]
    #[case::expr_mem_match(ValueAst::Expr(Box::new(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![0, 1]))), 1, Bindings::from([("h".to_string(), 1)]))]
    #[case::bind_in_set(ValueAst::bind("n", vec![1, 2, 3]), 2, Bindings::from([("n".to_string(), 2)]))]
    fn test_capture(#[case] pattern: ValueAst, #[case] value: i64, #[case] expected: Bindings) {
        assert_eq!(pattern.capture(value), Some(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_no_match(ValueAst::Lit(3), 4)]
    #[case::expr_lit_no_match(ValueAst::Expr(Box::new(Expr::Lit(3))), 4)]
    #[case::expr_rel_no_match(ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1))))), 0)]
    #[case::bind_not_in_set(ValueAst::bind("n", vec![1, 2]), 5)]
    #[case::reference(ValueAst::reference("n"), 3)]
    fn test_capture_no_match(#[case] pattern: ValueAst, #[case] value: i64) {
        assert_eq!(pattern.capture(value), None);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(Expr::Lit(5), Expr::Lit(5))]
    #[case::var(Expr::Var("x".into()), Expr::Var("x".into()))]
    #[case::neg_lit(Expr::Neg(Box::new(Expr::Lit(3))), Expr::Neg(Box::new(Expr::Lit(3))))]
    #[case::neg_neg_collapses(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Lit(3))))), Expr::Lit(3))]
    #[case::neg_neg_neg_collapses_to_one(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Var("x".into()))))))),
        Expr::Neg(Box::new(Expr::Var("x".into()))))]
    #[case::or_flattens_or_child(Expr::Or(vec![Expr::Var("a".into()), Expr::Or(vec![Expr::Var("b".into()), Expr::Var("c".into())])]),
        Expr::Or(vec![Expr::Var("a".into()), Expr::Var("b".into()), Expr::Var("c".into())]))]
    #[case::and_flattens_and_child(Expr::And(vec![Expr::And(vec![Expr::Var("a".into()), Expr::Var("b".into())]), Expr::Var("c".into())]),
        Expr::And(vec![Expr::Var("a".into()), Expr::Var("b".into()), Expr::Var("c".into())]))]
    #[case::or_does_not_flatten_and(Expr::Or(vec![Expr::Var("a".into()), Expr::And(vec![Expr::Var("b".into()), Expr::Var("c".into())])]),
        Expr::Or(vec![Expr::Var("a".into()), Expr::And(vec![Expr::Var("b".into()), Expr::Var("c".into())])]))]
    #[case::recursive_into_binop(Expr::BinOp(Box::new(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Lit(2)))))), ArithOp::Add, Box::new(Expr::Lit(3))),
        Expr::BinOp(Box::new(Expr::Lit(2)), ArithOp::Add, Box::new(Expr::Lit(3))))]
    #[case::recursive_into_rel(Expr::Rel(Box::new(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Var("h".into())))))), RelOp::Ge, Box::new(Expr::Lit(1))),
        Expr::Rel(Box::new(Expr::Var("h".into())), RelOp::Ge, Box::new(Expr::Lit(1))))]
    fn test_expr_simplify(#[case] input: Expr, #[case] expected: Expr) {
        assert_eq!(input.simplify(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ValueAst::Lit(5), ValueAst::Lit(5))]
    #[case::undetermined(ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::lit_set(ValueAst::LitSet(Box::new(vec![1, 2])), ValueAst::LitSet(Box::new(vec![1, 2])))]
    #[case::ref_stays(ValueAst::reference("h"), ValueAst::reference("h"))]
    #[case::bind_stays(ValueAst::bind("h", vec![1, 2]), ValueAst::bind("h", vec![1, 2]))]
    #[case::expr_lit_lifts(ValueAst::Expr(Box::new(Expr::Lit(5))), ValueAst::Lit(5))]
    #[case::expr_neg_lit_lifts(ValueAst::Expr(Box::new(Expr::Neg(Box::new(Expr::Lit(7))))), ValueAst::Lit(-7))]
    #[case::expr_neg_neg_lit_lifts(ValueAst::Expr(Box::new(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Lit(4))))))), ValueAst::Lit(4))]
    #[case::expr_var_lifts_to_ref(ValueAst::Expr(Box::new(Expr::Var("x".into()))), ValueAst::reference("x"))]
    #[case::expr_mem_var_lifts_to_bind(
        ValueAst::Expr(Box::new(Expr::Mem(Box::new(Expr::Var("h".into())), vec![1, 2, 3]))),
        ValueAst::bind("h", vec![1, 2, 3]),
    )]
    #[case::expr_neg_var_stays(ValueAst::Expr(Box::new(Expr::Neg(Box::new(Expr::Var("x".into()))))),
        ValueAst::Expr(Box::new(Expr::Neg(Box::new(Expr::Var("x".into()))))))]
    #[case::expr_binop_var_stays(
        ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Var("h".into())), ArithOp::Add, Box::new(Expr::Lit(1))))),
        ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Var("h".into())), ArithOp::Add, Box::new(Expr::Lit(1))))),
    )]
    #[case::expr_mem_compound_first_stays(
        ValueAst::Expr(Box::new(Expr::Mem(
            Box::new(Expr::BinOp(Box::new(Expr::Var("h".into())), ArithOp::Add, Box::new(Expr::Lit(1)))),
            vec![2, 3],
        ))),
        ValueAst::Expr(Box::new(Expr::Mem(
            Box::new(Expr::BinOp(Box::new(Expr::Var("h".into())), ArithOp::Add, Box::new(Expr::Lit(1)))),
            vec![2, 3],
        ))),
    )]
    fn test_value_ast_simplify(#[case] input: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::expr_var(ValueAst::Expr(Box::new(Expr::Var("h".into()))))]
    #[case::expr_mem_var(ValueAst::Expr(Box::new(Expr::Mem(Box::new(Expr::Var("h".into())), vec![1, 2]))))]
    #[case::ref_(ValueAst::reference("h"))]
    #[case::bind(ValueAst::bind("h", vec![1, 2, 3]))]
    #[case::expr_binop(ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Var("h".into())), ArithOp::Add, Box::new(Expr::Lit(1))))))]
    fn test_value_ast_simplify_idempotent(#[case] input: ValueAst) {
        let once = input.clone().simplify();
        let twice = once.clone().simplify();
        assert_eq!(once, twice);
    }

    #[rstest]
    #[case::neg_neg(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Lit(3))))))]
    #[case::nested_or(Expr::Or(vec![Expr::Or(vec![Expr::Var("a".into()), Expr::Var("b".into())]),
        Expr::Or(vec![Expr::Var("c".into()), Expr::Var("d".into())])]))]
    #[case::deep_neg(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Neg(Box::new(
        Expr::Neg(Box::new(Expr::Lit(1)))
    )))))))]
    fn test_expr_simplify_idempotent(#[case] input: Expr) {
        let once = input.clone().simplify();
        let twice = once.clone().simplify();
        assert_eq!(once, twice);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(ValueAst::Lit(2), ValueAst::Lit(3), ValueAst::Lit(5))]
    #[case::lit_negative(ValueAst::Lit(1), ValueAst::Lit(-4), ValueAst::Lit(-3))]
    #[case::lit_undetermined(ValueAst::Lit(2), ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::undetermined_lit(ValueAst::Undetermined, ValueAst::Lit(2), ValueAst::Undetermined)]
    #[case::litset_lit(ValueAst::LitSet(Box::new(vec![1, 2])), ValueAst::Lit(3), ValueAst::Undetermined)]
    #[case::lit_expr(ValueAst::Lit(2), ValueAst::Expr(Box::new(Expr::Var("x".into()))), ValueAst::Undetermined)]
    fn test_value_ast_add(#[case] lhs: ValueAst, #[case] rhs: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(lhs + rhs, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(ValueAst::Lit(5), ValueAst::Lit(3), ValueAst::Lit(2))]
    #[case::lit_negative_result(ValueAst::Lit(1), ValueAst::Lit(4), ValueAst::Lit(-3))]
    #[case::lit_undetermined(ValueAst::Lit(5), ValueAst::Undetermined, ValueAst::Undetermined)]
    fn test_value_ast_sub(#[case] lhs: ValueAst, #[case] rhs: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(lhs - rhs, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(ValueAst::Lit(4), ValueAst::Lit(3), ValueAst::Lit(12))]
    #[case::lit_undetermined(ValueAst::Lit(4), ValueAst::Undetermined, ValueAst::Undetermined)]
    fn test_value_ast_mul(#[case] lhs: ValueAst, #[case] rhs: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(lhs * rhs, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(ValueAst::Lit(10), ValueAst::Lit(3), ValueAst::Lit(3))]
    #[case::lit_undetermined(ValueAst::Lit(10), ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::undetermined_lit_zero(ValueAst::Undetermined, ValueAst::Lit(0), ValueAst::Undetermined)]
    fn test_value_ast_div(#[case] lhs: ValueAst, #[case] rhs: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(lhs / rhs, expected);
    }

    #[rstest]
    #[should_panic]
    fn test_value_ast_div_by_zero_panics() {
        let _ = ValueAst::Lit(5) / ValueAst::Lit(0);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_und(ValueAst::Undetermined, ValueAst::Undetermined, Some(ValueAst::Undetermined))]
    #[case::und_lit(ValueAst::Undetermined, ValueAst::Lit(3), Some(ValueAst::Lit(3)))]
    #[case::lit_und(ValueAst::Lit(3), ValueAst::Undetermined, Some(ValueAst::Lit(3)))]
    #[case::lit_lit_eq(ValueAst::Lit(3), ValueAst::Lit(3), Some(ValueAst::Lit(3)))]
    #[case::lit_lit_neq(ValueAst::Lit(3), ValueAst::Lit(4), None)]
    #[case::lit_litset_in(ValueAst::Lit(2), ValueAst::LitSet(Box::new(vec![1, 2, 3])), Some(ValueAst::Lit(2)))]
    #[case::lit_litset_out(ValueAst::Lit(5), ValueAst::LitSet(Box::new(vec![1, 2, 3])), None)]
    #[case::litset_lit_in(ValueAst::LitSet(Box::new(vec![1, 2, 3])), ValueAst::Lit(2), Some(ValueAst::Lit(2)))]
    #[case::litset_litset_multi(ValueAst::LitSet(Box::new(vec![1, 2, 3])), ValueAst::LitSet(Box::new(vec![2, 3, 4])),
        Some(ValueAst::LitSet(Box::new(vec![2, 3]))))]
    #[case::litset_litset_singleton(ValueAst::LitSet(Box::new(vec![1, 2])), ValueAst::LitSet(Box::new(vec![2, 3])), Some(ValueAst::Lit(2)))]
    #[case::litset_litset_empty(ValueAst::LitSet(Box::new(vec![1, 2])), ValueAst::LitSet(Box::new(vec![3, 4])), None)]
    #[case::expr_expr_eq(ValueAst::Expr(Box::new(Expr::Lit(5))), ValueAst::Expr(Box::new(Expr::Lit(5))), Some(ValueAst::Expr(Box::new(Expr::Lit(5)))))]
    #[case::expr_expr_neq(ValueAst::Expr(Box::new(Expr::Lit(5))), ValueAst::Expr(Box::new(Expr::Lit(6))), None)]
    #[case::expr_lit(ValueAst::Expr(Box::new(Expr::Var("x".into()))), ValueAst::Lit(5), None)]
    #[case::expr_und(ValueAst::Expr(Box::new(Expr::Var("x".into()))), ValueAst::Undetermined, Some(ValueAst::Expr(Box::new(Expr::Var("x".into())))))]
    #[case::bind_und(ValueAst::bind("n", vec![1, 2]), ValueAst::Undetermined, Some(ValueAst::bind("n", vec![1, 2])))]
    #[case::bind_bind_eq(ValueAst::bind("n", vec![1, 2]), ValueAst::bind("n", vec![1, 2]), Some(ValueAst::bind("n", vec![1, 2])))]
    #[case::bind_bind_id_neq(ValueAst::bind("n", vec![1, 2]), ValueAst::bind("m", vec![1, 2]), None)]
    #[case::bind_bind_set_neq(ValueAst::bind("n", vec![1, 2]), ValueAst::bind("n", vec![3, 4]), None)]
    #[case::bind_lit(ValueAst::bind("n", vec![1, 2]), ValueAst::Lit(1), None)]
    #[case::ref_ref_eq(ValueAst::reference("n"), ValueAst::reference("n"), Some(ValueAst::reference("n")))]
    #[case::ref_ref_neq(ValueAst::reference("n"), ValueAst::reference("m"), None)]
    #[case::ref_und(ValueAst::reference("n"), ValueAst::Undetermined, Some(ValueAst::reference("n")))]
    #[case::ref_lit(ValueAst::reference("n"), ValueAst::Lit(1), None)]
    #[case::bind_ref(ValueAst::bind("n", vec![1, 2]), ValueAst::reference("n"), None)]
    fn test_value_ast_meet(
        #[case] a: ValueAst,
        #[case] b: ValueAst,
        #[case] expected: Option<ValueAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_und(ValueAst::Undetermined, ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::und_lit(ValueAst::Undetermined, ValueAst::Lit(3), ValueAst::Undetermined)]
    #[case::lit_und(ValueAst::Lit(3), ValueAst::Undetermined, ValueAst::Undetermined)]
    #[case::lit_lit_eq(ValueAst::Lit(3), ValueAst::Lit(3), ValueAst::Lit(3))]
    #[case::lit_lit_neq(ValueAst::Lit(3), ValueAst::Lit(4), ValueAst::LitSet(Box::new(vec![3, 4])))]
    #[case::lit_litset_in(ValueAst::Lit(2), ValueAst::LitSet(Box::new(vec![1, 2, 3])), ValueAst::LitSet(Box::new(vec![2, 1, 3])))]
    #[case::lit_litset_out(ValueAst::Lit(5), ValueAst::LitSet(Box::new(vec![1, 2, 3])), ValueAst::LitSet(Box::new(vec![5, 1, 2, 3])))]
    #[case::litset_lit_in(ValueAst::LitSet(Box::new(vec![1, 2, 3])), ValueAst::Lit(2), ValueAst::LitSet(Box::new(vec![1, 2, 3])))]
    #[case::litset_lit_out(ValueAst::LitSet(Box::new(vec![1, 2, 3])), ValueAst::Lit(4), ValueAst::LitSet(Box::new(vec![1, 2, 3, 4])))]
    #[case::litset_litset_overlap(ValueAst::LitSet(Box::new(vec![1, 2])), ValueAst::LitSet(Box::new(vec![2, 3])),
        ValueAst::LitSet(Box::new(vec![1, 2, 3])))]
    #[case::expr_expr_eq(ValueAst::Expr(Box::new(Expr::Lit(5))), ValueAst::Expr(Box::new(Expr::Lit(5))), ValueAst::Expr(Box::new(Expr::Lit(5))))]
    #[case::expr_expr_neq(ValueAst::Expr(Box::new(Expr::Lit(5))), ValueAst::Expr(Box::new(Expr::Lit(6))), ValueAst::Undetermined)]
    #[case::expr_lit(ValueAst::Expr(Box::new(Expr::Var("x".into()))), ValueAst::Lit(5), ValueAst::Undetermined)]
    #[case::bind_bind_eq(ValueAst::bind("n", vec![1, 2]), ValueAst::bind("n", vec![1, 2]), ValueAst::bind("n", vec![1, 2]))]
    #[case::bind_bind_neq(ValueAst::bind("n", vec![1, 2]), ValueAst::bind("m", vec![1, 2]), ValueAst::Undetermined)]
    #[case::bind_lit(ValueAst::bind("n", vec![1, 2]), ValueAst::Lit(1), ValueAst::Undetermined)]
    #[case::ref_ref_eq(ValueAst::reference("n"), ValueAst::reference("n"), ValueAst::reference("n"))]
    #[case::ref_ref_neq(ValueAst::reference("n"), ValueAst::reference("m"), ValueAst::Undetermined)]
    #[case::ref_lit(ValueAst::reference("n"), ValueAst::Lit(1), ValueAst::Undetermined)]
    fn test_value_ast_join(
        #[case] a: ValueAst,
        #[case] b: ValueAst,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rstest]
    #[case::no_change(ValueAst::Lit(3), ValueAst::Lit(3), false, ValueAst::Lit(3))]
    #[case::tighten(ValueAst::Undetermined, ValueAst::Lit(3), true, ValueAst::Lit(3))]
    #[case::incompatible(ValueAst::Lit(3), ValueAst::Lit(4), false, ValueAst::Lit(3))]
    fn test_value_ast_narrow_from(
        #[case] mut target: ValueAst,
        #[case] source: ValueAst,
        #[case] expected_changed: bool,
        #[case] expected_after: ValueAst,
    ) {
        let changed = target.narrow_from(&source);
        assert_eq!(changed, expected_changed);
        assert_eq!(target, expected_after);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::no_change(ValueAst::Lit(3), ValueAst::Lit(3), false, ValueAst::Lit(3))]
    #[case::widen_to_set(ValueAst::Lit(3), ValueAst::Lit(4), true, ValueAst::LitSet(Box::new(vec![3, 4])))]
    #[case::widen_to_top(ValueAst::Lit(3), ValueAst::Undetermined, true, ValueAst::Undetermined)]
    fn test_value_ast_widen_with(
        #[case] mut target: ValueAst,
        #[case] source: ValueAst,
        #[case] expected_changed: bool,
        #[case] expected_after: ValueAst,
    ) {
        let changed = target.widen_with(&source);
        assert_eq!(changed, expected_changed);
        assert_eq!(target, expected_after);
    }
}
