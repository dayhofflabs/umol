//! Numeric form and its expression language.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::ops::{Add, Div, Mul, Sub};

use umol_chem::spin::SpinMultiplicity;

use super::error::{Contradiction, NoJoin};
use super::operators::{MemOp, RelOp};
use super::traits::{AsLit, Canonicalize, Lattice};

/// Integer-valued atom/bond field: undetermined (pattern wildcard), a literal,
/// a finite literal set, an arithmetic expression over variables, or a predicate
/// expression constraining the field. Used for charge, hydrogen count, isotope
/// mass, valence, bond order, etc.
///
/// Equality is **lazy**: derived `Eq`/`Hash`/`Ord` are structural ("same
/// tree"); semantic equality is `Canonicalize::canonical_eq`, which compares canonical
/// forms.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NumForm {
    #[default]
    Undetermined,
    Lit(i64),
    LitSet(Box<BTreeSet<i64>>),
    RangeFrom(i64),
    RangeTo(i64),
    ArithExpr(Box<ArithExpr>),
    PredExpr(Box<PredExpr>),
}

/// Arithmetic expression over `i64`, the value-sort half of the field grammar. `Sum`
/// and `Product` are n-ary (associative + commutative by construction);
/// subtraction lowers to `Sum([a, Neg(b)])`; `Div`/`Rem` stay binary. A ground
/// expression folds to a `Lit` (or `Neg(Lit)`), which `NumForm` lifts to `Lit`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ArithExpr {
    Lit(i64),
    Var(String),
    Neg(Box<ArithExpr>),
    Sum(Vec<ArithExpr>),
    Product(Vec<ArithExpr>),
    Div(Box<ArithExpr>, Box<ArithExpr>),
    Rem(Box<ArithExpr>, Box<ArithExpr>),
}

/// Predicate expression over arithmetic expressions, the constraint-sort half of the field grammar.
/// `Rel` and `Mem` operators are negation-closed (`RelOp` has `Ne`, `MemOp` has
/// `NotIn`), so canonicalization eliminates `Not` entirely — it survives only as
/// faithful parser input for `!`. ⊤/⊥ are not variants: a predicate that decides
/// is lifted by `NumForm` to `Undetermined` / `Err(Contradiction)`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PredExpr {
    Rel(ArithExpr, RelOp, ArithExpr),
    Mem(ArithExpr, MemOp, BTreeSet<i64>),
    Not(Box<PredExpr>),
    And(Vec<PredExpr>),
    Or(Vec<PredExpr>),
}

impl NumForm {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn lit(n: i64) -> Self {
        Self::Lit(n)
    }

    pub fn lit_set<I: IntoIterator<Item = i64>>(values: I) -> Self {
        Self::LitSet(Box::new(values.into_iter().collect()))
    }

    pub fn var(name: impl Into<String>) -> Self {
        Self::ArithExpr(Box::new(ArithExpr::Var(name.into())))
    }

    pub fn range_from(bound: i64) -> Self {
        Self::RangeFrom(bound)
    }

    pub fn range_to(bound: i64) -> Self {
        Self::RangeTo(bound)
    }

    pub fn arith_expr(expression: ArithExpr) -> Self {
        Self::ArithExpr(Box::new(expression))
    }

    pub fn pred_expr(expression: PredExpr) -> Self {
        Self::PredExpr(Box::new(expression))
    }
}

impl From<i64> for NumForm {
    fn from(value: i64) -> Self {
        Self::Lit(value)
    }
}

impl From<SpinMultiplicity> for NumForm {
    fn from(m: SpinMultiplicity) -> Self {
        Self::Lit(u8::from(m) as i64)
    }
}

impl From<Vec<i64>> for NumForm {
    fn from(values: Vec<i64>) -> Self {
        Self::lit_set(values)
    }
}

impl From<BTreeSet<i64>> for NumForm {
    fn from(values: BTreeSet<i64>) -> Self {
        Self::LitSet(Box::new(values))
    }
}

impl Canonicalize for NumForm {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            NumForm::Undetermined => NumForm::Undetermined,
            NumForm::Lit(n) => NumForm::Lit(n),
            NumForm::LitSet(set) => lift_set(*set)?,
            NumForm::RangeFrom(n) => NumForm::RangeFrom(n),
            NumForm::RangeTo(n) => NumForm::RangeTo(n),
            NumForm::ArithExpr(term) => lift_arith_expr(canon_arith_expr(*term)),
            NumForm::PredExpr(predicate) => match reduce_pred_expr(*predicate) {
                PredicateReduction::Top => NumForm::Undetermined,
                PredicateReduction::Bottom => return Err(Contradiction),
                PredicateReduction::Predicate(p) => NumForm::PredExpr(Box::new(p)),
            },
        })
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            NumForm::Undetermined | NumForm::Lit(_) => Ok(Cow::Borrowed(self)),
            _ => Ok(Cow::Owned(self.clone().canonicalize()?)),
        }
    }
}

impl Canonicalize for ArithExpr {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(canon_arith_expr(self))
    }
}

/// Lift a canonical set: empty is unsatisfiable, a singleton is a `Lit`.
fn lift_set(set: BTreeSet<i64>) -> Result<NumForm, Contradiction> {
    match set.len() {
        0 => Err(Contradiction),
        1 => Ok(NumForm::Lit(*set.iter().next().unwrap())),
        _ => Ok(NumForm::LitSet(Box::new(set))),
    }
}

/// Lift a canonical arithmetic expression: a ground expression becomes `Lit`.
fn lift_arith_expr(expression: ArithExpr) -> NumForm {
    match arith_expr_const(&expression) {
        Some(n) => NumForm::Lit(n),
        None => NumForm::ArithExpr(Box::new(expression)),
    }
}

/// Canonical literal expression: the sign lives in `Neg`, so `Lit` is always ≥ 0.
fn arith_expr_lit(n: i64) -> ArithExpr {
    if n < 0 {
        ArithExpr::Neg(Box::new(ArithExpr::Lit(-n)))
    } else {
        ArithExpr::Lit(n)
    }
}

/// The integer a canonical expression denotes, if ground.
fn arith_expr_const(expression: &ArithExpr) -> Option<i64> {
    match expression {
        ArithExpr::Lit(n) => Some(*n),
        ArithExpr::Neg(inner) => match inner.as_ref() {
            ArithExpr::Lit(n) => Some(-n),
            _ => None,
        },
        _ => None,
    }
}

fn canon_arith_expr(expression: ArithExpr) -> ArithExpr {
    match expression {
        ArithExpr::Lit(n) => arith_expr_lit(n),
        ArithExpr::Var(_) => expression,
        ArithExpr::Neg(inner) => canon_neg(canon_arith_expr(*inner)),
        ArithExpr::Sum(operands) => canon_sum(operands),
        ArithExpr::Product(operands) => canon_product(operands),
        ArithExpr::Div(a, b) => canon_div_rem(canon_arith_expr(*a), canon_arith_expr(*b), false),
        ArithExpr::Rem(a, b) => canon_div_rem(canon_arith_expr(*a), canon_arith_expr(*b), true),
    }
}

/// `inner` is already canonical.
fn canon_neg(inner: ArithExpr) -> ArithExpr {
    match inner {
        ArithExpr::Neg(grand) => *grand,
        ArithExpr::Lit(0) => ArithExpr::Lit(0),
        other => ArithExpr::Neg(Box::new(other)),
    }
}

fn canon_sum(operands: Vec<ArithExpr>) -> ArithExpr {
    let mut terms = Vec::new();
    let mut constant: i64 = 0;
    flatten_sum(operands, &mut terms, &mut constant);
    if constant != 0 {
        terms.push(arith_expr_lit(constant));
    }
    terms.sort();
    match terms.len() {
        0 => ArithExpr::Lit(0),
        1 => terms.pop().unwrap(),
        _ => ArithExpr::Sum(terms),
    }
}

fn flatten_sum(operands: Vec<ArithExpr>, terms: &mut Vec<ArithExpr>, constant: &mut i64) {
    for operand in operands {
        match canon_arith_expr(operand) {
            ArithExpr::Sum(inner) => flatten_sum(inner, terms, constant),
            other => match arith_expr_const(&other) {
                Some(n) => *constant += n,
                None => terms.push(other),
            },
        }
    }
}

fn canon_product(operands: Vec<ArithExpr>) -> ArithExpr {
    let mut terms = Vec::new();
    let mut constant: i64 = 1;
    flatten_product(operands, &mut terms, &mut constant);
    if constant == 0 {
        return ArithExpr::Lit(0);
    }
    if constant != 1 {
        terms.push(arith_expr_lit(constant));
    }
    terms.sort();
    match terms.len() {
        0 => ArithExpr::Lit(1),
        1 => terms.pop().unwrap(),
        _ => ArithExpr::Product(terms),
    }
}

fn flatten_product(operands: Vec<ArithExpr>, terms: &mut Vec<ArithExpr>, constant: &mut i64) {
    for operand in operands {
        match canon_arith_expr(operand) {
            ArithExpr::Product(inner) => flatten_product(inner, terms, constant),
            other => match arith_expr_const(&other) {
                Some(n) => *constant *= n,
                None => terms.push(other),
            },
        }
    }
}

/// `a`/`b` already canonical; folds a ground `(Lit) op (Lit)` when divisor ≠ 0.
fn canon_div_rem(a: ArithExpr, b: ArithExpr, is_rem: bool) -> ArithExpr {
    if let (Some(x), Some(y)) = (arith_expr_const(&a), arith_expr_const(&b)) {
        if y != 0 {
            return arith_expr_lit(if is_rem { x % y } else { x / y });
        }
    }
    let (a, b) = (Box::new(a), Box::new(b));
    if is_rem {
        ArithExpr::Rem(a, b)
    } else {
        ArithExpr::Div(a, b)
    }
}

/// Predicate canonical form, threading ⊤/⊥ that no `PredExpr` variant can
/// hold. Lifted to `NumForm` (⊤ → `Undetermined`, ⊥ → `Err`).
enum PredicateReduction {
    Top,
    Bottom,
    Predicate(PredExpr),
}

fn reduce_pred_expr(predicate: PredExpr) -> PredicateReduction {
    match predicate {
        PredExpr::Rel(a, op, b) => canon_rel(canon_arith_expr(a), op, canon_arith_expr(b)),
        PredExpr::Mem(e, op, set) => canon_mem(canon_arith_expr(e), op, set),
        PredExpr::Not(inner) => negate(reduce_pred_expr(*inner)),
        PredExpr::And(operands) => canon_junction(operands, true),
        PredExpr::Or(operands) => canon_junction(operands, false),
    }
}

/// `a`/`b` already canonical.
fn canon_rel(a: ArithExpr, op: RelOp, b: ArithExpr) -> PredicateReduction {
    if let (Some(x), Some(y)) = (arith_expr_const(&a), arith_expr_const(&b)) {
        let holds = match op {
            RelOp::Lt => x < y,
            RelOp::Le => x <= y,
            RelOp::Gt => x > y,
            RelOp::Ge => x >= y,
            RelOp::Eq => x == y,
            RelOp::Ne => x != y,
        };
        return if holds {
            PredicateReduction::Top
        } else {
            PredicateReduction::Bottom
        };
    }
    let rel = match op {
        RelOp::Gt => PredExpr::Rel(b, RelOp::Lt, a),
        RelOp::Ge => PredExpr::Rel(b, RelOp::Le, a),
        RelOp::Eq | RelOp::Ne => {
            let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
            PredExpr::Rel(lo, op, hi)
        }
        RelOp::Lt | RelOp::Le => PredExpr::Rel(a, op, b),
    };
    PredicateReduction::Predicate(rel)
}

fn neg_rel_op(op: RelOp) -> RelOp {
    match op {
        RelOp::Lt => RelOp::Ge,
        RelOp::Le => RelOp::Gt,
        RelOp::Gt => RelOp::Le,
        RelOp::Ge => RelOp::Lt,
        RelOp::Eq => RelOp::Ne,
        RelOp::Ne => RelOp::Eq,
    }
}

fn neg_mem_op(op: MemOp) -> MemOp {
    match op {
        MemOp::In => MemOp::NotIn,
        MemOp::NotIn => MemOp::In,
    }
}

/// `e` already canonical; `set` is sorted/deduped by type.
fn canon_mem(e: ArithExpr, op: MemOp, set: BTreeSet<i64>) -> PredicateReduction {
    if let Some(x) = arith_expr_const(&e) {
        let present = set.contains(&x);
        let holds = match op {
            MemOp::In => present,
            MemOp::NotIn => !present,
        };
        return if holds {
            PredicateReduction::Top
        } else {
            PredicateReduction::Bottom
        };
    }
    match set.len() {
        0 => match op {
            MemOp::In => PredicateReduction::Bottom,
            MemOp::NotIn => PredicateReduction::Top,
        },
        1 => {
            let only = *set.iter().next().unwrap();
            let rel_op = match op {
                MemOp::In => RelOp::Eq,
                MemOp::NotIn => RelOp::Ne,
            };
            canon_rel(e, rel_op, arith_expr_lit(only))
        }
        _ => PredicateReduction::Predicate(PredExpr::Mem(e, op, set)),
    }
}

fn negate(form: PredicateReduction) -> PredicateReduction {
    match form {
        PredicateReduction::Top => PredicateReduction::Bottom,
        PredicateReduction::Bottom => PredicateReduction::Top,
        PredicateReduction::Predicate(p) => negate_pred_expr(p),
    }
}

/// `predicate` is already canonical (so it carries no `Not`).
fn negate_pred_expr(predicate: PredExpr) -> PredicateReduction {
    match predicate {
        PredExpr::Rel(a, op, b) => canon_rel(a, neg_rel_op(op), b),
        PredExpr::Mem(e, op, set) => canon_mem(e, neg_mem_op(op), set),
        PredExpr::Not(inner) => reduce_pred_expr(*inner),
        PredExpr::And(operands) => {
            canon_junction_forms(operands.into_iter().map(negate_pred_expr).collect(), false)
        }
        PredExpr::Or(operands) => {
            canon_junction_forms(operands.into_iter().map(negate_pred_expr).collect(), true)
        }
    }
}

fn canon_junction(operands: Vec<PredExpr>, is_and: bool) -> PredicateReduction {
    canon_junction_forms(operands.into_iter().map(reduce_pred_expr).collect(), is_and)
}

fn canon_junction_forms(forms: Vec<PredicateReduction>, is_and: bool) -> PredicateReduction {
    let mut operands: Vec<PredExpr> = Vec::new();
    for form in forms {
        match form {
            PredicateReduction::Top => {
                if !is_and {
                    return PredicateReduction::Top;
                }
            }
            PredicateReduction::Bottom => {
                if is_and {
                    return PredicateReduction::Bottom;
                }
            }
            PredicateReduction::Predicate(p) => match p {
                PredExpr::And(inner) if is_and => operands.extend(inner),
                PredExpr::Or(inner) if !is_and => operands.extend(inner),
                other => operands.push(other),
            },
        }
    }
    operands.sort();
    operands.dedup();
    match operands.len() {
        0 => {
            if is_and {
                PredicateReduction::Top
            } else {
                PredicateReduction::Bottom
            }
        }
        1 => PredicateReduction::Predicate(operands.pop().unwrap()),
        _ => PredicateReduction::Predicate(if is_and {
            PredExpr::And(operands)
        } else {
            PredExpr::Or(operands)
        }),
    }
}

impl AsLit for NumForm {
    type Lit = i64;

    /// The single integer this value denotes, only when it is a literal.
    /// Non-canonicalizing: an `ArithExpr` or `LitSet` that would fold to a literal
    /// still returns `None` (canonicalize first if folding is wanted).
    #[inline]
    fn as_lit(&self) -> Option<i64> {
        match self {
            NumForm::Lit(n) => Some(*n),
            _ => None,
        }
    }
}

impl Lattice for NumForm {
    #[inline]
    fn is_undetermined(&self) -> bool {
        matches!(self, NumForm::Undetermined)
    }

    /// Bottom of the lattice — resolves to a single concrete integer. Aligned
    /// with `as_lit`: literal only, not canonicalizing.
    #[inline]
    fn is_ground(&self) -> bool {
        matches!(self, NumForm::Lit(_))
    }

    /// Greatest lower bound, canonicalizing both operands and the result.
    /// Distinct symbolic forms (`ArithExpr`/`PredExpr`) meet only when equal once
    /// canonical; symbolic versus concrete is rejected.
    fn meet(&self, other: &Self) -> Option<Self> {
        let a = self.canonical().ok()?;
        let b = other.canonical().ok()?;
        use NumForm::*;
        Some(match (a.as_ref(), b.as_ref()) {
            (Undetermined, _) => b.as_ref().clone(),
            (_, Undetermined) => a.as_ref().clone(),
            (Lit(x), Lit(y)) => {
                if x == y {
                    Lit(*x)
                } else {
                    return None;
                }
            }
            (Lit(x), LitSet(s)) | (LitSet(s), Lit(x)) => {
                if s.contains(x) {
                    Lit(*x)
                } else {
                    return None;
                }
            }
            (LitSet(s), LitSet(t)) => return lift_set(s.intersection(t).copied().collect()).ok(),
            (RangeFrom(i), Lit(n)) | (Lit(n), RangeFrom(i)) => {
                if n >= i {
                    Lit(*n)
                } else {
                    return None;
                }
            }
            (RangeTo(j), Lit(n)) | (Lit(n), RangeTo(j)) => {
                if n < j {
                    Lit(*n)
                } else {
                    return None;
                }
            }
            (RangeFrom(i), RangeFrom(j)) => RangeFrom((*i).max(*j)),
            (RangeTo(i), RangeTo(j)) => RangeTo((*i).min(*j)),
            (RangeFrom(i), RangeTo(j)) | (RangeTo(j), RangeFrom(i)) => {
                return lift_set((*i..*j).collect()).ok();
            }
            (RangeFrom(i), LitSet(s)) | (LitSet(s), RangeFrom(i)) => {
                return lift_set(s.iter().copied().filter(|&x| x >= *i).collect()).ok();
            }
            (RangeTo(j), LitSet(s)) | (LitSet(s), RangeTo(j)) => {
                return lift_set(s.iter().copied().filter(|&x| x < *j).collect()).ok();
            }
            (x, y) => {
                if x == y {
                    x.clone()
                } else {
                    return None;
                }
            }
        })
    }

    /// Least upper bound, canonicalizing both operands and the result.
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let a = self
            .canonical()
            .unwrap_or(Cow::Owned(NumForm::Undetermined));
        let b = other
            .canonical()
            .unwrap_or(Cow::Owned(NumForm::Undetermined));
        use NumForm::*;
        Ok(match (a.as_ref(), b.as_ref()) {
            (Undetermined, _) | (_, Undetermined) => Undetermined,
            (Lit(x), Lit(y)) => {
                if x == y {
                    Lit(*x)
                } else {
                    LitSet(Box::new([*x, *y].into_iter().collect()))
                }
            }
            (Lit(x), LitSet(s)) | (LitSet(s), Lit(x)) => {
                let mut union = s.as_ref().clone();
                union.insert(*x);
                litset_or_lit(union)
            }
            (LitSet(s), LitSet(t)) => litset_or_lit(s.union(t).copied().collect()),
            (RangeFrom(i), RangeFrom(j)) => RangeFrom((*i).min(*j)),
            (RangeTo(i), RangeTo(j)) => RangeTo((*i).max(*j)),
            (x, y) => {
                if x == y {
                    x.clone()
                } else {
                    Undetermined
                }
            }
        })
    }

    /// Partial-order check `target ⊑ self` (pattern admits every value target
    /// admits), specialized for the literal cases so the matcher's per-candidate
    /// path allocates nothing. `ArithExpr`/`PredExpr` on either side fall back to the
    /// canonicalizing `meet`-derived default, which this must equal.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::ArithExpr(_) | Self::PredExpr(_), _)
            | (_, Self::ArithExpr(_) | Self::PredExpr(_)) => {
                match (self.meet(target), target.canonical()) {
                    (Some(meet), Ok(target)) => meet == *target,
                    _ => false,
                }
            }
            (Self::Undetermined, Self::Undetermined | Self::Lit(_)) => true,
            (Self::Lit(_), Self::Undetermined) => false,
            (Self::Lit(p), Self::Lit(t)) => p == t,
            (Self::Undetermined, Self::LitSet(t)) => !t.is_empty(),
            (Self::LitSet(_), Self::Undetermined) => false,
            (Self::Lit(p), Self::LitSet(t)) => t.len() == 1 && t.contains(p),
            (Self::LitSet(p), Self::Lit(t)) => p.contains(t),
            (Self::LitSet(p), Self::LitSet(t)) => !t.is_empty() && t.iter().all(|x| p.contains(x)),
            (Self::RangeFrom(i), Self::Lit(t)) => t >= i,
            (Self::RangeFrom(i), Self::LitSet(t)) => !t.is_empty() && t.iter().all(|x| x >= i),
            (Self::RangeFrom(i), Self::RangeFrom(j)) => j >= i,
            (Self::RangeTo(j), Self::Lit(t)) => t < j,
            (Self::RangeTo(j), Self::LitSet(t)) => !t.is_empty() && t.iter().all(|x| x < j),
            (Self::RangeTo(j), Self::RangeTo(k)) => k <= j,
            (Self::Undetermined, Self::RangeFrom(_) | Self::RangeTo(_)) => true,
            (Self::RangeFrom(_) | Self::RangeTo(_), _)
            | (_, Self::RangeFrom(_) | Self::RangeTo(_)) => false,
        }
    }
}

/// A non-empty canonical set as a `NumForm`: a singleton collapses to `Lit`.
fn litset_or_lit(set: BTreeSet<i64>) -> NumForm {
    match set.len() {
        0 => NumForm::Undetermined,
        1 => NumForm::Lit(*set.iter().next().unwrap()),
        _ => NumForm::LitSet(Box::new(set)),
    }
}

// Arithmetic on `NumForm` propagates `Undetermined`: only `Lit op Lit` yields a
// `Lit`. Every binop has impls for all four `(owned|ref) × (owned|ref)`
// combinations delegating to the ref-ref form, plus a bare `i64` on either side.
macro_rules! impl_value_binop {
    ($Op:ident, $op:ident, $lit_op:tt) => {
        impl $Op<&NumForm> for &NumForm {
            type Output = NumForm;
            fn $op(self, rhs: &NumForm) -> NumForm {
                match (self, rhs) {
                    (NumForm::Lit(a), NumForm::Lit(b)) => NumForm::Lit(a $lit_op b),
                    _ => NumForm::Undetermined,
                }
            }
        }
        impl $Op<NumForm> for &NumForm {
            type Output = NumForm;
            fn $op(self, rhs: NumForm) -> NumForm { self.$op(&rhs) }
        }
        impl $Op<&NumForm> for NumForm {
            type Output = NumForm;
            fn $op(self, rhs: &NumForm) -> NumForm { (&self).$op(rhs) }
        }
        impl $Op<NumForm> for NumForm {
            type Output = NumForm;
            fn $op(self, rhs: NumForm) -> NumForm { (&self).$op(&rhs) }
        }
        impl $Op<i64> for &NumForm {
            type Output = NumForm;
            fn $op(self, rhs: i64) -> NumForm { self.$op(&NumForm::Lit(rhs)) }
        }
        impl $Op<i64> for NumForm {
            type Output = NumForm;
            fn $op(self, rhs: i64) -> NumForm { (&self).$op(&NumForm::Lit(rhs)) }
        }
        impl $Op<&NumForm> for i64 {
            type Output = NumForm;
            fn $op(self, rhs: &NumForm) -> NumForm { (&NumForm::Lit(self)).$op(rhs) }
        }
        impl $Op<NumForm> for i64 {
            type Output = NumForm;
            fn $op(self, rhs: NumForm) -> NumForm { (&NumForm::Lit(self)).$op(&rhs) }
        }
    };
}

impl_value_binop!(Add, add, +);
impl_value_binop!(Sub, sub, -);
impl_value_binop!(Mul, mul, *);
impl_value_binop!(Div, div, /);

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_set(NumForm::lit_set([2, 1, 2]), NumForm::LitSet(Box::new(BTreeSet::from([1, 2]))))]
    #[case::var(NumForm::var("x"), NumForm::ArithExpr(Box::new(ArithExpr::Var("x".to_string()))))]
    #[case::rel_predicate(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("r".to_string()), RelOp::Ge, ArithExpr::Lit(1))), NumForm::PredExpr(Box::new(PredExpr::Rel(ArithExpr::Var("r".to_string()), RelOp::Ge, ArithExpr::Lit(1)))))]
    #[case::term(NumForm::arith_expr(ArithExpr::Var("x".to_string())), NumForm::ArithExpr(Box::new(ArithExpr::Var("x".to_string()))))]
    #[case::predicate(NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2]))), NumForm::PredExpr(Box::new(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])))))]
    fn test_num_form_constructors(#[case] actual: NumForm, #[case] expected: NumForm) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::i64(NumForm::from(5_i64), NumForm::Lit(5))]
    #[case::spin_multiplicity(NumForm::from(SpinMultiplicity::TRIPLET), NumForm::Lit(3))]
    #[case::vec(NumForm::from(vec![2, 1, 2]), NumForm::LitSet(Box::new(BTreeSet::from([1, 2]))))]
    #[case::btreeset(NumForm::from(BTreeSet::from([3, 1])), NumForm::LitSet(Box::new(BTreeSet::from([1, 3]))))]
    fn test_num_form_from(#[case] actual: NumForm, #[case] expected: NumForm) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::litset_singleton(NumForm::lit_set([7]), Ok(NumForm::Lit(7)))]
    #[case::litset_empty(NumForm::lit_set([]), Err(Contradiction))]
    #[case::term_ground(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(2), ArithExpr::Lit(3)])), Ok(NumForm::Lit(5)))]
    #[case::term_neg_lit(NumForm::arith_expr(ArithExpr::Neg(Box::new(ArithExpr::Lit(4)))), Ok(NumForm::Lit(-4)))]
    #[case::term_neg_neg(NumForm::arith_expr(ArithExpr::Neg(Box::new(ArithExpr::Neg(Box::new(ArithExpr::Var("x".to_string())))))), Ok(NumForm::arith_expr(ArithExpr::Var("x".to_string()))))]
    #[case::term_sum_identity(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Var("x".to_string()), ArithExpr::Lit(0)])), Ok(NumForm::arith_expr(ArithExpr::Var("x".to_string()))))]
    #[case::term_sum_sorted_const_first(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Var("x".to_string()), ArithExpr::Lit(1)])), Ok(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(1), ArithExpr::Var("x".to_string())]))))]
    #[case::term_sum_flatten(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Sum(vec![ArithExpr::Var("x".to_string()), ArithExpr::Lit(1)]), ArithExpr::Lit(2)])), Ok(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(3), ArithExpr::Var("x".to_string())]))))]
    #[case::term_product_annihilator(NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Var("x".to_string()), ArithExpr::Lit(0)])), Ok(NumForm::Lit(0)))]
    #[case::term_product_identity(NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Var("x".to_string()), ArithExpr::Lit(1)])), Ok(NumForm::arith_expr(ArithExpr::Var("x".to_string()))))]
    #[case::term_div_fold(NumForm::arith_expr(ArithExpr::Div(Box::new(ArithExpr::Lit(10)), Box::new(ArithExpr::Lit(3)))), Ok(NumForm::Lit(3)))]
    #[case::term_rem_fold(NumForm::arith_expr(ArithExpr::Rem(Box::new(ArithExpr::Lit(10)), Box::new(ArithExpr::Lit(3)))), Ok(NumForm::Lit(1)))]
    #[case::pred_rel_true(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Lit(1), RelOp::Eq, ArithExpr::Lit(1))), Ok(NumForm::Undetermined))]
    #[case::pred_rel_false(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Lit(1), RelOp::Eq, ArithExpr::Lit(2))), Err(Contradiction))]
    #[case::pred_rel_orient_ge(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("x".to_string()), RelOp::Ge, ArithExpr::Lit(1))), Ok(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Lit(1), RelOp::Le, ArithExpr::Var("x".to_string())))))]
    #[case::pred_rel_eq_sorted(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("x".to_string()), RelOp::Eq, ArithExpr::Lit(0))), Ok(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Lit(0), RelOp::Eq, ArithExpr::Var("x".to_string())))))]
    #[case::pred_not_eq_to_ne(NumForm::pred_expr(PredExpr::Not(Box::new(PredExpr::Rel(ArithExpr::Var("x".to_string()), RelOp::Eq, ArithExpr::Lit(0))))), Ok(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Lit(0), RelOp::Ne, ArithExpr::Var("x".to_string())))))]
    #[case::pred_mem_singleton(NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([5]))), Ok(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Lit(5), RelOp::Eq, ArithExpr::Var("x".to_string())))))]
    #[case::pred_mem_notin_empty(NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::NotIn, BTreeSet::new())), Ok(NumForm::Undetermined))]
    #[case::pred_mem_in_empty(NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::new())), Err(Contradiction))]
    #[case::pred_not_mem_to_notin(NumForm::pred_expr(PredExpr::Not(Box::new(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2]))))), Ok(NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::NotIn, BTreeSet::from([1, 2])))))]
    #[case::pred_and_drops_top(NumForm::pred_expr(PredExpr::And(vec![PredExpr::Rel(ArithExpr::Var("x".to_string()), RelOp::Le, ArithExpr::Lit(3)), PredExpr::Rel(ArithExpr::Lit(1), RelOp::Eq, ArithExpr::Lit(1))])), Ok(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("x".to_string()), RelOp::Le, ArithExpr::Lit(3)))))]
    #[case::pred_demorgan(NumForm::pred_expr(PredExpr::Not(Box::new(PredExpr::And(vec![PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), PredExpr::Mem(ArithExpr::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4]))])))), Ok(NumForm::pred_expr(PredExpr::Or(vec![PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::NotIn, BTreeSet::from([1, 2])), PredExpr::Mem(ArithExpr::Var("y".to_string()), MemOp::NotIn, BTreeSet::from([3, 4]))]))))]
    #[case::term_sum_neg_const(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Var("x".to_string()), ArithExpr::Lit(-3)])), Ok(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Var("x".to_string()), ArithExpr::Neg(Box::new(ArithExpr::Lit(3)))]))))]
    #[case::term_neg_zero(NumForm::arith_expr(ArithExpr::Neg(Box::new(ArithExpr::Lit(0)))), Ok(NumForm::Lit(0)))]
    #[case::term_product_flatten(NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Product(vec![ArithExpr::Var("x".to_string()), ArithExpr::Var("y".to_string())]), ArithExpr::Var("z".to_string())])), Ok(NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Var("x".to_string()), ArithExpr::Var("y".to_string()), ArithExpr::Var("z".to_string())]))))]
    #[case::term_product_sort(NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Var("b".to_string()), ArithExpr::Var("a".to_string())])), Ok(NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Var("a".to_string()), ArithExpr::Var("b".to_string())]))))]
    #[case::term_product_const_fold(NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Lit(2), ArithExpr::Lit(3), ArithExpr::Var("x".to_string())])), Ok(NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Lit(6), ArithExpr::Var("x".to_string())]))))]
    #[case::term_sum_empty(NumForm::arith_expr(ArithExpr::Sum(vec![])), Ok(NumForm::Lit(0)))]
    #[case::term_product_empty(NumForm::arith_expr(ArithExpr::Product(vec![])), Ok(NumForm::Lit(1)))]
    #[case::term_div_by_zero(NumForm::arith_expr(ArithExpr::Div(Box::new(ArithExpr::Lit(10)), Box::new(ArithExpr::Lit(0)))), Ok(NumForm::arith_expr(ArithExpr::Div(Box::new(ArithExpr::Lit(10)), Box::new(ArithExpr::Lit(0))))))]
    #[case::pred_and_flatten(NumForm::pred_expr(PredExpr::And(vec![PredExpr::And(vec![PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), PredExpr::Mem(ArithExpr::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4]))]), PredExpr::Mem(ArithExpr::Var("z".to_string()), MemOp::In, BTreeSet::from([5, 6]))])), Ok(NumForm::pred_expr(PredExpr::And(vec![PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), PredExpr::Mem(ArithExpr::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4])), PredExpr::Mem(ArithExpr::Var("z".to_string()), MemOp::In, BTreeSet::from([5, 6]))]))))]
    #[case::pred_and_sort_dedup(NumForm::pred_expr(PredExpr::And(vec![PredExpr::Mem(ArithExpr::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4])), PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2]))])), Ok(NumForm::pred_expr(PredExpr::And(vec![PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), PredExpr::Mem(ArithExpr::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4]))]))))]
    #[case::pred_and_bottom(NumForm::pred_expr(PredExpr::And(vec![PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), PredExpr::Rel(ArithExpr::Lit(1), RelOp::Eq, ArithExpr::Lit(2))])), Err(Contradiction))]
    #[case::pred_or_drops_bottom(NumForm::pred_expr(PredExpr::Or(vec![PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), PredExpr::Rel(ArithExpr::Lit(1), RelOp::Eq, ArithExpr::Lit(2))])), Ok(NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])))))]
    #[case::pred_or_top(NumForm::pred_expr(PredExpr::Or(vec![PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), PredExpr::Rel(ArithExpr::Lit(1), RelOp::Eq, ArithExpr::Lit(1))])), Ok(NumForm::Undetermined))]
    #[case::pred_and_empty(NumForm::pred_expr(PredExpr::And(vec![])), Ok(NumForm::Undetermined))]
    #[case::pred_or_empty(NumForm::pred_expr(PredExpr::Or(vec![])), Err(Contradiction))]
    #[case::pred_not_not(NumForm::pred_expr(PredExpr::Not(Box::new(PredExpr::Not(Box::new(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2]))))))), Ok(NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])))))]
    #[case::pred_not_le(NumForm::pred_expr(PredExpr::Not(Box::new(PredExpr::Rel(ArithExpr::Var("x".to_string()), RelOp::Le, ArithExpr::Lit(3))))), Ok(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Lit(3), RelOp::Lt, ArithExpr::Var("x".to_string())))))]
    #[case::pred_rel_orient_gt(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("x".to_string()), RelOp::Gt, ArithExpr::Lit(1))), Ok(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Lit(1), RelOp::Lt, ArithExpr::Var("x".to_string())))))]
    #[case::pred_mem_notin_singleton(NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::NotIn, BTreeSet::from([5]))), Ok(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Lit(5), RelOp::Ne, ArithExpr::Var("x".to_string())))))]
    fn test_num_form_canonicalize(
        #[case] input: NumForm,
        #[case] expected: Result<NumForm, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(NumForm::Undetermined)]
    #[case::lit(NumForm::Lit(3))]
    #[case::litset(NumForm::lit_set([1, 2, 3]))]
    #[case::term_var(NumForm::arith_expr(ArithExpr::Var("x".to_string())))]
    fn test_num_form_canonicalize_identity(#[case] input: NumForm) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::sum(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Var("b".to_string()), ArithExpr::Lit(2), ArithExpr::Var("a".to_string()), ArithExpr::Lit(3)])))]
    #[case::product(NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Var("b".to_string()), ArithExpr::Var("a".to_string())])))]
    #[case::rel(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("x".to_string()), RelOp::Ge, ArithExpr::Lit(1))))]
    #[case::or(NumForm::pred_expr(PredExpr::Or(vec![PredExpr::Mem(ArithExpr::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4])), PredExpr::Mem(ArithExpr::Var("x".to_string()), MemOp::NotIn, BTreeSet::from([1, 2]))])))]
    fn test_num_form_canonicalize_idempotent(#[case] input: NumForm) {
        let once = input.canonicalize().unwrap();
        let twice = once.clone().canonicalize().unwrap();
        assert_eq!(once, twice);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(NumForm::Lit(3), Some(3))]
    #[case::lit_neg(NumForm::Lit(-5), Some(-5))]
    #[case::undetermined(NumForm::Undetermined, None)]
    #[case::litset(NumForm::lit_set([1, 2]), None)]
    #[case::term(NumForm::arith_expr(ArithExpr::Var("x".to_string())), None)]
    fn test_num_form_as_lit(#[case] ast: NumForm, #[case] expected: Option<i64>) {
        assert_eq!(ast.as_lit(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(NumForm::Undetermined, true)]
    #[case::lit(NumForm::Lit(3), false)]
    #[case::litset(NumForm::lit_set([1, 2]), false)]
    #[case::term(NumForm::var("x"), false)]
    #[case::predicate(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("r".to_string()), RelOp::Ge, ArithExpr::Lit(1))), false)]
    #[case::range_from(NumForm::RangeFrom(1), false)]
    fn test_num_form_is_undetermined(#[case] ast: NumForm, #[case] expected: bool) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(NumForm::Undetermined, NumForm::Lit(3), Some(NumForm::Lit(3)))]
    #[case::lit_und(NumForm::Lit(3), NumForm::Undetermined, Some(NumForm::Lit(3)))]
    #[case::lit_lit_eq(NumForm::Lit(3), NumForm::Lit(3), Some(NumForm::Lit(3)))]
    #[case::lit_lit_neq(NumForm::Lit(3), NumForm::Lit(4), None)]
    #[case::lit_set_in(NumForm::Lit(2), NumForm::lit_set([1, 2, 3]), Some(NumForm::Lit(2)))]
    #[case::lit_set_out(NumForm::Lit(5), NumForm::lit_set([1, 2, 3]), None)]
    #[case::set_set_multi(NumForm::lit_set([1, 2, 3]), NumForm::lit_set([2, 3, 4]), Some(NumForm::lit_set([2, 3])))]
    #[case::set_set_singleton(NumForm::lit_set([1, 2]), NumForm::lit_set([2, 3]), Some(NumForm::Lit(2)))]
    #[case::set_set_empty(NumForm::lit_set([1, 2]), NumForm::lit_set([3, 4]), None)]
    #[case::term_term_eq(NumForm::arith_expr(ArithExpr::Var("x".to_string())), NumForm::arith_expr(ArithExpr::Var("x".to_string())), Some(NumForm::arith_expr(ArithExpr::Var("x".to_string()))))]
    #[case::term_term_neq(NumForm::arith_expr(ArithExpr::Var("x".to_string())), NumForm::arith_expr(ArithExpr::Var("y".to_string())), None)]
    #[case::arith_expr_lit(NumForm::arith_expr(ArithExpr::Var("x".to_string())), NumForm::Lit(5), None)]
    #[case::rangefrom_lit_in(NumForm::RangeFrom(1), NumForm::Lit(2), Some(NumForm::Lit(2)))]
    #[case::rangefrom_lit_out(NumForm::RangeFrom(2), NumForm::Lit(1), None)]
    #[case::rangefrom_rangefrom(NumForm::RangeFrom(1), NumForm::RangeFrom(3), Some(NumForm::RangeFrom(3)))]
    #[case::rangeto_rangeto(NumForm::RangeTo(5), NumForm::RangeTo(3), Some(NumForm::RangeTo(3)))]
    #[case::rangefrom_rangeto_set(NumForm::RangeFrom(1), NumForm::RangeTo(4), Some(NumForm::lit_set([1, 2, 3])))]
    #[case::rangefrom_rangeto_empty(NumForm::RangeFrom(4), NumForm::RangeTo(2), None)]
    #[case::rangefrom_set(NumForm::RangeFrom(2), NumForm::lit_set([1, 2, 3]), Some(NumForm::lit_set([2, 3])))]
    fn test_num_form_meet(
        #[case] a: NumForm,
        #[case] b: NumForm,
        #[case] expected: Option<NumForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(NumForm::Undetermined, NumForm::Lit(3), NumForm::Undetermined)]
    #[case::lit_lit_eq(NumForm::Lit(3), NumForm::Lit(3), NumForm::Lit(3))]
    #[case::lit_lit_neq(NumForm::Lit(3), NumForm::Lit(4), NumForm::lit_set([3, 4]))]
    #[case::lit_set(NumForm::Lit(5), NumForm::lit_set([1, 2, 3]), NumForm::lit_set([1, 2, 3, 5]))]
    #[case::set_set(NumForm::lit_set([1, 2]), NumForm::lit_set([2, 3]), NumForm::lit_set([1, 2, 3]))]
    #[case::term_term_eq(NumForm::arith_expr(ArithExpr::Var("x".to_string())), NumForm::arith_expr(ArithExpr::Var("x".to_string())), NumForm::arith_expr(ArithExpr::Var("x".to_string())))]
    #[case::term_term_neq(NumForm::arith_expr(ArithExpr::Var("x".to_string())), NumForm::arith_expr(ArithExpr::Var("y".to_string())), NumForm::Undetermined)]
    #[case::rangefrom_rangefrom(NumForm::RangeFrom(3), NumForm::RangeFrom(1), NumForm::RangeFrom(1))]
    #[case::rangeto_rangeto(NumForm::RangeTo(3), NumForm::RangeTo(5), NumForm::RangeTo(5))]
    #[case::rangefrom_lit_overapprox(NumForm::RangeFrom(1), NumForm::Lit(5), NumForm::Undetermined)]
    fn test_num_form_join(#[case] a: NumForm, #[case] b: NumForm, #[case] expected: NumForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_und(NumForm::Undetermined, NumForm::Undetermined, true)]
    #[case::und_lit(NumForm::Undetermined, NumForm::Lit(3), true)]
    #[case::lit_und(NumForm::Lit(3), NumForm::Undetermined, false)]
    #[case::lit_lit(NumForm::Lit(3), NumForm::Lit(3), true)]
    #[case::lit_lit_neq(NumForm::Lit(3), NumForm::Lit(4), false)]
    #[case::set_lit_in(NumForm::lit_set([1, 2, 3]), NumForm::Lit(2), true)]
    #[case::set_lit_out(NumForm::lit_set([1, 2, 3]), NumForm::Lit(5), false)]
    #[case::set_set(NumForm::lit_set([1, 2, 3]), NumForm::lit_set([1, 2]), true)]
    #[case::rangefrom_lit_ge(NumForm::RangeFrom(1), NumForm::Lit(2), true)]
    #[case::rangefrom_lit_lt(NumForm::RangeFrom(2), NumForm::Lit(1), false)]
    #[case::rangefrom_rangefrom_wider(NumForm::RangeFrom(1), NumForm::RangeFrom(2), true)]
    #[case::rangefrom_rangefrom_narrower(NumForm::RangeFrom(2), NumForm::RangeFrom(1), false)]
    #[case::rangefrom_und(NumForm::RangeFrom(1), NumForm::Undetermined, false)]
    #[case::und_rangefrom(NumForm::Undetermined, NumForm::RangeFrom(1), true)]
    #[case::rangefrom_set_all_ge(NumForm::RangeFrom(1), NumForm::lit_set([2, 3]), true)]
    #[case::rangefrom_set_some_lt(NumForm::RangeFrom(2), NumForm::lit_set([1, 3]), false)]
    #[case::rangeto_lit_lt(NumForm::RangeTo(3), NumForm::Lit(2), true)]
    #[case::rangeto_lit_ge(NumForm::RangeTo(2), NumForm::Lit(3), false)]
    fn test_num_form_matches(
        #[case] pattern: NumForm,
        #[case] target: NumForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::no_change(NumForm::Lit(3), NumForm::Lit(3), false, NumForm::Lit(3))]
    #[case::tighten(NumForm::Undetermined, NumForm::Lit(3), true, NumForm::Lit(3))]
    #[case::incompatible(NumForm::Lit(3), NumForm::Lit(4), false, NumForm::Lit(3))]
    fn test_num_form_narrow_from(
        #[case] mut target: NumForm,
        #[case] source: NumForm,
        #[case] expected_changed: bool,
        #[case] expected_after: NumForm,
    ) {
        let changed = target.narrow_from(&source);
        assert_eq!(changed, expected_changed);
        assert_eq!(target, expected_after);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::no_change(NumForm::Lit(3), NumForm::Lit(3), false, NumForm::Lit(3))]
    #[case::widen_to_set(NumForm::Lit(3), NumForm::Lit(4), true, NumForm::lit_set([3, 4]))]
    #[case::widen_to_top(NumForm::Lit(3), NumForm::Undetermined, true, NumForm::Undetermined)]
    fn test_num_form_widen_with(
        #[case] mut target: NumForm,
        #[case] source: NumForm,
        #[case] expected_changed: bool,
        #[case] expected_after: NumForm,
    ) {
        let changed = target.widen_with(&source);
        assert_eq!(changed, Ok(expected_changed));
        assert_eq!(target, expected_after);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(NumForm::Lit(2), NumForm::Lit(3), NumForm::Lit(5))]
    #[case::lit_undetermined(NumForm::Lit(2), NumForm::Undetermined, NumForm::Undetermined)]
    fn test_num_form_add(#[case] lhs: NumForm, #[case] rhs: NumForm, #[case] expected: NumForm) {
        assert_eq!(lhs + rhs, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(NumForm::Lit(5), NumForm::Lit(3), NumForm::Lit(2))]
    #[case::lit_undetermined(NumForm::Lit(5), NumForm::Undetermined, NumForm::Undetermined)]
    fn test_num_form_sub(#[case] lhs: NumForm, #[case] rhs: NumForm, #[case] expected: NumForm) {
        assert_eq!(lhs - rhs, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(NumForm::Lit(4), NumForm::Lit(3), NumForm::Lit(12))]
    #[case::lit_undetermined(NumForm::Lit(4), NumForm::Undetermined, NumForm::Undetermined)]
    fn test_num_form_mul(#[case] lhs: NumForm, #[case] rhs: NumForm, #[case] expected: NumForm) {
        assert_eq!(lhs * rhs, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(NumForm::Lit(10), NumForm::Lit(3), NumForm::Lit(3))]
    #[case::lit_undetermined(NumForm::Lit(10), NumForm::Undetermined, NumForm::Undetermined)]
    fn test_num_form_div(#[case] lhs: NumForm, #[case] rhs: NumForm, #[case] expected: NumForm) {
        assert_eq!(lhs / rhs, expected);
    }

    #[rstest]
    #[should_panic]
    fn test_num_form_div_error() {
        let _ = NumForm::Lit(5) / NumForm::Lit(0);
    }
}
