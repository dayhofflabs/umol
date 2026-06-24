//! Value AST.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::ops::{Add, Div, Mul, Sub};

use umol_chem::spin::SpinMultiplicity;

use super::error::Contradiction;
use super::operators::{MemOp, RelOp};
use super::traits::{AsLit, Canonicalize, Lattice};

/// Integer-valued atom/bond field: undetermined (pattern wildcard), a literal,
/// a finite literal set, an arithmetic term over variables, or a boolean
/// predicate constraining the field. Used for charge, hydrogen count, isotope
/// mass, valence, bond order, etc.
///
/// Equality is **lazy**: derived `Eq`/`Hash`/`Ord` are structural ("same
/// tree"); semantic equality is `Canonicalize::equiv`, which compares canonical
/// forms.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ValueAst {
    #[default]
    Undetermined,
    Lit(i64),
    LitSet(Box<BTreeSet<i64>>),
    RangeFrom(i64),
    RangeTo(i64),
    Term(Box<ValueTerm>),
    Predicate(Box<ValuePredicate>),
}

/// Arithmetic term over `i64`, the value-sort half of the field grammar. `Sum`
/// and `Product` are n-ary (associative + commutative by construction);
/// subtraction lowers to `Sum([a, Neg(b)])`; `Div`/`Rem` stay binary. A ground
/// term folds to a `Lit` (or `Neg(Lit)`), which `ValueAst` lifts to `Lit`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ValueTerm {
    Lit(i64),
    Var(String),
    Neg(Box<ValueTerm>),
    Sum(Vec<ValueTerm>),
    Product(Vec<ValueTerm>),
    Div(Box<ValueTerm>, Box<ValueTerm>),
    Rem(Box<ValueTerm>, Box<ValueTerm>),
}

/// Boolean predicate over terms, the constraint-sort half of the field grammar.
/// `Rel` and `Mem` operators are negation-closed (`RelOp` has `Ne`, `MemOp` has
/// `NotIn`), so canonicalization eliminates `Not` entirely — it survives only as
/// faithful parser input for `!`. ⊤/⊥ are not variants: a predicate that decides
/// is lifted by `ValueAst` to `Undetermined` / `Err(Contradiction)`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ValuePredicate {
    Rel(ValueTerm, RelOp, ValueTerm),
    Mem(ValueTerm, MemOp, BTreeSet<i64>),
    Not(Box<ValuePredicate>),
    And(Vec<ValuePredicate>),
    Or(Vec<ValuePredicate>),
}

impl ValueAst {
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
        Self::Term(Box::new(ValueTerm::Var(name.into())))
    }

    pub fn range_from(bound: i64) -> Self {
        Self::RangeFrom(bound)
    }

    pub fn range_to(bound: i64) -> Self {
        Self::RangeTo(bound)
    }

    pub fn term(term: ValueTerm) -> Self {
        Self::Term(Box::new(term))
    }

    pub fn predicate(predicate: ValuePredicate) -> Self {
        Self::Predicate(Box::new(predicate))
    }
}

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
        Self::lit_set(values)
    }
}

impl From<BTreeSet<i64>> for ValueAst {
    fn from(values: BTreeSet<i64>) -> Self {
        Self::LitSet(Box::new(values))
    }
}

impl Canonicalize for ValueAst {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            ValueAst::Undetermined => ValueAst::Undetermined,
            ValueAst::Lit(n) => ValueAst::Lit(n),
            ValueAst::LitSet(set) => lift_set(*set)?,
            ValueAst::RangeFrom(n) => ValueAst::RangeFrom(n),
            ValueAst::RangeTo(n) => ValueAst::RangeTo(n),
            ValueAst::Term(term) => lift_term(canon_term(*term)),
            ValueAst::Predicate(predicate) => match canon_predicate(*predicate) {
                PredicateForm::Top => ValueAst::Undetermined,
                PredicateForm::Bottom => return Err(Contradiction),
                PredicateForm::Predicate(p) => ValueAst::Predicate(Box::new(p)),
            },
        })
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            ValueAst::Undetermined | ValueAst::Lit(_) => Ok(Cow::Borrowed(self)),
            _ => Ok(Cow::Owned(self.clone().canonicalize()?)),
        }
    }
}

impl Canonicalize for ValueTerm {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(canon_term(self))
    }
}

/// Lift a canonical set: empty is unsatisfiable, a singleton is a `Lit`.
fn lift_set(set: BTreeSet<i64>) -> Result<ValueAst, Contradiction> {
    match set.len() {
        0 => Err(Contradiction),
        1 => Ok(ValueAst::Lit(*set.iter().next().unwrap())),
        _ => Ok(ValueAst::LitSet(Box::new(set))),
    }
}

/// Lift a canonical term: a ground term (`Lit` or `Neg(Lit)`) becomes `Lit`.
fn lift_term(term: ValueTerm) -> ValueAst {
    match term_const(&term) {
        Some(n) => ValueAst::Lit(n),
        None => ValueAst::Term(Box::new(term)),
    }
}

/// Canonical literal term: the sign lives in `Neg`, so `Lit` is always ≥ 0.
fn term_lit(n: i64) -> ValueTerm {
    if n < 0 {
        ValueTerm::Neg(Box::new(ValueTerm::Lit(-n)))
    } else {
        ValueTerm::Lit(n)
    }
}

/// The integer a canonical term denotes, if ground (`Lit(n)` or `Neg(Lit(n))`).
fn term_const(term: &ValueTerm) -> Option<i64> {
    match term {
        ValueTerm::Lit(n) => Some(*n),
        ValueTerm::Neg(inner) => match inner.as_ref() {
            ValueTerm::Lit(n) => Some(-n),
            _ => None,
        },
        _ => None,
    }
}

fn canon_term(term: ValueTerm) -> ValueTerm {
    match term {
        ValueTerm::Lit(n) => term_lit(n),
        ValueTerm::Var(_) => term,
        ValueTerm::Neg(inner) => canon_neg(canon_term(*inner)),
        ValueTerm::Sum(operands) => canon_sum(operands),
        ValueTerm::Product(operands) => canon_product(operands),
        ValueTerm::Div(a, b) => canon_div_rem(canon_term(*a), canon_term(*b), false),
        ValueTerm::Rem(a, b) => canon_div_rem(canon_term(*a), canon_term(*b), true),
    }
}

/// `inner` is already canonical.
fn canon_neg(inner: ValueTerm) -> ValueTerm {
    match inner {
        ValueTerm::Neg(grand) => *grand,
        ValueTerm::Lit(0) => ValueTerm::Lit(0),
        other => ValueTerm::Neg(Box::new(other)),
    }
}

fn canon_sum(operands: Vec<ValueTerm>) -> ValueTerm {
    let mut terms = Vec::new();
    let mut constant: i64 = 0;
    flatten_sum(operands, &mut terms, &mut constant);
    if constant != 0 {
        terms.push(term_lit(constant));
    }
    terms.sort();
    match terms.len() {
        0 => ValueTerm::Lit(0),
        1 => terms.pop().unwrap(),
        _ => ValueTerm::Sum(terms),
    }
}

fn flatten_sum(operands: Vec<ValueTerm>, terms: &mut Vec<ValueTerm>, constant: &mut i64) {
    for operand in operands {
        match canon_term(operand) {
            ValueTerm::Sum(inner) => flatten_sum(inner, terms, constant),
            other => match term_const(&other) {
                Some(n) => *constant += n,
                None => terms.push(other),
            },
        }
    }
}

fn canon_product(operands: Vec<ValueTerm>) -> ValueTerm {
    let mut terms = Vec::new();
    let mut constant: i64 = 1;
    flatten_product(operands, &mut terms, &mut constant);
    if constant == 0 {
        return ValueTerm::Lit(0);
    }
    if constant != 1 {
        terms.push(term_lit(constant));
    }
    terms.sort();
    match terms.len() {
        0 => ValueTerm::Lit(1),
        1 => terms.pop().unwrap(),
        _ => ValueTerm::Product(terms),
    }
}

fn flatten_product(operands: Vec<ValueTerm>, terms: &mut Vec<ValueTerm>, constant: &mut i64) {
    for operand in operands {
        match canon_term(operand) {
            ValueTerm::Product(inner) => flatten_product(inner, terms, constant),
            other => match term_const(&other) {
                Some(n) => *constant *= n,
                None => terms.push(other),
            },
        }
    }
}

/// `a`/`b` already canonical; folds a ground `(Lit) op (Lit)` when divisor ≠ 0.
fn canon_div_rem(a: ValueTerm, b: ValueTerm, is_rem: bool) -> ValueTerm {
    if let (Some(x), Some(y)) = (term_const(&a), term_const(&b)) {
        if y != 0 {
            return term_lit(if is_rem { x % y } else { x / y });
        }
    }
    let (a, b) = (Box::new(a), Box::new(b));
    if is_rem {
        ValueTerm::Rem(a, b)
    } else {
        ValueTerm::Div(a, b)
    }
}

/// Predicate canonical form, threading ⊤/⊥ that no `ValuePredicate` variant can
/// hold. Lifted to `ValueAst` (⊤ → `Undetermined`, ⊥ → `Err`).
enum PredicateForm {
    Top,
    Bottom,
    Predicate(ValuePredicate),
}

fn canon_predicate(predicate: ValuePredicate) -> PredicateForm {
    match predicate {
        ValuePredicate::Rel(a, op, b) => canon_rel(canon_term(a), op, canon_term(b)),
        ValuePredicate::Mem(e, op, set) => canon_mem(canon_term(e), op, set),
        ValuePredicate::Not(inner) => negate(canon_predicate(*inner)),
        ValuePredicate::And(operands) => canon_junction(operands, true),
        ValuePredicate::Or(operands) => canon_junction(operands, false),
    }
}

/// `a`/`b` already canonical.
fn canon_rel(a: ValueTerm, op: RelOp, b: ValueTerm) -> PredicateForm {
    if let (Some(x), Some(y)) = (term_const(&a), term_const(&b)) {
        let holds = match op {
            RelOp::Lt => x < y,
            RelOp::Le => x <= y,
            RelOp::Gt => x > y,
            RelOp::Ge => x >= y,
            RelOp::Eq => x == y,
            RelOp::Ne => x != y,
        };
        return if holds {
            PredicateForm::Top
        } else {
            PredicateForm::Bottom
        };
    }
    let rel = match op {
        RelOp::Gt => ValuePredicate::Rel(b, RelOp::Lt, a),
        RelOp::Ge => ValuePredicate::Rel(b, RelOp::Le, a),
        RelOp::Eq | RelOp::Ne => {
            let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
            ValuePredicate::Rel(lo, op, hi)
        }
        RelOp::Lt | RelOp::Le => ValuePredicate::Rel(a, op, b),
    };
    PredicateForm::Predicate(rel)
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
fn canon_mem(e: ValueTerm, op: MemOp, set: BTreeSet<i64>) -> PredicateForm {
    if let Some(x) = term_const(&e) {
        let present = set.contains(&x);
        let holds = match op {
            MemOp::In => present,
            MemOp::NotIn => !present,
        };
        return if holds {
            PredicateForm::Top
        } else {
            PredicateForm::Bottom
        };
    }
    match set.len() {
        0 => match op {
            MemOp::In => PredicateForm::Bottom,
            MemOp::NotIn => PredicateForm::Top,
        },
        1 => {
            let only = *set.iter().next().unwrap();
            let rel_op = match op {
                MemOp::In => RelOp::Eq,
                MemOp::NotIn => RelOp::Ne,
            };
            canon_rel(e, rel_op, term_lit(only))
        }
        _ => PredicateForm::Predicate(ValuePredicate::Mem(e, op, set)),
    }
}

fn negate(form: PredicateForm) -> PredicateForm {
    match form {
        PredicateForm::Top => PredicateForm::Bottom,
        PredicateForm::Bottom => PredicateForm::Top,
        PredicateForm::Predicate(p) => negate_predicate(p),
    }
}

/// `predicate` is already canonical (so it carries no `Not`).
fn negate_predicate(predicate: ValuePredicate) -> PredicateForm {
    match predicate {
        ValuePredicate::Rel(a, op, b) => canon_rel(a, neg_rel_op(op), b),
        ValuePredicate::Mem(e, op, set) => canon_mem(e, neg_mem_op(op), set),
        ValuePredicate::Not(inner) => canon_predicate(*inner),
        ValuePredicate::And(operands) => {
            canon_junction_forms(operands.into_iter().map(negate_predicate).collect(), false)
        }
        ValuePredicate::Or(operands) => {
            canon_junction_forms(operands.into_iter().map(negate_predicate).collect(), true)
        }
    }
}

fn canon_junction(operands: Vec<ValuePredicate>, is_and: bool) -> PredicateForm {
    canon_junction_forms(operands.into_iter().map(canon_predicate).collect(), is_and)
}

fn canon_junction_forms(forms: Vec<PredicateForm>, is_and: bool) -> PredicateForm {
    let mut operands: Vec<ValuePredicate> = Vec::new();
    for form in forms {
        match form {
            PredicateForm::Top => {
                if !is_and {
                    return PredicateForm::Top;
                }
            }
            PredicateForm::Bottom => {
                if is_and {
                    return PredicateForm::Bottom;
                }
            }
            PredicateForm::Predicate(p) => match p {
                ValuePredicate::And(inner) if is_and => operands.extend(inner),
                ValuePredicate::Or(inner) if !is_and => operands.extend(inner),
                other => operands.push(other),
            },
        }
    }
    operands.sort();
    operands.dedup();
    match operands.len() {
        0 => {
            if is_and {
                PredicateForm::Top
            } else {
                PredicateForm::Bottom
            }
        }
        1 => PredicateForm::Predicate(operands.pop().unwrap()),
        _ => PredicateForm::Predicate(if is_and {
            ValuePredicate::And(operands)
        } else {
            ValuePredicate::Or(operands)
        }),
    }
}

impl AsLit for ValueAst {
    type Lit = i64;

    /// The single integer this value denotes, only when it is a literal.
    /// Non-canonicalizing: a `Term` or `LitSet` that would fold to a literal
    /// still returns `None` (canonicalize first if folding is wanted).
    #[inline]
    fn as_lit(&self) -> Option<i64> {
        match self {
            ValueAst::Lit(n) => Some(*n),
            _ => None,
        }
    }
}

impl Lattice for ValueAst {
    #[inline]
    fn is_undetermined(&self) -> bool {
        matches!(self, ValueAst::Undetermined)
    }

    /// Bottom of the lattice — resolves to a single concrete integer. Aligned
    /// with `as_lit`: literal only, not canonicalizing.
    #[inline]
    fn is_ground(&self) -> bool {
        matches!(self, ValueAst::Lit(_))
    }

    /// Greatest lower bound, canonicalizing both operands and the result.
    /// Distinct symbolic forms (`Term`/`Predicate`) meet only when equal once
    /// canonical; symbolic versus concrete is rejected.
    fn meet(&self, other: &Self) -> Option<Self> {
        let a = self.canonical().ok()?;
        let b = other.canonical().ok()?;
        use ValueAst::*;
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
    fn join(&self, other: &Self) -> Self {
        let a = self
            .canonical()
            .unwrap_or(Cow::Owned(ValueAst::Undetermined));
        let b = other
            .canonical()
            .unwrap_or(Cow::Owned(ValueAst::Undetermined));
        use ValueAst::*;
        match (a.as_ref(), b.as_ref()) {
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
        }
    }

    /// Partial-order check `target ⊑ self` (pattern admits every value target
    /// admits), specialized for the literal cases so the matcher's per-candidate
    /// path allocates nothing. `Term`/`Predicate` on either side fall back to the
    /// canonicalizing `meet`-derived default, which this must equal.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Term(_) | Self::Predicate(_), _) | (_, Self::Term(_) | Self::Predicate(_)) => {
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
            (Self::RangeFrom(_) | Self::RangeTo(_), _) | (_, Self::RangeFrom(_) | Self::RangeTo(_)) => {
                false
            }
        }
    }
}

/// A non-empty canonical set as a `ValueAst`: a singleton collapses to `Lit`.
fn litset_or_lit(set: BTreeSet<i64>) -> ValueAst {
    match set.len() {
        0 => ValueAst::Undetermined,
        1 => ValueAst::Lit(*set.iter().next().unwrap()),
        _ => ValueAst::LitSet(Box::new(set)),
    }
}

// Arithmetic on `ValueAst` propagates `Undetermined`: only `Lit op Lit` yields a
// `Lit`. Every binop has impls for all four `(owned|ref) × (owned|ref)`
// combinations delegating to the ref-ref form, plus a bare `i64` on either side.
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_set(ValueAst::lit_set([2, 1, 2]), ValueAst::LitSet(Box::new(BTreeSet::from([1, 2]))))]
    #[case::var(ValueAst::var("x"), ValueAst::Term(Box::new(ValueTerm::Var("x".to_string()))))]
    #[case::rel_predicate(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("r".to_string()), RelOp::Ge, ValueTerm::Lit(1))), ValueAst::Predicate(Box::new(ValuePredicate::Rel(ValueTerm::Var("r".to_string()), RelOp::Ge, ValueTerm::Lit(1)))))]
    #[case::term(ValueAst::term(ValueTerm::Var("x".to_string())), ValueAst::Term(Box::new(ValueTerm::Var("x".to_string()))))]
    #[case::predicate(ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2]))), ValueAst::Predicate(Box::new(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])))))]
    fn test_value_ast_constructors(#[case] actual: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::i64(ValueAst::from(5_i64), ValueAst::Lit(5))]
    #[case::spin_multiplicity(ValueAst::from(SpinMultiplicity::Triplet), ValueAst::Lit(3))]
    #[case::vec(ValueAst::from(vec![2, 1, 2]), ValueAst::LitSet(Box::new(BTreeSet::from([1, 2]))))]
    #[case::btreeset(ValueAst::from(BTreeSet::from([3, 1])), ValueAst::LitSet(Box::new(BTreeSet::from([1, 3]))))]
    fn test_value_ast_from(#[case] actual: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::litset_singleton(ValueAst::lit_set([7]), Ok(ValueAst::Lit(7)))]
    #[case::litset_empty(ValueAst::lit_set([]), Err(Contradiction))]
    #[case::term_ground(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Lit(2), ValueTerm::Lit(3)])), Ok(ValueAst::Lit(5)))]
    #[case::term_neg_lit(ValueAst::term(ValueTerm::Neg(Box::new(ValueTerm::Lit(4)))), Ok(ValueAst::Lit(-4)))]
    #[case::term_neg_neg(ValueAst::term(ValueTerm::Neg(Box::new(ValueTerm::Neg(Box::new(ValueTerm::Var("x".to_string())))))), Ok(ValueAst::term(ValueTerm::Var("x".to_string()))))]
    #[case::term_sum_identity(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Var("x".to_string()), ValueTerm::Lit(0)])), Ok(ValueAst::term(ValueTerm::Var("x".to_string()))))]
    #[case::term_sum_sorted_const_first(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Var("x".to_string()), ValueTerm::Lit(1)])), Ok(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Lit(1), ValueTerm::Var("x".to_string())]))))]
    #[case::term_sum_flatten(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Sum(vec![ValueTerm::Var("x".to_string()), ValueTerm::Lit(1)]), ValueTerm::Lit(2)])), Ok(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Lit(3), ValueTerm::Var("x".to_string())]))))]
    #[case::term_product_annihilator(ValueAst::term(ValueTerm::Product(vec![ValueTerm::Var("x".to_string()), ValueTerm::Lit(0)])), Ok(ValueAst::Lit(0)))]
    #[case::term_product_identity(ValueAst::term(ValueTerm::Product(vec![ValueTerm::Var("x".to_string()), ValueTerm::Lit(1)])), Ok(ValueAst::term(ValueTerm::Var("x".to_string()))))]
    #[case::term_div_fold(ValueAst::term(ValueTerm::Div(Box::new(ValueTerm::Lit(10)), Box::new(ValueTerm::Lit(3)))), Ok(ValueAst::Lit(3)))]
    #[case::term_rem_fold(ValueAst::term(ValueTerm::Rem(Box::new(ValueTerm::Lit(10)), Box::new(ValueTerm::Lit(3)))), Ok(ValueAst::Lit(1)))]
    #[case::pred_rel_true(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Lit(1), RelOp::Eq, ValueTerm::Lit(1))), Ok(ValueAst::Undetermined))]
    #[case::pred_rel_false(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Lit(1), RelOp::Eq, ValueTerm::Lit(2))), Err(Contradiction))]
    #[case::pred_rel_orient_ge(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("x".to_string()), RelOp::Ge, ValueTerm::Lit(1))), Ok(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Lit(1), RelOp::Le, ValueTerm::Var("x".to_string())))))]
    #[case::pred_rel_eq_sorted(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("x".to_string()), RelOp::Eq, ValueTerm::Lit(0))), Ok(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Lit(0), RelOp::Eq, ValueTerm::Var("x".to_string())))))]
    #[case::pred_not_eq_to_ne(ValueAst::predicate(ValuePredicate::Not(Box::new(ValuePredicate::Rel(ValueTerm::Var("x".to_string()), RelOp::Eq, ValueTerm::Lit(0))))), Ok(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Lit(0), RelOp::Ne, ValueTerm::Var("x".to_string())))))]
    #[case::pred_mem_singleton(ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([5]))), Ok(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Lit(5), RelOp::Eq, ValueTerm::Var("x".to_string())))))]
    #[case::pred_mem_notin_empty(ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::NotIn, BTreeSet::new())), Ok(ValueAst::Undetermined))]
    #[case::pred_mem_in_empty(ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::new())), Err(Contradiction))]
    #[case::pred_not_mem_to_notin(ValueAst::predicate(ValuePredicate::Not(Box::new(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2]))))), Ok(ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::NotIn, BTreeSet::from([1, 2])))))]
    #[case::pred_and_drops_top(ValueAst::predicate(ValuePredicate::And(vec![ValuePredicate::Rel(ValueTerm::Var("x".to_string()), RelOp::Le, ValueTerm::Lit(3)), ValuePredicate::Rel(ValueTerm::Lit(1), RelOp::Eq, ValueTerm::Lit(1))])), Ok(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("x".to_string()), RelOp::Le, ValueTerm::Lit(3)))))]
    #[case::pred_demorgan(ValueAst::predicate(ValuePredicate::Not(Box::new(ValuePredicate::And(vec![ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), ValuePredicate::Mem(ValueTerm::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4]))])))), Ok(ValueAst::predicate(ValuePredicate::Or(vec![ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::NotIn, BTreeSet::from([1, 2])), ValuePredicate::Mem(ValueTerm::Var("y".to_string()), MemOp::NotIn, BTreeSet::from([3, 4]))]))))]
    #[case::term_sum_neg_const(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Var("x".to_string()), ValueTerm::Lit(-3)])), Ok(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Var("x".to_string()), ValueTerm::Neg(Box::new(ValueTerm::Lit(3)))]))))]
    #[case::term_neg_zero(ValueAst::term(ValueTerm::Neg(Box::new(ValueTerm::Lit(0)))), Ok(ValueAst::Lit(0)))]
    #[case::term_product_flatten(ValueAst::term(ValueTerm::Product(vec![ValueTerm::Product(vec![ValueTerm::Var("x".to_string()), ValueTerm::Var("y".to_string())]), ValueTerm::Var("z".to_string())])), Ok(ValueAst::term(ValueTerm::Product(vec![ValueTerm::Var("x".to_string()), ValueTerm::Var("y".to_string()), ValueTerm::Var("z".to_string())]))))]
    #[case::term_product_sort(ValueAst::term(ValueTerm::Product(vec![ValueTerm::Var("b".to_string()), ValueTerm::Var("a".to_string())])), Ok(ValueAst::term(ValueTerm::Product(vec![ValueTerm::Var("a".to_string()), ValueTerm::Var("b".to_string())]))))]
    #[case::term_product_const_fold(ValueAst::term(ValueTerm::Product(vec![ValueTerm::Lit(2), ValueTerm::Lit(3), ValueTerm::Var("x".to_string())])), Ok(ValueAst::term(ValueTerm::Product(vec![ValueTerm::Lit(6), ValueTerm::Var("x".to_string())]))))]
    #[case::term_sum_empty(ValueAst::term(ValueTerm::Sum(vec![])), Ok(ValueAst::Lit(0)))]
    #[case::term_product_empty(ValueAst::term(ValueTerm::Product(vec![])), Ok(ValueAst::Lit(1)))]
    #[case::term_div_by_zero(ValueAst::term(ValueTerm::Div(Box::new(ValueTerm::Lit(10)), Box::new(ValueTerm::Lit(0)))), Ok(ValueAst::term(ValueTerm::Div(Box::new(ValueTerm::Lit(10)), Box::new(ValueTerm::Lit(0))))))]
    #[case::pred_and_flatten(ValueAst::predicate(ValuePredicate::And(vec![ValuePredicate::And(vec![ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), ValuePredicate::Mem(ValueTerm::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4]))]), ValuePredicate::Mem(ValueTerm::Var("z".to_string()), MemOp::In, BTreeSet::from([5, 6]))])), Ok(ValueAst::predicate(ValuePredicate::And(vec![ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), ValuePredicate::Mem(ValueTerm::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4])), ValuePredicate::Mem(ValueTerm::Var("z".to_string()), MemOp::In, BTreeSet::from([5, 6]))]))))]
    #[case::pred_and_sort_dedup(ValueAst::predicate(ValuePredicate::And(vec![ValuePredicate::Mem(ValueTerm::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4])), ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2]))])), Ok(ValueAst::predicate(ValuePredicate::And(vec![ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), ValuePredicate::Mem(ValueTerm::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4]))]))))]
    #[case::pred_and_bottom(ValueAst::predicate(ValuePredicate::And(vec![ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), ValuePredicate::Rel(ValueTerm::Lit(1), RelOp::Eq, ValueTerm::Lit(2))])), Err(Contradiction))]
    #[case::pred_or_drops_bottom(ValueAst::predicate(ValuePredicate::Or(vec![ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), ValuePredicate::Rel(ValueTerm::Lit(1), RelOp::Eq, ValueTerm::Lit(2))])), Ok(ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])))))]
    #[case::pred_or_top(ValueAst::predicate(ValuePredicate::Or(vec![ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])), ValuePredicate::Rel(ValueTerm::Lit(1), RelOp::Eq, ValueTerm::Lit(1))])), Ok(ValueAst::Undetermined))]
    #[case::pred_and_empty(ValueAst::predicate(ValuePredicate::And(vec![])), Ok(ValueAst::Undetermined))]
    #[case::pred_or_empty(ValueAst::predicate(ValuePredicate::Or(vec![])), Err(Contradiction))]
    #[case::pred_not_not(ValueAst::predicate(ValuePredicate::Not(Box::new(ValuePredicate::Not(Box::new(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2]))))))), Ok(ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::In, BTreeSet::from([1, 2])))))]
    #[case::pred_not_le(ValueAst::predicate(ValuePredicate::Not(Box::new(ValuePredicate::Rel(ValueTerm::Var("x".to_string()), RelOp::Le, ValueTerm::Lit(3))))), Ok(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Lit(3), RelOp::Lt, ValueTerm::Var("x".to_string())))))]
    #[case::pred_rel_orient_gt(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("x".to_string()), RelOp::Gt, ValueTerm::Lit(1))), Ok(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Lit(1), RelOp::Lt, ValueTerm::Var("x".to_string())))))]
    #[case::pred_mem_notin_singleton(ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::NotIn, BTreeSet::from([5]))), Ok(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Lit(5), RelOp::Ne, ValueTerm::Var("x".to_string())))))]
    fn test_value_ast_canonicalize(
        #[case] input: ValueAst,
        #[case] expected: Result<ValueAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(ValueAst::Undetermined)]
    #[case::lit(ValueAst::Lit(3))]
    #[case::litset(ValueAst::lit_set([1, 2, 3]))]
    #[case::term_var(ValueAst::term(ValueTerm::Var("x".to_string())))]
    fn test_value_ast_canonicalize_identity(#[case] input: ValueAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::sum(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Var("b".to_string()), ValueTerm::Lit(2), ValueTerm::Var("a".to_string()), ValueTerm::Lit(3)])))]
    #[case::product(ValueAst::term(ValueTerm::Product(vec![ValueTerm::Var("b".to_string()), ValueTerm::Var("a".to_string())])))]
    #[case::rel(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("x".to_string()), RelOp::Ge, ValueTerm::Lit(1))))]
    #[case::or(ValueAst::predicate(ValuePredicate::Or(vec![ValuePredicate::Mem(ValueTerm::Var("y".to_string()), MemOp::In, BTreeSet::from([3, 4])), ValuePredicate::Mem(ValueTerm::Var("x".to_string()), MemOp::NotIn, BTreeSet::from([1, 2]))])))]
    fn test_value_ast_canonicalize_idempotent(#[case] input: ValueAst) {
        let once = input.canonicalize().unwrap();
        let twice = once.clone().canonicalize().unwrap();
        assert_eq!(once, twice);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ValueAst::Lit(3), Some(3))]
    #[case::lit_neg(ValueAst::Lit(-5), Some(-5))]
    #[case::undetermined(ValueAst::Undetermined, None)]
    #[case::litset(ValueAst::lit_set([1, 2]), None)]
    #[case::term(ValueAst::term(ValueTerm::Var("x".to_string())), None)]
    fn test_value_ast_as_lit(#[case] ast: ValueAst, #[case] expected: Option<i64>) {
        assert_eq!(ast.as_lit(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_match(ValueAst::Lit(3), 3, true)]
    #[case::lit_mismatch(ValueAst::Lit(3), 4, false)]
    #[case::undetermined(ValueAst::Undetermined, 3, false)]
    #[case::litset(ValueAst::lit_set([1, 2]), 1, false)]
    fn test_value_ast_as_lit_matches(
        #[case] ast: ValueAst,
        #[case] value: i64,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.as_lit_matches(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(ValueAst::Undetermined, true)]
    #[case::lit(ValueAst::Lit(3), false)]
    #[case::litset(ValueAst::lit_set([1, 2]), false)]
    #[case::term(ValueAst::var("x"), false)]
    #[case::predicate(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("r".to_string()), RelOp::Ge, ValueTerm::Lit(1))), false)]
    #[case::range_from(ValueAst::RangeFrom(1), false)]
    fn test_value_ast_is_undetermined(#[case] ast: ValueAst, #[case] expected: bool) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ValueAst::Undetermined, ValueAst::Lit(3), Some(ValueAst::Lit(3)))]
    #[case::lit_und(ValueAst::Lit(3), ValueAst::Undetermined, Some(ValueAst::Lit(3)))]
    #[case::lit_lit_eq(ValueAst::Lit(3), ValueAst::Lit(3), Some(ValueAst::Lit(3)))]
    #[case::lit_lit_neq(ValueAst::Lit(3), ValueAst::Lit(4), None)]
    #[case::lit_set_in(ValueAst::Lit(2), ValueAst::lit_set([1, 2, 3]), Some(ValueAst::Lit(2)))]
    #[case::lit_set_out(ValueAst::Lit(5), ValueAst::lit_set([1, 2, 3]), None)]
    #[case::set_set_multi(ValueAst::lit_set([1, 2, 3]), ValueAst::lit_set([2, 3, 4]), Some(ValueAst::lit_set([2, 3])))]
    #[case::set_set_singleton(ValueAst::lit_set([1, 2]), ValueAst::lit_set([2, 3]), Some(ValueAst::Lit(2)))]
    #[case::set_set_empty(ValueAst::lit_set([1, 2]), ValueAst::lit_set([3, 4]), None)]
    #[case::term_term_eq(ValueAst::term(ValueTerm::Var("x".to_string())), ValueAst::term(ValueTerm::Var("x".to_string())), Some(ValueAst::term(ValueTerm::Var("x".to_string()))))]
    #[case::term_term_neq(ValueAst::term(ValueTerm::Var("x".to_string())), ValueAst::term(ValueTerm::Var("y".to_string())), None)]
    #[case::term_lit(ValueAst::term(ValueTerm::Var("x".to_string())), ValueAst::Lit(5), None)]
    #[case::rangefrom_lit_in(ValueAst::RangeFrom(1), ValueAst::Lit(2), Some(ValueAst::Lit(2)))]
    #[case::rangefrom_lit_out(ValueAst::RangeFrom(2), ValueAst::Lit(1), None)]
    #[case::rangefrom_rangefrom(ValueAst::RangeFrom(1), ValueAst::RangeFrom(3), Some(ValueAst::RangeFrom(3)))]
    #[case::rangeto_rangeto(ValueAst::RangeTo(5), ValueAst::RangeTo(3), Some(ValueAst::RangeTo(3)))]
    #[case::rangefrom_rangeto_set(ValueAst::RangeFrom(1), ValueAst::RangeTo(4), Some(ValueAst::lit_set([1, 2, 3])))]
    #[case::rangefrom_rangeto_empty(ValueAst::RangeFrom(4), ValueAst::RangeTo(2), None)]
    #[case::rangefrom_set(ValueAst::RangeFrom(2), ValueAst::lit_set([1, 2, 3]), Some(ValueAst::lit_set([2, 3])))]
    fn test_value_ast_meet(
        #[case] a: ValueAst,
        #[case] b: ValueAst,
        #[case] expected: Option<ValueAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ValueAst::Undetermined, ValueAst::Lit(3), ValueAst::Undetermined)]
    #[case::lit_lit_eq(ValueAst::Lit(3), ValueAst::Lit(3), ValueAst::Lit(3))]
    #[case::lit_lit_neq(ValueAst::Lit(3), ValueAst::Lit(4), ValueAst::lit_set([3, 4]))]
    #[case::lit_set(ValueAst::Lit(5), ValueAst::lit_set([1, 2, 3]), ValueAst::lit_set([1, 2, 3, 5]))]
    #[case::set_set(ValueAst::lit_set([1, 2]), ValueAst::lit_set([2, 3]), ValueAst::lit_set([1, 2, 3]))]
    #[case::term_term_eq(ValueAst::term(ValueTerm::Var("x".to_string())), ValueAst::term(ValueTerm::Var("x".to_string())), ValueAst::term(ValueTerm::Var("x".to_string())))]
    #[case::term_term_neq(ValueAst::term(ValueTerm::Var("x".to_string())), ValueAst::term(ValueTerm::Var("y".to_string())), ValueAst::Undetermined)]
    #[case::rangefrom_rangefrom(ValueAst::RangeFrom(3), ValueAst::RangeFrom(1), ValueAst::RangeFrom(1))]
    #[case::rangeto_rangeto(ValueAst::RangeTo(3), ValueAst::RangeTo(5), ValueAst::RangeTo(5))]
    #[case::rangefrom_lit_overapprox(ValueAst::RangeFrom(1), ValueAst::Lit(5), ValueAst::Undetermined)]
    fn test_value_ast_join(#[case] a: ValueAst, #[case] b: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_und(ValueAst::Undetermined, ValueAst::Undetermined, true)]
    #[case::und_lit(ValueAst::Undetermined, ValueAst::Lit(3), true)]
    #[case::lit_und(ValueAst::Lit(3), ValueAst::Undetermined, false)]
    #[case::lit_lit(ValueAst::Lit(3), ValueAst::Lit(3), true)]
    #[case::lit_lit_neq(ValueAst::Lit(3), ValueAst::Lit(4), false)]
    #[case::set_lit_in(ValueAst::lit_set([1, 2, 3]), ValueAst::Lit(2), true)]
    #[case::set_lit_out(ValueAst::lit_set([1, 2, 3]), ValueAst::Lit(5), false)]
    #[case::set_set(ValueAst::lit_set([1, 2, 3]), ValueAst::lit_set([1, 2]), true)]
    #[case::rangefrom_lit_ge(ValueAst::RangeFrom(1), ValueAst::Lit(2), true)]
    #[case::rangefrom_lit_lt(ValueAst::RangeFrom(2), ValueAst::Lit(1), false)]
    #[case::rangefrom_rangefrom_wider(ValueAst::RangeFrom(1), ValueAst::RangeFrom(2), true)]
    #[case::rangefrom_rangefrom_narrower(ValueAst::RangeFrom(2), ValueAst::RangeFrom(1), false)]
    #[case::rangefrom_und(ValueAst::RangeFrom(1), ValueAst::Undetermined, false)]
    #[case::und_rangefrom(ValueAst::Undetermined, ValueAst::RangeFrom(1), true)]
    #[case::rangefrom_set_all_ge(ValueAst::RangeFrom(1), ValueAst::lit_set([2, 3]), true)]
    #[case::rangefrom_set_some_lt(ValueAst::RangeFrom(2), ValueAst::lit_set([1, 3]), false)]
    #[case::rangeto_lit_lt(ValueAst::RangeTo(3), ValueAst::Lit(2), true)]
    #[case::rangeto_lit_ge(ValueAst::RangeTo(2), ValueAst::Lit(3), false)]
    fn test_value_ast_matches(
        #[case] pattern: ValueAst,
        #[case] target: ValueAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
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
    #[case::widen_to_set(ValueAst::Lit(3), ValueAst::Lit(4), true, ValueAst::lit_set([3, 4]))]
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

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(ValueAst::Lit(2), ValueAst::Lit(3), ValueAst::Lit(5))]
    #[case::lit_undetermined(ValueAst::Lit(2), ValueAst::Undetermined, ValueAst::Undetermined)]
    fn test_value_ast_add(#[case] lhs: ValueAst, #[case] rhs: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(lhs + rhs, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(ValueAst::Lit(5), ValueAst::Lit(3), ValueAst::Lit(2))]
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
    fn test_value_ast_div(#[case] lhs: ValueAst, #[case] rhs: ValueAst, #[case] expected: ValueAst) {
        assert_eq!(lhs / rhs, expected);
    }

    #[rstest]
    #[should_panic]
    fn test_value_ast_div_error() {
        let _ = ValueAst::Lit(5) / ValueAst::Lit(0);
    }
}
