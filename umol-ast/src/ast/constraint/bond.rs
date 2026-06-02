//! Localized bond constraints.  use std::mem;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use strum::EnumDiscriminants;

use super::super::remap::IdRemapping;
use super::super::traits::Lattice;
use super::super::value::ValueAst;

/// Localized bond constraint. Held inline on `BondAst` via
/// `BondConstraints`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(BondConstraintKind), derive(Hash))]
pub enum BondConstraint {
    Aromatic,
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl BondConstraint {
    pub fn ring_count(v: impl Into<ValueAst>) -> Self {
        Self::RingCount(v.into())
    }

    pub fn ring_size(v: impl Into<ValueAst>) -> Self {
        Self::RingSize(v.into())
    }

    pub fn kind(&self) -> BondConstraintKind {
        self.into()
    }

    /// `false` for variants that may legitimately appear multiple times on
    /// the same bond (currently only `RingSize`, where a bond shared between
    /// fused rings satisfies multiple ring-size assertions simultaneously).
    /// `true` for variants that are single-valued per bond.
    pub fn is_unique(&self) -> bool {
        self.kind() != BondConstraintKind::RingSize
    }

    /// `Aromatic` is a flag with no value. `RingCount` / `RingSize` are
    /// undetermined iff their inner value is undetermined.
    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Aromatic => false,
            Self::RingCount(v) | Self::RingSize(v) => v.is_undetermined(),
        }
    }

    /// Simplify the inner `ValueAst` of `RingCount` / `RingSize`. `Aromatic`
    /// is unchanged.
    pub fn simplify(self) -> Self {
        match self {
            Self::Aromatic => Self::Aromatic,
            Self::RingCount(v) => Self::RingCount(v.simplify()),
            Self::RingSize(v) => Self::RingSize(v.simplify()),
        }
    }
}

/// Per-bond constraint container. Enforces the per-variant cardinality policy
/// in [`BondConstraint::is_unique`] on insert: unique-kind variants replace any
/// existing entry of the same discriminant (last-wins); multi-kind variants
/// append.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct BondConstraints(Vec<BondConstraint>);

impl BondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[BondConstraint] {
        &self.0
    }

    pub fn contains(&self, kind: BondConstraintKind) -> bool {
        self.0.iter().any(|c| c.kind() == kind)
    }

    pub fn get(&self, kind: BondConstraintKind) -> Option<&BondConstraint> {
        self.0.iter().find(|c| c.kind() == kind)
    }

    pub fn get_mut(&mut self, kind: BondConstraintKind) -> Option<&mut BondConstraint> {
        self.0.iter_mut().find(|c| c.kind() == kind)
    }

    pub fn aromatic(&self) -> bool {
        self.contains(BondConstraintKind::Aromatic)
    }

    pub fn ring_count(&self) -> ValueAst {
        match self.get(BondConstraintKind::RingCount) {
            Some(BondConstraint::RingCount(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    /// Multi-valued ring-size assertions; a bond shared between fused rings
    /// may carry several. Iterator yields entries in insertion order; empty
    /// if none.
    pub fn ring_sizes(&self) -> impl Iterator<Item = &ValueAst> {
        self.get_all(BondConstraintKind::RingSize)
            .filter_map(|c| match c {
                BondConstraint::RingSize(v) => Some(v),
                _ => None,
            })
    }

    pub fn iter(&self) -> Iter<'_, BondConstraint> {
        self.0.iter()
    }

    /// Insert a constraint per the per-variant cardinality policy. Returns the
    /// replaced entry if `c.is_unique()` and a same-discriminant entry already
    /// existed; `None` otherwise.
    pub fn add(&mut self, c: BondConstraint) -> Option<BondConstraint> {
        if c.is_unique() {
            if let Some(pos) = self
                .0
                .iter()
                .position(|e| mem::discriminant(e) == mem::discriminant(&c))
            {
                return Some(mem::replace(&mut self.0[pos], c));
            }
        }
        self.0.push(c);
        None
    }

    /// Add multiple constraints at once, using semantics of `add`.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = BondConstraint>) {
        for constraint in constraints {
            self.add(constraint);
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&BondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = BondConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    /// Simplify each contained constraint's inner value in place.
    pub fn simplify_each(&mut self) {
        for c in self.0.iter_mut() {
            *c = mem::replace(c, BondConstraint::Aromatic).simplify();
        }
    }

    pub fn remove(&mut self, kind: BondConstraintKind) -> Option<BondConstraint> {
        let pos = self.0.iter().position(|c| c.kind() == kind)?;
        Some(self.0.remove(pos))
    }

    /// Remove the first entry exactly equal to `constraint`. Returns the
    /// removed entry if found; otherwise `None`.
    pub fn remove_entry(&mut self, constraint: &BondConstraint) -> Option<BondConstraint> {
        let pos = self.0.iter().position(|c| c == constraint)?;
        Some(self.0.remove(pos))
    }

    /// True if any entry exactly equals `constraint`.
    pub fn contains_entry(&self, constraint: &BondConstraint) -> bool {
        self.0.iter().any(|c| c == constraint)
    }

    /// Iterate over every entry of `kind`. Single-valued kinds yield at most
    /// one entry; multi-valued (`RingSize`) may yield several.
    pub fn get_all(&self, kind: BondConstraintKind) -> impl Iterator<Item = &BondConstraint> {
        self.0.iter().filter(move |c| c.kind() == kind)
    }

    /// Remove every entry of `kind`, returning them in insertion order.
    pub fn remove_all(&mut self, kind: BondConstraintKind) -> Vec<BondConstraint> {
        let mut out = Vec::new();
        self.0.retain(|c| {
            if c.kind() == kind {
                out.push(c.clone());
                false
            } else {
                true
            }
        });
        out
    }

    /// No-op: no `BondConstraint` variant carries an entity index.
    pub fn remap(self, _remap: &IdRemapping) -> Self {
        self
    }
}

impl Lattice for BondConstraints {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| match c {
            BondConstraint::Aromatic => false,
            BondConstraint::RingCount(v) | BondConstraint::RingSize(v) => v.is_undetermined(),
        })
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| match c {
            BondConstraint::Aromatic => true,
            BondConstraint::RingCount(v) | BondConstraint::RingSize(v) => v.is_ground(),
        })
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let mut result = Self::new();
        if self.aromatic() || other.aromatic() {
            result.add(BondConstraint::Aromatic);
        }
        let rc = self.ring_count().meet(&other.ring_count())?;
        if !rc.is_undetermined() {
            result.add(BondConstraint::RingCount(rc));
        }
        for v in self.ring_sizes().chain(other.ring_sizes()) {
            let entry = BondConstraint::RingSize(v.clone());
            if !v.is_undetermined() && !result.contains_entry(&entry) {
                result.add(entry);
            }
        }
        Some(result)
    }

    fn join(&self, other: &Self) -> Self {
        let mut result = Self::new();
        if self.aromatic() && other.aromatic() {
            result.add(BondConstraint::Aromatic);
        }
        if self.contains(BondConstraintKind::RingCount)
            && other.contains(BondConstraintKind::RingCount)
        {
            let joined = self.ring_count().join(&other.ring_count());
            if !joined.is_undetermined() {
                result.add(BondConstraint::RingCount(joined));
            }
        }
        for v in self.ring_sizes() {
            let entry = BondConstraint::RingSize(v.clone());
            if other
                .ring_sizes()
                .any(|o| BondConstraint::RingSize(o.clone()) == entry)
            {
                result.add(entry);
            }
        }
        result
    }

    /// `Aromatic` is a flag; pattern requires it iff target also has it.
    /// `RingCount` matches via `ValueAst::matches`. `RingSize` (multi-valued)
    /// requires every `self` assertion to be matchable in `target`.
    fn matches(&self, target: &Self) -> bool {
        (!self.aromatic() || target.aromatic())
            && self.ring_count().matches(&target.ring_count())
            && self
                .ring_sizes()
                .all(|p| target.ring_sizes().any(|t| p.matches(t)))
    }
}

impl FromIterator<BondConstraint> for BondConstraints {
    fn from_iter<I: IntoIterator<Item = BondConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

impl IntoIterator for BondConstraints {
    type Item = BondConstraint;
    type IntoIter = IntoIter<BondConstraint>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<BondConstraint> for BondConstraints {
    fn from(c: BondConstraint) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<BondConstraint>> for BondConstraints {
    fn from(cs: Vec<BondConstraint>) -> Self {
        Self::from_iter(cs)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Remapping;

    use super::*;
    use crate::ast::value::ValueExpr;

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_count(BondConstraint::ring_count(1), BondConstraint::RingCount(ValueAst::Lit(1)))]
    #[case::ring_size(BondConstraint::ring_size(6), BondConstraint::RingSize(ValueAst::Lit(6)))]
    fn test_bond_constraint_constructors(
        #[case] actual: BondConstraint,
        #[case] expected: BondConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic, BondConstraintKind::Aromatic)]
    #[case::ring_count(BondConstraint::ring_count(1), BondConstraintKind::RingCount)]
    #[case::ring_size(BondConstraint::ring_size(6), BondConstraintKind::RingSize)]
    fn test_bond_constraint_kind(#[case] c: BondConstraint, #[case] expected: BondConstraintKind) {
        assert_eq!(c.kind(), expected);
    }

    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic, true)]
    #[case::ring_count(BondConstraint::ring_count(1), true)]
    #[case::ring_size(BondConstraint::ring_size(6), false)]
    fn test_bond_constraint_is_unique(#[case] c: BondConstraint, #[case] expected: bool) {
        assert_eq!(c.is_unique(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic, false)]
    #[case::ring_count_lit(BondConstraint::ring_count(1), false)]
    #[case::ring_count_undetermined(BondConstraint::RingCount(ValueAst::Undetermined), true)]
    #[case::ring_size_lit(BondConstraint::ring_size(6), false)]
    #[case::ring_size_undetermined(BondConstraint::RingSize(ValueAst::Undetermined), true)]
    fn test_bond_constraint_is_undetermined(#[case] c: BondConstraint, #[case] expected: bool) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_count_folds_expr(
        BondConstraint::RingCount(ValueAst::Expr(Box::new(ValueExpr::Lit(2)))),
        BondConstraint::ring_count(2),
    )]
    #[case::ring_size_folds_expr(
        BondConstraint::RingSize(ValueAst::Expr(Box::new(ValueExpr::Lit(6)))),
        BondConstraint::ring_size(6),
    )]
    fn test_bond_constraint_simplify(
        #[case] input: BondConstraint,
        #[case] expected: BondConstraint,
    ) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic)]
    #[case::ring_count_lit(BondConstraint::ring_count(1))]
    #[case::ring_size_undetermined(BondConstraint::RingSize(ValueAst::Undetermined))]
    fn test_bond_constraint_simplify_identity(#[case] input: BondConstraint) {
        assert_eq!(input.clone().simplify(), input);
    }

    #[rstest]
    fn test_bond_constraints_new() {
        let cs = BondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[BondConstraint]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(BondConstraintKind::Aromatic, true)]
    #[case::ring_size_present(BondConstraintKind::RingSize, true)]
    #[case::ring_count_absent(BondConstraintKind::RingCount, false)]
    fn test_bond_constraints_contains(
        #[case] kind: BondConstraintKind,
        #[case] expected: bool,
    ) {
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_size(6),
        ]);
        assert_eq!(cs.contains(kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintKind::Aromatic, Some(BondConstraint::Aromatic))]
    #[case::ring_size(BondConstraintKind::RingSize, Some(BondConstraint::ring_size(6)))]
    #[case::ring_count_absent(BondConstraintKind::RingCount, None)]
    fn test_bond_constraints_get(
        #[case] kind: BondConstraintKind,
        #[case] expected: Option<BondConstraint>,
    ) {
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_size(6),
        ]);
        assert_eq!(cs.get(kind), expected.as_ref());
    }

    #[rstest]
    fn test_bond_constraints_get_mut() {
        let mut cs =
            BondConstraints::from_iter([BondConstraint::Aromatic, BondConstraint::ring_size(6)]);
        let entry = cs.get_mut(BondConstraintKind::RingSize).unwrap();
        *entry = BondConstraint::ring_size(5);
        assert_eq!(
            cs.as_slice(),
            &[BondConstraint::Aromatic, BondConstraint::ring_size(5),],
        );
    }

    #[rstest]
    fn test_bond_constraints_get_mut_absent() {
        let mut cs = BondConstraints::from_iter([BondConstraint::Aromatic]);
        assert!(cs.get_mut(BondConstraintKind::RingCount).is_none());
    }

    #[rstest]
    fn test_bond_constraints_iter() {
        let cs = BondConstraints::from_iter([
            BondConstraint::ring_size(6),
            BondConstraint::Aromatic,
            BondConstraint::ring_count(1),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                BondConstraint::ring_size(6),
                BondConstraint::Aromatic,
                BondConstraint::ring_count(1),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(
        vec![BondConstraint::Aromatic],
        vec![None],
        vec![BondConstraint::Aromatic],
    )]
    #[case::replace_same_kind(
        vec![
            BondConstraint::ring_count(1),
            BondConstraint::ring_count(2),
        ],
        vec![None, Some(BondConstraint::ring_count(1))],
        vec![BondConstraint::ring_count(2)],
    )]
    #[case::replace_unit_variant(
        vec![BondConstraint::Aromatic, BondConstraint::Aromatic],
        vec![None, Some(BondConstraint::Aromatic)],
        vec![BondConstraint::Aromatic],
    )]
    #[case::distinct_kinds(
        vec![
            BondConstraint::Aromatic,
            BondConstraint::ring_count(1),
            BondConstraint::ring_size(6),
        ],
        vec![None, None, None],
        vec![
            BondConstraint::Aromatic,
            BondConstraint::ring_count(1),
            BondConstraint::ring_size(6),
        ],
    )]
    fn test_bond_constraints_add(
        #[case] sequence: Vec<BondConstraint>,
        #[case] expected_returns: Vec<Option<BondConstraint>>,
        #[case] expected_state: Vec<BondConstraint>,
    ) {
        let mut cs = BondConstraints::new();
        let returns: Vec<_> = sequence.into_iter().map(|c| cs.add(c)).collect();
        assert_eq!(returns, expected_returns);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(
        |c: &BondConstraint| matches!(c, BondConstraint::Aromatic | BondConstraint::RingSize(_)),
        vec![
            BondConstraint::Aromatic,
            BondConstraint::ring_size(6),
        ],
    )]
    #[case::all_dropped(
        |_: &BondConstraint| false,
        vec![],
    )]
    fn test_bond_constraints_retain(
        #[case] predicate: impl FnMut(&BondConstraint) -> bool,
        #[case] expected: Vec<BondConstraint>,
    ) {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_count(1),
            BondConstraint::ring_size(6),
        ]);
        cs.retain(predicate);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_bond_constraints_clear() {
        let mut cs = BondConstraints::from_iter([BondConstraint::Aromatic]);
        cs.clear();
        assert_eq!(cs, BondConstraints::new());
    }

    #[rstest]
    fn test_bond_constraints_take() {
        let mut cs =
            BondConstraints::from_iter([BondConstraint::Aromatic, BondConstraint::ring_size(6)]);
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(
            drained,
            vec![BondConstraint::Aromatic, BondConstraint::ring_size(6),],
        );
        assert_eq!(cs, BondConstraints::new());
    }

    #[rstest]
    fn test_bond_constraints_simplify_each() {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::RingCount(ValueAst::Expr(Box::new(ValueExpr::Lit(1)))),
            BondConstraint::RingSize(ValueAst::Expr(Box::new(ValueExpr::Lit(6)))),
        ]);
        cs.simplify_each();
        assert_eq!(
            cs,
            BondConstraints::from_iter([
                BondConstraint::Aromatic,
                BondConstraint::ring_count(1),
                BondConstraint::ring_size(6),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(
        BondConstraintKind::Aromatic,
        Some(BondConstraint::Aromatic),
        vec![BondConstraint::ring_size(6)],
    )]
    #[case::ring_size_present(
        BondConstraintKind::RingSize,
        Some(BondConstraint::ring_size(6)),
        vec![BondConstraint::Aromatic],
    )]
    #[case::ring_count_absent(
        BondConstraintKind::RingCount,
        None,
        vec![
            BondConstraint::Aromatic,
            BondConstraint::ring_size(6),
        ],
    )]
    fn test_bond_constraints_remove(
        #[case] kind: BondConstraintKind,
        #[case] expected_returned: Option<BondConstraint>,
        #[case] expected_state: Vec<BondConstraint>,
    ) {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_size(6),
        ]);
        assert_eq!(cs.remove(kind), expected_returned);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rstest]
    fn test_bond_constraints_remap() {
        let cs =
            BondConstraints::from_iter([BondConstraint::Aromatic, BondConstraint::ring_size(6)]);
        let remap = IdRemapping::new(
            Remapping::new(vec![0, 1, 2], vec![0, 1]),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(cs.clone().remap(&remap), cs);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(
        BondConstraints::new(),
        BondConstraints::from_iter([BondConstraint::Aromatic]),
        true,
    )]
    #[case::aromatic_required_present(
        BondConstraints::from_iter([BondConstraint::Aromatic]),
        BondConstraints::from_iter([BondConstraint::Aromatic]),
        true,
    )]
    #[case::aromatic_required_absent(
        BondConstraints::from_iter([BondConstraint::Aromatic]),
        BondConstraints::new(),
        false,
    )]
    #[case::ring_count_wildcard_matches_lit(
        BondConstraints::from_iter([BondConstraint::RingCount(ValueAst::Undetermined)]),
        BondConstraints::from_iter([BondConstraint::ring_count(1)]),
        true,
    )]
    #[case::ring_count_lit_mismatch(
        BondConstraints::from_iter([BondConstraint::ring_count(1)]),
        BondConstraints::from_iter([BondConstraint::ring_count(2)]),
        false,
    )]
    #[case::ring_size_subset(
        BondConstraints::from_iter([BondConstraint::ring_size(5)]),
        BondConstraints::from_iter([BondConstraint::ring_size(5), BondConstraint::ring_size(6)]),
        true,
    )]
    #[case::ring_size_not_in_target(
        BondConstraints::from_iter([BondConstraint::ring_size(7)]),
        BondConstraints::from_iter([BondConstraint::ring_size(5)]),
        false,
    )]
    fn test_bond_constraints_matches(
        #[case] pattern: BondConstraints,
        #[case] target: BondConstraints,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(
        vec![
            BondConstraint::Aromatic,
            BondConstraint::ring_count(1),
        ],
        vec![
            BondConstraint::Aromatic,
            BondConstraint::ring_count(1),
        ],
    )]
    #[case::same_kind_last_wins(
        vec![
            BondConstraint::ring_count(1),
            BondConstraint::ring_count(2),
        ],
        vec![BondConstraint::ring_count(2)],
    )]
    #[case::empty(vec![], vec![])]
    fn test_bond_constraints_from_iter(
        #[case] input: Vec<BondConstraint>,
        #[case] expected: Vec<BondConstraint>,
    ) {
        let cs = BondConstraints::from_iter(input);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_bond_constraints_into_iter() {
        let cs =
            BondConstraints::from_iter([BondConstraint::Aromatic, BondConstraint::ring_size(6)]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![BondConstraint::Aromatic, BondConstraint::ring_size(6),],
        );
    }

    #[rstest]
    fn test_bond_constraints_from_bond_constraint() {
        let cs: BondConstraints = BondConstraint::Aromatic.into();
        assert_eq!(cs.as_slice(), &[BondConstraint::Aromatic]);
    }

    #[rstest]
    fn test_bond_constraints_from_vec() {
        let cs: BondConstraints =
            vec![BondConstraint::Aromatic, BondConstraint::ring_size(6)].into();
        assert_eq!(
            cs.as_slice(),
            &[BondConstraint::Aromatic, BondConstraint::ring_size(6),],
        );
    }
}
