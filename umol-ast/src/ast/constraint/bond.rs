//! Localized bond constraints.
use std::collections::BTreeSet;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use strum::{EnumDiscriminants, EnumIter};

use super::super::constraint::ring::{RingMembershipAst, RingScope};
use super::super::remap::IdRemapping;
use super::super::stereo::CisTransStereoAst;
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Localized bond constraint. Held inline on `BondAst` via
/// `BondConstraints`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(BondConstraintKind), derive(Hash, EnumIter))]
pub enum BondConstraint {
    Aromatic,
    RingMembership(RingMembershipAst),
    CisTransStereo(CisTransStereoAst),
}

impl BondConstraint {
    pub fn ring_membership(scope: RingScope, count: impl Into<ValueAst>) -> Self {
        Self::RingMembership(RingMembershipAst::new(scope, count))
    }

    pub fn cis_trans_stereo(c: CisTransStereoAst) -> Self {
        Self::CisTransStereo(c)
    }

    pub fn kind(&self) -> BondConstraintKind {
        self.into()
    }

    /// `false` for `RingMembership` (several per bond, one per `RingScope`); `true` otherwise.
    pub fn is_unique(&self) -> bool {
        self.kind() != BondConstraintKind::RingMembership
    }

    /// `Aromatic` is a flag with no value; the value-carrying variants are
    /// undetermined iff their inner value is undetermined.
    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Aromatic => false,
            Self::RingMembership(m) => m.count.is_undetermined(),
            Self::CisTransStereo(c) => c.is_undetermined(),
        }
    }

    /// Recursively simplify the contained value; the constraint kind is
    /// preserved.
    pub fn simplify(self) -> Self {
        match self {
            Self::Aromatic => Self::Aromatic,
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipAst::new(m.scope, m.count.simplify()))
            }
            Self::CisTransStereo(c) => Self::CisTransStereo(c.clone().canonicalize().unwrap_or(c)),
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

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &ValueAst)> {
        self.get_all(BondConstraintKind::RingMembership)
            .filter_map(|c| match c {
                BondConstraint::RingMembership(m) => Some((m.scope, &m.count)),
                _ => None,
            })
    }

    fn ring_value(&self, scope: RingScope) -> ValueAst {
        self.ring_memberships()
            .find(|(s, _)| *s == scope)
            .map(|(_, v)| v.clone())
            .unwrap_or(ValueAst::Undetermined)
    }

    pub fn ring_count(&self) -> ValueAst {
        self.ring_value(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> ValueAst {
        self.ring_value(RingScope::Size(s))
    }

    pub fn cis_trans_stereo(&self) -> CisTransStereoAst {
        match self.get(BondConstraintKind::CisTransStereo) {
            Some(BondConstraint::CisTransStereo(c)) => c.clone(),
            _ => CisTransStereoAst::Undetermined,
        }
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
    /// one entry; `RingMembership` may yield several (one per `RingScope`).
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
            BondConstraint::RingMembership(m) => m.count.is_undetermined(),
            BondConstraint::CisTransStereo(c) => c.is_undetermined(),
        })
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| match c {
            BondConstraint::Aromatic => true,
            BondConstraint::RingMembership(m) => m.count.is_ground(),
            BondConstraint::CisTransStereo(c) => c.is_ground(),
        })
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let mut result = Self::new();
        if self.aromatic() || other.aromatic() {
            result.add(BondConstraint::Aromatic);
        }
        let cts = self.cis_trans_stereo().meet(&other.cis_trans_stereo())?;
        if !cts.is_undetermined() {
            result.add(BondConstraint::CisTransStereo(cts));
        }
        let mut scopes: BTreeSet<RingScope> = self.ring_memberships().map(|(s, _)| s).collect();
        scopes.extend(other.ring_memberships().map(|(s, _)| s));
        for scope in scopes {
            let v = self.ring_value(scope).meet(&other.ring_value(scope))?;
            if !v.is_undetermined() {
                result.add(BondConstraint::RingMembership(RingMembershipAst::new(
                    scope, v,
                )));
            }
        }
        Some(result)
    }

    fn join(&self, other: &Self) -> Self {
        let mut result = Self::new();
        if self.aromatic() && other.aromatic() {
            result.add(BondConstraint::Aromatic);
        }
        if self.contains(BondConstraintKind::CisTransStereo)
            && other.contains(BondConstraintKind::CisTransStereo)
        {
            let joined = self.cis_trans_stereo().join(&other.cis_trans_stereo());
            if !joined.is_undetermined() {
                result.add(BondConstraint::CisTransStereo(joined));
            }
        }
        for (scope, v) in self.ring_memberships() {
            if other.ring_memberships().any(|(s, _)| s == scope) {
                let j = v.join(&other.ring_value(scope));
                if !j.is_undetermined() {
                    result.add(BondConstraint::RingMembership(RingMembershipAst::new(
                        scope, j,
                    )));
                }
            }
        }
        result
    }

    /// `Aromatic` is a flag; pattern requires it iff target also has it.
    fn matches(&self, target: &Self) -> bool {
        (!self.aromatic() || target.aromatic())
            && self.cis_trans_stereo().matches(&target.cis_trans_stereo())
            && self
                .ring_memberships()
                .all(|(scope, v)| v.matches(&target.ring_value(scope)))
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
    use crate::ast::stereo::{StereoCosetAst, StereoTerm};
    use crate::ast::value::ValueTerm;

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_membership_all(BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::All, ValueAst::Lit(1)))]
    #[case::ring_membership_size(BondConstraint::ring_membership(RingScope::Size(6), 1), BondConstraint::ring_membership(RingScope::Size(6), 1))]
    #[case::cis_trans_stereo(BondConstraint::cis_trans_stereo(CisTransStereoAst::NotStereo), BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo))]
    fn test_bond_constraint_constructors(
        #[case] actual: BondConstraint,
        #[case] expected: BondConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic, BondConstraintKind::Aromatic)]
    #[case::ring_membership_all(BondConstraint::ring_membership(RingScope::All, 1), BondConstraintKind::RingMembership)]
    #[case::ring_membership_size(BondConstraint::ring_membership(RingScope::Size(6), 1), BondConstraintKind::RingMembership)]
    #[case::cis_trans_stereo(BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo), BondConstraintKind::CisTransStereo)]
    fn test_bond_constraint_kind(#[case] c: BondConstraint, #[case] expected: BondConstraintKind) {
        assert_eq!(c.kind(), expected);
    }

    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic, true)]
    #[case::ring_membership(BondConstraint::ring_membership(RingScope::Size(6), 1), false)]
    #[case::cis_trans_stereo(BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo), true)]
    fn test_bond_constraint_is_unique(#[case] c: BondConstraint, #[case] expected: bool) {
        assert_eq!(c.is_unique(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic, false)]
    #[case::ring_membership_all_lit(BondConstraint::ring_membership(RingScope::All, 1), false)]
    #[case::ring_membership_all_undetermined(BondConstraint::ring_membership(RingScope::All, ValueAst::Undetermined), true)]
    #[case::ring_membership_size_lit(BondConstraint::ring_membership(RingScope::Size(6), 1), false)]
    #[case::ring_membership_size_undetermined(BondConstraint::ring_membership(RingScope::All, ValueAst::Undetermined), true)]
    #[case::cis_trans_not_stereo(BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo), false)]
    #[case::cis_trans_undetermined(BondConstraint::CisTransStereo(CisTransStereoAst::Undetermined), true)]
    fn test_bond_constraint_is_undetermined(#[case] c: BondConstraint, #[case] expected: bool) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_membership_all_folds_expr(BondConstraint::ring_membership(RingScope::All, ValueAst::term(ValueTerm::Lit(2))),
        BondConstraint::ring_membership(RingScope::All, 2))]
    #[case::ring_membership_size_folds_expr(BondConstraint::RingMembership(RingMembershipAst { scope: RingScope::Size(6), count: ValueAst::term(ValueTerm::Lit(1)) }),
        BondConstraint::ring_membership(RingScope::Size(6), 1))]
    #[case::cis_trans_lifts_term(BondConstraint::CisTransStereo(CisTransStereoAst::Stereo(StereoCosetAst::term(StereoTerm::Lit(1)))),
        BondConstraint::cis_trans_stereo(CisTransStereoAst::stereo(1_u32)))]
    fn test_bond_constraint_simplify(
        #[case] input: BondConstraint,
        #[case] expected: BondConstraint,
    ) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic)]
    #[case::ring_membership_all_lit(BondConstraint::ring_membership(RingScope::All, 1))]
    #[case::ring_membership_size_undetermined(BondConstraint::ring_membership(
        RingScope::All,
        ValueAst::Undetermined
    ))]
    #[case::cis_trans_not_stereo(BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo))]
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
    #[case::ring_membership_present(BondConstraintKind::RingMembership, true)]
    #[case::cis_trans_absent(BondConstraintKind::CisTransStereo, false)]
    fn test_bond_constraints_contains(
        #[case] kind: BondConstraintKind,
        #[case] expected: bool,
    ) {
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintKind::Aromatic, Some(BondConstraint::Aromatic))]
    #[case::ring_membership(BondConstraintKind::RingMembership, Some(BondConstraint::ring_membership(RingScope::Size(6), 1)))]
    #[case::cis_trans_absent(BondConstraintKind::CisTransStereo, None)]
    fn test_bond_constraints_get(
        #[case] kind: BondConstraintKind,
        #[case] expected: Option<BondConstraint>,
    ) {
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(kind), expected.as_ref());
    }

    #[rstest]
    fn test_bond_constraints_get_mut() {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let entry = cs.get_mut(BondConstraintKind::RingMembership).unwrap();
        *entry = BondConstraint::ring_membership(RingScope::Size(5), 1);
        assert_eq!(
            cs.as_slice(),
            &[
                BondConstraint::Aromatic,
                BondConstraint::ring_membership(RingScope::Size(5), 1),
            ],
        );
    }

    #[rstest]
    fn test_bond_constraints_get_mut_absent() {
        let mut cs = BondConstraints::from_iter([BondConstraint::Aromatic]);
        assert!(cs.get_mut(BondConstraintKind::RingMembership).is_none());
    }

    #[rstest]
    fn test_bond_constraints_iter() {
        let cs = BondConstraints::from_iter([
            BondConstraint::ring_membership(RingScope::Size(6), 1),
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::All, 1),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                BondConstraint::ring_membership(RingScope::Size(6), 1),
                BondConstraint::Aromatic,
                BondConstraint::ring_membership(RingScope::All, 1),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![BondConstraint::Aromatic], vec![None], vec![BondConstraint::Aromatic])]
    #[case::replace_value_kind(vec![BondConstraint::cis_trans_stereo(CisTransStereoAst::Undetermined), BondConstraint::cis_trans_stereo(CisTransStereoAst::NotStereo)],
        vec![None, Some(BondConstraint::cis_trans_stereo(CisTransStereoAst::Undetermined))], vec![BondConstraint::cis_trans_stereo(CisTransStereoAst::NotStereo)])]
    #[case::replace_unit_variant(vec![BondConstraint::Aromatic, BondConstraint::Aromatic],
        vec![None, Some(BondConstraint::Aromatic)], vec![BondConstraint::Aromatic])]
    #[case::append_multi_valued_ring(vec![BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::All, 2)],
        vec![None, None], vec![BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::All, 2)])]
    #[case::distinct_kinds(vec![BondConstraint::Aromatic, BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::Size(6), 1)],
        vec![None, None, None], vec![BondConstraint::Aromatic, BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::Size(6), 1)])]
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
    #[case::partial(|c: &BondConstraint| matches!(c, BondConstraint::Aromatic) || matches!(c, BondConstraint::RingMembership(m) if m.scope == RingScope::Size(6)), vec![
            BondConstraint::Aromatic, BondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::all_dropped(|_: &BondConstraint| false, vec![])]
    fn test_bond_constraints_retain(
        #[case] predicate: impl FnMut(&BondConstraint) -> bool,
        #[case] expected: Vec<BondConstraint>,
    ) {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::All, 1),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
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
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(
            drained,
            vec![
                BondConstraint::Aromatic,
                BondConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
        assert_eq!(cs, BondConstraints::new());
    }

    #[rstest]
    fn test_bond_constraints_simplify_each() {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::All, ValueAst::term(ValueTerm::Lit(1))),
            BondConstraint::RingMembership(RingMembershipAst {
                scope: RingScope::Size(6),
                count: ValueAst::term(ValueTerm::Lit(1)),
            }),
        ]);
        cs.simplify_each();
        assert_eq!(
            cs,
            BondConstraints::from_iter([
                BondConstraint::Aromatic,
                BondConstraint::ring_membership(RingScope::All, 1),
                BondConstraint::ring_membership(RingScope::Size(6), 1),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(BondConstraintKind::Aromatic, Some(BondConstraint::Aromatic), vec![BondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::ring_membership_present(BondConstraintKind::RingMembership, Some(BondConstraint::ring_membership(RingScope::Size(6), 1)), vec![BondConstraint::Aromatic])]
    #[case::cis_trans_absent(BondConstraintKind::CisTransStereo, None, vec![BondConstraint::Aromatic, BondConstraint::ring_membership(RingScope::Size(6), 1)])]
    fn test_bond_constraints_remove(
        #[case] kind: BondConstraintKind,
        #[case] expected_returned: Option<BondConstraint>,
        #[case] expected_state: Vec<BondConstraint>,
    ) {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.remove(kind), expected_returned);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rstest]
    fn test_bond_constraints_remap() {
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let remap = IdRemapping::new(
            Remapping::new(vec![0, 1, 2], vec![0, 1]),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(cs.clone().remap(&remap), cs);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(BondConstraints::new(), BondConstraints::from_iter([BondConstraint::Aromatic]), true)]
    #[case::aromatic_required_present(BondConstraints::from_iter([BondConstraint::Aromatic]),
        BondConstraints::from_iter([BondConstraint::Aromatic]), true)]
    #[case::aromatic_required_absent(BondConstraints::from_iter([BondConstraint::Aromatic]), BondConstraints::new(), false)]
    #[case::ring_membership_all_wildcard_matches_lit(BondConstraints::from_iter([BondConstraint::ring_membership(RingScope::All, ValueAst::Undetermined)]),
        BondConstraints::from_iter([BondConstraint::ring_membership(RingScope::All, 1)]), true)]
    #[case::ring_membership_all_lit_mismatch(BondConstraints::from_iter([BondConstraint::ring_membership(RingScope::All, 1)]),
        BondConstraints::from_iter([BondConstraint::ring_membership(RingScope::All, 2)]), false)]
    #[case::ring_membership_size_subset(BondConstraints::from_iter([BondConstraint::ring_membership(RingScope::Size(5), 1)]),
        BondConstraints::from_iter([BondConstraint::ring_membership(RingScope::Size(5), 1), BondConstraint::ring_membership(RingScope::Size(6), 1)]), true)]
    #[case::ring_membership_size_not_in_target(BondConstraints::from_iter([BondConstraint::ring_membership(RingScope::Size(7), 1)]),
        BondConstraints::from_iter([BondConstraint::ring_membership(RingScope::Size(5), 1)]), false)]
    #[case::cis_trans_match(BondConstraints::from_iter([BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo)]),
        BondConstraints::from_iter([BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo)]), true)]
    #[case::cis_trans_pattern_more_specific(BondConstraints::from_iter([BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo)]),
        BondConstraints::new(), false)]
    fn test_bond_constraints_matches(
        #[case] pattern: BondConstraints,
        #[case] target: BondConstraints,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(vec![BondConstraint::Aromatic, BondConstraint::ring_membership(RingScope::All, 1)],
        vec![BondConstraint::Aromatic, BondConstraint::ring_membership(RingScope::All, 1)])]
    #[case::unique_kind_last_wins(vec![BondConstraint::cis_trans_stereo(CisTransStereoAst::Undetermined), BondConstraint::cis_trans_stereo(CisTransStereoAst::NotStereo)],
        vec![BondConstraint::cis_trans_stereo(CisTransStereoAst::NotStereo)])]
    #[case::ring_appends(vec![BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::Size(6), 1)],
        vec![BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::Size(6), 1)])]
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
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                BondConstraint::Aromatic,
                BondConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rstest]
    fn test_bond_constraints_from_bond_constraint() {
        let cs: BondConstraints = BondConstraint::Aromatic.into();
        assert_eq!(cs.as_slice(), &[BondConstraint::Aromatic]);
    }

    #[rstest]
    fn test_bond_constraints_from_vec() {
        let cs: BondConstraints = vec![
            BondConstraint::Aromatic,
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]
        .into();
        assert_eq!(
            cs.as_slice(),
            &[
                BondConstraint::Aromatic,
                BondConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }
}
