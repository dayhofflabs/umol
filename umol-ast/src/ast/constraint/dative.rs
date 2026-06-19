//! Dative bond constraints.

use std::collections::BTreeSet;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use strum::EnumDiscriminants;

use super::super::constraint::ring::{RingMembershipAst, RingScope};
use super::super::remap::IdRemapping;
use super::super::traits::Lattice;
use super::super::value::ValueAst;

/// Dative-bond-scope constraint. Held inline on `DativeBondAst` via
/// `DativeBondConstraints`. `Aromatic` flags the dative bond as part of an
/// aromatic system (e.g. the N→B π-donation in borazine, O→B in boroxine,
/// or a C→M coordination spanning a metallaaromatic ring).
#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(DativeBondConstraintKind), derive(Hash))]
pub enum DativeBondConstraint {
    Aromatic,
    RingMembership(RingMembershipAst),
}

impl DativeBondConstraint {
    pub fn ring_membership(scope: RingScope, count: impl Into<ValueAst>) -> Self {
        Self::RingMembership(RingMembershipAst::new(scope, count))
    }

    pub fn kind(&self) -> DativeBondConstraintKind {
        self.into()
    }

    /// `false` for `RingMembership` (several per dative bond, one per `RingScope`); `true` for `Aromatic`.
    pub fn is_unique(&self) -> bool {
        self.kind() != DativeBondConstraintKind::RingMembership
    }

    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Aromatic => false,
            Self::RingMembership(m) => m.count.is_undetermined(),
        }
    }

    pub fn simplify(self) -> Self {
        match self {
            Self::Aromatic => Self::Aromatic,
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipAst::new(m.scope, m.count.simplify()))
            }
        }
    }

    pub fn remap(self, _remap: &IdRemapping) -> Option<Self> {
        // Value-only: no indices to remap.
        Some(self)
    }
}

/// Per-dative-bond constraint container. Enforces the per-variant cardinality
/// policy in [`DativeBondConstraint::is_unique`] on insert.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct DativeBondConstraints(Vec<DativeBondConstraint>);

impl DativeBondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[DativeBondConstraint] {
        &self.0
    }

    pub fn contains(&self, kind: DativeBondConstraintKind) -> bool {
        self.0.iter().any(|c| c.kind() == kind)
    }

    pub fn get(&self, kind: DativeBondConstraintKind) -> Option<&DativeBondConstraint> {
        self.0.iter().find(|c| c.kind() == kind)
    }

    pub fn get_mut(&mut self, kind: DativeBondConstraintKind) -> Option<&mut DativeBondConstraint> {
        self.0.iter_mut().find(|c| c.kind() == kind)
    }

    pub fn aromatic(&self) -> bool {
        self.contains(DativeBondConstraintKind::Aromatic)
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &ValueAst)> {
        self.get_all(DativeBondConstraintKind::RingMembership)
            .filter_map(|c| match c {
                DativeBondConstraint::RingMembership(m) => Some((m.scope, &m.count)),
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

    pub fn iter(&self) -> Iter<'_, DativeBondConstraint> {
        self.0.iter()
    }

    pub fn add(&mut self, c: DativeBondConstraint) -> Option<DativeBondConstraint> {
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
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = DativeBondConstraint>) {
        for constraint in constraints {
            self.add(constraint);
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&DativeBondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = DativeBondConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    /// Simplify each contained constraint's inner value in place.
    pub fn simplify_each(&mut self) {
        for c in self.0.iter_mut() {
            *c = mem::replace(c, DativeBondConstraint::Aromatic).simplify();
        }
    }

    pub fn remove(&mut self, kind: DativeBondConstraintKind) -> Option<DativeBondConstraint> {
        let pos = self.0.iter().position(|c| c.kind() == kind)?;
        Some(self.0.remove(pos))
    }

    /// Remove the first entry exactly equal to `constraint`.
    pub fn remove_entry(
        &mut self,
        constraint: &DativeBondConstraint,
    ) -> Option<DativeBondConstraint> {
        let pos = self.0.iter().position(|c| c == constraint)?;
        Some(self.0.remove(pos))
    }

    /// True if any entry exactly equals `constraint`.
    pub fn contains_entry(&self, constraint: &DativeBondConstraint) -> bool {
        self.0.iter().any(|c| c == constraint)
    }

    /// Iterate over every entry of `kind`. `RingMembership` may yield several
    /// (one per `RingScope`); other kinds at most one.
    pub fn get_all(
        &self,
        kind: DativeBondConstraintKind,
    ) -> impl Iterator<Item = &DativeBondConstraint> {
        self.0.iter().filter(move |c| c.kind() == kind)
    }

    /// Remove every entry of `kind`, returning them in insertion order.
    pub fn remove_all(&mut self, kind: DativeBondConstraintKind) -> Vec<DativeBondConstraint> {
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

    pub fn remap(self, remap: &IdRemapping) -> Self {
        Self(self.0.into_iter().filter_map(|c| c.remap(remap)).collect())
    }
}

impl Lattice for DativeBondConstraints {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| match c {
            DativeBondConstraint::Aromatic => false,
            DativeBondConstraint::RingMembership(m) => m.count.is_undetermined(),
        })
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| match c {
            DativeBondConstraint::Aromatic => true,
            DativeBondConstraint::RingMembership(m) => m.count.is_ground(),
        })
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let mut result = Self::new();
        if self.aromatic() || other.aromatic() {
            result.add(DativeBondConstraint::Aromatic);
        }
        let mut scopes: BTreeSet<RingScope> = self.ring_memberships().map(|(s, _)| s).collect();
        scopes.extend(other.ring_memberships().map(|(s, _)| s));
        for scope in scopes {
            let v = self.ring_value(scope).meet(&other.ring_value(scope))?;
            if !v.is_undetermined() {
                result.add(DativeBondConstraint::RingMembership(
                    RingMembershipAst::new(scope, v),
                ));
            }
        }
        Some(result)
    }

    fn join(&self, other: &Self) -> Self {
        let mut result = Self::new();
        if self.aromatic() && other.aromatic() {
            result.add(DativeBondConstraint::Aromatic);
        }
        for (scope, v) in self.ring_memberships() {
            if other.ring_memberships().any(|(s, _)| s == scope) {
                let j = v.join(&other.ring_value(scope));
                if !j.is_undetermined() {
                    result.add(DativeBondConstraint::RingMembership(
                        RingMembershipAst::new(scope, j),
                    ));
                }
            }
        }
        result
    }

    /// `Aromatic` is a flag; pattern requires it iff target also has it.
    fn matches(&self, target: &Self) -> bool {
        (!self.aromatic() || target.aromatic())
            && self
                .ring_memberships()
                .all(|(scope, v)| v.matches(&target.ring_value(scope)))
    }
}

impl FromIterator<DativeBondConstraint> for DativeBondConstraints {
    fn from_iter<I: IntoIterator<Item = DativeBondConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

impl IntoIterator for DativeBondConstraints {
    type Item = DativeBondConstraint;
    type IntoIter = IntoIter<DativeBondConstraint>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<DativeBondConstraint> for DativeBondConstraints {
    fn from(c: DativeBondConstraint) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<DativeBondConstraint>> for DativeBondConstraints {
    fn from(cs: Vec<DativeBondConstraint>) -> Self {
        Self::from_iter(cs)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Remapping;

    use super::*;
    use crate::ast::value::ValueTerm;

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_membership_all(DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::All, ValueAst::Lit(1)))]
    #[case::ring_membership_size(DativeBondConstraint::ring_membership(RingScope::Size(6), 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1))]
    fn test_dative_bond_constraint_constructors(
        #[case] actual: DativeBondConstraint,
        #[case] expected: DativeBondConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic, DativeBondConstraintKind::Aromatic)]
    #[case::ring_membership_all(DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraintKind::RingMembership)]
    #[case::ring_membership_size(DativeBondConstraint::ring_membership(RingScope::Size(6), 1), DativeBondConstraintKind::RingMembership)]
    fn test_dative_bond_constraint_kind(
        #[case] c: DativeBondConstraint,
        #[case] expected: DativeBondConstraintKind,
    ) {
        assert_eq!(c.kind(), expected);
    }

    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic, true)]
    #[case::ring_membership(DativeBondConstraint::ring_membership(RingScope::Size(6), 1), false)]
    fn test_dative_bond_constraint_is_unique(#[case] c: DativeBondConstraint, #[case] expected: bool) {
        assert_eq!(c.is_unique(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic, false)]
    #[case::ring_membership_all_lit(DativeBondConstraint::ring_membership(RingScope::All, 1), false)]
    #[case::ring_membership_all_undetermined(DativeBondConstraint::ring_membership(RingScope::All, ValueAst::Undetermined), true)]
    #[case::ring_membership_size_lit(DativeBondConstraint::ring_membership(RingScope::Size(6), 1), false)]
    #[case::ring_membership_size_undetermined(DativeBondConstraint::ring_membership(RingScope::All, ValueAst::Undetermined), true)]
    fn test_dative_bond_constraint_is_undetermined(
        #[case] c: DativeBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_membership_all_folds_expr(
        DativeBondConstraint::ring_membership(RingScope::All, ValueAst::term(ValueTerm::Lit(2))),
        DativeBondConstraint::ring_membership(RingScope::All, 2),
    )]
    #[case::ring_membership_size_folds_expr(
        DativeBondConstraint::RingMembership(RingMembershipAst { scope: RingScope::Size(6), count: ValueAst::term(ValueTerm::Lit(1)) }),
        DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
    )]
    fn test_dative_bond_constraint_simplify(
        #[case] input: DativeBondConstraint,
        #[case] expected: DativeBondConstraint,
    ) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic)]
    #[case::ring_membership_all_lit(DativeBondConstraint::ring_membership(RingScope::All, 1))]
    #[case::ring_membership_size_undetermined(DativeBondConstraint::ring_membership(
        RingScope::All,
        ValueAst::Undetermined
    ))]
    fn test_dative_bond_constraint_simplify_identity(#[case] input: DativeBondConstraint) {
        assert_eq!(input.clone().simplify(), input);
    }

    #[rstest]
    fn test_dative_bond_constraints_new() {
        let cs = DativeBondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[DativeBondConstraint]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(DativeBondConstraintKind::Aromatic, true)]
    #[case::ring_membership_present(DativeBondConstraintKind::RingMembership, true)]
    fn test_dative_bond_constraints_contains(
        #[case] kind: DativeBondConstraintKind,
        #[case] expected: bool,
    ) {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintKind::Aromatic, Some(DativeBondConstraint::Aromatic))]
    #[case::ring_membership(DativeBondConstraintKind::RingMembership, Some(DativeBondConstraint::ring_membership(RingScope::Size(6), 1)))]
    fn test_dative_bond_constraints_get(
        #[case] kind: DativeBondConstraintKind,
        #[case] expected: Option<DativeBondConstraint>,
    ) {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(kind), expected.as_ref());
    }

    #[rstest]
    fn test_dative_bond_constraints_get_mut() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let entry = cs
            .get_mut(DativeBondConstraintKind::RingMembership)
            .unwrap();
        *entry = DativeBondConstraint::ring_membership(RingScope::Size(5), 1);
        assert_eq!(
            cs.as_slice(),
            &[
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_membership(RingScope::Size(5), 1)
            ],
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_get_mut_absent() {
        let mut cs = DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic]);
        assert!(cs
            .get_mut(DativeBondConstraintKind::RingMembership)
            .is_none());
    }

    #[rstest]
    fn test_dative_bond_constraints_iter() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(RingScope::All, 1),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_membership(RingScope::All, 1),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![DativeBondConstraint::Aromatic], vec![None], vec![DativeBondConstraint::Aromatic])]
    #[case::append_multi_valued_ring(vec![DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::All, 2)],
        vec![None, None], vec![DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::All, 2)])]
    #[case::replace_unit_variant(vec![DativeBondConstraint::Aromatic, DativeBondConstraint::Aromatic],
        vec![None, Some(DativeBondConstraint::Aromatic)], vec![DativeBondConstraint::Aromatic])]
    #[case::distinct_kinds(vec![DativeBondConstraint::Aromatic, DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)],
        vec![None, None, None], vec![DativeBondConstraint::Aromatic, DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])]
    fn test_dative_bond_constraints_add(
        #[case] sequence: Vec<DativeBondConstraint>,
        #[case] expected_returns: Vec<Option<DativeBondConstraint>>,
        #[case] expected_state: Vec<DativeBondConstraint>,
    ) {
        let mut cs = DativeBondConstraints::new();
        let returns: Vec<_> = sequence.into_iter().map(|c| cs.add(c)).collect();
        assert_eq!(returns, expected_returns);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &DativeBondConstraint| matches!(c, DativeBondConstraint::Aromatic) || matches!(c, DativeBondConstraint::RingMembership(m) if m.scope == RingScope::Size(6)),
        vec![DativeBondConstraint::Aromatic, DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::all_dropped(|_: &DativeBondConstraint| false, vec![])]
    fn test_dative_bond_constraints_retain(
        #[case] predicate: impl FnMut(&DativeBondConstraint) -> bool,
        #[case] expected: Vec<DativeBondConstraint>,
    ) {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(RingScope::All, 1),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        cs.retain(predicate);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_dative_bond_constraints_clear() {
        let mut cs = DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic]);
        cs.clear();
        assert_eq!(cs, DativeBondConstraints::new());
    }

    #[rstest]
    fn test_dative_bond_constraints_take() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(
            drained,
            vec![
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1)
            ],
        );
        assert_eq!(cs, DativeBondConstraints::new());
    }

    #[rstest]
    fn test_dative_bond_constraints_simplify_each() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(
                RingScope::All,
                ValueAst::term(ValueTerm::Lit(1)),
            ),
            DativeBondConstraint::RingMembership(RingMembershipAst {
                scope: RingScope::Size(6),
                count: ValueAst::term(ValueTerm::Lit(1)),
            }),
        ]);
        cs.simplify_each();
        assert_eq!(
            cs,
            DativeBondConstraints::from_iter([
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_membership(RingScope::All, 1),
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(DativeBondConstraintKind::Aromatic, Some(DativeBondConstraint::Aromatic), vec![DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::ring_membership_present(DativeBondConstraintKind::RingMembership, Some(DativeBondConstraint::ring_membership(RingScope::Size(6), 1)), vec![DativeBondConstraint::Aromatic])]
    fn test_dative_bond_constraints_remove(
        #[case] kind: DativeBondConstraintKind,
        #[case] expected_returned: Option<DativeBondConstraint>,
        #[case] expected_state: Vec<DativeBondConstraint>,
    ) {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.remove(kind), expected_returned);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rstest]
    fn test_dative_bond_constraints_remap() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let remap = IdRemapping::new(
            Remapping::new(vec![1], vec![1]),
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
    #[case::distinct(vec![DativeBondConstraint::Aromatic, DativeBondConstraint::ring_membership(RingScope::All, 1)], vec![DativeBondConstraint::Aromatic, DativeBondConstraint::ring_membership(RingScope::All, 1)])]
    #[case::unit_kind_last_wins(vec![DativeBondConstraint::Aromatic, DativeBondConstraint::Aromatic], vec![DativeBondConstraint::Aromatic])]
    #[case::ring_appends(vec![DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)], vec![DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::empty(vec![], vec![])]
    fn test_dative_bond_constraints_from_iter(
        #[case] input: Vec<DativeBondConstraint>,
        #[case] expected: Vec<DativeBondConstraint>,
    ) {
        let cs = DativeBondConstraints::from_iter(input);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_dative_bond_constraints_into_iter() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1)
            ],
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_from_dative_bond_constraint() {
        let cs: DativeBondConstraints = DativeBondConstraint::Aromatic.into();
        assert_eq!(cs.as_slice(), &[DativeBondConstraint::Aromatic]);
    }

    #[rstest]
    fn test_dative_bond_constraints_from_vec() {
        let cs: DativeBondConstraints = vec![
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]
        .into();
        assert_eq!(
            cs.as_slice(),
            &[
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1)
            ],
        );
    }
}
