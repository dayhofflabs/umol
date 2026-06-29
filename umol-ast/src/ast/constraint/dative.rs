//! Dative bond constraints.

use std::collections::BTreeSet;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use strum::EnumDiscriminants;

use super::super::boolean::BooleanAst;
use super::super::constraint::ring::{RingMembershipAst, RingScope};
use super::super::error::Contradiction;
use super::super::remap::IdRemapping;
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Dative-bond-scope constraint. Held inline on `DativeBondAst` via
/// `DativeBondConstraints`. `Aromatic` flags the dative bond as part of an
/// aromatic system (e.g. the N→B π-donation in borazine, O→B in boroxine,
/// or a C→M coordination spanning a metallaaromatic ring).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, EnumDiscriminants)]
#[strum_discriminants(name(DativeBondConstraintKind), derive(Hash))]
pub enum DativeBondConstraint {
    Aromatic(BooleanAst),
    RingMembership(RingMembershipAst),
}

impl DativeBondConstraint {
    pub fn ring_membership(scope: RingScope, count: impl Into<ValueAst>) -> Self {
        Self::RingMembership(RingMembershipAst::new(scope, count))
    }

    pub fn kind(&self) -> DativeBondConstraintKind {
        self.into()
    }

    /// Entry identity for order/dedup: `kind()` plus `RingMembership`'s `RingScope`.
    pub fn key(&self) -> DativeBondConstraintKey {
        match self {
            Self::Aromatic(_) => DativeBondConstraintKey::Aromatic,
            Self::RingMembership(m) => DativeBondConstraintKey::RingMembership(m.scope),
        }
    }

    /// `false` for `RingMembership` (several per dative bond, one per `RingScope`); `true` for `Aromatic`.
    pub fn is_unique(&self) -> bool {
        self.kind() != DativeBondConstraintKind::RingMembership
    }

    /// Each variant is undetermined iff its inner value is undetermined.
    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Aromatic(b) => b.is_undetermined(),
            Self::RingMembership(m) => m.count.is_undetermined(),
        }
    }

    /// The same kind/sub-key with its value made `Undetermined` (the vacuous form).
    pub fn as_undetermined(&self) -> Self {
        match self {
            Self::Aromatic(_) => Self::Aromatic(BooleanAst::Undetermined),
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipAst::new(m.scope, ValueAst::Undetermined))
            }
        }
    }

    pub fn remap(self, _remap: &IdRemapping) -> Option<Self> {
        // Value-only: no indices to remap.
        Some(self)
    }
}

/// Entry identity: discriminant + sub-key. Variant order matches
/// `DativeBondConstraint`, so `Ord` agrees with `kind as u8`; the ring run
/// orders by `RingScope`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DativeBondConstraintKey {
    Aromatic,
    RingMembership(RingScope),
}

impl DativeBondConstraintKey {
    pub fn kind(self) -> DativeBondConstraintKind {
        match self {
            Self::Aromatic => DativeBondConstraintKind::Aromatic,
            Self::RingMembership(_) => DativeBondConstraintKind::RingMembership,
        }
    }
}

impl Canonicalize for DativeBondConstraint {
    /// Canonicalize the inner value; kind and sub-key are preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Aromatic(b) => Self::Aromatic(b.canonicalize()?),
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipAst::new(m.scope, m.count.canonicalize()?))
            }
        })
    }
}

/// Per-dative-bond constraint container, kept `key()`-sorted. On insert, unique
/// kinds replace the same-key entry (last-wins); ring appends at its scope.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

    /// The dative bond's aromatic value, or `Undetermined` when no `Aromatic` constraint is present.
    pub fn aromatic(&self) -> BooleanAst {
        match self.get(DativeBondConstraintKind::Aromatic) {
            Some(DativeBondConstraint::Aromatic(b)) => *b,
            _ => BooleanAst::Undetermined,
        }
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &ValueAst)> {
        self.get_all(DativeBondConstraintKind::RingMembership)
            .filter_map(|c| match c {
                DativeBondConstraint::RingMembership(m) => Some((m.scope, &m.count)),
                _ => None,
            })
    }

    fn ring_membership_value(&self, scope: RingScope) -> ValueAst {
        self.ring_memberships()
            .find(|(s, _)| *s == scope)
            .map(|(_, v)| v.clone())
            .unwrap_or(ValueAst::Undetermined)
    }

    pub fn ring_count(&self) -> ValueAst {
        self.ring_membership_value(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> ValueAst {
        self.ring_membership_value(RingScope::Size(s))
    }

    pub fn iter(&self) -> Iter<'_, DativeBondConstraint> {
        self.0.iter()
    }

    /// Insert at the `key()`-sorted position: unique kinds replace the same-key
    /// entry (returning it); ring appends, leaving duplicates for lazy dedup.
    pub fn add(&mut self, c: DativeBondConstraint) -> Option<DativeBondConstraint> {
        match self.find_by_key(c.key()) {
            Ok(i) if c.is_unique() => Some(mem::replace(&mut self.0[i], c)),
            Ok(i) => {
                let end = i + self.0[i..]
                    .iter()
                    .take_while(|e| e.key() == c.key())
                    .count();
                self.0.insert(end, c);
                None
            }
            Err(i) => {
                self.0.insert(i, c);
                None
            }
        }
    }

    fn find_by_key(&self, key: DativeBondConstraintKey) -> Result<usize, usize> {
        self.0.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains_key(&self, key: DativeBondConstraintKey) -> bool {
        self.find_by_key(key).is_ok()
    }

    pub fn get_by_key(&self, key: DativeBondConstraintKey) -> Option<&DativeBondConstraint> {
        self.find_by_key(key).ok().map(|i| &self.0[i])
    }

    pub fn get_by_key_mut(
        &mut self,
        key: DativeBondConstraintKey,
    ) -> Option<&mut DativeBondConstraint> {
        self.find_by_key(key).ok().map(|i| &mut self.0[i])
    }

    pub fn remove_by_key(&mut self, key: DativeBondConstraintKey) -> Option<DativeBondConstraint> {
        self.find_by_key(key).ok().map(|i| self.0.remove(i))
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

impl Canonicalize for DativeBondConstraints {
    /// Sort by `key()`, canonicalize each value, merge same-scope ring entries
    /// by value-`meet` (`Err` on contradiction), drop vacuous entries.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut input = self.0;
        input.sort_by_key(|c| c.key());
        let mut entries: Vec<DativeBondConstraint> = Vec::new();
        for c in input {
            let c = c.canonicalize()?;
            if let (
                Some(DativeBondConstraint::RingMembership(prev)),
                DativeBondConstraint::RingMembership(next),
            ) = (entries.last_mut(), &c)
            {
                if prev.scope == next.scope {
                    prev.count = prev.count.meet(&next.count).ok_or(Contradiction)?;
                    continue;
                }
            }
            entries.push(c);
        }
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl Lattice for DativeBondConstraints {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| match c {
            DativeBondConstraint::Aromatic(b) => b.is_undetermined(),
            DativeBondConstraint::RingMembership(m) => m.count.is_undetermined(),
        })
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| match c {
            DativeBondConstraint::Aromatic(b) => b.is_ground(),
            DativeBondConstraint::RingMembership(m) => m.count.is_ground(),
        })
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let mut result = Self::new();
        let aromatic = self.aromatic().meet(&other.aromatic())?;
        if !aromatic.is_undetermined() {
            result.add(DativeBondConstraint::Aromatic(aromatic));
        }
        let mut scopes: BTreeSet<RingScope> = self.ring_memberships().map(|(s, _)| s).collect();
        scopes.extend(other.ring_memberships().map(|(s, _)| s));
        for scope in scopes {
            let v = self
                .ring_membership_value(scope)
                .meet(&other.ring_membership_value(scope))?;
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
        let aromatic = self.aromatic().join(&other.aromatic());
        if !aromatic.is_undetermined() {
            result.add(DativeBondConstraint::Aromatic(aromatic));
        }
        for (scope, v) in self.ring_memberships() {
            if other.ring_memberships().any(|(s, _)| s == scope) {
                let j = v.join(&other.ring_membership_value(scope));
                if !j.is_undetermined() {
                    result.add(DativeBondConstraint::RingMembership(
                        RingMembershipAst::new(scope, j),
                    ));
                }
            }
        }
        result
    }

    /// Each value is matched on its own lattice; an empty pattern matches any target.
    fn matches(&self, target: &Self) -> bool {
        self.aromatic().matches(&target.aromatic())
            && self
                .ring_memberships()
                .all(|(scope, v)| v.matches(&target.ring_membership_value(scope)))
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
    #[case::aromatic(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintKind::Aromatic)]
    #[case::ring_membership_all(DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraintKind::RingMembership)]
    #[case::ring_membership_size(DativeBondConstraint::ring_membership(RingScope::Size(6), 1), DativeBondConstraintKind::RingMembership)]
    fn test_dative_bond_constraint_kind(
        #[case] c: DativeBondConstraint,
        #[case] expected: DativeBondConstraintKind,
    ) {
        assert_eq!(c.kind(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintKey::Aromatic)]
    #[case::ring_membership_all(DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraintKey::RingMembership(RingScope::All))]
    #[case::ring_membership_size(DativeBondConstraint::ring_membership(RingScope::Size(6), 1), DativeBondConstraintKey::RingMembership(RingScope::Size(6)))]
    fn test_dative_bond_constraint_key(
        #[case] c: DativeBondConstraint,
        #[case] expected: DativeBondConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintKey::Aromatic, DativeBondConstraintKind::Aromatic)]
    #[case::ring_membership_all(DativeBondConstraintKey::RingMembership(RingScope::All), DativeBondConstraintKind::RingMembership)]
    #[case::ring_membership_size(DativeBondConstraintKey::RingMembership(RingScope::Size(6)), DativeBondConstraintKind::RingMembership)]
    fn test_dative_bond_constraint_key_kind(
        #[case] key: DativeBondConstraintKey,
        #[case] expected: DativeBondConstraintKind,
    ) {
        assert_eq!(key.kind(), expected);
    }

    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), true)]
    #[case::ring_membership(DativeBondConstraint::ring_membership(RingScope::Size(6), 1), false)]
    fn test_dative_bond_constraint_is_unique(
        #[case] c: DativeBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_unique(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), false)]
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
    #[case::aromatic(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), Ok(DativeBondConstraint::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_count_litset_singleton(
        DativeBondConstraint::RingMembership(RingMembershipAst::new(RingScope::All, ValueAst::lit_set([2]))),
        Ok(DativeBondConstraint::ring_membership(RingScope::All, 2)))]
    #[case::empty_litset_contradiction(
        DativeBondConstraint::RingMembership(RingMembershipAst::new(RingScope::All, ValueAst::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_dative_bond_constraint_canonicalize(
        #[case] constraint: DativeBondConstraint,
        #[case] expected: Result<DativeBondConstraint, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
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
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintKind::Aromatic, Some(DativeBondConstraint::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_membership(DativeBondConstraintKind::RingMembership, Some(DativeBondConstraint::ring_membership(RingScope::Size(6), 1)))]
    fn test_dative_bond_constraints_get(
        #[case] kind: DativeBondConstraintKind,
        #[case] expected: Option<DativeBondConstraint>,
    ) {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(kind), expected.as_ref());
    }

    #[rstest]
    fn test_dative_bond_constraints_get_mut() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let entry = cs
            .get_mut(DativeBondConstraintKind::RingMembership)
            .unwrap();
        *entry = DativeBondConstraint::ring_membership(RingScope::Size(5), 1);
        assert_eq!(
            cs.as_slice(),
            &[
                DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraint::ring_membership(RingScope::Size(5), 1)
            ],
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_get_mut_absent() {
        let mut cs = DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(
            BooleanAst::Lit(true),
        )]);
        assert!(cs
            .get_mut(DativeBondConstraintKind::RingMembership)
            .is_none());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(DativeBondConstraintKey::Aromatic, true)]
    #[case::ring_all_present(DativeBondConstraintKey::RingMembership(RingScope::All), true)]
    #[case::ring_size_present(DativeBondConstraintKey::RingMembership(RingScope::Size(6)), true)]
    #[case::ring_size_absent(DativeBondConstraintKey::RingMembership(RingScope::Size(5)), false)]
    fn test_dative_bond_constraints_contains_key(
        #[case] key: DativeBondConstraintKey,
        #[case] expected: bool,
    ) {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::All, 2),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains_key(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintKey::Aromatic, Some(DativeBondConstraint::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_all(DativeBondConstraintKey::RingMembership(RingScope::All), Some(DativeBondConstraint::ring_membership(RingScope::All, 2)))]
    #[case::ring_size(DativeBondConstraintKey::RingMembership(RingScope::Size(6)), Some(DativeBondConstraint::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_size_absent(DativeBondConstraintKey::RingMembership(RingScope::Size(5)), None)]
    fn test_dative_bond_constraints_get_by_key(
        #[case] key: DativeBondConstraintKey,
        #[case] expected: Option<DativeBondConstraint>,
    ) {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::All, 2),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get_by_key(key), expected.as_ref());
    }

    #[rstest]
    fn test_dative_bond_constraints_get_by_key_mut() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::ring_membership(RingScope::All, 2),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let slot = cs
            .get_by_key_mut(DativeBondConstraintKey::RingMembership(RingScope::Size(6)))
            .unwrap();
        *slot = DativeBondConstraint::ring_membership(RingScope::Size(6), 2);
        assert_eq!(
            cs.get_by_key(DativeBondConstraintKey::RingMembership(RingScope::Size(6))),
            Some(&DativeBondConstraint::ring_membership(
                RingScope::Size(6),
                2
            )),
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_remove_by_key() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::All, 2),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let removed = cs.remove_by_key(DativeBondConstraintKey::RingMembership(RingScope::Size(6)));
        assert_eq!(
            removed,
            Some(DativeBondConstraint::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(
            cs.as_slice(),
            &[
                DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraint::ring_membership(RingScope::All, 2),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::merge_same_scope(
        DativeBondConstraints::from_iter([
            DativeBondConstraint::ring_membership(RingScope::All, ValueAst::lit_set([1, 2])),
            DativeBondConstraint::ring_membership(RingScope::All, ValueAst::lit_set([2, 3])),
        ]),
        Ok(DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, 2)])))]
    #[case::drop_vacuous(
        DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::All, ValueAst::Undetermined),
        ]),
        Ok(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))])))]
    #[case::contradiction_same_scope(
        DativeBondConstraints::from_iter([
            DativeBondConstraint::ring_membership(RingScope::All, 1),
            DativeBondConstraint::ring_membership(RingScope::All, 0),
        ]),
        Err(Contradiction))]
    fn test_dative_bond_constraints_canonicalize(
        #[case] constraints: DativeBondConstraints,
        #[case] expected: Result<DativeBondConstraints, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rstest]
    fn test_dative_bond_constraints_iter() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::All, 1),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraint::ring_membership(RingScope::All, 1),
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))], vec![None], vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))])]
    #[case::append_multi_valued_ring(vec![DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::All, 2)],
        vec![None, None], vec![DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::All, 2)])]
    #[case::replace_unit_variant(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::Aromatic(BooleanAst::Lit(true))],
        vec![None, Some(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)))], vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))])]
    #[case::distinct_kinds(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)],
        vec![None, None, None], vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])]
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
    #[case::partial(|c: &DativeBondConstraint| matches!(c, DativeBondConstraint::Aromatic(BooleanAst::Lit(true))) || matches!(c, DativeBondConstraint::RingMembership(m) if m.scope == RingScope::Size(6)),
        vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::all_dropped(|_: &DativeBondConstraint| false, vec![])]
    fn test_dative_bond_constraints_retain(
        #[case] predicate: impl FnMut(&DativeBondConstraint) -> bool,
        #[case] expected: Vec<DativeBondConstraint>,
    ) {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::All, 1),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        cs.retain(predicate);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_dative_bond_constraints_clear() {
        let mut cs = DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(
            BooleanAst::Lit(true),
        )]);
        cs.clear();
        assert_eq!(cs, DativeBondConstraints::new());
    }

    #[rstest]
    fn test_dative_bond_constraints_take() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(
            drained,
            vec![
                DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1)
            ],
        );
        assert_eq!(cs, DativeBondConstraints::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(DativeBondConstraintKind::Aromatic, Some(DativeBondConstraint::Aromatic(BooleanAst::Lit(true))), vec![DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::ring_membership_present(DativeBondConstraintKind::RingMembership, Some(DativeBondConstraint::ring_membership(RingScope::Size(6), 1)), vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))])]
    fn test_dative_bond_constraints_remove(
        #[case] kind: DativeBondConstraintKind,
        #[case] expected_returned: Option<DativeBondConstraint>,
        #[case] expected_state: Vec<DativeBondConstraint>,
    ) {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.remove(kind), expected_returned);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rstest]
    fn test_dative_bond_constraints_remap() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
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
    #[case::distinct(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1)], vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1)])]
    #[case::unit_kind_last_wins(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::Aromatic(BooleanAst::Lit(true))], vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))])]
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
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1)
            ],
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_from_dative_bond_constraint() {
        let cs: DativeBondConstraints =
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)).into();
        assert_eq!(
            cs.as_slice(),
            &[DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_from_vec() {
        let cs: DativeBondConstraints = vec![
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]
        .into();
        assert_eq!(
            cs.as_slice(),
            &[
                DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1)
            ],
        );
    }
}
