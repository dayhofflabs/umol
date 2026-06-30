//! Localized bond constraints.
use std::collections::BTreeSet;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use strum::{EnumDiscriminants, EnumIter};

use super::super::boolean::BooleanAst;
use super::super::constraint::ring::{RingMembershipAst, RingScope};
use super::super::error::Contradiction;
use super::super::remap::IdCompaction;
use super::super::stereo::CisTransStereoAst;
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Localized bond constraint. Held inline on `BondAst` via
/// `BondConstraints`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, EnumDiscriminants)]
#[strum_discriminants(name(BondConstraintKind), derive(Hash, EnumIter))]
pub enum BondConstraint {
    Aromatic(BooleanAst),
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

    /// Entry identity for order/dedup: `kind()` plus `RingMembership`'s `RingScope`.
    pub fn key(&self) -> BondConstraintKey {
        match self {
            Self::Aromatic(_) => BondConstraintKey::Aromatic,
            Self::RingMembership(m) => BondConstraintKey::RingMembership(m.scope),
            Self::CisTransStereo(_) => BondConstraintKey::CisTransStereo,
        }
    }

    /// `false` for `RingMembership` (several per bond, one per `RingScope`); `true` otherwise.
    pub fn is_unique(&self) -> bool {
        self.kind() != BondConstraintKind::RingMembership
    }

    /// Each variant is undetermined iff its inner value is undetermined.
    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Aromatic(b) => b.is_undetermined(),
            Self::RingMembership(m) => m.count.is_undetermined(),
            Self::CisTransStereo(c) => c.is_undetermined(),
        }
    }

    /// The same kind/sub-key with its value made `Undetermined` (the vacuous form). Renders as
    /// `#a*` / `#R*` / `#C*`; a reaction `:modify` uses it to mark a constraint for removal.
    pub fn as_undetermined(&self) -> Self {
        match self {
            Self::Aromatic(_) => Self::Aromatic(BooleanAst::Undetermined),
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipAst::new(m.scope, ValueAst::Undetermined))
            }
            Self::CisTransStereo(_) => Self::CisTransStereo(CisTransStereoAst::Undetermined),
        }
    }
}

/// Entry identity: discriminant + sub-key. Variant order matches `BondConstraint`,
/// so `Ord` agrees with `kind as u8`; the ring run orders by `RingScope`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondConstraintKey {
    Aromatic,
    RingMembership(RingScope),
    CisTransStereo,
}

impl BondConstraintKey {
    pub fn kind(self) -> BondConstraintKind {
        match self {
            Self::Aromatic => BondConstraintKind::Aromatic,
            Self::RingMembership(_) => BondConstraintKind::RingMembership,
            Self::CisTransStereo => BondConstraintKind::CisTransStereo,
        }
    }
}

impl Canonicalize for BondConstraint {
    /// Canonicalize the inner value; kind and sub-key are preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Aromatic(b) => Self::Aromatic(b.canonicalize()?),
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipAst::new(m.scope, m.count.canonicalize()?))
            }
            Self::CisTransStereo(c) => Self::CisTransStereo(c.canonicalize()?),
        })
    }
}

/// Per-bond constraint container, kept `key()`-sorted. On insert, unique kinds
/// replace the same-key entry (last-wins); ring appends at its scope position.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

    /// The bond's aromatic value, or `Undetermined` when no `Aromatic` constraint is present.
    pub fn aromatic(&self) -> BooleanAst {
        match self.get(BondConstraintKind::Aromatic) {
            Some(BondConstraint::Aromatic(b)) => *b,
            _ => BooleanAst::Undetermined,
        }
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &ValueAst)> {
        self.get_all(BondConstraintKind::RingMembership)
            .filter_map(|c| match c {
                BondConstraint::RingMembership(m) => Some((m.scope, &m.count)),
                _ => None,
            })
    }

    fn ring_membership_value(&self, scope: RingScope) -> Option<&ValueAst> {
        self.ring_memberships()
            .find(|(s, _)| *s == scope)
            .map(|(_, v)| v)
    }

    pub fn ring_count(&self) -> Option<&ValueAst> {
        self.ring_membership_value(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> Option<&ValueAst> {
        self.ring_membership_value(RingScope::Size(s))
    }

    pub fn cis_trans_stereo(&self) -> Option<&CisTransStereoAst> {
        match self.get(BondConstraintKind::CisTransStereo) {
            Some(BondConstraint::CisTransStereo(c)) => Some(c),
            _ => None,
        }
    }

    pub fn iter(&self) -> Iter<'_, BondConstraint> {
        self.0.iter()
    }

    /// Insert at the `key()`-sorted position: unique kinds replace the same-key
    /// entry (returning it); ring appends, leaving duplicates for lazy dedup.
    pub fn add(&mut self, c: BondConstraint) -> Option<BondConstraint> {
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

    fn find_by_key(&self, key: BondConstraintKey) -> Result<usize, usize> {
        self.0.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains_key(&self, key: BondConstraintKey) -> bool {
        self.find_by_key(key).is_ok()
    }

    pub fn get_by_key(&self, key: BondConstraintKey) -> Option<&BondConstraint> {
        self.find_by_key(key).ok().map(|i| &self.0[i])
    }

    pub fn get_by_key_mut(&mut self, key: BondConstraintKey) -> Option<&mut BondConstraint> {
        self.find_by_key(key).ok().map(|i| &mut self.0[i])
    }

    pub fn remove_by_key(&mut self, key: BondConstraintKey) -> Option<BondConstraint> {
        self.find_by_key(key).ok().map(|i| self.0.remove(i))
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
    pub fn compact(self, _compaction: &IdCompaction) -> Self {
        self
    }
}

impl Canonicalize for BondConstraints {
    /// Sort by `key()`, canonicalize each value, merge same-scope ring entries
    /// by value-`meet` (`Err` on contradiction), drop vacuous entries.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut input = self.0;
        input.sort_by_key(|c| c.key());
        let mut entries: Vec<BondConstraint> = Vec::new();
        for c in input {
            let c = c.canonicalize()?;
            if let (
                Some(BondConstraint::RingMembership(prev)),
                BondConstraint::RingMembership(next),
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

impl Lattice for BondConstraints {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| match c {
            BondConstraint::Aromatic(b) => b.is_undetermined(),
            BondConstraint::RingMembership(m) => m.count.is_undetermined(),
            BondConstraint::CisTransStereo(c) => c.is_undetermined(),
        })
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| match c {
            BondConstraint::Aromatic(b) => b.is_ground(),
            BondConstraint::RingMembership(m) => m.count.is_ground(),
            BondConstraint::CisTransStereo(c) => c.is_ground(),
        })
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let mut result = Self::new();
        let meet_val = |a: Option<&ValueAst>, b: Option<&ValueAst>| {
            a.unwrap_or(&ValueAst::Undetermined)
                .meet(b.unwrap_or(&ValueAst::Undetermined))
        };
        let aromatic = self.aromatic().meet(&other.aromatic())?;
        if !aromatic.is_undetermined() {
            result.add(BondConstraint::Aromatic(aromatic));
        }
        let cts = self
            .cis_trans_stereo()
            .unwrap_or(&CisTransStereoAst::Undetermined)
            .meet(
                other
                    .cis_trans_stereo()
                    .unwrap_or(&CisTransStereoAst::Undetermined),
            )?;
        if !cts.is_undetermined() {
            result.add(BondConstraint::CisTransStereo(cts));
        }
        let mut scopes: BTreeSet<RingScope> = self.ring_memberships().map(|(s, _)| s).collect();
        scopes.extend(other.ring_memberships().map(|(s, _)| s));
        for scope in scopes {
            let v = meet_val(
                self.ring_membership_value(scope),
                other.ring_membership_value(scope),
            )?;
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
        let aromatic = self.aromatic().join(&other.aromatic());
        if !aromatic.is_undetermined() {
            result.add(BondConstraint::Aromatic(aromatic));
        }
        if self.contains(BondConstraintKind::CisTransStereo)
            && other.contains(BondConstraintKind::CisTransStereo)
        {
            // both present (guarded above), so the accessors are `Some`
            let joined = self
                .cis_trans_stereo()
                .unwrap()
                .join(other.cis_trans_stereo().unwrap());
            if !joined.is_undetermined() {
                result.add(BondConstraint::CisTransStereo(joined));
            }
        }
        for (scope, v) in self.ring_memberships() {
            if other.ring_memberships().any(|(s, _)| s == scope) {
                let j = v.join(other.ring_membership_value(scope).unwrap());
                if !j.is_undetermined() {
                    result.add(BondConstraint::RingMembership(RingMembershipAst::new(
                        scope, j,
                    )));
                }
            }
        }
        result
    }

    /// Pattern-driven: every constraint the pattern carries must match the target,
    /// looked up by reference. Each value is matched on its own lattice; an empty
    /// pattern matches any target.
    fn matches(&self, target: &Self) -> bool {
        self.iter().all(|c| match c {
            BondConstraint::Aromatic(b) => b.matches(&target.aromatic()),
            BondConstraint::RingMembership(rm) => rm.count.matches(
                target
                    .ring_membership_value(rm.scope)
                    .unwrap_or(&ValueAst::Undetermined),
            ),
            BondConstraint::CisTransStereo(cts) => cts.matches(
                target
                    .cis_trans_stereo()
                    .unwrap_or(&CisTransStereoAst::Undetermined),
            ),
        })
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
    use umol_graph_core::Compaction;

    use super::*;
    use crate::ast::stereo::{StereoCosetAst, StereoTerm};
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
    #[case::aromatic(BondConstraint::Aromatic(BooleanAst::Lit(true)), BondConstraintKind::Aromatic)]
    #[case::ring_membership_all(BondConstraint::ring_membership(RingScope::All, 1), BondConstraintKind::RingMembership)]
    #[case::ring_membership_size(BondConstraint::ring_membership(RingScope::Size(6), 1), BondConstraintKind::RingMembership)]
    #[case::cis_trans_stereo(BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo), BondConstraintKind::CisTransStereo)]
    fn test_bond_constraint_kind(#[case] c: BondConstraint, #[case] expected: BondConstraintKind) {
        assert_eq!(c.kind(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic(BooleanAst::Lit(true)), BondConstraintKey::Aromatic)]
    #[case::ring_membership_all(BondConstraint::ring_membership(RingScope::All, 1), BondConstraintKey::RingMembership(RingScope::All))]
    #[case::ring_membership_size(BondConstraint::ring_membership(RingScope::Size(6), 1), BondConstraintKey::RingMembership(RingScope::Size(6)))]
    #[case::cis_trans_stereo(BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo), BondConstraintKey::CisTransStereo)]
    fn test_bond_constraint_key(#[case] c: BondConstraint, #[case] expected: BondConstraintKey) {
        assert_eq!(c.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintKey::Aromatic, BondConstraintKind::Aromatic)]
    #[case::ring_membership_all(BondConstraintKey::RingMembership(RingScope::All), BondConstraintKind::RingMembership)]
    #[case::ring_membership_size(BondConstraintKey::RingMembership(RingScope::Size(6)), BondConstraintKind::RingMembership)]
    #[case::cis_trans_stereo(BondConstraintKey::CisTransStereo, BondConstraintKind::CisTransStereo)]
    fn test_bond_constraint_key_kind(#[case] key: BondConstraintKey, #[case] expected: BondConstraintKind) {
        assert_eq!(key.kind(), expected);
    }

    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic(BooleanAst::Lit(true)), true)]
    #[case::ring_membership(BondConstraint::ring_membership(RingScope::Size(6), 1), false)]
    #[case::cis_trans_stereo(BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo), true)]
    fn test_bond_constraint_is_unique(#[case] c: BondConstraint, #[case] expected: bool) {
        assert_eq!(c.is_unique(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic(BooleanAst::Lit(true)), false)]
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
    #[case::aromatic(BondConstraint::Aromatic(BooleanAst::Lit(true)), BondConstraint::Aromatic(BooleanAst::Undetermined))]
    #[case::ring_membership_keeps_scope(BondConstraint::ring_membership(RingScope::Size(6), 1), BondConstraint::ring_membership(RingScope::Size(6), ValueAst::Undetermined))]
    #[case::cis_trans(BondConstraint::CisTransStereo(CisTransStereoAst::stereo(1_u32)), BondConstraint::CisTransStereo(CisTransStereoAst::Undetermined))]
    fn test_bond_constraint_as_undetermined(#[case] c: BondConstraint, #[case] expected: BondConstraint) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic(BooleanAst::Lit(true)), Ok(BondConstraint::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_count_litset_singleton(
        BondConstraint::RingMembership(RingMembershipAst::new(RingScope::All, ValueAst::lit_set([2]))),
        Ok(BondConstraint::ring_membership(RingScope::All, 2)))]
    #[case::cis_trans_lifts_term(
        BondConstraint::CisTransStereo(CisTransStereoAst::Stereo(StereoCosetAst::term(StereoTerm::Lit(1)))),
        Ok(BondConstraint::cis_trans_stereo(CisTransStereoAst::stereo(1_u32))))]
    #[case::empty_litset_contradiction(
        BondConstraint::RingMembership(RingMembershipAst::new(RingScope::All, ValueAst::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_bond_constraint_canonicalize(
        #[case] constraint: BondConstraint,
        #[case] expected: Result<BondConstraint, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
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
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintKind::Aromatic, Some(BondConstraint::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_membership(BondConstraintKind::RingMembership, Some(BondConstraint::ring_membership(RingScope::Size(6), 1)))]
    #[case::cis_trans_absent(BondConstraintKind::CisTransStereo, None)]
    fn test_bond_constraints_get(
        #[case] kind: BondConstraintKind,
        #[case] expected: Option<BondConstraint>,
    ) {
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(kind), expected.as_ref());
    }

    #[rstest]
    fn test_bond_constraints_get_mut() {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let entry = cs.get_mut(BondConstraintKind::RingMembership).unwrap();
        *entry = BondConstraint::ring_membership(RingScope::Size(5), 1);
        assert_eq!(
            cs.as_slice(),
            &[
                BondConstraint::Aromatic(BooleanAst::Lit(true)),
                BondConstraint::ring_membership(RingScope::Size(5), 1),
            ],
        );
    }

    #[rstest]
    fn test_bond_constraints_get_mut_absent() {
        let mut cs = BondConstraints::from_iter([BondConstraint::Aromatic(BooleanAst::Lit(true))]);
        assert!(cs.get_mut(BondConstraintKind::RingMembership).is_none());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(BondConstraintKey::Aromatic, true)]
    #[case::ring_all_present(BondConstraintKey::RingMembership(RingScope::All), true)]
    #[case::ring_size_present(BondConstraintKey::RingMembership(RingScope::Size(6)), true)]
    #[case::ring_size_absent(BondConstraintKey::RingMembership(RingScope::Size(5)), false)]
    #[case::cis_trans_absent(BondConstraintKey::CisTransStereo, false)]
    fn test_bond_constraints_contains_key(
        #[case] key: BondConstraintKey,
        #[case] expected: bool,
    ) {
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::All, 2),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains_key(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintKey::Aromatic, Some(BondConstraint::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_all(BondConstraintKey::RingMembership(RingScope::All), Some(BondConstraint::ring_membership(RingScope::All, 2)))]
    #[case::ring_size(BondConstraintKey::RingMembership(RingScope::Size(6)), Some(BondConstraint::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_size_absent(BondConstraintKey::RingMembership(RingScope::Size(5)), None)]
    fn test_bond_constraints_get_by_key(
        #[case] key: BondConstraintKey,
        #[case] expected: Option<BondConstraint>,
    ) {
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::All, 2),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get_by_key(key), expected.as_ref());
    }

    #[rstest]
    fn test_bond_constraints_get_by_key_mut() {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::ring_membership(RingScope::All, 2),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let slot = cs
            .get_by_key_mut(BondConstraintKey::RingMembership(RingScope::Size(6)))
            .unwrap();
        *slot = BondConstraint::ring_membership(RingScope::Size(6), 2);
        assert_eq!(
            cs.get_by_key(BondConstraintKey::RingMembership(RingScope::Size(6))),
            Some(&BondConstraint::ring_membership(RingScope::Size(6), 2)),
        );
    }

    #[rstest]
    fn test_bond_constraints_remove_by_key() {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::All, 2),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let removed = cs.remove_by_key(BondConstraintKey::RingMembership(RingScope::Size(6)));
        assert_eq!(
            removed,
            Some(BondConstraint::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(
            cs.as_slice(),
            &[
                BondConstraint::Aromatic(BooleanAst::Lit(true)),
                BondConstraint::ring_membership(RingScope::All, 2),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::merge_same_scope(
        BondConstraints::from_iter([
            BondConstraint::ring_membership(RingScope::All, ValueAst::lit_set([1, 2])),
            BondConstraint::ring_membership(RingScope::All, ValueAst::lit_set([2, 3])),
        ]),
        Ok(BondConstraints::from_iter([BondConstraint::ring_membership(RingScope::All, 2)])))]
    #[case::drop_vacuous(
        BondConstraints::from_iter([
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::All, ValueAst::Undetermined),
        ]),
        Ok(BondConstraints::from_iter([BondConstraint::Aromatic(BooleanAst::Lit(true))])))]
    #[case::canonicalizes_values(
        BondConstraints::from_iter([
            BondConstraint::CisTransStereo(CisTransStereoAst::Stereo(StereoCosetAst::term(StereoTerm::Lit(1)))),
        ]),
        Ok(BondConstraints::from_iter([BondConstraint::cis_trans_stereo(CisTransStereoAst::stereo(1_u32))])))]
    #[case::contradiction_same_scope(
        BondConstraints::from_iter([
            BondConstraint::ring_membership(RingScope::All, 1),
            BondConstraint::ring_membership(RingScope::All, 0),
        ]),
        Err(Contradiction))]
    fn test_bond_constraints_canonicalize(
        #[case] constraints: BondConstraints,
        #[case] expected: Result<BondConstraints, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rstest]
    fn test_bond_constraints_iter() {
        let cs = BondConstraints::from_iter([
            BondConstraint::ring_membership(RingScope::Size(6), 1),
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::All, 1),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                BondConstraint::Aromatic(BooleanAst::Lit(true)),
                BondConstraint::ring_membership(RingScope::All, 1),
                BondConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![BondConstraint::Aromatic(BooleanAst::Lit(true))], vec![None], vec![BondConstraint::Aromatic(BooleanAst::Lit(true))])]
    #[case::replace_value_kind(vec![BondConstraint::cis_trans_stereo(CisTransStereoAst::Undetermined), BondConstraint::cis_trans_stereo(CisTransStereoAst::NotStereo)],
        vec![None, Some(BondConstraint::cis_trans_stereo(CisTransStereoAst::Undetermined))], vec![BondConstraint::cis_trans_stereo(CisTransStereoAst::NotStereo)])]
    #[case::replace_unit_variant(vec![BondConstraint::Aromatic(BooleanAst::Lit(true)), BondConstraint::Aromatic(BooleanAst::Lit(true))],
        vec![None, Some(BondConstraint::Aromatic(BooleanAst::Lit(true)))], vec![BondConstraint::Aromatic(BooleanAst::Lit(true))])]
    #[case::append_multi_valued_ring(vec![BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::All, 2)],
        vec![None, None], vec![BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::All, 2)])]
    #[case::distinct_kinds(vec![BondConstraint::Aromatic(BooleanAst::Lit(true)), BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::Size(6), 1)],
        vec![None, None, None], vec![BondConstraint::Aromatic(BooleanAst::Lit(true)), BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::Size(6), 1)])]
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
    #[case::partial(|c: &BondConstraint| matches!(c, BondConstraint::Aromatic(BooleanAst::Lit(true))) || matches!(c, BondConstraint::RingMembership(m) if m.scope == RingScope::Size(6)), vec![
            BondConstraint::Aromatic(BooleanAst::Lit(true)), BondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::all_dropped(|_: &BondConstraint| false, vec![])]
    fn test_bond_constraints_retain(
        #[case] predicate: impl FnMut(&BondConstraint) -> bool,
        #[case] expected: Vec<BondConstraint>,
    ) {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::All, 1),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        cs.retain(predicate);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_bond_constraints_clear() {
        let mut cs = BondConstraints::from_iter([BondConstraint::Aromatic(BooleanAst::Lit(true))]);
        cs.clear();
        assert_eq!(cs, BondConstraints::new());
    }

    #[rstest]
    fn test_bond_constraints_take() {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(
            drained,
            vec![
                BondConstraint::Aromatic(BooleanAst::Lit(true)),
                BondConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
        assert_eq!(cs, BondConstraints::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(BondConstraintKind::Aromatic, Some(BondConstraint::Aromatic(BooleanAst::Lit(true))), vec![BondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::ring_membership_present(BondConstraintKind::RingMembership, Some(BondConstraint::ring_membership(RingScope::Size(6), 1)), vec![BondConstraint::Aromatic(BooleanAst::Lit(true))])]
    #[case::cis_trans_absent(BondConstraintKind::CisTransStereo, None, vec![BondConstraint::Aromatic(BooleanAst::Lit(true)), BondConstraint::ring_membership(RingScope::Size(6), 1)])]
    fn test_bond_constraints_remove(
        #[case] kind: BondConstraintKind,
        #[case] expected_returned: Option<BondConstraint>,
        #[case] expected_state: Vec<BondConstraint>,
    ) {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.remove(kind), expected_returned);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rstest]
    fn test_bond_constraints_compact() {
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let compaction = IdCompaction::new(
            Compaction::new(vec![0, 1, 2], vec![0, 1]),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(cs.clone().compact(&compaction), cs);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(BondConstraints::new(), BondConstraints::from_iter([BondConstraint::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::aromatic_required_present(BondConstraints::from_iter([BondConstraint::Aromatic(BooleanAst::Lit(true))]),
        BondConstraints::from_iter([BondConstraint::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::aromatic_required_absent(BondConstraints::from_iter([BondConstraint::Aromatic(BooleanAst::Lit(true))]), BondConstraints::new(), false)]
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
    #[case::distinct(vec![BondConstraint::Aromatic(BooleanAst::Lit(true)), BondConstraint::ring_membership(RingScope::All, 1)],
        vec![BondConstraint::Aromatic(BooleanAst::Lit(true)), BondConstraint::ring_membership(RingScope::All, 1)])]
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
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                BondConstraint::Aromatic(BooleanAst::Lit(true)),
                BondConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rstest]
    fn test_bond_constraints_from_bond_constraint() {
        let cs: BondConstraints = BondConstraint::Aromatic(BooleanAst::Lit(true)).into();
        assert_eq!(
            cs.as_slice(),
            &[BondConstraint::Aromatic(BooleanAst::Lit(true))]
        );
    }

    #[rstest]
    fn test_bond_constraints_from_vec() {
        let cs: BondConstraints = vec![
            BondConstraint::Aromatic(BooleanAst::Lit(true)),
            BondConstraint::ring_membership(RingScope::Size(6), 1),
        ]
        .into();
        assert_eq!(
            cs.as_slice(),
            &[
                BondConstraint::Aromatic(BooleanAst::Lit(true)),
                BondConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }
}
