//! Localized bond constraints.
use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use super::super::boolean::BooleanAst;
use super::super::constraint::ring::{RingMembershipAst, RingScope};
use super::super::error::{Contradiction, NoJoin};
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::stereo::CisTransStereoAst;
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Localized bond constraint.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondConstraintAst {
    Aromatic(BooleanAst),
    CisTransStereo(CisTransStereoAst),
    /// Ring count in the fixed Relevant ring projection, optionally restricted by size.
    RingMembership(RingMembershipAst),
}

impl BondConstraintAst {
    pub fn aromatic(b: impl Into<BooleanAst>) -> Self {
        Self::Aromatic(b.into())
    }

    pub fn cis_trans_stereo(c: impl Into<CisTransStereoAst>) -> Self {
        Self::CisTransStereo(c.into())
    }

    pub fn ring_membership(scope: RingScope, count: impl Into<ValueAst>) -> Self {
        Self::RingMembership(RingMembershipAst::new(scope, count.into()))
    }

    /// Bond constraint key, unique within a `BondConstraintsAst` container.
    pub fn key(&self) -> BondConstraintKey {
        match self {
            Self::Aromatic(_) => BondConstraintKey::Aromatic,
            Self::CisTransStereo(_) => BondConstraintKey::CisTransStereo,
            Self::RingMembership(m) => BondConstraintKey::RingMembership(m.scope),
        }
    }

    /// Vacuous form of constraint key, used for removal.
    pub fn as_undetermined(&self) -> Self {
        match self {
            Self::Aromatic(_) => Self::Aromatic(BooleanAst::Undetermined),
            Self::CisTransStereo(_) => Self::CisTransStereo(CisTransStereoAst::Undetermined),
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipAst::new(m.scope, ValueAst::Undetermined))
            }
        }
    }

    /// Value-only payload: no entity ids to compact, so this never drops.
    pub fn compact(self, _compaction: &IdCompaction) -> Option<Self> {
        Some(self)
    }

    /// Value-only payload: no entity ids to remap.
    pub(crate) fn remap(self, _map: &IdRemapping) -> Self {
        self
    }
}

impl Canonicalize for BondConstraintAst {
    /// Canonicalize the inner value; kind and sub-key are preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Aromatic(b) => Self::Aromatic(b.canonicalize()?),
            Self::CisTransStereo(c) => Self::CisTransStereo(c.canonicalize()?),
            Self::RingMembership(m) => Self::RingMembership(m.canonicalize()?),
        })
    }
}

impl Lattice for BondConstraintAst {
    fn is_undetermined(&self) -> bool {
        match self {
            Self::Aromatic(b) => b.is_undetermined(),
            Self::CisTransStereo(c) => c.is_undetermined(),
            Self::RingMembership(m) => m.is_undetermined(),
        }
    }

    fn is_ground(&self) -> bool {
        match self {
            Self::Aromatic(b) => b.is_ground(),
            Self::CisTransStereo(c) => c.is_ground(),
            Self::RingMembership(m) => m.is_ground(),
        }
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Aromatic(a), Self::Aromatic(b)) => a.meet(b).map(Self::Aromatic),
            (Self::CisTransStereo(a), Self::CisTransStereo(b)) => {
                a.meet(b).map(Self::CisTransStereo)
            }
            (Self::RingMembership(a), Self::RingMembership(b)) => {
                a.meet(b).map(Self::RingMembership)
            }
            _ => None,
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        match (self, other) {
            (Self::Aromatic(a), Self::Aromatic(b)) => Ok(Self::Aromatic(a.join(b)?)),
            (Self::CisTransStereo(a), Self::CisTransStereo(b)) => {
                Ok(Self::CisTransStereo(a.join(b)?))
            }
            (Self::RingMembership(a), Self::RingMembership(b)) => {
                a.join(b).map(Self::RingMembership)
            }
            _ => Err(NoJoin),
        }
    }

    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Aromatic(a), Self::Aromatic(b)) => a.matches(b),
            (Self::CisTransStereo(a), Self::CisTransStereo(b)) => a.matches(b),
            (Self::RingMembership(a), Self::RingMembership(b)) => a.matches(b),
            _ => false,
        }
    }

    fn is_compatible(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Aromatic(a), Self::Aromatic(b)) => a.is_compatible(b),
            (Self::CisTransStereo(a), Self::CisTransStereo(b)) => a.is_compatible(b),
            (Self::RingMembership(a), Self::RingMembership(b)) => a.is_compatible(b),
            _ => false,
        }
    }
}

/// Entry identity: discriminant + sub-key, BondConstraintsAst is ordered, unique by key.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondConstraintKey {
    Aromatic,
    CisTransStereo,
    RingMembership(RingScope),
}

/// Atom constraints container, ordered, unique by key, sorted flat vector storage.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondConstraintsAst(Vec<BondConstraintAst>);

impl BondConstraintsAst {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The bond's aromatic value, or `Undetermined` when no `Aromatic` constraint is present.
    pub fn aromatic(&self) -> BooleanAst {
        match self.get(BondConstraintKey::Aromatic) {
            Some(BondConstraintAst::Aromatic(b)) => *b,
            _ => BooleanAst::Undetermined,
        }
    }

    pub fn cis_trans_stereo(&self) -> Option<&CisTransStereoAst> {
        match self.get(BondConstraintKey::CisTransStereo) {
            Some(BondConstraintAst::CisTransStereo(c)) => Some(c),
            _ => None,
        }
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &ValueAst)> {
        self.iter().filter_map(|c| match c {
            BondConstraintAst::RingMembership(m) => Some((m.scope, &m.count)),
            _ => None,
        })
    }

    fn ring_membership(&self, scope: RingScope) -> Option<&ValueAst> {
        self.ring_memberships()
            .find(|(s, _)| *s == scope)
            .map(|(_, v)| v)
    }

    pub fn ring_count(&self) -> Option<&ValueAst> {
        self.ring_membership(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> Option<&ValueAst> {
        self.ring_membership(RingScope::Size(s))
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    fn find(&self, key: BondConstraintKey) -> Result<usize, usize> {
        self.0.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains(&self, key: BondConstraintKey) -> bool {
        self.find(key).is_ok()
    }

    pub fn get(&self, key: BondConstraintKey) -> Option<&BondConstraintAst> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: BondConstraintAst) {
        match self.find(c.key()) {
            Ok(i) => self.0[i] = c,
            Err(i) => self.0.insert(i, c),
        }
    }

    /// Transactional write at one key: verify the current value `canonical_eq` `old` (both absent
    /// matches), then apply `new` (`Some` sets, `None` removes). `old`/`new` address the same key.
    /// `Err` on a key or old-value mismatch; the store is unchanged when it errors. The delta
    /// apply/undo primitive.
    pub fn compare_and_set(
        &mut self,
        old: Option<BondConstraintAst>,
        new: Option<BondConstraintAst>,
    ) -> Result<(), Contradiction> {
        let key = match (&old, &new) {
            (Some(o), Some(n)) => {
                if o.key() != n.key() {
                    return Err(Contradiction);
                }
                o.key()
            }
            (Some(o), None) => o.key(),
            (None, Some(n)) => n.key(),
            (None, None) => return Ok(()),
        };
        let matches = match (self.get(key), old.as_ref()) {
            (None, None) => true,
            (Some(current), Some(old)) => current.canonical_eq(old),
            _ => false,
        };
        if !matches {
            return Err(Contradiction);
        }
        match new {
            Some(c) => self.set(c),
            None => {
                self.remove(key);
            }
        }
        Ok(())
    }

    pub fn remove(&mut self, key: BondConstraintKey) -> Option<BondConstraintAst> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = BondConstraintAst>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &BondConstraintsAst) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&BondConstraintAst) -> bool) {
        self.0.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = BondConstraintAst> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, BondConstraintAst> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &IdCompaction) -> Self {
        self
    }
}

impl Canonicalize for BondConstraintsAst {
    /// Canonicalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Canonicalize::canonicalize)
            .collect::<Result<Vec<BondConstraintAst>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl Lattice for BondConstraintsAst {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`BondConstraintAst::meet`; a `None` aborts the whole meet), an A-only /
    /// B-only key is kept (meet with the absent ⊤ is the value). Vacuous results are dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<BondConstraintAst> = Vec::new();
        let mut a = self.0.iter();
        let mut b = other.0.iter();
        let mut ca = a.next();
        let mut cb = b.next();
        loop {
            let (met, adv_a, adv_b) = match (ca, cb) {
                (Some(x), Some(y)) => match x.key().cmp(&y.key()) {
                    Ordering::Less => (x.clone(), true, false),
                    Ordering::Greater => (y.clone(), false, true),
                    Ordering::Equal => (x.meet(y)?, true, true),
                },
                (Some(x), None) => (x.clone(), true, false),
                (None, Some(y)) => (y.clone(), false, true),
                (None, None) => break,
            };
            if !met.is_undetermined() {
                entries.push(met);
            }
            if adv_a {
                ca = a.next();
            }
            if adv_b {
                cb = b.next();
            }
        }
        Some(Self(entries))
    }

    /// Least upper bound as a two-pointer merge: only keys present on *both* sides join
    /// (`BondConstraintAst::join`); a single-side key widens to the absent ⊤ and is dropped. The
    /// container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<BondConstraintAst> = Vec::new();
        let mut a = self.0.iter();
        let mut b = other.0.iter();
        let mut ca = a.next();
        let mut cb = b.next();
        while let (Some(x), Some(y)) = (ca, cb) {
            match x.key().cmp(&y.key()) {
                Ordering::Less => ca = a.next(),
                Ordering::Greater => cb = b.next(),
                Ordering::Equal => {
                    if let Ok(j) = x.join(y) {
                        if !j.is_undetermined() {
                            entries.push(j);
                        }
                    }
                    ca = a.next();
                    cb = b.next();
                }
            }
        }
        Ok(Self(entries))
    }

    /// Pattern-driven: every constraint the pattern carries must match the target,
    /// looked up by reference. Each value is matched on its own lattice; an empty
    /// pattern matches any target.
    fn matches(&self, target: &Self) -> bool {
        self.iter().all(|c| match c {
            BondConstraintAst::Aromatic(b) => b.matches(&target.aromatic()),
            BondConstraintAst::CisTransStereo(cts) => cts.matches(
                target
                    .cis_trans_stereo()
                    .unwrap_or(&CisTransStereoAst::Undetermined),
            ),
            BondConstraintAst::RingMembership(rm) => rm.count.matches(
                target
                    .ring_membership(rm.scope)
                    .unwrap_or(&ValueAst::Undetermined),
            ),
        })
    }

    /// Sorted merge, short-circuit: only shared keys can conflict; non-shared keys are always
    /// compatible. Cheaper than the `meet`-derived default — builds no result container.
    fn is_compatible(&self, other: &Self) -> bool {
        let mut a = self.0.iter();
        let mut b = other.0.iter();
        let mut ca = a.next();
        let mut cb = b.next();
        while let (Some(x), Some(y)) = (ca, cb) {
            match x.key().cmp(&y.key()) {
                Ordering::Less => ca = a.next(),
                Ordering::Greater => cb = b.next(),
                Ordering::Equal => {
                    if !x.is_compatible(y) {
                        return false;
                    }
                    ca = a.next();
                    cb = b.next();
                }
            }
        }
        true
    }
}

impl FromIterator<BondConstraintAst> for BondConstraintsAst {
    fn from_iter<I: IntoIterator<Item = BondConstraintAst>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for BondConstraintsAst {
    type Item = BondConstraintAst;
    type IntoIter = IntoIter<BondConstraintAst>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<BondConstraintAst> for BondConstraintsAst {
    fn from(c: BondConstraintAst) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<BondConstraintAst>> for BondConstraintsAst {
    fn from(cs: Vec<BondConstraintAst>) -> Self {
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
    #[case::ring_membership_all(BondConstraintAst::ring_membership(RingScope::All, 1), BondConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1)))]
    #[case::ring_membership_size(BondConstraintAst::ring_membership(RingScope::Size(6), 1), BondConstraintAst::ring_membership(RingScope::Size(6), 1))]
    #[case::cis_trans_stereo(BondConstraintAst::cis_trans_stereo(CisTransStereoAst::NotStereo), BondConstraintAst::CisTransStereo(CisTransStereoAst::NotStereo))]
    fn test_bond_constraint_ast_constructors(
        #[case] actual: BondConstraintAst,
        #[case] expected: BondConstraintAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintKey::Aromatic)]
    #[case::ring_membership_all(BondConstraintAst::ring_membership(RingScope::All, 1), BondConstraintKey::RingMembership(RingScope::All))]
    #[case::ring_membership_size(BondConstraintAst::ring_membership(RingScope::Size(6), 1), BondConstraintKey::RingMembership(RingScope::Size(6)))]
    #[case::cis_trans_stereo(BondConstraintAst::CisTransStereo(CisTransStereoAst::NotStereo), BondConstraintKey::CisTransStereo)]
    fn test_bond_constraint_ast_key(#[case] c: BondConstraintAst, #[case] expected: BondConstraintKey) {
        assert_eq!(c.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::Aromatic(BooleanAst::Undetermined))]
    #[case::ring_membership_keeps_scope(BondConstraintAst::ring_membership(RingScope::Size(6), 1), BondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined))]
    #[case::cis_trans(BondConstraintAst::CisTransStereo(CisTransStereoAst::stereo(1_u32)), BondConstraintAst::CisTransStereo(CisTransStereoAst::Undetermined))]
    fn test_bond_constraint_ast_as_undetermined(#[case] c: BondConstraintAst, #[case] expected: BondConstraintAst) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), Ok(BondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_count_litset_singleton(
        BondConstraintAst::RingMembership(RingMembershipAst::new(RingScope::All, ValueAst::lit_set([2]))),
        Ok(BondConstraintAst::ring_membership(RingScope::All, 2)))]
    #[case::cis_trans_lifts_term(
        BondConstraintAst::CisTransStereo(CisTransStereoAst::Stereo(StereoCosetAst::term(StereoTerm::Lit(1)))),
        Ok(BondConstraintAst::cis_trans_stereo(CisTransStereoAst::stereo(1_u32))))]
    #[case::empty_litset_contradiction(
        BondConstraintAst::RingMembership(RingMembershipAst::new(RingScope::All, ValueAst::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_bond_constraint_ast_canonicalize(
        #[case] constraint: BondConstraintAst,
        #[case] expected: Result<BondConstraintAst, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), false)]
    #[case::ring_membership_all_lit(BondConstraintAst::ring_membership(RingScope::All, 1), false)]
    #[case::ring_membership_all_undetermined(BondConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined), true)]
    #[case::ring_membership_size_lit(BondConstraintAst::ring_membership(RingScope::Size(6), 1), false)]
    #[case::ring_membership_size_undetermined(BondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined), true)]
    #[case::cis_trans_not_stereo(BondConstraintAst::CisTransStereo(CisTransStereoAst::NotStereo), false)]
    #[case::cis_trans_undetermined(BondConstraintAst::CisTransStereo(CisTransStereoAst::Undetermined), true)]
    fn test_bond_constraint_ast_is_undetermined(#[case] c: BondConstraintAst, #[case] expected: bool) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::Aromatic(BooleanAst::Undetermined), Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    #[case::same_key_incompatible(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::Aromatic(BooleanAst::Lit(false)), None)]
    #[case::different_key(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::All, 1), None)]
    fn test_bond_constraint_ast_meet(#[case] a: BondConstraintAst, #[case] b: BondConstraintAst, #[case] expected: Option<BondConstraintAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_widens(BondConstraintAst::ring_membership(RingScope::All, 1), BondConstraintAst::ring_membership(RingScope::All, 2), Ok(BondConstraintAst::ring_membership(RingScope::All, ValueAst::lit_set([1, 2]))))]
    #[case::different_key(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::All, 1), Err(NoJoin))]
    fn test_bond_constraint_ast_join(#[case] a: BondConstraintAst, #[case] b: BondConstraintAst, #[case] expected: Result<BondConstraintAst, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::Aromatic(BooleanAst::Lit(true)), true)]
    #[case::same_key_incompatible(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::Aromatic(BooleanAst::Lit(false)), false)]
    #[case::different_key(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::All, 1), false)]
    fn test_bond_constraint_ast_is_compatible(#[case] a: BondConstraintAst, #[case] b: BondConstraintAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_bond_constraints_ast_new() {
        let cs = BondConstraintsAst::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    fn test_bond_constraints_ast_iter() {
        let cs = BondConstraintsAst::from_iter([
            BondConstraintAst::ring_membership(RingScope::Size(6), 1),
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            BondConstraintAst::ring_membership(RingScope::All, 1),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                BondConstraintAst::ring_membership(RingScope::All, 1),
                BondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true))], vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true))])]
    #[case::overwrite_same_key(vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::Aromatic(BooleanAst::Lit(false))], vec![BondConstraintAst::Aromatic(BooleanAst::Lit(false))])]
    #[case::vacuous_stores(vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::Aromatic(BooleanAst::Undetermined)], vec![BondConstraintAst::Aromatic(BooleanAst::Undetermined)])]
    #[case::new_key_sorts(vec![BondConstraintAst::ring_membership(RingScope::Size(6), 1), BondConstraintAst::Aromatic(BooleanAst::Lit(true))], vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::Size(6), 1)])]
    #[case::ring_overwrite_scope(vec![BondConstraintAst::ring_membership(RingScope::Size(6), 1), BondConstraintAst::ring_membership(RingScope::Size(6), 2)], vec![BondConstraintAst::ring_membership(RingScope::Size(6), 2)])]
    fn test_bond_constraints_ast_set(#[case] sequence: Vec<BondConstraintAst>, #[case] expected: Vec<BondConstraintAst>) {
        let mut cs = BondConstraintsAst::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, BondConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite_shared(
        vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::All, 1)],
        vec![BondConstraintAst::Aromatic(BooleanAst::Lit(false))],
        vec![BondConstraintAst::Aromatic(BooleanAst::Lit(false)), BondConstraintAst::ring_membership(RingScope::All, 1)])]
    #[case::keeps_disjoint(
        vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true))],
        vec![BondConstraintAst::ring_membership(RingScope::All, 1)],
        vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::All, 1)])]
    #[case::vacuous_removes(
        vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::All, 1)],
        vec![BondConstraintAst::Aromatic(BooleanAst::Undetermined)],
        vec![BondConstraintAst::ring_membership(RingScope::All, 1)])]
    fn test_bond_constraints_ast_update(#[case] initial: Vec<BondConstraintAst>, #[case] other: Vec<BondConstraintAst>, #[case] expected: Vec<BondConstraintAst>) {
        let mut cs = BondConstraintsAst::from_iter(initial);
        cs.update(&BondConstraintsAst::from_iter(other));
        assert_eq!(cs, BondConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![BondConstraintAst::ring_membership(RingScope::All, 1)], Some(BondConstraintAst::ring_membership(RingScope::All, 1)), Some(BondConstraintAst::ring_membership(RingScope::All, 2)), Ok(()), vec![BondConstraintAst::ring_membership(RingScope::All, 2)])]
    #[case::remove(vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true))], Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))), Ok(()), vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true))])]
    #[case::old_mismatch(vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true))], Some(BondConstraintAst::Aromatic(BooleanAst::Lit(false))), None, Err(Contradiction), vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true))])]
    #[case::key_mismatch(vec![], Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))), Some(BondConstraintAst::ring_membership(RingScope::All, 1)), Err(Contradiction), vec![])]
    fn test_bond_constraints_ast_compare_and_set(
        #[case] initial: Vec<BondConstraintAst>,
        #[case] old: Option<BondConstraintAst>,
        #[case] new: Option<BondConstraintAst>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<BondConstraintAst>,
    ) {
        let mut cs = BondConstraintsAst::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, BondConstraintsAst::from_iter(expected_state));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(BondConstraintKey::Aromatic, true)]
    #[case::ring_all_present(BondConstraintKey::RingMembership(RingScope::All), true)]
    #[case::ring_size_present(BondConstraintKey::RingMembership(RingScope::Size(6)), true)]
    #[case::ring_size_absent(BondConstraintKey::RingMembership(RingScope::Size(5)), false)]
    #[case::cis_trans_absent(BondConstraintKey::CisTransStereo, false)]
    fn test_bond_constraints_ast_contains(
        #[case] key: BondConstraintKey,
        #[case] expected: bool,
    ) {
        let cs = BondConstraintsAst::from_iter([
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            BondConstraintAst::ring_membership(RingScope::All, 2),
            BondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintKey::Aromatic, Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_all(BondConstraintKey::RingMembership(RingScope::All), Some(BondConstraintAst::ring_membership(RingScope::All, 2)))]
    #[case::ring_size(BondConstraintKey::RingMembership(RingScope::Size(6)), Some(BondConstraintAst::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_size_absent(BondConstraintKey::RingMembership(RingScope::Size(5)), None)]
    fn test_bond_constraints_ast_get(
        #[case] key: BondConstraintKey,
        #[case] expected: Option<BondConstraintAst>,
    ) {
        let cs = BondConstraintsAst::from_iter([
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            BondConstraintAst::ring_membership(RingScope::All, 2),
            BondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(key), expected.as_ref());
    }

    #[rstest]
    fn test_bond_constraints_ast_remove() {
        let mut cs = BondConstraintsAst::from_iter([
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            BondConstraintAst::ring_membership(RingScope::All, 2),
            BondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        let removed = cs.remove(BondConstraintKey::RingMembership(RingScope::Size(6)));
        assert_eq!(
            removed,
            Some(BondConstraintAst::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(
            cs,
            BondConstraintsAst::from_iter([
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                BondConstraintAst::ring_membership(RingScope::All, 2),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &BondConstraintAst| matches!(c, BondConstraintAst::Aromatic(BooleanAst::Lit(true))) || matches!(c, BondConstraintAst::RingMembership(m) if m.scope == RingScope::Size(6)), vec![
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::Size(6), 1)])]
    #[case::all_dropped(|_: &BondConstraintAst| false, vec![])]
    fn test_bond_constraints_ast_retain(
        #[case] predicate: impl FnMut(&BondConstraintAst) -> bool,
        #[case] expected: Vec<BondConstraintAst>,
    ) {
        let mut cs = BondConstraintsAst::from_iter([
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            BondConstraintAst::ring_membership(RingScope::All, 1),
            BondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        cs.retain(predicate);
        assert_eq!(cs, BondConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_bond_constraints_ast_clear() {
        let mut cs =
            BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]);
        cs.clear();
        assert_eq!(cs, BondConstraintsAst::new());
    }

    #[rstest]
    fn test_bond_constraints_ast_take() {
        let mut cs = BondConstraintsAst::from_iter([
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            BondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(
            drained,
            vec![
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                BondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ],
        );
        assert_eq!(cs, BondConstraintsAst::new());
    }

    #[rstest]
    fn test_bond_constraints_ast_compact() {
        let cs = BondConstraintsAst::from_iter([
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            BondConstraintAst::ring_membership(RingScope::Size(6), 1),
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
    #[case::drop_vacuous(
        BondConstraintsAst::from_iter([
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            BondConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined),
        ]),
        Ok(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))])))]
    #[case::canonicalizes_values(
        BondConstraintsAst::from_iter([
            BondConstraintAst::CisTransStereo(CisTransStereoAst::Stereo(StereoCosetAst::term(StereoTerm::Lit(1)))),
        ]),
        Ok(BondConstraintsAst::from_iter([BondConstraintAst::cis_trans_stereo(CisTransStereoAst::stereo(1_u32))])))]
    fn test_bond_constraints_ast_canonicalize(
        #[case] constraints: BondConstraintsAst,
        #[case] expected: Result<BondConstraintsAst, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys_kept(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, 1)]),
        Some(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::All, 1)])))]
    #[case::shared_key_meets(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Undetermined)]),
        Some(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))])))]
    #[case::shared_key_contradicts(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(false))]), None)]
    #[case::ring_size_unions(BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::Size(5), 1)]), BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::Size(6), 1)]),
        Some(BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::Size(5), 1), BondConstraintAst::ring_membership(RingScope::Size(6), 1)])))]
    #[case::prunes_vacuous(BondConstraintsAst::new(), BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Undetermined)]), Some(BondConstraintsAst::new()))]
    fn test_bond_constraints_ast_meet(#[case] a: BondConstraintsAst, #[case] b: BondConstraintsAst, #[case] expected: Option<BondConstraintsAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::keeps_only_shared_keys(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::All, 1)]), BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]),
        BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]))]
    #[case::widens_value(BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, 1)]), BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, 2)]),
        BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, ValueAst::lit_set([1, 2]))]))]
    #[case::incompatible_drops_to_undetermined(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(false))]), BondConstraintsAst::new())]
    fn test_bond_constraints_ast_join(#[case] a: BondConstraintsAst, #[case] b: BondConstraintsAst, #[case] expected: BondConstraintsAst) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(BondConstraintsAst::new(), BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::aromatic_required_present(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]),
        BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::aromatic_required_absent(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), BondConstraintsAst::new(), false)]
    #[case::ring_membership_all_wildcard_matches_lit(BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined)]),
        BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, 1)]), true)]
    #[case::ring_membership_all_lit_mismatch(BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, 1)]),
        BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, 2)]), false)]
    #[case::ring_membership_size_subset(BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::Size(5), 1)]),
        BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::Size(5), 1), BondConstraintAst::ring_membership(RingScope::Size(6), 1)]), true)]
    #[case::ring_membership_size_not_in_target(BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::Size(7), 1)]),
        BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::Size(5), 1)]), false)]
    #[case::cis_trans_match(BondConstraintsAst::from_iter([BondConstraintAst::CisTransStereo(CisTransStereoAst::NotStereo)]),
        BondConstraintsAst::from_iter([BondConstraintAst::CisTransStereo(CisTransStereoAst::NotStereo)]), true)]
    #[case::cis_trans_pattern_more_specific(BondConstraintsAst::from_iter([BondConstraintAst::CisTransStereo(CisTransStereoAst::NotStereo)]),
        BondConstraintsAst::new(), false)]
    fn test_bond_constraints_ast_matches(
        #[case] pattern: BondConstraintsAst,
        #[case] target: BondConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, 1)]), true)]
    #[case::shared_key_compatible(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::shared_key_incompatible(BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]), BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(false))]), false)]
    fn test_bond_constraints_ast_is_compatible(#[case] a: BondConstraintsAst, #[case] b: BondConstraintsAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::All, 1)],
        vec![BondConstraintAst::Aromatic(BooleanAst::Lit(true)), BondConstraintAst::ring_membership(RingScope::All, 1)])]
    #[case::unique_kind_last_wins(vec![BondConstraintAst::cis_trans_stereo(CisTransStereoAst::Undetermined), BondConstraintAst::cis_trans_stereo(CisTransStereoAst::NotStereo)],
        vec![BondConstraintAst::cis_trans_stereo(CisTransStereoAst::NotStereo)])]
    #[case::ring_appends(vec![BondConstraintAst::ring_membership(RingScope::All, 1), BondConstraintAst::ring_membership(RingScope::Size(6), 1)],
        vec![BondConstraintAst::ring_membership(RingScope::All, 1), BondConstraintAst::ring_membership(RingScope::Size(6), 1)])]
    #[case::empty(vec![], vec![])]
    fn test_bond_constraints_ast_from_iter(
        #[case] input: Vec<BondConstraintAst>,
        #[case] expected: Vec<BondConstraintAst>,
    ) {
        let cs = BondConstraintsAst::from_iter(input);
        assert_eq!(cs, BondConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_bond_constraints_ast_into_iter() {
        let cs = BondConstraintsAst::from_iter([
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            BondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                BondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rstest]
    fn test_bond_constraints_ast_from_bond_constraint() {
        let cs: BondConstraintsAst = BondConstraintAst::Aromatic(BooleanAst::Lit(true)).into();
        assert_eq!(
            cs,
            BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))]),
        );
    }

    #[rstest]
    fn test_bond_constraints_ast_from_vec() {
        let cs: BondConstraintsAst = vec![
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            BondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]
        .into();
        assert_eq!(
            cs,
            BondConstraintsAst::from_iter([
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                BondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ]),
        );
    }
}
