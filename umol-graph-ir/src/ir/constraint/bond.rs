//! Localized bond constraints.
use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use super::super::boolean::BooleanForm;
use super::super::constraint::ring::{RingMembershipForm, RingScope};
use super::super::error::{Contradiction, NoJoin};
use super::super::num::NumForm;
use super::super::remap::{IdRemapping, MoleculeCompaction};
use super::super::stereo::CisTransStereoForm;
use super::super::traits::{Lattice, Normalize};

/// Localized bond constraint.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondConstraintForm {
    Aromatic(BooleanForm),
    CisTransStereo(CisTransStereoForm),
    /// Ring count in the fixed Relevant ring projection, optionally restricted by size.
    RingMembership(RingMembershipForm),
}

impl BondConstraintForm {
    pub fn aromatic(b: impl Into<BooleanForm>) -> Self {
        Self::Aromatic(b.into())
    }

    pub fn cis_trans_stereo(c: impl Into<CisTransStereoForm>) -> Self {
        Self::CisTransStereo(c.into())
    }

    pub fn ring_membership(scope: RingScope, count: impl Into<NumForm>) -> Self {
        Self::RingMembership(RingMembershipForm::new(scope, count.into()))
    }

    /// Bond constraint key, unique within a `BondConstraintsForm` container.
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
            Self::Aromatic(_) => Self::Aromatic(BooleanForm::Undetermined),
            Self::CisTransStereo(_) => Self::CisTransStereo(CisTransStereoForm::Undetermined),
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipForm::new(m.scope, NumForm::Undetermined))
            }
        }
    }

    /// Value-only payload: no entity ids to compact, so this never drops.
    pub fn compact(self, _compaction: &MoleculeCompaction) -> Option<Self> {
        Some(self)
    }

    /// Value-only payload: no entity ids to remap.
    pub(crate) fn remap(self, _map: &IdRemapping) -> Self {
        self
    }
}

impl Normalize for BondConstraintForm {
    /// Normalize the inner value; kind and sub-key are preserved.
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Aromatic(b) => Self::Aromatic(b.normalize()?),
            Self::CisTransStereo(c) => Self::CisTransStereo(c.normalize()?),
            Self::RingMembership(m) => Self::RingMembership(m.normalize()?),
        })
    }
}

impl Lattice for BondConstraintForm {
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

/// Entry identity: discriminant + sub-key, BondConstraintsForm is ordered, unique by key.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondConstraintKey {
    Aromatic,
    CisTransStereo,
    RingMembership(RingScope),
}

/// Atom constraints container, ordered, unique by key, sorted flat vector storage.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondConstraintsForm(Vec<BondConstraintForm>);

impl BondConstraintsForm {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The bond's aromatic value, or `Undetermined` when no `Aromatic` constraint is present.
    pub fn aromatic(&self) -> BooleanForm {
        match self.get(BondConstraintKey::Aromatic) {
            Some(BondConstraintForm::Aromatic(b)) => *b,
            _ => BooleanForm::Undetermined,
        }
    }

    pub fn cis_trans_stereo(&self) -> Option<&CisTransStereoForm> {
        match self.get(BondConstraintKey::CisTransStereo) {
            Some(BondConstraintForm::CisTransStereo(c)) => Some(c),
            _ => None,
        }
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &NumForm)> {
        self.iter().filter_map(|c| match c {
            BondConstraintForm::RingMembership(m) => Some((m.scope, &m.count)),
            _ => None,
        })
    }

    fn ring_membership(&self, scope: RingScope) -> Option<&NumForm> {
        self.ring_memberships()
            .find(|(s, _)| *s == scope)
            .map(|(_, v)| v)
    }

    pub fn ring_count(&self) -> Option<&NumForm> {
        self.ring_membership(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> Option<&NumForm> {
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

    pub fn get(&self, key: BondConstraintKey) -> Option<&BondConstraintForm> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: BondConstraintForm) {
        match self.find(c.key()) {
            Ok(i) => self.0[i] = c,
            Err(i) => self.0.insert(i, c),
        }
    }

    /// Transactional write at one key: verify the current value `normalized_eq` `old` (both absent
    /// matches), then apply `new` (`Some` sets, `None` removes). `old`/`new` address the same key.
    /// `Err` on a key or old-value mismatch; the store is unchanged when it errors. The delta
    /// apply/undo primitive.
    pub fn compare_and_set(
        &mut self,
        old: Option<BondConstraintForm>,
        new: Option<BondConstraintForm>,
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
            (Some(current), Some(old)) => current.normalized_eq(old),
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

    pub fn remove(&mut self, key: BondConstraintKey) -> Option<BondConstraintForm> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = BondConstraintForm>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &BondConstraintsForm) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&BondConstraintForm) -> bool) {
        self.0.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl ExactSizeIterator<Item = BondConstraintForm> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, BondConstraintForm> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &MoleculeCompaction) -> Self {
        self
    }
}

impl Normalize for BondConstraintsForm {
    /// Normalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn normalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Normalize::normalize)
            .collect::<Result<Vec<BondConstraintForm>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl Lattice for BondConstraintsForm {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`BondConstraintForm::meet`; a `None` aborts the whole meet), an A-only /
    /// B-only key is kept (meet with the absent ⊤ is the value). Vacuous results are dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<BondConstraintForm> = Vec::new();
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
    /// (`BondConstraintForm::join`); a single-side key widens to the absent ⊤ and is dropped. The
    /// container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<BondConstraintForm> = Vec::new();
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
            BondConstraintForm::Aromatic(b) => b.matches(&target.aromatic()),
            BondConstraintForm::CisTransStereo(cts) => cts.matches(
                target
                    .cis_trans_stereo()
                    .unwrap_or(&CisTransStereoForm::Undetermined),
            ),
            BondConstraintForm::RingMembership(rm) => rm.count.matches(
                target
                    .ring_membership(rm.scope)
                    .unwrap_or(&NumForm::Undetermined),
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

impl FromIterator<BondConstraintForm> for BondConstraintsForm {
    fn from_iter<I: IntoIterator<Item = BondConstraintForm>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for BondConstraintsForm {
    type Item = BondConstraintForm;
    type IntoIter = IntoIter<BondConstraintForm>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<BondConstraintForm> for BondConstraintsForm {
    fn from(c: BondConstraintForm) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<BondConstraintForm>> for BondConstraintsForm {
    fn from(cs: Vec<BondConstraintForm>) -> Self {
        Self::from_iter(cs)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::{EdgeId, GraphCompaction, NodeId};

    use super::*;
    use crate::ir::stereo::{StereoCoset, StereoTerm};
    #[rustfmt::skip]
    #[rstest]
    #[case::ring_membership_all(BondConstraintForm::ring_membership(RingScope::All, 1), BondConstraintForm::ring_membership(RingScope::All, NumForm::Lit(1)))]
    #[case::ring_membership_size(BondConstraintForm::ring_membership(RingScope::Size(6), 1), BondConstraintForm::ring_membership(RingScope::Size(6), 1))]
    #[case::cis_trans_stereo(BondConstraintForm::cis_trans_stereo(CisTransStereoForm::NotStereo), BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo))]
    fn test_bond_constraint_form_constructors(
        #[case] actual: BondConstraintForm,
        #[case] expected: BondConstraintForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintKey::Aromatic)]
    #[case::ring_membership_all(BondConstraintForm::ring_membership(RingScope::All, 1), BondConstraintKey::RingMembership(RingScope::All))]
    #[case::ring_membership_size(BondConstraintForm::ring_membership(RingScope::Size(6), 1), BondConstraintKey::RingMembership(RingScope::Size(6)))]
    #[case::cis_trans_stereo(BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo), BondConstraintKey::CisTransStereo)]
    fn test_bond_constraint_form_key(#[case] c: BondConstraintForm, #[case] expected: BondConstraintKey) {
        assert_eq!(c.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::Aromatic(BooleanForm::Undetermined))]
    #[case::ring_membership_keeps_scope(BondConstraintForm::ring_membership(RingScope::Size(6), 1), BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined))]
    #[case::cis_trans(BondConstraintForm::CisTransStereo(CisTransStereoForm::stereo(1_u32)), BondConstraintForm::CisTransStereo(CisTransStereoForm::Undetermined))]
    fn test_bond_constraint_form_as_undetermined(#[case] c: BondConstraintForm, #[case] expected: BondConstraintForm) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), Ok(BondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::ring_count_litset_singleton(
        BondConstraintForm::RingMembership(RingMembershipForm::new(RingScope::All, NumForm::lit_set([2]))),
        Ok(BondConstraintForm::ring_membership(RingScope::All, 2)))]
    #[case::cis_trans_lifts_term(
        BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::term(StereoTerm::Lit(1)))),
        Ok(BondConstraintForm::cis_trans_stereo(CisTransStereoForm::stereo(1_u32))))]
    #[case::empty_litset_contradiction(
        BondConstraintForm::RingMembership(RingMembershipForm::new(RingScope::All, NumForm::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_bond_constraint_form_normalize(
        #[case] constraint: BondConstraintForm,
        #[case] expected: Result<BondConstraintForm, Contradiction>,
    ) {
        assert_eq!(constraint.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), false)]
    #[case::ring_membership_all_lit(BondConstraintForm::ring_membership(RingScope::All, 1), false)]
    #[case::ring_membership_all_undetermined(BondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined), true)]
    #[case::ring_membership_size_lit(BondConstraintForm::ring_membership(RingScope::Size(6), 1), false)]
    #[case::ring_membership_size_undetermined(BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined), true)]
    #[case::cis_trans_not_stereo(BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo), false)]
    #[case::cis_trans_undetermined(BondConstraintForm::CisTransStereo(CisTransStereoForm::Undetermined), true)]
    fn test_bond_constraint_form_is_undetermined(#[case] c: BondConstraintForm, #[case] expected: bool) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::Aromatic(BooleanForm::Undetermined), Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::same_key_incompatible(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::Aromatic(BooleanForm::Lit(false)), None)]
    #[case::different_key(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::All, 1), None)]
    fn test_bond_constraint_form_meet(#[case] a: BondConstraintForm, #[case] b: BondConstraintForm, #[case] expected: Option<BondConstraintForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_widens(BondConstraintForm::ring_membership(RingScope::All, 1), BondConstraintForm::ring_membership(RingScope::All, 2), Ok(BondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([1, 2]))))]
    #[case::different_key(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::All, 1), Err(NoJoin))]
    fn test_bond_constraint_form_join(#[case] a: BondConstraintForm, #[case] b: BondConstraintForm, #[case] expected: Result<BondConstraintForm, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::Aromatic(BooleanForm::Lit(true)), true)]
    #[case::same_key_incompatible(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::Aromatic(BooleanForm::Lit(false)), false)]
    #[case::different_key(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::All, 1), false)]
    fn test_bond_constraint_form_is_compatible(#[case] a: BondConstraintForm, #[case] b: BondConstraintForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_bond_constraints_form_new() {
        let cs = BondConstraintsForm::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    fn test_bond_constraints_form_iter() {
        let cs = BondConstraintsForm::from_iter([
            BondConstraintForm::ring_membership(RingScope::Size(6), 1),
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::ring_membership(RingScope::All, 1),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                BondConstraintForm::ring_membership(RingScope::All, 1),
                BondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true))], vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true))])]
    #[case::overwrite_same_key(vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::Aromatic(BooleanForm::Lit(false))], vec![BondConstraintForm::Aromatic(BooleanForm::Lit(false))])]
    #[case::vacuous_stores(vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::Aromatic(BooleanForm::Undetermined)], vec![BondConstraintForm::Aromatic(BooleanForm::Undetermined)])]
    #[case::new_key_sorts(vec![BondConstraintForm::ring_membership(RingScope::Size(6), 1), BondConstraintForm::Aromatic(BooleanForm::Lit(true))], vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::Size(6), 1)])]
    #[case::ring_overwrite_scope(vec![BondConstraintForm::ring_membership(RingScope::Size(6), 1), BondConstraintForm::ring_membership(RingScope::Size(6), 2)], vec![BondConstraintForm::ring_membership(RingScope::Size(6), 2)])]
    fn test_bond_constraints_form_set(#[case] sequence: Vec<BondConstraintForm>, #[case] expected: Vec<BondConstraintForm>) {
        let mut cs = BondConstraintsForm::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, BondConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite_shared(
        vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::All, 1)],
        vec![BondConstraintForm::Aromatic(BooleanForm::Lit(false))],
        vec![BondConstraintForm::Aromatic(BooleanForm::Lit(false)), BondConstraintForm::ring_membership(RingScope::All, 1)])]
    #[case::keeps_disjoint(
        vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true))],
        vec![BondConstraintForm::ring_membership(RingScope::All, 1)],
        vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::All, 1)])]
    #[case::vacuous_removes(
        vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::All, 1)],
        vec![BondConstraintForm::Aromatic(BooleanForm::Undetermined)],
        vec![BondConstraintForm::ring_membership(RingScope::All, 1)])]
    fn test_bond_constraints_form_update(#[case] initial: Vec<BondConstraintForm>, #[case] other: Vec<BondConstraintForm>, #[case] expected: Vec<BondConstraintForm>) {
        let mut cs = BondConstraintsForm::from_iter(initial);
        cs.update(&BondConstraintsForm::from_iter(other));
        assert_eq!(cs, BondConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![BondConstraintForm::ring_membership(RingScope::All, 1)], Some(BondConstraintForm::ring_membership(RingScope::All, 1)), Some(BondConstraintForm::ring_membership(RingScope::All, 2)), Ok(()), vec![BondConstraintForm::ring_membership(RingScope::All, 2)])]
    #[case::remove(vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true))], Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))), Ok(()), vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true))])]
    #[case::old_mismatch(vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true))], Some(BondConstraintForm::Aromatic(BooleanForm::Lit(false))), None, Err(Contradiction), vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true))])]
    #[case::key_mismatch(vec![], Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))), Some(BondConstraintForm::ring_membership(RingScope::All, 1)), Err(Contradiction), vec![])]
    fn test_bond_constraints_form_compare_and_set(
        #[case] initial: Vec<BondConstraintForm>,
        #[case] old: Option<BondConstraintForm>,
        #[case] new: Option<BondConstraintForm>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<BondConstraintForm>,
    ) {
        let mut cs = BondConstraintsForm::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, BondConstraintsForm::from_iter(expected_state));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(BondConstraintKey::Aromatic, true)]
    #[case::ring_all_present(BondConstraintKey::RingMembership(RingScope::All), true)]
    #[case::ring_size_present(BondConstraintKey::RingMembership(RingScope::Size(6)), true)]
    #[case::ring_size_absent(BondConstraintKey::RingMembership(RingScope::Size(5)), false)]
    #[case::cis_trans_absent(BondConstraintKey::CisTransStereo, false)]
    fn test_bond_constraints_form_contains(
        #[case] key: BondConstraintKey,
        #[case] expected: bool,
    ) {
        let cs = BondConstraintsForm::from_iter([
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::ring_membership(RingScope::All, 2),
            BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintKey::Aromatic, Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::ring_all(BondConstraintKey::RingMembership(RingScope::All), Some(BondConstraintForm::ring_membership(RingScope::All, 2)))]
    #[case::ring_size(BondConstraintKey::RingMembership(RingScope::Size(6)), Some(BondConstraintForm::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_size_absent(BondConstraintKey::RingMembership(RingScope::Size(5)), None)]
    fn test_bond_constraints_form_get(
        #[case] key: BondConstraintKey,
        #[case] expected: Option<BondConstraintForm>,
    ) {
        let cs = BondConstraintsForm::from_iter([
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::ring_membership(RingScope::All, 2),
            BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(key), expected.as_ref());
    }

    #[rstest]
    fn test_bond_constraints_form_remove() {
        let mut cs = BondConstraintsForm::from_iter([
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::ring_membership(RingScope::All, 2),
            BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        let removed = cs.remove(BondConstraintKey::RingMembership(RingScope::Size(6)));
        assert_eq!(
            removed,
            Some(BondConstraintForm::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(
            cs,
            BondConstraintsForm::from_iter([
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                BondConstraintForm::ring_membership(RingScope::All, 2),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &BondConstraintForm| matches!(c, BondConstraintForm::Aromatic(BooleanForm::Lit(true))) || matches!(c, BondConstraintForm::RingMembership(m) if m.scope == RingScope::Size(6)), vec![
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::Size(6), 1)])]
    #[case::all_dropped(|_: &BondConstraintForm| false, vec![])]
    fn test_bond_constraints_form_retain(
        #[case] predicate: impl FnMut(&BondConstraintForm) -> bool,
        #[case] expected: Vec<BondConstraintForm>,
    ) {
        let mut cs = BondConstraintsForm::from_iter([
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::ring_membership(RingScope::All, 1),
            BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        cs.retain(predicate);
        assert_eq!(cs, BondConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_bond_constraints_form_clear() {
        let mut cs =
            BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]);
        cs.clear();
        assert_eq!(cs, BondConstraintsForm::new());
    }

    #[rstest]
    fn test_bond_constraints_form_take() {
        let mut empty = BondConstraintsForm::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let mut cs = BondConstraintsForm::from_iter([
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        let mut taken = cs.take();
        assert_eq!(taken.len(), 2);
        assert_eq!(taken.size_hint(), (2, Some(2)));
        assert_eq!(
            taken.next(),
            Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        );
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(
            taken.next(),
            Some(BondConstraintForm::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(cs, BondConstraintsForm::new());
    }

    #[rstest]
    fn test_bond_constraints_form_compact() {
        let cs = BondConstraintsForm::from_iter([
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                vec![NodeId(0), NodeId(1), NodeId(2)],
                vec![EdgeId(0), EdgeId(1)],
            ),
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
        BondConstraintsForm::from_iter([
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined),
        ]),
        Ok(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))])))]
    #[case::normalizes_values(
        BondConstraintsForm::from_iter([
            BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::term(StereoTerm::Lit(1)))),
        ]),
        Ok(BondConstraintsForm::from_iter([BondConstraintForm::cis_trans_stereo(CisTransStereoForm::stereo(1_u32))])))]
    fn test_bond_constraints_form_normalize(
        #[case] constraints: BondConstraintsForm,
        #[case] expected: Result<BondConstraintsForm, Contradiction>,
    ) {
        assert_eq!(constraints.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys_kept(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, 1)]),
        Some(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::All, 1)])))]
    #[case::shared_key_meets(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Undetermined)]),
        Some(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))])))]
    #[case::shared_key_contradicts(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(false))]), None)]
    #[case::ring_size_unions(BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::Size(5), 1)]), BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::Size(6), 1)]),
        Some(BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::Size(5), 1), BondConstraintForm::ring_membership(RingScope::Size(6), 1)])))]
    #[case::prunes_vacuous(BondConstraintsForm::new(), BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Undetermined)]), Some(BondConstraintsForm::new()))]
    fn test_bond_constraints_form_meet(#[case] a: BondConstraintsForm, #[case] b: BondConstraintsForm, #[case] expected: Option<BondConstraintsForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::keeps_only_shared_keys(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::All, 1)]), BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]),
        BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]))]
    #[case::widens_value(BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, 1)]), BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, 2)]),
        BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([1, 2]))]))]
    #[case::incompatible_drops_to_undetermined(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(false))]), BondConstraintsForm::new())]
    fn test_bond_constraints_form_join(#[case] a: BondConstraintsForm, #[case] b: BondConstraintsForm, #[case] expected: BondConstraintsForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(BondConstraintsForm::new(), BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), true)]
    #[case::aromatic_required_present(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]),
        BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), true)]
    #[case::aromatic_required_absent(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), BondConstraintsForm::new(), false)]
    #[case::ring_membership_all_wildcard_matches_lit(BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined)]),
        BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, 1)]), true)]
    #[case::ring_membership_all_lit_mismatch(BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, 1)]),
        BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, 2)]), false)]
    #[case::ring_membership_size_subset(BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::Size(5), 1)]),
        BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::Size(5), 1), BondConstraintForm::ring_membership(RingScope::Size(6), 1)]), true)]
    #[case::ring_membership_size_not_in_target(BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::Size(7), 1)]),
        BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::Size(5), 1)]), false)]
    #[case::cis_trans_match(BondConstraintsForm::from_iter([BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo)]),
        BondConstraintsForm::from_iter([BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo)]), true)]
    #[case::cis_trans_pattern_more_specific(BondConstraintsForm::from_iter([BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo)]),
        BondConstraintsForm::new(), false)]
    fn test_bond_constraints_form_matches(
        #[case] pattern: BondConstraintsForm,
        #[case] target: BondConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, 1)]), true)]
    #[case::shared_key_compatible(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), true)]
    #[case::shared_key_incompatible(BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]), BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(false))]), false)]
    fn test_bond_constraints_form_is_compatible(#[case] a: BondConstraintsForm, #[case] b: BondConstraintsForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::All, 1)],
        vec![BondConstraintForm::Aromatic(BooleanForm::Lit(true)), BondConstraintForm::ring_membership(RingScope::All, 1)])]
    #[case::unique_kind_last_wins(vec![BondConstraintForm::cis_trans_stereo(CisTransStereoForm::Undetermined), BondConstraintForm::cis_trans_stereo(CisTransStereoForm::NotStereo)],
        vec![BondConstraintForm::cis_trans_stereo(CisTransStereoForm::NotStereo)])]
    #[case::ring_appends(vec![BondConstraintForm::ring_membership(RingScope::All, 1), BondConstraintForm::ring_membership(RingScope::Size(6), 1)],
        vec![BondConstraintForm::ring_membership(RingScope::All, 1), BondConstraintForm::ring_membership(RingScope::Size(6), 1)])]
    #[case::empty(vec![], vec![])]
    fn test_bond_constraints_form_from_iter(
        #[case] input: Vec<BondConstraintForm>,
        #[case] expected: Vec<BondConstraintForm>,
    ) {
        let cs = BondConstraintsForm::from_iter(input);
        assert_eq!(cs, BondConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_bond_constraints_form_into_iter() {
        let cs = BondConstraintsForm::from_iter([
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                BondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rstest]
    fn test_bond_constraints_form_from_bond_constraint() {
        let cs: BondConstraintsForm = BondConstraintForm::Aromatic(BooleanForm::Lit(true)).into();
        assert_eq!(
            cs,
            BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]),
        );
    }

    #[rstest]
    fn test_bond_constraints_form_from_vec() {
        let cs: BondConstraintsForm = vec![
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]
        .into();
        assert_eq!(
            cs,
            BondConstraintsForm::from_iter([
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                BondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]),
        );
    }
}
