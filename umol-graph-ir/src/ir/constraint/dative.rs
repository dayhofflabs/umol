//! Dative bond constraints.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use umol_perm::DynPermutation;

use super::super::boolean::BooleanForm;
use super::super::compact::MoleculeCompaction;
use super::super::constraint::ring::{RingMembershipForm, RingScope};
use super::super::error::{Contradiction, NoJoin};
use super::super::num::NumForm;
use super::super::traits::{FrameTransport, Lattice, Normalize};

/// Dative-bond constraint.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DativeBondConstraintForm {
    Aromatic(BooleanForm),
    /// Asserted ring count, optionally restricted by size. Derivation from topology requires a
    /// ring model that includes dative overlays rather than the localized atom-bond projection.
    RingMembership(RingMembershipForm),
}

impl DativeBondConstraintForm {
    pub fn aromatic(value: impl Into<BooleanForm>) -> Self {
        Self::Aromatic(value.into())
    }

    pub fn ring_membership(scope: RingScope, count: impl Into<NumForm>) -> Self {
        Self::RingMembership(RingMembershipForm::new(scope, count))
    }

    /// Dative bond constraint key, unique within a `DativeBondConstraintsForm` container.
    pub fn key(&self) -> DativeBondConstraintKey {
        match self {
            Self::Aromatic(_) => DativeBondConstraintKey::Aromatic,
            Self::RingMembership(m) => DativeBondConstraintKey::RingMembership(m.scope),
        }
    }

    /// Vacuous form of constraint key, used for removal.
    pub fn as_undetermined(&self) -> Self {
        match self {
            Self::Aromatic(_) => Self::Aromatic(BooleanForm::Undetermined),
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipForm::new(m.scope, NumForm::Undetermined))
            }
        }
    }

    /// Value-only payload: no entity ids to compact, so this never drops.
    pub fn compact(self, _compaction: &MoleculeCompaction) -> Option<Self> {
        Some(self)
    }

    pub(crate) fn uses_participant_frame(&self) -> bool {
        match self {
            Self::Aromatic(_) | Self::RingMembership(_) => false,
        }
    }
}

impl Normalize for DativeBondConstraintForm {
    /// Normalize the inner value; kind and sub-key are preserved.
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Aromatic(b) => Self::Aromatic(b.normalize()?),
            Self::RingMembership(m) => Self::RingMembership(m.normalize()?),
        })
    }
}

impl FrameTransport for DativeBondConstraintForm {
    type Action = DynPermutation;

    fn reframe_by(self, _action: &Self::Action) -> Option<Self> {
        Some(match self {
            Self::Aromatic(value) => Self::Aromatic(value),
            Self::RingMembership(value) => Self::RingMembership(value),
        })
    }
}

impl Lattice for DativeBondConstraintForm {
    fn is_undetermined(&self) -> bool {
        match self {
            Self::Aromatic(b) => b.is_undetermined(),
            Self::RingMembership(m) => m.is_undetermined(),
        }
    }

    fn is_ground(&self) -> bool {
        match self {
            Self::Aromatic(b) => b.is_ground(),
            Self::RingMembership(m) => m.is_ground(),
        }
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Aromatic(a), Self::Aromatic(b)) => a.meet(b).map(Self::Aromatic),
            (Self::RingMembership(a), Self::RingMembership(b)) => {
                a.meet(b).map(Self::RingMembership)
            }
            _ => None,
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        match (self, other) {
            (Self::Aromatic(a), Self::Aromatic(b)) => Ok(Self::Aromatic(a.join(b)?)),
            (Self::RingMembership(a), Self::RingMembership(b)) => {
                a.join(b).map(Self::RingMembership)
            }
            _ => Err(NoJoin),
        }
    }

    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Aromatic(a), Self::Aromatic(b)) => a.matches(b),
            (Self::RingMembership(a), Self::RingMembership(b)) => a.matches(b),
            _ => false,
        }
    }

    fn is_compatible(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Aromatic(a), Self::Aromatic(b)) => a.is_compatible(b),
            (Self::RingMembership(a), Self::RingMembership(b)) => a.is_compatible(b),
            _ => false,
        }
    }
}

/// Entry identity: discriminant + sub-key. Dative bond constraints are ordered, unique by key.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DativeBondConstraintKey {
    Aromatic,
    RingMembership(RingScope),
}
/// Dative bond constraints container, ordered, unique by key, sorted flat vector storage.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DativeBondConstraintsForm(Vec<DativeBondConstraintForm>);

impl DativeBondConstraintsForm {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The dative bond's aromatic value, or `Undetermined` when no `Aromatic` constraint is present.
    pub fn aromatic(&self) -> BooleanForm {
        match self.get(DativeBondConstraintKey::Aromatic) {
            Some(DativeBondConstraintForm::Aromatic(b)) => *b,
            _ => BooleanForm::Undetermined,
        }
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &NumForm)> {
        self.iter().filter_map(|c| match c {
            DativeBondConstraintForm::RingMembership(m) => Some((m.scope, &m.count)),
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

    fn find(&self, key: DativeBondConstraintKey) -> Result<usize, usize> {
        self.0.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains(&self, key: DativeBondConstraintKey) -> bool {
        self.find(key).is_ok()
    }

    pub fn get(&self, key: DativeBondConstraintKey) -> Option<&DativeBondConstraintForm> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: DativeBondConstraintForm) {
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
        old: Option<DativeBondConstraintForm>,
        new: Option<DativeBondConstraintForm>,
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

    pub fn remove(&mut self, key: DativeBondConstraintKey) -> Option<DativeBondConstraintForm> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = DativeBondConstraintForm>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other`: for each entry, a vacuous (`Undetermined`) one `remove`s its key, else
    /// `set`. Disjoint keys are kept.
    pub fn update(&mut self, other: &DativeBondConstraintsForm) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&DativeBondConstraintForm) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl ExactSizeIterator<Item = DativeBondConstraintForm> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, DativeBondConstraintForm> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &MoleculeCompaction) -> Self {
        self
    }
}

impl Normalize for DativeBondConstraintsForm {
    /// Normalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn normalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Normalize::normalize)
            .collect::<Result<Vec<DativeBondConstraintForm>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl FrameTransport for DativeBondConstraintsForm {
    type Action = DynPermutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        if !self
            .iter()
            .any(DativeBondConstraintForm::uses_participant_frame)
        {
            return Some(self);
        }
        self.into_iter()
            .map(|constraint| constraint.reframe_by(action))
            .collect()
    }
}

impl Lattice for DativeBondConstraintsForm {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`DativeBondConstraintForm::meet`; a `None` aborts the whole meet), an
    /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<DativeBondConstraintForm> = Vec::new();
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
    /// (`DativeBondConstraintForm::join`); a single-side key widens to the absent ⊤ and is dropped.
    /// The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<DativeBondConstraintForm> = Vec::new();
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

    /// Each value is matched on its own lattice; an empty pattern matches any target.
    fn matches(&self, target: &Self) -> bool {
        self.iter().all(|c| match c {
            DativeBondConstraintForm::Aromatic(b) => b.matches(&target.aromatic()),
            DativeBondConstraintForm::RingMembership(rm) => rm.count.matches(
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

impl FromIterator<DativeBondConstraintForm> for DativeBondConstraintsForm {
    fn from_iter<I: IntoIterator<Item = DativeBondConstraintForm>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for DativeBondConstraintsForm {
    type Item = DativeBondConstraintForm;
    type IntoIter = IntoIter<DativeBondConstraintForm>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<DativeBondConstraintForm> for DativeBondConstraintsForm {
    fn from(c: DativeBondConstraintForm) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<DativeBondConstraintForm>> for DativeBondConstraintsForm {
    fn from(cs: Vec<DativeBondConstraintForm>) -> Self {
        Self::from_iter(cs)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::{Compaction, EdgeId, GraphCompaction, NodeId};

    use super::*;
    #[rustfmt::skip]
    #[rstest]
    #[case::ring_membership_all(DativeBondConstraintForm::ring_membership(RingScope::All, 1), DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::Lit(1)))]
    #[case::ring_membership_size(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1))]
    fn test_dative_bond_constraint_form_constructors(
        #[case] actual: DativeBondConstraintForm,
        #[case] expected: DativeBondConstraintForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintKey::Aromatic)]
    #[case::ring_membership_all(DativeBondConstraintForm::ring_membership(RingScope::All, 1), DativeBondConstraintKey::RingMembership(RingScope::All))]
    #[case::ring_membership_size(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1), DativeBondConstraintKey::RingMembership(RingScope::Size(6)))]
    fn test_dative_bond_constraint_form_key(
        #[case] c: DativeBondConstraintForm,
        #[case] expected: DativeBondConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::Aromatic(BooleanForm::Undetermined))]
    #[case::ring_membership_keeps_scope(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1), DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined))]
    fn test_dative_bond_constraint_form_as_undetermined(#[case] c: DativeBondConstraintForm, #[case] expected: DativeBondConstraintForm) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), Ok(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::ring_count_litset_singleton(
        DativeBondConstraintForm::RingMembership(RingMembershipForm::new(RingScope::All, NumForm::lit_set([2]))),
        Ok(DativeBondConstraintForm::ring_membership(RingScope::All, 2)))]
    #[case::empty_litset_contradiction(
        DativeBondConstraintForm::RingMembership(RingMembershipForm::new(RingScope::All, NumForm::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_dative_bond_constraint_form_normalize(
        #[case] constraint: DativeBondConstraintForm,
        #[case] expected: Result<DativeBondConstraintForm, Contradiction>,
    ) {
        assert_eq!(constraint.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintForm::aromatic(true), false)]
    #[case::ring_membership(DativeBondConstraintForm::ring_membership(RingScope::All, 1), false)]
    fn test_dative_bond_constraint_form_uses_participant_frame(
        #[case] constraint: DativeBondConstraintForm,
        #[case] expected: bool,
    ) {
        assert_eq!(constraint.uses_participant_frame(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintForm::aromatic(true))]
    #[case::ring_membership(DativeBondConstraintForm::ring_membership(RingScope::All, 1))]
    fn test_dative_bond_constraint_form_reframe_by(#[case] constraint: DativeBondConstraintForm) {
        let action = DynPermutation::try_from(vec![1, 0]).expect("the action is a permutation");

        assert_eq!(constraint.clone().reframe_by(&action), Some(constraint));
    }

    #[rstest]
    fn test_dative_bond_constraints_form_reframe_by() {
        let constraints = DativeBondConstraintsForm::from(vec![
            DativeBondConstraintForm::aromatic(true),
            DativeBondConstraintForm::ring_membership(RingScope::All, 1),
        ]);
        let action = DynPermutation::try_from(vec![1, 0]).expect("the action is a permutation");

        assert_eq!(constraints.clone().reframe_by(&action), Some(constraints),);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), false)]
    #[case::ring_membership_all_lit(DativeBondConstraintForm::ring_membership(RingScope::All, 1), false)]
    #[case::ring_membership_all_undetermined(DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined), true)]
    #[case::ring_membership_size_lit(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1), false)]
    #[case::ring_membership_size_undetermined(DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined), true)]
    fn test_dative_bond_constraint_form_is_undetermined(#[case] c: DativeBondConstraintForm, #[case] expected: bool) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::Aromatic(BooleanForm::Undetermined), Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::same_key_incompatible(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false)), None)]
    #[case::different_key(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::All, 1), None)]
    fn test_dative_bond_constraint_form_meet(#[case] a: DativeBondConstraintForm, #[case] b: DativeBondConstraintForm, #[case] expected: Option<DativeBondConstraintForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_widens(DativeBondConstraintForm::ring_membership(RingScope::All, 1), DativeBondConstraintForm::ring_membership(RingScope::All, 2), Ok(DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([1, 2]))))]
    #[case::different_key(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::All, 1), Err(NoJoin))]
    fn test_dative_bond_constraint_form_join(#[case] a: DativeBondConstraintForm, #[case] b: DativeBondConstraintForm, #[case] expected: Result<DativeBondConstraintForm, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), true)]
    #[case::same_key_incompatible(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false)), false)]
    #[case::different_key(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::All, 1), false)]
    fn test_dative_bond_constraint_form_is_compatible(#[case] a: DativeBondConstraintForm, #[case] b: DativeBondConstraintForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_dative_bond_constraints_form_new() {
        let cs = DativeBondConstraintsForm::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    fn test_dative_bond_constraints_form_iter() {
        let cs = DativeBondConstraintsForm::from_iter([
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::All, 1),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                DativeBondConstraintForm::ring_membership(RingScope::All, 1),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))], vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))])]
    #[case::overwrite_same_key(vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false))], vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false))])]
    #[case::vacuous_stores(vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::Aromatic(BooleanForm::Undetermined)], vec![DativeBondConstraintForm::Aromatic(BooleanForm::Undetermined)])]
    #[case::new_key_sorts(vec![DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))], vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)])]
    #[case::ring_overwrite_scope(vec![DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 2)], vec![DativeBondConstraintForm::ring_membership(RingScope::Size(6), 2)])]
    fn test_dative_bond_constraints_form_set(#[case] sequence: Vec<DativeBondConstraintForm>, #[case] expected: Vec<DativeBondConstraintForm>) {
        let mut cs = DativeBondConstraintsForm::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, DativeBondConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite_shared(
        vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false))],
        vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false)), DativeBondConstraintForm::ring_membership(RingScope::All, 1)])]
    #[case::keeps_disjoint(
        vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))],
        vec![DativeBondConstraintForm::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::All, 1)])]
    #[case::vacuous_removes(
        vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraintForm::Aromatic(BooleanForm::Undetermined)],
        vec![DativeBondConstraintForm::ring_membership(RingScope::All, 1)])]
    fn test_dative_bond_constraints_form_update(#[case] initial: Vec<DativeBondConstraintForm>, #[case] other: Vec<DativeBondConstraintForm>, #[case] expected: Vec<DativeBondConstraintForm>) {
        let mut cs = DativeBondConstraintsForm::from_iter(initial);
        cs.update(&DativeBondConstraintsForm::from_iter(other));
        assert_eq!(cs, DativeBondConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![DativeBondConstraintForm::ring_membership(RingScope::All, 1)], Some(DativeBondConstraintForm::ring_membership(RingScope::All, 1)), Some(DativeBondConstraintForm::ring_membership(RingScope::All, 2)), Ok(()), vec![DativeBondConstraintForm::ring_membership(RingScope::All, 2)])]
    #[case::remove(vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))], Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))), Ok(()), vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))])]
    #[case::old_mismatch(vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))], Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false))), None, Err(Contradiction), vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))])]
    #[case::key_mismatch(vec![], Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))), Some(DativeBondConstraintForm::ring_membership(RingScope::All, 1)), Err(Contradiction), vec![])]
    fn test_dative_bond_constraints_form_compare_and_set(
        #[case] initial: Vec<DativeBondConstraintForm>,
        #[case] old: Option<DativeBondConstraintForm>,
        #[case] new: Option<DativeBondConstraintForm>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<DativeBondConstraintForm>,
    ) {
        let mut cs = DativeBondConstraintsForm::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, DativeBondConstraintsForm::from_iter(expected_state));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(DativeBondConstraintKey::Aromatic, true)]
    #[case::ring_all_present(DativeBondConstraintKey::RingMembership(RingScope::All), true)]
    #[case::ring_size_present(DativeBondConstraintKey::RingMembership(RingScope::Size(6)), true)]
    #[case::ring_size_absent(DativeBondConstraintKey::RingMembership(RingScope::Size(5)), false)]
    fn test_dative_bond_constraints_form_contains(
        #[case] key: DativeBondConstraintKey,
        #[case] expected: bool,
    ) {
        let cs = DativeBondConstraintsForm::from_iter([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::All, 2),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintKey::Aromatic, Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::ring_all(DativeBondConstraintKey::RingMembership(RingScope::All), Some(DativeBondConstraintForm::ring_membership(RingScope::All, 2)))]
    #[case::ring_size(DativeBondConstraintKey::RingMembership(RingScope::Size(6)), Some(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_size_absent(DativeBondConstraintKey::RingMembership(RingScope::Size(5)), None)]
    fn test_dative_bond_constraints_form_get(
        #[case] key: DativeBondConstraintKey,
        #[case] expected: Option<DativeBondConstraintForm>,
    ) {
        let cs = DativeBondConstraintsForm::from_iter([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::All, 2),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(key), expected.as_ref());
    }

    #[rstest]
    fn test_dative_bond_constraints_form_remove() {
        let mut cs = DativeBondConstraintsForm::from_iter([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::All, 2),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        let removed = cs.remove(DativeBondConstraintKey::RingMembership(RingScope::Size(6)));
        assert_eq!(
            removed,
            Some(DativeBondConstraintForm::ring_membership(
                RingScope::Size(6),
                1
            )),
        );
        assert_eq!(
            cs,
            DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                DativeBondConstraintForm::ring_membership(RingScope::All, 2),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &DativeBondConstraintForm| matches!(c, DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))) || matches!(c, DativeBondConstraintForm::RingMembership(m) if m.scope == RingScope::Size(6)), vec![
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)])]
    #[case::all_dropped(|_: &DativeBondConstraintForm| false, vec![])]
    fn test_dative_bond_constraints_form_retain(
        #[case] predicate: impl FnMut(&DativeBondConstraintForm) -> bool,
        #[case] expected: Vec<DativeBondConstraintForm>,
    ) {
        let mut cs = DativeBondConstraintsForm::from_iter([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::All, 1),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        cs.retain(predicate);
        assert_eq!(cs, DativeBondConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_dative_bond_constraints_form_clear() {
        let mut cs = DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(
            BooleanForm::Lit(true),
        )]);
        cs.clear();
        assert_eq!(cs, DativeBondConstraintsForm::new());
    }

    #[rstest]
    fn test_dative_bond_constraints_form_take() {
        let mut empty = DativeBondConstraintsForm::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let mut cs = DativeBondConstraintsForm::from_iter([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        let mut taken = cs.take();
        assert_eq!(taken.len(), 2);
        assert_eq!(taken.size_hint(), (2, Some(2)));
        assert_eq!(
            taken.next(),
            Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        );
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(
            taken.next(),
            Some(DativeBondConstraintForm::ring_membership(
                RingScope::Size(6),
                1,
            )),
        );
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(cs, DativeBondConstraintsForm::new());
    }

    #[rstest]
    fn test_dative_bond_constraints_form_compact() {
        let cs = DativeBondConstraintsForm::from_iter([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::new(2, vec![NodeId(1)]).unwrap(),
                Compaction::new(2, vec![EdgeId(1)]).unwrap(),
            ),
            Compaction::identity(0),
            Compaction::identity(0),
            Compaction::identity(0),
            Compaction::identity(0),
            Compaction::identity(0),
            Compaction::identity(0),
        );
        assert_eq!(cs.clone().compact(&compaction), cs);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::drop_vacuous(
        DativeBondConstraintsForm::from_iter([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined),
        ]),
        Ok(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))])))]
    #[case::normalizes_values(
        DativeBondConstraintsForm::from_iter([
            DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([2])),
        ]),
        Ok(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::All, 2)])))]
    fn test_dative_bond_constraints_form_normalize(
        #[case] constraints: DativeBondConstraintsForm,
        #[case] expected: Result<DativeBondConstraintsForm, Contradiction>,
    ) {
        assert_eq!(constraints.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys_kept(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::All, 1)]),
        Some(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::All, 1)])))]
    #[case::shared_key_meets(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Undetermined)]),
        Some(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))])))]
    #[case::shared_key_contradicts(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false))]), None)]
    #[case::ring_size_unions(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1)]), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)]),
        Some(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)])))]
    #[case::prunes_vacuous(DativeBondConstraintsForm::new(), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Undetermined)]), Some(DativeBondConstraintsForm::new()))]
    fn test_dative_bond_constraints_form_meet(#[case] a: DativeBondConstraintsForm, #[case] b: DativeBondConstraintsForm, #[case] expected: Option<DativeBondConstraintsForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::keeps_only_shared_keys(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::All, 1)]), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]),
        DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]))]
    #[case::widens_value(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::All, 1)]), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::All, 2)]),
        DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([1, 2]))]))]
    #[case::incompatible_drops_to_undetermined(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false))]), DativeBondConstraintsForm::new())]
    fn test_dative_bond_constraints_form_join(#[case] a: DativeBondConstraintsForm, #[case] b: DativeBondConstraintsForm, #[case] expected: DativeBondConstraintsForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(DativeBondConstraintsForm::new(), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), true)]
    #[case::aromatic_required_present(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]),
        DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), true)]
    #[case::aromatic_required_absent(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), DativeBondConstraintsForm::new(), false)]
    #[case::ring_membership_all_wildcard_matches_lit(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined)]),
        DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::All, 1)]), true)]
    #[case::ring_membership_all_lit_mismatch(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::All, 1)]),
        DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::All, 2)]), false)]
    #[case::ring_membership_size_subset(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1)]),
        DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)]), true)]
    #[case::ring_membership_size_not_in_target(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::Size(7), 1)]),
        DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1)]), false)]
    fn test_dative_bond_constraints_form_matches(
        #[case] pattern: DativeBondConstraintsForm,
        #[case] target: DativeBondConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(RingScope::All, 1)]), true)]
    #[case::shared_key_compatible(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), true)]
    #[case::shared_key_incompatible(DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))]), DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false))]), false)]
    fn test_dative_bond_constraints_form_is_compatible(#[case] a: DativeBondConstraintsForm, #[case] b: DativeBondConstraintsForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::ring_membership(RingScope::All, 1)])]
    #[case::unique_kind_last_wins(vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false))],
        vec![DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false))])]
    #[case::ring_appends(vec![DativeBondConstraintForm::ring_membership(RingScope::All, 1), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)],
        vec![DativeBondConstraintForm::ring_membership(RingScope::All, 1), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)])]
    #[case::empty(vec![], vec![])]
    fn test_dative_bond_constraints_form_from_iter(
        #[case] input: Vec<DativeBondConstraintForm>,
        #[case] expected: Vec<DativeBondConstraintForm>,
    ) {
        let cs = DativeBondConstraintsForm::from_iter(input);
        assert_eq!(cs, DativeBondConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_dative_bond_constraints_form_into_iter() {
        let cs = DativeBondConstraintsForm::from_iter([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_form_from_dative_bond_constraint() {
        let cs: DativeBondConstraintsForm =
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)).into();
        assert_eq!(
            cs,
            DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::Aromatic(
                BooleanForm::Lit(true)
            )]),
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_form_from_vec() {
        let cs: DativeBondConstraintsForm = vec![
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]
        .into();
        assert_eq!(
            cs,
            DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]),
        );
    }
}
