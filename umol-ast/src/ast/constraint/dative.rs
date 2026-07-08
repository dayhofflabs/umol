//! Dative bond constraints.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use strum::EnumDiscriminants;

use super::super::boolean::BooleanAst;
use super::super::constraint::ring::{RingMembershipAst, RingScope};
use super::super::error::{Contradiction, NoJoin};
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Dative-bond constraint.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, EnumDiscriminants)]
#[strum_discriminants(name(DativeBondConstraintKind), derive(Hash))]
pub enum DativeBondConstraint {
    Aromatic(BooleanAst),
    RingMembership(RingMembershipAst),
}

impl DativeBondConstraint {
    pub fn aromatic(value: impl Into<BooleanAst>) -> Self {
        Self::Aromatic(value.into())
    }

    pub fn ring_membership(scope: RingScope, count: impl Into<ValueAst>) -> Self {
        Self::RingMembership(RingMembershipAst::new(scope, count))
    }

    pub fn kind(&self) -> DativeBondConstraintKind {
        self.into()
    }

    /// Dative bond constraint key, unique within a `DativeBondConstraints` container.
    pub fn key(&self) -> DativeBondConstraintKey {
        match self {
            Self::Aromatic(_) => DativeBondConstraintKey::Aromatic,
            Self::RingMembership(m) => DativeBondConstraintKey::RingMembership(m.scope),
        }
    }

    /// Vacuous form of constraint key, used for removal.
    pub fn as_undetermined(&self) -> Self {
        match self {
            Self::Aromatic(_) => Self::Aromatic(BooleanAst::Undetermined),
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

impl Lattice for DativeBondConstraint {
    fn is_undetermined(&self) -> bool {
        match self {
            Self::Aromatic(b) => b.is_undetermined(),
            Self::RingMembership(m) => m.count.is_undetermined(),
        }
    }

    fn is_ground(&self) -> bool {
        match self {
            Self::Aromatic(b) => b.is_ground(),
            Self::RingMembership(m) => m.count.is_ground(),
        }
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Aromatic(a), Self::Aromatic(b)) => a.meet(b).map(Self::Aromatic),
            (Self::RingMembership(a), Self::RingMembership(b)) => a
                .count
                .meet(&b.count)
                .map(|count| Self::RingMembership(RingMembershipAst::new(a.scope, count))),
            _ => None,
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        match (self, other) {
            (Self::Aromatic(a), Self::Aromatic(b)) => Ok(Self::Aromatic(a.join(b)?)),
            (Self::RingMembership(a), Self::RingMembership(b)) => Ok(Self::RingMembership(
                RingMembershipAst::new(a.scope, a.count.join(&b.count)?),
            )),
            _ => Err(NoJoin),
        }
    }

    fn is_compatible(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Aromatic(a), Self::Aromatic(b)) => a.is_compatible(b),
            (Self::RingMembership(a), Self::RingMembership(b)) => a.count.is_compatible(&b.count),
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
pub struct DativeBondConstraints(Vec<DativeBondConstraint>);

impl DativeBondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The dative bond's aromatic value, or `Undetermined` when no `Aromatic` constraint is present.
    pub fn aromatic(&self) -> BooleanAst {
        match self.get(DativeBondConstraintKey::Aromatic) {
            Some(DativeBondConstraint::Aromatic(b)) => *b,
            _ => BooleanAst::Undetermined,
        }
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &ValueAst)> {
        self.iter().filter_map(|c| match c {
            DativeBondConstraint::RingMembership(m) => Some((m.scope, &m.count)),
            _ => None,
        })
    }

    fn ring_membership(&self, scope: RingScope) -> ValueAst {
        self.ring_memberships()
            .find(|(s, _)| *s == scope)
            .map(|(_, v)| v.clone())
            .unwrap_or(ValueAst::Undetermined)
    }

    pub fn ring_count(&self) -> ValueAst {
        self.ring_membership(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> ValueAst {
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

    pub fn get(&self, key: DativeBondConstraintKey) -> Option<&DativeBondConstraint> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: DativeBondConstraint) {
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
        old: Option<DativeBondConstraint>,
        new: Option<DativeBondConstraint>,
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

    pub fn remove(&mut self, key: DativeBondConstraintKey) -> Option<DativeBondConstraint> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = DativeBondConstraint>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other`: for each entry, a vacuous (`Undetermined`) one `remove`s its key, else
    /// `set`. Disjoint keys are kept.
    pub fn update(&mut self, other: &DativeBondConstraints) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
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

    pub fn iter(&self) -> Iter<'_, DativeBondConstraint> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &IdCompaction) -> Self {
        self
    }
}

impl Canonicalize for DativeBondConstraints {
    /// Canonicalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Canonicalize::canonicalize)
            .collect::<Result<Vec<DativeBondConstraint>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl Lattice for DativeBondConstraints {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`DativeBondConstraint::meet`; a `None` aborts the whole meet), an
    /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<DativeBondConstraint> = Vec::new();
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
    /// (`DativeBondConstraint::join`); a single-side key widens to the absent ⊤ and is dropped.
    /// The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<DativeBondConstraint> = Vec::new();
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
        self.aromatic().matches(&target.aromatic())
            && self
                .ring_memberships()
                .all(|(scope, v)| v.matches(&target.ring_membership(scope)))
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

impl FromIterator<DativeBondConstraint> for DativeBondConstraints {
    fn from_iter<I: IntoIterator<Item = DativeBondConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
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
    use umol_graph_core::Compaction;

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
    #[case::aromatic(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::Aromatic(BooleanAst::Undetermined))]
    #[case::ring_membership_keeps_scope(DativeBondConstraint::ring_membership(RingScope::Size(6), 1), DativeBondConstraint::ring_membership(RingScope::Size(6), ValueAst::Undetermined))]
    fn test_dative_bond_constraint_as_undetermined(#[case] c: DativeBondConstraint, #[case] expected: DativeBondConstraint) {
        assert_eq!(c.as_undetermined(), expected);
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

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), false)]
    #[case::ring_membership_all_lit(DativeBondConstraint::ring_membership(RingScope::All, 1), false)]
    #[case::ring_membership_all_undetermined(DativeBondConstraint::ring_membership(RingScope::All, ValueAst::Undetermined), true)]
    #[case::ring_membership_size_lit(DativeBondConstraint::ring_membership(RingScope::Size(6), 1), false)]
    #[case::ring_membership_size_undetermined(DativeBondConstraint::ring_membership(RingScope::Size(6), ValueAst::Undetermined), true)]
    fn test_dative_bond_constraint_is_undetermined(#[case] c: DativeBondConstraint, #[case] expected: bool) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::Aromatic(BooleanAst::Undetermined), Some(DativeBondConstraint::Aromatic(BooleanAst::Lit(true))))]
    #[case::same_key_incompatible(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::Aromatic(BooleanAst::Lit(false)), None)]
    #[case::different_key(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1), None)]
    fn test_dative_bond_constraint_meet(#[case] a: DativeBondConstraint, #[case] b: DativeBondConstraint, #[case] expected: Option<DativeBondConstraint>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_widens(DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::All, 2), Ok(DativeBondConstraint::ring_membership(RingScope::All, ValueAst::lit_set([1, 2]))))]
    #[case::different_key(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1), Err(NoJoin))]
    fn test_dative_bond_constraint_join(#[case] a: DativeBondConstraint, #[case] b: DativeBondConstraint, #[case] expected: Result<DativeBondConstraint, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), true)]
    #[case::same_key_incompatible(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::Aromatic(BooleanAst::Lit(false)), false)]
    #[case::different_key(DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1), false)]
    fn test_dative_bond_constraint_is_compatible(#[case] a: DativeBondConstraint, #[case] b: DativeBondConstraint, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_dative_bond_constraints_new() {
        let cs = DativeBondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
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
    #[case::fresh(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))], vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))])]
    #[case::overwrite_same_key(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::Aromatic(BooleanAst::Lit(false))], vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(false))])]
    #[case::vacuous_stores(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::Aromatic(BooleanAst::Undetermined)], vec![DativeBondConstraint::Aromatic(BooleanAst::Undetermined)])]
    #[case::new_key_sorts(vec![DativeBondConstraint::ring_membership(RingScope::Size(6), 1), DativeBondConstraint::Aromatic(BooleanAst::Lit(true))], vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::ring_overwrite_scope(vec![DativeBondConstraint::ring_membership(RingScope::Size(6), 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 2)], vec![DativeBondConstraint::ring_membership(RingScope::Size(6), 2)])]
    fn test_dative_bond_constraints_set(#[case] sequence: Vec<DativeBondConstraint>, #[case] expected: Vec<DativeBondConstraint>) {
        let mut cs = DativeBondConstraints::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, DativeBondConstraints::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite_shared(
        vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(false))],
        vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(false)), DativeBondConstraint::ring_membership(RingScope::All, 1)])]
    #[case::keeps_disjoint(
        vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))],
        vec![DativeBondConstraint::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1)])]
    #[case::vacuous_removes(
        vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraint::Aromatic(BooleanAst::Undetermined)],
        vec![DativeBondConstraint::ring_membership(RingScope::All, 1)])]
    fn test_dative_bond_constraints_update(#[case] initial: Vec<DativeBondConstraint>, #[case] other: Vec<DativeBondConstraint>, #[case] expected: Vec<DativeBondConstraint>) {
        let mut cs = DativeBondConstraints::from_iter(initial);
        cs.update(&DativeBondConstraints::from_iter(other));
        assert_eq!(cs, DativeBondConstraints::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![DativeBondConstraint::ring_membership(RingScope::All, 1)], Some(DativeBondConstraint::ring_membership(RingScope::All, 1)), Some(DativeBondConstraint::ring_membership(RingScope::All, 2)), Ok(()), vec![DativeBondConstraint::ring_membership(RingScope::All, 2)])]
    #[case::remove(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))], Some(DativeBondConstraint::Aromatic(BooleanAst::Lit(true))), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(DativeBondConstraint::Aromatic(BooleanAst::Lit(true))), Ok(()), vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))])]
    #[case::old_mismatch(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))], Some(DativeBondConstraint::Aromatic(BooleanAst::Lit(false))), None, Err(Contradiction), vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true))])]
    #[case::key_mismatch(vec![], Some(DativeBondConstraint::Aromatic(BooleanAst::Lit(true))), Some(DativeBondConstraint::ring_membership(RingScope::All, 1)), Err(Contradiction), vec![])]
    fn test_dative_bond_constraints_compare_and_set(
        #[case] initial: Vec<DativeBondConstraint>,
        #[case] old: Option<DativeBondConstraint>,
        #[case] new: Option<DativeBondConstraint>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<DativeBondConstraint>,
    ) {
        let mut cs = DativeBondConstraints::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, DativeBondConstraints::from_iter(expected_state));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(DativeBondConstraintKey::Aromatic, true)]
    #[case::ring_all_present(DativeBondConstraintKey::RingMembership(RingScope::All), true)]
    #[case::ring_size_present(DativeBondConstraintKey::RingMembership(RingScope::Size(6)), true)]
    #[case::ring_size_absent(DativeBondConstraintKey::RingMembership(RingScope::Size(5)), false)]
    fn test_dative_bond_constraints_contains(
        #[case] key: DativeBondConstraintKey,
        #[case] expected: bool,
    ) {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::All, 2),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintKey::Aromatic, Some(DativeBondConstraint::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_all(DativeBondConstraintKey::RingMembership(RingScope::All), Some(DativeBondConstraint::ring_membership(RingScope::All, 2)))]
    #[case::ring_size(DativeBondConstraintKey::RingMembership(RingScope::Size(6)), Some(DativeBondConstraint::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_size_absent(DativeBondConstraintKey::RingMembership(RingScope::Size(5)), None)]
    fn test_dative_bond_constraints_get(
        #[case] key: DativeBondConstraintKey,
        #[case] expected: Option<DativeBondConstraint>,
    ) {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::All, 2),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(key), expected.as_ref());
    }

    #[rstest]
    fn test_dative_bond_constraints_remove() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::All, 2),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let removed = cs.remove(DativeBondConstraintKey::RingMembership(RingScope::Size(6)));
        assert_eq!(
            removed,
            Some(DativeBondConstraint::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(
            cs,
            DativeBondConstraints::from_iter([
                DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraint::ring_membership(RingScope::All, 2),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &DativeBondConstraint| matches!(c, DativeBondConstraint::Aromatic(BooleanAst::Lit(true))) || matches!(c, DativeBondConstraint::RingMembership(m) if m.scope == RingScope::Size(6)), vec![
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])]
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
        assert_eq!(cs, DativeBondConstraints::from_iter(expected));
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
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
        assert_eq!(cs, DativeBondConstraints::new());
    }

    #[rstest]
    fn test_dative_bond_constraints_compact() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let compaction = IdCompaction::new(
            Compaction::new(vec![1], vec![1]),
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
        DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraint::ring_membership(RingScope::All, ValueAst::Undetermined),
        ]),
        Ok(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))])))]
    #[case::canonicalizes_values(
        DativeBondConstraints::from_iter([
            DativeBondConstraint::ring_membership(RingScope::All, ValueAst::lit_set([2])),
        ]),
        Ok(DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, 2)])))]
    fn test_dative_bond_constraints_canonicalize(
        #[case] constraints: DativeBondConstraints,
        #[case] expected: Result<DativeBondConstraints, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys_kept(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, 1)]),
        Some(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1)])))]
    #[case::shared_key_meets(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Undetermined)]),
        Some(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))])))]
    #[case::shared_key_contradicts(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(false))]), None)]
    #[case::ring_size_unions(DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::Size(5), 1)]), DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::Size(6), 1)]),
        Some(DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::Size(5), 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])))]
    #[case::prunes_vacuous(DativeBondConstraints::new(), DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Undetermined)]), Some(DativeBondConstraints::new()))]
    fn test_dative_bond_constraints_meet(#[case] a: DativeBondConstraints, #[case] b: DativeBondConstraints, #[case] expected: Option<DativeBondConstraints>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::keeps_only_shared_keys(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1)]), DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]),
        DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]))]
    #[case::widens_value(DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, 1)]), DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, 2)]),
        DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, ValueAst::lit_set([1, 2]))]))]
    #[case::incompatible_drops_to_undetermined(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(false))]), DativeBondConstraints::new())]
    fn test_dative_bond_constraints_join(#[case] a: DativeBondConstraints, #[case] b: DativeBondConstraints, #[case] expected: DativeBondConstraints) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(DativeBondConstraints::new(), DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::aromatic_required_present(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]),
        DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::aromatic_required_absent(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraints::new(), false)]
    #[case::ring_membership_all_wildcard_matches_lit(DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, ValueAst::Undetermined)]),
        DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, 1)]), true)]
    #[case::ring_membership_all_lit_mismatch(DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, 1)]),
        DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, 2)]), false)]
    #[case::ring_membership_size_subset(DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::Size(5), 1)]),
        DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::Size(5), 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)]), true)]
    #[case::ring_membership_size_not_in_target(DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::Size(7), 1)]),
        DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::Size(5), 1)]), false)]
    fn test_dative_bond_constraints_matches(
        #[case] pattern: DativeBondConstraints,
        #[case] target: DativeBondConstraints,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(RingScope::All, 1)]), true)]
    #[case::shared_key_compatible(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::shared_key_incompatible(DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(false))]), false)]
    fn test_dative_bond_constraints_is_compatible(#[case] a: DativeBondConstraints, #[case] b: DativeBondConstraints, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::ring_membership(RingScope::All, 1)])]
    #[case::unique_kind_last_wins(vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(true)), DativeBondConstraint::Aromatic(BooleanAst::Lit(false))],
        vec![DativeBondConstraint::Aromatic(BooleanAst::Lit(false))])]
    #[case::ring_appends(vec![DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)],
        vec![DativeBondConstraint::ring_membership(RingScope::All, 1), DativeBondConstraint::ring_membership(RingScope::Size(6), 1)])]
    #[case::empty(vec![], vec![])]
    fn test_dative_bond_constraints_from_iter(
        #[case] input: Vec<DativeBondConstraint>,
        #[case] expected: Vec<DativeBondConstraint>,
    ) {
        let cs = DativeBondConstraints::from_iter(input);
        assert_eq!(cs, DativeBondConstraints::from_iter(expected));
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
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_from_dative_bond_constraint() {
        let cs: DativeBondConstraints =
            DativeBondConstraint::Aromatic(BooleanAst::Lit(true)).into();
        assert_eq!(
            cs,
            DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic(BooleanAst::Lit(
                true
            ))]),
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
            cs,
            DativeBondConstraints::from_iter([
                DativeBondConstraint::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraint::ring_membership(RingScope::Size(6), 1),
            ]),
        );
    }
}
