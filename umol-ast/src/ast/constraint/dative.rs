//! Dative bond constraints.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use super::super::boolean::BooleanAst;
use super::super::constraint::ring::{RingMembershipAst, RingScope};
use super::super::error::{Contradiction, NoJoin};
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Dative-bond constraint.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DativeBondConstraintAst {
    Aromatic(BooleanAst),
    /// Asserted ring count, optionally restricted by size. Derivation from topology requires a
    /// ring model that includes dative overlays rather than the localized atom-bond projection.
    RingMembership(RingMembershipAst),
}

impl DativeBondConstraintAst {
    pub fn aromatic(value: impl Into<BooleanAst>) -> Self {
        Self::Aromatic(value.into())
    }

    pub fn ring_membership(scope: RingScope, count: impl Into<ValueAst>) -> Self {
        Self::RingMembership(RingMembershipAst::new(scope, count))
    }

    /// Dative bond constraint key, unique within a `DativeBondConstraintsAst` container.
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

impl Canonicalize for DativeBondConstraintAst {
    /// Canonicalize the inner value; kind and sub-key are preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Aromatic(b) => Self::Aromatic(b.canonicalize()?),
            Self::RingMembership(m) => Self::RingMembership(m.canonicalize()?),
        })
    }
}

impl Lattice for DativeBondConstraintAst {
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
pub struct DativeBondConstraintsAst(Vec<DativeBondConstraintAst>);

impl DativeBondConstraintsAst {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The dative bond's aromatic value, or `Undetermined` when no `Aromatic` constraint is present.
    pub fn aromatic(&self) -> BooleanAst {
        match self.get(DativeBondConstraintKey::Aromatic) {
            Some(DativeBondConstraintAst::Aromatic(b)) => *b,
            _ => BooleanAst::Undetermined,
        }
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &ValueAst)> {
        self.iter().filter_map(|c| match c {
            DativeBondConstraintAst::RingMembership(m) => Some((m.scope, &m.count)),
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

    fn find(&self, key: DativeBondConstraintKey) -> Result<usize, usize> {
        self.0.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains(&self, key: DativeBondConstraintKey) -> bool {
        self.find(key).is_ok()
    }

    pub fn get(&self, key: DativeBondConstraintKey) -> Option<&DativeBondConstraintAst> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: DativeBondConstraintAst) {
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
        old: Option<DativeBondConstraintAst>,
        new: Option<DativeBondConstraintAst>,
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

    pub fn remove(&mut self, key: DativeBondConstraintKey) -> Option<DativeBondConstraintAst> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = DativeBondConstraintAst>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other`: for each entry, a vacuous (`Undetermined`) one `remove`s its key, else
    /// `set`. Disjoint keys are kept.
    pub fn update(&mut self, other: &DativeBondConstraintsAst) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&DativeBondConstraintAst) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = DativeBondConstraintAst> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, DativeBondConstraintAst> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &IdCompaction) -> Self {
        self
    }
}

impl Canonicalize for DativeBondConstraintsAst {
    /// Canonicalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Canonicalize::canonicalize)
            .collect::<Result<Vec<DativeBondConstraintAst>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl Lattice for DativeBondConstraintsAst {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`DativeBondConstraintAst::meet`; a `None` aborts the whole meet), an
    /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<DativeBondConstraintAst> = Vec::new();
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
    /// (`DativeBondConstraintAst::join`); a single-side key widens to the absent ⊤ and is dropped.
    /// The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<DativeBondConstraintAst> = Vec::new();
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
            DativeBondConstraintAst::Aromatic(b) => b.matches(&target.aromatic()),
            DativeBondConstraintAst::RingMembership(rm) => rm.count.matches(
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

impl FromIterator<DativeBondConstraintAst> for DativeBondConstraintsAst {
    fn from_iter<I: IntoIterator<Item = DativeBondConstraintAst>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for DativeBondConstraintsAst {
    type Item = DativeBondConstraintAst;
    type IntoIter = IntoIter<DativeBondConstraintAst>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<DativeBondConstraintAst> for DativeBondConstraintsAst {
    fn from(c: DativeBondConstraintAst) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<DativeBondConstraintAst>> for DativeBondConstraintsAst {
    fn from(cs: Vec<DativeBondConstraintAst>) -> Self {
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
    #[case::ring_membership_all(DativeBondConstraintAst::ring_membership(RingScope::All, 1), DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1)))]
    #[case::ring_membership_size(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1))]
    fn test_dative_bond_constraint_ast_ast_constructors(
        #[case] actual: DativeBondConstraintAst,
        #[case] expected: DativeBondConstraintAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintKey::Aromatic)]
    #[case::ring_membership_all(DativeBondConstraintAst::ring_membership(RingScope::All, 1), DativeBondConstraintKey::RingMembership(RingScope::All))]
    #[case::ring_membership_size(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1), DativeBondConstraintKey::RingMembership(RingScope::Size(6)))]
    fn test_dative_bond_constraint_ast_ast_key(
        #[case] c: DativeBondConstraintAst,
        #[case] expected: DativeBondConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined))]
    #[case::ring_membership_keeps_scope(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1), DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined))]
    fn test_dative_bond_constraint_ast_as_undetermined(#[case] c: DativeBondConstraintAst, #[case] expected: DativeBondConstraintAst) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), Ok(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_count_litset_singleton(
        DativeBondConstraintAst::RingMembership(RingMembershipAst::new(RingScope::All, ValueAst::lit_set([2]))),
        Ok(DativeBondConstraintAst::ring_membership(RingScope::All, 2)))]
    #[case::empty_litset_contradiction(
        DativeBondConstraintAst::RingMembership(RingMembershipAst::new(RingScope::All, ValueAst::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_dative_bond_constraint_ast_canonicalize(
        #[case] constraint: DativeBondConstraintAst,
        #[case] expected: Result<DativeBondConstraintAst, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), false)]
    #[case::ring_membership_all_lit(DativeBondConstraintAst::ring_membership(RingScope::All, 1), false)]
    #[case::ring_membership_all_undetermined(DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined), true)]
    #[case::ring_membership_size_lit(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1), false)]
    #[case::ring_membership_size_undetermined(DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined), true)]
    fn test_dative_bond_constraint_ast_is_undetermined(#[case] c: DativeBondConstraintAst, #[case] expected: bool) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined), Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    #[case::same_key_incompatible(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false)), None)]
    #[case::different_key(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::All, 1), None)]
    fn test_dative_bond_constraint_ast_meet(#[case] a: DativeBondConstraintAst, #[case] b: DativeBondConstraintAst, #[case] expected: Option<DativeBondConstraintAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_widens(DativeBondConstraintAst::ring_membership(RingScope::All, 1), DativeBondConstraintAst::ring_membership(RingScope::All, 2), Ok(DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::lit_set([1, 2]))))]
    #[case::different_key(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::All, 1), Err(NoJoin))]
    fn test_dative_bond_constraint_ast_join(#[case] a: DativeBondConstraintAst, #[case] b: DativeBondConstraintAst, #[case] expected: Result<DativeBondConstraintAst, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), true)]
    #[case::same_key_incompatible(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false)), false)]
    #[case::different_key(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::All, 1), false)]
    fn test_dative_bond_constraint_ast_is_compatible(#[case] a: DativeBondConstraintAst, #[case] b: DativeBondConstraintAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_new() {
        let cs = DativeBondConstraintsAst::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_iter() {
        let cs = DativeBondConstraintsAst::from_iter([
            DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::All, 1),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraintAst::ring_membership(RingScope::All, 1),
                DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))], vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))])]
    #[case::overwrite_same_key(vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))], vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))])]
    #[case::vacuous_stores(vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined)], vec![DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined)])]
    #[case::new_key_sorts(vec![DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))], vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)])]
    #[case::ring_overwrite_scope(vec![DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 2)], vec![DativeBondConstraintAst::ring_membership(RingScope::Size(6), 2)])]
    fn test_dative_bond_constraints_ast_set(#[case] sequence: Vec<DativeBondConstraintAst>, #[case] expected: Vec<DativeBondConstraintAst>) {
        let mut cs = DativeBondConstraintsAst::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, DativeBondConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite_shared(
        vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))],
        vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false)), DativeBondConstraintAst::ring_membership(RingScope::All, 1)])]
    #[case::keeps_disjoint(
        vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))],
        vec![DativeBondConstraintAst::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::All, 1)])]
    #[case::vacuous_removes(
        vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined)],
        vec![DativeBondConstraintAst::ring_membership(RingScope::All, 1)])]
    fn test_dative_bond_constraints_ast_update(#[case] initial: Vec<DativeBondConstraintAst>, #[case] other: Vec<DativeBondConstraintAst>, #[case] expected: Vec<DativeBondConstraintAst>) {
        let mut cs = DativeBondConstraintsAst::from_iter(initial);
        cs.update(&DativeBondConstraintsAst::from_iter(other));
        assert_eq!(cs, DativeBondConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![DativeBondConstraintAst::ring_membership(RingScope::All, 1)], Some(DativeBondConstraintAst::ring_membership(RingScope::All, 1)), Some(DativeBondConstraintAst::ring_membership(RingScope::All, 2)), Ok(()), vec![DativeBondConstraintAst::ring_membership(RingScope::All, 2)])]
    #[case::remove(vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))], Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))), Ok(()), vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))])]
    #[case::old_mismatch(vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))], Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))), None, Err(Contradiction), vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))])]
    #[case::key_mismatch(vec![], Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))), Some(DativeBondConstraintAst::ring_membership(RingScope::All, 1)), Err(Contradiction), vec![])]
    fn test_dative_bond_constraints_ast_compare_and_set(
        #[case] initial: Vec<DativeBondConstraintAst>,
        #[case] old: Option<DativeBondConstraintAst>,
        #[case] new: Option<DativeBondConstraintAst>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<DativeBondConstraintAst>,
    ) {
        let mut cs = DativeBondConstraintsAst::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, DativeBondConstraintsAst::from_iter(expected_state));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(DativeBondConstraintKey::Aromatic, true)]
    #[case::ring_all_present(DativeBondConstraintKey::RingMembership(RingScope::All), true)]
    #[case::ring_size_present(DativeBondConstraintKey::RingMembership(RingScope::Size(6)), true)]
    #[case::ring_size_absent(DativeBondConstraintKey::RingMembership(RingScope::Size(5)), false)]
    fn test_dative_bond_constraints_ast_contains(
        #[case] key: DativeBondConstraintKey,
        #[case] expected: bool,
    ) {
        let cs = DativeBondConstraintsAst::from_iter([
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::All, 2),
            DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintKey::Aromatic, Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    #[case::ring_all(DativeBondConstraintKey::RingMembership(RingScope::All), Some(DativeBondConstraintAst::ring_membership(RingScope::All, 2)))]
    #[case::ring_size(DativeBondConstraintKey::RingMembership(RingScope::Size(6)), Some(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_size_absent(DativeBondConstraintKey::RingMembership(RingScope::Size(5)), None)]
    fn test_dative_bond_constraints_ast_get(
        #[case] key: DativeBondConstraintKey,
        #[case] expected: Option<DativeBondConstraintAst>,
    ) {
        let cs = DativeBondConstraintsAst::from_iter([
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::All, 2),
            DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(key), expected.as_ref());
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_remove() {
        let mut cs = DativeBondConstraintsAst::from_iter([
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::All, 2),
            DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        let removed = cs.remove(DativeBondConstraintKey::RingMembership(RingScope::Size(6)));
        assert_eq!(
            removed,
            Some(DativeBondConstraintAst::ring_membership(
                RingScope::Size(6),
                1
            )),
        );
        assert_eq!(
            cs,
            DativeBondConstraintsAst::from_iter([
                DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraintAst::ring_membership(RingScope::All, 2),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &DativeBondConstraintAst| matches!(c, DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) || matches!(c, DativeBondConstraintAst::RingMembership(m) if m.scope == RingScope::Size(6)), vec![
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)])]
    #[case::all_dropped(|_: &DativeBondConstraintAst| false, vec![])]
    fn test_dative_bond_constraints_ast_retain(
        #[case] predicate: impl FnMut(&DativeBondConstraintAst) -> bool,
        #[case] expected: Vec<DativeBondConstraintAst>,
    ) {
        let mut cs = DativeBondConstraintsAst::from_iter([
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::All, 1),
            DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        cs.retain(predicate);
        assert_eq!(cs, DativeBondConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_clear() {
        let mut cs = DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(
            BooleanAst::Lit(true),
        )]);
        cs.clear();
        assert_eq!(cs, DativeBondConstraintsAst::new());
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_take() {
        let mut cs = DativeBondConstraintsAst::from_iter([
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(
            drained,
            vec![
                DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ],
        );
        assert_eq!(cs, DativeBondConstraintsAst::new());
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_compact() {
        let cs = DativeBondConstraintsAst::from_iter([
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
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
        DativeBondConstraintsAst::from_iter([
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined),
        ]),
        Ok(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))])))]
    #[case::canonicalizes_values(
        DativeBondConstraintsAst::from_iter([
            DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::lit_set([2])),
        ]),
        Ok(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, 2)])))]
    fn test_dative_bond_constraints_ast_canonicalize(
        #[case] constraints: DativeBondConstraintsAst,
        #[case] expected: Result<DativeBondConstraintsAst, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys_kept(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, 1)]),
        Some(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::All, 1)])))]
    #[case::shared_key_meets(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined)]),
        Some(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))])))]
    #[case::shared_key_contradicts(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))]), None)]
    #[case::ring_size_unions(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::Size(5), 1)]), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)]),
        Some(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::Size(5), 1), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)])))]
    #[case::prunes_vacuous(DativeBondConstraintsAst::new(), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined)]), Some(DativeBondConstraintsAst::new()))]
    fn test_dative_bond_constraints_ast_meet(#[case] a: DativeBondConstraintsAst, #[case] b: DativeBondConstraintsAst, #[case] expected: Option<DativeBondConstraintsAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::keeps_only_shared_keys(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::All, 1)]), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]),
        DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]))]
    #[case::widens_value(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, 1)]), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, 2)]),
        DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::lit_set([1, 2]))]))]
    #[case::incompatible_drops_to_undetermined(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))]), DativeBondConstraintsAst::new())]
    fn test_dative_bond_constraints_ast_join(#[case] a: DativeBondConstraintsAst, #[case] b: DativeBondConstraintsAst, #[case] expected: DativeBondConstraintsAst) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(DativeBondConstraintsAst::new(), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::aromatic_required_present(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]),
        DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::aromatic_required_absent(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraintsAst::new(), false)]
    #[case::ring_membership_all_wildcard_matches_lit(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined)]),
        DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, 1)]), true)]
    #[case::ring_membership_all_lit_mismatch(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, 1)]),
        DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, 2)]), false)]
    #[case::ring_membership_size_subset(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::Size(5), 1)]),
        DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::Size(5), 1), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)]), true)]
    #[case::ring_membership_size_not_in_target(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::Size(7), 1)]),
        DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::Size(5), 1)]), false)]
    fn test_dative_bond_constraints_ast_matches(
        #[case] pattern: DativeBondConstraintsAst,
        #[case] target: DativeBondConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, 1)]), true)]
    #[case::shared_key_compatible(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), true)]
    #[case::shared_key_incompatible(DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))]), DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))]), false)]
    fn test_dative_bond_constraints_ast_is_compatible(#[case] a: DativeBondConstraintsAst, #[case] b: DativeBondConstraintsAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::All, 1)],
        vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::All, 1)])]
    #[case::unique_kind_last_wins(vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))],
        vec![DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))])]
    #[case::ring_appends(vec![DativeBondConstraintAst::ring_membership(RingScope::All, 1), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)],
        vec![DativeBondConstraintAst::ring_membership(RingScope::All, 1), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)])]
    #[case::empty(vec![], vec![])]
    fn test_dative_bond_constraints_ast_from_iter(
        #[case] input: Vec<DativeBondConstraintAst>,
        #[case] expected: Vec<DativeBondConstraintAst>,
    ) {
        let cs = DativeBondConstraintsAst::from_iter(input);
        assert_eq!(cs, DativeBondConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_into_iter() {
        let cs = DativeBondConstraintsAst::from_iter([
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_from_dative_bond_constraint() {
        let cs: DativeBondConstraintsAst =
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)).into();
        assert_eq!(
            cs,
            DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(
                BooleanAst::Lit(true)
            )]),
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_from_vec() {
        let cs: DativeBondConstraintsAst = vec![
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]
        .into();
        assert_eq!(
            cs,
            DativeBondConstraintsAst::from_iter([
                DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ]),
        );
    }
}
