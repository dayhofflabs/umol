//! Per-multicenter-bond constraints.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use strum::EnumDiscriminants;

use super::super::error::{Contradiction, NoJoin};
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Multicenter-bond-scope constraint. Held inline on `MulticenterBondAst` via
/// `MulticenterBondConstraints`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, EnumDiscriminants)]
#[strum_discriminants(name(MulticenterBondConstraintKind), derive(Hash))]
pub enum MulticenterBondConstraint {
    /// Asserted total electron count for the multicenter bond. Cross-checked
    /// by the `ConsistencyValidator` against `sum(MulticenterBondAst::electrons)`.
    ElectronCount(ValueAst),
}

impl MulticenterBondConstraint {
    pub fn electron_count(v: impl Into<ValueAst>) -> Self {
        Self::ElectronCount(v.into())
    }

    pub fn kind(&self) -> MulticenterBondConstraintKind {
        self.into()
    }

    /// Multicenter-bond constraint key, unique within a `MulticenterBondConstraints` container.
    pub fn key(&self) -> MulticenterBondConstraintKey {
        match self {
            Self::ElectronCount(_) => MulticenterBondConstraintKey::ElectronCount,
        }
    }

    /// Vacuous form of constraint key, used for removal.
    pub fn as_undetermined(&self) -> Self {
        match self {
            Self::ElectronCount(_) => Self::ElectronCount(ValueAst::Undetermined),
        }
    }

    pub fn compact(self, _compaction: &IdCompaction) -> Option<Self> {
        // Value-only: no indices to compact.
        Some(self)
    }

    /// Value-only: no indices to remap.
    pub fn remap(self, _map: &IdRemapping) -> Self {
        self
    }
}

impl Canonicalize for MulticenterBondConstraint {
    /// Canonicalize the inner value; the kind is preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::ElectronCount(v) => Self::ElectronCount(v.canonicalize()?),
        })
    }
}

impl Lattice for MulticenterBondConstraint {
    fn is_undetermined(&self) -> bool {
        match self {
            Self::ElectronCount(v) => v.is_undetermined(),
        }
    }

    fn is_ground(&self) -> bool {
        match self {
            Self::ElectronCount(v) => v.is_ground(),
        }
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::ElectronCount(a), Self::ElectronCount(b)) => a.meet(b).map(Self::ElectronCount),
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        match (self, other) {
            (Self::ElectronCount(a), Self::ElectronCount(b)) => Ok(Self::ElectronCount(a.join(b)?)),
        }
    }

    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::ElectronCount(a), Self::ElectronCount(b)) => a.matches(b),
        }
    }

    fn is_compatible(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::ElectronCount(a), Self::ElectronCount(b)) => a.is_compatible(b),
        }
    }
}

/// Entry identity: discriminant only (every kind is single-valued, no sub-key).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MulticenterBondConstraintKey {
    ElectronCount,
}

/// Multicenter bond constraints container, ordered, unique by key, sorted flat vector storage.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterBondConstraints(Vec<MulticenterBondConstraint>);

impl MulticenterBondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The bond's electron count, or `Undetermined` when no `ElectronCount` constraint is present.
    pub fn electron_count(&self) -> ValueAst {
        match self.get(MulticenterBondConstraintKey::ElectronCount) {
            Some(MulticenterBondConstraint::ElectronCount(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    fn find(&self, key: MulticenterBondConstraintKey) -> Result<usize, usize> {
        self.0.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains(&self, key: MulticenterBondConstraintKey) -> bool {
        self.find(key).is_ok()
    }

    pub fn get(&self, key: MulticenterBondConstraintKey) -> Option<&MulticenterBondConstraint> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: MulticenterBondConstraint) {
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
        old: Option<MulticenterBondConstraint>,
        new: Option<MulticenterBondConstraint>,
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

    pub fn remove(
        &mut self,
        key: MulticenterBondConstraintKey,
    ) -> Option<MulticenterBondConstraint> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = MulticenterBondConstraint>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &MulticenterBondConstraints) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&MulticenterBondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = MulticenterBondConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, MulticenterBondConstraint> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &IdCompaction) -> Self {
        self
    }
}

impl Canonicalize for MulticenterBondConstraints {
    /// Canonicalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Canonicalize::canonicalize)
            .collect::<Result<Vec<MulticenterBondConstraint>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl Lattice for MulticenterBondConstraints {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`MulticenterBondConstraint::meet`; a `None` aborts the whole meet), an
    /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<MulticenterBondConstraint> = Vec::new();
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
    /// (`MulticenterBondConstraint::join`); a single-side key widens to the absent ⊤ and is dropped.
    /// The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<MulticenterBondConstraint> = Vec::new();
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

    /// Pattern-driven: the electron-count value is matched on its own lattice; an empty
    /// pattern matches any target.
    fn matches(&self, target: &Self) -> bool {
        self.electron_count().matches(&target.electron_count())
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

impl FromIterator<MulticenterBondConstraint> for MulticenterBondConstraints {
    fn from_iter<I: IntoIterator<Item = MulticenterBondConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for MulticenterBondConstraints {
    type Item = MulticenterBondConstraint;
    type IntoIter = IntoIter<MulticenterBondConstraint>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<MulticenterBondConstraint> for MulticenterBondConstraints {
    fn from(c: MulticenterBondConstraint) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<MulticenterBondConstraint>> for MulticenterBondConstraints {
    fn from(cs: Vec<MulticenterBondConstraint>) -> Self {
        Self::from_iter(cs)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Compaction;

    use super::*;

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraint::electron_count(6),
        MulticenterBondConstraint::ElectronCount(ValueAst::Lit(6))
    )]
    fn test_multicenter_bond_constraint_constructors(
        #[case] actual: MulticenterBondConstraint,
        #[case] expected: MulticenterBondConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraint::electron_count(6),
        MulticenterBondConstraintKind::ElectronCount
    )]
    fn test_multicenter_bond_constraint_kind(
        #[case] c: MulticenterBondConstraint,
        #[case] expected: MulticenterBondConstraintKind,
    ) {
        assert_eq!(c.kind(), expected);
    }

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraint::electron_count(6),
        MulticenterBondConstraintKey::ElectronCount
    )]
    fn test_multicenter_bond_constraint_key(
        #[case] c: MulticenterBondConstraint,
        #[case] expected: MulticenterBondConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraint::electron_count(6),
        MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined)
    )]
    fn test_multicenter_bond_constraint_as_undetermined(
        #[case] c: MulticenterBondConstraint,
        #[case] expected: MulticenterBondConstraint,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count_litset_singleton(MulticenterBondConstraint::ElectronCount(ValueAst::lit_set([6])), Ok(MulticenterBondConstraint::electron_count(6)))]
    #[case::empty_litset_contradiction(MulticenterBondConstraint::ElectronCount(ValueAst::lit_set(Vec::<i64>::new())), Err(Contradiction))]
    fn test_multicenter_bond_constraint_canonicalize(
        #[case] constraint: MulticenterBondConstraint,
        #[case] expected: Result<MulticenterBondConstraint, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
    }

    #[rstest]
    #[case::lit(MulticenterBondConstraint::electron_count(6), false)]
    #[case::undetermined(MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined), true)]
    fn test_multicenter_bond_constraint_is_undetermined(
        #[case] c: MulticenterBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(MulticenterBondConstraint::electron_count(6), MulticenterBondConstraint::electron_count(6), Some(MulticenterBondConstraint::electron_count(6)))]
    #[case::narrows_undetermined(MulticenterBondConstraint::electron_count(6), MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined), Some(MulticenterBondConstraint::electron_count(6)))]
    #[case::incompatible(MulticenterBondConstraint::electron_count(6), MulticenterBondConstraint::electron_count(2), None)]
    fn test_multicenter_bond_constraint_meet(#[case] a: MulticenterBondConstraint, #[case] b: MulticenterBondConstraint, #[case] expected: Option<MulticenterBondConstraint>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(MulticenterBondConstraint::electron_count(6), MulticenterBondConstraint::electron_count(6), Ok(MulticenterBondConstraint::electron_count(6)))]
    #[case::widens(MulticenterBondConstraint::electron_count(1), MulticenterBondConstraint::electron_count(2), Ok(MulticenterBondConstraint::ElectronCount(ValueAst::lit_set([1, 2]))))]
    fn test_multicenter_bond_constraint_join(#[case] a: MulticenterBondConstraint, #[case] b: MulticenterBondConstraint, #[case] expected: Result<MulticenterBondConstraint, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(MulticenterBondConstraint::electron_count(6), MulticenterBondConstraint::electron_count(6), true)]
    #[case::incompatible(MulticenterBondConstraint::electron_count(6), MulticenterBondConstraint::electron_count(2), false)]
    fn test_multicenter_bond_constraint_is_compatible(#[case] a: MulticenterBondConstraint, #[case] b: MulticenterBondConstraint, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_new() {
        let cs = MulticenterBondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_iter() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![MulticenterBondConstraint::electron_count(6)]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![MulticenterBondConstraint::electron_count(6)], vec![MulticenterBondConstraint::electron_count(6)])]
    #[case::overwrite_same_key(vec![MulticenterBondConstraint::electron_count(6), MulticenterBondConstraint::electron_count(10)], vec![MulticenterBondConstraint::electron_count(10)])]
    #[case::vacuous_stores(vec![MulticenterBondConstraint::electron_count(6), MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined)], vec![MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined)])]
    fn test_multicenter_bond_constraints_set(#[case] sequence: Vec<MulticenterBondConstraint>, #[case] expected: Vec<MulticenterBondConstraint>) {
        let mut cs = MulticenterBondConstraints::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, MulticenterBondConstraints::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite(
        vec![MulticenterBondConstraint::electron_count(6)],
        vec![MulticenterBondConstraint::electron_count(10)],
        vec![MulticenterBondConstraint::electron_count(10)])]
    #[case::adds_from_empty(
        vec![],
        vec![MulticenterBondConstraint::electron_count(6)],
        vec![MulticenterBondConstraint::electron_count(6)])]
    #[case::vacuous_removes(
        vec![MulticenterBondConstraint::electron_count(6)],
        vec![MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined)],
        vec![])]
    fn test_multicenter_bond_constraints_update(#[case] initial: Vec<MulticenterBondConstraint>, #[case] other: Vec<MulticenterBondConstraint>, #[case] expected: Vec<MulticenterBondConstraint>) {
        let mut cs = MulticenterBondConstraints::from_iter(initial);
        cs.update(&MulticenterBondConstraints::from_iter(other));
        assert_eq!(cs, MulticenterBondConstraints::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![MulticenterBondConstraint::electron_count(6)], Some(MulticenterBondConstraint::electron_count(6)), Some(MulticenterBondConstraint::electron_count(10)), Ok(()), vec![MulticenterBondConstraint::electron_count(10)])]
    #[case::remove(vec![MulticenterBondConstraint::electron_count(6)], Some(MulticenterBondConstraint::electron_count(6)), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(MulticenterBondConstraint::electron_count(6)), Ok(()), vec![MulticenterBondConstraint::electron_count(6)])]
    #[case::old_mismatch(vec![MulticenterBondConstraint::electron_count(6)], Some(MulticenterBondConstraint::electron_count(2)), None, Err(Contradiction), vec![MulticenterBondConstraint::electron_count(6)])]
    fn test_multicenter_bond_constraints_compare_and_set(
        #[case] initial: Vec<MulticenterBondConstraint>,
        #[case] old: Option<MulticenterBondConstraint>,
        #[case] new: Option<MulticenterBondConstraint>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<MulticenterBondConstraint>,
    ) {
        let mut cs = MulticenterBondConstraints::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, MulticenterBondConstraints::from_iter(expected_state));
    }

    #[rstest]
    #[case::present(
        MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)),
        true
    )]
    #[case::absent(MulticenterBondConstraints::new(), false)]
    fn test_multicenter_bond_constraints_contains(
        #[case] cs: MulticenterBondConstraints,
        #[case] expected: bool,
    ) {
        assert_eq!(
            cs.contains(MulticenterBondConstraintKey::ElectronCount),
            expected
        );
    }

    #[rstest]
    #[case::present(
        MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)),
        Some(MulticenterBondConstraint::electron_count(6))
    )]
    #[case::absent(MulticenterBondConstraints::new(), None)]
    fn test_multicenter_bond_constraints_get(
        #[case] cs: MulticenterBondConstraints,
        #[case] expected: Option<MulticenterBondConstraint>,
    ) {
        assert_eq!(
            cs.get(MulticenterBondConstraintKey::ElectronCount),
            expected.as_ref()
        );
    }

    #[rstest]
    #[case::present(
        MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)),
        Some(MulticenterBondConstraint::electron_count(6)),
        MulticenterBondConstraints::new()
    )]
    #[case::absent(
        MulticenterBondConstraints::new(),
        None,
        MulticenterBondConstraints::new()
    )]
    fn test_multicenter_bond_constraints_remove(
        #[case] mut cs: MulticenterBondConstraints,
        #[case] expected_removed: Option<MulticenterBondConstraint>,
        #[case] expected_state: MulticenterBondConstraints,
    ) {
        assert_eq!(
            cs.remove(MulticenterBondConstraintKey::ElectronCount),
            expected_removed
        );
        assert_eq!(cs, expected_state);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &MulticenterBondConstraint| matches!(c, MulticenterBondConstraint::ElectronCount(_)), vec![MulticenterBondConstraint::electron_count(6)])]
    #[case::all_dropped(|_: &MulticenterBondConstraint| false, vec![])]
    fn test_multicenter_bond_constraints_retain(
        #[case] predicate: impl FnMut(&MulticenterBondConstraint) -> bool,
        #[case] expected: Vec<MulticenterBondConstraint>,
    ) {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6));
        cs.retain(predicate);
        assert_eq!(cs, MulticenterBondConstraints::from_iter(expected));
    }

    #[rstest]
    fn test_multicenter_bond_constraints_clear() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6));
        cs.clear();
        assert_eq!(cs, MulticenterBondConstraints::new());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_take() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6));
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(drained, vec![MulticenterBondConstraint::electron_count(6)]);
        assert_eq!(cs, MulticenterBondConstraints::new());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_compact() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6));
        let compaction = IdCompaction::new(
            Compaction::new(vec![0, 1], vec![0]),
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
    #[case::canonicalizes_value(
        MulticenterBondConstraints::from(MulticenterBondConstraint::ElectronCount(ValueAst::lit_set([6]))),
        Ok(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6))))]
    #[case::drop_vacuous(
        MulticenterBondConstraints::from(MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined)),
        Ok(MulticenterBondConstraints::new()))]
    #[case::contradiction(
        MulticenterBondConstraints::from(MulticenterBondConstraint::ElectronCount(ValueAst::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_multicenter_bond_constraints_canonicalize(
        #[case] constraints: MulticenterBondConstraints,
        #[case] expected: Result<MulticenterBondConstraints, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::a_only_kept(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), MulticenterBondConstraints::new(),
        Some(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6))))]
    #[case::b_only_kept(MulticenterBondConstraints::new(), MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)),
        Some(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6))))]
    #[case::shared_key_meets(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), MulticenterBondConstraints::from(MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined)),
        Some(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6))))]
    #[case::shared_key_contradicts(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2)), None)]
    #[case::prunes_vacuous(MulticenterBondConstraints::new(), MulticenterBondConstraints::from(MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined)), Some(MulticenterBondConstraints::new()))]
    fn test_multicenter_bond_constraints_meet(#[case] a: MulticenterBondConstraints, #[case] b: MulticenterBondConstraints, #[case] expected: Option<MulticenterBondConstraints>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::widens_value(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(1)), MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2)),
        MulticenterBondConstraints::from(MulticenterBondConstraint::ElectronCount(ValueAst::lit_set([1, 2]))))]
    #[case::single_side_dropped(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), MulticenterBondConstraints::new(),
        MulticenterBondConstraints::new())]
    #[case::undetermined_drops(MulticenterBondConstraints::from(MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined)), MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)),
        MulticenterBondConstraints::new())]
    fn test_multicenter_bond_constraints_join(#[case] a: MulticenterBondConstraints, #[case] b: MulticenterBondConstraints, #[case] expected: MulticenterBondConstraints) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(MulticenterBondConstraints::new(), MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), true)]
    #[case::required_present(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), true)]
    #[case::required_absent(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), MulticenterBondConstraints::new(), false)]
    #[case::wildcard_matches_lit(MulticenterBondConstraints::from(MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined)), MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), true)]
    #[case::lit_mismatch(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2)), false)]
    fn test_multicenter_bond_constraints_matches(
        #[case] pattern: MulticenterBondConstraints,
        #[case] target: MulticenterBondConstraints,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_one_empty(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), MulticenterBondConstraints::new(), true)]
    #[case::shared_key_compatible(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), true)]
    #[case::shared_key_incompatible(MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6)), MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2)), false)]
    fn test_multicenter_bond_constraints_is_compatible(#[case] a: MulticenterBondConstraints, #[case] b: MulticenterBondConstraints, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(vec![MulticenterBondConstraint::electron_count(6)], vec![MulticenterBondConstraint::electron_count(6)])]
    #[case::same_key_last_wins(vec![MulticenterBondConstraint::electron_count(2), MulticenterBondConstraint::electron_count(6)], vec![MulticenterBondConstraint::electron_count(6)])]
    #[case::empty(vec![], vec![])]
    fn test_multicenter_bond_constraints_from_iter(
        #[case] input: Vec<MulticenterBondConstraint>,
        #[case] expected: Vec<MulticenterBondConstraint>,
    ) {
        let cs = MulticenterBondConstraints::from_iter(input);
        assert_eq!(cs, MulticenterBondConstraints::from_iter(expected));
    }

    #[rstest]
    fn test_multicenter_bond_constraints_into_iter() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(6));
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![MulticenterBondConstraint::electron_count(6)]
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_from_multicenter_bond_constraint() {
        let cs: MulticenterBondConstraints = MulticenterBondConstraint::electron_count(6).into();
        assert_eq!(
            cs,
            MulticenterBondConstraints::from_iter([MulticenterBondConstraint::electron_count(6)]),
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_from_vec() {
        let cs: MulticenterBondConstraints =
            vec![MulticenterBondConstraint::electron_count(6)].into();
        assert_eq!(
            cs,
            MulticenterBondConstraints::from_iter([MulticenterBondConstraint::electron_count(6)]),
        );
    }
}
