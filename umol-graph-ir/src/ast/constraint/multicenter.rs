//! Per-multicenter-bond constraints.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use super::super::error::{Contradiction, NoJoin};
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Multicenter-bond-scope constraint. Held inline on `MulticenterBondAst` via
/// `MulticenterBondConstraintsAst`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MulticenterBondConstraintAst {
    /// Asserted total electron count for the multicenter bond. Cross-checked
    /// by the `ConsistencyValidator` against `sum(MulticenterBondAst::electrons)`.
    ElectronCount(ValueAst),
}

impl MulticenterBondConstraintAst {
    pub fn electron_count(v: impl Into<ValueAst>) -> Self {
        Self::ElectronCount(v.into())
    }

    /// Multicenter-bond constraint key, unique within a `MulticenterBondConstraintsAst` container.
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

impl Canonicalize for MulticenterBondConstraintAst {
    /// Canonicalize the inner value; the kind is preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::ElectronCount(v) => Self::ElectronCount(v.canonicalize()?),
        })
    }
}

impl Lattice for MulticenterBondConstraintAst {
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
pub struct MulticenterBondConstraintsAst(Vec<MulticenterBondConstraintAst>);

impl MulticenterBondConstraintsAst {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The bond's electron count, or `Undetermined` when no `ElectronCount` constraint is present.
    pub fn electron_count(&self) -> ValueAst {
        match self.get(MulticenterBondConstraintKey::ElectronCount) {
            Some(MulticenterBondConstraintAst::ElectronCount(v)) => v.clone(),
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

    pub fn get(&self, key: MulticenterBondConstraintKey) -> Option<&MulticenterBondConstraintAst> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: MulticenterBondConstraintAst) {
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
        old: Option<MulticenterBondConstraintAst>,
        new: Option<MulticenterBondConstraintAst>,
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
    ) -> Option<MulticenterBondConstraintAst> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = MulticenterBondConstraintAst>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &MulticenterBondConstraintsAst) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&MulticenterBondConstraintAst) -> bool) {
        self.0.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl ExactSizeIterator<Item = MulticenterBondConstraintAst> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, MulticenterBondConstraintAst> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &IdCompaction) -> Self {
        self
    }
}

impl Canonicalize for MulticenterBondConstraintsAst {
    /// Canonicalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Canonicalize::canonicalize)
            .collect::<Result<Vec<MulticenterBondConstraintAst>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl Lattice for MulticenterBondConstraintsAst {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`MulticenterBondConstraintAst::meet`; a `None` aborts the whole meet), an
    /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<MulticenterBondConstraintAst> = Vec::new();
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
    /// (`MulticenterBondConstraintAst::join`); a single-side key widens to the absent ⊤ and is dropped.
    /// The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<MulticenterBondConstraintAst> = Vec::new();
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

impl FromIterator<MulticenterBondConstraintAst> for MulticenterBondConstraintsAst {
    fn from_iter<I: IntoIterator<Item = MulticenterBondConstraintAst>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for MulticenterBondConstraintsAst {
    type Item = MulticenterBondConstraintAst;
    type IntoIter = IntoIter<MulticenterBondConstraintAst>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<MulticenterBondConstraintAst> for MulticenterBondConstraintsAst {
    fn from(c: MulticenterBondConstraintAst) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<MulticenterBondConstraintAst>> for MulticenterBondConstraintsAst {
    fn from(cs: Vec<MulticenterBondConstraintAst>) -> Self {
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
        MulticenterBondConstraintAst::electron_count(6),
        MulticenterBondConstraintAst::ElectronCount(ValueAst::Lit(6))
    )]
    fn test_multicenter_bond_constraint_ast_constructors(
        #[case] actual: MulticenterBondConstraintAst,
        #[case] expected: MulticenterBondConstraintAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraintAst::electron_count(6),
        MulticenterBondConstraintKey::ElectronCount
    )]
    fn test_multicenter_bond_constraint_ast_key(
        #[case] c: MulticenterBondConstraintAst,
        #[case] expected: MulticenterBondConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraintAst::electron_count(6),
        MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined)
    )]
    fn test_multicenter_bond_constraint_ast_as_undetermined(
        #[case] c: MulticenterBondConstraintAst,
        #[case] expected: MulticenterBondConstraintAst,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count_litset_singleton(MulticenterBondConstraintAst::ElectronCount(ValueAst::lit_set([6])), Ok(MulticenterBondConstraintAst::electron_count(6)))]
    #[case::empty_litset_contradiction(MulticenterBondConstraintAst::ElectronCount(ValueAst::lit_set(Vec::<i64>::new())), Err(Contradiction))]
    fn test_multicenter_bond_constraint_ast_canonicalize(
        #[case] constraint: MulticenterBondConstraintAst,
        #[case] expected: Result<MulticenterBondConstraintAst, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
    }

    #[rstest]
    #[case::lit(MulticenterBondConstraintAst::electron_count(6), false)]
    #[case::undetermined(
        MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined),
        true
    )]
    fn test_multicenter_bond_constraint_ast_is_undetermined(
        #[case] c: MulticenterBondConstraintAst,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(MulticenterBondConstraintAst::electron_count(6), MulticenterBondConstraintAst::electron_count(6), Some(MulticenterBondConstraintAst::electron_count(6)))]
    #[case::narrows_undetermined(MulticenterBondConstraintAst::electron_count(6), MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined), Some(MulticenterBondConstraintAst::electron_count(6)))]
    #[case::incompatible(MulticenterBondConstraintAst::electron_count(6), MulticenterBondConstraintAst::electron_count(2), None)]
    fn test_multicenter_bond_constraint_ast_meet(#[case] a: MulticenterBondConstraintAst, #[case] b: MulticenterBondConstraintAst, #[case] expected: Option<MulticenterBondConstraintAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(MulticenterBondConstraintAst::electron_count(6), MulticenterBondConstraintAst::electron_count(6), Ok(MulticenterBondConstraintAst::electron_count(6)))]
    #[case::widens(MulticenterBondConstraintAst::electron_count(1), MulticenterBondConstraintAst::electron_count(2), Ok(MulticenterBondConstraintAst::ElectronCount(ValueAst::lit_set([1, 2]))))]
    fn test_multicenter_bond_constraint_ast_join(#[case] a: MulticenterBondConstraintAst, #[case] b: MulticenterBondConstraintAst, #[case] expected: Result<MulticenterBondConstraintAst, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(MulticenterBondConstraintAst::electron_count(6), MulticenterBondConstraintAst::electron_count(6), true)]
    #[case::incompatible(MulticenterBondConstraintAst::electron_count(6), MulticenterBondConstraintAst::electron_count(2), false)]
    fn test_multicenter_bond_constraint_ast_is_compatible(#[case] a: MulticenterBondConstraintAst, #[case] b: MulticenterBondConstraintAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_new() {
        let cs = MulticenterBondConstraintsAst::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_iter() {
        let cs =
            MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![MulticenterBondConstraintAst::electron_count(6)]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![MulticenterBondConstraintAst::electron_count(6)], vec![MulticenterBondConstraintAst::electron_count(6)])]
    #[case::overwrite_same_key(vec![MulticenterBondConstraintAst::electron_count(6), MulticenterBondConstraintAst::electron_count(10)], vec![MulticenterBondConstraintAst::electron_count(10)])]
    #[case::vacuous_stores(vec![MulticenterBondConstraintAst::electron_count(6), MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined)], vec![MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined)])]
    fn test_multicenter_bond_constraints_ast_set(#[case] sequence: Vec<MulticenterBondConstraintAst>, #[case] expected: Vec<MulticenterBondConstraintAst>) {
        let mut cs = MulticenterBondConstraintsAst::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, MulticenterBondConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite(
        vec![MulticenterBondConstraintAst::electron_count(6)],
        vec![MulticenterBondConstraintAst::electron_count(10)],
        vec![MulticenterBondConstraintAst::electron_count(10)])]
    #[case::adds_from_empty(
        vec![],
        vec![MulticenterBondConstraintAst::electron_count(6)],
        vec![MulticenterBondConstraintAst::electron_count(6)])]
    #[case::vacuous_removes(
        vec![MulticenterBondConstraintAst::electron_count(6)],
        vec![MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined)],
        vec![])]
    fn test_multicenter_bond_constraints_ast_update(#[case] initial: Vec<MulticenterBondConstraintAst>, #[case] other: Vec<MulticenterBondConstraintAst>, #[case] expected: Vec<MulticenterBondConstraintAst>) {
        let mut cs = MulticenterBondConstraintsAst::from_iter(initial);
        cs.update(&MulticenterBondConstraintsAst::from_iter(other));
        assert_eq!(cs, MulticenterBondConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![MulticenterBondConstraintAst::electron_count(6)], Some(MulticenterBondConstraintAst::electron_count(6)), Some(MulticenterBondConstraintAst::electron_count(10)), Ok(()), vec![MulticenterBondConstraintAst::electron_count(10)])]
    #[case::remove(vec![MulticenterBondConstraintAst::electron_count(6)], Some(MulticenterBondConstraintAst::electron_count(6)), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(MulticenterBondConstraintAst::electron_count(6)), Ok(()), vec![MulticenterBondConstraintAst::electron_count(6)])]
    #[case::old_mismatch(vec![MulticenterBondConstraintAst::electron_count(6)], Some(MulticenterBondConstraintAst::electron_count(2)), None, Err(Contradiction), vec![MulticenterBondConstraintAst::electron_count(6)])]
    fn test_multicenter_bond_constraints_ast_compare_and_set(
        #[case] initial: Vec<MulticenterBondConstraintAst>,
        #[case] old: Option<MulticenterBondConstraintAst>,
        #[case] new: Option<MulticenterBondConstraintAst>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<MulticenterBondConstraintAst>,
    ) {
        let mut cs = MulticenterBondConstraintsAst::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, MulticenterBondConstraintsAst::from_iter(expected_state));
    }

    #[rstest]
    #[case::present(
        MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)),
        true
    )]
    #[case::absent(MulticenterBondConstraintsAst::new(), false)]
    fn test_multicenter_bond_constraints_ast_contains(
        #[case] cs: MulticenterBondConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(
            cs.contains(MulticenterBondConstraintKey::ElectronCount),
            expected
        );
    }

    #[rstest]
    #[case::present(
        MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)),
        Some(MulticenterBondConstraintAst::electron_count(6))
    )]
    #[case::absent(MulticenterBondConstraintsAst::new(), None)]
    fn test_multicenter_bond_constraints_ast_get(
        #[case] cs: MulticenterBondConstraintsAst,
        #[case] expected: Option<MulticenterBondConstraintAst>,
    ) {
        assert_eq!(
            cs.get(MulticenterBondConstraintKey::ElectronCount),
            expected.as_ref()
        );
    }

    #[rstest]
    #[case::present(
        MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)),
        Some(MulticenterBondConstraintAst::electron_count(6)),
        MulticenterBondConstraintsAst::new()
    )]
    #[case::absent(
        MulticenterBondConstraintsAst::new(),
        None,
        MulticenterBondConstraintsAst::new()
    )]
    fn test_multicenter_bond_constraints_ast_remove(
        #[case] mut cs: MulticenterBondConstraintsAst,
        #[case] expected_removed: Option<MulticenterBondConstraintAst>,
        #[case] expected_state: MulticenterBondConstraintsAst,
    ) {
        assert_eq!(
            cs.remove(MulticenterBondConstraintKey::ElectronCount),
            expected_removed
        );
        assert_eq!(cs, expected_state);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &MulticenterBondConstraintAst| matches!(c, MulticenterBondConstraintAst::ElectronCount(_)), vec![MulticenterBondConstraintAst::electron_count(6)])]
    #[case::all_dropped(|_: &MulticenterBondConstraintAst| false, vec![])]
    fn test_multicenter_bond_constraints_ast_retain(
        #[case] predicate: impl FnMut(&MulticenterBondConstraintAst) -> bool,
        #[case] expected: Vec<MulticenterBondConstraintAst>,
    ) {
        let mut cs = MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6));
        cs.retain(predicate);
        assert_eq!(cs, MulticenterBondConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_clear() {
        let mut cs =
            MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6));
        cs.clear();
        assert_eq!(cs, MulticenterBondConstraintsAst::new());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_take() {
        let mut empty = MulticenterBondConstraintsAst::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let mut cs =
            MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6));
        let mut taken = cs.take();
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(
            taken.next(),
            Some(MulticenterBondConstraintAst::electron_count(6)),
        );
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.size_hint(), (0, Some(0)));
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(cs, MulticenterBondConstraintsAst::new());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_compact() {
        let cs =
            MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6));
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
        MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::ElectronCount(ValueAst::lit_set([6]))),
        Ok(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6))))]
    #[case::drop_vacuous(
        MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined)),
        Ok(MulticenterBondConstraintsAst::new()))]
    #[case::contradiction(
        MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::ElectronCount(ValueAst::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_multicenter_bond_constraints_ast_canonicalize(
        #[case] constraints: MulticenterBondConstraintsAst,
        #[case] expected: Result<MulticenterBondConstraintsAst, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::a_only_kept(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), MulticenterBondConstraintsAst::new(),
        Some(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6))))]
    #[case::b_only_kept(MulticenterBondConstraintsAst::new(), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)),
        Some(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6))))]
    #[case::shared_key_meets(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined)),
        Some(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6))))]
    #[case::shared_key_contradicts(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(2)), None)]
    #[case::prunes_vacuous(MulticenterBondConstraintsAst::new(), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined)), Some(MulticenterBondConstraintsAst::new()))]
    fn test_multicenter_bond_constraints_ast_meet(#[case] a: MulticenterBondConstraintsAst, #[case] b: MulticenterBondConstraintsAst, #[case] expected: Option<MulticenterBondConstraintsAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::widens_value(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(1)), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(2)),
        MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::ElectronCount(ValueAst::lit_set([1, 2]))))]
    #[case::single_side_dropped(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), MulticenterBondConstraintsAst::new(),
        MulticenterBondConstraintsAst::new())]
    #[case::undetermined_drops(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined)), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)),
        MulticenterBondConstraintsAst::new())]
    fn test_multicenter_bond_constraints_ast_join(#[case] a: MulticenterBondConstraintsAst, #[case] b: MulticenterBondConstraintsAst, #[case] expected: MulticenterBondConstraintsAst) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(MulticenterBondConstraintsAst::new(), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), true)]
    #[case::required_present(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), true)]
    #[case::required_absent(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), MulticenterBondConstraintsAst::new(), false)]
    #[case::wildcard_matches_lit(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined)), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), true)]
    #[case::lit_mismatch(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(2)), false)]
    fn test_multicenter_bond_constraints_ast_matches(
        #[case] pattern: MulticenterBondConstraintsAst,
        #[case] target: MulticenterBondConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_one_empty(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), MulticenterBondConstraintsAst::new(), true)]
    #[case::shared_key_compatible(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), true)]
    #[case::shared_key_incompatible(MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6)), MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(2)), false)]
    fn test_multicenter_bond_constraints_ast_is_compatible(#[case] a: MulticenterBondConstraintsAst, #[case] b: MulticenterBondConstraintsAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(vec![MulticenterBondConstraintAst::electron_count(6)], vec![MulticenterBondConstraintAst::electron_count(6)])]
    #[case::same_key_last_wins(vec![MulticenterBondConstraintAst::electron_count(2), MulticenterBondConstraintAst::electron_count(6)], vec![MulticenterBondConstraintAst::electron_count(6)])]
    #[case::empty(vec![], vec![])]
    fn test_multicenter_bond_constraints_ast_from_iter(
        #[case] input: Vec<MulticenterBondConstraintAst>,
        #[case] expected: Vec<MulticenterBondConstraintAst>,
    ) {
        let cs = MulticenterBondConstraintsAst::from_iter(input);
        assert_eq!(cs, MulticenterBondConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_into_iter() {
        let cs =
            MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(6));
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![MulticenterBondConstraintAst::electron_count(6)]
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_from_multicenter_bond_constraint() {
        let cs: MulticenterBondConstraintsAst =
            MulticenterBondConstraintAst::electron_count(6).into();
        assert_eq!(
            cs,
            MulticenterBondConstraintsAst::from_iter([
                MulticenterBondConstraintAst::electron_count(6)
            ]),
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_from_vec() {
        let cs: MulticenterBondConstraintsAst =
            vec![MulticenterBondConstraintAst::electron_count(6)].into();
        assert_eq!(
            cs,
            MulticenterBondConstraintsAst::from_iter([
                MulticenterBondConstraintAst::electron_count(6)
            ]),
        );
    }
}
