//! Per-aromatic-system constraints.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use super::super::error::{Contradiction, NoJoin};
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Aromatic-system-scope constraint. Held inline on `AromaticSystemAst` via
/// `AromaticSystemConstraintsAst`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticSystemConstraintAst {
    /// Asserted total π-electron count for the system. Cross-checked by the
    /// `ConsistencyValidator` against `sum(AromaticSystemAst::electrons)`.
    ElectronCount(ValueAst),
}

impl AromaticSystemConstraintAst {
    pub fn electron_count(v: impl Into<ValueAst>) -> Self {
        Self::ElectronCount(v.into())
    }

    /// Aromatic-system constraint key, unique within an `AromaticSystemConstraintsAst` container.
    pub fn key(&self) -> AromaticSystemConstraintKey {
        match self {
            Self::ElectronCount(_) => AromaticSystemConstraintKey::ElectronCount,
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

impl Canonicalize for AromaticSystemConstraintAst {
    /// Canonicalize the inner value; the kind is preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::ElectronCount(v) => Self::ElectronCount(v.canonicalize()?),
        })
    }
}

impl Lattice for AromaticSystemConstraintAst {
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
pub enum AromaticSystemConstraintKey {
    ElectronCount,
}

/// Aromatic system constraints container, ordered, unique by key, sorted flat vector storage.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AromaticSystemConstraintsAst(Vec<AromaticSystemConstraintAst>);

impl AromaticSystemConstraintsAst {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The system's electron count, or `Undetermined` when no `ElectronCount` constraint is present.
    pub fn electron_count(&self) -> ValueAst {
        match self.get(AromaticSystemConstraintKey::ElectronCount) {
            Some(AromaticSystemConstraintAst::ElectronCount(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    fn find(&self, key: AromaticSystemConstraintKey) -> Result<usize, usize> {
        self.0.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains(&self, key: AromaticSystemConstraintKey) -> bool {
        self.find(key).is_ok()
    }

    pub fn get(&self, key: AromaticSystemConstraintKey) -> Option<&AromaticSystemConstraintAst> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: AromaticSystemConstraintAst) {
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
        old: Option<AromaticSystemConstraintAst>,
        new: Option<AromaticSystemConstraintAst>,
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
        key: AromaticSystemConstraintKey,
    ) -> Option<AromaticSystemConstraintAst> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = AromaticSystemConstraintAst>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &AromaticSystemConstraintsAst) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&AromaticSystemConstraintAst) -> bool) {
        self.0.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = AromaticSystemConstraintAst> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, AromaticSystemConstraintAst> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &IdCompaction) -> Self {
        self
    }
}

impl Canonicalize for AromaticSystemConstraintsAst {
    /// Canonicalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Canonicalize::canonicalize)
            .collect::<Result<Vec<AromaticSystemConstraintAst>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl Lattice for AromaticSystemConstraintsAst {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`AromaticSystemConstraintAst::meet`; a `None` aborts the whole meet), an
    /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<AromaticSystemConstraintAst> = Vec::new();
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
    /// (`AromaticSystemConstraintAst::join`); a single-side key widens to the absent ⊤ and is dropped.
    /// The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<AromaticSystemConstraintAst> = Vec::new();
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

impl FromIterator<AromaticSystemConstraintAst> for AromaticSystemConstraintsAst {
    fn from_iter<I: IntoIterator<Item = AromaticSystemConstraintAst>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for AromaticSystemConstraintsAst {
    type Item = AromaticSystemConstraintAst;
    type IntoIter = IntoIter<AromaticSystemConstraintAst>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<AromaticSystemConstraintAst> for AromaticSystemConstraintsAst {
    fn from(c: AromaticSystemConstraintAst) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<AromaticSystemConstraintAst>> for AromaticSystemConstraintsAst {
    fn from(cs: Vec<AromaticSystemConstraintAst>) -> Self {
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
        AromaticSystemConstraintAst::electron_count(6),
        AromaticSystemConstraintAst::ElectronCount(ValueAst::Lit(6))
    )]
    fn test_aromatic_system_constraint_constructors(
        #[case] actual: AromaticSystemConstraintAst,
        #[case] expected: AromaticSystemConstraintAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::electron_count(
        AromaticSystemConstraintAst::electron_count(6),
        AromaticSystemConstraintKey::ElectronCount
    )]
    fn test_aromatic_system_constraint_key(
        #[case] c: AromaticSystemConstraintAst,
        #[case] expected: AromaticSystemConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rstest]
    #[case::electron_count(
        AromaticSystemConstraintAst::electron_count(6),
        AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)
    )]
    fn test_aromatic_system_constraint_as_undetermined(
        #[case] c: AromaticSystemConstraintAst,
        #[case] expected: AromaticSystemConstraintAst,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count_litset_singleton(AromaticSystemConstraintAst::ElectronCount(ValueAst::lit_set([6])), Ok(AromaticSystemConstraintAst::electron_count(6)))]
    #[case::empty_litset_contradiction(AromaticSystemConstraintAst::ElectronCount(ValueAst::lit_set(Vec::<i64>::new())), Err(Contradiction))]
    fn test_aromatic_system_constraint_canonicalize(
        #[case] constraint: AromaticSystemConstraintAst,
        #[case] expected: Result<AromaticSystemConstraintAst, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
    }

    #[rstest]
    #[case::lit(AromaticSystemConstraintAst::electron_count(6), false)]
    #[case::undetermined(
        AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined),
        true
    )]
    fn test_aromatic_system_constraint_is_undetermined(
        #[case] c: AromaticSystemConstraintAst,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(AromaticSystemConstraintAst::electron_count(6), AromaticSystemConstraintAst::electron_count(6), Some(AromaticSystemConstraintAst::electron_count(6)))]
    #[case::narrows_undetermined(AromaticSystemConstraintAst::electron_count(6), AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined), Some(AromaticSystemConstraintAst::electron_count(6)))]
    #[case::incompatible(AromaticSystemConstraintAst::electron_count(6), AromaticSystemConstraintAst::electron_count(2), None)]
    fn test_aromatic_system_constraint_meet(#[case] a: AromaticSystemConstraintAst, #[case] b: AromaticSystemConstraintAst, #[case] expected: Option<AromaticSystemConstraintAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(AromaticSystemConstraintAst::electron_count(6), AromaticSystemConstraintAst::electron_count(6), Ok(AromaticSystemConstraintAst::electron_count(6)))]
    #[case::widens(AromaticSystemConstraintAst::electron_count(1), AromaticSystemConstraintAst::electron_count(2), Ok(AromaticSystemConstraintAst::ElectronCount(ValueAst::lit_set([1, 2]))))]
    fn test_aromatic_system_constraint_join(#[case] a: AromaticSystemConstraintAst, #[case] b: AromaticSystemConstraintAst, #[case] expected: Result<AromaticSystemConstraintAst, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(AromaticSystemConstraintAst::electron_count(6), AromaticSystemConstraintAst::electron_count(6), true)]
    #[case::incompatible(AromaticSystemConstraintAst::electron_count(6), AromaticSystemConstraintAst::electron_count(2), false)]
    fn test_aromatic_system_constraint_is_compatible(#[case] a: AromaticSystemConstraintAst, #[case] b: AromaticSystemConstraintAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_aromatic_system_constraints_new() {
        let cs = AromaticSystemConstraintsAst::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    fn test_aromatic_system_constraints_iter() {
        let cs = AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![AromaticSystemConstraintAst::electron_count(6)]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![AromaticSystemConstraintAst::electron_count(6)], vec![AromaticSystemConstraintAst::electron_count(6)])]
    #[case::overwrite_same_key(vec![AromaticSystemConstraintAst::electron_count(6), AromaticSystemConstraintAst::electron_count(10)], vec![AromaticSystemConstraintAst::electron_count(10)])]
    #[case::vacuous_stores(vec![AromaticSystemConstraintAst::electron_count(6), AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)], vec![AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)])]
    fn test_aromatic_system_constraints_set(#[case] sequence: Vec<AromaticSystemConstraintAst>, #[case] expected: Vec<AromaticSystemConstraintAst>) {
        let mut cs = AromaticSystemConstraintsAst::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, AromaticSystemConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite(
        vec![AromaticSystemConstraintAst::electron_count(6)],
        vec![AromaticSystemConstraintAst::electron_count(10)],
        vec![AromaticSystemConstraintAst::electron_count(10)])]
    #[case::adds_from_empty(
        vec![],
        vec![AromaticSystemConstraintAst::electron_count(6)],
        vec![AromaticSystemConstraintAst::electron_count(6)])]
    #[case::vacuous_removes(
        vec![AromaticSystemConstraintAst::electron_count(6)],
        vec![AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)],
        vec![])]
    fn test_aromatic_system_constraints_update(#[case] initial: Vec<AromaticSystemConstraintAst>, #[case] other: Vec<AromaticSystemConstraintAst>, #[case] expected: Vec<AromaticSystemConstraintAst>) {
        let mut cs = AromaticSystemConstraintsAst::from_iter(initial);
        cs.update(&AromaticSystemConstraintsAst::from_iter(other));
        assert_eq!(cs, AromaticSystemConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![AromaticSystemConstraintAst::electron_count(6)], Some(AromaticSystemConstraintAst::electron_count(6)), Some(AromaticSystemConstraintAst::electron_count(10)), Ok(()), vec![AromaticSystemConstraintAst::electron_count(10)])]
    #[case::remove(vec![AromaticSystemConstraintAst::electron_count(6)], Some(AromaticSystemConstraintAst::electron_count(6)), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(AromaticSystemConstraintAst::electron_count(6)), Ok(()), vec![AromaticSystemConstraintAst::electron_count(6)])]
    #[case::old_mismatch(vec![AromaticSystemConstraintAst::electron_count(6)], Some(AromaticSystemConstraintAst::electron_count(2)), None, Err(Contradiction), vec![AromaticSystemConstraintAst::electron_count(6)])]
    fn test_aromatic_system_constraints_compare_and_set(
        #[case] initial: Vec<AromaticSystemConstraintAst>,
        #[case] old: Option<AromaticSystemConstraintAst>,
        #[case] new: Option<AromaticSystemConstraintAst>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<AromaticSystemConstraintAst>,
    ) {
        let mut cs = AromaticSystemConstraintsAst::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, AromaticSystemConstraintsAst::from_iter(expected_state));
    }

    #[rstest]
    #[case::present(
        AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)),
        true
    )]
    #[case::absent(AromaticSystemConstraintsAst::new(), false)]
    fn test_aromatic_system_constraints_contains(
        #[case] cs: AromaticSystemConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(
            cs.contains(AromaticSystemConstraintKey::ElectronCount),
            expected
        );
    }

    #[rstest]
    #[case::present(
        AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)),
        Some(AromaticSystemConstraintAst::electron_count(6))
    )]
    #[case::absent(AromaticSystemConstraintsAst::new(), None)]
    fn test_aromatic_system_constraints_get(
        #[case] cs: AromaticSystemConstraintsAst,
        #[case] expected: Option<AromaticSystemConstraintAst>,
    ) {
        assert_eq!(
            cs.get(AromaticSystemConstraintKey::ElectronCount),
            expected.as_ref()
        );
    }

    #[rstest]
    #[case::present(
        AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)),
        Some(AromaticSystemConstraintAst::electron_count(6)),
        AromaticSystemConstraintsAst::new()
    )]
    #[case::absent(
        AromaticSystemConstraintsAst::new(),
        None,
        AromaticSystemConstraintsAst::new()
    )]
    fn test_aromatic_system_constraints_remove(
        #[case] mut cs: AromaticSystemConstraintsAst,
        #[case] expected_removed: Option<AromaticSystemConstraintAst>,
        #[case] expected_state: AromaticSystemConstraintsAst,
    ) {
        assert_eq!(
            cs.remove(AromaticSystemConstraintKey::ElectronCount),
            expected_removed
        );
        assert_eq!(cs, expected_state);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &AromaticSystemConstraintAst| matches!(c, AromaticSystemConstraintAst::ElectronCount(_)), vec![AromaticSystemConstraintAst::electron_count(6)])]
    #[case::all_dropped(|_: &AromaticSystemConstraintAst| false, vec![])]
    fn test_aromatic_system_constraints_retain(
        #[case] predicate: impl FnMut(&AromaticSystemConstraintAst) -> bool,
        #[case] expected: Vec<AromaticSystemConstraintAst>,
    ) {
        let mut cs = AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6));
        cs.retain(predicate);
        assert_eq!(cs, AromaticSystemConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_aromatic_system_constraints_clear() {
        let mut cs =
            AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6));
        cs.clear();
        assert_eq!(cs, AromaticSystemConstraintsAst::new());
    }

    #[rstest]
    fn test_aromatic_system_constraints_take() {
        let mut cs =
            AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6));
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(
            drained,
            vec![AromaticSystemConstraintAst::electron_count(6)]
        );
        assert_eq!(cs, AromaticSystemConstraintsAst::new());
    }

    #[rstest]
    fn test_aromatic_system_constraints_compact() {
        let cs = AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6));
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
        AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::lit_set([6]))),
        Ok(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6))))]
    #[case::drop_vacuous(
        AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)),
        Ok(AromaticSystemConstraintsAst::new()))]
    #[case::contradiction(
        AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_aromatic_system_constraints_canonicalize(
        #[case] constraints: AromaticSystemConstraintsAst,
        #[case] expected: Result<AromaticSystemConstraintsAst, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::a_only_kept(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), AromaticSystemConstraintsAst::new(),
        Some(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6))))]
    #[case::b_only_kept(AromaticSystemConstraintsAst::new(), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)),
        Some(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6))))]
    #[case::shared_key_meets(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)),
        Some(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6))))]
    #[case::shared_key_contradicts(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(2)), None)]
    #[case::prunes_vacuous(AromaticSystemConstraintsAst::new(), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)), Some(AromaticSystemConstraintsAst::new()))]
    fn test_aromatic_system_constraints_meet(#[case] a: AromaticSystemConstraintsAst, #[case] b: AromaticSystemConstraintsAst, #[case] expected: Option<AromaticSystemConstraintsAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::widens_value(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(1)), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(2)),
        AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::lit_set([1, 2]))))]
    #[case::single_side_dropped(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), AromaticSystemConstraintsAst::new(),
        AromaticSystemConstraintsAst::new())]
    #[case::undetermined_drops(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)),
        AromaticSystemConstraintsAst::new())]
    fn test_aromatic_system_constraints_join(#[case] a: AromaticSystemConstraintsAst, #[case] b: AromaticSystemConstraintsAst, #[case] expected: AromaticSystemConstraintsAst) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(AromaticSystemConstraintsAst::new(), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), true)]
    #[case::required_present(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), true)]
    #[case::required_absent(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), AromaticSystemConstraintsAst::new(), false)]
    #[case::wildcard_matches_lit(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), true)]
    #[case::lit_mismatch(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(2)), false)]
    fn test_aromatic_system_constraints_matches(
        #[case] pattern: AromaticSystemConstraintsAst,
        #[case] target: AromaticSystemConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_one_empty(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), AromaticSystemConstraintsAst::new(), true)]
    #[case::shared_key_compatible(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), true)]
    #[case::shared_key_incompatible(AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)), AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(2)), false)]
    fn test_aromatic_system_constraints_is_compatible(#[case] a: AromaticSystemConstraintsAst, #[case] b: AromaticSystemConstraintsAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(vec![AromaticSystemConstraintAst::electron_count(6)], vec![AromaticSystemConstraintAst::electron_count(6)])]
    #[case::same_key_last_wins(vec![AromaticSystemConstraintAst::electron_count(2), AromaticSystemConstraintAst::electron_count(6)], vec![AromaticSystemConstraintAst::electron_count(6)])]
    #[case::empty(vec![], vec![])]
    fn test_aromatic_system_constraints_from_iter(
        #[case] input: Vec<AromaticSystemConstraintAst>,
        #[case] expected: Vec<AromaticSystemConstraintAst>,
    ) {
        let cs = AromaticSystemConstraintsAst::from_iter(input);
        assert_eq!(cs, AromaticSystemConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_aromatic_system_constraints_into_iter() {
        let cs = AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6));
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![AromaticSystemConstraintAst::electron_count(6)]
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_from_aromatic_system_constraint() {
        let cs: AromaticSystemConstraintsAst =
            AromaticSystemConstraintAst::electron_count(6).into();
        assert_eq!(
            cs,
            AromaticSystemConstraintsAst::from_iter([AromaticSystemConstraintAst::electron_count(
                6
            )]),
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_from_vec() {
        let cs: AromaticSystemConstraintsAst =
            vec![AromaticSystemConstraintAst::electron_count(6)].into();
        assert_eq!(
            cs,
            AromaticSystemConstraintsAst::from_iter([AromaticSystemConstraintAst::electron_count(
                6
            )]),
        );
    }
}
