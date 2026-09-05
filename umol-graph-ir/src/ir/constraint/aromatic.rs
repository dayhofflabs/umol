//! Per-aromatic-system constraints.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use umol_perm::DynPermutation;

use super::super::compact::MoleculeCompaction;
use super::super::error::{Contradiction, NoJoin};
use super::super::num::NumForm;
use super::super::traits::{FrameTransport, Lattice, Normalize};

/// Aromatic-system-scope constraint. Held inline on `AromaticSystemForm` via
/// `AromaticSystemConstraintsForm`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticSystemConstraintForm {
    /// Asserted total π-electron count for the system. Cross-checked by the
    /// `ConsistencyValidator` against `sum(AromaticSystemForm::electrons)`.
    ElectronCount(NumForm),
}

impl AromaticSystemConstraintForm {
    pub fn electron_count(v: impl Into<NumForm>) -> Self {
        Self::ElectronCount(v.into())
    }

    /// Aromatic-system constraint key, unique within an `AromaticSystemConstraintsForm` container.
    pub fn key(&self) -> AromaticSystemConstraintKey {
        match self {
            Self::ElectronCount(_) => AromaticSystemConstraintKey::ElectronCount,
        }
    }

    /// Vacuous form of constraint key, used for removal.
    pub fn as_undetermined(&self) -> Self {
        match self {
            Self::ElectronCount(_) => Self::ElectronCount(NumForm::Undetermined),
        }
    }

    pub fn compact(self, _compaction: &MoleculeCompaction) -> Option<Self> {
        // Value-only: no indices to compact.
        Some(self)
    }

    pub(crate) fn uses_participant_frame(&self) -> bool {
        match self {
            Self::ElectronCount(_) => false,
        }
    }
}

impl Normalize for AromaticSystemConstraintForm {
    /// Normalize the inner value; the kind is preserved.
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::ElectronCount(v) => Self::ElectronCount(v.normalize()?),
        })
    }
}

impl FrameTransport for AromaticSystemConstraintForm {
    type Action = DynPermutation;

    fn reframe_by(self, _action: &Self::Action) -> Option<Self> {
        Some(match self {
            Self::ElectronCount(value) => Self::ElectronCount(value),
        })
    }
}

impl Lattice for AromaticSystemConstraintForm {
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
pub struct AromaticSystemConstraintsForm(Vec<AromaticSystemConstraintForm>);

impl AromaticSystemConstraintsForm {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The system's electron count, or `Undetermined` when no `ElectronCount` constraint is present.
    pub fn electron_count(&self) -> NumForm {
        match self.get(AromaticSystemConstraintKey::ElectronCount) {
            Some(AromaticSystemConstraintForm::ElectronCount(v)) => v.clone(),
            _ => NumForm::Undetermined,
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

    pub fn get(&self, key: AromaticSystemConstraintKey) -> Option<&AromaticSystemConstraintForm> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: AromaticSystemConstraintForm) {
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
        old: Option<AromaticSystemConstraintForm>,
        new: Option<AromaticSystemConstraintForm>,
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

    pub fn remove(
        &mut self,
        key: AromaticSystemConstraintKey,
    ) -> Option<AromaticSystemConstraintForm> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = AromaticSystemConstraintForm>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &AromaticSystemConstraintsForm) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&AromaticSystemConstraintForm) -> bool) {
        self.0.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl ExactSizeIterator<Item = AromaticSystemConstraintForm> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, AromaticSystemConstraintForm> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &MoleculeCompaction) -> Self {
        self
    }
}

impl Normalize for AromaticSystemConstraintsForm {
    /// Normalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn normalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Normalize::normalize)
            .collect::<Result<Vec<AromaticSystemConstraintForm>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl FrameTransport for AromaticSystemConstraintsForm {
    type Action = DynPermutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        if !self
            .iter()
            .any(AromaticSystemConstraintForm::uses_participant_frame)
        {
            return Some(self);
        }
        self.into_iter()
            .map(|constraint| constraint.reframe_by(action))
            .collect()
    }
}

impl Lattice for AromaticSystemConstraintsForm {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`AromaticSystemConstraintForm::meet`; a `None` aborts the whole meet), an
    /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<AromaticSystemConstraintForm> = Vec::new();
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
    /// (`AromaticSystemConstraintForm::join`); a single-side key widens to the absent ⊤ and is dropped.
    /// The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<AromaticSystemConstraintForm> = Vec::new();
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

impl FromIterator<AromaticSystemConstraintForm> for AromaticSystemConstraintsForm {
    fn from_iter<I: IntoIterator<Item = AromaticSystemConstraintForm>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for AromaticSystemConstraintsForm {
    type Item = AromaticSystemConstraintForm;
    type IntoIter = IntoIter<AromaticSystemConstraintForm>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<AromaticSystemConstraintForm> for AromaticSystemConstraintsForm {
    fn from(c: AromaticSystemConstraintForm) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<AromaticSystemConstraintForm>> for AromaticSystemConstraintsForm {
    fn from(cs: Vec<AromaticSystemConstraintForm>) -> Self {
        Self::from_iter(cs)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::{Compaction, EdgeId, GraphCompaction, NodeId};

    use super::*;

    #[rstest]
    #[case::electron_count(
        AromaticSystemConstraintForm::electron_count(6),
        AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))
    )]
    fn test_aromatic_system_constraint_form_constructors(
        #[case] actual: AromaticSystemConstraintForm,
        #[case] expected: AromaticSystemConstraintForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::electron_count(
        AromaticSystemConstraintForm::electron_count(6),
        AromaticSystemConstraintKey::ElectronCount
    )]
    fn test_aromatic_system_constraint_form_key(
        #[case] c: AromaticSystemConstraintForm,
        #[case] expected: AromaticSystemConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rstest]
    #[case::electron_count(
        AromaticSystemConstraintForm::electron_count(6),
        AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)
    )]
    fn test_aromatic_system_constraint_form_as_undetermined(
        #[case] c: AromaticSystemConstraintForm,
        #[case] expected: AromaticSystemConstraintForm,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count_litset_singleton(AromaticSystemConstraintForm::ElectronCount(NumForm::lit_set([6])), Ok(AromaticSystemConstraintForm::electron_count(6)))]
    #[case::empty_litset_contradiction(AromaticSystemConstraintForm::ElectronCount(NumForm::lit_set(Vec::<i64>::new())), Err(Contradiction))]
    fn test_aromatic_system_constraint_form_normalize(
        #[case] constraint: AromaticSystemConstraintForm,
        #[case] expected: Result<AromaticSystemConstraintForm, Contradiction>,
    ) {
        assert_eq!(constraint.normalize(), expected);
    }

    #[rstest]
    #[case::electron_count(AromaticSystemConstraintForm::electron_count(6), false)]
    fn test_aromatic_system_constraint_form_uses_participant_frame(
        #[case] constraint: AromaticSystemConstraintForm,
        #[case] expected: bool,
    ) {
        assert_eq!(constraint.uses_participant_frame(), expected);
    }

    #[rstest]
    #[case::electron_count(AromaticSystemConstraintForm::electron_count(6))]
    fn test_aromatic_system_constraint_form_reframe_by(
        #[case] constraint: AromaticSystemConstraintForm,
    ) {
        let action = DynPermutation::try_from(vec![1, 0]).expect("the action is a permutation");

        assert_eq!(constraint.clone().reframe_by(&action), Some(constraint));
    }

    #[rstest]
    fn test_aromatic_system_constraints_form_reframe_by() {
        let constraints = AromaticSystemConstraintsForm::from(vec![
            AromaticSystemConstraintForm::electron_count(6),
        ]);
        let action = DynPermutation::try_from(vec![1, 0]).expect("the action is a permutation");

        assert_eq!(constraints.clone().reframe_by(&action), Some(constraints),);
    }

    #[rstest]
    #[case::lit(AromaticSystemConstraintForm::electron_count(6), false)]
    #[case::undetermined(
        AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined),
        true
    )]
    fn test_aromatic_system_constraint_form_is_undetermined(
        #[case] c: AromaticSystemConstraintForm,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(AromaticSystemConstraintForm::electron_count(6), AromaticSystemConstraintForm::electron_count(6), Some(AromaticSystemConstraintForm::electron_count(6)))]
    #[case::narrows_undetermined(AromaticSystemConstraintForm::electron_count(6), AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined), Some(AromaticSystemConstraintForm::electron_count(6)))]
    #[case::incompatible(AromaticSystemConstraintForm::electron_count(6), AromaticSystemConstraintForm::electron_count(2), None)]
    fn test_aromatic_system_constraint_form_meet(#[case] a: AromaticSystemConstraintForm, #[case] b: AromaticSystemConstraintForm, #[case] expected: Option<AromaticSystemConstraintForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(AromaticSystemConstraintForm::electron_count(6), AromaticSystemConstraintForm::electron_count(6), Ok(AromaticSystemConstraintForm::electron_count(6)))]
    #[case::widens(AromaticSystemConstraintForm::electron_count(1), AromaticSystemConstraintForm::electron_count(2), Ok(AromaticSystemConstraintForm::ElectronCount(NumForm::lit_set([1, 2]))))]
    fn test_aromatic_system_constraint_form_join(#[case] a: AromaticSystemConstraintForm, #[case] b: AromaticSystemConstraintForm, #[case] expected: Result<AromaticSystemConstraintForm, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(AromaticSystemConstraintForm::electron_count(6), AromaticSystemConstraintForm::electron_count(6), true)]
    #[case::incompatible(AromaticSystemConstraintForm::electron_count(6), AromaticSystemConstraintForm::electron_count(2), false)]
    fn test_aromatic_system_constraint_form_is_compatible(#[case] a: AromaticSystemConstraintForm, #[case] b: AromaticSystemConstraintForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_aromatic_system_constraints_form_new() {
        let cs = AromaticSystemConstraintsForm::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    fn test_aromatic_system_constraints_form_iter() {
        let cs =
            AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![AromaticSystemConstraintForm::electron_count(6)]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![AromaticSystemConstraintForm::electron_count(6)], vec![AromaticSystemConstraintForm::electron_count(6)])]
    #[case::overwrite_same_key(vec![AromaticSystemConstraintForm::electron_count(6), AromaticSystemConstraintForm::electron_count(10)], vec![AromaticSystemConstraintForm::electron_count(10)])]
    #[case::vacuous_stores(vec![AromaticSystemConstraintForm::electron_count(6), AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)], vec![AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)])]
    fn test_aromatic_system_constraints_form_set(#[case] sequence: Vec<AromaticSystemConstraintForm>, #[case] expected: Vec<AromaticSystemConstraintForm>) {
        let mut cs = AromaticSystemConstraintsForm::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, AromaticSystemConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite(
        vec![AromaticSystemConstraintForm::electron_count(6)],
        vec![AromaticSystemConstraintForm::electron_count(10)],
        vec![AromaticSystemConstraintForm::electron_count(10)])]
    #[case::adds_from_empty(
        vec![],
        vec![AromaticSystemConstraintForm::electron_count(6)],
        vec![AromaticSystemConstraintForm::electron_count(6)])]
    #[case::vacuous_removes(
        vec![AromaticSystemConstraintForm::electron_count(6)],
        vec![AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)],
        vec![])]
    fn test_aromatic_system_constraints_form_update(#[case] initial: Vec<AromaticSystemConstraintForm>, #[case] other: Vec<AromaticSystemConstraintForm>, #[case] expected: Vec<AromaticSystemConstraintForm>) {
        let mut cs = AromaticSystemConstraintsForm::from_iter(initial);
        cs.update(&AromaticSystemConstraintsForm::from_iter(other));
        assert_eq!(cs, AromaticSystemConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![AromaticSystemConstraintForm::electron_count(6)], Some(AromaticSystemConstraintForm::electron_count(6)), Some(AromaticSystemConstraintForm::electron_count(10)), Ok(()), vec![AromaticSystemConstraintForm::electron_count(10)])]
    #[case::remove(vec![AromaticSystemConstraintForm::electron_count(6)], Some(AromaticSystemConstraintForm::electron_count(6)), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(AromaticSystemConstraintForm::electron_count(6)), Ok(()), vec![AromaticSystemConstraintForm::electron_count(6)])]
    #[case::old_mismatch(vec![AromaticSystemConstraintForm::electron_count(6)], Some(AromaticSystemConstraintForm::electron_count(2)), None, Err(Contradiction), vec![AromaticSystemConstraintForm::electron_count(6)])]
    fn test_aromatic_system_constraints_form_compare_and_set(
        #[case] initial: Vec<AromaticSystemConstraintForm>,
        #[case] old: Option<AromaticSystemConstraintForm>,
        #[case] new: Option<AromaticSystemConstraintForm>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<AromaticSystemConstraintForm>,
    ) {
        let mut cs = AromaticSystemConstraintsForm::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, AromaticSystemConstraintsForm::from_iter(expected_state));
    }

    #[rstest]
    #[case::present(
        AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)),
        true
    )]
    #[case::absent(AromaticSystemConstraintsForm::new(), false)]
    fn test_aromatic_system_constraints_form_contains(
        #[case] cs: AromaticSystemConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(
            cs.contains(AromaticSystemConstraintKey::ElectronCount),
            expected
        );
    }

    #[rstest]
    #[case::present(
        AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)),
        Some(AromaticSystemConstraintForm::electron_count(6))
    )]
    #[case::absent(AromaticSystemConstraintsForm::new(), None)]
    fn test_aromatic_system_constraints_form_get(
        #[case] cs: AromaticSystemConstraintsForm,
        #[case] expected: Option<AromaticSystemConstraintForm>,
    ) {
        assert_eq!(
            cs.get(AromaticSystemConstraintKey::ElectronCount),
            expected.as_ref()
        );
    }

    #[rstest]
    #[case::present(
        AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)),
        Some(AromaticSystemConstraintForm::electron_count(6)),
        AromaticSystemConstraintsForm::new()
    )]
    #[case::absent(
        AromaticSystemConstraintsForm::new(),
        None,
        AromaticSystemConstraintsForm::new()
    )]
    fn test_aromatic_system_constraints_form_remove(
        #[case] mut cs: AromaticSystemConstraintsForm,
        #[case] expected_removed: Option<AromaticSystemConstraintForm>,
        #[case] expected_state: AromaticSystemConstraintsForm,
    ) {
        assert_eq!(
            cs.remove(AromaticSystemConstraintKey::ElectronCount),
            expected_removed
        );
        assert_eq!(cs, expected_state);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &AromaticSystemConstraintForm| matches!(c, AromaticSystemConstraintForm::ElectronCount(_)), vec![AromaticSystemConstraintForm::electron_count(6)])]
    #[case::all_dropped(|_: &AromaticSystemConstraintForm| false, vec![])]
    fn test_aromatic_system_constraints_form_retain(
        #[case] predicate: impl FnMut(&AromaticSystemConstraintForm) -> bool,
        #[case] expected: Vec<AromaticSystemConstraintForm>,
    ) {
        let mut cs = AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6));
        cs.retain(predicate);
        assert_eq!(cs, AromaticSystemConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_aromatic_system_constraints_form_clear() {
        let mut cs =
            AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6));
        cs.clear();
        assert_eq!(cs, AromaticSystemConstraintsForm::new());
    }

    #[rstest]
    fn test_aromatic_system_constraints_form_take() {
        let mut empty = AromaticSystemConstraintsForm::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let mut cs =
            AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6));
        let mut taken = cs.take();
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(
            taken.next(),
            Some(AromaticSystemConstraintForm::electron_count(6)),
        );
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.size_hint(), (0, Some(0)));
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(cs, AromaticSystemConstraintsForm::new());
    }

    #[rstest]
    fn test_aromatic_system_constraints_form_compact() {
        let cs =
            AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6));
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::new(2, vec![NodeId(0), NodeId(1)]).unwrap(),
                Compaction::new(1, vec![EdgeId(0)]).unwrap(),
            ),
            Compaction::empty(),
            Compaction::empty(),
            Compaction::empty(),
            Compaction::empty(),
            Compaction::empty(),
            Compaction::empty(),
        );
        assert_eq!(cs.clone().compact(&compaction), cs);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::normalizes_value(
        AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::lit_set([6]))),
        Ok(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6))))]
    #[case::drop_vacuous(
        AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)),
        Ok(AromaticSystemConstraintsForm::new()))]
    #[case::contradiction(
        AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_aromatic_system_constraints_form_normalize(
        #[case] constraints: AromaticSystemConstraintsForm,
        #[case] expected: Result<AromaticSystemConstraintsForm, Contradiction>,
    ) {
        assert_eq!(constraints.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::a_only_kept(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), AromaticSystemConstraintsForm::new(),
        Some(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6))))]
    #[case::b_only_kept(AromaticSystemConstraintsForm::new(), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)),
        Some(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6))))]
    #[case::shared_key_meets(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)),
        Some(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6))))]
    #[case::shared_key_contradicts(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(2)), None)]
    #[case::prunes_vacuous(AromaticSystemConstraintsForm::new(), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)), Some(AromaticSystemConstraintsForm::new()))]
    fn test_aromatic_system_constraints_form_meet(#[case] a: AromaticSystemConstraintsForm, #[case] b: AromaticSystemConstraintsForm, #[case] expected: Option<AromaticSystemConstraintsForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::widens_value(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(1)), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(2)),
        AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::lit_set([1, 2]))))]
    #[case::single_side_dropped(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), AromaticSystemConstraintsForm::new(),
        AromaticSystemConstraintsForm::new())]
    #[case::undetermined_drops(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)),
        AromaticSystemConstraintsForm::new())]
    fn test_aromatic_system_constraints_form_join(#[case] a: AromaticSystemConstraintsForm, #[case] b: AromaticSystemConstraintsForm, #[case] expected: AromaticSystemConstraintsForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(AromaticSystemConstraintsForm::new(), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), true)]
    #[case::required_present(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), true)]
    #[case::required_absent(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), AromaticSystemConstraintsForm::new(), false)]
    #[case::wildcard_matches_lit(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), true)]
    #[case::lit_mismatch(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(2)), false)]
    fn test_aromatic_system_constraints_form_matches(
        #[case] pattern: AromaticSystemConstraintsForm,
        #[case] target: AromaticSystemConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_one_empty(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), AromaticSystemConstraintsForm::new(), true)]
    #[case::shared_key_compatible(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), true)]
    #[case::shared_key_incompatible(AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)), AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(2)), false)]
    fn test_aromatic_system_constraints_form_is_compatible(#[case] a: AromaticSystemConstraintsForm, #[case] b: AromaticSystemConstraintsForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(vec![AromaticSystemConstraintForm::electron_count(6)], vec![AromaticSystemConstraintForm::electron_count(6)])]
    #[case::same_key_last_wins(vec![AromaticSystemConstraintForm::electron_count(2), AromaticSystemConstraintForm::electron_count(6)], vec![AromaticSystemConstraintForm::electron_count(6)])]
    #[case::empty(vec![], vec![])]
    fn test_aromatic_system_constraints_form_from_iter(
        #[case] input: Vec<AromaticSystemConstraintForm>,
        #[case] expected: Vec<AromaticSystemConstraintForm>,
    ) {
        let cs = AromaticSystemConstraintsForm::from_iter(input);
        assert_eq!(cs, AromaticSystemConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_aromatic_system_constraints_form_into_iter() {
        let cs =
            AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6));
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![AromaticSystemConstraintForm::electron_count(6)]
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_form_from_aromatic_system_constraint() {
        let cs: AromaticSystemConstraintsForm =
            AromaticSystemConstraintForm::electron_count(6).into();
        assert_eq!(
            cs,
            AromaticSystemConstraintsForm::from_iter([
                AromaticSystemConstraintForm::electron_count(6)
            ]),
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_form_from_vec() {
        let cs: AromaticSystemConstraintsForm =
            vec![AromaticSystemConstraintForm::electron_count(6)].into();
        assert_eq!(
            cs,
            AromaticSystemConstraintsForm::from_iter([
                AromaticSystemConstraintForm::electron_count(6)
            ]),
        );
    }
}
