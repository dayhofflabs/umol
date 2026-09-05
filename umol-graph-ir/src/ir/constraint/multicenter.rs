//! Per-multicenter-bond constraints.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use umol_perm::DynPermutation;

use super::super::compact::MoleculeCompaction;
use super::super::error::{Contradiction, NoJoin};
use super::super::num::NumForm;
use super::super::traits::{FrameTransport, Lattice, Normalize};

/// Multicenter-bond-scope constraint. Held inline on `MulticenterBondForm` via
/// `MulticenterBondConstraintsForm`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MulticenterBondConstraintForm {
    /// Asserted total electron count for the multicenter bond. Cross-checked
    /// by the `ConsistencyValidator` against `sum(MulticenterBondForm::electrons)`.
    ElectronCount(NumForm),
}

impl MulticenterBondConstraintForm {
    pub fn electron_count(v: impl Into<NumForm>) -> Self {
        Self::ElectronCount(v.into())
    }

    /// Multicenter-bond constraint key, unique within a `MulticenterBondConstraintsForm` container.
    pub fn key(&self) -> MulticenterBondConstraintKey {
        match self {
            Self::ElectronCount(_) => MulticenterBondConstraintKey::ElectronCount,
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

impl Normalize for MulticenterBondConstraintForm {
    /// Normalize the inner value; the kind is preserved.
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::ElectronCount(v) => Self::ElectronCount(v.normalize()?),
        })
    }
}

impl FrameTransport for MulticenterBondConstraintForm {
    type Action = DynPermutation;

    fn reframe_by(self, _action: &Self::Action) -> Option<Self> {
        Some(match self {
            Self::ElectronCount(value) => Self::ElectronCount(value),
        })
    }
}

impl Lattice for MulticenterBondConstraintForm {
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
pub struct MulticenterBondConstraintsForm(Vec<MulticenterBondConstraintForm>);

impl MulticenterBondConstraintsForm {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The bond's electron count, or `Undetermined` when no `ElectronCount` constraint is present.
    pub fn electron_count(&self) -> NumForm {
        match self.get(MulticenterBondConstraintKey::ElectronCount) {
            Some(MulticenterBondConstraintForm::ElectronCount(v)) => v.clone(),
            _ => NumForm::Undetermined,
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

    pub fn get(&self, key: MulticenterBondConstraintKey) -> Option<&MulticenterBondConstraintForm> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: MulticenterBondConstraintForm) {
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
        old: Option<MulticenterBondConstraintForm>,
        new: Option<MulticenterBondConstraintForm>,
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
        key: MulticenterBondConstraintKey,
    ) -> Option<MulticenterBondConstraintForm> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = MulticenterBondConstraintForm>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &MulticenterBondConstraintsForm) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&MulticenterBondConstraintForm) -> bool) {
        self.0.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl ExactSizeIterator<Item = MulticenterBondConstraintForm> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, MulticenterBondConstraintForm> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &MoleculeCompaction) -> Self {
        self
    }
}

impl Normalize for MulticenterBondConstraintsForm {
    /// Normalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn normalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Normalize::normalize)
            .collect::<Result<Vec<MulticenterBondConstraintForm>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl FrameTransport for MulticenterBondConstraintsForm {
    type Action = DynPermutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        if !self
            .iter()
            .any(MulticenterBondConstraintForm::uses_participant_frame)
        {
            return Some(self);
        }
        self.into_iter()
            .map(|constraint| constraint.reframe_by(action))
            .collect()
    }
}

impl Lattice for MulticenterBondConstraintsForm {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`MulticenterBondConstraintForm::meet`; a `None` aborts the whole meet), an
    /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<MulticenterBondConstraintForm> = Vec::new();
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
    /// (`MulticenterBondConstraintForm::join`); a single-side key widens to the absent ⊤ and is dropped.
    /// The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<MulticenterBondConstraintForm> = Vec::new();
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

impl FromIterator<MulticenterBondConstraintForm> for MulticenterBondConstraintsForm {
    fn from_iter<I: IntoIterator<Item = MulticenterBondConstraintForm>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for MulticenterBondConstraintsForm {
    type Item = MulticenterBondConstraintForm;
    type IntoIter = IntoIter<MulticenterBondConstraintForm>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<MulticenterBondConstraintForm> for MulticenterBondConstraintsForm {
    fn from(c: MulticenterBondConstraintForm) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<MulticenterBondConstraintForm>> for MulticenterBondConstraintsForm {
    fn from(cs: Vec<MulticenterBondConstraintForm>) -> Self {
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
        MulticenterBondConstraintForm::electron_count(6),
        MulticenterBondConstraintForm::ElectronCount(NumForm::Lit(6))
    )]
    fn test_multicenter_bond_constraint_form_constructors(
        #[case] actual: MulticenterBondConstraintForm,
        #[case] expected: MulticenterBondConstraintForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraintForm::electron_count(6),
        MulticenterBondConstraintKey::ElectronCount
    )]
    fn test_multicenter_bond_constraint_form_key(
        #[case] c: MulticenterBondConstraintForm,
        #[case] expected: MulticenterBondConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraintForm::electron_count(6),
        MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined)
    )]
    fn test_multicenter_bond_constraint_form_as_undetermined(
        #[case] c: MulticenterBondConstraintForm,
        #[case] expected: MulticenterBondConstraintForm,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count_litset_singleton(MulticenterBondConstraintForm::ElectronCount(NumForm::lit_set([6])), Ok(MulticenterBondConstraintForm::electron_count(6)))]
    #[case::empty_litset_contradiction(MulticenterBondConstraintForm::ElectronCount(NumForm::lit_set(Vec::<i64>::new())), Err(Contradiction))]
    fn test_multicenter_bond_constraint_form_normalize(
        #[case] constraint: MulticenterBondConstraintForm,
        #[case] expected: Result<MulticenterBondConstraintForm, Contradiction>,
    ) {
        assert_eq!(constraint.normalize(), expected);
    }

    #[rstest]
    #[case::electron_count(MulticenterBondConstraintForm::electron_count(6), false)]
    fn test_multicenter_bond_constraint_form_uses_participant_frame(
        #[case] constraint: MulticenterBondConstraintForm,
        #[case] expected: bool,
    ) {
        assert_eq!(constraint.uses_participant_frame(), expected);
    }

    #[rstest]
    #[case::electron_count(MulticenterBondConstraintForm::electron_count(6))]
    fn test_multicenter_bond_constraint_form_reframe_by(
        #[case] constraint: MulticenterBondConstraintForm,
    ) {
        let action = DynPermutation::try_from(vec![1, 0]).expect("the action is a permutation");

        assert_eq!(constraint.clone().reframe_by(&action), Some(constraint));
    }

    #[rstest]
    fn test_multicenter_bond_constraints_form_reframe_by() {
        let constraints = MulticenterBondConstraintsForm::from(vec![
            MulticenterBondConstraintForm::electron_count(6),
        ]);
        let action = DynPermutation::try_from(vec![1, 0]).expect("the action is a permutation");

        assert_eq!(constraints.clone().reframe_by(&action), Some(constraints),);
    }

    #[rstest]
    #[case::lit(MulticenterBondConstraintForm::electron_count(6), false)]
    #[case::undetermined(
        MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined),
        true
    )]
    fn test_multicenter_bond_constraint_form_is_undetermined(
        #[case] c: MulticenterBondConstraintForm,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(MulticenterBondConstraintForm::electron_count(6), MulticenterBondConstraintForm::electron_count(6), Some(MulticenterBondConstraintForm::electron_count(6)))]
    #[case::narrows_undetermined(MulticenterBondConstraintForm::electron_count(6), MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined), Some(MulticenterBondConstraintForm::electron_count(6)))]
    #[case::incompatible(MulticenterBondConstraintForm::electron_count(6), MulticenterBondConstraintForm::electron_count(2), None)]
    fn test_multicenter_bond_constraint_form_meet(#[case] a: MulticenterBondConstraintForm, #[case] b: MulticenterBondConstraintForm, #[case] expected: Option<MulticenterBondConstraintForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(MulticenterBondConstraintForm::electron_count(6), MulticenterBondConstraintForm::electron_count(6), Ok(MulticenterBondConstraintForm::electron_count(6)))]
    #[case::widens(MulticenterBondConstraintForm::electron_count(1), MulticenterBondConstraintForm::electron_count(2), Ok(MulticenterBondConstraintForm::ElectronCount(NumForm::lit_set([1, 2]))))]
    fn test_multicenter_bond_constraint_form_join(#[case] a: MulticenterBondConstraintForm, #[case] b: MulticenterBondConstraintForm, #[case] expected: Result<MulticenterBondConstraintForm, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(MulticenterBondConstraintForm::electron_count(6), MulticenterBondConstraintForm::electron_count(6), true)]
    #[case::incompatible(MulticenterBondConstraintForm::electron_count(6), MulticenterBondConstraintForm::electron_count(2), false)]
    fn test_multicenter_bond_constraint_form_is_compatible(#[case] a: MulticenterBondConstraintForm, #[case] b: MulticenterBondConstraintForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_form_new() {
        let cs = MulticenterBondConstraintsForm::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_form_iter() {
        let cs =
            MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![MulticenterBondConstraintForm::electron_count(6)]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![MulticenterBondConstraintForm::electron_count(6)], vec![MulticenterBondConstraintForm::electron_count(6)])]
    #[case::overwrite_same_key(vec![MulticenterBondConstraintForm::electron_count(6), MulticenterBondConstraintForm::electron_count(10)], vec![MulticenterBondConstraintForm::electron_count(10)])]
    #[case::vacuous_stores(vec![MulticenterBondConstraintForm::electron_count(6), MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined)], vec![MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined)])]
    fn test_multicenter_bond_constraints_form_set(#[case] sequence: Vec<MulticenterBondConstraintForm>, #[case] expected: Vec<MulticenterBondConstraintForm>) {
        let mut cs = MulticenterBondConstraintsForm::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, MulticenterBondConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite(
        vec![MulticenterBondConstraintForm::electron_count(6)],
        vec![MulticenterBondConstraintForm::electron_count(10)],
        vec![MulticenterBondConstraintForm::electron_count(10)])]
    #[case::adds_from_empty(
        vec![],
        vec![MulticenterBondConstraintForm::electron_count(6)],
        vec![MulticenterBondConstraintForm::electron_count(6)])]
    #[case::vacuous_removes(
        vec![MulticenterBondConstraintForm::electron_count(6)],
        vec![MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined)],
        vec![])]
    fn test_multicenter_bond_constraints_form_update(#[case] initial: Vec<MulticenterBondConstraintForm>, #[case] other: Vec<MulticenterBondConstraintForm>, #[case] expected: Vec<MulticenterBondConstraintForm>) {
        let mut cs = MulticenterBondConstraintsForm::from_iter(initial);
        cs.update(&MulticenterBondConstraintsForm::from_iter(other));
        assert_eq!(cs, MulticenterBondConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![MulticenterBondConstraintForm::electron_count(6)], Some(MulticenterBondConstraintForm::electron_count(6)), Some(MulticenterBondConstraintForm::electron_count(10)), Ok(()), vec![MulticenterBondConstraintForm::electron_count(10)])]
    #[case::remove(vec![MulticenterBondConstraintForm::electron_count(6)], Some(MulticenterBondConstraintForm::electron_count(6)), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(MulticenterBondConstraintForm::electron_count(6)), Ok(()), vec![MulticenterBondConstraintForm::electron_count(6)])]
    #[case::old_mismatch(vec![MulticenterBondConstraintForm::electron_count(6)], Some(MulticenterBondConstraintForm::electron_count(2)), None, Err(Contradiction), vec![MulticenterBondConstraintForm::electron_count(6)])]
    fn test_multicenter_bond_constraints_form_compare_and_set(
        #[case] initial: Vec<MulticenterBondConstraintForm>,
        #[case] old: Option<MulticenterBondConstraintForm>,
        #[case] new: Option<MulticenterBondConstraintForm>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<MulticenterBondConstraintForm>,
    ) {
        let mut cs = MulticenterBondConstraintsForm::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, MulticenterBondConstraintsForm::from_iter(expected_state));
    }

    #[rstest]
    #[case::present(
        MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)),
        true
    )]
    #[case::absent(MulticenterBondConstraintsForm::new(), false)]
    fn test_multicenter_bond_constraints_form_contains(
        #[case] cs: MulticenterBondConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(
            cs.contains(MulticenterBondConstraintKey::ElectronCount),
            expected
        );
    }

    #[rstest]
    #[case::present(
        MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)),
        Some(MulticenterBondConstraintForm::electron_count(6))
    )]
    #[case::absent(MulticenterBondConstraintsForm::new(), None)]
    fn test_multicenter_bond_constraints_form_get(
        #[case] cs: MulticenterBondConstraintsForm,
        #[case] expected: Option<MulticenterBondConstraintForm>,
    ) {
        assert_eq!(
            cs.get(MulticenterBondConstraintKey::ElectronCount),
            expected.as_ref()
        );
    }

    #[rstest]
    #[case::present(
        MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)),
        Some(MulticenterBondConstraintForm::electron_count(6)),
        MulticenterBondConstraintsForm::new()
    )]
    #[case::absent(
        MulticenterBondConstraintsForm::new(),
        None,
        MulticenterBondConstraintsForm::new()
    )]
    fn test_multicenter_bond_constraints_form_remove(
        #[case] mut cs: MulticenterBondConstraintsForm,
        #[case] expected_removed: Option<MulticenterBondConstraintForm>,
        #[case] expected_state: MulticenterBondConstraintsForm,
    ) {
        assert_eq!(
            cs.remove(MulticenterBondConstraintKey::ElectronCount),
            expected_removed
        );
        assert_eq!(cs, expected_state);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &MulticenterBondConstraintForm| matches!(c, MulticenterBondConstraintForm::ElectronCount(_)), vec![MulticenterBondConstraintForm::electron_count(6)])]
    #[case::all_dropped(|_: &MulticenterBondConstraintForm| false, vec![])]
    fn test_multicenter_bond_constraints_form_retain(
        #[case] predicate: impl FnMut(&MulticenterBondConstraintForm) -> bool,
        #[case] expected: Vec<MulticenterBondConstraintForm>,
    ) {
        let mut cs = MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6));
        cs.retain(predicate);
        assert_eq!(cs, MulticenterBondConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_multicenter_bond_constraints_form_clear() {
        let mut cs =
            MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6));
        cs.clear();
        assert_eq!(cs, MulticenterBondConstraintsForm::new());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_form_take() {
        let mut empty = MulticenterBondConstraintsForm::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let mut cs =
            MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6));
        let mut taken = cs.take();
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(
            taken.next(),
            Some(MulticenterBondConstraintForm::electron_count(6)),
        );
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.size_hint(), (0, Some(0)));
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(cs, MulticenterBondConstraintsForm::new());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_form_compact() {
        let cs =
            MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6));
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
        MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::ElectronCount(NumForm::lit_set([6]))),
        Ok(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6))))]
    #[case::drop_vacuous(
        MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined)),
        Ok(MulticenterBondConstraintsForm::new()))]
    #[case::contradiction(
        MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::ElectronCount(NumForm::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_multicenter_bond_constraints_form_normalize(
        #[case] constraints: MulticenterBondConstraintsForm,
        #[case] expected: Result<MulticenterBondConstraintsForm, Contradiction>,
    ) {
        assert_eq!(constraints.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::a_only_kept(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), MulticenterBondConstraintsForm::new(),
        Some(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6))))]
    #[case::b_only_kept(MulticenterBondConstraintsForm::new(), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)),
        Some(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6))))]
    #[case::shared_key_meets(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined)),
        Some(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6))))]
    #[case::shared_key_contradicts(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(2)), None)]
    #[case::prunes_vacuous(MulticenterBondConstraintsForm::new(), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined)), Some(MulticenterBondConstraintsForm::new()))]
    fn test_multicenter_bond_constraints_form_meet(#[case] a: MulticenterBondConstraintsForm, #[case] b: MulticenterBondConstraintsForm, #[case] expected: Option<MulticenterBondConstraintsForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::widens_value(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(1)), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(2)),
        MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::ElectronCount(NumForm::lit_set([1, 2]))))]
    #[case::single_side_dropped(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), MulticenterBondConstraintsForm::new(),
        MulticenterBondConstraintsForm::new())]
    #[case::undetermined_drops(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined)), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)),
        MulticenterBondConstraintsForm::new())]
    fn test_multicenter_bond_constraints_form_join(#[case] a: MulticenterBondConstraintsForm, #[case] b: MulticenterBondConstraintsForm, #[case] expected: MulticenterBondConstraintsForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(MulticenterBondConstraintsForm::new(), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), true)]
    #[case::required_present(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), true)]
    #[case::required_absent(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), MulticenterBondConstraintsForm::new(), false)]
    #[case::wildcard_matches_lit(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined)), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), true)]
    #[case::lit_mismatch(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(2)), false)]
    fn test_multicenter_bond_constraints_form_matches(
        #[case] pattern: MulticenterBondConstraintsForm,
        #[case] target: MulticenterBondConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_one_empty(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), MulticenterBondConstraintsForm::new(), true)]
    #[case::shared_key_compatible(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), true)]
    #[case::shared_key_incompatible(MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6)), MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(2)), false)]
    fn test_multicenter_bond_constraints_form_is_compatible(#[case] a: MulticenterBondConstraintsForm, #[case] b: MulticenterBondConstraintsForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(vec![MulticenterBondConstraintForm::electron_count(6)], vec![MulticenterBondConstraintForm::electron_count(6)])]
    #[case::same_key_last_wins(vec![MulticenterBondConstraintForm::electron_count(2), MulticenterBondConstraintForm::electron_count(6)], vec![MulticenterBondConstraintForm::electron_count(6)])]
    #[case::empty(vec![], vec![])]
    fn test_multicenter_bond_constraints_form_from_iter(
        #[case] input: Vec<MulticenterBondConstraintForm>,
        #[case] expected: Vec<MulticenterBondConstraintForm>,
    ) {
        let cs = MulticenterBondConstraintsForm::from_iter(input);
        assert_eq!(cs, MulticenterBondConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_multicenter_bond_constraints_form_into_iter() {
        let cs =
            MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6));
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![MulticenterBondConstraintForm::electron_count(6)]
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_form_from_multicenter_bond_constraint() {
        let cs: MulticenterBondConstraintsForm =
            MulticenterBondConstraintForm::electron_count(6).into();
        assert_eq!(
            cs,
            MulticenterBondConstraintsForm::from_iter([
                MulticenterBondConstraintForm::electron_count(6)
            ]),
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_form_from_vec() {
        let cs: MulticenterBondConstraintsForm =
            vec![MulticenterBondConstraintForm::electron_count(6)].into();
        assert_eq!(
            cs,
            MulticenterBondConstraintsForm::from_iter([
                MulticenterBondConstraintForm::electron_count(6)
            ]),
        );
    }
}
