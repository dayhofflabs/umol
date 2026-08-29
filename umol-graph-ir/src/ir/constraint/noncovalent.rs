//! Noncovalent bond constraints.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use super::super::boolean::BooleanForm;
use super::super::error::{Contradiction, NoJoin};
use super::super::remap::{IdRemapping, MoleculeCompaction};
use super::super::traits::{Lattice, Normalize};

/// Noncovalent-bond-scope constraint. Atom-ref and quantified-predicate forms
/// live at molecule scope via `RelationalConstraint`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NoncovalentBondConstraintForm {
    /// Whether the bond is intramolecular (`true`) or intermolecular (`false`);
    /// `Undetermined` when unspecified.
    Intramolecular(BooleanForm),
}

impl NoncovalentBondConstraintForm {
    pub fn intramolecular(b: impl Into<BooleanForm>) -> Self {
        Self::Intramolecular(b.into())
    }

    /// Noncovalent-bond constraint key, unique within a `NoncovalentBondConstraintsForm` container.
    pub fn key(&self) -> NoncovalentBondConstraintKey {
        match self {
            Self::Intramolecular(_) => NoncovalentBondConstraintKey::Intramolecular,
        }
    }

    /// Vacuous form of constraint key, used for removal.
    pub fn as_undetermined(&self) -> Self {
        match self {
            Self::Intramolecular(_) => Self::Intramolecular(BooleanForm::Undetermined),
        }
    }

    /// Value-only: no indices to compact.
    pub fn compact(self, _compaction: &MoleculeCompaction) -> Option<Self> {
        Some(self)
    }

    /// Value-only: no indices to remap.
    pub fn remap(self, _map: &IdRemapping) -> Self {
        self
    }
}

impl Normalize for NoncovalentBondConstraintForm {
    /// Normalize the inner value; the kind is preserved.
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Intramolecular(b) => Self::Intramolecular(b.normalize()?),
        })
    }
}

impl Lattice for NoncovalentBondConstraintForm {
    fn is_undetermined(&self) -> bool {
        match self {
            Self::Intramolecular(b) => b.is_undetermined(),
        }
    }

    fn is_ground(&self) -> bool {
        match self {
            Self::Intramolecular(b) => b.is_ground(),
        }
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Intramolecular(a), Self::Intramolecular(b)) => {
                a.meet(b).map(Self::Intramolecular)
            }
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        match (self, other) {
            (Self::Intramolecular(a), Self::Intramolecular(b)) => {
                Ok(Self::Intramolecular(a.join(b)?))
            }
        }
    }

    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Intramolecular(a), Self::Intramolecular(b)) => a.matches(b),
        }
    }

    fn is_compatible(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Intramolecular(a), Self::Intramolecular(b)) => a.is_compatible(b),
        }
    }
}

/// Entry identity: discriminant only (every kind is single-valued, no sub-key).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NoncovalentBondConstraintKey {
    Intramolecular,
}

/// Per-noncovalent-bond constraint container, ordered, unique by key, sorted flat vector storage.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NoncovalentBondConstraintsForm(Vec<NoncovalentBondConstraintForm>);

impl NoncovalentBondConstraintsForm {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The bond's intramolecular value, or `Undetermined` when no `Intramolecular` constraint is present.
    pub fn intramolecular(&self) -> BooleanForm {
        match self.get(NoncovalentBondConstraintKey::Intramolecular) {
            Some(NoncovalentBondConstraintForm::Intramolecular(b)) => *b,
            _ => BooleanForm::Undetermined,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    fn find(&self, key: NoncovalentBondConstraintKey) -> Result<usize, usize> {
        self.0.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains(&self, key: NoncovalentBondConstraintKey) -> bool {
        self.find(key).is_ok()
    }

    pub fn get(&self, key: NoncovalentBondConstraintKey) -> Option<&NoncovalentBondConstraintForm> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: NoncovalentBondConstraintForm) {
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
        old: Option<NoncovalentBondConstraintForm>,
        new: Option<NoncovalentBondConstraintForm>,
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
        key: NoncovalentBondConstraintKey,
    ) -> Option<NoncovalentBondConstraintForm> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = NoncovalentBondConstraintForm>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &NoncovalentBondConstraintsForm) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&NoncovalentBondConstraintForm) -> bool) {
        self.0.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl ExactSizeIterator<Item = NoncovalentBondConstraintForm> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, NoncovalentBondConstraintForm> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &MoleculeCompaction) -> Self {
        self
    }
}

impl Normalize for NoncovalentBondConstraintsForm {
    /// Normalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn normalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Normalize::normalize)
            .collect::<Result<Vec<NoncovalentBondConstraintForm>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl Lattice for NoncovalentBondConstraintsForm {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`NoncovalentBondConstraintForm::meet`; a `None` aborts the whole meet), an
    /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<NoncovalentBondConstraintForm> = Vec::new();
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
    /// (`NoncovalentBondConstraintForm::join`); a single-side key widens to the absent ⊤ and is dropped.
    /// The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<NoncovalentBondConstraintForm> = Vec::new();
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

    /// Pattern-driven: the intramolecular value is matched on its own lattice; an empty
    /// pattern matches any target.
    fn matches(&self, target: &Self) -> bool {
        self.intramolecular().matches(&target.intramolecular())
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

impl FromIterator<NoncovalentBondConstraintForm> for NoncovalentBondConstraintsForm {
    fn from_iter<I: IntoIterator<Item = NoncovalentBondConstraintForm>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for NoncovalentBondConstraintsForm {
    type Item = NoncovalentBondConstraintForm;
    type IntoIter = IntoIter<NoncovalentBondConstraintForm>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<NoncovalentBondConstraintForm> for NoncovalentBondConstraintsForm {
    fn from(c: NoncovalentBondConstraintForm) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<NoncovalentBondConstraintForm>> for NoncovalentBondConstraintsForm {
    fn from(cs: Vec<NoncovalentBondConstraintForm>) -> Self {
        Self::from_iter(cs)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::GraphCompaction;

    use super::*;

    #[rstest]
    #[case::intramolecular(
        NoncovalentBondConstraintForm::intramolecular(true),
        NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Lit(true))
    )]
    fn test_noncovalent_bond_constraint_form_constructors(
        #[case] actual: NoncovalentBondConstraintForm,
        #[case] expected: NoncovalentBondConstraintForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::intramolecular(
        NoncovalentBondConstraintForm::intramolecular(true),
        NoncovalentBondConstraintKey::Intramolecular
    )]
    fn test_noncovalent_bond_constraint_form_key(
        #[case] c: NoncovalentBondConstraintForm,
        #[case] expected: NoncovalentBondConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rstest]
    #[case::intramolecular(
        NoncovalentBondConstraintForm::intramolecular(true),
        NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined)
    )]
    fn test_noncovalent_bond_constraint_form_as_undetermined(
        #[case] c: NoncovalentBondConstraintForm,
        #[case] expected: NoncovalentBondConstraintForm,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rstest]
    #[case::intramolecular(
        NoncovalentBondConstraintForm::intramolecular(true),
        Ok(NoncovalentBondConstraintForm::intramolecular(true))
    )]
    fn test_noncovalent_bond_constraint_form_normalize(
        #[case] constraint: NoncovalentBondConstraintForm,
        #[case] expected: Result<NoncovalentBondConstraintForm, Contradiction>,
    ) {
        assert_eq!(constraint.normalize(), expected);
    }

    #[rstest]
    #[case::lit(NoncovalentBondConstraintForm::intramolecular(true), false)]
    #[case::undetermined(
        NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined),
        true
    )]
    fn test_noncovalent_bond_constraint_form_is_undetermined(
        #[case] c: NoncovalentBondConstraintForm,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::narrows_undetermined(NoncovalentBondConstraintForm::intramolecular(true), NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined), Some(NoncovalentBondConstraintForm::intramolecular(true)))]
    #[case::same_value(NoncovalentBondConstraintForm::intramolecular(true), NoncovalentBondConstraintForm::intramolecular(true), Some(NoncovalentBondConstraintForm::intramolecular(true)))]
    #[case::incompatible(NoncovalentBondConstraintForm::intramolecular(true), NoncovalentBondConstraintForm::intramolecular(false), None)]
    fn test_noncovalent_bond_constraint_form_meet(#[case] a: NoncovalentBondConstraintForm, #[case] b: NoncovalentBondConstraintForm, #[case] expected: Option<NoncovalentBondConstraintForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(NoncovalentBondConstraintForm::intramolecular(true), NoncovalentBondConstraintForm::intramolecular(true), Ok(NoncovalentBondConstraintForm::intramolecular(true)))]
    #[case::differ_widens(NoncovalentBondConstraintForm::intramolecular(true), NoncovalentBondConstraintForm::intramolecular(false), Ok(NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined)))]
    fn test_noncovalent_bond_constraint_form_join(#[case] a: NoncovalentBondConstraintForm, #[case] b: NoncovalentBondConstraintForm, #[case] expected: Result<NoncovalentBondConstraintForm, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(NoncovalentBondConstraintForm::intramolecular(true), NoncovalentBondConstraintForm::intramolecular(true), true)]
    #[case::incompatible(NoncovalentBondConstraintForm::intramolecular(true), NoncovalentBondConstraintForm::intramolecular(false), false)]
    fn test_noncovalent_bond_constraint_form_is_compatible(#[case] a: NoncovalentBondConstraintForm, #[case] b: NoncovalentBondConstraintForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_new() {
        let cs = NoncovalentBondConstraintsForm::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    #[case::present(
        NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)),
        BooleanForm::Lit(true)
    )]
    #[case::absent(NoncovalentBondConstraintsForm::new(), BooleanForm::Undetermined)]
    fn test_noncovalent_bond_constraints_form_intramolecular(
        #[case] cs: NoncovalentBondConstraintsForm,
        #[case] expected: BooleanForm,
    ) {
        assert_eq!(cs.intramolecular(), expected);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_iter() {
        let cs = NoncovalentBondConstraintsForm::from(
            NoncovalentBondConstraintForm::intramolecular(true),
        );
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![NoncovalentBondConstraintForm::intramolecular(true)]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![NoncovalentBondConstraintForm::intramolecular(true)], vec![NoncovalentBondConstraintForm::intramolecular(true)])]
    #[case::overwrite_same_key(vec![NoncovalentBondConstraintForm::intramolecular(true), NoncovalentBondConstraintForm::intramolecular(false)], vec![NoncovalentBondConstraintForm::intramolecular(false)])]
    #[case::vacuous_stores(vec![NoncovalentBondConstraintForm::intramolecular(true), NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined)], vec![NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined)])]
    fn test_noncovalent_bond_constraints_form_set(#[case] sequence: Vec<NoncovalentBondConstraintForm>, #[case] expected: Vec<NoncovalentBondConstraintForm>) {
        let mut cs = NoncovalentBondConstraintsForm::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, NoncovalentBondConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite(vec![NoncovalentBondConstraintForm::intramolecular(true)], vec![NoncovalentBondConstraintForm::intramolecular(false)], vec![NoncovalentBondConstraintForm::intramolecular(false)])]
    #[case::adds_from_empty(vec![], vec![NoncovalentBondConstraintForm::intramolecular(true)], vec![NoncovalentBondConstraintForm::intramolecular(true)])]
    #[case::vacuous_removes(vec![NoncovalentBondConstraintForm::intramolecular(true)], vec![NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined)], vec![])]
    fn test_noncovalent_bond_constraints_form_update(#[case] initial: Vec<NoncovalentBondConstraintForm>, #[case] other: Vec<NoncovalentBondConstraintForm>, #[case] expected: Vec<NoncovalentBondConstraintForm>) {
        let mut cs = NoncovalentBondConstraintsForm::from_iter(initial);
        cs.update(&NoncovalentBondConstraintsForm::from_iter(other));
        assert_eq!(cs, NoncovalentBondConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![NoncovalentBondConstraintForm::intramolecular(true)], Some(NoncovalentBondConstraintForm::intramolecular(true)), Some(NoncovalentBondConstraintForm::intramolecular(false)), Ok(()), vec![NoncovalentBondConstraintForm::intramolecular(false)])]
    #[case::remove(vec![NoncovalentBondConstraintForm::intramolecular(true)], Some(NoncovalentBondConstraintForm::intramolecular(true)), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(NoncovalentBondConstraintForm::intramolecular(true)), Ok(()), vec![NoncovalentBondConstraintForm::intramolecular(true)])]
    #[case::old_mismatch(vec![NoncovalentBondConstraintForm::intramolecular(true)], Some(NoncovalentBondConstraintForm::intramolecular(false)), None, Err(Contradiction), vec![NoncovalentBondConstraintForm::intramolecular(true)])]
    fn test_noncovalent_bond_constraints_form_compare_and_set(
        #[case] initial: Vec<NoncovalentBondConstraintForm>,
        #[case] old: Option<NoncovalentBondConstraintForm>,
        #[case] new: Option<NoncovalentBondConstraintForm>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<NoncovalentBondConstraintForm>,
    ) {
        let mut cs = NoncovalentBondConstraintsForm::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, NoncovalentBondConstraintsForm::from_iter(expected_state));
    }

    #[rstest]
    #[case::present(
        NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)),
        true
    )]
    #[case::absent(NoncovalentBondConstraintsForm::new(), false)]
    fn test_noncovalent_bond_constraints_form_contains(
        #[case] cs: NoncovalentBondConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(
            cs.contains(NoncovalentBondConstraintKey::Intramolecular),
            expected
        );
    }

    #[rstest]
    #[case::present(
        NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)),
        Some(NoncovalentBondConstraintForm::intramolecular(true))
    )]
    #[case::absent(NoncovalentBondConstraintsForm::new(), None)]
    fn test_noncovalent_bond_constraints_form_get(
        #[case] cs: NoncovalentBondConstraintsForm,
        #[case] expected: Option<NoncovalentBondConstraintForm>,
    ) {
        assert_eq!(
            cs.get(NoncovalentBondConstraintKey::Intramolecular),
            expected.as_ref()
        );
    }

    #[rstest]
    #[case::present(
        NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)),
        Some(NoncovalentBondConstraintForm::intramolecular(true)),
        NoncovalentBondConstraintsForm::new()
    )]
    #[case::absent(
        NoncovalentBondConstraintsForm::new(),
        None,
        NoncovalentBondConstraintsForm::new()
    )]
    fn test_noncovalent_bond_constraints_form_remove(
        #[case] mut cs: NoncovalentBondConstraintsForm,
        #[case] expected_removed: Option<NoncovalentBondConstraintForm>,
        #[case] expected_state: NoncovalentBondConstraintsForm,
    ) {
        assert_eq!(
            cs.remove(NoncovalentBondConstraintKey::Intramolecular),
            expected_removed
        );
        assert_eq!(cs, expected_state);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &NoncovalentBondConstraintForm| matches!(c, NoncovalentBondConstraintForm::Intramolecular(_)), vec![NoncovalentBondConstraintForm::intramolecular(true)])]
    #[case::all_dropped(|_: &NoncovalentBondConstraintForm| false, vec![])]
    fn test_noncovalent_bond_constraints_form_retain(
        #[case] predicate: impl FnMut(&NoncovalentBondConstraintForm) -> bool,
        #[case] expected: Vec<NoncovalentBondConstraintForm>,
    ) {
        let mut cs = NoncovalentBondConstraintsForm::from(
            NoncovalentBondConstraintForm::intramolecular(true),
        );
        cs.retain(predicate);
        assert_eq!(cs, NoncovalentBondConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_clear() {
        let mut cs = NoncovalentBondConstraintsForm::from(
            NoncovalentBondConstraintForm::intramolecular(true),
        );
        cs.clear();
        assert_eq!(cs, NoncovalentBondConstraintsForm::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_take() {
        let mut empty = NoncovalentBondConstraintsForm::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let mut cs = NoncovalentBondConstraintsForm::from(
            NoncovalentBondConstraintForm::intramolecular(true),
        );
        let mut taken = cs.take();
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(
            taken.next(),
            Some(NoncovalentBondConstraintForm::intramolecular(true)),
        );
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.size_hint(), (0, Some(0)));
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(cs, NoncovalentBondConstraintsForm::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_compact() {
        let cs = NoncovalentBondConstraintsForm::from(
            NoncovalentBondConstraintForm::intramolecular(true),
        );
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(Vec::new(), Vec::new()),
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
        NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined)),
        Ok(NoncovalentBondConstraintsForm::new()))]
    #[case::keeps_lit(
        NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)),
        Ok(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true))))]
    fn test_noncovalent_bond_constraints_form_normalize(
        #[case] constraints: NoncovalentBondConstraintsForm,
        #[case] expected: Result<NoncovalentBondConstraintsForm, Contradiction>,
    ) {
        assert_eq!(constraints.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::a_only_kept(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::new(),
        Some(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true))))]
    #[case::b_only_kept(NoncovalentBondConstraintsForm::new(), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)),
        Some(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true))))]
    #[case::shared_key_meets(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined)),
        Some(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true))))]
    #[case::shared_key_contradicts(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(false)), None)]
    fn test_noncovalent_bond_constraints_form_meet(#[case] a: NoncovalentBondConstraintsForm, #[case] b: NoncovalentBondConstraintsForm, #[case] expected: Option<NoncovalentBondConstraintsForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::differ_widens_to_undetermined(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(false)),
        NoncovalentBondConstraintsForm::new())]
    #[case::single_side_dropped(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::new(),
        NoncovalentBondConstraintsForm::new())]
    #[case::shared_same_kept(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)),
        NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)))]
    fn test_noncovalent_bond_constraints_form_join(#[case] a: NoncovalentBondConstraintsForm, #[case] b: NoncovalentBondConstraintsForm, #[case] expected: NoncovalentBondConstraintsForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(NoncovalentBondConstraintsForm::new(), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), true)]
    #[case::required_present(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), true)]
    #[case::required_absent(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::new(), false)]
    #[case::wildcard_matches_lit(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined)), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), true)]
    #[case::lit_mismatch(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(false)), false)]
    fn test_noncovalent_bond_constraints_form_matches(
        #[case] pattern: NoncovalentBondConstraintsForm,
        #[case] target: NoncovalentBondConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_one_empty(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::new(), true)]
    #[case::shared_key_compatible(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), true)]
    #[case::shared_key_incompatible(NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(false)), false)]
    fn test_noncovalent_bond_constraints_form_is_compatible(#[case] a: NoncovalentBondConstraintsForm, #[case] b: NoncovalentBondConstraintsForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(vec![NoncovalentBondConstraintForm::intramolecular(true)], vec![NoncovalentBondConstraintForm::intramolecular(true)])]
    #[case::same_key_last_wins(vec![NoncovalentBondConstraintForm::intramolecular(false), NoncovalentBondConstraintForm::intramolecular(true)], vec![NoncovalentBondConstraintForm::intramolecular(true)])]
    #[case::empty(vec![], vec![])]
    fn test_noncovalent_bond_constraints_form_from_iter(
        #[case] input: Vec<NoncovalentBondConstraintForm>,
        #[case] expected: Vec<NoncovalentBondConstraintForm>,
    ) {
        let cs = NoncovalentBondConstraintsForm::from_iter(input);
        assert_eq!(cs, NoncovalentBondConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_into_iter() {
        let cs = NoncovalentBondConstraintsForm::from(
            NoncovalentBondConstraintForm::intramolecular(true),
        );
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![NoncovalentBondConstraintForm::intramolecular(true)]
        );
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_from_noncovalent_bond_constraint() {
        let cs: NoncovalentBondConstraintsForm =
            NoncovalentBondConstraintForm::intramolecular(true).into();
        assert_eq!(
            cs,
            NoncovalentBondConstraintsForm::from_iter([
                NoncovalentBondConstraintForm::intramolecular(true)
            ]),
        );
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_from_vec() {
        let cs: NoncovalentBondConstraintsForm =
            vec![NoncovalentBondConstraintForm::intramolecular(true)].into();
        assert_eq!(
            cs,
            NoncovalentBondConstraintsForm::from_iter([
                NoncovalentBondConstraintForm::intramolecular(true)
            ]),
        );
    }
}
