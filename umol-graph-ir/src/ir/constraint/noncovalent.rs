//! Noncovalent bond constraints.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use super::super::boolean::BooleanForm;
use super::super::error::{Contradiction, NoJoin};
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::traits::{Canonicalize, Lattice};

/// Noncovalent-bond-scope constraint. Atom-ref and quantified-predicate forms
/// live at molecule scope via `RelationalConstraint`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NoncovalentBondConstraintAst {
    /// Whether the bond is intramolecular (`true`) or intermolecular (`false`);
    /// `Undetermined` when unspecified.
    Intramolecular(BooleanForm),
}

impl NoncovalentBondConstraintAst {
    pub fn intramolecular(b: impl Into<BooleanForm>) -> Self {
        Self::Intramolecular(b.into())
    }

    /// Noncovalent-bond constraint key, unique within a `NoncovalentBondConstraintsAst` container.
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
    pub fn compact(self, _compaction: &IdCompaction) -> Option<Self> {
        Some(self)
    }

    /// Value-only: no indices to remap.
    pub fn remap(self, _map: &IdRemapping) -> Self {
        self
    }
}

impl Canonicalize for NoncovalentBondConstraintAst {
    /// Canonicalize the inner value; the kind is preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Intramolecular(b) => Self::Intramolecular(b.canonicalize()?),
        })
    }
}

impl Lattice for NoncovalentBondConstraintAst {
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
pub struct NoncovalentBondConstraintsAst(Vec<NoncovalentBondConstraintAst>);

impl NoncovalentBondConstraintsAst {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// The bond's intramolecular value, or `Undetermined` when no `Intramolecular` constraint is present.
    pub fn intramolecular(&self) -> BooleanForm {
        match self.get(NoncovalentBondConstraintKey::Intramolecular) {
            Some(NoncovalentBondConstraintAst::Intramolecular(b)) => *b,
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

    pub fn get(&self, key: NoncovalentBondConstraintKey) -> Option<&NoncovalentBondConstraintAst> {
        self.find(key).ok().map(|i| &self.0[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: NoncovalentBondConstraintAst) {
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
        old: Option<NoncovalentBondConstraintAst>,
        new: Option<NoncovalentBondConstraintAst>,
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
        key: NoncovalentBondConstraintKey,
    ) -> Option<NoncovalentBondConstraintAst> {
        self.find(key).ok().map(|i| self.0.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = NoncovalentBondConstraintAst>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &NoncovalentBondConstraintsAst) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&NoncovalentBondConstraintAst) -> bool) {
        self.0.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl ExactSizeIterator<Item = NoncovalentBondConstraintAst> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn iter(&self) -> Iter<'_, NoncovalentBondConstraintAst> {
        self.0.iter()
    }

    pub fn compact(self, _compaction: &IdCompaction) -> Self {
        self
    }
}

impl Canonicalize for NoncovalentBondConstraintsAst {
    /// Canonicalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .0
            .into_iter()
            .map(Canonicalize::canonicalize)
            .collect::<Result<Vec<NoncovalentBondConstraintAst>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self(entries))
    }
}

impl Lattice for NoncovalentBondConstraintsAst {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`NoncovalentBondConstraintAst::meet`; a `None` aborts the whole meet), an
    /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: Vec<NoncovalentBondConstraintAst> = Vec::new();
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
    /// (`NoncovalentBondConstraintAst::join`); a single-side key widens to the absent ⊤ and is dropped.
    /// The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: Vec<NoncovalentBondConstraintAst> = Vec::new();
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

impl FromIterator<NoncovalentBondConstraintAst> for NoncovalentBondConstraintsAst {
    fn from_iter<I: IntoIterator<Item = NoncovalentBondConstraintAst>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for NoncovalentBondConstraintsAst {
    type Item = NoncovalentBondConstraintAst;
    type IntoIter = IntoIter<NoncovalentBondConstraintAst>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<NoncovalentBondConstraintAst> for NoncovalentBondConstraintsAst {
    fn from(c: NoncovalentBondConstraintAst) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<NoncovalentBondConstraintAst>> for NoncovalentBondConstraintsAst {
    fn from(cs: Vec<NoncovalentBondConstraintAst>) -> Self {
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
    #[case::intramolecular(
        NoncovalentBondConstraintAst::intramolecular(true),
        NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Lit(true))
    )]
    fn test_noncovalent_bond_constraint_ast_constructors(
        #[case] actual: NoncovalentBondConstraintAst,
        #[case] expected: NoncovalentBondConstraintAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::intramolecular(
        NoncovalentBondConstraintAst::intramolecular(true),
        NoncovalentBondConstraintKey::Intramolecular
    )]
    fn test_noncovalent_bond_constraint_ast_key(
        #[case] c: NoncovalentBondConstraintAst,
        #[case] expected: NoncovalentBondConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rstest]
    #[case::intramolecular(
        NoncovalentBondConstraintAst::intramolecular(true),
        NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined)
    )]
    fn test_noncovalent_bond_constraint_ast_as_undetermined(
        #[case] c: NoncovalentBondConstraintAst,
        #[case] expected: NoncovalentBondConstraintAst,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rstest]
    #[case::intramolecular(
        NoncovalentBondConstraintAst::intramolecular(true),
        Ok(NoncovalentBondConstraintAst::intramolecular(true))
    )]
    fn test_noncovalent_bond_constraint_ast_canonicalize(
        #[case] constraint: NoncovalentBondConstraintAst,
        #[case] expected: Result<NoncovalentBondConstraintAst, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
    }

    #[rstest]
    #[case::lit(NoncovalentBondConstraintAst::intramolecular(true), false)]
    #[case::undetermined(
        NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined),
        true
    )]
    fn test_noncovalent_bond_constraint_ast_is_undetermined(
        #[case] c: NoncovalentBondConstraintAst,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::narrows_undetermined(NoncovalentBondConstraintAst::intramolecular(true), NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined), Some(NoncovalentBondConstraintAst::intramolecular(true)))]
    #[case::same_value(NoncovalentBondConstraintAst::intramolecular(true), NoncovalentBondConstraintAst::intramolecular(true), Some(NoncovalentBondConstraintAst::intramolecular(true)))]
    #[case::incompatible(NoncovalentBondConstraintAst::intramolecular(true), NoncovalentBondConstraintAst::intramolecular(false), None)]
    fn test_noncovalent_bond_constraint_ast_meet(#[case] a: NoncovalentBondConstraintAst, #[case] b: NoncovalentBondConstraintAst, #[case] expected: Option<NoncovalentBondConstraintAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(NoncovalentBondConstraintAst::intramolecular(true), NoncovalentBondConstraintAst::intramolecular(true), Ok(NoncovalentBondConstraintAst::intramolecular(true)))]
    #[case::differ_widens(NoncovalentBondConstraintAst::intramolecular(true), NoncovalentBondConstraintAst::intramolecular(false), Ok(NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined)))]
    fn test_noncovalent_bond_constraint_ast_join(#[case] a: NoncovalentBondConstraintAst, #[case] b: NoncovalentBondConstraintAst, #[case] expected: Result<NoncovalentBondConstraintAst, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_value(NoncovalentBondConstraintAst::intramolecular(true), NoncovalentBondConstraintAst::intramolecular(true), true)]
    #[case::incompatible(NoncovalentBondConstraintAst::intramolecular(true), NoncovalentBondConstraintAst::intramolecular(false), false)]
    fn test_noncovalent_bond_constraint_ast_is_compatible(#[case] a: NoncovalentBondConstraintAst, #[case] b: NoncovalentBondConstraintAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_new() {
        let cs = NoncovalentBondConstraintsAst::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    #[case::present(
        NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)),
        BooleanForm::Lit(true)
    )]
    #[case::absent(NoncovalentBondConstraintsAst::new(), BooleanForm::Undetermined)]
    fn test_noncovalent_bond_constraints_ast_intramolecular(
        #[case] cs: NoncovalentBondConstraintsAst,
        #[case] expected: BooleanForm,
    ) {
        assert_eq!(cs.intramolecular(), expected);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_iter() {
        let cs =
            NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![NoncovalentBondConstraintAst::intramolecular(true)]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![NoncovalentBondConstraintAst::intramolecular(true)], vec![NoncovalentBondConstraintAst::intramolecular(true)])]
    #[case::overwrite_same_key(vec![NoncovalentBondConstraintAst::intramolecular(true), NoncovalentBondConstraintAst::intramolecular(false)], vec![NoncovalentBondConstraintAst::intramolecular(false)])]
    #[case::vacuous_stores(vec![NoncovalentBondConstraintAst::intramolecular(true), NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined)], vec![NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined)])]
    fn test_noncovalent_bond_constraints_ast_set(#[case] sequence: Vec<NoncovalentBondConstraintAst>, #[case] expected: Vec<NoncovalentBondConstraintAst>) {
        let mut cs = NoncovalentBondConstraintsAst::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, NoncovalentBondConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite(vec![NoncovalentBondConstraintAst::intramolecular(true)], vec![NoncovalentBondConstraintAst::intramolecular(false)], vec![NoncovalentBondConstraintAst::intramolecular(false)])]
    #[case::adds_from_empty(vec![], vec![NoncovalentBondConstraintAst::intramolecular(true)], vec![NoncovalentBondConstraintAst::intramolecular(true)])]
    #[case::vacuous_removes(vec![NoncovalentBondConstraintAst::intramolecular(true)], vec![NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined)], vec![])]
    fn test_noncovalent_bond_constraints_ast_update(#[case] initial: Vec<NoncovalentBondConstraintAst>, #[case] other: Vec<NoncovalentBondConstraintAst>, #[case] expected: Vec<NoncovalentBondConstraintAst>) {
        let mut cs = NoncovalentBondConstraintsAst::from_iter(initial);
        cs.update(&NoncovalentBondConstraintsAst::from_iter(other));
        assert_eq!(cs, NoncovalentBondConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![NoncovalentBondConstraintAst::intramolecular(true)], Some(NoncovalentBondConstraintAst::intramolecular(true)), Some(NoncovalentBondConstraintAst::intramolecular(false)), Ok(()), vec![NoncovalentBondConstraintAst::intramolecular(false)])]
    #[case::remove(vec![NoncovalentBondConstraintAst::intramolecular(true)], Some(NoncovalentBondConstraintAst::intramolecular(true)), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(NoncovalentBondConstraintAst::intramolecular(true)), Ok(()), vec![NoncovalentBondConstraintAst::intramolecular(true)])]
    #[case::old_mismatch(vec![NoncovalentBondConstraintAst::intramolecular(true)], Some(NoncovalentBondConstraintAst::intramolecular(false)), None, Err(Contradiction), vec![NoncovalentBondConstraintAst::intramolecular(true)])]
    fn test_noncovalent_bond_constraints_ast_compare_and_set(
        #[case] initial: Vec<NoncovalentBondConstraintAst>,
        #[case] old: Option<NoncovalentBondConstraintAst>,
        #[case] new: Option<NoncovalentBondConstraintAst>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<NoncovalentBondConstraintAst>,
    ) {
        let mut cs = NoncovalentBondConstraintsAst::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, NoncovalentBondConstraintsAst::from_iter(expected_state));
    }

    #[rstest]
    #[case::present(
        NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)),
        true
    )]
    #[case::absent(NoncovalentBondConstraintsAst::new(), false)]
    fn test_noncovalent_bond_constraints_ast_contains(
        #[case] cs: NoncovalentBondConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(
            cs.contains(NoncovalentBondConstraintKey::Intramolecular),
            expected
        );
    }

    #[rstest]
    #[case::present(
        NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)),
        Some(NoncovalentBondConstraintAst::intramolecular(true))
    )]
    #[case::absent(NoncovalentBondConstraintsAst::new(), None)]
    fn test_noncovalent_bond_constraints_ast_get(
        #[case] cs: NoncovalentBondConstraintsAst,
        #[case] expected: Option<NoncovalentBondConstraintAst>,
    ) {
        assert_eq!(
            cs.get(NoncovalentBondConstraintKey::Intramolecular),
            expected.as_ref()
        );
    }

    #[rstest]
    #[case::present(
        NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)),
        Some(NoncovalentBondConstraintAst::intramolecular(true)),
        NoncovalentBondConstraintsAst::new()
    )]
    #[case::absent(
        NoncovalentBondConstraintsAst::new(),
        None,
        NoncovalentBondConstraintsAst::new()
    )]
    fn test_noncovalent_bond_constraints_ast_remove(
        #[case] mut cs: NoncovalentBondConstraintsAst,
        #[case] expected_removed: Option<NoncovalentBondConstraintAst>,
        #[case] expected_state: NoncovalentBondConstraintsAst,
    ) {
        assert_eq!(
            cs.remove(NoncovalentBondConstraintKey::Intramolecular),
            expected_removed
        );
        assert_eq!(cs, expected_state);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &NoncovalentBondConstraintAst| matches!(c, NoncovalentBondConstraintAst::Intramolecular(_)), vec![NoncovalentBondConstraintAst::intramolecular(true)])]
    #[case::all_dropped(|_: &NoncovalentBondConstraintAst| false, vec![])]
    fn test_noncovalent_bond_constraints_ast_retain(
        #[case] predicate: impl FnMut(&NoncovalentBondConstraintAst) -> bool,
        #[case] expected: Vec<NoncovalentBondConstraintAst>,
    ) {
        let mut cs = NoncovalentBondConstraintsAst::from(
            NoncovalentBondConstraintAst::intramolecular(true),
        );
        cs.retain(predicate);
        assert_eq!(cs, NoncovalentBondConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_clear() {
        let mut cs =
            NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true));
        cs.clear();
        assert_eq!(cs, NoncovalentBondConstraintsAst::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_take() {
        let mut empty = NoncovalentBondConstraintsAst::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let mut cs =
            NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true));
        let mut taken = cs.take();
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(
            taken.next(),
            Some(NoncovalentBondConstraintAst::intramolecular(true)),
        );
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.size_hint(), (0, Some(0)));
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(cs, NoncovalentBondConstraintsAst::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_compact() {
        let cs =
            NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true));
        let compaction = IdCompaction::new(
            Compaction::new(Vec::new(), Vec::new()),
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
        NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined)),
        Ok(NoncovalentBondConstraintsAst::new()))]
    #[case::keeps_lit(
        NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)),
        Ok(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true))))]
    fn test_noncovalent_bond_constraints_ast_canonicalize(
        #[case] constraints: NoncovalentBondConstraintsAst,
        #[case] expected: Result<NoncovalentBondConstraintsAst, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::a_only_kept(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::new(),
        Some(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true))))]
    #[case::b_only_kept(NoncovalentBondConstraintsAst::new(), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)),
        Some(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true))))]
    #[case::shared_key_meets(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined)),
        Some(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true))))]
    #[case::shared_key_contradicts(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(false)), None)]
    fn test_noncovalent_bond_constraints_ast_meet(#[case] a: NoncovalentBondConstraintsAst, #[case] b: NoncovalentBondConstraintsAst, #[case] expected: Option<NoncovalentBondConstraintsAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::differ_widens_to_undetermined(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(false)),
        NoncovalentBondConstraintsAst::new())]
    #[case::single_side_dropped(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::new(),
        NoncovalentBondConstraintsAst::new())]
    #[case::shared_same_kept(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)),
        NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)))]
    fn test_noncovalent_bond_constraints_ast_join(#[case] a: NoncovalentBondConstraintsAst, #[case] b: NoncovalentBondConstraintsAst, #[case] expected: NoncovalentBondConstraintsAst) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_anything(NoncovalentBondConstraintsAst::new(), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), true)]
    #[case::required_present(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), true)]
    #[case::required_absent(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::new(), false)]
    #[case::wildcard_matches_lit(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined)), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), true)]
    #[case::lit_mismatch(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(false)), false)]
    fn test_noncovalent_bond_constraints_ast_matches(
        #[case] pattern: NoncovalentBondConstraintsAst,
        #[case] target: NoncovalentBondConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_one_empty(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::new(), true)]
    #[case::shared_key_compatible(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), true)]
    #[case::shared_key_incompatible(NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(false)), false)]
    fn test_noncovalent_bond_constraints_ast_is_compatible(#[case] a: NoncovalentBondConstraintsAst, #[case] b: NoncovalentBondConstraintsAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(vec![NoncovalentBondConstraintAst::intramolecular(true)], vec![NoncovalentBondConstraintAst::intramolecular(true)])]
    #[case::same_key_last_wins(vec![NoncovalentBondConstraintAst::intramolecular(false), NoncovalentBondConstraintAst::intramolecular(true)], vec![NoncovalentBondConstraintAst::intramolecular(true)])]
    #[case::empty(vec![], vec![])]
    fn test_noncovalent_bond_constraints_ast_from_iter(
        #[case] input: Vec<NoncovalentBondConstraintAst>,
        #[case] expected: Vec<NoncovalentBondConstraintAst>,
    ) {
        let cs = NoncovalentBondConstraintsAst::from_iter(input);
        assert_eq!(cs, NoncovalentBondConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_into_iter() {
        let cs =
            NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true));
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![NoncovalentBondConstraintAst::intramolecular(true)]
        );
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_from_noncovalent_bond_constraint() {
        let cs: NoncovalentBondConstraintsAst =
            NoncovalentBondConstraintAst::intramolecular(true).into();
        assert_eq!(
            cs,
            NoncovalentBondConstraintsAst::from_iter([
                NoncovalentBondConstraintAst::intramolecular(true)
            ]),
        );
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_ast_from_vec() {
        let cs: NoncovalentBondConstraintsAst =
            vec![NoncovalentBondConstraintAst::intramolecular(true)].into();
        assert_eq!(
            cs,
            NoncovalentBondConstraintsAst::from_iter([
                NoncovalentBondConstraintAst::intramolecular(true)
            ]),
        );
    }
}
