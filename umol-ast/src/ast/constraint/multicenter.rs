//! Multicenter bond constraints.

use std::mem::{self, replace};
use std::slice::Iter;
use std::vec::IntoIter;

use strum::EnumDiscriminants;

use super::super::error::Contradiction;
use super::super::remap::IdRemapping;
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

    /// Entry identity for order/dedup. Every kind is single-valued, so no sub-key.
    pub fn key(&self) -> MulticenterBondConstraintKey {
        match self {
            Self::ElectronCount(_) => MulticenterBondConstraintKey::ElectronCount,
        }
    }

    /// Every `MulticenterBondConstraint` variant is single-valued per bond.
    pub fn is_unique(&self) -> bool {
        true
    }

    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::ElectronCount(v) => v.is_undetermined(),
        }
    }

    pub fn remap(self, _remap: &IdRemapping) -> Option<Self> {
        // Value-only: no indices to remap.
        Some(self)
    }
}

/// Entry identity: discriminant only (every kind is single-valued, no sub-key).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MulticenterBondConstraintKey {
    ElectronCount,
}

impl MulticenterBondConstraintKey {
    pub fn kind(self) -> MulticenterBondConstraintKind {
        match self {
            Self::ElectronCount => MulticenterBondConstraintKind::ElectronCount,
        }
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

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterBondConstraints(Vec<MulticenterBondConstraint>);

impl MulticenterBondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[MulticenterBondConstraint] {
        &self.0
    }

    pub fn contains(&self, kind: MulticenterBondConstraintKind) -> bool {
        self.0.iter().any(|c| c.kind() == kind)
    }

    pub fn get(&self, kind: MulticenterBondConstraintKind) -> Option<&MulticenterBondConstraint> {
        self.0.iter().find(|c| c.kind() == kind)
    }

    pub fn get_mut(
        &mut self,
        kind: MulticenterBondConstraintKind,
    ) -> Option<&mut MulticenterBondConstraint> {
        self.0.iter_mut().find(|c| c.kind() == kind)
    }

    pub fn electron_count(&self) -> ValueAst {
        match self.get(MulticenterBondConstraintKind::ElectronCount) {
            Some(MulticenterBondConstraint::ElectronCount(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn iter(&self) -> Iter<'_, MulticenterBondConstraint> {
        self.0.iter()
    }

    /// Insert a constraint per the per-variant cardinality policy. Returns
    /// the replaced entry if `c.is_unique()` and a same-kind entry already
    /// existed; `None` otherwise.
    /// Insert at the `key()`-sorted position: unique kinds replace the same-key
    /// entry (returning it); non-unique kinds append after the same-key run.
    pub fn add(&mut self, c: MulticenterBondConstraint) -> Option<MulticenterBondConstraint> {
        match self.find_by_key(c.key()) {
            Ok(i) if c.is_unique() => Some(replace(&mut self.0[i], c)),
            Ok(i) => {
                let end = i + self.0[i..]
                    .iter()
                    .take_while(|e| e.key() == c.key())
                    .count();
                self.0.insert(end, c);
                None
            }
            Err(i) => {
                self.0.insert(i, c);
                None
            }
        }
    }

    fn find_by_key(&self, key: MulticenterBondConstraintKey) -> Result<usize, usize> {
        self.0.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains_key(&self, key: MulticenterBondConstraintKey) -> bool {
        self.find_by_key(key).is_ok()
    }

    pub fn get_by_key(
        &self,
        key: MulticenterBondConstraintKey,
    ) -> Option<&MulticenterBondConstraint> {
        self.find_by_key(key).ok().map(|i| &self.0[i])
    }

    pub fn get_by_key_mut(
        &mut self,
        key: MulticenterBondConstraintKey,
    ) -> Option<&mut MulticenterBondConstraint> {
        self.find_by_key(key).ok().map(|i| &mut self.0[i])
    }

    pub fn remove_by_key(
        &mut self,
        key: MulticenterBondConstraintKey,
    ) -> Option<MulticenterBondConstraint> {
        self.find_by_key(key).ok().map(|i| self.0.remove(i))
    }

    /// Add multiple constraints at once, using semantics of `add`.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = MulticenterBondConstraint>) {
        for constraint in constraints {
            self.add(constraint);
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&MulticenterBondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    pub fn take(&mut self) -> impl Iterator<Item = MulticenterBondConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn remove(
        &mut self,
        kind: MulticenterBondConstraintKind,
    ) -> Option<MulticenterBondConstraint> {
        let pos = self.0.iter().position(|c| c.kind() == kind)?;
        Some(self.0.remove(pos))
    }

    /// Iterate over every entry of `kind`. Currently every variant is
    /// single-valued so this yields at most one entry.
    pub fn get_all(
        &self,
        kind: MulticenterBondConstraintKind,
    ) -> impl Iterator<Item = &MulticenterBondConstraint> {
        self.0.iter().filter(move |c| c.kind() == kind)
    }

    /// Remove every entry of `kind`, returning them in insertion order.
    pub fn remove_all(
        &mut self,
        kind: MulticenterBondConstraintKind,
    ) -> Vec<MulticenterBondConstraint> {
        let mut out = Vec::new();
        self.0.retain(|c| {
            if c.kind() == kind {
                out.push(c.clone());
                false
            } else {
                true
            }
        });
        out
    }

    pub fn remap(self, remap: &IdRemapping) -> Self {
        Self(self.0.into_iter().filter_map(|c| c.remap(remap)).collect())
    }
}

impl Canonicalize for MulticenterBondConstraints {
    /// Sort by `key()`, canonicalize each value, drop vacuous entries. No merge
    /// clause — every kind is single-valued, so `add` admits no same-key duplicates.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut entries = self.0;
        entries.sort_by_key(|c| c.key());
        let mut out: Vec<MulticenterBondConstraint> = Vec::with_capacity(entries.len());
        for c in entries {
            out.push(c.canonicalize()?);
        }
        out.retain(|c| !c.is_undetermined());
        Ok(Self(out))
    }
}

impl Lattice for MulticenterBondConstraints {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| match c {
            MulticenterBondConstraint::ElectronCount(v) => v.is_undetermined(),
        })
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| match c {
            MulticenterBondConstraint::ElectronCount(v) => v.is_ground(),
        })
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let mut result = Self::new();
        let merged = self.electron_count().meet(&other.electron_count())?;
        if !merged.is_undetermined() {
            result.add(MulticenterBondConstraint::ElectronCount(merged));
        }
        Some(result)
    }

    fn join(&self, other: &Self) -> Self {
        let mut result = Self::new();
        let a_has = self.contains(MulticenterBondConstraintKind::ElectronCount);
        let b_has = other.contains(MulticenterBondConstraintKind::ElectronCount);
        if a_has && b_has {
            let joined = self.electron_count().join(&other.electron_count());
            if !joined.is_undetermined() {
                result.add(MulticenterBondConstraint::ElectronCount(joined));
            }
        }
        result
    }

    fn matches(&self, target: &Self) -> bool {
        self.electron_count().matches(&target.electron_count())
    }
}

impl FromIterator<MulticenterBondConstraint> for MulticenterBondConstraints {
    fn from_iter<I: IntoIterator<Item = MulticenterBondConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
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
    use umol_graph_core::Remapping;

    use super::*;
    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraint::electron_count(2),
        MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2))
    )]
    fn test_multicenter_bond_constraint_constructors(
        #[case] actual: MulticenterBondConstraint,
        #[case] expected: MulticenterBondConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraint::electron_count(2),
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
        MulticenterBondConstraintKey::ElectronCount,
        MulticenterBondConstraintKind::ElectronCount
    )]
    fn test_multicenter_bond_constraint_key_kind(
        #[case] key: MulticenterBondConstraintKey,
        #[case] expected: MulticenterBondConstraintKind,
    ) {
        assert_eq!(key.kind(), expected);
    }

    #[rstest]
    #[case::electron_count(MulticenterBondConstraint::electron_count(2))]
    fn test_multicenter_bond_constraint_is_unique(#[case] c: MulticenterBondConstraint) {
        assert!(c.is_unique());
    }

    #[rstest]
    #[case::lit(MulticenterBondConstraint::electron_count(2), false)]
    #[case::undetermined(MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined), true)]
    fn test_multicenter_bond_constraint_is_undetermined(
        #[case] c: MulticenterBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
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
    fn test_multicenter_bond_constraints_new() {
        let cs = MulticenterBondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[MulticenterBondConstraint]);
    }

    #[rstest]
    #[case::present(MulticenterBondConstraintKind::ElectronCount, true)]
    fn test_multicenter_bond_constraints_contains(
        #[case] kind: MulticenterBondConstraintKind,
        #[case] expected: bool,
    ) {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        assert_eq!(cs.contains(kind), expected);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_contains_absent() {
        let cs = MulticenterBondConstraints::new();
        assert!(!cs.contains(MulticenterBondConstraintKind::ElectronCount));
    }

    #[rstest]
    fn test_multicenter_bond_constraints_get() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        assert_eq!(
            cs.get(MulticenterBondConstraintKind::ElectronCount),
            Some(&MulticenterBondConstraint::electron_count(2)),
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_get_absent() {
        let cs = MulticenterBondConstraints::new();
        assert_eq!(cs.get(MulticenterBondConstraintKind::ElectronCount), None);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_get_mut() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        let entry = cs
            .get_mut(MulticenterBondConstraintKind::ElectronCount)
            .unwrap();
        *entry = MulticenterBondConstraint::electron_count(4);
        assert_eq!(
            cs.as_slice(),
            &[MulticenterBondConstraint::electron_count(4)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_get_mut_absent() {
        let mut cs = MulticenterBondConstraints::new();
        assert!(cs
            .get_mut(MulticenterBondConstraintKind::ElectronCount)
            .is_none());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_contains_key() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        assert!(cs.contains_key(MulticenterBondConstraintKey::ElectronCount));
        assert!(!MulticenterBondConstraints::new()
            .contains_key(MulticenterBondConstraintKey::ElectronCount));
    }

    #[rstest]
    fn test_multicenter_bond_constraints_get_by_key() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        assert_eq!(
            cs.get_by_key(MulticenterBondConstraintKey::ElectronCount),
            Some(&MulticenterBondConstraint::electron_count(2)),
        );
        assert_eq!(
            MulticenterBondConstraints::new()
                .get_by_key(MulticenterBondConstraintKey::ElectronCount),
            None,
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_get_by_key_mut() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        *cs.get_by_key_mut(MulticenterBondConstraintKey::ElectronCount)
            .unwrap() = MulticenterBondConstraint::electron_count(4);
        assert_eq!(
            cs.get_by_key(MulticenterBondConstraintKey::ElectronCount),
            Some(&MulticenterBondConstraint::electron_count(4)),
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_remove_by_key() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        assert_eq!(
            cs.remove_by_key(MulticenterBondConstraintKey::ElectronCount),
            Some(MulticenterBondConstraint::electron_count(2)),
        );
        assert_eq!(cs.as_slice(), &[] as &[MulticenterBondConstraint]);
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

    #[rstest]
    fn test_multicenter_bond_constraints_iter() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![MulticenterBondConstraint::electron_count(2)]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(
        vec![MulticenterBondConstraint::electron_count(2)],
        vec![None],
        vec![MulticenterBondConstraint::electron_count(2)],
    )]
    #[case::replace_same_kind(
        vec![
            MulticenterBondConstraint::electron_count(2),
            MulticenterBondConstraint::electron_count(4),
        ],
        vec![None, Some(MulticenterBondConstraint::electron_count(2))],
        vec![MulticenterBondConstraint::electron_count(4)],
    )]
    fn test_multicenter_bond_constraints_add(
        #[case] sequence: Vec<MulticenterBondConstraint>,
        #[case] expected_returns: Vec<Option<MulticenterBondConstraint>>,
        #[case] expected_state: Vec<MulticenterBondConstraint>,
    ) {
        let mut cs = MulticenterBondConstraints::new();
        let returns: Vec<_> = sequence.into_iter().map(|c| cs.add(c)).collect();
        assert_eq!(returns, expected_returns);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(
        |c: &MulticenterBondConstraint| matches!(c, MulticenterBondConstraint::ElectronCount(_)),
        vec![MulticenterBondConstraint::electron_count(2)],
    )]
    #[case::all_dropped(|_: &MulticenterBondConstraint| false, vec![])]
    fn test_multicenter_bond_constraints_retain(
        #[case] predicate: impl FnMut(&MulticenterBondConstraint) -> bool,
        #[case] expected: Vec<MulticenterBondConstraint>,
    ) {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        cs.retain(predicate);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_clear() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        cs.clear();
        assert_eq!(cs, MulticenterBondConstraints::new());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_take() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(drained, vec![MulticenterBondConstraint::electron_count(2)]);
        assert_eq!(cs, MulticenterBondConstraints::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::present(
        MulticenterBondConstraintKind::ElectronCount,
        Some(MulticenterBondConstraint::electron_count(2)),
        vec![],
    )]
    fn test_multicenter_bond_constraints_remove(
        #[case] kind: MulticenterBondConstraintKind,
        #[case] expected_returned: Option<MulticenterBondConstraint>,
        #[case] expected_state: Vec<MulticenterBondConstraint>,
    ) {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        assert_eq!(cs.remove(kind), expected_returned);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_remove_absent() {
        let mut cs = MulticenterBondConstraints::new();
        assert_eq!(
            cs.remove(MulticenterBondConstraintKind::ElectronCount),
            None,
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_remap() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        let remap = IdRemapping::new(
            Remapping::new(vec![0, 1], vec![0]),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(cs.clone().remap(&remap), cs);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(vec![MulticenterBondConstraint::electron_count(2)], vec![MulticenterBondConstraint::electron_count(2)])]
    #[case::same_kind_last_wins(vec![MulticenterBondConstraint::electron_count(2), MulticenterBondConstraint::electron_count(4)],
        vec![MulticenterBondConstraint::electron_count(4)])]
    #[case::empty(vec![], vec![])]
    fn test_multicenter_bond_constraints_from_iter(
        #[case] input: Vec<MulticenterBondConstraint>,
        #[case] expected: Vec<MulticenterBondConstraint>,
    ) {
        let cs = MulticenterBondConstraints::from_iter(input);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_into_iter() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![MulticenterBondConstraint::electron_count(2)]
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_from_multicenter_bond_constraint() {
        let cs: MulticenterBondConstraints = MulticenterBondConstraint::electron_count(2).into();
        assert_eq!(
            cs.as_slice(),
            &[MulticenterBondConstraint::electron_count(2)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_from_vec() {
        let cs: MulticenterBondConstraints =
            vec![MulticenterBondConstraint::electron_count(2)].into();
        assert_eq!(
            cs.as_slice(),
            &[MulticenterBondConstraint::electron_count(2)],
        );
    }
}
