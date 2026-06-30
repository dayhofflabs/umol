//! Per-aromatic-system constraints.

use std::mem::{self, replace};
use std::slice::Iter;
use std::vec::IntoIter;

use strum::EnumDiscriminants;

use super::super::error::Contradiction;
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Aromatic-system-scope constraint. Held inline on `AromaticSystemAst` via
/// `AromaticSystemConstraints`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, EnumDiscriminants)]
#[strum_discriminants(name(AromaticSystemConstraintKind), derive(Hash))]
pub enum AromaticSystemConstraint {
    /// Asserted total π-electron count for the system. Cross-checked by the
    /// `ConsistencyValidator` against `sum(AromaticSystemAst::electrons)`.
    ElectronCount(ValueAst),
}

impl AromaticSystemConstraint {
    pub fn electron_count(v: impl Into<ValueAst>) -> Self {
        Self::ElectronCount(v.into())
    }

    pub fn kind(&self) -> AromaticSystemConstraintKind {
        self.into()
    }

    /// Entry identity for order/dedup. Every kind is single-valued, so no sub-key.
    pub fn key(&self) -> AromaticSystemConstraintKey {
        match self {
            Self::ElectronCount(_) => AromaticSystemConstraintKey::ElectronCount,
        }
    }

    /// Every `AromaticSystemConstraint` variant is single-valued per system.
    pub fn is_unique(&self) -> bool {
        true
    }

    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::ElectronCount(v) => v.is_undetermined(),
        }
    }

    pub fn compact(self, _remap: &IdCompaction) -> Option<Self> {
        // Value-only: no indices to remap.
        Some(self)
    }

    /// Value-only: no indices to remap.
    pub(crate) fn remap(self, _map: &IdRemapping) -> Self {
        self
    }
}

/// Entry identity: discriminant only (every kind is single-valued, no sub-key).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticSystemConstraintKey {
    ElectronCount,
}

impl AromaticSystemConstraintKey {
    pub fn kind(self) -> AromaticSystemConstraintKind {
        match self {
            Self::ElectronCount => AromaticSystemConstraintKind::ElectronCount,
        }
    }
}

impl Canonicalize for AromaticSystemConstraint {
    /// Canonicalize the inner value; the kind is preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::ElectronCount(v) => Self::ElectronCount(v.canonicalize()?),
        })
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AromaticSystemConstraints(Vec<AromaticSystemConstraint>);

impl AromaticSystemConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[AromaticSystemConstraint] {
        &self.0
    }

    pub fn contains(&self, kind: AromaticSystemConstraintKind) -> bool {
        self.0.iter().any(|c| c.kind() == kind)
    }

    pub fn get(&self, kind: AromaticSystemConstraintKind) -> Option<&AromaticSystemConstraint> {
        self.0.iter().find(|c| c.kind() == kind)
    }

    pub fn get_mut(
        &mut self,
        kind: AromaticSystemConstraintKind,
    ) -> Option<&mut AromaticSystemConstraint> {
        self.0.iter_mut().find(|c| c.kind() == kind)
    }

    pub fn electron_count(&self) -> ValueAst {
        match self.get(AromaticSystemConstraintKind::ElectronCount) {
            Some(AromaticSystemConstraint::ElectronCount(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn iter(&self) -> Iter<'_, AromaticSystemConstraint> {
        self.0.iter()
    }

    /// Insert at the `key()`-sorted position: unique kinds replace the same-key
    /// entry (returning it); non-unique kinds append after the same-key run.
    pub fn add(&mut self, c: AromaticSystemConstraint) -> Option<AromaticSystemConstraint> {
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

    fn find_by_key(&self, key: AromaticSystemConstraintKey) -> Result<usize, usize> {
        self.0.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains_key(&self, key: AromaticSystemConstraintKey) -> bool {
        self.find_by_key(key).is_ok()
    }

    pub fn get_by_key(
        &self,
        key: AromaticSystemConstraintKey,
    ) -> Option<&AromaticSystemConstraint> {
        self.find_by_key(key).ok().map(|i| &self.0[i])
    }

    pub fn get_by_key_mut(
        &mut self,
        key: AromaticSystemConstraintKey,
    ) -> Option<&mut AromaticSystemConstraint> {
        self.find_by_key(key).ok().map(|i| &mut self.0[i])
    }

    pub fn remove_by_key(
        &mut self,
        key: AromaticSystemConstraintKey,
    ) -> Option<AromaticSystemConstraint> {
        self.find_by_key(key).ok().map(|i| self.0.remove(i))
    }

    /// Add multiple constraints at once, using semantics of `add`.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = AromaticSystemConstraint>) {
        for constraint in constraints {
            self.add(constraint);
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&AromaticSystemConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    pub fn take(&mut self) -> impl Iterator<Item = AromaticSystemConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn remove(
        &mut self,
        kind: AromaticSystemConstraintKind,
    ) -> Option<AromaticSystemConstraint> {
        let pos = self.0.iter().position(|c| c.kind() == kind)?;
        Some(self.0.remove(pos))
    }

    /// Iterate over every entry of `kind`. Currently every variant is
    /// single-valued so this yields at most one entry.
    pub fn get_all(
        &self,
        kind: AromaticSystemConstraintKind,
    ) -> impl Iterator<Item = &AromaticSystemConstraint> {
        self.0.iter().filter(move |c| c.kind() == kind)
    }

    /// Remove every entry of `kind`, returning them in insertion order.
    pub fn remove_all(
        &mut self,
        kind: AromaticSystemConstraintKind,
    ) -> Vec<AromaticSystemConstraint> {
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

    pub fn compact(self, remap: &IdCompaction) -> Self {
        Self(self.0.into_iter().filter_map(|c| c.compact(remap)).collect())
    }
}

impl Canonicalize for AromaticSystemConstraints {
    /// Sort by `key()`, canonicalize each value, drop vacuous entries. No merge
    /// clause — every kind is single-valued, so `add` admits no same-key duplicates.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut entries = self.0;
        entries.sort_by_key(|c| c.key());
        let mut out: Vec<AromaticSystemConstraint> = Vec::with_capacity(entries.len());
        for c in entries {
            out.push(c.canonicalize()?);
        }
        out.retain(|c| !c.is_undetermined());
        Ok(Self(out))
    }
}

impl Lattice for AromaticSystemConstraints {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| match c {
            AromaticSystemConstraint::ElectronCount(v) => v.is_undetermined(),
        })
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| match c {
            AromaticSystemConstraint::ElectronCount(v) => v.is_ground(),
        })
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let mut result = Self::new();
        let merged = self.electron_count().meet(&other.electron_count())?;
        if !merged.is_undetermined() {
            result.add(AromaticSystemConstraint::ElectronCount(merged));
        }
        Some(result)
    }

    fn join(&self, other: &Self) -> Self {
        let mut result = Self::new();
        let a_has = self.contains(AromaticSystemConstraintKind::ElectronCount);
        let b_has = other.contains(AromaticSystemConstraintKind::ElectronCount);
        if a_has && b_has {
            let joined = self.electron_count().join(&other.electron_count());
            if !joined.is_undetermined() {
                result.add(AromaticSystemConstraint::ElectronCount(joined));
            }
        }
        result
    }

    fn matches(&self, target: &Self) -> bool {
        self.electron_count().matches(&target.electron_count())
    }
}

impl FromIterator<AromaticSystemConstraint> for AromaticSystemConstraints {
    fn from_iter<I: IntoIterator<Item = AromaticSystemConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

impl IntoIterator for AromaticSystemConstraints {
    type Item = AromaticSystemConstraint;
    type IntoIter = IntoIter<AromaticSystemConstraint>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<AromaticSystemConstraint> for AromaticSystemConstraints {
    fn from(c: AromaticSystemConstraint) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<AromaticSystemConstraint>> for AromaticSystemConstraints {
    fn from(cs: Vec<AromaticSystemConstraint>) -> Self {
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
        AromaticSystemConstraint::electron_count(6),
        AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))
    )]
    fn test_aromatic_system_constraint_constructors(
        #[case] actual: AromaticSystemConstraint,
        #[case] expected: AromaticSystemConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::electron_count(
        AromaticSystemConstraint::electron_count(6),
        AromaticSystemConstraintKind::ElectronCount
    )]
    fn test_aromatic_system_constraint_kind(
        #[case] c: AromaticSystemConstraint,
        #[case] expected: AromaticSystemConstraintKind,
    ) {
        assert_eq!(c.kind(), expected);
    }

    #[rstest]
    #[case::electron_count(
        AromaticSystemConstraint::electron_count(6),
        AromaticSystemConstraintKey::ElectronCount
    )]
    fn test_aromatic_system_constraint_key(
        #[case] c: AromaticSystemConstraint,
        #[case] expected: AromaticSystemConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rstest]
    #[case::electron_count(
        AromaticSystemConstraintKey::ElectronCount,
        AromaticSystemConstraintKind::ElectronCount
    )]
    fn test_aromatic_system_constraint_key_kind(
        #[case] key: AromaticSystemConstraintKey,
        #[case] expected: AromaticSystemConstraintKind,
    ) {
        assert_eq!(key.kind(), expected);
    }

    #[rstest]
    #[case::electron_count(AromaticSystemConstraint::electron_count(6))]
    fn test_aromatic_system_constraint_is_unique(#[case] c: AromaticSystemConstraint) {
        assert!(c.is_unique());
    }

    #[rstest]
    #[case::lit(AromaticSystemConstraint::electron_count(6), false)]
    #[case::undetermined(AromaticSystemConstraint::ElectronCount(ValueAst::Undetermined), true)]
    fn test_aromatic_system_constraint_is_undetermined(
        #[case] c: AromaticSystemConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count_litset_singleton(AromaticSystemConstraint::ElectronCount(ValueAst::lit_set([6])), Ok(AromaticSystemConstraint::electron_count(6)))]
    #[case::empty_litset_contradiction(AromaticSystemConstraint::ElectronCount(ValueAst::lit_set(Vec::<i64>::new())), Err(Contradiction))]
    fn test_aromatic_system_constraint_canonicalize(
        #[case] constraint: AromaticSystemConstraint,
        #[case] expected: Result<AromaticSystemConstraint, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
    }

    #[rstest]
    fn test_aromatic_system_constraints_new() {
        let cs = AromaticSystemConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[AromaticSystemConstraint]);
    }

    #[rstest]
    #[case::present(AromaticSystemConstraintKind::ElectronCount, true)]
    fn test_aromatic_system_constraints_contains(
        #[case] kind: AromaticSystemConstraintKind,
        #[case] expected: bool,
    ) {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        assert_eq!(cs.contains(kind), expected);
    }

    #[rstest]
    fn test_aromatic_system_constraints_contains_absent() {
        let cs = AromaticSystemConstraints::new();
        assert!(!cs.contains(AromaticSystemConstraintKind::ElectronCount));
    }

    #[rstest]
    fn test_aromatic_system_constraints_get() {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        assert_eq!(
            cs.get(AromaticSystemConstraintKind::ElectronCount),
            Some(&AromaticSystemConstraint::electron_count(6)),
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_get_absent() {
        let cs = AromaticSystemConstraints::new();
        assert_eq!(cs.get(AromaticSystemConstraintKind::ElectronCount), None,);
    }

    #[rstest]
    fn test_aromatic_system_constraints_get_mut() {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        let entry = cs
            .get_mut(AromaticSystemConstraintKind::ElectronCount)
            .unwrap();
        *entry = AromaticSystemConstraint::electron_count(10);
        assert_eq!(
            cs.as_slice(),
            &[AromaticSystemConstraint::electron_count(10)],
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_get_mut_absent() {
        let mut cs = AromaticSystemConstraints::new();
        assert!(cs
            .get_mut(AromaticSystemConstraintKind::ElectronCount)
            .is_none());
    }

    #[rstest]
    fn test_aromatic_system_constraints_contains_key() {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        assert!(cs.contains_key(AromaticSystemConstraintKey::ElectronCount));
        assert!(!AromaticSystemConstraints::new()
            .contains_key(AromaticSystemConstraintKey::ElectronCount));
    }

    #[rstest]
    fn test_aromatic_system_constraints_get_by_key() {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        assert_eq!(
            cs.get_by_key(AromaticSystemConstraintKey::ElectronCount),
            Some(&AromaticSystemConstraint::electron_count(6)),
        );
        assert_eq!(
            AromaticSystemConstraints::new().get_by_key(AromaticSystemConstraintKey::ElectronCount),
            None,
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_get_by_key_mut() {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        *cs.get_by_key_mut(AromaticSystemConstraintKey::ElectronCount)
            .unwrap() = AromaticSystemConstraint::electron_count(10);
        assert_eq!(
            cs.get_by_key(AromaticSystemConstraintKey::ElectronCount),
            Some(&AromaticSystemConstraint::electron_count(10)),
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_remove_by_key() {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        assert_eq!(
            cs.remove_by_key(AromaticSystemConstraintKey::ElectronCount),
            Some(AromaticSystemConstraint::electron_count(6)),
        );
        assert_eq!(cs.as_slice(), &[] as &[AromaticSystemConstraint]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::canonicalizes_value(
        AromaticSystemConstraints::from(AromaticSystemConstraint::ElectronCount(ValueAst::lit_set([6]))),
        Ok(AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6))))]
    #[case::drop_vacuous(
        AromaticSystemConstraints::from(AromaticSystemConstraint::ElectronCount(ValueAst::Undetermined)),
        Ok(AromaticSystemConstraints::new()))]
    #[case::contradiction(
        AromaticSystemConstraints::from(AromaticSystemConstraint::ElectronCount(ValueAst::lit_set(Vec::<i64>::new()))),
        Err(Contradiction))]
    fn test_aromatic_system_constraints_canonicalize(
        #[case] constraints: AromaticSystemConstraints,
        #[case] expected: Result<AromaticSystemConstraints, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rstest]
    fn test_aromatic_system_constraints_iter() {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(collected, vec![AromaticSystemConstraint::electron_count(6)]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(
        vec![AromaticSystemConstraint::electron_count(6)],
        vec![None],
        vec![AromaticSystemConstraint::electron_count(6)],
    )]
    #[case::replace_same_kind(
        vec![
            AromaticSystemConstraint::electron_count(6),
            AromaticSystemConstraint::electron_count(10),
        ],
        vec![None, Some(AromaticSystemConstraint::electron_count(6))],
        vec![AromaticSystemConstraint::electron_count(10)],
    )]
    fn test_aromatic_system_constraints_add(
        #[case] sequence: Vec<AromaticSystemConstraint>,
        #[case] expected_returns: Vec<Option<AromaticSystemConstraint>>,
        #[case] expected_state: Vec<AromaticSystemConstraint>,
    ) {
        let mut cs = AromaticSystemConstraints::new();
        let returns: Vec<_> = sequence.into_iter().map(|c| cs.add(c)).collect();
        assert_eq!(returns, expected_returns);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(
        |c: &AromaticSystemConstraint| matches!(c, AromaticSystemConstraint::ElectronCount(_)),
        vec![AromaticSystemConstraint::electron_count(6)],
    )]
    #[case::all_dropped(|_: &AromaticSystemConstraint| false, vec![])]
    fn test_aromatic_system_constraints_retain(
        #[case] predicate: impl FnMut(&AromaticSystemConstraint) -> bool,
        #[case] expected: Vec<AromaticSystemConstraint>,
    ) {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        cs.retain(predicate);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_aromatic_system_constraints_clear() {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        cs.clear();
        assert_eq!(cs, AromaticSystemConstraints::new());
    }

    #[rstest]
    fn test_aromatic_system_constraints_take() {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(drained, vec![AromaticSystemConstraint::electron_count(6)]);
        assert_eq!(cs, AromaticSystemConstraints::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::present(
        AromaticSystemConstraintKind::ElectronCount,
        Some(AromaticSystemConstraint::electron_count(6)),
        vec![],
    )]
    fn test_aromatic_system_constraints_remove(
        #[case] kind: AromaticSystemConstraintKind,
        #[case] expected_returned: Option<AromaticSystemConstraint>,
        #[case] expected_state: Vec<AromaticSystemConstraint>,
    ) {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        assert_eq!(cs.remove(kind), expected_returned);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rstest]
    fn test_aromatic_system_constraints_remove_absent() {
        let mut cs = AromaticSystemConstraints::new();
        assert_eq!(cs.remove(AromaticSystemConstraintKind::ElectronCount), None,);
    }

    #[rstest]
    fn test_aromatic_system_constraints_remap() {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        let remap = IdCompaction::new(
            Compaction::new(vec![0, 1], vec![0]),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(cs.clone().compact(&remap), cs);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(
        vec![AromaticSystemConstraint::electron_count(6)],
        vec![AromaticSystemConstraint::electron_count(6)],
    )]
    #[case::same_kind_last_wins(
        vec![
            AromaticSystemConstraint::electron_count(2),
            AromaticSystemConstraint::electron_count(6),
        ],
        vec![AromaticSystemConstraint::electron_count(6)],
    )]
    #[case::empty(vec![], vec![])]
    fn test_aromatic_system_constraints_from_iter(
        #[case] input: Vec<AromaticSystemConstraint>,
        #[case] expected: Vec<AromaticSystemConstraint>,
    ) {
        let cs = AromaticSystemConstraints::from_iter(input);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_aromatic_system_constraints_into_iter() {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(collected, vec![AromaticSystemConstraint::electron_count(6)]);
    }

    #[rstest]
    fn test_aromatic_system_constraints_from_aromatic_system_constraint() {
        let cs: AromaticSystemConstraints = AromaticSystemConstraint::electron_count(6).into();
        assert_eq!(
            cs.as_slice(),
            &[AromaticSystemConstraint::electron_count(6)]
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_from_vec() {
        let cs: AromaticSystemConstraints =
            vec![AromaticSystemConstraint::electron_count(6)].into();
        assert_eq!(
            cs.as_slice(),
            &[AromaticSystemConstraint::electron_count(6)]
        );
    }
}
