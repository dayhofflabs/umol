//! Noncovalent bond constraints.

use std::mem;
use std::slice::Iter;

use super::super::error::Contradiction;
use super::super::remap::IdRemapping;
use super::super::traits::{Canonicalize, Lattice};

/// Noncovalent-bond-scope constraint. Atom-ref and quantified-predicate forms
/// live at molecule scope via `RelationalConstraint`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NoncovalentBondConstraint {}

impl NoncovalentBondConstraint {
    pub fn key(&self) -> NoncovalentBondConstraintKey {
        match *self {}
    }

    pub fn is_unique(&self) -> bool {
        match *self {}
    }

    pub fn is_undetermined(&self) -> bool {
        match *self {}
    }

    pub fn remap(self, _remap: &IdRemapping) -> Option<Self> {
        match self {}
    }
}

/// Entry identity for `NoncovalentBondConstraint`. Uninhabited, mirroring the
/// constraint enum; exists for parity with the other constraint families.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NoncovalentBondConstraintKey {}

impl Canonicalize for NoncovalentBondConstraint {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        match self {}
    }
}

/// Per-noncovalent-bond constraint container. Empty in practice until new
/// value-only variants land on `NoncovalentBondConstraint`.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NoncovalentBondConstraints(Vec<NoncovalentBondConstraint>);

impl NoncovalentBondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[NoncovalentBondConstraint] {
        &self.0
    }

    pub fn iter(&self) -> Iter<'_, NoncovalentBondConstraint> {
        self.0.iter()
    }

    pub fn add(&mut self, c: NoncovalentBondConstraint) -> Option<NoncovalentBondConstraint> {
        match c {}
    }

    fn find_by_key(&self, key: NoncovalentBondConstraintKey) -> Result<usize, usize> {
        match key {}
    }

    pub fn contains_key(&self, key: NoncovalentBondConstraintKey) -> bool {
        self.find_by_key(key).is_ok()
    }

    pub fn get_by_key(
        &self,
        key: NoncovalentBondConstraintKey,
    ) -> Option<&NoncovalentBondConstraint> {
        self.find_by_key(key).ok().map(|i| &self.0[i])
    }

    pub fn get_by_key_mut(
        &mut self,
        key: NoncovalentBondConstraintKey,
    ) -> Option<&mut NoncovalentBondConstraint> {
        self.find_by_key(key).ok().map(|i| &mut self.0[i])
    }

    pub fn remove_by_key(
        &mut self,
        key: NoncovalentBondConstraintKey,
    ) -> Option<NoncovalentBondConstraint> {
        self.find_by_key(key).ok().map(|i| self.0.remove(i))
    }

    /// Add multiple constraints at once, using semantics of `add`.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = NoncovalentBondConstraint>) {
        for constraint in constraints {
            self.add(constraint);
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&NoncovalentBondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = NoncovalentBondConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn remap(self, _remap: &IdRemapping) -> Self {
        self
    }
}

impl Canonicalize for NoncovalentBondConstraints {
    /// Always empty (uninhabited element), so canonicalization is the identity.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(self)
    }
}

impl Lattice for NoncovalentBondConstraints {
    fn is_undetermined(&self) -> bool {
        true
    }

    fn is_ground(&self) -> bool {
        true
    }

    fn meet(&self, _other: &Self) -> Option<Self> {
        Some(Self::new())
    }

    fn join(&self, _other: &Self) -> Self {
        Self::new()
    }

    fn matches(&self, _target: &Self) -> bool {
        true
    }
}

impl FromIterator<NoncovalentBondConstraint> for NoncovalentBondConstraints {
    fn from_iter<I: IntoIterator<Item = NoncovalentBondConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use std::iter::empty;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::RemovalRemapping;

    use super::*;

    #[rstest]
    fn test_noncovalent_bond_constraints_new() {
        let cs = NoncovalentBondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[NoncovalentBondConstraint]);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_iter() {
        let cs = NoncovalentBondConstraints::new();
        assert_eq!(cs.iter().count(), 0);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_retain() {
        let mut cs = NoncovalentBondConstraints::new();
        cs.retain(|_| true);
        assert_eq!(cs, NoncovalentBondConstraints::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_clear() {
        let mut cs = NoncovalentBondConstraints::new();
        cs.clear();
        assert_eq!(cs, NoncovalentBondConstraints::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_take() {
        let mut cs = NoncovalentBondConstraints::new();
        let drained: Vec<_> = cs.take().collect();
        assert!(drained.is_empty());
        assert_eq!(cs, NoncovalentBondConstraints::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_remap() {
        let cs = NoncovalentBondConstraints::new();
        let remap = IdRemapping::new(
            RemovalRemapping::new(Vec::new(), Vec::new()),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(cs.clone().remap(&remap), cs);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_from_iter() {
        let cs: NoncovalentBondConstraints = empty().collect();
        assert_eq!(cs, NoncovalentBondConstraints::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_canonicalize() {
        let cs = NoncovalentBondConstraints::new();
        assert_eq!(cs.clone().canonicalize(), Ok(cs));
    }
}
