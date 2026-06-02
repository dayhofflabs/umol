//! Stereo-element constraints — atom-centered and bond-centered.
//!
//! Both are empty in practice and split per site (atom vs bond), since their
//! projected structure will diverge as inhabited variants land. Each mirrors
//! `NoncovalentBondConstraint`: an uninhabited per-element enum plus a `Vec`
//! container exposing the uniform collection surface and a trivial `Lattice`.

use std::mem;
use std::slice::Iter;

use super::super::remap::IdRemapping;
use super::super::traits::Lattice;

/// Atom-centered stereo constraint. Uninhabited until value-only variants land.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum StereoAtomConstraint {}

impl StereoAtomConstraint {
    pub fn is_undetermined(&self) -> bool {
        match *self {}
    }

    pub fn remap(self, _remap: &IdRemapping) -> Option<Self> {
        match self {}
    }
}

/// Per-stereo-atom constraint container. Empty until variants land on
/// `StereoAtomConstraint`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct StereoAtomConstraints(Vec<StereoAtomConstraint>);

impl StereoAtomConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[StereoAtomConstraint] {
        &self.0
    }

    pub fn iter(&self) -> Iter<'_, StereoAtomConstraint> {
        self.0.iter()
    }

    pub fn add(&mut self, c: StereoAtomConstraint) -> Option<StereoAtomConstraint> {
        match c {}
    }

    /// Add multiple constraints at once, using the semantics of `add`.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = StereoAtomConstraint>) {
        for constraint in constraints {
            self.add(constraint);
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&StereoAtomConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = StereoAtomConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    /// No-op: the inner enum is uninhabited, so there are no values to simplify.
    pub fn simplify_each(&mut self) {}

    pub fn remap(self, _remap: &IdRemapping) -> Self {
        self
    }
}

impl Lattice for StereoAtomConstraints {
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

impl FromIterator<StereoAtomConstraint> for StereoAtomConstraints {
    fn from_iter<I: IntoIterator<Item = StereoAtomConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

/// Bond-centered stereo constraint. Uninhabited until value-only variants land.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum StereoBondConstraint {}

impl StereoBondConstraint {
    pub fn is_undetermined(&self) -> bool {
        match *self {}
    }

    pub fn remap(self, _remap: &IdRemapping) -> Option<Self> {
        match self {}
    }
}

/// Per-stereo-bond constraint container. Empty until variants land on
/// `StereoBondConstraint`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct StereoBondConstraints(Vec<StereoBondConstraint>);

impl StereoBondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[StereoBondConstraint] {
        &self.0
    }

    pub fn iter(&self) -> Iter<'_, StereoBondConstraint> {
        self.0.iter()
    }

    pub fn add(&mut self, c: StereoBondConstraint) -> Option<StereoBondConstraint> {
        match c {}
    }

    /// Add multiple constraints at once, using the semantics of `add`.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = StereoBondConstraint>) {
        for constraint in constraints {
            self.add(constraint);
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&StereoBondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = StereoBondConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    /// No-op: the inner enum is uninhabited, so there are no values to simplify.
    pub fn simplify_each(&mut self) {}

    pub fn remap(self, _remap: &IdRemapping) -> Self {
        self
    }
}

impl Lattice for StereoBondConstraints {
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

impl FromIterator<StereoBondConstraint> for StereoBondConstraints {
    fn from_iter<I: IntoIterator<Item = StereoBondConstraint>>(iter: I) -> Self {
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

    use super::*;

    #[rstest]
    fn test_stereo_atom_constraints_new() {
        let cs = StereoAtomConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[StereoAtomConstraint]);
    }

    #[rstest]
    fn test_stereo_atom_constraints_meet() {
        let cs = StereoAtomConstraints::new();
        assert_eq!(cs.meet(&cs), Some(StereoAtomConstraints::new()));
        assert!(cs.matches(&cs));
        assert!(cs.is_ground());
        assert!(cs.is_undetermined());
    }

    #[rstest]
    fn test_stereo_atom_constraints_from_iter() {
        let cs: StereoAtomConstraints = empty().collect();
        assert_eq!(cs, StereoAtomConstraints::new());
    }

    #[rstest]
    fn test_stereo_bond_constraints_new() {
        let cs = StereoBondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[StereoBondConstraint]);
    }

    #[rstest]
    fn test_stereo_bond_constraints_meet() {
        let cs = StereoBondConstraints::new();
        assert_eq!(cs.meet(&cs), Some(StereoBondConstraints::new()));
        assert!(cs.matches(&cs));
        assert!(cs.is_ground());
        assert!(cs.is_undetermined());
    }

    #[rstest]
    fn test_stereo_bond_constraints_from_iter() {
        let cs: StereoBondConstraints = empty().collect();
        assert_eq!(cs, StereoBondConstraints::new());
    }
}
