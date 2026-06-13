//! Stereo-element constraints — atom-centered and bond-centered.
//!
//! Both are empty in practice and split per site (atom vs bond), since their
//! projected structure will diverge as inhabited variants land. Each mirrors
//! `NoncovalentBondConstraint`: an uninhabited per-element enum plus a `Vec`
//! container exposing the uniform collection surface and a trivial `Lattice`.

use std::mem;
use std::slice::Iter;

use umol_perm::{Orientation, Permutation};

use super::super::ids::StereoLigandId;
use super::super::remap::IdRemapping;
use super::super::traits::Lattice;

/// A concrete permutation literal. A thin wrapper hosting AST-side impls that the
/// foreign `Permutation` cannot carry (matching, and EDN once the DSL lands). Not
/// a lattice — a single permutation has no top, so no `join`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PermutationAst(pub Permutation);

impl PermutationAst {
    /// A pattern matches a target iff they are the same permutation.
    pub fn matches(&self, target: &Self) -> bool {
        self.0 == target.0
    }
}

/// A concrete oriented-permutation literal: a permutation with a proper/improper
/// grade. Concrete (not a lattice); matchable.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct OrientedPermutationAst {
    pub perm: PermutationAst,
    pub orientation: Orientation,
}

impl OrientedPermutationAst {
    pub fn matches(&self, target: &Self) -> bool {
        self.perm.matches(&target.perm) && self.orientation == target.orientation
    }
}

/// Membership operator for a `±` ligand-symmetry literal: the permutation is
/// asserted present in (`In`) or absent from (`NotIn`) the ligand symmetry group.
/// Module-local (distinct from `value::MemOp`); consumed by `LigandSymmetryAst` (C3g.3).
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MemOp {
    In,
    NotIn,
}

/// The key of a per-pair topicity constraint: an unordered pair of ligand
/// positions, normalized so the lower position is `first`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LigandPairAst {
    first: StereoLigandId,
    second: StereoLigandId,
}

impl LigandPairAst {
    /// Normalizes the pair so the lower position is `first`.
    pub fn new(a: StereoLigandId, b: StereoLigandId) -> Self {
        if a.0 <= b.0 {
            Self { first: a, second: b }
        } else {
            Self { first: b, second: a }
        }
    }

    pub fn first(self) -> StereoLigandId {
        self.first
    }

    pub fn second(self) -> StereoLigandId {
        self.second
    }
}

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
    fn test_permutation_ast_matches() {
        let a = PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3]));
        let b = PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3]));
        let c = PermutationAst(Permutation::identity(4));
        assert!(a.matches(&b));
        assert!(!a.matches(&c));
    }

    #[rstest]
    fn test_oriented_permutation_ast_matches() {
        let perm = PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3]));
        let proper = OrientedPermutationAst { perm, orientation: Orientation::Proper };
        let same = OrientedPermutationAst { perm, orientation: Orientation::Proper };
        let flipped = OrientedPermutationAst { perm, orientation: Orientation::Improper };
        let other = OrientedPermutationAst {
            perm: PermutationAst(Permutation::identity(4)),
            orientation: Orientation::Proper,
        };
        assert!(proper.matches(&same));
        assert!(!proper.matches(&flipped));
        assert!(!proper.matches(&other));
    }

    #[rstest]
    #[case::ordered(StereoLigandId(1), StereoLigandId(2), StereoLigandId(1), StereoLigandId(2))]
    #[case::reversed(StereoLigandId(2), StereoLigandId(1), StereoLigandId(1), StereoLigandId(2))]
    #[case::equal(StereoLigandId(3), StereoLigandId(3), StereoLigandId(3), StereoLigandId(3))]
    fn test_ligand_pair_ast_new(
        #[case] a: StereoLigandId,
        #[case] b: StereoLigandId,
        #[case] first: StereoLigandId,
        #[case] second: StereoLigandId,
    ) {
        let pair = LigandPairAst::new(a, b);
        assert_eq!(pair.first(), first);
        assert_eq!(pair.second(), second);
    }

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
