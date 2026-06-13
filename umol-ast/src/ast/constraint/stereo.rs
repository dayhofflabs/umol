//! Stereo-element constraints — atom-centered and bond-centered.
//!
//! Both are empty in practice and split per site (atom vs bond), since their
//! projected structure will diverge as inhabited variants land. Each mirrors
//! `NoncovalentBondConstraint`: an uninhabited per-element enum plus a `Vec`
//! container exposing the uniform collection surface and a trivial `Lattice`.

use std::collections::BTreeSet;
use std::hash::{Hash, Hasher};
use std::mem;
use std::slice::Iter;

use strum::VariantArray;
use umol_perm::{Orientation, Permutation};

use super::super::ids::StereoLigandId;
use super::super::remap::IdRemapping;
use super::super::stereo::{Stereogenicity, Topicity};
use super::super::traits::{AsLit, Lattice};

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
/// Module-local (distinct from `value::MemOp`).
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
            Self {
                first: a,
                second: b,
            }
        } else {
            Self {
                first: b,
                second: a,
            }
        }
    }

    pub fn first(self) -> StereoLigandId {
        self.first
    }

    pub fn second(self) -> StereoLigandId {
        self.second
    }
}

/// Generates a subset-lattice over finite domain enum.
macro_rules! relation_ast {
    ($name:ident, $domain:ty) => {
        #[derive(Clone, Debug, Default)]
        pub enum $name {
            #[default]
            Undetermined,
            Lit($domain),
            LitSet(Vec<$domain>),
            NotSet(Vec<$domain>),
        }

        impl $name {
            fn to_set(&self) -> BTreeSet<$domain> {
                let domain = || <$domain as VariantArray>::VARIANTS.iter().copied();
                match self {
                    Self::Undetermined => domain().collect(),
                    Self::Lit(t) => BTreeSet::from([*t]),
                    Self::LitSet(values) => values.iter().copied().collect(),
                    Self::NotSet(values) => {
                        let excluded: BTreeSet<$domain> = values.iter().copied().collect();
                        domain().filter(|t| !excluded.contains(t)).collect()
                    }
                }
            }

            /// The canonical relation for a set of domain values; `None` if empty.
            fn from_set(set: BTreeSet<$domain>) -> Option<Self> {
                let domain: BTreeSet<$domain> = <$domain as VariantArray>::VARIANTS
                    .iter()
                    .copied()
                    .collect();
                if set.is_empty() {
                    None
                } else if set == domain {
                    Some(Self::Undetermined)
                } else if set.len() == 1 {
                    Some(Self::Lit(set.into_iter().next().unwrap()))
                } else {
                    let complement: Vec<$domain> = domain.difference(&set).copied().collect();
                    if set.len() <= complement.len() {
                        Some(Self::LitSet(set.into_iter().collect()))
                    } else {
                        Some(Self::NotSet(complement))
                    }
                }
            }
        }

        impl PartialEq for $name {
            fn eq(&self, other: &Self) -> bool {
                self.to_set() == other.to_set()
            }
        }

        impl Eq for $name {}

        impl Hash for $name {
            fn hash<H: Hasher>(&self, state: &mut H) {
                self.to_set().hash(state);
            }
        }

        impl Lattice for $name {
            fn is_undetermined(&self) -> bool {
                self.to_set().len() == <$domain as VariantArray>::VARIANTS.len()
            }

            fn is_ground(&self) -> bool {
                self.to_set().len() == 1
            }

            fn meet(&self, other: &Self) -> Option<Self> {
                Self::from_set(
                    self.to_set()
                        .intersection(&other.to_set())
                        .copied()
                        .collect(),
                )
            }

            fn join(&self, other: &Self) -> Self {
                Self::from_set(self.to_set().union(&other.to_set()).copied().collect())
                    .expect("union of two non-empty sets is non-empty")
            }

            fn matches(&self, target: &Self) -> bool {
                target.to_set().is_subset(&self.to_set())
            }
        }

        impl AsLit for $name {
            type Lit = $domain;

            fn as_lit(&self) -> Option<$domain> {
                let set = self.to_set();
                if set.len() == 1 {
                    set.into_iter().next()
                } else {
                    None
                }
            }
        }
    };
}

relation_ast! { TopicityRelationAst, Topicity }
relation_ast! { StereogenicityRelationAst, Stereogenicity }

/// `#p` — a signed ligand-symmetry literal: an oriented permutation asserted
/// present in (or absent from) the ligand symmetry group Π. Concrete (not a
/// lattice); matchable. Non-unique — a site may carry several.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct LigandSymmetryAst {
    pub perm: OrientedPermutationAst,
    pub mem: MemOp,
}

impl LigandSymmetryAst {
    pub fn matches(&self, target: &Self) -> bool {
        self.perm.matches(&target.perm) && self.mem == target.mem
    }
}

/// `#f` — a fluxionality move: a proper ligand-position permutation realized by
/// dynamics. Concrete (not a lattice); matchable. Non-unique.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FluxionalityAst {
    pub perm: PermutationAst,
}

impl FluxionalityAst {
    pub fn matches(&self, target: &Self) -> bool {
        self.perm.matches(&target.perm)
    }
}

/// `#o` — a per-pair topicity constraint: the relation holding between one ligand
/// pair. Keyed by `pair`; the lattice operates on the relation, per pair.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct TopicityAst {
    pub pair: LigandPairAst,
    pub rel: TopicityRelationAst,
}

impl Lattice for TopicityAst {
    fn is_undetermined(&self) -> bool {
        self.rel.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.rel.is_ground()
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        debug_assert_eq!(self.pair, other.pair, "topicity meet is per-pair");
        self.rel.meet(&other.rel).map(|rel| Self {
            pair: self.pair,
            rel,
        })
    }

    fn join(&self, other: &Self) -> Self {
        debug_assert_eq!(self.pair, other.pair, "topicity join is per-pair");
        Self {
            pair: self.pair,
            rel: self.rel.join(&other.rel),
        }
    }

    fn matches(&self, target: &Self) -> bool {
        self.pair == target.pair && self.rel.matches(&target.rel)
    }
}

/// `#g` — the stereogenicity classification constraint (unique per site).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct StereogenicityAst(pub StereogenicityRelationAst);

impl Lattice for StereogenicityAst {
    fn is_undetermined(&self) -> bool {
        self.0.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.0.is_ground()
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        self.0.meet(&other.0).map(Self)
    }

    fn join(&self, other: &Self) -> Self {
        Self(self.0.join(&other.0))
    }

    fn matches(&self, target: &Self) -> bool {
        self.0.matches(&target.0)
    }
}

impl AsLit for StereogenicityAst {
    type Lit = Stereogenicity;

    fn as_lit(&self) -> Option<Stereogenicity> {
        self.0.as_lit()
    }
}

/// Atom-centered stereo constraint.
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

/// Per-stereo-atom constraint container.
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

/// Bond-centered stereo constraint.
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

/// Per-stereo-bond constraint container.
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
        let proper = OrientedPermutationAst {
            perm,
            orientation: Orientation::Proper,
        };
        let same = OrientedPermutationAst {
            perm,
            orientation: Orientation::Proper,
        };
        let flipped = OrientedPermutationAst {
            perm,
            orientation: Orientation::Improper,
        };
        let other = OrientedPermutationAst {
            perm: PermutationAst(Permutation::identity(4)),
            orientation: Orientation::Proper,
        };
        assert!(proper.matches(&same));
        assert!(!proper.matches(&flipped));
        assert!(!proper.matches(&other));
    }

    #[rustfmt::skip]
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
    fn test_topicity_relation_ast_as_lit() {
        assert_eq!(
            TopicityRelationAst::Lit(Topicity::Homotopic).as_lit(),
            Some(Topicity::Homotopic),
        );
        assert_eq!(TopicityRelationAst::Undetermined.as_lit(), None);
        // NotSet([H]) denotes {E,D} — not a singleton.
        assert_eq!(
            TopicityRelationAst::NotSet(vec![Topicity::Homotopic]).as_lit(),
            None,
        );
    }

    #[rstest]
    fn test_topicity_relation_ast_eq_representation_independent() {
        // {H,E}: LitSet([H,E]) and NotSet([D]) denote the same set.
        assert_eq!(
            TopicityRelationAst::LitSet(vec![Topicity::Homotopic, Topicity::Enantiotopic]),
            TopicityRelationAst::NotSet(vec![Topicity::Diastereotopic]),
        );
        // The full set ≡ Undetermined.
        assert_eq!(
            TopicityRelationAst::LitSet(vec![
                Topicity::Homotopic,
                Topicity::Enantiotopic,
                Topicity::Diastereotopic,
            ]),
            TopicityRelationAst::Undetermined,
        );
    }

    #[rstest]
    fn test_topicity_relation_ast_lattice() {
        let h = TopicityRelationAst::Lit(Topicity::Homotopic);
        let e = TopicityRelationAst::Lit(Topicity::Enantiotopic);
        assert_eq!(h.meet(&h), Some(h.clone()));
        assert_eq!(h.meet(&e), None); // disjoint singletons ⇒ contradiction
        assert_eq!(TopicityRelationAst::Undetermined.meet(&h), Some(h.clone()));
        // {H} ∪ {E} = {H,E}, canonically NotSet([D]) on a 3-domain.
        assert_eq!(
            h.join(&e),
            TopicityRelationAst::NotSet(vec![Topicity::Diastereotopic])
        );
        assert!(TopicityRelationAst::Undetermined.matches(&h));
        assert!(h.matches(&h));
        assert!(!h.matches(&e));
        assert!(h.is_ground());
        assert!(!TopicityRelationAst::Undetermined.is_ground());
        assert!(TopicityRelationAst::Undetermined.is_undetermined());
        assert!(!h.is_undetermined());
    }

    #[rstest]
    fn test_ligand_symmetry_ast_matches() {
        let perm = OrientedPermutationAst {
            perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])),
            orientation: Orientation::Proper,
        };
        let present = LigandSymmetryAst { perm, mem: MemOp::In };
        let same = LigandSymmetryAst { perm, mem: MemOp::In };
        let absent = LigandSymmetryAst { perm, mem: MemOp::NotIn };
        let other = LigandSymmetryAst {
            perm: OrientedPermutationAst {
                perm: PermutationAst(Permutation::identity(4)),
                orientation: Orientation::Proper,
            },
            mem: MemOp::In,
        };
        assert!(present.matches(&same));
        assert!(!present.matches(&absent)); // different membership op
        assert!(!present.matches(&other)); // different permutation
    }

    #[rstest]
    fn test_fluxionality_ast_matches() {
        let a = FluxionalityAst {
            perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])),
        };
        let same = FluxionalityAst {
            perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])),
        };
        let other = FluxionalityAst {
            perm: PermutationAst(Permutation::identity(4)),
        };
        assert!(a.matches(&same));
        assert!(!a.matches(&other));
    }

    #[rstest]
    fn test_topicity_ast_lattice() {
        let pair = LigandPairAst::new(StereoLigandId(0), StereoLigandId(1));
        let h = TopicityAst { pair, rel: TopicityRelationAst::Lit(Topicity::Homotopic) };
        let e = TopicityAst { pair, rel: TopicityRelationAst::Lit(Topicity::Enantiotopic) };
        let open = TopicityAst { pair, rel: TopicityRelationAst::Undetermined };
        assert_eq!(h.meet(&e), None); // disjoint relations on the same pair
        assert_eq!(open.meet(&h), Some(h.clone()));
        assert!(open.matches(&h));
        assert!(h.matches(&h));
        assert!(!h.matches(&e));
        assert!(h.is_ground());
        assert!(open.is_undetermined());
        // A constraint on a different pair never matches.
        let elsewhere = TopicityAst {
            pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(2)),
            rel: TopicityRelationAst::Undetermined,
        };
        assert!(!elsewhere.matches(&h));
    }

    #[rstest]
    fn test_stereogenicity_ast_lattice_and_as_lit() {
        let g = StereogenicityAst(StereogenicityRelationAst::Lit(Stereogenicity::Stereogenic));
        let open = StereogenicityAst(StereogenicityRelationAst::Undetermined);
        assert_eq!(g.as_lit(), Some(Stereogenicity::Stereogenic));
        assert_eq!(open.as_lit(), None);
        assert_eq!(open.meet(&g), Some(g.clone()));
        assert!(open.matches(&g));
        assert!(g.is_ground());
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
