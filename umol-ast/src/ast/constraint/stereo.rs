//! Stereo-element constraints.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::hash::Hash;
use std::mem::{self, replace};
use std::slice::Iter;

use smallvec::SmallVec;
use strum::VariantArray;
use umol_perm::{Orientation, Permutation};

use super::super::error::Contradiction;
use super::super::ids::StereoLigandId;
use super::super::operators::MemOp;
use super::super::remap::IdRemapping;
use super::super::stereo::{Stereogenicity, Topicity};
use super::super::traits::{AsLit, Canonicalize, Lattice};

/// A concrete permutation literal. A thin wrapper hosting AST-side impls that the
/// foreign `Permutation` cannot carry. Not a lattice — a single permutation has no top, so no `join`.
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
        #[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
        pub enum $name {
            #[default]
            Undetermined,
            Lit($domain),
            LitSet(BTreeSet<$domain>),
            NotSet(BTreeSet<$domain>),
        }

        impl $name {
            /// The semantic set of admissible domain values.
            pub(crate) fn to_set(&self) -> BTreeSet<$domain> {
                let domain = || <$domain as VariantArray>::VARIANTS.iter().copied();
                match self {
                    Self::Undetermined => domain().collect(),
                    Self::Lit(t) => BTreeSet::from([*t]),
                    Self::LitSet(values) => values.clone(),
                    Self::NotSet(values) => domain().filter(|t| !values.contains(t)).collect(),
                }
            }

            /// The canonical relation for a set of admissible values; `None` if empty.
            /// Singleton → `Lit`, full → `Undetermined`, else the smaller of
            /// positive/complement (tiebreak positive).
            pub(crate) fn from_set(set: BTreeSet<$domain>) -> Option<Self> {
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
                    let complement: BTreeSet<$domain> = domain.difference(&set).copied().collect();
                    if set.len() <= complement.len() {
                        Some(Self::LitSet(set))
                    } else {
                        Some(Self::NotSet(complement))
                    }
                }
            }
        }

        impl Canonicalize for $name {
            /// Finite-domain canonical form: sort/dedup (`BTreeSet`), singleton → `Lit`,
            /// full → `Undetermined`, empty → `Err`, polarity (store the smaller side).
            fn canonicalize(self) -> Result<Self, Contradiction> {
                Self::from_set(self.to_set()).ok_or(Contradiction)
            }

            fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
                match self {
                    Self::Undetermined | Self::Lit(_) => Ok(Cow::Borrowed(self)),
                    _ => Ok(Cow::Owned(self.clone().canonicalize()?)),
                }
            }
        }

        impl Lattice for $name {
            fn is_undetermined(&self) -> bool {
                matches!(self, Self::Undetermined)
            }

            fn is_ground(&self) -> bool {
                matches!(self, Self::Lit(_))
            }

            /// Intersection of the admissible sets, canonicalized (∅ → `None`).
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
                match self {
                    Self::Lit(t) => Some(*t),
                    _ => None,
                }
            }
        }
    };
}

relation_ast! { TopicityRelationAst, Topicity }
relation_ast! { StereogenicityAst, Stereogenicity }

/// Ligand permutations with membership assertion. Non-unique.
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

/// Fluxionality move: proper ligand permutation realized by dynamics. Non-unique.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FluxionalityAst {
    pub perm: PermutationAst,
}

impl FluxionalityAst {
    pub fn matches(&self, target: &Self) -> bool {
        self.perm.matches(&target.perm)
    }
}

/// Per-pair topicity constraint: relation between pair of ligands.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct TopicityAst {
    pub pair: LigandPairAst,
    pub rel: TopicityRelationAst,
}

impl TopicityAst {
    /// Matches-only (not a `Lattice`): different pairs are incomparable, so there
    /// is no global top. `matches` = same `pair` and the per-pair `rel` matches;
    /// the per-pair lattice lives in [`TopicityRelationAst`].
    pub fn matches(&self, target: &Self) -> bool {
        self.pair == target.pair && self.rel.matches(&target.rel)
    }
}

/// Per-element stereo constraint for stereo atoms and bonds.
macro_rules! stereo_constraint {
    ($constraint:ident, $kind:ident, $constraints:ident) => {
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
        pub enum $kind {
            LigandSymmetry,
            Fluxionality,
            Topicity,
            Stereogenicity,
        }

        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub enum $constraint {
            LigandSymmetry(LigandSymmetryAst),
            Fluxionality(FluxionalityAst),
            Topicity(TopicityAst),
            Stereogenicity(StereogenicityAst),
        }

        impl $constraint {
            pub fn kind(&self) -> $kind {
                match self {
                    Self::LigandSymmetry(_) => $kind::LigandSymmetry,
                    Self::Fluxionality(_) => $kind::Fluxionality,
                    Self::Topicity(_) => $kind::Topicity,
                    Self::Stereogenicity(_) => $kind::Stereogenicity,
                }
            }

            /// Whether at most one constraint of this kind may be stored.
            pub fn is_unique(&self) -> bool {
                matches!(self, Self::Stereogenicity(_))
            }

            pub fn is_undetermined(&self) -> bool {
                match self {
                    Self::LigandSymmetry(_) | Self::Fluxionality(_) => false,
                    Self::Topicity(t) => t.rel.is_undetermined(),
                    Self::Stereogenicity(g) => g.is_undetermined(),
                }
            }

            /// Frame-relative ligand positions carry no atom ids, so remap is a no-op.
            pub fn remap(self, _remap: &IdRemapping) -> Option<Self> {
                Some(self)
            }
        }

        /// Stereo constraint container for stereo atoms and bonds.
        #[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
        pub struct $constraints {
            entries: SmallVec<[$constraint; 2]>,
        }

        impl $constraints {
            pub fn new() -> Self {
                Self::default()
            }

            pub fn is_empty(&self) -> bool {
                self.entries.is_empty()
            }

            pub fn len(&self) -> usize {
                self.entries.len()
            }

            pub fn contains(&self, kind: $kind) -> bool {
                self.find(kind).is_ok()
            }

            pub fn get(&self, kind: $kind) -> Option<&$constraint> {
                self.find(kind).ok().map(|i| &self.entries[i])
            }

            /// Ligand-symmetry literals (non-unique).
            pub fn ligand_symmetry(&self) -> impl Iterator<Item = &LigandSymmetryAst> {
                self.entries.iter().filter_map(|c| match c {
                    $constraint::LigandSymmetry(p) => Some(p),
                    _ => None,
                })
            }

            /// Fluxionality moves (non-unique).
            pub fn fluxionality(&self) -> impl Iterator<Item = &FluxionalityAst> {
                self.entries.iter().filter_map(|c| match c {
                    $constraint::Fluxionality(f) => Some(f),
                    _ => None,
                })
            }

            /// pair-specific topicity constraint.
            pub fn topicities(&self) -> impl Iterator<Item = &TopicityAst> {
                self.entries.iter().filter_map(|c| match c {
                    $constraint::Topicity(t) => Some(t),
                    _ => None,
                })
            }

            /// Topicity relation per ligand pair.
            pub fn topicity(&self, pair: LigandPairAst) -> TopicityRelationAst {
                self.topicities()
                    .find(|t| t.pair == pair)
                    .map(|t| t.rel.clone())
                    .unwrap_or_default()
            }

            /// Stored stereogenicity relation (unique); `Undetermined` if absent.
            pub fn stereogenicity(&self) -> StereogenicityAst {
                match self.get($kind::Stereogenicity) {
                    Some($constraint::Stereogenicity(g)) => g.clone(),
                    _ => StereogenicityAst::Undetermined,
                }
            }

            /// Add a constraint to the container. Unique constraint replace,
            /// non-unique constraint append.
            pub fn add(&mut self, c: $constraint) -> Option<$constraint> {
                let kind = c.kind();
                let start = self
                    .entries
                    .partition_point(|e| (e.kind() as u8) < (kind as u8));
                let end = start
                    + self.entries[start..]
                        .iter()
                        .take_while(|e| e.kind() == kind)
                        .count();
                if c.is_unique() {
                    if start < end {
                        return Some(replace(&mut self.entries[start], c));
                    }
                    self.entries.insert(start, c);
                    None
                } else if let $constraint::Topicity(ref t) = c {
                    let pair = t.pair;
                    let existing = self.entries[start..end]
                        .iter()
                        .position(|e| matches!(e, $constraint::Topicity(x) if x.pair == pair));
                    match existing {
                        Some(j) => Some(replace(&mut self.entries[start + j], c)),
                        None => {
                            self.entries.insert(end, c);
                            None
                        }
                    }
                } else {
                    self.entries.insert(end, c);
                    None
                }
            }

            /// Add multiple constraints at once, using the semantics of `add`.
            pub fn extend(&mut self, constraints: impl IntoIterator<Item = $constraint>) {
                for constraint in constraints {
                    self.add(constraint);
                }
            }

            pub fn retain(&mut self, mut f: impl FnMut(&$constraint) -> bool) {
                self.entries.retain(|c| f(c));
            }

            pub fn clear(&mut self) {
                self.entries.clear();
            }

            /// Move the entries out of the store, leaving it empty.
            pub fn take(&mut self) -> impl Iterator<Item = $constraint> {
                mem::take(&mut self.entries).into_iter()
            }

            /// No-op: stereo constraints hold no `ValueAst`; relations are already canonical.
            pub fn simplify_each(&mut self) {}

            pub fn remove(&mut self, kind: $kind) -> Option<$constraint> {
                self.find(kind).ok().map(|i| self.entries.remove(i))
            }

            /// Remove every entry of `kind`, returning them in store order.
            pub fn remove_all(&mut self, kind: $kind) -> Vec<$constraint> {
                let start = self
                    .entries
                    .partition_point(|e| (e.kind() as u8) < (kind as u8));
                let end = start
                    + self.entries[start..]
                        .iter()
                        .take_while(|e| e.kind() == kind)
                        .count();
                self.entries.drain(start..end).collect()
            }

            /// Iterate over every entry of `kind`.
            pub fn get_all(&self, kind: $kind) -> impl Iterator<Item = &$constraint> {
                let start = self
                    .entries
                    .partition_point(|c| (c.kind() as u8) < (kind as u8));
                self.entries[start..]
                    .iter()
                    .take_while(move |c| c.kind() == kind)
            }

            pub fn iter(&self) -> Iter<'_, $constraint> {
                self.entries.iter()
            }

            pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut $constraint> {
                self.entries.iter_mut()
            }

            /// True if any entry exactly equals `constraint`.
            pub fn contains_entry(&self, constraint: &$constraint) -> bool {
                self.entries.iter().any(|c| c == constraint)
            }

            /// No-op: frame-relative ligand positions carry no entity index.
            pub fn remap(self, _remap: &IdRemapping) -> Self {
                self
            }

            fn find(&self, kind: $kind) -> Result<usize, usize> {
                self.entries
                    .binary_search_by_key(&(kind as u8), |c| c.kind() as u8)
            }
        }

        impl Lattice for $constraints {
            fn is_undetermined(&self) -> bool {
                self.entries.iter().all(|c| c.is_undetermined())
            }

            fn is_ground(&self) -> bool {
                self.entries.iter().all(|c| match c {
                    $constraint::LigandSymmetry(_) | $constraint::Fluxionality(_) => true,
                    $constraint::Topicity(t) => t.rel.is_ground(),
                    $constraint::Stereogenicity(g) => g.is_ground(),
                })
            }

            fn meet(&self, other: &Self) -> Option<Self> {
                let mut result = Self::new();
                for p in self.ligand_symmetry().chain(other.ligand_symmetry()) {
                    let entry = $constraint::LigandSymmetry(*p);
                    if !result.contains_entry(&entry) {
                        result.add(entry);
                    }
                }
                for f in self.fluxionality().chain(other.fluxionality()) {
                    let entry = $constraint::Fluxionality(*f);
                    if !result.contains_entry(&entry) {
                        result.add(entry);
                    }
                }
                let pairs: BTreeSet<LigandPairAst> = self
                    .topicities()
                    .chain(other.topicities())
                    .map(|t| t.pair)
                    .collect();
                for pair in pairs {
                    let rel = self.topicity(pair).meet(&other.topicity(pair))?;
                    if !rel.is_undetermined() {
                        result.add($constraint::Topicity(TopicityAst { pair, rel }));
                    }
                }
                let g = self.stereogenicity().meet(&other.stereogenicity())?;
                if !g.is_undetermined() {
                    result.add($constraint::Stereogenicity(g));
                }
                Some(result)
            }

            fn join(&self, other: &Self) -> Self {
                let mut result = Self::new();
                for p in self.ligand_symmetry() {
                    if other.ligand_symmetry().any(|o| o == p) {
                        result.add($constraint::LigandSymmetry(*p));
                    }
                }
                for f in self.fluxionality() {
                    if other.fluxionality().any(|o| o == f) {
                        result.add($constraint::Fluxionality(*f));
                    }
                }
                for t in self.topicities() {
                    if other.topicities().any(|o| o.pair == t.pair) {
                        let rel = t.rel.join(&other.topicity(t.pair));
                        if !rel.is_undetermined() {
                            result.add($constraint::Topicity(TopicityAst {
                                pair: t.pair,
                                rel,
                            }));
                        }
                    }
                }
                if self.contains($kind::Stereogenicity) && other.contains($kind::Stereogenicity) {
                    let g = self.stereogenicity().join(&other.stereogenicity());
                    if !g.is_undetermined() {
                        result.add($constraint::Stereogenicity(g));
                    }
                }
                result
            }

            fn matches(&self, target: &Self) -> bool {
                self.ligand_symmetry()
                    .all(|p| target.ligand_symmetry().any(|t| p.matches(t)))
                    && self
                        .fluxionality()
                        .all(|p| target.fluxionality().any(|t| p.matches(t)))
                    && self
                        .topicities()
                        .all(|t| t.rel.matches(&target.topicity(t.pair)))
                    && self.stereogenicity().matches(&target.stereogenicity())
            }
        }

        impl FromIterator<$constraint> for $constraints {
            fn from_iter<I: IntoIterator<Item = $constraint>>(iter: I) -> Self {
                let mut out = Self::new();
                for c in iter {
                    out.add(c);
                }
                out
            }
        }

        impl IntoIterator for $constraints {
            type Item = $constraint;
            type IntoIter = smallvec::IntoIter<[$constraint; 2]>;

            fn into_iter(self) -> Self::IntoIter {
                self.entries.into_iter()
            }
        }

        impl From<$constraint> for $constraints {
            fn from(c: $constraint) -> Self {
                Self::from_iter([c])
            }
        }

        impl From<Vec<$constraint>> for $constraints {
            fn from(cs: Vec<$constraint>) -> Self {
                Self::from_iter(cs)
            }
        }
    };
}

stereo_constraint! { StereoAtomConstraint, StereoAtomConstraintKind, StereoAtomConstraints }
stereo_constraint! { StereoBondConstraint, StereoBondConstraintKind, StereoBondConstraints }

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

    #[rustfmt::skip]
    #[rstest]
    #[case::singleton_to_lit(TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic])), Ok(TopicityRelationAst::Lit(Topicity::Homotopic)))]
    #[case::litset_polarity_to_notset(TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])), Ok(TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic]))))]
    #[case::full_litset_to_undetermined(TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic, Topicity::Diastereotopic])), Ok(TopicityRelationAst::Undetermined))]
    #[case::empty_litset_err(TopicityRelationAst::LitSet(BTreeSet::new()), Err(Contradiction))]
    #[case::notset_complement_to_lit(TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])), Ok(TopicityRelationAst::Lit(Topicity::Diastereotopic)))]
    #[case::empty_notset_to_undetermined(TopicityRelationAst::NotSet(BTreeSet::new()), Ok(TopicityRelationAst::Undetermined))]
    #[case::full_notset_err(TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic, Topicity::Diastereotopic])), Err(Contradiction))]
    fn test_topicity_relation_ast_canonicalize(
        #[case] input: TopicityRelationAst,
        #[case] expected: Result<TopicityRelationAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rstest]
    #[case::undetermined(TopicityRelationAst::Undetermined)]
    #[case::lit(TopicityRelationAst::Lit(Topicity::Homotopic))]
    #[case::notset(TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])))]
    fn test_topicity_relation_ast_canonicalize_identity(#[case] input: TopicityRelationAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(TopicityRelationAst::Lit(Topicity::Homotopic), Some(Topicity::Homotopic))]
    #[case::undetermined(TopicityRelationAst::Undetermined, None)]
    #[case::notset(TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Homotopic])), None)]
    #[case::litset(TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])), None)]
    fn test_topicity_relation_ast_as_lit(#[case] r: TopicityRelationAst, #[case] expected: Option<Topicity>) {
        assert_eq!(r.as_lit(), expected);
        assert_eq!(r.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(TopicityRelationAst::Undetermined, true)]
    #[case::lit(TopicityRelationAst::Lit(Topicity::Homotopic), false)]
    #[case::notset(TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])), false)]
    fn test_topicity_relation_ast_is_undetermined(#[case] r: TopicityRelationAst, #[case] expected: bool) {
        assert_eq!(r.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(TopicityRelationAst::Undetermined, TopicityRelationAst::Lit(Topicity::Homotopic), Some(TopicityRelationAst::Lit(Topicity::Homotopic)))]
    #[case::lit_eq(TopicityRelationAst::Lit(Topicity::Homotopic), TopicityRelationAst::Lit(Topicity::Homotopic), Some(TopicityRelationAst::Lit(Topicity::Homotopic)))]
    #[case::lit_disjoint(TopicityRelationAst::Lit(Topicity::Homotopic), TopicityRelationAst::Lit(Topicity::Enantiotopic), None)]
    #[case::lit_notset_in(TopicityRelationAst::Lit(Topicity::Homotopic), TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])), Some(TopicityRelationAst::Lit(Topicity::Homotopic)))]
    #[case::lit_notset_out(TopicityRelationAst::Lit(Topicity::Diastereotopic), TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])), None)]
    #[case::notset_notset_to_lit(TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])), TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Enantiotopic])), Some(TopicityRelationAst::Lit(Topicity::Homotopic)))]
    fn test_topicity_relation_ast_meet(
        #[case] a: TopicityRelationAst,
        #[case] b: TopicityRelationAst,
        #[case] expected: Option<TopicityRelationAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und(TopicityRelationAst::Undetermined, TopicityRelationAst::Lit(Topicity::Homotopic), TopicityRelationAst::Undetermined)]
    #[case::lit_eq(TopicityRelationAst::Lit(Topicity::Homotopic), TopicityRelationAst::Lit(Topicity::Homotopic), TopicityRelationAst::Lit(Topicity::Homotopic))]
    #[case::lit_union_to_notset(TopicityRelationAst::Lit(Topicity::Homotopic), TopicityRelationAst::Lit(Topicity::Enantiotopic), TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])))]
    #[case::lit_notset_to_full(TopicityRelationAst::Lit(Topicity::Diastereotopic), TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])), TopicityRelationAst::Undetermined)]
    fn test_topicity_relation_ast_join(
        #[case] a: TopicityRelationAst,
        #[case] b: TopicityRelationAst,
        #[case] expected: TopicityRelationAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(TopicityRelationAst::Undetermined, TopicityRelationAst::Lit(Topicity::Homotopic), true)]
    #[case::lit_und(TopicityRelationAst::Lit(Topicity::Homotopic), TopicityRelationAst::Undetermined, false)]
    #[case::lit_eq(TopicityRelationAst::Lit(Topicity::Homotopic), TopicityRelationAst::Lit(Topicity::Homotopic), true)]
    #[case::lit_neq(TopicityRelationAst::Lit(Topicity::Homotopic), TopicityRelationAst::Lit(Topicity::Enantiotopic), false)]
    #[case::notset_in(TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])), TopicityRelationAst::Lit(Topicity::Homotopic), true)]
    #[case::notset_out(TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])), TopicityRelationAst::Lit(Topicity::Diastereotopic), false)]
    fn test_topicity_relation_ast_matches(
        #[case] pattern: TopicityRelationAst,
        #[case] target: TopicityRelationAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_ligand_symmetry_ast_matches() {
        let perm = OrientedPermutationAst {
            perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])),
            orientation: Orientation::Proper,
        };
        let present = LigandSymmetryAst {
            perm,
            mem: MemOp::In,
        };
        let same = LigandSymmetryAst {
            perm,
            mem: MemOp::In,
        };
        let absent = LigandSymmetryAst {
            perm,
            mem: MemOp::NotIn,
        };
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
    fn test_topicity_ast_matches() {
        let pair = LigandPairAst::new(StereoLigandId(0), StereoLigandId(1));
        let h = TopicityAst {
            pair,
            rel: TopicityRelationAst::Lit(Topicity::Homotopic),
        };
        let e = TopicityAst {
            pair,
            rel: TopicityRelationAst::Lit(Topicity::Enantiotopic),
        };
        let open = TopicityAst {
            pair,
            rel: TopicityRelationAst::Undetermined,
        };
        assert!(open.matches(&h));
        assert!(h.matches(&h));
        assert!(!h.matches(&e));
        // A constraint on a different pair never matches.
        let elsewhere = TopicityAst {
            pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(2)),
            rel: TopicityRelationAst::Undetermined,
        };
        assert!(!elsewhere.matches(&h));
    }

    // `StereogenicityAst` is the second `relation_ast!` instantiation; the macro's
    // full lattice/canonicalize logic is covered by the `TopicityRelationAst` tests
    // above and `test_stereogenicity_ast_lattice_laws` (property). These confirm the
    // instantiation over the `Stereogenicity` domain.
    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereogenicityAst::Lit(Stereogenicity::Stereogenic), Some(Stereogenicity::Stereogenic))]
    #[case::undetermined(StereogenicityAst::Undetermined, None)]
    #[case::notset(StereogenicityAst::NotSet(BTreeSet::from([Stereogenicity::Stereogenic])), None)]
    fn test_stereogenicity_ast_as_lit(#[case] g: StereogenicityAst, #[case] expected: Option<Stereogenicity>) {
        assert_eq!(g.as_lit(), expected);
        assert_eq!(g.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::singleton_to_lit(StereogenicityAst::LitSet(BTreeSet::from([Stereogenicity::Stereogenic])), Ok(StereogenicityAst::Lit(Stereogenicity::Stereogenic)))]
    #[case::full_to_undetermined(StereogenicityAst::LitSet(BTreeSet::from([Stereogenicity::Symmetric, Stereogenicity::Prochiral, Stereogenicity::Stereogenic])), Ok(StereogenicityAst::Undetermined))]
    #[case::empty_err(StereogenicityAst::LitSet(BTreeSet::new()), Err(Contradiction))]
    fn test_stereogenicity_ast_canonicalize(
        #[case] input: StereogenicityAst,
        #[case] expected: Result<StereogenicityAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rstest]
    fn test_stereo_atom_constraints_new() {
        let cs = StereoAtomConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.iter().count(), 0);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::append(
        vec![
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { perm: OrientedPermutationAst { perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, mem: MemOp::In }),
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { perm: OrientedPermutationAst { perm: PermutationAst(Permutation::from_image(4, &[0, 1, 3, 2])), orientation: Orientation::Proper }, mem: MemOp::In }),
        ],
        vec![
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { perm: OrientedPermutationAst { perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, mem: MemOp::In }),
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { perm: OrientedPermutationAst { perm: PermutationAst(Permutation::from_image(4, &[0, 1, 3, 2])), orientation: Orientation::Proper }, mem: MemOp::In }),
        ],
        None,
    )]
    #[case::unique_replace(
        vec![
            StereoAtomConstraint::Stereogenicity(StereogenicityAst::Undetermined),
            StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ],
        vec![StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))],
        Some(StereoAtomConstraint::Stereogenicity(StereogenicityAst::Undetermined)),
    )]
    #[case::keyed_replace(
        vec![
            StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)), rel: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)), rel: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        ],
        vec![StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)), rel: TopicityRelationAst::Lit(Topicity::Homotopic) })],
        Some(StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)), rel: TopicityRelationAst::Lit(Topicity::Enantiotopic) })),
    )]
    #[case::keyed_new_pair(
        vec![
            StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)), rel: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(2)), rel: TopicityRelationAst::Lit(Topicity::Diastereotopic) }),
        ],
        vec![
            StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)), rel: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(2)), rel: TopicityRelationAst::Lit(Topicity::Diastereotopic) }),
        ],
        None,
    )]
    #[case::kind_sorted(
        vec![
            StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
            StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)), rel: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraint::Fluxionality(FluxionalityAst { perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])) }),
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { perm: OrientedPermutationAst { perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, mem: MemOp::In }),
        ],
        vec![
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { perm: OrientedPermutationAst { perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, mem: MemOp::In }),
            StereoAtomConstraint::Fluxionality(FluxionalityAst { perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])) }),
            StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)), rel: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ],
        None,
    )]
    fn test_stereo_atom_constraints_add(
        #[case] adds: Vec<StereoAtomConstraint>,
        #[case] expected: Vec<StereoAtomConstraint>,
        #[case] last_return: Option<StereoAtomConstraint>,
    ) {
        let mut cs = StereoAtomConstraints::new();
        let mut returned = None;
        for c in adds {
            returned = cs.add(c);
        }
        assert_eq!(cs.iter().cloned().collect::<Vec<_>>(), expected);
        assert_eq!(returned, last_return);
    }

    #[rstest]
    fn test_stereo_atom_constraints_meet() {
        let pair = LigandPairAst::new(StereoLigandId(0), StereoLigandId(1));
        let p1 = LigandSymmetryAst {
            perm: OrientedPermutationAst {
                perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])),
                orientation: Orientation::Proper,
            },
            mem: MemOp::In,
        };
        let p2 = LigandSymmetryAst {
            perm: OrientedPermutationAst {
                perm: PermutationAst(Permutation::from_image(4, &[0, 1, 3, 2])),
                orientation: Orientation::Proper,
            },
            mem: MemOp::In,
        };

        let mut a = StereoAtomConstraints::new();
        a.add(StereoAtomConstraint::LigandSymmetry(p1));
        a.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            rel: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])),
        }));

        let mut b = StereoAtomConstraints::new();
        b.add(StereoAtomConstraint::LigandSymmetry(p1));
        b.add(StereoAtomConstraint::LigandSymmetry(p2));
        b.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            rel: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Homotopic])),
        }));
        b.add(StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)));

        let m = a.meet(&b).unwrap();
        // #p union+dedup.
        assert_eq!(
            m.ligand_symmetry().copied().collect::<Vec<_>>(),
            vec![p1, p2]
        );
        // #o per-pair value-meet: {H,E} ∩ {E,D} = {E}.
        assert_eq!(
            m.topicity(pair),
            TopicityRelationAst::Lit(Topicity::Enantiotopic)
        );
        // #g carried through from the side that has it.
        assert_eq!(
            m.stereogenicity(),
            StereogenicityAst::Lit(Stereogenicity::Stereogenic),
        );
    }

    #[rstest]
    fn test_stereo_atom_constraints_meet_error() {
        let pair = LigandPairAst::new(StereoLigandId(0), StereoLigandId(1));
        let mut a = StereoAtomConstraints::new();
        a.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            rel: TopicityRelationAst::Lit(Topicity::Homotopic),
        }));
        let mut b = StereoAtomConstraints::new();
        b.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            rel: TopicityRelationAst::Lit(Topicity::Enantiotopic),
        }));
        // Disjoint relations on the same pair contradict.
        assert_eq!(a.meet(&b), None);
    }

    #[rstest]
    fn test_stereo_atom_constraints_join() {
        let pair = LigandPairAst::new(StereoLigandId(0), StereoLigandId(1));
        let p1 = LigandSymmetryAst {
            perm: OrientedPermutationAst {
                perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])),
                orientation: Orientation::Proper,
            },
            mem: MemOp::In,
        };
        let p2 = LigandSymmetryAst {
            perm: OrientedPermutationAst {
                perm: PermutationAst(Permutation::from_image(4, &[0, 1, 3, 2])),
                orientation: Orientation::Proper,
            },
            mem: MemOp::In,
        };

        let mut a = StereoAtomConstraints::new();
        a.add(StereoAtomConstraint::LigandSymmetry(p1));
        a.add(StereoAtomConstraint::LigandSymmetry(p2));
        a.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            rel: TopicityRelationAst::Lit(Topicity::Homotopic),
        }));

        let mut b = StereoAtomConstraints::new();
        b.add(StereoAtomConstraint::LigandSymmetry(p1));
        b.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            rel: TopicityRelationAst::Lit(Topicity::Enantiotopic),
        }));

        let j = a.join(&b);
        // #p intersection.
        assert_eq!(j.ligand_symmetry().copied().collect::<Vec<_>>(), vec![p1]);
        // #o per-pair value-join: {H} ∪ {E} = {H,E}.
        assert_eq!(
            j.topicity(pair),
            TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])),
        );
    }

    #[rstest]
    fn test_stereo_atom_constraints_matches() {
        let pair = LigandPairAst::new(StereoLigandId(0), StereoLigandId(1));
        let p1 = LigandSymmetryAst {
            perm: OrientedPermutationAst {
                perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])),
                orientation: Orientation::Proper,
            },
            mem: MemOp::In,
        };

        let mut pattern = StereoAtomConstraints::new();
        pattern.add(StereoAtomConstraint::LigandSymmetry(p1));
        pattern.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            rel: TopicityRelationAst::Lit(Topicity::Enantiotopic),
        }));

        let mut target = StereoAtomConstraints::new();
        target.add(StereoAtomConstraint::LigandSymmetry(p1));
        target.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            rel: TopicityRelationAst::Lit(Topicity::Enantiotopic),
        }));
        assert!(pattern.matches(&target));

        // An empty target leaves both predicates unmatched.
        assert!(!pattern.matches(&StereoAtomConstraints::new()));
        // An empty pattern matches anything.
        assert!(StereoAtomConstraints::new().matches(&target));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoAtomConstraints::new(), true, true)]
    #[case::ligand_symmetry(
        StereoAtomConstraints::from(StereoAtomConstraint::LigandSymmetry(
            LigandSymmetryAst {
                perm: OrientedPermutationAst {
                    perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])),
                    orientation: Orientation::Proper,
                },
                mem: MemOp::In,
            },
        )),
        false,
        true,
    )]
    #[case::topicity_open(
        StereoAtomConstraints::from(StereoAtomConstraint::Topicity(TopicityAst {
            pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)),
            rel: TopicityRelationAst::Undetermined,
        })),
        true,
        false,
    )]
    #[case::stereogenicity(
        StereoAtomConstraints::from(StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        false,
        true,
    )]
    fn test_stereo_atom_constraints_is_undetermined_is_ground(
        #[case] cs: StereoAtomConstraints,
        #[case] is_undetermined: bool,
        #[case] is_ground: bool,
    ) {
        assert_eq!(cs.is_undetermined(), is_undetermined);
        assert_eq!(cs.is_ground(), is_ground);
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
        assert_eq!(cs.iter().count(), 0);
    }

    #[rstest]
    fn test_stereo_bond_constraints_add() {
        let mut cs = StereoBondConstraints::new();
        let f = FluxionalityAst {
            perm: PermutationAst(Permutation::from_image(4, &[1, 0, 2, 3])),
        };
        assert_eq!(cs.add(StereoBondConstraint::Fluxionality(f)), None);
        assert_eq!(cs.fluxionality().copied().collect::<Vec<_>>(), vec![f]);
    }

    #[rstest]
    fn test_stereo_bond_constraints_from_iter() {
        let cs: StereoBondConstraints = empty().collect();
        assert_eq!(cs, StereoBondConstraints::new());
    }
}
