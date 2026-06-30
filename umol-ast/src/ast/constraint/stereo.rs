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
use super::super::id::StereoLigandPosition;
use super::super::operators::MemOp;
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::stereo::{Stereogenicity, Topicity};
use super::super::traits::{AsLit, Canonicalize, Lattice};

/// Ligand permutation literal.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct LigandPermutation(pub Permutation);

impl LigandPermutation {
    /// A pattern matches a target iff they are the same permutation.
    pub fn matches(&self, target: &Self) -> bool {
        self.0 == target.0
    }
}

// Ligand permutation with proper/improper grade.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OrientedLigandPermutation {
    pub permutation: LigandPermutation,
    pub orientation: Orientation,
}

impl OrientedLigandPermutation {
    pub fn matches(&self, target: &Self) -> bool {
        self.permutation.matches(&target.permutation) && self.orientation == target.orientation
    }
}

/// The key of a per-pair topicity constraint: an unordered pair of ligand
/// positions, normalized so the lower position is `first`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StereoLigandPair {
    first: StereoLigandPosition,
    second: StereoLigandPosition,
}

impl StereoLigandPair {
    /// Normalizes the pair so the lower position is `first`.
    pub fn new(a: StereoLigandPosition, b: StereoLigandPosition) -> Self {
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

    pub fn first(self) -> StereoLigandPosition {
        self.first
    }

    pub fn second(self) -> StereoLigandPosition {
        self.second
    }
}

/// Generates a subset-lattice over finite domain enum.
macro_rules! relation_ast {
    ($name:ident, $domain:ty) => {
        #[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub enum $name {
            #[default]
            Undetermined,
            Lit($domain),
            LitSet(BTreeSet<$domain>),
            NotSet(BTreeSet<$domain>),
        }

        impl $name {
            pub fn undetermined() -> Self {
                Self::Undetermined
            }

            pub fn lit(value: $domain) -> Self {
                Self::Lit(value)
            }

            pub fn lit_set(values: impl IntoIterator<Item = $domain>) -> Self {
                Self::LitSet(values.into_iter().collect())
            }

            pub fn not(value: $domain) -> Self {
                Self::NotSet(BTreeSet::from([value]))
            }

            pub fn not_set(values: impl IntoIterator<Item = $domain>) -> Self {
                Self::NotSet(values.into_iter().collect())
            }

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
        }

        impl From<$domain> for $name {
            fn from(value: $domain) -> Self {
                Self::Lit(value)
            }
        }

        impl Canonicalize for $name {
            /// Finite-domain canonical form over `to_set`: empty → `Err`, full →
            /// `Undetermined`, singleton → `Lit`, else the smaller of positive /
            /// complement (tiebreak positive).
            fn canonicalize(self) -> Result<Self, Contradiction> {
                let set = self.to_set();
                let domain: BTreeSet<$domain> = <$domain as VariantArray>::VARIANTS
                    .iter()
                    .copied()
                    .collect();
                if set.is_empty() {
                    Err(Contradiction)
                } else if set == domain {
                    Ok(Self::Undetermined)
                } else if set.len() == 1 {
                    Ok(Self::Lit(set.into_iter().next().unwrap()))
                } else {
                    let complement: BTreeSet<$domain> = domain.difference(&set).copied().collect();
                    Ok(if set.len() <= complement.len() {
                        Self::LitSet(set)
                    } else {
                        Self::NotSet(complement)
                    })
                }
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

            /// Intersection of the admissible sets, folded by `canonicalize` (∅ → `None`).
            fn meet(&self, other: &Self) -> Option<Self> {
                Self::LitSet(
                    self.to_set()
                        .intersection(&other.to_set())
                        .copied()
                        .collect(),
                )
                .canonicalize()
                .ok()
            }

            fn join(&self, other: &Self) -> Self {
                Self::LitSet(self.to_set().union(&other.to_set()).copied().collect())
                    .canonicalize()
                    .unwrap_or(Self::Undetermined)
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
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LigandSymmetryAst {
    pub permutation: OrientedLigandPermutation,
    pub member: MemOp,
}

impl LigandSymmetryAst {
    pub fn matches(&self, target: &Self) -> bool {
        self.permutation.matches(&target.permutation) && self.member == target.member
    }
}

/// Fluxionality move: proper ligand permutation realized by dynamics. Non-unique.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FluxionalityAst {
    pub permutation: LigandPermutation,
}

impl FluxionalityAst {
    pub fn matches(&self, target: &Self) -> bool {
        self.permutation.matches(&target.permutation)
    }
}

/// Per-pair topicity constraint: relation between pair of ligands.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TopicityAst {
    pub pair: StereoLigandPair,
    pub relation: TopicityRelationAst,
}

impl TopicityAst {
    /// Matches-only (not a `Lattice`): different pairs are incomparable, so there
    /// is no global top. `matches` = same `pair` and the per-pair `rel` matches;
    /// the per-pair lattice lives in [`TopicityRelationAst`].
    pub fn matches(&self, target: &Self) -> bool {
        self.pair == target.pair && self.relation.matches(&target.relation)
    }
}

/// Per-element stereo constraint for stereo atoms and bonds.
macro_rules! stereo_constraint {
    ($constraint:ident, $kind:ident, $key:ident, $constraints:ident) => {
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
        pub enum $kind {
            LigandSymmetry,
            Fluxionality,
            Topicity,
            Stereogenicity,
        }

        /// Entry identity: discriminant + sub-key. Variant order matches `$kind`,
        /// so `Ord` agrees with `kind as u8`. `mem` is LigandSymmetry's value, not
        /// part of its key, so conflicting `In`/`NotIn` on one permutation contradict.
        #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub enum $key {
            LigandSymmetry(OrientedLigandPermutation),
            Fluxionality(LigandPermutation),
            Topicity(StereoLigandPair),
            Stereogenicity,
        }

        impl $key {
            pub fn kind(self) -> $kind {
                match self {
                    Self::LigandSymmetry(_) => $kind::LigandSymmetry,
                    Self::Fluxionality(_) => $kind::Fluxionality,
                    Self::Topicity(_) => $kind::Topicity,
                    Self::Stereogenicity => $kind::Stereogenicity,
                }
            }
        }

        #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

            /// Entry identity for order/dedup: `kind()` plus the sub-key
            /// (Topicity's pair, the permutation for `#f`/`#p`).
            pub fn key(&self) -> $key {
                match self {
                    Self::LigandSymmetry(ls) => $key::LigandSymmetry(ls.permutation),
                    Self::Fluxionality(f) => $key::Fluxionality(f.permutation),
                    Self::Topicity(t) => $key::Topicity(t.pair),
                    Self::Stereogenicity(_) => $key::Stereogenicity,
                }
            }

            /// Whether at most one constraint of this kind may be stored.
            pub fn is_unique(&self) -> bool {
                matches!(self, Self::Stereogenicity(_))
            }

            pub fn is_undetermined(&self) -> bool {
                match self {
                    Self::LigandSymmetry(_) | Self::Fluxionality(_) => false,
                    Self::Topicity(t) => t.relation.is_undetermined(),
                    Self::Stereogenicity(g) => g.is_undetermined(),
                }
            }

            /// Frame-relative ligand positions carry no atom ids, so compact is a no-op.
            pub fn compact(self, _compaction: &IdCompaction) -> Option<Self> {
                Some(self)
            }

            /// Frame-relative ligand positions carry no atom ids, so remap is a no-op.
            pub(crate) fn remap(self, _map: &IdRemapping) -> Self {
                self
            }
        }

        impl Canonicalize for $constraint {
            /// Canonicalize the inner relation value; `#f`/`#p` have no
            /// canonicalizable inner value (permutation/member are atomic).
            fn canonicalize(self) -> Result<Self, Contradiction> {
                Ok(match self {
                    Self::LigandSymmetry(ls) => Self::LigandSymmetry(ls),
                    Self::Fluxionality(f) => Self::Fluxionality(f),
                    Self::Topicity(t) => Self::Topicity(TopicityAst {
                        pair: t.pair,
                        relation: t.relation.canonicalize()?,
                    }),
                    Self::Stereogenicity(g) => Self::Stereogenicity(g.canonicalize()?),
                })
            }
        }

        /// Stereo constraint container for stereo atoms and bonds.
        #[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

            /// Membership polarity asserted for `permutation`, if any.
            fn ligand_mem(&self, permutation: OrientedLigandPermutation) -> Option<MemOp> {
                self.ligand_symmetry()
                    .find(|ls| ls.permutation == permutation)
                    .map(|ls| ls.member)
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
            pub fn topicity(&self, pair: StereoLigandPair) -> TopicityRelationAst {
                self.topicities()
                    .find(|t| t.pair == pair)
                    .map(|t| t.relation.clone())
                    .unwrap_or_default()
            }

            /// Stored stereogenicity relation (unique); `Undetermined` if absent.
            pub fn stereogenicity(&self) -> StereogenicityAst {
                match self.get($kind::Stereogenicity) {
                    Some($constraint::Stereogenicity(g)) => g.clone(),
                    _ => StereogenicityAst::Undetermined,
                }
            }

            /// Insert at the `key()`-sorted position: the unique kind
            /// (`Stereogenicity`) replaces the same-key entry (returning it);
            /// `#o`/`#f`/`#p` append, leaving same-key duplicates for
            /// `meet`/`canonicalize` to merge (lazy dedup).
            pub fn add(&mut self, c: $constraint) -> Option<$constraint> {
                match self.find_by_key(c.key()) {
                    Ok(i) if c.is_unique() => Some(replace(&mut self.entries[i], c)),
                    Ok(i) => {
                        let end = i + self.entries[i..]
                            .iter()
                            .take_while(|e| e.key() == c.key())
                            .count();
                        self.entries.insert(end, c);
                        None
                    }
                    Err(i) => {
                        self.entries.insert(i, c);
                        None
                    }
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
            pub fn compact(self, _compaction: &IdCompaction) -> Self {
                self
            }

            fn find(&self, kind: $kind) -> Result<usize, usize> {
                self.entries
                    .binary_search_by_key(&(kind as u8), |c| c.kind() as u8)
            }

            fn find_by_key(&self, key: $key) -> Result<usize, usize> {
                self.entries.binary_search_by(|c| c.key().cmp(&key))
            }

            pub fn contains_key(&self, key: $key) -> bool {
                self.find_by_key(key).is_ok()
            }

            pub fn get_by_key(&self, key: $key) -> Option<&$constraint> {
                self.find_by_key(key).ok().map(|i| &self.entries[i])
            }

            pub fn get_by_key_mut(&mut self, key: $key) -> Option<&mut $constraint> {
                self.find_by_key(key).ok().map(|i| &mut self.entries[i])
            }

            pub fn remove_by_key(&mut self, key: $key) -> Option<$constraint> {
                self.find_by_key(key).ok().map(|i| self.entries.remove(i))
            }

            /// Merge two same-key entries by value-`meet`: relations meet (`None`
            /// on contradiction), `#p` requires equal `mem` (else `None`), `#f`
            /// dedups. Caller guarantees `a.key() == b.key()`.
            fn merge_same_key(a: $constraint, b: $constraint) -> Option<$constraint> {
                match (a, b) {
                    ($constraint::Topicity(x), $constraint::Topicity(y)) => {
                        Some($constraint::Topicity(TopicityAst {
                            pair: x.pair,
                            relation: x.relation.meet(&y.relation)?,
                        }))
                    }
                    ($constraint::Stereogenicity(x), $constraint::Stereogenicity(y)) => {
                        Some($constraint::Stereogenicity(x.meet(&y)?))
                    }
                    ($constraint::LigandSymmetry(x), $constraint::LigandSymmetry(y)) => {
                        (x.member == y.member).then_some($constraint::LigandSymmetry(x))
                    }
                    ($constraint::Fluxionality(x), $constraint::Fluxionality(_)) => {
                        Some($constraint::Fluxionality(x))
                    }
                    _ => unreachable!("merge_same_key called with differing keys"),
                }
            }
        }

        impl Canonicalize for $constraints {
            /// Sort by `key()`, canonicalize each value, merge same-key entries
            /// by value-`meet` (`Err` on contradiction), drop vacuous entries.
            fn canonicalize(self) -> Result<Self, Contradiction> {
                let mut input = self.entries;
                input.sort_by_key(|c| c.key());
                let mut out: SmallVec<[$constraint; 2]> = SmallVec::new();
                for c in input {
                    let c = c.canonicalize()?;
                    if out.last().map(|p| p.key()) == Some(c.key()) {
                        let merged =
                            Self::merge_same_key(out.pop().unwrap(), c).ok_or(Contradiction)?;
                        out.push(merged);
                    } else {
                        out.push(c);
                    }
                }
                out.retain(|c| !c.is_undetermined());
                Ok(Self { entries: out })
            }
        }

        impl Lattice for $constraints {
            fn is_undetermined(&self) -> bool {
                self.entries.iter().all(|c| c.is_undetermined())
            }

            fn is_ground(&self) -> bool {
                self.entries.iter().all(|c| match c {
                    $constraint::LigandSymmetry(_) | $constraint::Fluxionality(_) => true,
                    $constraint::Topicity(t) => t.relation.is_ground(),
                    $constraint::Stereogenicity(g) => g.is_ground(),
                })
            }

            fn meet(&self, other: &Self) -> Option<Self> {
                let mut result = Self::new();
                let permutations: BTreeSet<OrientedLigandPermutation> = self
                    .ligand_symmetry()
                    .chain(other.ligand_symmetry())
                    .map(|ls| ls.permutation)
                    .collect();
                for permutation in permutations {
                    let member = match (self.ligand_mem(permutation), other.ligand_mem(permutation))
                    {
                        (Some(a), Some(b)) if a != b => return None,
                        (Some(m), _) | (_, Some(m)) => m,
                        (None, None) => continue,
                    };
                    result.add($constraint::LigandSymmetry(LigandSymmetryAst {
                        permutation,
                        member,
                    }));
                }
                for f in self.fluxionality().chain(other.fluxionality()) {
                    let entry = $constraint::Fluxionality(*f);
                    if !result.contains_entry(&entry) {
                        result.add(entry);
                    }
                }
                let pairs: BTreeSet<StereoLigandPair> = self
                    .topicities()
                    .chain(other.topicities())
                    .map(|t| t.pair)
                    .collect();
                for pair in pairs {
                    let relation = self.topicity(pair).meet(&other.topicity(pair))?;
                    if !relation.is_undetermined() {
                        result.add($constraint::Topicity(TopicityAst { pair, relation }));
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
                        let relation = t.relation.join(&other.topicity(t.pair));
                        if !relation.is_undetermined() {
                            result.add($constraint::Topicity(TopicityAst {
                                pair: t.pair,
                                relation,
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
                        .all(|t| t.relation.matches(&target.topicity(t.pair)))
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

stereo_constraint! { StereoAtomConstraint, StereoAtomConstraintKind, StereoAtomConstraintKey, StereoAtomConstraints }
stereo_constraint! { StereoBondConstraint, StereoBondConstraintKind, StereoBondConstraintKey, StereoBondConstraints }

#[cfg(test)]
mod tests {
    use std::iter::empty;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    fn test_permutation_ast_matches() {
        let a = LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3]));
        let b = LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3]));
        let c = LigandPermutation(Permutation::identity(4));
        assert!(a.matches(&b));
        assert!(!a.matches(&c));
    }

    #[rstest]
    fn test_oriented_permutation_ast_matches() {
        let permutation = LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3]));
        let proper = OrientedLigandPermutation {
            permutation,
            orientation: Orientation::Proper,
        };
        let same = OrientedLigandPermutation {
            permutation,
            orientation: Orientation::Proper,
        };
        let flipped = OrientedLigandPermutation {
            permutation,
            orientation: Orientation::Improper,
        };
        let other = OrientedLigandPermutation {
            permutation: LigandPermutation(Permutation::identity(4)),
            orientation: Orientation::Proper,
        };
        assert!(proper.matches(&same));
        assert!(!proper.matches(&flipped));
        assert!(!proper.matches(&other));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ordered(StereoLigandPosition(1), StereoLigandPosition(2), StereoLigandPosition(1), StereoLigandPosition(2))]
    #[case::reversed(StereoLigandPosition(2), StereoLigandPosition(1), StereoLigandPosition(1), StereoLigandPosition(2))]
    #[case::equal(StereoLigandPosition(3), StereoLigandPosition(3), StereoLigandPosition(3), StereoLigandPosition(3))]
    fn test_ligand_pair_ast_new(
        #[case] a: StereoLigandPosition,
        #[case] b: StereoLigandPosition,
        #[case] first: StereoLigandPosition,
        #[case] second: StereoLigandPosition,
    ) {
        let pair = StereoLigandPair::new(a, b);
        assert_eq!(pair.first(), first);
        assert_eq!(pair.second(), second);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(TopicityRelationAst::undetermined(), TopicityRelationAst::Undetermined)]
    #[case::lit(TopicityRelationAst::lit(Topicity::Homotopic), TopicityRelationAst::Lit(Topicity::Homotopic))]
    #[case::lit_set(TopicityRelationAst::lit_set([Topicity::Homotopic, Topicity::Enantiotopic]), TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])))]
    #[case::lit_set_singleton_raw(TopicityRelationAst::lit_set([Topicity::Homotopic]), TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic])))]
    #[case::not(TopicityRelationAst::not(Topicity::Homotopic), TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Homotopic])))]
    #[case::not_set(TopicityRelationAst::not_set([Topicity::Homotopic, Topicity::Enantiotopic]), TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])))]
    fn test_topicity_relation_ast_constructors(#[case] actual: TopicityRelationAst, #[case] expected: TopicityRelationAst) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::homotopic(Topicity::Homotopic, TopicityRelationAst::Lit(Topicity::Homotopic))]
    #[case::diastereotopic(
        Topicity::Diastereotopic,
        TopicityRelationAst::Lit(Topicity::Diastereotopic)
    )]
    fn test_topicity_relation_ast_from(
        #[case] value: Topicity,
        #[case] expected: TopicityRelationAst,
    ) {
        assert_eq!(TopicityRelationAst::from(value), expected);
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
        let permutation = OrientedLigandPermutation {
            permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
            orientation: Orientation::Proper,
        };
        let present = LigandSymmetryAst {
            permutation,
            member: MemOp::In,
        };
        let same = LigandSymmetryAst {
            permutation,
            member: MemOp::In,
        };
        let absent = LigandSymmetryAst {
            permutation,
            member: MemOp::NotIn,
        };
        let other = LigandSymmetryAst {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::identity(4)),
                orientation: Orientation::Proper,
            },
            member: MemOp::In,
        };
        assert!(present.matches(&same));
        assert!(!present.matches(&absent)); // different membership op
        assert!(!present.matches(&other)); // different permutation
    }

    #[rstest]
    fn test_fluxionality_ast_matches() {
        let a = FluxionalityAst {
            permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
        };
        let same = FluxionalityAst {
            permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
        };
        let other = FluxionalityAst {
            permutation: LigandPermutation(Permutation::identity(4)),
        };
        assert!(a.matches(&same));
        assert!(!a.matches(&other));
    }

    #[rstest]
    fn test_topicity_ast_matches() {
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1));
        let h = TopicityAst {
            pair,
            relation: TopicityRelationAst::Lit(Topicity::Homotopic),
        };
        let e = TopicityAst {
            pair,
            relation: TopicityRelationAst::Lit(Topicity::Enantiotopic),
        };
        let open = TopicityAst {
            pair,
            relation: TopicityRelationAst::Undetermined,
        };
        assert!(open.matches(&h));
        assert!(h.matches(&h));
        assert!(!h.matches(&e));
        // A constraint on a different pair never matches.
        let elsewhere = TopicityAst {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)),
            relation: TopicityRelationAst::Undetermined,
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
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, member: MemOp::In }),
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[0, 1, 3, 2])), orientation: Orientation::Proper }, member: MemOp::In }),
        ],
        vec![
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[0, 1, 3, 2])), orientation: Orientation::Proper }, member: MemOp::In }),
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, member: MemOp::In }),
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
    #[case::keyed_appends_lazy(
        vec![
            StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        ],
        vec![
            StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        ],
        None,
    )]
    #[case::keyed_new_pair(
        vec![
            StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationAst::Lit(Topicity::Diastereotopic) }),
        ],
        vec![
            StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationAst::Lit(Topicity::Diastereotopic) }),
        ],
        None,
    )]
    #[case::kind_sorted(
        vec![
            StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
            StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraint::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])) }),
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, member: MemOp::In }),
        ],
        vec![
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, member: MemOp::In }),
            StereoAtomConstraint::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])) }),
            StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
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

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(
        StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, member: MemOp::In }),
        StereoAtomConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }))]
    #[case::fluxionality(
        StereoAtomConstraint::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])) }),
        StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3]))))]
    #[case::topicity(
        StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))))]
    #[case::stereogenicity(
        StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintKey::Stereogenicity)]
    fn test_stereo_atom_constraint_key(
        #[case] c: StereoAtomConstraint,
        #[case] expected: StereoAtomConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(StereoAtomConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }), StereoAtomConstraintKind::LigandSymmetry)]
    #[case::fluxionality(StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::identity(4))), StereoAtomConstraintKind::Fluxionality)]
    #[case::topicity(StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))), StereoAtomConstraintKind::Topicity)]
    #[case::stereogenicity(StereoAtomConstraintKey::Stereogenicity, StereoAtomConstraintKind::Stereogenicity)]
    fn test_stereo_atom_constraint_key_kind(
        #[case] key: StereoAtomConstraintKey,
        #[case] expected: StereoAtomConstraintKind,
    ) {
        assert_eq!(key.kind(), expected);
    }

    #[rstest]
    fn test_stereo_atom_constraints_by_key() {
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1));
        let topicity = StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::Lit(Topicity::Homotopic),
        });
        let mut cs = StereoAtomConstraints::new();
        cs.add(topicity.clone());
        cs.add(StereoAtomConstraint::Stereogenicity(
            StereogenicityAst::Lit(Stereogenicity::Stereogenic),
        ));

        assert!(cs.contains_key(StereoAtomConstraintKey::Topicity(pair)));
        assert!(
            !cs.contains_key(StereoAtomConstraintKey::Topicity(StereoLigandPair::new(
                StereoLigandPosition(0),
                StereoLigandPosition(2)
            )))
        );
        assert_eq!(
            cs.get_by_key(StereoAtomConstraintKey::Topicity(pair)),
            Some(&topicity)
        );

        *cs.get_by_key_mut(StereoAtomConstraintKey::Stereogenicity)
            .unwrap() =
            StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Symmetric));
        assert_eq!(
            cs.stereogenicity(),
            StereogenicityAst::Lit(Stereogenicity::Symmetric)
        );

        assert_eq!(
            cs.remove_by_key(StereoAtomConstraintKey::Topicity(pair)),
            Some(topicity)
        );
        assert!(!cs.contains_key(StereoAtomConstraintKey::Topicity(pair)));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity_litset_singleton(
        StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic])) }),
        Ok(StereoAtomConstraint::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })))]
    #[case::fluxionality_identity(
        StereoAtomConstraint::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::identity(4)) }),
        Ok(StereoAtomConstraint::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::identity(4)) })))]
    fn test_stereo_atom_constraint_canonicalize(
        #[case] c: StereoAtomConstraint,
        #[case] expected: Result<StereoAtomConstraint, Contradiction>,
    ) {
        assert_eq!(c.canonicalize(), expected);
    }

    #[rstest]
    fn test_stereo_atom_constraints_canonicalize_merges_topicity() {
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1));
        let mut cs = StereoAtomConstraints::new();
        cs.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])),
        }));
        cs.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Homotopic])),
        }));
        assert_eq!(
            cs.canonicalize().unwrap(),
            StereoAtomConstraints::from_iter([StereoAtomConstraint::Topicity(TopicityAst {
                pair,
                relation: TopicityRelationAst::Lit(Topicity::Enantiotopic),
            })]),
        );
    }

    #[rstest]
    fn test_stereo_atom_constraints_canonicalize_ligand_mem_conflict() {
        let permutation = OrientedLigandPermutation {
            permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
            orientation: Orientation::Proper,
        };
        let mut cs = StereoAtomConstraints::new();
        cs.add(StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst {
            permutation,
            member: MemOp::In,
        }));
        cs.add(StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst {
            permutation,
            member: MemOp::NotIn,
        }));
        assert_eq!(cs.canonicalize(), Err(Contradiction));
    }

    #[rstest]
    fn test_stereo_atom_constraints_canonicalize_drops_vacuous() {
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1));
        let mut cs = StereoAtomConstraints::new();
        cs.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::Undetermined,
        }));
        cs.add(StereoAtomConstraint::Stereogenicity(
            StereogenicityAst::Lit(Stereogenicity::Stereogenic),
        ));
        assert_eq!(
            cs.canonicalize().unwrap(),
            StereoAtomConstraints::from_iter([StereoAtomConstraint::Stereogenicity(
                StereogenicityAst::Lit(Stereogenicity::Stereogenic)
            )]),
        );
    }

    #[rstest]
    fn test_stereo_atom_constraints_meet() {
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1));
        let p1 = LigandSymmetryAst {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
                orientation: Orientation::Proper,
            },
            member: MemOp::In,
        };
        let p2 = LigandSymmetryAst {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::from_image(4, &[0, 1, 3, 2])),
                orientation: Orientation::Proper,
            },
            member: MemOp::In,
        };

        let mut a = StereoAtomConstraints::new();
        a.add(StereoAtomConstraint::LigandSymmetry(p1));
        a.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])),
        }));

        let mut b = StereoAtomConstraints::new();
        b.add(StereoAtomConstraint::LigandSymmetry(p1));
        b.add(StereoAtomConstraint::LigandSymmetry(p2));
        b.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Homotopic])),
        }));
        b.add(StereoAtomConstraint::Stereogenicity(
            StereogenicityAst::Lit(Stereogenicity::Stereogenic),
        ));

        let m = a.meet(&b).unwrap();
        // #p union+dedup, key-sorted (p2's permutation sorts before p1's).
        assert_eq!(
            m.ligand_symmetry().copied().collect::<Vec<_>>(),
            vec![p2, p1]
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
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1));
        let mut a = StereoAtomConstraints::new();
        a.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::Lit(Topicity::Homotopic),
        }));
        let mut b = StereoAtomConstraints::new();
        b.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::Lit(Topicity::Enantiotopic),
        }));
        // Disjoint relations on the same pair contradict.
        assert_eq!(a.meet(&b), None);
    }

    #[rstest]
    fn test_stereo_atom_constraints_join() {
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1));
        let p1 = LigandSymmetryAst {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
                orientation: Orientation::Proper,
            },
            member: MemOp::In,
        };
        let p2 = LigandSymmetryAst {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::from_image(4, &[0, 1, 3, 2])),
                orientation: Orientation::Proper,
            },
            member: MemOp::In,
        };

        let mut a = StereoAtomConstraints::new();
        a.add(StereoAtomConstraint::LigandSymmetry(p1));
        a.add(StereoAtomConstraint::LigandSymmetry(p2));
        a.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::Lit(Topicity::Homotopic),
        }));

        let mut b = StereoAtomConstraints::new();
        b.add(StereoAtomConstraint::LigandSymmetry(p1));
        b.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::Lit(Topicity::Enantiotopic),
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
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1));
        let p1 = LigandSymmetryAst {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
                orientation: Orientation::Proper,
            },
            member: MemOp::In,
        };

        let mut pattern = StereoAtomConstraints::new();
        pattern.add(StereoAtomConstraint::LigandSymmetry(p1));
        pattern.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::Lit(Topicity::Enantiotopic),
        }));

        let mut target = StereoAtomConstraints::new();
        target.add(StereoAtomConstraint::LigandSymmetry(p1));
        target.add(StereoAtomConstraint::Topicity(TopicityAst {
            pair,
            relation: TopicityRelationAst::Lit(Topicity::Enantiotopic),
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
                permutation: OrientedLigandPermutation {
                    permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
                    orientation: Orientation::Proper,
                },
                member: MemOp::In,
            },
        )),
        false,
        true,
    )]
    #[case::topicity_open(
        StereoAtomConstraints::from(StereoAtomConstraint::Topicity(TopicityAst {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
            relation: TopicityRelationAst::Undetermined,
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
            permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
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
