//! Stereo-element constraints.

use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::hash::Hash;
use std::mem;
use std::slice::Iter;

use smallvec::SmallVec;
use strum::VariantArray;
use umol_perm::{Orientation, Permutation};

use super::super::boolean::BooleanAst;
use super::super::error::{Contradiction, NoJoin};
use super::super::id::StereoLigandPosition;
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::stereo::{Stereogenicity, Topicity};
use super::super::traits::{AsLit, Canonicalize, Lattice};

/// Stereo atom and bond constraints.
macro_rules! stereo_constraint {
    ($constraint:ident, $key:ident, $constraints:ident) => {
        /// Stereo constraint key: discriminant, unique within constraint container.
        #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub enum $key {
            LigandSymmetry(OrientedLigandPermutation),
            Fluxionality(LigandPermutation),
            Topicity(StereoLigandPair),
            Stereogenicity,
        }

        #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub enum $constraint {
            LigandSymmetry(LigandSymmetryAst),
            Fluxionality(FluxionalityAst),
            Topicity(TopicityAst),
            Stereogenicity(StereogenicityAst),
        }

        impl $constraint {
            /// Stereo constraint key, unique within a `StereoConstraints` container.
            pub fn key(&self) -> $key {
                match self {
                    Self::LigandSymmetry(ls) => $key::LigandSymmetry(ls.permutation),
                    Self::Fluxionality(f) => $key::Fluxionality(f.permutation),
                    Self::Topicity(t) => $key::Topicity(t.pair),
                    Self::Stereogenicity(_) => $key::Stereogenicity,
                }
            }

            /// Vacuous form of constraint key, used for removal.
            pub fn as_undetermined(&self) -> Self {
                match self {
                    Self::LigandSymmetry(ls) => Self::LigandSymmetry(LigandSymmetryAst {
                        permutation: ls.permutation,
                        present: BooleanAst::Undetermined,
                    }),
                    Self::Fluxionality(f) => Self::Fluxionality(FluxionalityAst {
                        permutation: f.permutation,
                        present: BooleanAst::Undetermined,
                    }),
                    Self::Topicity(t) => Self::Topicity(TopicityAst {
                        pair: t.pair,
                        relation: TopicityRelationAst::Undetermined,
                    }),
                    Self::Stereogenicity(_) => {
                        Self::Stereogenicity(StereogenicityAst::Undetermined)
                    }
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

        impl Lattice for $constraint {
            fn is_undetermined(&self) -> bool {
                match self {
                    Self::LigandSymmetry(ls) => ls.present.is_undetermined(),
                    Self::Fluxionality(f) => f.present.is_undetermined(),
                    Self::Topicity(t) => t.relation.is_undetermined(),
                    Self::Stereogenicity(g) => g.is_undetermined(),
                }
            }

            fn is_ground(&self) -> bool {
                match self {
                    Self::LigandSymmetry(ls) => ls.present.is_ground(),
                    Self::Fluxionality(f) => f.present.is_ground(),
                    Self::Topicity(t) => t.relation.is_ground(),
                    Self::Stereogenicity(g) => g.is_ground(),
                }
            }

            fn meet(&self, other: &Self) -> Option<Self> {
                match (self, other) {
                    (Self::LigandSymmetry(a), Self::LigandSymmetry(b)) => {
                        a.present.meet(&b.present).map(|present| {
                            Self::LigandSymmetry(LigandSymmetryAst {
                                permutation: a.permutation,
                                present,
                            })
                        })
                    }
                    (Self::Fluxionality(a), Self::Fluxionality(b)) => {
                        a.present.meet(&b.present).map(|present| {
                            Self::Fluxionality(FluxionalityAst {
                                permutation: a.permutation,
                                present,
                            })
                        })
                    }
                    (Self::Topicity(a), Self::Topicity(b)) => {
                        a.relation.meet(&b.relation).map(|relation| {
                            Self::Topicity(TopicityAst {
                                pair: a.pair,
                                relation,
                            })
                        })
                    }
                    (Self::Stereogenicity(a), Self::Stereogenicity(b)) => {
                        a.meet(b).map(Self::Stereogenicity)
                    }
                    _ => None,
                }
            }

            fn join(&self, other: &Self) -> Result<Self, NoJoin> {
                match (self, other) {
                    (Self::LigandSymmetry(a), Self::LigandSymmetry(b)) => {
                        Ok(Self::LigandSymmetry(LigandSymmetryAst {
                            permutation: a.permutation,
                            present: a.present.join(&b.present)?,
                        }))
                    }
                    (Self::Fluxionality(a), Self::Fluxionality(b)) => {
                        Ok(Self::Fluxionality(FluxionalityAst {
                            permutation: a.permutation,
                            present: a.present.join(&b.present)?,
                        }))
                    }
                    (Self::Topicity(a), Self::Topicity(b)) => Ok(Self::Topicity(TopicityAst {
                        pair: a.pair,
                        relation: a.relation.join(&b.relation)?,
                    })),
                    (Self::Stereogenicity(a), Self::Stereogenicity(b)) => {
                        Ok(Self::Stereogenicity(a.join(b)?))
                    }
                    _ => Err(NoJoin),
                }
            }

            fn matches(&self, target: &Self) -> bool {
                match (self, target) {
                    (Self::LigandSymmetry(a), Self::LigandSymmetry(b))
                        if a.permutation == b.permutation =>
                    {
                        a.present.matches(&b.present)
                    }
                    (Self::Fluxionality(a), Self::Fluxionality(b))
                        if a.permutation == b.permutation =>
                    {
                        a.present.matches(&b.present)
                    }
                    (Self::Topicity(a), Self::Topicity(b)) if a.pair == b.pair => {
                        a.relation.matches(&b.relation)
                    }
                    (Self::Stereogenicity(a), Self::Stereogenicity(b)) => a.matches(b),
                    _ => false,
                }
            }

            fn is_compatible(&self, other: &Self) -> bool {
                match (self, other) {
                    (Self::LigandSymmetry(a), Self::LigandSymmetry(b)) => {
                        a.present.is_compatible(&b.present)
                    }
                    (Self::Fluxionality(a), Self::Fluxionality(b)) => {
                        a.present.is_compatible(&b.present)
                    }
                    (Self::Topicity(a), Self::Topicity(b)) => a.relation.is_compatible(&b.relation),
                    (Self::Stereogenicity(a), Self::Stereogenicity(b)) => a.is_compatible(b),
                    _ => false,
                }
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

            /// Ligand-symmetry constraints
            pub fn ligand_symmetries(&self) -> impl Iterator<Item = &LigandSymmetryAst> {
                self.entries.iter().filter_map(|c| match c {
                    $constraint::LigandSymmetry(p) => Some(p),
                    _ => None,
                })
            }

            /// Ligand-symmetry constraint per ligand permutation.
            pub fn ligand_symmetry(
                &self,
                permutation: OrientedLigandPermutation,
            ) -> LigandSymmetryAst {
                self.ligand_symmetries()
                    .find(|ls| ls.permutation == permutation)
                    .map(|ls| ls.clone())
                    .unwrap_or(LigandSymmetryAst {
                        permutation,
                        present: BooleanAst::Undetermined,
                    })
            }

            /// Fluxionality constraints.
            pub fn fluxionalities(&self) -> impl Iterator<Item = &FluxionalityAst> {
                self.entries.iter().filter_map(|c| match c {
                    $constraint::Fluxionality(f) => Some(f),
                    _ => None,
                })
            }

            pub fn fluxionality(&self, permutation: LigandPermutation) -> FluxionalityAst {
                self.fluxionalities()
                    .find(|f| f.permutation == permutation)
                    .map(|f| f.clone())
                    .unwrap_or(FluxionalityAst {
                        permutation,
                        present: BooleanAst::Undetermined,
                    })
            }

            /// Topicity constraints.
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

            /// Stereogenicity constraint.
            pub fn stereogenicity(&self) -> StereogenicityAst {
                match self.get($key::Stereogenicity) {
                    Some($constraint::Stereogenicity(g)) => g.clone(),
                    _ => StereogenicityAst::Undetermined,
                }
            }

            pub fn is_empty(&self) -> bool {
                self.entries.is_empty()
            }

            pub fn len(&self) -> usize {
                self.entries.len()
            }

            fn find(&self, key: $key) -> Result<usize, usize> {
                self.entries.binary_search_by(|c| c.key().cmp(&key))
            }

            pub fn contains(&self, key: $key) -> bool {
                self.find(key).is_ok()
            }

            pub fn get(&self, key: $key) -> Option<&$constraint> {
                self.find(key).ok().map(|i| &self.entries[i])
            }

            /// Insert in sorted order by key, overwrite same key (last-wins).
            pub fn set(&mut self, c: $constraint) {
                match self.find(c.key()) {
                    Ok(i) => self.entries[i] = c,
                    Err(i) => self.entries.insert(i, c),
                }
            }

            /// Transactional write at one key: verify the current value `canonical_eq` `old` (both
            /// absent matches), then apply `new` (`Some` sets, `None` removes). `old`/`new` address
            /// the same key. `Err` on a key or old-value mismatch; the store is unchanged when it
            /// errors. The delta apply/undo primitive.
            pub fn compare_and_set(
                &mut self,
                old: Option<$constraint>,
                new: Option<$constraint>,
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

            pub fn remove(&mut self, key: $key) -> Option<$constraint> {
                self.find(key).ok().map(|i| self.entries.remove(i))
            }

            /// `set` each constraint in turn (last-wins), for bulk construction.
            pub fn extend(&mut self, constraints: impl IntoIterator<Item = $constraint>) {
                for constraint in constraints {
                    self.set(constraint);
                }
            }

            /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
            /// Undetermined entries in `other` remove.
            pub fn update(&mut self, other: &$constraints) {
                for c in other.iter() {
                    if c.is_undetermined() {
                        self.remove(c.key());
                    } else {
                        self.set(c.clone());
                    }
                }
            }

            /// Bulk-remove entries that don't satisfy the predicate.
            pub fn retain(&mut self, mut f: impl FnMut(&$constraint) -> bool) {
                self.entries.retain(|c| f(c));
            }

            /// Remove all entries.
            pub fn clear(&mut self) {
                self.entries.clear();
            }

            /// Move the entries out of the store, leaving it empty.
            pub fn take(&mut self) -> impl Iterator<Item = $constraint> {
                mem::take(&mut self.entries).into_iter()
            }

            pub fn iter(&self) -> Iter<'_, $constraint> {
                self.entries.iter()
            }

            /// No-op: frame-relative ligand positions carry no entity index.
            pub fn compact(self, _compaction: &IdCompaction) -> Self {
                self
            }
        }

        impl Canonicalize for $constraints {
            /// Canonicalize each value and drop the vacuous ones. Keys are already unique and
            /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
            /// canonicalizing a value never changes its `key()`.
            fn canonicalize(self) -> Result<Self, Contradiction> {
                let mut entries = self
                    .entries
                    .into_iter()
                    .map(Canonicalize::canonicalize)
                    .collect::<Result<SmallVec<[$constraint; 2]>, _>>()?;
                entries.retain(|c| !c.is_undetermined());
                Ok(Self { entries })
            }
        }

        impl Lattice for $constraints {
            fn is_undetermined(&self) -> bool {
                self.entries.iter().all(|c| c.is_undetermined())
            }

            fn is_ground(&self) -> bool {
                self.entries.iter().all(|c| c.is_ground())
            }

            /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared
            /// key meets its two values (`$constraint::meet`; a `None` aborts the whole meet), an
            /// A-only / B-only key is kept (meet with the absent ⊤ is the value). Vacuous results
            /// are dropped.
            fn meet(&self, other: &Self) -> Option<Self> {
                let mut entries: SmallVec<[$constraint; 2]> = SmallVec::new();
                let mut a = self.entries.iter();
                let mut b = other.entries.iter();
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
                Some(Self { entries })
            }

            /// Least upper bound as a two-pointer merge: only keys present on *both* sides join
            /// (`$constraint::join`); a single-side key widens to the absent ⊤ and is dropped. The
            /// container always has a top (the empty set), so this is total (`Ok`).
            fn join(&self, other: &Self) -> Result<Self, NoJoin> {
                let mut entries: SmallVec<[$constraint; 2]> = SmallVec::new();
                let mut a = self.entries.iter();
                let mut b = other.entries.iter();
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
                Ok(Self { entries })
            }

            fn matches(&self, target: &Self) -> bool {
                self.ligand_symmetries()
                    .all(|p| target.ligand_symmetries().any(|t| p.matches(t)))
                    && self
                        .fluxionalities()
                        .all(|p| target.fluxionalities().any(|t| p.matches(t)))
                    && self
                        .topicities()
                        .all(|t| t.relation.matches(&target.topicity(t.pair)))
                    && self.stereogenicity().matches(&target.stereogenicity())
            }

            /// Sorted merge, short-circuit: only shared keys can conflict; non-shared keys are
            /// always compatible.
            fn is_compatible(&self, other: &Self) -> bool {
                let mut a = self.entries.iter();
                let mut b = other.entries.iter();
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

        impl FromIterator<$constraint> for $constraints {
            fn from_iter<I: IntoIterator<Item = $constraint>>(iter: I) -> Self {
                let mut out = Self::new();
                for c in iter {
                    out.set(c);
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

stereo_constraint! { StereoAtomConstraintAst, StereoAtomConstraintKey, StereoAtomConstraintsAst }
stereo_constraint! { StereoBondConstraintAst, StereoBondConstraintKey, StereoBondConstraintsAst }

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

            fn join(&self, other: &Self) -> Result<Self, NoJoin> {
                Ok(
                    Self::LitSet(self.to_set().union(&other.to_set()).copied().collect())
                        .canonicalize()
                        .unwrap_or(Self::Undetermined),
                )
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

/// Ligand permutation with a presence assertion: whether the permutation is
/// (`present`) a ligand symmetry. Non-unique.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LigandSymmetryAst {
    pub permutation: OrientedLigandPermutation,
    pub present: BooleanAst,
}

impl LigandSymmetryAst {
    pub fn matches(&self, target: &Self) -> bool {
        self.permutation.matches(&target.permutation) && self.present.matches(&target.present)
    }
}

/// Fluxionality move: proper ligand permutation realized by dynamics, with a
/// presence assertion (whether the move is `present`). Non-unique.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FluxionalityAst {
    pub permutation: LigandPermutation,
    pub present: BooleanAst,
}

impl FluxionalityAst {
    pub fn matches(&self, target: &Self) -> bool {
        self.permutation.matches(&target.permutation) && self.present.matches(&target.present)
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::equal(
        LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
        LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
        true
    )]
    #[case::different(
        LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
        LigandPermutation(Permutation::identity(4)),
        false
    )]
    fn test_ligand_permutation_matches(
        #[case] a: LigandPermutation,
        #[case] b: LigandPermutation,
        #[case] expected: bool,
    ) {
        assert_eq!(a.matches(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::equal(
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper },
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper },
        true)]
    #[case::different_orientation(
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper },
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Improper },
        false)]
    #[case::different_permutation(
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper },
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper },
        false)]
    fn test_oriented_ligand_permutation_matches(
        #[case] a: OrientedLigandPermutation,
        #[case] b: OrientedLigandPermutation,
        #[case] expected: bool,
    ) {
        assert_eq!(a.matches(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ordered(StereoLigandPosition(1), StereoLigandPosition(2), StereoLigandPosition(1), StereoLigandPosition(2))]
    #[case::reversed(StereoLigandPosition(2), StereoLigandPosition(1), StereoLigandPosition(1), StereoLigandPosition(2))]
    #[case::equal(StereoLigandPosition(3), StereoLigandPosition(3), StereoLigandPosition(3), StereoLigandPosition(3))]
    fn test_stereo_ligand_pair_new(
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
        assert_eq!(a.join(&b), Ok(expected));
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

    // `StereogenicityAst` is the second `relation_ast!` instantiation; the macro's
    // full lattice/canonicalize logic is covered by the `TopicityRelationAst` tests
    // above. These confirm the instantiation over the `Stereogenicity` domain.
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

    #[rustfmt::skip]
    #[rstest]
    #[case::same(
        LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) },
        LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) },
        true)]
    #[case::different_presence(
        LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) },
        LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(false) },
        false)]
    #[case::different_permutation(
        LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) },
        LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) },
        false)]
    fn test_ligand_symmetry_ast_matches(
        #[case] pattern: LigandSymmetryAst,
        #[case] target: LigandSymmetryAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same(
        FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), present: BooleanAst::Lit(true) },
        FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), present: BooleanAst::Lit(true) },
        true)]
    #[case::different_permutation(
        FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), present: BooleanAst::Lit(true) },
        FluxionalityAst { permutation: LigandPermutation(Permutation::identity(4)), present: BooleanAst::Lit(true) },
        false)]
    #[case::different_presence(
        FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), present: BooleanAst::Lit(true) },
        FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), present: BooleanAst::Lit(false) },
        false)]
    fn test_fluxionality_ast_matches(
        #[case] pattern: FluxionalityAst,
        #[case] target: FluxionalityAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_pair_open_matches_lit(
        TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Undetermined },
        TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) },
        true)]
    #[case::same_pair_same_lit(
        TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) },
        TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) },
        true)]
    #[case::same_pair_different_lit(
        TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) },
        TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) },
        false)]
    #[case::different_pair(
        TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationAst::Undetermined },
        TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) },
        false)]
    fn test_topicity_ast_matches(
        #[case] pattern: TopicityAst,
        #[case] target: TopicityAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(
        StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
        StereoAtomConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }))]
    #[case::fluxionality(
        StereoAtomConstraintAst::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), present: BooleanAst::Lit(true) }),
        StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3]))))]
    #[case::topicity(
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))))]
    #[case::stereogenicity(
        StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintKey::Stereogenicity)]
    fn test_stereo_atom_constraint_ast_key(
        #[case] c: StereoAtomConstraintAst,
        #[case] expected: StereoAtomConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(
        StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
        StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Undetermined }))]
    #[case::fluxionality(
        StereoAtomConstraintAst::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::identity(4)), present: BooleanAst::Lit(true) }),
        StereoAtomConstraintAst::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::identity(4)), present: BooleanAst::Undetermined }))]
    #[case::topicity(
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Undetermined }))]
    #[case::stereogenicity(
        StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined))]
    fn test_stereo_atom_constraint_ast_as_undetermined(
        #[case] c: StereoAtomConstraintAst,
        #[case] expected: StereoAtomConstraintAst,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity_litset_singleton(
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic])) }),
        Ok(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })))]
    #[case::fluxionality_identity(
        StereoAtomConstraintAst::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::identity(4)), present: BooleanAst::Lit(true) }),
        Ok(StereoAtomConstraintAst::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::identity(4)), present: BooleanAst::Lit(true) })))]
    fn test_stereo_atom_constraint_ast_canonicalize(
        #[case] c: StereoAtomConstraintAst,
        #[case] expected: Result<StereoAtomConstraintAst, Contradiction>,
    ) {
        assert_eq!(c.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry_present(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }), false)]
    #[case::ligand_symmetry_undetermined(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Undetermined }), true)]
    #[case::topicity_lit(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }), false)]
    #[case::topicity_undetermined(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Undetermined }), true)]
    #[case::stereogenicity_lit(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)), false)]
    #[case::stereogenicity_undetermined(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined), true)]
    fn test_stereo_atom_constraint_ast_is_undetermined(#[case] c: StereoAtomConstraintAst, #[case] expected: bool) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry_narrows(
        StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Undetermined }),
        StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
        Some(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) })))]
    #[case::ligand_symmetry_conflict(
        StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
        StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(false) }),
        None)]
    #[case::topicity_disjoint(
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
        None)]
    #[case::different_key(
        StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        None)]
    fn test_stereo_atom_constraint_ast_meet(#[case] a: StereoAtomConstraintAst, #[case] b: StereoAtomConstraintAst, #[case] expected: Option<StereoAtomConstraintAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity_widens(
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
        Ok(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])) })))]
    #[case::different_key(
        StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        Err(NoJoin))]
    fn test_stereo_atom_constraint_ast_join(#[case] a: StereoAtomConstraintAst, #[case] b: StereoAtomConstraintAst, #[case] expected: Result<StereoAtomConstraintAst, NoJoin>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Undetermined }),
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        true)]
    #[case::same_key_incompatible(
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
        false)]
    #[case::different_key(
        StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
        false)]
    fn test_stereo_atom_constraint_ast_is_compatible(#[case] a: StereoAtomConstraintAst, #[case] b: StereoAtomConstraintAst, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_stereo_atom_constraints_ast_new() {
        let cs = StereoAtomConstraintsAst::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.iter().count(), 0);
    }

    #[rstest]
    #[case::present(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(
            StereogenicityAst::Lit(Stereogenicity::Stereogenic)
        )),
        StereogenicityAst::Lit(Stereogenicity::Stereogenic)
    )]
    #[case::absent(StereoAtomConstraintsAst::new(), StereogenicityAst::Undetermined)]
    fn test_stereo_atom_constraints_ast_stereogenicity(
        #[case] cs: StereoAtomConstraintsAst,
        #[case] expected: StereogenicityAst,
    ) {
        assert_eq!(cs.stereogenicity(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity_present(StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))), true)]
    #[case::topicity_absent(StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2))), false)]
    #[case::stereogenicity_present(StereoAtomConstraintKey::Stereogenicity, true)]
    #[case::fluxionality_absent(StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::identity(4))), false)]
    fn test_stereo_atom_constraints_ast_contains(#[case] key: StereoAtomConstraintKey, #[case] expected: bool) {
        let cs = StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ]);
        assert_eq!(cs.contains(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity(
        StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))),
        Some(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })))]
    #[case::stereogenicity(
        StereoAtomConstraintKey::Stereogenicity,
        Some(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))))]
    #[case::absent(
        StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2))),
        None)]
    fn test_stereo_atom_constraints_ast_get(#[case] key: StereoAtomConstraintKey, #[case] expected: Option<StereoAtomConstraintAst>) {
        let cs = StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ]);
        assert_eq!(cs.get(key), expected.as_ref());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))], vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))])]
    #[case::overwrite_unique(
        vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined), StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))],
        vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))])]
    #[case::overwrite_same_ligand_permutation(
        vec![
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(false) }),
        ],
        vec![StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(false) })])]
    #[case::kind_sorted(
        vec![
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
        ],
        vec![
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ])]
    fn test_stereo_atom_constraints_ast_set(
        #[case] sequence: Vec<StereoAtomConstraintAst>,
        #[case] expected: Vec<StereoAtomConstraintAst>,
    ) {
        let mut cs = StereoAtomConstraintsAst::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, StereoAtomConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(
        vec![StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })],
        Some(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        Some(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })),
        Ok(()),
        vec![StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })])]
    #[case::remove(
        vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))],
        Some(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        None,
        Ok(()),
        vec![])]
    #[case::add_from_absent(
        vec![],
        None,
        Some(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        Ok(()),
        vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))])]
    #[case::old_mismatch(
        vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))],
        Some(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Symmetric))),
        None,
        Err(Contradiction),
        vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))])]
    #[case::key_mismatch(
        vec![],
        Some(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        Some(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        Err(Contradiction),
        vec![])]
    fn test_stereo_atom_constraints_ast_compare_and_set(
        #[case] initial: Vec<StereoAtomConstraintAst>,
        #[case] old: Option<StereoAtomConstraintAst>,
        #[case] new: Option<StereoAtomConstraintAst>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<StereoAtomConstraintAst>,
    ) {
        let mut cs = StereoAtomConstraintsAst::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, StereoAtomConstraintsAst::from_iter(expected_state));
    }

    #[rstest]
    fn test_stereo_atom_constraints_ast_remove() {
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1));
        let mut cs = StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::Topicity(TopicityAst {
                pair,
                relation: TopicityRelationAst::Lit(Topicity::Homotopic),
            }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(
                Stereogenicity::Stereogenic,
            )),
        ]);
        assert_eq!(
            cs.remove(StereoAtomConstraintKey::Topicity(pair)),
            Some(StereoAtomConstraintAst::Topicity(TopicityAst {
                pair,
                relation: TopicityRelationAst::Lit(Topicity::Homotopic),
            })),
        );
        assert_eq!(
            cs,
            StereoAtomConstraintsAst::from_iter([StereoAtomConstraintAst::Stereogenicity(
                StereogenicityAst::Lit(Stereogenicity::Stereogenic)
            )]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite_shared(
        vec![StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })],
        vec![StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })],
        vec![StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })])]
    #[case::keeps_disjoint(
        vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))],
        vec![StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })],
        vec![
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ])]
    #[case::vacuous_removes(
        vec![
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ],
        vec![StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Undetermined })],
        vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))])]
    fn test_stereo_atom_constraints_ast_update(
        #[case] initial: Vec<StereoAtomConstraintAst>,
        #[case] other: Vec<StereoAtomConstraintAst>,
        #[case] expected: Vec<StereoAtomConstraintAst>,
    ) {
        let mut cs = StereoAtomConstraintsAst::from_iter(initial);
        cs.update(&StereoAtomConstraintsAst::from_iter(other));
        assert_eq!(cs, StereoAtomConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::drop_vacuous(
        StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Undetermined }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ]),
        Ok(StereoAtomConstraintsAst::from_iter([StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))])))]
    #[case::canonicalizes_values(
        StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic])) }),
        ]),
        Ok(StereoAtomConstraintsAst::from_iter([StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })])))]
    fn test_stereo_atom_constraints_ast_canonicalize(
        #[case] constraints: StereoAtomConstraintsAst,
        #[case] expected: Result<StereoAtomConstraintsAst, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoAtomConstraintsAst::new(), true)]
    #[case::ligand_symmetry_present(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) })),
        false)]
    #[case::topicity_open(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Undetermined })),
        true)]
    #[case::stereogenicity_lit(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        false)]
    fn test_stereo_atom_constraints_ast_is_undetermined(
        #[case] cs: StereoAtomConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(cs.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoAtomConstraintsAst::new(), true)]
    #[case::ligand_symmetry_present(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) })),
        true)]
    #[case::topicity_open(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Undetermined })),
        false)]
    #[case::stereogenicity_lit(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        true)]
    fn test_stereo_atom_constraints_ast_is_ground(
        #[case] cs: StereoAtomConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(cs.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys_kept(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        Some(StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ])))]
    #[case::shared_topicity_value_meet(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])) })),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Homotopic])) })),
        Some(StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }))))]
    #[case::ligand_symmetry_union(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) })),
        StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[0, 1, 3, 2])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
        ]),
        Some(StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[0, 1, 3, 2])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
        ])))]
    #[case::stereogenicity_carried_through(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ]),
        Some(StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ])))]
    #[case::incompatible_same_key_none(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })),
        None)]
    fn test_stereo_atom_constraints_ast_meet(
        #[case] a: StereoAtomConstraintsAst,
        #[case] b: StereoAtomConstraintsAst,
        #[case] expected: Option<StereoAtomConstraintsAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::shared_topicity_widens(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })),
        Ok(StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])) }))))]
    #[case::ligand_symmetry_intersection(
        StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[0, 1, 3, 2])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
        ]),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) })),
        Ok(StereoAtomConstraintsAst::from(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }))))]
    #[case::disjoint_keys_drop_to_empty(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        Ok(StereoAtomConstraintsAst::new()))]
    fn test_stereo_atom_constraints_ast_join(
        #[case] a: StereoAtomConstraintsAst,
        #[case] b: StereoAtomConstraintsAst,
        #[case] expected: Result<StereoAtomConstraintsAst, NoJoin>,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_any(
        StereoAtomConstraintsAst::new(),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) })),
        true)]
    #[case::specific_pattern_absent_in_target(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })),
        StereoAtomConstraintsAst::new(),
        false)]
    #[case::same_ligand_symmetry_and_topicity(
        StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
        ]),
        StereoAtomConstraintsAst::from_iter([
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) }),
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }),
        ]),
        true)]
    #[case::ligand_symmetry_missing_in_target(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) })),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(4, &[0, 1, 3, 2])), orientation: Orientation::Proper }, present: BooleanAst::Lit(true) })),
        false)]
    #[case::topicity_subset(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])) })),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        true)]
    #[case::stereogenicity_mismatch(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Symmetric))),
        false)]
    fn test_stereo_atom_constraints_ast_matches(
        #[case] pattern: StereoAtomConstraintsAst,
        #[case] target: StereoAtomConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        true)]
    #[case::shared_key_compatible(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Undetermined })),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        true)]
    #[case::shared_key_incompatible(
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })),
        false)]
    fn test_stereo_atom_constraints_ast_is_compatible(
        #[case] a: StereoAtomConstraintsAst,
        #[case] b: StereoAtomConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(
        vec![
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ],
        vec![
            StereoAtomConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ])]
    #[case::overwrite_same_key(
        vec![
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ],
        vec![StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))])]
    #[case::empty(vec![], vec![])]
    fn test_stereo_atom_constraints_ast_from_iter(
        #[case] input: Vec<StereoAtomConstraintAst>,
        #[case] expected: Vec<StereoAtomConstraintAst>,
    ) {
        assert_eq!(
            StereoAtomConstraintsAst::from_iter(input),
            StereoAtomConstraintsAst::from_iter(expected),
        );
    }

    #[rstest]
    fn test_stereo_bond_constraints_ast_new() {
        let cs = StereoBondConstraintsAst::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.iter().count(), 0);
    }

    #[rstest]
    fn test_stereo_bond_constraints_ast_set() {
        let mut cs = StereoBondConstraintsAst::new();
        let f = FluxionalityAst {
            permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
            present: BooleanAst::Lit(true),
        };
        cs.set(StereoBondConstraintAst::Fluxionality(f));
        assert_eq!(cs.fluxionalities().copied().collect::<Vec<_>>(), vec![f]);
    }

    // `StereoBondConstraintsAst` is the second `stereo_constraint!` instantiation; the shared macro
    // logic is exercised by the `StereoAtomConstraintsAst` tests above. These confirm the bond
    // instantiation's transactional write and lattice operations independently.
    #[rustfmt::skip]
    #[rstest]
    #[case::modify(
        vec![StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })],
        Some(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        Some(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })),
        Ok(()),
        vec![StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })])]
    #[case::remove(
        vec![StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))],
        Some(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        None,
        Ok(()),
        vec![])]
    #[case::key_mismatch(
        vec![],
        Some(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        Some(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        Err(Contradiction),
        vec![])]
    fn test_stereo_bond_constraints_ast_compare_and_set(
        #[case] initial: Vec<StereoBondConstraintAst>,
        #[case] old: Option<StereoBondConstraintAst>,
        #[case] new: Option<StereoBondConstraintAst>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<StereoBondConstraintAst>,
    ) {
        let mut cs = StereoBondConstraintsAst::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, StereoBondConstraintsAst::from_iter(expected_state));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys_kept(
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        Some(StereoBondConstraintsAst::from_iter([
            StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ])))]
    #[case::shared_topicity_value_meet(
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])) })),
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Homotopic])) })),
        Some(StereoBondConstraintsAst::from(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) }))))]
    #[case::incompatible_same_key_none(
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })),
        None)]
    fn test_stereo_bond_constraints_ast_meet(
        #[case] a: StereoBondConstraintsAst,
        #[case] b: StereoBondConstraintsAst,
        #[case] expected: Option<StereoBondConstraintsAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_any(
        StereoBondConstraintsAst::new(),
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        true)]
    #[case::topicity_subset(
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::NotSet(BTreeSet::from([Topicity::Diastereotopic])) })),
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) })),
        true)]
    #[case::stereogenicity_mismatch(
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        StereoBondConstraintsAst::from(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Symmetric))),
        false)]
    fn test_stereo_bond_constraints_ast_matches(
        #[case] pattern: StereoBondConstraintsAst,
        #[case] target: StereoBondConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(
        vec![
            StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ],
        vec![
            StereoBondConstraintAst::Topicity(TopicityAst { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationAst::Lit(Topicity::Homotopic) }),
            StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        ])]
    #[case::empty(vec![], vec![])]
    fn test_stereo_bond_constraints_ast_from_iter(
        #[case] input: Vec<StereoBondConstraintAst>,
        #[case] expected: Vec<StereoBondConstraintAst>,
    ) {
        assert_eq!(
            StereoBondConstraintsAst::from_iter(input),
            StereoBondConstraintsAst::from_iter(expected),
        );
    }
}
