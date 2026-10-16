//! Stereo-element constraints.

use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::hash::Hash;
use std::mem;
use std::slice::Iter;

use smallvec::SmallVec;
use strum::VariantArray;
use umol_perm::{Orientation, Permutation};

use super::super::boolean::BooleanForm;
use super::super::compact::MoleculeCompaction;
use super::super::error::{Contradiction, NoJoinError};
use super::super::id::StereoLigandPosition;
use super::super::stereo::{Stereogenicity, Topicity};
use super::super::traits::{AsLit, Lattice, Normalize};

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
            LigandSymmetry(LigandSymmetryForm),
            Fluxionality(FluxionalityForm),
            Topicity(TopicityForm),
            Stereogenicity(StereogenicityForm),
        }

        impl $constraint {
            /// Vacuous form at `key`: the top of the key's fiber.
            pub fn undetermined(key: $key) -> Self {
                match key {
                    $key::LigandSymmetry(permutation) => Self::LigandSymmetry(LigandSymmetryForm {
                        permutation,
                        invariant: BooleanForm::Undetermined,
                    }),
                    $key::Fluxionality(permutation) => Self::Fluxionality(FluxionalityForm {
                        permutation,
                        active: BooleanForm::Undetermined,
                    }),
                    $key::Topicity(pair) => Self::Topicity(TopicityForm {
                        pair,
                        relation: TopicityRelationForm::Undetermined,
                    }),
                    $key::Stereogenicity => Self::Stereogenicity(StereogenicityForm::Undetermined),
                }
            }

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
                Self::undetermined(self.key())
            }

            /// Frame-relative ligand positions carry no atom ids, so compact is a no-op.
            pub fn compact(self, _compaction: &MoleculeCompaction) -> Option<Self> {
                Some(self)
            }

            pub(crate) fn uses_participant_frame(&self) -> bool {
                match self {
                    Self::LigandSymmetry(_) | Self::Fluxionality(_) | Self::Topicity(_) => true,
                    Self::Stereogenicity(_) => false,
                }
            }
        }

        impl Normalize for $constraint {
            /// Normalize each value in its own fiber; the key is preserved.
            fn normalize(self) -> Result<Self, Contradiction> {
                Ok(match self {
                    Self::LigandSymmetry(ls) => Self::LigandSymmetry(ls.normalize()?),
                    Self::Fluxionality(f) => Self::Fluxionality(f.normalize()?),
                    Self::Topicity(t) => Self::Topicity(t.normalize()?),
                    Self::Stereogenicity(g) => Self::Stereogenicity(g.normalize()?),
                })
            }
        }

        impl Lattice for $constraint {
            fn is_undetermined(&self) -> bool {
                match self {
                    Self::LigandSymmetry(ls) => ls.is_undetermined(),
                    Self::Fluxionality(f) => f.is_undetermined(),
                    Self::Topicity(t) => t.is_undetermined(),
                    Self::Stereogenicity(g) => g.is_undetermined(),
                }
            }

            fn is_ground(&self) -> bool {
                match self {
                    Self::LigandSymmetry(ls) => ls.is_ground(),
                    Self::Fluxionality(f) => f.is_ground(),
                    Self::Topicity(t) => t.is_ground(),
                    Self::Stereogenicity(g) => g.is_ground(),
                }
            }

            fn meet(&self, other: &Self) -> Option<Self> {
                match (self, other) {
                    (Self::LigandSymmetry(a), Self::LigandSymmetry(b)) => {
                        a.meet(b).map(Self::LigandSymmetry)
                    }
                    (Self::Fluxionality(a), Self::Fluxionality(b)) => {
                        a.meet(b).map(Self::Fluxionality)
                    }
                    (Self::Topicity(a), Self::Topicity(b)) => a.meet(b).map(Self::Topicity),
                    (Self::Stereogenicity(a), Self::Stereogenicity(b)) => {
                        a.meet(b).map(Self::Stereogenicity)
                    }
                    _ => None,
                }
            }

            fn join(&self, other: &Self) -> Result<Self, NoJoinError> {
                match (self, other) {
                    (Self::LigandSymmetry(a), Self::LigandSymmetry(b)) => {
                        a.join(b).map(Self::LigandSymmetry)
                    }
                    (Self::Fluxionality(a), Self::Fluxionality(b)) => {
                        a.join(b).map(Self::Fluxionality)
                    }
                    (Self::Topicity(a), Self::Topicity(b)) => a.join(b).map(Self::Topicity),
                    (Self::Stereogenicity(a), Self::Stereogenicity(b)) => {
                        a.join(b).map(Self::Stereogenicity)
                    }
                    _ => Err(NoJoinError),
                }
            }

            fn matches(&self, target: &Self) -> bool {
                match (self, target) {
                    (Self::LigandSymmetry(a), Self::LigandSymmetry(b)) => a.matches(b),
                    (Self::Fluxionality(a), Self::Fluxionality(b)) => a.matches(b),
                    (Self::Topicity(a), Self::Topicity(b)) => a.matches(b),
                    (Self::Stereogenicity(a), Self::Stereogenicity(b)) => a.matches(b),
                    _ => false,
                }
            }

            fn is_compatible(&self, other: &Self) -> bool {
                match (self, other) {
                    (Self::LigandSymmetry(a), Self::LigandSymmetry(b)) => a.is_compatible(b),
                    (Self::Fluxionality(a), Self::Fluxionality(b)) => a.is_compatible(b),
                    (Self::Topicity(a), Self::Topicity(b)) => a.is_compatible(b),
                    (Self::Stereogenicity(a), Self::Stereogenicity(b)) => a.is_compatible(b),
                    _ => false,
                }
            }
        }

        /// Stereo constraint container for stereo atoms and bonds.
        /// Absent keys contribute top to matching, as stored undetermined values do.
        #[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $constraints {
            entries: SmallVec<[$constraint; 2]>,
        }

        impl $constraints {
            pub fn new() -> Self {
                Self::default()
            }

            /// The normalized own assertion: absence is top; contradiction is `None`.
            /// Reads only this key and leaves storage unchanged.
            pub fn asserted(&self, key: $key) -> Option<$constraint> {
                match self.get(key) {
                    Some(value) => value.clone().normalize().ok(),
                    None => Some($constraint::undetermined(key)),
                }
            }

            /// The normalized own assertion. This entity kind has no closure negative,
            /// so this agrees with [`Self::asserted`]. Reads no molecule context.
            pub fn asserted_complete(&self, key: $key) -> Option<$constraint> {
                self.asserted(key)
            }

            /// Ligand-symmetry constraints
            pub fn ligand_symmetries(&self) -> impl Iterator<Item = &LigandSymmetryForm> {
                self.entries.iter().filter_map(|c| match c {
                    $constraint::LigandSymmetry(p) => Some(p),
                    _ => None,
                })
            }

            /// Ligand-symmetry invariant per ligand permutation.
            pub fn ligand_symmetry(&self, permutation: OrientedLigandPermutation) -> &BooleanForm {
                self.ligand_symmetries()
                    .find(|ls| ls.permutation == permutation)
                    .map_or(&BooleanForm::Undetermined, |ls| &ls.invariant)
            }

            /// Fluxionality constraints.
            pub fn fluxionalities(&self) -> impl Iterator<Item = &FluxionalityForm> {
                self.entries.iter().filter_map(|c| match c {
                    $constraint::Fluxionality(f) => Some(f),
                    _ => None,
                })
            }

            pub fn fluxionality(&self, permutation: LigandPermutation) -> &BooleanForm {
                self.fluxionalities()
                    .find(|f| f.permutation == permutation)
                    .map_or(&BooleanForm::Undetermined, |f| &f.active)
            }

            /// Topicity constraints.
            pub fn topicities(&self) -> impl Iterator<Item = &TopicityForm> {
                self.entries.iter().filter_map(|c| match c {
                    $constraint::Topicity(t) => Some(t),
                    _ => None,
                })
            }

            /// Topicity relation per ligand pair.
            pub fn topicity(&self, pair: StereoLigandPair) -> &TopicityRelationForm {
                self.topicities()
                    .find(|t| t.pair == pair)
                    .map_or(&TopicityRelationForm::Undetermined, |t| &t.relation)
            }

            /// Stereogenicity constraint.
            pub fn stereogenicity(&self) -> &StereogenicityForm {
                match self.get($key::Stereogenicity) {
                    Some($constraint::Stereogenicity(g)) => g,
                    _ => &StereogenicityForm::Undetermined,
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

            /// Transactional write at one key: verify the current value `normalized_eq` `old` (both
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
                    (Some(current), Some(old)) => current.normalized_eq(old),
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
            pub fn take(&mut self) -> impl ExactSizeIterator<Item = $constraint> {
                mem::take(&mut self.entries).into_iter()
            }

            pub fn iter(&self) -> Iter<'_, $constraint> {
                self.entries.iter()
            }

            /// No-op: frame-relative ligand positions carry no entity index.
            pub fn compact(self, _compaction: &MoleculeCompaction) -> Self {
                self
            }
        }

        impl Normalize for $constraints {
            /// Normalize each value and drop the vacuous ones. Keys are already unique and
            /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
            /// canonicalizing a value never changes its `key()`.
            fn normalize(self) -> Result<Self, Contradiction> {
                let mut entries = self
                    .entries
                    .into_iter()
                    .map(Normalize::normalize)
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
            fn join(&self, other: &Self) -> Result<Self, NoJoinError> {
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
                    .all(|p| p.invariant.matches(target.ligand_symmetry(p.permutation)))
                    && self
                        .fluxionalities()
                        .all(|p| p.active.matches(target.fluxionality(p.permutation)))
                    && self
                        .topicities()
                        .all(|t| t.relation.matches(target.topicity(t.pair)))
                    && self.stereogenicity().matches(target.stereogenicity())
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

stereo_constraint! { StereoAtomConstraintForm, StereoAtomConstraintKey, StereoAtomConstraintsForm }
stereo_constraint! { StereoBondConstraintForm, StereoBondConstraintKey, StereoBondConstraintsForm }

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
macro_rules! relation_form {
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

        impl Normalize for $name {
            /// Finite-domain normal form over `to_set`: empty → `Err`, full →
            /// `Undetermined`, singleton → `Lit`, else the smaller of positive /
            /// complement (tiebreak positive).
            fn normalize(self) -> Result<Self, Contradiction> {
                if matches!(self, Self::Undetermined | Self::Lit(_)) {
                    return Ok(self);
                }
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

            /// Structural when both values are trivially normal; otherwise compares normal
            /// forms.
            fn normalized_eq(&self, other: &Self) -> bool {
                match (self, other) {
                    (Self::Undetermined | Self::Lit(_), Self::Undetermined | Self::Lit(_)) => {
                        self == other
                    }
                    _ => self == other || self.clone().normalize() == other.clone().normalize(),
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

            fn is_compatible(&self, other: &Self) -> bool {
                match (self, other) {
                    (Self::Undetermined, Self::Undetermined | Self::Lit(_))
                    | (Self::Lit(_), Self::Undetermined) => true,
                    (Self::Undetermined, Self::LitSet(s))
                    | (Self::LitSet(s), Self::Undetermined) => !s.is_empty(),
                    (Self::Undetermined, Self::NotSet(s))
                    | (Self::NotSet(s), Self::Undetermined) => {
                        s.len() < <$domain as VariantArray>::VARIANTS.len()
                    }
                    (Self::Lit(a), Self::Lit(b)) => a == b,
                    (Self::Lit(value), Self::LitSet(s)) | (Self::LitSet(s), Self::Lit(value)) => {
                        s.contains(value)
                    }
                    (Self::Lit(value), Self::NotSet(s)) | (Self::NotSet(s), Self::Lit(value)) => {
                        !s.contains(value)
                    }
                    (Self::LitSet(a), Self::LitSet(b)) => !a.is_disjoint(b),
                    (Self::LitSet(a), Self::NotSet(b)) | (Self::NotSet(b), Self::LitSet(a)) => {
                        a.iter().any(|value| !b.contains(value))
                    }
                    (Self::NotSet(a), Self::NotSet(b)) => <$domain as VariantArray>::VARIANTS
                        .iter()
                        .any(|value| !a.contains(value) && !b.contains(value)),
                }
            }

            /// Intersection of the admissible sets, folded by `normalize` (∅ → `None`).
            fn meet(&self, other: &Self) -> Option<Self> {
                Self::LitSet(
                    self.to_set()
                        .intersection(&other.to_set())
                        .copied()
                        .collect(),
                )
                .normalize()
                .ok()
            }

            fn join(&self, other: &Self) -> Result<Self, NoJoinError> {
                Ok(
                    Self::LitSet(self.to_set().union(&other.to_set()).copied().collect())
                        .normalize()
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

relation_form! { TopicityRelationForm, Topicity }
relation_form! { StereogenicityForm, Stereogenicity }

/// Ligand permutation with a presence assertion: whether the permutation is
/// (`invariant`) a ligand symmetry. Non-unique.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LigandSymmetryForm {
    pub permutation: OrientedLigandPermutation,
    pub invariant: BooleanForm,
}

impl Normalize for LigandSymmetryForm {
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(Self {
            permutation: self.permutation,
            invariant: self.invariant.normalize()?,
        })
    }
}

/// Meet-semilattice keyed by `permutation`: same permutation delegates to the
/// `invariant` boolean lattice, different permutations lie in different fibers
/// (`meet` → `None`, `join` → `Err(NoJoinError)`).
impl Lattice for LigandSymmetryForm {
    fn is_undetermined(&self) -> bool {
        self.invariant.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.invariant.is_ground()
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        if self.permutation != other.permutation {
            return None;
        }
        self.invariant.meet(&other.invariant).map(|invariant| Self {
            permutation: self.permutation,
            invariant,
        })
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoinError> {
        if self.permutation != other.permutation {
            return Err(NoJoinError);
        }
        Ok(Self {
            permutation: self.permutation,
            invariant: self.invariant.join(&other.invariant)?,
        })
    }

    fn matches(&self, target: &Self) -> bool {
        self.permutation == target.permutation && self.invariant.matches(&target.invariant)
    }

    fn is_compatible(&self, other: &Self) -> bool {
        self.permutation == other.permutation && self.invariant.is_compatible(&other.invariant)
    }
}

/// Fluxionality move: proper ligand permutation realized by dynamics, with a
/// presence assertion (whether the move is `active`). Non-unique.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FluxionalityForm {
    pub permutation: LigandPermutation,
    pub active: BooleanForm,
}

impl Normalize for FluxionalityForm {
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(Self {
            permutation: self.permutation,
            active: self.active.normalize()?,
        })
    }
}

/// Meet-semilattice keyed by `permutation`: same permutation delegates to the
/// `active` boolean lattice, different permutations lie in different fibers
/// (`meet` → `None`, `join` → `Err(NoJoinError)`).
impl Lattice for FluxionalityForm {
    fn is_undetermined(&self) -> bool {
        self.active.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.active.is_ground()
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        if self.permutation != other.permutation {
            return None;
        }
        self.active.meet(&other.active).map(|active| Self {
            permutation: self.permutation,
            active,
        })
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoinError> {
        if self.permutation != other.permutation {
            return Err(NoJoinError);
        }
        Ok(Self {
            permutation: self.permutation,
            active: self.active.join(&other.active)?,
        })
    }

    fn matches(&self, target: &Self) -> bool {
        self.permutation == target.permutation && self.active.matches(&target.active)
    }

    fn is_compatible(&self, other: &Self) -> bool {
        self.permutation == other.permutation && self.active.is_compatible(&other.active)
    }
}

/// Per-pair topicity constraint: relation between pair of ligands.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TopicityForm {
    pub pair: StereoLigandPair,
    pub relation: TopicityRelationForm,
}

impl Normalize for TopicityForm {
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(Self {
            pair: self.pair,
            relation: self.relation.normalize()?,
        })
    }
}

/// Meet-semilattice keyed by `pair`: same pair delegates to the per-pair
/// `relation` lattice, different pairs lie in different fibers (`meet` → `None`,
/// `join` → `Err(NoJoinError)`) — they are incomparable, there is no global top.
impl Lattice for TopicityForm {
    fn is_undetermined(&self) -> bool {
        self.relation.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.relation.is_ground()
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        if self.pair != other.pair {
            return None;
        }
        self.relation.meet(&other.relation).map(|relation| Self {
            pair: self.pair,
            relation,
        })
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoinError> {
        if self.pair != other.pair {
            return Err(NoJoinError);
        }
        Ok(Self {
            pair: self.pair,
            relation: self.relation.join(&other.relation)?,
        })
    }

    fn matches(&self, target: &Self) -> bool {
        self.pair == target.pair && self.relation.matches(&target.relation)
    }

    fn is_compatible(&self, other: &Self) -> bool {
        self.pair == other.pair && self.relation.is_compatible(&other.relation)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
                orientation: Orientation::Proper,
            },
            invariant: BooleanForm::Lit(true),
        }),
        true,
    )]
    #[case::fluxionality(
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
            permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
            active: BooleanForm::Lit(true),
        }),
        true,
    )]
    #[case::topicity(
        StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        }),
        true,
    )]
    #[case::stereogenicity(
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        false,
    )]
    fn test_stereo_atom_constraint_form_uses_participant_frame(
        #[case] constraint: StereoAtomConstraintForm,
        #[case] expected: bool,
    ) {
        assert_eq!(constraint.uses_participant_frame(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(
        StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
                orientation: Orientation::Proper,
            },
            invariant: BooleanForm::Lit(true),
        }),
        true,
    )]
    #[case::fluxionality(
        StereoBondConstraintForm::Fluxionality(FluxionalityForm {
            permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
            active: BooleanForm::Lit(true),
        }),
        true,
    )]
    #[case::topicity(
        StereoBondConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        }),
        true,
    )]
    #[case::stereogenicity(
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        false,
    )]
    fn test_stereo_bond_constraint_form_uses_participant_frame(
        #[case] constraint: StereoBondConstraintForm,
        #[case] expected: bool,
    ) {
        assert_eq!(constraint.uses_participant_frame(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(
        StereoBondConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }),
        StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined }))]
    #[case::fluxionality(
        StereoBondConstraintKey::Fluxionality(LigandPermutation(Permutation::identity(4))),
        StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Undetermined }))]
    #[case::topicity(
        StereoBondConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))),
        StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }))]
    #[case::stereogenicity(StereoBondConstraintKey::Stereogenicity, StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined))]
    fn test_stereo_bond_constraint_form_undetermined(
        #[case] key: StereoBondConstraintKey,
        #[case] expected: StereoBondConstraintForm,
    ) {
        assert_eq!(StereoBondConstraintForm::undetermined(key), expected);
    }

    #[rstest]
    #[case::equal(
        LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
        LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
        true
    )]
    #[case::different(
        LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
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
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper },
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper },
        true)]
    #[case::different_orientation(
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper },
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Improper },
        false)]
    #[case::different_permutation(
        OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper },
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
    #[case::undetermined(TopicityRelationForm::undetermined(), TopicityRelationForm::Undetermined)]
    #[case::lit(TopicityRelationForm::lit(Topicity::Homotopic), TopicityRelationForm::Lit(Topicity::Homotopic))]
    #[case::lit_set(TopicityRelationForm::lit_set([Topicity::Homotopic, Topicity::Enantiotopic]), TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])))]
    #[case::lit_set_singleton_raw(TopicityRelationForm::lit_set([Topicity::Homotopic]), TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic])))]
    #[case::not(TopicityRelationForm::not(Topicity::Homotopic), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic])))]
    #[case::not_set(TopicityRelationForm::not_set([Topicity::Homotopic, Topicity::Enantiotopic]), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])))]
    fn test_topicity_relation_form_constructors(#[case] actual: TopicityRelationForm, #[case] expected: TopicityRelationForm) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::homotopic(Topicity::Homotopic, TopicityRelationForm::Lit(Topicity::Homotopic))]
    #[case::diastereotopic(
        Topicity::Diastereotopic,
        TopicityRelationForm::Lit(Topicity::Diastereotopic)
    )]
    fn test_topicity_relation_form_from(
        #[case] value: Topicity,
        #[case] expected: TopicityRelationForm,
    ) {
        assert_eq!(TopicityRelationForm::from(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::singleton_to_lit(TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic])), Ok(TopicityRelationForm::Lit(Topicity::Homotopic)))]
    #[case::litset_polarity_to_notset(TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])), Ok(TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic]))))]
    #[case::full_litset_to_undetermined(TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic, Topicity::Diastereotopic])), Ok(TopicityRelationForm::Undetermined))]
    #[case::empty_litset_err(TopicityRelationForm::LitSet(BTreeSet::new()), Err(Contradiction))]
    #[case::notset_complement_to_lit(TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])), Ok(TopicityRelationForm::Lit(Topicity::Diastereotopic)))]
    #[case::empty_notset_to_undetermined(TopicityRelationForm::NotSet(BTreeSet::new()), Ok(TopicityRelationForm::Undetermined))]
    #[case::full_notset_err(TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic, Topicity::Diastereotopic])), Err(Contradiction))]
    fn test_topicity_relation_form_normalize(
        #[case] input: TopicityRelationForm,
        #[case] expected: Result<TopicityRelationForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rstest]
    #[case::undetermined(TopicityRelationForm::Undetermined)]
    #[case::lit(TopicityRelationForm::Lit(Topicity::Homotopic))]
    #[case::notset(TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])))]
    fn test_topicity_relation_form_normalize_identity(#[case] input: TopicityRelationForm) {
        assert_eq!(input.clone().normalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    fn test_topicity_relation_form_normalized_eq_default(
        #[values(TopicityRelationForm::Undetermined, TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Lit(Topicity::Diastereotopic), TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic])), TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic, Topicity::Diastereotopic])), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])), TopicityRelationForm::LitSet(BTreeSet::new()))] left: TopicityRelationForm,
        #[values(TopicityRelationForm::Undetermined, TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Lit(Topicity::Diastereotopic), TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic])), TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic, Topicity::Diastereotopic])), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])), TopicityRelationForm::LitSet(BTreeSet::new()))] right: TopicityRelationForm,
    ) {
        assert_eq!(left.normalized_eq(&right), left.clone().normalize() == right.clone().normalize());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(TopicityRelationForm::Lit(Topicity::Homotopic), Some(Topicity::Homotopic))]
    #[case::undetermined(TopicityRelationForm::Undetermined, None)]
    #[case::notset(TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic])), None)]
    #[case::litset(TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])), None)]
    fn test_topicity_relation_form_as_lit(#[case] r: TopicityRelationForm, #[case] expected: Option<Topicity>) {
        assert_eq!(r.as_lit(), expected);
        assert_eq!(r.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(TopicityRelationForm::Undetermined, true)]
    #[case::lit(TopicityRelationForm::Lit(Topicity::Homotopic), false)]
    #[case::notset(TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])), false)]
    fn test_topicity_relation_form_is_undetermined(#[case] r: TopicityRelationForm, #[case] expected: bool) {
        assert_eq!(r.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(TopicityRelationForm::Undetermined, TopicityRelationForm::Lit(Topicity::Homotopic), Some(TopicityRelationForm::Lit(Topicity::Homotopic)))]
    #[case::lit_eq(TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Lit(Topicity::Homotopic), Some(TopicityRelationForm::Lit(Topicity::Homotopic)))]
    #[case::lit_disjoint(TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Lit(Topicity::Enantiotopic), None)]
    #[case::lit_notset_in(TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])), Some(TopicityRelationForm::Lit(Topicity::Homotopic)))]
    #[case::lit_notset_out(TopicityRelationForm::Lit(Topicity::Diastereotopic), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])), None)]
    #[case::notset_notset_to_lit(TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Enantiotopic])), Some(TopicityRelationForm::Lit(Topicity::Homotopic)))]
    fn test_topicity_relation_form_meet(
        #[case] a: TopicityRelationForm,
        #[case] b: TopicityRelationForm,
        #[case] expected: Option<TopicityRelationForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und(TopicityRelationForm::Undetermined, TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Undetermined)]
    #[case::lit_eq(TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Lit(Topicity::Homotopic))]
    #[case::lit_union_to_notset(TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Lit(Topicity::Enantiotopic), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])))]
    #[case::lit_notset_to_full(TopicityRelationForm::Lit(Topicity::Diastereotopic), TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])), TopicityRelationForm::Undetermined)]
    fn test_topicity_relation_form_join(
        #[case] a: TopicityRelationForm,
        #[case] b: TopicityRelationForm,
        #[case] expected: TopicityRelationForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(TopicityRelationForm::Undetermined, TopicityRelationForm::Lit(Topicity::Homotopic), true)]
    #[case::lit_und(TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Undetermined, false)]
    #[case::lit_eq(TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Lit(Topicity::Homotopic), true)]
    #[case::lit_neq(TopicityRelationForm::Lit(Topicity::Homotopic), TopicityRelationForm::Lit(Topicity::Enantiotopic), false)]
    #[case::notset_in(TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])), TopicityRelationForm::Lit(Topicity::Homotopic), true)]
    #[case::notset_out(TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])), TopicityRelationForm::Lit(Topicity::Diastereotopic), false)]
    fn test_topicity_relation_form_matches(
        #[case] pattern: TopicityRelationForm,
        #[case] target: TopicityRelationForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    // `StereogenicityForm` is the second `relation_form!` instantiation; the macro's
    // full lattice/normalize logic is covered by the `TopicityRelationForm` tests
    // above. These confirm the instantiation over the `Stereogenicity` domain.
    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereogenicityForm::Lit(Stereogenicity::Stereogenic), Some(Stereogenicity::Stereogenic))]
    #[case::undetermined(StereogenicityForm::Undetermined, None)]
    #[case::notset(StereogenicityForm::NotSet(BTreeSet::from([Stereogenicity::Stereogenic])), None)]
    fn test_stereogenicity_form_as_lit(#[case] g: StereogenicityForm, #[case] expected: Option<Stereogenicity>) {
        assert_eq!(g.as_lit(), expected);
        assert_eq!(g.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::singleton_to_lit(StereogenicityForm::LitSet(BTreeSet::from([Stereogenicity::Stereogenic])), Ok(StereogenicityForm::Lit(Stereogenicity::Stereogenic)))]
    #[case::full_to_undetermined(StereogenicityForm::LitSet(BTreeSet::from([Stereogenicity::Symmetric, Stereogenicity::Prochiral, Stereogenicity::Stereogenic])), Ok(StereogenicityForm::Undetermined))]
    #[case::empty_err(StereogenicityForm::LitSet(BTreeSet::new()), Err(Contradiction))]
    fn test_stereogenicity_form_normalize(
        #[case] input: StereogenicityForm,
        #[case] expected: Result<StereogenicityForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    fn test_stereogenicity_form_normalized_eq_default(
        #[values(StereogenicityForm::Undetermined, StereogenicityForm::Lit(Stereogenicity::Stereogenic), StereogenicityForm::Lit(Stereogenicity::Symmetric), StereogenicityForm::LitSet(BTreeSet::from([Stereogenicity::Stereogenic])), StereogenicityForm::LitSet(BTreeSet::from([Stereogenicity::Symmetric, Stereogenicity::Prochiral, Stereogenicity::Stereogenic])), StereogenicityForm::NotSet(BTreeSet::from([Stereogenicity::Symmetric])), StereogenicityForm::LitSet(BTreeSet::new()))] left: StereogenicityForm,
        #[values(StereogenicityForm::Undetermined, StereogenicityForm::Lit(Stereogenicity::Stereogenic), StereogenicityForm::Lit(Stereogenicity::Symmetric), StereogenicityForm::LitSet(BTreeSet::from([Stereogenicity::Stereogenic])), StereogenicityForm::LitSet(BTreeSet::from([Stereogenicity::Symmetric, Stereogenicity::Prochiral, Stereogenicity::Stereogenic])), StereogenicityForm::NotSet(BTreeSet::from([Stereogenicity::Symmetric])), StereogenicityForm::LitSet(BTreeSet::new()))] right: StereogenicityForm,
    ) {
        assert_eq!(left.normalized_eq(&right), left.clone().normalize() == right.clone().normalize());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same(
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        true)]
    #[case::different_presence(
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) },
        false)]
    #[case::different_permutation(
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        false)]
    fn test_ligand_symmetry_form_matches(
        #[case] pattern: LigandSymmetryForm,
        #[case] target: LigandSymmetryForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same(
        FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) },
        FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) },
        true)]
    #[case::different_permutation(
        FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) },
        FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(true) },
        false)]
    #[case::different_presence(
        FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) },
        FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(false) },
        false)]
    fn test_fluxionality_form_matches(
        #[case] pattern: FluxionalityForm,
        #[case] target: FluxionalityForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_pair_open_matches_lit(
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined },
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        true)]
    #[case::same_pair_same_lit(
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        true)]
    #[case::same_pair_different_lit(
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) },
        false)]
    #[case::different_pair(
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationForm::Undetermined },
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        false)]
    fn test_topicity_form_matches(
        #[case] pattern: TopicityForm,
        #[case] target: TopicityForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_permutation_narrows(
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined },
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        Some(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }))]
    #[case::same_permutation_conflict(
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) },
        None)]
    #[case::different_permutation(
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        None)]
    fn test_ligand_symmetry_form_meet(
        #[case] a: LigandSymmetryForm,
        #[case] b: LigandSymmetryForm,
        #[case] expected: Option<LigandSymmetryForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_permutation(
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        Ok(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }))]
    #[case::different_permutation(
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) },
        Err(NoJoinError))]
    fn test_ligand_symmetry_form_join(
        #[case] a: LigandSymmetryForm,
        #[case] b: LigandSymmetryForm,
        #[case] expected: Result<LigandSymmetryForm, NoJoinError>,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_permutation_narrows(
        FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined },
        FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) },
        Some(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) }))]
    #[case::different_permutation(
        FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) },
        FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(true) },
        None)]
    fn test_fluxionality_form_meet(
        #[case] a: FluxionalityForm,
        #[case] b: FluxionalityForm,
        #[case] expected: Option<FluxionalityForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_permutation(
        FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(true) },
        FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(true) },
        Ok(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(true) }))]
    #[case::different_permutation(
        FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) },
        FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(true) },
        Err(NoJoinError))]
    fn test_fluxionality_form_join(
        #[case] a: FluxionalityForm,
        #[case] b: FluxionalityForm,
        #[case] expected: Result<FluxionalityForm, NoJoinError>,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_pair_narrows(
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined },
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        Some(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }))]
    #[case::same_pair_conflict(
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) },
        None)]
    #[case::different_pair(
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        None)]
    fn test_topicity_form_meet(
        #[case] a: TopicityForm,
        #[case] b: TopicityForm,
        #[case] expected: Option<TopicityForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_pair(
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        Ok(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }))]
    #[case::different_pair(
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) },
        Err(NoJoinError))]
    fn test_topicity_form_join(
        #[case] a: TopicityForm,
        #[case] b: TopicityForm,
        #[case] expected: Result<TopicityForm, NoJoinError>,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(
        StereoAtomConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }),
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined }))]
    #[case::fluxionality(
        StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::identity(4))),
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Undetermined }))]
    #[case::topicity(
        StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }))]
    #[case::stereogenicity(StereoAtomConstraintKey::Stereogenicity, StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined))]
    fn test_stereo_atom_constraint_form_undetermined(
        #[case] key: StereoAtomConstraintKey,
        #[case] expected: StereoAtomConstraintForm,
    ) {
        assert_eq!(StereoAtomConstraintForm::undetermined(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        StereoAtomConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }))]
    #[case::fluxionality(
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) }),
        StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::from_image(&[1, 0, 2, 3]))))]
    #[case::topicity(
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))))]
    #[case::stereogenicity(
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintKey::Stereogenicity)]
    fn test_stereo_atom_constraint_form_key(
        #[case] c: StereoAtomConstraintForm,
        #[case] expected: StereoAtomConstraintKey,
    ) {
        assert_eq!(c.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry(
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined }))]
    #[case::fluxionality(
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(true) }),
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Undetermined }))]
    #[case::topicity(
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }))]
    #[case::stereogenicity(
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined))]
    fn test_stereo_atom_constraint_form_as_undetermined(
        #[case] c: StereoAtomConstraintForm,
        #[case] expected: StereoAtomConstraintForm,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity_litset_singleton(
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic])) }),
        Ok(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })))]
    #[case::fluxionality_identity(
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(true) }),
        Ok(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(true) })))]
    fn test_stereo_atom_constraint_form_normalize(
        #[case] c: StereoAtomConstraintForm,
        #[case] expected: Result<StereoAtomConstraintForm, Contradiction>,
    ) {
        assert_eq!(c.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry_present(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }), false)]
    #[case::ligand_symmetry_undetermined(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined }), true)]
    #[case::topicity_lit(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }), false)]
    #[case::topicity_undetermined(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }), true)]
    #[case::stereogenicity_lit(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)), false)]
    #[case::stereogenicity_undetermined(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined), true)]
    fn test_stereo_atom_constraint_form_is_undetermined(#[case] c: StereoAtomConstraintForm, #[case] expected: bool) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ligand_symmetry_narrows(
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined }),
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        Some(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })))]
    #[case::ligand_symmetry_conflict(
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) }),
        None)]
    #[case::topicity_disjoint(
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }),
        None)]
    #[case::different_key(
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        None)]
    fn test_stereo_atom_constraint_form_meet(#[case] a: StereoAtomConstraintForm, #[case] b: StereoAtomConstraintForm, #[case] expected: Option<StereoAtomConstraintForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity_widens(
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }),
        Ok(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])) })))]
    #[case::different_key(
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        Err(NoJoinError))]
    fn test_stereo_atom_constraint_form_join(#[case] a: StereoAtomConstraintForm, #[case] b: StereoAtomConstraintForm, #[case] expected: Result<StereoAtomConstraintForm, NoJoinError>) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_compatible(
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        true)]
    #[case::same_key_incompatible(
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }),
        false)]
    #[case::different_key(
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        false)]
    fn test_stereo_atom_constraint_form_is_compatible(#[case] a: StereoAtomConstraintForm, #[case] b: StereoAtomConstraintForm, #[case] expected: bool) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rstest]
    fn test_stereo_atom_constraints_form_new() {
        let cs = StereoAtomConstraintsForm::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.iter().count(), 0);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::absent(None, StereoAtomConstraintKey::Stereogenicity, Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)))]
    #[case::top(Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), StereoAtomConstraintKey::Stereogenicity, Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)))]
    #[case::literal(Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), StereoAtomConstraintKey::Stereogenicity, Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))))]
    #[case::singleton(Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic]))), StereoAtomConstraintKey::Stereogenicity, Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))))]
    #[case::bottom(Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([]))), StereoAtomConstraintKey::Stereogenicity, None, None)]
    #[case::topicity_absent(None, StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))), Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })), Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })))]
    #[case::topicity_top(Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })), StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))), Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })), Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })))]
    #[case::ligand_symmetry_absent(None, StereoAtomConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }), Some(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })), Some(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })))]
    #[case::ligand_symmetry_top(Some(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })), StereoAtomConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }), Some(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })), Some(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })))]
    #[case::fluxionality_absent(None, StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::from_image(&[1, 0, 2, 3]))), Some(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })), Some(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })))]
    #[case::fluxionality_top(Some(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })), StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::from_image(&[1, 0, 2, 3]))), Some(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })), Some(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })))]
    fn test_stereo_atom_constraints_form_asserted(
        #[case] source: Option<StereoAtomConstraintForm>,
        #[case] key: StereoAtomConstraintKey,
        #[case] expected: Option<StereoAtomConstraintForm>,
        #[case] complete: Option<StereoAtomConstraintForm>,
    ) {
        let constraints: StereoAtomConstraintsForm = source.into_iter().collect();
        let before = constraints.clone();
        assert_eq!(constraints.asserted(key), expected);
        assert_eq!(constraints.asserted_complete(key), complete);
        assert_eq!(constraints, before);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, BooleanForm::Lit(true))]
    #[case::absent(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Improper }, BooleanForm::Undetermined)]
    fn test_stereo_atom_constraints_form_ligand_symmetry(#[case] permutation: OrientedLigandPermutation, #[case] expected: BooleanForm) {
        let cs = StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper },
            invariant: BooleanForm::Lit(true),
        }));
        assert_eq!(cs.ligand_symmetry(permutation), &expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored(LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), BooleanForm::Lit(false))]
    #[case::absent(LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), BooleanForm::Undetermined)]
    fn test_stereo_atom_constraints_form_fluxionality(#[case] permutation: LigandPermutation, #[case] expected: BooleanForm) {
        let cs = StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
            permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
            active: BooleanForm::Lit(false),
        }));
        assert_eq!(cs.fluxionality(permutation), &expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), TopicityRelationForm::Lit(Topicity::Homotopic))]
    #[case::absent(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), TopicityRelationForm::Undetermined)]
    fn test_stereo_atom_constraints_form_topicity(#[case] pair: StereoLigandPair, #[case] expected: TopicityRelationForm) {
        let cs = StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        }));
        assert_eq!(cs.topicity(pair), &expected);
    }

    #[rstest]
    #[case::present(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(Stereogenicity::Stereogenic)
        )),
        StereogenicityForm::Lit(Stereogenicity::Stereogenic)
    )]
    #[case::absent(StereoAtomConstraintsForm::new(), StereogenicityForm::Undetermined)]
    fn test_stereo_atom_constraints_form_stereogenicity(
        #[case] cs: StereoAtomConstraintsForm,
        #[case] expected: StereogenicityForm,
    ) {
        assert_eq!(cs.stereogenicity(), &expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity_present(StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))), true)]
    #[case::topicity_absent(StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2))), false)]
    #[case::stereogenicity_present(StereoAtomConstraintKey::Stereogenicity, true)]
    #[case::fluxionality_absent(StereoAtomConstraintKey::Fluxionality(LigandPermutation(Permutation::identity(4))), false)]
    fn test_stereo_atom_constraints_form_contains(#[case] key: StereoAtomConstraintKey, #[case] expected: bool) {
        let cs = StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ]);
        assert_eq!(cs.contains(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity(
        StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))),
        Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })))]
    #[case::stereogenicity(
        StereoAtomConstraintKey::Stereogenicity,
        Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))))]
    #[case::absent(
        StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2))),
        None)]
    fn test_stereo_atom_constraints_form_get(#[case] key: StereoAtomConstraintKey, #[case] expected: Option<StereoAtomConstraintForm>) {
        let cs = StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ]);
        assert_eq!(cs.get(key), expected.as_ref());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))], vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))])]
    #[case::overwrite_unique(
        vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined), StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))],
        vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))])]
    #[case::overwrite_same_ligand_permutation(
        vec![
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) }),
        ],
        vec![StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) })])]
    #[case::kind_sorted(
        vec![
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        ],
        vec![
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ])]
    fn test_stereo_atom_constraints_form_set(
        #[case] sequence: Vec<StereoAtomConstraintForm>,
        #[case] expected: Vec<StereoAtomConstraintForm>,
    ) {
        let mut cs = StereoAtomConstraintsForm::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, StereoAtomConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(
        vec![StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })],
        Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })),
        Ok(()),
        vec![StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })])]
    #[case::remove(
        vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))],
        Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        None,
        Ok(()),
        vec![])]
    #[case::add_from_absent(
        vec![],
        None,
        Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        Ok(()),
        vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))])]
    #[case::old_mismatch(
        vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))],
        Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric))),
        None,
        Err(Contradiction),
        vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))])]
    #[case::key_mismatch(
        vec![],
        Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        Err(Contradiction),
        vec![])]
    fn test_stereo_atom_constraints_form_compare_and_set(
        #[case] initial: Vec<StereoAtomConstraintForm>,
        #[case] old: Option<StereoAtomConstraintForm>,
        #[case] new: Option<StereoAtomConstraintForm>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<StereoAtomConstraintForm>,
    ) {
        let mut cs = StereoAtomConstraintsForm::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, StereoAtomConstraintsForm::from_iter(expected_state));
    }

    #[rstest]
    fn test_stereo_atom_constraints_form_remove() {
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1));
        let mut cs = StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::Topicity(TopicityForm {
                pair,
                relation: TopicityRelationForm::Lit(Topicity::Homotopic),
            }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            )),
        ]);
        assert_eq!(
            cs.remove(StereoAtomConstraintKey::Topicity(pair)),
            Some(StereoAtomConstraintForm::Topicity(TopicityForm {
                pair,
                relation: TopicityRelationForm::Lit(Topicity::Homotopic),
            })),
        );
        assert_eq!(
            cs,
            StereoAtomConstraintsForm::from_iter([StereoAtomConstraintForm::Stereogenicity(
                StereogenicityForm::Lit(Stereogenicity::Stereogenic)
            )]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite_shared(
        vec![StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })],
        vec![StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })],
        vec![StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })])]
    #[case::keeps_disjoint(
        vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))],
        vec![StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })],
        vec![
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ])]
    #[case::vacuous_removes(
        vec![
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ],
        vec![StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })],
        vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))])]
    fn test_stereo_atom_constraints_form_update(
        #[case] initial: Vec<StereoAtomConstraintForm>,
        #[case] other: Vec<StereoAtomConstraintForm>,
        #[case] expected: Vec<StereoAtomConstraintForm>,
    ) {
        let mut cs = StereoAtomConstraintsForm::from_iter(initial);
        cs.update(&StereoAtomConstraintsForm::from_iter(other));
        assert_eq!(cs, StereoAtomConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_stereo_atom_constraints_form_take() {
        let mut empty = StereoAtomConstraintsForm::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let constraint = StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        ));
        let mut constraints = StereoAtomConstraintsForm::from(constraint.clone());
        let mut taken = constraints.take();
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(taken.next(), Some(constraint));
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.size_hint(), (0, Some(0)));
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(constraints, StereoAtomConstraintsForm::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::drop_vacuous(
        StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ]),
        Ok(StereoAtomConstraintsForm::from_iter([StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))])))]
    #[case::normalizes_values(
        StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic])) }),
        ]),
        Ok(StereoAtomConstraintsForm::from_iter([StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })])))]
    fn test_stereo_atom_constraints_form_normalize(
        #[case] constraints: StereoAtomConstraintsForm,
        #[case] expected: Result<StereoAtomConstraintsForm, Contradiction>,
    ) {
        assert_eq!(constraints.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoAtomConstraintsForm::new(), true)]
    #[case::ligand_symmetry_present(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        false)]
    #[case::topicity_open(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })),
        true)]
    #[case::stereogenicity_lit(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        false)]
    fn test_stereo_atom_constraints_form_is_undetermined(
        #[case] cs: StereoAtomConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(cs.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoAtomConstraintsForm::new(), true)]
    #[case::ligand_symmetry_present(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        true)]
    #[case::topicity_open(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })),
        false)]
    #[case::stereogenicity_lit(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        true)]
    fn test_stereo_atom_constraints_form_is_ground(
        #[case] cs: StereoAtomConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(cs.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys_kept(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        Some(StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ])))]
    #[case::shared_topicity_value_meet(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic])) })),
        Some(StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }))))]
    #[case::ligand_symmetry_union(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        ]),
        Some(StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        ])))]
    #[case::stereogenicity_carried_through(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ]),
        Some(StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ])))]
    #[case::incompatible_same_key_none(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })),
        None)]
    fn test_stereo_atom_constraints_form_meet(
        #[case] a: StereoAtomConstraintsForm,
        #[case] b: StereoAtomConstraintsForm,
        #[case] expected: Option<StereoAtomConstraintsForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::shared_topicity_widens(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })),
        Ok(StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])) }))))]
    #[case::ligand_symmetry_intersection(
        StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        ]),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        Ok(StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }))))]
    #[case::disjoint_keys_drop_to_empty(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        Ok(StereoAtomConstraintsForm::new()))]
    fn test_stereo_atom_constraints_form_join(
        #[case] a: StereoAtomConstraintsForm,
        #[case] b: StereoAtomConstraintsForm,
        #[case] expected: Result<StereoAtomConstraintsForm, NoJoinError>,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_any(
        StereoAtomConstraintsForm::new(),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        true)]
    #[case::specific_pattern_absent_in_target(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })),
        StereoAtomConstraintsForm::new(),
        false)]
    #[case::same_ligand_symmetry_and_topicity(
        StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }),
        ]),
        StereoAtomConstraintsForm::from_iter([
            StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }),
        ]),
        true)]
    #[case::ligand_symmetry_missing_in_target(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        false)]
    #[case::topicity_subset(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        true)]
    #[case::stereogenicity_mismatch(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric))),
        false)]
    #[case::ligand_symmetry_vacuous_absent(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })),
        StereoAtomConstraintsForm::new(),
        true)]
    #[case::ligand_symmetry_determined_absent(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::new(),
        false)]
    #[case::ligand_symmetry_determined_vacuous(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })),
        false)]
    #[case::ligand_symmetry_determined_equal(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        true)]
    #[case::ligand_symmetry_determined_conflicting(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) })),
        false)]
    #[case::ligand_symmetry_vacuous_other_permutation(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        true)]
    #[case::ligand_symmetry_determined_other_permutation(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        false)]
    #[case::ligand_symmetry_vacuous_other_orientation(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Improper }, invariant: BooleanForm::Lit(true) })),
        true)]
    #[case::ligand_symmetry_determined_other_orientation(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Improper }, invariant: BooleanForm::Lit(true) })),
        false)]
    #[case::fluxionality_vacuous_absent(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })),
        StereoAtomConstraintsForm::new(),
        true)]
    #[case::fluxionality_determined_absent(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::new(),
        false)]
    #[case::fluxionality_determined_vacuous(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })),
        false)]
    #[case::fluxionality_determined_equal(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        true)]
    #[case::fluxionality_determined_conflicting(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(false) })),
        false)]
    #[case::fluxionality_vacuous_other_permutation(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), active: BooleanForm::Lit(true) })),
        true)]
    #[case::fluxionality_determined_other_permutation(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), active: BooleanForm::Lit(true) })),
        false)]
    #[case::topicity_vacuous_absent(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })),
        StereoAtomConstraintsForm::new(),
        true)]
    #[case::topicity_determined_absent(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsForm::new(),
        false)]
    #[case::topicity_determined_vacuous(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })),
        false)]
    #[case::topicity_determined_equal(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        true)]
    #[case::topicity_determined_conflicting(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })),
        false)]
    #[case::topicity_vacuous_other_pair(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        true)]
    #[case::topicity_determined_other_pair(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        false)]
    #[case::stereogenicity_vacuous_absent(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
        StereoAtomConstraintsForm::new(),
        true)]
    #[case::stereogenicity_determined_absent(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsForm::new(),
        false)]
    #[case::stereogenicity_determined_vacuous(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
        false)]
    #[case::stereogenicity_determined_equal(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        true)]
    #[case::stereogenicity_determined_conflicting(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric))),
        false)]
    fn test_stereo_atom_constraints_form_matches(
        #[case] pattern: StereoAtomConstraintsForm,
        #[case] target: StereoAtomConstraintsForm,
        #[case] expected: bool,
    ) {
        let by_meet = match (pattern.meet(&target), target.clone().normalize()) {
            (Some(meet), Ok(target)) => meet == target,
            _ => false,
        };
        assert_eq!(by_meet, expected);
        assert_eq!(pattern.matches(&target), expected);
        assert_eq!(
            pattern.normalize().unwrap().matches(&target.normalize().unwrap()),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        true)]
    #[case::shared_key_compatible(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        true)]
    #[case::shared_key_incompatible(
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })),
        false)]
    fn test_stereo_atom_constraints_form_is_compatible(
        #[case] a: StereoAtomConstraintsForm,
        #[case] b: StereoAtomConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(a.is_compatible(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(
        vec![
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ],
        vec![
            StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ])]
    #[case::overwrite_same_key(
        vec![
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ],
        vec![StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))])]
    #[case::empty(vec![], vec![])]
    fn test_stereo_atom_constraints_form_from_iter(
        #[case] input: Vec<StereoAtomConstraintForm>,
        #[case] expected: Vec<StereoAtomConstraintForm>,
    ) {
        assert_eq!(
            StereoAtomConstraintsForm::from_iter(input),
            StereoAtomConstraintsForm::from_iter(expected),
        );
    }

    #[rstest]
    fn test_stereo_bond_constraints_form_new() {
        let cs = StereoBondConstraintsForm::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.iter().count(), 0);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::absent(None, StereoBondConstraintKey::Stereogenicity, Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)))]
    #[case::top(Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), StereoBondConstraintKey::Stereogenicity, Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)))]
    #[case::literal(Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), StereoBondConstraintKey::Stereogenicity, Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))))]
    #[case::singleton(Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic]))), StereoBondConstraintKey::Stereogenicity, Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))))]
    #[case::bottom(Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([]))), StereoBondConstraintKey::Stereogenicity, None, None)]
    #[case::topicity_absent(None, StereoBondConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))), Some(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })), Some(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })))]
    #[case::topicity_top(Some(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })), StereoBondConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))), Some(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })), Some(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })))]
    #[case::ligand_symmetry_absent(None, StereoBondConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }), Some(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })), Some(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })))]
    #[case::ligand_symmetry_top(Some(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })), StereoBondConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }), Some(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })), Some(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })))]
    #[case::fluxionality_absent(None, StereoBondConstraintKey::Fluxionality(LigandPermutation(Permutation::from_image(&[1, 0, 2, 3]))), Some(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })), Some(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })))]
    #[case::fluxionality_top(Some(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })), StereoBondConstraintKey::Fluxionality(LigandPermutation(Permutation::from_image(&[1, 0, 2, 3]))), Some(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })), Some(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })))]
    fn test_stereo_bond_constraints_form_asserted(
        #[case] source: Option<StereoBondConstraintForm>,
        #[case] key: StereoBondConstraintKey,
        #[case] expected: Option<StereoBondConstraintForm>,
        #[case] complete: Option<StereoBondConstraintForm>,
    ) {
        let constraints: StereoBondConstraintsForm = source.into_iter().collect();
        let before = constraints.clone();
        assert_eq!(constraints.asserted(key), expected);
        assert_eq!(constraints.asserted_complete(key), complete);
        assert_eq!(constraints, before);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, BooleanForm::Lit(true))]
    #[case::absent(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Improper }, BooleanForm::Undetermined)]
    fn test_stereo_bond_constraints_form_ligand_symmetry(#[case] permutation: OrientedLigandPermutation, #[case] expected: BooleanForm) {
        let cs = StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper },
            invariant: BooleanForm::Lit(true),
        }));
        assert_eq!(cs.ligand_symmetry(permutation), &expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored(LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), BooleanForm::Lit(false))]
    #[case::absent(LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), BooleanForm::Undetermined)]
    fn test_stereo_bond_constraints_form_fluxionality(#[case] permutation: LigandPermutation, #[case] expected: BooleanForm) {
        let cs = StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm {
            permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
            active: BooleanForm::Lit(false),
        }));
        assert_eq!(cs.fluxionality(permutation), &expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), TopicityRelationForm::Lit(Topicity::Homotopic))]
    #[case::absent(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), TopicityRelationForm::Undetermined)]
    fn test_stereo_bond_constraints_form_topicity(#[case] pair: StereoLigandPair, #[case] expected: TopicityRelationForm) {
        let cs = StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        }));
        assert_eq!(cs.topicity(pair), &expected);
    }

    #[rstest]
    #[case::stored(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(Stereogenicity::Stereogenic)
        )),
        StereogenicityForm::Lit(Stereogenicity::Stereogenic)
    )]
    #[case::absent(StereoBondConstraintsForm::new(), StereogenicityForm::Undetermined)]
    fn test_stereo_bond_constraints_form_stereogenicity(
        #[case] cs: StereoBondConstraintsForm,
        #[case] expected: StereogenicityForm,
    ) {
        assert_eq!(cs.stereogenicity(), &expected);
    }

    #[rstest]
    fn test_stereo_bond_constraints_form_set() {
        let mut cs = StereoBondConstraintsForm::new();
        let f = FluxionalityForm {
            permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
            active: BooleanForm::Lit(true),
        };
        cs.set(StereoBondConstraintForm::Fluxionality(f));
        assert_eq!(cs.fluxionalities().copied().collect::<Vec<_>>(), vec![f]);
    }

    // `StereoBondConstraintsForm` is the second `stereo_constraint!` instantiation; the shared macro
    // logic is exercised by the `StereoAtomConstraintsForm` tests above. These confirm the bond
    // instantiation's transactional write and lattice operations independently.
    #[rustfmt::skip]
    #[rstest]
    #[case::modify(
        vec![StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })],
        Some(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        Some(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })),
        Ok(()),
        vec![StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })])]
    #[case::remove(
        vec![StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))],
        Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        None,
        Ok(()),
        vec![])]
    #[case::key_mismatch(
        vec![],
        Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        Some(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        Err(Contradiction),
        vec![])]
    fn test_stereo_bond_constraints_form_compare_and_set(
        #[case] initial: Vec<StereoBondConstraintForm>,
        #[case] old: Option<StereoBondConstraintForm>,
        #[case] new: Option<StereoBondConstraintForm>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<StereoBondConstraintForm>,
    ) {
        let mut cs = StereoBondConstraintsForm::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, StereoBondConstraintsForm::from_iter(expected_state));
    }

    #[rstest]
    fn test_stereo_bond_constraints_form_take() {
        let mut empty = StereoBondConstraintsForm::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let constraint = StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
            Stereogenicity::Stereogenic,
        ));
        let mut constraints = StereoBondConstraintsForm::from(constraint.clone());
        let mut taken = constraints.take();
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(taken.next(), Some(constraint));
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.size_hint(), (0, Some(0)));
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(constraints, StereoBondConstraintsForm::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::disjoint_keys_kept(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        Some(StereoBondConstraintsForm::from_iter([
            StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ])))]
    #[case::shared_topicity_value_meet(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic])) })),
        Some(StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }))))]
    #[case::incompatible_same_key_none(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })),
        None)]
    fn test_stereo_bond_constraints_form_meet(
        #[case] a: StereoBondConstraintsForm,
        #[case] b: StereoBondConstraintsForm,
        #[case] expected: Option<StereoBondConstraintsForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_pattern_matches_any(
        StereoBondConstraintsForm::new(),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        true)]
    #[case::topicity_subset(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        true)]
    #[case::stereogenicity_mismatch(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric))),
        false)]
    #[case::ligand_symmetry_vacuous_absent(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })),
        StereoBondConstraintsForm::new(),
        true)]
    #[case::ligand_symmetry_determined_absent(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::new(),
        false)]
    #[case::ligand_symmetry_determined_vacuous(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })),
        false)]
    #[case::ligand_symmetry_determined_equal(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        true)]
    #[case::ligand_symmetry_determined_conflicting(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) })),
        false)]
    #[case::ligand_symmetry_vacuous_other_permutation(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        true)]
    #[case::ligand_symmetry_determined_other_permutation(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        false)]
    #[case::ligand_symmetry_vacuous_other_orientation(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Improper }, invariant: BooleanForm::Lit(true) })),
        true)]
    #[case::ligand_symmetry_determined_other_orientation(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Improper }, invariant: BooleanForm::Lit(true) })),
        false)]
    #[case::fluxionality_vacuous_absent(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })),
        StereoBondConstraintsForm::new(),
        true)]
    #[case::fluxionality_determined_absent(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::new(),
        false)]
    #[case::fluxionality_determined_vacuous(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })),
        false)]
    #[case::fluxionality_determined_equal(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        true)]
    #[case::fluxionality_determined_conflicting(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(false) })),
        false)]
    #[case::fluxionality_vacuous_other_permutation(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Undetermined })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), active: BooleanForm::Lit(true) })),
        true)]
    #[case::fluxionality_determined_other_permutation(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(true) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])), active: BooleanForm::Lit(true) })),
        false)]
    #[case::topicity_vacuous_absent(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })),
        StereoBondConstraintsForm::new(),
        true)]
    #[case::topicity_determined_absent(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoBondConstraintsForm::new(),
        false)]
    #[case::topicity_determined_vacuous(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })),
        false)]
    #[case::topicity_determined_equal(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        true)]
    #[case::topicity_determined_conflicting(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })),
        false)]
    #[case::topicity_vacuous_other_pair(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        true)]
    #[case::topicity_determined_other_pair(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
        false)]
    #[case::stereogenicity_vacuous_absent(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
        StereoBondConstraintsForm::new(),
        true)]
    #[case::stereogenicity_determined_absent(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoBondConstraintsForm::new(),
        false)]
    #[case::stereogenicity_determined_vacuous(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
        false)]
    #[case::stereogenicity_determined_equal(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        true)]
    #[case::stereogenicity_determined_conflicting(
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric))),
        false)]
    fn test_stereo_bond_constraints_form_matches(
        #[case] pattern: StereoBondConstraintsForm,
        #[case] target: StereoBondConstraintsForm,
        #[case] expected: bool,
    ) {
        let by_meet = match (pattern.meet(&target), target.clone().normalize()) {
            (Some(meet), Ok(target)) => meet == target,
            _ => false,
        };
        assert_eq!(by_meet, expected);
        assert_eq!(pattern.matches(&target), expected);
        assert_eq!(
            pattern.normalize().unwrap().matches(&target.normalize().unwrap()),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(
        vec![
            StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ],
        vec![
            StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        ])]
    #[case::empty(vec![], vec![])]
    fn test_stereo_bond_constraints_form_from_iter(
        #[case] input: Vec<StereoBondConstraintForm>,
        #[case] expected: Vec<StereoBondConstraintForm>,
    ) {
        assert_eq!(
            StereoBondConstraintsForm::from_iter(input),
            StereoBondConstraintsForm::from_iter(expected),
        );
    }
}
