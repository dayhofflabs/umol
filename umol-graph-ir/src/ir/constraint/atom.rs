//! Atom constraints.

use std::borrow::Cow;
use std::cmp::Ordering;
use std::mem;

use smallvec::SmallVec;

use super::super::compact::MoleculeCompaction;
use super::super::constraint::ring::{RingMembershipForm, RingScope};
use super::super::error::{Contradiction, NoJoin};
use super::super::num::NumForm;
use super::super::remap::IdRemapping;
use super::super::stereo::TetrahedralStereoForm;
use super::super::traits::{AsLit, Lattice, Normalize};

/// Atom-scope constraint.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AtomConstraintForm {
    Valence(NumForm),
    DonatedPairs(NumForm),
    AcceptedPairs(NumForm),
    AromaticValence(AromaticValenceForm),
    MulticenterValence(MulticenterValenceForm),
    TetrahedralStereo(TetrahedralStereoForm),
    Degree(NumForm),
    TotalDegree(NumForm),
    TotalValence(NumForm),
    /// Incident bonds belonging to the fixed Relevant ring projection.
    RingDegree(NumForm),
    /// Sum of incident bond orders in the fixed Relevant ring projection.
    RingValence(NumForm),
    TotalHydrogens(NumForm),
    /// Ring count in the fixed Relevant ring projection, optionally restricted by size.
    RingMembership(RingMembershipForm),
}

impl AtomConstraintForm {
    pub fn valence(v: impl Into<NumForm>) -> Self {
        Self::Valence(v.into())
    }

    pub fn donated_pairs(v: impl Into<NumForm>) -> Self {
        Self::DonatedPairs(v.into())
    }

    pub fn accepted_pairs(v: impl Into<NumForm>) -> Self {
        Self::AcceptedPairs(v.into())
    }

    pub fn aromatic_valence(v: impl Into<AromaticValenceForm>) -> Self {
        Self::AromaticValence(v.into())
    }

    pub fn multicenter_valence(v: impl Into<MulticenterValenceForm>) -> Self {
        Self::MulticenterValence(v.into())
    }

    pub fn tetrahedral_stereo(c: impl Into<TetrahedralStereoForm>) -> Self {
        Self::TetrahedralStereo(c.into())
    }

    pub fn degree(v: impl Into<NumForm>) -> Self {
        Self::Degree(v.into())
    }

    pub fn total_degree(v: impl Into<NumForm>) -> Self {
        Self::TotalDegree(v.into())
    }

    pub fn total_valence(v: impl Into<NumForm>) -> Self {
        Self::TotalValence(v.into())
    }

    pub fn ring_degree(v: impl Into<NumForm>) -> Self {
        Self::RingDegree(v.into())
    }

    pub fn ring_valence(v: impl Into<NumForm>) -> Self {
        Self::RingValence(v.into())
    }

    pub fn total_hydrogens(v: impl Into<NumForm>) -> Self {
        Self::TotalHydrogens(v.into())
    }

    pub fn ring_membership(scope: RingScope, count: impl Into<NumForm>) -> Self {
        Self::RingMembership(RingMembershipForm::new(scope, count))
    }

    /// Atom constraint key, unique within an `AtomConstraintsForm` container.
    pub fn key(&self) -> AtomConstraintKey {
        match self {
            Self::Valence(_) => AtomConstraintKey::Valence,
            Self::DonatedPairs(_) => AtomConstraintKey::DonatedPairs,
            Self::AcceptedPairs(_) => AtomConstraintKey::AcceptedPairs,
            Self::AromaticValence(_) => AtomConstraintKey::AromaticValence,
            Self::MulticenterValence(_) => AtomConstraintKey::MulticenterValence,
            Self::TetrahedralStereo(_) => AtomConstraintKey::TetrahedralStereo,
            Self::Degree(_) => AtomConstraintKey::Degree,
            Self::TotalDegree(_) => AtomConstraintKey::TotalDegree,
            Self::TotalValence(_) => AtomConstraintKey::TotalValence,
            Self::RingDegree(_) => AtomConstraintKey::RingDegree,
            Self::RingValence(_) => AtomConstraintKey::RingValence,
            Self::TotalHydrogens(_) => AtomConstraintKey::TotalHydrogens,
            Self::RingMembership(m) => AtomConstraintKey::RingMembership(m.scope),
        }
    }

    /// Vacuous form of constraint key, used for removal.
    pub fn as_undetermined(&self) -> Self {
        match self {
            Self::Valence(_) => Self::Valence(NumForm::Undetermined),
            Self::DonatedPairs(_) => Self::DonatedPairs(NumForm::Undetermined),
            Self::AcceptedPairs(_) => Self::AcceptedPairs(NumForm::Undetermined),
            Self::AromaticValence(_) => Self::AromaticValence(AromaticValenceForm::Undetermined),
            Self::MulticenterValence(_) => {
                Self::MulticenterValence(MulticenterValenceForm::Undetermined)
            }
            Self::TetrahedralStereo(_) => {
                Self::TetrahedralStereo(TetrahedralStereoForm::Undetermined)
            }
            Self::Degree(_) => Self::Degree(NumForm::Undetermined),
            Self::TotalDegree(_) => Self::TotalDegree(NumForm::Undetermined),
            Self::TotalValence(_) => Self::TotalValence(NumForm::Undetermined),
            Self::RingDegree(_) => Self::RingDegree(NumForm::Undetermined),
            Self::RingValence(_) => Self::RingValence(NumForm::Undetermined),
            Self::TotalHydrogens(_) => Self::TotalHydrogens(NumForm::Undetermined),
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipForm::new(m.scope, NumForm::Undetermined))
            }
        }
    }

    /// Value-only payload: no entity ids to compact, so this never drops.
    pub fn compact(self, _compaction: &MoleculeCompaction) -> Option<Self> {
        Some(self)
    }

    /// Value-only payload: no entity ids to remap.
    pub fn remap(self, _map: &IdRemapping) -> Self {
        self
    }
}

impl Normalize for AtomConstraintForm {
    /// Normalize the inner value; kind and sub-key are preserved.
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Valence(v) => Self::Valence(v.normalize()?),
            Self::DonatedPairs(v) => Self::DonatedPairs(v.normalize()?),
            Self::AcceptedPairs(v) => Self::AcceptedPairs(v.normalize()?),
            Self::AromaticValence(c) => Self::AromaticValence(c.normalize()?),
            Self::MulticenterValence(c) => Self::MulticenterValence(c.normalize()?),
            Self::TetrahedralStereo(c) => Self::TetrahedralStereo(c.normalize()?),
            Self::Degree(v) => Self::Degree(v.normalize()?),
            Self::TotalDegree(v) => Self::TotalDegree(v.normalize()?),
            Self::TotalValence(v) => Self::TotalValence(v.normalize()?),
            Self::RingDegree(v) => Self::RingDegree(v.normalize()?),
            Self::RingValence(v) => Self::RingValence(v.normalize()?),
            Self::TotalHydrogens(v) => Self::TotalHydrogens(v.normalize()?),
            Self::RingMembership(m) => Self::RingMembership(m.normalize()?),
        })
    }
}

impl Lattice for AtomConstraintForm {
    fn is_undetermined(&self) -> bool {
        match self {
            Self::Valence(v)
            | Self::DonatedPairs(v)
            | Self::AcceptedPairs(v)
            | Self::Degree(v)
            | Self::TotalDegree(v)
            | Self::TotalValence(v)
            | Self::RingDegree(v)
            | Self::RingValence(v)
            | Self::TotalHydrogens(v) => v.is_undetermined(),
            Self::AromaticValence(c) => c.is_undetermined(),
            Self::MulticenterValence(c) => c.is_undetermined(),
            Self::TetrahedralStereo(c) => c.is_undetermined(),
            Self::RingMembership(m) => m.is_undetermined(),
        }
    }

    fn is_ground(&self) -> bool {
        match self {
            Self::Valence(v)
            | Self::DonatedPairs(v)
            | Self::AcceptedPairs(v)
            | Self::Degree(v)
            | Self::TotalDegree(v)
            | Self::TotalValence(v)
            | Self::RingDegree(v)
            | Self::RingValence(v)
            | Self::TotalHydrogens(v) => v.is_ground(),
            Self::AromaticValence(c) => c.is_ground(),
            Self::MulticenterValence(c) => c.is_ground(),
            Self::TetrahedralStereo(c) => c.is_ground(),
            Self::RingMembership(m) => m.is_ground(),
        }
    }

    /// Meet per-key, None on an incompatible meet.
    fn meet(&self, other: &AtomConstraintForm) -> Option<AtomConstraintForm> {
        match (self, other) {
            (Self::Valence(a), Self::Valence(b)) => a.meet(b).map(Self::Valence),
            (Self::DonatedPairs(a), Self::DonatedPairs(b)) => a.meet(b).map(Self::DonatedPairs),
            (Self::AcceptedPairs(a), Self::AcceptedPairs(b)) => a.meet(b).map(Self::AcceptedPairs),
            (Self::AromaticValence(a), Self::AromaticValence(b)) => {
                a.meet(b).map(Self::AromaticValence)
            }
            (Self::MulticenterValence(a), Self::MulticenterValence(b)) => {
                a.meet(b).map(Self::MulticenterValence)
            }
            (Self::TetrahedralStereo(a), Self::TetrahedralStereo(b)) => {
                a.meet(b).map(Self::TetrahedralStereo)
            }
            (Self::Degree(a), Self::Degree(b)) => a.meet(b).map(Self::Degree),
            (Self::TotalDegree(a), Self::TotalDegree(b)) => a.meet(b).map(Self::TotalDegree),
            (Self::TotalValence(a), Self::TotalValence(b)) => a.meet(b).map(Self::TotalValence),
            (Self::RingDegree(a), Self::RingDegree(b)) => a.meet(b).map(Self::RingDegree),
            (Self::RingValence(a), Self::RingValence(b)) => a.meet(b).map(Self::RingValence),
            (Self::TotalHydrogens(a), Self::TotalHydrogens(b)) => {
                a.meet(b).map(Self::TotalHydrogens)
            }
            (Self::RingMembership(a), Self::RingMembership(b)) => {
                a.meet(b).map(Self::RingMembership)
            }
            _ => None,
        }
    }

    /// Join per-key, Err(NoJoin) when operands are in different fibers (different keys).
    fn join(&self, other: &AtomConstraintForm) -> Result<AtomConstraintForm, NoJoin> {
        match (self, other) {
            (Self::Valence(a), Self::Valence(b)) => Ok(Self::Valence(a.join(b)?)),
            (Self::DonatedPairs(a), Self::DonatedPairs(b)) => Ok(Self::DonatedPairs(a.join(b)?)),
            (Self::AcceptedPairs(a), Self::AcceptedPairs(b)) => Ok(Self::AcceptedPairs(a.join(b)?)),
            (Self::AromaticValence(a), Self::AromaticValence(b)) => {
                Ok(Self::AromaticValence(a.join(b)?))
            }
            (Self::MulticenterValence(a), Self::MulticenterValence(b)) => {
                Ok(Self::MulticenterValence(a.join(b)?))
            }
            (Self::TetrahedralStereo(a), Self::TetrahedralStereo(b)) => {
                Ok(Self::TetrahedralStereo(a.join(b)?))
            }
            (Self::Degree(a), Self::Degree(b)) => Ok(Self::Degree(a.join(b)?)),
            (Self::TotalDegree(a), Self::TotalDegree(b)) => Ok(Self::TotalDegree(a.join(b)?)),
            (Self::TotalValence(a), Self::TotalValence(b)) => Ok(Self::TotalValence(a.join(b)?)),
            (Self::RingDegree(a), Self::RingDegree(b)) => Ok(Self::RingDegree(a.join(b)?)),
            (Self::RingValence(a), Self::RingValence(b)) => Ok(Self::RingValence(a.join(b)?)),
            (Self::TotalHydrogens(a), Self::TotalHydrogens(b)) => {
                Ok(Self::TotalHydrogens(a.join(b)?))
            }
            (Self::RingMembership(a), Self::RingMembership(b)) => {
                a.join(b).map(Self::RingMembership)
            }
            _ => Err(NoJoin),
        }
    }

    /// Pattern-driven per-variant match: same key with a matching payload; a mismatched sub-key
    /// (ring scope) or a different kind never matches. Overrides the `meet`-derived default.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Valence(a), Self::Valence(b)) => a.matches(b),
            (Self::DonatedPairs(a), Self::DonatedPairs(b)) => a.matches(b),
            (Self::AcceptedPairs(a), Self::AcceptedPairs(b)) => a.matches(b),
            (Self::AromaticValence(a), Self::AromaticValence(b)) => a.matches(b),
            (Self::MulticenterValence(a), Self::MulticenterValence(b)) => a.matches(b),
            (Self::TetrahedralStereo(a), Self::TetrahedralStereo(b)) => a.matches(b),
            (Self::Degree(a), Self::Degree(b)) => a.matches(b),
            (Self::TotalDegree(a), Self::TotalDegree(b)) => a.matches(b),
            (Self::TotalValence(a), Self::TotalValence(b)) => a.matches(b),
            (Self::RingDegree(a), Self::RingDegree(b)) => a.matches(b),
            (Self::RingValence(a), Self::RingValence(b)) => a.matches(b),
            (Self::TotalHydrogens(a), Self::TotalHydrogens(b)) => a.matches(b),
            (Self::RingMembership(a), Self::RingMembership(b)) => a.matches(b),
            _ => false,
        }
    }

    /// Compatible iff same key with compatible payloads; different keys are incompatible (one
    /// constraint can't be two kinds). Overrides the `meet`-derived default to skip building the
    /// `AtomConstraintForm`.
    fn is_compatible(&self, other: &AtomConstraintForm) -> bool {
        match (self, other) {
            (Self::Valence(a), Self::Valence(b)) => a.is_compatible(b),
            (Self::DonatedPairs(a), Self::DonatedPairs(b)) => a.is_compatible(b),
            (Self::AcceptedPairs(a), Self::AcceptedPairs(b)) => a.is_compatible(b),
            (Self::AromaticValence(a), Self::AromaticValence(b)) => a.is_compatible(b),
            (Self::MulticenterValence(a), Self::MulticenterValence(b)) => a.is_compatible(b),
            (Self::TetrahedralStereo(a), Self::TetrahedralStereo(b)) => a.is_compatible(b),
            (Self::Degree(a), Self::Degree(b)) => a.is_compatible(b),
            (Self::TotalDegree(a), Self::TotalDegree(b)) => a.is_compatible(b),
            (Self::TotalValence(a), Self::TotalValence(b)) => a.is_compatible(b),
            (Self::RingDegree(a), Self::RingDegree(b)) => a.is_compatible(b),
            (Self::RingValence(a), Self::RingValence(b)) => a.is_compatible(b),
            (Self::TotalHydrogens(a), Self::TotalHydrogens(b)) => a.is_compatible(b),
            (Self::RingMembership(a), Self::RingMembership(b)) => a.is_compatible(b),
            _ => false,
        }
    }
}

/// Entry identity: discriminant + sub-key, AtomConstraintsForm is ordered, unique by key.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AtomConstraintKey {
    Valence,
    DonatedPairs,
    AcceptedPairs,
    AromaticValence,
    MulticenterValence,
    TetrahedralStereo,
    Degree,
    TotalDegree,
    TotalValence,
    RingDegree,
    RingValence,
    TotalHydrogens,
    RingMembership(RingScope),
}

/// Atom constraint container, unique by key, sorted flat vector storage.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomConstraintsForm {
    entries: SmallVec<[AtomConstraintForm; 2]>,
}

impl AtomConstraintsForm {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn valence(&self) -> Option<&NumForm> {
        match self.get(AtomConstraintKey::Valence) {
            Some(AtomConstraintForm::Valence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn aromatic_valence(&self) -> Option<&AromaticValenceForm> {
        match self.get(AtomConstraintKey::AromaticValence) {
            Some(AtomConstraintForm::AromaticValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn multicenter_valence(&self) -> Option<&MulticenterValenceForm> {
        match self.get(AtomConstraintKey::MulticenterValence) {
            Some(AtomConstraintForm::MulticenterValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn tetrahedral_stereo(&self) -> Option<&TetrahedralStereoForm> {
        match self.get(AtomConstraintKey::TetrahedralStereo) {
            Some(AtomConstraintForm::TetrahedralStereo(c)) => Some(c),
            _ => None,
        }
    }

    pub fn degree(&self) -> Option<&NumForm> {
        match self.get(AtomConstraintKey::Degree) {
            Some(AtomConstraintForm::Degree(v)) => Some(v),
            _ => None,
        }
    }

    pub fn total_degree(&self) -> Option<&NumForm> {
        match self.get(AtomConstraintKey::TotalDegree) {
            Some(AtomConstraintForm::TotalDegree(v)) => Some(v),
            _ => None,
        }
    }

    pub fn total_valence(&self) -> Option<&NumForm> {
        match self.get(AtomConstraintKey::TotalValence) {
            Some(AtomConstraintForm::TotalValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn ring_degree(&self) -> Option<&NumForm> {
        match self.get(AtomConstraintKey::RingDegree) {
            Some(AtomConstraintForm::RingDegree(v)) => Some(v),
            _ => None,
        }
    }

    pub fn ring_valence(&self) -> Option<&NumForm> {
        match self.get(AtomConstraintKey::RingValence) {
            Some(AtomConstraintForm::RingValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn total_hydrogens(&self) -> Option<&NumForm> {
        match self.get(AtomConstraintKey::TotalHydrogens) {
            Some(AtomConstraintForm::TotalHydrogens(v)) => Some(v),
            _ => None,
        }
    }

    pub fn donated_pairs(&self) -> Option<&NumForm> {
        match self.get(AtomConstraintKey::DonatedPairs) {
            Some(AtomConstraintForm::DonatedPairs(v)) => Some(v),
            _ => None,
        }
    }

    pub fn accepted_pairs(&self) -> Option<&NumForm> {
        match self.get(AtomConstraintKey::AcceptedPairs) {
            Some(AtomConstraintForm::AcceptedPairs(v)) => Some(v),
            _ => None,
        }
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &NumForm)> {
        self.iter().filter_map(|c| match c {
            AtomConstraintForm::RingMembership(m) => Some((m.scope, &m.count)),
            _ => None,
        })
    }

    fn ring_membership(&self, scope: RingScope) -> Option<&NumForm> {
        self.ring_memberships()
            .find(|(s, _)| *s == scope)
            .map(|(_, v)| v)
    }

    pub fn ring_count(&self) -> Option<&NumForm> {
        self.ring_membership(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> Option<&NumForm> {
        self.ring_membership(RingScope::Size(s))
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    fn find(&self, key: AtomConstraintKey) -> Result<usize, usize> {
        self.entries.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains(&self, key: AtomConstraintKey) -> bool {
        self.find(key).is_ok()
    }

    pub fn get(&self, key: AtomConstraintKey) -> Option<&AtomConstraintForm> {
        self.find(key).ok().map(|i| &self.entries[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: AtomConstraintForm) {
        match self.find(c.key()) {
            Ok(i) => self.entries[i] = c,
            Err(i) => self.entries.insert(i, c),
        }
    }

    /// Transactional write at one key: verify the current value equals `old` (by `normalized_eq`;
    /// both absent is a match), then apply `new` (`Some` sets, `None` removes).
    /// `Err` on a key or old-value mismatch; the store is unchanged when it errors.
    pub fn compare_and_set(
        &mut self,
        old: Option<AtomConstraintForm>,
        new: Option<AtomConstraintForm>,
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

    pub fn remove(&mut self, key: AtomConstraintKey) -> Option<AtomConstraintForm> {
        self.find(key).ok().map(|i| self.entries.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = AtomConstraintForm>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &AtomConstraintsForm) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&AtomConstraintForm) -> bool) {
        self.entries.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl ExactSizeIterator<Item = AtomConstraintForm> {
        mem::take(&mut self.entries).into_iter()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = &AtomConstraintForm> {
        self.entries.iter()
    }

    /// No-op: no `AtomConstraintForm` variant carries an entity index.
    pub fn compact(self, _compaction: &MoleculeCompaction) -> Self {
        self
    }
}

impl Normalize for AtomConstraintsForm {
    /// Normalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn normalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .entries
            .into_iter()
            .map(Normalize::normalize)
            .collect::<Result<SmallVec<[AtomConstraintForm; 2]>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self { entries })
    }
}

impl Lattice for AtomConstraintsForm {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`AtomConstraintForm::meet`; a `None` aborts the whole meet), an A-only /
    /// B-only key is kept (meet with the absent ⊤ is the value). Vacuous results are dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: SmallVec<[AtomConstraintForm; 2]> = SmallVec::new();
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
    /// (`AtomConstraintForm::join`); a single-side key widens to the absent ⊤ and is dropped. A
    /// same-key join never returns `Err(NoJoin)`, but if it did the key would simply drop
    /// (widen to ⊤). The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: SmallVec<[AtomConstraintForm; 2]> = SmallVec::new();
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

    /// Pattern-driven: every constraint the pattern carries must match the target's
    /// corresponding value, looked up by reference (an absent target value is
    /// `Undetermined`, which matches). Cost is proportional to the pattern's
    /// constraints, not the field count; an empty pattern matches any target.
    fn matches(&self, target: &Self) -> bool {
        self.iter().all(|c| match c {
            AtomConstraintForm::Valence(v) => {
                v.matches(target.valence().unwrap_or(&NumForm::Undetermined))
            }
            AtomConstraintForm::DonatedPairs(v) => {
                v.matches(target.donated_pairs().unwrap_or(&NumForm::Undetermined))
            }
            AtomConstraintForm::AcceptedPairs(v) => {
                v.matches(target.accepted_pairs().unwrap_or(&NumForm::Undetermined))
            }
            AtomConstraintForm::AromaticValence(av) => av.matches(
                target
                    .aromatic_valence()
                    .unwrap_or(&AromaticValenceForm::Undetermined),
            ),
            AtomConstraintForm::MulticenterValence(mv) => mv.matches(
                target
                    .multicenter_valence()
                    .unwrap_or(&MulticenterValenceForm::Undetermined),
            ),
            AtomConstraintForm::TetrahedralStereo(ts) => ts.matches(
                target
                    .tetrahedral_stereo()
                    .unwrap_or(&TetrahedralStereoForm::Undetermined),
            ),
            AtomConstraintForm::Degree(v) => {
                v.matches(target.degree().unwrap_or(&NumForm::Undetermined))
            }
            AtomConstraintForm::TotalDegree(v) => {
                v.matches(target.total_degree().unwrap_or(&NumForm::Undetermined))
            }
            AtomConstraintForm::TotalValence(v) => {
                v.matches(target.total_valence().unwrap_or(&NumForm::Undetermined))
            }
            AtomConstraintForm::RingDegree(v) => {
                v.matches(target.ring_degree().unwrap_or(&NumForm::Undetermined))
            }
            AtomConstraintForm::RingValence(v) => {
                v.matches(target.ring_valence().unwrap_or(&NumForm::Undetermined))
            }
            AtomConstraintForm::TotalHydrogens(v) => {
                v.matches(target.total_hydrogens().unwrap_or(&NumForm::Undetermined))
            }
            AtomConstraintForm::RingMembership(rm) => rm.count.matches(
                target
                    .ring_membership(rm.scope)
                    .unwrap_or(&NumForm::Undetermined),
            ),
        })
    }

    /// Sorted merge, short-circuit: only shared keys can conflict; non-shared keys are always
    /// compatible. Cheaper than the `meet`-derived default — builds no result container.
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

impl FromIterator<AtomConstraintForm> for AtomConstraintsForm {
    fn from_iter<I: IntoIterator<Item = AtomConstraintForm>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for AtomConstraintsForm {
    type Item = AtomConstraintForm;
    type IntoIter = smallvec::IntoIter<[AtomConstraintForm; 2]>;

    fn into_iter(self) -> Self::IntoIter {
        self.entries.into_iter()
    }
}

impl From<AtomConstraintForm> for AtomConstraintsForm {
    fn from(c: AtomConstraintForm) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<AtomConstraintForm>> for AtomConstraintsForm {
    fn from(cs: Vec<AtomConstraintForm>) -> Self {
        Self::from_iter(cs)
    }
}

/// Aromatic-valence state of an atom: `Undetermined`, explicitly
/// `NotAromatic`, or participating in an aromatic system with the given
/// aromatic-valence count.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticValence {
    NotAromatic,
    Aromatic(i64),
}

impl AromaticValence {
    /// Aromatic-valence count for calculations that identify absence with zero.
    pub const fn valence_count(self) -> i64 {
        match self {
            Self::NotAromatic => 0,
            Self::Aromatic(valence) => valence,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticValenceForm {
    #[default]
    Undetermined,
    NotAromatic,
    Aromatic(NumForm),
}

impl AromaticValenceForm {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn not_aromatic() -> Self {
        Self::NotAromatic
    }

    pub fn aromatic(v: impl Into<NumForm>) -> Self {
        Self::Aromatic(v.into())
    }

    pub fn is_aromatic(&self) -> bool {
        matches!(self, Self::Aromatic(_))
    }

    pub fn aromatic_covalence(&self) -> NumForm {
        self.as_lit()
            .map(AromaticValence::valence_count)
            .map(aromatic_covalence)
            .map(NumForm::Lit)
            .unwrap_or(NumForm::Undetermined)
    }

    /// Pattern matches value.
    pub fn matches_value(&self, value: i64) -> bool {
        match self {
            Self::Aromatic(v) => v.matches(&NumForm::Lit(value)),
            Self::NotAromatic => value == 0,
            Self::Undetermined => true,
        }
    }
}

impl From<AromaticValence> for AromaticValenceForm {
    fn from(valence: AromaticValence) -> Self {
        match valence {
            AromaticValence::NotAromatic => Self::NotAromatic,
            AromaticValence::Aromatic(valence) => Self::Aromatic(NumForm::Lit(valence)),
        }
    }
}

impl Normalize for AromaticValenceForm {
    /// Delegate to the inner `NumForm`; `NotAromatic`/`Undetermined` identity.
    /// No cross-variant fold (`Aromatic(Lit(0))` stays distinct from `NotAromatic`).
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Aromatic(v) => Self::Aromatic(v.normalize()?),
            other => other,
        })
    }

    fn normalized(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            Self::Aromatic(v) => Ok(match v.normalized()? {
                Cow::Borrowed(_) => Cow::Borrowed(self),
                Cow::Owned(cv) => Cow::Owned(Self::Aromatic(cv)),
            }),
            _ => Ok(Cow::Borrowed(self)),
        }
    }
}

impl Lattice for AromaticValenceForm {
    #[inline]
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    fn is_ground(&self) -> bool {
        match self {
            Self::Undetermined => false,
            Self::NotAromatic => true,
            Self::Aromatic(v) => v.is_ground(),
        }
    }

    /// Greatest lower bound, canonicalizing operands and output.
    fn meet(&self, other: &Self) -> Option<Self> {
        let a = self.normalized().ok()?;
        let b = other.normalized().ok()?;
        match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::NotAromatic, Self::NotAromatic) => Some(Self::NotAromatic),
            (Self::NotAromatic, Self::Aromatic(_)) | (Self::Aromatic(_), Self::NotAromatic) => None,
            (Self::Aromatic(p), Self::Aromatic(q)) => p.meet(q).map(Self::Aromatic),
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let a = self.normalized().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.normalized().unwrap_or(Cow::Owned(Self::Undetermined));
        Ok(match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::NotAromatic, Self::NotAromatic) => Self::NotAromatic,
            (Self::NotAromatic, Self::Aromatic(_)) | (Self::Aromatic(_), Self::NotAromatic) => {
                Self::Undetermined
            }
            (Self::Aromatic(p), Self::Aromatic(q)) => Self::Aromatic(p.join(q)?),
        })
    }

    /// Partial-order check `target ⊑ self`, allocation-free — defers to the inner
    /// `NumForm::matches` for the `Aromatic` value, never building a `meet`.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, Self::Undetermined | Self::NotAromatic) => true,
            (Self::Undetermined, Self::Aromatic(v)) => v.normalized().is_ok(),
            (Self::NotAromatic, Self::NotAromatic) => true,
            (Self::Aromatic(p), Self::Aromatic(q)) => p.matches(q),
            _ => false,
        }
    }
}

impl AsLit for AromaticValenceForm {
    type Lit = AromaticValence;

    /// The exact absence or aromatic-valence value when ground.
    #[inline]
    fn as_lit(&self) -> Option<AromaticValence> {
        match self {
            Self::NotAromatic => Some(AromaticValence::NotAromatic),
            Self::Aromatic(value) => value.as_lit().map(AromaticValence::Aromatic),
            Self::Undetermined => None,
        }
    }
}

/// Covalence supplied by aromatic bonding at the given aromatic valence.
pub const fn aromatic_covalence(aromatic_valence: i64) -> i64 {
    match aromatic_valence {
        1 => 1,
        _ => 0,
    }
}

/// Multicenter-valence state of an atom: `Undetermined`, explicitly
/// `NotMulticenter`, or participating in a multicenter bond with the given
/// multicenter-valence count.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MulticenterValence {
    NotMulticenter,
    Multicenter(i64),
}

impl MulticenterValence {
    /// Multicenter-valence count for calculations that identify absence with zero.
    pub const fn valence_count(self) -> i64 {
        match self {
            Self::NotMulticenter => 0,
            Self::Multicenter(valence) => valence,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MulticenterValenceForm {
    #[default]
    Undetermined,
    NotMulticenter,
    Multicenter(NumForm),
}

impl MulticenterValenceForm {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn not_multicenter() -> Self {
        Self::NotMulticenter
    }

    pub fn multicenter(v: impl Into<NumForm>) -> Self {
        Self::Multicenter(v.into())
    }

    pub fn is_multicenter(&self) -> bool {
        matches!(self, Self::Multicenter(_))
    }

    /// Pattern matches value.
    pub fn matches_value(&self, value: i64) -> bool {
        match self {
            Self::Multicenter(v) => v.matches(&NumForm::Lit(value)),
            Self::NotMulticenter => value == 0,
            Self::Undetermined => true,
        }
    }
}

impl From<MulticenterValence> for MulticenterValenceForm {
    fn from(valence: MulticenterValence) -> Self {
        match valence {
            MulticenterValence::NotMulticenter => Self::NotMulticenter,
            MulticenterValence::Multicenter(valence) => Self::Multicenter(NumForm::Lit(valence)),
        }
    }
}

impl Normalize for MulticenterValenceForm {
    /// Delegate to the inner `NumForm`; `NotMulticenter`/`Undetermined` identity.
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Multicenter(v) => Self::Multicenter(v.normalize()?),
            other => other,
        })
    }

    fn normalized(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            Self::Multicenter(v) => Ok(match v.normalized()? {
                Cow::Borrowed(_) => Cow::Borrowed(self),
                Cow::Owned(cv) => Cow::Owned(Self::Multicenter(cv)),
            }),
            _ => Ok(Cow::Borrowed(self)),
        }
    }
}

impl Lattice for MulticenterValenceForm {
    #[inline]
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    fn is_ground(&self) -> bool {
        match self {
            Self::Undetermined => false,
            Self::NotMulticenter => true,
            Self::Multicenter(v) => v.is_ground(),
        }
    }

    /// Greatest lower bound, canonicalizing operands and output.
    fn meet(&self, other: &Self) -> Option<Self> {
        let a = self.normalized().ok()?;
        let b = other.normalized().ok()?;
        match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::NotMulticenter, Self::NotMulticenter) => Some(Self::NotMulticenter),
            (Self::NotMulticenter, Self::Multicenter(_))
            | (Self::Multicenter(_), Self::NotMulticenter) => None,
            (Self::Multicenter(p), Self::Multicenter(q)) => p.meet(q).map(Self::Multicenter),
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let a = self.normalized().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.normalized().unwrap_or(Cow::Owned(Self::Undetermined));
        Ok(match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::NotMulticenter, Self::NotMulticenter) => Self::NotMulticenter,
            (Self::NotMulticenter, Self::Multicenter(_))
            | (Self::Multicenter(_), Self::NotMulticenter) => Self::Undetermined,
            (Self::Multicenter(p), Self::Multicenter(q)) => Self::Multicenter(p.join(q)?),
        })
    }

    /// Partial-order check `target ⊑ self`, allocation-free — defers to the inner
    /// `NumForm::matches` for the `Multicenter` value, never building a `meet`.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, Self::Undetermined | Self::NotMulticenter) => true,
            (Self::Undetermined, Self::Multicenter(v)) => v.normalized().is_ok(),
            (Self::NotMulticenter, Self::NotMulticenter) => true,
            (Self::Multicenter(p), Self::Multicenter(q)) => p.matches(q),
            _ => false,
        }
    }
}

impl AsLit for MulticenterValenceForm {
    type Lit = MulticenterValence;

    /// The exact absence or multicenter-valence value when ground.
    #[inline]
    fn as_lit(&self) -> Option<MulticenterValence> {
        match self {
            Self::NotMulticenter => Some(MulticenterValence::NotMulticenter),
            Self::Multicenter(value) => value.as_lit().map(MulticenterValence::Multicenter),
            Self::Undetermined => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::{EdgeId, GraphCompaction, NodeId};

    use super::*;
    use crate::ir::num::ArithExpr;

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraintForm::valence(4), AtomConstraintForm::Valence(NumForm::Lit(4)))]
    #[case::total_valence(AtomConstraintForm::total_valence(5), AtomConstraintForm::TotalValence(NumForm::Lit(5)))]
    #[case::donated_pairs(AtomConstraintForm::donated_pairs(1), AtomConstraintForm::DonatedPairs(NumForm::Lit(1)))]
    #[case::accepted_pairs(AtomConstraintForm::accepted_pairs(2), AtomConstraintForm::AcceptedPairs(NumForm::Lit(2)))]
    #[case::degree(AtomConstraintForm::degree(3), AtomConstraintForm::Degree(NumForm::Lit(3)))]
    #[case::total_degree(AtomConstraintForm::total_degree(4), AtomConstraintForm::TotalDegree(NumForm::Lit(4)))]
    #[case::ring_degree(AtomConstraintForm::ring_degree(2), AtomConstraintForm::RingDegree(NumForm::Lit(2)))]
    #[case::ring_valence(AtomConstraintForm::ring_valence(3), AtomConstraintForm::RingValence(NumForm::Lit(3)))]
    #[case::total_hydrogens(AtomConstraintForm::total_hydrogens(3), AtomConstraintForm::TotalHydrogens(NumForm::Lit(3)))]
    #[case::ring_membership_all(AtomConstraintForm::ring_membership(RingScope::All, 1), AtomConstraintForm::RingMembership(RingMembershipForm { scope: RingScope::All, count: NumForm::Lit(1) }))]
    #[case::ring_membership_size(AtomConstraintForm::ring_membership(RingScope::Size(6), 1), AtomConstraintForm::RingMembership(RingMembershipForm { scope: RingScope::Size(6), count: NumForm::Lit(1) }))]
    #[case::aromatic_valence(
        AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic),
        AtomConstraintForm::AromaticValence(AromaticValenceForm::NotAromatic),
    )]
    #[case::multicenter_valence(
        AtomConstraintForm::multicenter_valence(MulticenterValenceForm::NotMulticenter),
        AtomConstraintForm::MulticenterValence(MulticenterValenceForm::NotMulticenter),
    )]
    #[case::tetrahedral_stereo(
        AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::NotStereo),
        AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo),
    )]
    fn test_atom_constraint_form_constructors(
        #[case] actual: AtomConstraintForm,
        #[case] expected: AtomConstraintForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraintForm::valence(4), AtomConstraintKey::Valence)]
    #[case::total_valence(AtomConstraintForm::total_valence(5), AtomConstraintKey::TotalValence)]
    #[case::aromatic_valence(AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic), AtomConstraintKey::AromaticValence)]
    #[case::multicenter_valence(AtomConstraintForm::multicenter_valence(MulticenterValenceForm::Undetermined), AtomConstraintKey::MulticenterValence)]
    #[case::donated_pairs(AtomConstraintForm::donated_pairs(1), AtomConstraintKey::DonatedPairs)]
    #[case::accepted_pairs(AtomConstraintForm::accepted_pairs(2), AtomConstraintKey::AcceptedPairs)]
    #[case::degree(AtomConstraintForm::degree(3), AtomConstraintKey::Degree)]
    #[case::total_degree(AtomConstraintForm::total_degree(4), AtomConstraintKey::TotalDegree)]
    #[case::ring_degree(AtomConstraintForm::ring_degree(2), AtomConstraintKey::RingDegree)]
    #[case::ring_valence(AtomConstraintForm::ring_valence(3), AtomConstraintKey::RingValence)]
    #[case::total_hydrogens(AtomConstraintForm::total_hydrogens(3), AtomConstraintKey::TotalHydrogens)]
    #[case::ring_membership_all(AtomConstraintForm::ring_membership(RingScope::All, 1), AtomConstraintKey::RingMembership(RingScope::All))]
    #[case::ring_membership_size(AtomConstraintForm::ring_membership(RingScope::Size(6), 1), AtomConstraintKey::RingMembership(RingScope::Size(6)))]
    #[case::tetrahedral_stereo(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo), AtomConstraintKey::TetrahedralStereo)]
    fn test_atom_constraint_form_key(
        #[case] constraint: AtomConstraintForm,
        #[case] expected: AtomConstraintKey,
    ) {
        assert_eq!(constraint.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_litset_singleton(AtomConstraintForm::Valence(NumForm::lit_set([4])), Ok(AtomConstraintForm::valence(4)))]
    #[case::ring_count_litset_singleton(
        AtomConstraintForm::RingMembership(RingMembershipForm::new(RingScope::Size(6), NumForm::lit_set([2]))),
        Ok(AtomConstraintForm::ring_membership(RingScope::Size(6), 2)))]
    #[case::empty_litset_contradiction(AtomConstraintForm::Valence(NumForm::lit_set(Vec::<i64>::new())), Err(Contradiction))]
    fn test_atom_constraint_form_normalize(
        #[case] constraint: AtomConstraintForm,
        #[case] expected: Result<AtomConstraintForm, Contradiction>,
    ) {
        assert_eq!(constraint.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_lit(AtomConstraintForm::valence(4), false)]
    #[case::valence_undetermined(AtomConstraintForm::Valence(NumForm::Undetermined), true)]
    #[case::degree_undetermined(AtomConstraintForm::Degree(NumForm::Undetermined), true)]
    #[case::ring_membership_undetermined(AtomConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined), true)]
    #[case::aromatic_undetermined(AtomConstraintForm::aromatic_valence(AromaticValenceForm::Undetermined), true)]
    #[case::aromatic_not_aromatic(AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic), false)]
    #[case::aromatic_with_value(AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(1)), false)]
    #[case::multicenter_undetermined(AtomConstraintForm::multicenter_valence(MulticenterValenceForm::Undetermined), true)]
    #[case::multicenter_not(AtomConstraintForm::multicenter_valence(MulticenterValenceForm::NotMulticenter), false)]
    #[case::multicenter_with_value(AtomConstraintForm::multicenter_valence(MulticenterValenceForm::multicenter(1)), false)]
    #[case::tetrahedral_not_stereo(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo), false)]
    #[case::tetrahedral_undetermined(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Undetermined), true)]
    fn test_atom_constraint_form_is_undetermined(
        #[case] c: AtomConstraintForm,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraintForm::valence(4), AtomConstraintForm::Valence(NumForm::Undetermined))]
    #[case::degree(AtomConstraintForm::Degree(NumForm::Lit(2)), AtomConstraintForm::Degree(NumForm::Undetermined))]
    #[case::ring_membership_keeps_scope(AtomConstraintForm::ring_membership(RingScope::Size(6), 1), AtomConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined))]
    #[case::aromatic(AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(1)), AtomConstraintForm::aromatic_valence(AromaticValenceForm::Undetermined))]
    #[case::multicenter(AtomConstraintForm::multicenter_valence(MulticenterValenceForm::multicenter(1)), AtomConstraintForm::multicenter_valence(MulticenterValenceForm::Undetermined))]
    #[case::tetrahedral(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo), AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Undetermined))]
    fn test_atom_constraint_form_as_undetermined(
        #[case] c: AtomConstraintForm,
        #[case] expected: AtomConstraintForm,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_widens(AtomConstraintForm::valence(4), AtomConstraintForm::valence(3), Ok(AtomConstraintForm::Valence(NumForm::lit_set([3, 4]))))]
    #[case::different_key(AtomConstraintForm::valence(4), AtomConstraintForm::degree(3), Err(NoJoin))]
    fn test_atom_constraint_form_join(
        #[case] a: AtomConstraintForm,
        #[case] b: AtomConstraintForm,
        #[case] expected: Result<AtomConstraintForm, NoJoin>,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rstest]
    #[case::not_aromatic(AromaticValence::NotAromatic, 0)]
    #[case::aromatic_zero(AromaticValence::Aromatic(0), 0)]
    #[case::aromatic_value(AromaticValence::Aromatic(3), 3)]
    fn test_aromatic_valence_valence_count(
        #[case] valence: AromaticValence,
        #[case] expected: i64,
    ) {
        assert_eq!(valence.valence_count(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(AromaticValenceForm::Undetermined, AromaticValenceForm::Undetermined)]
    #[case::not_aromatic(AromaticValenceForm::NotAromatic, AromaticValenceForm::NotAromatic)]
    #[case::aromatic(AromaticValenceForm::aromatic(1), AromaticValenceForm::Aromatic(NumForm::Lit(1)))]
    fn test_aromatic_valence_form_constructors(
        #[case] actual: AromaticValenceForm,
        #[case] expected: AromaticValenceForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::not_aromatic(AromaticValence::NotAromatic, AromaticValenceForm::NotAromatic)]
    #[case::aromatic(
        AromaticValence::Aromatic(1),
        AromaticValenceForm::Aromatic(NumForm::Lit(1))
    )]
    fn test_aromatic_valence_form_from(
        #[case] valence: AromaticValence,
        #[case] expected: AromaticValenceForm,
    ) {
        assert_eq!(AromaticValenceForm::from(valence), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceForm::Undetermined, false)]
    #[case::not_aromatic(AromaticValenceForm::NotAromatic, false)]
    #[case::aromatic_undetermined(AromaticValenceForm::Aromatic(NumForm::Undetermined), true)]
    #[case::aromatic_lit(AromaticValenceForm::aromatic(1), true)]
    fn test_aromatic_valence_form_is_aromatic(
        #[case] v: AromaticValenceForm,
        #[case] expected: bool,
    ) {
        assert_eq!(v.is_aromatic(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceForm::Undetermined, NumForm::Undetermined)]
    #[case::not_aromatic(AromaticValenceForm::NotAromatic, NumForm::Lit(0))]
    #[case::aromatic_undetermined(
        AromaticValenceForm::Aromatic(NumForm::Undetermined),
        NumForm::Undetermined
    )]
    #[case::aromatic_one(AromaticValenceForm::aromatic(1), NumForm::Lit(1))]
    #[case::aromatic_zero(AromaticValenceForm::aromatic(0), NumForm::Lit(0))]
    #[case::aromatic_two(AromaticValenceForm::aromatic(2), NumForm::Lit(0))]
    fn test_aromatic_valence_form_aromatic_covalence(
        #[case] v: AromaticValenceForm,
        #[case] expected: NumForm,
    ) {
        assert_eq!(v.aromatic_covalence(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceForm::Undetermined, true)]
    #[case::not_aromatic(AromaticValenceForm::NotAromatic, false)]
    #[case::aromatic_lit(AromaticValenceForm::aromatic(1), false)]
    #[case::aromatic_inner_undetermined(
        AromaticValenceForm::Aromatic(NumForm::Undetermined),
        false
    )]
    fn test_aromatic_valence_form_is_undetermined(
        #[case] v: AromaticValenceForm,
        #[case] expected: bool,
    ) {
        assert_eq!(v.is_undetermined(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceForm::Undetermined, None)]
    #[case::not_aromatic(AromaticValenceForm::NotAromatic, Some(AromaticValence::NotAromatic))]
    #[case::aromatic_undetermined(AromaticValenceForm::Aromatic(NumForm::Undetermined), None)]
    #[case::aromatic_zero(AromaticValenceForm::aromatic(0), Some(AromaticValence::Aromatic(0)))]
    #[case::aromatic_lit(AromaticValenceForm::aromatic(3), Some(AromaticValence::Aromatic(3)))]
    #[case::aromatic_term_unresolved(
        AromaticValenceForm::Aromatic(NumForm::arith_expr(ArithExpr::Lit(2))),
        None
    )]
    fn test_aromatic_valence_form_as_lit(
        #[case] v: AromaticValenceForm,
        #[case] expected: Option<AromaticValence>,
    ) {
        assert_eq!(v.as_lit(), expected);
        assert_eq!(v.is_ground(), expected.is_some());
    }

    #[rstest]
    #[case::aromatic_folds_inner(
        AromaticValenceForm::Aromatic(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(1), ArithExpr::Lit(1)]))),
        Ok(AromaticValenceForm::aromatic(2)),
    )]
    #[case::aromatic_zero_not_collapsed(
        AromaticValenceForm::aromatic(0),
        Ok(AromaticValenceForm::aromatic(0))
    )]
    fn test_aromatic_valence_form_normalize(
        #[case] input: AromaticValenceForm,
        #[case] expected: Result<AromaticValenceForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceForm::Undetermined)]
    #[case::not_aromatic(AromaticValenceForm::NotAromatic)]
    #[case::aromatic_lit(AromaticValenceForm::aromatic(1))]
    #[case::aromatic_zero(AromaticValenceForm::aromatic(0))]
    fn test_aromatic_valence_form_normalize_identity(#[case] input: AromaticValenceForm) {
        assert_eq!(input.clone().normalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_aromatic(AromaticValenceForm::Undetermined, AromaticValenceForm::aromatic(1), Some(AromaticValenceForm::aromatic(1)))]
    #[case::not_not(AromaticValenceForm::NotAromatic, AromaticValenceForm::NotAromatic, Some(AromaticValenceForm::NotAromatic))]
    #[case::not_aromatic(AromaticValenceForm::NotAromatic, AromaticValenceForm::aromatic(1), None)]
    #[case::aromatic_eq(AromaticValenceForm::aromatic(1), AromaticValenceForm::aromatic(1), Some(AromaticValenceForm::aromatic(1)))]
    #[case::aromatic_neq(AromaticValenceForm::aromatic(1), AromaticValenceForm::aromatic(2), None)]
    #[case::aromatic_inner_wildcard(AromaticValenceForm::Aromatic(NumForm::Undetermined), AromaticValenceForm::aromatic(2), Some(AromaticValenceForm::aromatic(2)))]
    fn test_aromatic_valence_form_meet(
        #[case] a: AromaticValenceForm,
        #[case] b: AromaticValenceForm,
        #[case] expected: Option<AromaticValenceForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und(AromaticValenceForm::Undetermined, AromaticValenceForm::aromatic(1), AromaticValenceForm::Undetermined)]
    #[case::not_not(AromaticValenceForm::NotAromatic, AromaticValenceForm::NotAromatic, AromaticValenceForm::NotAromatic)]
    #[case::not_aromatic(AromaticValenceForm::NotAromatic, AromaticValenceForm::aromatic(1), AromaticValenceForm::Undetermined)]
    #[case::aromatic_eq(AromaticValenceForm::aromatic(1), AromaticValenceForm::aromatic(1), AromaticValenceForm::aromatic(1))]
    #[case::aromatic_inner_wildcard(AromaticValenceForm::Aromatic(NumForm::Undetermined), AromaticValenceForm::aromatic(1), AromaticValenceForm::Aromatic(NumForm::Undetermined))]
    fn test_aromatic_valence_form_join(
        #[case] a: AromaticValenceForm,
        #[case] b: AromaticValenceForm,
        #[case] expected: AromaticValenceForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rstest]
    #[case::wildcard_vs_not_aromatic(
        AromaticValenceForm::Undetermined,
        AromaticValenceForm::NotAromatic,
        true
    )]
    #[case::wildcard_vs_aromatic_lit(
        AromaticValenceForm::Undetermined,
        AromaticValenceForm::aromatic(1),
        true
    )]
    #[case::not_aromatic_vs_aromatic(
        AromaticValenceForm::NotAromatic,
        AromaticValenceForm::aromatic(1),
        false
    )]
    #[case::aromatic_vs_not_aromatic(
        AromaticValenceForm::aromatic(1),
        AromaticValenceForm::NotAromatic,
        false
    )]
    #[case::not_aromatic_vs_not_aromatic(
        AromaticValenceForm::NotAromatic,
        AromaticValenceForm::NotAromatic,
        true
    )]
    #[case::aromatic_lit_match(
        AromaticValenceForm::aromatic(1),
        AromaticValenceForm::aromatic(1),
        true
    )]
    #[case::aromatic_lit_mismatch(
        AromaticValenceForm::aromatic(1),
        AromaticValenceForm::aromatic(2),
        false
    )]
    #[case::aromatic_wildcard_inner(
        AromaticValenceForm::Aromatic(NumForm::Undetermined),
        AromaticValenceForm::aromatic(2),
        true
    )]
    #[case::specific_vs_undetermined(
        AromaticValenceForm::aromatic(1),
        AromaticValenceForm::Undetermined,
        false
    )]
    fn test_aromatic_valence_form_matches(
        #[case] pattern: AromaticValenceForm,
        #[case] target: AromaticValenceForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::not_multicenter(MulticenterValence::NotMulticenter, 0)]
    #[case::multicenter_zero(MulticenterValence::Multicenter(0), 0)]
    #[case::multicenter_value(MulticenterValence::Multicenter(2), 2)]
    fn test_multicenter_valence_valence_count(
        #[case] valence: MulticenterValence,
        #[case] expected: i64,
    ) {
        assert_eq!(valence.valence_count(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterValenceForm::Undetermined, MulticenterValenceForm::Undetermined)]
    #[case::not_multicenter(MulticenterValenceForm::NotMulticenter, MulticenterValenceForm::NotMulticenter)]
    #[case::multicenter(MulticenterValenceForm::multicenter(2), MulticenterValenceForm::Multicenter(NumForm::Lit(2)))]
    fn test_multicenter_valence_form_constructors(
        #[case] actual: MulticenterValenceForm,
        #[case] expected: MulticenterValenceForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::not_multicenter(
        MulticenterValence::NotMulticenter,
        MulticenterValenceForm::NotMulticenter
    )]
    #[case::multicenter(
        MulticenterValence::Multicenter(2),
        MulticenterValenceForm::Multicenter(NumForm::Lit(2))
    )]
    fn test_multicenter_valence_form_from(
        #[case] valence: MulticenterValence,
        #[case] expected: MulticenterValenceForm,
    ) {
        assert_eq!(MulticenterValenceForm::from(valence), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterValenceForm::Undetermined, false)]
    #[case::not_multicenter(MulticenterValenceForm::NotMulticenter, false)]
    #[case::multicenter_undetermined(MulticenterValenceForm::Multicenter(NumForm::Undetermined), true)]
    #[case::multicenter_lit(MulticenterValenceForm::multicenter(1), true)]
    fn test_multicenter_valence_form_is_multicenter(
        #[case] v: MulticenterValenceForm,
        #[case] expected: bool,
    ) {
        assert_eq!(v.is_multicenter(), expected);
    }

    #[rstest]
    #[case::undetermined(MulticenterValenceForm::Undetermined, true)]
    #[case::not_multicenter(MulticenterValenceForm::NotMulticenter, false)]
    #[case::multicenter_lit(MulticenterValenceForm::multicenter(1), false)]
    fn test_multicenter_valence_form_is_undetermined(
        #[case] v: MulticenterValenceForm,
        #[case] expected: bool,
    ) {
        assert_eq!(v.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterValenceForm::Undetermined, None)]
    #[case::not_multicenter(MulticenterValenceForm::NotMulticenter, Some(MulticenterValence::NotMulticenter))]
    #[case::multicenter_undetermined(MulticenterValenceForm::Multicenter(NumForm::Undetermined), None)]
    #[case::multicenter_zero(MulticenterValenceForm::multicenter(0), Some(MulticenterValence::Multicenter(0)))]
    #[case::multicenter_lit(MulticenterValenceForm::multicenter(2), Some(MulticenterValence::Multicenter(2)))]
    #[case::multicenter_term_unresolved(MulticenterValenceForm::Multicenter(NumForm::arith_expr(ArithExpr::Lit(3))), None)]
    fn test_multicenter_valence_form_as_lit(
        #[case] v: MulticenterValenceForm,
        #[case] expected: Option<MulticenterValence>,
    ) {
        assert_eq!(v.as_lit(), expected);
        assert_eq!(v.is_ground(), expected.is_some());
    }

    #[rstest]
    #[case::multicenter_folds_inner(
        MulticenterValenceForm::Multicenter(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(1), ArithExpr::Lit(2)]))),
        Ok(MulticenterValenceForm::multicenter(3)),
    )]
    #[case::multicenter_zero_not_collapsed(
        MulticenterValenceForm::multicenter(0),
        Ok(MulticenterValenceForm::multicenter(0))
    )]
    fn test_multicenter_valence_form_normalize(
        #[case] input: MulticenterValenceForm,
        #[case] expected: Result<MulticenterValenceForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rstest]
    #[case::undetermined(MulticenterValenceForm::Undetermined)]
    #[case::not_multicenter(MulticenterValenceForm::NotMulticenter)]
    #[case::multicenter_lit(MulticenterValenceForm::multicenter(1))]
    #[case::multicenter_zero(MulticenterValenceForm::multicenter(0))]
    fn test_multicenter_valence_form_normalize_identity(#[case] input: MulticenterValenceForm) {
        assert_eq!(input.clone().normalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_multicenter(MulticenterValenceForm::Undetermined, MulticenterValenceForm::multicenter(2), Some(MulticenterValenceForm::multicenter(2)))]
    #[case::not_not(MulticenterValenceForm::NotMulticenter, MulticenterValenceForm::NotMulticenter, Some(MulticenterValenceForm::NotMulticenter))]
    #[case::not_multicenter(MulticenterValenceForm::NotMulticenter, MulticenterValenceForm::multicenter(2), None)]
    #[case::multicenter_eq(MulticenterValenceForm::multicenter(2), MulticenterValenceForm::multicenter(2), Some(MulticenterValenceForm::multicenter(2)))]
    #[case::multicenter_neq(MulticenterValenceForm::multicenter(2), MulticenterValenceForm::multicenter(3), None)]
    #[case::multicenter_inner_wildcard(MulticenterValenceForm::Multicenter(NumForm::Undetermined), MulticenterValenceForm::multicenter(3), Some(MulticenterValenceForm::multicenter(3)))]
    fn test_multicenter_valence_form_meet(
        #[case] a: MulticenterValenceForm,
        #[case] b: MulticenterValenceForm,
        #[case] expected: Option<MulticenterValenceForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und(MulticenterValenceForm::Undetermined, MulticenterValenceForm::multicenter(2), MulticenterValenceForm::Undetermined)]
    #[case::not_not(MulticenterValenceForm::NotMulticenter, MulticenterValenceForm::NotMulticenter, MulticenterValenceForm::NotMulticenter)]
    #[case::not_multicenter(MulticenterValenceForm::NotMulticenter, MulticenterValenceForm::multicenter(2), MulticenterValenceForm::Undetermined)]
    #[case::multicenter_eq(MulticenterValenceForm::multicenter(2), MulticenterValenceForm::multicenter(2), MulticenterValenceForm::multicenter(2))]
    #[case::multicenter_inner_wildcard(MulticenterValenceForm::Multicenter(NumForm::Undetermined), MulticenterValenceForm::multicenter(2), MulticenterValenceForm::Multicenter(NumForm::Undetermined))]
    fn test_multicenter_valence_form_join(
        #[case] a: MulticenterValenceForm,
        #[case] b: MulticenterValenceForm,
        #[case] expected: MulticenterValenceForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard_vs_not_multicenter(MulticenterValenceForm::Undetermined, MulticenterValenceForm::NotMulticenter, true)]
    #[case::wildcard_vs_multicenter_lit(MulticenterValenceForm::Undetermined, MulticenterValenceForm::multicenter(2), true)]
    #[case::not_multicenter_vs_multicenter(MulticenterValenceForm::NotMulticenter, MulticenterValenceForm::multicenter(2), false)]
    #[case::multicenter_vs_not_multicenter(MulticenterValenceForm::multicenter(2), MulticenterValenceForm::NotMulticenter, false)]
    #[case::not_multicenter_vs_not_multicenter(MulticenterValenceForm::NotMulticenter, MulticenterValenceForm::NotMulticenter, true)]
    #[case::multicenter_lit_match(MulticenterValenceForm::multicenter(2), MulticenterValenceForm::multicenter(2), true)]
    #[case::multicenter_lit_mismatch(MulticenterValenceForm::multicenter(2), MulticenterValenceForm::multicenter(3), false)]
    #[case::specific_vs_undetermined(MulticenterValenceForm::multicenter(2), MulticenterValenceForm::Undetermined, false)]
    fn test_multicenter_valence_form_matches(
        #[case] pattern: MulticenterValenceForm,
        #[case] target: MulticenterValenceForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_atom_constraints_form_new() {
        let cs = AtomConstraintsForm::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_present(AtomConstraintKey::Valence, true)]
    #[case::ring_all_present(AtomConstraintKey::RingMembership(RingScope::All), true)]
    #[case::ring_size_present(AtomConstraintKey::RingMembership(RingScope::Size(6)), true)]
    #[case::ring_size_absent(AtomConstraintKey::RingMembership(RingScope::Size(5)), false)]
    #[case::degree_absent(AtomConstraintKey::Degree, false)]
    fn test_atom_constraints_form_contains(
        #[case] key: AtomConstraintKey,
        #[case] expected: bool,
    ) {
        let cs = AtomConstraintsForm::from_iter([
            AtomConstraintForm::valence(4),
            AtomConstraintForm::ring_membership(RingScope::All, 2),
            AtomConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_all(AtomConstraintKey::RingMembership(RingScope::All), Some(AtomConstraintForm::ring_membership(RingScope::All, 2)))]
    #[case::ring_size(AtomConstraintKey::RingMembership(RingScope::Size(6)), Some(AtomConstraintForm::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_size_absent(AtomConstraintKey::RingMembership(RingScope::Size(5)), None)]
    #[case::valence(AtomConstraintKey::Valence, Some(AtomConstraintForm::valence(4)))]
    fn test_atom_constraints_form_get(
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        let cs = AtomConstraintsForm::from_iter([
            AtomConstraintForm::valence(4),
            AtomConstraintForm::ring_membership(RingScope::All, 2),
            AtomConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(key), expected.as_ref());
    }

    #[rstest]
    fn test_atom_constraints_form_remove() {
        let mut cs = AtomConstraintsForm::from_iter([
            AtomConstraintForm::valence(4),
            AtomConstraintForm::ring_membership(RingScope::All, 2),
            AtomConstraintForm::ring_membership(RingScope::Size(6), 1),
        ]);
        let removed = cs.remove(AtomConstraintKey::RingMembership(RingScope::Size(6)));
        assert_eq!(
            removed,
            Some(AtomConstraintForm::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(
            cs.iter().cloned().collect::<Vec<_>>(),
            vec![
                AtomConstraintForm::valence(4),
                AtomConstraintForm::ring_membership(RingScope::All, 2),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::drop_vacuous(
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::Valence(NumForm::Undetermined),
            AtomConstraintForm::degree(3),
        ]),
        Ok(AtomConstraintsForm::from_iter([AtomConstraintForm::degree(3)])))]
    #[case::normalizes_values(
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::Degree(NumForm::lit_set([3])),
            AtomConstraintForm::Valence(NumForm::lit_set([4])),
        ]),
        Ok(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4), AtomConstraintForm::degree(3)])))]
    fn test_atom_constraints_form_normalize(
        #[case] constraints: AtomConstraintsForm,
        #[case] expected: Result<AtomConstraintsForm, Contradiction>,
    ) {
        assert_eq!(constraints.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![AtomConstraintForm::valence(4)], vec![AtomConstraintForm::valence(4)])]
    #[case::overwrite_same_key(vec![AtomConstraintForm::valence(3), AtomConstraintForm::valence(4)], vec![AtomConstraintForm::valence(4)])]
    #[case::vacuous_overwrites(vec![AtomConstraintForm::valence(4), AtomConstraintForm::Valence(NumForm::Undetermined)], vec![AtomConstraintForm::Valence(NumForm::Undetermined)])]
    #[case::vacuous_absent_inserts(vec![AtomConstraintForm::Valence(NumForm::Undetermined)], vec![AtomConstraintForm::Valence(NumForm::Undetermined)])]
    #[case::new_key_sorts(vec![AtomConstraintForm::degree(3), AtomConstraintForm::valence(4)], vec![AtomConstraintForm::valence(4), AtomConstraintForm::degree(3)])]
    #[case::ring_overwrite_scope(vec![AtomConstraintForm::ring_membership(RingScope::Size(6), 1), AtomConstraintForm::ring_membership(RingScope::Size(6), 2)], vec![AtomConstraintForm::ring_membership(RingScope::Size(6), 2)])]
    #[case::ring_vacuous_overwrites_scope(vec![AtomConstraintForm::ring_membership(RingScope::All, 2), AtomConstraintForm::ring_membership(RingScope::Size(6), 1), AtomConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)], vec![AtomConstraintForm::ring_membership(RingScope::All, 2), AtomConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)])]
    fn test_atom_constraints_form_set(
        #[case] sequence: Vec<AtomConstraintForm>,
        #[case] expected_state: Vec<AtomConstraintForm>,
    ) {
        let mut cs = AtomConstraintsForm::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, AtomConstraintsForm::from_iter(expected_state));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite_shared(
        vec![AtomConstraintForm::valence(3), AtomConstraintForm::degree(3)],
        vec![AtomConstraintForm::valence(4)],
        vec![AtomConstraintForm::valence(4), AtomConstraintForm::degree(3)])]
    #[case::keeps_disjoint(
        vec![AtomConstraintForm::valence(4)],
        vec![AtomConstraintForm::degree(3)],
        vec![AtomConstraintForm::valence(4), AtomConstraintForm::degree(3)])]
    #[case::vacuous_removes(
        vec![AtomConstraintForm::valence(4), AtomConstraintForm::degree(3)],
        vec![AtomConstraintForm::Valence(NumForm::Undetermined)],
        vec![AtomConstraintForm::degree(3)])]
    #[case::empty_other(
        vec![AtomConstraintForm::valence(4)],
        vec![],
        vec![AtomConstraintForm::valence(4)])]
    #[case::ring_overwrite_scope(
        vec![AtomConstraintForm::ring_membership(RingScope::All, 2), AtomConstraintForm::ring_membership(RingScope::Size(6), 1)],
        vec![AtomConstraintForm::ring_membership(RingScope::Size(6), 3)],
        vec![AtomConstraintForm::ring_membership(RingScope::All, 2), AtomConstraintForm::ring_membership(RingScope::Size(6), 3)])]
    fn test_atom_constraints_form_update(
        #[case] initial: Vec<AtomConstraintForm>,
        #[case] other: Vec<AtomConstraintForm>,
        #[case] expected: Vec<AtomConstraintForm>,
    ) {
        let mut cs = AtomConstraintsForm::from_iter(initial);
        let overlay = AtomConstraintsForm::from_iter(other);
        cs.update(&overlay);
        assert_eq!(cs, AtomConstraintsForm::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![AtomConstraintForm::valence(3)], Some(AtomConstraintForm::valence(3)), Some(AtomConstraintForm::valence(4)), Ok(()), vec![AtomConstraintForm::valence(4)])]
    #[case::remove(vec![AtomConstraintForm::valence(4)], Some(AtomConstraintForm::valence(4)), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(AtomConstraintForm::valence(4)), Ok(()), vec![AtomConstraintForm::valence(4)])]
    #[case::normalized_match(vec![AtomConstraintForm::Valence(NumForm::lit_set([4]))], Some(AtomConstraintForm::valence(4)), Some(AtomConstraintForm::valence(5)), Ok(()), vec![AtomConstraintForm::valence(5)])]
    #[case::old_mismatch(vec![AtomConstraintForm::valence(3)], Some(AtomConstraintForm::valence(4)), Some(AtomConstraintForm::valence(5)), Err(Contradiction), vec![AtomConstraintForm::valence(3)])]
    #[case::old_absent_mismatch(vec![AtomConstraintForm::valence(3)], None, Some(AtomConstraintForm::valence(4)), Err(Contradiction), vec![AtomConstraintForm::valence(3)])]
    #[case::key_mismatch(vec![], Some(AtomConstraintForm::valence(3)), Some(AtomConstraintForm::degree(4)), Err(Contradiction), vec![])]
    #[case::noop(vec![AtomConstraintForm::valence(4)], None, None, Ok(()), vec![AtomConstraintForm::valence(4)])]
    fn test_atom_constraints_form_compare_and_set(
        #[case] initial: Vec<AtomConstraintForm>,
        #[case] old: Option<AtomConstraintForm>,
        #[case] new: Option<AtomConstraintForm>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<AtomConstraintForm>,
    ) {
        let mut cs = AtomConstraintsForm::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, AtomConstraintsForm::from_iter(expected_state));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &AtomConstraintForm| matches!(c, AtomConstraintForm::Valence(_) | AtomConstraintForm::RingMembership(_)), vec![AtomConstraintForm::valence(4), AtomConstraintForm::ring_membership(RingScope::All, 2)])]
    #[case::all_dropped(|_: &AtomConstraintForm| false, vec![])]
    fn test_atom_constraints_form_retain(
        #[case] predicate: impl FnMut(&AtomConstraintForm) -> bool,
        #[case] expected: Vec<AtomConstraintForm>,
    ) {
        let mut cs = AtomConstraintsForm::from_iter([
            AtomConstraintForm::valence(4),
            AtomConstraintForm::degree(3),
            AtomConstraintForm::ring_membership(RingScope::All, 2),
        ]);
        cs.retain(predicate);
        assert_eq!(cs, AtomConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_atom_constraints_form_clear() {
        let mut cs = AtomConstraintsForm::from_iter([
            AtomConstraintForm::valence(4),
            AtomConstraintForm::degree(3),
        ]);
        cs.clear();
        assert_eq!(cs, AtomConstraintsForm::new());
    }

    #[rstest]
    fn test_atom_constraints_form_take() {
        let mut empty = AtomConstraintsForm::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let mut cs = AtomConstraintsForm::from_iter([
            AtomConstraintForm::valence(4),
            AtomConstraintForm::degree(3),
        ]);
        let mut taken = cs.take();
        assert_eq!(taken.len(), 2);
        assert_eq!(taken.size_hint(), (2, Some(2)));
        assert_eq!(taken.next(), Some(AtomConstraintForm::valence(4)));
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(taken.next(), Some(AtomConstraintForm::degree(3)));
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(cs, AtomConstraintsForm::new());
    }

    #[rstest]
    fn test_atom_constraints_form_iter() {
        let empty = AtomConstraintsForm::new();
        let mut empty_iter = empty.iter();
        assert_eq!(empty_iter.len(), 0);
        assert_eq!(empty_iter.size_hint(), (0, Some(0)));
        assert_eq!(empty_iter.next(), None);

        let cs = AtomConstraintsForm::from_iter([
            AtomConstraintForm::ring_membership(RingScope::Size(6), 1),
            AtomConstraintForm::valence(4),
            AtomConstraintForm::degree(3),
        ]);
        let mut iter = cs.iter();
        assert_eq!(iter.len(), 3);
        assert_eq!(iter.size_hint(), (3, Some(3)));
        assert_eq!(iter.next(), Some(&AtomConstraintForm::valence(4)));
        assert_eq!(iter.len(), 2);
        assert_eq!(iter.size_hint(), (2, Some(2)));
        assert_eq!(iter.next(), Some(&AtomConstraintForm::degree(3)));
        assert_eq!(
            iter.next(),
            Some(&AtomConstraintForm::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(iter.len(), 0);
        assert_eq!(iter.next(), None);
    }

    #[rstest]
    fn test_atom_constraints_form_compact() {
        let cs = AtomConstraintsForm::from_iter([
            AtomConstraintForm::valence(4),
            AtomConstraintForm::degree(3),
        ]);
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(vec![NodeId(0), NodeId(1), NodeId(2)], vec![EdgeId(0)]),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(cs.clone().compact(&compaction), cs);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(vec![AtomConstraintForm::valence(4), AtomConstraintForm::degree(3)], vec![AtomConstraintForm::valence(4), AtomConstraintForm::degree(3)])]
    #[case::same_kind_last_wins(vec![AtomConstraintForm::valence(3), AtomConstraintForm::valence(4)], vec![AtomConstraintForm::valence(4)])]
    #[case::empty(vec![], vec![])]
    fn test_atom_constraints_form_from_iter(
        #[case] input: Vec<AtomConstraintForm>,
        #[case] expected: Vec<AtomConstraintForm>,
    ) {
        let cs = AtomConstraintsForm::from_iter(input);
        assert_eq!(cs, AtomConstraintsForm::from_iter(expected));
    }

    #[rstest]
    fn test_atom_constraints_form_into_iter() {
        let cs = AtomConstraintsForm::from_iter([
            AtomConstraintForm::valence(4),
            AtomConstraintForm::degree(3),
        ]);
        let collected: Vec<AtomConstraintForm> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                AtomConstraintForm::valence(4),
                AtomConstraintForm::degree(3)
            ],
        );
    }

    #[rstest]
    fn test_atom_constraints_form_from_atom_constraint() {
        let cs: AtomConstraintsForm = AtomConstraintForm::valence(4).into();
        assert_eq!(
            cs,
            AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)])
        );
    }

    #[rstest]
    fn test_atom_constraints_form_from_vec() {
        let cs: AtomConstraintsForm = vec![
            AtomConstraintForm::valence(4),
            AtomConstraintForm::donated_pairs(1),
        ]
        .into();
        assert_eq!(
            cs,
            AtomConstraintsForm::from_iter([
                AtomConstraintForm::valence(4),
                AtomConstraintForm::donated_pairs(1),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_empty(AtomConstraintsForm::new(), AtomConstraintsForm::new(), Some(AtomConstraintsForm::new()))]
    #[case::adds_kind_from_other(AtomConstraintsForm::new(), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), Some(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)])))]
    #[case::narrows_undetermined_to_lit(AtomConstraintsForm::from_iter([AtomConstraintForm::Valence(NumForm::Undetermined)]), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]),
        Some(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)])))]
    #[case::lit_lit_match_preserved(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]),
        Some(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)])))]
    #[case::lit_lit_mismatch_none(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(3)]), None)]
    #[case::multi_kind_combines(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), AtomConstraintsForm::from_iter([AtomConstraintForm::degree(3)]),
        Some(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4), AtomConstraintForm::degree(3)])))]
    #[case::aromatic_valence_narrows(AtomConstraintsForm::from_iter([AtomConstraintForm::aromatic_valence(AromaticValenceForm::Undetermined)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(1))]),
        Some(AtomConstraintsForm::from_iter([AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(1))])))]
    #[case::aromatic_valence_not_vs_aromatic_none(AtomConstraintsForm::from_iter([AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(1))]), None)]
    #[case::ring_membership_size_unions(AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(5), 1)]), AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(6), 1)]),
        Some(AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(5), 1), AtomConstraintForm::ring_membership(RingScope::Size(6), 1)])))]
    #[case::ring_membership_size_dedup(AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(5), 1)]), AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(5), 1)]),
        Some(AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(5), 1)])))]
    #[case::prunes_vacuous(AtomConstraintsForm::new(), AtomConstraintsForm::from_iter([AtomConstraintForm::Valence(NumForm::Undetermined)]), Some(AtomConstraintsForm::new()))]
    #[case::tetrahedral_narrows_from_absent(AtomConstraintsForm::new(),
        AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)]),
        Some(AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)])))]
    #[case::tetrahedral_not_stereo_vs_stereo_contradicts(AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::stereo(0_u32))]), None)]
    fn test_atom_constraints_form_meet(
        #[case] a: AtomConstraintsForm,
        #[case] b: AtomConstraintsForm,
        #[case] expected: Option<AtomConstraintsForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::extends_self(AtomConstraintsForm::new(), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), true, AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]))]
    #[case::no_change(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), false,
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]))]
    #[case::contradiction_leaves_self_unchanged(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(3)]), false,
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]))]
    fn test_atom_constraints_form_narrow_from(
        #[case] mut target: AtomConstraintsForm,
        #[case] source: AtomConstraintsForm,
        #[case] expected_changed: bool,
        #[case] expected_after: AtomConstraintsForm,
    ) {
        let changed = target.narrow_from(&source);
        assert_eq!(changed, expected_changed);
        assert_eq!(target, expected_after);
    }

    #[rstest]
    #[case::keeps_only_shared_kinds(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4), AtomConstraintForm::degree(2)]), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]))]
    #[case::widens_value(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(3)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::Valence(NumForm::lit_set([4, 3]))]))]
    #[case::tetrahedral_same(AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)]))]
    #[case::tetrahedral_incompatible_drops_to_undetermined(AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::stereo(0_u32))]), AtomConstraintsForm::new())]
    fn test_atom_constraints_form_join(
        #[case] a: AtomConstraintsForm,
        #[case] b: AtomConstraintsForm,
        #[case] expected: AtomConstraintsForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rstest]
    #[case::empty_pattern_matches_anything(AtomConstraintsForm::new(), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), true)]
    #[case::missing_in_target_when_pattern_specific(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), AtomConstraintsForm::new(), false)]
    #[case::same_lit(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), true)]
    #[case::lit_lit_mismatch(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]), AtomConstraintsForm::from_iter([AtomConstraintForm::valence(3)]), false)]
    #[case::aromatic_wildcard_matches_aromatic(AtomConstraintsForm::from_iter([AtomConstraintForm::aromatic_valence(AromaticValenceForm::Undetermined)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(1))]), true)]
    #[case::aromatic_not_vs_aromatic_mismatch(AtomConstraintsForm::from_iter([AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(1))]), false)]
    #[case::ring_membership_size_subset(AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(5), 1)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(5), 1), AtomConstraintForm::ring_membership(RingScope::Size(6), 1)]), true)]
    #[case::ring_membership_size_not_present_in_target(AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(7), 1)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(5), 1), AtomConstraintForm::ring_membership(RingScope::Size(6), 1)]), false)]
    #[case::multi_kind_all_must_match(AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4), AtomConstraintForm::degree(3)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4), AtomConstraintForm::degree(2)]), false)]
    #[case::tetrahedral_same(AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)]),
        AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)]), true)]
    #[case::tetrahedral_pattern_specific_vs_absent(AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)]),
        AtomConstraintsForm::new(), false)]
    fn test_atom_constraints_form_matches(
        #[case] pattern: AtomConstraintsForm,
        #[case] target: AtomConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
