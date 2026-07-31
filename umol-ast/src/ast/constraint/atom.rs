//! Atom constraints.

use std::borrow::Cow;
use std::cmp::Ordering;
use std::mem;

use smallvec::SmallVec;

use super::super::constraint::ring::{RingMembershipAst, RingScope};
use super::super::error::{Contradiction, NoJoin};
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::stereo::TetrahedralStereoAst;
use super::super::traits::{AsLit, Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Atom-scope constraint.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AtomConstraintAst {
    Valence(ValueAst),
    DonatedPairs(ValueAst),
    AcceptedPairs(ValueAst),
    AromaticValence(AromaticValenceAst),
    MulticenterValence(MulticenterValenceAst),
    TetrahedralStereo(TetrahedralStereoAst),
    Degree(ValueAst),
    TotalDegree(ValueAst),
    TotalValence(ValueAst),
    /// Incident bonds belonging to the fixed Relevant ring projection.
    RingDegree(ValueAst),
    /// Sum of incident bond orders in the fixed Relevant ring projection.
    RingValence(ValueAst),
    TotalHydrogens(ValueAst),
    /// Ring count in the fixed Relevant ring projection, optionally restricted by size.
    RingMembership(RingMembershipAst),
}

impl AtomConstraintAst {
    pub fn valence(v: impl Into<ValueAst>) -> Self {
        Self::Valence(v.into())
    }

    pub fn donated_pairs(v: impl Into<ValueAst>) -> Self {
        Self::DonatedPairs(v.into())
    }

    pub fn accepted_pairs(v: impl Into<ValueAst>) -> Self {
        Self::AcceptedPairs(v.into())
    }

    pub fn aromatic_valence(v: impl Into<AromaticValenceAst>) -> Self {
        Self::AromaticValence(v.into())
    }

    pub fn multicenter_valence(v: impl Into<MulticenterValenceAst>) -> Self {
        Self::MulticenterValence(v.into())
    }

    pub fn tetrahedral_stereo(c: impl Into<TetrahedralStereoAst>) -> Self {
        Self::TetrahedralStereo(c.into())
    }

    pub fn degree(v: impl Into<ValueAst>) -> Self {
        Self::Degree(v.into())
    }

    pub fn total_degree(v: impl Into<ValueAst>) -> Self {
        Self::TotalDegree(v.into())
    }

    pub fn total_valence(v: impl Into<ValueAst>) -> Self {
        Self::TotalValence(v.into())
    }

    pub fn ring_degree(v: impl Into<ValueAst>) -> Self {
        Self::RingDegree(v.into())
    }

    pub fn ring_valence(v: impl Into<ValueAst>) -> Self {
        Self::RingValence(v.into())
    }

    pub fn total_hydrogens(v: impl Into<ValueAst>) -> Self {
        Self::TotalHydrogens(v.into())
    }

    pub fn ring_membership(scope: RingScope, count: impl Into<ValueAst>) -> Self {
        Self::RingMembership(RingMembershipAst::new(scope, count))
    }

    /// Atom constraint key, unique within an `AtomConstraintsAst` container.
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
            Self::Valence(_) => Self::Valence(ValueAst::Undetermined),
            Self::DonatedPairs(_) => Self::DonatedPairs(ValueAst::Undetermined),
            Self::AcceptedPairs(_) => Self::AcceptedPairs(ValueAst::Undetermined),
            Self::AromaticValence(_) => Self::AromaticValence(AromaticValenceAst::Undetermined),
            Self::MulticenterValence(_) => {
                Self::MulticenterValence(MulticenterValenceAst::Undetermined)
            }
            Self::TetrahedralStereo(_) => {
                Self::TetrahedralStereo(TetrahedralStereoAst::Undetermined)
            }
            Self::Degree(_) => Self::Degree(ValueAst::Undetermined),
            Self::TotalDegree(_) => Self::TotalDegree(ValueAst::Undetermined),
            Self::TotalValence(_) => Self::TotalValence(ValueAst::Undetermined),
            Self::RingDegree(_) => Self::RingDegree(ValueAst::Undetermined),
            Self::RingValence(_) => Self::RingValence(ValueAst::Undetermined),
            Self::TotalHydrogens(_) => Self::TotalHydrogens(ValueAst::Undetermined),
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipAst::new(m.scope, ValueAst::Undetermined))
            }
        }
    }

    /// Value-only payload: no entity ids to compact, so this never drops.
    pub fn compact(self, _compaction: &IdCompaction) -> Option<Self> {
        Some(self)
    }

    /// Value-only payload: no entity ids to remap.
    pub fn remap(self, _map: &IdRemapping) -> Self {
        self
    }
}

impl Canonicalize for AtomConstraintAst {
    /// Canonicalize the inner value; kind and sub-key are preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Valence(v) => Self::Valence(v.canonicalize()?),
            Self::DonatedPairs(v) => Self::DonatedPairs(v.canonicalize()?),
            Self::AcceptedPairs(v) => Self::AcceptedPairs(v.canonicalize()?),
            Self::AromaticValence(c) => Self::AromaticValence(c.canonicalize()?),
            Self::MulticenterValence(c) => Self::MulticenterValence(c.canonicalize()?),
            Self::TetrahedralStereo(c) => Self::TetrahedralStereo(c.canonicalize()?),
            Self::Degree(v) => Self::Degree(v.canonicalize()?),
            Self::TotalDegree(v) => Self::TotalDegree(v.canonicalize()?),
            Self::TotalValence(v) => Self::TotalValence(v.canonicalize()?),
            Self::RingDegree(v) => Self::RingDegree(v.canonicalize()?),
            Self::RingValence(v) => Self::RingValence(v.canonicalize()?),
            Self::TotalHydrogens(v) => Self::TotalHydrogens(v.canonicalize()?),
            Self::RingMembership(m) => Self::RingMembership(m.canonicalize()?),
        })
    }
}

impl Lattice for AtomConstraintAst {
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
    fn meet(&self, other: &AtomConstraintAst) -> Option<AtomConstraintAst> {
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
    fn join(&self, other: &AtomConstraintAst) -> Result<AtomConstraintAst, NoJoin> {
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
    /// `AtomConstraintAst`.
    fn is_compatible(&self, other: &AtomConstraintAst) -> bool {
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

/// Entry identity: discriminant + sub-key, AtomConstraintsAst is ordered, unique by key.
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
pub struct AtomConstraintsAst {
    entries: SmallVec<[AtomConstraintAst; 2]>,
}

impl AtomConstraintsAst {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn valence(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKey::Valence) {
            Some(AtomConstraintAst::Valence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn aromatic_valence(&self) -> Option<&AromaticValenceAst> {
        match self.get(AtomConstraintKey::AromaticValence) {
            Some(AtomConstraintAst::AromaticValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn multicenter_valence(&self) -> Option<&MulticenterValenceAst> {
        match self.get(AtomConstraintKey::MulticenterValence) {
            Some(AtomConstraintAst::MulticenterValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn tetrahedral_stereo(&self) -> Option<&TetrahedralStereoAst> {
        match self.get(AtomConstraintKey::TetrahedralStereo) {
            Some(AtomConstraintAst::TetrahedralStereo(c)) => Some(c),
            _ => None,
        }
    }

    pub fn degree(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKey::Degree) {
            Some(AtomConstraintAst::Degree(v)) => Some(v),
            _ => None,
        }
    }

    pub fn total_degree(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKey::TotalDegree) {
            Some(AtomConstraintAst::TotalDegree(v)) => Some(v),
            _ => None,
        }
    }

    pub fn total_valence(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKey::TotalValence) {
            Some(AtomConstraintAst::TotalValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn ring_degree(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKey::RingDegree) {
            Some(AtomConstraintAst::RingDegree(v)) => Some(v),
            _ => None,
        }
    }

    pub fn ring_valence(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKey::RingValence) {
            Some(AtomConstraintAst::RingValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn total_hydrogens(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKey::TotalHydrogens) {
            Some(AtomConstraintAst::TotalHydrogens(v)) => Some(v),
            _ => None,
        }
    }

    pub fn donated_pairs(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKey::DonatedPairs) {
            Some(AtomConstraintAst::DonatedPairs(v)) => Some(v),
            _ => None,
        }
    }

    pub fn accepted_pairs(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKey::AcceptedPairs) {
            Some(AtomConstraintAst::AcceptedPairs(v)) => Some(v),
            _ => None,
        }
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &ValueAst)> {
        self.iter().filter_map(|c| match c {
            AtomConstraintAst::RingMembership(m) => Some((m.scope, &m.count)),
            _ => None,
        })
    }

    fn ring_membership(&self, scope: RingScope) -> Option<&ValueAst> {
        self.ring_memberships()
            .find(|(s, _)| *s == scope)
            .map(|(_, v)| v)
    }

    pub fn ring_count(&self) -> Option<&ValueAst> {
        self.ring_membership(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> Option<&ValueAst> {
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

    pub fn get(&self, key: AtomConstraintKey) -> Option<&AtomConstraintAst> {
        self.find(key).ok().map(|i| &self.entries[i])
    }

    /// Insert in sorted order by key, overwrite same key (last-wins).
    pub fn set(&mut self, c: AtomConstraintAst) {
        match self.find(c.key()) {
            Ok(i) => self.entries[i] = c,
            Err(i) => self.entries.insert(i, c),
        }
    }

    /// Transactional write at one key: verify the current value equals `old` (by `canonical_eq`;
    /// both absent is a match), then apply `new` (`Some` sets, `None` removes).
    /// `Err` on a key or old-value mismatch; the store is unchanged when it errors.
    pub fn compare_and_set(
        &mut self,
        old: Option<AtomConstraintAst>,
        new: Option<AtomConstraintAst>,
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

    pub fn remove(&mut self, key: AtomConstraintKey) -> Option<AtomConstraintAst> {
        self.find(key).ok().map(|i| self.entries.remove(i))
    }

    /// `set` each constraint in turn (last-wins), for bulk construction.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = AtomConstraintAst>) {
        for constraint in constraints {
            self.set(constraint);
        }
    }

    /// Overlay `other` onto self by `set`-ing each of its entries (last-wins).
    /// Undetermined entries in `other` remove.
    pub fn update(&mut self, other: &AtomConstraintsAst) {
        for c in other.iter() {
            if c.is_undetermined() {
                self.remove(c.key());
            } else {
                self.set(c.clone());
            }
        }
    }

    /// Bulk-remove entries that don't satisfy the predicate.
    pub fn retain(&mut self, mut f: impl FnMut(&AtomConstraintAst) -> bool) {
        self.entries.retain(|c| f(c));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl ExactSizeIterator<Item = AtomConstraintAst> {
        mem::take(&mut self.entries).into_iter()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = &AtomConstraintAst> {
        self.entries.iter()
    }

    /// No-op: no `AtomConstraintAst` variant carries an entity index.
    pub fn compact(self, _compaction: &IdCompaction) -> Self {
        self
    }
}

impl Canonicalize for AtomConstraintsAst {
    /// Canonicalize each value and drop the vacuous ones. Keys are already unique and
    /// key-sorted (every write goes through `set`), so no dedup or re-sort is needed —
    /// canonicalizing a value never changes its `key()`.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut entries = self
            .entries
            .into_iter()
            .map(Canonicalize::canonicalize)
            .collect::<Result<SmallVec<[AtomConstraintAst; 2]>, _>>()?;
        entries.retain(|c| !c.is_undetermined());
        Ok(Self { entries })
    }
}

impl Lattice for AtomConstraintsAst {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| c.is_ground())
    }

    /// Greatest lower bound as a two-pointer merge over the key-sorted entries: a shared key
    /// meets its two values (`AtomConstraintAst::meet`; a `None` aborts the whole meet), an A-only /
    /// B-only key is kept (meet with the absent ⊤ is the value). Vacuous results are dropped.
    fn meet(&self, other: &Self) -> Option<Self> {
        let mut entries: SmallVec<[AtomConstraintAst; 2]> = SmallVec::new();
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
    /// (`AtomConstraintAst::join`); a single-side key widens to the absent ⊤ and is dropped. A
    /// same-key join never returns `Err(NoJoin)`, but if it did the key would simply drop
    /// (widen to ⊤). The container always has a top (the empty set), so this is total (`Ok`).
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let mut entries: SmallVec<[AtomConstraintAst; 2]> = SmallVec::new();
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
            AtomConstraintAst::Valence(v) => {
                v.matches(target.valence().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraintAst::DonatedPairs(v) => {
                v.matches(target.donated_pairs().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraintAst::AcceptedPairs(v) => {
                v.matches(target.accepted_pairs().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraintAst::AromaticValence(av) => av.matches(
                target
                    .aromatic_valence()
                    .unwrap_or(&AromaticValenceAst::Undetermined),
            ),
            AtomConstraintAst::MulticenterValence(mv) => mv.matches(
                target
                    .multicenter_valence()
                    .unwrap_or(&MulticenterValenceAst::Undetermined),
            ),
            AtomConstraintAst::TetrahedralStereo(ts) => ts.matches(
                target
                    .tetrahedral_stereo()
                    .unwrap_or(&TetrahedralStereoAst::Undetermined),
            ),
            AtomConstraintAst::Degree(v) => {
                v.matches(target.degree().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraintAst::TotalDegree(v) => {
                v.matches(target.total_degree().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraintAst::TotalValence(v) => {
                v.matches(target.total_valence().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraintAst::RingDegree(v) => {
                v.matches(target.ring_degree().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraintAst::RingValence(v) => {
                v.matches(target.ring_valence().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraintAst::TotalHydrogens(v) => {
                v.matches(target.total_hydrogens().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraintAst::RingMembership(rm) => rm.count.matches(
                target
                    .ring_membership(rm.scope)
                    .unwrap_or(&ValueAst::Undetermined),
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

impl FromIterator<AtomConstraintAst> for AtomConstraintsAst {
    fn from_iter<I: IntoIterator<Item = AtomConstraintAst>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}

impl IntoIterator for AtomConstraintsAst {
    type Item = AtomConstraintAst;
    type IntoIter = smallvec::IntoIter<[AtomConstraintAst; 2]>;

    fn into_iter(self) -> Self::IntoIter {
        self.entries.into_iter()
    }
}

impl From<AtomConstraintAst> for AtomConstraintsAst {
    fn from(c: AtomConstraintAst) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<AtomConstraintAst>> for AtomConstraintsAst {
    fn from(cs: Vec<AtomConstraintAst>) -> Self {
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

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticValenceAst {
    #[default]
    Undetermined,
    NotAromatic,
    Aromatic(ValueAst),
}

impl AromaticValenceAst {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn not_aromatic() -> Self {
        Self::NotAromatic
    }

    pub fn aromatic(v: impl Into<ValueAst>) -> Self {
        Self::Aromatic(v.into())
    }

    pub fn is_aromatic(&self) -> bool {
        matches!(self, Self::Aromatic(_))
    }

    pub fn aromatic_increment(&self) -> ValueAst {
        match self {
            Self::Aromatic(v) => match v.as_lit() {
                Some(a) => ValueAst::Lit(aromatic_increment(a)),
                None => ValueAst::Undetermined,
            },
            Self::NotAromatic => ValueAst::Lit(0),
            Self::Undetermined => ValueAst::Undetermined,
        }
    }

    /// Pattern matches value.
    pub fn matches_value(&self, value: i64) -> bool {
        match self {
            Self::Aromatic(v) => v.matches(&ValueAst::Lit(value)),
            Self::NotAromatic => value == 0,
            Self::Undetermined => true,
        }
    }
}

impl From<AromaticValence> for AromaticValenceAst {
    fn from(valence: AromaticValence) -> Self {
        match valence {
            AromaticValence::NotAromatic => Self::NotAromatic,
            AromaticValence::Aromatic(valence) => Self::Aromatic(ValueAst::Lit(valence)),
        }
    }
}

impl Canonicalize for AromaticValenceAst {
    /// Delegate to the inner `ValueAst`; `NotAromatic`/`Undetermined` identity.
    /// No cross-variant fold (`Aromatic(Lit(0))` stays distinct from `NotAromatic`).
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Aromatic(v) => Self::Aromatic(v.canonicalize()?),
            other => other,
        })
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            Self::Aromatic(v) => Ok(match v.canonical()? {
                Cow::Borrowed(_) => Cow::Borrowed(self),
                Cow::Owned(cv) => Cow::Owned(Self::Aromatic(cv)),
            }),
            _ => Ok(Cow::Borrowed(self)),
        }
    }
}

impl Lattice for AromaticValenceAst {
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
        let a = self.canonical().ok()?;
        let b = other.canonical().ok()?;
        match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::NotAromatic, Self::NotAromatic) => Some(Self::NotAromatic),
            (Self::NotAromatic, Self::Aromatic(_)) | (Self::Aromatic(_), Self::NotAromatic) => None,
            (Self::Aromatic(p), Self::Aromatic(q)) => p.meet(q).map(Self::Aromatic),
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
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
    /// `ValueAst::matches` for the `Aromatic` value, never building a `meet`.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, Self::Undetermined | Self::NotAromatic) => true,
            (Self::Undetermined, Self::Aromatic(v)) => v.canonical().is_ok(),
            (Self::NotAromatic, Self::NotAromatic) => true,
            (Self::Aromatic(p), Self::Aromatic(q)) => p.matches(q),
            _ => false,
        }
    }
}

impl AsLit for AromaticValenceAst {
    type Lit = i64;

    /// Inner literal π count when `Aromatic(Lit(n))`; `None` for
    /// `Undetermined` or `Aromatic` wrapping a non-literal.
    #[inline]
    fn as_lit(&self) -> Option<i64> {
        match self {
            Self::Aromatic(v) => v.as_lit(),
            Self::NotAromatic => Some(0),
            _ => None,
        }
    }
}

/// Compute aromatic increment from aromatic valence.
pub fn aromatic_increment(aromatic_valence: i64) -> i64 {
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

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MulticenterValenceAst {
    #[default]
    Undetermined,
    NotMulticenter,
    Multicenter(ValueAst),
}

impl MulticenterValenceAst {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn not_multicenter() -> Self {
        Self::NotMulticenter
    }

    pub fn multicenter(v: impl Into<ValueAst>) -> Self {
        Self::Multicenter(v.into())
    }

    pub fn is_multicenter(&self) -> bool {
        matches!(self, Self::Multicenter(_))
    }

    /// Pattern matches value.
    pub fn matches_value(&self, value: i64) -> bool {
        match self {
            Self::Multicenter(v) => v.matches(&ValueAst::Lit(value)),
            Self::NotMulticenter => value == 0,
            Self::Undetermined => true,
        }
    }
}

impl From<MulticenterValence> for MulticenterValenceAst {
    fn from(valence: MulticenterValence) -> Self {
        match valence {
            MulticenterValence::NotMulticenter => Self::NotMulticenter,
            MulticenterValence::Multicenter(valence) => Self::Multicenter(ValueAst::Lit(valence)),
        }
    }
}

impl Canonicalize for MulticenterValenceAst {
    /// Delegate to the inner `ValueAst`; `NotMulticenter`/`Undetermined` identity.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Multicenter(v) => Self::Multicenter(v.canonicalize()?),
            other => other,
        })
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            Self::Multicenter(v) => Ok(match v.canonical()? {
                Cow::Borrowed(_) => Cow::Borrowed(self),
                Cow::Owned(cv) => Cow::Owned(Self::Multicenter(cv)),
            }),
            _ => Ok(Cow::Borrowed(self)),
        }
    }
}

impl Lattice for MulticenterValenceAst {
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
        let a = self.canonical().ok()?;
        let b = other.canonical().ok()?;
        match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::NotMulticenter, Self::NotMulticenter) => Some(Self::NotMulticenter),
            (Self::NotMulticenter, Self::Multicenter(_))
            | (Self::Multicenter(_), Self::NotMulticenter) => None,
            (Self::Multicenter(p), Self::Multicenter(q)) => p.meet(q).map(Self::Multicenter),
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        Ok(match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::NotMulticenter, Self::NotMulticenter) => Self::NotMulticenter,
            (Self::NotMulticenter, Self::Multicenter(_))
            | (Self::Multicenter(_), Self::NotMulticenter) => Self::Undetermined,
            (Self::Multicenter(p), Self::Multicenter(q)) => Self::Multicenter(p.join(q)?),
        })
    }

    /// Partial-order check `target ⊑ self`, allocation-free — defers to the inner
    /// `ValueAst::matches` for the `Multicenter` value, never building a `meet`.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, Self::Undetermined | Self::NotMulticenter) => true,
            (Self::Undetermined, Self::Multicenter(v)) => v.canonical().is_ok(),
            (Self::NotMulticenter, Self::NotMulticenter) => true,
            (Self::Multicenter(p), Self::Multicenter(q)) => p.matches(q),
            _ => false,
        }
    }
}

impl AsLit for MulticenterValenceAst {
    type Lit = i64;

    /// Inner literal multicenter valence when `Multicenter(Lit(n))`; `None`
    /// for `Undetermined` or non-literal inner.
    #[inline]
    fn as_lit(&self) -> Option<i64> {
        match self {
            Self::Multicenter(v) => v.as_lit(),
            Self::NotMulticenter => Some(0),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Compaction;

    use super::*;
    use crate::ast::value::ValueTerm;

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraintAst::valence(4), AtomConstraintAst::Valence(ValueAst::Lit(4)))]
    #[case::total_valence(AtomConstraintAst::total_valence(5), AtomConstraintAst::TotalValence(ValueAst::Lit(5)))]
    #[case::donated_pairs(AtomConstraintAst::donated_pairs(1), AtomConstraintAst::DonatedPairs(ValueAst::Lit(1)))]
    #[case::accepted_pairs(AtomConstraintAst::accepted_pairs(2), AtomConstraintAst::AcceptedPairs(ValueAst::Lit(2)))]
    #[case::degree(AtomConstraintAst::degree(3), AtomConstraintAst::Degree(ValueAst::Lit(3)))]
    #[case::total_degree(AtomConstraintAst::total_degree(4), AtomConstraintAst::TotalDegree(ValueAst::Lit(4)))]
    #[case::ring_degree(AtomConstraintAst::ring_degree(2), AtomConstraintAst::RingDegree(ValueAst::Lit(2)))]
    #[case::ring_valence(AtomConstraintAst::ring_valence(3), AtomConstraintAst::RingValence(ValueAst::Lit(3)))]
    #[case::total_hydrogens(AtomConstraintAst::total_hydrogens(3), AtomConstraintAst::TotalHydrogens(ValueAst::Lit(3)))]
    #[case::ring_membership_all(AtomConstraintAst::ring_membership(RingScope::All, 1), AtomConstraintAst::RingMembership(RingMembershipAst { scope: RingScope::All, count: ValueAst::Lit(1) }))]
    #[case::ring_membership_size(AtomConstraintAst::ring_membership(RingScope::Size(6), 1), AtomConstraintAst::RingMembership(RingMembershipAst { scope: RingScope::Size(6), count: ValueAst::Lit(1) }))]
    #[case::aromatic_valence(
        AtomConstraintAst::aromatic_valence(AromaticValenceAst::NotAromatic),
        AtomConstraintAst::AromaticValence(AromaticValenceAst::NotAromatic),
    )]
    #[case::multicenter_valence(
        AtomConstraintAst::multicenter_valence(MulticenterValenceAst::NotMulticenter),
        AtomConstraintAst::MulticenterValence(MulticenterValenceAst::NotMulticenter),
    )]
    #[case::tetrahedral_stereo(
        AtomConstraintAst::tetrahedral_stereo(TetrahedralStereoAst::NotStereo),
        AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo),
    )]
    fn test_atom_constraint_ast_constructors(
        #[case] actual: AtomConstraintAst,
        #[case] expected: AtomConstraintAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraintAst::valence(4), AtomConstraintKey::Valence)]
    #[case::total_valence(AtomConstraintAst::total_valence(5), AtomConstraintKey::TotalValence)]
    #[case::aromatic_valence(AtomConstraintAst::aromatic_valence(AromaticValenceAst::NotAromatic), AtomConstraintKey::AromaticValence)]
    #[case::multicenter_valence(AtomConstraintAst::multicenter_valence(MulticenterValenceAst::Undetermined), AtomConstraintKey::MulticenterValence)]
    #[case::donated_pairs(AtomConstraintAst::donated_pairs(1), AtomConstraintKey::DonatedPairs)]
    #[case::accepted_pairs(AtomConstraintAst::accepted_pairs(2), AtomConstraintKey::AcceptedPairs)]
    #[case::degree(AtomConstraintAst::degree(3), AtomConstraintKey::Degree)]
    #[case::total_degree(AtomConstraintAst::total_degree(4), AtomConstraintKey::TotalDegree)]
    #[case::ring_degree(AtomConstraintAst::ring_degree(2), AtomConstraintKey::RingDegree)]
    #[case::ring_valence(AtomConstraintAst::ring_valence(3), AtomConstraintKey::RingValence)]
    #[case::total_hydrogens(AtomConstraintAst::total_hydrogens(3), AtomConstraintKey::TotalHydrogens)]
    #[case::ring_membership_all(AtomConstraintAst::ring_membership(RingScope::All, 1), AtomConstraintKey::RingMembership(RingScope::All))]
    #[case::ring_membership_size(AtomConstraintAst::ring_membership(RingScope::Size(6), 1), AtomConstraintKey::RingMembership(RingScope::Size(6)))]
    #[case::tetrahedral_stereo(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo), AtomConstraintKey::TetrahedralStereo)]
    fn test_atom_constraint_ast_key(
        #[case] constraint: AtomConstraintAst,
        #[case] expected: AtomConstraintKey,
    ) {
        assert_eq!(constraint.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_litset_singleton(AtomConstraintAst::Valence(ValueAst::lit_set([4])), Ok(AtomConstraintAst::valence(4)))]
    #[case::ring_count_litset_singleton(
        AtomConstraintAst::RingMembership(RingMembershipAst::new(RingScope::Size(6), ValueAst::lit_set([2]))),
        Ok(AtomConstraintAst::ring_membership(RingScope::Size(6), 2)))]
    #[case::empty_litset_contradiction(AtomConstraintAst::Valence(ValueAst::lit_set(Vec::<i64>::new())), Err(Contradiction))]
    fn test_atom_constraint_ast_canonicalize(
        #[case] constraint: AtomConstraintAst,
        #[case] expected: Result<AtomConstraintAst, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_lit(AtomConstraintAst::valence(4), false)]
    #[case::valence_undetermined(AtomConstraintAst::Valence(ValueAst::Undetermined), true)]
    #[case::degree_undetermined(AtomConstraintAst::Degree(ValueAst::Undetermined), true)]
    #[case::ring_membership_undetermined(AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined), true)]
    #[case::aromatic_undetermined(AtomConstraintAst::aromatic_valence(AromaticValenceAst::Undetermined), true)]
    #[case::aromatic_not_aromatic(AtomConstraintAst::aromatic_valence(AromaticValenceAst::NotAromatic), false)]
    #[case::aromatic_with_value(AtomConstraintAst::aromatic_valence(AromaticValenceAst::aromatic(1)), false)]
    #[case::multicenter_undetermined(AtomConstraintAst::multicenter_valence(MulticenterValenceAst::Undetermined), true)]
    #[case::multicenter_not(AtomConstraintAst::multicenter_valence(MulticenterValenceAst::NotMulticenter), false)]
    #[case::multicenter_with_value(AtomConstraintAst::multicenter_valence(MulticenterValenceAst::multicenter(1)), false)]
    #[case::tetrahedral_not_stereo(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo), false)]
    #[case::tetrahedral_undetermined(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Undetermined), true)]
    fn test_atom_constraint_ast_is_undetermined(
        #[case] c: AtomConstraintAst,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraintAst::valence(4), AtomConstraintAst::Valence(ValueAst::Undetermined))]
    #[case::degree(AtomConstraintAst::Degree(ValueAst::Lit(2)), AtomConstraintAst::Degree(ValueAst::Undetermined))]
    #[case::ring_membership_keeps_scope(AtomConstraintAst::ring_membership(RingScope::Size(6), 1), AtomConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined))]
    #[case::aromatic(AtomConstraintAst::aromatic_valence(AromaticValenceAst::aromatic(1)), AtomConstraintAst::aromatic_valence(AromaticValenceAst::Undetermined))]
    #[case::multicenter(AtomConstraintAst::multicenter_valence(MulticenterValenceAst::multicenter(1)), AtomConstraintAst::multicenter_valence(MulticenterValenceAst::Undetermined))]
    #[case::tetrahedral(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo), AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Undetermined))]
    fn test_atom_constraint_ast_as_undetermined(
        #[case] c: AtomConstraintAst,
        #[case] expected: AtomConstraintAst,
    ) {
        assert_eq!(c.as_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_key_widens(AtomConstraintAst::valence(4), AtomConstraintAst::valence(3), Ok(AtomConstraintAst::Valence(ValueAst::lit_set([3, 4]))))]
    #[case::different_key(AtomConstraintAst::valence(4), AtomConstraintAst::degree(3), Err(NoJoin))]
    fn test_atom_constraint_ast_join(
        #[case] a: AtomConstraintAst,
        #[case] b: AtomConstraintAst,
        #[case] expected: Result<AtomConstraintAst, NoJoin>,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined, AromaticValenceAst::Undetermined)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, AromaticValenceAst::NotAromatic)]
    #[case::aromatic(AromaticValenceAst::aromatic(1), AromaticValenceAst::Aromatic(ValueAst::Lit(1)))]
    fn test_aromatic_valence_ast_constructors(
        #[case] actual: AromaticValenceAst,
        #[case] expected: AromaticValenceAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::not_aromatic(AromaticValence::NotAromatic, AromaticValenceAst::NotAromatic)]
    #[case::aromatic(
        AromaticValence::Aromatic(1),
        AromaticValenceAst::Aromatic(ValueAst::Lit(1))
    )]
    fn test_aromatic_valence_ast_from(
        #[case] valence: AromaticValence,
        #[case] expected: AromaticValenceAst,
    ) {
        assert_eq!(AromaticValenceAst::from(valence), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined, false)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, false)]
    #[case::aromatic_undetermined(AromaticValenceAst::Aromatic(ValueAst::Undetermined), true)]
    #[case::aromatic_lit(AromaticValenceAst::aromatic(1), true)]
    fn test_aromatic_valence_ast_is_aromatic(
        #[case] v: AromaticValenceAst,
        #[case] expected: bool,
    ) {
        assert_eq!(v.is_aromatic(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined, ValueAst::Undetermined)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, ValueAst::Lit(0))]
    #[case::aromatic_undetermined(
        AromaticValenceAst::Aromatic(ValueAst::Undetermined),
        ValueAst::Undetermined
    )]
    #[case::aromatic_one(AromaticValenceAst::aromatic(1), ValueAst::Lit(1))]
    #[case::aromatic_zero(AromaticValenceAst::aromatic(0), ValueAst::Lit(0))]
    #[case::aromatic_two(AromaticValenceAst::aromatic(2), ValueAst::Lit(0))]
    fn test_aromatic_valence_ast_aromatic_increment(
        #[case] v: AromaticValenceAst,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(v.aromatic_increment(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined, true)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, false)]
    #[case::aromatic_lit(AromaticValenceAst::aromatic(1), false)]
    #[case::aromatic_inner_undetermined(
        AromaticValenceAst::Aromatic(ValueAst::Undetermined),
        false
    )]
    fn test_aromatic_valence_ast_is_undetermined(
        #[case] v: AromaticValenceAst,
        #[case] expected: bool,
    ) {
        assert_eq!(v.is_undetermined(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined, None)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, Some(0))]
    #[case::aromatic_undetermined(AromaticValenceAst::Aromatic(ValueAst::Undetermined), None)]
    #[case::aromatic_lit(AromaticValenceAst::aromatic(3), Some(3))]
    #[case::aromatic_term_unresolved(
        AromaticValenceAst::Aromatic(ValueAst::term(ValueTerm::Lit(2))),
        None
    )]
    fn test_aromatic_valence_ast_as_lit(
        #[case] v: AromaticValenceAst,
        #[case] expected: Option<i64>,
    ) {
        assert_eq!(v.as_lit(), expected);
        assert_eq!(v.is_ground(), expected.is_some());
    }

    #[rstest]
    #[case::aromatic_folds_inner(
        AromaticValenceAst::Aromatic(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Lit(1), ValueTerm::Lit(1)]))),
        Ok(AromaticValenceAst::aromatic(2)),
    )]
    #[case::aromatic_zero_not_collapsed(
        AromaticValenceAst::aromatic(0),
        Ok(AromaticValenceAst::aromatic(0))
    )]
    fn test_aromatic_valence_ast_canonicalize(
        #[case] input: AromaticValenceAst,
        #[case] expected: Result<AromaticValenceAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic)]
    #[case::aromatic_lit(AromaticValenceAst::aromatic(1))]
    #[case::aromatic_zero(AromaticValenceAst::aromatic(0))]
    fn test_aromatic_valence_ast_canonicalize_identity(#[case] input: AromaticValenceAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_aromatic(AromaticValenceAst::Undetermined, AromaticValenceAst::aromatic(1), Some(AromaticValenceAst::aromatic(1)))]
    #[case::not_not(AromaticValenceAst::NotAromatic, AromaticValenceAst::NotAromatic, Some(AromaticValenceAst::NotAromatic))]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, AromaticValenceAst::aromatic(1), None)]
    #[case::aromatic_eq(AromaticValenceAst::aromatic(1), AromaticValenceAst::aromatic(1), Some(AromaticValenceAst::aromatic(1)))]
    #[case::aromatic_neq(AromaticValenceAst::aromatic(1), AromaticValenceAst::aromatic(2), None)]
    #[case::aromatic_inner_wildcard(AromaticValenceAst::Aromatic(ValueAst::Undetermined), AromaticValenceAst::aromatic(2), Some(AromaticValenceAst::aromatic(2)))]
    fn test_aromatic_valence_ast_meet(
        #[case] a: AromaticValenceAst,
        #[case] b: AromaticValenceAst,
        #[case] expected: Option<AromaticValenceAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und(AromaticValenceAst::Undetermined, AromaticValenceAst::aromatic(1), AromaticValenceAst::Undetermined)]
    #[case::not_not(AromaticValenceAst::NotAromatic, AromaticValenceAst::NotAromatic, AromaticValenceAst::NotAromatic)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, AromaticValenceAst::aromatic(1), AromaticValenceAst::Undetermined)]
    #[case::aromatic_eq(AromaticValenceAst::aromatic(1), AromaticValenceAst::aromatic(1), AromaticValenceAst::aromatic(1))]
    #[case::aromatic_inner_wildcard(AromaticValenceAst::Aromatic(ValueAst::Undetermined), AromaticValenceAst::aromatic(1), AromaticValenceAst::Aromatic(ValueAst::Undetermined))]
    fn test_aromatic_valence_ast_join(
        #[case] a: AromaticValenceAst,
        #[case] b: AromaticValenceAst,
        #[case] expected: AromaticValenceAst,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rstest]
    #[case::wildcard_vs_not_aromatic(
        AromaticValenceAst::Undetermined,
        AromaticValenceAst::NotAromatic,
        true
    )]
    #[case::wildcard_vs_aromatic_lit(
        AromaticValenceAst::Undetermined,
        AromaticValenceAst::aromatic(1),
        true
    )]
    #[case::not_aromatic_vs_aromatic(
        AromaticValenceAst::NotAromatic,
        AromaticValenceAst::aromatic(1),
        false
    )]
    #[case::aromatic_vs_not_aromatic(
        AromaticValenceAst::aromatic(1),
        AromaticValenceAst::NotAromatic,
        false
    )]
    #[case::not_aromatic_vs_not_aromatic(
        AromaticValenceAst::NotAromatic,
        AromaticValenceAst::NotAromatic,
        true
    )]
    #[case::aromatic_lit_match(
        AromaticValenceAst::aromatic(1),
        AromaticValenceAst::aromatic(1),
        true
    )]
    #[case::aromatic_lit_mismatch(
        AromaticValenceAst::aromatic(1),
        AromaticValenceAst::aromatic(2),
        false
    )]
    #[case::aromatic_wildcard_inner(
        AromaticValenceAst::Aromatic(ValueAst::Undetermined),
        AromaticValenceAst::aromatic(2),
        true
    )]
    #[case::specific_vs_undetermined(
        AromaticValenceAst::aromatic(1),
        AromaticValenceAst::Undetermined,
        false
    )]
    fn test_aromatic_valence_ast_matches(
        #[case] pattern: AromaticValenceAst,
        #[case] target: AromaticValenceAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined, MulticenterValenceAst::Undetermined)]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, MulticenterValenceAst::NotMulticenter)]
    #[case::multicenter(MulticenterValenceAst::multicenter(2), MulticenterValenceAst::Multicenter(ValueAst::Lit(2)))]
    fn test_multicenter_valence_ast_constructors(
        #[case] actual: MulticenterValenceAst,
        #[case] expected: MulticenterValenceAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::not_multicenter(
        MulticenterValence::NotMulticenter,
        MulticenterValenceAst::NotMulticenter
    )]
    #[case::multicenter(
        MulticenterValence::Multicenter(2),
        MulticenterValenceAst::Multicenter(ValueAst::Lit(2))
    )]
    fn test_multicenter_valence_ast_from(
        #[case] valence: MulticenterValence,
        #[case] expected: MulticenterValenceAst,
    ) {
        assert_eq!(MulticenterValenceAst::from(valence), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined, false)]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, false)]
    #[case::multicenter_undetermined(MulticenterValenceAst::Multicenter(ValueAst::Undetermined), true)]
    #[case::multicenter_lit(MulticenterValenceAst::multicenter(1), true)]
    fn test_multicenter_valence_ast_is_multicenter(
        #[case] v: MulticenterValenceAst,
        #[case] expected: bool,
    ) {
        assert_eq!(v.is_multicenter(), expected);
    }

    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined, true)]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, false)]
    #[case::multicenter_lit(MulticenterValenceAst::multicenter(1), false)]
    fn test_multicenter_valence_ast_is_undetermined(
        #[case] v: MulticenterValenceAst,
        #[case] expected: bool,
    ) {
        assert_eq!(v.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined, None)]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, Some(0))]
    #[case::multicenter_undetermined(MulticenterValenceAst::Multicenter(ValueAst::Undetermined), None)]
    #[case::multicenter_lit(MulticenterValenceAst::multicenter(2), Some(2))]
    #[case::multicenter_term_unresolved(MulticenterValenceAst::Multicenter(ValueAst::term(ValueTerm::Lit(3))), None)]
    fn test_multicenter_valence_ast_as_lit(
        #[case] v: MulticenterValenceAst,
        #[case] expected: Option<i64>,
    ) {
        assert_eq!(v.as_lit(), expected);
        assert_eq!(v.is_ground(), expected.is_some());
    }

    #[rstest]
    #[case::multicenter_folds_inner(
        MulticenterValenceAst::Multicenter(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Lit(1), ValueTerm::Lit(2)]))),
        Ok(MulticenterValenceAst::multicenter(3)),
    )]
    #[case::multicenter_zero_not_collapsed(
        MulticenterValenceAst::multicenter(0),
        Ok(MulticenterValenceAst::multicenter(0))
    )]
    fn test_multicenter_valence_ast_canonicalize(
        #[case] input: MulticenterValenceAst,
        #[case] expected: Result<MulticenterValenceAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined)]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter)]
    #[case::multicenter_lit(MulticenterValenceAst::multicenter(1))]
    #[case::multicenter_zero(MulticenterValenceAst::multicenter(0))]
    fn test_multicenter_valence_ast_canonicalize_identity(#[case] input: MulticenterValenceAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_multicenter(MulticenterValenceAst::Undetermined, MulticenterValenceAst::multicenter(2), Some(MulticenterValenceAst::multicenter(2)))]
    #[case::not_not(MulticenterValenceAst::NotMulticenter, MulticenterValenceAst::NotMulticenter, Some(MulticenterValenceAst::NotMulticenter))]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, MulticenterValenceAst::multicenter(2), None)]
    #[case::multicenter_eq(MulticenterValenceAst::multicenter(2), MulticenterValenceAst::multicenter(2), Some(MulticenterValenceAst::multicenter(2)))]
    #[case::multicenter_neq(MulticenterValenceAst::multicenter(2), MulticenterValenceAst::multicenter(3), None)]
    #[case::multicenter_inner_wildcard(MulticenterValenceAst::Multicenter(ValueAst::Undetermined), MulticenterValenceAst::multicenter(3), Some(MulticenterValenceAst::multicenter(3)))]
    fn test_multicenter_valence_ast_meet(
        #[case] a: MulticenterValenceAst,
        #[case] b: MulticenterValenceAst,
        #[case] expected: Option<MulticenterValenceAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und(MulticenterValenceAst::Undetermined, MulticenterValenceAst::multicenter(2), MulticenterValenceAst::Undetermined)]
    #[case::not_not(MulticenterValenceAst::NotMulticenter, MulticenterValenceAst::NotMulticenter, MulticenterValenceAst::NotMulticenter)]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, MulticenterValenceAst::multicenter(2), MulticenterValenceAst::Undetermined)]
    #[case::multicenter_eq(MulticenterValenceAst::multicenter(2), MulticenterValenceAst::multicenter(2), MulticenterValenceAst::multicenter(2))]
    #[case::multicenter_inner_wildcard(MulticenterValenceAst::Multicenter(ValueAst::Undetermined), MulticenterValenceAst::multicenter(2), MulticenterValenceAst::Multicenter(ValueAst::Undetermined))]
    fn test_multicenter_valence_ast_join(
        #[case] a: MulticenterValenceAst,
        #[case] b: MulticenterValenceAst,
        #[case] expected: MulticenterValenceAst,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard_vs_not_multicenter(MulticenterValenceAst::Undetermined, MulticenterValenceAst::NotMulticenter, true)]
    #[case::wildcard_vs_multicenter_lit(MulticenterValenceAst::Undetermined, MulticenterValenceAst::multicenter(2), true)]
    #[case::not_multicenter_vs_multicenter(MulticenterValenceAst::NotMulticenter, MulticenterValenceAst::multicenter(2), false)]
    #[case::multicenter_vs_not_multicenter(MulticenterValenceAst::multicenter(2), MulticenterValenceAst::NotMulticenter, false)]
    #[case::not_multicenter_vs_not_multicenter(MulticenterValenceAst::NotMulticenter, MulticenterValenceAst::NotMulticenter, true)]
    #[case::multicenter_lit_match(MulticenterValenceAst::multicenter(2), MulticenterValenceAst::multicenter(2), true)]
    #[case::multicenter_lit_mismatch(MulticenterValenceAst::multicenter(2), MulticenterValenceAst::multicenter(3), false)]
    #[case::specific_vs_undetermined(MulticenterValenceAst::multicenter(2), MulticenterValenceAst::Undetermined, false)]
    fn test_multicenter_valence_ast_matches(
        #[case] pattern: MulticenterValenceAst,
        #[case] target: MulticenterValenceAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_atom_constraints_ast_new() {
        let cs = AtomConstraintsAst::new();
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
    fn test_atom_constraints_ast_contains(
        #[case] key: AtomConstraintKey,
        #[case] expected: bool,
    ) {
        let cs = AtomConstraintsAst::from_iter([
            AtomConstraintAst::valence(4),
            AtomConstraintAst::ring_membership(RingScope::All, 2),
            AtomConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_all(AtomConstraintKey::RingMembership(RingScope::All), Some(AtomConstraintAst::ring_membership(RingScope::All, 2)))]
    #[case::ring_size(AtomConstraintKey::RingMembership(RingScope::Size(6)), Some(AtomConstraintAst::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_size_absent(AtomConstraintKey::RingMembership(RingScope::Size(5)), None)]
    #[case::valence(AtomConstraintKey::Valence, Some(AtomConstraintAst::valence(4)))]
    fn test_atom_constraints_ast_get(
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintAst>,
    ) {
        let cs = AtomConstraintsAst::from_iter([
            AtomConstraintAst::valence(4),
            AtomConstraintAst::ring_membership(RingScope::All, 2),
            AtomConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get(key), expected.as_ref());
    }

    #[rstest]
    fn test_atom_constraints_ast_remove() {
        let mut cs = AtomConstraintsAst::from_iter([
            AtomConstraintAst::valence(4),
            AtomConstraintAst::ring_membership(RingScope::All, 2),
            AtomConstraintAst::ring_membership(RingScope::Size(6), 1),
        ]);
        let removed = cs.remove(AtomConstraintKey::RingMembership(RingScope::Size(6)));
        assert_eq!(
            removed,
            Some(AtomConstraintAst::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(
            cs.iter().cloned().collect::<Vec<_>>(),
            vec![
                AtomConstraintAst::valence(4),
                AtomConstraintAst::ring_membership(RingScope::All, 2),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::drop_vacuous(
        AtomConstraintsAst::from_iter([
            AtomConstraintAst::Valence(ValueAst::Undetermined),
            AtomConstraintAst::degree(3),
        ]),
        Ok(AtomConstraintsAst::from_iter([AtomConstraintAst::degree(3)])))]
    #[case::canonicalizes_values(
        AtomConstraintsAst::from_iter([
            AtomConstraintAst::Degree(ValueAst::lit_set([3])),
            AtomConstraintAst::Valence(ValueAst::lit_set([4])),
        ]),
        Ok(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4), AtomConstraintAst::degree(3)])))]
    fn test_atom_constraints_ast_canonicalize(
        #[case] constraints: AtomConstraintsAst,
        #[case] expected: Result<AtomConstraintsAst, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![AtomConstraintAst::valence(4)], vec![AtomConstraintAst::valence(4)])]
    #[case::overwrite_same_key(vec![AtomConstraintAst::valence(3), AtomConstraintAst::valence(4)], vec![AtomConstraintAst::valence(4)])]
    #[case::vacuous_overwrites(vec![AtomConstraintAst::valence(4), AtomConstraintAst::Valence(ValueAst::Undetermined)], vec![AtomConstraintAst::Valence(ValueAst::Undetermined)])]
    #[case::vacuous_absent_inserts(vec![AtomConstraintAst::Valence(ValueAst::Undetermined)], vec![AtomConstraintAst::Valence(ValueAst::Undetermined)])]
    #[case::new_key_sorts(vec![AtomConstraintAst::degree(3), AtomConstraintAst::valence(4)], vec![AtomConstraintAst::valence(4), AtomConstraintAst::degree(3)])]
    #[case::ring_overwrite_scope(vec![AtomConstraintAst::ring_membership(RingScope::Size(6), 1), AtomConstraintAst::ring_membership(RingScope::Size(6), 2)], vec![AtomConstraintAst::ring_membership(RingScope::Size(6), 2)])]
    #[case::ring_vacuous_overwrites_scope(vec![AtomConstraintAst::ring_membership(RingScope::All, 2), AtomConstraintAst::ring_membership(RingScope::Size(6), 1), AtomConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined)], vec![AtomConstraintAst::ring_membership(RingScope::All, 2), AtomConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined)])]
    fn test_atom_constraints_ast_set(
        #[case] sequence: Vec<AtomConstraintAst>,
        #[case] expected_state: Vec<AtomConstraintAst>,
    ) {
        let mut cs = AtomConstraintsAst::new();
        for c in sequence {
            cs.set(c);
        }
        assert_eq!(cs, AtomConstraintsAst::from_iter(expected_state));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite_shared(
        vec![AtomConstraintAst::valence(3), AtomConstraintAst::degree(3)],
        vec![AtomConstraintAst::valence(4)],
        vec![AtomConstraintAst::valence(4), AtomConstraintAst::degree(3)])]
    #[case::keeps_disjoint(
        vec![AtomConstraintAst::valence(4)],
        vec![AtomConstraintAst::degree(3)],
        vec![AtomConstraintAst::valence(4), AtomConstraintAst::degree(3)])]
    #[case::vacuous_removes(
        vec![AtomConstraintAst::valence(4), AtomConstraintAst::degree(3)],
        vec![AtomConstraintAst::Valence(ValueAst::Undetermined)],
        vec![AtomConstraintAst::degree(3)])]
    #[case::empty_other(
        vec![AtomConstraintAst::valence(4)],
        vec![],
        vec![AtomConstraintAst::valence(4)])]
    #[case::ring_overwrite_scope(
        vec![AtomConstraintAst::ring_membership(RingScope::All, 2), AtomConstraintAst::ring_membership(RingScope::Size(6), 1)],
        vec![AtomConstraintAst::ring_membership(RingScope::Size(6), 3)],
        vec![AtomConstraintAst::ring_membership(RingScope::All, 2), AtomConstraintAst::ring_membership(RingScope::Size(6), 3)])]
    fn test_atom_constraints_ast_update(
        #[case] initial: Vec<AtomConstraintAst>,
        #[case] other: Vec<AtomConstraintAst>,
        #[case] expected: Vec<AtomConstraintAst>,
    ) {
        let mut cs = AtomConstraintsAst::from_iter(initial);
        let overlay = AtomConstraintsAst::from_iter(other);
        cs.update(&overlay);
        assert_eq!(cs, AtomConstraintsAst::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(vec![AtomConstraintAst::valence(3)], Some(AtomConstraintAst::valence(3)), Some(AtomConstraintAst::valence(4)), Ok(()), vec![AtomConstraintAst::valence(4)])]
    #[case::remove(vec![AtomConstraintAst::valence(4)], Some(AtomConstraintAst::valence(4)), None, Ok(()), vec![])]
    #[case::add_from_absent(vec![], None, Some(AtomConstraintAst::valence(4)), Ok(()), vec![AtomConstraintAst::valence(4)])]
    #[case::canonical_match(vec![AtomConstraintAst::Valence(ValueAst::lit_set([4]))], Some(AtomConstraintAst::valence(4)), Some(AtomConstraintAst::valence(5)), Ok(()), vec![AtomConstraintAst::valence(5)])]
    #[case::old_mismatch(vec![AtomConstraintAst::valence(3)], Some(AtomConstraintAst::valence(4)), Some(AtomConstraintAst::valence(5)), Err(Contradiction), vec![AtomConstraintAst::valence(3)])]
    #[case::old_absent_mismatch(vec![AtomConstraintAst::valence(3)], None, Some(AtomConstraintAst::valence(4)), Err(Contradiction), vec![AtomConstraintAst::valence(3)])]
    #[case::key_mismatch(vec![], Some(AtomConstraintAst::valence(3)), Some(AtomConstraintAst::degree(4)), Err(Contradiction), vec![])]
    #[case::noop(vec![AtomConstraintAst::valence(4)], None, None, Ok(()), vec![AtomConstraintAst::valence(4)])]
    fn test_atom_constraints_ast_compare_and_set(
        #[case] initial: Vec<AtomConstraintAst>,
        #[case] old: Option<AtomConstraintAst>,
        #[case] new: Option<AtomConstraintAst>,
        #[case] expected_result: Result<(), Contradiction>,
        #[case] expected_state: Vec<AtomConstraintAst>,
    ) {
        let mut cs = AtomConstraintsAst::from_iter(initial);
        assert_eq!(cs.compare_and_set(old, new), expected_result);
        assert_eq!(cs, AtomConstraintsAst::from_iter(expected_state));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &AtomConstraintAst| matches!(c, AtomConstraintAst::Valence(_) | AtomConstraintAst::RingMembership(_)), vec![AtomConstraintAst::valence(4), AtomConstraintAst::ring_membership(RingScope::All, 2)])]
    #[case::all_dropped(|_: &AtomConstraintAst| false, vec![])]
    fn test_atom_constraints_ast_retain(
        #[case] predicate: impl FnMut(&AtomConstraintAst) -> bool,
        #[case] expected: Vec<AtomConstraintAst>,
    ) {
        let mut cs = AtomConstraintsAst::from_iter([
            AtomConstraintAst::valence(4),
            AtomConstraintAst::degree(3),
            AtomConstraintAst::ring_membership(RingScope::All, 2),
        ]);
        cs.retain(predicate);
        assert_eq!(cs, AtomConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_atom_constraints_ast_clear() {
        let mut cs = AtomConstraintsAst::from_iter([
            AtomConstraintAst::valence(4),
            AtomConstraintAst::degree(3),
        ]);
        cs.clear();
        assert_eq!(cs, AtomConstraintsAst::new());
    }

    #[rstest]
    fn test_atom_constraints_ast_take() {
        let mut empty = AtomConstraintsAst::new();
        let mut empty_taken = empty.take();
        assert_eq!(empty_taken.len(), 0);
        assert_eq!(empty_taken.size_hint(), (0, Some(0)));
        assert_eq!(empty_taken.next(), None);

        let mut cs = AtomConstraintsAst::from_iter([
            AtomConstraintAst::valence(4),
            AtomConstraintAst::degree(3),
        ]);
        let mut taken = cs.take();
        assert_eq!(taken.len(), 2);
        assert_eq!(taken.size_hint(), (2, Some(2)));
        assert_eq!(taken.next(), Some(AtomConstraintAst::valence(4)));
        assert_eq!(taken.len(), 1);
        assert_eq!(taken.size_hint(), (1, Some(1)));
        assert_eq!(taken.next(), Some(AtomConstraintAst::degree(3)));
        assert_eq!(taken.len(), 0);
        assert_eq!(taken.next(), None);
        drop(taken);
        assert_eq!(cs, AtomConstraintsAst::new());
    }

    #[rstest]
    fn test_atom_constraints_ast_iter() {
        let empty = AtomConstraintsAst::new();
        let mut empty_iter = empty.iter();
        assert_eq!(empty_iter.len(), 0);
        assert_eq!(empty_iter.size_hint(), (0, Some(0)));
        assert_eq!(empty_iter.next(), None);

        let cs = AtomConstraintsAst::from_iter([
            AtomConstraintAst::ring_membership(RingScope::Size(6), 1),
            AtomConstraintAst::valence(4),
            AtomConstraintAst::degree(3),
        ]);
        let mut iter = cs.iter();
        assert_eq!(iter.len(), 3);
        assert_eq!(iter.size_hint(), (3, Some(3)));
        assert_eq!(iter.next(), Some(&AtomConstraintAst::valence(4)));
        assert_eq!(iter.len(), 2);
        assert_eq!(iter.size_hint(), (2, Some(2)));
        assert_eq!(iter.next(), Some(&AtomConstraintAst::degree(3)));
        assert_eq!(
            iter.next(),
            Some(&AtomConstraintAst::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(iter.len(), 0);
        assert_eq!(iter.next(), None);
    }

    #[rstest]
    fn test_atom_constraints_ast_compact() {
        let cs = AtomConstraintsAst::from_iter([
            AtomConstraintAst::valence(4),
            AtomConstraintAst::degree(3),
        ]);
        let compaction = IdCompaction::new(
            Compaction::new(vec![0, 1, 2], vec![0]),
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
    #[case::distinct(vec![AtomConstraintAst::valence(4), AtomConstraintAst::degree(3)], vec![AtomConstraintAst::valence(4), AtomConstraintAst::degree(3)])]
    #[case::same_kind_last_wins(vec![AtomConstraintAst::valence(3), AtomConstraintAst::valence(4)], vec![AtomConstraintAst::valence(4)])]
    #[case::empty(vec![], vec![])]
    fn test_atom_constraints_ast_from_iter(
        #[case] input: Vec<AtomConstraintAst>,
        #[case] expected: Vec<AtomConstraintAst>,
    ) {
        let cs = AtomConstraintsAst::from_iter(input);
        assert_eq!(cs, AtomConstraintsAst::from_iter(expected));
    }

    #[rstest]
    fn test_atom_constraints_ast_into_iter() {
        let cs = AtomConstraintsAst::from_iter([
            AtomConstraintAst::valence(4),
            AtomConstraintAst::degree(3),
        ]);
        let collected: Vec<AtomConstraintAst> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![AtomConstraintAst::valence(4), AtomConstraintAst::degree(3)],
        );
    }

    #[rstest]
    fn test_atom_constraints_ast_from_atom_constraint() {
        let cs: AtomConstraintsAst = AtomConstraintAst::valence(4).into();
        assert_eq!(
            cs,
            AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)])
        );
    }

    #[rstest]
    fn test_atom_constraints_ast_from_vec() {
        let cs: AtomConstraintsAst = vec![
            AtomConstraintAst::valence(4),
            AtomConstraintAst::donated_pairs(1),
        ]
        .into();
        assert_eq!(
            cs,
            AtomConstraintsAst::from_iter([
                AtomConstraintAst::valence(4),
                AtomConstraintAst::donated_pairs(1),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_empty(AtomConstraintsAst::new(), AtomConstraintsAst::new(), Some(AtomConstraintsAst::new()))]
    #[case::adds_kind_from_other(AtomConstraintsAst::new(), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), Some(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)])))]
    #[case::narrows_undetermined_to_lit(AtomConstraintsAst::from_iter([AtomConstraintAst::Valence(ValueAst::Undetermined)]), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]),
        Some(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)])))]
    #[case::lit_lit_match_preserved(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]),
        Some(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)])))]
    #[case::lit_lit_mismatch_none(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(3)]), None)]
    #[case::multi_kind_combines(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), AtomConstraintsAst::from_iter([AtomConstraintAst::degree(3)]),
        Some(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4), AtomConstraintAst::degree(3)])))]
    #[case::aromatic_valence_narrows(AtomConstraintsAst::from_iter([AtomConstraintAst::aromatic_valence(AromaticValenceAst::Undetermined)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::aromatic_valence(AromaticValenceAst::aromatic(1))]),
        Some(AtomConstraintsAst::from_iter([AtomConstraintAst::aromatic_valence(AromaticValenceAst::aromatic(1))])))]
    #[case::aromatic_valence_not_vs_aromatic_none(AtomConstraintsAst::from_iter([AtomConstraintAst::aromatic_valence(AromaticValenceAst::NotAromatic)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::aromatic_valence(AromaticValenceAst::aromatic(1))]), None)]
    #[case::ring_membership_size_unions(AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(5), 1)]), AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(6), 1)]),
        Some(AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(5), 1), AtomConstraintAst::ring_membership(RingScope::Size(6), 1)])))]
    #[case::ring_membership_size_dedup(AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(5), 1)]), AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(5), 1)]),
        Some(AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(5), 1)])))]
    #[case::prunes_vacuous(AtomConstraintsAst::new(), AtomConstraintsAst::from_iter([AtomConstraintAst::Valence(ValueAst::Undetermined)]), Some(AtomConstraintsAst::new()))]
    #[case::tetrahedral_narrows_from_absent(AtomConstraintsAst::new(),
        AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        Some(AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)])))]
    #[case::tetrahedral_not_stereo_vs_stereo_contradicts(AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::stereo(0_u32))]), None)]
    fn test_atom_constraints_ast_meet(
        #[case] a: AtomConstraintsAst,
        #[case] b: AtomConstraintsAst,
        #[case] expected: Option<AtomConstraintsAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::extends_self(AtomConstraintsAst::new(), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), true, AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]))]
    #[case::no_change(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), false,
        AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]))]
    #[case::contradiction_leaves_self_unchanged(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(3)]), false,
        AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]))]
    fn test_atom_constraints_ast_narrow_from(
        #[case] mut target: AtomConstraintsAst,
        #[case] source: AtomConstraintsAst,
        #[case] expected_changed: bool,
        #[case] expected_after: AtomConstraintsAst,
    ) {
        let changed = target.narrow_from(&source);
        assert_eq!(changed, expected_changed);
        assert_eq!(target, expected_after);
    }

    #[rstest]
    #[case::keeps_only_shared_kinds(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4), AtomConstraintAst::degree(2)]), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]))]
    #[case::widens_value(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(3)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::Valence(ValueAst::lit_set([4, 3]))]))]
    #[case::tetrahedral_same(AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]))]
    #[case::tetrahedral_incompatible_drops_to_undetermined(AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::stereo(0_u32))]), AtomConstraintsAst::new())]
    fn test_atom_constraints_ast_join(
        #[case] a: AtomConstraintsAst,
        #[case] b: AtomConstraintsAst,
        #[case] expected: AtomConstraintsAst,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rstest]
    #[case::empty_pattern_matches_anything(AtomConstraintsAst::new(), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), true)]
    #[case::missing_in_target_when_pattern_specific(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), AtomConstraintsAst::new(), false)]
    #[case::same_lit(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), true)]
    #[case::lit_lit_mismatch(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4)]), AtomConstraintsAst::from_iter([AtomConstraintAst::valence(3)]), false)]
    #[case::aromatic_wildcard_matches_aromatic(AtomConstraintsAst::from_iter([AtomConstraintAst::aromatic_valence(AromaticValenceAst::Undetermined)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::aromatic_valence(AromaticValenceAst::aromatic(1))]), true)]
    #[case::aromatic_not_vs_aromatic_mismatch(AtomConstraintsAst::from_iter([AtomConstraintAst::aromatic_valence(AromaticValenceAst::NotAromatic)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::aromatic_valence(AromaticValenceAst::aromatic(1))]), false)]
    #[case::ring_membership_size_subset(AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(5), 1)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(5), 1), AtomConstraintAst::ring_membership(RingScope::Size(6), 1)]), true)]
    #[case::ring_membership_size_not_present_in_target(AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(7), 1)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(5), 1), AtomConstraintAst::ring_membership(RingScope::Size(6), 1)]), false)]
    #[case::multi_kind_all_must_match(AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4), AtomConstraintAst::degree(3)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::valence(4), AtomConstraintAst::degree(2)]), false)]
    #[case::tetrahedral_same(AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]), true)]
    #[case::tetrahedral_pattern_specific_vs_absent(AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraintsAst::new(), false)]
    fn test_atom_constraints_ast_matches(
        #[case] pattern: AtomConstraintsAst,
        #[case] target: AtomConstraintsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
