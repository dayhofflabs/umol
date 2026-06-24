//! Atom constraints.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::mem::{self, replace};

use smallvec::SmallVec;
use strum::{EnumCount, EnumDiscriminants, EnumIter};

use super::super::constraint::ring::{RingMembershipAst, RingScope};
use super::super::error::Contradiction;
use super::super::remap::IdRemapping;
use super::super::stereo::TetrahedralStereoAst;
use super::super::traits::{AsLit, Canonicalize, Lattice};
use super::super::value::ValueAst;

/// Atom-scope constraint: a predicate that pattern-matches a single atom
/// on a topological or valence property (valence, degree, ring membership,
/// etc.). Held inline on `AtomAst` via `AtomConstraints`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, EnumDiscriminants)]
#[strum_discriminants(name(AtomConstraintKind), derive(Hash, EnumCount, EnumIter))]
#[repr(u8)]
pub enum AtomConstraint {
    Valence(ValueAst),
    TotalValence(ValueAst),
    AromaticValence(AromaticValenceAst),
    MulticenterValence(MulticenterValenceAst),
    DonatedPairs(ValueAst),
    AcceptedPairs(ValueAst),
    Degree(ValueAst),
    TotalDegree(ValueAst),
    RingDegree(ValueAst),
    RingValence(ValueAst),
    TotalHydrogens(ValueAst),
    RingMembership(RingMembershipAst),
    TetrahedralStereo(TetrahedralStereoAst),
}

impl AtomConstraint {
    pub fn valence(v: impl Into<ValueAst>) -> Self {
        Self::Valence(v.into())
    }

    pub fn total_valence(v: impl Into<ValueAst>) -> Self {
        Self::TotalValence(v.into())
    }

    pub fn aromatic_valence(v: AromaticValenceAst) -> Self {
        Self::AromaticValence(v)
    }

    pub fn multicenter_valence(v: MulticenterValenceAst) -> Self {
        Self::MulticenterValence(v)
    }

    pub fn donated_pairs(v: impl Into<ValueAst>) -> Self {
        Self::DonatedPairs(v.into())
    }

    pub fn accepted_pairs(v: impl Into<ValueAst>) -> Self {
        Self::AcceptedPairs(v.into())
    }

    pub fn degree(v: impl Into<ValueAst>) -> Self {
        Self::Degree(v.into())
    }

    pub fn total_degree(v: impl Into<ValueAst>) -> Self {
        Self::TotalDegree(v.into())
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

    pub fn tetrahedral_stereo(c: TetrahedralStereoAst) -> Self {
        Self::TetrahedralStereo(c)
    }

    pub fn kind(&self) -> AtomConstraintKind {
        self.into()
    }

    /// Entry identity for order/dedup: `kind()` plus `RingMembership`'s `RingScope`.
    pub fn key(&self) -> AtomConstraintKey {
        match self {
            Self::Valence(_) => AtomConstraintKey::Valence,
            Self::TotalValence(_) => AtomConstraintKey::TotalValence,
            Self::AromaticValence(_) => AtomConstraintKey::AromaticValence,
            Self::MulticenterValence(_) => AtomConstraintKey::MulticenterValence,
            Self::DonatedPairs(_) => AtomConstraintKey::DonatedPairs,
            Self::AcceptedPairs(_) => AtomConstraintKey::AcceptedPairs,
            Self::Degree(_) => AtomConstraintKey::Degree,
            Self::TotalDegree(_) => AtomConstraintKey::TotalDegree,
            Self::RingDegree(_) => AtomConstraintKey::RingDegree,
            Self::RingValence(_) => AtomConstraintKey::RingValence,
            Self::TotalHydrogens(_) => AtomConstraintKey::TotalHydrogens,
            Self::RingMembership(m) => AtomConstraintKey::RingMembership(m.scope),
            Self::TetrahedralStereo(_) => AtomConstraintKey::TetrahedralStereo,
        }
    }

    /// `false` for `RingMembership` (several per atom, one per `RingScope`); `true` otherwise.
    pub fn is_unique(&self) -> bool {
        !matches!(self.kind(), AtomConstraintKind::RingMembership)
    }

    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Valence(v)
            | Self::TotalValence(v)
            | Self::DonatedPairs(v)
            | Self::AcceptedPairs(v)
            | Self::Degree(v)
            | Self::TotalDegree(v)
            | Self::RingDegree(v)
            | Self::RingValence(v)
            | Self::TotalHydrogens(v) => v.is_undetermined(),
            Self::RingMembership(m) => m.count.is_undetermined(),
            Self::AromaticValence(c) => c.is_undetermined(),
            Self::MulticenterValence(c) => c.is_undetermined(),
            Self::TetrahedralStereo(c) => c.is_undetermined(),
        }
    }
}

impl Canonicalize for AtomConstraint {
    /// Canonicalize the inner value; kind and sub-key are preserved.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Valence(v) => Self::Valence(v.canonicalize()?),
            Self::TotalValence(v) => Self::TotalValence(v.canonicalize()?),
            Self::DonatedPairs(v) => Self::DonatedPairs(v.canonicalize()?),
            Self::AcceptedPairs(v) => Self::AcceptedPairs(v.canonicalize()?),
            Self::Degree(v) => Self::Degree(v.canonicalize()?),
            Self::TotalDegree(v) => Self::TotalDegree(v.canonicalize()?),
            Self::RingDegree(v) => Self::RingDegree(v.canonicalize()?),
            Self::RingValence(v) => Self::RingValence(v.canonicalize()?),
            Self::TotalHydrogens(v) => Self::TotalHydrogens(v.canonicalize()?),
            Self::RingMembership(m) => {
                Self::RingMembership(RingMembershipAst::new(m.scope, m.count.canonicalize()?))
            }
            Self::AromaticValence(c) => Self::AromaticValence(c.canonicalize()?),
            Self::MulticenterValence(c) => Self::MulticenterValence(c.canonicalize()?),
            Self::TetrahedralStereo(c) => Self::TetrahedralStereo(c.canonicalize()?),
        })
    }
}

/// Entry identity: discriminant + sub-key. Variant order matches `AtomConstraint`,
/// so `Ord` agrees with `kind as u8`; the ring run orders by `RingScope`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AtomConstraintKey {
    Valence,
    TotalValence,
    AromaticValence,
    MulticenterValence,
    DonatedPairs,
    AcceptedPairs,
    Degree,
    TotalDegree,
    RingDegree,
    RingValence,
    TotalHydrogens,
    RingMembership(RingScope),
    TetrahedralStereo,
}

impl AtomConstraintKey {
    pub fn kind(self) -> AtomConstraintKind {
        match self {
            Self::Valence => AtomConstraintKind::Valence,
            Self::TotalValence => AtomConstraintKind::TotalValence,
            Self::AromaticValence => AtomConstraintKind::AromaticValence,
            Self::MulticenterValence => AtomConstraintKind::MulticenterValence,
            Self::DonatedPairs => AtomConstraintKind::DonatedPairs,
            Self::AcceptedPairs => AtomConstraintKind::AcceptedPairs,
            Self::Degree => AtomConstraintKind::Degree,
            Self::TotalDegree => AtomConstraintKind::TotalDegree,
            Self::RingDegree => AtomConstraintKind::RingDegree,
            Self::RingValence => AtomConstraintKind::RingValence,
            Self::TotalHydrogens => AtomConstraintKind::TotalHydrogens,
            Self::RingMembership(_) => AtomConstraintKind::RingMembership,
            Self::TetrahedralStereo => AtomConstraintKind::TetrahedralStereo,
        }
    }
}

/// Aromatic-valence state of an atom: `Undetermined`, explicitly
/// `NotAromatic`, or participating in an aromatic system with the given
/// aromatic-valence count.
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

    fn join(&self, other: &Self) -> Self {
        let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::NotAromatic, Self::NotAromatic) => Self::NotAromatic,
            (Self::NotAromatic, Self::Aromatic(_)) | (Self::Aromatic(_), Self::NotAromatic) => {
                Self::Undetermined
            }
            (Self::Aromatic(p), Self::Aromatic(q)) => Self::Aromatic(p.join(q)),
        }
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

    fn join(&self, other: &Self) -> Self {
        let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::NotMulticenter, Self::NotMulticenter) => Self::NotMulticenter,
            (Self::NotMulticenter, Self::Multicenter(_))
            | (Self::Multicenter(_), Self::NotMulticenter) => Self::Undetermined,
            (Self::Multicenter(p), Self::Multicenter(q)) => Self::Multicenter(p.join(q)),
        }
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

/// Per-atom constraints.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomConstraints {
    entries: SmallVec<[AtomConstraint; 2]>,
}

impl AtomConstraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn contains(&self, kind: AtomConstraintKind) -> bool {
        self.find(kind).is_ok()
    }

    pub fn get(&self, kind: AtomConstraintKind) -> Option<&AtomConstraint> {
        self.find(kind).ok().map(|i| &self.entries[i])
    }

    pub fn get_mut(&mut self, kind: AtomConstraintKind) -> Option<&mut AtomConstraint> {
        match self.find(kind) {
            Ok(i) => Some(&mut self.entries[i]),
            Err(_) => None,
        }
    }

    pub fn valence(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKind::Valence) {
            Some(AtomConstraint::Valence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn total_valence(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKind::TotalValence) {
            Some(AtomConstraint::TotalValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn degree(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKind::Degree) {
            Some(AtomConstraint::Degree(v)) => Some(v),
            _ => None,
        }
    }

    pub fn total_degree(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKind::TotalDegree) {
            Some(AtomConstraint::TotalDegree(v)) => Some(v),
            _ => None,
        }
    }

    pub fn ring_degree(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKind::RingDegree) {
            Some(AtomConstraint::RingDegree(v)) => Some(v),
            _ => None,
        }
    }

    pub fn ring_valence(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKind::RingValence) {
            Some(AtomConstraint::RingValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn total_hydrogens(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKind::TotalHydrogens) {
            Some(AtomConstraint::TotalHydrogens(v)) => Some(v),
            _ => None,
        }
    }

    pub fn donated_pairs(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKind::DonatedPairs) {
            Some(AtomConstraint::DonatedPairs(v)) => Some(v),
            _ => None,
        }
    }

    pub fn accepted_pairs(&self) -> Option<&ValueAst> {
        match self.get(AtomConstraintKind::AcceptedPairs) {
            Some(AtomConstraint::AcceptedPairs(v)) => Some(v),
            _ => None,
        }
    }

    fn ring_memberships(&self) -> impl Iterator<Item = (RingScope, &ValueAst)> {
        self.get_all(AtomConstraintKind::RingMembership)
            .filter_map(|c| match c {
                AtomConstraint::RingMembership(m) => Some((m.scope, &m.count)),
                _ => None,
            })
    }

    fn ring_membership_value(&self, scope: RingScope) -> Option<&ValueAst> {
        self.ring_memberships()
            .find(|(s, _)| *s == scope)
            .map(|(_, v)| v)
    }

    pub fn ring_count(&self) -> Option<&ValueAst> {
        self.ring_membership_value(RingScope::All)
    }

    pub fn ring_size_count(&self, s: u8) -> Option<&ValueAst> {
        self.ring_membership_value(RingScope::Size(s))
    }

    pub fn aromatic_valence(&self) -> Option<&AromaticValenceAst> {
        match self.get(AtomConstraintKind::AromaticValence) {
            Some(AtomConstraint::AromaticValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn multicenter_valence(&self) -> Option<&MulticenterValenceAst> {
        match self.get(AtomConstraintKind::MulticenterValence) {
            Some(AtomConstraint::MulticenterValence(v)) => Some(v),
            _ => None,
        }
    }

    pub fn tetrahedral_stereo(&self) -> Option<&TetrahedralStereoAst> {
        match self.get(AtomConstraintKind::TetrahedralStereo) {
            Some(AtomConstraint::TetrahedralStereo(c)) => Some(c),
            _ => None,
        }
    }

    /// Insert at the `key()`-sorted position: unique kinds replace the same-key
    /// entry (returning it); ring appends, leaving duplicates for lazy dedup.
    pub fn add(&mut self, constraint: AtomConstraint) -> Option<AtomConstraint> {
        match self.find_by_key(constraint.key()) {
            Ok(i) if constraint.is_unique() => Some(replace(&mut self.entries[i], constraint)),
            Ok(i) => {
                let end = i + self.entries[i..]
                    .iter()
                    .take_while(|c| c.key() == constraint.key())
                    .count();
                self.entries.insert(end, constraint);
                None
            }
            Err(i) => {
                self.entries.insert(i, constraint);
                None
            }
        }
    }

    /// Add multiple constraints at once, using semantics of `add`.
    pub fn extend(&mut self, constraints: impl IntoIterator<Item = AtomConstraint>) {
        for constraint in constraints {
            self.add(constraint);
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&AtomConstraint) -> bool) {
        self.entries.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// Move the entries out of the store, leaving it empty. Returned items
    /// are in the store's internal sorted-by-kind order.
    pub fn take(&mut self) -> impl Iterator<Item = AtomConstraint> {
        mem::take(&mut self.entries).into_iter()
    }

    pub fn remove(&mut self, kind: AtomConstraintKind) -> Option<AtomConstraint> {
        self.find(kind).ok().map(|i| self.entries.remove(i))
    }

    /// Remove the first entry exactly equal to `constraint`. Returns the
    /// removed entry if found; otherwise `None`.
    pub fn remove_entry(&mut self, constraint: &AtomConstraint) -> Option<AtomConstraint> {
        let pos = self.entries.iter().position(|c| c == constraint)?;
        Some(self.entries.remove(pos))
    }

    /// True if any entry exactly equals `constraint`.
    pub fn contains_entry(&self, constraint: &AtomConstraint) -> bool {
        self.entries.iter().any(|c| c == constraint)
    }

    /// Iterate over every entry of `kind`. Single-valued kinds yield at most
    /// one entry; `RingMembership` may yield several (one per scope).
    pub fn get_all(&self, kind: AtomConstraintKind) -> impl Iterator<Item = &AtomConstraint> {
        let start = self
            .entries
            .partition_point(|c| (c.kind() as u8) < (kind as u8));
        self.entries[start..]
            .iter()
            .take_while(move |c| c.kind() == kind)
    }

    /// Remove every entry of `kind`, returning them in store order. Single-
    /// valued kinds drain at most one entry.
    pub fn remove_all(&mut self, kind: AtomConstraintKind) -> Vec<AtomConstraint> {
        let start = self
            .entries
            .partition_point(|c| (c.kind() as u8) < (kind as u8));
        let end = start
            + self.entries[start..]
                .iter()
                .take_while(|c| c.kind() == kind)
                .count();
        self.entries.drain(start..end).collect()
    }

    pub fn iter(&self) -> impl Iterator<Item = &AtomConstraint> {
        self.entries.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut AtomConstraint> {
        self.entries.iter_mut()
    }

    /// No-op: no `AtomConstraint` variant carries an entity index.
    pub fn remap(self, _remap: &IdRemapping) -> Self {
        self
    }

    fn find(&self, kind: AtomConstraintKind) -> Result<usize, usize> {
        self.entries
            .binary_search_by_key(&(kind as u8), |c| c.kind() as u8)
    }

    fn find_by_key(&self, key: AtomConstraintKey) -> Result<usize, usize> {
        self.entries.binary_search_by(|c| c.key().cmp(&key))
    }

    pub fn contains_key(&self, key: AtomConstraintKey) -> bool {
        self.find_by_key(key).is_ok()
    }

    pub fn get_by_key(&self, key: AtomConstraintKey) -> Option<&AtomConstraint> {
        self.find_by_key(key).ok().map(|i| &self.entries[i])
    }

    pub fn get_by_key_mut(&mut self, key: AtomConstraintKey) -> Option<&mut AtomConstraint> {
        self.find_by_key(key).ok().map(|i| &mut self.entries[i])
    }

    pub fn remove_by_key(&mut self, key: AtomConstraintKey) -> Option<AtomConstraint> {
        self.find_by_key(key).ok().map(|i| self.entries.remove(i))
    }
}

impl Canonicalize for AtomConstraints {
    /// Sort by `key()`, canonicalize each value, merge same-scope ring entries
    /// by value-`meet` (`Err` on contradiction), drop vacuous entries.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut input = self.entries;
        input.sort_by_key(|c| c.key());
        let mut entries: SmallVec<[AtomConstraint; 2]> = SmallVec::new();
        for c in input {
            let c = c.canonicalize()?;
            if let (
                Some(AtomConstraint::RingMembership(prev)),
                AtomConstraint::RingMembership(next),
            ) = (entries.last_mut(), &c)
            {
                if prev.scope == next.scope {
                    prev.count = prev.count.meet(&next.count).ok_or(Contradiction)?;
                    continue;
                }
            }
            entries.push(c);
        }
        entries.retain(|c| !c.is_undetermined());
        Ok(Self { entries })
    }
}

impl Lattice for AtomConstraints {
    fn is_undetermined(&self) -> bool {
        self.iter().all(|c| c.is_undetermined())
    }

    fn is_ground(&self) -> bool {
        self.iter().all(|c| match c {
            AtomConstraint::Valence(v)
            | AtomConstraint::TotalValence(v)
            | AtomConstraint::DonatedPairs(v)
            | AtomConstraint::AcceptedPairs(v)
            | AtomConstraint::Degree(v)
            | AtomConstraint::TotalDegree(v)
            | AtomConstraint::RingDegree(v)
            | AtomConstraint::RingValence(v)
            | AtomConstraint::TotalHydrogens(v) => v.is_ground(),
            AtomConstraint::RingMembership(m) => m.count.is_ground(),
            AtomConstraint::AromaticValence(c) => c.is_ground(),
            AtomConstraint::MulticenterValence(c) => c.is_ground(),
            AtomConstraint::TetrahedralStereo(c) => c.is_ground(),
        })
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let mut result = Self::new();
        let meet_val = |a: Option<&ValueAst>, b: Option<&ValueAst>| {
            a.unwrap_or(&ValueAst::Undetermined)
                .meet(b.unwrap_or(&ValueAst::Undetermined))
        };
        let v = meet_val(self.valence(), other.valence())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::Valence(v));
        }
        let v = meet_val(self.total_valence(), other.total_valence())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::TotalValence(v));
        }
        let v = self
            .aromatic_valence()
            .unwrap_or(&AromaticValenceAst::Undetermined)
            .meet(
                other
                    .aromatic_valence()
                    .unwrap_or(&AromaticValenceAst::Undetermined),
            )?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::AromaticValence(v));
        }
        let v = self
            .multicenter_valence()
            .unwrap_or(&MulticenterValenceAst::Undetermined)
            .meet(
                other
                    .multicenter_valence()
                    .unwrap_or(&MulticenterValenceAst::Undetermined),
            )?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::MulticenterValence(v));
        }
        let v = meet_val(self.donated_pairs(), other.donated_pairs())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::DonatedPairs(v));
        }
        let v = meet_val(self.accepted_pairs(), other.accepted_pairs())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::AcceptedPairs(v));
        }
        let v = meet_val(self.degree(), other.degree())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::Degree(v));
        }
        let v = meet_val(self.total_degree(), other.total_degree())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::TotalDegree(v));
        }
        let v = meet_val(self.ring_degree(), other.ring_degree())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::RingDegree(v));
        }
        let v = meet_val(self.ring_valence(), other.ring_valence())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::RingValence(v));
        }
        let v = meet_val(self.total_hydrogens(), other.total_hydrogens())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::TotalHydrogens(v));
        }
        let v = self
            .tetrahedral_stereo()
            .unwrap_or(&TetrahedralStereoAst::Undetermined)
            .meet(
                other
                    .tetrahedral_stereo()
                    .unwrap_or(&TetrahedralStereoAst::Undetermined),
            )?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::TetrahedralStereo(v));
        }
        let mut scopes: BTreeSet<RingScope> = self.ring_memberships().map(|(s, _)| s).collect();
        scopes.extend(other.ring_memberships().map(|(s, _)| s));
        for scope in scopes {
            let v = meet_val(
                self.ring_membership_value(scope),
                other.ring_membership_value(scope),
            )?;
            if !v.is_undetermined() {
                result.add(AtomConstraint::RingMembership(RingMembershipAst::new(
                    scope, v,
                )));
            }
        }
        Some(result)
    }

    fn join(&self, other: &Self) -> Self {
        let mut result = Self::new();
        macro_rules! join_unique_value {
            ($kind:ident, $accessor:ident, $variant:ident) => {
                if self.contains(AtomConstraintKind::$kind)
                    && other.contains(AtomConstraintKind::$kind)
                {
                    // both present (guarded above), so the accessors are `Some`
                    let joined = self.$accessor().unwrap().join(other.$accessor().unwrap());
                    if !joined.is_undetermined() {
                        result.add(AtomConstraint::$variant(joined));
                    }
                }
            };
        }
        join_unique_value!(Valence, valence, Valence);
        join_unique_value!(TotalValence, total_valence, TotalValence);
        join_unique_value!(AromaticValence, aromatic_valence, AromaticValence);
        join_unique_value!(MulticenterValence, multicenter_valence, MulticenterValence);
        join_unique_value!(DonatedPairs, donated_pairs, DonatedPairs);
        join_unique_value!(AcceptedPairs, accepted_pairs, AcceptedPairs);
        join_unique_value!(Degree, degree, Degree);
        join_unique_value!(TotalDegree, total_degree, TotalDegree);
        join_unique_value!(RingDegree, ring_degree, RingDegree);
        join_unique_value!(RingValence, ring_valence, RingValence);
        join_unique_value!(TotalHydrogens, total_hydrogens, TotalHydrogens);
        join_unique_value!(TetrahedralStereo, tetrahedral_stereo, TetrahedralStereo);
        for (scope, v) in self.ring_memberships() {
            if other.ring_memberships().any(|(s, _)| s == scope) {
                let j = v.join(other.ring_membership_value(scope).unwrap());
                if !j.is_undetermined() {
                    result.add(AtomConstraint::RingMembership(RingMembershipAst::new(
                        scope, j,
                    )));
                }
            }
        }
        result
    }

    /// Pattern-driven: every constraint the pattern carries must match the target's
    /// corresponding value, looked up by reference (an absent target value is
    /// `Undetermined`, which matches). Cost is proportional to the pattern's
    /// constraints, not the field count; an empty pattern matches any target.
    fn matches(&self, target: &Self) -> bool {
        self.iter().all(|c| match c {
            AtomConstraint::Valence(v) => {
                v.matches(target.valence().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraint::TotalValence(v) => {
                v.matches(target.total_valence().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraint::AromaticValence(av) => av.matches(
                target
                    .aromatic_valence()
                    .unwrap_or(&AromaticValenceAst::Undetermined),
            ),
            AtomConstraint::MulticenterValence(mv) => mv.matches(
                target
                    .multicenter_valence()
                    .unwrap_or(&MulticenterValenceAst::Undetermined),
            ),
            AtomConstraint::DonatedPairs(v) => {
                v.matches(target.donated_pairs().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraint::AcceptedPairs(v) => {
                v.matches(target.accepted_pairs().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraint::Degree(v) => {
                v.matches(target.degree().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraint::TotalDegree(v) => {
                v.matches(target.total_degree().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraint::RingDegree(v) => {
                v.matches(target.ring_degree().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraint::RingValence(v) => {
                v.matches(target.ring_valence().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraint::TotalHydrogens(v) => {
                v.matches(target.total_hydrogens().unwrap_or(&ValueAst::Undetermined))
            }
            AtomConstraint::RingMembership(rm) => rm.count.matches(
                target
                    .ring_membership_value(rm.scope)
                    .unwrap_or(&ValueAst::Undetermined),
            ),
            AtomConstraint::TetrahedralStereo(ts) => ts.matches(
                target
                    .tetrahedral_stereo()
                    .unwrap_or(&TetrahedralStereoAst::Undetermined),
            ),
        })
    }
}

impl FromIterator<AtomConstraint> for AtomConstraints {
    fn from_iter<I: IntoIterator<Item = AtomConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

impl IntoIterator for AtomConstraints {
    type Item = AtomConstraint;
    type IntoIter = smallvec::IntoIter<[AtomConstraint; 2]>;

    fn into_iter(self) -> Self::IntoIter {
        self.entries.into_iter()
    }
}

impl From<AtomConstraint> for AtomConstraints {
    fn from(c: AtomConstraint) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<AtomConstraint>> for AtomConstraints {
    fn from(cs: Vec<AtomConstraint>) -> Self {
        Self::from_iter(cs)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Remapping;

    use super::*;
    use crate::ast::value::ValueTerm;

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraint::valence(4), AtomConstraint::Valence(ValueAst::Lit(4)))]
    #[case::total_valence(AtomConstraint::total_valence(5), AtomConstraint::TotalValence(ValueAst::Lit(5)))]
    #[case::donated_pairs(AtomConstraint::donated_pairs(1), AtomConstraint::DonatedPairs(ValueAst::Lit(1)))]
    #[case::accepted_pairs(AtomConstraint::accepted_pairs(2), AtomConstraint::AcceptedPairs(ValueAst::Lit(2)))]
    #[case::degree(AtomConstraint::degree(3), AtomConstraint::Degree(ValueAst::Lit(3)))]
    #[case::total_degree(AtomConstraint::total_degree(4), AtomConstraint::TotalDegree(ValueAst::Lit(4)))]
    #[case::ring_degree(AtomConstraint::ring_degree(2), AtomConstraint::RingDegree(ValueAst::Lit(2)))]
    #[case::ring_valence(AtomConstraint::ring_valence(3), AtomConstraint::RingValence(ValueAst::Lit(3)))]
    #[case::total_hydrogens(AtomConstraint::total_hydrogens(3), AtomConstraint::TotalHydrogens(ValueAst::Lit(3)))]
    #[case::ring_membership_all(AtomConstraint::ring_membership(RingScope::All, 1), AtomConstraint::RingMembership(RingMembershipAst { scope: RingScope::All, count: ValueAst::Lit(1) }))]
    #[case::ring_membership_size(AtomConstraint::ring_membership(RingScope::Size(6), 1), AtomConstraint::RingMembership(RingMembershipAst { scope: RingScope::Size(6), count: ValueAst::Lit(1) }))]
    #[case::aromatic_valence(
        AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic),
        AtomConstraint::AromaticValence(AromaticValenceAst::NotAromatic),
    )]
    #[case::multicenter_valence(
        AtomConstraint::multicenter_valence(MulticenterValenceAst::NotMulticenter),
        AtomConstraint::MulticenterValence(MulticenterValenceAst::NotMulticenter),
    )]
    #[case::tetrahedral_stereo(
        AtomConstraint::tetrahedral_stereo(TetrahedralStereoAst::NotStereo),
        AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo),
    )]
    fn test_atom_constraint_constructors(
        #[case] actual: AtomConstraint,
        #[case] expected: AtomConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraint::valence(4), AtomConstraintKind::Valence)]
    #[case::total_valence(AtomConstraint::total_valence(5), AtomConstraintKind::TotalValence)]
    #[case::aromatic_valence(AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic), AtomConstraintKind::AromaticValence)]
    #[case::multicenter_valence(AtomConstraint::multicenter_valence(MulticenterValenceAst::Undetermined), AtomConstraintKind::MulticenterValence)]
    #[case::donated_pairs(AtomConstraint::donated_pairs(1), AtomConstraintKind::DonatedPairs)]
    #[case::accepted_pairs(AtomConstraint::accepted_pairs(2), AtomConstraintKind::AcceptedPairs)]
    #[case::degree(AtomConstraint::degree(3), AtomConstraintKind::Degree)]
    #[case::total_degree(AtomConstraint::total_degree(4), AtomConstraintKind::TotalDegree)]
    #[case::ring_degree(AtomConstraint::ring_degree(2), AtomConstraintKind::RingDegree)]
    #[case::ring_valence(AtomConstraint::ring_valence(3), AtomConstraintKind::RingValence)]
    #[case::total_hydrogens(AtomConstraint::total_hydrogens(3), AtomConstraintKind::TotalHydrogens)]
    #[case::ring_membership_all(AtomConstraint::ring_membership(RingScope::All, 1), AtomConstraintKind::RingMembership)]
    #[case::ring_membership_size(AtomConstraint::ring_membership(RingScope::Size(6), 1), AtomConstraintKind::RingMembership)]
    #[case::tetrahedral_stereo(AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo), AtomConstraintKind::TetrahedralStereo)]
    fn test_atom_constraint_kind(
        #[case] constraint: AtomConstraint,
        #[case] expected: AtomConstraintKind,
    ) {
        assert_eq!(constraint.kind(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraint::valence(4), AtomConstraintKey::Valence)]
    #[case::total_valence(AtomConstraint::total_valence(5), AtomConstraintKey::TotalValence)]
    #[case::aromatic_valence(AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic), AtomConstraintKey::AromaticValence)]
    #[case::multicenter_valence(AtomConstraint::multicenter_valence(MulticenterValenceAst::Undetermined), AtomConstraintKey::MulticenterValence)]
    #[case::donated_pairs(AtomConstraint::donated_pairs(1), AtomConstraintKey::DonatedPairs)]
    #[case::accepted_pairs(AtomConstraint::accepted_pairs(2), AtomConstraintKey::AcceptedPairs)]
    #[case::degree(AtomConstraint::degree(3), AtomConstraintKey::Degree)]
    #[case::total_degree(AtomConstraint::total_degree(4), AtomConstraintKey::TotalDegree)]
    #[case::ring_degree(AtomConstraint::ring_degree(2), AtomConstraintKey::RingDegree)]
    #[case::ring_valence(AtomConstraint::ring_valence(3), AtomConstraintKey::RingValence)]
    #[case::total_hydrogens(AtomConstraint::total_hydrogens(3), AtomConstraintKey::TotalHydrogens)]
    #[case::ring_membership_all(AtomConstraint::ring_membership(RingScope::All, 1), AtomConstraintKey::RingMembership(RingScope::All))]
    #[case::ring_membership_size(AtomConstraint::ring_membership(RingScope::Size(6), 1), AtomConstraintKey::RingMembership(RingScope::Size(6)))]
    #[case::tetrahedral_stereo(AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo), AtomConstraintKey::TetrahedralStereo)]
    fn test_atom_constraint_key(
        #[case] constraint: AtomConstraint,
        #[case] expected: AtomConstraintKey,
    ) {
        assert_eq!(constraint.key(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraintKey::Valence, AtomConstraintKind::Valence)]
    #[case::ring_membership_all(AtomConstraintKey::RingMembership(RingScope::All), AtomConstraintKind::RingMembership)]
    #[case::ring_membership_size(AtomConstraintKey::RingMembership(RingScope::Size(6)), AtomConstraintKind::RingMembership)]
    #[case::tetrahedral_stereo(AtomConstraintKey::TetrahedralStereo, AtomConstraintKind::TetrahedralStereo)]
    fn test_atom_constraint_key_kind(
        #[case] key: AtomConstraintKey,
        #[case] expected: AtomConstraintKind,
    ) {
        assert_eq!(key.kind(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_litset_singleton(AtomConstraint::Valence(ValueAst::lit_set([4])), Ok(AtomConstraint::valence(4)))]
    #[case::ring_count_litset_singleton(
        AtomConstraint::RingMembership(RingMembershipAst::new(RingScope::Size(6), ValueAst::lit_set([2]))),
        Ok(AtomConstraint::ring_membership(RingScope::Size(6), 2)))]
    #[case::empty_litset_contradiction(AtomConstraint::Valence(ValueAst::lit_set(Vec::<i64>::new())), Err(Contradiction))]
    fn test_atom_constraint_canonicalize(
        #[case] constraint: AtomConstraint,
        #[case] expected: Result<AtomConstraint, Contradiction>,
    ) {
        assert_eq!(constraint.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_lit(AtomConstraint::valence(4), false)]
    #[case::valence_undetermined(AtomConstraint::Valence(ValueAst::Undetermined), true)]
    #[case::degree_undetermined(AtomConstraint::Degree(ValueAst::Undetermined), true)]
    #[case::ring_membership_undetermined(AtomConstraint::ring_membership(RingScope::All, ValueAst::Undetermined), true)]
    #[case::aromatic_undetermined(AtomConstraint::aromatic_valence(AromaticValenceAst::Undetermined), true)]
    #[case::aromatic_not_aromatic(AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic), false)]
    #[case::aromatic_with_value(AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1)), false)]
    #[case::multicenter_undetermined(AtomConstraint::multicenter_valence(MulticenterValenceAst::Undetermined), true)]
    #[case::multicenter_not(AtomConstraint::multicenter_valence(MulticenterValenceAst::NotMulticenter), false)]
    #[case::multicenter_with_value(AtomConstraint::multicenter_valence(MulticenterValenceAst::multicenter(1)), false)]
    #[case::tetrahedral_not_stereo(AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo), false)]
    #[case::tetrahedral_undetermined(AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::Undetermined), true)]
    fn test_atom_constraint_is_undetermined(
        #[case] c: AtomConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
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
        assert_eq!(a.join(&b), expected);
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
        assert_eq!(a.join(&b), expected);
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
    fn test_atom_constraints_new() {
        let cs = AtomConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_present(AtomConstraintKind::Valence, true)]
    #[case::aromatic_present(AtomConstraintKind::AromaticValence, true)]
    #[case::degree_absent(AtomConstraintKind::Degree, false)]
    fn test_atom_constraints_contains(
        #[case] kind: AtomConstraintKind,
        #[case] expected: bool,
    ) {
        let cs = AtomConstraints::from_iter([
            AtomConstraint::valence(4),
            AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic),
        ]);
        assert_eq!(cs.contains(kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_present(AtomConstraintKind::Valence, Some(AtomConstraint::valence(4)))]
    #[case::aromatic_present(AtomConstraintKind::AromaticValence, Some(AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic)))]
    #[case::degree_absent(AtomConstraintKind::Degree, None)]
    fn test_atom_constraints_get(
        #[case] kind: AtomConstraintKind,
        #[case] expected: Option<AtomConstraint>,
    ) {
        let cs = AtomConstraints::from_iter([
            AtomConstraint::valence(4),
            AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic),
        ]);
        assert_eq!(cs.get(kind), expected.as_ref());
    }

    #[rstest]
    fn test_atom_constraints_get_mut() {
        let mut cs = AtomConstraints::from_iter([AtomConstraint::valence(3)]);
        let slot = cs.get_mut(AtomConstraintKind::Valence).unwrap();
        *slot = AtomConstraint::valence(5);
        assert_eq!(
            cs.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::valence(5)),
        );
    }

    #[rstest]
    fn test_atom_constraints_get_mut_absent() {
        let mut cs = AtomConstraints::new();
        assert!(cs.get_mut(AtomConstraintKind::Valence).is_none());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_present(AtomConstraintKey::Valence, true)]
    #[case::ring_all_present(AtomConstraintKey::RingMembership(RingScope::All), true)]
    #[case::ring_size_present(AtomConstraintKey::RingMembership(RingScope::Size(6)), true)]
    #[case::ring_size_absent(AtomConstraintKey::RingMembership(RingScope::Size(5)), false)]
    #[case::degree_absent(AtomConstraintKey::Degree, false)]
    fn test_atom_constraints_contains_key(
        #[case] key: AtomConstraintKey,
        #[case] expected: bool,
    ) {
        let cs = AtomConstraints::from_iter([
            AtomConstraint::valence(4),
            AtomConstraint::ring_membership(RingScope::All, 2),
            AtomConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.contains_key(key), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_all(AtomConstraintKey::RingMembership(RingScope::All), Some(AtomConstraint::ring_membership(RingScope::All, 2)))]
    #[case::ring_size(AtomConstraintKey::RingMembership(RingScope::Size(6)), Some(AtomConstraint::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_size_absent(AtomConstraintKey::RingMembership(RingScope::Size(5)), None)]
    #[case::valence(AtomConstraintKey::Valence, Some(AtomConstraint::valence(4)))]
    fn test_atom_constraints_get_by_key(
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraint>,
    ) {
        let cs = AtomConstraints::from_iter([
            AtomConstraint::valence(4),
            AtomConstraint::ring_membership(RingScope::All, 2),
            AtomConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        assert_eq!(cs.get_by_key(key), expected.as_ref());
    }

    #[rstest]
    fn test_atom_constraints_get_by_key_mut() {
        let mut cs = AtomConstraints::from_iter([
            AtomConstraint::ring_membership(RingScope::All, 2),
            AtomConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let slot = cs
            .get_by_key_mut(AtomConstraintKey::RingMembership(RingScope::Size(6)))
            .unwrap();
        *slot = AtomConstraint::ring_membership(RingScope::Size(6), 2);
        assert_eq!(
            cs.get_by_key(AtomConstraintKey::RingMembership(RingScope::Size(6))),
            Some(&AtomConstraint::ring_membership(RingScope::Size(6), 2)),
        );
    }

    #[rstest]
    fn test_atom_constraints_remove_by_key() {
        let mut cs = AtomConstraints::from_iter([
            AtomConstraint::valence(4),
            AtomConstraint::ring_membership(RingScope::All, 2),
            AtomConstraint::ring_membership(RingScope::Size(6), 1),
        ]);
        let removed = cs.remove_by_key(AtomConstraintKey::RingMembership(RingScope::Size(6)));
        assert_eq!(
            removed,
            Some(AtomConstraint::ring_membership(RingScope::Size(6), 1)),
        );
        assert_eq!(
            cs.iter().cloned().collect::<Vec<_>>(),
            vec![
                AtomConstraint::valence(4),
                AtomConstraint::ring_membership(RingScope::All, 2),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::merge_same_scope(
        AtomConstraints::from_iter([
            AtomConstraint::ring_membership(RingScope::All, ValueAst::lit_set([1, 2])),
            AtomConstraint::ring_membership(RingScope::All, ValueAst::lit_set([2, 3])),
        ]),
        Ok(AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::All, 2)])))]
    #[case::drop_vacuous(
        AtomConstraints::from_iter([
            AtomConstraint::Valence(ValueAst::Undetermined),
            AtomConstraint::degree(3),
        ]),
        Ok(AtomConstraints::from_iter([AtomConstraint::degree(3)])))]
    #[case::canonicalizes_values(
        AtomConstraints::from_iter([
            AtomConstraint::Degree(ValueAst::lit_set([3])),
            AtomConstraint::Valence(ValueAst::lit_set([4])),
        ]),
        Ok(AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(3)])))]
    #[case::contradiction_same_scope(
        AtomConstraints::from_iter([
            AtomConstraint::ring_membership(RingScope::All, 1),
            AtomConstraint::ring_membership(RingScope::All, 0),
        ]),
        Err(Contradiction))]
    fn test_atom_constraints_canonicalize(
        #[case] constraints: AtomConstraints,
        #[case] expected: Result<AtomConstraints, Contradiction>,
    ) {
        assert_eq!(constraints.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(vec![AtomConstraint::valence(4)], vec![None], vec![AtomConstraint::valence(4)])]
    #[case::replace_same_kind(vec![AtomConstraint::valence(3), AtomConstraint::valence(4)], vec![None, Some(AtomConstraint::valence(3))], vec![AtomConstraint::valence(4)])]
    #[case::distinct_kinds(vec![AtomConstraint::valence(4), AtomConstraint::degree(3), AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic)],
        vec![None, None, None], vec![AtomConstraint::valence(4), AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic), AtomConstraint::degree(3)])]
    fn test_atom_constraints_add(
        #[case] sequence: Vec<AtomConstraint>,
        #[case] expected_returns: Vec<Option<AtomConstraint>>,
        #[case] expected_state: Vec<AtomConstraint>,
    ) {
        let mut cs = AtomConstraints::new();
        let returns: Vec<_> = sequence.into_iter().map(|c| cs.add(c)).collect();
        assert_eq!(returns, expected_returns);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(collected, expected_state);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(|c: &AtomConstraint| matches!(c, AtomConstraint::Valence(_) | AtomConstraint::RingMembership(_)), vec![AtomConstraint::valence(4), AtomConstraint::ring_membership(RingScope::All, 2)])]
    #[case::all_dropped(|_: &AtomConstraint| false, vec![])]
    fn test_atom_constraints_retain(
        #[case] predicate: impl FnMut(&AtomConstraint) -> bool,
        #[case] expected: Vec<AtomConstraint>,
    ) {
        let mut cs = AtomConstraints::from_iter([
            AtomConstraint::valence(4),
            AtomConstraint::degree(3),
            AtomConstraint::ring_membership(RingScope::All, 2),
        ]);
        cs.retain(predicate);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(collected, expected);
    }

    #[rstest]
    fn test_atom_constraints_clear() {
        let mut cs =
            AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(3)]);
        cs.clear();
        assert_eq!(cs, AtomConstraints::new());
    }

    #[rstest]
    fn test_atom_constraints_take() {
        let mut cs =
            AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(3)]);
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(
            drained,
            vec![AtomConstraint::valence(4), AtomConstraint::degree(3)],
        );
        assert_eq!(cs, AtomConstraints::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_present(AtomConstraintKind::Valence, Some(AtomConstraint::valence(4)), vec![AtomConstraint::degree(3)])]
    #[case::degree_present(AtomConstraintKind::Degree, Some(AtomConstraint::degree(3)), vec![AtomConstraint::valence(4)])]
    #[case::absent(AtomConstraintKind::RingMembership, None, vec![AtomConstraint::valence(4), AtomConstraint::degree(3)])]
    fn test_atom_constraints_remove(
        #[case] kind: AtomConstraintKind,
        #[case] expected_returned: Option<AtomConstraint>,
        #[case] expected_state: Vec<AtomConstraint>,
    ) {
        let mut cs = AtomConstraints::from_iter([
            AtomConstraint::valence(4),
            AtomConstraint::degree(3),
        ]);
        assert_eq!(cs.remove(kind), expected_returned);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(collected, expected_state);
    }

    #[rstest]
    fn test_atom_constraints_iter() {
        let cs = AtomConstraints::from_iter([
            AtomConstraint::ring_membership(RingScope::Size(6), 1),
            AtomConstraint::valence(4),
            AtomConstraint::degree(3),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                AtomConstraint::valence(4),
                AtomConstraint::degree(3),
                AtomConstraint::ring_membership(RingScope::Size(6), 1),
            ],
        );
    }

    #[rstest]
    fn test_atom_constraints_iter_mut() {
        let mut cs =
            AtomConstraints::from_iter([AtomConstraint::valence(3), AtomConstraint::degree(2)]);
        for c in cs.iter_mut() {
            if let AtomConstraint::Valence(v) = c {
                *v = ValueAst::Lit(7);
            }
        }
        assert_eq!(
            cs,
            AtomConstraints::from_iter([AtomConstraint::valence(7), AtomConstraint::degree(2),]),
        );
    }

    #[rstest]
    fn test_atom_constraints_remap() {
        let cs =
            AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(3)]);
        let remap = IdRemapping::new(
            Remapping::new(vec![0, 1, 2], vec![0]),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(cs.clone().remap(&remap), cs);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(vec![AtomConstraint::valence(4), AtomConstraint::degree(3)], vec![AtomConstraint::valence(4), AtomConstraint::degree(3)])]
    #[case::same_kind_last_wins(vec![AtomConstraint::valence(3), AtomConstraint::valence(4)], vec![AtomConstraint::valence(4)])]
    #[case::empty(vec![], vec![])]
    fn test_atom_constraints_from_iter(
        #[case] input: Vec<AtomConstraint>,
        #[case] expected: Vec<AtomConstraint>,
    ) {
        let cs = AtomConstraints::from_iter(input);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(collected, expected);
    }

    #[rstest]
    fn test_atom_constraints_into_iter() {
        let cs =
            AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(3)]);
        let collected: Vec<AtomConstraint> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![AtomConstraint::valence(4), AtomConstraint::degree(3)],
        );
    }

    #[rstest]
    fn test_atom_constraints_from_atom_constraint() {
        let cs: AtomConstraints = AtomConstraint::valence(4).into();
        assert_eq!(cs, AtomConstraints::from_iter([AtomConstraint::valence(4)]));
    }

    #[rstest]
    fn test_atom_constraints_from_vec() {
        let cs: AtomConstraints =
            vec![AtomConstraint::valence(4), AtomConstraint::donated_pairs(1)].into();
        assert_eq!(
            cs,
            AtomConstraints::from_iter([
                AtomConstraint::valence(4),
                AtomConstraint::donated_pairs(1),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty_empty(AtomConstraints::new(), AtomConstraints::new(), Some(AtomConstraints::new()))]
    #[case::adds_kind_from_other(AtomConstraints::new(), AtomConstraints::from_iter([AtomConstraint::valence(4)]), Some(AtomConstraints::from_iter([AtomConstraint::valence(4)])))]
    #[case::narrows_undetermined_to_lit(AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Undetermined)]), AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        Some(AtomConstraints::from_iter([AtomConstraint::valence(4)])))]
    #[case::lit_lit_match_preserved(AtomConstraints::from_iter([AtomConstraint::valence(4)]), AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        Some(AtomConstraints::from_iter([AtomConstraint::valence(4)])))]
    #[case::lit_lit_mismatch_none(AtomConstraints::from_iter([AtomConstraint::valence(4)]), AtomConstraints::from_iter([AtomConstraint::valence(3)]), None)]
    #[case::multi_kind_combines(AtomConstraints::from_iter([AtomConstraint::valence(4)]), AtomConstraints::from_iter([AtomConstraint::degree(3)]),
        Some(AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(3)])))]
    #[case::aromatic_valence_narrows(AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::Undetermined)]),
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1))]),
        Some(AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1))])))]
    #[case::aromatic_valence_not_vs_aromatic_none(AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic)]),
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1))]), None)]
    #[case::ring_membership_size_unions(AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::Size(5), 1)]), AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::Size(6), 1)]),
        Some(AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::Size(5), 1), AtomConstraint::ring_membership(RingScope::Size(6), 1)])))]
    #[case::ring_membership_size_dedup(AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::Size(5), 1)]), AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::Size(5), 1)]),
        Some(AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::Size(5), 1)])))]
    #[case::prunes_vacuous(AtomConstraints::new(), AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Undetermined)]), Some(AtomConstraints::new()))]
    #[case::tetrahedral_narrows_from_absent(AtomConstraints::new(),
        AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        Some(AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo)])))]
    #[case::tetrahedral_not_stereo_vs_stereo_contradicts(AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::stereo(0_u32))]), None)]
    fn test_atom_constraints_meet(
        #[case] a: AtomConstraints,
        #[case] b: AtomConstraints,
        #[case] expected: Option<AtomConstraints>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::extends_self(AtomConstraints::new(), AtomConstraints::from_iter([AtomConstraint::valence(4)]), true, AtomConstraints::from_iter([AtomConstraint::valence(4)]))]
    #[case::no_change(AtomConstraints::from_iter([AtomConstraint::valence(4)]), AtomConstraints::from_iter([AtomConstraint::valence(4)]), false,
        AtomConstraints::from_iter([AtomConstraint::valence(4)]))]
    #[case::contradiction_leaves_self_unchanged(AtomConstraints::from_iter([AtomConstraint::valence(4)]), AtomConstraints::from_iter([AtomConstraint::valence(3)]), false,
        AtomConstraints::from_iter([AtomConstraint::valence(4)]))]
    fn test_atom_constraints_narrow_from(
        #[case] mut target: AtomConstraints,
        #[case] source: AtomConstraints,
        #[case] expected_changed: bool,
        #[case] expected_after: AtomConstraints,
    ) {
        let changed = target.narrow_from(&source);
        assert_eq!(changed, expected_changed);
        assert_eq!(target, expected_after);
    }

    #[rstest]
    #[case::keeps_only_shared_kinds(AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(2)]), AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::from_iter([AtomConstraint::valence(4)]))]
    #[case::widens_value(AtomConstraints::from_iter([AtomConstraint::valence(4)]), AtomConstraints::from_iter([AtomConstraint::valence(3)]),
        AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::lit_set([4, 3]))]))]
    #[case::tetrahedral_same(AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]))]
    #[case::tetrahedral_incompatible_drops_to_undetermined(AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::stereo(0_u32))]), AtomConstraints::new())]
    fn test_atom_constraints_join(
        #[case] a: AtomConstraints,
        #[case] b: AtomConstraints,
        #[case] expected: AtomConstraints,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rstest]
    #[case::empty_pattern_matches_anything(AtomConstraints::new(), AtomConstraints::from_iter([AtomConstraint::valence(4)]), true)]
    #[case::missing_in_target_when_pattern_specific(AtomConstraints::from_iter([AtomConstraint::valence(4)]), AtomConstraints::new(), false)]
    #[case::same_lit(AtomConstraints::from_iter([AtomConstraint::valence(4)]), AtomConstraints::from_iter([AtomConstraint::valence(4)]), true)]
    #[case::lit_lit_mismatch(AtomConstraints::from_iter([AtomConstraint::valence(4)]), AtomConstraints::from_iter([AtomConstraint::valence(3)]), false)]
    #[case::aromatic_wildcard_matches_aromatic(AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::Undetermined)]),
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1))]), true)]
    #[case::aromatic_not_vs_aromatic_mismatch(AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic)]),
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1))]), false)]
    #[case::ring_membership_size_subset(AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::Size(5), 1)]),
        AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::Size(5), 1), AtomConstraint::ring_membership(RingScope::Size(6), 1)]), true)]
    #[case::ring_membership_size_not_present_in_target(AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::Size(7), 1)]),
        AtomConstraints::from_iter([AtomConstraint::ring_membership(RingScope::Size(5), 1), AtomConstraint::ring_membership(RingScope::Size(6), 1)]), false)]
    #[case::multi_kind_all_must_match(AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(3)]),
        AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(2)]), false)]
    #[case::tetrahedral_same(AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]), true)]
    #[case::tetrahedral_pattern_specific_vs_absent(AtomConstraints::from_iter([AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]),
        AtomConstraints::new(), false)]
    fn test_atom_constraints_matches(
        #[case] pattern: AtomConstraints,
        #[case] target: AtomConstraints,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
