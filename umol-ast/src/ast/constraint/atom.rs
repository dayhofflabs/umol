//! Atom constraints.

use std::mem::{self, replace};

use smallvec::SmallVec;
use strum::{EnumCount, EnumDiscriminants, EnumIter};

use super::super::remap::IdRemapping;
use super::super::traits::{AsLit, Lattice};
use super::super::value::ValueAst;
use super::joint_domain::JointDomainAst;

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
    RingCount(ValueAst),
    RingSize(ValueAst),
    JointDomain(JointDomainAst),
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

    pub fn ring_count(v: impl Into<ValueAst>) -> Self {
        Self::RingCount(v.into())
    }

    pub fn ring_size(v: impl Into<ValueAst>) -> Self {
        Self::RingSize(v.into())
    }

    pub fn joint_domain(jd: JointDomainAst) -> Self {
        Self::JointDomain(jd)
    }

    pub fn kind(&self) -> AtomConstraintKind {
        self.into()
    }

    /// `false` for variants that may legitimately appear multiple times on
    /// the same atom: `RingSize` (an atom in fused rings satisfies multiple
    /// ring-size assertions) and `JointDomainAst` (each entry constrains a
    /// distinct var-tuple). `true` for variants that are single-valued.
    pub fn is_unique(&self) -> bool {
        !matches!(
            self.kind(),
            AtomConstraintKind::RingSize | AtomConstraintKind::JointDomain
        )
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
            | Self::TotalHydrogens(v)
            | Self::RingCount(v)
            | Self::RingSize(v) => v.is_undetermined(),
            Self::AromaticValence(c) => c.is_undetermined(),
            Self::MulticenterValence(c) => c.is_undetermined(),
            Self::JointDomain(jd) => jd.is_undetermined(),
        }
    }

    /// Recursively simplify the contained value. The constraint kind is
    /// preserved.
    pub fn simplify(self) -> Self {
        match self {
            Self::Valence(v) => Self::Valence(v.simplify()),
            Self::TotalValence(v) => Self::TotalValence(v.simplify()),
            Self::AromaticValence(c) => Self::AromaticValence(c.simplify()),
            Self::MulticenterValence(c) => Self::MulticenterValence(c.simplify()),
            Self::DonatedPairs(v) => Self::DonatedPairs(v.simplify()),
            Self::AcceptedPairs(v) => Self::AcceptedPairs(v.simplify()),
            Self::Degree(v) => Self::Degree(v.simplify()),
            Self::TotalDegree(v) => Self::TotalDegree(v.simplify()),
            Self::RingDegree(v) => Self::RingDegree(v.simplify()),
            Self::RingValence(v) => Self::RingValence(v.simplify()),
            Self::TotalHydrogens(v) => Self::TotalHydrogens(v.simplify()),
            Self::RingCount(v) => Self::RingCount(v.simplify()),
            Self::RingSize(v) => Self::RingSize(v.simplify()),
            Self::JointDomain(jd) => Self::JointDomain(jd.simplify()),
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

    /// Simplify the inner `ValueAst` of `Aromatic(_)`. Other variants are
    /// already canonical.
    pub fn simplify(self) -> Self {
        match self {
            Self::Aromatic(v) => Self::Aromatic(v.simplify()),
            other => other,
        }
    }

    /// Pattern matches value.
    pub fn matches_value(&self, value: i64) -> bool {
        match self {
            Self::Aromatic(v) => v.matches_value(value),
            Self::NotAromatic => value == 0,
            Self::Undetermined => true,
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

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::NotAromatic, Self::NotAromatic) => Some(Self::NotAromatic),
            (Self::NotAromatic, Self::Aromatic(_)) | (Self::Aromatic(_), Self::NotAromatic) => None,
            (Self::Aromatic(a), Self::Aromatic(b)) => a.meet(b).map(Self::Aromatic),
        }
    }

    fn join(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::NotAromatic, Self::NotAromatic) => Self::NotAromatic,
            (Self::NotAromatic, Self::Aromatic(_)) | (Self::Aromatic(_), Self::NotAromatic) => {
                Self::Undetermined
            }
            (Self::Aromatic(a), Self::Aromatic(b)) => Self::Aromatic(a.join(b)),
        }
    }

    /// Pattern matches target. `Undetermined` is a wildcard pattern.
    /// `NotAromatic` and `Aromatic(_)` are mutually exclusive; an `Aromatic`
    /// pattern recurses on the inner `ValueAst::matches`.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::NotAromatic, Self::NotAromatic) => true,
            (Self::NotAromatic, Self::Aromatic(_)) | (Self::Aromatic(_), Self::NotAromatic) => {
                false
            }
            (Self::Aromatic(p), Self::Aromatic(t)) => p.matches(t),
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
    pub fn multicenter(v: impl Into<ValueAst>) -> Self {
        Self::Multicenter(v.into())
    }

    pub fn is_multicenter(&self) -> bool {
        matches!(self, Self::Multicenter(_))
    }

    /// Simplify the inner `ValueAst` of `Multicenter(_)`. Other variants
    /// are already canonical.
    pub fn simplify(self) -> Self {
        match self {
            Self::Multicenter(v) => Self::Multicenter(v.simplify()),
            other => other,
        }
    }

    /// Pattern matches value.
    pub fn matches_value(&self, value: i64) -> bool {
        match self {
            Self::Multicenter(v) => v.matches_value(value),
            Self::NotMulticenter => value == 0,
            Self::Undetermined => true,
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

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::NotMulticenter, Self::NotMulticenter) => Some(Self::NotMulticenter),
            (Self::NotMulticenter, Self::Multicenter(_))
            | (Self::Multicenter(_), Self::NotMulticenter) => None,
            (Self::Multicenter(a), Self::Multicenter(b)) => a.meet(b).map(Self::Multicenter),
        }
    }

    fn join(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::NotMulticenter, Self::NotMulticenter) => Self::NotMulticenter,
            (Self::NotMulticenter, Self::Multicenter(_))
            | (Self::Multicenter(_), Self::NotMulticenter) => Self::Undetermined,
            (Self::Multicenter(a), Self::Multicenter(b)) => Self::Multicenter(a.join(b)),
        }
    }

    /// Pattern matches target. `Undetermined` is a wildcard pattern.
    /// `NotMulticenter` and `Multicenter(_)` are mutually exclusive; a
    /// `Multicenter` pattern recurses on the inner `ValueAst::matches`.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::NotMulticenter, Self::NotMulticenter) => true,
            (Self::NotMulticenter, Self::Multicenter(_))
            | (Self::Multicenter(_), Self::NotMulticenter) => false,
            (Self::Multicenter(p), Self::Multicenter(t)) => p.matches(t),
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

/// Per-atom constraints: at most one entry per [`AtomConstraintKind`].
/// Stored kind-sorted in an inline-capacity-2 `SmallVec`; the common cases
/// after resolution (0–2 constraints) never touch the heap.
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

    pub fn valence(&self) -> ValueAst {
        match self.get(AtomConstraintKind::Valence) {
            Some(AtomConstraint::Valence(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn total_valence(&self) -> ValueAst {
        match self.get(AtomConstraintKind::TotalValence) {
            Some(AtomConstraint::TotalValence(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn degree(&self) -> ValueAst {
        match self.get(AtomConstraintKind::Degree) {
            Some(AtomConstraint::Degree(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn total_degree(&self) -> ValueAst {
        match self.get(AtomConstraintKind::TotalDegree) {
            Some(AtomConstraint::TotalDegree(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn ring_degree(&self) -> ValueAst {
        match self.get(AtomConstraintKind::RingDegree) {
            Some(AtomConstraint::RingDegree(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn ring_valence(&self) -> ValueAst {
        match self.get(AtomConstraintKind::RingValence) {
            Some(AtomConstraint::RingValence(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn total_hydrogens(&self) -> ValueAst {
        match self.get(AtomConstraintKind::TotalHydrogens) {
            Some(AtomConstraint::TotalHydrogens(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn donated_pairs(&self) -> ValueAst {
        match self.get(AtomConstraintKind::DonatedPairs) {
            Some(AtomConstraint::DonatedPairs(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn accepted_pairs(&self) -> ValueAst {
        match self.get(AtomConstraintKind::AcceptedPairs) {
            Some(AtomConstraint::AcceptedPairs(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn ring_count(&self) -> ValueAst {
        match self.get(AtomConstraintKind::RingCount) {
            Some(AtomConstraint::RingCount(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn aromatic_valence(&self) -> AromaticValenceAst {
        match self.get(AtomConstraintKind::AromaticValence) {
            Some(AtomConstraint::AromaticValence(v)) => v.clone(),
            _ => AromaticValenceAst::Undetermined,
        }
    }

    pub fn multicenter_valence(&self) -> MulticenterValenceAst {
        match self.get(AtomConstraintKind::MulticenterValence) {
            Some(AtomConstraint::MulticenterValence(v)) => v.clone(),
            _ => MulticenterValenceAst::Undetermined,
        }
    }

    /// Multi-valued ring-size assertions; an atom in fused rings may carry
    /// several. Iterator yields entries in store order; empty if none.
    pub fn ring_sizes(&self) -> impl Iterator<Item = &ValueAst> {
        self.get_all(AtomConstraintKind::RingSize)
            .filter_map(|c| match c {
                AtomConstraint::RingSize(v) => Some(v),
                _ => None,
            })
    }

    /// Joint-domain (relational) constraints. Multiple entries may coexist —
    /// each constrains a distinct (or overlapping) tuple of atom-level vars.
    pub fn joint_domains(&self) -> impl Iterator<Item = &JointDomainAst> {
        self.get_all(AtomConstraintKind::JointDomain)
            .filter_map(|c| match c {
                AtomConstraint::JointDomain(jd) => Some(jd),
                _ => None,
            })
    }

    /// Insert a constraint per the per-variant cardinality policy. Single-
    /// valued kinds (`AtomConstraint::is_unique` → true) replace any
    /// existing entry of the same kind, last-wins; multi-valued kinds append
    /// after the existing cluster of that kind, preserving the sorted-by-
    /// kind invariant. Returns the replaced entry if a unique-kind same-kind
    /// entry existed; `None` for the append path.
    pub fn add(&mut self, constraint: AtomConstraint) -> Option<AtomConstraint> {
        match self.find(constraint.kind()) {
            Ok(i) => {
                if constraint.is_unique() {
                    Some(replace(&mut self.entries[i], constraint))
                } else {
                    let after = self.entries[i..]
                        .iter()
                        .take_while(|c| c.kind() == constraint.kind())
                        .count();
                    self.entries.insert(i + after, constraint);
                    None
                }
            }
            Err(i) => {
                self.entries.insert(i, constraint);
                None
            }
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

    /// Simplify each contained constraint's value in place. Kind is
    /// preserved by `AtomConstraint::simplify`, so the sorted-by-kind
    /// invariant holds without re-sorting.
    pub fn simplify_each(&mut self) {
        for c in self.entries.iter_mut() {
            *c = mem::replace(c, AtomConstraint::Valence(ValueAst::Undetermined)).simplify();
        }
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
    /// one entry; multi-valued (`RingSize`) may yield several.
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
            | AtomConstraint::TotalHydrogens(v)
            | AtomConstraint::RingCount(v)
            | AtomConstraint::RingSize(v) => v.is_ground(),
            AtomConstraint::AromaticValence(c) => c.is_ground(),
            AtomConstraint::MulticenterValence(c) => c.is_ground(),
            AtomConstraint::JointDomain(jd) => jd.is_ground(),
        })
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let mut result = Self::new();
        let v = self.valence().meet(&other.valence())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::Valence(v));
        }
        let v = self.total_valence().meet(&other.total_valence())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::TotalValence(v));
        }
        let v = self.aromatic_valence().meet(&other.aromatic_valence())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::AromaticValence(v));
        }
        let v = self
            .multicenter_valence()
            .meet(&other.multicenter_valence())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::MulticenterValence(v));
        }
        let v = self.donated_pairs().meet(&other.donated_pairs())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::DonatedPairs(v));
        }
        let v = self.accepted_pairs().meet(&other.accepted_pairs())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::AcceptedPairs(v));
        }
        let v = self.degree().meet(&other.degree())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::Degree(v));
        }
        let v = self.total_degree().meet(&other.total_degree())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::TotalDegree(v));
        }
        let v = self.ring_degree().meet(&other.ring_degree())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::RingDegree(v));
        }
        let v = self.ring_valence().meet(&other.ring_valence())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::RingValence(v));
        }
        let v = self.total_hydrogens().meet(&other.total_hydrogens())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::TotalHydrogens(v));
        }
        let v = self.ring_count().meet(&other.ring_count())?;
        if !v.is_undetermined() {
            result.add(AtomConstraint::RingCount(v));
        }
        for v in self.ring_sizes().chain(other.ring_sizes()) {
            if v.is_undetermined() {
                continue;
            }
            let entry = AtomConstraint::RingSize(v.clone());
            if !result.contains_entry(&entry) {
                result.add(entry);
            }
        }
        for jd in self.joint_domains().chain(other.joint_domains()) {
            let entry = AtomConstraint::JointDomain(jd.clone());
            if !result.contains_entry(&entry) {
                result.add(entry);
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
                    let joined = self.$accessor().join(&other.$accessor());
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
        join_unique_value!(RingCount, ring_count, RingCount);
        for v in self.ring_sizes() {
            if v.is_undetermined() {
                continue;
            }
            let entry = AtomConstraint::RingSize(v.clone());
            if other
                .ring_sizes()
                .any(|o| AtomConstraint::RingSize(o.clone()) == entry)
            {
                result.add(entry);
            }
        }
        for jd in self.joint_domains() {
            if other.joint_domains().any(|o| o == jd) {
                result.add(AtomConstraint::JointDomain(jd.clone()));
            }
        }
        result
    }

    /// Field-wise per-kind: single-valued kinds match via the corresponding
    /// `Lattice::matches`; `RingSize` (multi-valued) requires every `self`
    /// assertion to be matchable by some `target` assertion.
    fn matches(&self, target: &Self) -> bool {
        self.valence().matches(&target.valence())
            && self.total_valence().matches(&target.total_valence())
            && self.aromatic_valence().matches(&target.aromatic_valence())
            && self
                .multicenter_valence()
                .matches(&target.multicenter_valence())
            && self.donated_pairs().matches(&target.donated_pairs())
            && self.accepted_pairs().matches(&target.accepted_pairs())
            && self.degree().matches(&target.degree())
            && self.total_degree().matches(&target.total_degree())
            && self.ring_degree().matches(&target.ring_degree())
            && self.ring_valence().matches(&target.ring_valence())
            && self.total_hydrogens().matches(&target.total_hydrogens())
            && self.ring_count().matches(&target.ring_count())
            && self
                .ring_sizes()
                .all(|p| target.ring_sizes().any(|t| p.matches(t)))
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

    use super::*;
    use crate::ast::constraint::joint_domain::{JointDomainAst, JointVar};
    use crate::ast::value::Expr;

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
    #[case::ring_count(AtomConstraint::ring_count(1), AtomConstraint::RingCount(ValueAst::Lit(1)))]
    #[case::ring_size(AtomConstraint::ring_size(6), AtomConstraint::RingSize(ValueAst::Lit(6)))]
    #[case::aromatic_valence(
        AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic),
        AtomConstraint::AromaticValence(AromaticValenceAst::NotAromatic),
    )]
    #[case::multicenter_valence(
        AtomConstraint::multicenter_valence(MulticenterValenceAst::NotMulticenter),
        AtomConstraint::MulticenterValence(MulticenterValenceAst::NotMulticenter),
    )]
    #[case::joint_domain(
        AtomConstraint::joint_domain(
            JointDomainAst::from_ints(
                vec![JointVar::LonePairs, JointVar::UnpairedElectrons],
                vec![vec![3, 0], vec![1, 4]],
            ).unwrap(),
        ),
        AtomConstraint::JointDomain(
            JointDomainAst::from_ints(
                vec![JointVar::LonePairs, JointVar::UnpairedElectrons],
                vec![vec![3, 0], vec![1, 4]],
            ).unwrap(),
        ),
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
    #[case::ring_count(AtomConstraint::ring_count(1), AtomConstraintKind::RingCount)]
    #[case::ring_size(AtomConstraint::ring_size(6), AtomConstraintKind::RingSize)]
    #[case::joint_domain(
        AtomConstraint::joint_domain(JointDomainAst::from_ints(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens],
            vec![vec![0, 1], vec![1, 2]],
        ).unwrap()),
        AtomConstraintKind::JointDomain,
    )]
    fn test_atom_constraint_kind(
        #[case] constraint: AtomConstraint,
        #[case] expected: AtomConstraintKind,
    ) {
        assert_eq!(constraint.kind(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_lit(AtomConstraint::valence(4), false)]
    #[case::valence_undetermined(AtomConstraint::Valence(ValueAst::Undetermined), true)]
    #[case::degree_undetermined(AtomConstraint::Degree(ValueAst::Undetermined), true)]
    #[case::ring_size_undetermined(AtomConstraint::RingSize(ValueAst::Undetermined), true)]
    #[case::aromatic_undetermined(AtomConstraint::aromatic_valence(AromaticValenceAst::Undetermined), true)]
    #[case::aromatic_not_aromatic(AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic), false)]
    #[case::aromatic_with_value(AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1)), false)]
    #[case::multicenter_undetermined(AtomConstraint::multicenter_valence(MulticenterValenceAst::Undetermined), true)]
    #[case::multicenter_not(AtomConstraint::multicenter_valence(MulticenterValenceAst::NotMulticenter), false)]
    #[case::multicenter_with_value(AtomConstraint::multicenter_valence(MulticenterValenceAst::multicenter(1)), false)]
    #[case::joint_domain(AtomConstraint::joint_domain(JointDomainAst::from_ints(
        vec![JointVar::Charge, JointVar::ImplicitHydrogens],
        vec![vec![0, 1], vec![1, 2]],
    ).unwrap()), false)]
    fn test_atom_constraint_is_undetermined(
        #[case] c: AtomConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_folds_expr(
        AtomConstraint::Valence(ValueAst::Expr(Box::new(Expr::Lit(4)))),
        AtomConstraint::valence(4),
    )]
    #[case::degree_folds_expr(
        AtomConstraint::Degree(ValueAst::Expr(Box::new(Expr::Lit(3)))),
        AtomConstraint::degree(3),
    )]
    #[case::aromatic_valence_folds_inner(
        AtomConstraint::aromatic_valence(AromaticValenceAst::Aromatic(ValueAst::Expr(Box::new(Expr::Lit(2))))),
        AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(2)),
    )]
    #[case::multicenter_valence_folds_inner(
        AtomConstraint::multicenter_valence(MulticenterValenceAst::Multicenter(ValueAst::Expr(Box::new(Expr::Lit(3))))),
        AtomConstraint::multicenter_valence(MulticenterValenceAst::multicenter(3)),
    )]
    fn test_atom_constraint_simplify(
        #[case] input: AtomConstraint,
        #[case] expected: AtomConstraint,
    ) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::valence_lit(AtomConstraint::valence(4))]
    #[case::aromatic_not_aromatic(AtomConstraint::aromatic_valence(
        AromaticValenceAst::NotAromatic
    ))]
    #[case::multicenter_undetermined(AtomConstraint::multicenter_valence(
        MulticenterValenceAst::Undetermined
    ))]
    #[case::joint_domain(AtomConstraint::joint_domain(JointDomainAst::from_ints(
        vec![JointVar::Charge, JointVar::ImplicitHydrogens],
        vec![vec![0, 1], vec![1, 2]],
    ).unwrap()))]
    fn test_atom_constraint_simplify_identity(#[case] input: AtomConstraint) {
        assert_eq!(input.clone().simplify(), input);
    }

    #[rstest]
    fn test_atom_constraint_is_unique_joint_domain() {
        let jd = AtomConstraint::joint_domain(
            JointDomainAst::from_ints(
                vec![JointVar::Charge, JointVar::ImplicitHydrogens],
                vec![vec![0, 1], vec![1, 2]],
            )
            .unwrap(),
        );
        assert!(!jd.is_unique());
    }

    #[rstest]
    fn test_atom_constraints_joint_domains_accessor() {
        let jd1 = JointDomainAst::from_ints(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens],
            vec![vec![0, 1], vec![1, 2]],
        )
        .unwrap();
        let jd2 = JointDomainAst::from_ints(
            vec![JointVar::LonePairs, JointVar::UnpairedElectrons],
            vec![vec![3, 0], vec![1, 4]],
        )
        .unwrap();
        let mut cs = AtomConstraints::new();
        cs.add(AtomConstraint::joint_domain(jd1.clone()));
        cs.add(AtomConstraint::joint_domain(jd2.clone()));
        let collected: Vec<&JointDomainAst> = cs.joint_domains().collect();
        assert_eq!(collected, vec![&jd1, &jd2]);
    }

    #[rstest]
    fn test_atom_constraints_is_ground_false_with_joint_domain() {
        let mut cs = AtomConstraints::new();
        cs.add(AtomConstraint::Valence(ValueAst::Lit(4)));
        assert!(cs.is_ground());
        cs.add(AtomConstraint::joint_domain(
            JointDomainAst::from_ints(
                vec![JointVar::Charge, JointVar::ImplicitHydrogens],
                vec![vec![0, 1], vec![1, 2]],
            )
            .unwrap(),
        ));
        assert!(!cs.is_ground());
    }

    #[rstest]
    #[case::lit(
        AromaticValenceAst::aromatic(1),
        AromaticValenceAst::Aromatic(ValueAst::Lit(1))
    )]
    fn test_aromatic_valence_ast_aromatic(
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
    #[case::aromatic_expr_folds(
        AromaticValenceAst::Aromatic(ValueAst::Expr(Box::new(Expr::Lit(2)))),
        Some(2)
    )]
    fn test_aromatic_valence_ast_as_lit(
        #[case] v: AromaticValenceAst,
        #[case] expected: Option<i64>,
    ) {
        assert_eq!(v.as_lit(), expected);
    }

    #[rstest]
    #[case::aromatic_folds_expr(
        AromaticValenceAst::Aromatic(ValueAst::Expr(Box::new(Expr::Lit(2)))),
        AromaticValenceAst::aromatic(2)
    )]
    fn test_aromatic_valence_ast_simplify(
        #[case] input: AromaticValenceAst,
        #[case] expected: AromaticValenceAst,
    ) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic)]
    #[case::aromatic_lit(AromaticValenceAst::aromatic(1))]
    fn test_aromatic_valence_ast_simplify_identity(#[case] input: AromaticValenceAst) {
        assert_eq!(input.clone().simplify(), input);
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

    #[rstest]
    #[case::lit(
        MulticenterValenceAst::multicenter(2),
        MulticenterValenceAst::Multicenter(ValueAst::Lit(2))
    )]
    fn test_multicenter_valence_ast_multicenter(
        #[case] actual: MulticenterValenceAst,
        #[case] expected: MulticenterValenceAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined, false)]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, false)]
    #[case::multicenter_undetermined(
        MulticenterValenceAst::Multicenter(ValueAst::Undetermined),
        true
    )]
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

    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined, None)]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, Some(0))]
    #[case::multicenter_undetermined(
        MulticenterValenceAst::Multicenter(ValueAst::Undetermined),
        None
    )]
    #[case::multicenter_lit(MulticenterValenceAst::multicenter(2), Some(2))]
    #[case::multicenter_expr_folds(
        MulticenterValenceAst::Multicenter(ValueAst::Expr(Box::new(Expr::Lit(3)))),
        Some(3)
    )]
    fn test_multicenter_valence_ast_as_lit(
        #[case] v: MulticenterValenceAst,
        #[case] expected: Option<i64>,
    ) {
        assert_eq!(v.as_lit(), expected);
    }

    #[rstest]
    #[case::multicenter_folds_expr(
        MulticenterValenceAst::Multicenter(ValueAst::Expr(Box::new(Expr::Lit(3)))),
        MulticenterValenceAst::multicenter(3)
    )]
    fn test_multicenter_valence_ast_simplify(
        #[case] input: MulticenterValenceAst,
        #[case] expected: MulticenterValenceAst,
    ) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined)]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter)]
    #[case::multicenter_lit(MulticenterValenceAst::multicenter(1))]
    fn test_multicenter_valence_ast_simplify_identity(#[case] input: MulticenterValenceAst) {
        assert_eq!(input.clone().simplify(), input);
    }

    #[rstest]
    #[case::wildcard_vs_not_multicenter(
        MulticenterValenceAst::Undetermined,
        MulticenterValenceAst::NotMulticenter,
        true
    )]
    #[case::wildcard_vs_multicenter_lit(
        MulticenterValenceAst::Undetermined,
        MulticenterValenceAst::multicenter(2),
        true
    )]
    #[case::not_multicenter_vs_multicenter(
        MulticenterValenceAst::NotMulticenter,
        MulticenterValenceAst::multicenter(2),
        false
    )]
    #[case::multicenter_vs_not_multicenter(
        MulticenterValenceAst::multicenter(2),
        MulticenterValenceAst::NotMulticenter,
        false
    )]
    #[case::not_multicenter_vs_not_multicenter(
        MulticenterValenceAst::NotMulticenter,
        MulticenterValenceAst::NotMulticenter,
        true
    )]
    #[case::multicenter_lit_match(
        MulticenterValenceAst::multicenter(2),
        MulticenterValenceAst::multicenter(2),
        true
    )]
    #[case::multicenter_lit_mismatch(
        MulticenterValenceAst::multicenter(2),
        MulticenterValenceAst::multicenter(3),
        false
    )]
    #[case::specific_vs_undetermined(
        MulticenterValenceAst::multicenter(2),
        MulticenterValenceAst::Undetermined,
        false
    )]
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
    #[case::aromatic_present(AtomConstraintKind::AromaticValence,
        Some(AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic)))]
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
    #[case::fresh(
        vec![AtomConstraint::valence(4)],
        vec![None],
        vec![AtomConstraint::valence(4)],
    )]
    #[case::replace_same_kind(
        vec![AtomConstraint::valence(3), AtomConstraint::valence(4)],
        vec![None, Some(AtomConstraint::valence(3))],
        vec![AtomConstraint::valence(4)],
    )]
    #[case::distinct_kinds(
        vec![
            AtomConstraint::valence(4),
            AtomConstraint::degree(3),
            AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic),
        ],
        vec![None, None, None],
        vec![
            AtomConstraint::valence(4),
            AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic),
            AtomConstraint::degree(3),
        ],
    )]
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
    #[case::partial(
        |c: &AtomConstraint| matches!(c, AtomConstraint::Valence(_) | AtomConstraint::RingCount(_)),
        vec![AtomConstraint::valence(4), AtomConstraint::ring_count(2)],
    )]
    #[case::all_dropped(|_: &AtomConstraint| false, vec![])]
    fn test_atom_constraints_retain(
        #[case] predicate: impl FnMut(&AtomConstraint) -> bool,
        #[case] expected: Vec<AtomConstraint>,
    ) {
        let mut cs = AtomConstraints::from_iter([
            AtomConstraint::valence(4),
            AtomConstraint::degree(3),
            AtomConstraint::ring_count(2),
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

    #[rstest]
    fn test_atom_constraints_simplify_each() {
        let mut cs = AtomConstraints::from_iter([
            AtomConstraint::Valence(ValueAst::Expr(Box::new(Expr::Lit(4)))),
            AtomConstraint::Degree(ValueAst::Expr(Box::new(Expr::Lit(3)))),
            AtomConstraint::aromatic_valence(AromaticValenceAst::Aromatic(ValueAst::Expr(
                Box::new(Expr::Lit(2)),
            ))),
        ]);
        cs.simplify_each();
        assert_eq!(
            cs,
            AtomConstraints::from_iter([
                AtomConstraint::valence(4),
                AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(2)),
                AtomConstraint::degree(3),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_present(
        AtomConstraintKind::Valence,
        Some(AtomConstraint::valence(4)),
        vec![AtomConstraint::degree(3)],
    )]
    #[case::degree_present(
        AtomConstraintKind::Degree,
        Some(AtomConstraint::degree(3)),
        vec![AtomConstraint::valence(4)],
    )]
    #[case::absent(
        AtomConstraintKind::RingCount,
        None,
        vec![AtomConstraint::valence(4), AtomConstraint::degree(3)],
    )]
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
            AtomConstraint::ring_size(6),
            AtomConstraint::valence(4),
            AtomConstraint::degree(3),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                AtomConstraint::valence(4),
                AtomConstraint::degree(3),
                AtomConstraint::ring_size(6),
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
            umol_graph_core::Remapping {
                removed_nodes: vec![0, 1, 2],
                removed_edges: vec![0],
            },
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(cs.clone().remap(&remap), cs);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::distinct(
        vec![AtomConstraint::valence(4), AtomConstraint::degree(3)],
        vec![AtomConstraint::valence(4), AtomConstraint::degree(3)],
    )]
    #[case::same_kind_last_wins(
        vec![AtomConstraint::valence(3), AtomConstraint::valence(4)],
        vec![AtomConstraint::valence(4)],
    )]
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

    #[rstest]
    #[case::empty_empty(
        AtomConstraints::new(),
        AtomConstraints::new(),
        Some(AtomConstraints::new())
    )]
    #[case::adds_kind_from_other(
        AtomConstraints::new(),
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        Some(AtomConstraints::from_iter([AtomConstraint::valence(4)])),
    )]
    #[case::narrows_undetermined_to_lit(
        AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Undetermined)]),
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        Some(AtomConstraints::from_iter([AtomConstraint::valence(4)])),
    )]
    #[case::lit_lit_match_preserved(
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        Some(AtomConstraints::from_iter([AtomConstraint::valence(4)])),
    )]
    #[case::lit_lit_mismatch_none(
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::from_iter([AtomConstraint::valence(3)]),
        None,
    )]
    #[case::multi_kind_combines(
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::from_iter([AtomConstraint::degree(3)]),
        Some(AtomConstraints::from_iter([
            AtomConstraint::valence(4),
            AtomConstraint::degree(3),
        ])),
    )]
    #[case::aromatic_valence_narrows(
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::Undetermined)]),
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1))]),
        Some(AtomConstraints::from_iter([
            AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1)),
        ])),
    )]
    #[case::aromatic_valence_not_vs_aromatic_none(
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic)]),
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1))]),
        None,
    )]
    #[case::ring_size_unions(
        AtomConstraints::from_iter([AtomConstraint::ring_size(5)]),
        AtomConstraints::from_iter([AtomConstraint::ring_size(6)]),
        Some(AtomConstraints::from_iter([
            AtomConstraint::ring_size(5),
            AtomConstraint::ring_size(6),
        ])),
    )]
    #[case::ring_size_dedup(
        AtomConstraints::from_iter([AtomConstraint::ring_size(5)]),
        AtomConstraints::from_iter([AtomConstraint::ring_size(5)]),
        Some(AtomConstraints::from_iter([AtomConstraint::ring_size(5)])),
    )]
    #[case::prunes_vacuous(
        AtomConstraints::new(),
        AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Undetermined)]),
        Some(AtomConstraints::new()),
    )]
    fn test_atom_constraints_meet(
        #[case] a: AtomConstraints,
        #[case] b: AtomConstraints,
        #[case] expected: Option<AtomConstraints>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::extends_self(
        AtomConstraints::new(),
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        true,
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
    )]
    #[case::no_change(
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        false,
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
    )]
    #[case::contradiction_leaves_self_unchanged(
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::from_iter([AtomConstraint::valence(3)]),
        false,
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
    )]
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
    #[case::keeps_only_shared_kinds(
        AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(2)]),
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
    )]
    #[case::widens_value(
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::from_iter([AtomConstraint::valence(3)]),
        AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Set(Box::new(vec![4, 3])))]),
    )]
    fn test_atom_constraints_join(
        #[case] a: AtomConstraints,
        #[case] b: AtomConstraints,
        #[case] expected: AtomConstraints,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rstest]
    #[case::empty_pattern_matches_anything(
        AtomConstraints::new(),
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        true,
    )]
    #[case::missing_in_target_when_pattern_specific(
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::new(),
        false,
    )]
    #[case::same_lit(
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        true,
    )]
    #[case::lit_lit_mismatch(
        AtomConstraints::from_iter([AtomConstraint::valence(4)]),
        AtomConstraints::from_iter([AtomConstraint::valence(3)]),
        false,
    )]
    #[case::aromatic_wildcard_matches_aromatic(
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::Undetermined)]),
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1))]),
        true,
    )]
    #[case::aromatic_not_vs_aromatic_mismatch(
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic)]),
        AtomConstraints::from_iter([AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(1))]),
        false,
    )]
    #[case::ring_size_subset(
        AtomConstraints::from_iter([AtomConstraint::ring_size(5)]),
        AtomConstraints::from_iter([AtomConstraint::ring_size(5), AtomConstraint::ring_size(6)]),
        true,
    )]
    #[case::ring_size_not_present_in_target(
        AtomConstraints::from_iter([AtomConstraint::ring_size(7)]),
        AtomConstraints::from_iter([AtomConstraint::ring_size(5), AtomConstraint::ring_size(6)]),
        false,
    )]
    #[case::multi_kind_all_must_match(
        AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(3)]),
        AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::degree(2)]),
        false,
    )]
    fn test_atom_constraints_matches(
        #[case] pattern: AtomConstraints,
        #[case] target: AtomConstraints,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
