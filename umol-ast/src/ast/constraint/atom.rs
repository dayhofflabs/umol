//! Atom constraints.

use std::mem::{self, replace};

use smallvec::SmallVec;
use strum::{EnumCount, EnumDiscriminants, EnumIter};

use super::super::remap::IdxRemapping;
use super::super::value::ValueAst;

/// Atom-scope constraint: a predicate that pattern-matches a single atom
/// on a topological or valence property (valence, degree, ring membership,
/// etc.). Held inline on `AtomAst` via `AtomConstraints`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, EnumDiscriminants)]
#[strum_discriminants(name(AtomConstraintKind), derive(Hash, EnumCount, EnumIter))]
#[repr(u8)]
pub enum AtomConstraint {
    Valence(ValueAst),
    AromaticValence(AromaticValenceAst),
    MulticenterValence(MulticenterValenceAst),
    DonatedPairs(ValueAst),
    AcceptedPairs(ValueAst),
    Degree(ValueAst),
    Connectivity(ValueAst),
    RingConnectivity(ValueAst),
    TotalHydrogens(ValueAst),
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl AtomConstraint {
    pub fn valence(v: impl Into<ValueAst>) -> Self {
        Self::Valence(v.into())
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

    pub fn connectivity(v: impl Into<ValueAst>) -> Self {
        Self::Connectivity(v.into())
    }

    pub fn ring_connectivity(v: impl Into<ValueAst>) -> Self {
        Self::RingConnectivity(v.into())
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

    pub fn kind(&self) -> AtomConstraintKind {
        self.into()
    }

    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Valence(v)
            | Self::DonatedPairs(v)
            | Self::AcceptedPairs(v)
            | Self::Degree(v)
            | Self::Connectivity(v)
            | Self::RingConnectivity(v)
            | Self::TotalHydrogens(v)
            | Self::RingCount(v)
            | Self::RingSize(v) => v.is_undetermined(),
            Self::AromaticValence(c) => c.is_undetermined(),
            Self::MulticenterValence(c) => c.is_undetermined(),
        }
    }

    /// Recursively simplify the contained value. The constraint kind is
    /// preserved.
    pub fn simplify(self) -> Self {
        match self {
            Self::Valence(v) => Self::Valence(v.simplify()),
            Self::AromaticValence(c) => Self::AromaticValence(c.simplify()),
            Self::MulticenterValence(c) => Self::MulticenterValence(c.simplify()),
            Self::DonatedPairs(v) => Self::DonatedPairs(v.simplify()),
            Self::AcceptedPairs(v) => Self::AcceptedPairs(v.simplify()),
            Self::Degree(v) => Self::Degree(v.simplify()),
            Self::Connectivity(v) => Self::Connectivity(v.simplify()),
            Self::RingConnectivity(v) => Self::RingConnectivity(v.simplify()),
            Self::TotalHydrogens(v) => Self::TotalHydrogens(v.simplify()),
            Self::RingCount(v) => Self::RingCount(v.simplify()),
            Self::RingSize(v) => Self::RingSize(v.simplify()),
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

    pub fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    /// Simplify the inner `ValueAst` of `Aromatic(_)`. Other variants are
    /// already canonical.
    pub fn simplify(self) -> Self {
        match self {
            Self::Aromatic(v) => Self::Aromatic(v.simplify()),
            other => other,
        }
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

    pub fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    /// Simplify the inner `ValueAst` of `Multicenter(_)`. Other variants
    /// are already canonical.
    pub fn simplify(self) -> Self {
        match self {
            Self::Multicenter(v) => Self::Multicenter(v.simplify()),
            other => other,
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

    /// Insert a constraint at its kind's sorted position, returning the
    /// previous entry of the same kind if any. Every `AtomConstraintKind` is
    /// single-valued per atom, so `add` always replaces same-kind entries
    /// (last-wins).
    pub fn add(&mut self, constraint: AtomConstraint) -> Option<AtomConstraint> {
        match self.find(constraint.kind()) {
            Ok(i) => Some(replace(&mut self.entries[i], constraint)),
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

    pub fn iter(&self) -> impl Iterator<Item = &AtomConstraint> {
        self.entries.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut AtomConstraint> {
        self.entries.iter_mut()
    }

    /// No-op: no `AtomConstraint` variant carries an entity index.
    pub fn remap(self, _remap: &IdxRemapping) -> Self {
        self
    }

    fn find(&self, kind: AtomConstraintKind) -> Result<usize, usize> {
        self.entries
            .binary_search_by_key(&(kind as u8), |c| c.kind() as u8)
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
    use crate::ast::value::Expr;

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraint::valence(4), AtomConstraint::Valence(ValueAst::Lit(4)))]
    #[case::donated_pairs(AtomConstraint::donated_pairs(1), AtomConstraint::DonatedPairs(ValueAst::Lit(1)))]
    #[case::accepted_pairs(AtomConstraint::accepted_pairs(2), AtomConstraint::AcceptedPairs(ValueAst::Lit(2)))]
    #[case::degree(AtomConstraint::degree(3), AtomConstraint::Degree(ValueAst::Lit(3)))]
    #[case::connectivity(AtomConstraint::connectivity(4), AtomConstraint::Connectivity(ValueAst::Lit(4)))]
    #[case::ring_connectivity(AtomConstraint::ring_connectivity(2), AtomConstraint::RingConnectivity(ValueAst::Lit(2)))]
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
    fn test_atom_constraint_constructors(
        #[case] actual: AtomConstraint,
        #[case] expected: AtomConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraint::valence(4), AtomConstraintKind::Valence)]
    #[case::aromatic_valence(AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic), AtomConstraintKind::AromaticValence)]
    #[case::multicenter_valence(AtomConstraint::multicenter_valence(MulticenterValenceAst::Undetermined), AtomConstraintKind::MulticenterValence)]
    #[case::donated_pairs(AtomConstraint::donated_pairs(1), AtomConstraintKind::DonatedPairs)]
    #[case::accepted_pairs(AtomConstraint::accepted_pairs(2), AtomConstraintKind::AcceptedPairs)]
    #[case::degree(AtomConstraint::degree(3), AtomConstraintKind::Degree)]
    #[case::connectivity(AtomConstraint::connectivity(4), AtomConstraintKind::Connectivity)]
    #[case::ring_connectivity(AtomConstraint::ring_connectivity(2), AtomConstraintKind::RingConnectivity)]
    #[case::total_hydrogens(AtomConstraint::total_hydrogens(3), AtomConstraintKind::TotalHydrogens)]
    #[case::ring_count(AtomConstraint::ring_count(1), AtomConstraintKind::RingCount)]
    #[case::ring_size(AtomConstraint::ring_size(6), AtomConstraintKind::RingSize)]
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
    fn test_atom_constraint_is_undetermined(
        #[case] c: AtomConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_folds_expr(
        AtomConstraint::Valence(ValueAst::Expr(Expr::Lit(4))),
        AtomConstraint::valence(4),
    )]
    #[case::degree_folds_expr(
        AtomConstraint::Degree(ValueAst::Expr(Expr::Lit(3))),
        AtomConstraint::degree(3),
    )]
    #[case::aromatic_valence_folds_inner(
        AtomConstraint::aromatic_valence(AromaticValenceAst::Aromatic(ValueAst::Expr(Expr::Lit(2)))),
        AtomConstraint::aromatic_valence(AromaticValenceAst::aromatic(2)),
    )]
    #[case::multicenter_valence_folds_inner(
        AtomConstraint::multicenter_valence(MulticenterValenceAst::Multicenter(ValueAst::Expr(Expr::Lit(3)))),
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
    fn test_atom_constraint_simplify_identity(#[case] input: AtomConstraint) {
        assert_eq!(input.clone().simplify(), input);
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
    #[case::aromatic_folds_expr(
        AromaticValenceAst::Aromatic(ValueAst::Expr(Expr::Lit(2))),
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
    #[case::multicenter_folds_expr(
        MulticenterValenceAst::Multicenter(ValueAst::Expr(Expr::Lit(3))),
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
            AtomConstraint::Valence(ValueAst::Expr(Expr::Lit(4))),
            AtomConstraint::Degree(ValueAst::Expr(Expr::Lit(3))),
            AtomConstraint::aromatic_valence(AromaticValenceAst::Aromatic(ValueAst::Expr(
                Expr::Lit(2),
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
        let remap = IdxRemapping::new(
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
}
