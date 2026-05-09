//! Dative bond constraints.

use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use strum::EnumDiscriminants;

use super::super::remap::IdxRemapping;
use super::super::value::ValueAst;

/// Dative-bond-scope constraint. Held inline on `DativeBondAst` via
/// `DativeBondConstraints`. `Aromatic` flags the dative bond as part of an
/// aromatic system (e.g. the N→B π-donation in borazine, O→B in boroxine,
/// or a C→M coordination spanning a metallaaromatic ring).
#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(DativeBondConstraintKind), derive(Hash))]
pub enum DativeBondConstraint {
    Aromatic,
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl DativeBondConstraint {
    pub fn ring_count(v: impl Into<ValueAst>) -> Self {
        Self::RingCount(v.into())
    }

    pub fn ring_size(v: impl Into<ValueAst>) -> Self {
        Self::RingSize(v.into())
    }

    pub fn kind(&self) -> DativeBondConstraintKind {
        self.into()
    }

    /// Every `DativeBondConstraint` variant is single-valued per dative bond.
    pub fn is_unique(&self) -> bool {
        true
    }

    /// `Aromatic` is a flag with no value. `RingCount`/`RingSize` are
    /// undetermined iff their inner value is.
    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Aromatic => false,
            Self::RingCount(v) | Self::RingSize(v) => v.is_undetermined(),
        }
    }

    /// Simplify the inner `ValueAst` of `RingCount` / `RingSize`. `Aromatic`
    /// is unchanged.
    pub fn simplify(self) -> Self {
        match self {
            Self::Aromatic => Self::Aromatic,
            Self::RingCount(v) => Self::RingCount(v.simplify()),
            Self::RingSize(v) => Self::RingSize(v.simplify()),
        }
    }

    pub fn remap(self, _remap: &IdxRemapping) -> Option<Self> {
        // Value-only: no indices to remap.
        Some(self)
    }
}

/// Per-dative-bond constraint container. Enforces the per-variant cardinality
/// policy in [`DativeBondConstraint::is_unique`] on insert.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct DativeBondConstraints(Vec<DativeBondConstraint>);

impl DativeBondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[DativeBondConstraint] {
        &self.0
    }

    pub fn contains(&self, kind: DativeBondConstraintKind) -> bool {
        self.0.iter().any(|c| c.kind() == kind)
    }

    pub fn get(&self, kind: DativeBondConstraintKind) -> Option<&DativeBondConstraint> {
        self.0.iter().find(|c| c.kind() == kind)
    }

    pub fn get_mut(&mut self, kind: DativeBondConstraintKind) -> Option<&mut DativeBondConstraint> {
        self.0.iter_mut().find(|c| c.kind() == kind)
    }

    pub fn iter(&self) -> Iter<'_, DativeBondConstraint> {
        self.0.iter()
    }

    pub fn add(&mut self, c: DativeBondConstraint) -> Option<DativeBondConstraint> {
        if c.is_unique() {
            if let Some(pos) = self
                .0
                .iter()
                .position(|e| mem::discriminant(e) == mem::discriminant(&c))
            {
                return Some(mem::replace(&mut self.0[pos], c));
            }
        }
        self.0.push(c);
        None
    }

    pub fn retain(&mut self, mut f: impl FnMut(&DativeBondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = DativeBondConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    /// Simplify each contained constraint's inner value in place.
    pub fn simplify_each(&mut self) {
        for c in self.0.iter_mut() {
            *c = mem::replace(c, DativeBondConstraint::Aromatic).simplify();
        }
    }

    pub fn remove(&mut self, kind: DativeBondConstraintKind) -> Option<DativeBondConstraint> {
        let pos = self.0.iter().position(|c| c.kind() == kind)?;
        Some(self.0.remove(pos))
    }

    pub fn remap(self, remap: &IdxRemapping) -> Self {
        Self(self.0.into_iter().filter_map(|c| c.remap(remap)).collect())
    }
}

impl FromIterator<DativeBondConstraint> for DativeBondConstraints {
    fn from_iter<I: IntoIterator<Item = DativeBondConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

impl IntoIterator for DativeBondConstraints {
    type Item = DativeBondConstraint;
    type IntoIter = IntoIter<DativeBondConstraint>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<DativeBondConstraint> for DativeBondConstraints {
    fn from(c: DativeBondConstraint) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<DativeBondConstraint>> for DativeBondConstraints {
    fn from(cs: Vec<DativeBondConstraint>) -> Self {
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
    #[case::ring_count(DativeBondConstraint::ring_count(1), DativeBondConstraint::RingCount(ValueAst::Lit(1)))]
    #[case::ring_size(DativeBondConstraint::ring_size(6), DativeBondConstraint::RingSize(ValueAst::Lit(6)))]
    fn test_dative_bond_constraint_constructors(
        #[case] actual: DativeBondConstraint,
        #[case] expected: DativeBondConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic, DativeBondConstraintKind::Aromatic)]
    #[case::ring_count(DativeBondConstraint::ring_count(1), DativeBondConstraintKind::RingCount)]
    #[case::ring_size(DativeBondConstraint::ring_size(6), DativeBondConstraintKind::RingSize)]
    fn test_dative_bond_constraint_kind(
        #[case] c: DativeBondConstraint,
        #[case] expected: DativeBondConstraintKind,
    ) {
        assert_eq!(c.kind(), expected);
    }

    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic)]
    #[case::ring_count(DativeBondConstraint::ring_count(1))]
    #[case::ring_size(DativeBondConstraint::ring_size(6))]
    fn test_dative_bond_constraint_is_unique(#[case] c: DativeBondConstraint) {
        assert!(c.is_unique());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic, false)]
    #[case::ring_count_lit(DativeBondConstraint::ring_count(1), false)]
    #[case::ring_count_undetermined(DativeBondConstraint::RingCount(ValueAst::Undetermined), true)]
    #[case::ring_size_lit(DativeBondConstraint::ring_size(6), false)]
    #[case::ring_size_undetermined(DativeBondConstraint::RingSize(ValueAst::Undetermined), true)]
    fn test_dative_bond_constraint_is_undetermined(
        #[case] c: DativeBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_count_folds_expr(
        DativeBondConstraint::RingCount(ValueAst::Expr(Expr::Lit(2))),
        DativeBondConstraint::ring_count(2),
    )]
    #[case::ring_size_folds_expr(
        DativeBondConstraint::RingSize(ValueAst::Expr(Expr::Lit(6))),
        DativeBondConstraint::ring_size(6),
    )]
    fn test_dative_bond_constraint_simplify(
        #[case] input: DativeBondConstraint,
        #[case] expected: DativeBondConstraint,
    ) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::aromatic(DativeBondConstraint::Aromatic)]
    #[case::ring_count_lit(DativeBondConstraint::ring_count(1))]
    #[case::ring_size_undetermined(DativeBondConstraint::RingSize(ValueAst::Undetermined))]
    fn test_dative_bond_constraint_simplify_identity(#[case] input: DativeBondConstraint) {
        assert_eq!(input.clone().simplify(), input);
    }

    #[rstest]
    fn test_dative_bond_constraints_new() {
        let cs = DativeBondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[DativeBondConstraint]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(DativeBondConstraintKind::Aromatic, true)]
    #[case::ring_size_present(DativeBondConstraintKind::RingSize, true)]
    #[case::ring_count_absent(DativeBondConstraintKind::RingCount, false)]
    fn test_dative_bond_constraints_contains(
        #[case] kind: DativeBondConstraintKind,
        #[case] expected: bool,
    ) {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_size(6),
        ]);
        assert_eq!(cs.contains(kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintKind::Aromatic, Some(DativeBondConstraint::Aromatic))]
    #[case::ring_size(DativeBondConstraintKind::RingSize, Some(DativeBondConstraint::ring_size(6)))]
    #[case::ring_count_absent(DativeBondConstraintKind::RingCount, None)]
    fn test_dative_bond_constraints_get(
        #[case] kind: DativeBondConstraintKind,
        #[case] expected: Option<DativeBondConstraint>,
    ) {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_size(6),
        ]);
        assert_eq!(cs.get(kind), expected.as_ref());
    }

    #[rstest]
    fn test_dative_bond_constraints_get_mut() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_size(6),
        ]);
        let entry = cs.get_mut(DativeBondConstraintKind::RingSize).unwrap();
        *entry = DativeBondConstraint::ring_size(5);
        assert_eq!(
            cs.as_slice(),
            &[
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_size(5)
            ],
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_get_mut_absent() {
        let mut cs = DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic]);
        assert!(cs.get_mut(DativeBondConstraintKind::RingCount).is_none());
    }

    #[rstest]
    fn test_dative_bond_constraints_iter() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::ring_size(6),
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_count(1),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                DativeBondConstraint::ring_size(6),
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_count(1),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(
        vec![DativeBondConstraint::Aromatic],
        vec![None],
        vec![DativeBondConstraint::Aromatic],
    )]
    #[case::replace_same_kind(
        vec![DativeBondConstraint::ring_count(1), DativeBondConstraint::ring_count(2)],
        vec![None, Some(DativeBondConstraint::ring_count(1))],
        vec![DativeBondConstraint::ring_count(2)],
    )]
    #[case::replace_unit_variant(
        vec![DativeBondConstraint::Aromatic, DativeBondConstraint::Aromatic],
        vec![None, Some(DativeBondConstraint::Aromatic)],
        vec![DativeBondConstraint::Aromatic],
    )]
    #[case::distinct_kinds(
        vec![
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_count(1),
            DativeBondConstraint::ring_size(6),
        ],
        vec![None, None, None],
        vec![
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_count(1),
            DativeBondConstraint::ring_size(6),
        ],
    )]
    fn test_dative_bond_constraints_add(
        #[case] sequence: Vec<DativeBondConstraint>,
        #[case] expected_returns: Vec<Option<DativeBondConstraint>>,
        #[case] expected_state: Vec<DativeBondConstraint>,
    ) {
        let mut cs = DativeBondConstraints::new();
        let returns: Vec<_> = sequence.into_iter().map(|c| cs.add(c)).collect();
        assert_eq!(returns, expected_returns);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(
        |c: &DativeBondConstraint| matches!(c, DativeBondConstraint::Aromatic | DativeBondConstraint::RingSize(_)),
        vec![DativeBondConstraint::Aromatic, DativeBondConstraint::ring_size(6)],
    )]
    #[case::all_dropped(|_: &DativeBondConstraint| false, vec![])]
    fn test_dative_bond_constraints_retain(
        #[case] predicate: impl FnMut(&DativeBondConstraint) -> bool,
        #[case] expected: Vec<DativeBondConstraint>,
    ) {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_count(1),
            DativeBondConstraint::ring_size(6),
        ]);
        cs.retain(predicate);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_dative_bond_constraints_clear() {
        let mut cs = DativeBondConstraints::from_iter([DativeBondConstraint::Aromatic]);
        cs.clear();
        assert_eq!(cs, DativeBondConstraints::new());
    }

    #[rstest]
    fn test_dative_bond_constraints_take() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_size(6),
        ]);
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(
            drained,
            vec![
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_size(6)
            ],
        );
        assert_eq!(cs, DativeBondConstraints::new());
    }

    #[rstest]
    fn test_dative_bond_constraints_simplify_each() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::RingCount(ValueAst::Expr(Expr::Lit(1))),
            DativeBondConstraint::RingSize(ValueAst::Expr(Expr::Lit(6))),
        ]);
        cs.simplify_each();
        assert_eq!(
            cs,
            DativeBondConstraints::from_iter([
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_count(1),
                DativeBondConstraint::ring_size(6),
            ]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic_present(
        DativeBondConstraintKind::Aromatic,
        Some(DativeBondConstraint::Aromatic),
        vec![DativeBondConstraint::ring_size(6)],
    )]
    #[case::ring_size_present(
        DativeBondConstraintKind::RingSize,
        Some(DativeBondConstraint::ring_size(6)),
        vec![DativeBondConstraint::Aromatic],
    )]
    #[case::ring_count_absent(
        DativeBondConstraintKind::RingCount,
        None,
        vec![DativeBondConstraint::Aromatic, DativeBondConstraint::ring_size(6)],
    )]
    fn test_dative_bond_constraints_remove(
        #[case] kind: DativeBondConstraintKind,
        #[case] expected_returned: Option<DativeBondConstraint>,
        #[case] expected_state: Vec<DativeBondConstraint>,
    ) {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_size(6),
        ]);
        assert_eq!(cs.remove(kind), expected_returned);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rstest]
    fn test_dative_bond_constraints_remap() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_size(6),
        ]);
        let remap = IdxRemapping::new(
            umol_graph_core::Remapping {
                removed_nodes: vec![1],
                removed_edges: vec![1],
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
        vec![DativeBondConstraint::Aromatic, DativeBondConstraint::ring_count(1)],
        vec![DativeBondConstraint::Aromatic, DativeBondConstraint::ring_count(1)],
    )]
    #[case::same_kind_last_wins(
        vec![DativeBondConstraint::ring_count(1), DativeBondConstraint::ring_count(2)],
        vec![DativeBondConstraint::ring_count(2)],
    )]
    #[case::empty(vec![], vec![])]
    fn test_dative_bond_constraints_from_iter(
        #[case] input: Vec<DativeBondConstraint>,
        #[case] expected: Vec<DativeBondConstraint>,
    ) {
        let cs = DativeBondConstraints::from_iter(input);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_dative_bond_constraints_into_iter() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_size(6),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_size(6)
            ],
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_from_dative_bond_constraint() {
        let cs: DativeBondConstraints = DativeBondConstraint::Aromatic.into();
        assert_eq!(cs.as_slice(), &[DativeBondConstraint::Aromatic]);
    }

    #[rstest]
    fn test_dative_bond_constraints_from_vec() {
        let cs: DativeBondConstraints = vec![
            DativeBondConstraint::Aromatic,
            DativeBondConstraint::ring_size(6),
        ]
        .into();
        assert_eq!(
            cs.as_slice(),
            &[
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_size(6)
            ],
        );
    }
}
