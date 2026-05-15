//! Multicenter bond constraints.

use std::mem::{self, replace};
use std::slice::Iter;
use std::vec::IntoIter;

use strum::EnumDiscriminants;

use super::super::remap::IdRemapping;
use super::super::value::ValueAst;

/// Multicenter-bond-scope constraint. Held inline on `MulticenterBondAst` via
/// `MulticenterBondConstraints`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(MulticenterBondConstraintKind), derive(Hash))]
pub enum MulticenterBondConstraint {
    /// Asserted total electron count for the multicenter bond. Cross-checked
    /// by the `ConsistencyValidator` against `sum(MulticenterBondAst::electrons)`.
    ElectronCount(ValueAst),
}

impl MulticenterBondConstraint {
    pub fn electron_count(v: impl Into<ValueAst>) -> Self {
        Self::ElectronCount(v.into())
    }

    pub fn kind(&self) -> MulticenterBondConstraintKind {
        self.into()
    }

    /// Every `MulticenterBondConstraint` variant is single-valued per bond.
    pub fn is_unique(&self) -> bool {
        true
    }

    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::ElectronCount(v) => v.is_undetermined(),
        }
    }

    pub fn simplify(self) -> Self {
        match self {
            Self::ElectronCount(v) => Self::ElectronCount(v.simplify()),
        }
    }

    pub fn remap(self, _remap: &IdRemapping) -> Option<Self> {
        // Value-only: no indices to remap.
        Some(self)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MulticenterBondConstraints(Vec<MulticenterBondConstraint>);

impl MulticenterBondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[MulticenterBondConstraint] {
        &self.0
    }

    pub fn contains(&self, kind: MulticenterBondConstraintKind) -> bool {
        self.0.iter().any(|c| c.kind() == kind)
    }

    pub fn get(&self, kind: MulticenterBondConstraintKind) -> Option<&MulticenterBondConstraint> {
        self.0.iter().find(|c| c.kind() == kind)
    }

    pub fn get_mut(
        &mut self,
        kind: MulticenterBondConstraintKind,
    ) -> Option<&mut MulticenterBondConstraint> {
        self.0.iter_mut().find(|c| c.kind() == kind)
    }

    pub fn electron_count(&self) -> ValueAst {
        match self.get(MulticenterBondConstraintKind::ElectronCount) {
            Some(MulticenterBondConstraint::ElectronCount(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn iter(&self) -> Iter<'_, MulticenterBondConstraint> {
        self.0.iter()
    }

    /// Insert a constraint per the per-variant cardinality policy. Returns
    /// the replaced entry if `c.is_unique()` and a same-kind entry already
    /// existed; `None` otherwise.
    pub fn add(&mut self, c: MulticenterBondConstraint) -> Option<MulticenterBondConstraint> {
        if c.is_unique() {
            if let Some(i) = self.0.iter().position(|x| x.kind() == c.kind()) {
                return Some(replace(&mut self.0[i], c));
            }
        }
        self.0.push(c);
        None
    }

    pub fn retain(&mut self, mut f: impl FnMut(&MulticenterBondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    pub fn take(&mut self) -> impl Iterator<Item = MulticenterBondConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn simplify_each(&mut self) {
        for c in self.0.iter_mut() {
            *c = mem::replace(
                c,
                MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined),
            )
            .simplify();
        }
    }

    pub fn remove(
        &mut self,
        kind: MulticenterBondConstraintKind,
    ) -> Option<MulticenterBondConstraint> {
        let pos = self.0.iter().position(|c| c.kind() == kind)?;
        Some(self.0.remove(pos))
    }

    /// Iterate over every entry of `kind`. Currently every variant is
    /// single-valued so this yields at most one entry.
    pub fn get_all(
        &self,
        kind: MulticenterBondConstraintKind,
    ) -> impl Iterator<Item = &MulticenterBondConstraint> {
        self.0.iter().filter(move |c| c.kind() == kind)
    }

    /// Remove every entry of `kind`, returning them in insertion order.
    pub fn remove_all(
        &mut self,
        kind: MulticenterBondConstraintKind,
    ) -> Vec<MulticenterBondConstraint> {
        let mut out = Vec::new();
        self.0.retain(|c| {
            if c.kind() == kind {
                out.push(c.clone());
                false
            } else {
                true
            }
        });
        out
    }

    pub fn remap(self, remap: &IdRemapping) -> Self {
        Self(self.0.into_iter().filter_map(|c| c.remap(remap)).collect())
    }
}

impl FromIterator<MulticenterBondConstraint> for MulticenterBondConstraints {
    fn from_iter<I: IntoIterator<Item = MulticenterBondConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

impl IntoIterator for MulticenterBondConstraints {
    type Item = MulticenterBondConstraint;
    type IntoIter = IntoIter<MulticenterBondConstraint>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<MulticenterBondConstraint> for MulticenterBondConstraints {
    fn from(c: MulticenterBondConstraint) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<MulticenterBondConstraint>> for MulticenterBondConstraints {
    fn from(cs: Vec<MulticenterBondConstraint>) -> Self {
        Self::from_iter(cs)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Remapping;

    use super::*;
    use crate::ast::value::Expr;

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraint::electron_count(2),
        MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2))
    )]
    fn test_multicenter_bond_constraint_constructors(
        #[case] actual: MulticenterBondConstraint,
        #[case] expected: MulticenterBondConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::electron_count(
        MulticenterBondConstraint::electron_count(2),
        MulticenterBondConstraintKind::ElectronCount
    )]
    fn test_multicenter_bond_constraint_kind(
        #[case] c: MulticenterBondConstraint,
        #[case] expected: MulticenterBondConstraintKind,
    ) {
        assert_eq!(c.kind(), expected);
    }

    #[rstest]
    #[case::electron_count(MulticenterBondConstraint::electron_count(2))]
    fn test_multicenter_bond_constraint_is_unique(#[case] c: MulticenterBondConstraint) {
        assert!(c.is_unique());
    }

    #[rstest]
    #[case::lit(MulticenterBondConstraint::electron_count(2), false)]
    #[case::undetermined(MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined), true)]
    fn test_multicenter_bond_constraint_is_undetermined(
        #[case] c: MulticenterBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rstest]
    #[case::folds_expr(
        MulticenterBondConstraint::ElectronCount(ValueAst::Expr(Box::new(Expr::Lit(2)))),
        MulticenterBondConstraint::electron_count(2)
    )]
    fn test_multicenter_bond_constraint_simplify(
        #[case] input: MulticenterBondConstraint,
        #[case] expected: MulticenterBondConstraint,
    ) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::lit(MulticenterBondConstraint::electron_count(2))]
    #[case::undetermined(MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined))]
    fn test_multicenter_bond_constraint_simplify_identity(
        #[case] input: MulticenterBondConstraint,
    ) {
        assert_eq!(input.clone().simplify(), input);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_new() {
        let cs = MulticenterBondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[MulticenterBondConstraint]);
    }

    #[rstest]
    #[case::present(MulticenterBondConstraintKind::ElectronCount, true)]
    fn test_multicenter_bond_constraints_contains(
        #[case] kind: MulticenterBondConstraintKind,
        #[case] expected: bool,
    ) {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        assert_eq!(cs.contains(kind), expected);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_contains_absent() {
        let cs = MulticenterBondConstraints::new();
        assert!(!cs.contains(MulticenterBondConstraintKind::ElectronCount));
    }

    #[rstest]
    fn test_multicenter_bond_constraints_get() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        assert_eq!(
            cs.get(MulticenterBondConstraintKind::ElectronCount),
            Some(&MulticenterBondConstraint::electron_count(2)),
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_get_absent() {
        let cs = MulticenterBondConstraints::new();
        assert_eq!(cs.get(MulticenterBondConstraintKind::ElectronCount), None);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_get_mut() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        let entry = cs
            .get_mut(MulticenterBondConstraintKind::ElectronCount)
            .unwrap();
        *entry = MulticenterBondConstraint::electron_count(4);
        assert_eq!(
            cs.as_slice(),
            &[MulticenterBondConstraint::electron_count(4)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_get_mut_absent() {
        let mut cs = MulticenterBondConstraints::new();
        assert!(cs
            .get_mut(MulticenterBondConstraintKind::ElectronCount)
            .is_none());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_iter() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![MulticenterBondConstraint::electron_count(2)]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(
        vec![MulticenterBondConstraint::electron_count(2)],
        vec![None],
        vec![MulticenterBondConstraint::electron_count(2)],
    )]
    #[case::replace_same_kind(
        vec![
            MulticenterBondConstraint::electron_count(2),
            MulticenterBondConstraint::electron_count(4),
        ],
        vec![None, Some(MulticenterBondConstraint::electron_count(2))],
        vec![MulticenterBondConstraint::electron_count(4)],
    )]
    fn test_multicenter_bond_constraints_add(
        #[case] sequence: Vec<MulticenterBondConstraint>,
        #[case] expected_returns: Vec<Option<MulticenterBondConstraint>>,
        #[case] expected_state: Vec<MulticenterBondConstraint>,
    ) {
        let mut cs = MulticenterBondConstraints::new();
        let returns: Vec<_> = sequence.into_iter().map(|c| cs.add(c)).collect();
        assert_eq!(returns, expected_returns);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(
        |c: &MulticenterBondConstraint| matches!(c, MulticenterBondConstraint::ElectronCount(_)),
        vec![MulticenterBondConstraint::electron_count(2)],
    )]
    #[case::all_dropped(|_: &MulticenterBondConstraint| false, vec![])]
    fn test_multicenter_bond_constraints_retain(
        #[case] predicate: impl FnMut(&MulticenterBondConstraint) -> bool,
        #[case] expected: Vec<MulticenterBondConstraint>,
    ) {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        cs.retain(predicate);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_clear() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        cs.clear();
        assert_eq!(cs, MulticenterBondConstraints::new());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_take() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(drained, vec![MulticenterBondConstraint::electron_count(2)]);
        assert_eq!(cs, MulticenterBondConstraints::new());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_simplify_each() {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::ElectronCount(
            ValueAst::Expr(Box::new(Expr::Lit(2))),
        ));
        cs.simplify_each();
        assert_eq!(
            cs,
            MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2)),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::present(
        MulticenterBondConstraintKind::ElectronCount,
        Some(MulticenterBondConstraint::electron_count(2)),
        vec![],
    )]
    fn test_multicenter_bond_constraints_remove(
        #[case] kind: MulticenterBondConstraintKind,
        #[case] expected_returned: Option<MulticenterBondConstraint>,
        #[case] expected_state: Vec<MulticenterBondConstraint>,
    ) {
        let mut cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        assert_eq!(cs.remove(kind), expected_returned);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_remove_absent() {
        let mut cs = MulticenterBondConstraints::new();
        assert_eq!(
            cs.remove(MulticenterBondConstraintKind::ElectronCount),
            None,
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_remap() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        let remap = IdRemapping::new(
            Remapping {
                removed_nodes: vec![0, 1],
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
    #[case::single(
        vec![MulticenterBondConstraint::electron_count(2)],
        vec![MulticenterBondConstraint::electron_count(2)],
    )]
    #[case::same_kind_last_wins(
        vec![
            MulticenterBondConstraint::electron_count(2),
            MulticenterBondConstraint::electron_count(4),
        ],
        vec![MulticenterBondConstraint::electron_count(4)],
    )]
    #[case::empty(vec![], vec![])]
    fn test_multicenter_bond_constraints_from_iter(
        #[case] input: Vec<MulticenterBondConstraint>,
        #[case] expected: Vec<MulticenterBondConstraint>,
    ) {
        let cs = MulticenterBondConstraints::from_iter(input);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_into_iter() {
        let cs = MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2));
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![MulticenterBondConstraint::electron_count(2)]
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_from_multicenter_bond_constraint() {
        let cs: MulticenterBondConstraints = MulticenterBondConstraint::electron_count(2).into();
        assert_eq!(
            cs.as_slice(),
            &[MulticenterBondConstraint::electron_count(2)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_from_vec() {
        let cs: MulticenterBondConstraints =
            vec![MulticenterBondConstraint::electron_count(2)].into();
        assert_eq!(
            cs.as_slice(),
            &[MulticenterBondConstraint::electron_count(2)],
        );
    }
}
