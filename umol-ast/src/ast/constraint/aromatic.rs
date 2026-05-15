//! Per-aromatic-system constraints.

use std::mem::{self, replace};
use std::slice::Iter;
use std::vec::IntoIter;

use strum::EnumDiscriminants;

use super::super::remap::IdRemapping;
use super::super::value::ValueAst;

/// Aromatic-system-scope constraint. Held inline on `AromaticSystemAst` via
/// `AromaticSystemConstraints`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(AromaticSystemConstraintKind), derive(Hash))]
pub enum AromaticSystemConstraint {
    /// Asserted total π-electron count for the system. Cross-checked by the
    /// `ConsistencyValidator` against `sum(AromaticSystemAst::electrons)`.
    ElectronCount(ValueAst),
}

impl AromaticSystemConstraint {
    pub fn electron_count(v: impl Into<ValueAst>) -> Self {
        Self::ElectronCount(v.into())
    }

    pub fn kind(&self) -> AromaticSystemConstraintKind {
        self.into()
    }

    /// Every `AromaticSystemConstraint` variant is single-valued per system.
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
pub struct AromaticSystemConstraints(Vec<AromaticSystemConstraint>);

impl AromaticSystemConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[AromaticSystemConstraint] {
        &self.0
    }

    pub fn contains(&self, kind: AromaticSystemConstraintKind) -> bool {
        self.0.iter().any(|c| c.kind() == kind)
    }

    pub fn get(&self, kind: AromaticSystemConstraintKind) -> Option<&AromaticSystemConstraint> {
        self.0.iter().find(|c| c.kind() == kind)
    }

    pub fn get_mut(
        &mut self,
        kind: AromaticSystemConstraintKind,
    ) -> Option<&mut AromaticSystemConstraint> {
        self.0.iter_mut().find(|c| c.kind() == kind)
    }

    pub fn electron_count(&self) -> ValueAst {
        match self.get(AromaticSystemConstraintKind::ElectronCount) {
            Some(AromaticSystemConstraint::ElectronCount(v)) => v.clone(),
            _ => ValueAst::Undetermined,
        }
    }

    pub fn iter(&self) -> Iter<'_, AromaticSystemConstraint> {
        self.0.iter()
    }

    /// Insert a constraint per the per-variant cardinality policy. Returns
    /// the replaced entry if `c.is_unique()` and a same-kind entry already
    /// existed; `None` otherwise.
    pub fn add(&mut self, c: AromaticSystemConstraint) -> Option<AromaticSystemConstraint> {
        if c.is_unique() {
            if let Some(i) = self.0.iter().position(|x| x.kind() == c.kind()) {
                return Some(replace(&mut self.0[i], c));
            }
        }
        self.0.push(c);
        None
    }

    pub fn retain(&mut self, mut f: impl FnMut(&AromaticSystemConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    pub fn take(&mut self) -> impl Iterator<Item = AromaticSystemConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn simplify_each(&mut self) {
        for c in self.0.iter_mut() {
            *c = mem::replace(
                c,
                AromaticSystemConstraint::ElectronCount(ValueAst::Undetermined),
            )
            .simplify();
        }
    }

    pub fn remove(
        &mut self,
        kind: AromaticSystemConstraintKind,
    ) -> Option<AromaticSystemConstraint> {
        let pos = self.0.iter().position(|c| c.kind() == kind)?;
        Some(self.0.remove(pos))
    }

    /// Iterate over every entry of `kind`. Currently every variant is
    /// single-valued so this yields at most one entry.
    pub fn get_all(
        &self,
        kind: AromaticSystemConstraintKind,
    ) -> impl Iterator<Item = &AromaticSystemConstraint> {
        self.0.iter().filter(move |c| c.kind() == kind)
    }

    /// Remove every entry of `kind`, returning them in insertion order.
    pub fn remove_all(
        &mut self,
        kind: AromaticSystemConstraintKind,
    ) -> Vec<AromaticSystemConstraint> {
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

impl FromIterator<AromaticSystemConstraint> for AromaticSystemConstraints {
    fn from_iter<I: IntoIterator<Item = AromaticSystemConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

impl IntoIterator for AromaticSystemConstraints {
    type Item = AromaticSystemConstraint;
    type IntoIter = IntoIter<AromaticSystemConstraint>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<AromaticSystemConstraint> for AromaticSystemConstraints {
    fn from(c: AromaticSystemConstraint) -> Self {
        Self::from_iter([c])
    }
}

impl From<Vec<AromaticSystemConstraint>> for AromaticSystemConstraints {
    fn from(cs: Vec<AromaticSystemConstraint>) -> Self {
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
        AromaticSystemConstraint::electron_count(6),
        AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))
    )]
    fn test_aromatic_system_constraint_constructors(
        #[case] actual: AromaticSystemConstraint,
        #[case] expected: AromaticSystemConstraint,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::electron_count(
        AromaticSystemConstraint::electron_count(6),
        AromaticSystemConstraintKind::ElectronCount
    )]
    fn test_aromatic_system_constraint_kind(
        #[case] c: AromaticSystemConstraint,
        #[case] expected: AromaticSystemConstraintKind,
    ) {
        assert_eq!(c.kind(), expected);
    }

    #[rstest]
    #[case::electron_count(AromaticSystemConstraint::electron_count(6))]
    fn test_aromatic_system_constraint_is_unique(#[case] c: AromaticSystemConstraint) {
        assert!(c.is_unique());
    }

    #[rstest]
    #[case::lit(AromaticSystemConstraint::electron_count(6), false)]
    #[case::undetermined(AromaticSystemConstraint::ElectronCount(ValueAst::Undetermined), true)]
    fn test_aromatic_system_constraint_is_undetermined(
        #[case] c: AromaticSystemConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rstest]
    #[case::folds_expr(
        AromaticSystemConstraint::ElectronCount(ValueAst::Expr(Box::new(Expr::Lit(6)))),
        AromaticSystemConstraint::electron_count(6)
    )]
    fn test_aromatic_system_constraint_simplify(
        #[case] input: AromaticSystemConstraint,
        #[case] expected: AromaticSystemConstraint,
    ) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::lit(AromaticSystemConstraint::electron_count(6))]
    #[case::undetermined(AromaticSystemConstraint::ElectronCount(ValueAst::Undetermined))]
    fn test_aromatic_system_constraint_simplify_identity(#[case] input: AromaticSystemConstraint) {
        assert_eq!(input.clone().simplify(), input);
    }

    #[rstest]
    fn test_aromatic_system_constraints_new() {
        let cs = AromaticSystemConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[AromaticSystemConstraint]);
    }

    #[rstest]
    #[case::present(AromaticSystemConstraintKind::ElectronCount, true)]
    fn test_aromatic_system_constraints_contains(
        #[case] kind: AromaticSystemConstraintKind,
        #[case] expected: bool,
    ) {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        assert_eq!(cs.contains(kind), expected);
    }

    #[rstest]
    fn test_aromatic_system_constraints_contains_absent() {
        let cs = AromaticSystemConstraints::new();
        assert!(!cs.contains(AromaticSystemConstraintKind::ElectronCount));
    }

    #[rstest]
    fn test_aromatic_system_constraints_get() {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        assert_eq!(
            cs.get(AromaticSystemConstraintKind::ElectronCount),
            Some(&AromaticSystemConstraint::electron_count(6)),
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_get_absent() {
        let cs = AromaticSystemConstraints::new();
        assert_eq!(cs.get(AromaticSystemConstraintKind::ElectronCount), None,);
    }

    #[rstest]
    fn test_aromatic_system_constraints_get_mut() {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        let entry = cs
            .get_mut(AromaticSystemConstraintKind::ElectronCount)
            .unwrap();
        *entry = AromaticSystemConstraint::electron_count(10);
        assert_eq!(
            cs.as_slice(),
            &[AromaticSystemConstraint::electron_count(10)],
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_get_mut_absent() {
        let mut cs = AromaticSystemConstraints::new();
        assert!(cs
            .get_mut(AromaticSystemConstraintKind::ElectronCount)
            .is_none());
    }

    #[rstest]
    fn test_aromatic_system_constraints_iter() {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(collected, vec![AromaticSystemConstraint::electron_count(6)]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fresh(
        vec![AromaticSystemConstraint::electron_count(6)],
        vec![None],
        vec![AromaticSystemConstraint::electron_count(6)],
    )]
    #[case::replace_same_kind(
        vec![
            AromaticSystemConstraint::electron_count(6),
            AromaticSystemConstraint::electron_count(10),
        ],
        vec![None, Some(AromaticSystemConstraint::electron_count(6))],
        vec![AromaticSystemConstraint::electron_count(10)],
    )]
    fn test_aromatic_system_constraints_add(
        #[case] sequence: Vec<AromaticSystemConstraint>,
        #[case] expected_returns: Vec<Option<AromaticSystemConstraint>>,
        #[case] expected_state: Vec<AromaticSystemConstraint>,
    ) {
        let mut cs = AromaticSystemConstraints::new();
        let returns: Vec<_> = sequence.into_iter().map(|c| cs.add(c)).collect();
        assert_eq!(returns, expected_returns);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::partial(
        |c: &AromaticSystemConstraint| matches!(c, AromaticSystemConstraint::ElectronCount(_)),
        vec![AromaticSystemConstraint::electron_count(6)],
    )]
    #[case::all_dropped(|_: &AromaticSystemConstraint| false, vec![])]
    fn test_aromatic_system_constraints_retain(
        #[case] predicate: impl FnMut(&AromaticSystemConstraint) -> bool,
        #[case] expected: Vec<AromaticSystemConstraint>,
    ) {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        cs.retain(predicate);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_aromatic_system_constraints_clear() {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        cs.clear();
        assert_eq!(cs, AromaticSystemConstraints::new());
    }

    #[rstest]
    fn test_aromatic_system_constraints_take() {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        let drained: Vec<_> = cs.take().collect();
        assert_eq!(drained, vec![AromaticSystemConstraint::electron_count(6)]);
        assert_eq!(cs, AromaticSystemConstraints::new());
    }

    #[rstest]
    fn test_aromatic_system_constraints_simplify_each() {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::ElectronCount(
            ValueAst::Expr(Box::new(Expr::Lit(6))),
        ));
        cs.simplify_each();
        assert_eq!(
            cs,
            AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6)),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::present(
        AromaticSystemConstraintKind::ElectronCount,
        Some(AromaticSystemConstraint::electron_count(6)),
        vec![],
    )]
    fn test_aromatic_system_constraints_remove(
        #[case] kind: AromaticSystemConstraintKind,
        #[case] expected_returned: Option<AromaticSystemConstraint>,
        #[case] expected_state: Vec<AromaticSystemConstraint>,
    ) {
        let mut cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        assert_eq!(cs.remove(kind), expected_returned);
        assert_eq!(cs.as_slice(), expected_state.as_slice());
    }

    #[rstest]
    fn test_aromatic_system_constraints_remove_absent() {
        let mut cs = AromaticSystemConstraints::new();
        assert_eq!(cs.remove(AromaticSystemConstraintKind::ElectronCount), None,);
    }

    #[rstest]
    fn test_aromatic_system_constraints_remap() {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
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
        vec![AromaticSystemConstraint::electron_count(6)],
        vec![AromaticSystemConstraint::electron_count(6)],
    )]
    #[case::same_kind_last_wins(
        vec![
            AromaticSystemConstraint::electron_count(2),
            AromaticSystemConstraint::electron_count(6),
        ],
        vec![AromaticSystemConstraint::electron_count(6)],
    )]
    #[case::empty(vec![], vec![])]
    fn test_aromatic_system_constraints_from_iter(
        #[case] input: Vec<AromaticSystemConstraint>,
        #[case] expected: Vec<AromaticSystemConstraint>,
    ) {
        let cs = AromaticSystemConstraints::from_iter(input);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_aromatic_system_constraints_into_iter() {
        let cs = AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6));
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(collected, vec![AromaticSystemConstraint::electron_count(6)]);
    }

    #[rstest]
    fn test_aromatic_system_constraints_from_aromatic_system_constraint() {
        let cs: AromaticSystemConstraints = AromaticSystemConstraint::electron_count(6).into();
        assert_eq!(
            cs.as_slice(),
            &[AromaticSystemConstraint::electron_count(6)]
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_from_vec() {
        let cs: AromaticSystemConstraints =
            vec![AromaticSystemConstraint::electron_count(6)].into();
        assert_eq!(
            cs.as_slice(),
            &[AromaticSystemConstraint::electron_count(6)]
        );
    }
}
