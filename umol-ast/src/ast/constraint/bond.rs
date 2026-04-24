//! Per-covalent-bond constraints.

use std::mem;
use std::slice::Iter;

use strum::EnumDiscriminants;

use super::super::remap::IdxRemapping;
use super::super::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(BondConstraintKind), derive(Hash))]
pub enum BondConstraint {
    Aromatic,
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl BondConstraint {
    pub fn kind(&self) -> BondConstraintKind {
        self.into()
    }

    /// Every `BondConstraint` variant is single-valued per bond.
    pub fn is_unique(&self) -> bool {
        true
    }

    /// `Aromatic` is a flag with no value. `RingCount` / `RingSize` are
    /// undetermined iff their inner value is undetermined.
    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Aromatic => false,
            Self::RingCount(v) | Self::RingSize(v) => v.is_undetermined(),
        }
    }
}

/// Per-bond constraint container. Enforces the per-variant cardinality policy
/// in [`BondConstraint::is_unique`] on insert: unique-kind variants replace any
/// existing entry of the same discriminant (last-wins); multi-kind variants
/// append.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct BondConstraints(Vec<BondConstraint>);

impl BondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[BondConstraint] {
        &self.0
    }

    pub fn iter(&self) -> Iter<'_, BondConstraint> {
        self.0.iter()
    }

    /// Insert a constraint per the per-variant cardinality policy. Returns the
    /// replaced entry if `c.is_unique()` and a same-discriminant entry already
    /// existed; `None` otherwise.
    pub fn add(&mut self, c: BondConstraint) -> Option<BondConstraint> {
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

    pub fn retain(&mut self, mut f: impl FnMut(&BondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// No-op: no `BondConstraint` variant carries an entity index.
    pub fn remap(self, _remap: &IdxRemapping) -> Self {
        self
    }
}

impl FromIterator<BondConstraint> for BondConstraints {
    fn from_iter<I: IntoIterator<Item = BondConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Remapping;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic, BondConstraintKind::Aromatic)]
    #[case::ring_count(BondConstraint::RingCount(ValueAst::Lit(1)), BondConstraintKind::RingCount)]
    #[case::ring_size(BondConstraint::RingSize(ValueAst::Lit(6)), BondConstraintKind::RingSize)]
    fn test_bond_constraint_kind(#[case] c: BondConstraint, #[case] expected: BondConstraintKind) {
        assert_eq!(c.kind(), expected);
    }

    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic)]
    #[case::ring_count(BondConstraint::RingCount(ValueAst::Lit(1)))]
    #[case::ring_size(BondConstraint::RingSize(ValueAst::Lit(6)))]
    fn test_bond_constraint_is_unique_true_everywhere(#[case] c: BondConstraint) {
        assert!(c.is_unique());
    }

    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic, false)]
    #[case::ring_count_lit(BondConstraint::RingCount(ValueAst::Lit(1)), false)]
    #[case::ring_count_undetermined(BondConstraint::RingCount(ValueAst::Undetermined), true)]
    #[case::ring_size_lit(BondConstraint::RingSize(ValueAst::Lit(6)), false)]
    #[case::ring_size_undetermined(BondConstraint::RingSize(ValueAst::Undetermined), true)]
    fn test_bond_constraint_is_undetermined(#[case] c: BondConstraint, #[case] expected: bool) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rstest]
    fn test_bond_constraints_new_is_empty() {
        let cs = BondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[BondConstraint]);
    }

    #[rstest]
    fn test_bond_constraints_add_fresh_returns_none() {
        let mut cs = BondConstraints::new();
        let prev = cs.add(BondConstraint::Aromatic);
        assert_eq!(prev, None);
        assert_eq!(cs.as_slice(), &[BondConstraint::Aromatic]);
    }

    #[rstest]
    fn test_bond_constraints_add_same_kind_replaces() {
        let mut cs = BondConstraints::new();
        cs.add(BondConstraint::RingCount(ValueAst::Lit(1)));
        let prev = cs.add(BondConstraint::RingCount(ValueAst::Lit(2)));
        assert_eq!(prev, Some(BondConstraint::RingCount(ValueAst::Lit(1))));
        assert_eq!(
            cs.as_slice(),
            &[BondConstraint::RingCount(ValueAst::Lit(2))]
        );
        assert_eq!(cs.len(), 1);
    }

    #[rstest]
    fn test_bond_constraints_add_distinct_kinds_coexist() {
        let mut cs = BondConstraints::new();
        assert_eq!(cs.add(BondConstraint::Aromatic), None);
        assert_eq!(cs.add(BondConstraint::RingCount(ValueAst::Lit(1))), None);
        assert_eq!(cs.add(BondConstraint::RingSize(ValueAst::Lit(6))), None);
        assert_eq!(
            cs.as_slice(),
            &[
                BondConstraint::Aromatic,
                BondConstraint::RingCount(ValueAst::Lit(1)),
                BondConstraint::RingSize(ValueAst::Lit(6)),
            ]
        );
    }

    #[rstest]
    fn test_bond_constraints_retain_keeps_matching() {
        let mut cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::RingCount(ValueAst::Lit(1)),
            BondConstraint::RingSize(ValueAst::Lit(6)),
        ]);
        cs.retain(|c| matches!(c, BondConstraint::Aromatic | BondConstraint::RingSize(_)));
        assert_eq!(
            cs.as_slice(),
            &[
                BondConstraint::Aromatic,
                BondConstraint::RingSize(ValueAst::Lit(6))
            ]
        );
    }

    #[rstest]
    fn test_bond_constraints_clear() {
        let mut cs = BondConstraints::from_iter([BondConstraint::Aromatic]);
        cs.clear();
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_bond_constraints_iter_preserves_insertion_order() {
        let cs = BondConstraints::from_iter([
            BondConstraint::RingSize(ValueAst::Lit(6)),
            BondConstraint::Aromatic,
            BondConstraint::RingCount(ValueAst::Lit(1)),
        ]);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                BondConstraint::RingSize(ValueAst::Lit(6)),
                BondConstraint::Aromatic,
                BondConstraint::RingCount(ValueAst::Lit(1)),
            ]
        );
    }

    #[rstest]
    fn test_bond_constraints_from_iter_last_wins_same_kind() {
        let cs = BondConstraints::from_iter([
            BondConstraint::RingCount(ValueAst::Lit(1)),
            BondConstraint::RingCount(ValueAst::Lit(2)),
        ]);
        assert_eq!(
            cs.as_slice(),
            &[BondConstraint::RingCount(ValueAst::Lit(2))]
        );
    }

    #[rstest]
    fn test_bond_constraints_remap_is_noop() {
        let cs = BondConstraints::from_iter([
            BondConstraint::Aromatic,
            BondConstraint::RingSize(ValueAst::Lit(6)),
        ]);
        let remap = IdxRemapping::new(
            Remapping {
                removed_nodes: vec![0, 1, 2],
                removed_edges: vec![0, 1],
            },
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        let after = cs.clone().remap(&remap);
        assert_eq!(after, cs);
    }
}
