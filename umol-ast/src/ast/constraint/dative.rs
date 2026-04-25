//! Per-dative-bond constraints.

use std::mem;
use std::slice::Iter;

use strum::EnumDiscriminants;

use super::super::remap::IdxRemapping;
use super::super::value::ValueAst;

/// Dative-bond-scope constraint. Held inline on `DativeBondAst` via
/// `DativeBondConstraints`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(DativeBondConstraintKind), derive(Hash))]
pub enum DativeBondConstraint {
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl DativeBondConstraint {
    pub fn kind(&self) -> DativeBondConstraintKind {
        self.into()
    }

    /// Single-valued per dative bond: both variants pin a single integer
    /// shape per kind.
    pub fn is_unique(&self) -> bool {
        true
    }

    /// `RingCount`/`RingSize` are undetermined iff their inner value is.
    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::RingCount(v) | Self::RingSize(v) => v.is_undetermined(),
        }
    }

    /// Simplify the inner `ValueAst`.
    pub fn simplify(self) -> Self {
        match self {
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
            *c =
                mem::replace(c, DativeBondConstraint::RingCount(ValueAst::Undetermined)).simplify();
        }
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_count(DativeBondConstraint::RingCount(ValueAst::Lit(1)), DativeBondConstraintKind::RingCount)]
    #[case::ring_size(DativeBondConstraint::RingSize(ValueAst::Lit(6)), DativeBondConstraintKind::RingSize)]
    fn test_dative_bond_constraint_kind(
        #[case] c: DativeBondConstraint,
        #[case] expected: DativeBondConstraintKind,
    ) {
        assert_eq!(c.kind(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_count(DativeBondConstraint::RingCount(ValueAst::Lit(1)), true)]
    #[case::ring_size(DativeBondConstraint::RingSize(ValueAst::Lit(6)), true)]
    fn test_dative_bond_constraint_is_unique(
        #[case] c: DativeBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_unique(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_count_lit(DativeBondConstraint::RingCount(ValueAst::Lit(1)), false)]
    #[case::ring_count_undetermined(DativeBondConstraint::RingCount(ValueAst::Undetermined), true)]
    #[case::ring_size_undetermined(DativeBondConstraint::RingSize(ValueAst::Undetermined), true)]
    fn test_dative_bond_constraint_is_undetermined(
        #[case] c: DativeBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rstest]
    fn test_dative_bond_constraints_add_unique_replaces() {
        let mut cs = DativeBondConstraints::new();
        cs.add(DativeBondConstraint::RingCount(ValueAst::Lit(1)));
        let prev = cs.add(DativeBondConstraint::RingCount(ValueAst::Lit(2)));
        assert_eq!(
            prev,
            Some(DativeBondConstraint::RingCount(ValueAst::Lit(1)))
        );
        assert_eq!(
            cs.as_slice(),
            &[DativeBondConstraint::RingCount(ValueAst::Lit(2))]
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_add_distinct_kinds_coexist() {
        let mut cs = DativeBondConstraints::new();
        cs.add(DativeBondConstraint::RingCount(ValueAst::Lit(1)));
        cs.add(DativeBondConstraint::RingSize(ValueAst::Lit(6)));
        assert_eq!(cs.len(), 2);
    }

    #[rstest]
    fn test_dative_bond_constraints_retain() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::RingCount(ValueAst::Lit(2)),
            DativeBondConstraint::RingSize(ValueAst::Lit(6)),
        ]);
        cs.retain(|c| matches!(c, DativeBondConstraint::RingCount(_)));
        assert_eq!(
            cs.as_slice(),
            &[DativeBondConstraint::RingCount(ValueAst::Lit(2))]
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_clear() {
        let mut cs =
            DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(1))]);
        cs.clear();
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_dative_bond_constraints_remap_is_identity_for_value_only() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::RingCount(ValueAst::Lit(1)),
            DativeBondConstraint::RingSize(ValueAst::Lit(6)),
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
        let after = cs.clone().remap(&remap);
        assert_eq!(after, cs);
    }
}
