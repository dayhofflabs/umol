//! Per-dative-bond constraints.

use std::mem;
use std::slice::Iter;

use strum::EnumDiscriminants;

use super::atom::AtomConstraint;
use crate::ast::idx::{AtomIdx, BondIdx};
use crate::ast::remap::IdxRemapping;
use crate::ast::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(DativeBondConstraintKind), derive(Hash))]
pub enum DativeBondConstraint {
    RingCount(ValueAst),
    RingSize(ValueAst),
    Donor(AtomIdx),
    Acceptor(AtomIdx),
    DonorSatisfies(Box<AtomConstraint>),
    AcceptorSatisfies(Box<AtomConstraint>),
    Parallels(BondIdx),
}

impl DativeBondConstraint {
    pub fn kind(&self) -> DativeBondConstraintKind {
        self.into()
    }

    /// Single-valued per dative bond: `RingCount`, `RingSize`, `Donor`,
    /// `Acceptor`, `Parallels` (a dative bond has at most one donor, one
    /// acceptor, and one parallel covalent bond). Multi-valued:
    /// `DonorSatisfies`, `AcceptorSatisfies` (independent predicates AND
    /// together).
    pub fn is_unique(&self) -> bool {
        matches!(
            self,
            Self::RingCount(_)
                | Self::RingSize(_)
                | Self::Donor(_)
                | Self::Acceptor(_)
                | Self::Parallels(_)
        )
    }

    /// `Donor` / `Acceptor` / `Parallels` carry entity refs and are never
    /// undetermined. `RingCount` / `RingSize` are undetermined iff their
    /// inner value is. `DonorSatisfies` / `AcceptorSatisfies` delegate to
    /// the inner atom constraint.
    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Donor(_) | Self::Acceptor(_) | Self::Parallels(_) => false,
            Self::RingCount(v) | Self::RingSize(v) => v.is_undetermined(),
            Self::DonorSatisfies(c) | Self::AcceptorSatisfies(c) => c.is_undetermined(),
        }
    }

    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        match self {
            Self::RingCount(v) => Some(Self::RingCount(v)),
            Self::RingSize(v) => Some(Self::RingSize(v)),
            Self::Donor(a) => remap.atom(a).map(Self::Donor),
            Self::Acceptor(a) => remap.atom(a).map(Self::Acceptor),
            Self::DonorSatisfies(c) => Some(Self::DonorSatisfies(c)),
            Self::AcceptorSatisfies(c) => Some(Self::AcceptorSatisfies(c)),
            Self::Parallels(b) => remap.bond(b).map(Self::Parallels),
        }
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

    fn idx_remapping(removed_nodes: Vec<u32>, removed_edges: Vec<u32>) -> IdxRemapping {
        IdxRemapping::new(
            umol_graph_core::Remapping {
                removed_nodes,
                removed_edges,
            },
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_count(DativeBondConstraint::RingCount(ValueAst::Lit(1)), DativeBondConstraintKind::RingCount)]
    #[case::donor(DativeBondConstraint::Donor(AtomIdx(3)), DativeBondConstraintKind::Donor)]
    #[case::parallels(DativeBondConstraint::Parallels(BondIdx(2)), DativeBondConstraintKind::Parallels)]
    #[case::donor_satisfies(DativeBondConstraint::DonorSatisfies(Box::new(AtomConstraint::Valence(ValueAst::Lit(3)))), DativeBondConstraintKind::DonorSatisfies)]
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
    #[case::donor(DativeBondConstraint::Donor(AtomIdx(3)), true)]
    #[case::acceptor(DativeBondConstraint::Acceptor(AtomIdx(4)), true)]
    #[case::parallels(DativeBondConstraint::Parallels(BondIdx(2)), true)]
    #[case::donor_satisfies(DativeBondConstraint::DonorSatisfies(Box::new(AtomConstraint::Valence(ValueAst::Lit(3)))), false)]
    #[case::acceptor_satisfies(DativeBondConstraint::AcceptorSatisfies(Box::new(AtomConstraint::Degree(ValueAst::Lit(4)))), false)]
    fn test_dative_bond_constraint_is_unique(
        #[case] c: DativeBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_unique(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::donor(DativeBondConstraint::Donor(AtomIdx(3)), false)]
    #[case::acceptor(DativeBondConstraint::Acceptor(AtomIdx(4)), false)]
    #[case::parallels(DativeBondConstraint::Parallels(BondIdx(2)), false)]
    #[case::ring_count_lit(DativeBondConstraint::RingCount(ValueAst::Lit(1)), false)]
    #[case::ring_count_undetermined(DativeBondConstraint::RingCount(ValueAst::Undetermined), true)]
    #[case::ring_size_undetermined(DativeBondConstraint::RingSize(ValueAst::Undetermined), true)]
    #[case::donor_satisfies_lit(DativeBondConstraint::DonorSatisfies(Box::new(AtomConstraint::Valence(ValueAst::Lit(3)))), false)]
    #[case::donor_satisfies_undetermined(DativeBondConstraint::DonorSatisfies(Box::new(AtomConstraint::Valence(ValueAst::Undetermined))), true)]
    fn test_dative_bond_constraint_is_undetermined(
        #[case] c: DativeBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rstest]
    fn test_dative_bond_constraints_add_unique_replaces() {
        let mut cs = DativeBondConstraints::new();
        cs.add(DativeBondConstraint::Donor(AtomIdx(1)));
        let prev = cs.add(DativeBondConstraint::Donor(AtomIdx(2)));
        assert_eq!(prev, Some(DativeBondConstraint::Donor(AtomIdx(1))));
        assert_eq!(cs.as_slice(), &[DativeBondConstraint::Donor(AtomIdx(2))]);
    }

    #[rstest]
    fn test_dative_bond_constraints_add_multi_appends() {
        let mut cs = DativeBondConstraints::new();
        let c1 = DativeBondConstraint::DonorSatisfies(Box::new(AtomConstraint::Valence(
            ValueAst::Lit(3),
        )));
        let c2 = DativeBondConstraint::DonorSatisfies(Box::new(AtomConstraint::Degree(
            ValueAst::Lit(2),
        )));
        assert_eq!(cs.add(c1.clone()), None);
        assert_eq!(cs.add(c2.clone()), None);
        assert_eq!(cs.as_slice(), &[c1, c2]);
    }

    #[rstest]
    fn test_dative_bond_constraints_add_distinct_kinds_coexist() {
        let mut cs = DativeBondConstraints::new();
        cs.add(DativeBondConstraint::Donor(AtomIdx(1)));
        cs.add(DativeBondConstraint::Acceptor(AtomIdx(2)));
        cs.add(DativeBondConstraint::RingCount(ValueAst::Lit(1)));
        assert_eq!(cs.len(), 3);
    }

    #[rstest]
    fn test_dative_bond_constraints_retain() {
        let mut cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Donor(AtomIdx(1)),
            DativeBondConstraint::RingCount(ValueAst::Lit(2)),
            DativeBondConstraint::RingSize(ValueAst::Lit(6)),
        ]);
        cs.retain(|c| {
            matches!(
                c,
                DativeBondConstraint::RingCount(_) | DativeBondConstraint::RingSize(_)
            )
        });
        assert_eq!(
            cs.as_slice(),
            &[
                DativeBondConstraint::RingCount(ValueAst::Lit(2)),
                DativeBondConstraint::RingSize(ValueAst::Lit(6)),
            ]
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_clear() {
        let mut cs = DativeBondConstraints::from_iter([DativeBondConstraint::Donor(AtomIdx(0))]);
        cs.clear();
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_dative_bond_constraints_remap_shifts_atom_ref() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Donor(AtomIdx(3)),
            DativeBondConstraint::RingCount(ValueAst::Lit(1)),
        ]);
        let remap = idx_remapping(vec![1], vec![]);
        let after = cs.remap(&remap);
        assert_eq!(
            after.as_slice(),
            &[
                DativeBondConstraint::Donor(AtomIdx(2)),
                DativeBondConstraint::RingCount(ValueAst::Lit(1)),
            ]
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_remap_drops_removed_atom() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Donor(AtomIdx(1)),
            DativeBondConstraint::RingCount(ValueAst::Lit(2)),
        ]);
        let remap = idx_remapping(vec![1], vec![]);
        let after = cs.remap(&remap);
        assert_eq!(
            after.as_slice(),
            &[DativeBondConstraint::RingCount(ValueAst::Lit(2))]
        );
    }

    #[rstest]
    fn test_dative_bond_constraints_remap_drops_removed_bond() {
        let cs = DativeBondConstraints::from_iter([
            DativeBondConstraint::Parallels(BondIdx(2)),
            DativeBondConstraint::RingCount(ValueAst::Lit(1)),
        ]);
        let remap = idx_remapping(vec![], vec![2]);
        let after = cs.remap(&remap);
        assert_eq!(
            after.as_slice(),
            &[DativeBondConstraint::RingCount(ValueAst::Lit(1))]
        );
    }
}
