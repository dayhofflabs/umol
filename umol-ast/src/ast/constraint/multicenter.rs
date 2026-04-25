//! Per-multicenter-bond constraints.

use std::mem::{self, replace};
use std::slice::Iter;

use super::super::remap::IdxRemapping;
use super::super::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum MulticenterBondConstraint {
    /// Asserted total electron count for the multicenter bond. Cross-checked
    /// by the `ConsistencyValidator` against `sum(MulticenterBondAst::electrons)`.
    ElectronCount(ValueAst),
}

impl MulticenterBondConstraint {
    pub fn is_unique(&self) -> bool {
        match self {
            Self::ElectronCount(_) => true,
        }
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

    pub fn remap(self, _remap: &IdxRemapping) -> Option<Self> {
        match self {
            Self::ElectronCount(_) => Some(self),
        }
    }

    fn matches_kind(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (Self::ElectronCount(_), Self::ElectronCount(_)),
        )
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

    pub fn iter(&self) -> Iter<'_, MulticenterBondConstraint> {
        self.0.iter()
    }

    pub fn add(&mut self, c: MulticenterBondConstraint) -> Option<MulticenterBondConstraint> {
        if c.is_unique() {
            if let Some(i) = self.0.iter().position(|x| x.matches_kind(&c)) {
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

    pub fn remap(self, remap: &IdxRemapping) -> Self {
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

#[cfg(test)]
mod tests {
    use std::iter::empty;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Remapping;

    use super::*;

    fn empty_remapping() -> IdxRemapping {
        IdxRemapping::new(
            Remapping {
                removed_nodes: Vec::new(),
                removed_edges: Vec::new(),
            },
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    }

    #[rstest]
    fn test_multicenter_bond_constraints_empty_methods() {
        let mut cs = MulticenterBondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert!(cs.as_slice().is_empty());
        assert_eq!(cs.iter().count(), 0);
        cs.retain(|_| true);
        assert!(cs.is_empty());
        let drained: Vec<_> = cs.take().collect();
        assert!(drained.is_empty());
        cs.clear();
        assert!(cs.is_empty());
        let from_empty: MulticenterBondConstraints = empty().collect();
        assert!(from_empty.is_empty());
        let remapped = MulticenterBondConstraints::new().remap(&empty_remapping());
        assert!(remapped.is_empty());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_add_inserts_new() {
        let mut cs = MulticenterBondConstraints::new();
        let displaced = cs.add(MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2)));
        assert_eq!(displaced, None);
        assert_eq!(cs.len(), 1);
        assert_eq!(
            cs.as_slice()[0],
            MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2)),
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_add_replaces_existing_unique_kind() {
        let mut cs = MulticenterBondConstraints::new();
        cs.add(MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2)));
        let displaced =
            cs.add(MulticenterBondConstraint::ElectronCount(ValueAst::Lit(4)));
        assert_eq!(
            displaced,
            Some(MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2))),
        );
        assert_eq!(cs.len(), 1);
        assert_eq!(
            cs.as_slice()[0],
            MulticenterBondConstraint::ElectronCount(ValueAst::Lit(4)),
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraint_is_undetermined() {
        assert!(MulticenterBondConstraint::ElectronCount(ValueAst::Undetermined).is_undetermined());
        assert!(!MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2)).is_undetermined());
    }

    #[rstest]
    fn test_multicenter_bond_constraints_simplify_each() {
        use crate::ast::value::Expr;
        let mut cs = MulticenterBondConstraints::new();
        cs.add(MulticenterBondConstraint::ElectronCount(ValueAst::Expr(
            Expr::Lit(2),
        )));
        cs.simplify_each();
        assert_eq!(
            cs.as_slice()[0],
            MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2)),
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraints_remap_passes_through() {
        let mut cs = MulticenterBondConstraints::new();
        cs.add(MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2)));
        let remapped = cs.remap(&empty_remapping());
        assert_eq!(remapped.len(), 1);
    }

    #[rstest]
    fn test_multicenter_bond_constraints_from_iter_uses_add() {
        let cs: MulticenterBondConstraints = [
            MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2)),
            MulticenterBondConstraint::ElectronCount(ValueAst::Lit(4)),
        ]
        .into_iter()
        .collect();
        assert_eq!(cs.len(), 1);
        assert_eq!(
            cs.as_slice()[0],
            MulticenterBondConstraint::ElectronCount(ValueAst::Lit(4)),
        );
    }
}
