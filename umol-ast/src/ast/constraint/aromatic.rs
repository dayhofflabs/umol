//! Per-aromatic-system constraints.

use std::mem::{self, replace};
use std::slice::Iter;

use super::super::remap::IdxRemapping;
use super::super::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticSystemConstraint {
    /// Asserted total π-electron count for the system. Cross-checked by the
    /// `ConsistencyValidator` against `sum(AromaticSystemAst::electrons)`.
    ElectronCount(ValueAst),
}

impl AromaticSystemConstraint {
    /// Each kind admits at most one entry per system. Used by the container
    /// to decide replace-vs-insert on `add`.
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

    pub fn iter(&self) -> Iter<'_, AromaticSystemConstraint> {
        self.0.iter()
    }

    /// Insert a constraint, replacing any existing entry whose kind shadows
    /// it (single-valued kinds report `is_unique`). Returns the displaced
    /// entry.
    pub fn add(&mut self, c: AromaticSystemConstraint) -> Option<AromaticSystemConstraint> {
        if c.is_unique() {
            if let Some(i) = self.0.iter().position(|x| x.matches_kind(&c)) {
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
            *c = mem::replace(c, AromaticSystemConstraint::ElectronCount(ValueAst::Undetermined))
                .simplify();
        }
    }

    pub fn remap(self, remap: &IdxRemapping) -> Self {
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
    fn test_aromatic_system_constraints_empty_methods() {
        let mut cs = AromaticSystemConstraints::new();
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
        let from_empty: AromaticSystemConstraints = empty().collect();
        assert!(from_empty.is_empty());
        let remapped = AromaticSystemConstraints::new().remap(&empty_remapping());
        assert!(remapped.is_empty());
    }

    #[rstest]
    fn test_aromatic_system_constraints_add_inserts_new() {
        let mut cs = AromaticSystemConstraints::new();
        let displaced = cs.add(AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6)));
        assert_eq!(displaced, None);
        assert_eq!(cs.len(), 1);
        assert_eq!(
            cs.as_slice()[0],
            AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6)),
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_add_replaces_existing_unique_kind() {
        let mut cs = AromaticSystemConstraints::new();
        cs.add(AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6)));
        let displaced = cs.add(AromaticSystemConstraint::ElectronCount(ValueAst::Lit(10)));
        assert_eq!(
            displaced,
            Some(AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))),
        );
        assert_eq!(cs.len(), 1);
        assert_eq!(
            cs.as_slice()[0],
            AromaticSystemConstraint::ElectronCount(ValueAst::Lit(10)),
        );
    }

    #[rstest]
    fn test_aromatic_system_constraint_is_undetermined() {
        assert!(AromaticSystemConstraint::ElectronCount(ValueAst::Undetermined).is_undetermined());
        assert!(!AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6)).is_undetermined());
    }

    #[rstest]
    fn test_aromatic_system_constraints_simplify_each() {
        use crate::ast::value::Expr;
        let mut cs = AromaticSystemConstraints::new();
        cs.add(AromaticSystemConstraint::ElectronCount(ValueAst::Expr(
            Expr::Lit(6),
        )));
        cs.simplify_each();
        assert_eq!(
            cs.as_slice()[0],
            AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6)),
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_remap_passes_through() {
        let mut cs = AromaticSystemConstraints::new();
        cs.add(AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6)));
        let remapped = cs.remap(&empty_remapping());
        assert_eq!(remapped.len(), 1);
        assert_eq!(
            remapped.as_slice()[0],
            AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6)),
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_from_iter_uses_add() {
        let cs: AromaticSystemConstraints = [
            AromaticSystemConstraint::ElectronCount(ValueAst::Lit(2)),
            AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6)),
        ]
        .into_iter()
        .collect();
        assert_eq!(cs.len(), 1);
        assert_eq!(
            cs.as_slice()[0],
            AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6)),
        );
    }
}
