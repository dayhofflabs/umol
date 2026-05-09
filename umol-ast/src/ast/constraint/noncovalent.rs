//! Noncovalent bond constraints.

use std::mem;
use std::slice::Iter;

use super::super::remap::IdxRemapping;

/// Noncovalent-bond-scope constraint. Atom-ref and quantified-predicate forms
/// live at molecule scope via `RelationalConstraint`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondConstraint {}

impl NoncovalentBondConstraint {
    pub fn is_unique(&self) -> bool {
        match *self {}
    }

    pub fn is_undetermined(&self) -> bool {
        match *self {}
    }

    pub fn remap(self, _remap: &IdxRemapping) -> Option<Self> {
        match self {}
    }
}

/// Per-noncovalent-bond constraint container. Empty in practice until new
/// value-only variants land on `NoncovalentBondConstraint`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct NoncovalentBondConstraints(Vec<NoncovalentBondConstraint>);

impl NoncovalentBondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[NoncovalentBondConstraint] {
        &self.0
    }

    pub fn iter(&self) -> Iter<'_, NoncovalentBondConstraint> {
        self.0.iter()
    }

    pub fn add(&mut self, c: NoncovalentBondConstraint) -> Option<NoncovalentBondConstraint> {
        match c {}
    }

    pub fn retain(&mut self, mut f: impl FnMut(&NoncovalentBondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = NoncovalentBondConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    /// No-op: the inner enum is uninhabited so the store has no values to
    /// simplify. Kept for API symmetry with the inhabited containers.
    pub fn simplify_each(&mut self) {}

    pub fn remap(self, _remap: &IdxRemapping) -> Self {
        self
    }
}

impl FromIterator<NoncovalentBondConstraint> for NoncovalentBondConstraints {
    fn from_iter<I: IntoIterator<Item = NoncovalentBondConstraint>>(iter: I) -> Self {
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

    #[rstest]
    fn test_noncovalent_bond_constraints_new() {
        let cs = NoncovalentBondConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[NoncovalentBondConstraint]);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_iter() {
        let cs = NoncovalentBondConstraints::new();
        assert_eq!(cs.iter().count(), 0);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_retain() {
        let mut cs = NoncovalentBondConstraints::new();
        cs.retain(|_| true);
        assert_eq!(cs, NoncovalentBondConstraints::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_clear() {
        let mut cs = NoncovalentBondConstraints::new();
        cs.clear();
        assert_eq!(cs, NoncovalentBondConstraints::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_take() {
        let mut cs = NoncovalentBondConstraints::new();
        let drained: Vec<_> = cs.take().collect();
        assert!(drained.is_empty());
        assert_eq!(cs, NoncovalentBondConstraints::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_simplify_each() {
        let mut cs = NoncovalentBondConstraints::new();
        cs.simplify_each();
        assert_eq!(cs, NoncovalentBondConstraints::new());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_remap() {
        let cs = NoncovalentBondConstraints::new();
        let remap = IdxRemapping::new(
            Remapping {
                removed_nodes: Vec::new(),
                removed_edges: Vec::new(),
            },
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        assert_eq!(cs.clone().remap(&remap), cs);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_from_iter() {
        let cs: NoncovalentBondConstraints = empty().collect();
        assert_eq!(cs, NoncovalentBondConstraints::new());
    }
}
