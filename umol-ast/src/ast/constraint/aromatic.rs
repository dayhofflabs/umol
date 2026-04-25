//! Per-aromatic-system constraints.
//!
//! All previous variants (`Atoms`, `Contains`, `ContainsAll`, `AllAtoms`,
//! `AnyAtom`) were atom-ref-bearing or carried a delegated atom predicate;
//! those moved to `RelationalConstraint` at molecule scope. The enum is
//! kept (empty for now) so future value-only aromatic-system constraints
//! can be added here without reshaping the AST or DSL surface.

use std::mem;
use std::slice::Iter;

use super::super::remap::IdxRemapping;

/// Aromatic-system-scope constraint. Currently uninhabited — placeholder
/// for future value-only variants. Atom-ref and quantified-predicate forms
/// live at molecule scope via `RelationalConstraint`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticSystemConstraint {}

impl AromaticSystemConstraint {
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

/// Per-aromatic-system constraint container. Empty in practice until new
/// value-only variants land on `AromaticSystemConstraint`.
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

    pub fn add(&mut self, c: AromaticSystemConstraint) -> Option<AromaticSystemConstraint> {
        match c {}
    }

    pub fn retain(&mut self, mut f: impl FnMut(&AromaticSystemConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = AromaticSystemConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    /// No-op: the inner enum is uninhabited so the store has no values to
    /// simplify. Kept for API symmetry with the inhabited containers.
    pub fn simplify_each(&mut self) {}

    pub fn remap(self, _remap: &IdxRemapping) -> Self {
        // No inhabitants → vec is always empty → no-op.
        self
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

    /// Exercise the container methods on an empty store. The constraint
    /// enum is uninhabited so every entry-bearing method (`add`,
    /// `from_iter` body, etc.) is structurally unreachable at runtime; the
    /// container's bookkeeping methods are still callable and must work.
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
}
