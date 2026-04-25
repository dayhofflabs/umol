//! Per-aromatic-system constraints.
//!
//! All previous variants (`Atoms`, `Contains`, `ContainsAll`, `AllAtoms`,
//! `AnyAtom`) were atom-ref-bearing or carried a delegated atom predicate;
//! those moved to `RelationalConstraint` at molecule scope. The enum is
//! kept (empty for now) so future value-only aromatic-system constraints
//! can be added here without reshaping the AST or DSL surface.

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
        std::mem::take(&mut self.0).into_iter()
    }

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
