//! Per-multicenter-bond constraints.
//!
//! All previous variants (`Atoms`, `Contains`, `ContainsAll`, `AllAtoms`,
//! `AnyAtom`) were atom-ref-bearing or carried a delegated atom predicate;
//! those moved to `RelationalConstraint` at molecule scope. The enum is
//! kept (empty for now) so future value-only multicenter-bond constraints
//! can be added here without reshaping the AST or DSL surface.

use std::mem;
use std::slice::Iter;

use super::super::remap::IdxRemapping;

/// Multicenter-bond-scope constraint. Currently uninhabited — placeholder
/// for future value-only variants. Atom-ref and quantified-predicate forms
/// live at molecule scope via `RelationalConstraint`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum MulticenterBondConstraint {}

impl MulticenterBondConstraint {
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

/// Per-multicenter-bond constraint container. Empty in practice until new
/// value-only variants land on `MulticenterBondConstraint`.
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
        match c {}
    }

    pub fn retain(&mut self, mut f: impl FnMut(&MulticenterBondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Move the entries out of the store, leaving it empty.
    pub fn take(&mut self) -> impl Iterator<Item = MulticenterBondConstraint> {
        mem::take(&mut self.0).into_iter()
    }

    pub fn remap(self, _remap: &IdxRemapping) -> Self {
        self
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
