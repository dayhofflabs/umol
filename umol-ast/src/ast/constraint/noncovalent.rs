//! Per-noncovalent-bond constraints.
//!
//! All previous variants (`Ends`, `Contains`, `EndsSatisfy`) were atom-ref-
//! bearing or carried a delegated atom predicate; those moved to
//! `RelationalConstraint` at molecule scope. The enum is kept (empty for
//! now) so future value-only noncovalent-bond constraints can be added here
//! without reshaping the AST or DSL surface.

use std::mem;
use std::slice::Iter;

use super::super::remap::IdxRemapping;

/// Noncovalent-bond-scope constraint. Currently uninhabited — placeholder
/// for future value-only variants. Atom-ref and quantified-predicate forms
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
