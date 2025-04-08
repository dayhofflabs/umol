// Traits and structs for pattern matching

use crate::core::types::{AtomIndex, BondIndex};
use super::{Atom, Bond, Molecule, Query};

/// A pattern matcher for finding subgraphs in molecules
pub struct PatternMatcher {
    query: Query,
}

impl PatternMatcher {
    /// Create a new pattern matcher from a query
    pub fn new(query: Query) -> Self {
        Self { query }
    }

    /// Find all matches of the query in a molecule
    pub fn find_matches(&self, molecule: &Molecule) -> Vec<Vec<AtomIndex>> {
        // TODO: Implement subgraph isomorphism algorithm
        Vec::new()
    }
}
