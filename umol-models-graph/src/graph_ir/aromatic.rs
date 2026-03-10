//! Aromatic system representation for GraphIR.
//!
//! An aromatic system is a molecule-level object describing a set of atoms
//! participating in a delocalized π system, each contributing a fixed number
//! of electrons (`aromatic_valence`). The system as a whole is subject to
//! validation (ring membership, Kekulé feasibility, optional Hückel 4n+2).

use super::molecule::AtomIndex;

/// Per-atom contribution to an aromatic system.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AromaticContribution {
    atom: AtomIndex,
    aromatic_valence: u8,
}

impl AromaticContribution {
    pub fn new(atom: AtomIndex, aromatic_valence: u8) -> Self {
        Self {
            atom,
            aromatic_valence,
        }
    }

    pub fn atom(&self) -> AtomIndex {
        self.atom
    }

    pub fn aromatic_valence(&self) -> u8 {
        self.aromatic_valence
    }
}

// TODO: Add charge and multiplicity fields for delocalized charge and spin, not assignable to individual atoms.
/// An aromatic system: a set of atoms each contributing electrons to a
/// delocalized π system. Contributions are canonicalized by atom index;
/// each atom appears at most once.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AromaticSystem {
    contributions: Vec<AromaticContribution>,
}

impl AromaticSystem {
    pub fn new<I>(contributions: I) -> Self
    where
        I: IntoIterator<Item = AromaticContribution>,
    {
        let mut contributions: Vec<AromaticContribution> = contributions.into_iter().collect();
        contributions.sort_unstable();
        contributions.dedup_by_key(|c| c.atom);
        Self { contributions }
    }

    pub fn contributions(&self) -> &[AromaticContribution] {
        &self.contributions
    }

    pub fn atom_count(&self) -> usize {
        self.contributions.len()
    }

    pub fn electron_count(&self) -> u8 {
        self.contributions.iter().map(|c| c.aromatic_valence).sum()
    }

    pub fn contains_atom(&self, atom: AtomIndex) -> bool {
        self.contributions
            .binary_search_by_key(&atom, |c| c.atom)
            .is_ok()
    }

    pub fn atoms(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.contributions.iter().map(|c| c.atom)
    }
}
