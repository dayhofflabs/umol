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

// TODO: Add multiplicity field
/// An aromatic system consisting of atoms and number of electrons contributed
/// Charge is delocalized charge, not assignable to any individual atom
/// Each atom can participate in at most one aromatic system, appears only once
/// in the contributions list
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AromaticSystem {
    contributions: Vec<AromaticContribution>,
    /// Delocalized charge
    charge: i8,
    /// List of atomic indices per ring
    rings: Vec<Vec<AtomIndex>>,
}

impl AromaticSystem {
    pub fn new<I>(contributions: I) -> Self
    where
        I: IntoIterator<Item = AromaticContribution>,
    {
        let mut contributions: Vec<AromaticContribution> = contributions.into_iter().collect();
        contributions.sort_unstable();
        contributions.dedup_by_key(|c| c.atom);
        Self {
            contributions,
            charge: 0,
            rings: Vec::new(),
        }
    }

    pub fn with_rings<I>(contributions: I, rings: Vec<Vec<AtomIndex>>) -> Self
    where
        I: IntoIterator<Item = AromaticContribution>,
    {
        let mut contributions: Vec<AromaticContribution> = contributions.into_iter().collect();
        contributions.sort_unstable();
        contributions.dedup_by_key(|c| c.atom);
        Self {
            contributions,
            charge: 0,
            rings,
        }
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

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn set_charge(&mut self, charge: i8) {
        self.charge = charge;
    }

    pub fn rings(&self) -> &[Vec<AtomIndex>] {
        &self.rings
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
