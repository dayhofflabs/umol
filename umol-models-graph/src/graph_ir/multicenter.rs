//! Multi-center bond representation for GraphIR.
//!
//! A multicenter bond consists of one or more multicenter sets. Each set contains per-atom
//! contributions (atom index + optional electron count) and a system-level charge.
//! The two-level structure preserves the set grouping from input formats (e.g. center vs ligand
//! sets in CXSMILES) while allowing per-atom electron data to be populated during resolution.
//!
//! Electron counts are `Option<u8>`: `None` when not yet resolved (topology-only),
//! `Some(n)` after the valence phase populates them.

use smallvec::SmallVec;

use super::molecule::AtomIndex;

/// Per-atom contribution to a multicenter set.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterContribution {
    atom: AtomIndex,
    electrons: Option<u8>,
}

impl MulticenterContribution {
    pub fn new(atom: AtomIndex, electrons: Option<u8>) -> Self {
        Self { atom, electrons }
    }

    pub fn topology_only(atom: AtomIndex) -> Self {
        Self {
            atom,
            electrons: None,
        }
    }

    pub fn atom(&self) -> AtomIndex {
        self.atom
    }

    pub fn electrons(&self) -> Option<u8> {
        self.electrons
    }

    pub fn is_resolved(&self) -> bool {
        self.electrons.is_some()
    }
}

/// A set of atoms participating in a multicenter bond, with a system-level charge.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterSet {
    contributions: SmallVec<[MulticenterContribution; 8]>,
    charge: i8,
}

impl MulticenterSet {
    pub fn new<I>(contributions: I, charge: i8) -> Self
    where
        I: IntoIterator<Item = MulticenterContribution>,
    {
        let mut contributions: SmallVec<[MulticenterContribution; 8]> =
            contributions.into_iter().collect();
        contributions.sort_unstable();
        contributions.dedup_by_key(|c| c.atom);
        let result = Self {
            contributions,
            charge,
        };
        debug_assert!(
            !result.is_resolved() || result.electron_sum() as i16 - charge as i16 >= 0,
            "negative total electron count"
        );
        result
    }

    pub fn topology_only<I>(atoms: I) -> Self
    where
        I: IntoIterator<Item = AtomIndex>,
    {
        Self::new(
            atoms
                .into_iter()
                .map(MulticenterContribution::topology_only),
            0,
        )
    }

    pub fn contributions(&self) -> &[MulticenterContribution] {
        self.contributions.as_slice()
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn atom_count(&self) -> usize {
        self.contributions.len()
    }

    pub fn is_resolved(&self) -> bool {
        self.contributions.iter().all(|c| c.is_resolved())
    }

    pub fn contains_atom(&self, atom: AtomIndex) -> bool {
        self.contributions
            .binary_search_by_key(&atom, |c| c.atom)
            .is_ok()
    }

    pub fn atoms(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.contributions.iter().map(|c| c.atom)
    }

    fn electron_sum(&self) -> u8 {
        self.contributions
            .iter()
            .filter_map(|c| c.electrons)
            .sum()
    }

    pub fn electron_count(&self) -> u8 {
        let sum = self.electron_sum() as i16 - self.charge as i16;
        sum as u8
    }
}

/// A multi-center bond consisting of one or more multicenter sets.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterBond {
    sets: SmallVec<[MulticenterSet; 4]>,
}

impl MulticenterBond {
    pub fn new<I>(sets: I) -> Self
    where
        I: IntoIterator<Item = MulticenterSet>,
    {
        let mut sets: SmallVec<[MulticenterSet; 4]> = sets.into_iter().collect();
        sets.sort_unstable();
        sets.dedup();
        Self { sets }
    }

    pub fn sets(&self) -> &[MulticenterSet] {
        self.sets.as_slice()
    }

    pub fn is_resolved(&self) -> bool {
        self.sets.iter().all(|s| s.is_resolved())
    }

    pub fn electron_count(&self) -> u8 {
        self.sets.iter().map(MulticenterSet::electron_count).sum()
    }

    pub fn atom_count(&self) -> usize {
        self.sets.iter().map(MulticenterSet::atom_count).sum()
    }

    pub fn contains_atom(&self, atom: AtomIndex) -> bool {
        self.sets.iter().any(|s| s.contains_atom(atom))
    }

    pub fn all_atoms(&self) -> Vec<AtomIndex> {
        self.sets
            .iter()
            .flat_map(|s| s.atoms())
            .collect()
    }
}
