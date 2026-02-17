//! Multi-center bond representation for GraphIR.
//!
//! Multi-center bonds generarlize bonds to involve more than two atoms, where electrons are shared
//! or donated across multiple centers simultaneously. Each multicenter bond consists of one or more
//! multicenter sets, each contributing a fixed number of electrons. Multicenter sets are sorted and
//! deduplicated by atom index. Atom indices use `AtomIndex` from the molecule module.


use smallvec::SmallVec;

use super::molecule::AtomIndex;

/// A set of atoms contributing to a multi-center bond.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterSet {
    atoms: SmallVec<[AtomIndex; 8]>,
    electrons: u8,
}

impl MulticenterSet {
    pub fn new<I>(atoms: I, electrons: u8) -> Self
    where
        I: IntoIterator<Item = AtomIndex>,
    {
        let mut atoms: SmallVec<[AtomIndex; 8]> = atoms.into_iter().collect();
        atoms.sort_unstable();
        atoms.dedup();
        Self { atoms, electrons }
    }

    pub fn single(atom: AtomIndex, electrons: u8) -> Self {
        Self {
            atoms: SmallVec::from_elem(atom, 1),
            electrons,
        }
    }

    pub fn atoms(&self) -> &[AtomIndex] {
        self.atoms.as_slice()
    }

    pub fn electrons(&self) -> u8 {
        self.electrons
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn electron_count(&self) -> u8 {
        self.electrons
    }

    pub fn contains_atom(&self, atom: AtomIndex) -> bool {
        self.atoms.binary_search(&atom).is_ok()
    }
}

/// A multi-center bond involving three or more atoms.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterBond {
    contributions: SmallVec<[MulticenterSet; 4]>,
}

impl MulticenterBond {
    pub fn new<I>(contributions: I) -> Self
    where
        I: IntoIterator<Item = MulticenterSet>,
    {
        let mut contributions: SmallVec<[MulticenterSet; 4]> =
            contributions.into_iter().collect();
        contributions.sort_unstable();
        contributions.dedup();
        Self { contributions }
    }

    pub fn contributions(&self) -> &[MulticenterSet] {
        self.contributions.as_slice()
    }

    pub fn electron_count(&self) -> u8 {
        self.contributions
            .iter()
            .map(MulticenterSet::electrons)
            .sum()
    }

    pub fn atom_count(&self) -> usize {
        self.contributions
            .iter()
            .map(MulticenterSet::atom_count)
            .sum()
    }

    pub fn all_atoms(&self) -> Vec<AtomIndex> {
        self.contributions
            .iter()
            .flat_map(|c| c.atoms().iter().copied())
            .collect()
    }
}
