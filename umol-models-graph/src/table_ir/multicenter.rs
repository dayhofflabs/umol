//! Multi-center bond representation for TableIR.
//!
//! Multi-center bonds describe bonding interactions involving more than two atoms,
//! where electrons are shared or donated across multiple centers simultaneously.
//!
//! # Examples
//!
//! **Haptic bonding (ferrocene η⁵-Cp→Fe):**
//! One Cp ring donates 6 electrons to Fe, which contributes 0:
//! ```text
//! contributions: [([0,1,2,3,4], 6), ([10], 0)]
//! ```
//!
//! **Electron-deficient 3c-2e (diborane B-H-B bridge):**
//! Two B atoms and one bridging H share 2 electrons total:
//! ```text
//! contributions: [([0], 1), ([1], 0), ([5], 1)]
//! ```
//!
//! **Electron-rich 3c-4e (XeF₂):**
//! Xe contributes a lone pair, each F contributes one electron:
//! ```text
//! contributions: [([0], 2), ([1], 1), ([2], 1)]
//! ```

/// A set of atoms contributing to a multi-center bond.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MulticenterSet {
    /// Atom indices participating in this contribution.
    pub atoms: Vec<u32>,

    /// Number of electrons contributed by this group.
    pub electrons: u8,
}

impl MulticenterSet {
    pub fn new(atoms: Vec<u32>, electrons: u8) -> Self {
        Self { atoms, electrons }
    }

    pub fn single(atom: u32, electrons: u8) -> Self {
        Self {
            atoms: vec![atom],
            electrons,
        }
    }
}

/// A multi-center bond involving three or more atoms.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MulticenterBond {
    pub contributions: Vec<MulticenterSet>,
}

impl MulticenterBond {
    pub fn new(contributions: Vec<MulticenterSet>) -> Self {
        Self { contributions }
    }

    /// Number of electrons in this multi-center bond.
    pub fn electron_count(&self) -> u8 {
        self.contributions.iter().map(|c| c.electrons).sum()
    }

    /// Total number of atoms involved.
    pub fn atom_count(&self) -> usize {
        self.contributions.iter().map(|c| c.atoms.len()).sum()
    }

    /// All atom indices involved in this bond (flattened).
    pub fn all_atoms(&self) -> Vec<u32> {
        self.contributions
            .iter()
            .flat_map(|c| c.atoms.iter().copied())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::ferrocene(MulticenterBond::new(vec![MulticenterSet::new(vec![0, 1, 2, 3, 4], 6), MulticenterSet::single(9, 0)]), 6, 6, vec![0, 1, 2, 3, 4, 9])]
    #[case::diborane(MulticenterBond::new(vec![MulticenterSet::single(0, 1), MulticenterSet::single(1, 0), MulticenterSet::single(5, 1)]), 2, 3, vec![0, 1, 5])]
    #[case::xef2(MulticenterBond::new(vec![MulticenterSet::single(0, 2), MulticenterSet::single(1, 1), MulticenterSet::single(2, 1)]), 4, 3, vec![0, 1, 2])]
    fn test_multicenter_bond(
        #[case] bond: MulticenterBond,
        #[case] expected_electron_count: u8,
        #[case] expected_atom_count: usize,
        #[case] expected_atoms: Vec<u32>,
    ) {
        assert_eq!(bond.electron_count(), expected_electron_count);
        assert_eq!(bond.atom_count(), expected_atom_count);
        assert_eq!(bond.all_atoms(), expected_atoms);
    }
}
