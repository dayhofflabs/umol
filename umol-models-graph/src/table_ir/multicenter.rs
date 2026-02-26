//! Multi-center bond representation for table-based molecular models.
//!
//! Multi-center bonds describe bonding interactions involving more than two atoms,
//! where electrons are shared or donated across multiple centers simultaneously.
//! An ordinary bond of order `n` is equivalent to a multi-center bond with two
//! contributions, each consisting of a single atom and `n` electrons.
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

use std::collections::HashMap;

use smallvec::SmallVec;

/// A set of atoms participating in a multi-center bond.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterSet {
    atoms: SmallVec<[u32; 8]>,
}

impl MulticenterSet {
    pub fn new(atoms: Vec<u32>) -> Self {
        Self::from_atoms(atoms)
    }

    pub fn from_atoms<I>(atoms: I) -> Self
    where
        I: IntoIterator<Item = u32>,
    {
        let mut atoms: SmallVec<[u32; 8]> = atoms.into_iter().collect();
        atoms.sort_unstable();
        atoms.dedup();
        Self { atoms }
    }

    pub fn single(atom: u32) -> Self {
        Self {
            atoms: SmallVec::from_elem(atom, 1),
        }
    }

    pub fn atoms(&self) -> &[u32] {
        self.atoms.as_slice()
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn set_atoms<I>(&mut self, atoms: I)
    where
        I: IntoIterator<Item = u32>,
    {
        self.atoms = atoms.into_iter().collect();
        self.atoms.sort_unstable();
        self.atoms.dedup();
    }

    pub fn insert_atom(&mut self, atom: u32) {
        match self.atoms.binary_search(&atom) {
            Ok(_) => {}
            Err(pos) => self.atoms.insert(pos, atom),
        }
    }

    pub fn remove_atom(&mut self, atom: u32) -> bool {
        match self.atoms.binary_search(&atom) {
            Ok(pos) => {
                self.atoms.remove(pos);
                true
            }
            Err(_) => false,
        }
    }
}

/// A multi-center bond involving three or more atoms.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterBond {
    contributions: SmallVec<[MulticenterSet; 4]>,
}

impl MulticenterBond {
    pub fn new(contributions: Vec<MulticenterSet>) -> Self {
        Self::from_contributions(contributions)
    }

    pub fn from_contributions<I>(contributions: I) -> Self
    where
        I: IntoIterator<Item = MulticenterSet>,
    {
        let mut contributions: SmallVec<[MulticenterSet; 4]> = contributions.into_iter().collect();
        contributions.sort_unstable();
        contributions.dedup();
        Self { contributions }
    }

    pub fn contributions(&self) -> &[MulticenterSet] {
        self.contributions.as_slice()
    }

    pub fn set_contributions<I>(&mut self, contributions: I)
    where
        I: IntoIterator<Item = MulticenterSet>,
    {
        self.contributions = contributions.into_iter().collect();
        self.contributions.sort_unstable();
        self.contributions.dedup();
    }

    pub fn insert_contribution(&mut self, contribution: MulticenterSet) {
        match self.contributions.binary_search(&contribution) {
            Ok(_) => {}
            Err(pos) => self.contributions.insert(pos, contribution),
        }
    }

    pub fn remove_contribution(&mut self, contribution: &MulticenterSet) -> bool {
        match self.contributions.binary_search(contribution) {
            Ok(pos) => {
                self.contributions.remove(pos);
                true
            }
            Err(_) => false,
        }
    }

    pub fn atom_count(&self) -> usize {
        self.contributions
            .iter()
            .map(MulticenterSet::atom_count)
            .sum()
    }

    pub fn all_atoms(&self) -> Vec<u32> {
        self.contributions
            .iter()
            .flat_map(|c| c.atoms().iter().copied())
            .collect()
    }

    /// Return a multicenter bond with remapped atom indices.
    /// Returns `None` when at least one atom has no mapping.
    pub fn update_atoms(&self, index_map: &HashMap<u32, u32>) -> Option<Self> {
        let contributions = self
            .contributions
            .iter()
            .map(|contribution| {
                contribution
                    .atoms()
                    .iter()
                    .map(|old_idx| index_map.get(old_idx).copied())
                    .collect::<Option<Vec<u32>>>()
                    .map(|atoms| MulticenterSet::new(atoms))
            })
            .collect::<Option<Vec<MulticenterSet>>>()?;

        Some(Self::new(contributions))
    }
}

#[cfg(test)]
mod tests {
    use map_macro::hash_map;
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::ferrocene(MulticenterSet::new(vec![0, 1, 2, 3, 4]), 5, vec![0, 1, 2, 3, 4])]
    #[case::single(MulticenterSet::single(0), 1, vec![0])]
    fn test_multicenter_set_new(
        #[case] set: MulticenterSet,
        #[case] expected_atom_count: usize,
        #[case] expected_atoms: Vec<u32>,
    ) {
        assert_eq!(set.atom_count(), expected_atom_count);
        assert_eq!(set.atoms.into_iter().collect::<Vec<u32>>(), expected_atoms);
    }

    #[rstest]
    #[case::sorted_unique(vec![3, 1, 2, 2], vec![1, 2, 3])]
    #[case::single(vec![7], vec![7])]
    fn test_multicenter_set_set_atoms(#[case] input: Vec<u32>, #[case] expected: Vec<u32>) {
        let mut set = MulticenterSet::single(0);
        set.set_atoms(input);
        assert_eq!(set.atoms.as_slice(), expected.as_slice());
    }

    #[rstest]
    #[case::insert_new(vec![1, 3], 2, vec![1, 2, 3])]
    #[case::insert_existing(vec![1, 2, 3], 2, vec![1, 2, 3])]
    fn test_multicenter_set_insert_atom(
        #[case] atoms: Vec<u32>,
        #[case] atom: u32,
        #[case] expected: Vec<u32>,
    ) {
        let mut set = MulticenterSet::new(atoms);
        set.insert_atom(atom);
        assert_eq!(set.atoms.as_slice(), expected.as_slice());
    }

    #[rstest]
    #[case::remove_present(vec![1, 2, 3], 2, true, vec![1, 3])]
    #[case::remove_missing(vec![1, 2, 3], 4, false, vec![1, 2, 3])]
    fn test_multicenter_set_remove_atom(
        #[case] atoms: Vec<u32>,
        #[case] atom: u32,
        #[case] expected_removed: bool,
        #[case] expected_atoms: Vec<u32>,
    ) {
        let mut set = MulticenterSet::new(atoms);
        let removed = set.remove_atom(atom);
        assert_eq!(removed, expected_removed);
        assert_eq!(set.atoms.as_slice(), expected_atoms.as_slice());
    }

    #[rstest]
    #[case::ferrocene(MulticenterBond::new(vec![MulticenterSet::new(vec![0, 1, 2, 3, 4]), MulticenterSet::single(9)]), 6, vec![0, 1, 2, 3, 4, 9])]
    #[case::diborane(MulticenterBond::new(vec![MulticenterSet::single(0), MulticenterSet::single(1), MulticenterSet::single(5)]), 3, vec![0, 1, 5])]
    #[case::xef2(MulticenterBond::new(vec![MulticenterSet::single(0), MulticenterSet::single(1), MulticenterSet::single(2)]), 3, vec![0, 1, 2])]
    #[case::one_set(MulticenterBond::new(vec![MulticenterSet::new(vec![0, 1, 2])]), 3, vec![0, 1, 2])]
    #[case::permuted(MulticenterBond::new(vec![MulticenterSet::new(vec![2, 1, 0])]), 3, vec![0, 1, 2])]
    #[case::non_contiguous(MulticenterBond::new(vec![MulticenterSet::new(vec![0, 2, 5])]), 3, vec![0, 2, 5])]
    #[case::duplicate_atoms(MulticenterBond::new(vec![MulticenterSet::new(vec![0, 1, 1, 2])]), 3, vec![0, 1, 2])]
    #[case::multiple_sets(MulticenterBond::new(vec![MulticenterSet::new(vec![4, 5, 6]),
           MulticenterSet::new(vec![0, 1, 2])]), 6, vec![0, 1, 2, 4, 5, 6])]
    #[case::multiple_sets_overlapping(MulticenterBond::new(vec![MulticenterSet::new(vec![0, 1, 2]),
           MulticenterSet::new(vec![1, 2, 3])]), 6, vec![0, 1, 2, 1, 2, 3])]
    fn test_multicenter_bond_new(
        #[case] bond: MulticenterBond,
        #[case] expected_atom_count: usize,
        #[case] expected_atoms: Vec<u32>,
    ) {
        assert_eq!(bond.atom_count(), expected_atom_count);
        assert_eq!(bond.all_atoms(), expected_atoms);
    }

    #[rstest]
    #[case::unsorted_duplicate(vec![MulticenterSet::new(vec![5, 4]), MulticenterSet::new(vec![1, 0]), MulticenterSet::new(vec![5, 4])],
        vec![MulticenterSet::new(vec![0, 1]), MulticenterSet::new(vec![4, 5])],)]
    fn test_multicenter_bond_set_contributions(
        #[case] input: Vec<MulticenterSet>,
        #[case] expected: Vec<MulticenterSet>,
    ) {
        let mut bond = MulticenterBond::new(vec![MulticenterSet::single(9)]);
        bond.set_contributions(input);
        assert_eq!(bond.contributions.as_slice(), expected.as_slice());
    }

    #[rstest]
    #[case::new(vec![MulticenterSet::new(vec![0, 1])], MulticenterSet::new(vec![2, 3]),
        vec![MulticenterSet::new(vec![0, 1]), MulticenterSet::new(vec![2, 3])])]
    #[case::existing(vec![MulticenterSet::new(vec![0, 1]), MulticenterSet::new(vec![2, 3])],
        MulticenterSet::new(vec![3, 2]), vec![MulticenterSet::new(vec![0, 1]), MulticenterSet::new(vec![2, 3])])]
    fn test_multicenter_bond_insert_contribution(
        #[case] initial: Vec<MulticenterSet>,
        #[case] contribution: MulticenterSet,
        #[case] expected: Vec<MulticenterSet>,
    ) {
        let mut bond = MulticenterBond::new(initial);
        bond.insert_contribution(contribution);
        assert_eq!(bond.contributions.as_slice(), expected.as_slice());
    }

    #[rstest]
    #[case::remove_present(vec![ MulticenterSet::new(vec![0, 1]), MulticenterSet::new(vec![2, 3])],
        MulticenterSet::new(vec![3, 2]), true, vec![MulticenterSet::new(vec![0, 1])])]
    #[case::remove_missing(vec![MulticenterSet::new(vec![0, 1])], MulticenterSet::new(vec![2, 3]),
        false, vec![MulticenterSet::new(vec![0, 1])])]
    fn test_multicenter_bond_remove_contribution(
        #[case] initial: Vec<MulticenterSet>,
        #[case] contribution: MulticenterSet,
        #[case] expected_removed: bool,
        #[case] expected: Vec<MulticenterSet>,
    ) {
        let mut bond = MulticenterBond::new(initial);
        let removed = bond.remove_contribution(&contribution);
        assert_eq!(removed, expected_removed);
        assert_eq!(bond.contributions.as_slice(), expected.as_slice());
    }

    #[rstest]
    #[case::complete_mapping(MulticenterBond::new(vec![ MulticenterSet::new(vec![2, 3]), MulticenterSet::single(5)]),
        hash_map! { 2u32 => 0u32, 3u32 => 1u32, 5u32 => 2u32 },
        Some(MulticenterBond::new(vec![MulticenterSet::new(vec![0, 1]), MulticenterSet::single(2)])))]
    #[case::missing_mapping(MulticenterBond::new(vec![MulticenterSet::new(vec![2, 3])]),
        hash_map! { 2u32 => 0u32 }, None)]
    fn test_multicenter_bond_update_atoms(
        #[case] bond: MulticenterBond,
        #[case] map: HashMap<u32, u32>,
        #[case] expected: Option<MulticenterBond>,
    ) {
        assert_eq!(bond.update_atoms(&map), expected);
    }
}
