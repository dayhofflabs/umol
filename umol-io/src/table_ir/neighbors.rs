//! Neighbor tables of TableIR records.

use super::bond::AtomPair;

/// One bond at an atom: the atom at its other end and the bond's index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Neighbor {
    pub atom: u32,
    pub bond: u32,
}

/// Bonds at each atom of a TableIR record, in bond order, with the atom at the other end.
/// Computed from the bond list on request; the record's fields stay authoritative, so a value is
/// held only across one operation on an unchanged record. Every bond is listed, duplicates
/// included; a bond joining an atom to itself is listed once at that atom; a bond naming an atom
/// outside the record contributes no entry, and the raise's integrity gate rejects it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomNeighbors {
    neighbors: Vec<Vec<Neighbor>>,
}

impl AtomNeighbors {
    /// Table for `atom_count` atoms joined by `bonds`, taken in bond order.
    pub fn new(atom_count: usize, bonds: impl IntoIterator<Item = AtomPair>) -> Self {
        let mut neighbors = vec![Vec::new(); atom_count];
        for (bond, pair) in (0u32..).zip(bonds) {
            let (first, second) = pair.as_tuple();
            if let Some(entries) = neighbors.get_mut(first as usize) {
                if (second as usize) < atom_count {
                    entries.push(Neighbor { atom: second, bond });
                }
            }
            if first != second {
                if let Some(entries) = neighbors.get_mut(second as usize) {
                    if (first as usize) < atom_count {
                        entries.push(Neighbor { atom: first, bond });
                    }
                }
            }
        }
        Self { neighbors }
    }

    /// Neighbors of `atom` in bond order.
    pub fn neighbors(&self, atom: u32) -> &[Neighbor] {
        &self.neighbors[atom as usize]
    }

    /// Number of distinct atoms bonded to `atom`.
    pub fn degree(&self, atom: u32) -> usize {
        let mut atoms: Vec<u32> = self.neighbors(atom).iter().map(|n| n.atom).collect();
        atoms.sort_unstable();
        atoms.dedup();
        atoms.len()
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::chain(3, vec![AtomPair::new(0, 1), AtomPair::new(1, 2)], vec![
        vec![Neighbor { atom: 1, bond: 0 }],
        vec![Neighbor { atom: 0, bond: 0 }, Neighbor { atom: 2, bond: 1 }],
        vec![Neighbor { atom: 1, bond: 1 }],
    ])]
    #[case::ring_closure_last(3, vec![AtomPair::new(0, 1), AtomPair::new(1, 2), AtomPair::new(0, 2)], vec![
        vec![Neighbor { atom: 1, bond: 0 }, Neighbor { atom: 2, bond: 2 }],
        vec![Neighbor { atom: 0, bond: 0 }, Neighbor { atom: 2, bond: 1 }],
        vec![Neighbor { atom: 1, bond: 1 }, Neighbor { atom: 0, bond: 2 }],
    ])]
    #[case::duplicate_bond(2, vec![AtomPair::new(0, 1), AtomPair::new(1, 0)], vec![
        vec![Neighbor { atom: 1, bond: 0 }, Neighbor { atom: 1, bond: 1 }],
        vec![Neighbor { atom: 0, bond: 0 }, Neighbor { atom: 0, bond: 1 }],
    ])]
    #[case::self_bond(1, vec![AtomPair::new(0, 0)], vec![vec![Neighbor { atom: 0, bond: 0 }]])]
    #[case::endpoint_outside_record(2, vec![AtomPair::new(0, 1), AtomPair::new(1, 5)], vec![
        vec![Neighbor { atom: 1, bond: 0 }],
        vec![Neighbor { atom: 0, bond: 0 }],
    ])]
    #[case::no_bonds(2, vec![], vec![vec![], vec![]])]
    fn test_atom_neighbors_new(
        #[case] atom_count: usize,
        #[case] bonds: Vec<AtomPair>,
        #[case] expected: Vec<Vec<Neighbor>>,
    ) {
        let table = AtomNeighbors::new(atom_count, bonds);
        let actual: Vec<Vec<Neighbor>> = (0..atom_count as u32)
            .map(|atom| table.neighbors(atom).to_vec())
            .collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::two(AtomNeighbors::new(3, [AtomPair::new(0, 1), AtomPair::new(1, 2)]), 1, 2)]
    #[case::duplicate_bond(AtomNeighbors::new(2, [AtomPair::new(0, 1), AtomPair::new(1, 0)]), 0, 1)]
    #[case::isolated(AtomNeighbors::new(2, [AtomPair::new(0, 1)]), 1, 1)]
    #[case::no_bonds(AtomNeighbors::new(1, []), 0, 0)]
    fn test_atom_neighbors_degree(
        #[case] table: AtomNeighbors,
        #[case] atom: u32,
        #[case] expected: usize,
    ) {
        assert_eq!(table.degree(atom), expected);
    }
}
