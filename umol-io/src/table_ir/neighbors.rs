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
/// held only while atom count, bond-table indices, and endpoints stay unchanged.
///
/// # Semantic properties
///
/// Every bond is listed, duplicates included; a bond joining an atom to itself is listed once
/// at that atom. A bond naming an atom outside the record contributes no entry, and raise's
/// integrity gate rejects it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomNeighbors {
    offsets: Vec<usize>,
    neighbors: Vec<Neighbor>,
}

impl AtomNeighbors {
    /// Table for `atom_count` atoms joined by `bonds`, taken in bond order.
    pub fn new(atom_count: usize, bonds: impl IntoIterator<Item = AtomPair>) -> Self {
        let bonds: Vec<_> = bonds.into_iter().collect();
        let mut offsets = vec![0; atom_count + 1];
        for pair in &bonds {
            let (first, second) = pair.as_tuple();
            if first as usize >= atom_count || second as usize >= atom_count {
                continue;
            }
            offsets[first as usize] += 1;
            if first != second {
                offsets[second as usize] += 1;
            }
        }
        let mut total = 0;
        for offset in &mut offsets {
            total += *offset;
            *offset = total;
        }
        let mut neighbors = vec![Neighbor { atom: 0, bond: 0 }; total];
        for (bond, pair) in bonds.into_iter().enumerate().rev() {
            let (first, second) = pair.as_tuple();
            if first as usize >= atom_count || second as usize >= atom_count {
                continue;
            }
            offsets[first as usize] -= 1;
            neighbors[offsets[first as usize]] = Neighbor {
                atom: second,
                bond: bond as u32,
            };
            if first != second {
                offsets[second as usize] -= 1;
                neighbors[offsets[second as usize]] = Neighbor {
                    atom: first,
                    bond: bond as u32,
                };
            }
        }
        Self { offsets, neighbors }
    }

    /// Neighbors of `atom` in bond order.
    pub fn neighbors(&self, atom: u32) -> &[Neighbor] {
        &self.neighbors[self.offsets[atom as usize]..self.offsets[atom as usize + 1]]
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
    #[case::empty(0, vec![], vec![])]
    #[case::invalid_before_valid(3, vec![AtomPair::new(0,9), AtomPair::new(0,2), AtomPair::new(1,1)], vec![
        vec![Neighbor { atom: 2, bond: 1 }],
        vec![Neighbor { atom: 1, bond: 2 }],
        vec![Neighbor { atom: 0, bond: 1 }],
    ])]
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

#[cfg(all(test, feature = "proptest"))]
mod properties {
    use std::collections::BTreeSet;

    use proptest::prelude::*;

    use super::{AtomNeighbors, AtomPair, Neighbor};

    proptest! {
        // Incident-list semantics against a direct scan, including duplicate and invalid pairs.
        #[test]
        fn test_atom_neighbors_new_incidence(
            atom_count in 0usize..24,
            pairs in prop::collection::vec((0u32..28, 0u32..28), 0..96),
        ) {
            let table = AtomNeighbors::new(atom_count, pairs.iter().map(|&(a,b)| AtomPair::new(a,b)));
            for atom in 0..atom_count as u32 {
                let expected: Vec<_> = pairs.iter().enumerate().filter_map(|(bond, &(a,b))| {
                    if a as usize >= atom_count || b as usize >= atom_count {
                        None
                    } else if a == atom {
                        Some(Neighbor { atom: b, bond: bond as u32 })
                    } else if b == atom {
                        Some(Neighbor { atom: a, bond: bond as u32 })
                    } else {
                        None
                    }
                }).collect();
                let degree = expected.iter().map(|neighbor| neighbor.atom).collect::<BTreeSet<_>>().len();
                prop_assert_eq!(table.neighbors(atom), expected.as_slice());
                prop_assert_eq!(table.degree(atom), degree);
            }
        }
    }
}
