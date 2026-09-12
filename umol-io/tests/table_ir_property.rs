//! Bond incidences and degrees compared with the supplied bond rows.

use std::collections::BTreeSet;

use proptest::prelude::*;

use umol_io::table_ir::{AtomNeighbors, AtomPair, Neighbor};

proptest! {
    // Incident-list semantics against a direct scan, including duplicate and invalid pairs.
    #[test]
    fn test_atom_neighbors_new_incidence(
        atom_count in 0usize..24,
        pairs in prop::collection::vec((0u32..28, 0u32..28), 0..96),
    ) {
        let table = AtomNeighbors::new(atom_count, pairs.iter().map(|&(a,b)| AtomPair::new(a,b)));
        for atom in 0..atom_count as u32 {
            let mut expected: Vec<_> = pairs.iter().enumerate().filter_map(|(bond, &(a,b))| {
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
            expected.sort_unstable_by_key(|neighbor| (neighbor.atom, neighbor.bond));
            let degree = expected.iter().map(|neighbor| neighbor.atom).collect::<BTreeSet<_>>().len();
            prop_assert_eq!(table.neighbors(atom), expected.as_slice());
            prop_assert_eq!(table.degree(atom), degree);
        }
    }
}
