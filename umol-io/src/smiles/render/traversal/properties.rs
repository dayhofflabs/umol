//! Compare preorder with recursive bond-table scanning, and reconstruct complete edge incidence.
//! Generated tables include disconnected multigraphs, loops, and out-of-range endpoints: topology
//! traversal is total on these, while format acceptance remains a separate boundary concern.
//! Ring-label tests use a direct set definition of availability, including delayed reuse.

use std::collections::{BTreeMap, BTreeSet};

use proptest::prelude::*;

use super::Traversal;
use crate::table_ir::{Atom, AtomPair, Bond, BondOrder, Molecule, Neighbor};

fn visit(
    atom: u32,
    parent: Option<Neighbor>,
    count: usize,
    bonds: &[Bond],
    output: &mut Vec<(u32, Option<Neighbor>, usize)>,
) {
    let position = output.len();
    output.push((atom, parent, 0));
    for (bond, entry) in bonds.iter().enumerate() {
        if let Some(other) = entry.atoms.other(atom) {
            if (other as usize) < count && !output.iter().any(|&(seen, _, _)| seen == other) {
                visit(
                    other,
                    Some(Neighbor {
                        atom,
                        bond: bond as u32,
                    }),
                    count,
                    bonds,
                    output,
                );
            }
        }
    }
    output[position].2 = output.len();
}

proptest! {
    #[test]
    fn test_traversal_new(
        count in 0usize..16,
        endpoints in prop::collection::vec((0u32..18, 0u32..18), 0..70),
    ) {
        let molecule = Molecule {
            atoms: vec![Atom::wildcard(); count],
            bonds: endpoints.into_iter().map(|(a,b)| Bond::new(a,b,BondOrder::Single)).collect(),
            ..Molecule::empty()
        };
        let mut expected = Vec::new();
        for atom in 0..count as u32 {
            if !expected.iter().any(|&(seen, _, _)| seen == atom) {
                visit(atom, None, count, &molecule.bonds, &mut expected);
            }
        }
        let traversal = Traversal::new(&molecule);
        let actual: Vec<_> = traversal.atoms.iter().map(|a| (a.atom, a.parent, a.subtree_end)).collect();
        prop_assert_eq!(&actual, &expected);

        let mut emitted: Vec<_> = traversal.bonds().map(|(bond, [a,b])| (bond, AtomPair::new(a,b))).collect();
        emitted.sort_unstable();
        let edges: Vec<_> = molecule.bonds.iter().enumerate()
            .filter(|(_, b)| (b.atoms.second() as usize) < count)
            .map(|(i, b)| (i as u32, b.atoms)).collect();
        prop_assert_eq!(emitted, edges);

        let mut active = BTreeMap::new();
        let mut cursor = 0;
        for (position, atom) in traversal.atoms.iter().enumerate() {
            prop_assert_eq!(atom.rings.start, cursor);
            cursor = atom.rings.end;
            let children: Vec<_> = traversal.children(position).map(|a| a.atom).collect();
            let expected_children: Vec<_> = expected.iter()
                .filter(|(_, parent, _)| parent.is_some_and(|p| p.atom == atom.atom))
                .map(|&(child, _, _)| child).collect();
            prop_assert_eq!(children, expected_children);

            let mut actual_neighbors: Vec<_> = traversal.neighbors(position).map(|n| (n.atom,n.bond)).collect();
            actual_neighbors.sort_unstable();
            let mut expected_neighbors = Vec::new();
            for (bond, entry) in molecule.bonds.iter().enumerate() {
                let pair = entry.atoms;
                if (pair.second() as usize) >= count { continue; }
                if pair.first() == atom.atom { expected_neighbors.push((pair.second(), bond as u32)); }
                if pair.second() == atom.atom { expected_neighbors.push((pair.first(), bond as u32)); }
            }
            expected_neighbors.sort_unstable();
            prop_assert_eq!(actual_neighbors, expected_neighbors);

            let mut closed = BTreeSet::new();
            let rings = &traversal.rings[atom.rings.clone()];
            prop_assert!(rings.windows(2).all(|pair| pair[0].neighbor.bond <= pair[1].neighbor.bond));
            for ring in rings {
                if ring.opening {
                    let smallest = (1usize..).find(|label| !active.contains_key(label) && !closed.contains(label)).unwrap();
                    prop_assert_eq!(ring.label, smallest);
                    prop_assert_eq!(active.insert(ring.label, ring.neighbor.bond), None);
                } else {
                    prop_assert_eq!(active.remove(&ring.label), Some(ring.neighbor.bond));
                    closed.insert(ring.label);
                }
            }
        }
        prop_assert_eq!(cursor, traversal.rings.len());
        prop_assert_eq!(active, BTreeMap::new());
    }
}
