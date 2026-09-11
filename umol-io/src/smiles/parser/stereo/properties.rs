use std::cell::OnceCell;

use proptest::prelude::*;

use super::{derive_stereo_bonds, DirectionMarker};
use crate::table_ir::BondDirection::{Falling, Rising};
use crate::table_ir::BondOrder::{Double, Single};
use crate::table_ir::BondRelation::{OppositeSide, SameSide};
use crate::table_ir::{AtomPair, BondConfiguration, StereoBond};

proptest! {
    // Renumbering and global marker reversal preserve physical sides; the producer selects
    // references afresh in the target numbering. Partial evidence remains unasserted.
    #[test]
    fn test_derive_stereo_bonds_remapping(
        atom_keys in any::<[u16; 6]>(),
        bond_keys in any::<[u16; 5]>(),
        same_side in any::<bool>(),
        marked in any::<[bool; 4]>(),
        reverse in any::<bool>(),
    ) {
        let mut atom_order = [0, 1, 2, 3, 4, 5];
        atom_order.sort_by_key(|&atom| atom_keys[atom]);
        let mut image = [0; 6];
        for (new, old) in atom_order.into_iter().enumerate() {
            image[old] = new as u32;
        }
        let markers = [Falling, Falling,
            if same_side { Falling } else { Rising },
            if same_side { Falling } else { Rising }];
        let markers = markers.map(|marker| if reverse { marker.flip() } else { marker });
        let mut bonds = [
            (1, 4, Double, None),
            (0, 1, Single, marked[0].then_some(markers[0])),
            (1, 2, Single, marked[1].then_some(markers[1])),
            (3, 4, Single, marked[2].then_some(markers[2])),
            (4, 5, Single, marked[3].then_some(markers[3])),
        ].map(|(a, b, order, marker)| {
            let (a, b) = (image[a], image[b]);
            (AtomPair::new(a, b), order, if a > b { marker.map(|m| m.flip()) } else { marker })
        }).into_iter().enumerate().collect::<Vec<_>>();
        bonds.sort_by_key(|&(old, _)| bond_keys[old]);
        let site = bonds.iter().position(|&(old, _)| old == 0).unwrap() as u32;
        let bonds = bonds.into_iter().map(|(_, (atoms, order, marker))| (atoms, order, marker.map(DirectionMarker::new))).collect::<Vec<_>>();

        let expected = if (marked[0] || marked[1]) && (marked[2] || marked[3]) {
            let first = if image[0] < image[2] { 0 } else { 2 };
            let second = if image[3] < image[5] { 3 } else { 5 };
            let mut references = [image[first], image[second]];
            if image[1] > image[4] {
                references.swap(0, 1);
            }
            let first_above = first == 0;
            let second_above = if second == 3 { same_side } else { !same_side };
            vec![StereoBond { bond: site, configuration: BondConfiguration::Framed {
                references, relation: if first_above == second_above { SameSide } else { OppositeSide },
            }}]
        } else {
            vec![]
        };
        prop_assert_eq!(derive_stereo_bonds(6, &bonds, |bond| (bond.0, bond.1, bond.2.as_ref()), &OnceCell::new()), Ok(expected));
    }
}
