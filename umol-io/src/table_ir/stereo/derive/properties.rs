use std::array;

use proptest::prelude::*;
use umol_geometric_core::Point3D;

use super::derive_stereo_bonds;
use crate::table_ir::{
    AtomPair, BondConfiguration, BondOrder, BondRelation, BondStereo, StereoBond,
};

proptest! {
    // Signed coordinate permutations, translations, and positive scaling preserve relative sides.
    #[test]
    fn test_derive_stereo_bonds_similarity(
        offset in any::<[i16;3]>(),
        signs in any::<[bool;3]>(),
        cycle in 0_usize..3,
        scale in 1_u16..1000,
        same in any::<bool>(),
        reverse_atoms in any::<bool>(),
        reverse_bonds in any::<bool>(),
    ) {
        let points = [[0.,1.,0.],[0.,0.,0.],[2.,0.,0.],[2.,if same {1.} else {-1.},0.]];
        let mut points = points.map(|p| {
            let q: [f64;3] = array::from_fn(|i| {
                f64::from(offset[i]) + f64::from(scale) * p[(i+cycle)%3] * if signs[i] {1.} else {-1.}
            });
            Point3D::new(q[0],q[1],q[2])
        });
        let mut bonds = [(0,1,BondOrder::Single),(1,2,BondOrder::Double),(2,3,BondOrder::Single)]
            .map(|(a,b,o)| (AtomPair::new(if reverse_atoms {3-a} else {a}, if reverse_atoms {3-b} else {b}),o,None));
        if reverse_atoms { points.reverse(); }
        if reverse_bonds { bonds.reverse(); }
        let references = [0,3];
        prop_assert_eq!(derive_stereo_bonds(4, &bonds, |bond| *bond, Some(&points), Vec::new(), Vec::new()), Ok(vec![StereoBond {
            bond:1, configuration:BondConfiguration::Framed {
                references, relation:if same { BondRelation::SameSide } else { BondRelation::OppositeSide },
            },
        }]));
    }
    // Redundant source assertions preserve an existing frame in either choice of references.
    #[test]
    fn test_derive_stereo_bonds_frame_preservation(
        same in any::<bool>(),
        complements in any::<[bool;2]>(),
        copies in 0_usize..32,
    ) {
        let bonds = [
            (AtomPair::new(1,4), BondOrder::Double, None),
            (AtomPair::new(0,1), BondOrder::Single, None), (AtomPair::new(1,2), BondOrder::Single, None),
            (AtomPair::new(3,4), BondOrder::Single, None), (AtomPair::new(4,5), BondOrder::Single, None),
        ];
        let first_above = !complements[0];
        let second_above = if same { !complements[1] } else { complements[1] };
        let frames = vec![StereoBond { bond:0, configuration:BondConfiguration::Framed {
            references:[if complements[0] {2} else {0}, if complements[1] {5} else {3}],
            relation:if first_above == second_above {BondRelation::SameSide} else {BondRelation::OppositeSide},
        }}];
        let expected = frames.clone();
        let code = if same {BondStereo::Cis} else {BondStereo::Trans};
        prop_assert_eq!(derive_stereo_bonds(6, &bonds, |bond| *bond, None, frames, vec![(0,code);copies]), Ok(expected));
    }

}
