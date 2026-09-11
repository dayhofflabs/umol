//! Coordinate changes and row ordering preserve parsed cis/trans stereo.

use std::array;

use proptest::prelude::*;
use umol_io::ctfile::parser::parse_mol_to_table_ir;
use umol_io::table_ir::{BondConfiguration, BondRelation, StereoBond};

proptest! {
    #[test]
    fn test_parse_mol_to_table_ir_stereo(
        offset in prop::array::uniform3(-100_i16..100), signs in any::<[bool; 3]>(),
        cycle in 0_usize..3, scale in 1_u16..100, same in any::<bool>(),
        reverse_atoms in any::<bool>(), reverse_bonds in any::<bool>(),
    ) {
        let points = [[0.,1.,0.], [0.,0.,0.], [2.,0.,0.], [2.,if same {1.} else {-1.},0.]];
        let mut points: [[f64; 3]; 4] = points.map(|p| array::from_fn(|i|
            f64::from(offset[i]) + f64::from(scale) * p[(i + cycle) % 3]
                * if signs[i] {1.} else {-1.}));
        if reverse_atoms { points.reverse(); }
        let mut input = String::from("stereo\n  umol          3D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n");
        for (i, [x, y, z]) in points.into_iter().enumerate() {
            input.push_str(&format!("{x:10.4}{y:10.4}{z:10.4} {:3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
                if i == 0 || i == 3 { "F" } else { "C" }));
        }
        let mut bonds = [(1, 2, 1), (2, 3, 2), (3, 4, 1)];
        if reverse_bonds { bonds.reverse(); }
        for (a, b, order) in bonds {
            let (a, b) = if reverse_atoms { (5 - a, 5 - b) } else { (a, b) };
            input.push_str(&format!("{a:3}{b:3}{order:3}  0  0  0  0\n"));
        }
        input.push_str("M  END\n");
        let table = parse_mol_to_table_ir(&input).unwrap();
        prop_assert_eq!(table.stereo_bonds, vec![StereoBond { bond: 1,
            configuration: BondConfiguration::Framed { references: [0, 3],
                relation: if same { BondRelation::SameSide } else { BondRelation::OppositeSide },
            },
        }]);
    }
}
