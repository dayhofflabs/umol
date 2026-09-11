use rstest::rstest;

use super::{derive_stereo_bonds, DirectionError};
use crate::smiles::{parse_extended_smiles_bytes, ParseError, Smiles};
use crate::table_ir::BondDirection::{Falling, Rising};
use crate::table_ir::BondOrder::{Double, Single};
use crate::table_ir::BondRelation::{OppositeSide, SameSide};
use crate::table_ir::{
    AtomPair, BondConfiguration, BondDirection, BondOrder, BondRelation, StereoBond,
};

#[rstest]
#[case::absent("FC=CF", vec![])]
#[case::first_partial("F/C=CF", vec![])]
#[case::first_partial_reversed("F\\C=CF", vec![])]
#[case::second_partial("FC=C/F", vec![])]
#[case::opposite("F/C=C/F", vec![(1, [0, 3], OppositeSide)])]
#[case::reversed("F\\C=C\\F", vec![(1, [0, 3], OppositeSide)])]
#[case::same("F/C=C\\F", vec![(1, [0, 3], SameSide)])]
#[case::same_reversed("F\\C=C/F", vec![(1, [0, 3], SameSide)])]
#[case::branch("C(/F)=C/F", vec![(1, [1, 3], SameSide)])]
#[case::four_substituents("F/C(Cl)=C(/Br)I", vec![(2, [0, 4], OppositeSide)])]
#[case::redundant("F/C(/Cl)=C(/Br)\\I", vec![(2, [0, 4], OppositeSide)])]
#[case::mixed_glyphs("C(/F)(\\Cl)=C(/Br)\\I", vec![(2, [1, 4], SameSide)])]
#[case::complement("FC(/Cl)=C/Br", vec![(2, [0, 4], OppositeSide)])]
#[case::complement_partial("FC(/Cl)=CBr", vec![])]
#[case::redundant_partial("F/C(/Cl)=CBr", vec![])]
#[case::shared("C/C=C/C=C/C", vec![(1, [0, 3], OppositeSide), (3, [2, 5], OppositeSide)])]
#[case::shared_opposite("C/C=C/C=C\\C", vec![(1, [0, 3], OppositeSide), (3, [2, 5], SameSide)])]
#[case::shared_partial("C/C=C/C=CC", vec![(1, [0, 3], OppositeSide)])]
#[case::shared_only("CC=C/C=CC", vec![])]
#[case::outer_only("C/C=CC=C/C", vec![])]
#[case::spacer("F/C=C/CC=CC", vec![(1, [0, 3], OppositeSide)])]
#[case::triene("C/C=C/C=C/C=C/C", vec![
    (1, [0, 3], OppositeSide), (3, [2, 5], OppositeSide), (5, [4, 7], OppositeSide),
])]
#[case::triene_partial("C/C=C/C=CC=C/C", vec![(1, [0, 3], OppositeSide)])]
#[case::explicit_h("C/C=C/C=CC(/[H])=C/C", vec![(1, [0, 3], OppositeSide), (6, [4, 8], OppositeSide)])]
#[case::explicit_h_reversed("C/C=C/C=CC(\\[H])=C/C", vec![(1, [0, 3], OppositeSide), (6, [4, 8], SameSide)])]
#[case::ring_open("C/C=C/1CO1", vec![(1, [0, 3], SameSide)])]
#[case::ring_close("C/C=C1CO\\1", vec![(1, [0, 3], SameSide)])]
#[case::ring_both("C/C=C/1CO\\1", vec![(1, [0, 3], SameSide)])]
#[case::shared_branch("F/C=C(/C=C/F)C=C", vec![(1, [0, 3], OppositeSide), (3, [2, 5], OppositeSide)])]
#[case::shared_cycle("C1/C=C/C=C/CCC1", vec![(2, [0, 3], OppositeSide), (4, [2, 5], OppositeSide)])]
#[case::tetrahedral_ring("C[C@H]1/C=C/CCO1", vec![(3, [1, 4], OppositeSide)])]
#[case::components("F/C=C/F.CC.F/C=C\\Cl", vec![(1, [0, 3], OppositeSide), (5, [6, 9], SameSide)])]
#[case::reverse_branch(r"C(\F)=C/F", vec![(1,[1,3],OppositeSide)])]
#[case::substituted_cis_1(r"C/C(/F)=C(\F)/C", vec![(2,[0,4],SameSide)])]
#[case::substituted_cis_2(r"C/C(/F)=C(/C)\F", vec![(2,[0,4],OppositeSide)])]
#[case::substituted_cis_3(r"C/C(F)=C(/C)F", vec![(2,[0,4],OppositeSide)])]
#[case::substituted_cis_4(r"CC(/F)=C(/C)F", vec![(2,[0,4],OppositeSide)])]
#[case::substituted_cis_5(r"C/C(F)=C(C)\F", vec![(2,[0,4],OppositeSide)])]
#[case::substituted_cis_6(r"CC(/F)=C(C)\F", vec![(2,[0,4],OppositeSide)])]
#[case::separate_double_bonds(r"F/C=C/C/C=C\C", vec![(1,[0,3],OppositeSide),(4,[3,6],SameSide)])]
fn test_parse_molecule_stereo_bonds(
    #[case] input: &str,
    #[case] expected: Vec<(u32, [u32; 2], BondRelation)>,
    #[values(false, true)] extended: bool,
) {
    let frames = if extended {
        parse_extended_smiles_bytes(input.as_bytes())
            .unwrap()
            .stereo_bonds
    } else {
        Smiles::parse(input).unwrap().into_table_ir().stereo_bonds
    };
    assert_eq!(
        frames,
        expected
            .into_iter()
            .map(|(bond, references, relation)| StereoBond {
                bond,
                configuration: BondConfiguration::Framed {
                    references,
                    relation
                },
            })
            .collect::<Vec<_>>()
    );
}

#[rstest]
#[case::complete("F/C(\\Cl)=C/Br", DirectionError::CisTransConflict { atom: 1 })]
#[case::partial("F/C(\\Cl)=CBr", DirectionError::CisTransConflict { atom: 1 })]
#[case::second_partial("FC=C(/Br)/I", DirectionError::CisTransConflict { atom: 2 })]
#[case::dangling("C/C", DirectionError::DanglingBondDirection { bond: 0 })]
#[case::missing_reference("F/C=C", DirectionError::DanglingBondDirection { bond: 0 })]
#[case::excess_ligands("F/C(Cl)(Br)=C/I", DirectionError::UnsupportedSite { bond: 3 })]
#[case::shared_ligand("C1/C=C1", DirectionError::UnsupportedSite { bond: 2 })]
#[case::odd_cumulene("F/C=C=C=C/F", DirectionError::UnsupportedSite { bond: 1 })]
#[case::even_cumulene("F/C=C=C/F", DirectionError::UnsupportedSite { bond: 1 })]
#[case::partial_cumulene("FC=C=C=C/F", DirectionError::UnsupportedSite { bond: 3 })]
#[case::long_cumulene(r"F/C=C=C=C=C=C/F", DirectionError::UnsupportedSite { bond:1 })]
#[case::long_cumulene_opposite(r"F/C=C=C=C=C=C\F", DirectionError::UnsupportedSite { bond:1 })]
#[case::cumulene_opposite(r"F/C=C=C=C\F", DirectionError::UnsupportedSite { bond:1 })]
#[case::shared_cis_trans_ligand(r"SSC=S1CC1\2C=112", DirectionError::UnsupportedSite { bond:8 })]
fn test_parse_molecule_stereo_bonds_error(
    #[case] input: &str,
    #[case] expected: DirectionError,
    #[values(false, true)] extended: bool,
) {
    let result = if extended {
        parse_extended_smiles_bytes(input.as_bytes()).map(|_| ())
    } else {
        Smiles::parse(input).map(|_| ())
    };
    assert_eq!(result, Err(ParseError::from(expected)));
}

#[rstest]
#[case::opening_order(vec![
    (AtomPair::new(0, 1), Single, Some(Rising)),
    (AtomPair::new(1, 2), Double, None),
    (AtomPair::new(2, 4), Single, Some(Rising)),
    (AtomPair::new(2, 3), Single, None),
    (AtomPair::new(3, 4), Single, None),
], 1)]
#[case::completion_order(vec![
    (AtomPair::new(0, 1), Single, Some(Rising)),
    (AtomPair::new(1, 2), Double, None),
    (AtomPair::new(2, 3), Single, None),
    (AtomPair::new(3, 4), Single, None),
    (AtomPair::new(2, 4), Single, Some(Rising)),
], 1)]
#[case::reversed_table(vec![
    (AtomPair::new(3, 4), Single, None),
    (AtomPair::new(2, 3), Single, None),
    (AtomPair::new(2, 4), Single, Some(Rising)),
    (AtomPair::new(1, 2), Double, None),
    (AtomPair::new(0, 1), Single, Some(Rising)),
], 3)]
fn test_derive_stereo_bonds_order(
    #[case] bonds: Vec<(AtomPair, BondOrder, Option<BondDirection>)>,
    #[case] bond: u32,
) {
    assert_eq!(
        derive_stereo_bonds(5, &bonds),
        Ok(vec![StereoBond {
            bond,
            configuration: BondConfiguration::Framed {
                references: [0, 3],
                relation: SameSide
            },
        }])
    );
}

#[rstest]
#[case::first_endpoint(2, vec![(AtomPair::new(2, 3), Single, Some(Rising))], DirectionError::AtomIndexOutOfBounds { atom: 2 })]
#[case::second_endpoint(2, vec![(AtomPair::new(0, 3), Single, Some(Rising))], DirectionError::AtomIndexOutOfBounds { atom: 3 })]
#[case::site_as_reference(4, vec![
    (AtomPair::new(1, 2), Double, None),
    (AtomPair::new(1, 1), Single, Some(Rising)),
    (AtomPair::new(2, 3), Single, None),
], DirectionError::UnsupportedSite { bond: 0 })]
#[case::parallel_conflict(4, vec![
    (AtomPair::new(0, 1), Single, Some(Rising)),
    (AtomPair::new(0, 1), Single, Some(Falling)),
    (AtomPair::new(1, 2), Double, None),
    (AtomPair::new(2, 3), Single, None),
], DirectionError::CisTransConflict { atom: 1 })]
fn test_derive_stereo_bonds_incidence_error(
    #[case] atom_count: usize,
    #[case] bonds: Vec<(AtomPair, BondOrder, Option<BondDirection>)>,
    #[case] expected: DirectionError,
) {
    assert_eq!(derive_stereo_bonds(atom_count, &bonds), Err(expected));
}

#[rstest]
#[case::both_ends("C/C=C/1CO/1", ParseError::MismatchedRingBondDirections { pos: 10, open_pos: 6 })]
fn test_derive_stereo_bonds_ring_error(#[case] input: &str, #[case] expected: ParseError) {
    assert_eq!(Smiles::parse(input), Err(expected.clone()));
    assert_eq!(parse_extended_smiles_bytes(input.as_bytes()), Err(expected));
}

#[rstest]
fn test_derive_stereo_bonds_exhaustive() {
    // Enumerate physical half-plane assignments, independently of marker propagation.
    let arrangements = [
        [true, false, true, false],
        [true, false, false, true],
        [false, true, true, false],
        [false, true, false, true],
    ];
    for code in 0..81 {
        let mut digits = code;
        let markers = [0; 4].map(|_| {
            let marker = [None, Some(Rising), Some(Falling)][digits % 3];
            digits /= 3;
            marker
        });
        let bonds = [
            (AtomPair::new(1, 4), Double, None),
            (AtomPair::new(0, 1), Single, markers[0]),
            (AtomPair::new(1, 2), Single, markers[1]),
            (AtomPair::new(3, 4), Single, markers[2]),
            (AtomPair::new(4, 5), Single, markers[3]),
        ];
        let compatible = |arrangement: [bool; 4], indices: &[usize]| {
            let rising = [
                !arrangement[0],
                arrangement[1],
                !arrangement[2],
                arrangement[3],
            ];
            indices
                .iter()
                .all(|&i| markers[i].is_none_or(|marker| (marker == Rising) == rising[i]))
        };
        let possible: Vec<_> = arrangements
            .into_iter()
            .filter(|&a| compatible(a, &[0, 1, 2, 3]))
            .collect();
        let expected = if possible.is_empty() {
            let atom = if arrangements.into_iter().any(|a| compatible(a, &[0, 1])) {
                4
            } else {
                1
            };
            Err(DirectionError::CisTransConflict { atom })
        } else if possible
            .iter()
            .all(|a| (a[0] == a[2]) == (possible[0][0] == possible[0][2]))
        {
            Ok(vec![StereoBond {
                bond: 0,
                configuration: BondConfiguration::Framed {
                    references: [0, 3],
                    relation: if possible[0][0] == possible[0][2] {
                        SameSide
                    } else {
                        OppositeSide
                    },
                },
            }])
        } else {
            Ok(vec![])
        };
        assert_eq!(
            derive_stereo_bonds(6, &bonds),
            expected,
            "marker assignment {code}"
        );
    }
}
