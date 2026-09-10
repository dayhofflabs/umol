use rstest::rstest;

use super::super::*;
use crate::table_ir::{BondDirection, BondOrder, StereoAtom, StereoLigand, Winding};

#[rstest]
#[case::opening("C[C@H]1CCCCO1", vec![StereoAtom { atom: 1, ligands: vec![StereoLigand::Atom(0), StereoLigand::ImplicitHydrogen, StereoLigand::Atom(6), StereoLigand::Atom(2)], winding: Winding::CounterClockwise }])]
#[case::closing("O1CCCC[C@@H]1C", vec![StereoAtom { atom: 5, ligands: vec![StereoLigand::Atom(4), StereoLigand::ImplicitHydrogen, StereoLigand::Atom(0), StereoLigand::Atom(6)], winding: Winding::Clockwise }])]
#[case::closing_mirror("O1CCCC[C@H]1C", vec![StereoAtom { atom: 5, ligands: vec![StereoLigand::Atom(4), StereoLigand::ImplicitHydrogen, StereoLigand::Atom(0), StereoLigand::Atom(6)], winding: Winding::CounterClockwise }])]
#[case::explicit_hydrogen("C[C@]1([H])CCCCO1", vec![StereoAtom { atom: 1, ligands: vec![StereoLigand::Atom(0), StereoLigand::Atom(7), StereoLigand::Atom(2), StereoLigand::Atom(3)], winding: Winding::CounterClockwise }])]
#[case::acyclic("N[C@H](F)Cl", vec![StereoAtom { atom: 1, ligands: vec![StereoLigand::Atom(0), StereoLigand::ImplicitHydrogen, StereoLigand::Atom(2), StereoLigand::Atom(3)], winding: Winding::CounterClockwise }])]
#[case::root("[C@H](F)(Cl)Br", vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::ImplicitHydrogen, StereoLigand::Atom(1), StereoLigand::Atom(2), StereoLigand::Atom(3)], winding: Winding::CounterClockwise }])]
#[case::later_root("C.[C@H](F)(Cl)Br", vec![StereoAtom { atom: 1, ligands: vec![StereoLigand::ImplicitHydrogen, StereoLigand::Atom(2), StereoLigand::Atom(3), StereoLigand::Atom(4)], winding: Winding::CounterClockwise }])]
#[case::lone_pair("C[S@](=O)CC", vec![StereoAtom { atom: 1, ligands: vec![StereoLigand::Atom(0), StereoLigand::LonePair, StereoLigand::Atom(2), StereoLigand::Atom(3)], winding: Winding::CounterClockwise }])]
#[case::mixed_open_close("O1CCC[C@]21CCCC2", vec![StereoAtom { atom: 4, ligands: vec![StereoLigand::Atom(3), StereoLigand::Atom(8), StereoLigand::Atom(0), StereoLigand::Atom(5)], winding: Winding::CounterClockwise }])]
#[case::mixed_close_open("O1CCC[C@@]12CCCC2", vec![StereoAtom { atom: 4, ligands: vec![StereoLigand::Atom(3), StereoLigand::Atom(0), StereoLigand::Atom(8), StereoLigand::Atom(5)], winding: Winding::Clockwise }])]
#[case::two_closures("O1CCC2CC[C@]12F", vec![StereoAtom { atom: 6, ligands: vec![StereoLigand::Atom(5), StereoLigand::Atom(0), StereoLigand::Atom(3), StereoLigand::Atom(7)], winding: Winding::CounterClockwise }])]
#[case::branch_before_close("O1CCC[C@](F)1Cl", vec![StereoAtom { atom: 4, ligands: vec![StereoLigand::Atom(3), StereoLigand::Atom(5), StereoLigand::Atom(0), StereoLigand::Atom(6)], winding: Winding::CounterClockwise }])]
#[case::branch_return("[C@](CCC1)(F)(Cl)1", vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::Atom(1), StereoLigand::Atom(4), StereoLigand::Atom(5), StereoLigand::Atom(3)], winding: Winding::CounterClockwise }])]
#[case::large_label("O%99CCCC[C@@H]%99C", vec![StereoAtom { atom: 5, ligands: vec![StereoLigand::Atom(4), StereoLigand::ImplicitHydrogen, StereoLigand::Atom(0), StereoLigand::Atom(6)], winding: Winding::Clockwise }])]
#[case::label_reuse("C1CC1.O1CCCC[C@@H]1C", vec![StereoAtom { atom: 8, ligands: vec![StereoLigand::Atom(7), StereoLigand::ImplicitHydrogen, StereoLigand::Atom(3), StereoLigand::Atom(9)], winding: Winding::Clockwise }])]
#[case::incomplete("[C@]", vec![StereoAtom { atom: 0, ligands: vec![], winding: Winding::CounterClockwise }])]
#[case::repeated_hydrogen("[C@H2](F)Cl", vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::ImplicitHydrogen, StereoLigand::ImplicitHydrogen, StereoLigand::Atom(1), StereoLigand::Atom(2)], winding: Winding::CounterClockwise }])]
#[case::unsupported("[C@SP1](F)(Cl)(Br)I", vec![])]
#[case::unmarked_ring("C1CCCCC1", vec![])]
#[case::th1("N[C@TH1H](F)Cl", vec![StereoAtom { atom: 1, ligands: vec![StereoLigand::Atom(0), StereoLigand::ImplicitHydrogen, StereoLigand::Atom(2), StereoLigand::Atom(3)], winding: Winding::CounterClockwise }])]
#[case::th2("N[C@TH2H](F)Cl", vec![StereoAtom { atom: 1, ligands: vec![StereoLigand::Atom(0), StereoLigand::ImplicitHydrogen, StereoLigand::Atom(2), StereoLigand::Atom(3)], winding: Winding::Clockwise }])]
#[case::two_roots("[C@H](F)(Cl)Br.[C@@H](F)(Cl)Br", vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::ImplicitHydrogen, StereoLigand::Atom(1), StereoLigand::Atom(2), StereoLigand::Atom(3)], winding: Winding::CounterClockwise }, StereoAtom { atom: 4, ligands: vec![StereoLigand::ImplicitHydrogen, StereoLigand::Atom(5), StereoLigand::Atom(6), StereoLigand::Atom(7)], winding: Winding::Clockwise }])]
#[case::marked_both_ends("[C@]1(F)(Cl)CC[C@@]1(Br)I", vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::Atom(5), StereoLigand::Atom(1), StereoLigand::Atom(2), StereoLigand::Atom(3)], winding: Winding::CounterClockwise }, StereoAtom { atom: 5, ligands: vec![StereoLigand::Atom(4), StereoLigand::Atom(0), StereoLigand::Atom(6), StereoLigand::Atom(7)], winding: Winding::Clockwise }])]
#[case::two_openings("[C@]12(CCC1)CCC2", vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::Atom(3), StereoLigand::Atom(6), StereoLigand::Atom(1), StereoLigand::Atom(4)], winding: Winding::CounterClockwise }])]
#[case::allene("[C@AL1](F)(Cl)(Br)I", vec![])]
#[case::tb("[C@TB1](F)(Cl)(Br)I", vec![])]
#[case::oh("[C@OH1](F)(Cl)(Br)I", vec![])]
#[case::root_lone_pair("[S@](=O)(C)CC", vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::LonePair, StereoLigand::Atom(1), StereoLigand::Atom(2), StereoLigand::Atom(3)], winding: Winding::CounterClockwise }])]
#[case::self_loop("[C@]11", vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::Atom(0), StereoLigand::Atom(0)], winding: Winding::CounterClockwise }])]
#[case::parallel_bonds("C[C]2[C@@]2[C-]", vec![StereoAtom { atom: 2, ligands: vec![StereoLigand::Atom(1), StereoLigand::LonePair, StereoLigand::Atom(1), StereoLigand::Atom(3)], winding: Winding::Clockwise }])]
fn test_parse_molecule_stereo(#[case] input: &str, #[case] expected: Vec<StereoAtom>) {
    let basic = parse_molecule(input.as_bytes(), &SmilesIoConfig::default()).unwrap();
    let extended =
        parse_extended_smiles_bytes_with(input.as_bytes(), &SmilesIoConfig::default()).unwrap();
    assert_eq!(basic.stereo_atoms, expected);
    assert_eq!(extended.stereo_atoms, expected);
    assert_eq!(ExtendedMolecule::from(basic).stereo_atoms, expected);
    assert_eq!(Molecule::try_from(extended).unwrap().stereo_atoms, expected);
}

#[rstest]
fn test_parse_reaction_stereo() {
    let input = b"N[C@H](F)Cl>[C@@H](F)(Cl)Br>C.[C@H](F)(Cl)Br";
    let expected = [
        vec![StereoAtom {
            atom: 1,
            ligands: vec![
                StereoLigand::Atom(0),
                StereoLigand::ImplicitHydrogen,
                StereoLigand::Atom(2),
                StereoLigand::Atom(3),
            ],
            winding: Winding::CounterClockwise,
        }],
        vec![StereoAtom {
            atom: 0,
            ligands: vec![
                StereoLigand::ImplicitHydrogen,
                StereoLigand::Atom(1),
                StereoLigand::Atom(2),
                StereoLigand::Atom(3),
            ],
            winding: Winding::Clockwise,
        }],
        vec![StereoAtom {
            atom: 1,
            ligands: vec![
                StereoLigand::ImplicitHydrogen,
                StereoLigand::Atom(2),
                StereoLigand::Atom(3),
                StereoLigand::Atom(4),
            ],
            winding: Winding::CounterClockwise,
        }],
    ];
    let basic = parse_reaction(input, &SmilesIoConfig::default()).unwrap();
    let extended =
        parse_extended_reaction_smiles_bytes_with(input, &SmilesIoConfig::default()).unwrap();
    assert_eq!(
        [
            basic.reactants.stereo_atoms,
            basic.agents.stereo_atoms,
            basic.products.stereo_atoms
        ],
        expected
    );
    assert_eq!(
        [
            extended.reactants.stereo_atoms,
            extended.agents.stereo_atoms,
            extended.products.stereo_atoms
        ],
        expected
    );
}

#[rstest]
#[case::reverse_closures(
    "[C@]12(CCC2)CCC1",
    vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::Atom(6), StereoLigand::Atom(3), StereoLigand::Atom(1), StereoLigand::Atom(4)], winding: Winding::CounterClockwise }],
    vec![(0, 6), (0, 3), (0, 1), (1, 2), (2, 3), (0, 4), (4, 5), (5, 6)]
)]
#[case::restored_label_reuse(
    "[C@]1(CC1)1CC1",
    vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::Atom(2), StereoLigand::Atom(1), StereoLigand::Atom(4), StereoLigand::Atom(3)], winding: Winding::CounterClockwise }],
    vec![(0, 2), (0, 1), (1, 2), (0, 4), (0, 3), (3, 4)]
)]
#[case::connected_traversal_roots(
    "[C@H]1(F)Cl.[C@@H]1(Br)I",
    vec![
        StereoAtom { atom: 0, ligands: vec![StereoLigand::ImplicitHydrogen, StereoLigand::Atom(3), StereoLigand::Atom(1), StereoLigand::Atom(2)], winding: Winding::CounterClockwise },
        StereoAtom { atom: 3, ligands: vec![StereoLigand::ImplicitHydrogen, StereoLigand::Atom(0), StereoLigand::Atom(4), StereoLigand::Atom(5)], winding: Winding::Clockwise },
    ],
    vec![(0, 3), (0, 1), (0, 2), (3, 4), (3, 5)]
)]
#[case::nested_marked_branch(
    "[C@H](F)([C@@H](Cl)Br)I",
    vec![
        StereoAtom { atom: 0, ligands: vec![StereoLigand::ImplicitHydrogen, StereoLigand::Atom(1), StereoLigand::Atom(2), StereoLigand::Atom(5)], winding: Winding::CounterClockwise },
        StereoAtom { atom: 2, ligands: vec![StereoLigand::Atom(0), StereoLigand::ImplicitHydrogen, StereoLigand::Atom(3), StereoLigand::Atom(4)], winding: Winding::Clockwise },
    ],
    vec![(0, 1), (0, 2), (2, 3), (2, 4), (0, 5)]
)]
#[case::terminal_hydrogen(
    "C[C@H]",
    vec![StereoAtom { atom: 1, ligands: vec![StereoLigand::Atom(0), StereoLigand::ImplicitHydrogen], winding: Winding::CounterClockwise }],
    vec![(0, 1)]
)]
#[case::terminal_repeated_hydrogen(
    "C[C@H2]",
    vec![StereoAtom { atom: 1, ligands: vec![StereoLigand::Atom(0), StereoLigand::ImplicitHydrogen, StereoLigand::ImplicitHydrogen], winding: Winding::CounterClockwise }],
    vec![(0, 1)]
)]
#[case::isolated_hydrogen(
    "[C@H]",
    vec![StereoAtom { atom: 0, ligands: vec![StereoLigand::ImplicitHydrogen], winding: Winding::CounterClockwise }],
    vec![]
)]
fn test_parse_molecule_assembly(
    #[case] input: &str,
    #[case] expected_stereo: Vec<StereoAtom>,
    #[case] expected_bonds: Vec<(u32, u32)>,
) {
    let basic = parse_molecule(input.as_bytes(), &SmilesIoConfig::default()).unwrap();
    let extended =
        parse_extended_smiles_bytes_with(input.as_bytes(), &SmilesIoConfig::default()).unwrap();
    assert_eq!(basic.stereo_atoms, expected_stereo);
    assert_eq!(extended.stereo_atoms, expected_stereo);
    assert_eq!(
        basic
            .bonds
            .iter()
            .map(|bond| (bond.atoms.first(), bond.atoms.second()))
            .collect::<Vec<_>>(),
        expected_bonds,
    );
    assert_eq!(
        extended
            .bonds
            .iter()
            .map(|bond| (bond.atoms.first(), bond.atoms.second()))
            .collect::<Vec<_>>(),
        expected_bonds,
    );
}

#[rstest]
#[case::forward(vec!["F/C=C/Cl", "Cl/C=C/F"], vec![(0, 1, BondOrder::Single, Some(BondDirection::Rising)), (1, 2, BondOrder::Double, None), (2, 3, BondOrder::Single, Some(BondDirection::Rising))])]
#[case::both_signs(vec!["F\\C=C\\Cl"], vec![(0, 1, BondOrder::Single, Some(BondDirection::Falling)), (1, 2, BondOrder::Double, None), (2, 3, BondOrder::Single, Some(BondDirection::Falling))])]
#[case::opposite(vec!["F/C=C\\Cl"], vec![(0, 1, BondOrder::Single, Some(BondDirection::Rising)), (1, 2, BondOrder::Double, None), (2, 3, BondOrder::Single, Some(BondDirection::Falling))])]
#[case::ring_marker(vec!["C/C=C1CO\\1", "C/C=C/1CO1"], vec![(0, 1, BondOrder::Single, Some(BondDirection::Rising)), (1, 2, BondOrder::Double, None), (2, 4, BondOrder::Single, Some(BondDirection::Rising)), (2, 3, BondOrder::Single, None), (3, 4, BondOrder::Single, None)])]
#[case::one_sided(vec!["F/C=CCl"], vec![(0, 1, BondOrder::Single, Some(BondDirection::Rising)), (1, 2, BondOrder::Double, None), (2, 3, BondOrder::Single, None)])]
fn test_parse_molecule_bond_stereo(
    #[case] inputs: Vec<&str>,
    #[case] expected: Vec<(u32, u32, BondOrder, Option<BondDirection>)>,
) {
    for input in inputs {
        let basic = parse_molecule(input.as_bytes(), &SmilesIoConfig::default()).unwrap();
        let extended = parse_extended_smiles_bytes(input.as_bytes()).unwrap();
        assert_eq!(ExtendedMolecule::from(basic.clone()), extended);
        let converted = Molecule::try_from(extended).unwrap();
        assert_eq!(converted, basic);
        for table in [basic, converted] {
            assert_eq!(
                table
                    .bonds
                    .iter()
                    .map(|bond| (
                        bond.atoms.first(),
                        bond.atoms.second(),
                        bond.order,
                        bond.direction
                    ))
                    .collect::<Vec<_>>(),
                expected
            );
        }
    }
}

#[rstest]
#[case::rising("C/1CC/1", ParseError::MismatchedRingBondDirections { pos: 6, open_pos: 2 })]
#[case::falling("C\\1CC\\1", ParseError::MismatchedRingBondDirections { pos: 6, open_pos: 2 })]
#[case::large_label("C/%99CC/%99", ParseError::MismatchedRingBondDirections { pos: 8, open_pos: 2 })]
fn test_parse_molecule_bond_stereo_error(#[case] input: &str, #[case] expected: ParseError) {
    assert_eq!(
        parse_molecule(input.as_bytes(), &SmilesIoConfig::default()),
        Err(expected.clone())
    );
    assert_eq!(parse_extended_smiles_bytes(input.as_bytes()), Err(expected));
}
