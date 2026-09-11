use rstest::rstest;

use super::super::*;
use crate::table_ir::{
    BondConfiguration, BondRelation, StereoAtom, StereoBond, StereoLigand, Winding,
};

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
#[case::forward(vec!["F/C=C/Cl", "Cl/C=C/F", "F\\C=C\\Cl"], vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide } }])]
#[case::opposite(vec!["F/C=C\\Cl"], vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::SameSide } }])]
#[case::ring(vec!["C/C=C1CO\\1", "C/C=C/1CO1"], vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::SameSide } }])]
#[case::partial(vec!["F/C=CCl"], vec![])]
fn test_parse_molecule_bond_stereo(#[case] inputs: Vec<&str>, #[case] expected: Vec<StereoBond>) {
    for input in inputs {
        let basic = parse_molecule(input.as_bytes(), &SmilesIoConfig::default()).unwrap();
        let extended = parse_extended_smiles_bytes(input.as_bytes()).unwrap();
        assert_eq!(ExtendedMolecule::from(basic.clone()), extended);
        assert_eq!(Molecule::try_from(extended).unwrap(), basic);
        assert_eq!(basic.stereo_bonds, expected);
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

#[rstest]
#[case::directions_cis(r"F/C=C\F |c:1|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::SameSide } }])]
#[case::directions_trans("F/C=C/F |t:1|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide } }])]
#[case::duplicate("FC=CF |c:1,c:1|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::SameSide } }])]
#[case::partial_either("F/C=CF |ctu:1|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Either }])]
#[case::terminal_either("C=C |ctu:0|", vec![StereoBond { bond: 0, configuration: BondConfiguration::Either }])]
#[case::wavy_either("FC=CF |w:1.0,ctu:1|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Either }])]
#[case::labels("F/C=C/F |$first$|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide } }])]
#[case::new_before_existing("FC=CF.F/C=C/F |c:1,c:1|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::SameSide } }, StereoBond { bond: 4, configuration: BondConfiguration::Framed { references: [4,7], relation: BondRelation::OppositeSide } }])]
#[case::new_after_existing("F/C=C/F.FC=CF |c:4,c:4|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide } }, StereoBond { bond: 4, configuration: BondConfiguration::Framed { references: [4,7], relation: BondRelation::SameSide } }])]
#[case::unsorted_duplicates("FC=CF.FC=CF |t:4,c:1,t:4,c:1|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::SameSide } }, StereoBond { bond: 4, configuration: BondConfiguration::Framed { references: [4,7], relation: BondRelation::OppositeSide } }])]
#[case::overwritten_wavy("FC=CF |w:1.0,w:0.0|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Either }])]
#[case::existing_geometry("F/C=C/F |(0,1,;0,0,;2,0,;2,-1,)|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::OppositeSide } }])]
#[case::geometry_only("FC=CF |(0,1,;0,0,;2,0,;2,1,)|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Framed { references: [0,3], relation: BondRelation::SameSide } }])]
#[case::either_geometry("FC=CF |ctu:1,(0,1,;0,0,;2,0,;2,1,)|", vec![StereoBond { bond: 1, configuration: BondConfiguration::Either }])]
fn test_parse_molecule_cx_frames(#[case] input: &str, #[case] expected: Vec<StereoBond>) {
    let config = SmilesIoConfig::chemaxon();
    let basic = parse_molecule(input.as_bytes(), &config).unwrap();
    let extended = parse_extended_smiles_bytes_with(input.as_bytes(), &config).unwrap();
    assert_eq!(basic.stereo_bonds, expected);
    assert_eq!(ExtendedMolecule::from(basic.clone()), extended);
    assert_eq!(Molecule::try_from(extended).unwrap(), basic);
}

#[rstest]
#[case::direction_code("F/C=C/F |c:1|", ParseError::ConflictingBondConfiguration { bond: 1 })]
#[case::direction_either("F/C=C/F |ctu:1|", ParseError::ConflictingBondConfiguration { bond: 1 })]
#[case::codes("FC=CF |c:1,t:1|", ParseError::ConflictingBondConfiguration { bond: 1 })]
#[case::code_either("FC=CF |ctu:1,c:1|", ParseError::ConflictingBondConfiguration { bond: 1 })]
#[case::wavy_code("FC=CF |w:1.0,c:1|", ParseError::ConflictingBondConfiguration { bond: 1 })]
#[case::overwritten_wavy("FC=CF |w:1.0,w:0.0,c:1|", ParseError::ConflictingBondConfiguration { bond: 1 })]
#[case::geometry_code("FC=CF |(0,1,;0,0,;2,0,;2,-1,),c:1|", ParseError::ConflictingBondConfiguration { bond: 1 })]
#[case::single("CC |ctu:0|", ParseError::UnsupportedStereoBond { bond: 0 })]
#[case::existing_geometry_conflict("F/C=C/F |(0,1,;0,0,;2,0,;2,1,)|", ParseError::ConflictingBondConfiguration { bond: 1 })]
#[case::existing_short_positions("F/C=C/F |(0,1,;0,0,;2,0,)|", ParseError::MissingPosition { atom: 3 })]
#[case::existing_site_order_change("F/C=C/F |H:1.1|", ParseError::UnsupportedStereoBond { bond: 1 })]
#[case::unsorted_conflict("FC=CF.FC=CF |t:4,c:1,c:4|", ParseError::ConflictingBondConfiguration { bond: 4 })]
fn test_parse_molecule_cx_frames_error(#[case] input: &str, #[case] expected: ParseError) {
    let config = SmilesIoConfig::chemaxon();
    assert_eq!(
        parse_molecule(input.as_bytes(), &config),
        Err(expected.clone())
    );
    assert_eq!(
        parse_extended_smiles_bytes_with(input.as_bytes(), &config),
        Err(expected)
    );
}

#[rstest]
#[case::ordinary(b"F/C=C/F>C=C>FC=CF |ctu:3,c:5|", 1, [0,3])]
#[case::uneven(b"CC.F/C=C/F>C=C>FC=CF |ctu:4,c:6|", 2, [2,5])]
fn test_parse_reaction_bond_frames(
    #[case] input: &[u8],
    #[case] bond: u32,
    #[case] references: [u32; 2],
) {
    let config = SmilesIoConfig::chemaxon();
    let basic = parse_reaction(input, &config).unwrap();
    let extended = parse_extended_reaction_smiles_bytes_with(input, &config).unwrap();
    let expected = [
        vec![StereoBond {
            bond,
            configuration: BondConfiguration::Framed {
                references,
                relation: BondRelation::OppositeSide,
            },
        }],
        vec![StereoBond {
            bond: 0,
            configuration: BondConfiguration::Either,
        }],
        vec![StereoBond {
            bond: 1,
            configuration: BondConfiguration::Framed {
                references: [0, 3],
                relation: BondRelation::SameSide,
            },
        }],
    ];
    assert_eq!(
        [
            basic.reactants.stereo_bonds,
            basic.agents.stereo_bonds,
            basic.products.stereo_bonds
        ],
        expected
    );
    assert_eq!(
        [
            extended.reactants.stereo_bonds,
            extended.agents.stereo_bonds,
            extended.products.stereo_bonds
        ],
        expected
    );
}

#[rstest]
#[case::code("|t:1|")]
#[case::geometry("|(0,1,;0,0,;2,0,;2,-1,)|")]
#[case::geometry_conflict("|(0,1,;0,0,;2,0,;2,1,)|")]
#[case::site_order_change("|H:1.1|")]
#[case::substituent_order_change("|C:0.0|")]
fn test_update_molecule_lookup(#[case] annotations: &str) {
    let config = SmilesIoConfig::chemaxon();
    let neighbors = OnceCell::new();
    let (_, (mut molecule, _, _)) = parse_smiles_inner(
        b"F/C=C/F",
        0,
        false,
        true,
        config.syntax_flags,
        None,
        &neighbors,
    )
    .unwrap();
    let entries = parse_cx_annotations(annotations.as_bytes(), config.syntax_flags).unwrap();
    let mut fresh = molecule.clone();
    let expected = update_molecule(&mut fresh, entries.clone(), &OnceCell::new()).map(|()| fresh);
    let actual = update_molecule(&mut molecule, entries, &neighbors).map(|()| molecule);
    assert_eq!(actual, expected);
}
