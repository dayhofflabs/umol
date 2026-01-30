use bstr::ByteSlice;
use pretty_assertions::assert_eq;
use rstest::*;
use umol_data::SpinMultiplicity;

use super::super::*;
use crate::position::Point3D;
use crate::table_ir::{
    BondDonation, BondNoncovalent, BondOrder, BondStereo, BondWedge, StereoInterpretation,
    StereoSet, StereoSetMode, UnpairedElectrons,
};

fn parse_basic_cxsmiles(input: &[u8]) -> Result<Molecule, ParseError> {
    parse_smiles_bytes_with(input, &SmilesIoConfig::basic_chemaxon())
}

fn parse_extended_cxsmiles(input: &[u8]) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_smiles_bytes_with(input, &SmilesIoConfig::chemaxon())
}

#[rustfmt::skip]
#[rstest]
#[case::single_atom_3d(b"C |(1,2,3)|", Some(vec![Point3D::new(1.0, 2.0, 3.0)]))]
#[case::two_atoms_2d(b"CC |(1.5,2.5;3.5,4.5)|", Some(vec![Point3D::new(1.5, 2.5, 0.0), Point3D::new(3.5, 4.5, 0.0)]))]
#[case::empty_coordinates(b"C |()|", Some(vec![]))]
fn cx_coordinates(#[case] input: &[u8], #[case] expected: Option<Vec<Point3D>>) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.positions, expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.positions, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::two_labels(b"CC |$C1;N1$|", vec![Some("C1".to_string()), Some("N1".to_string())])]
#[case::single_label(b"CC |$C1$|", vec![Some("C1".to_string()), None])]
fn cx_atom_labels(#[case] input: &[u8], #[case] expected: Vec<Option<String>>) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    let labels: Vec<_> = mol.atoms.iter().map(|a| a.label.clone()).collect();
    assert_eq!(labels, expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    let labels: Vec<_> = mol.atoms.iter().map(|a| a.label.clone()).collect();
    assert_eq!(labels, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::two_values(b"CC |$_AV:v0;v1$|", vec![Some("v0".to_string()), Some("v1".to_string())])]
#[case::single_value(b"CC |$_AV:v0$|", vec![Some("v0".to_string()), None])]
fn cx_atom_values(#[case] input: &[u8], #[case] expected: Vec<Option<String>>) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    let values: Vec<_> = mol.atoms.iter().map(|a| a.value.clone()).collect();
    assert_eq!(values, expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    let values: Vec<_> = mol.atoms.iter().map(|a| a.value.clone()).collect();
    assert_eq!(values, expected);
}

#[rustfmt::skip]
#[rstest]
#[case::monovalent(b"C |^1:0|", UnpairedElectrons { count: 1, multiplicity: None })]
#[case::divalent_triplet(b"C |^4:0|", UnpairedElectrons { count: 2, multiplicity: Some(SpinMultiplicity::Triplet)})]
fn cx_radicals(#[case] input: &[u8], #[case] expected: UnpairedElectrons) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.atoms[0].unpaired_electrons, Some(expected));

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.atoms[0].unpaired_electrons, Some(expected));
}

#[rstest]
#[case::either(b"CCC |w:0.0|", 0usize, BondWedge::Either)]
#[case::either_up(b"CCC |wU:1.0|", 1usize, BondWedge::EitherUp)]
#[case::either_down(b"CCC |wD:0.0|", 0usize, BondWedge::EitherDown)]
fn cx_wiggly_bonds(#[case] input: &[u8], #[case] bond_idx: usize, #[case] wedge: BondWedge) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds[bond_idx].wedge, Some(wedge));

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds[bond_idx].wedge, Some(wedge));
}

#[rstest]
#[case::cis(b"C=C |c:0|", 0usize)]
fn cx_cis_bonds(#[case] input: &[u8], #[case] bond_idx: usize) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds[bond_idx].stereo, Some(BondStereo::Cis));

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds[bond_idx].stereo, Some(BondStereo::Cis));
}

#[rstest]
#[case::trans(b"C=C |t:0|", 0usize)]
fn cx_trans_bonds(#[case] input: &[u8], #[case] bond_idx: usize) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds[bond_idx].stereo, Some(BondStereo::Trans));

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds[bond_idx].stereo, Some(BondStereo::Trans));
}

#[rstest]
#[case::dative_bond(b"CCC |C:0.2|", (0u32, 2u32))]
fn cx_coordinate_bonds(#[case] input: &[u8], #[case] atoms: (u32, u32)) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds.len(), 3);
    let b = mol.bonds.last().expect("dative bond appended");
    assert_eq!(b.order, BondOrder::Single);
    assert_eq!(b.donation, Some(BondDonation::Donating));
    assert_eq!(b.atoms.as_tuple(), atoms);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds.len(), 3);
    let b = mol.bonds.last().expect("dative bond appended");
    assert_eq!(b.order, BondOrder::Single);
    assert_eq!(b.donation, Some(BondDonation::Donating));
    assert_eq!(b.atoms.as_tuple(), atoms);
}

#[rstest]
#[case::hbond(b"CCC |H:0.2|", (0u32, 2u32))]
fn cx_hydrogen_bonds(#[case] input: &[u8], #[case] atoms: (u32, u32)) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds.len(), 3);
    let b = mol.bonds.last().expect("hydrogen bond appended");
    assert_eq!(b.order, BondOrder::Zero);
    assert_eq!(b.noncovalent, Some(BondNoncovalent::Hydrogen));
    assert_eq!(b.atoms.as_tuple(), atoms);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds.len(), 3);
    let b = mol.bonds.last().expect("hydrogen bond appended");
    assert_eq!(b.order, BondOrder::Zero);
    assert_eq!(b.noncovalent, Some(BondNoncovalent::Hydrogen));
    assert_eq!(b.atoms.as_tuple(), atoms);
}

#[rstest]
#[case::fragment_groups(b"CCC |f:0.1,2|", vec![vec![0u32, 1u32], vec![2u32]])]
fn cx_fragment_groups(#[case] input: &[u8], #[case] expected: Vec<Vec<u32>>) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxProperty { pos: 0 });

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    let components = mol
        .cx_data
        .as_ref()
        .and_then(|d| d.components.as_ref())
        .cloned();
    assert_eq!(components, Some(expected));
}

#[rustfmt::skip]
#[rstest]
#[case::absolute(b"CC |a:0,1|", Some(StereoInterpretation::Absolute), None)]
#[case::or_group(b"CC |o1:0,1|", None, Some(( 1u32, StereoSet { atoms: vec![0u32, 1u32], mode: StereoSetMode::Correlated})))]
#[case::and_group(b"CC |&2:1|", None, Some(( 2u32, StereoSet { atoms: vec![1u32], mode: StereoSetMode::Independent})))]
fn cx_stereo_groups(
    #[case] input: &[u8],
    #[case] expected_interpretation: Option<StereoInterpretation>,
    #[case] expected_group: Option<(u32, StereoSet)>,
) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxProperty { pos: 0 });

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.stereo_interpretation, expected_interpretation);

    let cx_data = mol.cx_data.as_ref().expect("cx_data should be present");
    assert_eq!(cx_data.stereo_interpretation, expected_interpretation);
    if let Some((idx, set)) = expected_group {
        assert_eq!(cx_data.stereo_groups.get(&idx), Some(&set));
    } else {
        assert!(cx_data.stereo_groups.is_empty());
    }
}

#[rstest]
#[case::relative(b"C |r|")]
#[case::relative_with_components(b"C |r:0|")]
fn cx_relative_stereo(#[case] input: &[u8]) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxProperty { pos: 0 });

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(
        mol.stereo_interpretation,
        Some(StereoInterpretation::Relative)
    );
    assert_eq!(
        mol.cx_data.as_ref().and_then(|d| d.stereo_interpretation),
        Some(StereoInterpretation::Relative)
    );
}

#[rstest]
#[case::atom_properties(b"C |atomProp:0.key.value|", "key", "value")]
fn cx_atom_properties(#[case] input: &[u8], #[case] key: &str, #[case] value: &str) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxProperty { pos: 0 });

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.cx_data, None, "atomProp does not populate cx_data");
    assert_eq!(mol.atoms[0].properties.get(key), Some(&value.to_string()));
}
