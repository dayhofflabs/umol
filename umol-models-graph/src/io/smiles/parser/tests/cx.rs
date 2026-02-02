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

#[rstest]
#[case::too_many_coords(
    b"C |(0,0;1,1)|",
    ParseError::AtomIndexOutOfBounds { atom_idx: 1 }
)]
fn cx_coordinates_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
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

#[rstest]
#[case::label_out_of_range(b"C |$a;b$|", ParseError::AtomIndexOutOfBounds { atom_idx: 1 })]
fn cx_atom_labels_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
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
#[case::atom_index_out_of_range(
    b"C |^1:1|",
    ParseError::AtomIndexOutOfBounds { atom_idx: 1 }
)]
fn cx_radicals_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
}

#[rstest]
#[case::either(b"CCC |w:0.0|", 0usize, BondWedge::Either)]
#[case::either_up(b"CCC |wU:1.0|", 0usize, BondWedge::EitherUp)]
#[case::either_down(b"CCC |wD:2.1|", 1usize, BondWedge::EitherDown)]
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
#[case::atom_not_in_bond(b"CCC |w:0.1|", ParseError::MismatchedAtomBondIndices { atom_idx: 0, bond_idx: 1 })]
fn cx_wiggly_bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
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
#[case::bond_index_out_of_range(
    b"C=C |c:1|",
    ParseError::BondIndexOutOfBounds { bond_idx: 1 }
)]
fn cx_cis_bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
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
#[case::unspec(b"C=C |ctu:0|", 0usize)]
fn cx_unspec_bonds(#[case] input: &[u8], #[case] bond_idx: usize) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds[bond_idx].stereo, Some(BondStereo::Either));

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.bonds[bond_idx].stereo, Some(BondStereo::Either));
}

#[rstest]
#[case::donating(b"CC |C:0.0|", Ok(Some(BondDonation::Donating)))]
#[case::accepting(b"CC |C:1.0|", Ok(Some(BondDonation::Accepting)))]
#[case::bond_out_of_range(b"CC |C:0.1|", Err(ParseError::BondIndexOutOfBounds { bond_idx: 1 }))]
fn cx_coordinate_bonds(
    #[case] input: &[u8],
    #[case] expected: Result<Option<BondDonation>, ParseError>,
) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    match expected.clone() {
        Ok(expected_donation) => {
            assert!(
                res.is_ok(),
                "{:?} should have succeeded: {:?}",
                input_str,
                res
            );
            let mol = res.unwrap();
            assert_eq!(mol.bonds.len(), 1);
            assert_eq!(mol.bonds[0].order, BondOrder::Single);
            assert_eq!(mol.bonds[0].donation, expected_donation);
        }
        Err(expected_err) => {
            assert!(
                res.is_err(),
                "{:?} should have failed: {:?}",
                input_str,
                res
            );
            assert_eq!(res.unwrap_err(), expected_err);
        }
    }

    let res = parse_extended_cxsmiles(input);
    match expected {
        Ok(expected_donation) => {
            assert!(
                res.is_ok(),
                "{:?} should have succeeded: {:?}",
                input_str,
                res
            );
            let mol = res.unwrap();
            assert_eq!(mol.bonds.len(), 1);
            assert_eq!(mol.bonds[0].order, BondOrder::Single);
            assert_eq!(mol.bonds[0].donation, expected_donation);
        }
        Err(expected_err) => {
            assert!(
                res.is_err(),
                "{:?} should have failed: {:?}",
                input_str,
                res
            );
            assert_eq!(res.unwrap_err(), expected_err);
        }
    }
}

#[rstest]
#[case::hbond(b"CC |H:0.0|", Ok((BondOrder::Zero, Some(BondNoncovalent::Hydrogen))))]
#[case::bond_out_of_range(b"CC |H:0.1|", Err(ParseError::BondIndexOutOfBounds { bond_idx: 1 }))]
fn cx_hydrogen_bonds(
    #[case] input: &[u8],
    #[case] expected: Result<(BondOrder, Option<BondNoncovalent>), ParseError>,
) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    match expected.clone() {
        Ok((expected_order, expected_noncovalent)) => {
            assert!(
                res.is_ok(),
                "{:?} should have succeeded: {:?}",
                input_str,
                res
            );
            let mol = res.unwrap();
            assert_eq!(mol.bonds.len(), 1);
            assert_eq!(mol.bonds[0].order, expected_order);
            assert_eq!(mol.bonds[0].noncovalent, expected_noncovalent);
        }
        Err(expected_err) => {
            assert!(
                res.is_err(),
                "{:?} should have failed: {:?}",
                input_str,
                res
            );
            assert_eq!(res.unwrap_err(), expected_err);
        }
    }

    let res = parse_extended_cxsmiles(input);
    match expected {
        Ok((expected_order, expected_noncovalent)) => {
            assert!(
                res.is_ok(),
                "{:?} should have succeeded: {:?}",
                input_str,
                res
            );
            let mol = res.unwrap();
            assert_eq!(mol.bonds.len(), 1);
            assert_eq!(mol.bonds[0].order, expected_order);
            assert_eq!(mol.bonds[0].noncovalent, expected_noncovalent);
        }
        Err(expected_err) => {
            assert!(
                res.is_err(),
                "{:?} should have failed: {:?}",
                input_str,
                res
            );
            assert_eq!(res.unwrap_err(), expected_err);
        }
    }
}

#[rstest]
#[case::atom_not_in_bond(
    b"CCC |C:0.1|",
    ParseError::MismatchedAtomBondIndices { atom_idx: 0, bond_idx: 1 }
)]
#[case::atom_not_in_hbond(
    b"CCC |H:0.1|",
    ParseError::MismatchedAtomBondIndices { atom_idx: 0, bond_idx: 1 }
)]
fn cx_bond_indexed_tags_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
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
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxTag { pos: 0 });

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

#[rstest]
#[case::atom_index_out_of_range(
    b"CCC |f:0.3|",
    ParseError::AtomIndexOutOfBounds { atom_idx: 3 }
)]
fn cx_fragment_groups_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxTag { pos: 0 });

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
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
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxTag { pos: 0 });

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.stereo_interpretation, expected_interpretation);

    if let Some((idx, set)) = expected_group {
        let cx_data = mol.cx_data.as_ref().expect("cx_data should be present");
        assert_eq!(cx_data.stereo_groups.get(&idx), Some(&set));
    } else {
        assert!(mol.cx_data.is_none());
    }
}

#[rstest]
#[case::absolute_atom_index_out_of_range(
    b"CC |a:0,2|",
    ParseError::AtomIndexOutOfBounds { atom_idx: 2 }
)]
#[case::or_group_atom_index_out_of_range(
    b"CC |o1:2|",
    ParseError::AtomIndexOutOfBounds { atom_idx: 2 }
)]
#[case::and_group_atom_index_out_of_range(
    b"CC |&1:2|",
    ParseError::AtomIndexOutOfBounds { atom_idx: 2 }
)]
fn cx_stereo_groups_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxTag { pos: 0 });

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
}

#[rstest]
#[case::relative(b"C |r|")]
fn cx_relative_stereo(#[case] input: &[u8]) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxTag { pos: 0 });

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
    assert_eq!(mol.cx_data, None);
}

#[rstest]
#[case::with_component_list(b"C |r:0|", ParseError::InvalidCxTag { pos: 0 })]
fn cx_relative_stereo_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
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
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxTag { pos: 0 });

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

#[rstest]
#[case::atom_index_out_of_range(
    b"C |atomProp:1.key.value|",
    ParseError::AtomIndexOutOfBounds { atom_idx: 1 }
)]
fn cx_atom_properties_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxTag { pos: 0 });

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::lp_explicit_count(b"CC |lp:0:2,1:1|", vec![Some(2), Some(1)])]
#[case::lp_implicit_count(b"CC |LP:0,1|", vec![Some(1), Some(1)])]
#[case::lp_single(b"CC |lp:0:3|", vec![Some(3), None])]
#[case::lp_uppercase_single(b"CC |LP:1|", vec![None, Some(1)])]
fn cx_lone_pairs(#[case] input: &[u8], #[case] expected: Vec<Option<u8>>) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    for (i, &exp) in expected.iter().enumerate() {
        assert_eq!(mol.atoms[i].lone_pairs, exp, "atom {} lone_pairs", i);
    }

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    for (i, &exp) in expected.iter().enumerate() {
        assert_eq!(mol.atoms[i].lone_pairs, exp, "atom {} lone_pairs", i);
    }
}

#[rstest]
#[case::atom_index_out_of_range(b"C |lp:1:2|", ParseError::AtomIndexOutOfBounds { atom_idx: 1 })]
fn cx_lone_pairs_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::ferrocene(b"[Fe]c1cccc1.c1cccc1 |m:0:1.2.3.4.5,0:6.7.8.9.10|", 2, vec![(0, vec![1, 2, 3, 4, 5]), (0, vec![6, 7, 8, 9, 10])])]
#[case::single_multicenter(b"CC[Zr]CC |m:2:0.1|", 1, vec![(2, vec![0, 1])])]
fn cx_multicenter_bonds(
    #[case] input: &[u8],
    #[case] expected_count: usize,
    #[case] expected_bonds: Vec<(u32, Vec<u32>)>,
) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.multicenter_bond_count(), expected_count);
    for (i, (center, ligands)) in expected_bonds.iter().enumerate() {
        let bond = &mol.multicenter_bonds[i];
        assert_eq!(bond.contributions.len(), 2, "bond {} contributions", i);
        assert_eq!(bond.contributions[0].atoms, *ligands, "ligands for bond {}", i);
        assert_eq!(bond.contributions[1].atoms, vec![*center], "center for bond {}", i);
    }

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.multicenter_bond_count(), expected_count);
}

#[rustfmt::skip]
#[rstest]
#[case::center_out_of_range(b"CC |m:5:0.1|", ParseError::AtomIndexOutOfBounds { atom_idx: 5 })]
#[case::ligand_out_of_range(b"CC |m:0:1.5|", ParseError::AtomIndexOutOfBounds { atom_idx: 5 })]
fn cx_multicenter_bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
}
