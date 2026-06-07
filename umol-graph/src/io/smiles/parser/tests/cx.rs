use bstr::ByteSlice;
use pretty_assertions::assert_eq;
use rstest::*;
use umol_shared::spin::SpinMultiplicity;

use super::super::*;
use crate::position::Point3D;
use crate::table_ir::bond::BondNoncovalent;
use crate::table_ir::{
    BondDonation, BondOrder, BondStereo, BondWedge, ConfigurationScope, LinkAtom, RingBondCount,
    SGroupBracketCoords, SGroupBracketOrientation, SGroupBracketStyle, SGroupConnectivity,
    SGroupDataType, SGroupType, StereoSet, StereoSetRelation, SubstitutionCount, UnsaturatedAtom,
};

fn parse_basic_cxsmiles(input: &[u8]) -> Result<Molecule, ParseError> {
    parse_smiles_bytes_to_table_ir_with(input, &SmilesIoConfig::basic_chemaxon())
}

fn parse_extended_cxsmiles(input: &[u8]) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_smiles_bytes_with(input, &SmilesIoConfig::chemaxon())
}

#[rustfmt::skip]
#[rstest]
#[case::single_atom_3d(b"C |(1,2,3)|", Some(vec![Point3D::new(1.0, 2.0, 3.0)]))]
#[case::two_atoms_2d(b"CC |(1.5,2.5;3.5,4.5)|", Some(vec![Point3D::new(1.5, 2.5, 0.0), Point3D::new(3.5, 4.5, 0.0)]))]
#[case::empty_coordinates(b"C |()|", Some(vec![]))]
fn test_cx_coordinates(#[case] input: &[u8], #[case] expected: Option<Vec<Point3D>>) {
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
#[case::too_many_coords(b"C |(0,0;1,1)|", ParseError::AtomIndexOutOfBounds { atom_idx: 1 })]
fn test_cx_coordinates_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
fn test_cx_atom_labels(#[case] input: &[u8], #[case] expected: Vec<Option<String>>) {
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
fn test_cx_atom_labels_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
fn test_cx_atom_values(#[case] input: &[u8], #[case] expected: Vec<Option<String>>) {
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
#[case::monovalent(b"C |^1:0|", 1, None)]
#[case::divalent_triplet(b"C |^4:0|", 2, Some(SpinMultiplicity::Triplet))]
fn test_cx_radicals(
    #[case] input: &[u8],
    #[case] expected_unpaired: u8,
    #[case] expected_multiplicity: Option<SpinMultiplicity>,
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
    assert_eq!(mol.atoms[0].unpaired_electrons, Some(expected_unpaired));
    assert_eq!(mol.atoms[0].multiplicity, expected_multiplicity);

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should have succeeded: {:?}",
        input_str,
        res
    );
    let mol = res.unwrap();
    assert_eq!(mol.atoms[0].unpaired_electrons, Some(expected_unpaired));
    assert_eq!(mol.atoms[0].multiplicity, expected_multiplicity);
}

#[rstest]
#[case::atom_index_out_of_range(b"C |^1:1|", ParseError::AtomIndexOutOfBounds { atom_idx: 1 })]
fn test_cx_radicals_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
fn test_cx_wiggly_bonds(#[case] input: &[u8], #[case] bond_idx: usize, #[case] wedge: BondWedge) {
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
fn test_cx_wiggly_bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
fn test_cx_cis_bonds(#[case] input: &[u8], #[case] bond_idx: usize) {
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
#[case::bond_index_out_of_range(b"C=C |c:1|", ParseError::BondIndexOutOfBounds { bond_idx: 1 })]
fn test_cx_cis_bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
fn test_cx_trans_bonds(#[case] input: &[u8], #[case] bond_idx: usize) {
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
fn test_cx_unspec_bonds(#[case] input: &[u8], #[case] bond_idx: usize) {
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
fn test_cx_coordinate_bonds(
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
fn test_cx_hydrogen_bonds(
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
#[case::atom_not_in_bond(b"CCC |C:0.1|", ParseError::MismatchedAtomBondIndices { atom_idx: 0, bond_idx: 1 })]
#[case::atom_not_in_hbond(b"CCC |H:0.1|", ParseError::MismatchedAtomBondIndices { atom_idx: 0, bond_idx: 1 })]
fn test_cx_bond_indexed_tags_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
fn test_cx_fragment_groups(#[case] input: &[u8], #[case] expected: Vec<Vec<u32>>) {
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
#[case::atom_index_out_of_range(b"CCC |f:0.3|", ParseError::AtomIndexOutOfBounds { atom_idx: 3 })]
fn test_cx_fragment_groups_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
#[case::absolute(b"CC |a:0,1|", Some(ConfigurationScope::Absolute), None)]
#[case::or_group(b"CC |o1:0,1|", None, Some(( 1u32, StereoSet { atoms: vec![0u32, 1u32], relation: StereoSetRelation::Correlated})))]
#[case::and_group(b"CC |&2:1|", None, Some(( 2u32, StereoSet { atoms: vec![1u32], relation: StereoSetRelation::Independent})))]
fn test_cx_stereo_groups(
    #[case] input: &[u8],
    #[case] expected_interpretation: Option<ConfigurationScope>,
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
    assert_eq!(mol.configuration_scope, expected_interpretation);

    if let Some((idx, set)) = expected_group {
        let cx_data = mol.cx_data.as_ref().expect("cx_data should be present");
        assert_eq!(cx_data.stereo_groups.get(&idx), Some(&set));
    } else {
        assert!(mol.cx_data.is_none());
    }
}

#[rstest]
#[case::absolute_atom_index_out_of_range(b"CC |a:0,2|", ParseError::AtomIndexOutOfBounds { atom_idx: 2 })]
#[case::or_group_atom_index_out_of_range(b"CC |o1:2|", ParseError::AtomIndexOutOfBounds { atom_idx: 2 })]
#[case::and_group_atom_index_out_of_range(b"CC |&1:2|", ParseError::AtomIndexOutOfBounds { atom_idx: 2 })]
fn test_cx_stereo_groups_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
fn test_cx_relative_stereo(#[case] input: &[u8]) {
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
    assert_eq!(mol.configuration_scope, Some(ConfigurationScope::Relative));
    assert_eq!(mol.cx_data, None);
}

#[rstest]
#[case::with_component_list(b"C |r:0|", ParseError::InvalidCxTag { pos: 0 })]
fn test_cx_relative_stereo_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
fn test_cx_atom_properties(#[case] input: &[u8], #[case] key: &str, #[case] value: &str) {
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

#[rustfmt::skip]
#[rstest]
#[case::atom_index_out_of_range(b"C |atomProp:1.key.value|", ParseError::AtomIndexOutOfBounds { atom_idx: 1 })]
fn test_cx_atom_properties_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
fn test_cx_lone_pairs(#[case] input: &[u8], #[case] expected: Vec<Option<u8>>) {
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
fn test_cx_lone_pairs_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
fn test_cx_multicenter_bonds(
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
        assert_eq!(bond.contributions().len(), 2, "bond {} contributions", i);
        let ligand_contribution = bond
            .contributions()
            .iter()
            .find(|contribution| contribution.atoms().len() > 1)
            .unwrap();
        let center_contribution = bond
            .contributions()
            .iter()
            .find(|contribution| contribution.atoms().len() == 1)
            .unwrap();
        assert_eq!(
            ligand_contribution.atoms(),
            ligands.as_slice(),
            "ligands for bond {}",
            i
        );
        assert_eq!(
            center_contribution.atoms(),
            [*center].as_slice(),
            "center for bond {}",
            i
        );
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
fn test_cx_multicenter_bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
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
#[case::rb_single(b"CCC |rb:1:2|", vec![None, Some(RingBondCount::R2), None])]
#[case::rb_multiple(b"CCCC |rb:0:3,2:4,3:*|", vec![Some(RingBondCount::R3), None, Some(RingBondCount::R4Plus), Some(RingBondCount::AsDrawn)])]
#[case::rb_as_drawn(b"CC |rb:0:*|", vec![Some(RingBondCount::AsDrawn), None])]
#[case::rb_no_ring(b"CC |rb:1:-1|", vec![None, Some(RingBondCount::NoRingBonds)])]
fn test_cx_ring_bond_count(#[case] input: &[u8], #[case] expected: Vec<Option<RingBondCount>>) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed (rb is extended only): {:?}",
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
    for (i, &exp) in expected.iter().enumerate() {
        assert_eq!(mol.atoms[i].ring_bond_count, exp, "atom {} ring_bond_count", i);
    }
}

#[rustfmt::skip]
#[rstest]
#[case::rb_atom_out_of_range(b"CC |rb:2:2|", ParseError::AtomIndexOutOfBounds { atom_idx: 2 })]
fn test_cx_ring_bond_count_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_extended_cxsmiles(input);
    assert!(res.is_err(), "{:?} should have failed: {:?}", input_str, res);
    assert_eq!(res.unwrap_err(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::s_single(b"CCC |s:1:2|", vec![None, Some(SubstitutionCount::S2), None])]
#[case::s_multiple(b"CCCC |s:0:1,2:2,3:*|", vec![Some(SubstitutionCount::S1), None, Some(SubstitutionCount::S2), Some(SubstitutionCount::AsDrawn)])]
#[case::s_as_drawn(b"CC |s:0:*|", vec![Some(SubstitutionCount::AsDrawn), None])]
#[case::s_no_sub(b"CC |s:1:-1|", vec![None, Some(SubstitutionCount::NoSubstitution)])]
#[case::s6_plus(b"CCCCCC |s:5:6|", vec![None, None, None, None, None, Some(SubstitutionCount::S6Plus)])]
fn test_cx_substitution_count(#[case] input: &[u8], #[case] expected: Vec<Option<SubstitutionCount>>) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed (s is extended only): {:?}",
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
    for (i, &exp) in expected.iter().enumerate() {
        assert_eq!(mol.atoms[i].substitution_count, exp, "atom {} substitution_count", i);
    }
}

#[rustfmt::skip]
#[rstest]
#[case::s_atom_out_of_range(b"CC |s:3:2|", ParseError::AtomIndexOutOfBounds { atom_idx: 3 })]
fn test_cx_substitution_count_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();
    let res = parse_extended_cxsmiles(input);
    assert!(res.is_err(), "{:?} should have failed: {:?}", input_str, res);
    assert_eq!(res.unwrap_err(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::u_single(b"CCC |u:1|", vec![None, Some(UnsaturatedAtom), None])]
#[case::u_multiple(b"CCCC |u:0,2,3|", vec![Some(UnsaturatedAtom), None, Some(UnsaturatedAtom), Some(UnsaturatedAtom)])]
fn test_cx_unsaturated(#[case] input: &[u8], #[case] expected: Vec<Option<UnsaturatedAtom>>) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed (u is extended only): {:?}",
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
    for (i, &exp) in expected.iter().enumerate() {
        assert_eq!(mol.atoms[i].unsaturated.is_some(), exp.is_some(), "atom {} unsaturated", i);
    }
}

#[rustfmt::skip]
#[rstest]
#[case::u_atom_out_of_range(b"CC |u:3|", ParseError::AtomIndexOutOfBounds { atom_idx: 3 })]
fn test_cx_unsaturated_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();
    let res = parse_extended_cxsmiles(input);
    assert!(res.is_err(), "{:?} should have failed: {:?}", input_str, res);
    assert_eq!(res.unwrap_err(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::lo_single_center(b"CCCC |LO:1:0.2.3|", 1, vec![(1u32, vec![(0u32, 1u8), (2u32, 2u8), (3u32, 3u8)])])]
#[case::lo_two_centers(b"CCCCC |LO:1:0.2,3:2.4|", 2, vec![(1u32, vec![(0u32, 1u8), (2u32, 2u8)]), (3u32, vec![(2u32, 1u8), (4u32, 2u8)])])]
fn test_cx_ligand_order(
    #[case] input: &[u8],
    #[case] center_count: usize,
    #[case] expected: Vec<(u32, Vec<(u32, u8)>)>,
) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed (LO is extended only): {:?}",
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
    let mut centers_with_lo = 0;
    for (i, atom) in mol.atoms.iter().enumerate() {
        if let Some(ref lo) = atom.ligand_order {
            centers_with_lo += 1;
            let exp = expected.iter().find(|(c, _)| *c == i as u32).expect("unexpected center");
            assert_eq!(lo, &exp.1, "atom {} ligand_order", i);
        }
    }
    assert_eq!(centers_with_lo, center_count);
}

#[rustfmt::skip]
#[rstest]
#[case::lo_center_out_of_range(b"CC |LO:5:0.1|", ParseError::AtomIndexOutOfBounds { atom_idx: 5 })]
#[case::lo_neighbor_out_of_range(b"CC |LO:0:1.5|", ParseError::AtomIndexOutOfBounds { atom_idx: 5 })]
fn test_cx_ligand_order_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

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
#[case::ln_min_max(b"CC |LN:1:1.5|", vec![(1u32, LinkAtom { min_repeat: 1, repeat_count: 5, subs_index1: 0, subs_index2: None })])]
#[case::ln_with_outer(b"CCC |LN:1:1.2.0.2|", vec![(1u32, LinkAtom { min_repeat: 1, repeat_count: 2, subs_index1: 0, subs_index2: Some(2) })])]
#[case::ln_two_entries(b"CCCC |LN:1:0.3,3:1.2.2.3|", vec![
    (1u32, LinkAtom { min_repeat: 0, repeat_count: 3, subs_index1: 0, subs_index2: None }),
    (3u32, LinkAtom { min_repeat: 1, repeat_count: 2, subs_index1: 2, subs_index2: Some(3) }),
])]
fn test_cx_link_nodes(#[case] input: &[u8], #[case] expected: Vec<(u32, LinkAtom)>) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed (LN is extended only): {:?}",
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
    for (idx, exp_link) in expected {
        assert_eq!(mol.atoms[idx as usize].link_atom, Some(exp_link));
    }
}

#[rustfmt::skip]
#[rstest]
#[case::ln_atom_out_of_range(b"CC |LN:5:1.2|", ParseError::AtomIndexOutOfBounds { atom_idx: 5 })]
#[case::bad_value_count(b"CCC |LN:1:1.2.3|", ParseError::InvalidCxTag { pos: 0 })]
fn test_cx_link_nodes_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

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
#[case::sgd_minimal(b"CC |SgD:0,1:MW:150:::::|", 0, vec![0u32, 1u32], "MW", "150")]
#[case::sgd_single_atom(b"C |SgD:0:Name:value:::::|", 0, vec![0u32], "Name", "value")]
fn test_cx_sgroup_data(
    #[case] input: &[u8],
    #[case] sgroup_idx: u32,
    #[case] expected_atoms: Vec<u32>,
    #[case] expected_name: &str,
    #[case] expected_data: &str,
) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed (SgD is extended only): {:?}",
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
    let sgroups = mol.sgroups();
    assert!(!sgroups.is_empty());
    let sg = sgroups.get(&sgroup_idx).unwrap();
    assert_eq!(sg.group_type, SGroupType::Data);
    assert_eq!(sg.atom_indices, expected_atoms);
    let data = sg.data.as_ref().unwrap();
    assert_eq!(data.field_type, SGroupDataType::Text);
    assert_eq!(data.field_name, expected_name);
    assert_eq!(
        data.data_content,
        Some(vec![expected_data.to_string()])
    );
}

#[rustfmt::skip]
#[rstest]
#[case::sg_minimal(b"CC |Sg:n:0,1|", 0, SGroupType::RepeatingUnit, vec![0u32, 1u32], None, None)]
#[case::sg_with_subscript(b"CCC |Sg:n:0,1,2:n|", 0, SGroupType::RepeatingUnit, vec![0u32, 1u32, 2u32], Some("n".to_string()), None)]
#[case::sg_with_connectivity(b"CC |Sg:n:0,1:n:ht|", 0, SGroupType::RepeatingUnit, vec![0u32, 1u32], Some("n".to_string()), Some(SGroupConnectivity::HeadToTail))]
#[case::sg_copolymer_alt(b"CCCC |Sg:co:alt:0,1,2,3:n|", 0, SGroupType::Copolymer, vec![0u32, 1u32, 2u32, 3u32], Some("n".to_string()), None)]
#[case::sg_followed_by_entry(b"CCC |Sg:n:0,1:n,C:1.0|", 0, SGroupType::RepeatingUnit, vec![0u32, 1u32], Some("n".to_string()), None)]
fn test_cx_sgroup(
    #[case] input: &[u8],
    #[case] sgroup_idx: u32,
    #[case] expected_type: SGroupType,
    #[case] expected_atoms: Vec<u32>,
    #[case] expected_subscript: Option<String>,
    #[case] expected_connectivity: Option<SGroupConnectivity>,
) {
    let input_str = input.to_str_lossy();

    let res = parse_basic_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed (Sg is extended only): {:?}",
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
    let sgroups = mol.sgroups();
    assert!(!sgroups.is_empty());
    let sg = sgroups.get(&sgroup_idx).unwrap();
    assert_eq!(sg.group_type, expected_type);
    assert_eq!(sg.atom_indices, expected_atoms);
    assert_eq!(sg.subscript, expected_subscript);
    assert_eq!(sg.connectivity, expected_connectivity);
}

#[rustfmt::skip]
#[rstest]
#[case::sg_with_bracket(b"CC |Sg:n:0,1:n:ht:::s,b,1,2,3,4|", Some(SGroupBracketOrientation::Straight), Some(SGroupBracketStyle::Default),
    Some(SGroupBracketCoords { bracket1: (1.0, 2.0), bracket2: (3.0, 4.0), bracket3: None, bracket4: None }))]
#[case::sg_with_bracket_4coords(b"CC |Sg:n:0,1:n:::d,c,1,2,3,4,5,6,7,8|", Some(SGroupBracketOrientation::Down), Some(SGroupBracketStyle::Curved),
    Some(SGroupBracketCoords { bracket1: (1.0, 2.0), bracket2: (3.0, 4.0), bracket3: Some((5.0, 6.0)), bracket4: Some((7.0, 8.0)) }))]
fn test_cx_sgroup_bracket(
    #[case] input: &[u8],
    #[case] expected_orientation: Option<SGroupBracketOrientation>,
    #[case] expected_style: Option<SGroupBracketStyle>,
    #[case] expected_coords: Option<SGroupBracketCoords>,
) {
    let res = parse_extended_cxsmiles(input);
    assert!(res.is_ok(), "{:?} should succeed: {:?}", input.to_str_lossy(), res);
    let mol = res.unwrap();
    let sg = mol.sgroups().get(&0).unwrap();
    assert_eq!(sg.bracket_orientation, expected_orientation);
    assert_eq!(sg.bracket_style, expected_style);
    assert_eq!(sg.bracket_coords, expected_coords);
}

#[rustfmt::skip]
#[rstest]
#[case::sg_with_flip(b"CC |Sg:n:0,1:n:ht,1|", Some(true))]
#[case::sg_with_flip_false(b"CC |Sg:n:0,1:n:ht,0|", Some(false))]
fn test_cx_sgroup_connectivity_flip(#[case] input: &[u8], #[case] expected_flip: Option<bool>) {
    let res = parse_extended_cxsmiles(input);
    assert!(res.is_ok(), "{:?} should succeed: {:?}", input.to_str_lossy(), res);
    let mol = res.unwrap();
    let sg = mol.sgroups().get(&0).unwrap();
    assert_eq!(sg.connectivity_flip, expected_flip);
}

#[rustfmt::skip]
#[rstest]
#[case::sgd_with_coords(b"CC |SgD:0,1:MW:150:::::1.5,2.5,3.5,4.5|", Some(SGroupBracketCoords { bracket1: (1.5, 2.5), bracket2: (3.5, 4.5), bracket3: None, bracket4: None }))]
#[case::sgd_atom_attached(b"CC |SgD:0,1:MW:150:::::(-1)|", None)]
fn test_cx_sgroup_data_coords(
    #[case] input: &[u8],
    #[case] expected_coords: Option<SGroupBracketCoords>,
) {
    let res = parse_extended_cxsmiles(input);
    assert!(res.is_ok(), "{:?} should succeed: {:?}", input.to_str_lossy(), res);
    let mol = res.unwrap();
    let sg = mol.sgroups().get(&0).unwrap();
    assert_eq!(sg.bracket_coords, expected_coords);
}

#[rustfmt::skip]
#[rstest]
#[case::sg_atom_out_of_range(b"CC |Sg:n:0,5|", ParseError::AtomIndexOutOfBounds { atom_idx: 5 })]
#[case::sg_unknown_type(b"CC |Sg:xyz:0,1|", ParseError::InvalidCxTag { pos: 0 })]
fn test_cx_sgroup_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

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
#[case::sgh_parent_child(b"CC |Sg:n:0,1:n,SgD:0:Name:val:::::,SgH:0:1|", 1u32, Some(0u32))]
#[case::sgh_multiple_children(b"CCC |Sg:n:0,1:n,SgD:0:MW:150:::::,SgD:1,2:MW:200:::::,SgH:0:1.2|", 1u32, Some(0u32))]
fn test_cx_sgroup_hierarchy(
    #[case] input: &[u8],
    #[case] child_idx: u32,
    #[case] expected_parent: Option<u32>,
) {
    let res = parse_basic_cxsmiles(input);
    assert!(res.is_err(), "SgH is extended only");
    assert_eq!(res.unwrap_err(), ParseError::InvalidCxTag { pos: 0 });

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_ok(),
        "{:?} should succeed: {:?}",
        input.to_str_lossy(),
        res
    );
    let mol = res.unwrap();
    let sg = mol.sgroups().get(&child_idx).unwrap();
    assert_eq!(sg.hierarchy_parent, expected_parent);
}

#[rustfmt::skip]
#[rstest]
#[case::sgh_parent_out_of_range(b"CC |Sg:n:0,1,SgH:5:0|", ParseError::SgroupIndexOutOfBounds { sgroup_idx: 5 })]
#[case::sgh_child_out_of_range(b"CC |Sg:n:0,1,SgH:0:5|", ParseError::SgroupIndexOutOfBounds { sgroup_idx: 5 })]
fn test_cx_sgroup_hierarchy_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let res = parse_extended_cxsmiles(input);
    assert!(res.is_err(), "{:?} should have failed: {:?}", input.to_str_lossy(), res);
    assert_eq!(res.unwrap_err(), expected);
}

#[rustfmt::skip]
#[rstest]
#[case::sgd_atom_out_of_range(b"CC |SgD:0,5:MW:150:::::|", ParseError::AtomIndexOutOfBounds { atom_idx: 5 })]
fn test_cx_sgroup_data_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
    let input_str = input.to_str_lossy();

    let res = parse_extended_cxsmiles(input);
    assert!(
        res.is_err(),
        "{:?} should have failed: {:?}",
        input_str,
        res
    );
    assert_eq!(res.unwrap_err(), expected);
}
