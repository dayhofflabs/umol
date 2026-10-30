use rstest::rstest;

use super::*;

#[rstest]
#[case::atom(GraphStereoSite::Atom(GraphIrAtomId(3)), StereoSite::Atom(3))]
#[case::bond(GraphStereoSite::Bond(GraphIrBondId(2)), StereoSite::Bond(2))]
fn test_stereo_site_from_rust(#[case] value: GraphStereoSite, #[case] expected: StereoSite) {
    assert_eq!(StereoSite::from_rust(value), expected);
    assert_eq!(expected.to_rust(), value);
}

#[rstest]
#[case::tetrahedral(GraphStereoInconsistency::TetrahedralStereoFailure { atom: GraphIrAtomId(3) }, StereoInconsistency::TetrahedralStereoFailure { atom: 3 })]
#[case::atom(GraphStereoInconsistency::StereoAtomFailure { stereo_atom: GraphIrStereoAtomId(2) }, StereoInconsistency::StereoAtomFailure { stereo_atom: 2 })]
#[case::atom_mismatch(GraphStereoInconsistency::TetrahedralStereoMismatch { atom: GraphIrAtomId(3), stereo_atom: GraphIrStereoAtomId(2) }, StereoInconsistency::TetrahedralStereoMismatch { atom: 3, stereo_atom: 2 })]
#[case::cis_trans(GraphStereoInconsistency::CisTransStereoFailure { bond: GraphIrBondId(3) }, StereoInconsistency::CisTransStereoFailure { bond: 3 })]
#[case::bond(GraphStereoInconsistency::StereoBondFailure { stereo_bond: GraphIrStereoBondId(2) }, StereoInconsistency::StereoBondFailure { stereo_bond: 2 })]
#[case::bond_mismatch(GraphStereoInconsistency::CisTransStereoMismatch { bond: GraphIrBondId(3), stereo_bond: GraphIrStereoBondId(2) }, StereoInconsistency::CisTransStereoMismatch { bond: 3, stereo_bond: 2 })]
fn test_stereo_inconsistency_from_rust(
    #[case] value: GraphStereoInconsistency,
    #[case] expected: StereoInconsistency,
) {
    assert_eq!(StereoInconsistency::from_rust(value), expected);
    assert_eq!(expected.to_rust(), value);
}
