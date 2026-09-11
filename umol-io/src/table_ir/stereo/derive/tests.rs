use rstest::{fixture, rstest};
use umol_geometric_core::Point3D;

use super::{derive_stereo_bonds, StereoDerivationError};
use crate::ctfile::parser::{parse_extended_mol_bytes, parse_mol_bytes_to_table_ir};
use crate::smiles::{parse_extended_smiles_bytes_with, Smiles, SmilesIoConfig};
use crate::table_ir::BondConfiguration::{Either, Framed};
use crate::table_ir::BondOrder::{Double, Single};
use crate::table_ir::BondRelation::{OppositeSide, SameSide};
use crate::table_ir::{
    AtomPair, BondConfiguration, BondOrder, BondOrientation, BondStereo, BondTaper, BondWedge,
    StereoBond,
};

#[fixture]
fn alkene() -> Vec<(AtomPair, BondOrder, Option<BondWedge>)> {
    vec![
        (AtomPair::new(0, 1), Single, None),
        (AtomPair::new(1, 2), Double, None),
        (AtomPair::new(2, 3), Single, None),
    ]
}

#[rstest]
#[case::absent(None, None)]
#[case::same_2d(Some(vec![[0.,1.,0.], [0.,0.,0.], [2.,0.,0.], [2.,1.,0.]]), Some(Framed { references: [0,3], relation: SameSide }))]
#[case::opposite_2d(Some(vec![[0.,1.,0.], [0.,0.,0.], [2.,0.,0.], [2.,-1.,0.]]), Some(Framed { references: [0,3], relation: OppositeSide }))]
#[case::same_3d(Some(vec![[0.,0.,1.], [0.,0.,0.], [0.,2.,0.], [0.,2.,1.]]), Some(Framed { references: [0,3], relation: SameSide }))]
#[case::opposite_3d(Some(vec![[0.,0.,1.], [0.,0.,0.], [0.,2.,0.], [0.,2.,-1.]]), Some(Framed { references: [0,3], relation: OppositeSide }))]
#[case::perpendicular_3d(Some(vec![[0.,1.,0.], [0.,0.,0.], [2.,0.,0.], [2.,0.,1.]]), Some(Framed { references: [0,3], relation: OppositeSide }))]
#[case::all_zero(Some(vec![[0.,0.,0.]; 4]), None)]
#[case::zero_axis(Some(vec![[0.,1.,0.], [0.,0.,0.], [0.,0.,0.], [2.,1.,0.]]), None)]
#[case::on_axis(Some(vec![[1.,0.,0.], [0.,0.,0.], [2.,0.,0.], [2.,1.,0.]]), None)]
#[case::below_tolerance(Some(vec![[0.,1e-7,0.], [0.,0.,0.], [2.,0.,0.], [2.,1.,0.]]), None)]
#[case::nan(Some(vec![[0.,f64::NAN,0.], [0.,0.,0.], [2.,0.,0.], [2.,1.,0.]]), None)]
#[case::infinite(Some(vec![[0.,f64::INFINITY,0.], [0.,0.,0.], [2.,0.,0.], [2.,1.,0.]]), None)]
#[case::large_scale(Some(vec![[0.,1e300,0.], [0.,0.,0.], [2e300,0.,0.], [2e300,-1e300,0.]]), Some(Framed { references: [0,3], relation: OppositeSide }))]
fn test_derive_stereo_bonds_geometry(
    alkene: Vec<(AtomPair, BondOrder, Option<BondWedge>)>,
    #[case] positions: Option<Vec<[f64; 3]>>,
    #[case] expected: Option<BondConfiguration>,
) {
    let positions = positions.map(|p| {
        p.into_iter()
            .map(|[x, y, z]| Point3D::new(x, y, z))
            .collect::<Vec<_>>()
    });
    assert_eq!(
        derive_stereo_bonds(4, &alkene, positions.as_deref(), &[]),
        Ok(expected
            .into_iter()
            .map(|configuration| StereoBond {
                bond: 1,
                configuration
            })
            .collect())
    );
}

#[rstest]
#[case::four(1., -1., true, Some(Framed { references: [0,3], relation: SameSide }))]
#[case::three(-1., 1., false, Some(Framed { references: [0,3], relation: OppositeSide }))]
#[case::folded(1., 1., true, None)]
#[case::complement_on_axis(1., 0., true, None)]
fn test_derive_stereo_bonds_substituents(
    #[case] reference_y: f64,
    #[case] complement_y: f64,
    #[case] fourth: bool,
    #[case] expected: Option<BondConfiguration>,
) {
    let mut bonds = vec![
        (AtomPair::new(1, 4), Double, None),
        (AtomPair::new(1, 2), Single, None),
        (AtomPair::new(0, 1), Single, None),
        (AtomPair::new(3, 4), Single, None),
    ];
    if fourth {
        bonds.push((AtomPair::new(4, 5), Single, None));
    }
    let positions = [
        [0., 1., 0.],
        [0., 0., 0.],
        [0., -1., 0.],
        [2., reference_y, 0.],
        [2., 0., 0.],
        [2., complement_y, 0.],
    ]
    .map(|[x, y, z]| Point3D::new(x, y, z));
    assert_eq!(
        derive_stereo_bonds(6, &bonds, Some(&positions), &[]),
        Ok(expected
            .into_iter()
            .map(|configuration| StereoBond {
                bond: 0,
                configuration
            })
            .collect())
    );
}

#[rstest]
#[case::cis(vec![(1, BondStereo::Cis)], None, Ok(vec![StereoBond { bond: 1, configuration: Framed { references: [0,3], relation: SameSide }}]))]
#[case::trans(vec![(1, BondStereo::Trans)], None, Ok(vec![StereoBond { bond: 1, configuration: Framed { references: [0,3], relation: OppositeSide }}]))]
#[case::either(vec![(1, BondStereo::Either)], Some(vec![]), Ok(vec![StereoBond { bond: 1, configuration: Either }]))]
#[case::duplicate(vec![(1, BondStereo::Cis), (1, BondStereo::Cis)], None, Ok(vec![StereoBond { bond: 1, configuration: Framed { references: [0,3], relation: SameSide }}]))]
#[case::geometry_agrees(vec![(1, BondStereo::Cis)], Some(vec![[0.,1.,0.],[0.,0.,0.],[2.,0.,0.],[2.,1.,0.]]), Ok(vec![StereoBond { bond: 1, configuration: Framed { references: [0,3], relation: SameSide }}]))]
#[case::degenerate_with_code(vec![(1, BondStereo::Cis)], Some(vec![[0.,0.,0.];4]), Ok(vec![StereoBond { bond: 1, configuration: Framed { references: [0,3], relation: SameSide }}]))]
#[case::either_with_geometry(vec![(1, BondStereo::Either)], Some(vec![[0.,1.,0.],[0.,0.,0.],[2.,0.,0.],[2.,1.,0.]]), Ok(vec![StereoBond { bond: 1, configuration: Either }]))]
fn test_derive_stereo_bonds_annotations(
    alkene: Vec<(AtomPair, BondOrder, Option<BondWedge>)>,
    #[case] annotations: Vec<(u32, BondStereo)>,
    #[case] positions: Option<Vec<[f64; 3]>>,
    #[case] expected: Result<Vec<StereoBond>, StereoDerivationError>,
) {
    let positions = positions.map(|p| {
        p.into_iter()
            .map(|[x, y, z]| Point3D::new(x, y, z))
            .collect::<Vec<_>>()
    });
    assert_eq!(
        derive_stereo_bonds(4, &alkene, positions.as_deref(), &annotations),
        expected
    );
}

#[rstest]
#[case::code_conflict(vec![(1, BondStereo::Cis), (1, BondStereo::Trans)], None, Err(StereoDerivationError::ConflictingConfiguration { bond: 1 }))]
#[case::either_conflict(vec![(1, BondStereo::Either), (1, BondStereo::Cis)], None, Err(StereoDerivationError::ConflictingConfiguration { bond: 1 }))]
#[case::invalid_bond(vec![(3, BondStereo::Either)], None, Err(StereoDerivationError::BondIndexOutOfBounds { bond: 3 }))]
#[case::single_bond(vec![(0, BondStereo::Cis)], None, Err(StereoDerivationError::UnsupportedSite { bond: 0 }))]
#[case::geometry_conflicts(vec![(1, BondStereo::Trans)], Some(vec![[0.,1.,0.],[0.,0.,0.],[2.,0.,0.],[2.,1.,0.]]), Err(StereoDerivationError::ConflictingConfiguration { bond: 1 }))]
fn test_derive_stereo_bonds_annotations_error(
    alkene: Vec<(AtomPair, BondOrder, Option<BondWedge>)>,
    #[case] annotations: Vec<(u32, BondStereo)>,
    #[case] positions: Option<Vec<[f64; 3]>>,
    #[case] expected: Result<Vec<StereoBond>, StereoDerivationError>,
) {
    let positions = positions.map(|p| {
        p.into_iter()
            .map(|[x, y, z]| Point3D::new(x, y, z))
            .collect::<Vec<_>>()
    });
    assert_eq!(
        derive_stereo_bonds(4, &alkene, positions.as_deref(), &annotations),
        expected
    );
}

#[rstest]
#[case::narrow(BondOrientation::Either, BondTaper::Narrowing, vec![StereoBond { bond: 1, configuration: Either }])]
#[case::wide(BondOrientation::Either, BondTaper::Widening, vec![])]
#[case::up_unknown(BondOrientation::EitherUp, BondTaper::Narrowing, vec![StereoBond { bond: 1, configuration: Either }])]
#[case::down_unknown(BondOrientation::EitherDown, BondTaper::Narrowing, vec![StereoBond { bond: 1, configuration: Either }])]
#[case::up(BondOrientation::Up, BondTaper::Narrowing, vec![])]
#[case::down(BondOrientation::Down, BondTaper::Narrowing, vec![])]
fn test_derive_stereo_bonds_wavy(
    mut alkene: Vec<(AtomPair, BondOrder, Option<BondWedge>)>,
    #[case] orientation: BondOrientation,
    #[case] taper: BondTaper,
    #[case] expected: Vec<StereoBond>,
) {
    alkene[0].2 = Some(BondWedge { orientation, taper });
    assert_eq!(derive_stereo_bonds(4, &alkene, None, &[]), Ok(expected));
}

#[rstest]
#[case::cumulene(vec![
    (AtomPair::new(0,1), Single, Some(BondWedge { orientation: BondOrientation::Either, taper: BondTaper::Narrowing })),
    (AtomPair::new(1,2), Double, None), (AtomPair::new(1,3), Double, None),
], vec![], Ok(vec![]))]
#[case::terminal(vec![(AtomPair::new(0,1), Double, None)], vec![(0,BondStereo::Either)], Ok(vec![StereoBond { bond:0, configuration:Either }]))]
#[case::wavy_redundant(vec![
    (AtomPair::new(0,1), Single, Some(BondWedge { orientation: BondOrientation::Either, taper: BondTaper::Narrowing })),
    (AtomPair::new(1,2), Double, None), (AtomPair::new(2,3), Single, None),
], vec![(1,BondStereo::Either)], Ok(vec![StereoBond { bond:1, configuration:Either }]))]
fn test_derive_stereo_bonds_context(
    #[case] bonds: Vec<(AtomPair, BondOrder, Option<BondWedge>)>,
    #[case] annotations: Vec<(u32, BondStereo)>,
    #[case] expected: Result<Vec<StereoBond>, StereoDerivationError>,
) {
    assert_eq!(derive_stereo_bonds(4, &bonds, None, &annotations), expected);
}

#[rstest]
#[case::terminal_definite(vec![(AtomPair::new(0,1), Double, None)], vec![(0,BondStereo::Cis)], Err(StereoDerivationError::UnsupportedSite { bond:0 }))]
#[case::wavy_conflict(vec![
    (AtomPair::new(0,1), Single, Some(BondWedge { orientation: BondOrientation::Either, taper: BondTaper::Narrowing })),
    (AtomPair::new(1,2), Double, None), (AtomPair::new(2,3), Single, None),
], vec![(1,BondStereo::Cis)], Err(StereoDerivationError::ConflictingConfiguration { bond:1 }))]
fn test_derive_stereo_bonds_context_error(
    #[case] bonds: Vec<(AtomPair, BondOrder, Option<BondWedge>)>,
    #[case] annotations: Vec<(u32, BondStereo)>,
    #[case] expected: Result<Vec<StereoBond>, StereoDerivationError>,
) {
    assert_eq!(derive_stereo_bonds(4, &bonds, None, &annotations), expected);
}

#[rstest]
fn test_derive_stereo_bonds_positions_error(alkene: Vec<(AtomPair, BondOrder, Option<BondWedge>)>) {
    let points = [
        Point3D::new(0., 1., 0.),
        Point3D::zero(),
        Point3D::new(2., 0., 0.),
    ];
    assert_eq!(
        derive_stereo_bonds(4, &alkene, Some(&points), &[]),
        Err(StereoDerivationError::MissingPosition { atom: 3 })
    );
}

#[rstest]
fn test_derive_stereo_bonds_index_error() {
    assert_eq!(
        derive_stereo_bonds(2, &[(AtomPair::new(0, 2), Double, None)], Some(&[]), &[]),
        Err(StereoDerivationError::AtomIndexOutOfBounds { atom: 2 })
    );
}

#[rstest]
#[case::no_references("C=C |ctu:0|", vec![StereoBond { bond: 0, configuration: Either }])]
#[case::acyclic_cis("CC(F)=C(C)F |c:2|", vec![StereoBond { bond: 2, configuration: Framed { references: [0,4], relation: SameSide }}])]
#[case::acyclic_trans("CC(F)=C(C)F |t:2|", vec![StereoBond { bond: 2, configuration: Framed { references: [0,4], relation: OppositeSide }}])]
#[case::ring("FC1(=C(F)CCCCCCCCCC1) |c:1|", vec![StereoBond { bond: 2, configuration: Framed { references: [0,3], relation: SameSide }}])]
#[case::closure_site("C1CCCCCCCCCC=1 |t:10|", vec![StereoBond { bond: 0, configuration: Framed { references: [1,9], relation: OppositeSide }}])]
#[case::wavy("CC=CC |w:1.0|", vec![StereoBond { bond: 1, configuration: Either }])]
fn test_derive_stereo_bonds_cx(
    #[case] input: &str,
    #[case] expected: Vec<StereoBond>,
    #[values(false, true)] extended: bool,
) {
    let config = SmilesIoConfig::chemaxon();
    let (count, bonds, positions, annotations) = if extended {
        let mol = parse_extended_smiles_bytes_with(input.as_bytes(), &config).unwrap();
        let annotations = mol
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(i, b)| b.stereo.map(|s| (i as u32, s)))
            .collect::<Vec<_>>();
        (
            mol.atoms.len(),
            mol.bonds
                .iter()
                .map(|b| (b.atoms, b.order, b.wedge))
                .collect::<Vec<_>>(),
            mol.positions,
            annotations,
        )
    } else {
        let mol = Smiles::parse_with(input, &config).unwrap().into_table_ir();
        let annotations = mol
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(i, b)| b.stereo.map(|s| (i as u32, s)))
            .collect::<Vec<_>>();
        (
            mol.atoms.len(),
            mol.bonds
                .iter()
                .map(|b| (b.atoms, b.order, b.wedge))
                .collect::<Vec<_>>(),
            mol.positions,
            annotations,
        )
    };
    assert_eq!(
        derive_stereo_bonds(count, &bonds, positions.as_deref(), &annotations),
        Ok(expected)
    );
}

#[rstest]
#[case::same(1., false, Framed { references: [0,3], relation: SameSide })]
#[case::opposite(-1., false, Framed { references: [0,3], relation: OppositeSide })]
#[case::either(1., true, Either)]
fn test_derive_stereo_bonds_ctfile(
    #[case] y: f64,
    #[case] either: bool,
    #[case] configuration: BondConfiguration,
    #[values(false, true)] reversed: bool,
    #[values(false, true)] extended: bool,
) {
    let (a, b) = if reversed { (3, 2) } else { (2, 3) };
    let input = {
        let code = if either { 3 } else { 0 };
        format!("\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    1.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000{y:10.4}    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n{a:3}{b:3}  2{code:3}        0\n  3  4  1  0        0\nM  END\n")
    };
    let (count, bonds, positions, annotations) = if extended {
        let mol = parse_extended_mol_bytes(input.as_bytes()).unwrap();
        let annotations = mol
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(i, b)| b.stereo.map(|s| (i as u32, s)))
            .collect::<Vec<_>>();
        (
            mol.atoms.len(),
            mol.bonds
                .iter()
                .map(|b| (b.atoms, b.order, b.wedge))
                .collect::<Vec<_>>(),
            mol.positions,
            annotations,
        )
    } else {
        let mol = parse_mol_bytes_to_table_ir(input.as_bytes()).unwrap();
        let annotations = mol
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(i, b)| b.stereo.map(|s| (i as u32, s)))
            .collect::<Vec<_>>();
        (
            mol.atoms.len(),
            mol.bonds
                .iter()
                .map(|b| (b.atoms, b.order, b.wedge))
                .collect::<Vec<_>>(),
            mol.positions,
            annotations,
        )
    };
    assert_eq!(
        derive_stereo_bonds(count, &bonds, positions.as_deref(), &annotations),
        Ok(vec![StereoBond {
            bond: 1,
            configuration
        }])
    );
}
