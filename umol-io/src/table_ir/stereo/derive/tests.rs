use std::cell::OnceCell;

use rstest::{fixture, rstest};
use umol_geometric_core::Point3D;

use super::{derive_stereo_bonds, StereoDerivationError};
use crate::ctfile::config::{CtabParseFlags, CtfileIoConfig};
use crate::ctfile::error::ParseError as CtfileParseError;
use crate::ctfile::parser::{
    parse_extended_mol_bytes, parse_extended_mol_bytes_with, parse_mol_bytes_to_table_ir,
    parse_mol_bytes_to_table_ir_with,
};
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
        derive_stereo_bonds(
            4,
            &alkene,
            |bond| *bond,
            positions.as_deref(),
            Vec::new(),
            Vec::new(),
            &OnceCell::new()
        ),
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
        derive_stereo_bonds(
            6,
            &bonds,
            |bond| *bond,
            Some(&positions),
            Vec::new(),
            Vec::new(),
            &OnceCell::new()
        ),
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
        derive_stereo_bonds(
            4,
            &alkene,
            |bond| *bond,
            positions.as_deref(),
            Vec::new(),
            annotations,
            &OnceCell::new()
        ),
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
        derive_stereo_bonds(
            4,
            &alkene,
            |bond| *bond,
            positions.as_deref(),
            Vec::new(),
            annotations,
            &OnceCell::new()
        ),
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
    assert_eq!(
        derive_stereo_bonds(
            4,
            &alkene,
            |bond| *bond,
            None,
            Vec::new(),
            Vec::new(),
            &OnceCell::new()
        ),
        Ok(expected)
    );
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
    assert_eq!(
        derive_stereo_bonds(
            4,
            &bonds,
            |bond| *bond,
            None,
            Vec::new(),
            annotations,
            &OnceCell::new()
        ),
        expected
    );
}

#[rstest]
#[case::cumulated_definite(vec![
    (AtomPair::new(0,1), Double, None),
    (AtomPair::new(0,2), Double, None),
    (AtomPair::new(1,3), Single, None),
], vec![(0,BondStereo::Cis)], Err(StereoDerivationError::UnsupportedSite { bond:0 }))]
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
    assert_eq!(
        derive_stereo_bonds(
            4,
            &bonds,
            |bond| *bond,
            None,
            Vec::new(),
            annotations,
            &OnceCell::new()
        ),
        expected
    );
}

#[rstest]
#[case::one(&[2], &[4], true)]
#[case::two(&[3,2], &[5,4], true)]
#[case::duplicates(&[3,3,2,3,2], &[5,4,5,4], true)]
#[case::zero(&[], &[4], false)]
#[case::excess(&[2,3,5], &[4], false)]
#[case::second_excess(&[2], &[3,4,5], false)]
#[case::shared(&[2,3], &[3,4], false)]
fn test_derive_stereo_bonds_incidences(
    #[case] first: &[u32],
    #[case] second: &[u32],
    #[case] supported: bool,
    #[values(false, true)] explicit: bool,
) {
    let mut bonds = vec![(AtomPair::new(0, 1), Double, None)];
    for (endpoint, ligands) in [(0, first), (1, second)] {
        bonds.extend(
            ligands
                .iter()
                .map(|&atom| (AtomPair::new(endpoint, atom), Single, None)),
        );
    }
    let positions = [
        Point3D::new(0., 0., 0.),
        Point3D::new(2., 0., 0.),
        Point3D::new(0., 1., 0.),
        Point3D::new(0., -1., 0.),
        Point3D::new(2., 1., 0.),
        Point3D::new(2., -1., 0.),
    ];
    let expected = if supported {
        Ok(vec![StereoBond {
            bond: 0,
            configuration: Framed {
                references: [2, 4],
                relation: SameSide,
            },
        }])
    } else if explicit {
        Err(StereoDerivationError::UnsupportedSite { bond: 0 })
    } else {
        Ok(vec![])
    };
    assert_eq!(
        derive_stereo_bonds(
            6,
            &bonds,
            |bond| *bond,
            Some(&positions),
            vec![],
            if explicit {
                vec![(0, BondStereo::Cis)]
            } else {
                vec![]
            },
            &OnceCell::new()
        ),
        expected
    );
}

#[rstest]
fn test_derive_stereo_bonds_positions_error(alkene: Vec<(AtomPair, BondOrder, Option<BondWedge>)>) {
    let points = [
        Point3D::new(0., 1., 0.),
        Point3D::zero(),
        Point3D::new(2., 0., 0.),
    ];
    assert_eq!(
        derive_stereo_bonds(
            4,
            &alkene,
            |bond| *bond,
            Some(&points),
            Vec::new(),
            Vec::new(),
            &OnceCell::new()
        ),
        Err(StereoDerivationError::MissingPosition { atom: 3 })
    );
}

#[rstest]
fn test_derive_stereo_bonds_index_error() {
    assert_eq!(
        derive_stereo_bonds(
            2,
            &[(AtomPair::new(0, 2), Double, None)],
            |bond| *bond,
            Some(&[]),
            Vec::new(),
            Vec::new(),
            &OnceCell::new()
        ),
        Err(StereoDerivationError::AtomIndexOutOfBounds { atom: 2 })
    );
}

#[rstest]
#[case::minimum(Framed { references: [0,3], relation: SameSide }, BondStereo::Cis)]
#[case::first_complement(Framed { references: [2,3], relation: OppositeSide }, BondStereo::Cis)]
#[case::second_complement(Framed { references: [0,5], relation: OppositeSide }, BondStereo::Cis)]
#[case::both_complements(Framed { references: [2,5], relation: SameSide }, BondStereo::Cis)]
#[case::either(Either, BondStereo::Either)]
fn test_derive_stereo_bonds_frames_identity(
    #[case] configuration: BondConfiguration,
    #[case] code: BondStereo,
    #[values(false, true)] geometry: bool,
) {
    let bonds = [
        (AtomPair::new(1, 4), Double, None),
        (AtomPair::new(0, 1), Single, None),
        (AtomPair::new(1, 2), Single, None),
        (AtomPair::new(3, 4), Single, None),
        (AtomPair::new(4, 5), Single, None),
    ];
    let positions = [
        [0., 1., 0.],
        [0., 0., 0.],
        [0., -1., 0.],
        [2., 1., 0.],
        [2., 0., 0.],
        [2., -1., 0.],
    ]
    .map(|[x, y, z]| Point3D::new(x, y, z));
    let frames = vec![StereoBond {
        bond: 0,
        configuration,
    }];
    let expected = frames.clone();
    let allocation = frames.as_ptr();
    let actual = derive_stereo_bonds(
        6,
        &bonds,
        |bond| *bond,
        geometry.then_some(positions.as_slice()),
        frames,
        vec![(0, code), (0, code)],
        &OnceCell::new(),
    )
    .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(actual.as_ptr(), allocation);
}

#[rstest]
#[case::complement_conflict(Framed { references: [2,3], relation: OppositeSide }, BondStereo::Trans, StereoDerivationError::ConflictingConfiguration { bond: 0 })]
#[case::either_conflict(Either, BondStereo::Cis, StereoDerivationError::ConflictingConfiguration { bond: 0 })]
#[case::definite_either_conflict(Framed { references: [0,3], relation: SameSide }, BondStereo::Either, StereoDerivationError::ConflictingConfiguration { bond: 0 })]
#[case::invalid_reference(Framed { references: [4,3], relation: SameSide }, BondStereo::Cis, StereoDerivationError::UnsupportedSite { bond: 0 })]
fn test_derive_stereo_bonds_frames_error(
    #[case] configuration: BondConfiguration,
    #[case] code: BondStereo,
    #[case] expected: StereoDerivationError,
) {
    let bonds = [
        (AtomPair::new(1, 4), Double, None),
        (AtomPair::new(0, 1), Single, None),
        (AtomPair::new(1, 2), Single, None),
        (AtomPair::new(3, 4), Single, None),
        (AtomPair::new(4, 5), Single, None),
    ];
    assert_eq!(
        derive_stereo_bonds(
            6,
            &bonds,
            |bond| *bond,
            None,
            vec![StereoBond {
                bond: 0,
                configuration
            }],
            vec![(0, code)],
            &OnceCell::new()
        ),
        Err(expected)
    );
}

#[rstest]
#[case::no_references("C=C |ctu:0|", vec![StereoBond { bond: 0, configuration: Either }])]
#[case::acyclic_cis("CC(F)=C(C)F |c:2|", vec![StereoBond { bond: 2, configuration: Framed { references: [0,4], relation: SameSide }}])]
#[case::acyclic_trans("CC(F)=C(C)F |t:2|", vec![StereoBond { bond: 2, configuration: Framed { references: [0,4], relation: OppositeSide }}])]
#[case::ring("FC1(=C(F)CCCCCCCCCC1) |c:1|", vec![StereoBond { bond: 2, configuration: Framed { references: [0,3], relation: SameSide }}])]
#[case::closure_site("C1CCCCCCCCCC=1 |t:10|", vec![StereoBond { bond: 0, configuration: Framed { references: [1,9], relation: OppositeSide }}])]
#[case::wavy("CC=CC |w:1.0|", vec![StereoBond { bond: 1, configuration: Either }])]
fn test_parse_molecule_cx_stereo_bonds(
    #[case] input: &str,
    #[case] expected: Vec<StereoBond>,
    #[values(false, true)] extended: bool,
) {
    let config = SmilesIoConfig::chemaxon();
    let actual = if extended {
        parse_extended_smiles_bytes_with(input.as_bytes(), &config)
            .unwrap()
            .stereo_bonds
    } else {
        Smiles::parse_with(input, &config)
            .unwrap()
            .into_table_ir()
            .stereo_bonds
    };
    assert_eq!(actual, expected);
}

#[rstest]
#[case::same(1., false, Framed { references: [0,3], relation: SameSide })]
#[case::opposite(-1., false, Framed { references: [0,3], relation: OppositeSide })]
#[case::either(1., true, Either)]
fn test_parse_mol_bytes_to_table_ir_stereo_bonds(
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
    let actual = if extended {
        parse_extended_mol_bytes(input.as_bytes())
            .unwrap()
            .stereo_bonds
    } else {
        parse_mol_bytes_to_table_ir(input.as_bytes())
            .unwrap()
            .stereo_bonds
    };
    assert_eq!(
        actual,
        vec![StereoBond {
            bond: 1,
            configuration
        }]
    );
}

#[rstest]
#[case::promote(1, 0, "M  ZBO  1   2   2\n", Double, vec![StereoBond { bond: 1, configuration: Framed { references: [0,3], relation: SameSide }}])]
#[case::demote(2, 0, "M  ZBO  1   2   1\n", Single, vec![])]
#[case::retain_either(2, 3, "M  ZBO  1   2   2\n", Double, vec![StereoBond { bond: 1, configuration: Either }])]
fn test_parse_mol_bytes_to_table_ir_stereo_properties(
    #[case] order: u8,
    #[case] code: u8,
    #[case] property: &str,
    #[case] expected_order: BondOrder,
    #[case] expected_frames: Vec<StereoBond>,
    #[values(false, true)] extended: bool,
) {
    let input = format!("\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1\n  2  3{order:3}{code:3}\n  3  4  1\n{property}M  CHG  1   1  -1\nM  ISO  1   4  13\nM  END\n");
    let mut config = if extended {
        CtfileIoConfig::extended()
    } else {
        CtfileIoConfig::basic()
    };
    config.parse_flags |= CtabParseFlags::CLARK_EXTENSIONS;
    let actual = if extended {
        let mol = parse_extended_mol_bytes_with(input.as_bytes(), &config).unwrap();
        (
            mol.bonds[1].order,
            mol.stereo_bonds,
            mol.atoms[0].charge,
            mol.atoms[3].isotope_mass,
        )
    } else {
        let mol = parse_mol_bytes_to_table_ir_with(input.as_bytes(), &config).unwrap();
        (
            mol.bonds[1].order,
            mol.stereo_bonds,
            mol.atoms[0].charge,
            mol.atoms[3].isotope_mass,
        )
    };
    assert_eq!(
        actual,
        (expected_order, expected_frames, Some(-1), Some(13))
    );
}

#[rstest]
fn test_parse_mol_bytes_to_table_ir_stereo_properties_error(#[values(false, true)] extended: bool) {
    let input = b"\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  3\nM  ZBO  1   1   1\nM  END\n";
    let mut config = if extended {
        CtfileIoConfig::extended()
    } else {
        CtfileIoConfig::basic()
    };
    config.parse_flags |= CtabParseFlags::CLARK_EXTENSIONS;
    let actual = if extended {
        parse_extended_mol_bytes_with(input, &config).map(|_| ())
    } else {
        parse_mol_bytes_to_table_ir_with(input, &config).map(|_| ())
    };
    assert_eq!(
        actual,
        Err(CtfileParseError::UnsupportedStereoBond { bond: 0 })
    );
}

#[rstest]
#[case::none(vec![], vec![])]
#[case::either(vec![(1, BondStereo::Either)], vec![StereoBond { bond: 1, configuration: Either }])]
fn test_derive_stereo_bonds_lookup(
    alkene: Vec<(AtomPair, BondOrder, Option<BondWedge>)>,
    #[case] assertions: Vec<(u32, BondStereo)>,
    #[case] expected: Vec<StereoBond>,
) {
    let neighbors = OnceCell::new();
    assert_eq!(
        derive_stereo_bonds(
            4,
            &alkene,
            |bond| *bond,
            None,
            vec![],
            assertions,
            &neighbors
        ),
        Ok(expected)
    );
    assert_eq!(neighbors.get(), None);
}
