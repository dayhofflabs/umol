use rstest::rstest;
use umol_geometric_core::Point2D;
#[cfg(feature = "coordgen")]
use umol_graph_ir::ir::Molecule;
use umol_graph_ir::ir::{AtomId, BondId};
use umol_graph_ir::mol;
#[cfg(feature = "coordgen")]
use umol_graph_ir::mol_dsl;
#[cfg(feature = "coordgen")]
use umol_io::layout::{layout_molecule, MoleculeLayoutAlgorithm};
use umol_io::layout::{MoleculeLayout, MoleculeLayoutError};

#[cfg(feature = "coordgen")]
#[rstest]
#[case::coordgen(MoleculeLayoutAlgorithm::CoordGen)]
fn test_layout_molecule_frame(#[case] algorithm: MoleculeLayoutAlgorithm) {
    let molecule = mol! {
        (oxygen: O), (nitrogen: N), (fluorine: F), (carbon: C),
        (carbon) - (oxygen),
        (carbon) = (nitrogen),
        (carbon) - (fluorine),
    };
    let layout = layout_molecule(&molecule, algorithm).expect("layout succeeds");

    assert_eq!(layout.atom_count(), molecule.atoms().count());
    for atom_id in molecule.atoms().ids() {
        assert!(layout.position(atom_id).is_some());
    }
    for terminal in [AtomId(0), AtomId(1), AtomId(2)] {
        let distance = layout_distance(&layout, AtomId(3), terminal);
        assert!((distance - 1.0).abs() < 1e-3, "bond length was {distance}");
    }
    assert!(layout_distance(&layout, AtomId(0), AtomId(1)) > 1.0);
}

#[cfg(feature = "coordgen")]
#[rstest]
#[case::coordgen(MoleculeLayoutAlgorithm::CoordGen)]
fn test_layout_molecule_empty_frame(#[case] algorithm: MoleculeLayoutAlgorithm) {
    let layout = layout_molecule(&Molecule::new(), algorithm).expect("layout succeeds");

    assert_eq!(layout.positions(), &[]);
}

#[cfg(feature = "coordgen")]
#[rstest]
#[case::coordgen(MoleculeLayoutAlgorithm::CoordGen)]
fn test_layout_molecule_scale(#[case] algorithm: MoleculeLayoutAlgorithm) {
    let molecule = mol! {
        (carbon_0: C) - (carbon_1: C),
        (carbon_1) - (oxygen: O),
    };
    let layout = layout_molecule(&molecule, algorithm).expect("layout succeeds");

    assert!(layout
        .positions()
        .iter()
        .all(|position| position.x.is_finite() && position.y.is_finite()));
    for bond in molecule.bonds().iter() {
        let [atom_0, atom_1] = bond.atom_ids();
        let distance = layout_distance(&layout, atom_0, atom_1);
        assert!((distance - 1.0).abs() < 1e-3, "bond length was {distance}");
    }
}

#[cfg(feature = "coordgen")]
#[rstest]
#[case::coordgen(MoleculeLayoutAlgorithm::CoordGen)]
fn test_layout_molecule_projection(#[case] algorithm: MoleculeLayoutAlgorithm) {
    let molecule = mol! {
        (generic: "*") -[ "*" ]- (oxygen: O),
    };

    let layout = layout_molecule(&molecule, algorithm).expect("layout succeeds");

    assert_eq!(layout.atom_count(), 2);
    assert!(layout
        .positions()
        .iter()
        .all(|point| point.x.is_finite() && point.y.is_finite()));
    let distance = layout_distance(&layout, AtomId(0), AtomId(1));
    assert!((distance - 1.0).abs() < 1e-3, "bond length was {distance}");
}

#[cfg(feature = "coordgen")]
#[rstest]
#[case::coordgen(MoleculeLayoutAlgorithm::CoordGen)]
fn test_layout_molecule_determinism(#[case] algorithm: MoleculeLayoutAlgorithm) {
    let molecule = mol! {
        (carbon_0: C) - (carbon_1: C) = (oxygen: O),
        (carbon_1) - (nitrogen: N),
    };

    let first = layout_molecule(&molecule, algorithm).expect("first layout succeeds");
    let second = layout_molecule(&molecule, algorithm).expect("second layout succeeds");

    assert_eq!(first, second);
}

#[cfg(feature = "coordgen")]
#[rstest]
#[case::z_implicit_hydrogen(
    r#"{:atoms ["C" "C" "C" "C"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
        :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]}"#,
    BondId(1),
    AtomId(0),
    AtomId(3),
    true
)]
#[case::e_implicit_hydrogen(
    r#"{:atoms ["C" "C" "C" "C"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
        :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#,
    BondId(1),
    AtomId(0),
    AtomId(3),
    false
)]
#[case::reversed_endpoint_frame(
    r#"{:atoms ["C" "C" "C" "C"]
        :bonds [[0 1 "1"] [2 1 "2"] [2 3 "1"]]
        :stereo-bonds [{:site 1 :ligands [3 [:h 2] 0 [:h 1]] :attrs "Ct0"}]}"#,
    BondId(1),
    AtomId(3),
    AtomId(0),
    true
)]
#[case::reframed_implicit_hydrogen(
    r#"{:atoms ["C" "C" "C" "C"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
        :stereo-bonds [{:site 1 :ligands [[:h 1] 0 3 [:h 2]] :attrs "Ct0"}]}"#,
    BondId(1),
    AtomId(0),
    AtomId(3),
    false
)]
fn test_layout_molecule_cis_trans(
    #[case] input: &str,
    #[case] site: BondId,
    #[case] first_ligand: AtomId,
    #[case] second_ligand: AtomId,
    #[case] expected_same_side: bool,
) {
    let molecule = mol_dsl!(input);

    let layout = layout_molecule(&molecule, MoleculeLayoutAlgorithm::CoordGen)
        .expect("cis/trans layout succeeds");
    let [site_0, site_1] = molecule.bond(site).atom_ids();

    assert_eq!(
        layout_same_side(&layout, site_0, site_1, first_ligand, second_ligand),
        expected_same_side
    );
}

#[rstest]
#[case::nan(Point2D::new(0.0, f64::NAN))]
#[case::positive_infinity(Point2D::new(f64::INFINITY, 0.0))]
#[case::negative_infinity(Point2D::new(0.0, f64::NEG_INFINITY))]
fn test_molecule_layout_try_new(#[case] position: Point2D) {
    let error = MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), position])
        .expect_err("non-finite position must be rejected");

    assert!(matches!(
        error,
        MoleculeLayoutError::NonFinitePosition {
            atom_id: AtomId(1),
            position: actual,
        } if actual.x.to_bits() == position.x.to_bits()
            && actual.y.to_bits() == position.y.to_bits()
    ));
}

#[rstest]
fn test_molecule_layout_accessors() {
    let layout = MoleculeLayout::try_new(Vec::new()).expect("empty layout is valid");

    assert_eq!(layout.atom_count(), 0);
    assert!(layout.is_empty());
    assert!(layout.positions().is_empty());
    assert_eq!(layout.position(AtomId(0)), None);
}

#[rstest]
fn test_molecule_layout_set_position() {
    let initial = Point2D::new(1.0, 2.0);
    let updated = Point2D::new(-3.0, 4.0);
    let mut layout = MoleculeLayout::try_new(vec![initial]).expect("finite layout");

    layout
        .set_position(AtomId(0), updated)
        .expect("in-frame finite edit");
    assert_eq!(layout.position(AtomId(0)), Some(&updated));

    let nonfinite = Point2D::new(f64::NAN, 5.0);
    let error = layout
        .set_position(AtomId(0), nonfinite)
        .expect_err("non-finite edit must be rejected");
    assert!(matches!(
        error,
        MoleculeLayoutError::NonFinitePosition {
            atom_id: AtomId(0),
            position,
        } if position.x.is_nan() && position.y == 5.0
    ));
    assert_eq!(layout.position(AtomId(0)), Some(&updated));

    assert_eq!(
        layout.set_position(AtomId(1), Point2D::new(0.0, 0.0)),
        Err(MoleculeLayoutError::AtomOutOfFrame {
            atom_id: AtomId(1),
            frame_size: 1,
        })
    );
    assert_eq!(layout.positions(), &[updated]);
}

#[rstest]
fn test_molecule_layout_check_frame() {
    let molecule = mol! {
        (carbon: C), (oxygen: O),
    };
    let matching = MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)])
        .expect("finite layout");
    let mismatched = MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).expect("finite layout");

    assert_eq!(matching.check_frame(&molecule), Ok(()));
    assert_eq!(
        mismatched.check_frame(&molecule),
        Err(MoleculeLayoutError::FrameSizeMismatch {
            molecule_atom_count: 2,
            layout_atom_count: 1,
        })
    );
}

#[cfg(feature = "coordgen")]
fn layout_distance(layout: &MoleculeLayout, atom_0: AtomId, atom_1: AtomId) -> f64 {
    let point_0 = layout
        .position(atom_0)
        .expect("first atom is in the layout frame");
    let point_1 = layout
        .position(atom_1)
        .expect("second atom is in the layout frame");
    (point_1.x - point_0.x).hypot(point_1.y - point_0.y)
}

#[cfg(feature = "coordgen")]
fn layout_same_side(
    layout: &MoleculeLayout,
    site_0: AtomId,
    site_1: AtomId,
    first_ligand: AtomId,
    second_ligand: AtomId,
) -> bool {
    let first = layout_half_plane(layout, site_0, site_1, first_ligand);
    let second = layout_half_plane(layout, site_0, site_1, second_ligand);
    assert!(first.abs() > 1e-6);
    assert!(second.abs() > 1e-6);
    first.is_sign_positive() == second.is_sign_positive()
}

#[cfg(feature = "coordgen")]
fn layout_half_plane(
    layout: &MoleculeLayout,
    site_0: AtomId,
    site_1: AtomId,
    ligand: AtomId,
) -> f64 {
    let site_0 = layout.position(site_0).expect("site atom is in frame");
    let site_1 = layout.position(site_1).expect("site atom is in frame");
    let ligand = layout.position(ligand).expect("ligand atom is in frame");
    (site_1.x - site_0.x) * (ligand.y - site_0.y) - (site_1.y - site_0.y) * (ligand.x - site_0.x)
}
