use rstest::rstest;
use umol_geometric_core::Point2D;
use umol_graph_ir::ir::AtomId;
#[cfg(feature = "coordgen")]
use umol_graph_ir::ir::{BondId, Molecule};
use umol_graph_ir::mol;
#[cfg(feature = "coordgen")]
use umol_graph_ir::mol_dsl;
#[cfg(feature = "coordgen")]
use umol_io::depict::{Depict, DepictConfig};
#[cfg(feature = "coordgen")]
use umol_io::layout::MoleculeLayoutAlgorithm;
use umol_io::layout::{MoleculeLayout, MoleculeLayoutError, ReactionLayout, ReactionLayoutError};

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
    let layout = molecule
        .layout_with(&DepictConfig {
            layout_algorithm: algorithm,
        })
        .expect("layout succeeds");

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
    let layout = Molecule::new()
        .layout_with(&DepictConfig {
            layout_algorithm: algorithm,
        })
        .expect("layout succeeds");

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
    let layout = molecule
        .layout_with(&DepictConfig {
            layout_algorithm: algorithm,
        })
        .expect("layout succeeds");

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

    let layout = molecule
        .layout_with(&DepictConfig {
            layout_algorithm: algorithm,
        })
        .expect("layout succeeds");

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

    let first = molecule
        .layout_with(&DepictConfig {
            layout_algorithm: algorithm,
        })
        .expect("first layout succeeds");
    let second = molecule
        .layout_with(&DepictConfig {
            layout_algorithm: algorithm,
        })
        .expect("second layout succeeds");

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

    let layout = molecule
        .layout_with(&DepictConfig {
            layout_algorithm: MoleculeLayoutAlgorithm::CoordGen,
        })
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

#[rstest]
fn test_reaction_layout_arrange() {
    let lhs = MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)])
        .expect("finite layout");
    let rhs = MoleculeLayout::try_new(vec![Point2D::new(0.0, 1.0), Point2D::new(1.0, 3.0)])
        .expect("finite layout");

    let layout = ReactionLayout::arrange(lhs, rhs).expect("finite sides arrange");

    assert_eq!(
        layout.lhs().positions(),
        [Point2D::new(-2.75, 0.0), Point2D::new(-1.75, 0.0)]
    );
    assert_eq!(
        layout.rhs().positions(),
        [Point2D::new(1.75, -1.0), Point2D::new(2.75, 1.0)]
    );
    assert_eq!(layout.arrow_start(), Point2D::new(-0.75, 0.0));
    assert_eq!(layout.arrow_end(), Point2D::new(0.75, 0.0));
}

#[rstest]
fn test_reaction_layout_arrange_empty_sides() {
    let empty = MoleculeLayout::try_new(Vec::new()).expect("empty layout is valid");

    let layout =
        ReactionLayout::arrange(empty.clone(), empty.clone()).expect("empty sides arrange");

    assert_eq!(layout.lhs(), &empty);
    assert_eq!(layout.rhs(), &empty);
    assert_eq!(layout.arrow_start(), Point2D::new(-0.75, 0.0));
    assert_eq!(layout.arrow_end(), Point2D::new(0.75, 0.0));
}

#[rstest]
#[case::lhs(true)]
#[case::rhs(false)]
fn test_reaction_layout_arrange_overflow(#[case] lhs: bool) {
    let far = MoleculeLayout::try_new(vec![Point2D::new(0.0, f64::MAX)]).expect("finite layout");
    let near = MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).expect("finite layout");
    let (lhs_layout, rhs_layout) = if lhs { (far, near) } else { (near, far) };

    let error = ReactionLayout::arrange(lhs_layout, rhs_layout).expect_err("translation overflows");

    let translation = match error {
        ReactionLayoutError::LhsTranslation(source) if lhs => source,
        ReactionLayoutError::RhsTranslation(source) if !lhs => source,
        other => panic!("unexpected error {other:?}"),
    };
    assert!(matches!(
        translation,
        MoleculeLayoutError::NonFinitePosition {
            atom_id: AtomId(0),
            ..
        }
    ));
}

#[rstest]
fn test_reaction_layout_try_new() {
    let lhs = MoleculeLayout::try_new(vec![Point2D::new(-2.0, 0.0)]).expect("finite layout");
    let rhs = MoleculeLayout::try_new(vec![Point2D::new(2.0, 0.0)]).expect("finite layout");
    let start = Point2D::new(-1.0, 0.5);
    let end = Point2D::new(1.0, 0.5);

    let layout = ReactionLayout::try_new(lhs.clone(), rhs.clone(), start, end)
        .expect("finite distinct arrow endpoints");

    assert_eq!(layout.lhs(), &lhs);
    assert_eq!(layout.rhs(), &rhs);
    assert_eq!(layout.arrow_start(), start);
    assert_eq!(layout.arrow_end(), end);
}

#[rstest]
#[case::nan_start(Point2D::new(f64::NAN, 0.0), Point2D::new(1.0, 0.0))]
#[case::infinite_end(Point2D::new(-1.0, 0.0), Point2D::new(f64::INFINITY, 0.0))]
fn test_reaction_layout_non_finite_arrow(#[case] start: Point2D, #[case] end: Point2D) {
    let side = MoleculeLayout::try_new(Vec::new()).expect("empty layout is valid");

    let error = ReactionLayout::try_new(side.clone(), side.clone(), start, end)
        .expect_err("non-finite arrow must be rejected");

    assert!(matches!(
        error,
        ReactionLayoutError::NonFiniteArrow { start: actual_start, end: actual_end }
            if actual_start.x.to_bits() == start.x.to_bits()
                && actual_end.x.to_bits() == end.x.to_bits()
    ));
}

#[rstest]
fn test_reaction_layout_degenerate_arrow() {
    let side = MoleculeLayout::try_new(Vec::new()).expect("empty layout is valid");
    let position = Point2D::new(0.5, -0.5);

    assert_eq!(
        ReactionLayout::try_new(side.clone(), side, position, position),
        Err(ReactionLayoutError::DegenerateArrow { position })
    );
}

#[rstest]
fn test_reaction_layout_set_arrow() {
    let side = MoleculeLayout::try_new(Vec::new()).expect("empty layout is valid");
    let mut layout = ReactionLayout::arrange(side.clone(), side).expect("empty sides arrange");
    let unchanged = layout.clone();

    assert_eq!(
        layout.set_arrow(Point2D::new(0.0, 0.0), Point2D::new(0.0, 0.0)),
        Err(ReactionLayoutError::DegenerateArrow {
            position: Point2D::new(0.0, 0.0),
        })
    );
    assert_eq!(layout, unchanged);
    assert!(matches!(
        layout.set_arrow(Point2D::new(0.0, 0.0), Point2D::new(f64::NAN, 0.0)),
        Err(ReactionLayoutError::NonFiniteArrow { .. })
    ));
    assert_eq!(layout, unchanged);

    layout
        .set_arrow(Point2D::new(-2.0, 1.0), Point2D::new(2.0, 1.0))
        .expect("finite distinct arrow endpoints");
    assert_eq!(layout.arrow_start(), Point2D::new(-2.0, 1.0));
    assert_eq!(layout.arrow_end(), Point2D::new(2.0, 1.0));
}

#[rstest]
fn test_reaction_layout_side_mutation() {
    let side = MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).expect("finite layout");
    let mut layout = ReactionLayout::try_new(
        side.clone(),
        side,
        Point2D::new(-1.0, 0.0),
        Point2D::new(1.0, 0.0),
    )
    .expect("finite distinct arrow endpoints");

    layout
        .lhs_mut()
        .set_position(AtomId(0), Point2D::new(-3.0, 0.0))
        .expect("in-frame finite edit");
    layout
        .rhs_mut()
        .set_position(AtomId(0), Point2D::new(3.0, 0.0))
        .expect("in-frame finite edit");

    assert_eq!(layout.lhs().positions(), [Point2D::new(-3.0, 0.0)]);
    assert_eq!(layout.rhs().positions(), [Point2D::new(3.0, 0.0)]);
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
