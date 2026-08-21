use rstest::rstest;
use umol_graph_ir::ir::AtomId;
use umol_graph_ir::mol;
use umol_io::layout::{MoleculeLayout, MoleculeLayoutError, Point2D};

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
