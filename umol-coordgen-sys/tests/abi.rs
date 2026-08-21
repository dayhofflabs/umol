use rstest::rstest;
use umol_coordgen_sys::{generate_coordinates, Bond, CoordgenError};

#[rstest]
fn test_generate_coordinates_empty_graph() {
    let points = generate_coordinates(&[], &[]).unwrap();

    assert!(points.is_empty());
}

#[rstest]
fn test_generate_coordinates_one_atom() {
    let points = generate_coordinates(&[6], &[]).unwrap();

    assert_eq!(points.len(), 1);
    assert!(points[0].x.is_finite());
    assert!(points[0].y.is_finite());
}

#[rstest]
fn test_generate_coordinates_bonded_graph() {
    let points = generate_coordinates(
        &[6, 8],
        &[Bond {
            atom_0: 0,
            atom_1: 1,
            order: 1,
        }],
    )
    .unwrap();

    assert_eq!(points.len(), 2);
    assert!(points
        .iter()
        .all(|point| point.x.is_finite() && point.y.is_finite()));
    assert_ne!(points[0], points[1]);
}

#[rstest]
#[case::first_endpoint(
    Bond {
        atom_0: 2,
        atom_1: 0,
        order: 1,
    },
    2
)]
#[case::second_endpoint(
    Bond {
        atom_0: 0,
        atom_1: 2,
        order: 1,
    },
    2
)]
fn test_generate_coordinates_bond_atom_out_of_bounds(
    #[case] bond: Bond,
    #[case] atom_index: usize,
) {
    assert_eq!(
        generate_coordinates(&[6], &[bond]),
        Err(CoordgenError::BondAtomOutOfBounds {
            bond_index: 0,
            atom_index,
            atom_count: 1,
        })
    );
}
