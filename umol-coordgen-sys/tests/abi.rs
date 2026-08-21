use rstest::rstest;
use umol_coordgen_sys::{generate_coordinates, Bond, CoordgenError};

#[rstest]
#[case::empty_result(Vec::new(), Vec::new())]
fn test_generate_coordinates_empty_graph(
    #[case] atomic_numbers: Vec<u16>,
    #[case] bonds: Vec<Bond>,
) {
    let points = generate_coordinates(&atomic_numbers, &bonds).expect("coordinate generation");

    assert_eq!(points, Vec::new());
}

#[rstest]
#[case::finite_point(vec![6], Vec::new())]
fn test_generate_coordinates_one_atom(#[case] atomic_numbers: Vec<u16>, #[case] bonds: Vec<Bond>) {
    let points = generate_coordinates(&atomic_numbers, &bonds).expect("coordinate generation");

    assert_eq!(points.len(), 1);
    assert!(points[0].x.is_finite());
    assert!(points[0].y.is_finite());
}

#[rstest]
#[case::preserved_frame(
    vec![8, 7, 9, 6],
    vec![
        Bond {
            atom_0: 3,
            atom_1: 0,
            order: 1,
        },
        Bond {
            atom_0: 3,
            atom_1: 1,
            order: 2,
        },
        Bond {
            atom_0: 3,
            atom_1: 2,
            order: 1,
        },
    ]
)]
fn test_generate_coordinates_bonded_graph(
    #[case] atomic_numbers: Vec<u16>,
    #[case] bonds: Vec<Bond>,
) {
    let points = generate_coordinates(&atomic_numbers, &bonds).expect("coordinate generation");

    assert_eq!(points.len(), 4);
    assert!(points
        .iter()
        .all(|point| point.x.is_finite() && point.y.is_finite()));
    for terminal in 0..3 {
        let distance = (points[3].x - points[terminal].x).hypot(points[3].y - points[terminal].y);
        assert!((distance - 50.0).abs() < 0.1, "bond length was {distance}");
    }
    let terminal_distance = (points[0].x - points[1].x).hypot(points[0].y - points[1].y);
    assert!(terminal_distance > 50.0);
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
fn test_generate_coordinates_error(#[case] bond: Bond, #[case] atom_index: usize) {
    assert_eq!(
        generate_coordinates(&[6], &[bond]),
        Err(CoordgenError::BondAtomOutOfBounds {
            bond_index: 0,
            atom_index,
            atom_count: 1,
        })
    );
}
