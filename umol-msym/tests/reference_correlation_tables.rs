use nalgebra::Vector3;
use rstest::rstest;
use umol_msym::{lower_symmetry, SchoenfliesSymbol, SymmetryCenter, Thresholds};

fn water() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter {
            atomic_number: 8,
            mass: 15.999,
            position: Vector3::new(0.0, 0.0, 0.117_370_3),
            name: "O".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: Vector3::new(0.0, 0.757_160_4, -0.469_481_2),
            name: "H".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: Vector3::new(0.0, -0.757_160_4, -0.469_481_2),
            name: "H".into(),
        },
    ]
}

fn methane() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter {
            atomic_number: 6,
            mass: 12.011,
            position: Vector3::new(0.0, 0.0, 0.0),
            name: "C".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: Vector3::new(0.629_118_5, 0.629_118_5, 0.629_118_5),
            name: "H".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: Vector3::new(-0.629_118_5, -0.629_118_5, 0.629_118_5),
            name: "H".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: Vector3::new(-0.629_118_5, 0.629_118_5, -0.629_118_5),
            name: "H".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: Vector3::new(0.629_118_5, -0.629_118_5, -0.629_118_5),
            name: "H".into(),
        },
    ]
}

fn sf6() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter {
            atomic_number: 16,
            mass: 32.06,
            position: Vector3::new(0.0, 0.0, 0.0),
            name: "S".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: Vector3::new(1.564, 0.0, 0.0),
            name: "F".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: Vector3::new(-1.564, 0.0, 0.0),
            name: "F".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: Vector3::new(0.0, 1.564, 0.0),
            name: "F".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: Vector3::new(0.0, -1.564, 0.0),
            name: "F".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: Vector3::new(0.0, 0.0, 1.564),
            name: "F".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: Vector3::new(0.0, 0.0, -1.564),
            name: "F".into(),
        },
    ]
}

fn hcl() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter {
            atomic_number: 17,
            mass: 35.453,
            position: Vector3::new(0.0, 0.0, 0.0),
            name: "Cl".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: Vector3::new(0.0, 0.0, 1.275),
            name: "H".into(),
        },
    ]
}

fn co2() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter {
            atomic_number: 8,
            mass: 15.999,
            position: Vector3::new(0.0, 0.0, -1.16),
            name: "O".into(),
        },
        SymmetryCenter {
            atomic_number: 6,
            mass: 12.011,
            position: Vector3::new(0.0, 0.0, 0.0),
            name: "C".into(),
        },
        SymmetryCenter {
            atomic_number: 8,
            mass: 15.999,
            position: Vector3::new(0.0, 0.0, 1.16),
            name: "O".into(),
        },
    ]
}

fn assert_orthogonal_transform(t: nalgebra::Matrix3<f64>) {
    for i in 0..3 {
        let norm_sq = t.row(i).norm_squared();
        assert!(
            (norm_sq - 1.0).abs() < 1e-10,
            "transform row {i} not normalized: {norm_sq}"
        );
    }
}

#[rstest]
#[case::to_cs(SchoenfliesSymbol::Cs)]
#[case::to_c2(SchoenfliesSymbol::Cn(2))]
fn test_lower_symmetry_water(#[case] target: SchoenfliesSymbol) {
    let result = lower_symmetry(&water(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.symbol(), SchoenfliesSymbol::Cnv(2));
    assert_eq!(result.child_group.symbol(), target);
    assert_orthogonal_transform(result.transform);
}

#[rstest]
#[case::to_c3v(SchoenfliesSymbol::Cnv(3))]
#[case::to_c2v(SchoenfliesSymbol::Cnv(2))]
#[case::to_d2d(SchoenfliesSymbol::Dnd(2))]
fn test_lower_symmetry_methane(#[case] target: SchoenfliesSymbol) {
    let result = lower_symmetry(&methane(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.symbol(), SchoenfliesSymbol::Td);
    assert_eq!(result.child_group.symbol(), target);
    assert_orthogonal_transform(result.transform);
}

#[rstest]
#[case::to_td(SchoenfliesSymbol::Td)]
#[case::to_o(SchoenfliesSymbol::O)]
#[case::to_d4h(SchoenfliesSymbol::Dnh(4))]
#[case::to_d3d(SchoenfliesSymbol::Dnd(3))]
#[case::to_d2h(SchoenfliesSymbol::Dnh(2))]
fn test_lower_symmetry_sf6(#[case] target: SchoenfliesSymbol) {
    let result = lower_symmetry(&sf6(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.symbol(), SchoenfliesSymbol::Oh);
    assert_eq!(result.child_group.symbol(), target);
    assert_orthogonal_transform(result.transform);
}

#[rstest]
#[case(SchoenfliesSymbol::Cnv(4))]
#[case(SchoenfliesSymbol::Cnv(2))]
#[case(SchoenfliesSymbol::Cn(1))]
fn test_lower_symmetry_linear(#[case] target: SchoenfliesSymbol) {
    let result = lower_symmetry(&hcl(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.symbol(), SchoenfliesSymbol::Coov);
    assert_eq!(result.child_group.symbol(), target);
    assert_eq!(result.centers.len(), 2);
    assert_orthogonal_transform(result.transform);
}

#[rstest]
#[case(SchoenfliesSymbol::Dnh(2))]
#[case(SchoenfliesSymbol::Cnv(2))]
#[case(SchoenfliesSymbol::Cn(1))]
fn test_lower_symmetry_co2(#[case] target: SchoenfliesSymbol) {
    let result = lower_symmetry(&co2(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.symbol(), SchoenfliesSymbol::Dooh);
    assert_eq!(result.child_group.symbol(), target);
    assert_eq!(result.centers.len(), 3);
    assert_orthogonal_transform(result.transform);
}

#[rstest]
fn test_lower_symmetry_invalid_subgroup() {
    let result = lower_symmetry(&water(), SchoenfliesSymbol::Ci, Thresholds::default());
    assert!(result.is_err());
}
