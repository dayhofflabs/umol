use rstest::rstest;
use umol_msym::{lower_symmetry, SchoenfliesLabel, SymmetryCenter, Thresholds};

fn water() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter {
            atomic_number: 8,
            mass: 15.999,
            position: [0.0, 0.0, 0.117_370_3],
            name: "O".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: [0.0, 0.757_160_4, -0.469_481_2],
            name: "H".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: [0.0, -0.757_160_4, -0.469_481_2],
            name: "H".into(),
        },
    ]
}

fn methane() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter {
            atomic_number: 6,
            mass: 12.011,
            position: [0.0, 0.0, 0.0],
            name: "C".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: [0.629_118_5, 0.629_118_5, 0.629_118_5],
            name: "H".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: [-0.629_118_5, -0.629_118_5, 0.629_118_5],
            name: "H".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: [-0.629_118_5, 0.629_118_5, -0.629_118_5],
            name: "H".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: [0.629_118_5, -0.629_118_5, -0.629_118_5],
            name: "H".into(),
        },
    ]
}

fn assert_correlation_valid(result: &umol_msym::SymmetryDescentResult) {
    let ct = match &result.correlation {
        Some(ct) => ct,
        None => return, // no correlation for infinite → finite
    };
    let parent_irreps = result.parent_group.irreps();
    assert_eq!(ct.rows.len(), parent_irreps.len());

    for (i, row) in ct.rows.iter().enumerate() {
        assert!(
            !row.is_empty(),
            "empty correlation row for {} in {} → {}",
            parent_irreps[i].symbol(),
            result.parent_group.label(),
            result.child_group.label()
        );

        for (ir, n) in row {
            assert!(
                *n > 0,
                "zero multiplicity for {} in row {} of {} → {}",
                ir.symbol(),
                parent_irreps[i].symbol(),
                result.parent_group.label(),
                result.child_group.label()
            );
        }
    }

    if !result.child_group.has_complex_irreps() {
        for (i, row) in ct.rows.iter().enumerate() {
            let parent_dim = parent_irreps[i].dimension();
            let child_sum: i32 = row.iter().map(|(ir, n)| ir.dimension() * *n as i32).sum();
            assert_eq!(
                parent_dim, child_sum,
                "dimension mismatch for {} → {}: {} has dim {} but child sums to {}",
                result.parent_group.label(),
                result.child_group.label(),
                parent_irreps[i].symbol(),
                parent_dim,
                child_sum
            );
        }
    }
}

fn assert_orthogonal_transform(t: [[f64; 3]; 3]) {
    for i in 0..3 {
        let norm: f64 = (0..3).map(|j| t[i][j] * t[i][j]).sum();
        assert!(
            (norm - 1.0).abs() < 1e-10,
            "transform row {i} not normalized: {norm}"
        );
    }
}

#[rstest]
#[case(SchoenfliesLabel::Cnv(2))]
#[case(SchoenfliesLabel::Cs)]
#[case(SchoenfliesLabel::Cn(2))]
#[case(SchoenfliesLabel::Cn(1))]
fn test_lower_symmetry_water(#[case] target: SchoenfliesLabel) {
    let result = lower_symmetry(&water(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.label(), SchoenfliesLabel::Cnv(2));
    assert_eq!(result.child_group.label(), target);
    assert_eq!(result.centers.len(), 3);
    assert_correlation_valid(&result);
    assert_orthogonal_transform(result.transform);
}

#[rstest]
#[case(SchoenfliesLabel::T)]
#[case(SchoenfliesLabel::Cnv(3))]
#[case(SchoenfliesLabel::Cnv(2))]
#[case(SchoenfliesLabel::Dnd(2))]
fn test_lower_symmetry_methane(#[case] target: SchoenfliesLabel) {
    let result = lower_symmetry(&methane(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.label(), SchoenfliesLabel::Td);
    assert_eq!(result.child_group.label(), target);
    assert_eq!(result.centers.len(), 5);
    assert_correlation_valid(&result);
    assert_orthogonal_transform(result.transform);
}

fn hcl() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter { atomic_number: 17, mass: 35.453, position: [0.0, 0.0, 0.0], name: "Cl".into() },
        SymmetryCenter { atomic_number: 1, mass: 1.008, position: [0.0, 0.0, 1.275], name: "H".into() },
    ]
}

#[rstest]
#[case(SchoenfliesLabel::Cnv(4))]
#[case(SchoenfliesLabel::Cnv(2))]
#[case(SchoenfliesLabel::Cn(1))]
fn test_lower_symmetry_linear(#[case] target: SchoenfliesLabel) {
    let result = lower_symmetry(&hcl(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.label(), SchoenfliesLabel::Coov);
    assert_eq!(result.child_group.label(), target);
    assert_eq!(result.centers.len(), 2);
    if target == SchoenfliesLabel::Cn(1) {
        assert!(result.correlation.is_none());
    } else {
        assert!(result.correlation.is_none());
    }
    assert_orthogonal_transform(result.transform);
}

fn co2() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter { atomic_number: 8, mass: 15.999, position: [0.0, 0.0, -1.16], name: "O".into() },
        SymmetryCenter { atomic_number: 6, mass: 12.011, position: [0.0, 0.0, 0.0], name: "C".into() },
        SymmetryCenter { atomic_number: 8, mass: 15.999, position: [0.0, 0.0, 1.16], name: "O".into() },
    ]
}

#[rstest]
#[case(SchoenfliesLabel::Dnh(2))]
#[case(SchoenfliesLabel::Cnv(2))]
#[case(SchoenfliesLabel::Cn(1))]
fn test_lower_symmetry_co2(#[case] target: SchoenfliesLabel) {
    let result = lower_symmetry(&co2(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.label(), SchoenfliesLabel::Dooh);
    assert_eq!(result.child_group.label(), target);
    assert_eq!(result.centers.len(), 3);
    assert!(result.correlation.is_none());
    assert_orthogonal_transform(result.transform);
}

#[rstest]
fn test_lower_symmetry_invalid_subgroup() {
    let result = lower_symmetry(&water(), SchoenfliesLabel::Ci, Thresholds::default());
    assert!(result.is_err());
}
