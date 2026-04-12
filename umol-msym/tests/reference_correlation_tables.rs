use std::collections::BTreeMap;

use rstest::rstest;
use umol_msym::{
    lower_symmetry, Irrep, SchoenfliesLabel, SymmetryCenter, SymmetryDescentResult, Thresholds,
};

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

fn sf6() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter {
            atomic_number: 16,
            mass: 32.06,
            position: [0.0, 0.0, 0.0],
            name: "S".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: [1.564, 0.0, 0.0],
            name: "F".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: [-1.564, 0.0, 0.0],
            name: "F".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: [0.0, 1.564, 0.0],
            name: "F".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: [0.0, -1.564, 0.0],
            name: "F".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: [0.0, 0.0, 1.564],
            name: "F".into(),
        },
        SymmetryCenter {
            atomic_number: 9,
            mass: 18.998,
            position: [0.0, 0.0, -1.564],
            name: "F".into(),
        },
    ]
}

fn hcl() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter {
            atomic_number: 17,
            mass: 35.453,
            position: [0.0, 0.0, 0.0],
            name: "Cl".into(),
        },
        SymmetryCenter {
            atomic_number: 1,
            mass: 1.008,
            position: [0.0, 0.0, 1.275],
            name: "H".into(),
        },
    ]
}

fn co2() -> Vec<SymmetryCenter> {
    vec![
        SymmetryCenter {
            atomic_number: 8,
            mass: 15.999,
            position: [0.0, 0.0, -1.16],
            name: "O".into(),
        },
        SymmetryCenter {
            atomic_number: 6,
            mass: 12.011,
            position: [0.0, 0.0, 0.0],
            name: "C".into(),
        },
        SymmetryCenter {
            atomic_number: 8,
            mass: 15.999,
            position: [0.0, 0.0, 1.16],
            name: "O".into(),
        },
    ]
}

// Altmann & Herzig, "Point-Group Theory Tables" (2nd ed.) subduction columns.
// Each AltmannColumn is (parent_irrep, child_decomposition_multiset) from T n.9
// tables. Multiset comparison absorbs outer-automorphism label permutations
// (e.g. B1 ↔ B2); multiple columns accept non-conjugate orientations.

type AltmannRow = (&'static str, &'static [&'static str]);
type AltmannColumn = &'static [AltmannRow];

// T 50.9 (p. 483): C2v → Cs with σv = σx plane
const C2V_TO_CS_SIGMA_X: AltmannColumn = &[
    ("A1", &["A'"]),
    ("A2", &["A''"]),
    ("B1", &["A''"]),
    ("B2", &["A'"]),
];

// T 50.9 (p. 483): C2v → Cs with σv = σy plane
const C2V_TO_CS_SIGMA_Y: AltmannColumn = &[
    ("A1", &["A'"]),
    ("A2", &["A''"]),
    ("B1", &["A'"]),
    ("B2", &["A''"]),
];

// T 50.9 (p. 483): C2v → C2
const C2V_TO_C2: AltmannColumn = &[
    ("A1", &["A"]),
    ("A2", &["A"]),
    ("B1", &["B"]),
    ("B2", &["B"]),
];

// T 73.9 (p. 639): Td → (C3v). Cyclic E in libmsym is split-labeled "E1".
const TD_TO_C3V: AltmannColumn = &[
    ("A1", &["A1"]),
    ("A2", &["A2"]),
    ("E", &["E1"]),
    ("T1", &["A2", "E1"]),
    ("T2", &["A1", "E1"]),
];

// T 73.9 (p. 639): Td → (C2v). Altmann's column lists A2 → A1 by reassigning
// labels under the outer automorphism; libmsym yields the direct character-
// restriction result A2 → A2, which is mathematically equivalent.
const TD_TO_C2V: AltmannColumn = &[
    ("A1", &["A1"]),
    ("A2", &["A2"]),
    ("E", &["A1", "A2"]),
    ("T1", &["A2", "B1", "B2"]),
    ("T2", &["A1", "B1", "B2"]),
];

// T 73.9 (p. 639): Td → (D2d). Cyclic E is "E1" in libmsym.
const TD_TO_D2D: AltmannColumn = &[
    ("A1", &["A1"]),
    ("A2", &["B1"]),
    ("E", &["A1", "B1"]),
    ("T1", &["A2", "E1"]),
    ("T2", &["B2", "E1"]),
];

// T 71.9 (p. 629): Oh → (Td). Td is not a canonical subgroup of Oh; this
// column applies when the Td axes are chosen aligned with Oh's cube.
const OH_TO_TD: AltmannColumn = &[
    ("A1g", &["A1"]),
    ("A2g", &["A2"]),
    ("Eg", &["E"]),
    ("T1g", &["T1"]),
    ("T2g", &["T2"]),
    ("A1u", &["A2"]),
    ("A2u", &["A1"]),
    ("Eu", &["E"]),
    ("T1u", &["T2"]),
    ("T2u", &["T1"]),
];

// T 71.9 (p. 629): Oh → O
const OH_TO_O: AltmannColumn = &[
    ("A1g", &["A1"]),
    ("A2g", &["A2"]),
    ("Eg", &["E"]),
    ("T1g", &["T1"]),
    ("T2g", &["T2"]),
    ("A1u", &["A1"]),
    ("A2u", &["A2"]),
    ("Eu", &["E"]),
    ("T1u", &["T1"]),
    ("T2u", &["T2"]),
];

// T 71.9 (p. 630): Oh → (D4h), C4 along z. Altmann's labeling; libmsym may
// emit a B1↔B2 swap (an outer-automorphism orientation) which is accepted via
// OH_TO_D4H_SWAPPED. Cyclic Eg/Eu are "E1g/E1u" in libmsym.
const OH_TO_D4H: AltmannColumn = &[
    ("A1g", &["A1g"]),
    ("A2g", &["B1g"]),
    ("Eg", &["A1g", "B1g"]),
    ("T1g", &["A2g", "E1g"]),
    ("T2g", &["B2g", "E1g"]),
    ("A1u", &["A1u"]),
    ("A2u", &["B1u"]),
    ("Eu", &["A1u", "B1u"]),
    ("T1u", &["A2u", "E1u"]),
    ("T2u", &["B2u", "E1u"]),
];

// T 71.9 (p. 630): Oh → (D4h), B1 ↔ B2 swap of OH_TO_D4H.
const OH_TO_D4H_SWAPPED: AltmannColumn = &[
    ("A1g", &["A1g"]),
    ("A2g", &["B2g"]),
    ("Eg", &["A1g", "B2g"]),
    ("T1g", &["A2g", "E1g"]),
    ("T2g", &["B1g", "E1g"]),
    ("A1u", &["A1u"]),
    ("A2u", &["B2u"]),
    ("Eu", &["A1u", "B2u"]),
    ("T1u", &["A2u", "E1u"]),
    ("T2u", &["B1u", "E1u"]),
];

// T 71.9 (p. 630): Oh → (D3d), C3 along a body diagonal. Cyclic Eg/Eu are
// "E1g/E1u" in libmsym.
const OH_TO_D3D: AltmannColumn = &[
    ("A1g", &["A1g"]),
    ("A2g", &["A2g"]),
    ("Eg", &["E1g"]),
    ("T1g", &["A2g", "E1g"]),
    ("T2g", &["A1g", "E1g"]),
    ("A1u", &["A1u"]),
    ("A2u", &["A2u"]),
    ("Eu", &["E1u"]),
    ("T1u", &["A2u", "E1u"]),
    ("T2u", &["A1u", "E1u"]),
];

// T 71.9 (p. 630): Oh → D2h, C2 orientation (three orthogonal Oh C2 axes).
// libmsym names the totally-symmetric D2h irrep "A1g" rather than Altmann's "Ag".
const OH_TO_D2H_C2: AltmannColumn = &[
    ("A1g", &["A1g"]),
    ("A2g", &["A1g"]),
    ("Eg", &["A1g", "A1g"]),
    ("T1g", &["B1g", "B2g", "B3g"]),
    ("T2g", &["B1g", "B2g", "B3g"]),
    ("A1u", &["A1u"]),
    ("A2u", &["A1u"]),
    ("Eu", &["A1u", "A1u"]),
    ("T1u", &["B1u", "B2u", "B3u"]),
    ("T2u", &["B1u", "B2u", "B3u"]),
];

// T 71.9 (p. 630): Oh → (D2h), C2' orientation (face-diagonal C2 axes).
// The S3 outer automorphism of D2h permutes {B1,B2,B3}, so libmsym may emit
// any of three distinct orbit representatives: B1 fixed, B2 fixed, or B3 fixed.
const OH_TO_D2H_C2_PRIME_B1: AltmannColumn = &[
    ("A1g", &["A1g"]),
    ("A2g", &["B1g"]),
    ("Eg", &["A1g", "B1g"]),
    ("T1g", &["B1g", "B2g", "B3g"]),
    ("T2g", &["A1g", "B2g", "B3g"]),
    ("A1u", &["A1u"]),
    ("A2u", &["B1u"]),
    ("Eu", &["A1u", "B1u"]),
    ("T1u", &["B1u", "B2u", "B3u"]),
    ("T2u", &["A1u", "B2u", "B3u"]),
];

const OH_TO_D2H_C2_PRIME_B2: AltmannColumn = &[
    ("A1g", &["A1g"]),
    ("A2g", &["B2g"]),
    ("Eg", &["A1g", "B2g"]),
    ("T1g", &["B1g", "B2g", "B3g"]),
    ("T2g", &["A1g", "B1g", "B3g"]),
    ("A1u", &["A1u"]),
    ("A2u", &["B2u"]),
    ("Eu", &["A1u", "B2u"]),
    ("T1u", &["B1u", "B2u", "B3u"]),
    ("T2u", &["A1u", "B1u", "B3u"]),
];

const OH_TO_D2H_C2_PRIME_B3: AltmannColumn = &[
    ("A1g", &["A1g"]),
    ("A2g", &["B3g"]),
    ("Eg", &["A1g", "B3g"]),
    ("T1g", &["B1g", "B2g", "B3g"]),
    ("T2g", &["A1g", "B1g", "B2g"]),
    ("A1u", &["A1u"]),
    ("A2u", &["B3u"]),
    ("Eu", &["A1u", "B3u"]),
    ("T1u", &["B1u", "B2u", "B3u"]),
    ("T2u", &["A1u", "B1u", "B2u"]),
];

fn decomp_multiset(row: &[(Irrep, u32)]) -> Vec<String> {
    let mut symbols: Vec<String> = row
        .iter()
        .flat_map(|(ir, n)| std::iter::repeat(ir.symbol().to_string()).take(*n as usize))
        .collect();
    symbols.sort();
    symbols
}

fn build_actual(result: &SymmetryDescentResult) -> BTreeMap<String, Vec<String>> {
    let ct = result
        .correlation
        .as_ref()
        .expect("correlation table missing");
    result
        .parent_group
        .irreps()
        .iter()
        .zip(ct.rows.iter())
        .map(|(ir, row)| (ir.symbol().to_string(), decomp_multiset(row)))
        .collect()
}

fn build_expected(column: AltmannColumn) -> BTreeMap<String, Vec<String>> {
    column
        .iter()
        .map(|(parent, children)| {
            let mut sorted: Vec<String> = children.iter().map(|s| (*s).to_string()).collect();
            sorted.sort();
            ((*parent).to_string(), sorted)
        })
        .collect()
}

fn format_decomp(ms: &BTreeMap<String, Vec<String>>) -> String {
    ms.iter()
        .map(|(p, c)| format!("  {} → {}", p, c.join(" + ")))
        .collect::<Vec<_>>()
        .join("\n")
}

fn assert_descent_matches_altmann(
    result: &SymmetryDescentResult,
    acceptable: &[AltmannColumn],
) {
    let actual = build_actual(result);
    for column in acceptable {
        if actual == build_expected(column) {
            return;
        }
    }
    let expected_dump: String = acceptable
        .iter()
        .enumerate()
        .map(|(i, col)| {
            format!(
                "--- accepted column {i} ---\n{}",
                format_decomp(&build_expected(col))
            )
        })
        .collect::<Vec<_>>()
        .join("\n");
    panic!(
        "{} → {} subduction does not match Altmann.\n--- actual ---\n{}\n{}",
        result.parent_group.label(),
        result.child_group.label(),
        format_decomp(&actual),
        expected_dump,
    );
}

fn assert_orthogonal_transform(t: [[f64; 3]; 3]) {
    (0..3).for_each(|i| {
        let norm: f64 = (0..3).map(|j| t[i][j] * t[i][j]).sum();
        assert!(
            (norm - 1.0).abs() < 1e-10,
            "transform row {i} not normalized: {norm}"
        );
    });
}

#[rstest]
#[case::to_cs(SchoenfliesLabel::Cs, &[C2V_TO_CS_SIGMA_X, C2V_TO_CS_SIGMA_Y])]
#[case::to_c2(SchoenfliesLabel::Cn(2), &[C2V_TO_C2])]
fn test_lower_symmetry_water(
    #[case] target: SchoenfliesLabel,
    #[case] acceptable: &[AltmannColumn],
) {
    let result = lower_symmetry(&water(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.label(), SchoenfliesLabel::Cnv(2));
    assert_eq!(result.child_group.label(), target);
    assert_descent_matches_altmann(&result, acceptable);
    assert_orthogonal_transform(result.transform);
}

#[rstest]
#[case::to_c3v(SchoenfliesLabel::Cnv(3), &[TD_TO_C3V])]
#[case::to_c2v(SchoenfliesLabel::Cnv(2), &[TD_TO_C2V])]
#[case::to_d2d(SchoenfliesLabel::Dnd(2), &[TD_TO_D2D])]
fn test_lower_symmetry_methane(
    #[case] target: SchoenfliesLabel,
    #[case] acceptable: &[AltmannColumn],
) {
    let result = lower_symmetry(&methane(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.label(), SchoenfliesLabel::Td);
    assert_eq!(result.child_group.label(), target);
    assert_descent_matches_altmann(&result, acceptable);
    assert_orthogonal_transform(result.transform);
}

#[rstest]
#[case::to_td(SchoenfliesLabel::Td, &[OH_TO_TD])]
#[case::to_o(SchoenfliesLabel::O, &[OH_TO_O])]
#[case::to_d4h(SchoenfliesLabel::Dnh(4), &[OH_TO_D4H, OH_TO_D4H_SWAPPED])]
#[case::to_d3d(SchoenfliesLabel::Dnd(3), &[OH_TO_D3D])]
#[case::to_d2h(
    SchoenfliesLabel::Dnh(2),
    &[
        OH_TO_D2H_C2,
        OH_TO_D2H_C2_PRIME_B1,
        OH_TO_D2H_C2_PRIME_B2,
        OH_TO_D2H_C2_PRIME_B3,
    ],
)]
fn test_lower_symmetry_sf6(
    #[case] target: SchoenfliesLabel,
    #[case] acceptable: &[AltmannColumn],
) {
    let result = lower_symmetry(&sf6(), target, Thresholds::default()).unwrap();
    assert_eq!(result.parent_group.label(), SchoenfliesLabel::Oh);
    assert_eq!(result.child_group.label(), target);
    assert_descent_matches_altmann(&result, acceptable);
    assert_orthogonal_transform(result.transform);
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
    assert!(result.correlation.is_none());
    assert_orthogonal_transform(result.transform);
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
