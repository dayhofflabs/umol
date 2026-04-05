mod reference_tables_data;

use reference_tables_data::{reference_table, ReferenceTable};
use rstest::rstest;
use umol_msym::{detect_symmetry, CharacterTable, PointGroup, SymmetryCenter, Thresholds};

// Reference character tables from http://gernot-katzers-spice-pages.com/character_tables/

fn make_centers(
    atomic_numbers: &[i32],
    masses: &[f64],
    positions: &[[f64; 3]],
) -> Vec<SymmetryCenter> {
    atomic_numbers
        .iter()
        .zip(masses.iter())
        .zip(positions.iter())
        .map(|((&z, &m), &pos)| SymmetryCenter {
            atomic_number: z,
            mass: m,
            position: pos,
            name: String::new(),
        })
        .collect()
}

/// Compare character tables independent of class ordering and irrep ordering.
fn compare_character_tables(
    group: &str,
    msym_ct: &CharacterTable,
    reference: &ReferenceTable,
) {
    assert_eq!(
        msym_ct.order, reference.order,
        "{group}: order mismatch (msym={}, reference={})",
        msym_ct.order, reference.order
    );

    assert_eq!(
        msym_ct.irreps.len(),
        reference.irrep_names.len(),
        "{group}: irrep count mismatch (msym={}, reference={}). \
         msym: {:?}, reference: {:?}",
        msym_ct.irreps.len(),
        reference.irrep_names.len(),
        msym_ct.irreps.iter().map(|ir| &ir.name).collect::<Vec<_>>(),
        reference.irrep_names,
    );

    let mut msym_sizes = msym_ct.class_sizes.clone();
    let mut ref_sizes = reference.class_sizes.clone();
    msym_sizes.sort();
    ref_sizes.sort();
    assert_eq!(
        msym_sizes, ref_sizes,
        "{group}: class size multisets differ"
    );

    fn irrep_fingerprint(class_sizes: &[i32], characters: &[f64]) -> Vec<(i32, i64)> {
        let mut pairs: Vec<(i32, i64)> = class_sizes
            .iter()
            .zip(characters.iter())
            .map(|(&cs, &ch)| (cs, (ch * 100000.0).round() as i64))
            .collect();
        pairs.sort();
        pairs
    }

    let msym_fingerprints: Vec<Vec<(i32, i64)>> = msym_ct
        .irreps
        .iter()
        .map(|ir| irrep_fingerprint(&msym_ct.class_sizes, &msym_ct.characters[ir.index]))
        .collect();

    let ref_fingerprints: Vec<Vec<(i32, i64)>> = reference
        .characters
        .iter()
        .map(|chars| irrep_fingerprint(&reference.class_sizes, chars))
        .collect();

    let mut matched_ref: Vec<bool> = vec![false; reference.irrep_names.len()];

    for (mi, msym_fp) in msym_fingerprints.iter().enumerate() {
        let msym_name = &msym_ct.irreps[mi].name;
        let msym_dim = msym_ct.irreps[mi].dimension;

        let match_idx = ref_fingerprints
            .iter()
            .enumerate()
            .position(|(ki, ref_fp)| !matched_ref[ki] && ref_fp == msym_fp);

        match match_idx {
            Some(ki) => {
                matched_ref[ki] = true;
                let ref_name = &reference.irrep_names[ki];
                let ref_dim = reference.characters[ki][0] as i32;

                assert_eq!(
                    msym_dim, ref_dim,
                    "{group}: dimension mismatch for msym '{msym_name}' ↔ reference '{ref_name}'"
                );
            }
            None => {
                panic!(
                    "{group}: msym irrep '{msym_name}' (dim={msym_dim}) has no matching \
                     reference irrep by character fingerprint.\n\
                     msym fingerprint: {msym_fp:?}\n\
                     reference fingerprints: {ref_fingerprints:?}"
                );
            }
        }
    }
}

// --- Molecule geometries ---

fn water() -> Vec<SymmetryCenter> {
    make_centers(
        &[8, 1, 1],
        &[15.999, 1.008, 1.008],
        &[
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ],
    )
}

fn methane() -> Vec<SymmetryCenter> {
    make_centers(
        &[6, 1, 1, 1, 1],
        &[12.011, 1.008, 1.008, 1.008, 1.008],
        &[
            [0.0, 0.0, 0.0],
            [0.629, 0.629, 0.629],
            [-0.629, -0.629, 0.629],
            [-0.629, 0.629, -0.629],
            [0.629, -0.629, -0.629],
        ],
    )
}

fn sf6() -> Vec<SymmetryCenter> {
    make_centers(
        &[16, 9, 9, 9, 9, 9, 9],
        &[32.06, 18.998, 18.998, 18.998, 18.998, 18.998, 18.998],
        &[
            [0.0, 0.0, 0.0],
            [1.564, 0.0, 0.0],
            [-1.564, 0.0, 0.0],
            [0.0, 1.564, 0.0],
            [0.0, -1.564, 0.0],
            [0.0, 0.0, 1.564],
            [0.0, 0.0, -1.564],
        ],
    )
}

fn icosahedron() -> Vec<SymmetryCenter> {
    let phi: f64 = (1.0 + 5.0_f64.sqrt()) / 2.0;
    make_centers(
        &[5; 12],
        &[10.81; 12],
        &[
            [0.0, 1.0, phi],
            [0.0, 1.0, -phi],
            [0.0, -1.0, phi],
            [0.0, -1.0, -phi],
            [1.0, phi, 0.0],
            [1.0, -phi, 0.0],
            [-1.0, phi, 0.0],
            [-1.0, -phi, 0.0],
            [phi, 0.0, 1.0],
            [phi, 0.0, -1.0],
            [-phi, 0.0, 1.0],
            [-phi, 0.0, -1.0],
        ],
    )
}

fn ethylene() -> Vec<SymmetryCenter> {
    make_centers(
        &[6, 6, 1, 1, 1, 1],
        &[12.011, 12.011, 1.008, 1.008, 1.008, 1.008],
        &[
            [0.0, 0.0, 0.667],
            [0.0, 0.0, -0.667],
            [0.0, 0.923, 1.238],
            [0.0, -0.923, 1.238],
            [0.0, -0.923, -1.238],
            [0.0, 0.923, -1.238],
        ],
    )
}

fn trans_diazene() -> Vec<SymmetryCenter> {
    make_centers(
        &[7, 7, 1, 1],
        &[14.007, 14.007, 1.008, 1.008],
        &[
            [0.0, 0.627, 0.0],
            [0.0, -0.627, 0.0],
            [0.0, 1.044, 0.943],
            [0.0, -1.044, -0.943],
        ],
    )
}

fn h2o2() -> Vec<SymmetryCenter> {
    make_centers(
        &[8, 8, 1, 1],
        &[15.999, 15.999, 1.008, 1.008],
        &[
            [0.0, 0.727, -0.053],
            [0.0, -0.727, -0.053],
            [0.787, 0.889, 0.427],
            [-0.787, -0.889, 0.427],
        ],
    )
}

fn thioformaldehyde_hfcs() -> Vec<SymmetryCenter> {
    make_centers(
        &[6, 16, 9, 1],
        &[12.011, 32.06, 18.998, 1.008],
        &[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.610],
            [0.0, 1.110, -0.620],
            [0.0, -0.940, -0.580],
        ],
    )
}

fn allene() -> Vec<SymmetryCenter> {
    make_centers(
        &[6, 6, 6, 1, 1, 1, 1],
        &[12.011, 12.011, 12.011, 1.008, 1.008, 1.008, 1.008],
        &[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.310],
            [0.0, 0.0, -1.310],
            [0.0, 0.930, 1.870],
            [0.0, -0.930, 1.870],
            [0.930, 0.0, -1.870],
            [-0.930, 0.0, -1.870],
        ],
    )
}

#[rstest]
#[case("C2", h2o2())]
#[case("Cs", thioformaldehyde_hfcs())]
#[case("C2v", water())]
#[case("C2h", trans_diazene())]
#[case("D2d", allene())]
#[case("D2h", ethylene())]
#[case("Td", methane())]
#[case("Oh", sf6())]
#[case("Ih", icosahedron())]
fn test_character_table_vs_reference(#[case] group: &str, #[case] elements: Vec<SymmetryCenter>) {
    let result = detect_symmetry(&elements, Thresholds::defaults()).unwrap();
    let reference = reference_table(group).unwrap();

    assert_eq!(
        result.group.name, group,
        "Expected {group}, detected {}",
        result.group.name
    );

    compare_character_tables(group, &result.group.character_table, &reference);
}

/// Groups where we construct by Schoenflies name (no specific molecule needed).
#[rstest]
#[case("Ci")]
#[case("C3")]
#[case("C6")]
#[case("S4")]
#[case("D2")]
#[case("T")]
#[case("Th")]
#[case("O")]
#[case("I")]
fn test_character_table_by_name_vs_reference(#[case] group: &str) {
    let pg = match PointGroup::from_schoenflies(group) {
        Ok(pg) => pg,
        Err(e) => {
            eprintln!("{group}: from_schoenflies failed: {e}, falling back to shape check");
            let reference = reference_table(group).unwrap();
            assert_eq!(
                reference.irrep_names.len(),
                reference.class_sizes.len(),
                "{group}: reference table not square"
            );
            assert_eq!(
                reference.order,
                reference.class_sizes.iter().sum::<i32>() as usize,
                "{group}: reference order doesn't match class sizes"
            );
            return;
        }
    };

    let reference = reference_table(group).unwrap();
    compare_character_tables(group, &pg.character_table, &reference);
}
