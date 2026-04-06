use crate::basis::{BasisFunction, IrrepBasis, Salc, SalcBasis};
use crate::context::Context;
use crate::error::Error;
use crate::point_group::PointGroup;
use crate::types::{EquivalenceSet, SchoenfliesLabel, SymmetryCenter, Thresholds};

/// Result of symmetry detection or symmetrization.
#[derive(Debug)]
pub struct SymmetryResult {
    pub group: &'static PointGroup,
    pub equivalence_sets: Vec<EquivalenceSet>,
    /// Atom positions after processing. Same as input for `detect_symmetry`,
    /// snapped to exact symmetry for `symmetrize`.
    pub centers: Vec<SymmetryCenter>,
}

fn c1_result(centers: &[SymmetryCenter]) -> SymmetryResult {
    let equivalence_sets = centers
        .iter()
        .map(|c| EquivalenceSet {
            centers: vec![c.clone()],
            max_error: 0.0,
        })
        .collect();
    SymmetryResult {
        group: PointGroup::c1(),
        equivalence_sets,
        centers: centers.to_vec(),
    }
}

/// Detect point group symmetry of a set of atoms.
///
/// Centers must have positions in Angstroms (libmsym convention).
pub fn detect_symmetry(
    centers: &[SymmetryCenter],
    thresholds: Thresholds,
) -> Result<SymmetryResult, Error> {
    let mut ctx = Context::new()?;
    ctx.set_centers(centers)?;
    ctx.set_thresholds(&thresholds)?;
    match ctx.find_symmetry() {
        Ok(()) => {}
        Err(e) if e.code == umol_msym_sys::MSYM_POINT_GROUP_ERROR => {
            return Ok(c1_result(centers));
        }
        Err(e) => return Err(e),
    }

    let group = PointGroup::from_context(&ctx)?;
    let equivalence_sets = ctx.equivalence_sets()?;
    let centers = ctx.centers()?;

    Ok(SymmetryResult {
        group,
        equivalence_sets,
        centers,
    })
}

/// Detect symmetry and snap atom positions to exact symmetry.
///
/// Centers must have positions in Angstroms (libmsym convention).
/// Returned centers have symmetrized positions.
pub fn symmetrize(
    centers: &[SymmetryCenter],
    thresholds: Thresholds,
) -> Result<SymmetryResult, Error> {
    let mut ctx = Context::new()?;
    ctx.set_centers(centers)?;
    ctx.set_thresholds(&thresholds)?;
    match ctx.find_symmetry() {
        Ok(()) => {}
        Err(e) if e.code == umol_msym_sys::MSYM_POINT_GROUP_ERROR => {
            return Ok(c1_result(centers));
        }
        Err(e) => return Err(e),
    }
    let group = PointGroup::from_context(&ctx)?;
    let equivalence_sets = ctx.equivalence_sets()?;
    ctx.symmetrize_centers()?;
    let centers = ctx.centers()?;

    Ok(SymmetryResult {
        group,
        equivalence_sets,
        centers,
    })
}

/// Generate a full molecule from an asymmetric unit and a target point group.
///
/// The asymmetric unit atoms are replicated by the group operations.
/// Centers must have positions in Angstroms (libmsym convention).
/// Generate a full molecule from an asymmetric unit and a target point group.
///
/// Atoms in the asymmetric unit must be positioned relative to the molecular
/// center of mass (the origin). On-axis atoms will generate fewer copies than
/// general-position atoms. Centers must have positions in Angstroms.
pub fn symmetrize_to(
    label: SchoenfliesLabel,
    asymmetric_unit: &[SymmetryCenter],
    thresholds: Thresholds,
) -> Result<SymmetryResult, Error> {
    let mut ctx = Context::new()?;
    ctx.set_thresholds(&thresholds)?;
    // Set asymmetric unit to establish center of mass in the context.
    // Then override center of mass to origin so generation works correctly.
    ctx.set_centers(asymmetric_unit)?;
    ctx.set_center_of_mass([0.0, 0.0, 0.0])?;
    ctx.set_point_group(label)?;
    ctx.generate_centers(asymmetric_unit)?;
    ctx.find_symmetry()?;

    let group = PointGroup::from_context(&ctx)?;
    let equivalence_sets = ctx.equivalence_sets()?;
    let centers = ctx.centers()?;

    Ok(SymmetryResult {
        group,
        equivalence_sets,
        centers,
    })
}

/// Compute symmetry-adapted linear combinations (SALCs) for a set of basis functions.
///
/// Centers and basis functions must be consistent: each `BasisFunction::atom_index`
/// must be a valid index into `centers`. Symmetry is detected first, then basis
/// functions are projected onto the irreps of the detected point group.
pub fn compute_salcs(
    centers: &[SymmetryCenter],
    basis: &[BasisFunction],
    thresholds: Thresholds,
) -> Result<SalcBasis, Error> {
    let mut ctx = Context::new()?;
    ctx.set_centers(centers)?;
    ctx.set_thresholds(&thresholds)?;
    ctx.find_symmetry()?;

    let group = PointGroup::from_context(&ctx)?;
    ctx.set_basis_functions(basis)?;

    let l = basis.len();
    let (coefficients, species) = ctx.salcs(l)?;

    let irreps = group.irreps();
    let zero_thresh = thresholds.zero;

    // Group SALCs by irrep
    let mut irrep_salcs: Vec<Vec<Salc>> = vec![Vec::new(); irreps.len()];

    for (salc_idx, &species_idx) in species.iter().enumerate() {
        let row = &coefficients[salc_idx * l..(salc_idx + 1) * l];
        let sparse: Vec<(usize, f64)> = row
            .iter()
            .enumerate()
            .filter(|(_, &c)| c.abs() > zero_thresh)
            .map(|(j, &c)| (j, c))
            .collect();
        if sparse.is_empty() {
            continue;
        }
        irrep_salcs[species_idx as usize].push(Salc { coefficients: sparse });
    }

    let irrep_bases: Vec<IrrepBasis> = irreps
        .into_iter()
        .zip(irrep_salcs)
        .filter(|(_, salcs)| !salcs.is_empty())
        .map(|(irrep, salcs)| IrrepBasis { irrep, salcs })
        .collect();

    Ok(SalcBasis {
        basis_functions: basis.to_vec(),
        irreps: irrep_bases,
    })
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::basis::{BasisKind, CartesianAxis};

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

    #[rstest]
    #[case(
        &[8, 1, 1],
        &[15.999, 1.008, 1.008],
        &[[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
        "C2v", 4, 2
    )]
    #[case(
        &[6, 1, 1, 1, 1],
        &[12.011, 1.008, 1.008, 1.008, 1.008],
        &[[0.0, 0.0, 0.0], [0.629, 0.629, 0.629], [-0.629, -0.629, 0.629],
          [-0.629, 0.629, -0.629], [0.629, -0.629, -0.629]],
        "Td", 24, 2
    )]
    fn test_detect_symmetry(
        #[case] zs: &[i32],
        #[case] masses: &[f64],
        #[case] positions: &[[f64; 3]],
        #[case] expected_group: &str,
        #[case] expected_order: usize,
        #[case] expected_eq_sets: usize,
    ) {
        let centers = make_centers(zs, masses, positions);
        let result = detect_symmetry(&centers, Thresholds::default()).unwrap();
        assert_eq!(result.group.to_string(), expected_group);
        assert_eq!(result.group.order(), expected_order);
        assert_eq!(result.equivalence_sets.len(), expected_eq_sets);
        assert_eq!(result.centers.len(), centers.len());
    }

    #[rstest]
    fn test_symmetrize() {
        let centers = make_centers(
            &[8, 1, 1],
            &[15.999, 1.008, 1.008],
            &[[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
        );
        let result = symmetrize(&centers, Thresholds::default()).unwrap();
        assert_eq!(result.group.to_string(), "C2v");
        assert_eq!(result.centers.len(), 3);

        // H atoms should have exactly opposite y coordinates after symmetrization
        let y1 = result.centers[1].position[1];
        let y2 = result.centers[2].position[1];
        assert!((y1 + y2).abs() < 1e-10, "H atoms not symmetric: {y1} vs {y2}");
    }

    #[rstest]
    #[case(
        SchoenfliesLabel::Cnv(2),
        &[8], &[15.999], &[[0.0, 0.0, 0.117]],
        &[1],  &[1.008], &[[0.0, 0.757, -0.469]],
        "C2v", 3
    )]
    #[case(
        SchoenfliesLabel::Td,
        &[6], &[12.011], &[[0.0, 0.0, 0.0]],
        &[1], &[1.008],  &[[0.629, 0.629, 0.629]],
        "Td", 5
    )]
    #[case(
        SchoenfliesLabel::Dnh(6),
        &[6], &[12.011], &[[1.4, 0.0, 0.0]],
        &[1], &[1.008],  &[[2.48, 0.0, 0.0]],
        "D6h", 12
    )]
    fn test_symmetrize_to(
        #[case] label: SchoenfliesLabel,
        #[case] zs1: &[i32],
        #[case] ms1: &[f64],
        #[case] ps1: &[[f64; 3]],
        #[case] zs2: &[i32],
        #[case] ms2: &[f64],
        #[case] ps2: &[[f64; 3]],
        #[case] expected_group: &str,
        #[case] expected_atoms: usize,
    ) {
        let mut asym = make_centers(zs1, ms1, ps1);
        asym.extend(make_centers(zs2, ms2, ps2));
        let result = symmetrize_to(label, &asym, Thresholds::default()).unwrap();
        assert_eq!(result.group.to_string(), expected_group);
        assert_eq!(result.centers.len(), expected_atoms);
    }

    fn s_basis(n_atoms: usize) -> Vec<BasisFunction> {
        (0..n_atoms)
            .map(|i| BasisFunction {
                atom_index: i,
                kind: BasisKind::RealSphericalHarmonic,
                n: 1,
                l: 0,
                m: 0,
            })
            .collect()
    }

    fn displacement_basis(n_atoms: usize) -> Vec<BasisFunction> {
        let axes = [
            (CartesianAxis::X, 1),
            (CartesianAxis::Y, -1),
            (CartesianAxis::Z, 0),
        ];
        (0..n_atoms)
            .flat_map(|i| {
                axes.iter().map(move |&(axis, m)| BasisFunction {
                    atom_index: i,
                    kind: BasisKind::Displacement(axis),
                    n: 2,
                    l: 1,
                    m,
                })
            })
            .collect()
    }

    #[rstest]
    #[case(
        &[8, 1, 1],
        &[15.999, 1.008, 1.008],
        &[[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
        3, &["A1", "A1", "B1"]
    )]
    #[case(
        &[6, 1, 1, 1, 1],
        &[12.011, 1.008, 1.008, 1.008, 1.008],
        &[[0.0, 0.0, 0.0], [0.629, 0.629, 0.629], [-0.629, -0.629, 0.629],
          [-0.629, 0.629, -0.629], [0.629, -0.629, -0.629]],
        5, &["A1", "A1", "T2", "T2", "T2"]
    )]
    fn test_compute_salcs_s_basis(
        #[case] zs: &[i32],
        #[case] masses: &[f64],
        #[case] positions: &[[f64; 3]],
        #[case] expected_total: usize,
        #[case] expected_salc_irreps: &[&str],
    ) {
        let centers = make_centers(zs, masses, positions);
        let basis = s_basis(centers.len());
        let result = compute_salcs(&centers, &basis, Thresholds::default()).unwrap();

        let total: usize = result.irreps.iter().map(|ib| ib.salcs.len()).sum();
        assert_eq!(total, expected_total);

        let mut salc_symbols: Vec<String> = result
            .irreps
            .iter()
            .flat_map(|ib| {
                std::iter::repeat_n(ib.irrep.symbol().to_owned(), ib.salcs.len())
            })
            .collect();
        salc_symbols.sort();
        let mut expected: Vec<&str> = expected_salc_irreps.to_vec();
        expected.sort();
        assert_eq!(salc_symbols, expected);
    }

    #[rstest]
    fn test_compute_salcs_displacement_basis() {
        let centers = make_centers(
            &[8, 1, 1],
            &[15.999, 1.008, 1.008],
            &[[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
        );
        let basis = displacement_basis(centers.len());
        let result = compute_salcs(&centers, &basis, Thresholds::default()).unwrap();

        let total: usize = result.irreps.iter().map(|ib| ib.salcs.len()).sum();
        assert_eq!(total, 9); // 3N = 9

        // Verify orthogonality: SALC rows should be orthonormal
        for ib in &result.irreps {
            for (i, s1) in ib.salcs.iter().enumerate() {
                for (j, s2) in ib.salcs.iter().enumerate() {
                    let dot: f64 = s1.coefficients.iter().map(|&(k, c1)| {
                        s2.coefficients.iter()
                            .find(|&&(k2, _)| k2 == k)
                            .map_or(0.0, |&(_, c2)| c1 * c2)
                    }).sum();
                    if i == j {
                        assert!((dot - 1.0).abs() < 1e-8, "SALC {i} not normalized: {dot}");
                    } else {
                        assert!(dot.abs() < 1e-8, "SALCs {i},{j} not orthogonal: {dot}");
                    }
                }
            }
        }
    }
}
