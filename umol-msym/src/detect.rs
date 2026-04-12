use crate::basis::{BasisFunction, IrrepBasis, Salc, SalcBasis};
use crate::context::Context;
use crate::error::MsymError;
use crate::linear;
use crate::matrix_rep::MatrixRep;
use crate::point_group::PointGroup;
use crate::subgroup::{correlation_table, CorrelationTable};
use crate::thresholds::Thresholds;
use crate::types::{EquivalenceSet, SchoenfliesSymbol, SymmetryCenter};

/// Result of symmetry detection or symmetrization.
#[derive(Debug)]
pub struct SymmetryResult {
    pub group: &'static PointGroup,
    /// 3×3 matrices for each symmetry operation, placed in the molecule's
    /// coordinate frame. Use this (not `group.ops()` alone) for atom permutations
    /// and coordinate transforms.
    pub representation: MatrixRep,
    pub equivalence_sets: Vec<EquivalenceSet>,
    /// Atom positions after processing. Same as input for `detect_symmetry`,
    /// snapped to exact symmetry for `symmetrize`.
    pub centers: Vec<SymmetryCenter>,
}

fn c1_result(centers: &[SymmetryCenter]) -> SymmetryResult {
    let group = PointGroup::c1();
    let equivalence_sets = centers
        .iter()
        .map(|c| EquivalenceSet {
            centers: vec![c.clone()],
        })
        .collect();
    SymmetryResult {
        group,
        representation: MatrixRep::identity_only(group),
        equivalence_sets,
        centers: centers.to_vec(),
    }
}

fn group_from_ctx(ctx: &Context) -> Result<&'static PointGroup, MsymError> {
    PointGroup::from_symbol(ctx.point_group()?)
}

/// Detect point group symmetry of a set of atoms.
///
/// Centers must have positions in Angstroms (libmsym convention).
pub fn detect_symmetry(
    centers: &[SymmetryCenter],
    thresholds: Thresholds,
) -> Result<SymmetryResult, MsymError> {
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

    let group = group_from_ctx(&ctx)?;
    let representation = ctx.symmetry_representation(group)?;
    let equivalence_sets = ctx.equivalence_sets()?;
    let centers = ctx.centers()?;

    Ok(SymmetryResult {
        group,
        representation,
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
) -> Result<SymmetryResult, MsymError> {
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
    let group = group_from_ctx(&ctx)?;
    let representation = ctx.symmetry_representation(group)?;
    let equivalence_sets = ctx.equivalence_sets()?;
    ctx.symmetrize_centers()?;
    let centers = ctx.centers()?;

    Ok(SymmetryResult {
        group,
        representation,
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
pub fn generate_symmetry_images(
    label: SchoenfliesSymbol,
    asymmetric_unit: &[SymmetryCenter],
    thresholds: Thresholds,
) -> Result<SymmetryResult, MsymError> {
    let mut ctx = Context::new()?;
    ctx.set_thresholds(&thresholds)?;
    // Set asymmetric unit to establish center of mass in the context.
    // Then override center of mass to origin so generation works correctly.
    ctx.set_centers(asymmetric_unit)?;
    ctx.set_center_of_mass([0.0, 0.0, 0.0])?;
    ctx.set_point_group(label)?;
    ctx.generate_centers(asymmetric_unit)?;
    ctx.find_symmetry()?;

    let group = group_from_ctx(&ctx)?;
    let representation = ctx.symmetry_representation(group)?;
    let equivalence_sets = ctx.equivalence_sets()?;
    let centers = ctx.centers()?;

    Ok(SymmetryResult {
        group,
        representation,
        equivalence_sets,
        centers,
    })
}

#[derive(Debug)]
pub struct SymmetryDescentResult {
    pub parent_group: &'static PointGroup,
    pub child_group: &'static PointGroup,
    pub child_representation: MatrixRep,
    pub transform: [[f64; 3]; 3],
    pub correlation: Option<CorrelationTable>,
    pub centers: Vec<SymmetryCenter>,
    pub equivalence_sets: Vec<EquivalenceSet>,
}

fn build_correlation(
    parent: &'static PointGroup,
    child: &'static PointGroup,
    class_map: &[usize],
) -> Result<CorrelationTable, MsymError> {
    correlation_table(parent, child, class_map).map_err(|e| MsymError {
        code: umol_msym_sys::MSYM_INVALID_CHARACTER_TABLE,
        message: format!("correlation table computation failed: {e}"),
    })
}

/// Lower the symmetry of a molecule to a specified subgroup.
pub fn lower_symmetry(
    centers: &[SymmetryCenter],
    target: SchoenfliesSymbol,
    thresholds: Thresholds,
) -> Result<SymmetryDescentResult, MsymError> {
    let mut ctx = Context::new()?;
    ctx.set_centers(centers)?;
    ctx.set_thresholds(&thresholds)?;
    ctx.find_symmetry()?;

    let parent_group = group_from_ctx(&ctx)?;
    let centers_out = ctx.centers()?;

    let identity_transform = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    // Identity descent: target is the same group.
    if target == parent_group.symbol() {
        let equivalence_sets = ctx.equivalence_sets()?;
        let child_representation = ctx.symmetry_representation(parent_group)?;
        let correlation = if parent_group.is_linear() {
            None
        } else {
            let n_classes = parent_group.class_reps().len();
            let class_map: Vec<usize> = (0..n_classes).collect();
            Some(build_correlation(parent_group, parent_group, &class_map)?)
        };
        return Ok(SymmetryDescentResult {
            parent_group,
            child_group: parent_group,
            child_representation,
            transform: identity_transform,
            correlation,
            centers: centers_out,
            equivalence_sets,
        });
    }

    // C1 is always a subgroup but libmsym doesn't list it.
    if target == SchoenfliesSymbol::Cn(1) {
        let child_group = PointGroup::c1();
        let class_map = vec![0]; // E → E
        let equivalence_sets = centers_out
            .iter()
            .map(|c| EquivalenceSet {
                centers: vec![c.clone()],
                })
            .collect();
        let correlation = if parent_group.is_linear() {
            None
        } else {
            Some(build_correlation(parent_group, child_group, &class_map)?)
        };
        return Ok(SymmetryDescentResult {
            parent_group,
            child_group,
            child_representation: MatrixRep::identity_only(child_group),
            transform: identity_transform,
            correlation,
            centers: centers_out,
            equivalence_sets,
        });
    }

    // Infinite → finite: re-perceive the molecule under the target finite group.
    // No correlation table (infinite parent has no finite character table).
    if parent_group.is_linear() {
        let child_name = target.to_string();
        let mut ctx2 = Context::new()?;
        ctx2.set_centers(centers)?;
        ctx2.set_thresholds(&thresholds)?;
        ctx2.set_point_group_by_name(&child_name)?;
        ctx2.find_symmetry()?;

        let detected = ctx2.point_group()?;
        if detected != target {
            return Err(MsymError {
                code: umol_msym_sys::MSYM_INVALID_SUBGROUPS,
                message: format!(
                    "{target} is not a subgroup of {}, or the molecule cannot be perceived under it",
                    parent_group.symbol()
                ),
            });
        }

        let child_group = PointGroup::parse(&child_name)?;
        let child_representation = ctx2.symmetry_representation(child_group)?;
        let centers_out = ctx2.centers()?;
        let equivalence_sets = ctx2.equivalence_sets()?;

        return Ok(SymmetryDescentResult {
            parent_group,
            child_group,
            child_representation,
            transform: identity_transform,
            correlation: None,
            centers: centers_out,
            equivalence_sets,
        });
    }

    let subgroups = ctx.subgroups()?;
    let sg = subgroups
        .iter()
        .find(|sg| sg.symbol == target)
        .ok_or_else(|| MsymError {
            code: umol_msym_sys::MSYM_INVALID_SUBGROUPS,
            message: format!("{target} is not a subgroup of {}", parent_group.symbol()),
        })?
        .clone();

    let parent_ops_with_class = sg.parent_ops.clone();

    ctx.select_subgroup(&sg)?;

    let child_name = ctx.point_group_name()?;
    let child_group = PointGroup::parse(&child_name)?;
    let transform = ctx.alignment_transform()?;
    let child_representation = ctx.symmetry_representation(child_group)?;
    let centers_out = ctx.centers()?;
    let equivalence_sets = ctx.equivalence_sets()?;

    // Build class_map: child class → parent class.
    // Match each child op to the closest parent subgroup op by 3×3 matrix.
    let child_n_classes = child_group.class_reps().len();
    let mut class_map = vec![0usize; child_n_classes];
    let mut class_assigned = vec![false; child_n_classes];

    for child_op in child_group.ops() {
        let child_class = child_op.class();
        if class_assigned[child_class] {
            continue;
        }
        let child_matrix = child_representation.matrix(child_op);
        let (_, parent_class) = parent_ops_with_class
            .iter()
            .map(|(pmat, pcls)| {
                let diff = pmat - child_matrix;
                let dist: f64 = diff.iter().map(|x| x * x).sum();
                (dist, *pcls)
            })
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
            .unwrap();
        class_map[child_class] = parent_class;
        class_assigned[child_class] = true;
    }

    let correlation = Some(build_correlation(parent_group, child_group, &class_map)?);

    Ok(SymmetryDescentResult {
        parent_group,
        child_group,
        child_representation,
        transform,
        correlation,
        centers: centers_out,
        equivalence_sets,
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
) -> Result<SalcBasis, MsymError> {
    let mut ctx = Context::new()?;
    ctx.set_centers(centers)?;
    ctx.set_thresholds(&thresholds)?;
    ctx.find_symmetry()?;

    let group = group_from_ctx(&ctx)?;

    if group.is_linear() {
        let aligned = ctx.centers()?;
        return Ok(linear::compute_salcs(
            &aligned,
            basis,
            group,
            thresholds.equivalence,
        ));
    }

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
        irrep_salcs[species_idx as usize].push(Salc {
            coefficients: sparse,
        });
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
#[allow(clippy::too_many_arguments)]
mod tests {
    use std::iter::repeat_n;

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
    // HCl (C∞v)
    #[case(
        &[17, 1],
        &[35.453, 1.008],
        &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.275]],
        "C∞v", 0, 2
    )]
    // CO₂ (D∞h)
    #[case(
        &[8, 6, 8],
        &[15.999, 12.011, 15.999],
        &[[0.0, 0.0, -1.16], [0.0, 0.0, 0.0], [0.0, 0.0, 1.16]],
        "D∞h", 0, 2
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
            &[
                [0.0, 0.0, 0.117],
                [0.0, 0.757, -0.469],
                [0.0, -0.757, -0.469],
            ],
        );
        let result = symmetrize(&centers, Thresholds::default()).unwrap();
        assert_eq!(result.group.to_string(), "C2v");
        assert_eq!(result.centers.len(), 3);

        // H atoms should have exactly opposite y coordinates after symmetrization
        let y1 = result.centers[1].position[1];
        let y2 = result.centers[2].position[1];
        assert!(
            (y1 + y2).abs() < 1e-10,
            "H atoms not symmetric: {y1} vs {y2}"
        );
    }

    #[rstest]
    #[case(
        SchoenfliesSymbol::Cnv(2),
        &[8], &[15.999], &[[0.0, 0.0, 0.117]],
        &[1],  &[1.008], &[[0.0, 0.757, -0.469]],
        "C2v", 3
    )]
    #[case(
        SchoenfliesSymbol::Td,
        &[6], &[12.011], &[[0.0, 0.0, 0.0]],
        &[1], &[1.008],  &[[0.629, 0.629, 0.629]],
        "Td", 5
    )]
    #[case(
        SchoenfliesSymbol::Dnh(6),
        &[6], &[12.011], &[[1.4, 0.0, 0.0]],
        &[1], &[1.008],  &[[2.48, 0.0, 0.0]],
        "D6h", 12
    )]
    fn test_generate_symmetry_images(
        #[case] label: SchoenfliesSymbol,
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
        let result = generate_symmetry_images(label, &asym, Thresholds::default()).unwrap();
        assert_eq!(result.group.to_string(), expected_group);
        assert_eq!(result.centers.len(), expected_atoms);
    }

    fn s_basis(n_atoms: usize) -> Vec<BasisFunction> {
        (0..n_atoms)
            .map(|i| BasisFunction {
                atom_index: i,
                kind: BasisKind::RealSphericalHarmonic,
                shell_index: 1,
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
                    shell_index: 2,
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
    // HCl (C∞v): 2 atoms → 2 Σ+ SALCs
    #[case(
        &[17, 1],
        &[35.453, 1.008],
        &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.275]],
        2, &["Σ+", "Σ+"]
    )]
    // CO₂ (D∞h): 3 atoms → 2 Σ+g + 1 Σ+u
    #[case(
        &[8, 6, 8],
        &[15.999, 12.011, 15.999],
        &[[0.0, 0.0, -1.16], [0.0, 0.0, 0.0], [0.0, 0.0, 1.16]],
        3, &["Σ+g", "Σ+g", "Σ+u"]
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
            .flat_map(|ib| repeat_n(ib.irrep.symbol().to_owned(), ib.salcs.len()))
            .collect();
        salc_symbols.sort();
        let mut expected: Vec<&str> = expected_salc_irreps.to_vec();
        expected.sort();
        assert_eq!(salc_symbols, expected);
    }

    fn assert_salcs_orthonormal(result: &SalcBasis) {
        for ib in &result.irreps {
            for (i, s1) in ib.salcs.iter().enumerate() {
                for (j, s2) in ib.salcs.iter().enumerate() {
                    let dot: f64 = s1
                        .coefficients
                        .iter()
                        .map(|&(k, c1)| {
                            s2.coefficients
                                .iter()
                                .find(|&&(k2, _)| k2 == k)
                                .map_or(0.0, |&(_, c2)| c1 * c2)
                        })
                        .sum();
                    if i == j {
                        assert!((dot - 1.0).abs() < 1e-8, "SALC {i} not normalized: {dot}");
                    } else {
                        assert!(dot.abs() < 1e-8, "SALCs {i},{j} not orthogonal: {dot}");
                    }
                }
            }
        }
    }

    #[rstest]
    // Water (C2v): 3N=9 displacement SALCs
    #[case(
        &[8, 1, 1],
        &[15.999, 1.008, 1.008],
        &[[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
        9
    )]
    // HCl (C∞v): 3N=6 displacement SALCs: 2Σ+ (z) + 2Π (x,y)
    #[case(
        &[17, 1],
        &[35.453, 1.008],
        &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.275]],
        6
    )]
    // CO₂ (D∞h): 3N=9 displacement SALCs
    #[case(
        &[8, 6, 8],
        &[15.999, 12.011, 15.999],
        &[[0.0, 0.0, -1.16], [0.0, 0.0, 0.0], [0.0, 0.0, 1.16]],
        9
    )]
    fn test_compute_salcs_displacement_basis(
        #[case] zs: &[i32],
        #[case] masses: &[f64],
        #[case] positions: &[[f64; 3]],
        #[case] expected_total: usize,
    ) {
        let centers = make_centers(zs, masses, positions);
        let basis = displacement_basis(centers.len());
        let result = compute_salcs(&centers, &basis, Thresholds::default()).unwrap();

        let total: usize = result.irreps.iter().map(|ib| ib.salcs.len()).sum();
        assert_eq!(total, expected_total);
        assert_salcs_orthonormal(&result);
    }
}
