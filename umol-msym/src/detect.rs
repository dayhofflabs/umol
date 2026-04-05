use crate::context::Context;
use crate::error::Error;
use crate::point_group::PointGroup;
use crate::types::{EquivalenceSet, SymmetryCenter, Thresholds};

/// Result of symmetry detection or symmetrization.
#[derive(Debug)]
pub struct SymmetryResult {
    pub group: &'static PointGroup,
    pub equivalence_sets: Vec<EquivalenceSet>,
    /// Atom positions after processing. Same as input for `detect_symmetry`,
    /// snapped to exact symmetry for `symmetrize`.
    pub centers: Vec<SymmetryCenter>,
}

/// Detect point group symmetry of a set of atoms.
///
/// Centers must have positions in Angstroms (libmsym convention).
pub fn detect_symmetry(
    centers: &[SymmetryCenter],
    thresholds: Thresholds,
) -> Result<SymmetryResult, Error> {
    let mut ctx = Context::new()?;
    ctx.set_elements(centers)?;
    ctx.set_thresholds(&thresholds)?;
    ctx.find_symmetry()?;

    let group = PointGroup::from_context(&ctx)?;
    let equivalence_sets = ctx.equivalence_sets()?;
    let centers = ctx.elements()?;

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
    ctx.set_elements(centers)?;
    ctx.set_thresholds(&thresholds)?;
    ctx.find_symmetry()?;
    ctx.symmetrize_elements()?;
    // Re-detect to repopulate the character table (cleared by symmetrize_elements)
    ctx.find_symmetry()?;

    let group = PointGroup::from_context(&ctx)?;
    let equivalence_sets = ctx.equivalence_sets()?;
    let centers = ctx.elements()?;

    Ok(SymmetryResult {
        group,
        equivalence_sets,
        centers,
    })
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

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
}
