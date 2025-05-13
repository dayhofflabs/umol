//! Conformer parsing functions for MOL files.

use crate::Point3D;

/// Check if the conformer is 3D based on the presence of non-zero coordinates.
pub(crate) fn is_3d(positions: &[Point3D]) -> bool {
    positions
        .iter()
        .any(|pos| pos.x * pos.x + pos.y * pos.y + pos.z * pos.z > f64::EPSILON)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_3d() {
        let positions = vec![Point3D::new(0.0, 0.0, 0.0)];
        assert!(!is_3d(&positions));

        let positions = vec![Point3D::new(1.0, 0.0, 0.0)];
        assert!(is_3d(&positions));
    }
}
