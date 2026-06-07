//! Orientation (handedness) of ordered point tuples.

use crate::Point3D;

/// Signed volume of the tetrahedron `(a, b, c, d)` — the 3×3 determinant of the
/// edge vectors from `a`, i.e. six times the geometric volume. Its **sign** is the
/// orientation / handedness of the ordered tuple; callers take the sign.
pub fn signed_volume(a: Point3D, b: Point3D, c: Point3D, d: Point3D) -> f64 {
    let u = [b.x - a.x, b.y - a.y, b.z - a.z];
    let v = [c.x - a.x, c.y - a.y, c.z - a.z];
    let w = [d.x - a.x, d.y - a.y, d.z - a.z];
    u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0])
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::right_handed(
        Point3D::new(1.0, 0.0, 0.0),
        Point3D::new(0.0, 1.0, 0.0),
        Point3D::new(0.0, 0.0, 1.0),
        1.0
    )]
    #[case::left_handed(
        Point3D::new(0.0, 1.0, 0.0),
        Point3D::new(1.0, 0.0, 0.0),
        Point3D::new(0.0, 0.0, 1.0),
        -1.0
    )]
    #[case::degenerate_coplanar(
        Point3D::new(1.0, 0.0, 0.0),
        Point3D::new(0.0, 1.0, 0.0),
        Point3D::new(1.0, 1.0, 0.0),
        0.0
    )]
    fn test_signed_volume(
        #[case] b: Point3D,
        #[case] c: Point3D,
        #[case] d: Point3D,
        #[case] expected: f64,
    ) {
        assert_eq!(signed_volume(Point3D::zero(), b, c, d), expected);
    }
}
