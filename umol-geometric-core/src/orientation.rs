//! Orientation (handedness) of ordered point tuples.

use crate::point::Point3D;

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

/// Relative tolerance below which a component perpendicular to an axis counts as absent: the
/// perpendicular length is compared with `AXIS_SIDE_TOLERANCE` times the axis length.
pub const AXIS_SIDE_TOLERANCE: f64 = 1e-6;

/// Whether `first` and `second` lie on the same side of the line through `axis_start` and
/// `axis_end`, judged by their components perpendicular to it. `None` when the axis has no
/// length or either perpendicular component is within [`AXIS_SIDE_TOLERANCE`] of the axis
/// length, so all-zero and collinear positions are never read as a side. One formula serves
/// planar and spatial coordinates.
pub fn same_side_of_axis(
    axis_start: Point3D,
    axis_end: Point3D,
    first: Point3D,
    second: Point3D,
) -> Option<bool> {
    let axis = [
        axis_end.x - axis_start.x,
        axis_end.y - axis_start.y,
        axis_end.z - axis_start.z,
    ];
    let axis_norm_sq = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
    if axis_norm_sq == 0.0 {
        return None;
    }
    let perpendicular = |point: Point3D| {
        let r = [
            point.x - axis_start.x,
            point.y - axis_start.y,
            point.z - axis_start.z,
        ];
        let t = (r[0] * axis[0] + r[1] * axis[1] + r[2] * axis[2]) / axis_norm_sq;
        [r[0] - t * axis[0], r[1] - t * axis[1], r[2] - t * axis[2]]
    };
    let p = perpendicular(first);
    let q = perpendicular(second);
    let limit = AXIS_SIDE_TOLERANCE * AXIS_SIDE_TOLERANCE * axis_norm_sq;
    if p[0] * p[0] + p[1] * p[1] + p[2] * p[2] <= limit
        || q[0] * q[0] + q[1] * q[1] + q[2] * q[2] <= limit
    {
        return None;
    }
    Some(p[0] * q[0] + p[1] * q[1] + p[2] * q[2] > 0.0)
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

    #[rstest]
    #[case::same_side_planar(Point3D::new(0.3, 1.0, 0.0), Point3D::new(0.7, 0.5, 0.0), Some(true))]
    #[case::opposite_side_planar(Point3D::new(0.3, 1.0, 0.0), Point3D::new(0.7, -0.5, 0.0), Some(false))]
    #[case::same_side_spatial(Point3D::new(0.3, 0.0, 1.0), Point3D::new(0.7, 0.4, 0.9), Some(true))]
    #[case::opposite_side_spatial(Point3D::new(0.3, 0.0, 1.0), Point3D::new(0.7, 0.0, -1.0), Some(false))]
    #[case::perpendicular_in_space(
        Point3D::new(0.3, 1.0, 0.0),
        Point3D::new(0.7, 0.0, 1.0),
        Some(false)
    )]
    #[case::first_on_axis(Point3D::new(0.3, 0.0, 0.0), Point3D::new(0.7, 1.0, 0.0), None)]
    #[case::second_on_axis(Point3D::new(0.3, 1.0, 0.0), Point3D::new(1.5, 0.0, 0.0), None)]
    #[case::at_tolerance(Point3D::new(0.3, 1e-6, 0.0), Point3D::new(0.7, 1.0, 0.0), None)]
    #[case::above_tolerance(Point3D::new(0.3, 2e-6, 0.0), Point3D::new(0.7, 1.0, 0.0), Some(true))]
    fn test_same_side_of_axis(
        #[case] first: Point3D,
        #[case] second: Point3D,
        #[case] expected: Option<bool>,
    ) {
        assert_eq!(
            same_side_of_axis(Point3D::zero(), Point3D::new(1.0, 0.0, 0.0), first, second),
            expected
        );
    }

    #[rstest]
    #[case::zero_length_axis(Point3D::zero(), Point3D::zero())]
    #[case::all_zero_positions(Point3D::zero(), Point3D::zero())]
    fn test_same_side_of_axis_degenerate_axis(
        #[case] axis_start: Point3D,
        #[case] axis_end: Point3D,
    ) {
        assert_eq!(
            same_side_of_axis(
                axis_start,
                axis_end,
                Point3D::new(0.0, 1.0, 0.0),
                Point3D::zero()
            ),
            None
        );
    }
}
