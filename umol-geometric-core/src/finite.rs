//! Finiteness of quantities derived from planar coordinates.
//!
//! Finite coordinates do not guarantee finite derived geometry: the difference of two finite points
//! can overflow, and coincident points have no direction. These predicates report whether a
//! difference, a length, or a direction derived from two points exists as a finite value.

use crate::point::Point2D;

/// The vector from `start` to `end` when both components are finite.
pub fn finite_difference(start: Point2D, end: Point2D) -> Option<Point2D> {
    let difference = Point2D::new(end.x - start.x, end.y - start.y);
    difference.is_finite().then_some(difference)
}

/// The Euclidean length of `vector` when it is finite and strictly positive.
pub fn nonzero_length(vector: Point2D) -> Option<f64> {
    let length = vector.x.hypot(vector.y);
    (length.is_finite() && length > 0.0).then_some(length)
}

/// The unit vector from `start` to `end` when the difference is finite with non-zero length.
pub fn normalizable_direction(start: Point2D, end: Point2D) -> Option<Point2D> {
    let difference = finite_difference(start, end)?;
    let length = nonzero_length(difference)?;
    Some(Point2D::new(difference.x / length, difference.y / length))
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::finite(Point2D::new(1.0, 2.0), Point2D::new(4.0, -2.0), Some(Point2D::new(3.0, -4.0)))]
    #[case::coincident(Point2D::new(1.0, 2.0), Point2D::new(1.0, 2.0), Some(Point2D::zero()))]
    #[case::overflow_x(Point2D::new(-1e308, 0.0), Point2D::new(1e308, 0.0), None)]
    #[case::overflow_y(Point2D::new(0.0, 1e308), Point2D::new(0.0, -1e308), None)]
    fn test_finite_difference(
        #[case] start: Point2D,
        #[case] end: Point2D,
        #[case] expected: Option<Point2D>,
    ) {
        assert_eq!(finite_difference(start, end), expected);
    }

    #[rstest]
    #[case::unit(Point2D::new(3.0, 4.0), Some(5.0))]
    #[case::subnormal(Point2D::new(1e-320, 0.0), Some(1e-320))]
    #[case::zero(Point2D::zero(), None)]
    #[case::negative_zero(Point2D::new(-0.0, 0.0), None)]
    #[case::overflowing_hypot(Point2D::new(1.5e308, 1.5e308), None)]
    #[case::infinite(Point2D::new(f64::INFINITY, 0.0), None)]
    #[case::nan(Point2D::new(f64::NAN, 0.0), None)]
    fn test_nonzero_length(#[case] vector: Point2D, #[case] expected: Option<f64>) {
        assert_eq!(nonzero_length(vector), expected);
    }

    #[rstest]
    #[case::axis(
        Point2D::new(1.0, 1.0),
        Point2D::new(1.0, 4.0),
        Some(Point2D::new(0.0, 1.0))
    )]
    #[case::diagonal(Point2D::zero(), Point2D::new(-3.0, 4.0), Some(Point2D::new(-0.6, 0.8)))]
    #[case::coincident(Point2D::new(2.0, 2.0), Point2D::new(2.0, 2.0), None)]
    #[case::overflow(Point2D::new(1e308, 0.0), Point2D::new(-1e308, 0.0), None)]
    #[case::overflowing_hypot(Point2D::new(-1.5e308, -1.5e308), Point2D::new(0.0, 0.0), None)]
    fn test_normalizable_direction(
        #[case] start: Point2D,
        #[case] end: Point2D,
        #[case] expected: Option<Point2D>,
    ) {
        assert_eq!(normalizable_direction(start, end), expected);
    }
}
