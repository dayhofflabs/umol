//! In-plane direction helpers for reading 2D depictions.

use crate::Point3D;

/// The open in-plane (z = 0) direction at `center`: opposite the mean of the unit
/// vectors to `neighbors`. This is where a virtual ligand sits when a wedge is read
/// against a 2D depiction.
pub fn complementary_direction(center: Point3D, neighbors: &[Point3D]) -> Point3D {
    let (mut sum_x, mut sum_y) = (0.0, 0.0);
    for neighbor in neighbors {
        let (dx, dy) = (neighbor.x - center.x, neighbor.y - center.y);
        let length = (dx * dx + dy * dy).sqrt();
        if length > 0.0 {
            sum_x += dx / length;
            sum_y += dy / length;
        }
    }
    Point3D::new(center.x - sum_x, center.y - sum_y, 0.0)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::two_orthogonal(
        vec![Point3D::new(1.0, 0.0, 0.0), Point3D::new(0.0, 1.0, 0.0)],
        Point3D::new(-1.0, -1.0, 0.0)
    )]
    #[case::opposite_pair_cancels(
        vec![Point3D::new(1.0, 0.0, 0.0), Point3D::new(-1.0, 0.0, 0.0)],
        Point3D::new(0.0, 0.0, 0.0)
    )]
    fn test_complementary_direction(#[case] neighbors: Vec<Point3D>, #[case] expected: Point3D) {
        assert_eq!(
            complementary_direction(Point3D::zero(), &neighbors),
            expected
        );
    }
}
