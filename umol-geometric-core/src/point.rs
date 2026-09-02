/// Dimensionless 2D and 3D coordinate types

// 2D coordinate type (Cartesian coordinates)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

impl Point2D {
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub const fn zero() -> Self {
        Self { x: 0.0, y: 0.0 }
    }

    pub const fn is_finite(&self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }

    pub const fn is_zero(&self) -> bool {
        self.x.to_bits() == 0 && self.y.to_bits() == 0
    }
}

/// 3D coordinate type (Cartesian coordinates)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3D {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub const fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    pub const fn is_zero(&self) -> bool {
        self.x.to_bits() == 0 && self.y.to_bits() == 0 && self.z.to_bits() == 0
    }

    pub const fn is_finite(&self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    pub fn all_zero(positions: &[Point3D]) -> bool {
        positions.iter().all(|p| p.is_zero())
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::finite(Point2D::new(1.0, 2.0), true)]
    #[case::infinite(Point2D::new(f64::INFINITY, 2.0), false)]
    #[case::nan(Point2D::new(f64::NAN, 2.0), false)]
    fn test_point2d_is_finite(#[case] point: Point2D, #[case] expected: bool) {
        assert_eq!(point.is_finite(), expected);
    }

    #[rstest]
    #[case::zero(Point2D::zero(), true)]
    #[case::negative_zero(Point2D::new(-0.0, 0.0), false)]
    #[case::nonzero(Point2D::new(0.0, 1.0), false)]
    fn test_point2d_is_zero(#[case] point: Point2D, #[case] expected: bool) {
        assert_eq!(point.is_zero(), expected);
    }

    #[rstest]
    #[case::zero(Point3D::zero(), true)]
    #[case::negative_zero(Point3D::new(-0.0, 0.0, 0.0), false)]
    #[case::nonzero(Point3D::new(0.0, 1.0, 0.0), false)]
    fn test_point3d_is_zero(#[case] point: Point3D, #[case] expected: bool) {
        assert_eq!(point.is_zero(), expected);
    }

    #[rstest]
    #[case::finite(Point3D::new(1.0, 2.0, 3.0), true)]
    #[case::infinite(Point3D::new(f64::INFINITY, 2.0, 3.0), false)]
    #[case::nan(Point3D::new(f64::NAN, 2.0, 3.0), false)]
    fn test_point3d_is_finite(#[case] point: Point3D, #[case] expected: bool) {
        assert_eq!(point.is_finite(), expected);
    }

    #[rstest]
    #[case::zero(vec![Point3D::zero()], true)]
    #[case::negative_zero(vec![Point3D::new(-0.0, 0.0, 0.0)], false)]
    #[case::nonzero(vec![Point3D::new(0.0, 1.0, 0.0)], false)]
    fn test_point3d_all_zero(#[case] positions: Vec<Point3D>, #[case] expected: bool) {
        assert_eq!(Point3D::all_zero(&positions), expected);
    }
}
