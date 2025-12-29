//! 3D coordinate type.

/// 3D coordinate type (Cartesian coordinates)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3D {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub fn zero() -> Self {
        Self { x: 0.0, y: 0.0, z: 0.0 }
    }

    pub fn is_zero(&self) -> bool {
        self.x.to_bits() == 0 && self.y.to_bits() == 0 && self.z.to_bits() == 0
    }
}

/// Check if all positions are zero
pub fn all_zero(positions: &[Point3D]) -> bool {
    positions.iter().all(|p| p.is_zero())
}