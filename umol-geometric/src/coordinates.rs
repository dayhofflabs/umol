//! Coordinate data structures for the Born-Oppenheimer model.

use nalgebra::DMatrix;

/// Coordinate representation for a molecular geometry.
pub enum Coordinates {
    /// Cartesian coordinates: 3×N matrix, each column is (x, y, z) in Bohr.
    Cartesian(DMatrix<f64>),
}
