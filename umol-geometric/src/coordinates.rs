//! Coordinate data structures for the Born-Oppenheimer model

use nalgebra::DMatrix;

use crate::point_group::PointGroup;

/// Coordinate representation for a molecular geometry.
pub enum Coordinates<G: PointGroup> {
    /// Cartesian coordinates: 3×N matrix, each column is (x, y, z) in Angstroms.
    Cartesian(DMatrix<f64>),

    /// Symmetry-adapted representation (future).
    #[allow(dead_code)]
    Symmetric {
        group: G,
        unique_atoms: Vec<usize>,
        full_coords: DMatrix<f64>,
    },
}
