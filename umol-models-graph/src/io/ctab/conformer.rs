//! Conformer type for CTab format.

/// Type alias for 3D coordinates using nalgebra.
pub type Point3D = nalgebra::Point3<f64>;

/// Represents a single conformation (set of 3D coordinates) for a molecule.
#[derive(Debug, Clone)]
pub struct Conformer {
    pub positions: Vec<Point3D>,
}

impl Conformer {
    /// Create new conformer with all positions at the origin
    pub fn new(num_atoms: usize) -> Self {
        Self {
            positions: vec![Point3D::origin(); num_atoms],
        }
    }

    /// Create new conformer from given positions
    pub fn from_positions(positions: Vec<Point3D>) -> Self {
        Self { positions }
    }

    /// Get number of atoms
    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// Set 3D position by index
    ///
    /// Panic if the index is out of bounds
    pub fn set_position(&mut self, idx: usize, pos: Point3D) {
        if let Some(p) = self.positions.get_mut(idx) {
            *p = pos;
        } else {
            panic!(
                "Attempted to set position for out-of-bounds atom index {} (conformer size {})",
                idx,
                self.positions.len()
            );
        }
    }

    /// Get the 3D position by index
    pub fn get_position(&self, idx: usize) -> Option<&Point3D> {
        self.positions.get(idx)
    }
}
