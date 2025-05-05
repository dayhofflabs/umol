//! Conformer type for the molecular graph model.

use crate::AtomIndex;
use std::collections::HashMap;

/// Type alias for 3D coordinates using nalgebra.
pub type Point3D = nalgebra::Point3<f64>;

/// Represents a single conformation (set of 3D coordinates) for a molecule.
#[derive(Debug, Clone)]
pub struct Conformer {
    /// Optional identifier for the conformer.
    pub id: Option<usize>,
    /// Atomic coordinates for this conformation. The vector index corresponds
    /// to the atom's index *in the graph* (which should align with MOL file order
    /// if constructed correctly).
    pub positions: Vec<Point3D>,
    /// Flag indicating if the conformer represents 3D coordinates.
    pub is_3d: bool,
    /// Generic string-based properties associated with this conformer.
    pub properties: HashMap<String, String>,
}

impl Conformer {
    /// Creates a new Conformer for a given number of atoms,
    /// initializing all positions to the origin.
    pub fn new(num_atoms: usize, is_3d: bool) -> Self {
        Self {
            id: None,
            // Initialize positions with default points (e.g., origin)
            positions: vec![Point3D::origin(); num_atoms],
            is_3d,
            properties: HashMap::new(),
        }
    }

    /// Sets the 3D position for a specific atom graph index within this conformer.
    ///
    /// Panics if the index derived from `idx` is out of bounds
    /// for the `positions` vector. Assumes the conformer was correctly sized by the Molecule.
    pub fn set_position(&mut self, idx: AtomIndex, pos: Point3D) {
        if let Some(p) = self.positions.get_mut(idx.index()) {
            *p = pos;
        } else {
            // This shouldn't happen if Molecule manages conformer size correctly.
            panic!(
                "Attempted to set position for out-of-bounds atom index {} (conformer size {})",
                idx.index(),
                self.positions.len()
            );
        }
    }

    /// Gets the 3D position for a specific atom graph index within this conformer.
    pub fn get_position(&self, idx: AtomIndex) -> Option<&Point3D> {
        self.positions.get(idx.index())
    }
}
