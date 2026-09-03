//! Two-dimensional coordinate assignments coupled to graph-IR atom frames.

use std::any::Any;

use thiserror::Error;
use umol_geometric_core::Point2D;
use umol_graph_ir::ir::{AtomId, Molecule};
use umol_utils::error::UmolError;

#[cfg(feature = "coordgen")]
mod coordgen;

/// Algorithm used to generate a two-dimensional molecule layout.
///
/// This selector has no default: callers choose the operational backend explicitly.
#[cfg(feature = "coordgen")]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MoleculeLayoutAlgorithm {
    /// The vendored CoordGen 2D coordinate generator.
    CoordGen,
}

/// Generates a two-dimensional layout for `molecule` with the selected algorithm.
///
/// The result preserves the supplied dense [`AtomId`] frame. Backend projection uses localized
/// topology and literal element and bond-order hints, with generic fallback values for forms the
/// backend cannot represent. It does not canonicalize or change the graph-IR molecule.
///
/// # Errors
///
/// Returns [`LayoutError::CoordGen`] if the selected backend cannot generate coordinates.
#[cfg(feature = "coordgen")]
pub fn layout_molecule(
    molecule: &Molecule,
    algorithm: MoleculeLayoutAlgorithm,
) -> Result<MoleculeLayout, LayoutError> {
    match algorithm {
        MoleculeLayoutAlgorithm::CoordGen => coordgen::layout(molecule).map_err(LayoutError::from),
    }
}

/// Failure while generating a molecule layout.
#[cfg(feature = "coordgen")]
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum LayoutError {
    #[error("CoordGen layout failed: {0}")]
    CoordGen(#[from] umol_coordgen_sys::CoordgenError),
}

#[cfg(feature = "coordgen")]
impl UmolError for LayoutError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// An editable two-dimensional coordinate assignment in a dense [`AtomId`] frame.
///
/// Position `i` belongs to `AtomId(i)`. The layout carries no chemical attributes and does not
/// canonicalize or otherwise change the supplied atom frame.
#[derive(Clone, Debug, PartialEq)]
pub struct MoleculeLayout {
    positions: Vec<Point2D>,
}

impl MoleculeLayout {
    /// Constructs a layout from positions in ascending atom-id order.
    ///
    /// # Errors
    ///
    /// Returns [`MoleculeLayoutError::NonFinitePosition`] for the first point containing a NaN or
    /// infinity.
    pub fn try_new(positions: Vec<Point2D>) -> Result<Self, MoleculeLayoutError> {
        if let Some((index, position)) = positions
            .iter()
            .copied()
            .enumerate()
            .find(|(_, position)| !position.is_finite())
        {
            return Err(MoleculeLayoutError::NonFinitePosition {
                atom_id: AtomId::from(index),
                position,
            });
        }

        Ok(Self { positions })
    }

    /// Number of atoms in this layout's frame.
    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// Whether this layout has an empty atom frame.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Positions in ascending atom-id order.
    pub fn positions(&self) -> &[Point2D] {
        &self.positions
    }

    /// Returns the position assigned to `atom_id`, if it belongs to this layout's frame.
    pub fn position(&self, atom_id: AtomId) -> Option<&Point2D> {
        self.positions.get(atom_id.index())
    }

    /// Replaces the position assigned to `atom_id`.
    ///
    /// A failed edit leaves the layout unchanged.
    ///
    /// # Errors
    ///
    /// Returns [`MoleculeLayoutError::AtomOutOfFrame`] if `atom_id` is outside the layout frame,
    /// or [`MoleculeLayoutError::NonFinitePosition`] if `position` contains a NaN or infinity.
    pub fn set_position(
        &mut self,
        atom_id: AtomId,
        position: Point2D,
    ) -> Result<(), MoleculeLayoutError> {
        let frame_size = self.positions.len();
        let stored =
            self.positions
                .get_mut(atom_id.index())
                .ok_or(MoleculeLayoutError::AtomOutOfFrame {
                    atom_id,
                    frame_size,
                })?;
        if !position.is_finite() {
            return Err(MoleculeLayoutError::NonFinitePosition { atom_id, position });
        }
        *stored = position;
        Ok(())
    }

    /// Checks whether `molecule` uses the same dense atom-frame size as this layout.
    ///
    /// This check establishes only frame agreement. It does not validate molecular chemistry or
    /// interpret the coordinates.
    ///
    /// # Errors
    ///
    /// Returns [`MoleculeLayoutError::FrameSizeMismatch`] when the atom counts differ.
    pub fn check_frame(&self, molecule: &Molecule) -> Result<(), MoleculeLayoutError> {
        let molecule_atom_count = molecule.atoms().count();
        let layout_atom_count = self.positions.len();
        if molecule_atom_count != layout_atom_count {
            return Err(MoleculeLayoutError::FrameSizeMismatch {
                molecule_atom_count,
                layout_atom_count,
            });
        }
        Ok(())
    }
}

/// Failures while constructing, editing, or contextually combining a [`MoleculeLayout`].
#[derive(Clone, Debug, Error, PartialEq)]
pub enum MoleculeLayoutError {
    #[error("atom {atom_id} has non-finite position {position:?}")]
    NonFinitePosition { atom_id: AtomId, position: Point2D },
    #[error("atom {atom_id} is outside layout frame of size {frame_size}")]
    AtomOutOfFrame { atom_id: AtomId, frame_size: usize },
    #[error(
        "molecule atom count {molecule_atom_count} does not match layout atom count {layout_atom_count}"
    )]
    FrameSizeMismatch {
        molecule_atom_count: usize,
        layout_atom_count: usize,
    },
}

impl UmolError for MoleculeLayoutError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}
