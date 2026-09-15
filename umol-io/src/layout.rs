//! Two-dimensional coordinate assignments coupled to graph-IR atom frames.

use std::any::Any;

use thiserror::Error;
use umol_geometric_core::Point2D;
use umol_graph_ir::ir::{AtomId, Molecule};
use umol_utils::error::UmolError;

#[cfg(feature = "coordgen")]
mod coordgen;
#[cfg(feature = "coordgen")]
pub(crate) mod stereo;

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
pub(crate) fn layout_molecule(
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

const ARROW_HALF_LENGTH: f64 = 0.75;
const SIDE_ARROW_GAP: f64 = 1.0;

/// Two side layouts and a reaction arrow in one reaction coordinate system.
///
/// Every side position and both arrow endpoints share one coordinate system. Like
/// [`MoleculeLayout`], the reaction layout carries no chemical attributes and is not bound to a
/// reaction: its side frames are checked against the materialized sides only when it is depicted.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionLayout {
    lhs: MoleculeLayout,
    rhs: MoleculeLayout,
    arrow_start: Point2D,
    arrow_end: Point2D,
}

impl ReactionLayout {
    /// Constructs a layout from two side layouts and explicit arrow endpoints.
    ///
    /// # Errors
    ///
    /// Returns [`ReactionLayoutError::NonFiniteArrow`] if an arrow endpoint contains a NaN or
    /// infinity, and [`ReactionLayoutError::DegenerateArrow`] if the endpoints coincide.
    pub fn try_new(
        lhs: MoleculeLayout,
        rhs: MoleculeLayout,
        arrow_start: Point2D,
        arrow_end: Point2D,
    ) -> Result<Self, ReactionLayoutError> {
        check_arrow(arrow_start, arrow_end)?;
        Ok(Self {
            lhs,
            rhs,
            arrow_start,
            arrow_end,
        })
    }

    /// Places `lhs` left of a horizontal arrow at the origin and `rhs` right of it.
    ///
    /// Each side is translated so that it is vertically centered on the arrow and separated from
    /// the arrow tip by a one-bond gap; the arrow runs from `(-0.75, 0)` to `(0.75, 0)`. An empty
    /// side is left where it is.
    ///
    /// # Errors
    ///
    /// Returns [`ReactionLayoutError::LhsTranslation`] or [`ReactionLayoutError::RhsTranslation`]
    /// if translating a side produces a non-finite position.
    pub fn arrange(lhs: MoleculeLayout, rhs: MoleculeLayout) -> Result<Self, ReactionLayoutError> {
        let lhs = translate_layout(&lhs, side_offset(&lhs, Side::Lhs))
            .map_err(ReactionLayoutError::LhsTranslation)?;
        let rhs = translate_layout(&rhs, side_offset(&rhs, Side::Rhs))
            .map_err(ReactionLayoutError::RhsTranslation)?;
        Ok(Self {
            lhs,
            rhs,
            arrow_start: Point2D::new(-ARROW_HALF_LENGTH, 0.0),
            arrow_end: Point2D::new(ARROW_HALF_LENGTH, 0.0),
        })
    }

    /// The left-hand side layout.
    pub fn lhs(&self) -> &MoleculeLayout {
        &self.lhs
    }

    /// The right-hand side layout.
    pub fn rhs(&self) -> &MoleculeLayout {
        &self.rhs
    }

    /// Mutable access to the left-hand side layout.
    pub fn lhs_mut(&mut self) -> &mut MoleculeLayout {
        &mut self.lhs
    }

    /// Mutable access to the right-hand side layout.
    pub fn rhs_mut(&mut self) -> &mut MoleculeLayout {
        &mut self.rhs
    }

    /// Tail of the reaction arrow.
    pub fn arrow_start(&self) -> Point2D {
        self.arrow_start
    }

    /// Tip of the reaction arrow.
    pub fn arrow_end(&self) -> Point2D {
        self.arrow_end
    }

    /// Replaces both arrow endpoints.
    ///
    /// A failed edit leaves the layout unchanged.
    ///
    /// # Errors
    ///
    /// Returns [`ReactionLayoutError::NonFiniteArrow`] if an endpoint contains a NaN or infinity,
    /// and [`ReactionLayoutError::DegenerateArrow`] if the endpoints coincide.
    pub fn set_arrow(&mut self, start: Point2D, end: Point2D) -> Result<(), ReactionLayoutError> {
        check_arrow(start, end)?;
        self.arrow_start = start;
        self.arrow_end = end;
        Ok(())
    }
}

pub(crate) fn check_arrow(start: Point2D, end: Point2D) -> Result<(), ReactionLayoutError> {
    if !start.is_finite() || !end.is_finite() {
        return Err(ReactionLayoutError::NonFiniteArrow { start, end });
    }
    if start == end {
        return Err(ReactionLayoutError::DegenerateArrow { position: start });
    }
    Ok(())
}

#[derive(Clone, Copy)]
enum Side {
    Lhs,
    Rhs,
}

fn side_offset(layout: &MoleculeLayout, side: Side) -> Point2D {
    let Some(first) = layout.positions().first() else {
        return Point2D::new(0.0, 0.0);
    };
    let mut min = *first;
    let mut max = *first;
    for &position in &layout.positions()[1..] {
        min.x = min.x.min(position.x);
        min.y = min.y.min(position.y);
        max.x = max.x.max(position.x);
        max.y = max.y.max(position.y);
    }
    let x = match side {
        Side::Lhs => -ARROW_HALF_LENGTH - SIDE_ARROW_GAP - max.x,
        Side::Rhs => ARROW_HALF_LENGTH + SIDE_ARROW_GAP - min.x,
    };
    Point2D::new(x, -(min.y + max.y) / 2.0)
}

fn translate_layout(
    layout: &MoleculeLayout,
    offset: Point2D,
) -> Result<MoleculeLayout, MoleculeLayoutError> {
    MoleculeLayout::try_new(
        layout
            .positions()
            .iter()
            .map(|position| Point2D::new(position.x + offset.x, position.y + offset.y))
            .collect(),
    )
}

/// Failures while constructing or editing a [`ReactionLayout`].
#[derive(Clone, Debug, Error, PartialEq)]
pub enum ReactionLayoutError {
    #[error("reaction arrow from {start:?} to {end:?} has a non-finite endpoint")]
    NonFiniteArrow { start: Point2D, end: Point2D },
    #[error("reaction arrow starts and ends at {position:?}")]
    DegenerateArrow { position: Point2D },
    #[error("lhs translation: {0}")]
    LhsTranslation(#[source] MoleculeLayoutError),
    #[error("rhs translation: {0}")]
    RhsTranslation(#[source] MoleculeLayoutError),
}

impl UmolError for ReactionLayoutError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}
