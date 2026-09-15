//! Two-dimensional molecule and reaction layout bindings, enabled by the `depiction` feature.

use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use umol_geometric_core::Point2D;
use umol_graph_ir::ir::AtomId as GraphIrAtomId;
use umol_io::layout::{
    MoleculeLayout as IoMoleculeLayout, MoleculeLayoutError as IoMoleculeLayoutError,
    ReactionLayout as IoReactionLayout, ReactionLayoutError as IoReactionLayoutError,
};

/// A two-dimensional coordinate assignment in a dense atom-id frame.
///
/// Position `i` belongs to atom `i`. The layout carries no chemical attributes and is not bound
/// to a molecule: `set_position` accepts any finite in-frame point, and the coordinates are
/// checked against a molecule only by `Molecule.verify_layout` and `Molecule.depict_layout`.
/// Mutable, value-equal, and unhashable.
#[pyclass(eq, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq)]
pub struct MoleculeLayout(IoMoleculeLayout);

#[pymethods]
impl MoleculeLayout {
    /// Construct from `(x, y)` positions in ascending atom-id order.
    ///
    /// Raises `ValueError` for a NaN or infinite coordinate.
    #[new]
    fn new(positions: Vec<(f64, f64)>) -> PyResult<Self> {
        IoMoleculeLayout::try_new(positions.into_iter().map(point).collect())
            .map(Self)
            .map_err(molecule_layout_error)
    }

    /// `(x, y)` positions in ascending atom-id order, as a tuple.
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.0.positions().iter().copied().map(coordinates))
    }

    /// Return the `(x, y)` position of `atom_id`.
    ///
    /// Raises `IndexError` if `atom_id` is outside the layout frame.
    fn position(&self, atom_id: u32) -> PyResult<(f64, f64)> {
        self.0
            .position(GraphIrAtomId(atom_id))
            .copied()
            .map(coordinates)
            .ok_or_else(|| {
                PyIndexError::new_err(
                    IoMoleculeLayoutError::AtomOutOfFrame {
                        atom_id: GraphIrAtomId(atom_id),
                        frame_size: self.0.atom_count(),
                    }
                    .to_string(),
                )
            })
    }

    /// Move `atom_id` to `position` in place.
    ///
    /// Raises `ValueError` if `atom_id` is outside the layout frame or `position` is not finite;
    /// the layout is unchanged on failure.
    fn set_position(&mut self, atom_id: u32, position: (f64, f64)) -> PyResult<()> {
        self.0
            .set_position(GraphIrAtomId(atom_id), point(position))
            .map_err(molecule_layout_error)
    }

    fn __len__(&self) -> usize {
        self.0.atom_count()
    }
}

impl MoleculeLayout {
    pub(crate) fn from_rust(layout: IoMoleculeLayout) -> Self {
        Self(layout)
    }

    pub(crate) fn to_rust(&self) -> &IoMoleculeLayout {
        &self.0
    }
}

/// Two side layouts and a reaction arrow in one reaction coordinate system.
///
/// `lhs` and `rhs` return copies; a side position is edited through `set_lhs_position` or
/// `set_rhs_position`, so every write goes to this layout and no read aliases it. Like
/// `MoleculeLayout`, the reaction layout is not bound to a reaction: its side frames are checked
/// against the materialized sides only by `Reaction.verify_layout` and `Reaction.depict_layout`.
/// Mutable, value-equal, and unhashable.
#[pyclass(eq, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionLayout(IoReactionLayout);

#[pymethods]
impl ReactionLayout {
    /// Construct from two side layouts and `(x, y)` arrow endpoints.
    ///
    /// Raises `ValueError` if an arrow endpoint is not finite or the endpoints coincide.
    #[new]
    fn new(
        lhs: &MoleculeLayout,
        rhs: &MoleculeLayout,
        arrow_start: (f64, f64),
        arrow_end: (f64, f64),
    ) -> PyResult<Self> {
        IoReactionLayout::try_new(
            lhs.0.clone(),
            rhs.0.clone(),
            point(arrow_start),
            point(arrow_end),
        )
        .map(Self)
        .map_err(reaction_layout_error)
    }

    /// Place `lhs` left of a horizontal arrow at the origin and `rhs` right of it.
    ///
    /// Raises `ValueError` if translating a side produces a non-finite position.
    #[staticmethod]
    fn arrange(lhs: &MoleculeLayout, rhs: &MoleculeLayout) -> PyResult<Self> {
        IoReactionLayout::arrange(lhs.0.clone(), rhs.0.clone())
            .map(Self)
            .map_err(reaction_layout_error)
    }

    /// A copy of the left-hand side layout.
    #[getter]
    fn lhs(&self) -> MoleculeLayout {
        MoleculeLayout(self.0.lhs().clone())
    }

    /// A copy of the right-hand side layout.
    #[getter]
    fn rhs(&self) -> MoleculeLayout {
        MoleculeLayout(self.0.rhs().clone())
    }

    /// Move left-hand side atom `atom_id` to `position` in place.
    ///
    /// Raises `ValueError` if `atom_id` is outside the side's frame or `position` is not finite;
    /// the layout is unchanged on failure.
    fn set_lhs_position(&mut self, atom_id: u32, position: (f64, f64)) -> PyResult<()> {
        self.0
            .lhs_mut()
            .set_position(GraphIrAtomId(atom_id), point(position))
            .map_err(molecule_layout_error)
    }

    /// Move right-hand side atom `atom_id` to `position` in place.
    ///
    /// Raises `ValueError` if `atom_id` is outside the side's frame or `position` is not finite;
    /// the layout is unchanged on failure.
    fn set_rhs_position(&mut self, atom_id: u32, position: (f64, f64)) -> PyResult<()> {
        self.0
            .rhs_mut()
            .set_position(GraphIrAtomId(atom_id), point(position))
            .map_err(molecule_layout_error)
    }

    /// `(x, y)` tail of the reaction arrow.
    #[getter]
    fn arrow_start(&self) -> (f64, f64) {
        coordinates(self.0.arrow_start())
    }

    /// `(x, y)` tip of the reaction arrow.
    #[getter]
    fn arrow_end(&self) -> (f64, f64) {
        coordinates(self.0.arrow_end())
    }

    /// Replace both arrow endpoints in place.
    ///
    /// Raises `ValueError` if an endpoint is not finite or the endpoints coincide; the layout is
    /// unchanged on failure.
    fn set_arrow(&mut self, start: (f64, f64), end: (f64, f64)) -> PyResult<()> {
        self.0
            .set_arrow(point(start), point(end))
            .map_err(reaction_layout_error)
    }
}

impl ReactionLayout {
    pub(crate) fn from_rust(layout: IoReactionLayout) -> Self {
        Self(layout)
    }

    pub(crate) fn to_rust(&self) -> &IoReactionLayout {
        &self.0
    }
}

fn point((x, y): (f64, f64)) -> Point2D {
    Point2D::new(x, y)
}

fn coordinates(position: Point2D) -> (f64, f64) {
    (position.x, position.y)
}

fn molecule_layout_error(error: IoMoleculeLayoutError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn reaction_layout_error(error: IoReactionLayoutError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::{PyIndexError, PyValueError};
    use rstest::rstest;

    use super::*;

    fn message(py: Python<'_>, error: &PyErr) -> String {
        error.value(py).str().unwrap().extract::<String>().unwrap()
    }

    fn layout(positions: Vec<(f64, f64)>) -> MoleculeLayout {
        MoleculeLayout::new(positions).unwrap()
    }

    #[rstest]
    fn test_molecule_layout_new() {
        Python::attach(|py| {
            let layout = layout(vec![(0.0, 1.0), (2.0, -1.0)]);

            let positions = layout.positions(py).unwrap();

            assert_eq!(
                positions.extract::<Vec<(f64, f64)>>().unwrap(),
                vec![(0.0, 1.0), (2.0, -1.0)]
            );
            assert_eq!(layout.__len__(), 2);
            assert_eq!(
                layout.to_rust(),
                &IoMoleculeLayout::try_new(vec![Point2D::new(0.0, 1.0), Point2D::new(2.0, -1.0)])
                    .unwrap()
            );
        });
    }

    #[rstest]
    #[case::nan(
        vec![(0.0, 0.0), (f64::NAN, 1.0)],
        "atom 1 has non-finite position Point2D { x: NaN, y: 1.0 }"
    )]
    #[case::infinity(
        vec![(f64::INFINITY, 0.0)],
        "atom 0 has non-finite position Point2D { x: inf, y: 0.0 }"
    )]
    fn test_molecule_layout_new_error(#[case] positions: Vec<(f64, f64)>, #[case] expected: &str) {
        Python::attach(|py| {
            let error = MoleculeLayout::new(positions).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(message(py, &error), expected);
        });
    }

    #[rstest]
    fn test_molecule_layout_position() {
        let layout = layout(vec![(0.0, 1.0), (2.0, -1.0)]);

        assert_eq!(layout.position(1).unwrap(), (2.0, -1.0));
    }

    #[rstest]
    fn test_molecule_layout_position_error() {
        Python::attach(|py| {
            let layout = layout(vec![(0.0, 1.0), (2.0, -1.0)]);

            let error = layout.position(2).unwrap_err();

            assert!(error.is_instance_of::<PyIndexError>(py));
            assert_eq!(
                message(py, &error),
                "atom 2 is outside layout frame of size 2"
            );
        });
    }

    #[rstest]
    fn test_molecule_layout_set_position() {
        let mut moved = layout(vec![(0.0, 1.0), (2.0, -1.0)]);

        moved.set_position(1, (3.5, 4.0)).unwrap();

        assert_eq!(moved, layout(vec![(0.0, 1.0), (3.5, 4.0)]));
    }

    #[rstest]
    #[case::frame(2, (0.0, 0.0), "atom 2 is outside layout frame of size 2")]
    #[case::nan(
        0,
        (f64::NAN, 0.0),
        "atom 0 has non-finite position Point2D { x: NaN, y: 0.0 }"
    )]
    fn test_molecule_layout_set_position_error(
        #[case] atom_id: u32,
        #[case] position: (f64, f64),
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let mut unchanged = layout(vec![(0.0, 1.0), (2.0, -1.0)]);

            let error = unchanged.set_position(atom_id, position).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(message(py, &error), expected);
            assert_eq!(unchanged, layout(vec![(0.0, 1.0), (2.0, -1.0)]));
        });
    }

    #[rstest]
    fn test_reaction_layout_new() {
        let lhs = layout(vec![(0.0, 0.0), (1.0, 0.0)]);
        let rhs = layout(vec![(4.0, 0.0)]);

        let reaction = ReactionLayout::new(&lhs, &rhs, (2.0, 0.0), (3.0, 0.0)).unwrap();

        assert_eq!(reaction.lhs(), lhs);
        assert_eq!(reaction.rhs(), rhs);
        assert_eq!(reaction.arrow_start(), (2.0, 0.0));
        assert_eq!(reaction.arrow_end(), (3.0, 0.0));
        assert_eq!(
            reaction.to_rust(),
            &IoReactionLayout::try_new(
                lhs.to_rust().clone(),
                rhs.to_rust().clone(),
                Point2D::new(2.0, 0.0),
                Point2D::new(3.0, 0.0),
            )
            .unwrap()
        );
    }

    #[rstest]
    #[case::degenerate(
        (1.0, 1.0),
        (1.0, 1.0),
        "reaction arrow starts and ends at Point2D { x: 1.0, y: 1.0 }"
    )]
    #[case::non_finite(
        (0.0, 0.0),
        (f64::NAN, 0.0),
        "reaction arrow from Point2D { x: 0.0, y: 0.0 } to Point2D { x: NaN, y: 0.0 } has a non-finite endpoint"
    )]
    fn test_reaction_layout_new_error(
        #[case] start: (f64, f64),
        #[case] end: (f64, f64),
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let side = layout(vec![(0.0, 0.0)]);

            let error = ReactionLayout::new(&side, &side, start, end).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(message(py, &error), expected);
        });
    }

    #[rstest]
    fn test_reaction_layout_arrange() {
        let lhs = layout(vec![(0.0, 0.0), (1.0, 0.0)]);
        let rhs = layout(vec![(4.0, 0.0)]);

        let reaction = ReactionLayout::arrange(&lhs, &rhs).unwrap();

        assert_eq!(
            reaction.to_rust(),
            &IoReactionLayout::arrange(lhs.to_rust().clone(), rhs.to_rust().clone()).unwrap()
        );
        assert_eq!(reaction.arrow_start(), (-0.75, 0.0));
        assert_eq!(reaction.arrow_end(), (0.75, 0.0));
    }

    #[rstest]
    fn test_reaction_layout_side_copies() {
        let side = layout(vec![(0.0, 0.0), (1.0, 0.0)]);
        let mut reaction = ReactionLayout::new(&side, &side, (2.0, 0.0), (3.0, 0.0)).unwrap();

        let mut lhs = reaction.lhs();
        lhs.set_position(0, (5.0, 5.0)).unwrap();

        assert_eq!(reaction.lhs(), side);

        reaction.set_lhs_position(0, (5.0, 5.0)).unwrap();
        reaction.set_rhs_position(1, (6.0, 6.0)).unwrap();

        assert_eq!(reaction.lhs(), lhs);
        assert_eq!(reaction.rhs(), layout(vec![(0.0, 0.0), (6.0, 6.0)]));
    }

    #[rstest]
    fn test_reaction_layout_set_side_position_error() {
        Python::attach(|py| {
            let side = layout(vec![(0.0, 0.0)]);
            let mut reaction = ReactionLayout::new(&side, &side, (2.0, 0.0), (3.0, 0.0)).unwrap();
            let unchanged = reaction.clone();

            let error = reaction.set_rhs_position(1, (0.0, 0.0)).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                message(py, &error),
                "atom 1 is outside layout frame of size 1"
            );
            assert_eq!(reaction, unchanged);
        });
    }

    #[rstest]
    fn test_reaction_layout_set_arrow() {
        Python::attach(|py| {
            let side = layout(vec![(0.0, 0.0)]);
            let mut reaction = ReactionLayout::new(&side, &side, (2.0, 0.0), (3.0, 0.0)).unwrap();

            reaction.set_arrow((0.0, -1.0), (0.0, 1.0)).unwrap();

            assert_eq!(reaction.arrow_start(), (0.0, -1.0));
            assert_eq!(reaction.arrow_end(), (0.0, 1.0));

            let error = reaction.set_arrow((0.0, 0.0), (0.0, 0.0)).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(reaction.arrow_start(), (0.0, -1.0));
            assert_eq!(reaction.arrow_end(), (0.0, 1.0));
        });
    }
}
