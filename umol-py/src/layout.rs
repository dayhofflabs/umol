//! Two-dimensional molecule layout bindings, enabled by the `depiction` feature.

use pyo3::exceptions::{PyIndexError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use umol_geometric_core::Point2D;
use umol_graph_ir::ir::AtomId as GraphIrAtomId;
use umol_io::layout::{
    layout_molecule, LayoutError as IoLayoutError, MoleculeLayout as IoMoleculeLayout,
    MoleculeLayoutError as IoMoleculeLayoutError,
};

use crate::depict::MoleculeLayoutAlgorithm;
use crate::molecule::Molecule;

/// An immutable two-dimensional coordinate assignment in a dense atom-id frame.
///
/// Position `i` belongs to atom `i`. The layout carries no chemical attributes; editing returns a
/// new layout and leaves this one unchanged.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq)]
pub struct MoleculeLayout(IoMoleculeLayout);

#[pymethods]
impl MoleculeLayout {
    /// Construct from `(x, y)` positions in ascending atom-id order.
    ///
    /// Raises `ValueError` for a NaN or infinite coordinate.
    #[new]
    fn new(positions: Vec<(f64, f64)>) -> PyResult<Self> {
        IoMoleculeLayout::try_new(
            positions
                .into_iter()
                .map(|(x, y)| Point2D::new(x, y))
                .collect(),
        )
        .map(Self)
        .map_err(molecule_layout_error)
    }

    /// `(x, y)` positions in ascending atom-id order.
    #[getter]
    fn positions(&self) -> Vec<(f64, f64)> {
        self.0
            .positions()
            .iter()
            .map(|position| (position.x, position.y))
            .collect()
    }

    /// Return the `(x, y)` position of `atom_id`.
    ///
    /// Raises `IndexError` if `atom_id` is outside the layout frame.
    fn position(&self, atom_id: u32) -> PyResult<(f64, f64)> {
        self.0
            .position(GraphIrAtomId(atom_id))
            .map(|position| (position.x, position.y))
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

    /// Return a new layout with `atom_id` at `position`; this layout is unchanged.
    ///
    /// Raises `ValueError` if `atom_id` is outside the layout frame or `position` is not finite.
    fn with_position(&self, atom_id: u32, position: (f64, f64)) -> PyResult<Self> {
        let mut layout = self.0.clone();
        layout
            .set_position(GraphIrAtomId(atom_id), Point2D::new(position.0, position.1))
            .map_err(molecule_layout_error)?;
        Ok(Self(layout))
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

#[pymethods]
impl Molecule {
    /// Generate a two-dimensional layout of this molecule with `algorithm`.
    ///
    /// Raises `RuntimeError` if the layout backend fails.
    #[pyo3(signature = (*, algorithm=MoleculeLayoutAlgorithm::CoordGen()))]
    fn layout(&self, algorithm: MoleculeLayoutAlgorithm) -> PyResult<MoleculeLayout> {
        layout_molecule(self.to_rust(), algorithm.to_rust())
            .map(MoleculeLayout::from_rust)
            .map_err(layout_error)
    }
}

fn molecule_layout_error(error: IoMoleculeLayoutError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn layout_error(error: IoLayoutError) -> PyErr {
    PyRuntimeError::new_err(error.to_string())
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::{PyIndexError, PyValueError};
    use rstest::rstest;
    use umol_graph_ir::mol_dsl;
    use umol_io::layout::MoleculeLayoutAlgorithm as IoMoleculeLayoutAlgorithm;

    use super::*;

    #[rstest]
    fn test_molecule_layout_new() {
        let layout = MoleculeLayout::new(vec![(0.0, 1.0), (2.0, -1.0)]).unwrap();

        assert_eq!(layout.positions(), vec![(0.0, 1.0), (2.0, -1.0)]);
        assert_eq!(layout.__len__(), 2);
        assert_eq!(
            layout.to_rust(),
            &IoMoleculeLayout::try_new(vec![Point2D::new(0.0, 1.0), Point2D::new(2.0, -1.0)])
                .unwrap()
        );
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
    fn test_molecule_layout_new_error(#[case] positions: Vec<(f64, f64)>, #[case] message: &str) {
        Python::attach(|py| {
            let error = MoleculeLayout::new(positions).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                message
            );
        });
    }

    #[rstest]
    fn test_molecule_layout_position() {
        let layout = MoleculeLayout::new(vec![(0.0, 1.0), (2.0, -1.0)]).unwrap();

        assert_eq!(layout.position(1).unwrap(), (2.0, -1.0));
    }

    #[rstest]
    fn test_molecule_layout_position_error() {
        Python::attach(|py| {
            let layout = MoleculeLayout::new(vec![(0.0, 1.0), (2.0, -1.0)]).unwrap();

            let error = layout.position(2).unwrap_err();

            assert!(error.is_instance_of::<PyIndexError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "atom 2 is outside layout frame of size 2"
            );
        });
    }

    #[rstest]
    fn test_molecule_layout_with_position() {
        let layout = MoleculeLayout::new(vec![(0.0, 1.0), (2.0, -1.0)]).unwrap();

        let moved = layout.with_position(1, (3.5, 4.0)).unwrap();

        assert_eq!(
            moved,
            MoleculeLayout::new(vec![(0.0, 1.0), (3.5, 4.0)]).unwrap()
        );
        assert_eq!(
            layout,
            MoleculeLayout::new(vec![(0.0, 1.0), (2.0, -1.0)]).unwrap()
        );
    }

    #[rstest]
    #[case::frame(2, (0.0, 0.0), "atom 2 is outside layout frame of size 2")]
    #[case::nan(
        0,
        (f64::NAN, 0.0),
        "atom 0 has non-finite position Point2D { x: NaN, y: 0.0 }"
    )]
    fn test_molecule_layout_with_position_error(
        #[case] atom_id: u32,
        #[case] position: (f64, f64),
        #[case] message: &str,
    ) {
        Python::attach(|py| {
            let layout = MoleculeLayout::new(vec![(0.0, 1.0), (2.0, -1.0)]).unwrap();

            let error = layout.with_position(atom_id, position).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                message
            );
        });
    }

    #[rstest]
    fn test_molecule_layout() {
        let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#));
        let expected =
            layout_molecule(molecule.to_rust(), IoMoleculeLayoutAlgorithm::CoordGen).unwrap();

        let layout = molecule
            .layout(MoleculeLayoutAlgorithm::CoordGen())
            .unwrap();

        assert_eq!(layout.to_rust(), &expected);
    }
}
