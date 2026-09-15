//! Explicit layout and SVG depiction bindings, enabled by the `depiction` feature.

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use umol_io::depict::{
    Depict as IoDepict, DepictConfig as IoDepictConfig, Depiction as IoDepiction,
    MoleculeDepictionError as IoMoleculeDepictionError,
    ReactionDepictionError as IoReactionDepictionError,
};
use umol_io::layout::MoleculeLayoutAlgorithm as IoMoleculeLayoutAlgorithm;

use crate::error::contradiction_error;
use crate::layout::{MoleculeLayout, ReactionLayout};
use crate::molecule::Molecule;
use crate::reaction::Reaction;

/// Algorithm used to generate two-dimensional molecule layouts for depiction.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MoleculeLayoutAlgorithm {
    CoordGen(),
}

#[pymethods]
impl MoleculeLayoutAlgorithm {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __repr__(&self) -> &'static str {
        match self {
            Self::CoordGen() => "MoleculeLayoutAlgorithm.CoordGen()",
        }
    }
}

impl MoleculeLayoutAlgorithm {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for molecule-layout algorithms"
    )]
    pub(crate) fn from_rust(algorithm: IoMoleculeLayoutAlgorithm) -> Self {
        match algorithm {
            IoMoleculeLayoutAlgorithm::CoordGen => Self::CoordGen(),
        }
    }

    pub(crate) fn to_rust(self) -> IoMoleculeLayoutAlgorithm {
        match self {
            Self::CoordGen() => IoMoleculeLayoutAlgorithm::CoordGen,
        }
    }
}

/// Operational configuration for molecule and reaction depiction.
///
/// The default configuration selects CoordGen, currently the only layout algorithm.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DepictConfig {
    layout_algorithm: MoleculeLayoutAlgorithm,
}

impl Default for DepictConfig {
    fn default() -> Self {
        Self {
            layout_algorithm: MoleculeLayoutAlgorithm::CoordGen(),
        }
    }
}

#[pymethods]
impl DepictConfig {
    #[new]
    #[pyo3(signature = (*, layout_algorithm=MoleculeLayoutAlgorithm::CoordGen()))]
    fn new(layout_algorithm: MoleculeLayoutAlgorithm) -> Self {
        Self { layout_algorithm }
    }

    #[staticmethod]
    fn default() -> Self {
        Default::default()
    }

    #[getter]
    fn layout_algorithm(&self) -> MoleculeLayoutAlgorithm {
        self.layout_algorithm
    }

    fn __repr__(&self) -> &'static str {
        "DepictConfig.default()"
    }
}

impl DepictConfig {
    pub(crate) fn to_rust(self) -> IoDepictConfig {
        IoDepictConfig {
            layout_algorithm: self.layout_algorithm.to_rust(),
        }
    }
}

/// An opaque, format-neutral molecule or reaction depiction.
#[pyclass(frozen, skip_from_py_object)]
pub struct Depiction(IoDepiction);

#[pymethods]
impl Depiction {
    /// Render this depiction as a complete SVG document suitable for writing to an SVG file.
    fn render_svg(&self) -> String {
        self.0.render_svg()
    }

    /// Return the complete SVG document through Jupyter's rich-display protocol.
    fn _repr_svg_(&self) -> String {
        self.render_svg()
    }
}

impl Depiction {
    pub(crate) fn from_rust(depiction: IoDepiction) -> Self {
        Self(depiction)
    }
}

#[pymethods]
impl Molecule {
    /// Generate a two-dimensional layout of this molecule with `algorithm`.
    ///
    /// Raises `RuntimeError` if the layout backend fails.
    #[pyo3(signature = (*, algorithm=MoleculeLayoutAlgorithm::CoordGen()))]
    fn layout(&self, algorithm: MoleculeLayoutAlgorithm) -> PyResult<MoleculeLayout> {
        let config = IoDepictConfig {
            layout_algorithm: algorithm.to_rust(),
        };
        self.to_rust()
            .layout_with(&config)
            .map(MoleculeLayout::from_rust)
            .map_err(generated_molecule_error)
    }

    /// Construct a format-neutral depiction using the default configuration.
    ///
    /// Raises `RuntimeError` if layout or tetrahedral depiction fails.
    fn depict(&self) -> PyResult<Depiction> {
        self.to_rust()
            .depict()
            .map(Depiction::from_rust)
            .map_err(generated_molecule_error)
    }

    /// Construct a format-neutral depiction using `config`.
    ///
    /// Raises `RuntimeError` if layout or tetrahedral depiction fails.
    fn depict_with(&self, config: DepictConfig) -> PyResult<Depiction> {
        let config = config.to_rust();
        self.to_rust()
            .depict_with(&config)
            .map(Depiction::from_rust)
            .map_err(generated_molecule_error)
    }

    /// Check that `layout` can depict this molecule without moving any coordinate.
    ///
    /// Raises `ValueError` naming the first entity whose supplied geometry cannot be depicted:
    /// a frame mismatch, a definite cis/trans bond drawn in the other configuration or with a
    /// ligand on its axis, non-finite derived geometry, or a tetrahedral centre admitting no wedge.
    fn verify_layout(&self, layout: &MoleculeLayout) -> PyResult<()> {
        self.to_rust()
            .verify_layout(layout.to_rust())
            .map_err(supplied_molecule_error)
    }

    /// Construct a format-neutral depiction in a supplied `layout`, verifying it first.
    ///
    /// Raises `ValueError` as `verify_layout` does.
    fn depict_layout(&self, layout: &MoleculeLayout) -> PyResult<Depiction> {
        self.to_rust()
            .depict_layout(layout.to_rust())
            .map(Depiction::from_rust)
            .map_err(supplied_molecule_error)
    }
}

#[pymethods]
impl Reaction {
    /// Generate a two-dimensional layout of this reaction with `algorithm`.
    ///
    /// Raises `ContradictionError` if the reaction cannot be materialized and `RuntimeError` if
    /// layout of either materialized side or their arrangement fails.
    #[pyo3(signature = (*, algorithm=MoleculeLayoutAlgorithm::CoordGen()))]
    fn layout(
        &self,
        py: Python<'_>,
        algorithm: MoleculeLayoutAlgorithm,
    ) -> PyResult<ReactionLayout> {
        let config = IoDepictConfig {
            layout_algorithm: algorithm.to_rust(),
        };
        self.to_rust(py)?
            .layout_with(&config)
            .map(ReactionLayout::from_rust)
            .map_err(generated_reaction_error)
    }

    /// Construct a format-neutral depiction using the default configuration.
    ///
    /// Raises `ContradictionError` if the reaction cannot be materialized and `RuntimeError` if
    /// layout or depiction of either materialized side fails.
    fn depict(&self, py: Python<'_>) -> PyResult<Depiction> {
        self.to_rust(py)?
            .depict()
            .map(Depiction::from_rust)
            .map_err(generated_reaction_error)
    }

    /// Construct a format-neutral depiction using `config`.
    ///
    /// Raises `ContradictionError` if the reaction cannot be materialized and `RuntimeError` if
    /// layout or depiction of either materialized side fails.
    fn depict_with(&self, py: Python<'_>, config: DepictConfig) -> PyResult<Depiction> {
        let config = config.to_rust();
        self.to_rust(py)?
            .depict_with(&config)
            .map(Depiction::from_rust)
            .map_err(generated_reaction_error)
    }

    /// Check that `layout` can depict this reaction without moving any coordinate.
    ///
    /// Raises `ContradictionError` if the reaction cannot be materialized and `ValueError` if
    /// either side rejects its layout as `Molecule.verify_layout` does or the arrow cannot be
    /// drawn.
    fn verify_layout(&self, py: Python<'_>, layout: &ReactionLayout) -> PyResult<()> {
        self.to_rust(py)?
            .verify_layout(layout.to_rust())
            .map_err(supplied_reaction_error)
    }

    /// Construct a format-neutral depiction in a supplied `layout`, verifying it first.
    ///
    /// Raises `ContradictionError` and `ValueError` as `verify_layout` does.
    fn depict_layout(&self, py: Python<'_>, layout: &ReactionLayout) -> PyResult<Depiction> {
        self.to_rust(py)?
            .depict_layout(layout.to_rust())
            .map(Depiction::from_rust)
            .map_err(supplied_reaction_error)
    }
}

/// Maps a failure of `layout`, `depict`, or `depict_with`: the caller supplied no coordinates, so
/// every failure is operational.
fn generated_molecule_error(error: IoMoleculeDepictionError) -> PyErr {
    PyRuntimeError::new_err(error.to_string())
}

/// Maps a failure of `verify_layout` or `depict_layout`: the supplied coordinates do not fit the
/// molecule, except for a backend failure, which the supplied path cannot produce.
fn supplied_molecule_error(error: IoMoleculeDepictionError) -> PyErr {
    match error {
        error @ IoMoleculeDepictionError::Layout(_) => PyRuntimeError::new_err(error.to_string()),
        error => PyValueError::new_err(error.to_string()),
    }
}

fn generated_reaction_error(error: IoReactionDepictionError) -> PyErr {
    match error {
        IoReactionDepictionError::Materialization(error) => contradiction_error(error),
        error => PyRuntimeError::new_err(error.to_string()),
    }
}

fn supplied_reaction_error(error: IoReactionDepictionError) -> PyErr {
    match error {
        IoReactionDepictionError::Materialization(error) => contradiction_error(error),
        error @ (IoReactionDepictionError::LhsDepiction(IoMoleculeDepictionError::Layout(_))
        | IoReactionDepictionError::RhsDepiction(IoMoleculeDepictionError::Layout(_))) => {
            PyRuntimeError::new_err(error.to_string())
        }
        error => PyValueError::new_err(error.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::{PyRuntimeError, PyTypeError, PyValueError};
    use rstest::rstest;
    use umol_geometric_core::Point2D;
    use umol_graph_ir::ir::{
        AtomId, BondDelta, BondFieldChange, BondId, Delta, Deltas, Molecule as GraphIrMolecule,
        NumForm, Reaction as GraphIrReaction, StereoAtomId,
    };
    use umol_graph_ir::mol_dsl;
    use umol_io::layout::{
        MoleculeLayout as IoMoleculeLayout, MoleculeLayoutError as IoMoleculeLayoutError,
        ReactionLayout as IoReactionLayout, ReactionLayoutError as IoReactionLayoutError,
    };

    use super::*;
    use crate::error::ContradictionError;

    const TRANS_BUTENE: &str = r#"{:atoms ["C" "C" "C" "C"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
        :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#;
    const TETRAHEDRAL: &str = r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#;

    fn message(py: Python<'_>, error: &PyErr) -> String {
        error.value(py).str().unwrap().extract::<String>().unwrap()
    }

    fn layout(positions: &[(f64, f64)]) -> MoleculeLayout {
        MoleculeLayout::from_rust(
            IoMoleculeLayout::try_new(positions.iter().map(|&(x, y)| Point2D::new(x, y)).collect())
                .unwrap(),
        )
    }

    fn contradictory_reaction(py: Python<'_>) -> Reaction {
        let reaction = GraphIrReaction::new(
            mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(2),
                    new: NumForm::Lit(3),
                },
            })]),
        );
        Reaction::from_rust(py, reaction).unwrap()
    }

    #[rstest]
    fn test_molecule_layout_algorithm_conversion() {
        let rust = IoMoleculeLayoutAlgorithm::CoordGen;
        let python = MoleculeLayoutAlgorithm::from_rust(rust);

        assert_eq!(python.to_rust(), rust);
        assert_eq!(python.__repr__(), "MoleculeLayoutAlgorithm.CoordGen()");
    }

    #[rstest]
    fn test_depict_config_new() {
        let config = DepictConfig::new(MoleculeLayoutAlgorithm::CoordGen());

        assert_eq!(config, DepictConfig::default());
        assert_eq!(
            config.layout_algorithm(),
            MoleculeLayoutAlgorithm::CoordGen()
        );
        assert_eq!(config.__repr__(), "DepictConfig.default()");
        assert_eq!(
            config.to_rust().layout_algorithm,
            IoMoleculeLayoutAlgorithm::CoordGen
        );
    }

    #[rstest]
    fn test_depiction_render_svg() {
        let molecule = Molecule::from_rust(GraphIrMolecule::new());
        let rust = molecule.to_rust().depict().unwrap();
        let expected = rust.render_svg();
        let depiction = Depiction::from_rust(rust);

        assert_eq!(depiction.render_svg(), expected);
        assert_eq!(depiction._repr_svg_(), expected);
    }

    #[rstest]
    fn test_depiction_constructor_error() {
        Python::attach(|py| {
            let error = py.get_type::<Depiction>().call0().unwrap_err();

            assert!(error.is_instance_of::<PyTypeError>(py));
        });
    }

    #[rstest]
    fn test_molecule_layout() {
        let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#));
        let expected = molecule.to_rust().layout().unwrap();

        let layout = molecule
            .layout(MoleculeLayoutAlgorithm::CoordGen())
            .unwrap();

        assert_eq!(layout.to_rust(), &expected);
    }

    #[rstest]
    fn test_molecule_depict() {
        let molecule = Molecule::from_rust(GraphIrMolecule::new());
        let expected = molecule.to_rust().depict().unwrap().render_svg();

        let depiction = molecule.depict().unwrap();

        assert_eq!(depiction.render_svg(), expected);
        assert_eq!(depiction._repr_svg_(), expected);
        assert_eq!(
            depiction.render_svg(),
            r#"<svg xmlns="http://www.w3.org/2000/svg" class="umol-depiction" viewBox="-0.5 -0.5 1 1">
</svg>"#
        );
    }

    #[rstest]
    fn test_molecule_depict_with() {
        let molecule = Molecule::from_rust(GraphIrMolecule::new());
        let config = DepictConfig::default();
        let expected = molecule
            .to_rust()
            .depict_with(&config.to_rust())
            .unwrap()
            .render_svg();

        let depiction = molecule.depict_with(config).unwrap();

        assert_eq!(depiction.render_svg(), expected);
    }

    #[rstest]
    fn test_molecule_verify_layout() {
        let molecule = Molecule::from_rust(mol_dsl!(TRANS_BUTENE));
        let generated = molecule
            .layout(MoleculeLayoutAlgorithm::CoordGen())
            .unwrap();

        assert!(molecule.verify_layout(&generated).is_ok());
    }

    #[rstest]
    fn test_molecule_verify_layout_cis_trans_error() {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(TRANS_BUTENE));
            let drawn_as_cis = layout(&[(0.0, 1.0), (1.0, 0.0), (2.0, 0.0), (3.0, 1.0)]);

            let error = molecule.verify_layout(&drawn_as_cis).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                message(py, &error),
                "bond 1 is stored as E but its supplied coordinates draw Z"
            );
        });
    }

    #[rstest]
    fn test_molecule_depict_layout() {
        let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#));
        let generated = molecule
            .layout(MoleculeLayoutAlgorithm::CoordGen())
            .unwrap();
        let expected = molecule.to_rust().depict().unwrap().render_svg();

        let depiction = molecule.depict_layout(&generated).unwrap();

        assert_eq!(depiction.render_svg(), expected);
    }

    #[rstest]
    fn test_molecule_depict_layout_frame_error() {
        Python::attach(|py| {
            let molecule =
                Molecule::from_rust(mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#));

            let error = molecule
                .depict_layout(&layout(&[(0.0, 0.0)]))
                .err()
                .unwrap();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                message(py, &error),
                "layout frame: molecule atom count 2 does not match layout atom count 1"
            );
        });
    }

    #[rstest]
    fn test_molecule_depict_layout_tetrahedral_error() {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(TETRAHEDRAL));
            let collinear = layout(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)]);

            let error = molecule.depict_layout(&collinear).err().unwrap();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                message(py, &error),
                "tetrahedral geometry cannot establish a display wedge for stereo atom 0"
            );
        });
    }

    #[rstest]
    fn test_generated_molecule_error() {
        Python::attach(|py| {
            let error = generated_molecule_error(IoMoleculeDepictionError::TetrahedralGeometry {
                stereo_atom: StereoAtomId(3),
            });

            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                message(py, &error),
                "tetrahedral geometry cannot establish a display wedge for stereo atom 3"
            );
        });
    }

    #[rstest]
    #[case::layout_frame(
        IoMoleculeDepictionError::LayoutFrame(IoMoleculeLayoutError::FrameSizeMismatch {
            molecule_atom_count: 2,
            layout_atom_count: 1,
        }),
        "layout frame: molecule atom count 2 does not match layout atom count 1"
    )]
    #[case::tetrahedral_geometry(
        IoMoleculeDepictionError::TetrahedralGeometry {
            stereo_atom: StereoAtomId(3),
        },
        "tetrahedral geometry cannot establish a display wedge for stereo atom 3"
    )]
    fn test_supplied_molecule_error(
        #[case] input: IoMoleculeDepictionError,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let error = supplied_molecule_error(input);

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(message(py, &error), expected);
        });
    }

    #[rstest]
    #[case::lhs(IoReactionDepictionError::LhsDepiction(
        IoMoleculeDepictionError::TetrahedralGeometry {
            stereo_atom: StereoAtomId(1),
        }
    ), "lhs depiction: tetrahedral geometry cannot establish a display wedge for stereo atom 1")]
    #[case::rhs(IoReactionDepictionError::RhsDepiction(
        IoMoleculeDepictionError::TetrahedralGeometry {
            stereo_atom: StereoAtomId(2),
        }
    ), "rhs depiction: tetrahedral geometry cannot establish a display wedge for stereo atom 2")]
    #[case::layout(IoReactionDepictionError::Layout(
        IoReactionLayoutError::LhsTranslation(IoMoleculeLayoutError::NonFinitePosition {
            atom_id: AtomId(0),
            position: Point2D::new(f64::INFINITY, 0.0),
        })
    ), "reaction layout: lhs translation: atom 0 has non-finite position Point2D { x: inf, y: 0.0 }")]
    fn test_generated_reaction_error(
        #[case] input: IoReactionDepictionError,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let error = generated_reaction_error(input);

            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(message(py, &error), expected);
        });
    }

    #[rstest]
    #[case::lhs(IoReactionDepictionError::LhsDepiction(
        IoMoleculeDepictionError::TetrahedralGeometry {
            stereo_atom: StereoAtomId(1),
        }
    ), "lhs depiction: tetrahedral geometry cannot establish a display wedge for stereo atom 1")]
    #[case::rhs(IoReactionDepictionError::RhsDepiction(
        IoMoleculeDepictionError::LayoutFrame(IoMoleculeLayoutError::FrameSizeMismatch {
            molecule_atom_count: 2,
            layout_atom_count: 1,
        })
    ), "rhs depiction: layout frame: molecule atom count 2 does not match layout atom count 1")]
    #[case::layout(IoReactionDepictionError::Layout(
        IoReactionLayoutError::DegenerateArrow {
            position: Point2D::new(0.0, 0.0),
        }
    ), "reaction layout: reaction arrow starts and ends at Point2D { x: 0.0, y: 0.0 }")]
    fn test_supplied_reaction_error(
        #[case] input: IoReactionDepictionError,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let error = supplied_reaction_error(input);

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(message(py, &error), expected);
        });
    }

    #[rstest]
    fn test_reaction_layout() {
        Python::attach(|py| {
            let reaction = GraphIrReaction::new(
                mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
                Deltas::new(),
            );
            let expected = reaction.layout().unwrap();
            let reaction = Reaction::from_rust(py, reaction).unwrap();

            let layout = reaction
                .layout(py, MoleculeLayoutAlgorithm::CoordGen())
                .unwrap();

            assert_eq!(layout.to_rust(), &expected);
        });
    }

    #[rstest]
    fn test_reaction_layout_error() {
        Python::attach(|py| {
            let reaction = contradictory_reaction(py);

            let error = reaction
                .layout(py, MoleculeLayoutAlgorithm::CoordGen())
                .unwrap_err();

            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(message(py, &error), "reached a contradiction");
        });
    }

    #[rstest]
    fn test_reaction_depict() {
        Python::attach(|py| {
            let reaction =
                GraphIrReaction::new(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#), Deltas::new());
            let expected = reaction.depict().unwrap().render_svg();
            let reaction = Reaction::from_rust(py, reaction).unwrap();

            let depiction = reaction.depict(py).unwrap();

            assert_eq!(depiction.render_svg(), expected);
            assert_eq!(depiction._repr_svg_(), expected);
        });
    }

    #[rstest]
    fn test_reaction_depict_with() {
        Python::attach(|py| {
            let reaction =
                GraphIrReaction::new(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#), Deltas::new());
            let config = DepictConfig::default();
            let expected = reaction
                .depict_with(&config.to_rust())
                .unwrap()
                .render_svg();
            let reaction = Reaction::from_rust(py, reaction).unwrap();

            let depiction = reaction.depict_with(py, config).unwrap();

            assert_eq!(depiction.render_svg(), expected);
        });
    }

    #[rstest]
    fn test_reaction_depict_with_error() {
        Python::attach(|py| {
            let reaction = contradictory_reaction(py);

            let error = reaction
                .depict_with(py, DepictConfig::default())
                .err()
                .unwrap();

            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(message(py, &error), "reached a contradiction");
        });
    }

    #[rstest]
    fn test_reaction_verify_and_depict_layout() {
        Python::attach(|py| {
            let reaction = GraphIrReaction::new(
                mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
                Deltas::new(),
            );
            let expected = reaction.depict().unwrap().render_svg();
            let reaction = Reaction::from_rust(py, reaction).unwrap();
            let generated = reaction
                .layout(py, MoleculeLayoutAlgorithm::CoordGen())
                .unwrap();

            assert!(reaction.verify_layout(py, &generated).is_ok());
            assert_eq!(
                reaction.depict_layout(py, &generated).unwrap().render_svg(),
                expected
            );
        });
    }

    #[rstest]
    fn test_reaction_depict_layout_frame_error() {
        Python::attach(|py| {
            let reaction = GraphIrReaction::new(
                mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
                Deltas::new(),
            );
            let reaction = Reaction::from_rust(py, reaction).unwrap();
            let side = layout(&[(0.0, 0.0)]);
            let mismatched = ReactionLayout::from_rust(
                IoReactionLayout::try_new(
                    side.to_rust().clone(),
                    side.to_rust().clone(),
                    Point2D::new(-1.0, 0.0),
                    Point2D::new(1.0, 0.0),
                )
                .unwrap(),
            );

            let error = reaction.depict_layout(py, &mismatched).err().unwrap();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                message(py, &error),
                "lhs depiction: layout frame: molecule atom count 2 does not match layout atom count 1"
            );
        });
    }
}
