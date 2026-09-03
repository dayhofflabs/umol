//! Explicit layout and SVG display bindings.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use umol_io::depict::molecule::MoleculeDepictionError as IoMoleculeDepictionError;
use umol_io::depict::reaction::ReactionDepictionError as IoReactionDepictionError;
use umol_io::depict::{
    Depict as IoDepict, DepictConfig as IoDepictConfig, Depiction as IoDepiction,
};
use umol_io::layout::MoleculeLayoutAlgorithm as IoMoleculeLayoutAlgorithm;

use crate::error::contradiction_error;
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
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion boundary for depiction configuration"
    )]
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
    /// Render this depiction as a complete SVG document fragment.
    fn render_svg(&self) -> String {
        self.0.render_svg()
    }

    /// Return the complete SVG fragment through Jupyter's rich-display protocol.
    fn _repr_svg_(&self) -> String {
        self.render_svg()
    }
}

impl Depiction {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion boundary for format-neutral depictions"
    )]
    pub(crate) fn from_rust(depiction: IoDepiction) -> Self {
        Self(depiction)
    }
}

/// An already rendered SVG value for notebook display.
#[pyclass(frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Svg(String);

#[pymethods]
impl Svg {
    /// Return the complete SVG fragment through Jupyter's rich-display protocol.
    fn _repr_svg_(&self) -> &str {
        &self.0
    }
}

impl Svg {
    pub(crate) fn from_rust(depiction: &IoDepiction) -> Self {
        Self(depiction.render_svg())
    }
}

#[pymethods]
impl Molecule {
    /// Generate and render an SVG depiction with an explicitly selected layout algorithm.
    fn depict_with(&self, layout_algorithm: MoleculeLayoutAlgorithm) -> PyResult<Svg> {
        self.to_rust()
            .depict_with(layout_algorithm.to_rust())
            .map(|depiction| Svg::from_rust(&depiction))
            .map_err(molecule_depiction_error)
    }
}

#[pymethods]
impl Reaction {
    /// Generate and render an SVG depiction with an explicitly selected layout algorithm.
    fn depict_with(
        &self,
        py: Python<'_>,
        layout_algorithm: MoleculeLayoutAlgorithm,
    ) -> PyResult<Svg> {
        self.to_rust(py)?
            .depict_with(layout_algorithm.to_rust())
            .map(|depiction| Svg::from_rust(&depiction))
            .map_err(reaction_depiction_error)
    }
}

fn molecule_depiction_error(error: IoMoleculeDepictionError) -> PyErr {
    PyRuntimeError::new_err(error.to_string())
}

fn reaction_depiction_error(error: IoReactionDepictionError) -> PyErr {
    match error {
        IoReactionDepictionError::Materialization(error) => contradiction_error(error),
        error @ (IoReactionDepictionError::LhsDepiction(_)
        | IoReactionDepictionError::RhsDepiction(_)) => PyRuntimeError::new_err(error.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::{PyRuntimeError, PyTypeError};
    use rstest::rstest;
    use umol_graph_ir::ir::{
        BondDelta, BondFieldChange, BondId, Delta, Deltas, Molecule as GraphIrMolecule, NumForm,
        Reaction as GraphIrReaction, StereoAtomId,
    };

    use super::*;
    use crate::error::ContradictionError;

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
    fn test_svg_repr_svg() {
        let molecule = Molecule::from_rust(GraphIrMolecule::new());
        let depiction = molecule
            .to_rust()
            .depict_with(IoMoleculeLayoutAlgorithm::CoordGen)
            .unwrap();
        let expected = depiction.render_svg();

        let svg = molecule
            .depict_with(MoleculeLayoutAlgorithm::CoordGen())
            .unwrap();

        assert_eq!(svg._repr_svg_(), expected);
        assert_eq!(
            svg._repr_svg_(),
            r#"<svg xmlns="http://www.w3.org/2000/svg" class="umol-depiction" viewBox="-0.5 -0.5 1 1">
</svg>"#
        );
    }

    #[rstest]
    fn test_svg_constructor_error() {
        Python::attach(|py| {
            let error = py.get_type::<Svg>().call0().unwrap_err();

            assert!(error.is_instance_of::<PyTypeError>(py));
        });
    }

    #[rstest]
    fn test_molecule_depiction_error() {
        Python::attach(|py| {
            let error = molecule_depiction_error(IoMoleculeDepictionError::TetrahedralGeometry {
                stereo_atom: StereoAtomId(3),
            });

            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "tetrahedral geometry cannot establish a display wedge for stereo atom 3"
            );
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
    fn test_reaction_depiction_error(
        #[case] input: IoReactionDepictionError,
        #[case] message: &str,
    ) {
        Python::attach(|py| {
            let error = reaction_depiction_error(input);

            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                message
            );
        });
    }

    #[rstest]
    fn test_reaction_depict_with() {
        Python::attach(|py| {
            let reaction = GraphIrReaction::new(
                umol_graph_ir::mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
                Deltas::new(),
            );
            let expected = reaction
                .depict_with(IoMoleculeLayoutAlgorithm::CoordGen)
                .unwrap()
                .render_svg();
            let reaction = Reaction::from_rust(py, reaction).unwrap();

            let svg = reaction
                .depict_with(py, MoleculeLayoutAlgorithm::CoordGen())
                .unwrap();

            assert_eq!(svg._repr_svg_(), expected);
        });
    }

    #[rstest]
    fn test_reaction_depict_with_error() {
        Python::attach(|py| {
            let reaction = GraphIrReaction::new(
                umol_graph_ir::mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
                Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order {
                        old: NumForm::Lit(2),
                        new: NumForm::Lit(3),
                    },
                })]),
            );
            let reaction = Reaction::from_rust(py, reaction).unwrap();

            let error = reaction
                .depict_with(py, MoleculeLayoutAlgorithm::CoordGen())
                .unwrap_err();

            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reached a contradiction"
            );
        });
    }
}
