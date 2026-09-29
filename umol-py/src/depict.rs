//! Explicit layout and SVG depiction bindings, enabled by the `depiction` feature.

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use umol_io::depict::{
    Depict as IoDepict, DepictConfig as IoDepictConfig, Depiction as IoDepiction,
    MoleculeDepictionError as IoMoleculeDepictionError,
    ReactionDepictionError as IoReactionDepictionError, SvgConfig as IoSvgConfig,
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

/// Rendering configuration for `Depiction.render_svg_with`.
///
/// `mapping_index_text_size` is the font size of correspondence-pair numbers in bond-length units;
/// `mask_id` identifies the atom-label mask definition, which must differ between depictions
/// inlined in one document. The default renders exactly as `Depiction.render_svg`.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SvgConfig(IoSvgConfig);

#[pymethods]
impl SvgConfig {
    /// Raises `ValueError` if `mapping_index_text_size` is not a finite positive number or
    /// `mask_id` is empty.
    #[new]
    #[pyo3(signature = (*, mapping_index_text_size=None, mask_id=None))]
    fn new(mapping_index_text_size: Option<f64>, mask_id: Option<String>) -> PyResult<Self> {
        let defaults = IoSvgConfig::default();
        let mapping_index_text_size =
            mapping_index_text_size.unwrap_or(defaults.mapping_index_text_size);
        if !mapping_index_text_size.is_finite() || mapping_index_text_size <= 0.0 {
            return Err(PyValueError::new_err(format!(
                "mapping_index_text_size must be a finite positive number, got {mapping_index_text_size}"
            )));
        }
        let mask_id = mask_id.unwrap_or(defaults.mask_id);
        if mask_id.is_empty() {
            return Err(PyValueError::new_err("mask_id must not be empty"));
        }
        Ok(Self(IoSvgConfig {
            mapping_index_text_size,
            mask_id,
        }))
    }

    #[staticmethod]
    fn default() -> Self {
        Default::default()
    }

    #[getter]
    fn mapping_index_text_size(&self) -> f64 {
        self.0.mapping_index_text_size
    }

    #[getter]
    fn mask_id(&self) -> &str {
        &self.0.mask_id
    }

    fn __repr__(&self) -> String {
        if self == &Self::default() {
            return "SvgConfig.default()".to_owned();
        }
        format!(
            "SvgConfig(mapping_index_text_size={:?}, mask_id={:?})",
            self.0.mapping_index_text_size, self.0.mask_id
        )
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

    /// Render this depiction as a complete SVG document using `config`.
    fn render_svg_with(&self, config: &SvgConfig) -> String {
        self.0.render_svg_with(&config.0)
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
    /// Construct a format-neutral depiction using the default configuration.
    ///
    /// Raises `RuntimeError` if layout or tetrahedral depiction fails.
    fn depict(&self) -> PyResult<Depiction> {
        self.to_rust()
            .depict()
            .map(Depiction::from_rust)
            .map_err(molecule_depiction_error)
    }

    /// Construct a format-neutral depiction using `config`.
    ///
    /// Raises `RuntimeError` if layout or tetrahedral depiction fails.
    fn depict_with(&self, config: DepictConfig) -> PyResult<Depiction> {
        let config = config.to_rust();
        self.to_rust()
            .depict_with(&config)
            .map(Depiction::from_rust)
            .map_err(molecule_depiction_error)
    }
}

#[pymethods]
impl Reaction {
    /// Construct a format-neutral depiction using the default configuration.
    ///
    /// Raises `ContradictionError` if the reaction cannot be materialized and `RuntimeError` if
    /// layout or depiction of either materialized side fails.
    fn depict(&self, py: Python<'_>) -> PyResult<Depiction> {
        self.to_rust(py)?
            .depict()
            .map(Depiction::from_rust)
            .map_err(reaction_depiction_error)
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
    use pyo3::exceptions::{PyRuntimeError, PyTypeError, PyValueError};
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
        assert_eq!(depiction.render_svg_with(&SvgConfig::default()), expected);
    }

    #[rstest]
    fn test_svg_config_new() {
        let config = SvgConfig::new(None, None).unwrap();

        assert_eq!(config, SvgConfig::default());
        assert_eq!(config.mapping_index_text_size(), 0.3825);
        assert_eq!(config.mask_id(), "umol-atom-label-mask");
        assert_eq!(config.__repr__(), "SvgConfig.default()");
    }

    #[rstest]
    fn test_svg_config_new_custom() {
        let config = SvgConfig::new(Some(0.25), Some("lhs-mask".to_owned())).unwrap();

        assert_ne!(config, SvgConfig::default());
        assert_eq!(config.mapping_index_text_size(), 0.25);
        assert_eq!(config.mask_id(), "lhs-mask");
        assert_eq!(
            config.__repr__(),
            r#"SvgConfig(mapping_index_text_size=0.25, mask_id="lhs-mask")"#
        );
    }

    #[rstest]
    #[case::zero(
        Some(0.0),
        None,
        "mapping_index_text_size must be a finite positive number, got 0"
    )]
    #[case::negative(Some(-0.1), None, "mapping_index_text_size must be a finite positive number, got -0.1")]
    #[case::nan(
        Some(f64::NAN),
        None,
        "mapping_index_text_size must be a finite positive number, got NaN"
    )]
    #[case::infinite(
        Some(f64::INFINITY),
        None,
        "mapping_index_text_size must be a finite positive number, got inf"
    )]
    #[case::empty_mask_id(None, Some(String::new()), "mask_id must not be empty")]
    fn test_svg_config_new_error(
        #[case] mapping_index_text_size: Option<f64>,
        #[case] mask_id: Option<String>,
        #[case] message: &str,
    ) {
        Python::attach(|py| {
            let error = SvgConfig::new(mapping_index_text_size, mask_id).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                message
            );
        });
    }

    #[rstest]
    fn test_depiction_render_svg_with() {
        let reaction = GraphIrReaction::new(
            umol_graph_ir::mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            Deltas::new(),
        );
        let config = SvgConfig::new(Some(0.25), Some("lhs-mask".to_owned())).unwrap();
        let rust = reaction.depict().unwrap();
        let expected = rust.render_svg_with(&config.0);
        let depiction = Depiction::from_rust(rust);

        let svg = depiction.render_svg_with(&config);

        assert_eq!(svg, expected);
        assert!(svg.contains(r#"<mask id="lhs-mask""#));
        assert!(svg.contains(r#"font-size="0.25""#));
        assert_ne!(svg, depiction.render_svg());
    }

    #[rstest]
    fn test_depiction_constructor_error() {
        Python::attach(|py| {
            let error = py.get_type::<Depiction>().call0().unwrap_err();

            assert!(error.is_instance_of::<PyTypeError>(py));
        });
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
    fn test_reaction_depict() {
        Python::attach(|py| {
            let reaction = GraphIrReaction::new(
                umol_graph_ir::mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
                Deltas::new(),
            );
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
            let reaction = GraphIrReaction::new(
                umol_graph_ir::mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
                Deltas::new(),
            );
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
                .depict_with(py, DepictConfig::default())
                .err()
                .unwrap();

            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reached a contradiction"
            );
        });
    }
}
