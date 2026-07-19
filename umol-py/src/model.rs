//! Python bindings for chemistry-model values.

pub(crate) mod aromaticity;
pub(crate) mod stereo;
pub(crate) mod valence;

use pyo3::prelude::*;
use umol_chem::element::Element as ChemElement;
use umol_graph::ops::model::{
    ChemistryModel as GraphChemistryModel, ElementScope as GraphElementScope,
};

use self::aromaticity::AromaticityModel;
use self::stereo::StereoModel;
use self::valence::ValenceModel;
use crate::element::Element;

/// Semantic configuration for valence, aromaticity, and stereo perception.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq)]
pub struct ChemistryModel {
    valence: ValenceModel,
    aromaticity: AromaticityModel,
    stereo: StereoModel,
}

#[pymethods]
impl ChemistryModel {
    #[new]
    #[pyo3(signature = (*, valence, aromaticity, stereo))]
    fn new(valence: ValenceModel, aromaticity: AromaticityModel, stereo: StereoModel) -> Self {
        Self {
            valence,
            aromaticity,
            stereo,
        }
    }

    #[staticmethod]
    fn default() -> Self {
        Self::from_rust(&GraphChemistryModel::default())
    }

    #[getter]
    fn valence(&self) -> ValenceModel {
        self.valence.clone()
    }

    #[getter]
    fn aromaticity(&self) -> AromaticityModel {
        self.aromaticity.clone()
    }

    #[getter]
    fn stereo(&self) -> StereoModel {
        self.stereo.clone()
    }

    fn __repr__(&self) -> String {
        if self == &Self::from_rust(&GraphChemistryModel::default()) {
            return "ChemistryModel.default()".to_owned();
        }
        format!(
            "ChemistryModel(valence={}, aromaticity={}, stereo={})",
            self.valence.__repr__(),
            self.aromaticity.__repr__(),
            self.stereo.__repr__(),
        )
    }
}

impl ChemistryModel {
    pub(crate) fn from_rust(model: &GraphChemistryModel) -> Self {
        Self {
            valence: ValenceModel::from_rust(&model.valence),
            aromaticity: AromaticityModel::from_rust(&model.aromaticity),
            stereo: StereoModel::from_rust(&model.stereo),
        }
    }

    pub(crate) fn to_rust(&self) -> GraphChemistryModel {
        GraphChemistryModel {
            valence: self.valence.to_rust(),
            aromaticity: self.aromaticity.to_rust(),
            stereo: self.stereo.to_rust(),
        }
    }
}

/// Elements eligible for a chemistry-model operation.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ElementScope {
    /// Every element is eligible.
    Any {},
    /// Only the listed elements are eligible.
    #[pyo3(constructor = (elements))]
    AllowList { elements: Vec<Element> },
}

#[pymethods]
impl ElementScope {
    fn __repr__(&self) -> String {
        match self {
            Self::Any {} => "ElementScope.Any()".to_owned(),
            Self::AllowList { elements } => {
                let elements = elements
                    .iter()
                    .map(|element| format!("Element('{}')", ChemElement::from(element).symbol()))
                    .collect::<Vec<_>>()
                    .join(", ");
                format!("ElementScope.AllowList([{elements}])")
            }
        }
    }
}

impl ElementScope {
    pub(crate) fn from_rust(scope: &GraphElementScope) -> Self {
        match scope {
            GraphElementScope::Any => Self::Any {},
            GraphElementScope::AllowList(elements) => Self::AllowList {
                elements: elements.iter().copied().map(Element::from).collect(),
            },
        }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for aggregate model configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphElementScope {
        match self {
            Self::Any {} => GraphElementScope::Any,
            Self::AllowList { elements } => {
                GraphElementScope::AllowList(elements.iter().map(ChemElement::from).collect())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::rstest;
    use umol_graph::ops::model::{
        AromaticityModel as GraphAromaticityModel, InconsistencyPolicy as GraphInconsistencyPolicy,
        RingLimits as GraphRingLimits, StereoModel as GraphStereoModel,
        ValenceModel as GraphValenceModel,
    };
    use umol_graph::valence_table;

    use super::*;

    #[rstest]
    fn test_chemistry_model_default() {
        let model = GraphChemistryModel::default();

        assert_eq!(
            ChemistryModel::default(),
            ChemistryModel {
                valence: ValenceModel::from_rust(&model.valence),
                aromaticity: AromaticityModel::from_rust(&model.aromaticity),
                stereo: StereoModel::from_rust(&model.stereo),
            }
        );
    }

    #[rstest]
    #[case::default(ChemistryModel::default(), "ChemistryModel.default()")]
    #[case::configured(
        ChemistryModel::new(
            ValenceModel::from_rust(&GraphValenceModel::Counts {
                table: Cow::Owned(valence_table![C => [4]]),
            }),
            AromaticityModel::from_rust(&GraphAromaticityModel::Hmo {
                scope: GraphElementScope::Any,
                stabilization_threshold: 0.375,
            }),
            StereoModel::from_rust(&GraphStereoModel {
                para_stereo: true,
                max_iterations: 8,
                inconsistency: GraphInconsistencyPolicy::Keep,
                ..GraphStereoModel::default()
            }),
        ),
        "ChemistryModel(valence=ValenceModel.Counts(table=ValenceTable(entries={Element('C'): ValenceEntry(target_covalences=[4], aromatic_valences=[])})), aromaticity=AromaticityModel.Hmo(scope=ElementScope.Any(), stabilization_threshold=0.375), stereo=StereoModel(kind_models={StereoKind.Tetrahedral: StereoKindModel(scope=ElementScope.Any(), fluxionality=False), StereoKind.CisTrans: StereoKindModel(scope=ElementScope.Any(), fluxionality=False)}, para_stereo=True, max_iterations=8, inconsistency=InconsistencyPolicy.Keep))"
    )]
    fn test_chemistry_model_repr(#[case] model: ChemistryModel, #[case] expected: &str) {
        assert_eq!(model.__repr__(), expected);
    }

    #[rstest]
    fn test_chemistry_model_from_rust() {
        let model = GraphChemistryModel {
            valence: GraphValenceModel::Counts {
                table: Cow::Owned(valence_table![C => [4], O => [2]]),
            },
            aromaticity: GraphAromaticityModel::Clar {
                scope: GraphElementScope::AllowList(vec![ChemElement::C]),
                ring_limits: GraphRingLimits {
                    min_ring_size: 6,
                    ..GraphRingLimits::default()
                },
            },
            stereo: GraphStereoModel {
                para_stereo: true,
                max_iterations: 8,
                inconsistency: GraphInconsistencyPolicy::Strip,
                ..GraphStereoModel::default()
            },
        };

        assert_eq!(
            ChemistryModel::from_rust(&model),
            ChemistryModel {
                valence: ValenceModel::from_rust(&model.valence),
                aromaticity: AromaticityModel::from_rust(&model.aromaticity),
                stereo: StereoModel::from_rust(&model.stereo),
            }
        );
    }

    #[rstest]
    fn test_chemistry_model_to_rust() {
        let expected = GraphChemistryModel {
            valence: GraphValenceModel::Counts {
                table: Cow::Owned(valence_table![C => [4], O => [2]]),
            },
            aromaticity: GraphAromaticityModel::Clar {
                scope: GraphElementScope::AllowList(vec![ChemElement::C]),
                ring_limits: GraphRingLimits {
                    min_ring_size: 6,
                    ..GraphRingLimits::default()
                },
            },
            stereo: GraphStereoModel {
                para_stereo: true,
                max_iterations: 8,
                inconsistency: GraphInconsistencyPolicy::Strip,
                ..GraphStereoModel::default()
            },
        };
        let model = ChemistryModel {
            valence: ValenceModel::from_rust(&expected.valence),
            aromaticity: AromaticityModel::from_rust(&expected.aromaticity),
            stereo: StereoModel::from_rust(&expected.stereo),
        };

        assert_eq!(model.to_rust(), expected);
    }

    #[rstest]
    #[case::any(ElementScope::Any {}, "ElementScope.Any()")]
    #[case::empty(
        ElementScope::AllowList {
            elements: Vec::new(),
        },
        "ElementScope.AllowList([])"
    )]
    #[case::populated(
        ElementScope::AllowList {
            elements: vec![Element::from(ChemElement::C), Element::from(ChemElement::N)],
        },
        "ElementScope.AllowList([Element('C'), Element('N')])"
    )]
    fn test_element_scope_repr(#[case] scope: ElementScope, #[case] expected: &str) {
        assert_eq!(scope.__repr__(), expected);
    }

    #[rstest]
    #[case::any(GraphElementScope::Any, ElementScope::Any {})]
    #[case::empty(
        GraphElementScope::AllowList(Vec::new()),
        ElementScope::AllowList {
            elements: Vec::new(),
        }
    )]
    #[case::populated(
        GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N, ChemElement::C]),
        ElementScope::AllowList {
            elements: vec![
                Element::from(ChemElement::C),
                Element::from(ChemElement::N),
                Element::from(ChemElement::C),
            ],
        }
    )]
    fn test_element_scope_from_rust(
        #[case] scope: GraphElementScope,
        #[case] expected: ElementScope,
    ) {
        assert_eq!(ElementScope::from_rust(&scope), expected);
    }

    #[rstest]
    #[case::any(ElementScope::Any {}, GraphElementScope::Any)]
    #[case::empty(
        ElementScope::AllowList {
            elements: Vec::new(),
        },
        GraphElementScope::AllowList(Vec::new())
    )]
    #[case::populated(
        ElementScope::AllowList {
            elements: vec![
                Element::from(ChemElement::C),
                Element::from(ChemElement::N),
                Element::from(ChemElement::C),
            ],
        },
        GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N, ChemElement::C])
    )]
    fn test_element_scope_to_rust(
        #[case] scope: ElementScope,
        #[case] expected: GraphElementScope,
    ) {
        assert_eq!(scope.to_rust(), expected);
    }
}
