//! Python bindings for chemistry-model values.

pub(crate) mod aromaticity;
pub(crate) mod connectivity;
pub(crate) mod stereo;
pub(crate) mod valence;

use pyo3::prelude::*;
use umol_chem::element::Element as ChemElement;
use umol_graph::ops::model::{
    ChemistryModel as GraphChemistryModel, ElementScope as GraphElementScope,
};

use self::aromaticity::AromaticityModel;
use self::connectivity::ConnectivityModel;
use self::stereo::StereoModel;
use self::valence::ValenceModel;
use crate::convert::TupleItems;
use crate::element::Element;

/// Semantic configuration for connectivity, valence, aromaticity, and stereo perception.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq)]
pub struct ChemistryModel(GraphChemistryModel);

#[pymethods]
impl ChemistryModel {
    #[new]
    #[pyo3(signature = (*, connectivity, valence, aromaticity, stereo))]
    fn new(
        connectivity: PyRef<'_, ConnectivityModel>,
        valence: PyRef<'_, ValenceModel>,
        aromaticity: PyRef<'_, AromaticityModel>,
        stereo: PyRef<'_, StereoModel>,
    ) -> Self {
        Self::from_rust(GraphChemistryModel {
            connectivity: connectivity.to_rust(),
            valence: valence.to_rust(),
            aromaticity: aromaticity.to_rust(),
            stereo: stereo.to_rust(),
        })
    }

    #[staticmethod]
    fn default() -> Self {
        Self::from_rust(GraphChemistryModel::default())
    }

    #[getter]
    fn connectivity(&self) -> ConnectivityModel {
        ConnectivityModel::from_rust(&self.0.connectivity)
    }

    #[getter]
    fn valence(&self) -> ValenceModel {
        ValenceModel::from_rust(&self.0.valence)
    }

    #[getter]
    fn aromaticity(&self) -> AromaticityModel {
        AromaticityModel::from_rust(&self.0.aromaticity)
    }

    #[getter]
    fn stereo(&self) -> StereoModel {
        StereoModel::from_rust(&self.0.stereo)
    }

    fn __repr__(&self) -> String {
        if self.0 == GraphChemistryModel::default() {
            return "ChemistryModel.default()".to_owned();
        }
        format!(
            "ChemistryModel(connectivity={}, valence={}, aromaticity={}, stereo={})",
            self.connectivity().__repr__(),
            self.valence().__repr__(),
            self.aromaticity().__repr__(),
            self.stereo().__repr__(),
        )
    }
}

impl ChemistryModel {
    pub(crate) fn from_rust(model: GraphChemistryModel) -> Self {
        Self(model)
    }

    pub(crate) fn to_rust(&self) -> &GraphChemistryModel {
        &self.0
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
    AllowList { elements: TupleItems<Vec<Element>> },
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
        AromaticityModel as GraphAromaticityModel, AromaticityRule as GraphAromaticityRule,
        AromaticityTieBreak as GraphAromaticityTieBreak, StereoModel as GraphStereoModel,
        ValenceCandidateSource as GraphValenceCandidateSource, ValenceModel as GraphValenceModel,
    };
    use umol_graph::ops::validate::ConnectivityModel as GraphConnectivityModel;
    use umol_graph::valence_table;

    use super::*;

    #[rstest]
    fn test_chemistry_model_default() {
        let model = ChemistryModel::default();

        assert_eq!(model, ChemistryModel(GraphChemistryModel::default()));
        assert!(matches!(
            model.to_rust().valence.candidates,
            GraphValenceCandidateSource::AtomTyping {
                registry: Cow::Borrowed(_)
            }
        ));
    }

    #[rstest]
    #[case::default(ChemistryModel::default(), "ChemistryModel.default()")]
    #[case::configured(
        ChemistryModel(GraphChemistryModel {
            connectivity: GraphConnectivityModel::default(),
            valence: GraphValenceModel::counts(Cow::Owned(valence_table![C => [4]])),
            aromaticity: GraphAromaticityModel { scope: GraphElementScope::Any, rule: GraphAromaticityRule::Hmo { stabilization_threshold: 0.375 }, tie_break: GraphAromaticityTieBreak::Strict },
            stereo: GraphStereoModel {
                para_stereo: true,
                ..GraphStereoModel::default()
            },
        }),
        "ChemistryModel(connectivity=ConnectivityModel.default(), valence=ValenceModel(candidates=ValenceCandidateSource.Counts(table=ValenceTable(entries={Element('C'): ValenceEntry(target_covalences=[4], aromatic_valences=[], fallback_aromatic_valences=[])})), tie_break=ValenceTieBreak.Strict), aromaticity=AromaticityModel(scope=ElementScope.Any(), rule=AromaticityRule.Hmo(stabilization_threshold=0.375), tie_break=AromaticityTieBreak.Strict), stereo=StereoModel(kind_models={StereoKind.Tetrahedral: StereoKindModel(scope=ElementScope.Any(), fluxionality=False), StereoKind.CisTrans: StereoKindModel(scope=ElementScope.Any(), fluxionality=False)}, para_stereo=True, stereo_bond_minimum_ring_size=8))"
    )]
    fn test_chemistry_model_repr(#[case] model: ChemistryModel, #[case] expected: &str) {
        assert_eq!(model.__repr__(), expected);
    }

    #[rstest]
    fn test_chemistry_model_from_rust() {
        let model = GraphChemistryModel {
            connectivity: GraphConnectivityModel {
                allow_disconnected: false,
                ..GraphConnectivityModel::default()
            },
            valence: GraphValenceModel::counts(Cow::Owned(valence_table![C => [4], O => [2]])),
            aromaticity: GraphAromaticityModel {
                scope: GraphElementScope::AllowList(vec![ChemElement::C]),
                rule: GraphAromaticityRule::Clar,
                tie_break: GraphAromaticityTieBreak::Strict,
            },
            stereo: GraphStereoModel {
                para_stereo: true,
                ..GraphStereoModel::default()
            },
        };

        assert_eq!(
            ChemistryModel::from_rust(model.clone()),
            ChemistryModel(model)
        );
    }

    #[rstest]
    fn test_chemistry_model_to_rust() {
        let expected = GraphChemistryModel {
            connectivity: GraphConnectivityModel {
                allow_disconnected: false,
                ..GraphConnectivityModel::default()
            },
            valence: GraphValenceModel::counts(Cow::Owned(valence_table![C => [4], O => [2]])),
            aromaticity: GraphAromaticityModel {
                scope: GraphElementScope::AllowList(vec![ChemElement::C]),
                rule: GraphAromaticityRule::Clar,
                tie_break: GraphAromaticityTieBreak::Strict,
            },
            stereo: GraphStereoModel {
                para_stereo: true,
                ..GraphStereoModel::default()
            },
        };
        let model = ChemistryModel(expected.clone());

        assert_eq!(model.to_rust(), &expected);
    }

    #[rstest]
    #[case::any(ElementScope::Any {}, "ElementScope.Any()")]
    #[case::empty(
        ElementScope::AllowList {
            elements: TupleItems(Vec::new()),
        },
        "ElementScope.AllowList([])"
    )]
    #[case::populated(
        ElementScope::AllowList {
            elements: TupleItems(vec![Element::from(ChemElement::C), Element::from(ChemElement::N)]),
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
            elements: TupleItems(Vec::new()),
        }
    )]
    #[case::populated(
        GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N, ChemElement::C]),
        ElementScope::AllowList {
            elements: TupleItems(vec![
                Element::from(ChemElement::C),
                Element::from(ChemElement::N),
                Element::from(ChemElement::C),
            ]),
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
            elements: TupleItems(Vec::new()),
        },
        GraphElementScope::AllowList(Vec::new())
    )]
    #[case::populated(
        ElementScope::AllowList {
            elements: TupleItems(vec![
                Element::from(ChemElement::C),
                Element::from(ChemElement::N),
                Element::from(ChemElement::C),
            ]),
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
