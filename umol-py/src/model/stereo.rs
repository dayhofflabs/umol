//! Python bindings for stereo-model values.

use std::array;
use std::collections::BTreeMap;

use pyo3::prelude::*;
use umol_graph::ops::model::{
    StereoKindModel as GraphStereoKindModel, StereoModel as GraphStereoModel,
};

use super::ElementScope;
use crate::stereo::StereoKind;

/// Per-kind element eligibility and fluxionality settings.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoKindModel {
    scope: ElementScope,
    fluxionality: bool,
}

#[pymethods]
impl StereoKindModel {
    #[new]
    #[pyo3(signature = (*, scope, fluxionality))]
    fn new(scope: ElementScope, fluxionality: bool) -> Self {
        Self {
            scope,
            fluxionality,
        }
    }

    #[getter]
    fn scope(&self) -> ElementScope {
        self.scope.clone()
    }

    #[getter]
    fn fluxionality(&self) -> bool {
        self.fluxionality
    }

    fn __repr__(&self) -> String {
        format!(
            "StereoKindModel(scope={}, fluxionality={})",
            self.scope.__repr__(),
            if self.fluxionality { "True" } else { "False" },
        )
    }
}

impl StereoKindModel {
    pub(crate) fn from_rust(model: &GraphStereoKindModel) -> Self {
        Self {
            scope: ElementScope::from_rust(&model.scope),
            fluxionality: model.fluxionality,
        }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for StereoModel configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphStereoKindModel {
        GraphStereoKindModel {
            scope: self.scope.to_rust(),
            fluxionality: self.fluxionality,
        }
    }
}

const STEREO_KINDS: &[StereoKind] = &[
    StereoKind::Tetrahedral,
    StereoKind::CisTrans,
    StereoKind::Axial,
    StereoKind::SquarePlanar,
    StereoKind::TrigonalBipyramidal,
    StereoKind::Octahedral,
];

/// Stereo perception model with per-kind settings.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoModel {
    kind_models: BTreeMap<StereoKind, StereoKindModel>,
    para_stereo: bool,
}

#[pymethods]
impl StereoModel {
    #[new]
    #[pyo3(signature = (*, kind_models, para_stereo))]
    fn new(kind_models: BTreeMap<StereoKind, StereoKindModel>, para_stereo: bool) -> Self {
        Self {
            kind_models,
            para_stereo,
        }
    }

    #[staticmethod]
    fn default() -> Self {
        Self::from_rust(&GraphStereoModel::default())
    }

    #[getter]
    fn kind_models(&self) -> BTreeMap<StereoKind, StereoKindModel> {
        self.kind_models.clone()
    }

    #[getter]
    fn para_stereo(&self) -> bool {
        self.para_stereo
    }

    pub(crate) fn __repr__(&self) -> String {
        if self == &Self::from_rust(&GraphStereoModel::default()) {
            return "StereoModel.default()".to_owned();
        }

        let kind_models = self
            .kind_models
            .iter()
            .map(|(kind, model)| format!("StereoKind.{kind:?}: {}", model.__repr__()))
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "StereoModel(kind_models={{{kind_models}}}, para_stereo={})",
            if self.para_stereo { "True" } else { "False" },
        )
    }
}

impl StereoModel {
    pub(crate) fn from_rust(model: &GraphStereoModel) -> Self {
        let kind_models = STEREO_KINDS
            .iter()
            .copied()
            .filter_map(|kind| {
                model
                    .kind_model(kind.to_rust())
                    .map(|kind_model| (kind, StereoKindModel::from_rust(kind_model)))
            })
            .collect();
        Self {
            kind_models,
            para_stereo: model.para_stereo,
        }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for ChemistryModel configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphStereoModel {
        let mut kind_models = array::from_fn(|_| None);
        for (kind, model) in &self.kind_models {
            kind_models[kind.to_rust() as usize] = Some(model.to_rust());
        }
        GraphStereoModel {
            kind_models,
            para_stereo: self.para_stereo,
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element as ChemElement;
    use umol_graph::ops::model::ElementScope as GraphElementScope;

    use super::*;

    #[rstest]
    #[case::any(
        ElementScope::Any {},
        false,
        StereoKindModel {
            scope: ElementScope::Any {},
            fluxionality: false,
        }
    )]
    #[case::allow_list(
        ElementScope::AllowList {
            elements: vec![ChemElement::C.into(), ChemElement::N.into()],
        },
        true,
        StereoKindModel {
            scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into(), ChemElement::N.into()],
            },
            fluxionality: true,
        }
    )]
    fn test_stereo_kind_model_new(
        #[case] scope: ElementScope,
        #[case] fluxionality: bool,
        #[case] expected: StereoKindModel,
    ) {
        assert_eq!(StereoKindModel::new(scope, fluxionality), expected);
    }

    #[rstest]
    #[case::any(
        StereoKindModel::new(ElementScope::Any {}, false),
        "StereoKindModel(scope=ElementScope.Any(), fluxionality=False)"
    )]
    #[case::allow_list(
        StereoKindModel::new(
            ElementScope::AllowList {
                elements: vec![ChemElement::C.into(), ChemElement::N.into()],
            },
            true,
        ),
        "StereoKindModel(scope=ElementScope.AllowList([Element('C'), Element('N')]), fluxionality=True)"
    )]
    fn test_stereo_kind_model_repr(#[case] model: StereoKindModel, #[case] expected: &str) {
        assert_eq!(model.__repr__(), expected);
    }

    #[rstest]
    #[case::any(
        GraphStereoKindModel {
            scope: GraphElementScope::Any,
            fluxionality: false,
        },
        StereoKindModel::new(ElementScope::Any {}, false)
    )]
    #[case::allow_list(
        GraphStereoKindModel {
            scope: GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N]),
            fluxionality: true,
        },
        StereoKindModel::new(
            ElementScope::AllowList {
                elements: vec![ChemElement::C.into(), ChemElement::N.into()],
            },
            true,
        )
    )]
    fn test_stereo_kind_model_from_rust(
        #[case] model: GraphStereoKindModel,
        #[case] expected: StereoKindModel,
    ) {
        assert_eq!(StereoKindModel::from_rust(&model), expected);
    }

    #[rstest]
    #[case::any(
        StereoKindModel::new(ElementScope::Any {}, false),
        GraphStereoKindModel {
            scope: GraphElementScope::Any,
            fluxionality: false,
        }
    )]
    #[case::allow_list(
        StereoKindModel::new(
            ElementScope::AllowList {
                elements: vec![ChemElement::C.into(), ChemElement::N.into()],
            },
            true,
        ),
        GraphStereoKindModel {
            scope: GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N]),
            fluxionality: true,
        }
    )]
    fn test_stereo_kind_model_to_rust(
        #[case] model: StereoKindModel,
        #[case] expected: GraphStereoKindModel,
    ) {
        assert_eq!(model.to_rust(), expected);
    }

    #[rstest]
    fn test_stereo_model_default() {
        assert_eq!(
            StereoModel::default(),
            StereoModel {
                kind_models: BTreeMap::from([
                    (
                        StereoKind::Tetrahedral,
                        StereoKindModel::new(ElementScope::Any {}, false),
                    ),
                    (
                        StereoKind::CisTrans,
                        StereoKindModel::new(ElementScope::Any {}, false),
                    ),
                ]),
                para_stereo: false,
            }
        );
    }

    #[rstest]
    #[case::default(StereoModel::default(), "StereoModel.default()")]
    #[case::configured(
        StereoModel::new(
            BTreeMap::from([
                (
                    StereoKind::Tetrahedral,
                    StereoKindModel::new(ElementScope::Any {}, false),
                ),
                (
                    StereoKind::Octahedral,
                    StereoKindModel::new(
                        ElementScope::AllowList {
                            elements: vec![ChemElement::Fe.into()],
                        },
                        true,
                    ),
                ),
            ]),
            true,
        ),
        "StereoModel(kind_models={StereoKind.Tetrahedral: StereoKindModel(scope=ElementScope.Any(), fluxionality=False), StereoKind.Octahedral: StereoKindModel(scope=ElementScope.AllowList([Element('Fe')]), fluxionality=True)}, para_stereo=True)"
    )]
    fn test_stereo_model_repr(#[case] model: StereoModel, #[case] expected: &str) {
        assert_eq!(model.__repr__(), expected);
    }

    #[rstest]
    fn test_stereo_model_from_rust() {
        let model = GraphStereoModel {
            kind_models: [
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::Any,
                    fluxionality: false,
                }),
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::AllowList(vec![ChemElement::C]),
                    fluxionality: true,
                }),
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::AllowList(vec![ChemElement::N]),
                    fluxionality: false,
                }),
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::AllowList(vec![ChemElement::O]),
                    fluxionality: true,
                }),
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::AllowList(vec![ChemElement::S]),
                    fluxionality: false,
                }),
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::AllowList(vec![ChemElement::Fe]),
                    fluxionality: true,
                }),
            ],
            para_stereo: true,
        };

        assert_eq!(
            StereoModel::from_rust(&model),
            StereoModel {
                kind_models: BTreeMap::from([
                    (
                        StereoKind::Tetrahedral,
                        StereoKindModel::new(ElementScope::Any {}, false),
                    ),
                    (
                        StereoKind::CisTrans,
                        StereoKindModel::new(
                            ElementScope::AllowList {
                                elements: vec![ChemElement::C.into()],
                            },
                            true,
                        ),
                    ),
                    (
                        StereoKind::Axial,
                        StereoKindModel::new(
                            ElementScope::AllowList {
                                elements: vec![ChemElement::N.into()],
                            },
                            false,
                        ),
                    ),
                    (
                        StereoKind::SquarePlanar,
                        StereoKindModel::new(
                            ElementScope::AllowList {
                                elements: vec![ChemElement::O.into()],
                            },
                            true,
                        ),
                    ),
                    (
                        StereoKind::TrigonalBipyramidal,
                        StereoKindModel::new(
                            ElementScope::AllowList {
                                elements: vec![ChemElement::S.into()],
                            },
                            false,
                        ),
                    ),
                    (
                        StereoKind::Octahedral,
                        StereoKindModel::new(
                            ElementScope::AllowList {
                                elements: vec![ChemElement::Fe.into()],
                            },
                            true,
                        ),
                    ),
                ]),
                para_stereo: true,
            }
        );
    }

    #[rstest]
    fn test_stereo_model_to_rust() {
        let model = StereoModel::new(
            BTreeMap::from([
                (
                    StereoKind::Tetrahedral,
                    StereoKindModel::new(ElementScope::Any {}, false),
                ),
                (
                    StereoKind::CisTrans,
                    StereoKindModel::new(
                        ElementScope::AllowList {
                            elements: vec![ChemElement::C.into()],
                        },
                        true,
                    ),
                ),
                (
                    StereoKind::Axial,
                    StereoKindModel::new(
                        ElementScope::AllowList {
                            elements: vec![ChemElement::N.into()],
                        },
                        false,
                    ),
                ),
                (
                    StereoKind::SquarePlanar,
                    StereoKindModel::new(
                        ElementScope::AllowList {
                            elements: vec![ChemElement::O.into()],
                        },
                        true,
                    ),
                ),
                (
                    StereoKind::TrigonalBipyramidal,
                    StereoKindModel::new(
                        ElementScope::AllowList {
                            elements: vec![ChemElement::S.into()],
                        },
                        false,
                    ),
                ),
                (
                    StereoKind::Octahedral,
                    StereoKindModel::new(
                        ElementScope::AllowList {
                            elements: vec![ChemElement::Fe.into()],
                        },
                        true,
                    ),
                ),
            ]),
            true,
        );

        assert_eq!(
            model.to_rust(),
            GraphStereoModel {
                kind_models: [
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::Any,
                        fluxionality: false,
                    }),
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::AllowList(vec![ChemElement::C]),
                        fluxionality: true,
                    }),
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::AllowList(vec![ChemElement::N]),
                        fluxionality: false,
                    }),
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::AllowList(vec![ChemElement::O]),
                        fluxionality: true,
                    }),
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::AllowList(vec![ChemElement::S]),
                        fluxionality: false,
                    }),
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::AllowList(vec![ChemElement::Fe]),
                        fluxionality: true,
                    }),
                ],
                para_stereo: true,
            }
        );
    }
}
