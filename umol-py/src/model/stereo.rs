//! Python bindings for stereo-model values.

use std::array;
use std::collections::BTreeMap;

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyMapping, PyMappingProxy};
use umol_graph::ops::model::{
    StereoKindModel as GraphStereoKindModel, StereoModel as GraphStereoModel,
};

use super::ElementScope;
use crate::stereo::StereoKind;

/// Per-kind element eligibility settings.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoKindModel {
    scope: ElementScope,
}

#[pymethods]
impl StereoKindModel {
    #[new]
    #[pyo3(signature = (*, scope))]
    fn new(scope: ElementScope) -> Self {
        Self { scope }
    }

    #[getter]
    fn scope(&self) -> ElementScope {
        self.scope.clone()
    }

    fn __repr__(&self) -> String {
        format!("StereoKindModel(scope={})", self.scope.__repr__())
    }
}

impl StereoKindModel {
    pub(crate) fn from_rust(model: &GraphStereoKindModel) -> Self {
        Self {
            scope: ElementScope::from_rust(&model.scope),
        }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for StereoModel configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphStereoKindModel {
        GraphStereoKindModel {
            scope: self.scope.to_rust(),
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
    stereo_bond_min_ring_size: u32,
    stereo_bond_e_min_ring_size: u32,
}

#[pymethods]
impl StereoModel {
    #[new]
    #[pyo3(signature = (*, kind_models, stereo_bond_min_ring_size, stereo_bond_e_min_ring_size))]
    fn new(
        kind_models: &Bound<'_, PyMapping>,
        stereo_bond_min_ring_size: u32,
        stereo_bond_e_min_ring_size: u32,
    ) -> PyResult<Self> {
        Ok(Self {
            kind_models: kind_models
                .try_iter()?
                .map(|key| {
                    let key = key?;
                    let value = kind_models.get_item(&key)?;
                    Ok((key.extract()?, value.extract()?))
                })
                .collect::<PyResult<_>>()?,
            stereo_bond_min_ring_size,
            stereo_bond_e_min_ring_size,
        })
    }

    #[staticmethod]
    fn default() -> Self {
        Self::from_rust(&GraphStereoModel::default())
    }

    #[getter]
    fn kind_models<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyMappingProxy>> {
        let values = PyDict::new(py);
        for (kind, model) in &self.kind_models {
            values.set_item(*kind, model.clone())?;
        }
        Ok(PyMappingProxy::new(py, values.as_mapping()))
    }

    #[getter]
    fn stereo_bond_min_ring_size(&self) -> u32 {
        self.stereo_bond_min_ring_size
    }

    #[getter]
    fn stereo_bond_e_min_ring_size(&self) -> u32 {
        self.stereo_bond_e_min_ring_size
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
            "StereoModel(kind_models={{{kind_models}}}, stereo_bond_min_ring_size={}, stereo_bond_e_min_ring_size={})",
            self.stereo_bond_min_ring_size,
            self.stereo_bond_e_min_ring_size,
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
            stereo_bond_min_ring_size: model.stereo_bond_min_ring_size,
            stereo_bond_e_min_ring_size: model.stereo_bond_e_min_ring_size,
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
            stereo_bond_min_ring_size: self.stereo_bond_min_ring_size,
            stereo_bond_e_min_ring_size: self.stereo_bond_e_min_ring_size,
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element as ChemElement;
    use umol_graph::ops::model::ElementScope as GraphElementScope;

    use super::*;
    use crate::convert::TupleItems;

    #[rstest]
    #[case::any(
        ElementScope::Any {},
        StereoKindModel {
            scope: ElementScope::Any {},
        }
    )]
    #[case::allow_list(
        ElementScope::AllowList {
            elements: TupleItems(vec![ChemElement::C.into(), ChemElement::N.into()]),
        },
        StereoKindModel {
            scope: ElementScope::AllowList {
                elements: TupleItems(vec![ChemElement::C.into(), ChemElement::N.into()]),
            },
        }
    )]
    fn test_stereo_kind_model_new(#[case] scope: ElementScope, #[case] expected: StereoKindModel) {
        assert_eq!(StereoKindModel::new(scope), expected);
    }

    #[rstest]
    #[case::any(
        StereoKindModel::new(ElementScope::Any {}),
        "StereoKindModel(scope=ElementScope.Any())"
    )]
    #[case::allow_list(
        StereoKindModel::new(
            ElementScope::AllowList {
                elements: TupleItems(vec![ChemElement::C.into(), ChemElement::N.into()]),
            }),
        "StereoKindModel(scope=ElementScope.AllowList([Element('C'), Element('N')]))"
    )]
    fn test_stereo_kind_model_repr(#[case] model: StereoKindModel, #[case] expected: &str) {
        assert_eq!(model.__repr__(), expected);
    }

    #[rstest]
    #[case::any(
        GraphStereoKindModel {
            scope: GraphElementScope::Any,
        },
        StereoKindModel::new(ElementScope::Any {})
    )]
    #[case::allow_list(
        GraphStereoKindModel {
            scope: GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N]),
        },
        StereoKindModel::new(
            ElementScope::AllowList {
                elements: TupleItems(vec![ChemElement::C.into(), ChemElement::N.into()]),
            })
    )]
    fn test_stereo_kind_model_from_rust(
        #[case] model: GraphStereoKindModel,
        #[case] expected: StereoKindModel,
    ) {
        assert_eq!(StereoKindModel::from_rust(&model), expected);
    }

    #[rstest]
    #[case::any(
        StereoKindModel::new(ElementScope::Any {}),
        GraphStereoKindModel {
            scope: GraphElementScope::Any,
        }
    )]
    #[case::allow_list(
        StereoKindModel::new(
            ElementScope::AllowList {
                elements: TupleItems(vec![ChemElement::C.into(), ChemElement::N.into()]),
            }),
        GraphStereoKindModel {
            scope: GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N]),
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
                        StereoKindModel::new(ElementScope::Any {}),
                    ),
                    (
                        StereoKind::CisTrans,
                        StereoKindModel::new(ElementScope::Any {}),
                    ),
                ]),
                stereo_bond_min_ring_size: 8,
                stereo_bond_e_min_ring_size: 8,
            }
        );
    }

    #[rstest]
    #[case::default(StereoModel::default(), "StereoModel.default()")]
    #[case::configured(
        StereoModel { kind_models: BTreeMap::from([
                (
                    StereoKind::Tetrahedral,
                    StereoKindModel::new(ElementScope::Any {}),
                ),
                (
                    StereoKind::Octahedral,
                    StereoKindModel::new(
                        ElementScope::AllowList {
                            elements: TupleItems(vec![ChemElement::Fe.into()]),
                        }),
                ),
            ]), stereo_bond_min_ring_size: 8, stereo_bond_e_min_ring_size: 8 },
        "StereoModel(kind_models={StereoKind.Tetrahedral: StereoKindModel(scope=ElementScope.Any()), StereoKind.Octahedral: StereoKindModel(scope=ElementScope.AllowList([Element('Fe')]))}, stereo_bond_min_ring_size=8, stereo_bond_e_min_ring_size=8)"
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
                }),
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::AllowList(vec![ChemElement::C]),
                }),
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::AllowList(vec![ChemElement::N]),
                }),
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::AllowList(vec![ChemElement::O]),
                }),
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::AllowList(vec![ChemElement::S]),
                }),
                Some(GraphStereoKindModel {
                    scope: GraphElementScope::AllowList(vec![ChemElement::Fe]),
                }),
            ],
            stereo_bond_min_ring_size: 10,
            stereo_bond_e_min_ring_size: 12,
        };

        assert_eq!(
            StereoModel::from_rust(&model),
            StereoModel {
                kind_models: BTreeMap::from([
                    (
                        StereoKind::Tetrahedral,
                        StereoKindModel::new(ElementScope::Any {}),
                    ),
                    (
                        StereoKind::CisTrans,
                        StereoKindModel::new(ElementScope::AllowList {
                            elements: TupleItems(vec![ChemElement::C.into()]),
                        }),
                    ),
                    (
                        StereoKind::Axial,
                        StereoKindModel::new(ElementScope::AllowList {
                            elements: TupleItems(vec![ChemElement::N.into()]),
                        }),
                    ),
                    (
                        StereoKind::SquarePlanar,
                        StereoKindModel::new(ElementScope::AllowList {
                            elements: TupleItems(vec![ChemElement::O.into()]),
                        }),
                    ),
                    (
                        StereoKind::TrigonalBipyramidal,
                        StereoKindModel::new(ElementScope::AllowList {
                            elements: TupleItems(vec![ChemElement::S.into()]),
                        }),
                    ),
                    (
                        StereoKind::Octahedral,
                        StereoKindModel::new(ElementScope::AllowList {
                            elements: TupleItems(vec![ChemElement::Fe.into()]),
                        }),
                    ),
                ]),
                stereo_bond_min_ring_size: 10,
                stereo_bond_e_min_ring_size: 12,
            }
        );
    }

    #[rstest]
    fn test_stereo_model_to_rust() {
        let model = StereoModel {
            kind_models: BTreeMap::from([
                (
                    StereoKind::Tetrahedral,
                    StereoKindModel::new(ElementScope::Any {}),
                ),
                (
                    StereoKind::CisTrans,
                    StereoKindModel::new(ElementScope::AllowList {
                        elements: TupleItems(vec![ChemElement::C.into()]),
                    }),
                ),
                (
                    StereoKind::Axial,
                    StereoKindModel::new(ElementScope::AllowList {
                        elements: TupleItems(vec![ChemElement::N.into()]),
                    }),
                ),
                (
                    StereoKind::SquarePlanar,
                    StereoKindModel::new(ElementScope::AllowList {
                        elements: TupleItems(vec![ChemElement::O.into()]),
                    }),
                ),
                (
                    StereoKind::TrigonalBipyramidal,
                    StereoKindModel::new(ElementScope::AllowList {
                        elements: TupleItems(vec![ChemElement::S.into()]),
                    }),
                ),
                (
                    StereoKind::Octahedral,
                    StereoKindModel::new(ElementScope::AllowList {
                        elements: TupleItems(vec![ChemElement::Fe.into()]),
                    }),
                ),
            ]),
            stereo_bond_min_ring_size: 10,
            stereo_bond_e_min_ring_size: 12,
        };

        assert_eq!(
            model.to_rust(),
            GraphStereoModel {
                kind_models: [
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::Any,
                    }),
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::AllowList(vec![ChemElement::C]),
                    }),
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::AllowList(vec![ChemElement::N]),
                    }),
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::AllowList(vec![ChemElement::O]),
                    }),
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::AllowList(vec![ChemElement::S]),
                    }),
                    Some(GraphStereoKindModel {
                        scope: GraphElementScope::AllowList(vec![ChemElement::Fe]),
                    }),
                ],
                stereo_bond_min_ring_size: 10,
                stereo_bond_e_min_ring_size: 12,
            }
        );
    }
}
