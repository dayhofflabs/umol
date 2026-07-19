//! Python bindings for chemistry-model values.

pub(crate) mod aromaticity;
pub(crate) mod valence;

use pyo3::prelude::*;
use umol_chem::element::Element as ChemElement;
use umol_graph::ops::model::ElementScope as GraphElementScope;

use crate::element::Element;

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
    use rstest::rstest;

    use super::*;

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
