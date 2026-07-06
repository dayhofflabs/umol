//! Atom-field value types mirroring `umol_ast::ast` (`ElementAst`; `IsotopeMassAst`
//! and `SpinStateAst` follow at S2b/c, `AtomAst` at S3).

use std::collections::BTreeSet;

use pyo3::prelude::*;
use umol_ast::ast::{
    ElementAst as AstElementAst, IsotopeMassAst as AstIsotopeMassAst,
    SpinStateAst as AstSpinStateAst,
};
use umol_chem::element::Element as ChemElement;

use crate::element::Element;
use crate::value::{MemOp, ValueAst};

/// Element expression: undetermined, a single element, a finite element set, a
/// complement set (`!{…}`), or a variable with an optional membership restriction.
#[pyclass]
pub enum ElementAst {
    Undetermined(),
    Lit(Element),
    LitSet(BTreeSet<Element>),
    NotSet(BTreeSet<Element>),
    Var(String, Option<(MemOp, BTreeSet<Element>)>),
}

// The AST bridge; consumed by `AtomAst` at S3 (unused in the lib until then).
#[allow(dead_code)]
impl ElementAst {
    pub(crate) fn from_ast(ast: &AstElementAst) -> ElementAst {
        match ast {
            AstElementAst::Undetermined => ElementAst::Undetermined(),
            AstElementAst::Lit(e) => ElementAst::Lit(Element::from(*e)),
            AstElementAst::LitSet(members) => {
                ElementAst::LitSet(members.iter().copied().map(Element::from).collect())
            }
            AstElementAst::NotSet(members) => {
                ElementAst::NotSet(members.iter().copied().map(Element::from).collect())
            }
            AstElementAst::Var(boxed) => {
                let (name, restriction) = &**boxed;
                ElementAst::Var(
                    name.clone(),
                    restriction.as_ref().map(|(op, members)| {
                        (
                            MemOp::from_ast(*op),
                            members.iter().copied().map(Element::from).collect(),
                        )
                    }),
                )
            }
        }
    }

    pub(crate) fn to_ast(&self) -> AstElementAst {
        match self {
            ElementAst::Undetermined() => AstElementAst::Undetermined,
            ElementAst::Lit(e) => AstElementAst::Lit(ChemElement::from(e)),
            ElementAst::LitSet(members) => {
                AstElementAst::LitSet(Box::new(members.iter().map(ChemElement::from).collect()))
            }
            ElementAst::NotSet(members) => {
                AstElementAst::NotSet(Box::new(members.iter().map(ChemElement::from).collect()))
            }
            ElementAst::Var(name, restriction) => AstElementAst::Var(Box::new((
                name.clone(),
                restriction.as_ref().map(|(op, members)| {
                    (op.to_ast(), members.iter().map(ChemElement::from).collect())
                }),
            ))),
        }
    }
}

/// Isotope-mass expression: undetermined, the natural isotopic mixture, a single
/// mass number, a finite mass set, or a variable with an optional mass-set restriction.
#[pyclass]
pub enum IsotopeMassAst {
    Undetermined(),
    Natural(),
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Var(String, Option<BTreeSet<u32>>),
}

// The AST bridge; consumed by `AtomAst` at S3 (unused in the lib until then).
#[allow(dead_code)]
impl IsotopeMassAst {
    pub(crate) fn from_ast(ast: &AstIsotopeMassAst) -> IsotopeMassAst {
        match ast {
            AstIsotopeMassAst::Undetermined => IsotopeMassAst::Undetermined(),
            AstIsotopeMassAst::Natural => IsotopeMassAst::Natural(),
            AstIsotopeMassAst::Lit(mass) => IsotopeMassAst::Lit(*mass),
            AstIsotopeMassAst::LitSet(masses) => IsotopeMassAst::LitSet((**masses).clone()),
            AstIsotopeMassAst::Var(boxed) => {
                let (name, restriction) = &**boxed;
                IsotopeMassAst::Var(name.clone(), restriction.clone())
            }
        }
    }

    pub(crate) fn to_ast(&self) -> AstIsotopeMassAst {
        match self {
            IsotopeMassAst::Undetermined() => AstIsotopeMassAst::Undetermined,
            IsotopeMassAst::Natural() => AstIsotopeMassAst::Natural,
            IsotopeMassAst::Lit(mass) => AstIsotopeMassAst::Lit(*mass),
            IsotopeMassAst::LitSet(masses) => AstIsotopeMassAst::LitSet(Box::new(masses.clone())),
            IsotopeMassAst::Var(name, restriction) => {
                AstIsotopeMassAst::Var(Box::new((name.clone(), restriction.clone())))
            }
        }
    }
}

/// Spin state: unpaired-electron count and multiplicity as independent value fields.
#[pyclass]
pub struct SpinStateAst {
    #[pyo3(get)]
    unpaired: Py<ValueAst>,
    #[pyo3(get)]
    multiplicity: Py<ValueAst>,
}

#[pymethods]
impl SpinStateAst {
    #[new]
    fn new(unpaired: Py<ValueAst>, multiplicity: Py<ValueAst>) -> Self {
        SpinStateAst {
            unpaired,
            multiplicity,
        }
    }
}

// The AST bridge; consumed by `AtomAst` at S3 (unused in the lib until then).
#[allow(dead_code)]
impl SpinStateAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstSpinStateAst) -> PyResult<SpinStateAst> {
        Ok(SpinStateAst {
            unpaired: Py::new(py, ValueAst::from_ast(py, &ast.unpaired)?)?,
            multiplicity: Py::new(py, ValueAst::from_ast(py, &ast.multiplicity)?)?,
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstSpinStateAst {
        AstSpinStateAst {
            unpaired: self.unpaired.bind(py).borrow().to_ast(py),
            multiplicity: self.multiplicity.bind(py).borrow().to_ast(py),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use umol_ast::ast::{MemOp as AstMemOp, ValueAst as AstValueAst};

    #[rstest]
    #[case(AstElementAst::Undetermined)]
    #[case(AstElementAst::Lit(ChemElement::C))]
    #[case(AstElementAst::LitSet(Box::new(BTreeSet::from([ChemElement::C, ChemElement::N]))))]
    #[case(AstElementAst::NotSet(Box::new(BTreeSet::from([ChemElement::O]))))]
    #[case(AstElementAst::Var(Box::new(("x".to_string(), None))))]
    #[case(AstElementAst::Var(Box::new((
        "y".to_string(),
        Some((AstMemOp::In, BTreeSet::from([ChemElement::C, ChemElement::N]))),
    ))))]
    fn test_element_ast_roundtrip(#[case] ast: AstElementAst) {
        assert_eq!(ElementAst::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstIsotopeMassAst::Undetermined)]
    #[case(AstIsotopeMassAst::Natural)]
    #[case(AstIsotopeMassAst::Lit(13))]
    #[case(AstIsotopeMassAst::LitSet(Box::new(BTreeSet::from([12, 13, 14]))))]
    #[case(AstIsotopeMassAst::Var(Box::new(("x".to_string(), None))))]
    #[case(AstIsotopeMassAst::Var(Box::new((
        "y".to_string(),
        Some(BTreeSet::from([12, 13])),
    ))))]
    fn test_isotope_mass_ast_roundtrip(#[case] ast: AstIsotopeMassAst) {
        assert_eq!(IsotopeMassAst::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstSpinStateAst { unpaired: AstValueAst::Lit(1), multiplicity: AstValueAst::Lit(2) })]
    #[case(AstSpinStateAst {
        unpaired: AstValueAst::Undetermined,
        multiplicity: AstValueAst::Undetermined,
    })]
    fn test_spin_state_ast_roundtrip(#[case] ast: AstSpinStateAst) {
        Python::attach(|py| {
            assert_eq!(SpinStateAst::from_ast(py, &ast).unwrap().to_ast(py), ast);
        });
    }
}
