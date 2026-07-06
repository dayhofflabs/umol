//! Atom-field value types mirroring `umol_ast::ast` (`ElementAst`; `IsotopeMassAst`
//! and `SpinStateAst` follow at S2b/c, `AtomAst` at S3).

use std::collections::BTreeSet;

use pyo3::prelude::*;
use umol_ast::ast::{
    AtomAst as AstAtomAst, ElementAst as AstElementAst, IsotopeMassAst as AstIsotopeMassAst,
    SpinStateAst as AstSpinStateAst,
};
use umol_chem::element::Element as ChemElement;

use crate::convert::into_py_variant;
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

impl SpinStateAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstSpinStateAst) -> PyResult<SpinStateAst> {
        Ok(SpinStateAst {
            unpaired: into_py_variant(py, ValueAst::from_ast(py, &ast.unpaired)?)?,
            multiplicity: into_py_variant(py, ValueAst::from_ast(py, &ast.multiplicity)?)?,
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstSpinStateAst {
        AstSpinStateAst {
            unpaired: self.unpaired.bind(py).borrow().to_ast(py),
            multiplicity: self.multiplicity.bind(py).borrow().to_ast(py),
        }
    }
}

/// An atom: element, isotope, charge, implicit hydrogens, lone pairs, and spin.
/// (Constraints are not yet exposed — S5.)
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct AtomAst(AstAtomAst);

#[pymethods]
impl AtomAst {
    /// Construct from an element expression.
    #[new]
    fn new(element: PyRef<'_, ElementAst>) -> Self {
        AtomAst(AstAtomAst::new(element.to_ast()))
    }

    /// Construct from a single element.
    #[staticmethod]
    fn from_element(element: Element) -> Self {
        AtomAst(AstAtomAst::from_element(ChemElement::from(&element)))
    }

    #[getter]
    fn element(&self) -> ElementAst {
        ElementAst::from_ast(&self.0.element)
    }

    #[getter]
    fn isotope_mass(&self) -> IsotopeMassAst {
        IsotopeMassAst::from_ast(&self.0.isotope_mass)
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.charge)
    }

    #[getter]
    fn implicit_hydrogens(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.implicit_hydrogens)
    }

    #[getter]
    fn lone_pairs(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.lone_pairs)
    }

    #[getter]
    fn spin(&self, py: Python<'_>) -> PyResult<SpinStateAst> {
        SpinStateAst::from_ast(py, &self.0.spin)
    }

    /// A copy with the element replaced.
    fn with_element(&self, element: PyRef<'_, ElementAst>) -> AtomAst {
        AtomAst(self.0.clone().with_element(element.to_ast()))
    }

    /// A copy with the isotope mass replaced.
    fn with_isotope_mass(&self, isotope_mass: PyRef<'_, IsotopeMassAst>) -> AtomAst {
        AtomAst(self.0.clone().with_isotope_mass(isotope_mass.to_ast()))
    }

    /// A copy with the charge replaced.
    fn with_charge(&self, py: Python<'_>, charge: PyRef<'_, ValueAst>) -> AtomAst {
        AtomAst(self.0.clone().with_charge(charge.to_ast(py)))
    }

    /// A copy with the implicit-hydrogen count replaced.
    fn with_implicit_hydrogens(
        &self,
        py: Python<'_>,
        implicit_hydrogens: PyRef<'_, ValueAst>,
    ) -> AtomAst {
        AtomAst(
            self.0
                .clone()
                .with_implicit_hydrogens(implicit_hydrogens.to_ast(py)),
        )
    }

    /// A copy with the lone-pair count replaced.
    fn with_lone_pairs(&self, py: Python<'_>, lone_pairs: PyRef<'_, ValueAst>) -> AtomAst {
        AtomAst(self.0.clone().with_lone_pairs(lone_pairs.to_ast(py)))
    }

    /// A copy with the spin state replaced.
    fn with_spin(&self, py: Python<'_>, spin: PyRef<'_, SpinStateAst>) -> AtomAst {
        AtomAst(self.0.clone().with_spin(spin.to_ast(py)))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{MemOp as AstMemOp, ValueAst as AstValueAst};

    use super::*;

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

    #[rstest]
    fn test_atom_ast_from_element() {
        let atom = AtomAst::from_element(Element::from(ChemElement::C));
        assert_eq!(atom.0, AstAtomAst::from_element(ChemElement::C));
    }
}
