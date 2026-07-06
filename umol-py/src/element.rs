//! `Element` — a periodic-table element, wrapping `umol_chem::Element`.
// Module complete: blanket-allow the `absolute_paths` false positives from pyo3's
// `hash` derive (hygienic `::std::…` paths). Hand-written code here imports at top.
#![allow(clippy::absolute_paths)]

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_chem::element::Element as ChemElement;

/// A chemical element (periodic-table entry).
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Element(ChemElement);

impl From<ChemElement> for Element {
    fn from(element: ChemElement) -> Self {
        Element(element)
    }
}

impl From<&Element> for ChemElement {
    fn from(element: &Element) -> Self {
        element.0
    }
}

#[pymethods]
impl Element {
    /// Construct from an IUPAC element symbol, e.g. `"C"`, `"Cl"`.
    #[new]
    fn new(symbol: &str) -> PyResult<Self> {
        ChemElement::try_from(symbol)
            .map(Self)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Construct from atomic number (1–118).
    #[staticmethod]
    fn from_atomic_number(number: u8) -> PyResult<Self> {
        ChemElement::from_atomic_number(number)
            .map(Self)
            .ok_or_else(|| PyValueError::new_err(format!("invalid atomic number: {number}")))
    }

    /// IUPAC element symbol.
    #[getter]
    fn symbol(&self) -> &'static str {
        self.0.symbol()
    }

    /// Atomic number (proton count).
    #[getter]
    fn atomic_number(&self) -> u8 {
        self.0.atomic_number()
    }

    fn __repr__(&self) -> String {
        format!("Element('{}')", self.0.symbol())
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case("H", 1)]
    #[case("C", 6)]
    #[case("Cl", 17)]
    #[case("Og", 118)]
    fn test_element_new(#[case] symbol: &str, #[case] atomic_number: u8) {
        let element = Element::new(symbol).unwrap();
        assert_eq!(element.symbol(), symbol);
        assert_eq!(element.atomic_number(), atomic_number);
    }

    #[rstest]
    #[case("X")]
    #[case("")]
    #[case("carbon")]
    fn test_element_new_error(#[case] symbol: &str) {
        assert!(Element::new(symbol).is_err());
    }

    #[rstest]
    #[case(1, "H")]
    #[case(6, "C")]
    #[case(118, "Og")]
    fn test_element_from_atomic_number(#[case] number: u8, #[case] symbol: &str) {
        let element = Element::from_atomic_number(number).unwrap();
        assert_eq!(element.symbol(), symbol);
        assert_eq!(element.atomic_number(), number);
    }

    #[rstest]
    #[case(0)]
    #[case(119)]
    #[case(255)]
    fn test_element_from_atomic_number_error(#[case] number: u8) {
        assert!(Element::from_atomic_number(number).is_err());
    }

    #[rstest]
    fn test_element_repr() {
        assert_eq!(Element::new("C").unwrap().__repr__(), "Element('C')");
    }
}
