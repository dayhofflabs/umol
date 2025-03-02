// GraphAtom implementation

use crate::atom::{AtomSite, Element};
use crate::error::MoleculeError;
use std::fmt::{self, Display};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GraphAtom {
    element: Option<Element>,
    charge: i8,
    unpaired_electrons: u8,
    implicit_hydrogens: u8,
}

impl AtomSite for GraphAtom {
    fn element(&self) -> Option<Element> {
        self.element
    }
}

impl GraphAtom {
    pub fn new<T: Into<GraphAtom>>(value: T) -> Self {
        value.into()
    }

    pub fn with_charge(mut self, charge: i8) -> Result<Self, MoleculeError> {
        if let Some(element) = self.element {
            element.validate_charge(charge)?;
        }
        self.charge = charge;
        Ok(self)
    }

    pub fn with_unpaired_electrons(
        mut self,
        unpaired_electrons: u8,
    ) -> Result<Self, MoleculeError> {
        if let Some(element) = self.element {
            element.validate_unpaired_electrons(unpaired_electrons)?;
        }
        self.unpaired_electrons = unpaired_electrons;
        Ok(self)
    }

    pub fn with_implicit_hydrogens(mut self, implicit_hydrogens: u8) -> Self {
        self.implicit_hydrogens = implicit_hydrogens;
        self
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.unpaired_electrons
    }

    pub fn implicit_hydrogens(&self) -> u8 {
        self.implicit_hydrogens
    }
}

pub struct GraphAtomBuilder {
    element: Option<Element>,
    charge: Option<i8>,
    unpaired_electrons: Option<u8>,
    implicit_hydrogens: Option<u8>,
}

impl GraphAtomBuilder {
    pub fn new(element: Element) -> Self {
        Self {
            element: Some(element),
            charge: None,
            unpaired_electrons: None,
            implicit_hydrogens: None,
        }
    }

    pub fn charge(mut self, charge: i8) -> Self {
        self.charge = Some(charge);
        self
    }

    pub fn unpaired_electrons(mut self, unpaired_electrons: u8) -> Self {
        self.unpaired_electrons = Some(unpaired_electrons);
        self
    }

    pub fn implicit_hydrogens(mut self, implicit_hydrogens: u8) -> Self {
        self.implicit_hydrogens = Some(implicit_hydrogens);
        self
    }

    pub fn build(self) -> Result<GraphAtom, MoleculeError> {
        let element = self.element.unwrap();
        if let Some(charge) = self.charge {
            element.validate_charge(charge)?;
        }
        if let Some(unpaired) = self.unpaired_electrons {
            element.validate_unpaired_electrons(unpaired)?;
        }

        Ok(GraphAtom {
            element: self.element,
            charge: self.charge.unwrap_or(0),
            unpaired_electrons: self.unpaired_electrons.unwrap_or(0),
            implicit_hydrogens: self.implicit_hydrogens.unwrap_or(0),
        })
    }
}

impl From<Element> for GraphAtom {
    fn from(element: Element) -> Self {
        GraphAtom {
            element: Some(element),
            charge: 0,
            unpaired_electrons: 0,
            implicit_hydrogens: 0,
        }
    }
}

impl TryFrom<&str> for GraphAtom {
    type Error = MoleculeError;

    fn try_from(symbol: &str) -> Result<Self, Self::Error> {
        match Element::from_symbol(symbol) {
            Some(element) => Ok(element.into()),
            None => Err(MoleculeError::InvalidElementSymbol(symbol.to_string())),
        }
    }
}

impl Display for GraphAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let element = self.element.unwrap();
        let symbol = element.symbol();

        // Always use brackets if:
        // - Symbol is 2 characters
        // - Has charge
        // - Has unpaired electrons
        let needs_brackets = symbol.len() == 2 || self.charge != 0 || self.unpaired_electrons > 0;

        if needs_brackets {
            write!(f, "[{}", symbol)?;

            if self.charge > 0 {
                write!(f, "+{}", self.charge)?;
            } else if self.charge < 0 {
                write!(f, "{}", self.charge)?; // Negative sign included automatically
            }

            if self.unpaired_electrons > 0 {
                write!(f, "^{}", self.unpaired_electrons)?;
            }

            write!(f, "]")
        } else {
            write!(f, "{}", symbol)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Element;
    use rstest::rstest;

    #[test]
    fn test_atom_new() {
        let atom = GraphAtom::new(Element::C);
        assert_eq!(atom.element(), Some(Element::C));
        assert_eq!(atom.charge(), 0);
        assert_eq!(atom.unpaired_electrons(), 0);
        assert_eq!(atom.implicit_hydrogens(), 0);
    }

    #[test]
    fn test_atom_with_charge() {
        let atom = GraphAtom::new(Element::C);
        let atom = atom.with_charge(1).unwrap();
        assert_eq!(atom.charge(), 1);

        let atom = atom.with_charge(-1).unwrap();
        assert_eq!(atom.charge(), -1);

        // Test validation
        let atom = GraphAtom::new(Element::H);
        assert!(atom.with_charge(2).is_err());
        assert!(atom.with_charge(-2).is_err());
    }

    #[test]
    fn test_atom_with_unpaired_electrons() {
        let atom = GraphAtom::new(Element::C);
        let atom = atom.with_unpaired_electrons(2).unwrap();
        assert_eq!(atom.unpaired_electrons(), 2);

        // Test validation
        let atom = GraphAtom::new(Element::H);
        assert!(atom.with_unpaired_electrons(2).is_err());
        assert!(atom.with_unpaired_electrons(0).is_ok());
        assert!(atom.with_unpaired_electrons(1).is_ok());
    }

    #[test]
    fn test_atom_with_implicit_hydrogens() {
        let atom = GraphAtom::new(Element::C);
        let atom = atom.with_implicit_hydrogens(4);
        assert_eq!(atom.implicit_hydrogens(), 4);
    }

    #[test]
    fn test_atom_builder() {
        let atom = GraphAtomBuilder::new(Element::N)
            .charge(-1)
            .unpaired_electrons(2)
            .implicit_hydrogens(1)
            .build()
            .unwrap();

        assert_eq!(atom.element(), Some(Element::N));
        assert_eq!(atom.charge(), -1);
        assert_eq!(atom.unpaired_electrons(), 2);
        assert_eq!(atom.implicit_hydrogens(), 1);
    }

    #[test]
    fn test_atom_builder_validation() {
        // Invalid charge
        let result = GraphAtomBuilder::new(Element::H)
            .charge(2)
            .build();
        assert!(result.is_err());

        // Invalid unpaired electrons
        let result = GraphAtomBuilder::new(Element::H)
            .unpaired_electrons(2)
            .build();
        assert!(result.is_err());

        // Valid combination
        let result = GraphAtomBuilder::new(Element::C)
            .charge(-1)
            .unpaired_electrons(3)
            .build();
        assert!(result.is_ok());
    }

    #[test]
    fn test_atom_from_element() {
        let atom: GraphAtom = Element::O.into();
        assert_eq!(atom.element(), Some(Element::O));
        assert_eq!(atom.charge(), 0);
        assert_eq!(atom.unpaired_electrons(), 0);
    }

    #[test]
    fn test_atom_try_from_str() {
        let atom = GraphAtom::try_from("C").unwrap();
        assert_eq!(atom.element(), Some(Element::C));

        let atom = GraphAtom::try_from("N").unwrap();
        assert_eq!(atom.element(), Some(Element::N));

        assert!(GraphAtom::try_from("Xx").is_err());
    }

    #[rstest]
    #[case(Element::C, 0, 0, "C")]
    #[case(Element::C, 1, 0, "[C+1]")]
    #[case(Element::C, -1, 0, "[C-1]")]
    #[case(Element::C, 0, 2, "[C^2]")]
    #[case(Element::C, 1, 2, "[C+1^2]")]
    #[case(Element::He, 0, 0, "[He]")]  // Two-letter element always gets brackets
    fn test_atom_display(
        #[case] element: Element,
        #[case] charge: i8,
        #[case] unpaired: u8,
        #[case] expected: &str,
    ) {
        let mut atom = GraphAtom::new(element);
        if charge != 0 {
            atom = atom.with_charge(charge).unwrap();
        }
        if unpaired != 0 {
            atom = atom.with_unpaired_electrons(unpaired).unwrap();
        }
        assert_eq!(format!("{}", atom), expected);
    }

    #[test]
    fn test_atom_site_trait() {
        let atom = GraphAtom::new(Element::F);
        assert_eq!(atom.element(), Some(Element::F));
    }
}
