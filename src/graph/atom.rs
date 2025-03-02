// GraphAtom implementation

use crate::atom::{AtomSite, Element};
use crate::error::MoleculeError;
use std::fmt::{self, Display};

#[derive(Debug, Clone)]
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
