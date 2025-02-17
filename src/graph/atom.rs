// GraphAtom implementation

use crate::atom::{AtomSite, Element};
use crate::error::MoleculeError;

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

    pub fn with_charge(mut self, charge: i8) -> Self {
        self.charge = charge;
        self
    }

    pub fn with_unpaired_electrons(mut self, unpaired_electrons: u8) -> Self {
        self.unpaired_electrons = unpaired_electrons;
        self
    }

    pub fn with_implicit_hydrogens(mut self, implicit_hydrogens: u8) -> Self {
        self.implicit_hydrogens = implicit_hydrogens;
        self
    }
}

impl From<Element> for GraphAtom {
    fn from(element: Element) -> Self {
        Self {
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
            Some(element) => Ok(GraphAtom::new(element)),
            None => Err(MoleculeError::InvalidElementSymbol(symbol.to_string())),
        }
    }
}
